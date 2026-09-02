//! COUP-3 §3.1 — the **one deterministic SDC-IMEX step**, plus the COUP-2
//! machinery it arms every step: the wall-exchange (Robin-Robin) coupling
//! (COUP-2 §3.5) and the flux-telescoping conservation audit (COUP-2 §3.1).
//! This module retires every declared-scaffolding integrator (session-5
//! explicit conduction, session-7 MOL SSP-RK2, session-10/11 flux-matched
//! coupled splitting): the production advance is THIS step, everywhere.
//!
//! ## The step (2 Lobatto nodes, fixed sweeps — Minion IMEX-SDC)
//!
//! Write the system `U_t = A(U) + D(U)`: `A` = the explicit hyperbolic
//! class (the Euler operator's `eval_rhs`), `D` = the spatially-coupled
//! implicit diffusion class (the conduction operator, solved per sweep by
//! a deterministic fixed-structure preconditioned-CG solve — §3.1's
//! "fixed-cycle multigrid or CG" with CG as the first occupant; a
//! geometric-multigrid occupant may supersede it behind the same seam).
//! With node states `U₀ = Uⁿ` and `U₁ ≈ Uⁿ⁺¹`:
//!
//! - **Predictor** (IMEX Euler): `U₁⁽⁰⁾ = U₀ + Δt·A(U₀) + Δt·D(U₁⁽⁰⁾)`.
//! - **Correction sweep k** (trapezoid quadrature):
//!   `U₁⁽ᵏ⁾ = U₀ + Δt/2·[A(U₀) + A(U₁⁽ᵏ⁻¹⁾)] + Δt/2·D(U₀) − Δt/2·D(U₁⁽ᵏ⁻¹⁾)
//!    + Δt·D(U₁⁽ᵏ⁾)`, for a **fixed** [`N_SDC_CORRECTIONS`] sweeps.
//!
//! The fixed point is the trapezoidal rule — 2nd order in both classes;
//! the truncated sweep count keeps |R| < 1 for stiff `D` (the sweeps damp
//! toward trapezoid from the L-stable BE predictor) and gives the explicit
//! part the stability polynomial `1 + z + z²/2 + z³/4` (imaginary-axis
//! stable to |z| ≤ 2 — the fluid CFL alone sets Δt; class `D` has **no**
//! step limit, which is the point: fine wall cells stop binding Δt).
//! On cut-geometry worlds every node-state composition ends with the
//! State-Redistribution pass (Berger–Giuliani apply SRD stagewise).
//!
//! ## Robin-Robin wall exchange (COUP-2 §3.5, placed here)
//!
//! Gas↔solid heat exchange enters each sweep's class-`D` solve as **Robin
//! interface conditions with the wall-function h as the Robin coefficient**
//! (linear in the solid T — unconditionally stable at Biot > 1), iterated
//! by a fixed count of Picard sweeps ([`N_ROBIN_SWEEPS`]) with clamped
//! Aitken relaxation over the gas-operand refresh. The gas debits exactly
//! the per-face heats the accepted solid solve received (one assembly is
//! the single owner of those numbers — conservation by construction), and
//! the post-sweep converged wall-heat integral is the jacket enthalpy rise
//! COUP-3 §3.5's expander solve consumes.
//!
//! ## The every-step audit (COUP-2 §3.1)
//!
//! `Δ(stored) = Σ(port fluxes) + Σ(sources)` per conserved quantity, where
//! `Σ(sources)` is the **integrator-applied increments** (the exact
//! SDC-node-weighted composition below — never a rate×Δt recomputation):
//! interior fluxes telescope, so only run-boundary ports (domain BCs, wall
//! faces incl. declared stair transpiration), applied volumetric sources
//! (geometric/closure/external), the exchange pair (cancelling), and the
//! solid's non-interior heats are ledgered. Tolerance = the derived
//! `TOL_AUDIT[q]` of §3.1.1 (`K_AUDIT·ε·√N·S[q]` + declared floor). A
//! violation is a **halt with diagnosis** (META-1 P6), never a warning.
//!
//! Named constants below are manifest-recording candidates; until FND-6's
//! results-bundle wave lands they are compile-time constants pinned by the
//! recorded build fingerprint (META-1 §2.2) and printed in certificates.
//!
//! Determinism (COUP-3 §3.7): fixed sweep counts, fixed CG iteration
//! structure (absolute tolerance + fixed max-iters, per META-1 §2.2),
//! fixed clamped Aitken relaxation, fixed-order reductions everywhere —
//! bit-identical at any thread count (asserted by the S2 battery).

use std::collections::BTreeMap;

use crate::conduction::{
    AssembleMode, Conduction, Domain, ExchangeKey, GasFaceRobin, HeatLedger, SolverError,
};
use crate::euler::{
    Combustion, Cons, EosLaw, Euler, EulerFields, EulerWorkspace, FlowError, FlowLedger, I_EN,
    I_MR, I_MT, I_MZ, I_RB, I_RC, NCOMP, Prim, srd_neighborhood,
};
use crate::gas_diffusion::{
    GasComp, GasDiffError, GasDiffusion, GasOperands, GasTransportField, GasWork,
};
use crate::transport::TransportProps;
use crate::wall_heat::{NearWallGas, WallHeatError, WallLaw};
use crucible_grid::{BRICK, BRICK_CELLS, FaceDir, FieldId, Grid, InterfaceFace, tree_combine};

// --- Named constants (COUP-3 §3.7 / COUP-2 §3.1.1) ---------------------------

/// Fixed SDC correction-sweep count after the IMEX-Euler predictor
/// (COUP-3 §3.1: "2–3 sweeps to 2nd order" — predictor + 2 corrections).
/// The audit's final-composition weights assume the last sweep is a
/// trapezoid correction — compile-time-guarded below.
pub const N_SDC_CORRECTIONS: usize = 2;
const _: () = assert!(N_SDC_CORRECTIONS >= 1);

/// Fixed Picard sweep count of the Robin-Robin wall-exchange solve inside
/// each SDC sweep's class-`D` solve (COUP-2 §3.5).
///
/// **Raised 3 → 5 at S4 (review finding).** Before the FND-7 spine, the
/// wall law's `k` and `c_p` were config constants, so `h` did not depend on
/// the near-wall gas temperature and the map these sweeps relax,
/// `q(T_gas) = h·(T_aw(T_gas) − T_w)`, was **affine** — three sweeps
/// converged it to round-off. With per-cell transport `h` is a function of
/// the operand the sweeps move, and the map is nonlinear; the contraction
/// survives (it is still geometric) but the *margin* on
/// [`EPS_ROBIN_RESID`] does not. Measured on a four-class duct with a
/// near-wall gradient of ~1.7e5 K/m — milder than an RL10 chamber wall:
///
/// | sweeps | 3 | 4 | 5 | 6 |
/// |---|---|---|---|---|
/// | residual | 3.8e-6 (**halt**) | 5.5e-10 | 1.6e-11 | 1.2e-13 |
///
/// Four would clear the acceptance; five is chosen because the quantity
/// that sets the required count — `dh/dT_gas` through the equilibrium
/// conductivity — grows with dissociation, and the shipped surface reaches
/// `k_eff/k_frozen` = 22.7 at its hottest, thinnest corner. The extra sweep
/// costs a class-`D` solve per SDC sweep; a halt on a legal config costs
/// the run.
pub const N_ROBIN_SWEEPS: usize = 5;

/// Residual acceptance of the Robin-Robin solve: max relative change of
/// any patch's exchange heat between the last two Picard sweeps. Failure ⇒
/// `COUPLING_RESIDUAL` (numerical — a solver defect, never a verdict).
/// Sized on the COUP-3 §3.5 principle (orders below the wall law's ±20–30%
/// band, so the floor never contributes to the physics error) — NOT at
/// machine zero: the residual measures only the gas-operand staleness of
/// the last sweep (the exchange pair itself is exactly conservative at any
/// residual — both sides use the accepted assembly's numbers), and during
/// violent start transients the operands legitimately move ~1e-8/sweep
/// (measured, RL10 smoke march).
///
/// **S4 note:** that ~1e-8/sweep was measured when `h` was independent of
/// the operands. It no longer is (see [`N_ROBIN_SWEEPS`]), so the sweep
/// count — not this constant — carries the margin.
///
/// **S7 note (1e-6 → 1e-4):** the startup march's real pre-spark state is
/// a NEAR-VACUUM cold fill (~10² Pa), where the wall-adjacent gas cell's
/// thermal mass is ~10⁴× smaller than at the stations' dense fills — the
/// exchange map's contraction factor genuinely weakens (the operand swings
/// per debited joule grow as 1/ρ) and the fixed five sweeps land at a
/// measured ~3e-5 relative on a ~10 W exchange (0.3 mW of staleness —
/// physically nothing). The acceptance is a HALT GATE, not a solution
/// modifier: relaxing it changes no accepted number anywhere (the stations'
/// residuals remain ≪ 1e-6 and their certificates byte-identical); it only
/// stops a legal near-vacuum start from being declared a solver defect.
/// 1e-4 still sits 3+ orders below the wall law's ±20–30% band — the
/// COUP-3 §3.5 sizing principle — and a genuinely non-contracting map
/// (O(1) residual) is still caught.
pub const EPS_ROBIN_RESID: f64 = 1e-4;

/// Clamped Aitken relaxation bounds of the Picard sweeps (COUP-3 §3.5's
/// deterministic clamp discipline, applied to the exchange iteration).
pub const ROBIN_OMEGA_MIN: f64 = 0.1;
pub const ROBIN_OMEGA_MAX: f64 = 2.0;

/// Residual acceptance of the gas class-`D` cross-term Picard (S3): the
/// STATE effect over this step of the rate change between the last two
/// Picard sweeps, relative to the conserved components' magnitudes (see
/// `rate_resid`). This is a **contraction guard, not an
/// accuracy floor**: the truncated Picard is the same kind of fixed-count
/// iteration as the SDC sweeps themselves — its remainder is a temporal-
/// truncation term of the integrator (verified by the dt-Richardson and
/// MMS order gates), legitimately ~1e-2–1e-4 relative during violent
/// transients and ~roundoff near steady state. What must NEVER happen is
/// non-contraction: the lagged remainder's structural gain is ≲ 1/12 at
/// any Δt (AM-GM over the implicit diagonals), so a final-sweep relative
/// change above this ceiling means the cross terms are not contracting —
/// a broken assembly or a regime outside the bound. Failure ⇒
/// `COUPLING_RESIDUAL` (numerical — a solver defect, never a verdict).
///
/// The practical margin assumes [`N_ROBIN_SWEEPS`] ≥ 2: on an impulsive
/// start from rest the composed-state momentum scale is O(Δt·rate), so the
/// FIRST iterate's self-change is O(1) relative. The check measures the
/// LAST iterate (post-contraction), which is why it passes there; dropping
/// the Picard count to 1 would make cold starts trip this spuriously —
/// loudly, never silently (S3 review finding).
pub const EPS_GAS_DIFF_RESID: f64 = 0.25;

/// Fixed iteration cap of the class-`D` preconditioned-CG solve. The stop
/// is the absolute+deterministic rule of META-1 §2.2: residual below
/// [`EPS_CG_RESID`]·‖b‖ or the cap, whichever first — a pure function of
/// the data. The residual-acceptance check after the loop raises
/// `COUPLING_RESIDUAL` if the cap was hit unconverged.
pub const N_CG_ITERS_MAX: usize = 512;

/// Class-`D` CG residual acceptance, relative to the right-hand side norm.
/// Sized so the solve's conservation leakage sits below `TOL_AUDIT` (the
/// audit closes THROUGH the implicit solve, not around it).
pub const EPS_CG_RESID: f64 = 1e-12;

/// COUP-2 §3.1.1 audit safety factor (headroom for flux-aggregation and
/// ledger arithmetic beyond the bare `ε·√N` reduction bound).
pub const K_AUDIT: f64 = 100.0;

// --- Wall patches (SOLV-1 §3.5 exchange surface; single owner) ---------------

/// One gas-cell wall patch — the SOLV-1 §3.5 exchange surface. On
/// cut-geometry worlds a patch is a gas cell's embedded interface: area
/// |W| and normal from the grid's closure identity, spanning its
/// gas↔solid grid faces (the flux carriers), with the SRD-neighborhood
/// debit set (a sliver cell cannot absorb its own wall debit). On box/
/// stair worlds a patch is one interface face with its own area and
/// axis-aligned normal — the same law, the geometry decides the shape
/// (data, not a code branch).
pub struct WallPatch {
    pub gas: (usize, usize),
    pub faces: Vec<InterfaceFace>,
    /// Full grid area of each spanned face (per full ring, N_θ = 1).
    pub face_areas: Vec<f64>,
    /// Exchange area: |W| (cut) or the face area (box).
    pub area: f64,
    /// `area / Σ face_areas` — maps stair-face carriers onto the exchange
    /// area (exactly 1.0 on box worlds).
    pub area_scale: f64,
    /// Outward (gas→wall) unit normal (n_r, n_z).
    pub n_hat: (f64, f64),
    /// Energy-debit cells (the gas cell's SRD neighborhood; `[self]` for
    /// regular cells) and Σ κV over them.
    pub debit_cells: Vec<(usize, usize)>,
    pub debit_kv_sum: f64,
}

fn opposite(d: FaceDir) -> FaceDir {
    match d {
        FaceDir::RMinus => FaceDir::RPlus,
        FaceDir::RPlus => FaceDir::RMinus,
        FaceDir::ZMinus => FaceDir::ZPlus,
        FaceDir::ZPlus => FaceDir::ZMinus,
        // S9 compile completeness: wall faces come from `gas_solid_faces`,
        // which enumerates only r/z faces (regions are (r,z)-shaped;
        // per-θ wall patches ride plan S11).
        FaceDir::ThetaMinus => FaceDir::ThetaPlus,
        FaceDir::ThetaPlus => FaceDir::ThetaMinus,
    }
}

/// The [`ExchangeKey`] of a patch face — the solid cell and the face
/// direction as seen from it.
fn face_key(face: &InterfaceFace) -> ExchangeKey {
    (face.solid.0, face.solid.1, opposite(face.dir).index() as u8)
}

/// Enumerate wall patches in deterministic order (gas-cell lexicographic,
/// then face direction). See [`WallPatch`] for the two world classes. The
/// single-valued contour class cannot produce a multi-sided wall (slot)
/// inside one cell; per-side reconstruction for genuine slots is the FND-3
/// PLIC/CSG wave.
pub fn build_wall_patches(g: &Grid) -> Result<Vec<WallPatch>, String> {
    if g.bricks().iter().any(|b| b.n_theta() != 1) {
        // The patch areas, closure vectors, debit κV sums, AND the
        // exchange-heat keying are per-full-ring (N_θ = 1) here; per-θ
        // patches arrive with the coarse-3-D RL10 wave (plan S11) — refuse rather than
        // undercount (S2 review finding).
        return Err(
            "wall patches at N_θ > 1 arrive with the coarse-3-D RL10 wave (plan S11); \
             refusing rather than guessing"
                .to_string(),
        );
    }
    let cut = g.has_cut_geometry();
    let mut patches: Vec<WallPatch> = Vec::new();
    for face in g.gas_solid_faces() {
        let start_new = if cut {
            patches.last().map(|p| p.gas) != Some(face.gas)
        } else {
            true // box worlds: one patch per face
        };
        if start_new {
            patches.push(WallPatch {
                gas: face.gas,
                faces: Vec::new(),
                face_areas: Vec::new(),
                area: 0.0,
                area_scale: 1.0,
                n_hat: (0.0, 0.0),
                debit_cells: Vec::new(),
                debit_kv_sum: 0.0,
            });
        }
        let p = patches.last_mut().expect("just pushed");
        p.face_areas.push(g.interface_area_per_theta(&face, 1));
        p.faces.push(face);
    }
    for p in &mut patches {
        let (i_r, i_z) = p.gas;
        if cut {
            let (w_r, w_z) = g.wall_closure(i_r, i_z, 1);
            let area = (w_r * w_r + w_z * w_z).sqrt();
            if !area.is_finite() || area <= 0.0 {
                return Err(format!(
                    "wall cell ({i_r}, {i_z}): zero closure interface area yet gas↔solid \
                     faces exist — geometry incoherent (a slot-class wall? per-side \
                     interface reconstruction is the FND-3 PLIC wave); refusing"
                ));
            }
            p.area = area;
            p.area_scale = area / p.face_areas.iter().sum::<f64>();
            p.n_hat = (-w_r / area, -w_z / area);
        } else {
            p.area = p.face_areas[0];
            p.area_scale = 1.0;
            p.n_hat = match p.faces[0].dir {
                FaceDir::RMinus => (-1.0, 0.0),
                FaceDir::RPlus => (1.0, 0.0),
                FaceDir::ZMinus => (0.0, -1.0),
                FaceDir::ZPlus => (0.0, 1.0),
                FaceDir::ThetaMinus | FaceDir::ThetaPlus => unreachable!(
                    "gas_solid_faces enumerates only r/z faces — regions are \
                     (r,z)-shaped (per-θ wall patches ride plan S11)"
                ),
            };
        }
        match srd_neighborhood(g, i_r, i_z).map_err(|e| format!("wall patch: {e}"))? {
            Some(members) => {
                p.debit_kv_sum = members.iter().map(|(_, kv)| kv).sum();
                p.debit_cells = members.into_iter().map(|(c, _)| c).collect();
            }
            None => {
                p.debit_cells = vec![p.gas];
                p.debit_kv_sum = g.kappa(i_r, i_z) * g.cell_volume(i_r, 1);
            }
        }
    }
    Ok(patches)
}

// --- Operator classes (the per-step schedule arguments) ----------------------

/// Class `A` — the explicit hyperbolic operator and its fields.
pub struct FlowClass<'a, 'b, E: EosLaw> {
    pub op: &'a Euler<'b, E>,
    pub fields: &'a EulerFields,
}

/// Class `D` — the spatially-coupled implicit diffusion operator: the
/// conduction spatial discretization, its temperature field, and a scratch
/// field for the CG direction vector (the retired explicit-rate slot).
pub struct DiffusionClass<'a, 'b> {
    pub op: &'a Conduction<'b>,
    pub t_field: FieldId,
    pub scratch_field: FieldId,
}

/// The Robin-Robin wall-exchange coupling (COUP-2 §3.5): the config-time
/// patch list, the one wall-function law, and the FND-7 temperature query
/// for the gas operands (an explicit closure until the spine's transport
/// wave folds temperature into the `EosLaw` trait — OFFL-5).
pub struct ExchangeClass<'a> {
    pub patches: &'a [WallPatch],
    pub law: &'a WallLaw,
    pub temperature: &'a (dyn Fn(&Prim) -> Result<f64, &'static str> + Sync),
    /// The FND-7 §3.3 spine query at a cell's own state — the wall law's
    /// transport operands (SOLV-1 §3.5, 0.4.1). The **same** closure the
    /// gas class-`D` occupant reads, so the law and the resolved `F_visc`
    /// beside it can never disagree about the medium.
    pub transport: &'a (dyn Fn(&Prim) -> Result<TransportProps, &'static str> + Sync),
}

/// The gas-phase class-`D` occupant (S3): `F_visc` — compressible viscous
/// stress + Fourier conduction + species diffusion (`crate::gas_diffusion`),
/// solved inside the same fixed Picard sweeps as the wall exchange. The
/// temperature query is the same FND-7 seam closure as [`ExchangeClass`].
pub struct GasDiffusionClass<'a> {
    pub op: &'a GasDiffusion<'a>,
    pub temperature: &'a (dyn Fn(&Prim) -> Result<f64, &'static str> + Sync),
    /// The FND-7 §3.3 spine query (see [`ExchangeClass::transport`]).
    pub transport: &'a (dyn Fn(&Prim) -> Result<TransportProps, &'static str> + Sync),
}

/// Class `R` — the cell-local implicit stiff-reaction occupant (COUP-3
/// §3.3, S7): the SOLV-4 §3.6 auto-ignition term, advanced per sweep by the
/// operator's own fixed-structure implicit node solve
/// ([`Combustion::implicit_auto_update`]) with the trapezoid quadrature
/// carried on **realized** rates (COUP-2: applied increments, never
/// rate×Δt). The propagation/diffusion limbs of SOLV-4.4 stay in class `A`
/// (they ride `Euler::eval_rhs`); this class holds only the term whose
/// `τ_ign` collapses below the acoustic Δt at chamber conditions. The
/// `Combustion` operator here must be the SAME one the flow class's `Euler`
/// carries (one blend, one closure surface — the caller wires both from one
/// build; a mismatch would be a config defect, not a runtime branch).
pub struct ReactionClass<'a> {
    pub op: &'a Combustion<'a>,
}

/// COUP-2 §3.1.1 — the audit's declared reference scales (the absolute
/// floor ingredient `TOL_AUDIT_FLOOR[q] = K_AUDIT·ε·√N·ref[q]`). Zero is
/// legal: the throughput term of `S[q]` already scales every quantity that
/// moved at all; a floor guards quantities that are identically zero AND
/// unmoved (then both sides of the identity are zero too).
#[derive(Debug, Clone, Copy)]
pub struct AuditSpec {
    pub k_audit: f64,
    /// Reference scales: [mass, momentum, energy, composition].
    pub ref_scale: [f64; 4],
}

impl Default for AuditSpec {
    fn default() -> Self {
        AuditSpec {
            k_audit: K_AUDIT,
            ref_scale: [0.0; 4],
        }
    }
}

// --- Errors & reports ---------------------------------------------------------

#[derive(Debug)]
pub enum SdcError {
    Flow(FlowError),
    Solid(SolverError),
    Wall(WallHeatError),
    /// The gas class-`D` operator refused (bad coefficient, non-finite
    /// operand, N_θ > 1) — `crate::gas_diffusion`.
    Gas(GasDiffError),
    /// A fixed-sweep coupling solve missed its named residual acceptance —
    /// COUP-4's `COUPLING_RESIDUAL` halt class (numerical: a solver defect,
    /// never an engine verdict).
    CouplingResidual {
        solve: &'static str,
        resid: f64,
        eps: f64,
    },
    /// COUP-2 §2: the conservation identity failed beyond `TOL_AUDIT` —
    /// halt with diagnosis (quantity, magnitudes), never a tolerance to
    /// live with.
    AuditViolation {
        quantity: &'static str,
        delta: f64,
        applied: f64,
        tol: f64,
    },
    /// A step-schedule contradiction (e.g. exchange without both classes).
    Config(&'static str),
}

impl std::fmt::Display for SdcError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Flow(e) => write!(f, "hyperbolic class: {e}"),
            Self::Solid(e) => write!(f, "diffusion class: {e}"),
            Self::Wall(e) => write!(f, "wall law: {e}"),
            Self::Gas(e) => write!(f, "gas diffusion class: {e}"),
            Self::CouplingResidual { solve, resid, eps } => write!(
                f,
                "COUPLING_RESIDUAL: {solve} residual {resid:.3e} exceeds its acceptance \
                 {eps:.1e} after the fixed sweeps (numerical — a solver defect, never an \
                 engine verdict; COUP-3 §2)"
            ),
            Self::AuditViolation {
                quantity,
                delta,
                applied,
                tol,
            } => write!(
                f,
                "CONSERVATION AUDIT VIOLATION ({quantity}): Δstored = {delta:.6e} but the \
                 integrator applied {applied:.6e} (|gap| > TOL_AUDIT = {tol:.3e}) — an \
                 increment escaped the ledger; halt with diagnosis (COUP-2 §2, META-1 P6)"
            ),
            Self::Config(m) => write!(f, "SDC step configuration: {m}"),
        }
    }
}

impl std::error::Error for SdcError {}

impl From<FlowError> for SdcError {
    fn from(e: FlowError) -> Self {
        Self::Flow(e)
    }
}
impl From<SolverError> for SdcError {
    fn from(e: SolverError) -> Self {
        Self::Solid(e)
    }
}
impl From<WallHeatError> for SdcError {
    fn from(e: WallHeatError) -> Self {
        Self::Wall(e)
    }
}

impl From<GasDiffError> for SdcError {
    fn from(e: GasDiffError) -> Self {
        // A missed CG acceptance is the COUP-4 `COUPLING_RESIDUAL` class,
        // reported uniformly with the other fixed-sweep solves.
        if let GasDiffError::CgUnconverged { resid, .. } = e {
            return Self::CouplingResidual {
                solve: "class-D gas-diffusion CG (COUP-3 §3.1, S3)",
                resid,
                eps: EPS_CG_RESID,
            };
        }
        Self::Gas(e)
    }
}

/// One audited quantity's per-step closure record.
#[derive(Debug, Clone, Copy)]
pub struct AuditRow {
    pub quantity: &'static str,
    pub delta: f64,
    pub applied: f64,
    pub tol: f64,
}

/// One accepted step's report.
#[derive(Debug, Clone, Default)]
pub struct StepReport {
    /// Per-quantity audit closure (empty only if nothing was audited).
    pub audit: Vec<AuditRow>,
    /// Exchange readout at the accepted solve (None ⇔ no exchange class).
    pub exchange: Option<ExchangeStepReport>,
    /// Last SOLID class-`D` CG solve's iteration count and residual.
    pub cg_iters: usize,
    pub cg_resid: f64,
    /// Gas class-`D` readout (S3): worst per-component CG iterations/
    /// residual of the accepted solve, and the accepted cross-term Picard
    /// residual (≤ [`EPS_GAS_DIFF_RESID`] by acceptance). Zero ⇔ no gas
    /// diffusion scheduled.
    pub gas_cg_iters: usize,
    pub gas_cg_resid: f64,
    pub gas_picard_resid: f64,
}

/// The Robin-Robin exchange readout of one step.
#[derive(Debug, Clone, Default)]
pub struct ExchangeStepReport {
    /// Per-patch exchange heat [W] at the accepted solve (gas→solid > 0).
    pub q_w: Vec<f64>,
    /// Per-patch film coefficient and recovery temperature (final operands).
    pub h: Vec<f64>,
    pub t_aw: Vec<f64>,
    /// Σ q_w — the jacket enthalpy-rise integrand COUP-3 §3.5 consumes.
    pub jacket_w: f64,
    /// The step's composed exchange energy INTO the solid [J] (the gas was
    /// debited exactly −this), composed exterior-face heat [J], and
    /// composed domain-edge BC heat [J] — a duct-class liner that reaches
    /// the domain edge takes its coolant Robin THERE (S2 review finding:
    /// a coolant ledger reading only the exterior line records zero).
    pub applied_exchange_j: f64,
    pub applied_exterior_j: f64,
    pub applied_bc_j: f64,
    /// Final Picard residual (≤ [`EPS_ROBIN_RESID`] by acceptance).
    pub robin_resid: f64,
}

// --- The integrator -----------------------------------------------------------

/// Per-brick f64 buffers shaped like grid fields.
type BufF = Vec<Vec<f64>>;

struct SolidBufs {
    t0: BufF,
    heat0: BufF,
    heat_prev: BufF,
    heat_cur: BufF,
    apply: BufF,
    diag: BufF,
    b: BufF,
    r: BufF,
    z: BufF,
    q: BufF,
    delta: BufF,
}

/// Gas class-`D` state (S3): the solved/lagged operand sets, the operator
/// scratch, and the SDC-node rate records (the solid's heat0/heat_prev/
/// heat_cur pattern, five components wide; `dlag` is the Picard trial —
/// the q_trial analogue).
struct GasBufs {
    sol: GasOperands,
    lag: GasOperands,
    work: GasWork,
    /// Per-cell spine transport, refreshed with the operands each Picard
    /// iterate (COUP-3 §3.1, 0.4.3).
    tr: GasTransportField,
    d0: Vec<Vec<Cons>>,
    dprev: Vec<Vec<Cons>>,
    dcur: Vec<Vec<Cons>>,
    dlag: Vec<Vec<Cons>>,
    dstage: Vec<Vec<Cons>>,
}

/// Class-`R` node-rate records (S7): the auto-ignition term's node-0 rate,
/// the previous sweep's realized rate, this sweep's realized rate, and a
/// staging buffer for the two-pass (read-solve, then write) application.
struct ReactionBufs {
    r0: BufF,
    r_prev: BufF,
    r_trial: BufF,
    x_new: BufF,
}

/// The one deterministic integrator (persistent workspaces; one instance
/// per grid). Class arguments arrive per step so callers may evolve
/// boundary schedules between steps (COUP-7 declared schedules).
pub struct Sdc {
    ws: Option<EulerWorkspace>,
    rate_e0: Vec<Vec<Cons>>,
    sb: Option<SolidBufs>,
    gb: Option<GasBufs>,
    rb: Option<ReactionBufs>,
    pub audit_spec: AuditSpec,
}

impl Default for Sdc {
    fn default() -> Self {
        Self::new()
    }
}

impl Sdc {
    pub fn new() -> Self {
        Sdc {
            ws: None,
            rate_e0: Vec::new(),
            sb: None,
            gb: None,
            rb: None,
            audit_spec: AuditSpec::default(),
        }
    }

    pub fn with_audit(audit: AuditSpec) -> Self {
        Sdc {
            audit_spec: audit,
            ..Self::new()
        }
    }

    /// The COUP-3 Δt rule: fixed CFL over the hyperbolic wave speeds (the
    /// stiff classes need no CFL — implicit). Warm-started through the
    /// persistent workspace.
    pub fn stable_dt<E: EosLaw + Sync>(
        &mut self,
        g: &Grid,
        flow: &FlowClass<'_, '_, E>,
        cfl: f64,
    ) -> Result<f64, SdcError> {
        let ws = self.ensure_ws(g, flow)?;
        Ok(flow.op.stable_dt_ws(g, flow.fields, ws, cfl)?)
    }

    fn ensure_ws<E: EosLaw + Sync>(
        &mut self,
        g: &Grid,
        flow: &FlowClass<'_, '_, E>,
    ) -> Result<&mut EulerWorkspace, SdcError> {
        // Staleness guard (S2 review): a θ-refined or swapped grid must
        // rebuild the flow workspaces exactly as `ensure_sb` rebuilds the
        // solid ones — a stale `u0`/`rate_e0` would silently mis-snapshot.
        // Per-brick plane shapes (S8): every brick's own N_θ is compared,
        // so a single brick's θ-collapse/refine rebuilds.
        let plane = |bi: usize| g.brick(bi).n_theta() as usize * BRICK_CELLS;
        let stale = self.rate_e0.len() != g.n_bricks()
            || self
                .rate_e0
                .iter()
                .enumerate()
                .any(|(bi, v)| v.len() != plane(bi));
        if self.ws.is_none() || stale {
            self.ws = Some(flow.op.workspace(g)?);
            self.rate_e0 = (0..g.n_bricks())
                .map(|bi| vec![[0.0; NCOMP]; plane(bi)])
                .collect();
        }
        Ok(self.ws.as_mut().expect("just ensured"))
    }

    fn ensure_sb(&mut self, g: &Grid) -> &mut SolidBufs {
        let plane = g.brick(0).n_theta() as usize * BRICK_CELLS;
        let nb = g.n_bricks();
        let stale = self
            .sb
            .as_ref()
            .is_none_or(|s| s.t0.len() != nb || s.t0.first().is_none_or(|v| v.len() != plane));
        if stale {
            let mk = || vec![vec![0.0f64; plane]; nb];
            self.sb = Some(SolidBufs {
                t0: mk(),
                heat0: mk(),
                heat_prev: mk(),
                heat_cur: mk(),
                apply: mk(),
                diag: mk(),
                b: mk(),
                r: mk(),
                z: mk(),
                q: mk(),
                delta: mk(),
            });
        }
        self.sb.as_mut().expect("just ensured")
    }

    /// Gas class-`D` buffers (per-brick θ-plane sizes since S9). Shape
    /// staleness (brick count + per-brick plane), like the reaction
    /// buffers: this is sound ONLY because every buffer is fully rewritten
    /// before it is read each step — nothing carries across steps. A
    /// future warm start (persisting `sol` between steps) must strengthen
    /// this test first, or a same-shape grid swap would silently read
    /// stale operands (S3 review finding).
    fn ensure_gb(&mut self, g: &Grid) -> &mut GasBufs {
        let nb = g.n_bricks();
        let plane = |bi: usize| g.brick(bi).n_theta() as usize * BRICK_CELLS;
        let stale = self.gb.as_ref().is_none_or(|s| {
            s.d0.len() != nb || s.d0.iter().enumerate().any(|(bi, v)| v.len() != plane(bi))
        });
        if stale {
            let mkc = || -> Vec<Vec<Cons>> {
                (0..nb).map(|bi| vec![[0.0f64; NCOMP]; plane(bi)]).collect()
            };
            self.gb = Some(GasBufs {
                sol: GasOperands::alloc(g),
                lag: GasOperands::alloc(g),
                work: GasWork::alloc(g),
                tr: GasTransportField::alloc(g),
                d0: mkc(),
                dprev: mkc(),
                dcur: mkc(),
                dlag: mkc(),
                dstage: mkc(),
            });
        }
        self.gb.as_mut().expect("just ensured")
    }

    /// Class-`R` buffers (S7; per-brick θ-plane sizes since S8). Shape
    /// staleness only, like the gas buffers: every buffer is fully
    /// rewritten before it is read each step.
    fn ensure_rb(&mut self, g: &Grid) -> &mut ReactionBufs {
        let nb = g.n_bricks();
        let plane = |bi: usize| g.brick(bi).n_theta() as usize * BRICK_CELLS;
        let stale = self.rb.as_ref().is_none_or(|s| {
            s.r0.len() != nb || s.r0.iter().enumerate().any(|(bi, v)| v.len() != plane(bi))
        });
        if stale {
            let mk = || -> Vec<Vec<f64>> { (0..nb).map(|bi| vec![0.0f64; plane(bi)]).collect() };
            self.rb = Some(ReactionBufs {
                r0: mk(),
                r_prev: mk(),
                r_trial: mk(),
                x_new: mk(),
            });
        }
        self.rb.as_mut().expect("just ensured")
    }

    /// Pure-diffusion convenience (the Goal-A studies): the same step with
    /// only class `D` scheduled.
    pub fn step_diffusion(
        &mut self,
        g: &mut Grid,
        diffusion: &DiffusionClass<'_, '_>,
        t: f64,
        dt: f64,
    ) -> Result<StepReport, SdcError> {
        self.step::<crate::euler::GammaLaw>(g, None, Some(diffusion), None, None, None, t, dt)
    }

    /// March `n_steps` of the pure-diffusion step from `t0`; returns the
    /// final time.
    pub fn advance_diffusion(
        &mut self,
        g: &mut Grid,
        diffusion: &DiffusionClass<'_, '_>,
        t0: f64,
        dt: f64,
        n_steps: usize,
    ) -> Result<f64, SdcError> {
        let mut t = t0;
        for _ in 0..n_steps {
            self.step_diffusion(g, diffusion, t, dt)?;
            t += dt;
        }
        Ok(t)
    }

    /// Flow-only convenience: one step with only class `A` scheduled.
    pub fn step_flow<E: EosLaw + Sync>(
        &mut self,
        g: &mut Grid,
        flow: &FlowClass<'_, '_, E>,
        t: f64,
        dt: f64,
    ) -> Result<StepReport, SdcError> {
        self.step(g, Some(flow), None, None, None, None, t, dt)
    }

    /// CFL-paced flow-only march from `t0` to `t_final`; returns steps taken.
    pub fn march_flow<E: EosLaw + Sync>(
        &mut self,
        g: &mut Grid,
        flow: &FlowClass<'_, '_, E>,
        cfl: f64,
        t0: f64,
        t_final: f64,
    ) -> Result<usize, SdcError> {
        let mut t = t0;
        let mut steps = 0usize;
        while t < t_final {
            let dt = self.stable_dt(g, flow, cfl)?.min(t_final - t);
            if t + dt == t {
                // dt fell below one ulp of t: the march can no longer
                // advance (CFL collapse or a mis-scaled t0) — refuse loudly
                // rather than spin forever (S2 review finding).
                return Err(SdcError::Config(
                    "march stalled: dt below one ulp of t (CFL collapse or mis-scaled t0)",
                ));
            }
            self.step_flow(g, flow, t, dt)?;
            t += dt;
            steps += 1;
        }
        Ok(steps)
    }

    /// One SDC-IMEX step of size `dt` at time `t` (module doc). The
    /// schedule is the fixed source order: class `A` evaluation, then the
    /// class-`D` solves (gas F_visc, then solid conduction with the
    /// Robin-Robin exchange) inside the fixed Picard sweeps, then the
    /// class-`R` cell-local implicit reaction on the accepted composition,
    /// per SDC sweep; the COUP-2 audit closes the step.
    #[allow(clippy::too_many_lines, clippy::too_many_arguments)]
    pub fn step<E: EosLaw + Sync>(
        &mut self,
        g: &mut Grid,
        flow: Option<&FlowClass<'_, '_, E>>,
        diffusion: Option<&DiffusionClass<'_, '_>>,
        gas: Option<&GasDiffusionClass<'_>>,
        exchange: Option<&ExchangeClass<'_>>,
        reaction: Option<&ReactionClass<'_>>,
        t: f64,
        dt: f64,
    ) -> Result<StepReport, SdcError> {
        // --- Schedule validation (loud, META-1 P6) -----------------------
        if flow.is_none() && diffusion.is_none() {
            return Err(SdcError::Config("no operator class scheduled"));
        }
        if reaction.is_some() && flow.is_none() {
            return Err(SdcError::Config(
                "the reaction class rides the flow state; schedule it with \
                 the hyperbolic class",
            ));
        }
        if exchange.is_some() && (flow.is_none() || diffusion.is_none()) {
            return Err(SdcError::Config(
                "wall exchange needs both the hyperbolic and diffusion classes",
            ));
        }
        if gas.is_some() && flow.is_none() {
            return Err(SdcError::Config(
                "the gas diffusion class rides the flow state; schedule it with \
                 the hyperbolic class",
            ));
        }
        if gas.is_some() {
            // Uniform N_θ ≥ 1 is legal since plan S9 (the full θ-stress
            // tensor, COUP-3 0.4.5 as built).
            let nt0 = g.brick(0).n_theta();
            if g.bricks().iter().any(|b| b.n_theta() != nt0) {
                return Err(SdcError::Config(
                    "gas diffusion (class D) at MIXED N_θ: the ring-interface coupling \
                     of an implicit operator rides plan S11 (COUP-3 0.4.5); uniform N_θ \
                     only — refusing rather than guessing",
                ));
            }
            if nt0 != 1 && g.has_cut_geometry() {
                return Err(SdcError::Config(
                    "gas diffusion (class D) with cut geometry at N_θ > 1 rides plan \
                     S11 (θ-face apertures + per-θ wall patches); refusing rather than \
                     guessing",
                ));
            }
        }
        if let Some(dc) = diffusion
            && dc.t_field == dc.scratch_field
        {
            return Err(SdcError::Config(
                "the class-D scratch field must be distinct from the temperature \
                 field (the CG direction vector would clobber the solve's warm start)",
            ));
        }
        if let (Some(dc), Some(_)) = (diffusion, flow) {
            if dc.op.domain != Domain::Solid {
                return Err(SdcError::Config(
                    "a scalar gas-domain conduction class alongside the flow class is \
                     superseded by the S3 gas-diffusion class (which owns the gas \
                     energy); this schedule couples flow to SOLID conduction only",
                ));
            }
            if g.bricks().iter().any(|b| b.n_theta() != 1) {
                // All-brick scan (S8 review): brick 0 is not representative
                // on a disconnected mixed world — a partial check panics in
                // the buffer sizing instead of refusing.
                return Err(SdcError::Config(
                    "coupled flow+conduction at N_θ > 1 rides the coarse-3-D RL10 wave \
                     (plan S11: per-θ wall patches + exchange keying); refusing rather \
                     than guessing",
                ));
            }
        }

        // --- Stored totals BEFORE (audit LHS) ----------------------------
        let stored_before = self.audit_stored(g, flow, diffusion);

        // --- Node-0 snapshots and evaluations ----------------------------
        if let Some(fc) = flow {
            let ws = self.ensure_ws(g, fc)?;
            let ids = fc.fields.ids();
            for (bi, u0) in ws.0.u0.iter_mut().enumerate() {
                let b = g.brick(bi);
                for k in 0..NCOMP {
                    let src = b.field(ids[k]);
                    for (cell, u) in u0.iter_mut().enumerate() {
                        u[k] = src[cell];
                    }
                }
            }
        }
        if let Some(dc) = diffusion {
            let sb = self.ensure_sb(g);
            for (bi, t0) in sb.t0.iter_mut().enumerate() {
                t0.copy_from_slice(g.brick(bi).field(dc.t_field));
            }
        }

        // Exchange operands at node 0 (pre-step gas state) + heatW⁰.
        let mut ex_ops: Vec<(f64, f64)> = Vec::new(); // (h_film, t_aw) per patch
        let mut ex_map: BTreeMap<ExchangeKey, GasFaceRobin> = BTreeMap::new();
        let mut q0: Vec<f64> = Vec::new();
        let mut hl0 = HeatLedger::default();
        if let (Some(xc), Some(fc)) = (exchange, flow) {
            ex_ops = self.exchange_operands(g, fc, xc)?;
            ex_map = build_ex_map(xc, &ex_ops);
        }
        if let Some(dc) = diffusion {
            let sb = self.sb.as_mut().expect("ensured");
            let mut rec: Vec<(ExchangeKey, f64)> = Vec::new();
            dc.op.assemble_heat(
                g,
                dc.t_field,
                t,
                AssembleMode::Affine,
                exchange.map(|_| &ex_map),
                &mut sb.heat0,
                None,
                Some(&mut hl0),
                exchange.map(|_| &mut rec),
            )?;
            if let Some(xc) = exchange {
                q0 = patch_heats(xc.patches, &rec);
            }
        }
        let mut l0 = FlowLedger::default();
        if let Some(fc) = flow {
            let ws = self.ws.as_mut().expect("ensured");
            fc.op.eval_rhs(g, fc.fields, ws, t)?;
            l0 = *ws.ledger();
            for (bi, dst) in self.rate_e0.iter_mut().enumerate() {
                dst.copy_from_slice(&ws.0.rate[bi]);
            }
        }

        // Class-R node-0 rates R(U⁰) (S7): the auto-ignition term's explicit
        // evaluation at the pre-step state, from the primitive cache eval_rhs
        // just filled — the trapezoid quadrature's +Δt/2·R⁰ operand.
        if let (Some(rc), Some(_)) = (reaction, flow) {
            self.ensure_rb(g);
            let (rb, ws) = (
                self.rb.as_mut().expect("ensured"),
                self.ws.as_ref().expect("ensured"),
            );
            for bi in 0..g.n_bricks() {
                let brick = g.brick(bi);
                let mask = brick.mask();
                let nt = brick.n_theta();
                let (r0, prim) = (&mut rb.r0[bi], &ws.0.prim[bi]);
                for j in 0..nt {
                    for local in 0..BRICK_CELLS {
                        let idx = j as usize * BRICK_CELLS + local;
                        r0[idx] = 0.0;
                        if mask & (1u64 << local) == 0 {
                            continue;
                        }
                        let (i_r, i_z) = brick.global_rz(local);
                        // The node-0 rate is the REALIZABLE one (capped at the
                        // guarded parking point over dt) — see `auto_rate_node0`.
                        r0[idx] = rc.op.auto_rate_node0(&prim[idx], dt).map_err(|what| {
                            FlowError::NonPhysicalState {
                                i_r,
                                i_z,
                                i_theta: j,
                                what,
                            }
                        })?;
                    }
                }
                // r_prev seeds as R⁰ (the predictor's wrprev weight is zero,
                // so it is never read then; sweep 1 rolls in the predictor's
                // realized rates first) — the heat_prev/dprev pattern.
                rb.r_prev[bi].copy_from_slice(&rb.r0[bi]);
            }
        }

        // Gas class-D at node 0: D(U⁰) rates + ledger (eval_rhs above just
        // filled the primitive cache from U⁰ — the operand source).
        let mut gd0 = FlowLedger::default();
        if let (Some(gc), Some(_)) = (gas, flow) {
            self.ensure_gb(g);
            let (gb, ws) = (
                self.gb.as_mut().expect("ensured"),
                self.ws.as_ref().expect("ensured"),
            );
            derive_gas_operands(
                g,
                &ws.0.prim,
                gc.temperature,
                gc.transport,
                &mut gb.sol,
                &mut gb.tr,
                &mut gb.work.ke_base,
            )?;
            gb.lag.clone_from(&gb.sol);
            gc.op.fill_lag_gradients(g, &gb.lag, &mut gb.work)?;
            gc.op.assemble_rates(
                g,
                &gb.sol,
                &gb.lag,
                &gb.work,
                &gb.tr,
                t,
                &mut gb.d0,
                Some(&mut gd0),
            )?;
            // dprev starts as D(U⁰) (the predictor's weights are zero —
            // the buffer must hold SOMETHING shaped right; sweep 1 reads
            // the predictor's accepted rates), dlag likewise as the first
            // trial. The heat_prev/q_prev pattern, five components wide.
            for (bi, src) in gb.d0.iter().enumerate() {
                gb.dprev[bi].copy_from_slice(src);
                gb.dlag[bi].copy_from_slice(src);
            }
        }

        // --- The fixed sweeps --------------------------------------------
        let mut l_last = l0;
        let mut hl_prev = hl0;
        let mut hl_last = HeatLedger::default();
        let mut gd_prev = gd0;
        let mut gd_last = FlowLedger::default();
        let mut gd_iter = FlowLedger::default();
        let mut q_prev = q0.clone();
        let mut debit_applied = 0.0f64; // κV-weighted gas energy debited (J)
        let mut burn_applied = 0.0f64; // κV-weighted class-R ρb applied (kg)
        let mut burn_gross = 0.0f64; // gross magnitude for the audit tolerance
        let mut report = StepReport::default();
        let mut ex_report = ExchangeStepReport::default();

        if diffusion.is_some() {
            // heat_prev starts as heatW⁰ (the predictor doesn't read it —
            // its weights are zero — but the buffer must hold SOMETHING
            // shaped right; sweep 1 reads the predictor's accepted value).
            let sb = self.sb.as_mut().expect("ensured");
            for (bi, dst) in sb.heat_prev.iter_mut().enumerate() {
                dst.copy_from_slice(&sb.heat0[bi]);
            }
        }

        for sweep in 0..=N_SDC_CORRECTIONS {
            let predictor = sweep == 0;
            // Weights of this sweep's composition (module doc).
            let (we0, we1) = if predictor {
                (dt, 0.0)
            } else {
                (0.5 * dt, 0.5 * dt)
            };
            let (wq0, wqprev) = if predictor {
                (0.0, 0.0)
            } else {
                (0.5 * dt, -0.5 * dt)
            };
            let wqnew = dt;

            // Class A at the current node-1 iterate (corrections only; the
            // predictor reuses F_E⁰ — the grid still holds U⁰ then).
            if !predictor && let Some(fc) = flow {
                let ws = self.ws.as_mut().expect("ensured");
                fc.op.eval_rhs(g, fc.fields, ws, t + dt)?;
                l_last = *ws.ledger();
            }

            // The fixed Picard loop (single pass when neither coupled
            // solve — wall exchange or gas cross terms — is scheduled).
            let n_picard = if exchange.is_some() || gas.is_some() {
                N_ROBIN_SWEEPS
            } else {
                1
            };
            // Per-sweep trials: the previous sweep's accepted values.
            let mut q_trial = q_prev.clone();
            if gas.is_some() {
                let gb = self.gb.as_mut().expect("ensured");
                for bi in 0..gb.dlag.len() {
                    let (dlag, dprev) = (&mut gb.dlag[bi], &gb.dprev[bi]);
                    dlag.copy_from_slice(dprev);
                }
            }
            let mut omega = 1.0f64;
            let mut resid_hist: Option<Vec<f64>> = None; // previous residual vector
            let mut robin_resid = 0.0f64;
            let mut gas_resid = 0.0f64;
            let mut q_new: Vec<f64> = Vec::new();
            let mut hl_tmp = HeatLedger::default();

            for ps in 0..n_picard {
                if let Some(fc) = flow {
                    self.compose_gas(
                        g,
                        fc,
                        we0,
                        we1,
                        exchange,
                        wq0,
                        &q0,
                        wqprev,
                        &q_prev,
                        wqnew,
                        &q_trial,
                        gas.is_some(),
                    );
                    let ws = self.ws.as_ref().expect("ensured");
                    fc.op.apply_srd(g, fc.fields, ws);
                }
                // Gas class-D solves (S3): operate on module buffers only
                // (the state advances through the composition above, next
                // iterate, with the freshly-accepted rates in dlag).
                if let (Some(gc), Some(fc)) = (gas, flow) {
                    {
                        let ws = self.ws.as_mut().expect("ensured");
                        fc.op.refresh_prims(g, fc.fields, ws)?;
                    }
                    let gb = self.gb.as_mut().expect("ensured");
                    let ws = self.ws.as_ref().expect("ensured");
                    // Cross-term lag = the previous iterate's solutions
                    // (iterate 0: the composed base itself).
                    if ps > 0 {
                        std::mem::swap(&mut gb.lag, &mut gb.sol);
                    }
                    derive_gas_operands(
                        g,
                        &ws.0.prim,
                        gc.temperature,
                        gc.transport,
                        &mut gb.sol,
                        &mut gb.tr,
                        &mut gb.work.ke_base,
                    )?;
                    if ps == 0 {
                        gb.lag.clone_from(&gb.sol);
                    }
                    gc.op.fill_lag_gradients(g, &gb.lag, &mut gb.work)?;
                    // Stage rates at the base operands (the affine RHS of
                    // the velocity solves), then the fixed solve order:
                    // u_r, u_z, ω — then T (work fluxes at the accepted
                    // velocities) and C off a re-staged assembly.
                    gc.op.assemble_rates(
                        g,
                        &gb.sol,
                        &gb.lag,
                        &gb.work,
                        &gb.tr,
                        t + dt,
                        &mut gb.dstage,
                        None,
                    )?;
                    let mut it_max = 0usize;
                    let mut rs_max = 0.0f64;
                    for comp in [GasComp::Ur, GasComp::Uz, GasComp::Om] {
                        fill_gas_rhs(
                            g,
                            comp,
                            wqnew,
                            &gb.dstage,
                            &gb.dlag,
                            &gb.sol,
                            &gb.work.ke_base,
                            &mut gb.work.b,
                        );
                        gc.op.fill_mass(g, comp, &gb.sol.rho, &gb.tr, &mut gb.work);
                        let x = match comp {
                            GasComp::Ur => &mut gb.sol.ur,
                            GasComp::Uz => &mut gb.sol.uz,
                            _ => &mut gb.sol.om,
                        };
                        let (it, rs) = gc.op.cg_solve(g, comp, &gb.tr, wqnew, x, &mut gb.work)?;
                        it_max = it_max.max(it);
                        rs_max = rs_max.max(rs);
                    }
                    gc.op.assemble_rates(
                        g,
                        &gb.sol,
                        &gb.lag,
                        &gb.work,
                        &gb.tr,
                        t + dt,
                        &mut gb.dstage,
                        None,
                    )?;
                    for comp in [GasComp::T, GasComp::C] {
                        fill_gas_rhs(
                            g,
                            comp,
                            wqnew,
                            &gb.dstage,
                            &gb.dlag,
                            &gb.sol,
                            &gb.work.ke_base,
                            &mut gb.work.b,
                        );
                        gc.op.fill_mass(g, comp, &gb.sol.rho, &gb.tr, &mut gb.work);
                        let x = match comp {
                            GasComp::T => &mut gb.sol.tt,
                            _ => &mut gb.sol.cc,
                        };
                        let (it, rs) = gc.op.cg_solve(g, comp, &gb.tr, wqnew, x, &mut gb.work)?;
                        it_max = it_max.max(it);
                        rs_max = rs_max.max(rs);
                    }
                    // The accepted assembly at the solutions — the rates
                    // the composition applies and the audit ledgers.
                    gd_iter = FlowLedger::default();
                    gc.op.assemble_rates(
                        g,
                        &gb.sol,
                        &gb.lag,
                        &gb.work,
                        &gb.tr,
                        t + dt,
                        &mut gb.dcur,
                        Some(&mut gd_iter),
                    )?;
                    gas_resid = rate_resid(g, &fc.fields.ids(), &gb.dcur, &gb.dlag, wqnew);
                    std::mem::swap(&mut gb.dlag, &mut gb.dcur);
                    report.gas_cg_iters = it_max;
                    report.gas_cg_resid = rs_max;
                }
                if let (Some(xc), Some(fc)) = (exchange, flow) {
                    ex_ops = self.exchange_operands(g, fc, xc)?;
                    ex_map = build_ex_map(xc, &ex_ops);
                }
                if let Some(dc) = diffusion {
                    let (iters, resid) = self.class_d_solve(
                        g,
                        dc,
                        exchange.map(|_| &ex_map),
                        t + dt,
                        dt,
                        we0 != dt,
                    )?;
                    report.cg_iters = iters;
                    report.cg_resid = resid;
                    // The accepted assembly at the solution — the single
                    // owner of this sweep's F_I record AND the per-face
                    // exchange heats the gas must debit (conservation by
                    // construction).
                    let sb = self.sb.as_mut().expect("ensured");
                    let mut rec: Vec<(ExchangeKey, f64)> = Vec::new();
                    hl_tmp = HeatLedger::default();
                    dc.op.assemble_heat(
                        g,
                        dc.t_field,
                        t + dt,
                        AssembleMode::Affine,
                        exchange.map(|_| &ex_map),
                        &mut sb.heat_cur,
                        None,
                        Some(&mut hl_tmp),
                        exchange.map(|_| &mut rec),
                    )?;
                    if let Some(xc) = exchange {
                        q_new = patch_heats(xc.patches, &rec);
                    }
                }
                if exchange.is_some() {
                    // Picard residual + clamped Aitken relaxation on the
                    // exchange-heat vector (COUP-2 §3.5's accelerated fixed
                    // sweeps; scalar secant ω over the residual vector).
                    let resid_vec: Vec<f64> =
                        q_new.iter().zip(&q_trial).map(|(n, o)| n - o).collect();
                    robin_resid = rel_resid(&q_new, &resid_vec);
                    let last = ps + 1 == n_picard;
                    if last {
                        q_trial = q_new.clone();
                    } else {
                        if let Some(prev) = &resid_hist {
                            let num: f64 =
                                prev.iter().zip(&resid_vec).map(|(a, b)| a * (b - a)).sum();
                            let den: f64 = prev
                                .iter()
                                .zip(&resid_vec)
                                .map(|(a, b)| (b - a) * (b - a))
                                .sum();
                            if den > 0.0 {
                                omega =
                                    (-omega * num / den).clamp(ROBIN_OMEGA_MIN, ROBIN_OMEGA_MAX);
                            }
                        }
                        for (qt, dr) in q_trial.iter_mut().zip(&resid_vec) {
                            *qt += omega * dr;
                        }
                        resid_hist = Some(resid_vec);
                    }
                }
            }
            if exchange.is_some() {
                // NaN-safe acceptance: a NaN residual is a violation too.
                if robin_resid.is_nan() || robin_resid > EPS_ROBIN_RESID {
                    return Err(SdcError::CouplingResidual {
                        solve: "Robin-Robin wall exchange (COUP-2 §3.5)",
                        resid: robin_resid,
                        eps: EPS_ROBIN_RESID,
                    });
                }
                ex_report.robin_resid = robin_resid;
            }
            if gas.is_some() {
                // NaN-safe acceptance of the gas cross-term Picard (S3).
                if gas_resid.is_nan() || gas_resid > EPS_GAS_DIFF_RESID {
                    return Err(SdcError::CouplingResidual {
                        solve: "gas class-D cross-term Picard (COUP-3 §3.1, S3)",
                        resid: gas_resid,
                        eps: EPS_GAS_DIFF_RESID,
                    });
                }
                report.gas_picard_resid = gas_resid;
            }

            // Accept the sweep: final gas composition with the ACCEPTED
            // exchange heats (exactly what the solid solve received) and
            // the ACCEPTED gas diffusion rates (in dlag after the swap).
            if let Some(fc) = flow {
                self.compose_gas(
                    g,
                    fc,
                    we0,
                    we1,
                    exchange,
                    wq0,
                    &q0,
                    wqprev,
                    &q_prev,
                    wqnew,
                    &q_trial,
                    gas.is_some(),
                );
                let ws = self.ws.as_ref().expect("ensured");
                fc.op.apply_srd(g, fc.fields, ws);
            }
            // Class-R on the accepted composition (S7): the cell-local
            // implicit auto-ignition node solve, applied once per sweep so
            // the next sweep's class-A evaluation (the propagation reaction)
            // sees the auto-ignited b — the seed→propagate coupling. The
            // realized rates roll into the trapezoid quadrature exactly as
            // the exchange heats do (wr0·R⁰ + wrprev·R_prev added to the
            // base; wrnew·R_new realized by the solve itself).
            if let (Some(rc), Some(fc)) = (reaction, flow) {
                let (net, gross) = self.apply_reaction(g, fc, rc, wq0, wqprev, wqnew)?;
                burn_applied = net;
                burn_gross = gross;
            }
            let last_sweep = sweep == N_SDC_CORRECTIONS;
            if last_sweep {
                // The audit's applied-increment records (final composition).
                if exchange.is_some() {
                    let sum = |v: &[f64]| -> f64 { v.iter().sum() };
                    let applied_exchange_j =
                        wq0 * sum(&q0) + wqprev * sum(&q_prev) + wqnew * sum(&q_trial);
                    debit_applied = -applied_exchange_j;
                    let applied_exterior_j = wq0 * hl0.exterior_w
                        + wqprev * hl_prev.exterior_w
                        + wqnew * hl_tmp.exterior_w;
                    let applied_bc_j = wq0 * hl0.bc_w + wqprev * hl_prev.bc_w + wqnew * hl_tmp.bc_w;
                    ex_report.q_w = q_trial.clone();
                    ex_report.h = ex_ops.iter().map(|(h, _)| *h).collect();
                    ex_report.t_aw = ex_ops.iter().map(|(_, ta)| *ta).collect();
                    ex_report.jacket_w = q_trial.iter().sum();
                    ex_report.applied_exchange_j = applied_exchange_j;
                    ex_report.applied_exterior_j = applied_exterior_j;
                    ex_report.applied_bc_j = applied_bc_j;
                }
                hl_last = hl_tmp;
                gd_last = gd_iter;
            } else {
                if diffusion.is_some() {
                    // Roll the F_I record: heat_prev ← this sweep's
                    // accepted assembly; q_prev/hl_prev likewise.
                    let sb = self.sb.as_mut().expect("ensured");
                    std::mem::swap(&mut sb.heat_prev, &mut sb.heat_cur);
                    hl_prev = hl_tmp;
                    q_prev = q_trial.clone();
                }
                if gas.is_some() {
                    // Roll the gas D record: dprev ← this sweep's accepted
                    // rates (sitting in dlag after the iterate swap).
                    let gb = self.gb.as_mut().expect("ensured");
                    std::mem::swap(&mut gb.dprev, &mut gb.dlag);
                    gd_prev = gd_iter;
                }
                if reaction.is_some() {
                    // Roll the class-R record: r_prev ← this sweep's
                    // realized rates.
                    let rb = self.rb.as_mut().expect("ensured");
                    std::mem::swap(&mut rb.r_prev, &mut rb.r_trial);
                }
            }
        }

        // --- The COUP-2 audit (module doc) --------------------------------
        let stored_after = self.audit_stored(g, flow, diffusion);
        self.audit_check(
            g,
            flow.is_some(),
            diffusion.is_some(),
            gas.is_some(),
            dt,
            &stored_before,
            &stored_after,
            &l0,
            &l_last,
            &hl0,
            &hl_prev,
            &hl_last,
            &gd0,
            &gd_prev,
            &gd_last,
            debit_applied,
            burn_applied,
            burn_gross,
            &mut report,
        )?;
        if exchange.is_some() {
            report.exchange = Some(ex_report);
        }
        Ok(report)
    }

    /// The class-`R` application (S7): per active gas cell, the implicit
    /// auto-ignition node solve on the accepted composition — two passes
    /// (read + solve into staging, then write) so the borrow of the grid's
    /// conserved fields stays clean. Returns the κV-weighted (net, gross)
    /// applied-ρb record for the audit's `burn_progress` row. Serial,
    /// fixed order, bit-deterministic.
    fn apply_reaction<E: EosLaw + Sync>(
        &mut self,
        g: &mut Grid,
        fc: &FlowClass<'_, '_, E>,
        rc: &ReactionClass<'_>,
        wr0: f64,
        wrprev: f64,
        wrnew: f64,
    ) -> Result<(f64, f64), SdcError> {
        let ids = fc.fields.ids();
        let rb = self.rb.as_mut().expect("ensured");
        let mut net = 0.0f64;
        let mut gross = 0.0f64;
        for bi in 0..g.n_bricks() {
            let brick = g.brick(bi);
            let mask = brick.mask();
            let nt = brick.n_theta();
            let fields: [&[f64]; NCOMP] = std::array::from_fn(|k| brick.field(ids[k]));
            let (r0, r_prev, r_trial, x_new) = (
                &rb.r0[bi],
                &rb.r_prev[bi],
                &mut rb.r_trial[bi],
                &mut rb.x_new[bi],
            );
            for j in 0..nt {
                for local in 0..BRICK_CELLS {
                    let idx = j as usize * BRICK_CELLS + local;
                    r_trial[idx] = 0.0;
                    x_new[idx] = 0.0;
                    if mask & (1u64 << local) == 0 {
                        continue;
                    }
                    let (i_r, i_z) = brick.global_rz(local);
                    let kappa = brick.kappa_rz(local);
                    if kappa <= 0.0 {
                        x_new[idx] = fields[I_RB][idx];
                        continue;
                    }
                    let u: Cons = std::array::from_fn(|k| fields[k][idx]);
                    let base = u[I_RB] + wr0 * r0[idx] + wrprev * r_prev[idx];
                    let (x, r) = rc
                        .op
                        .implicit_auto_update(&u, base, wrnew)
                        .map_err(|what| FlowError::NonPhysicalState {
                            i_r,
                            i_z,
                            i_theta: j,
                            what,
                        })?;
                    x_new[idx] = x;
                    r_trial[idx] = r;
                    let kv = kappa * g.cell_volume(i_r, nt);
                    net += kv * (x - u[I_RB]);
                    gross += kv
                        * ((wr0 * r0[idx]).abs()
                            + (wrprev * r_prev[idx]).abs()
                            + (wrnew * r).abs());
                }
            }
        }
        for bi in 0..g.n_bricks() {
            let (mask, nt) = (g.brick(bi).mask(), g.brick(bi).n_theta());
            let x_new = &rb.x_new[bi];
            let dst = g.brick_field_mut(bi, ids[I_RB]);
            for j in 0..nt {
                for local in 0..BRICK_CELLS {
                    if mask & (1u64 << local) != 0 {
                        // κ ≤ 0 cells staged their unchanged value (a bit-exact
                        // no-op write), so the masked write is unconditional.
                        let idx = j as usize * BRICK_CELLS + local;
                        dst[idx] = x_new[idx];
                    }
                }
            }
        }
        Ok((net, gross))
    }

    /// Gas node-1 composition: `U = U⁰ + we0·F_E⁰ + we1·F_E(cur)` plus the
    /// gas class-D quadrature `wq0·D⁰ + wqprev·D_prev + wqnew·D_trial`
    /// (S3; dlag holds the trial) plus the exchange debit
    /// `−(wq0·q0 + wqprev·q_prev + wqnew·q_trial)` per patch distributed
    /// uniformly over its SRD debit set. Fixed order; serial.
    #[allow(clippy::too_many_arguments)]
    fn compose_gas<E: EosLaw + Sync>(
        &self,
        g: &mut Grid,
        fc: &FlowClass<'_, '_, E>,
        we0: f64,
        we1: f64,
        exchange: Option<&ExchangeClass<'_>>,
        wq0: f64,
        q0: &[f64],
        wqprev: f64,
        q_prev: &[f64],
        wqnew: f64,
        q_trial: &[f64],
        gas: bool,
    ) {
        let ids = fc.fields.ids();
        let ws = self.ws.as_ref().expect("ensured");
        for bi in 0..g.n_bricks() {
            for (k, &id) in ids.iter().enumerate() {
                let u0 = &ws.0.u0[bi];
                let re0 = &self.rate_e0[bi];
                let rl = &ws.0.rate[bi];
                let dst = g.brick_field_mut(bi, id);
                if we1 == 0.0 {
                    for (cell, v) in dst.iter_mut().enumerate() {
                        *v = u0[cell][k] + we0 * re0[cell][k];
                    }
                } else {
                    for (cell, v) in dst.iter_mut().enumerate() {
                        *v = u0[cell][k] + we0 * re0[cell][k] + we1 * rl[cell][k];
                    }
                }
            }
        }
        if gas {
            // The gas diffusion rates (mass slot carries none — skipped).
            // This writes EVERY slot, masked cells included, relying on
            // `assemble_rates`' contract that non-gas slots are exact
            // zeros. If that contract ever broke, solid/exterior conserved
            // fields would drift SILENTLY (κ = 0 hides it from both the
            // audit and the stored reductions) — so assert it in debug
            // builds rather than trust it (S3 review finding).
            let gb = self.gb.as_ref().expect("gas buffers ensured");
            debug_assert!(
                (0..g.n_bricks()).all(|bi| {
                    let mask = g.brick(bi).mask();
                    (0..gb.d0[bi].len()).all(|c| {
                        mask & (1u64 << (c % BRICK_CELLS)) != 0
                            || (gb.d0[bi][c].iter().all(|v| *v == 0.0)
                                && gb.dprev[bi][c].iter().all(|v| *v == 0.0)
                                && gb.dlag[bi][c].iter().all(|v| *v == 0.0))
                    })
                }),
                "gas diffusion rates must be exactly zero outside the gas mask"
            );
            for bi in 0..g.n_bricks() {
                for (k, &id) in ids.iter().enumerate().skip(1) {
                    let (d0, dp, dl) = (&gb.d0[bi], &gb.dprev[bi], &gb.dlag[bi]);
                    let dst = g.brick_field_mut(bi, id);
                    for (cell, v) in dst.iter_mut().enumerate() {
                        *v += wq0 * d0[cell][k] + wqprev * dp[cell][k] + wqnew * dl[cell][k];
                    }
                }
            }
        }
        if let Some(xc) = exchange {
            let en = ids[I_EN];
            for (p, patch) in xc.patches.iter().enumerate() {
                let comp = wq0 * q0[p] + wqprev * q_prev[p] + wqnew * q_trial[p];
                let de = -comp / patch.debit_kv_sum;
                for &(i_r, i_z) in &patch.debit_cells {
                    cell_add(g, en, i_r, i_z, de);
                }
            }
        }
    }

    /// The gas-side exchange operands per patch: (h_film, T_aw) from the
    /// patch cell's CURRENT state projected on the interface normal
    /// (SOLV-1 §3.5 — local near-wall operands only).
    fn exchange_operands<E: EosLaw + Sync>(
        &self,
        g: &Grid,
        fc: &FlowClass<'_, '_, E>,
        xc: &ExchangeClass<'_>,
    ) -> Result<Vec<(f64, f64)>, SdcError> {
        let ids = fc.fields.ids();
        let (dr, dz) = (g.spec().dr, g.spec().dz);
        let mut out = Vec::with_capacity(xc.patches.len());
        for patch in xc.patches {
            let mut u = [0.0f64; NCOMP];
            for (k, id) in ids.iter().enumerate() {
                u[k] = cell_value(g, *id, patch.gas.0, patch.gas.1);
            }
            let w = fc
                .op
                .eos
                .prim_checked(&u)
                .map_err(|what| FlowError::NonPhysicalState {
                    i_r: patch.gas.0,
                    i_z: patch.gas.1,
                    i_theta: 0,
                    what,
                })?;
            let temperature = (xc.temperature)(&w).map_err(|what| FlowError::NonPhysicalState {
                i_r: patch.gas.0,
                i_z: patch.gas.1,
                i_theta: 0,
                what,
            })?;
            // The spine at the SAME near-wall cell state the temperature
            // came from — SOLV-1 §3.5, 0.4.1: the law states no transport
            // constant of its own.
            let tr = (xc.transport)(&w).map_err(|what| FlowError::NonPhysicalState {
                i_r: patch.gas.0,
                i_z: patch.gas.1,
                i_theta: 0,
                what,
            })?;
            let (n_r, n_z) = patch.n_hat;
            let v_n = w[1] * n_r + w[3] * n_z;
            let u_t = ((w[1] * w[1] + w[3] * w[3] - v_n * v_n).max(0.0) + w[2] * w[2]).sqrt();
            let y = 0.5 * (n_r.abs() * dr + n_z.abs() * dz);
            let gas = NearWallGas {
                rho: w[0],
                u_t,
                temperature,
                y,
                tr,
            };
            let h = xc.law.film_coefficient(&gas)?;
            let t_aw = xc.law.adiabatic_wall_temperature(&gas)?;
            out.push((h, t_aw));
        }
        Ok(out)
    }

    /// The class-`D` solve of one sweep (module doc): backward-Euler-form
    /// correction solved matrix-free by Jacobi-preconditioned CG with the
    /// fixed iteration structure. `trapezoid` selects the correction-sweep
    /// right-hand side (predictor otherwise). Solves in δ = x − x₀ (x₀ =
    /// the field's current content — the warm start), so a converged
    /// steady state costs zero iterations.
    fn class_d_solve(
        &mut self,
        g: &mut Grid,
        dc: &DiffusionClass<'_, '_>,
        ex_map: Option<&BTreeMap<ExchangeKey, GasFaceRobin>>,
        time: f64,
        dt: f64,
        trapezoid: bool,
    ) -> Result<(usize, f64), SdcError> {
        let nt = g.brick(0).n_theta();
        let rho_cp = dc.op.rho_cp;
        let sb = self.ensure_sb(g);

        // RHS assembly: heat_affine at x₀ (+ diag for the preconditioner).
        dc.op.assemble_heat(
            g,
            dc.t_field,
            time,
            AssembleMode::Affine,
            ex_map,
            &mut sb.heat_cur,
            Some(&mut sb.diag),
            None,
            None,
        )?;
        // b = ρc_pV·(t0 − x0) + [dt/2·heat0 − dt/2·heat_prev]_trap + dt·heat_affine(x0).
        let domain_solid = dc.op.domain == Domain::Solid;
        for bi in 0..g.n_bricks() {
            let brick = g.brick(bi);
            let dmask = if domain_solid {
                brick.solid_mask()
            } else {
                brick.mask()
            };
            let x0 = brick.field(dc.t_field);
            for j in 0..nt {
                for local in 0..BRICK_CELLS {
                    let idx = j as usize * BRICK_CELLS + local;
                    if dmask & (1u64 << local) == 0 {
                        sb.b[bi][idx] = 0.0;
                        continue;
                    }
                    let (i_r, _) = brick.global_rz(local);
                    let m = rho_cp * g.cell_volume(i_r, nt);
                    let mut b = m * (sb.t0[bi][idx] - x0[idx]) + dt * sb.heat_cur[bi][idx];
                    if trapezoid {
                        b += 0.5 * dt * sb.heat0[bi][idx] - 0.5 * dt * sb.heat_prev[bi][idx];
                    }
                    sb.b[bi][idx] = b;
                }
            }
        }

        // Jacobi diagonal of A_cg = ρc_pV − dt·Ã (diag_c ≤ 0 ⇒ positive).
        for bi in 0..g.n_bricks() {
            let brick = g.brick(bi);
            let dmask = if domain_solid {
                brick.solid_mask()
            } else {
                brick.mask()
            };
            for j in 0..nt {
                for local in 0..BRICK_CELLS {
                    let idx = j as usize * BRICK_CELLS + local;
                    if dmask & (1u64 << local) == 0 {
                        sb.diag[bi][idx] = 1.0; // never read; keeps 1/diag finite
                        continue;
                    }
                    let (i_r, _) = brick.global_rz(local);
                    let m = rho_cp * g.cell_volume(i_r, nt);
                    sb.diag[bi][idx] = m - dt * sb.diag[bi][idx];
                }
            }
        }

        // CG on δ (init 0), direction p in the scratch grid field.
        for v in sb.delta.iter_mut() {
            v.fill(0.0);
        }
        for (bi, rb) in sb.r.iter_mut().enumerate() {
            rb.copy_from_slice(&sb.b[bi]);
        }
        let masked = |g: &Grid, bi: usize, solid: bool| -> u64 {
            let b = g.brick(bi);
            if solid { b.solid_mask() } else { b.mask() }
        };
        let dot = |g: &Grid, a: &BufF, b: &BufF, solid: bool| -> f64 {
            let partials: Vec<f64> = (0..g.n_bricks())
                .map(|bi| {
                    let dmask = masked(g, bi, solid);
                    let mut acc = 0.0f64;
                    for j in 0..nt {
                        for local in 0..BRICK_CELLS {
                            if dmask & (1u64 << local) != 0 {
                                let idx = j as usize * BRICK_CELLS + local;
                                acc += a[bi][idx] * b[bi][idx];
                            }
                        }
                    }
                    acc
                })
                .collect();
            tree_combine(&partials)
        };

        let b_norm2 = dot(g, &sb.b, &sb.b, domain_solid);
        let eps2 = EPS_CG_RESID * EPS_CG_RESID * b_norm2;
        // z = r/diag; p = z (into the scratch field, zeroed outside domain).
        for bi in 0..g.n_bricks() {
            let (zb, rb, db) = (&mut sb.z[bi], &sb.r[bi], &sb.diag[bi]);
            for ((zv, rv), dv) in zb.iter_mut().zip(rb).zip(db) {
                *zv = rv / dv;
            }
        }
        for bi in 0..g.n_bricks() {
            // p = z everywhere: z is exactly 0 outside the domain (r = b = 0
            // there, diag = 1), so no stale scratch content survives.
            g.brick_field_mut(bi, dc.scratch_field)
                .copy_from_slice(&sb.z[bi]);
        }
        let mut rz = dot(g, &sb.r, &sb.z, domain_solid);
        let mut r_norm2 = dot(g, &sb.r, &sb.r, domain_solid);
        let mut iters = 0usize;
        while iters < N_CG_ITERS_MAX && r_norm2 > eps2 && rz > 0.0 {
            // q = A_cg·p = ρc_pV∘p − dt·Ã·p (Linear assembly: exact, no
            // affine constants to cancel).
            dc.op.assemble_heat(
                g,
                dc.scratch_field,
                time,
                AssembleMode::Linear,
                ex_map,
                &mut sb.apply,
                None,
                None,
                None,
            )?;
            for bi in 0..g.n_bricks() {
                let brick = g.brick(bi);
                let dmask = if domain_solid {
                    brick.solid_mask()
                } else {
                    brick.mask()
                };
                let p = brick.field(dc.scratch_field);
                for j in 0..nt {
                    for local in 0..BRICK_CELLS {
                        let idx = j as usize * BRICK_CELLS + local;
                        if dmask & (1u64 << local) == 0 {
                            sb.q[bi][idx] = 0.0;
                            continue;
                        }
                        let (i_r, _) = brick.global_rz(local);
                        let m = rho_cp * g.cell_volume(i_r, nt);
                        sb.q[bi][idx] = m * p[idx] - dt * sb.apply[bi][idx];
                    }
                }
            }
            let pq = {
                let partials: Vec<f64> = (0..g.n_bricks())
                    .map(|bi| {
                        let dmask = masked(g, bi, domain_solid);
                        let p = g.brick(bi).field(dc.scratch_field);
                        let mut acc = 0.0f64;
                        for j in 0..nt {
                            for local in 0..BRICK_CELLS {
                                if dmask & (1u64 << local) != 0 {
                                    let idx = j as usize * BRICK_CELLS + local;
                                    acc += p[idx] * sb.q[bi][idx];
                                }
                            }
                        }
                        acc
                    })
                    .collect();
                tree_combine(&partials)
            };
            if pq.is_nan() || pq <= 0.0 {
                break; // SPD breakdown ⇒ the acceptance check below decides
            }
            let alpha = rz / pq;
            for bi in 0..g.n_bricks() {
                let p = g.brick(bi).field(dc.scratch_field);
                for (d, pv) in sb.delta[bi].iter_mut().zip(p) {
                    *d += alpha * pv;
                }
            }
            for bi in 0..g.n_bricks() {
                let (rb, qb) = (&mut sb.r[bi], &sb.q[bi]);
                for (rv, qv) in rb.iter_mut().zip(qb) {
                    *rv -= alpha * qv;
                }
            }
            for bi in 0..g.n_bricks() {
                let (zb, rb, db) = (&mut sb.z[bi], &sb.r[bi], &sb.diag[bi]);
                for ((zv, rv), dv) in zb.iter_mut().zip(rb).zip(db) {
                    *zv = rv / dv;
                }
            }
            let rz_new = dot(g, &sb.r, &sb.z, domain_solid);
            let beta = rz_new / rz;
            rz = rz_new;
            for bi in 0..g.n_bricks() {
                let dst = g.brick_field_mut(bi, dc.scratch_field);
                for (pv, zv) in dst.iter_mut().zip(&sb.z[bi]) {
                    *pv = zv + beta * *pv;
                }
            }
            r_norm2 = dot(g, &sb.r, &sb.r, domain_solid);
            iters += 1;
        }
        let resid = (r_norm2 / b_norm2.max(f64::MIN_POSITIVE)).sqrt();
        // NaN-safe acceptance (a NaN right-hand side or residual is a
        // violation — `b_norm2 > 0.0` alone is false for NaN and would
        // silently accept; S2 review finding).
        if b_norm2.is_nan() || (b_norm2 > 0.0 && (resid.is_nan() || resid > EPS_CG_RESID)) {
            return Err(SdcError::CouplingResidual {
                solve: "class-D fixed-cycle CG (COUP-3 §3.1)",
                resid,
                eps: EPS_CG_RESID,
            });
        }
        // x = x₀ + δ over domain cells.
        for bi in 0..g.n_bricks() {
            let dmask = masked(g, bi, domain_solid);
            let dst = g.brick_field_mut(bi, dc.t_field);
            for j in 0..nt {
                for local in 0..BRICK_CELLS {
                    if dmask & (1u64 << local) != 0 {
                        let idx = j as usize * BRICK_CELLS + local;
                        dst[idx] += sb.delta[bi][idx];
                    }
                }
            }
        }
        Ok((iters, resid))
    }

    /// Stored totals for the audit rows: gas κV·U per component (when flow
    /// is scheduled) and the diffusion domain's ρc_p·V·T (when scheduled).
    /// The middle array is the GROSS magnitude `Σ κV·|q|` — the §3.1.1
    /// scale of what the reduction actually summed (S8: a cancelling
    /// stored total, e.g. mirror-symmetric θ-momentum, must not collapse
    /// the tolerance below the reduction's own rounding).
    fn audit_stored<E: EosLaw>(
        &self,
        g: &Grid,
        flow: Option<&FlowClass<'_, '_, E>>,
        diffusion: Option<&DiffusionClass<'_, '_>>,
    ) -> ([f64; NCOMP], [f64; NCOMP], f64) {
        let mut gas = [0.0f64; NCOMP];
        let mut gas_gross = [0.0f64; NCOMP];
        if let Some(fc) = flow {
            for (k, &id) in fc.fields.ids().iter().enumerate() {
                gas[k] = g.reduce_kappa_volume_weighted(id);
                gas_gross[k] = g.reduce_kappa_volume_weighted_abs(id);
            }
        }
        let mut solid = 0.0f64;
        if let Some(dc) = diffusion {
            solid = dc.op.rho_cp
                * match dc.op.domain {
                    Domain::Solid => g.reduce_solid_volume_weighted(dc.t_field),
                    Domain::FlowActive => g.reduce_volume_weighted(dc.t_field),
                };
        }
        (gas, gas_gross, solid)
    }

    /// The COUP-2 §3.1 identity per audited quantity, at the §3.1.1
    /// derived tolerance. The applied side composes the SAME node weights
    /// the final state was built with.
    #[allow(clippy::too_many_arguments)]
    fn audit_check(
        &self,
        g: &Grid,
        has_flow: bool,
        has_diffusion: bool,
        has_gas: bool,
        dt: f64,
        before: &([f64; NCOMP], [f64; NCOMP], f64),
        after: &([f64; NCOMP], [f64; NCOMP], f64),
        l0: &FlowLedger,
        l_last: &FlowLedger,
        hl0: &HeatLedger,
        hl_prev: &HeatLedger,
        hl_last: &HeatLedger,
        gd0: &FlowLedger,
        gd_prev: &FlowLedger,
        gd_last: &FlowLedger,
        debit_applied: f64,
        burn_applied: f64,
        burn_gross: f64,
        report: &mut StepReport,
    ) -> Result<(), SdcError> {
        let spec = &self.audit_spec;
        let eps = f64::EPSILON;
        let mut n_cells = 0usize;
        for b in g.bricks() {
            n_cells += (b.mask().count_ones() + b.solid_mask().count_ones()) as usize
                * b.n_theta() as usize;
        }
        let sqrt_n = (n_cells as f64).sqrt();
        // Final-composition weights (the fixed sweep structure guarantees
        // the last sweep is a trapezoid correction).
        let (we0, we1) = (0.5 * dt, 0.5 * dt);
        let (wh0, whprev, whnew) = (0.5 * dt, -0.5 * dt, dt);

        let solid_applied = if has_diffusion {
            wh0 * hl0.applied_w() + whprev * hl_prev.applied_w() + whnew * hl_last.applied_w()
        } else {
            0.0
        };
        let solid_gross = if has_diffusion {
            wh0.abs() * hl0.gross_w + whprev.abs() * hl_prev.gross_w + whnew.abs() * hl_last.gross_w
        } else {
            0.0
        };

        let mut rows: Vec<AuditRow> = Vec::new();
        let names = [
            "mass",
            "momentum_r",
            "momentum_theta",
            "momentum_z",
            "energy",
            "composition",
            "burn_progress",
        ];
        // Burn progress (k = I_RB) shares the composition floor scale
        // (ref_scale[3]): both are density × a [0,1] fraction. The per-step
        // throughput term S[q] already scales it independently.
        let ref_of = |k: usize| -> f64 {
            match k {
                0 => spec.ref_scale[0],
                1..=3 => spec.ref_scale[1],
                4 => spec.ref_scale[2],
                _ => spec.ref_scale[3],
            }
        };
        if has_flow {
            #[allow(clippy::needless_range_loop)] // k indexes five parallel arrays
            for k in 0..NCOMP {
                let mut delta = after.0[k] - before.0[k];
                let mut applied = we0 * (l0.port_net[k] + l0.src_net[k])
                    + we1 * (l_last.port_net[k] + l_last.src_net[k]);
                let mut gross = we0 * (l0.port_abs[k] + l0.src_abs[k])
                    + we1 * (l_last.port_abs[k] + l_last.src_abs[k]);
                if has_gas {
                    // The gas class-D quadrature (S3): same weights the
                    // final composition applied its rates with.
                    applied += wh0 * (gd0.port_net[k] + gd0.src_net[k])
                        + whprev * (gd_prev.port_net[k] + gd_prev.src_net[k])
                        + whnew * (gd_last.port_net[k] + gd_last.src_net[k]);
                    gross += wh0.abs() * (gd0.port_abs[k] + gd0.src_abs[k])
                        + whprev.abs() * (gd_prev.port_abs[k] + gd_prev.src_abs[k])
                        + whnew.abs() * (gd_last.port_abs[k] + gd_last.src_abs[k]);
                }
                // S8: the GROSS stored magnitude, both endpoints — a
                // cancelling net total must not collapse the tolerance
                // below the reduction's (and the per-cell composition's)
                // own rounding (COUP-2 §3.1.1: scaled by what was actually
                // summed; interior telescoping fluxes — the θ sweep — are
                // never ledgered, so their composed magnitudes appear only
                // through the after-state's gross content).
                let mut stored_scale = before.1[k] + after.1[k];
                if k == I_RB {
                    // The class-R applied increments (S7): the implicit
                    // auto-ignition node solve's realized composition,
                    // recorded exactly as applied (COUP-2).
                    applied += burn_applied;
                    gross += burn_gross;
                }
                if k == I_EN {
                    // The combined-energy row: the exchange pair cancels
                    // between the gas debit and the solid's received heats.
                    applied += debit_applied;
                    if has_diffusion {
                        delta += after.2 - before.2;
                        applied += solid_applied;
                        gross += debit_applied.abs() + solid_gross;
                        stored_scale += before.2.abs();
                    }
                }
                let s = stored_scale + gross;
                let tol =
                    (spec.k_audit * eps * sqrt_n * s).max(spec.k_audit * eps * sqrt_n * ref_of(k));
                rows.push(AuditRow {
                    quantity: names[k],
                    delta,
                    applied,
                    tol,
                });
            }
        } else if has_diffusion {
            let delta = after.2 - before.2;
            let s = before.2.abs() + solid_gross;
            let tol = (spec.k_audit * eps * sqrt_n * s)
                .max(spec.k_audit * eps * sqrt_n * spec.ref_scale[2]);
            rows.push(AuditRow {
                quantity: "energy",
                delta,
                applied: solid_applied,
                tol,
            });
        }
        for row in &rows {
            let gap = (row.delta - row.applied).abs();
            // NaN-safe: a NaN anywhere in the ledger is a violation.
            if gap.is_nan() || gap > row.tol {
                return Err(SdcError::AuditViolation {
                    quantity: row.quantity,
                    delta: row.delta,
                    applied: row.applied,
                    tol: row.tol,
                });
            }
        }
        report.audit = rows;
        Ok(())
    }
}

/// The exchange-residual acceptance's absolute heat floor [W] (S7): below
/// this per-patch magnitude the relative acceptance is meaningless — a
/// NEAR-VACUUM cold fill (the startup march's real pre-spark state, ~10²
/// Pa) exchanges micro-watts, where iteration noise is a large FRACTION of
/// a physically-nothing number. A centiwatt-class disagreement is orders
/// below the wall-function band's resolution at any exchange this
/// instrument certifies (station fixtures 10²–10⁴ W, the RL10 jacket
/// ~10⁷ W). One owner, named; a configuration whose REAL exchange lives at
/// this floor (micro-thruster class) must revisit it, loudly.
pub const EPS_ROBIN_Q_FLOOR_W: f64 = 1.0e-2;

/// Max relative residual of the exchange-heat vector, floored at
/// [`EPS_ROBIN_Q_FLOOR_W`] so a ~zero-heat state (near-vacuum cold start)
/// cannot fail on noise over nothing. NaN anywhere is
/// returned as NaN (S2 review finding: `f64::max` silently DROPS NaN, so
/// a fold alone would defeat the caller's NaN-safe acceptance).
fn rel_resid(q_new: &[f64], resid: &[f64]) -> f64 {
    if q_new.iter().chain(resid).any(|v| v.is_nan()) {
        return f64::NAN;
    }
    let scale = q_new
        .iter()
        .fold(0.0f64, |m, q| m.max(q.abs()))
        .max(EPS_ROBIN_Q_FLOOR_W);
    resid.iter().fold(0.0f64, |m, r| m.max(r.abs())) / scale
}

/// Per-patch heat sums from an assembly's recorded per-face exchange
/// heats (fixed patch/face order — the visit order is keyed, not assumed).
/// N_θ = 1 only (guarded upstream): at N_θ > 1 the per-θ records would
/// collapse under one key — the S8 wave re-keys this by (face, θ).
fn patch_heats(patches: &[WallPatch], rec: &[(ExchangeKey, f64)]) -> Vec<f64> {
    let map: BTreeMap<ExchangeKey, f64> = rec.iter().copied().collect();
    patches
        .iter()
        .map(|p| p.faces.iter().map(|f| map[&face_key(f)]).sum())
        .collect()
}

/// The exchange Robin map from per-patch operands.
fn build_ex_map(xc: &ExchangeClass<'_>, ops: &[(f64, f64)]) -> BTreeMap<ExchangeKey, GasFaceRobin> {
    let mut map = BTreeMap::new();
    for (patch, &(h, t_aw)) in xc.patches.iter().zip(ops) {
        for face in &patch.faces {
            map.insert(
                face_key(face),
                GasFaceRobin {
                    h_film: h,
                    t_aw,
                    area_scale: patch.area_scale,
                },
            );
        }
    }
    map
}

// --- Gas class-D helpers (S3) --------------------------------------------------

/// Derive the gas class-D operands from the primitive cache (the current
/// composed state): ρ, u_r, ω = u_θ/r̄, u_z, T (the FND-7 seam closure),
/// C — plus the per-cell base kinetic energy ½|u|² the T-solve's
/// dissipation bookkeeping needs. N_θ = 1 (validated upstream).
#[allow(clippy::too_many_arguments)]
fn derive_gas_operands(
    g: &Grid,
    prim: &[Vec<Prim>],
    temperature: &(dyn Fn(&Prim) -> Result<f64, &'static str> + Sync),
    transport: &(dyn Fn(&Prim) -> Result<TransportProps, &'static str> + Sync),
    sol: &mut GasOperands,
    tr: &mut GasTransportField,
    ke_base: &mut [Vec<f64>],
) -> Result<(), SdcError> {
    for bi in 0..g.n_bricks() {
        let b = g.brick(bi);
        let nt = b.n_theta() as usize;
        for jt in 0..nt {
            for local in 0..BRICK_CELLS {
                let idx = jt * BRICK_CELLS + local;
                if b.mask() & (1u64 << local) == 0 {
                    sol.rho[bi][idx] = 1.0; // never read; keeps masses finite
                    sol.ur[bi][idx] = 0.0;
                    sol.om[bi][idx] = 0.0;
                    sol.uz[bi][idx] = 0.0;
                    sol.tt[bi][idx] = 0.0;
                    sol.cc[bi][idx] = 0.0;
                    tr.set_masked(bi, idx);
                    ke_base[bi][idx] = 0.0;
                    continue;
                }
                let (i_r, i_z) = b.global_rz(local);
                let w = &prim[bi][idx];
                let tt = temperature(w).map_err(|what| {
                    SdcError::Flow(FlowError::NonPhysicalState {
                        i_r,
                        i_z,
                        i_theta: jt as u32,
                        what,
                    })
                })?;
                sol.rho[bi][idx] = w[0];
                sol.ur[bi][idx] = w[1];
                sol.om[bi][idx] = w[2] / g.r_center(i_r);
                sol.uz[bi][idx] = w[3];
                sol.tt[bi][idx] = tt;
                sol.cc[bi][idx] = w[5];
                // The spine reading at this cell's own state — the COUP-3
                // 0.4.3 refresh point. Its refusal carries the cell.
                let props = transport(w).map_err(|what| {
                    SdcError::Flow(FlowError::NonPhysicalState {
                        i_r,
                        i_z,
                        i_theta: jt as u32,
                        what,
                    })
                })?;
                tr.set(bi, idx, i_r, i_z, &props).map_err(SdcError::Gas)?;
                ke_base[bi][idx] = 0.5 * (w[1] * w[1] + w[2] * w[2] + w[3] * w[3]);
            }
        }
    }
    Ok(())
}

/// The right-hand side of one gas component's implicit solve (module doc
/// of `crate::gas_diffusion`): `b = wqnew·(unit·(R_affine(x₀) − d_lag))`,
/// in the component's solve units (κV for u_r/u_z/C; λ-units r̄κV for ω;
/// J for T, which also carries the dissipation bookkeeping
/// `ρκV·(ke_base − ke_new)` — the KE the momentum solves moved into
/// internal energy). Zero outside the gas mask.
#[allow(clippy::too_many_arguments)]
fn fill_gas_rhs(
    g: &Grid,
    comp: GasComp,
    wqnew: f64,
    dstage: &[Vec<Cons>],
    dlag: &[Vec<Cons>],
    sol: &GasOperands,
    ke_base: &[Vec<f64>],
    b: &mut [Vec<f64>],
) {
    let k = match comp {
        GasComp::Ur => I_MR,
        GasComp::Uz => I_MZ,
        GasComp::Om => I_MT,
        GasComp::T => I_EN,
        GasComp::C => I_RC,
    };
    for bi in 0..g.n_bricks() {
        let brick = g.brick(bi);
        let nt = brick.n_theta();
        for idx in 0..b[bi].len() {
            let local = idx % BRICK_CELLS;
            if brick.mask() & (1u64 << local) == 0 {
                b[bi][idx] = 0.0;
                continue;
            }
            let (i_r, _) = brick.global_rz(local);
            let kv = brick.kappa_rz(local) * g.cell_volume(i_r, nt);
            let ddiff = dstage[bi][idx][k] - dlag[bi][idx][k];
            let mut val = match comp {
                GasComp::Om => wqnew * g.r_center(i_r) * kv * ddiff,
                _ => wqnew * kv * ddiff,
            };
            if comp == GasComp::T {
                let r = g.r_center(i_r);
                let ut = sol.om[bi][idx] * r;
                let ke_new = 0.5
                    * (sol.ur[bi][idx] * sol.ur[bi][idx]
                        + ut * ut
                        + sol.uz[bi][idx] * sol.uz[bi][idx]);
                val += sol.rho[bi][idx] * kv * (ke_base[bi][idx] - ke_new);
            }
            b[bi][idx] = val;
        }
    }
}

/// S13c GPU cross-check for the FULL RESIDENT class-D gas Picard iterate
/// (doc-hidden; **additive** — runs the real `fill_lag_gradients` /
/// `assemble_rates` / `fill_gas_rhs` / `fill_mass` / `cg_solve` in the exact
/// SDC-inner order, changes no production number). Mirrors the sweep the
/// device `gpu_class_d_iterate` performs: assemble(sol,lag) → {Ur,Uz,Om}
/// solves → re-assemble → {T,C} solves, and returns the updated `sol`
/// (dense) plus every dense operand for the cross-check. N_θ = 1 only.
/// `c = i_r*n_z + i_z`.
#[doc(hidden)]
#[allow(clippy::too_many_arguments)]
pub fn xcheck_class_d_iterate_dense(
    g: &Grid,
    op: &GasDiffusion<'_>,
    tr: &GasTransportField,
    wqnew: f64,
    sol_rho: impl Fn(usize, usize) -> f64,
    sol_ur: impl Fn(usize, usize) -> f64,
    sol_om: impl Fn(usize, usize) -> f64,
    sol_uz: impl Fn(usize, usize) -> f64,
    sol_tt: impl Fn(usize, usize) -> f64,
    sol_cc: impl Fn(usize, usize) -> f64,
    lag_ur: impl Fn(usize, usize) -> f64,
    lag_uz: impl Fn(usize, usize) -> f64,
    dlag_of: impl Fn(usize, usize) -> [f64; NCOMP],
) -> ClassDIterateDense {
    let mut sol = GasOperands::alloc(g);
    let mut lag = GasOperands::alloc(g);
    let mut dlag = vec![vec![[0.0f64; NCOMP]; BRICK_CELLS]; g.n_bricks()];
    let mut dstage = vec![vec![[0.0f64; NCOMP]; BRICK_CELLS]; g.n_bricks()];
    let mut work = GasWork::alloc(g);
    g.for_each_active_cell(|cell| {
        let (i_r, i_z, bi, idx) = (cell.i_r, cell.i_z, cell.bi, cell.idx);
        sol.rho[bi][idx] = sol_rho(i_r, i_z);
        sol.ur[bi][idx] = sol_ur(i_r, i_z);
        sol.om[bi][idx] = sol_om(i_r, i_z);
        sol.uz[bi][idx] = sol_uz(i_r, i_z);
        sol.tt[bi][idx] = sol_tt(i_r, i_z);
        sol.cc[bi][idx] = sol_cc(i_r, i_z);
        lag.ur[bi][idx] = lag_ur(i_r, i_z);
        lag.uz[bi][idx] = lag_uz(i_r, i_z);
        dlag[bi][idx] = dlag_of(i_r, i_z);
        // ke_base at the ORIGINAL operands (derive_gas_operands' capture).
        let r = g.r_center(i_r);
        let ut = sol.om[bi][idx] * r;
        work.ke_base[bi][idx] =
            0.5 * (sol.ur[bi][idx] * sol.ur[bi][idx] + ut * ut + sol.uz[bi][idx] * sol.uz[bi][idx]);
    });
    // Capture the initial sol for the dense-input echo (the GPU gets the same).
    let sol0 = {
        let mut s = GasOperands::alloc(g);
        s.clone_from(&sol);
        s
    };
    op.fill_lag_gradients(g, &lag, &mut work)
        .expect("xcheck fill_lag_gradients");
    // assemble(sol,lag) → dstage; velocity block {Ur,Uz,Om}.
    op.assemble_rates(g, &sol, &lag, &work, tr, 0.0, &mut dstage, None)
        .expect("xcheck assemble (velocity stage)");
    for comp in [GasComp::Ur, GasComp::Uz, GasComp::Om] {
        fill_gas_rhs(
            g,
            comp,
            wqnew,
            &dstage,
            &dlag,
            &sol,
            &work.ke_base,
            &mut work.b,
        );
        op.fill_mass(g, comp, &sol.rho, tr, &mut work);
        let x = match comp {
            GasComp::Ur => &mut sol.ur,
            GasComp::Uz => &mut sol.uz,
            _ => &mut sol.om,
        };
        op.cg_solve(g, comp, tr, wqnew, x, &mut work)
            .expect("xcheck cg (velocity)");
    }
    // re-assemble at the new velocities; {T,C}.
    op.assemble_rates(g, &sol, &lag, &work, tr, 0.0, &mut dstage, None)
        .expect("xcheck assemble (scalar stage)");
    for comp in [GasComp::T, GasComp::C] {
        fill_gas_rhs(
            g,
            comp,
            wqnew,
            &dstage,
            &dlag,
            &sol,
            &work.ke_base,
            &mut work.b,
        );
        op.fill_mass(g, comp, &sol.rho, tr, &mut work);
        let x = match comp {
            GasComp::T => &mut sol.tt,
            _ => &mut sol.cc,
        };
        op.cg_solve(g, comp, tr, wqnew, x, &mut work)
            .expect("xcheck cg (scalar)");
    }
    let (n_r, n_z) = (g.spec().n_r, g.spec().n_z);
    let ncell = n_r * n_z;
    let z = || vec![0.0f64; ncell];
    let mut out = ClassDIterateDense {
        ur: z(),
        om: z(),
        uz: z(),
        tt: z(),
        cc: z(),
        rho: z(),
        lag_ur: z(),
        lag_uz: z(),
        mu: z(),
        k: z(),
        rhod: z(),
        dhdz: z(),
        cv: z(),
        gas: z(),
        dlag: vec![[0.0f64; NCOMP]; ncell],
        n_r,
        n_z,
        ..Default::default()
    };
    g.for_each_active_cell(|cell| {
        let c = cell.i_r * n_z + cell.i_z;
        let (bi, idx) = (cell.bi, cell.idx);
        out.ur[c] = sol.ur[bi][idx];
        out.om[c] = sol.om[bi][idx];
        out.uz[c] = sol.uz[bi][idx];
        out.tt[c] = sol.tt[bi][idx];
        out.cc[c] = sol.cc[bi][idx];
        out.rho[c] = sol0.rho[bi][idx];
        out.lag_ur[c] = lag.ur[bi][idx];
        out.lag_uz[c] = lag.uz[bi][idx];
        out.mu[c] = tr.mu[bi][idx];
        out.k[c] = tr.k[bi][idx];
        out.rhod[c] = tr.rho_d[bi][idx];
        out.dhdz[c] = tr.dh_dz[bi][idx];
        out.cv[c] = tr.cv[bi][idx];
        out.gas[c] = 1.0;
        out.dlag[c] = dlag[bi][idx];
        // Seed the sol echo with the INITIAL operands (what the GPU is given).
        // (ur/om/uz/tt/cc above are the FINAL — the CPU oracle to compare.)
    });
    // The initial operands the GPU iterate consumes (separate from the oracle).
    out.init_ur = g_dense(g, n_z, &sol0.ur);
    out.init_om = g_dense(g, n_z, &sol0.om);
    out.init_uz = g_dense(g, n_z, &sol0.uz);
    out.init_tt = g_dense(g, n_z, &sol0.tt);
    out.init_cc = g_dense(g, n_z, &sol0.cc);
    out
}

fn g_dense(g: &Grid, n_z: usize, f: &[Vec<f64>]) -> Vec<f64> {
    let mut v = vec![0.0f64; g.spec().n_r * n_z];
    g.for_each_active_cell(|cell| {
        v[cell.i_r * n_z + cell.i_z] = f[cell.bi][cell.idx];
    });
    v
}

/// Dense result of the S13c resident-iterate cross-check (doc-hidden). The
/// `ur..cc` fields are the CPU **oracle** (final, post-iterate); `init_*` are
/// the **initial** operands the GPU iterate consumes. `c = i_r*n_z + i_z`.
#[doc(hidden)]
#[derive(Default)]
pub struct ClassDIterateDense {
    pub ur: Vec<f64>,
    pub om: Vec<f64>,
    pub uz: Vec<f64>,
    pub tt: Vec<f64>,
    pub cc: Vec<f64>,
    pub rho: Vec<f64>,
    pub lag_ur: Vec<f64>,
    pub lag_uz: Vec<f64>,
    pub mu: Vec<f64>,
    pub k: Vec<f64>,
    pub rhod: Vec<f64>,
    pub dhdz: Vec<f64>,
    pub cv: Vec<f64>,
    pub gas: Vec<f64>,
    pub dlag: Vec<[f64; NCOMP]>,
    pub init_ur: Vec<f64>,
    pub init_om: Vec<f64>,
    pub init_uz: Vec<f64>,
    pub init_tt: Vec<f64>,
    pub init_cc: Vec<f64>,
    pub n_r: usize,
    pub n_z: usize,
}

/// The gas Picard residual: the STATE effect of the rate staleness over
/// this step (`wqnew·|Δrate|`) relative to each conserved component's own
/// magnitude on the current composed state — the three momentum
/// components share one scale (same units; the lagged cross terms couple
/// exactly the momentum block, so momentum-joint is the contraction
/// claim). A component's rate change is NEVER divided by its own possibly-
/// degenerate rate scale (a quiescent component's fp-noise rates would
/// read as divergence). NaN anywhere returns NaN (the caller's NaN-safe
/// acceptance). Serial fixed order.
fn rate_resid(
    g: &Grid,
    ids: &[FieldId; NCOMP],
    dcur: &[Vec<Cons>],
    dlag: &[Vec<Cons>],
    wqnew: f64,
) -> f64 {
    let mut s_mom = 0.0f64;
    let mut s_en = 0.0f64;
    let mut s_rc = 0.0f64;
    for bi in 0..g.n_bricks() {
        let b = g.brick(bi);
        let plane = b.n_theta() as usize * BRICK_CELLS;
        for idx in 0..plane {
            if b.mask() & (1u64 << (idx % BRICK_CELLS)) == 0 {
                continue;
            }
            for k in [I_MR, I_MT, I_MZ, I_EN, I_RC] {
                // `f64::max` DROPS NaN, so a NaN scale would silently
                // shrink to a finite one and let the acceptance pass on a
                // poisoned state — the S2 `rel_resid` finding class, and
                // the reason this is an explicit test (S3 review finding).
                let v = b.field(ids[k])[idx];
                if v.is_nan() {
                    return f64::NAN;
                }
                match k {
                    I_EN => s_en = s_en.max(v.abs()),
                    I_RC => s_rc = s_rc.max(v.abs()),
                    _ => s_mom = s_mom.max(v.abs()),
                }
            }
        }
    }
    let mut resid = 0.0f64;
    for k in 1..NCOMP {
        let scale = match k {
            I_EN => s_en,
            I_RC => s_rc,
            _ => s_mom,
        };
        let mut delta = 0.0f64;
        for (a, b) in dcur.iter().zip(dlag) {
            for (ca, cb) in a.iter().zip(b) {
                let (x, y) = (ca[k], cb[k]);
                if x.is_nan() || y.is_nan() {
                    return f64::NAN;
                }
                delta = delta.max((x - y).abs());
            }
        }
        resid = resid.max(wqnew * delta / scale.max(f64::MIN_POSITIVE));
    }
    resid
}

// --- Small field utilities (random access at coupler rate) --------------------

fn cell_value(g: &Grid, f: FieldId, i_r: usize, i_z: usize) -> f64 {
    let bi = g
        .brick_index_by_coords((i_r / BRICK) as u32, (i_z / BRICK) as u32)
        .expect("cell in an allocated brick");
    let b = g.brick(bi);
    b.field(f)[b.cell_index(0, (i_r % BRICK) * BRICK + i_z % BRICK)]
}

fn cell_add(g: &mut Grid, f: FieldId, i_r: usize, i_z: usize, dv: f64) {
    let bi = g
        .brick_index_by_coords((i_r / BRICK) as u32, (i_z / BRICK) as u32)
        .expect("cell in an allocated brick");
    let idx = {
        let b = g.brick(bi);
        b.cell_index(0, (i_r % BRICK) * BRICK + i_z % BRICK)
    };
    g.brick_field_mut(bi, f)[idx] += dv;
}
