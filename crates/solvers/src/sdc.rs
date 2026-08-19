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
    Cons, EosLaw, Euler, EulerFields, EulerWorkspace, FlowError, FlowLedger, I_EN, NCOMP, Prim,
    srd_neighborhood,
};
use crate::wall_heat::{NearWallGas, WallHeatError, WallLaw};
use crucible_grid::{BRICK, BRICK_CELLS, FaceDir, FieldId, Grid, InterfaceFace, tree_combine};

// --- Named constants (COUP-3 §3.7 / COUP-2 §3.1.1) ---------------------------

/// Fixed SDC correction-sweep count after the IMEX-Euler predictor
/// (COUP-3 §3.1: "2–3 sweeps to 2nd order" — predictor + 2 corrections).
pub const N_SDC_CORRECTIONS: usize = 2;

/// Fixed Picard sweep count of the Robin-Robin wall-exchange solve inside
/// each SDC sweep's class-`D` solve (COUP-2 §3.5).
pub const N_ROBIN_SWEEPS: usize = 3;

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
pub const EPS_ROBIN_RESID: f64 = 1e-6;

/// Clamped Aitken relaxation bounds of the Picard sweeps (COUP-3 §3.5's
/// deterministic clamp discipline, applied to the exchange iteration).
pub const ROBIN_OMEGA_MIN: f64 = 0.1;
pub const ROBIN_OMEGA_MAX: f64 = 2.0;

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
    /// Last class-`D` CG solve's iteration count and terminal residual.
    pub cg_iters: usize,
    pub cg_resid: f64,
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
    /// debited exactly −this) and composed exterior (coolant) heat [J].
    pub applied_exchange_j: f64,
    pub applied_exterior_j: f64,
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

/// The one deterministic integrator (persistent workspaces; one instance
/// per grid). Class arguments arrive per step so callers may evolve
/// boundary schedules between steps (COUP-7 declared schedules).
pub struct Sdc {
    ws: Option<EulerWorkspace>,
    rate_e0: Vec<Vec<Cons>>,
    sb: Option<SolidBufs>,
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
        if self.ws.is_none() {
            self.ws = Some(flow.op.workspace(g)?);
            let plane = g.brick(0).n_theta() as usize * BRICK_CELLS;
            self.rate_e0 = vec![vec![[0.0; NCOMP]; plane]; g.n_bricks()];
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

    /// Pure-diffusion convenience (the Goal-A studies): the same step with
    /// only class `D` scheduled.
    pub fn step_diffusion(
        &mut self,
        g: &mut Grid,
        diffusion: &DiffusionClass<'_, '_>,
        t: f64,
        dt: f64,
    ) -> Result<StepReport, SdcError> {
        self.step::<crate::euler::GammaLaw>(g, None, Some(diffusion), None, t, dt)
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
        self.step(g, Some(flow), None, None, t, dt)
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
            self.step_flow(g, flow, t, dt)?;
            t += dt;
            steps += 1;
        }
        Ok(steps)
    }

    /// One SDC-IMEX step of size `dt` at time `t` (module doc). The
    /// schedule is the fixed source order: class `A` evaluation, then the
    /// class-`D` solve with the Robin-Robin exchange inside, per sweep;
    /// the COUP-2 audit closes the step.
    pub fn step<E: EosLaw + Sync>(
        &mut self,
        g: &mut Grid,
        flow: Option<&FlowClass<'_, '_, E>>,
        diffusion: Option<&DiffusionClass<'_, '_>>,
        exchange: Option<&ExchangeClass<'_>>,
        t: f64,
        dt: f64,
    ) -> Result<StepReport, SdcError> {
        // --- Schedule validation (loud, META-1 P6) -----------------------
        if flow.is_none() && diffusion.is_none() {
            return Err(SdcError::Config("no operator class scheduled"));
        }
        if exchange.is_some() && (flow.is_none() || diffusion.is_none()) {
            return Err(SdcError::Config(
                "wall exchange needs both the hyperbolic and diffusion classes",
            ));
        }
        if let (Some(dc), Some(_)) = (diffusion, flow) {
            if dc.op.domain != Domain::Solid {
                return Err(SdcError::Config(
                    "a gas-domain diffusion class alongside the flow class is the S3 \
                     viscous wave; this schedule couples flow to SOLID conduction only",
                ));
            }
            if g.brick(0).n_theta() != 1 {
                return Err(SdcError::Config(
                    "coupled flow+conduction at N_θ > 1 arrives with the 3-D wave; \
                     refusing rather than guessing",
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

        // --- The fixed sweeps --------------------------------------------
        let mut l_last = l0;
        let mut hl_prev = hl0;
        let mut hl_last = HeatLedger::default();
        let mut q_prev = q0.clone();
        let mut debit_applied = 0.0f64; // κV-weighted gas energy debited (J)
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

            // Robin-Picard loop (single pass when no exchange is scheduled).
            let n_picard = if exchange.is_some() {
                N_ROBIN_SWEEPS
            } else {
                1
            };
            let mut q_trial = q_prev.clone();
            let mut omega = 1.0f64;
            let mut resid_hist: Option<Vec<f64>> = None; // previous residual vector
            let mut robin_resid = 0.0f64;
            let mut q_new: Vec<f64> = Vec::new();
            let mut hl_tmp = HeatLedger::default();

            for ps in 0..n_picard {
                if let Some(fc) = flow {
                    self.compose_gas(
                        g, fc, we0, we1, exchange, wq0, &q0, wqprev, &q_prev, wqnew, &q_trial,
                    );
                    let ws = self.ws.as_ref().expect("ensured");
                    fc.op.apply_srd(g, fc.fields, ws);
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

            // Accept the sweep: final gas composition with the ACCEPTED
            // exchange heats (exactly what the solid solve received).
            if let Some(fc) = flow {
                self.compose_gas(
                    g, fc, we0, we1, exchange, wq0, &q0, wqprev, &q_prev, wqnew, &q_trial,
                );
                let ws = self.ws.as_ref().expect("ensured");
                fc.op.apply_srd(g, fc.fields, ws);
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
                    ex_report.q_w = q_trial.clone();
                    ex_report.h = ex_ops.iter().map(|(h, _)| *h).collect();
                    ex_report.t_aw = ex_ops.iter().map(|(_, ta)| *ta).collect();
                    ex_report.jacket_w = q_trial.iter().sum();
                    ex_report.applied_exchange_j = applied_exchange_j;
                    ex_report.applied_exterior_j = applied_exterior_j;
                }
                hl_last = hl_tmp;
            } else if diffusion.is_some() {
                // Roll the F_I record: heat_prev ← this sweep's accepted
                // assembly; q_prev/hl_prev likewise.
                let sb = self.sb.as_mut().expect("ensured");
                std::mem::swap(&mut sb.heat_prev, &mut sb.heat_cur);
                hl_prev = hl_tmp;
                q_prev = q_trial.clone();
            }
        }

        // --- The COUP-2 audit (module doc) --------------------------------
        let stored_after = self.audit_stored(g, flow, diffusion);
        self.audit_check(
            g,
            flow.is_some(),
            diffusion.is_some(),
            dt,
            &stored_before,
            &stored_after,
            &l0,
            &l_last,
            &hl0,
            &hl_prev,
            &hl_last,
            debit_applied,
            &mut report,
        )?;
        if exchange.is_some() {
            report.exchange = Some(ex_report);
        }
        Ok(report)
    }

    /// Gas node-1 composition: `U = U⁰ + we0·F_E⁰ + we1·F_E(cur)` plus the
    /// exchange debit `−(wq0·q0 + wqprev·q_prev + wqnew·q_trial)` per patch
    /// distributed uniformly over its SRD debit set. Fixed order; serial.
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
            let (n_r, n_z) = patch.n_hat;
            let v_n = w[1] * n_r + w[3] * n_z;
            let u_t = ((w[1] * w[1] + w[3] * w[3] - v_n * v_n).max(0.0) + w[2] * w[2]).sqrt();
            let y = 0.5 * (n_r.abs() * dr + n_z.abs() * dz);
            let gas = NearWallGas {
                rho: w[0],
                u_t,
                temperature,
                y,
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
        // NaN-safe acceptance (a NaN residual is a violation).
        if b_norm2 > 0.0 && (resid.is_nan() || resid > EPS_CG_RESID) {
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
    fn audit_stored<E: EosLaw>(
        &self,
        g: &Grid,
        flow: Option<&FlowClass<'_, '_, E>>,
        diffusion: Option<&DiffusionClass<'_, '_>>,
    ) -> ([f64; NCOMP], f64) {
        let mut gas = [0.0f64; NCOMP];
        if let Some(fc) = flow {
            for (k, &id) in fc.fields.ids().iter().enumerate() {
                gas[k] = g.reduce_kappa_volume_weighted(id);
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
        (gas, solid)
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
        dt: f64,
        before: &([f64; NCOMP], f64),
        after: &([f64; NCOMP], f64),
        l0: &FlowLedger,
        l_last: &FlowLedger,
        hl0: &HeatLedger,
        hl_prev: &HeatLedger,
        hl_last: &HeatLedger,
        debit_applied: f64,
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
        ];
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
                let mut stored_scale = before.0[k].abs();
                if k == I_EN {
                    // The combined-energy row: the exchange pair cancels
                    // between the gas debit and the solid's received heats.
                    applied += debit_applied;
                    if has_diffusion {
                        delta += after.1 - before.1;
                        applied += solid_applied;
                        gross += debit_applied.abs() + solid_gross;
                        stored_scale += before.1.abs();
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
            let delta = after.1 - before.1;
            let s = before.1.abs() + solid_gross;
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

/// Max relative residual of the exchange-heat vector, floored so a
/// zero-heat patch (cold start) cannot divide by zero.
fn rel_resid(q_new: &[f64], resid: &[f64]) -> f64 {
    let scale = q_new
        .iter()
        .fold(0.0f64, |m, q| m.max(q.abs()))
        .max(f64::MIN_POSITIVE);
    resid.iter().fold(0.0f64, |m, r| m.max(r.abs())) / scale
}

/// Per-patch heat sums from an assembly's recorded per-face exchange
/// heats (fixed patch/face order — the visit order is keyed, not assumed).
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
