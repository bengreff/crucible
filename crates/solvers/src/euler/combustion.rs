//! SOLV-4 §3.6 — the **burn-progress source operator** (the doc's `c` is the
//! code's `b`/`I_RB`). It accumulates the one continuous rate law (SOLV-4.4)
//! into the class-`A` explicit rate at every SDC node, riding inside
//! [`super::Euler::eval_rhs`] so the SDC step composes and audits it exactly
//! like the flow's own geometric source — no new SDC-step plumbing.
//!
//! **The reaction-diffusion (thickened-flame) form (SOLV-4.4):**
//! `∂(ρb)/∂t|_source = ∇·(ρ D_c ∇b) + ρ_u·K·b(1−b)(b−a) + ρ·(1−b)/τ_ign`, with
//! the **bistable (Nagumo) matched coefficients** `D_c = w·S_T/(1−2a)`,
//! `K = 2·S_T/((1−2a)·w)`, `w = Θ·Δ` (`S_T = S_L·wrinkling`, `Δ =` the local
//! cell length, `a =` the ignition threshold). The Nagumo pushed front then
//! travels at exactly `√(D_c·K/2)·(1−2a) = S_T` with width `√(2D_c/K) = w =
//! Θ·Δ` — a fixed `Θ` cells at **every** resolution, so the front *speed* is
//! closure-set (the `flame_1d` gate) while the *width* is the declared
//! front-carrier device. The bistable (not monostable-KPP) form makes unburnt
//! `b = 0` metastable — a **pushed** front (speed set by the core, robust and
//! grid-independent), and no spurious self-ignition of the unburnt domain.
//! Only `ρb` is sourced — the heat release is implicit in the blended EOS
//! (SOLV-4 §3.6): as `b` rises the same conserved enthalpy reads hotter off
//! the burnt branch.
//!
//! **Integration scope (S7 — the class split of record):** the propagation
//! reaction and the front-thickening diffusion stay **explicit class A**
//! (non-stiff by construction: the matched `D_c`/`K` clock is `Δ/S_T`, far
//! above the acoustic Δt since `S_T ≪ a`) behind the **loud positivity
//! guard** — `b` leaving `[0,1]` beyond [`EPS_BURN_BOUND`] is a halt (under-
//! resolution), never a clamp (META-1 P6). The **auto-ignition term is the
//! cell-local implicit class-`R` occupant** (COUP-3 §3.3, landed at S7 —
//! the S6 explicit tier retired): at chamber pressure and `T_u ≳ 1300 K`
//! `τ_ign` collapses below the acoustic Δt (`dt/τ → O(1)`), where an
//! explicit advance overshoots the guard; the [`Self::implicit_auto_update`]
//! node solve is unconditionally stable, exact in the stiff limit (the
//! guarded law parks at its `1 − BURN_COMPLETE` fixed point), and reduces to
//! the explicit value when mild — ONE treatment, no regime branch (Rule 12).
//!
//! **The cold-side non-reactive floor (S7, OFFL-3 §3.3):** where the cell's
//! `p` or `T_u` sits below the ignition surface's own declared envelope
//! floor, both rate terms read **zero by declaration** — the cold analogue
//! of the [`BURN_COMPLETE`] domain guard. Near-vacuum fill (~10² Pa) and
//! sub-cryo corners are declared no-burn (cm-scale flames genuinely cannot
//! propagate there), never a refusal mid-march and never a clamp of a
//! mid-range value; the surface's cold rows (ignition v0.3.0) carry the real
//! `S_L` down to the unburnt branch's own validity floor, so the front's
//! consumption of injector-cold reactants is data, not this guard.
//!
//! The serial cell loop is bit-deterministic at any thread count. **S8:**
//! the operator runs per θ-plane at the grid's (uniform) N_θ — the
//! front-thickening diffusion gains its θ-direction faces (periodic ring
//! stencil, within-brick; face `ρD_c` = the two-cell arithmetic mean, like
//! the r/z faces), so a spark kernel that is a point in θ propagates
//! azimuthally. Mixed per-brick N_θ refuses upstream (`Euler::validate` —
//! the ring-interface c-diffusion exchange rides plan S11). **The S_T-CFL
//! (SOLV-4 0.4.6):** [`Combustion::front_carrier_signal`] is the carrier's
//! explicit stability rate `σ_front = 2·D_c·Σ_d 1/Δ_d² +
//! C_NAGUMO_SLOPE·(ρ_u/ρ)·K` — a member of the Δt reduction (COUP-3 §3.4)
//! wherever the rate law is live, with the [`S_T_MACH_LIMIT`] model-form
//! scale-separation refusal at the sonic end (S_T vs the LOCAL sound speed
//! — deliberately velocities, not mesh rates), so a declared wrinkling > 1
//! marches honestly at the carrier's own Δt instead of tripping the
//! positivity guard.

use super::{BurnBlendEos, Cons, EosLaw, FlowError, I_RB, I_RC, I_RHO, NCOMP, Prim, Scratch};
use crucible_grid::{BRICK, BRICK_CELLS, FaceDir, FieldId, Grid};
use crucible_tables::{BoundColumn, Table};

/// The declared front width parameter in cells `Θ` (SOLV-4.4): the Nagumo
/// wave scale `w = Θ·Δ`. ~1.5 gives a ~9-cell 2–98% front — smooth enough for
/// the central-difference Laplacian while staying local; a front-carrier
/// device, not a flame-thickness claim.
pub const THETA_CELLS: f64 = 1.5;

/// A cell with `c` within this of 1 is treated as **fully burnt**: the
/// reaction has vanished there and the unburnt-reactant ignition closures no
/// longer describe its (post-flame) gas, so the operator skips their queries
/// (a domain guard on the burnt fixed point — see the `accumulate` note).
/// **One owner (SOLV-4 0.4.3):** this constant also derives the blend's
/// pure-burnt threshold (`EPS_B_PURE_BURNT = 2·BURN_COMPLETE` — sized to
/// *strictly contain* the pinning attractor — the S6 explicit tier
/// asymptoted just *under* `1 − BURN_COMPLETE`; the S7 class-`R` solve
/// parks exactly AT it (SOLV-4 0.4.4); an equal threshold is a
/// knife edge) and is the diagnostics' guard
/// ([`reacting_measure`]/[`consumption_rate`]). A tighter blend threshold
/// would leave burnt cells interrogating the unburnt branch at ~10⁻³ weight
/// forever, refusing on hot states that branch cannot describe.
pub const BURN_COMPLETE: f64 = 1.0e-3;

/// The declared **sub-cell ignition threshold** `a ∈ (0, ½)` of the bistable
/// (Nagumo) propagation reaction `c(1−c)(c−a)` (SOLV-4.4). It makes unburnt
/// `c = 0` **metastable** (`c < a` decays — reactants do not burn without a
/// trigger) and the front a robust **pushed** front. 0.2 keeps a clear
/// metastable well while leaving the speed factor `(1−2a) = 0.6` healthy.
pub const IGN_THRESHOLD_A: f64 = 0.2;

/// Max slope `|d/dc[c(1−c)(c−a)]|` of the bistable shape on [0, 1] at
/// `a = IGN_THRESHOLD_A` (attained at c = 1: `|−3 + 2(1+a) − a| = 1 − a·…`
/// = 0.8 for a = 0.2) — the propagation reaction's explicit stability
/// scale in [`Combustion::front_carrier_signal`] (SOLV-4 §3.6, S8). Tied
/// to [`IGN_THRESHOLD_A`]: recompute if that constant ever moves.
pub const C_NAGUMO_SLOPE: f64 = 0.8;

/// SOLV-4 §3.6 (S8, restated on MODEL-FORM quantities at 0.4.7) — the
/// scale-separation guard: a cell whose `S_T` exceeds this fraction of its
/// own sound speed refuses loudly (`FlowError::FrontCarrierScaleSeparation`)
/// — the quasi-isobaric flamelet premises are broken (fast-deflagration/
/// DDT class, outside the declared model form). Deliberately a ratio of
/// VELOCITIES, not of mesh signal rates: the S8 review measured the
/// rate-ratio form resolution-dependent (the θ-arc's 1/arc² carrier-
/// diffusion rate outgrows the 1/arc acoustic rate, tightening the
/// threshold linearly in N_θ at the innermost rings until it refuses mild
/// flames the certified tier marches). The S7-certified laminar march's
/// worst crossover-band cells sit near S_T/c ≈ 0.27, comfortably inside;
/// 2/3 refuses only the genuinely sonic-class front.
pub const S_T_MACH_LIMIT: f64 = 2.0 / 3.0;

/// The matched Nagumo `(D_c, K)` for a front of speed `S_T` and width `Θ·Δ`
/// (SOLV-4.4): `D_c = w·S_T/(1−2a)`, `K = 2·S_T/((1−2a)·w)`, `w = Θ·Δ`. Then
/// `√(D_c·K/2)·(1−2a) = S_T` (speed) and `√(2D_c/K) = w` (width), both scaling
/// with `Δ` ⇒ grid-independent. `S_T = 0` ⇒ `(0, 0)` (extinct: no propagation).
#[inline]
fn front_coeffs(s_t: f64, delta: f64, theta: f64) -> (f64, f64) {
    if s_t <= 0.0 {
        return (0.0, 0.0);
    }
    let w = theta * delta;
    let f = 1.0 - 2.0 * IGN_THRESHOLD_A;
    (w * s_t / f, 2.0 * s_t / (f * w))
}

/// The bistable propagation reaction rate into `ρc` [kg/m³/s]: `ρ_u·K·b(1−b)(b−a)`.
/// Negative for `b < a` (the metastable pull back to unburnt).
#[inline]
fn propagation_rate(rho_u: f64, k: f64, b_c: f64) -> f64 {
    rho_u * k * b_c * (1.0 - b_c) * (b_c - IGN_THRESHOLD_A)
}

/// Loud-guard tolerance on the burn progress `b ∈ [0,1]`: a value beyond
/// `[−EPS_BURN_BOUND, 1+EPS_BURN_BOUND]` after a step is under-resolved
/// ignition and halts (never a silent clamp). Widened 1e-6 → 1e-4 at S7:
/// an engine-scale ignition BLAST advects `ρc` and `ρ` through the same
/// captured front with independently flux-limited reconstructions, and the
/// ratio `b = ρc/ρ` legitimately excurses at the 1e-5 scale there (the
/// scalar-consistency property of conservative component-wise advection —
/// a numerical fact about the shock, not an ignition pathology). 1e-4
/// stays an order under the partition's percent-class weights and three under any physical
/// excursion, so a genuine source-side blowup still halts loudly.
pub const EPS_BURN_BOUND: f64 = 1.0e-4;

/// The **reacting-measure ignition floor** (SOLV-4 §3.6, the single owner of
/// the COUP-4 `NEVER_IGNITED`/`FLAMEOUT` thresholds): a self-sustaining burn
/// has [`reacting_measure`] `R = ∫ b(1−b)·rate dV` above this, a dead one
/// below it. `NEVER_IGNITED` = the igniter schedule exhausted with `R` never
/// exceeding it; `FLAMEOUT` = `R` collapsing below it after having exceeded.
/// Units are kg/s (a mass-burning rate). **Scaling:** at the RL10 tier (S7)
/// the run scales this by its injected mass rate (a fixed fraction — a flame
/// that consumes < `EPS_IGNITED` is not established); at the mini-sim tier the
/// absolute floor below distinguishes a lit kernel from a dying one.
pub const EPS_IGNITED: f64 = 1.0e-6;

/// The OFFL-3 §3.3 ignition closure surface, bound once: `S_L(p, T_u, Z)` and
/// `τ_ign(p, T_u, Z)` (SOLV-4 §3.6). Extinction is the columns' own values
/// (`S_L → 0`, `τ_ign → τ_max`), never a runtime branch.
#[derive(Debug, Clone)]
pub struct IgnitionColumns<'t> {
    flame_speed: BoundColumn<'t>,
    ignition_delay: BoundColumn<'t>,
    /// The surface's own declared envelope floors in `p` [Pa] and `T_u` [K],
    /// cached at bind — the single owner of the cold-side non-reactive floor
    /// (module header): below either, the rate law is zero by declaration.
    p_floor: f64,
    tu_floor: f64,
}

impl<'t> IgnitionColumns<'t> {
    /// Bind to a loaded `(p, T_u, Z)` ignition surface; refuses (never
    /// guesses) if it is not that schema or a column is missing/relabeled.
    pub fn bind(table: &'t Table) -> Result<Self, String> {
        let flame_speed = table
            .bind("laminar_flame_speed", "m/s")
            .map_err(|e| format!("ignition surface: {e}"))?;
        let ignition_delay = table
            .bind("ignition_delay", "s")
            .map_err(|e| format!("ignition surface: {e}"))?;
        let expected = ["p", "T_u", "Z"];
        for (i, want) in expected.iter().enumerate() {
            let got = flame_speed.axis_name(i);
            if got != *want {
                return Err(format!(
                    "ignition surface axis {i} is {got:?}, expected {want:?} — \
                     not a (p, T_u, Z) closure surface (OFFL-3 §3.3); refusing the bind"
                ));
            }
        }
        Ok(Self {
            p_floor: flame_speed.axis_envelope(0).0,
            tu_floor: flame_speed.axis_envelope(1).0,
            flame_speed,
            ignition_delay,
        })
    }

    /// The cold-side **non-reactive floor** (module header, OFFL-3 §3.3):
    /// `true` when the state sits below the surface's own declared envelope
    /// floor in `p` or `T_u` — declared no-burn, both rate terms zero. The
    /// hot/rich edges stay hard refusals (the envelope-consistency contract
    /// guarantees a legal blend state is covered there; an excursion is a
    /// real defect, never declared away).
    pub fn non_reactive_floor(&self, p: f64, t_u: f64) -> bool {
        p < self.p_floor || t_u < self.tu_floor
    }

    fn s_l(&self, p: f64, t_u: f64, z: f64) -> Result<f64, &'static str> {
        self.flame_speed
            .interpolate(&[p, t_u, z])
            .map_err(|_| "flame-speed query outside the ignition-surface envelope")
    }

    fn tau_ign(&self, p: f64, t_u: f64, z: f64) -> Result<f64, &'static str> {
        self.ignition_delay
            .interpolate(&[p, t_u, z])
            .map_err(|_| "ignition-delay query outside the ignition-surface envelope")
    }

    /// Public `S_L(p, T_u, Z)` [m/s] — the VAL-2 anchor / diagnostic accessor.
    pub fn laminar_flame_speed(&self, p: f64, t_u: f64, z: f64) -> Result<f64, &'static str> {
        self.s_l(p, t_u, z)
    }

    /// Public `τ_ign(p, T_u, Z)` [s] — the VAL-2 anchor / diagnostic accessor.
    pub fn induction_time(&self, p: f64, t_u: f64, z: f64) -> Result<f64, &'static str> {
        self.tau_ign(p, t_u, z)
    }
}

/// The burn-progress source operator (SOLV-4 §3.6). Held by `Euler` as an
/// optional occupant; a shifting-only run schedules none and the `ρb` slot
/// stays inert.
pub struct Combustion<'a> {
    /// The blended EOS — the `T_u`/`ρ_u` provider (the rate-law coordinate).
    pub blend: &'a BurnBlendEos<'a>,
    /// The `(p, T_u, Z)` ignition closures.
    pub ignition: IgnitionColumns<'a>,
    /// Turbulent wrinkling `S_T/S_L` (`turbulent-flame-speed`, banded); 1.0 =
    /// laminar (the mini-sim tier). Pure data.
    pub wrinkling: f64,
    /// Front thickness in cells `Θ` (see [`THETA_CELLS`]).
    pub theta: f64,
}

impl Combustion<'_> {
    /// Accumulate SOLV-4.4 into the class-`A` rate (`s.rate[·][I_RB]`) + the
    /// COUP-2 source ledger, from the sweep-current primitive cache (`s.prim`)
    /// and grid state (for the neighbour `b`-gradient / diffusion). Serial +
    /// deterministic. Only `ρb` is sourced.
    pub(crate) fn accumulate(
        &self,
        g: &Grid,
        ids: &[FieldId; NCOMP],
        s: &mut Scratch,
    ) -> Result<(), FlowError> {
        // Uniform N_θ is guaranteed upstream (`Euler::validate`, S8: mixed
        // N_θ + combustion refuses with the S11 owner named).
        let nt = g.brick(0).n_theta();
        debug_assert!(g.bricks().iter().all(|b| b.n_theta() == nt));
        let spec = g.spec();
        let (dr, dz) = (spec.dr, spec.dz);
        let delta = (dr * dz).sqrt(); // local cell length Δ (isotropic on a square grid)
        let (rho_id, rhob_id) = (ids[I_RHO], ids[I_RB]);

        // b at a global cell + θ-plane (ρb/ρ); returns None outside the
        // domain or on a non-gas cell (a wall face ⇒ zero-flux for b,
        // SOLV-4 §3.6).
        let b_at = |i_r: usize, i_z: usize, j: u32| -> Option<f64> {
            if i_r >= spec.n_r
                || i_z >= spec.n_z
                || !g.is_active(i_r, i_z)
                || g.kappa(i_r, i_z) <= 0.0
            {
                return None;
            }
            let bi = g.brick_index_by_coords((i_r / BRICK) as u32, (i_z / BRICK) as u32)?;
            let b = g.brick(bi);
            let cell = b.cell_index(j, (i_r % BRICK) * BRICK + i_z % BRICK);
            Some(b.field(rhob_id)[cell] / b.field(rho_id)[cell])
        };

        let Scratch {
            prim, rate, ledger, ..
        } = s;
        for bi in 0..g.n_bricks() {
            let brick = g.brick(bi);
            let mask = brick.mask();
            for local in 0..BRICK_CELLS {
                if mask & (1u64 << local) == 0 {
                    continue;
                }
                let (i_r, i_z) = brick.global_rz(local);
                for j in 0..nt {
                    let idx = j as usize * BRICK_CELLS + local;
                    let w: &Prim = &prim[bi][idx];
                    let (rho, p) = (w[I_RHO], w[4]);
                    let (z, b) = (w[I_RC], w[I_RB]);
                    // Loud positivity guard (META-1 P6) — never a clamp.
                    if !(-EPS_BURN_BOUND..=1.0 + EPS_BURN_BOUND).contains(&b) {
                        return Err(FlowError::NonPhysicalState {
                            i_r,
                            i_z,
                            i_theta: j,
                            what: "burn progress left [0,1] — under-resolved ignition \
                               (refine Δt or gentle the igniter ramp)",
                        });
                    }
                    let b_c = b.clamp(0.0, 1.0);
                    // The ignition closures `S_L`/`τ_ign` are unburnt-REACTANT
                    // properties (the OFFL-3 surface's `T_u` coordinate); a
                    // fully-burnt cell (`c → 1`) is a stable reaction fixed point
                    // (`c(1−c) → 0`, `1−c → 0`) and its post-flame gas is not a
                    // reactant, so the unburnt branch's hot `T_u` extrapolation
                    // legitimately leaves the ignition envelope. Skip its queries
                    // (no reaction, no front-carrier) rather than fail on a
                    // quantity that does not apply — a domain guard on the burnt
                    // fixed point, not a physics regime branch (cf. the blend's
                    // pure-limit skip).
                    // The auto-ignition term is the class-`R` implicit occupant
                    // (module header, S7) — accumulated by the SDC step's
                    // reaction class, NOT here; this class-A pass carries the
                    // propagation reaction + front diffusion only.
                    let (react, rho_dc_here) =
                        if b_c >= 1.0 - BURN_COMPLETE || self.blend.below_unburnt_floor(w) {
                            // Burnt fixed point, or colder than any representable
                            // reactant (the floor's second face) — non-reactive.
                            (0.0, 0.0)
                        } else {
                            let t_u = self
                                .blend
                                .unburnt_temperature(w)
                                .map_err(|what| nonphys(i_r, i_z, j, what))?;
                            if self.ignition.non_reactive_floor(p, t_u) {
                                // Cold-side declared no-burn floor (module header).
                                (0.0, 0.0)
                            } else {
                                let rho_u = self
                                    .blend
                                    .unburnt_density(w)
                                    .map_err(|what| nonphys(i_r, i_z, j, what))?;
                                let s_l = self
                                    .ignition
                                    .s_l(p, t_u, z)
                                    .map_err(|what| nonphys(i_r, i_z, j, what))?;
                                let s_t = (s_l * self.wrinkling).max(0.0);
                                let (d_c, k) = front_coeffs(s_t, delta, self.theta);
                                (propagation_rate(rho_u, k, b_c), rho * d_c)
                            }
                        };

                    // Front-thickening diffusion ∇·(ρ D_c ∇b) — symmetric two-point
                    // face fluxes so interior faces telescope (conservative). Each
                    // face's ρD_c is the two-cell arithmetic mean; a wall/edge face
                    // is zero-flux (skipped).
                    let vol = g.cell_volume(i_r, nt);
                    // Per-sector κ (S11): on a revolved cut wall every sector's
                    // κ equals plane 0's, so this is bitwise the S9 form at
                    // N_θ = 1 and on box worlds; it becomes load-bearing only
                    // when the θ-varying wave lands a genuinely per-sector κ.
                    let kv = brick.kappa_cell(j, local) * vol;
                    let mut diff = 0.0;
                    let faces = [
                        (
                            FaceDir::RMinus,
                            i_r.checked_sub(1).map(|r| (r, i_z)),
                            g.face_area_r(i_r, false, nt),
                            dr,
                        ),
                        (
                            FaceDir::RPlus,
                            Some((i_r + 1, i_z)),
                            g.face_area_r(i_r, true, nt),
                            dr,
                        ),
                        (
                            FaceDir::ZMinus,
                            i_z.checked_sub(1).map(|zz| (i_r, zz)),
                            g.face_area_z(i_r, nt),
                            dz,
                        ),
                        (
                            FaceDir::ZPlus,
                            Some((i_r, i_z + 1)),
                            g.face_area_z(i_r, nt),
                            dz,
                        ),
                    ];
                    for (dir, nbr, area, d) in faces {
                        // Per-sector meridional aperture (S11): θ-uniform on a
                        // revolved wall (bitwise the S9 `g.aperture` view at
                        // N_θ = 1 / box worlds), per-sector when it matters.
                        let ap = g.aperture_at(i_r, j, i_z, dir);
                        if ap <= 0.0 {
                            continue;
                        }
                        let Some((nr, nz)) = nbr else { continue };
                        let Some(b_nbr) = b_at(nr, nz, j) else {
                            continue;
                        };
                        // Neighbour ρD_c (its own S_T); the face value is the mean.
                        let bi_n = g
                            .brick_index_by_coords((nr / BRICK) as u32, (nz / BRICK) as u32)
                            .expect("active neighbour's brick");
                        let idx_n = j as usize * BRICK_CELLS + (nr % BRICK) * BRICK + nz % BRICK;
                        let rho_dc_nbr = match self.rho_dc_of(&prim[bi_n][idx_n], delta) {
                            Ok(v) => v,
                            Err(what) => return Err(nonphys(nr, nz, j, what)),
                        };
                        let rho_dc_face = 0.5 * (rho_dc_here + rho_dc_nbr);
                        diff += ap * area * rho_dc_face * (b_nbr - b) / d;
                    }
                    // θ-direction faces (S8): the periodic ring stencil, within
                    // this brick at this local cell. Skipped at N_θ = 1 (the two
                    // faces are the same face — zero net flux, zero cost; the
                    // pre-S8 arithmetic is untouched bitwise). The θ-face
                    // aperture (S11) weights the flux — a revolved cut cell's
                    // constant-θ plane is partly blocked, and full dr·dz would
                    // over-diffuse the front across a wall-clipped sector. On
                    // box worlds the aperture is 1.0 (×1 exact ⇒ the S8
                    // arithmetic is bitwise untouched).
                    if nt > 1 {
                        let a_th = g.face_area_theta();
                        let arc = g.r_center(i_r) * (std::f64::consts::TAU / f64::from(nt));
                        for (jn, dir) in [
                            ((j + nt - 1) % nt, FaceDir::ThetaMinus),
                            ((j + 1) % nt, FaceDir::ThetaPlus),
                        ] {
                            let ap_th = g.aperture_at(i_r, j, i_z, dir);
                            if ap_th <= 0.0 {
                                continue;
                            }
                            let idx_n = jn as usize * BRICK_CELLS + local;
                            let b_nbr = brick.field(rhob_id)[idx_n] / brick.field(rho_id)[idx_n];
                            let rho_dc_nbr = match self.rho_dc_of(&prim[bi][idx_n], delta) {
                                Ok(v) => v,
                                Err(what) => return Err(nonphys(i_r, i_z, jn, what)),
                            };
                            let rho_dc_face = 0.5 * (rho_dc_here + rho_dc_nbr);
                            diff += ap_th * a_th * rho_dc_face * (b_nbr - b) / arc;
                        }
                    }
                    let diff_div = if kv > 0.0 { diff / kv } else { 0.0 };

                    let src = react + diff_div;
                    rate[bi][idx][I_RB] += src;
                    // Reaction + diffusion both go to the source ledger; the
                    // diffusion telescopes to the domain boundary (zero-flux ⇒ 0)
                    // in the global sum to within the ulp-level asymmetry of the
                    // two sides' b reconstruction (prim-cache ρb·(1/ρ) vs direct
                    // ρb/ρ), inside TOL_AUDIT — the audit balances against the
                    // stored Δ(ρb) (COUP-2 §3.1.1).
                    ledger.src_net[I_RB] += kv * src;
                    ledger.src_abs[I_RB] += (kv * src).abs();
                }
            }
        }
        Ok(())
    }

    /// A cell's `ρ·D_c` from its primitive — the face-mean diffusion
    /// coefficient's per-side operand (its own `ρ`, `S_T` from its `T_u`,
    /// and the matched `D_c`); 0 for non-reactive cells (burnt fixed point,
    /// sub-floor, cold-side declared no-burn).
    fn rho_dc_of(&self, w: &Prim, delta: f64) -> Result<f64, &'static str> {
        let (rho, p, z, b) = (w[I_RHO], w[4], w[I_RC], w[I_RB]);
        if b.clamp(0.0, 1.0) >= 1.0 - BURN_COMPLETE || self.blend.below_unburnt_floor(w) {
            return Ok(0.0); // burnt or sub-floor neighbour: no front-carrier
        }
        let t_u = self.blend.unburnt_temperature(w)?;
        if self.ignition.non_reactive_floor(p, t_u) {
            return Ok(0.0); // cold-side declared no-burn floor: no front-carrier
        }
        let s_l = self.ignition.s_l(p, t_u, z)?;
        let s_t = (s_l * self.wrinkling).max(0.0);
        let (d_c, _k) = front_coeffs(s_t, delta, self.theta);
        Ok(rho * d_c)
    }

    /// SOLV-4 §3.6 (S8) — the front-carrier's explicit stability rate
    /// `σ_front = 2·D_c·Σ_d 1/Δ_d² + C_NAGUMO_SLOPE·(ρ_u/ρ)·K` [1/s], the
    /// Δt-reduction member (COUP-3 §3.4), returned WITH the cell's `S_T`
    /// (the [`S_T_MACH_LIMIT`] guard's model-form operand — the caller
    /// compares it to the local sound speed). `inv_sq` is the caller's
    /// Σ_d 1/Δ_d² over the open directions (θ included at N_θ > 1); `delta`
    /// the same Δ the rate law's matched coefficients use. (0, 0) for
    /// non-reactive cells — the member exists exactly where the rate law
    /// is live.
    pub fn front_carrier_signal(
        &self,
        w: &Prim,
        delta: f64,
        inv_sq: f64,
    ) -> Result<(f64, f64), &'static str> {
        let (rho, p, z, b) = (w[I_RHO], w[4], w[I_RC], w[I_RB]);
        let b_c = b.clamp(0.0, 1.0);
        if b_c >= 1.0 - BURN_COMPLETE || self.blend.below_unburnt_floor(w) {
            return Ok((0.0, 0.0));
        }
        let t_u = self.blend.unburnt_temperature(w)?;
        if self.ignition.non_reactive_floor(p, t_u) {
            return Ok((0.0, 0.0));
        }
        let s_l = self.ignition.s_l(p, t_u, z)?;
        let s_t = (s_l * self.wrinkling).max(0.0);
        if s_t <= 0.0 {
            return Ok((0.0, 0.0));
        }
        let rho_u = self.blend.unburnt_density(w)?;
        let (d_c, k) = front_coeffs(s_t, delta, self.theta);
        Ok((2.0 * d_c * inv_sq + C_NAGUMO_SLOPE * (rho_u / rho) * k, s_t))
    }

    /// The auto-ignition rate `R = (ρ − ρb)/τ_ign` [kg/m³/s] at a projected
    /// primitive — the class-`R` term's explicit evaluation (the SDC node-0
    /// rate `R(U⁰)` of the trapezoid quadrature). Guarded exactly as the
    /// implicit solve: zero at the burnt fixed point ([`BURN_COMPLETE`]) and
    /// below the cold non-reactive floor. `(ρ − ρb)/τ` equals `ρ(1−b)/τ` on
    /// the physical range and is the linear-in-`ρb` form the implicit node
    /// solve inverts in closed form.
    pub fn auto_rate(&self, w: &Prim) -> Result<f64, &'static str> {
        let (rho, p, z, b) = (w[I_RHO], w[4], w[I_RC], w[I_RB]);
        let b_c = b.clamp(0.0, 1.0);
        if b_c >= 1.0 - BURN_COMPLETE || self.blend.below_unburnt_floor(w) {
            return Ok(0.0);
        }
        let t_u = self.blend.unburnt_temperature(w)?;
        if self.ignition.non_reactive_floor(p, t_u) {
            return Ok(0.0);
        }
        let tau = self.ignition.tau_ign(p, t_u, z)?;
        Ok(rho * (1.0 - b_c) / tau)
    }

    /// The class-`R` **node-0 quadrature rate**: [`Self::auto_rate`] capped at
    /// the rate that exactly reaches the guarded law's parking point over the
    /// step, `(ρ(1 − BURN_COMPLETE) − ρb)/Δt`. The trapezoid correction
    /// composes `+Δt/2·R⁰` explicitly, and at extreme stiffness the raw rate
    /// `ρ/τ` is orders beyond what the guarded trajectory can realize —
    /// composing it would overshoot `ρb` past every bound the implicit solve
    /// protects. The exact trajectory from `U⁰` parks within the step, so its
    /// realizable node-0 rate IS the capped value — the applied-increments
    /// doctrine (COUP-2) extended to the quadrature's explicit operand.
    /// Inactive (bit-identical to `auto_rate`) in the mild regime.
    pub fn auto_rate_node0(&self, w: &Prim, dt: f64) -> Result<f64, &'static str> {
        let r = self.auto_rate(w)?;
        if r <= 0.0 || dt <= 0.0 {
            return Ok(r);
        }
        let rho = w[I_RHO];
        let cap = rho * (1.0 - BURN_COMPLETE);
        let rb = rho * w[I_RB].clamp(0.0, 1.0);
        Ok(r.min(((cap - rb) / dt).max(0.0)))
    }

    /// The **class-`R` implicit node solve** (COUP-3 §3.3's cell-local slot,
    /// chemical occupant, S7): advance `ρb` through the auto-ignition term
    /// over a node weight `w_new`, backward-Euler —
    /// `x = base + w_new·(ρ − x)/τ(state(x))` — by the **fixed-structure
    /// τ-refreeze scheme**: the equation is linear in `x` at frozen `τ`
    /// (closed form, unconditionally stable, no iteration to diverge), and
    /// `τ`'s weak dependence on `x` (through the blend's projected `p`) is
    /// converged by a fixed [`N_TAU_REFREEZE`] re-evaluations of `τ` at the
    /// current iterate — a pure function of the cell state, bit-deterministic.
    /// The guarded law's burnt fixed point is honored **exactly**: the source
    /// vanishes at `ρb ≥ ρ(1 − BURN_COMPLETE)`, so when the unguarded update
    /// crosses it the trajectory parks there — that cap is the exact integral
    /// of the declared discontinuous rate law, not a clamp of an overshoot.
    /// In the mild limit (`w_new ≪ τ`) the update reduces to the explicit
    /// value to `O((w_new/τ)²)` — one treatment across the whole regime, no
    /// stiffness branch (Rule 12). Returns `(ρb_new, realized rate)` where
    /// the realized rate `(ρb_new − base)/w_new` is what the ledger and the
    /// SDC quadrature carry (COUP-2: applied increments, never rate×Δt).
    ///
    /// `u` is the cell's composed conserved state EXCLUDING this term's new-
    /// node contribution; `base` is its `ρb` slot plus the quadrature's
    /// node-0/previous-sweep terms.
    pub fn implicit_auto_update(
        &self,
        u: &Cons,
        base: f64,
        w_new: f64,
    ) -> Result<(f64, f64), &'static str> {
        if !base.is_finite() {
            return Err("non-finite composed burn progress entering the class-R solve");
        }
        let rho = u[I_RHO];
        let cap = rho * (1.0 - BURN_COMPLETE);
        // A cell whose ADVECTED ρc already sits at/past the cap (the
        // EPS_BURN_BOUND-scale shock-consistency excursion of b = ρc/ρ):
        // the declared law's source is zero there, so the exact trajectory
        // leaves ρc unchanged — shaving it here would be a silent clamp of
        // advected state (S8 review). The loud positivity guard owns any
        // real excursion.
        if u[I_RB] >= cap {
            return Ok((u[I_RB], (u[I_RB] - base) / w_new));
        }
        // The node result is PROJECTED onto the guarded law's invariant set
        // [0, cap] (S7, the ◆C2 shake-out finding): [0, cap] is forward-
        // invariant for the declared law (source ≥ 0, zero at the cap), so
        // every exact trajectory from an admissible state stays inside —
        // but the trapezoid QUADRATURE base mixes the node-0 rate (the
        // pre-runaway state, small) with the previous sweep's realized rate
        // (mid-runaway, large), and during a single-step thermal runaway
        // that mismatch legitimately drives the composed base ~0.1·ρ
        // OUTSIDE the set (measured: −0.08ρ at the ◆C2 spark kernel). The
        // projection is the same exact-integral statement as the parked
        // cap — the applied (realized) increment is what the ledger and
        // the audit carry, so conservation bookkeeping is exact either
        // way; a base inside the set is untouched.
        let project = |x: f64| x.clamp(0.0, cap);
        // SYMMETRIC base projection (S8 review — the ◆C2 finding made
        // two-sided): the trapezoid quadrature base legitimately leaves the
        // law's invariant set [0, cap] on EITHER side during a single-step
        // runaway (measured −0.08ρ below on the spark kernel; above on the
        // burnt fixed point). The exact trajectory from the admissible
        // state stays inside, so the ARTIFACT — never the state — is
        // discarded, on both sides alike (the S7 code projected only the
        // high side before the solve; a negative base was integrated
        // through the BE and leaked into the accepted ρc at O(1/(1+w/τ))).
        let base_raw = base;
        let base = project(base);
        if base == cap {
            // Parked: the solve from the cap stays at the cap exactly.
            return Ok((cap, (cap - base_raw) / w_new));
        }
        let mut x = base;
        for _ in 0..N_TAU_REFREEZE {
            let mut ut = *u;
            ut[I_RB] = x;
            let w = self.blend.prim_checked(&ut)?;
            if self.blend.below_unburnt_floor(&w) {
                // Colder than any representable reactant: source zero — the
                // exact trajectory from the projected base is the projected
                // base; the realized rate carries the FULL discrepancy vs
                // the raw quadrature base (ledger identity).
                return Ok((base, (base - base_raw) / w_new));
            }
            let (p, z) = (w[4], w[I_RC]);
            let t_u = self.blend.unburnt_temperature(&w)?;
            if self.ignition.non_reactive_floor(p, t_u) {
                return Ok((base, (base - base_raw) / w_new));
            }
            let tau = self.ignition.tau_ign(p, t_u, z)?;
            // BE at frozen τ: x = (base + w·ρ/τ)/(1 + w/τ), parked at the cap.
            x = ((base + w_new * rho / tau) / (1.0 + w_new / tau)).min(cap);
            if x == cap {
                break; // parked: further τ refreshes cannot move it
            }
        }
        let x = project(x);
        Ok((x, (x - base_raw) / w_new))
    }
}

/// Fixed τ-refreeze count of the class-`R` implicit node solve
/// ([`Combustion::implicit_auto_update`]). **Contract (restated at S9,
/// SOLV-4 §3.6 v0.4.8):** the fixed count is STRUCTURE, not a convergence
/// claim — the S9 witness measured the original "weak τ-dependence"
/// premise false in the reaction-driven-compression band (within one node
/// solve the burn's constant-volume compression heating can drive τ down
/// ~×1/100; measured node lag vs the converged frozen-τ fixed point:
/// ~3.0e-1 at w/τ = 0.3, 2.2e-2 at 1, 2.1e-3 at 3). The residual per-node
/// lag is a temporal-truncation-class term the SDC sweeps' own
/// re-evaluations absorb — measured composed temporal order 1.8–2.1
/// through that same band (the dt-Richardson gate owns accuracy; the
/// envelope pin catches lag growth) — the same split as the truncated
/// gas-diffusion Picard (COUP-3 §3.1). Fixed structure, never adaptive
/// (COUP-3 §3.7).
pub const N_TAU_REFREEZE: usize = 2;

/// The **reacting measure** `R = ∫ b(1−b)·(the live burn rate) dV` (SOLV-4
/// §3.6) — the run's flame content, the single owner of the COUP-4
/// `NEVER_IGNITED`/`FLAMEOUT` thresholds. Computed from the committed grid
/// state (not the sweep cache) for the post-step halt check.
pub fn reacting_measure(
    g: &Grid,
    ids: &[FieldId; NCOMP],
    blend: &BurnBlendEos<'_>,
    ignition: &IgnitionColumns<'_>,
    wrinkling: f64,
    theta: f64,
) -> Result<f64, FlowError> {
    let spec = g.spec();
    let delta = (spec.dr * spec.dz).sqrt();
    let mut r = 0.0;
    for bi in 0..g.n_bricks() {
        let brick = g.brick(bi);
        let mask = brick.mask();
        let nt = brick.n_theta();
        let ids_local: [&[f64]; NCOMP] = std::array::from_fn(|k| brick.field(ids[k]));
        for local in 0..BRICK_CELLS {
            if mask & (1u64 << local) == 0 {
                continue;
            }
            let (i_r, i_z) = brick.global_rz(local);
            for j in 0..nt {
                let cell = brick.cell_index(j, local);
                let u: [f64; NCOMP] = std::array::from_fn(|k| ids_local[k][cell]);
                let w = blend
                    .prim_checked(&u)
                    .map_err(|what| nonphys(i_r, i_z, j, what))?;
                let (rho, p, z, b) = (w[I_RHO], w[4], w[I_RC], w[I_RB]);
                let b_c = b.clamp(0.0, 1.0);
                // The BURN_COMPLETE domain guard (as `accumulate`): a fully-burnt
                // cell's rate is 0 and its post-flame gas is not a reactant, so
                // the unburnt/ignition closures are not queried on it (its
                // c(1−c)-weighted contribution here is ≤ 1e-3 and exactly 0 now).
                if b_c >= 1.0 - BURN_COMPLETE || blend.below_unburnt_floor(&w) {
                    continue;
                }
                let t_u = blend
                    .unburnt_temperature(&w)
                    .map_err(|what| nonphys(i_r, i_z, j, what))?;
                if ignition.non_reactive_floor(p, t_u) {
                    continue; // cold-side declared no-burn floor: rate is 0 there
                }
                let s_l = ignition
                    .s_l(p, t_u, z)
                    .map_err(|what| nonphys(i_r, i_z, j, what))?;
                let tau = ignition
                    .tau_ign(p, t_u, z)
                    .map_err(|what| nonphys(i_r, i_z, j, what))?;
                let s_t = (s_l * wrinkling).max(0.0);
                let (_d_c, k) = front_coeffs(s_t, delta, theta);
                let rho_u = blend
                    .unburnt_density(&w)
                    .map_err(|what| nonphys(i_r, i_z, j, what))?;
                // The live burn rate at the cell (reaction + auto-ignition), the
                // same terms the source applies; weighted by b(1−b) (flame content)
                // and the cell volume. The propagation reaction's positive part
                // (the c > a core drives it; the metastable c < a fringe pulls
                // back) keeps R a clean flame-activity indicator.
                let rate = propagation_rate(rho_u, k, b_c).max(0.0) + rho * (1.0 - b_c) / tau;
                // Per-sector κ (S11), matching `accumulate`'s source weighting;
                // bit-identical to the plane-0 `kappa_rz` on θ-uniform worlds.
                r += b_c * (1.0 - b_c) * rate * brick.kappa_cell(j, local) * g.cell_volume(i_r, nt);
            }
        }
    }
    Ok(r)
}

/// The **consumption rate** `∫(ρ_u·K·b(1−b) + ρ(1−b)/τ_ign) dV` — the total
/// mass-burning rate [kg/s]. For a resolved premixed front it equals
/// `ρ_u·S_L·A_flame`, so `S_L = consumption_rate/(ρ_u·A)` is the closure-set
/// **consumption speed** (frame- and volume-independent — the clean
/// grid-independence measure). Computed from the committed grid state.
pub fn consumption_rate(
    g: &Grid,
    ids: &[FieldId; NCOMP],
    blend: &BurnBlendEos<'_>,
    ignition: &IgnitionColumns<'_>,
    wrinkling: f64,
    theta: f64,
) -> Result<f64, FlowError> {
    let spec = g.spec();
    let delta = (spec.dr * spec.dz).sqrt();
    let mut acc = 0.0;
    for bi in 0..g.n_bricks() {
        let brick = g.brick(bi);
        let mask = brick.mask();
        let nt = brick.n_theta();
        let ids_local: [&[f64]; NCOMP] = std::array::from_fn(|k| brick.field(ids[k]));
        for local in 0..BRICK_CELLS {
            if mask & (1u64 << local) == 0 {
                continue;
            }
            let (i_r, i_z) = brick.global_rz(local);
            for j in 0..nt {
                let cell = brick.cell_index(j, local);
                let u: [f64; NCOMP] = std::array::from_fn(|k| ids_local[k][cell]);
                let w = blend
                    .prim_checked(&u)
                    .map_err(|what| nonphys(i_r, i_z, j, what))?;
                let (rho, p, z, b) = (w[I_RHO], w[4], w[I_RC], w[I_RB]);
                let b_c = b.clamp(0.0, 1.0);
                // The BURN_COMPLETE domain guard (as `accumulate`): consumption of
                // the last ≤ 1e-3 is complete — no closure queries on post-flame gas.
                if b_c >= 1.0 - BURN_COMPLETE || blend.below_unburnt_floor(&w) {
                    continue;
                }
                let t_u = blend
                    .unburnt_temperature(&w)
                    .map_err(|w2| nonphys(i_r, i_z, j, w2))?;
                if ignition.non_reactive_floor(p, t_u) {
                    continue; // cold-side declared no-burn floor: rate is 0 there
                }
                let rho_u = blend
                    .unburnt_density(&w)
                    .map_err(|w2| nonphys(i_r, i_z, j, w2))?;
                let s_l = ignition
                    .s_l(p, t_u, z)
                    .map_err(|w2| nonphys(i_r, i_z, j, w2))?;
                let tau = ignition
                    .tau_ign(p, t_u, z)
                    .map_err(|w2| nonphys(i_r, i_z, j, w2))?;
                let s_t = (s_l * wrinkling).max(0.0);
                let (_d_c, k) = front_coeffs(s_t, delta, theta);
                let rate = propagation_rate(rho_u, k, b_c) + rho * (1.0 - b_c) / tau;
                // Per-sector κ (S11), matching `accumulate`; bit-identical to
                // plane-0 `kappa_rz` on θ-uniform worlds.
                acc += rate * brick.kappa_cell(j, local) * g.cell_volume(i_r, nt);
            }
        }
    }
    Ok(acc)
}

fn nonphys(i_r: usize, i_z: usize, i_theta: u32, what: &'static str) -> FlowError {
    FlowError::NonPhysicalState {
        i_r,
        i_z,
        i_theta,
        what,
    }
}
