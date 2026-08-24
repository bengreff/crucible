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
//! **Integration scope (S6, mini-sim tier):** all three terms are advanced
//! **explicitly** (the matched `D_c` is non-stiff since `S_T ≪ a`, deflagration
//! `τ_ign` is long) behind a **loud positivity guard** — `b` leaving `[0,1]`
//! beyond [`EPS_BURN_BOUND`] is a halt (under-resolution: refine Δt or gentle
//! the igniter ramp), never a clamp (META-1 P6). The stiff class-`R` implicit
//! auto-ignition (RL10 knock-class end-gas) is a declared S7 hardening. The
//! serial cell loop is bit-deterministic at any thread count; N_θ = 1 only
//! (the 3-D re-key rides plan S8, like `gas_diffusion`).

use super::{BurnBlendEos, EosLaw, FlowError, I_RB, I_RC, I_RHO, NCOMP, Prim, Scratch};
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
/// *strictly contain* the pinning attractor, since `c` asymptotes to just
/// *under* `1 − BURN_COMPLETE` and never crosses; an equal threshold is a
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
/// ignition and halts (never a silent clamp). Sized well above passive-scalar
/// round-off and well below any physical excursion.
pub const EPS_BURN_BOUND: f64 = 1.0e-6;

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
            flame_speed,
            ignition_delay,
        })
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
        let nt = g.brick(0).n_theta();
        if nt != 1 {
            return Err(FlowError::NonPhysicalState {
                i_r: usize::MAX,
                i_z: usize::MAX,
                i_theta: 0,
                what: "combustion at N_θ > 1 arrives with the 3-D wave (plan S8)",
            });
        }
        let spec = g.spec();
        let (dr, dz) = (spec.dr, spec.dz);
        let delta = (dr * dz).sqrt(); // local cell length Δ (isotropic on a square grid)
        let (rho_id, rhob_id) = (ids[I_RHO], ids[I_RB]);

        // b at a global cell (ρb/ρ); returns None outside the domain or on a
        // non-gas cell (a wall face ⇒ zero-flux for b, SOLV-4 §3.6).
        let b_at = |i_r: usize, i_z: usize| -> Option<f64> {
            if i_r >= spec.n_r
                || i_z >= spec.n_z
                || !g.is_active(i_r, i_z)
                || g.kappa(i_r, i_z) <= 0.0
            {
                return None;
            }
            let bi = g.brick_index_by_coords((i_r / BRICK) as u32, (i_z / BRICK) as u32)?;
            let b = g.brick(bi);
            let cell = b.cell_index(0, (i_r % BRICK) * BRICK + i_z % BRICK);
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
                let w: &Prim = &prim[bi][local];
                let (rho, p) = (w[I_RHO], w[4]);
                let (z, b) = (w[I_RC], w[I_RB]);
                // Loud positivity guard (META-1 P6) — never a clamp.
                if !(-EPS_BURN_BOUND..=1.0 + EPS_BURN_BOUND).contains(&b) {
                    return Err(FlowError::NonPhysicalState {
                        i_r,
                        i_z,
                        i_theta: 0,
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
                let (react, auto, rho_dc_here) = if b_c >= 1.0 - BURN_COMPLETE {
                    (0.0, 0.0, 0.0)
                } else {
                    let t_u = self
                        .blend
                        .unburnt_temperature(w)
                        .map_err(|what| nonphys(i_r, i_z, what))?;
                    let rho_u = self
                        .blend
                        .unburnt_density(w)
                        .map_err(|what| nonphys(i_r, i_z, what))?;
                    let s_l = self
                        .ignition
                        .s_l(p, t_u, z)
                        .map_err(|what| nonphys(i_r, i_z, what))?;
                    let tau = self
                        .ignition
                        .tau_ign(p, t_u, z)
                        .map_err(|what| nonphys(i_r, i_z, what))?;
                    let s_t = (s_l * self.wrinkling).max(0.0);
                    let (d_c, k) = front_coeffs(s_t, delta, self.theta);
                    (
                        propagation_rate(rho_u, k, b_c),
                        rho * (1.0 - b_c) / tau,
                        rho * d_c,
                    )
                };

                // Front-thickening diffusion ∇·(ρ D_c ∇b) — symmetric two-point
                // face fluxes so interior faces telescope (conservative). Each
                // face's ρD_c is the two-cell arithmetic mean; a wall/edge face
                // is zero-flux (skipped).
                let vol = g.cell_volume(i_r, nt);
                let kv = brick.kappa_rz(local) * vol;
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
                    let ap = g.aperture(i_r, i_z, dir);
                    if ap <= 0.0 {
                        continue;
                    }
                    let Some((nr, nz)) = nbr else { continue };
                    let Some(b_nbr) = b_at(nr, nz) else { continue };
                    // Neighbour ρD_c (its own S_T); the face value is the mean.
                    let rho_dc_nbr = match self.face_rho_dc(g, prim, nr, nz, z, delta) {
                        Ok(v) => v,
                        Err(what) => return Err(nonphys(nr, nz, what)),
                    };
                    let rho_dc_face = 0.5 * (rho_dc_here + rho_dc_nbr);
                    diff += ap * area * rho_dc_face * (b_nbr - b) / d;
                }
                let diff_div = if kv > 0.0 { diff / kv } else { 0.0 };

                let src = react + auto + diff_div;
                rate[bi][local][I_RB] += src;
                // Reaction + diffusion both go to the source ledger; the
                // diffusion telescopes to the domain boundary (zero-flux ⇒ 0)
                // in the global sum, so the audit balances against the stored
                // Δ(ρb) exactly (COUP-2 §3.1.1).
                ledger.src_net[I_RB] += kv * src;
                ledger.src_abs[I_RB] += (kv * src).abs();
            }
        }
        Ok(())
    }

    /// A neighbour cell's `ρ·D_c` for the face-mean diffusion coefficient: its
    /// own `ρ`, `S_T` (from its `T_u`) and the matched `D_c`. Reads the
    /// neighbour's primitive from the cache (same sweep state).
    fn face_rho_dc(
        &self,
        g: &Grid,
        prim: &[Vec<Prim>],
        i_r: usize,
        i_z: usize,
        _z_self: f64,
        delta: f64,
    ) -> Result<f64, &'static str> {
        let bi = g
            .brick_index_by_coords((i_r / BRICK) as u32, (i_z / BRICK) as u32)
            .ok_or("neighbour cell not in an allocated brick")?;
        let local = (i_r % BRICK) * BRICK + i_z % BRICK;
        let w = &prim[bi][local];
        let (rho, p, z, b) = (w[I_RHO], w[4], w[I_RC], w[I_RB]);
        if b.clamp(0.0, 1.0) >= 1.0 - BURN_COMPLETE {
            return Ok(0.0); // burnt neighbour: no front-carrier (see accumulate)
        }
        let t_u = self.blend.unburnt_temperature(w)?;
        let s_l = self.ignition.s_l(p, t_u, z)?;
        let s_t = (s_l * self.wrinkling).max(0.0);
        let (d_c, _k) = front_coeffs(s_t, delta, self.theta);
        Ok(rho * d_c)
    }
}

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
    let nt = g.brick(0).n_theta();
    let spec = g.spec();
    let delta = (spec.dr * spec.dz).sqrt();
    let mut r = 0.0;
    for bi in 0..g.n_bricks() {
        let brick = g.brick(bi);
        let mask = brick.mask();
        let ids_local: [&[f64]; NCOMP] = std::array::from_fn(|k| brick.field(ids[k]));
        for local in 0..BRICK_CELLS {
            if mask & (1u64 << local) == 0 {
                continue;
            }
            let (i_r, i_z) = brick.global_rz(local);
            let cell = brick.cell_index(0, local);
            let u: [f64; NCOMP] = std::array::from_fn(|k| ids_local[k][cell]);
            let w = blend
                .prim_checked(&u)
                .map_err(|what| nonphys(i_r, i_z, what))?;
            let (rho, p, z, b) = (w[I_RHO], w[4], w[I_RC], w[I_RB]);
            let b_c = b.clamp(0.0, 1.0);
            // The BURN_COMPLETE domain guard (as `accumulate`): a fully-burnt
            // cell's rate is 0 and its post-flame gas is not a reactant, so
            // the unburnt/ignition closures are not queried on it (its
            // c(1−c)-weighted contribution here is ≤ 1e-3 and exactly 0 now).
            if b_c >= 1.0 - BURN_COMPLETE {
                continue;
            }
            let t_u = blend
                .unburnt_temperature(&w)
                .map_err(|what| nonphys(i_r, i_z, what))?;
            let s_l = ignition
                .s_l(p, t_u, z)
                .map_err(|what| nonphys(i_r, i_z, what))?;
            let tau = ignition
                .tau_ign(p, t_u, z)
                .map_err(|what| nonphys(i_r, i_z, what))?;
            let s_t = (s_l * wrinkling).max(0.0);
            let (_d_c, k) = front_coeffs(s_t, delta, theta);
            let rho_u = blend
                .unburnt_density(&w)
                .map_err(|what| nonphys(i_r, i_z, what))?;
            // The live burn rate at the cell (reaction + auto-ignition), the
            // same terms the source applies; weighted by b(1−b) (flame content)
            // and the cell volume. The propagation reaction's positive part
            // (the c > a core drives it; the metastable c < a fringe pulls
            // back) keeps R a clean flame-activity indicator.
            let rate = propagation_rate(rho_u, k, b_c).max(0.0) + rho * (1.0 - b_c) / tau;
            r += b_c * (1.0 - b_c) * rate * brick.kappa_rz(local) * g.cell_volume(i_r, nt);
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
    let nt = g.brick(0).n_theta();
    let spec = g.spec();
    let delta = (spec.dr * spec.dz).sqrt();
    let mut acc = 0.0;
    for bi in 0..g.n_bricks() {
        let brick = g.brick(bi);
        let mask = brick.mask();
        let ids_local: [&[f64]; NCOMP] = std::array::from_fn(|k| brick.field(ids[k]));
        for local in 0..BRICK_CELLS {
            if mask & (1u64 << local) == 0 {
                continue;
            }
            let (i_r, i_z) = brick.global_rz(local);
            let cell = brick.cell_index(0, local);
            let u: [f64; NCOMP] = std::array::from_fn(|k| ids_local[k][cell]);
            let w = blend
                .prim_checked(&u)
                .map_err(|what| nonphys(i_r, i_z, what))?;
            let (rho, p, z, b) = (w[I_RHO], w[4], w[I_RC], w[I_RB]);
            let b_c = b.clamp(0.0, 1.0);
            // The BURN_COMPLETE domain guard (as `accumulate`): consumption of
            // the last ≤ 1e-3 is complete — no closure queries on post-flame gas.
            if b_c >= 1.0 - BURN_COMPLETE {
                continue;
            }
            let t_u = blend
                .unburnt_temperature(&w)
                .map_err(|w2| nonphys(i_r, i_z, w2))?;
            let rho_u = blend
                .unburnt_density(&w)
                .map_err(|w2| nonphys(i_r, i_z, w2))?;
            let s_l = ignition
                .s_l(p, t_u, z)
                .map_err(|w2| nonphys(i_r, i_z, w2))?;
            let tau = ignition
                .tau_ign(p, t_u, z)
                .map_err(|w2| nonphys(i_r, i_z, w2))?;
            let s_t = (s_l * wrinkling).max(0.0);
            let (_d_c, k) = front_coeffs(s_t, delta, theta);
            let rate = propagation_rate(rho_u, k, b_c) + rho * (1.0 - b_c) / tau;
            acc += rate * brick.kappa_rz(local) * g.cell_volume(i_r, nt);
        }
    }
    Ok(acc)
}

fn nonphys(i_r: usize, i_z: usize, what: &'static str) -> FlowError {
    FlowError::NonPhysicalState {
        i_r,
        i_z,
        i_theta: 0,
        what,
    }
}
