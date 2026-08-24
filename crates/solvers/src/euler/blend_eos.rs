//! SOLV-4 §3.6 — the **burn-progress blended EOS occupant**. The chemical
//! regime's thermochemistry is the `b`-blend (`b` = the doc's burn-progress
//! `c`, the code's `I_RB` slot) of two OFFL-3 branches, both bound as
//! [`TableEos`] on the same `(p, h, Z)` schema: **unburnt** (the S5 frozen-
//! reactant surface) and **burnt** (the shifting-equilibrium surface).
//!
//! **The energy-conserving flamelet closure (SOLV-4 §3.6).** Both sub-states
//! share the cell pressure `p` **and** the cell's actual static specific
//! enthalpy `h = e + p/ρ` — the adiabatic constant-`p` identity
//! `h_b = h_u = h` on S5's shared CEA formation reference (burning conserves
//! enthalpy), so the two branches are interrogated at the **same** `(p, h, Z)`
//! and differ only in composition. The cell's one mean density is partitioned
//! by **mass-weighted specific volume**
//! `1/ρ = (1−b)/ρ_u(p, h, Z) + b/ρ_b(p, h + h_offset, Z)`, and the pressure is
//! the root of that identity — the [`TableEos`] equilibrium projection
//! generalized to interrogate both branches (one deterministic Illinois
//! root find). It recovers the limits exactly: `b = 1` is the pure burnt
//! projection **bit-for-bit** (shifting mode), `b = 0` the pure unburnt one.
//!
//! `T_u = T_unburnt(p, h, Z)` (the unburnt branch's temperature at the cell's
//! actual `h`) is the reaction coordinate the SOLV-4 §3.6 rate law keys on:
//! cold at an energy-conserving flame front (`h` ≈ uniform), and **hot where
//! an igniter deposit or compression raised `h`** — which is exactly why the
//! igniter fires `τ_ign` (`T_u` sees the deposited enthalpy).

use super::{
    Cons, EosLaw, FlowError, I_EI, I_EN, I_G1, I_MR, I_MT, I_MZ, I_RB, I_RC, I_RHO, Prim, TableEos,
};

/// **Pure-BURNT threshold** on the `[0,1]` progress variable `b`: at
/// `b ≥ 1 − EPS_B_PURE_BURNT` the cell is read as pure burnt (the unburnt
/// branch is skipped). **Derived from the one owner,
/// [`super::BURN_COMPLETE`] (SOLV-4 0.4.3):** the combustion source zeroes
/// its rate (and its closure queries: post-flame gas is not a reactant, the
/// §3.6 domain guard) at `b ≥ 1 − BURN_COMPLETE`, so `b` **asymptotes
/// toward that fixed point from below and never crosses it** (measured:
/// cells park ~2e-7 under it). The blend's pure-burnt region must therefore
/// *strictly contain* that pinning attractor — an equal threshold is a
/// knife edge on whose wrong side every burnt cell keeps interrogating the
/// unburnt branch at ~10⁻³ weight forever, and hot burnt gas (whose `h`
/// legitimately exceeds the unburnt surface's ceiling) refuses on a branch
/// describing 0.1 % of its mass. Factor 2 gives the attractor ~100× margin
/// over the largest per-step creep (`Δb ~ (1−b)·Δt/τ`). The crossing step
/// this skip declares is `EPS·(v_u/v_b)` ≈ **2.6e-4** relative at flame
/// states (the dropped unburnt branch is the *dense* one, `v_u ≪ v_b`) —
/// far inside the density column's own interp bound. **Not** a physics
/// threshold; the physical `b ∈ [0,1]` bound is owned by the combustion
/// source's loud positivity guard (SOLV-4 §3.6), never clamped here.
pub const EPS_B_PURE_BURNT: f64 = 2.0 * super::BURN_COMPLETE;

/// **Pure-UNBURNT threshold**: at `b ≤ EPS_B_PURE_UNBURNT` the cell is read
/// as pure unburnt (the burnt branch is skipped — the envelope-robustness
/// motive: a cold `b = 0`⁺ cell projects without a spurious burnt-branch
/// query). **Deliberately asymmetric to [`EPS_B_PURE_BURNT`] (review
/// finding, S6 close):** the reaction pins an attractor only at the BURNT
/// end — on the cold side the metastable fringe (`b < a`) decays smoothly
/// through every value to 0, parking nowhere — and the cold-side crossing
/// step scales as `EPS·(v_b/v_u)` ≈ **7.8·EPS** at flame states (here the
/// dropped burnt branch is the *light* one), so a 2e-3 threshold would step
/// the specific volume ~1.5 % — *outside* the density column's own ~0.8 %
/// bound. At 1e-9 the step is ~8e-9: a true regularization.
pub const EPS_B_PURE_UNBURNT: f64 = 1.0e-9;

/// Inward relative shrink on the closed-form h-derived pressure bracket so
/// endpoint queries cannot fall an ulp outside either surface's envelope.
const H_BRACKET_MARGIN: f64 = 1.0e-12;

/// Fixed log-scan resolution of the non-monotone fallback (mirrors
/// [`TableEos`]'s `N_P_SCAN`); sized to match the tables' own p-axis density.
const N_P_SCAN: usize = 64;

/// The SOLV-4 §3.6 blended occupant: two `(p, h, Z)` branches + the S18
/// `h_offset` knockdown (applied to the **burnt** interrogation only — the
/// unburnt reactants carry no combustion-completeness deficit).
#[derive(Debug, Clone)]
pub struct BurnBlendEos<'t> {
    unburnt: TableEos<'t>,
    burnt: TableEos<'t>,
    /// S18 η_c\* knockdown on the burnt branch (J/kg); 0 = full equilibrium.
    pub h_offset: f64,
}

impl<'t> BurnBlendEos<'t> {
    /// Bind the blend to its two already-bound branches. The unburnt branch
    /// must be the S5 frozen-reactant surface, the burnt branch the
    /// equilibrium surface — both `(p, h, Z)` (checked by `TableEos::bind`).
    pub fn new(unburnt: TableEos<'t>, burnt: TableEos<'t>) -> Self {
        Self {
            unburnt,
            burnt,
            h_offset: 0.0,
        }
    }

    /// The burnt branch (S18 knockdown carrier + shifting-mode fast paths).
    pub fn burnt(&self) -> &TableEos<'t> {
        &self.burnt
    }

    /// Conserved state from `(p, h, Z, b)` and a velocity — the forward blend
    /// (density = the mass-weighted specific volume at the common `(p, h, Z)`),
    /// for constructor initial conditions. `b = 0` is the pure unburnt state,
    /// `b = 1` the pure burnt — both at the *same* `h` (the adiabatic-flamelet
    /// identity), so a flame front is initialized by varying only `b`.
    pub fn cons_from_phzb(
        &self,
        p: f64,
        h: f64,
        z: f64,
        b: f64,
        vel: [f64; 3],
    ) -> Result<Cons, &'static str> {
        let b_c = b.clamp(0.0, 1.0);
        let rho = 1.0 / self.blend_inv_rho(p, h, z, b_c)?;
        let e = h - p / rho;
        let ke = 0.5 * (vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]);
        Ok([
            rho,
            rho * vel[0],
            rho * vel[1],
            rho * vel[2],
            rho * (e + ke),
            rho * z,
            rho * b,
        ])
    }

    /// `1/ρ` of the blend at a trial `(p, h, Z)` and burn fraction `b`: the
    /// mass-weighted specific volume, skipping a branch of ~zero weight (so
    /// the pure limits never query the other surface's envelope).
    fn blend_inv_rho(&self, p: f64, h: f64, z: f64, b: f64) -> Result<f64, &'static str> {
        let mut v = 0.0;
        if b < 1.0 - EPS_B_PURE_BURNT {
            let ru = self
                .unburnt
                .rho_at(p, h, z)
                .map_err(|_| "unburnt branch query failed inside the blend projection bracket")?;
            v += (1.0 - b) / ru;
        }
        if b > EPS_B_PURE_UNBURNT {
            let rb = self
                .burnt
                .rho_at(p, h + self.h_offset, z)
                .map_err(|_| "burnt branch query failed inside the blend projection bracket")?;
            v += b / rb;
        }
        Ok(v)
    }

    /// Project the pressure: root of `blend_inv_rho(p) − 1/ρ = 0` over the
    /// admissible bracket (the intersection of the *required* branches'
    /// envelopes). Deterministic Illinois; a mild non-monotone corner falls
    /// back to a fixed log-scan for a sign change, and a genuinely off-surface
    /// state refuses (META-1 P6 — the near-vacuum tangency-acceptance the pure
    /// `TableEos` carries is an RL10-plume feature, deferred to S7).
    fn project_pressure(&self, rho: f64, e: f64, z: f64, b: f64) -> Result<f64, &'static str> {
        let inv = 1.0 / rho;
        let need_u = b < 1.0 - EPS_B_PURE_BURNT;
        let need_b = b > EPS_B_PURE_UNBURNT;
        let [(pu_lo, pu_hi), (hu_lo, hu_hi), _] = self.unburnt.envelopes();
        let [(pb_lo, pb_hi), (hb_lo, hb_hi), _] = self.burnt.envelopes();

        // h-window (of h = e + p/ρ) and p-window, intersected over the
        // required branches; the burnt window is shifted by −h_offset because
        // that branch is queried at h + h_offset.
        let (mut h_lo, mut h_hi) = (f64::NEG_INFINITY, f64::INFINITY);
        let (mut p_lo_env, mut p_hi_env) = (0.0_f64, f64::INFINITY);
        if need_u {
            h_lo = h_lo.max(hu_lo);
            h_hi = h_hi.min(hu_hi);
            p_lo_env = p_lo_env.max(pu_lo);
            p_hi_env = p_hi_env.min(pu_hi);
        }
        if need_b {
            h_lo = h_lo.max(hb_lo - self.h_offset);
            h_hi = h_hi.min(hb_hi - self.h_offset);
            p_lo_env = p_lo_env.max(pb_lo);
            p_hi_env = p_hi_env.min(pb_hi);
        }
        let p_from_h_lo = rho * (h_lo - e);
        let p_from_h_hi = rho * (h_hi - e);
        let mut lo = p_lo_env
            .max(p_from_h_lo * (1.0 + H_BRACKET_MARGIN))
            .max(0.0);
        let hi = p_hi_env.min(p_from_h_hi * (1.0 - H_BRACKET_MARGIN));
        if !(lo.is_finite() && hi.is_finite()) || lo >= hi || hi <= 0.0 {
            return Err("no admissible pressure bracket: blended state outside the table envelope");
        }
        lo = lo.max(0.0);

        let g = |p: f64| -> Result<f64, &'static str> {
            Ok(self.blend_inv_rho(p, e + p * inv, z, b)? - inv)
        };
        let ga0 = g(lo)?;
        let gb0 = g(hi)?;
        if ga0 == 0.0 {
            return Ok(lo);
        }
        if gb0 == 0.0 {
            return Ok(hi);
        }
        if ga0 * gb0 < 0.0 {
            return TableEos::illinois_root(lo, hi, ga0, gb0, &g);
        }
        // Mild non-monotone corner: fixed log-scan for the first sign change.
        let ratio = hi / lo;
        let mut prev_p = lo;
        let mut prev_g = ga0;
        for k in 1..=N_P_SCAN {
            let pk = lo * ratio.powf(k as f64 / N_P_SCAN as f64);
            let gk = g(pk)?;
            if gk == 0.0 {
                return Ok(pk);
            }
            if prev_g * gk < 0.0 {
                return TableEos::illinois_root(prev_p, pk, prev_g, gk, &g);
            }
            prev_p = pk;
            prev_g = gk;
        }
        Err(
            "blended state off both surfaces over the admissible bracket — \
             no sign change (S7 adds the near-vacuum tangency acceptance)",
        )
    }

    /// The unburnt-branch temperature `T_u(p, h, Z)` at a projected primitive
    /// — the SOLV-4 §3.6 rate-law coordinate (`S_L`, `τ_ign`). `h` is the
    /// cell's actual static enthalpy `e + p/ρ`, so `T_u` rises with an igniter
    /// deposit.
    pub fn unburnt_temperature(&self, w: &Prim) -> Result<f64, &'static str> {
        let (rho, p, z) = (w[I_RHO], w[4], w[I_RC]);
        let h = w[I_EI] + p / rho;
        self.unburnt
            .temp_at(p, h, z)
            .map_err(|_| "unburnt temperature query failed (T_u off the ignition-surface envelope)")
    }

    /// The unburnt-branch density `ρ_u(p, h, Z)` — the `ρ_u` factor of the
    /// propagation source `ρ_u·S_T·|∇c|` (SOLV-4.4).
    pub fn unburnt_density(&self, w: &Prim) -> Result<f64, &'static str> {
        let (rho, p, z) = (w[I_RHO], w[4], w[I_RC]);
        let h = w[I_EI] + p / rho;
        self.unburnt
            .rho_at(p, h, z)
            .map_err(|_| "unburnt density query failed inside the blend")
    }

    fn prim_checked_impl(&self, u: &Cons) -> Result<Prim, &'static str> {
        let rho = u[I_RHO];
        if !rho.is_finite() || rho <= 0.0 {
            return Err("non-positive or non-finite density");
        }
        let inv = 1.0 / rho;
        let (ur, ut, uz) = (u[I_MR] * inv, u[I_MT] * inv, u[I_MZ] * inv);
        let e = u[I_EN] * inv - 0.5 * (ur * ur + ut * ut + uz * uz);
        if !e.is_finite() {
            return Err("non-finite internal energy");
        }
        let z = u[I_RC] * inv;
        if !z.is_finite() {
            return Err("non-finite elemental mixture fraction");
        }
        let b_raw = u[I_RB] * inv;
        if !b_raw.is_finite() {
            return Err("non-finite burn progress");
        }
        // The blend math needs b ∈ [0,1]; a tiny advective over/undershoot is
        // a numerical fact about a [0,1] fraction, so it is bounded HERE for
        // the projection only — the STORED b (and the audited ρb) stay the
        // true advected value, and the physical [0,1] violation guard belongs
        // to the combustion source (SOLV-4 §3.6), not this projection.
        let b = b_raw.clamp(0.0, 1.0);
        let p = self.project_pressure(rho, e, z, b)?;
        let h = e + p * inv;
        // Mass-weighted blended sound speed a² = (1−b)a_u² + b·a_b² (declared
        // model-form; skips a ~zero-weight branch, matching the projection).
        let mut a2 = 0.0;
        if b < 1.0 - EPS_B_PURE_BURNT {
            let au = self
                .unburnt
                .sound_at(p, h, z)
                .map_err(|_| "unburnt sound-speed query failed at the projected state")?;
            a2 += (1.0 - b) * au * au;
        }
        if b > EPS_B_PURE_UNBURNT {
            let ab = self
                .burnt
                .sound_at(p, h + self.h_offset, z)
                .map_err(|_| "burnt sound-speed query failed at the projected state")?;
            a2 += b * ab * ab;
        }
        let g1 = rho * a2 / p;
        Ok([rho, ur, ut, uz, p, z, b_raw, e, g1])
    }
}

impl EosLaw for BurnBlendEos<'_> {
    fn prim_checked(&self, u: &Cons) -> Result<Prim, &'static str> {
        self.prim_checked_impl(u)
    }

    fn prim_to_cons(&self, w: &Prim) -> Cons {
        let rho = w[I_RHO];
        let ke = 0.5 * (w[1] * w[1] + w[2] * w[2] + w[3] * w[3]);
        [
            rho,
            rho * w[1],
            rho * w[2],
            rho * w[3],
            rho * (w[I_EI] + ke),
            rho * w[I_RC],
            rho * w[I_RB],
        ]
    }

    fn total_energy(&self, w: &Prim) -> f64 {
        let rho = w[I_RHO];
        rho * (w[I_EI] + 0.5 * (w[1] * w[1] + w[2] * w[2] + w[3] * w[3]))
    }

    fn sound_speed_w(&self, w: &Prim) -> f64 {
        (w[I_G1] * w[4] / w[I_RHO]).sqrt()
    }

    fn roe_sound_speed(
        &self,
        wl: &Prim,
        wr: &Prim,
        _h_roe: f64,
        _q2_roe: f64,
        sql: f64,
        sqr: f64,
        inv: f64,
    ) -> f64 {
        // The datum-free general-EOS Roe sound speed (as `TableEos`): average
        // the squared sound speeds `c² = Γ₁·p/ρ` from each side's own
        // projected (blended) state — always positive, exact for identical
        // states, and formation-reference-safe.
        let cl2 = wl[I_G1] * wl[4] / wl[I_RHO];
        let cr2 = wr[I_G1] * wr[4] / wr[I_RHO];
        ((sql * cl2 + sqr * cr2) * inv).sqrt()
    }

    fn stagnation_ghost(
        &self,
        _p0: f64,
        _rho0: f64,
        _c_frac: f64,
        _u_n: f64,
        _normal: usize,
    ) -> Result<Prim, FlowError> {
        // No closed-form isentrope on tabulated surfaces; the COUP-7 injector
        // object owns inflow. Refuse, never approximate (META-1 P6).
        Err(FlowError::BcUnsupportedByEos {
            bc: "StagnationInflow",
        })
    }

    /// Injector inflow: the injected propellant is **unburnt** (`b = 0`), so
    /// the face state is the unburnt branch's mass-flow ghost (which already
    /// sets the burn slot to 0). The reactants burn in the chamber.
    fn mass_flow_inflow_ghost(
        &self,
        mdot_per_area: f64,
        h_total: f64,
        c_frac: f64,
        p_int: f64,
        sign: f64,
        normal: usize,
    ) -> Result<Prim, FlowError> {
        self.unburnt
            .mass_flow_inflow_ghost(mdot_per_area, h_total, c_frac, p_int, sign, normal)
    }
}
