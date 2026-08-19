//! SOLV-1 §3.4 — **shifting-equilibrium mode**: the [`EosLaw`] occupant
//! backed by the OFFL-3 `(p, h, Z)` equilibrium surface. The composition
//! slot of `U` carries the **elemental mixture fraction Z** (advected,
//! source-free — element conservation exact by construction); every
//! conserved→primitive conversion performs the **per-cell equilibrium
//! projection**: solve for the pressure `p` at which the tabulated
//! equilibrium density at `(p, h = e + p/ρ, Z)` matches the cell's conserved
//! density. Species, temperature, and heat release are whatever the surface
//! says at that state — combustion happens *in the EOS*, never as a separate
//! `if(burning)` branch (Rule 12).
//!
//! **The p-iteration** (SOLV-1 §3.4): `h = e + p/ρ` is linear in `p` at
//! fixed conserved state, so the admissible pressure bracket (the p-range
//! keeping both `p` and `h` inside the declared table envelope) is closed
//! form; the root of `ρ_tab(p, h(p), Z) − ρ` is then located by a
//! **deterministic Illinois regula-falsi** with a fixed relative tolerance
//! [`EPS_P_PROJECTION`] and fixed max iteration count [`N_P_ITER_MAX`] —
//! the META-1 §2.2 sanctioned convergence form (absolute + deterministic;
//! never wall-clock). Failure to bracket or converge is a halt diagnosis,
//! never a clamp: a state off the tabulated surface refuses (FND-5 §3.5).
//!
//! **η_c\* knockdown hook (S18, COUP-7 §3.2.1):** [`TableEos::h_offset`] is
//! the source-term-level combustion-completeness knockdown — an enthalpy
//! offset added to the `h`-coordinate of every equilibrium interrogation
//! (δh < 0 ⇒ the surface is read as if that much reaction enthalpy had not
//! been released, so T, a, Γ₁ and the emergent p_c all carry the deficit
//! self-consistently). The conserved energy is untouched (the deficit is
//! unreleased chemical energy riding in `U`); output-side multiplication of
//! c\*/Isp is forbidden (SOLV-1 §3.4). The COUP-7 injector prior tier sets
//! δh (calibrated so delivered c\* = η_c\*·c\*_ideal at the anchor state);
//! 0.0 = full equilibrium.

use super::{Cons, EosLaw, FlowError, I_EI, I_EN, I_G1, I_MR, I_MT, I_MZ, I_RC, I_RHO, Prim};
use crucible_tables::{BoundColumn, Table, TableError};

/// Fixed maximum iteration count of the equilibrium pressure projection.
/// Illinois regula-falsi typically converges in < 15; plain bisection would
/// need ~38 to reach [`EPS_P_PROJECTION`] over the full table p-envelope, so
/// 48 is headroom, not a tunable.
pub const N_P_ITER_MAX: usize = 48;

/// Relative bracket-width tolerance of the pressure projection — orders
/// below the table's own interpolation-error bound, so the projection never
/// contributes to the physics error budget.
pub const EPS_P_PROJECTION: f64 = 1e-11;

/// Inward relative shrink applied to the closed-form h-derived pressure
/// bracket so endpoint queries cannot fall an ulp outside the envelope.
const H_BRACKET_MARGIN: f64 = 1e-12;

/// Fixed sign-change scan resolution of the projection's slow path (the
/// near-vacuum non-monotone corner; see `project_pressure`).
pub const N_P_SCAN: usize = 16;

/// Slow-path acceptance: when the constraint line `h = e + p/ρ` does not
/// cross the *interpolated* ρ-contour, the closest on-line state is
/// accepted iff its density mismatch is within the **table's own declared
/// interpolation-error bound** (the producer-measured `interp_error_bound`
/// of the density column, cached at bind). Rationale (SOLV-1 §3.4): the
/// advected (ρ, e, Z) is off-surface by construction — the per-step
/// equilibrium projection re-equilibrates it, and a crossing miss within
/// the surface's own declared error IS the surface, within its band. A
/// mismatch beyond the declared bound is a genuinely off-surface state and
/// refuses (META-1 P6) — the bound is data, never a tunable.
const _EPS_PROJ_ACCEPT_DOC: () = ();

/// Fixed iteration count of the mass-flow-inflow face solve (module-header
/// determinism form): `ρ = ρ_tab(p, h_total − ½(j/ρ)², Z)` converges
/// geometrically at injector-plane Mach numbers (~0.02 for the RL10 — the
/// KE correction is ~10⁻⁴ of h_total); 8 is far past machine convergence.
pub const N_INFLOW_ITER: usize = 8;

/// The shifting-equilibrium EOS occupant: bound columns of one OFFL-3
/// `(p, h, Z)` equilibrium surface (columns resolved once — `BoundColumn`,
/// FND-5 §3.3 — for the ~10⁸ projections of an anchor run).
#[derive(Debug, Clone)]
pub struct TableEos<'t> {
    rho: BoundColumn<'t>,
    sound: BoundColumn<'t>,
    temperature: BoundColumn<'t>,
    /// The density column's producer-measured interpolation-error bound
    /// [kg/m³] — the slow-path acceptance (see module consts).
    rho_err_bound: f64,
    /// Declared envelopes of the (p, h, Z) axes, cached at bind.
    p_env: (f64, f64),
    h_env: (f64, f64),
    z_env: (f64, f64),
    /// S18 knockdown: enthalpy offset (J/kg) added to the h-coordinate of
    /// every equilibrium interrogation. See the module header.
    pub h_offset: f64,
}

impl<'t> TableEos<'t> {
    /// Bind the occupant to a loaded equilibrium surface. Refuses (never
    /// guesses) if the table is not the `(p, h, Z)` schema or a column is
    /// missing/relabeled (META-2 §4 ★ units gate runs per column).
    pub fn bind(table: &'t Table) -> Result<Self, String> {
        let rho = table
            .bind("density", "kg/m^3")
            .map_err(|e| format!("equilibrium surface: {e}"))?;
        let sound = table
            .bind("sound_speed", "m/s")
            .map_err(|e| format!("equilibrium surface: {e}"))?;
        let temperature = table
            .bind("temperature", "K")
            .map_err(|e| format!("equilibrium surface: {e}"))?;
        let expected = ["p", "h", "Z"];
        for (i, want) in expected.iter().enumerate() {
            let got = rho.axis_name(i);
            if got != *want {
                return Err(format!(
                    "equilibrium surface axis {i} is {got:?}, expected {want:?} — \
                     not a (p, h, Z) local-state surface (OFFL-3 §3.3); refusing the bind"
                ));
            }
        }
        Ok(Self {
            p_env: rho.axis_envelope(0),
            h_env: rho.axis_envelope(1),
            z_env: rho.axis_envelope(2),
            rho_err_bound: rho.interp_error_bound(),
            rho,
            sound,
            temperature,
            h_offset: 0.0,
        })
    }

    /// Warm-started projection: try a tight bracket around a hint pressure
    /// (the cell's previous-stage projection — the state moves a CFL-limited
    /// fraction per stage, so the root almost always sits within a few
    /// percent). Falls back to the full-envelope [`Self::project_pressure`]
    /// when the tight bracket does not straddle the root. A pure
    /// data-dependent path (no schedule dependence): bit-reproducible at
    /// any thread count. Measured session 12: the full-bracket Illinois
    /// dominated the entire march (~70% of wall clock in the surface
    /// interpolation it drives).
    fn project_pressure_hinted(
        &self,
        rho: f64,
        e_q: f64,
        z: f64,
        p_hint: f64,
    ) -> Result<f64, &'static str> {
        const HINT_SPREAD: f64 = 1.05;
        let inv = 1.0 / rho;
        // The same closed-form admissible bounds as the cold path.
        let p_from_h_lo = rho * (self.h_env.0 - e_q);
        let p_from_h_hi = rho * (self.h_env.1 - e_q);
        let lo_adm = self
            .p_env
            .0
            .max(p_from_h_lo * (1.0 + H_BRACKET_MARGIN))
            .max(0.0);
        let hi_adm = self.p_env.1.min(p_from_h_hi * (1.0 - H_BRACKET_MARGIN));
        let a = (p_hint / HINT_SPREAD).max(lo_adm);
        let b = (p_hint * HINT_SPREAD).min(hi_adm);
        if !(a.is_finite() && b.is_finite()) || a >= b {
            return self.project_pressure(rho, e_q, z);
        }
        let g = |p: f64| -> Result<f64, &'static str> {
            self.rho
                .interpolate(&[p, e_q + p * inv, z])
                .map(|r| r - rho)
                .map_err(|_| "equilibrium surface query failed inside the projection bracket")
        };
        let ga = g(a)?;
        if ga == 0.0 {
            return Ok(a);
        }
        let gb = g(b)?;
        if gb == 0.0 {
            return Ok(b);
        }
        if ga * gb < 0.0 {
            return Self::illinois_root(a, b, ga, gb, &g);
        }
        self.project_pressure(rho, e_q, z)
    }

    /// The deterministic Illinois regula-falsi over a sign-changing bracket
    /// (shared by the cold and warm-started projections — one root finder).
    fn illinois_root(
        mut a: f64,
        mut b: f64,
        mut ga: f64,
        mut gb: f64,
        g: &impl Fn(f64) -> Result<f64, &'static str>,
    ) -> Result<f64, &'static str> {
        for _ in 0..N_P_ITER_MAX {
            if (b - a).abs() <= EPS_P_PROJECTION * a.abs().max(b.abs()) {
                return Ok(0.5 * (a + b));
            }
            // Regula-falsi iterate, bisection fallback if it degenerates.
            let denom = gb - ga;
            let mut p = if denom != 0.0 {
                b - gb * (b - a) / denom
            } else {
                0.5 * (a + b)
            };
            let (loe, hie) = (a.min(b), a.max(b));
            if !(p > loe && p < hie) {
                p = 0.5 * (a + b);
            }
            let gp = g(p)?;
            if gp == 0.0 {
                return Ok(p);
            }
            if gp * gb < 0.0 {
                a = b;
                ga = gb;
            } else {
                // Illinois: halve the stale endpoint's residual so the
                // bracket cannot stagnate on one side.
                ga *= 0.5;
            }
            b = p;
            gb = gp;
        }
        if (b - a).abs() <= EPS_P_PROJECTION * a.abs().max(b.abs()) * 10.0 {
            return Ok(0.5 * (a + b));
        }
        Err("equilibrium pressure projection did not converge in the fixed iteration budget")
    }

    /// The equilibrium pressure projection: root of
    /// `ρ_tab(p, h(p), Z) − ρ` over the admissible bracket, where
    /// `h(p) = e_q + p/ρ`, `e_q = e + h_offset`. Deterministic Illinois
    /// regula-falsi (module header). Returns the located pressure.
    fn project_pressure(&self, rho: f64, e_q: f64, z: f64) -> Result<f64, &'static str> {
        let inv = 1.0 / rho;
        // Closed-form admissible bracket: h(p) = e_q + p/ρ must lie in the
        // h envelope AND p in the p envelope; shrink the h-derived bounds
        // inward so endpoint rounding cannot leave the envelope.
        let p_from_h_lo = rho * (self.h_env.0 - e_q);
        let p_from_h_hi = rho * (self.h_env.1 - e_q);
        let mut lo = self.p_env.0.max(p_from_h_lo * (1.0 + H_BRACKET_MARGIN));
        let hi = self.p_env.1.min(p_from_h_hi * (1.0 - H_BRACKET_MARGIN));
        if !(lo.is_finite() && hi.is_finite()) || lo >= hi || hi <= 0.0 {
            return Err("no admissible pressure bracket: state outside the table envelope");
        }
        lo = lo.max(0.0);

        let g = |p: f64| -> Result<f64, &'static str> {
            self.rho
                .interpolate(&[p, e_q + p * inv, z])
                .map(|r| r - rho)
                .map_err(|_| "equilibrium surface query failed inside the projection bracket")
        };

        let ga0 = g(lo)?;
        let gb0 = g(hi)?;
        if ga0 == 0.0 {
            return Ok(lo);
        }
        if gb0 == 0.0 {
            return Ok(hi);
        }
        let (a, b, ga, gb);
        if ga0 * gb0 < 0.0 {
            // Fast path (the dense interior of the envelope): the endpoints
            // bracket — ∂ρ/∂p dominates and g is effectively monotone.
            (a, b, ga, gb) = (lo, hi, ga0, gb0);
        } else {
            // Near-vacuum corner: the constraint line h = e + p/ρ runs
            // almost parallel to the ρ-contour of the surface (the p/ρ term
            // dominates h), so g is non-monotone and the crossing is
            // shallow. Fixed log-spaced scan for a sign change; failing
            // that, a fixed-count golden-section on |g| accepts a tangency
            // root to EPS_PROJ_ACCEPT (all deterministic, fixed order).
            let mut prev_p = lo;
            let mut prev_g = ga0;
            let mut found: Option<(f64, f64, f64, f64)> = None;
            let mut best = (prev_p, prev_g.abs());
            let ratio = hi / lo;
            for k in 1..=N_P_SCAN {
                let pk = lo * ratio.powf(k as f64 / N_P_SCAN as f64);
                let gk = g(pk)?;
                if gk == 0.0 {
                    return Ok(pk);
                }
                if gk.abs() < best.1 {
                    best = (pk, gk.abs());
                }
                if prev_g * gk < 0.0 && found.is_none() {
                    found = Some((prev_p, pk, prev_g, gk));
                }
                prev_p = pk;
                prev_g = gk;
            }
            match found {
                Some((pa, pb, gaa, gbb)) => (a, b, ga, gb) = (pa, pb, gaa, gbb),
                None => {
                    // Tangency: golden-section minimize |g| around the best
                    // sample, then accept iff the surface is met to
                    // EPS_PROJ_ACCEPT relative in ρ.
                    let phi = 0.618_033_988_749_894_9_f64;
                    let (mut x0, mut x3) = (
                        (best.0 / ratio.powf(1.0 / N_P_SCAN as f64)).max(lo),
                        (best.0 * ratio.powf(1.0 / N_P_SCAN as f64)).min(hi),
                    );
                    let mut x1 = x3 - phi * (x3 - x0);
                    let mut x2 = x0 + phi * (x3 - x0);
                    let mut f1 = g(x1)?.abs();
                    let mut f2 = g(x2)?.abs();
                    for _ in 0..N_P_ITER_MAX {
                        if f1 < f2 {
                            x3 = x2;
                            x2 = x1;
                            f2 = f1;
                            x1 = x3 - phi * (x3 - x0);
                            f1 = g(x1)?.abs();
                        } else {
                            x0 = x1;
                            x1 = x2;
                            f1 = f2;
                            x2 = x0 + phi * (x3 - x0);
                            f2 = g(x2)?.abs();
                        }
                    }
                    let (p_best, g_best) = if f1 < f2 { (x1, f1) } else { (x2, f2) };
                    if g_best <= self.rho_err_bound {
                        return Ok(p_best);
                    }
                    return Err("state off the equilibrium surface beyond its declared \
                         interpolation-error bound — no admissible projection");
                }
            }
        }
        Self::illinois_root(a, b, ga, gb, &g)
    }

    /// Equilibrium temperature at a primitive state produced by this
    /// occupant (diagnostics + the wall-law operand feed).
    pub fn temperature_w(&self, w: &Prim) -> Result<f64, TableError> {
        let (rho, p, z) = (w[I_RHO], w[4], w[I_RC]);
        let h = w[I_EI] + self.h_offset + p / rho;
        self.temperature.interpolate(&[p, h, z])
    }

    /// Conserved state on the equilibrium surface from `(p, h, Z)` — true
    /// enthalpy — and a velocity (constructor initial conditions and the
    /// COUP-7 injector inflow). **S18 knockdown semantics (fixed session
    /// 12):** `h_offset` shifts the equilibrium *interrogation* coordinate
    /// only — the realized state is the surface at `h + h_offset` — while
    /// the STORED energy is always the true `h − p/ρ`. The deficit
    /// enthalpy is conserved but equilibrium-invisible (sequestered as the
    /// incomplete-combustion energy η_c\* models); a symmetric shift on
    /// both sides is a pure gauge relabeling that changes nothing (the
    /// session-11 wiring — found by the calibration trial coming back
    /// bit-identical).
    pub fn cons_from_phz(&self, p: f64, h: f64, z: f64, vel: [f64; 3]) -> Result<Cons, TableError> {
        let rho = self.rho.interpolate(&[p, h + self.h_offset, z])?;
        let e_true = h - p / rho;
        let ke = 0.5 * (vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]);
        Ok([
            rho,
            rho * vel[0],
            rho * vel[1],
            rho * vel[2],
            rho * (e_true + ke),
            rho * z,
        ])
    }

    /// The declared (p, h, Z) envelopes (assembly-time sizing/refusals).
    pub fn envelopes(&self) -> [(f64, f64); 3] {
        [self.p_env, self.h_env, self.z_env]
    }
}

impl TableEos<'_> {
    fn prim_checked_impl(&self, u: &Cons, hint: Option<f64>) -> Result<Prim, &'static str> {
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
        if z < self.z_env.0 || z > self.z_env.1 {
            return Err("elemental mixture fraction outside the table envelope");
        }
        let e_q = e + self.h_offset;
        let p = match hint {
            Some(ph) if ph.is_finite() && ph > 0.0 => {
                self.project_pressure_hinted(rho, e_q, z, ph)?
            }
            _ => self.project_pressure(rho, e_q, z)?,
        };
        let h = e_q + p * inv;
        let a = self
            .sound
            .interpolate(&[p, h, z])
            .map_err(|_| "sound-speed query failed at the projected state")?;
        let g1 = rho * a * a / p;
        Ok([rho, ur, ut, uz, p, z, e, g1])
    }
}

impl EosLaw for TableEos<'_> {
    fn prim_checked(&self, u: &Cons) -> Result<Prim, &'static str> {
        self.prim_checked_impl(u, None)
    }

    fn prim_checked_hinted(&self, u: &Cons, hint: Option<f64>) -> Result<Prim, &'static str> {
        self.prim_checked_impl(u, hint)
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
        h_roe: f64,
        q2_roe: f64,
        sql: f64,
        sqr: f64,
        inv: f64,
    ) -> f64 {
        // Roe-averaged Γ₁ — the standard general-convex-EOS extension of the
        // Batten bounds (the aux slot is reconstructed like every primitive).
        let g1_roe = (sql * wl[I_G1] + sqr * wr[I_G1]) * inv;
        ((g1_roe - 1.0) * (h_roe - 0.5 * q2_roe)).max(0.0).sqrt()
    }

    fn stagnation_ghost(
        &self,
        _p0: f64,
        _rho0: f64,
        _c_frac: f64,
        _u_n: f64,
        _normal: usize,
    ) -> Result<Prim, FlowError> {
        // No closed-form isentrope on a tabulated surface; the COUP-7
        // injector object owns inflow in shifting mode. Refuse, never
        // approximate (META-1 P6).
        Err(FlowError::BcUnsupportedByEos {
            bc: "StagnationInflow",
        })
    }

    /// COUP-7 §3.2.1 prior-tier inflow: declared (ṁ/A, h_total, Z), interior
    /// static pressure. Fixed-count solve of `ρ = ρ_tab(p, h_total − ½u², Z)`
    /// with `u = (ṁ/A)/ρ` (deterministic, [`N_INFLOW_ITER`]); the resulting
    /// ghost is the premixed equilibrium injection state — combustion
    /// completes at the plane by construction of the equilibrium surface
    /// (the prior tier's declared meaning, SOLV-1 §3.4).
    ///
    /// **Startup regularization (declared):** the face velocity is capped at
    /// the local sound speed — a physical injector face cannot exceed
    /// sonic injection, so while the chamber is filling from near-vacuum
    /// the delivered ṁ is the choked-face value and grows with the interior
    /// pressure. The cap is inactive at any established operating point
    /// (chamber-face Mach ~10⁻², orders below unity), so it never touches a
    /// converged state — a startup path device, not a physics closure.
    fn mass_flow_inflow_ghost(
        &self,
        mdot_per_area: f64,
        h_total: f64,
        c_frac: f64,
        p_int: f64,
        sign: f64,
        normal: usize,
    ) -> Result<Prim, FlowError> {
        let off = |what: &'static str| FlowError::NonPhysicalState {
            i_r: usize::MAX,
            i_z: usize::MAX,
            i_theta: 0,
            what,
        };
        if !(p_int.is_finite() && p_int > 0.0) {
            return Err(off("non-physical interior pressure at the injector face"));
        }
        let mut h_s = h_total;
        let (mut rho, mut u) = (f64::NAN, 0.0f64);
        // Surface queries at the knocked coordinate h + h_offset; stored
        // energy true (see cons_from_phz — the S18 asymmetry).
        for _ in 0..N_INFLOW_ITER {
            rho = self
                .rho
                .interpolate(&[p_int, h_s + self.h_offset, c_frac])
                .map_err(|_| off("injector inflow state outside the table envelope"))?;
            let a = self
                .sound
                .interpolate(&[p_int, h_s + self.h_offset, c_frac])
                .map_err(|_| off("injector inflow sound speed outside the table envelope"))?;
            u = (mdot_per_area / rho).min(a);
            h_s = h_total - 0.5 * u * u;
        }
        let a = self
            .sound
            .interpolate(&[p_int, h_s + self.h_offset, c_frac])
            .map_err(|_| off("injector inflow sound speed outside the table envelope"))?;
        let e_true = h_s - p_int / rho;
        let g1 = rho * a * a / p_int;
        let mut m = [rho, 0.0, 0.0, 0.0, p_int, c_frac, e_true, g1];
        m[normal] = sign * u;
        Ok(m)
    }
}
