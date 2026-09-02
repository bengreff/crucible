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

use super::{Cons, EosLaw, FlowError, I_EI, I_EN, I_G1, I_MR, I_MT, I_MZ, I_RB, I_RC, I_RHO, Prim};
use crucible_tables::{BoundColumn, ColumnMarshal, Table, TableError};

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

/// Fixed bisection budget of the root finder's guaranteed closer (S7):
/// enough to reach `EPS_P_PROJECTION` relative from any admissible
/// bracket (2⁻⁶⁴ of even a seven-decade span is below one ulp).
pub const N_P_BISECT: usize = 64;

/// Fixed sign-change scan resolution of the projection's slow path (the
/// near-vacuum non-monotone corner; see `project_pressure`). Sized by
/// need, not tuned: 16 missed the shallow crossing of a dial-12 RL10
/// establishment shear-layer state whose root demonstrably sits inside
/// the bracket (session 12 — the halt's (ρ, e) solved by hand); 64 gives
/// ~11 points/decade over the widest admissible bracket, matching the
/// table's own p-axis density, so the scan cannot under-resolve a
/// contour crossing the surface itself represents.
pub const N_P_SCAN: usize = 64;

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
    /// The density column's rule-space (|Δln ρ|, ≈ relative) bound —
    /// present on 0.3.2+ artifacts; the acceptance uses it when available
    /// (an absolute bound attained at chamber densities is vacuous at
    /// fringe densities — session-12 review). `None` ⇒ absolute fallback,
    /// preserving pre-0.3.2 artifact behavior bit-for-bit.
    rho_err_bound_log: Option<f64>,
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
            rho_err_bound_log: rho.interp_error_bound_log(),
            rho,
            sound,
            temperature,
            h_offset: 0.0,
        })
    }

    /// Warm-started projection: try a tight bracket around a hint pressure
    /// (the cell's previous-stage projection — the state moves a CFL-limited
    /// fraction per stage, so the root almost always sits within a few
    /// percent). **Root-uniqueness guard (session-12 review):** the tight
    /// bracket is used ONLY when the full admissible bracket's endpoints
    /// straddle the root — the cold fast path's own precondition, under
    /// which the module declares g effectively monotone (unique root), so
    /// the hinted result is the SAME root the cold path finds and the hint
    /// stays a pure acceleration. In the non-monotone near-vacuum corner
    /// (endpoints don't straddle) a tight bracket could lock onto a
    /// different branch than the cold scan — history-dependent physics —
    /// so that case falls through to the cold slow path unconditionally.
    /// Deterministic (data-dependent branch only); bit-reproducible at any
    /// thread count. Measured session 12: the full-bracket Illinois
    /// dominated the march (~70% of wall clock in surface interpolation).
    /// `pub(crate)`: the SOLV-4 §3.6 blend's pure-limit hinted projections
    /// delegate here (S7), mirroring the cold-path delegation.
    pub(crate) fn project_pressure_hinted(
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
        // Uniqueness precondition: the full bracket must straddle. (NaN
        // products fall through to the cold path — never the warm one.)
        let ga_full = g(lo_adm)?;
        let gb_full = g(hi_adm)?;
        let straddles = ga_full * gb_full < 0.0;
        if !straddles {
            return self.project_pressure(rho, e_q, z);
        }
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
        // Straddling full bracket but not the tight one: run the fast path
        // on the full bracket (endpoint values already in hand).
        Self::illinois_root(lo_adm, hi_adm, ga_full, gb_full, &g)
    }

    /// The deterministic Illinois regula-falsi over a sign-changing bracket
    /// (shared by the cold and warm-started projections, and the SOLV-4 §3.6
    /// blend projection in `blend_eos.rs` — one root finder).
    pub(crate) fn illinois_root(
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
        // Guaranteed closer (S7): the bracket is sign-changed by
        // construction, so a fixed bisection CANNOT fail to reach the
        // tolerance — regula-falsi is the fast path, this the deterministic
        // backstop (measured need: a wide blend purge-state bracket where
        // the secant iterate hugged the flat endpoint through the whole
        // budget). Any projection that converged before is bit-unchanged
        // (this path only runs where the old code refused).
        let (mut a, mut b, mut ga) = if a < b { (a, b, ga) } else { (b, a, gb) };
        for _ in 0..N_P_BISECT {
            let m = 0.5 * (a + b);
            if (b - a).abs() <= EPS_P_PROJECTION * a.abs().max(b.abs()) {
                return Ok(m);
            }
            let gm = g(m)?;
            if gm == 0.0 {
                return Ok(m);
            }
            if ga * gm < 0.0 {
                b = m;
            } else {
                a = m;
                ga = gm;
            }
        }
        Ok(0.5 * (a + b))
    }

    /// The equilibrium pressure projection: root of
    /// `ρ_tab(p, h(p), Z) − ρ` over the admissible bracket, where
    /// `h(p) = e_q + p/ρ`, `e_q = e + h_offset`. Deterministic Illinois
    /// regula-falsi (module header). Returns the located pressure.
    /// `pub(crate)`: the SOLV-4 §3.6 blend delegates its pure-limit
    /// projections here (S7), so `b = 0`/`b = 1` cells run THIS routine —
    /// near-vacuum tangency acceptance included — exactly as a pure
    /// single-surface occupant would.
    pub(crate) fn project_pressure(&self, rho: f64, e_q: f64, z: f64) -> Result<f64, &'static str> {
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
                    let mut g1 = g(x1)?; // signed: ρ_tab − ρ
                    let mut g2 = g(x2)?;
                    for _ in 0..N_P_ITER_MAX {
                        if g1.abs() < g2.abs() {
                            x3 = x2;
                            x2 = x1;
                            g2 = g1;
                            x1 = x3 - phi * (x3 - x0);
                            g1 = g(x1)?;
                        } else {
                            x0 = x1;
                            x1 = x2;
                            g1 = g2;
                            x2 = x0 + phi * (x3 - x0);
                            g2 = g(x2)?;
                        }
                    }
                    let (p_best, g_best) = if g1.abs() < g2.abs() {
                        (x1, g1)
                    } else {
                        (x2, g2)
                    };
                    // Acceptance in the density column's own rule space
                    // when the producer stamped it (session-12 review: an
                    // absolute bound attained at chamber ρ accepts states
                    // 10⁴× the local ρ off-surface at the fringe — a
                    // de-facto clamp); absolute fallback preserves
                    // pre-0.3.2 artifacts.
                    let accepted = match self.rho_err_bound_log {
                        Some(lb) => {
                            let ratio = 1.0 + g_best / rho;
                            ratio > 0.0 && ratio.ln().abs() <= lb
                        }
                        None => g_best.abs() <= self.rho_err_bound,
                    };
                    if accepted {
                        return Ok(p_best);
                    }
                    return Err("state off the equilibrium surface beyond its declared \
                         interpolation-error bound — no admissible projection");
                }
            }
        }
        Self::illinois_root(a, b, ga, gb, &g)
    }

    /// **The one interrogation coordinate** of the pinned surfaces at a
    /// primitive state: `(p, h, Z)` with `h` the surface's own enthalpy
    /// coordinate, `h_offset` included (S18 knockdown semantics — see
    /// [`TableEos::cons_from_phz`]). Every consumer that interrogates a
    /// (p, h, Z)-keyed surface at a cell — this occupant's temperature, and
    /// the FND-7 §3.3 transport spine — forms the triple **here**, so the
    /// EOS and the transport surface can never be read at different states
    /// (OFFL-5 §3.1a's envelope-identity requirement is pointless if the
    /// two coordinates drift).
    pub fn interrogation_php(&self, w: &Prim) -> crate::transport::MediumState {
        let (rho, p, z) = (w[I_RHO], w[4], w[I_RC]);
        crate::transport::MediumState {
            p,
            h: w[I_EI] + self.h_offset + p / rho,
            z,
        }
    }

    /// Equilibrium temperature at a primitive state produced by this
    /// occupant (diagnostics + the wall-law operand feed).
    pub fn temperature_w(&self, w: &Prim) -> Result<f64, TableError> {
        let st = self.interrogation_php(w);
        self.temperature.interpolate(&[st.p, st.h, st.z])
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
            0.0, // burn progress: shifting-mode init is the inert b ≡ 0 corner
        ])
    }

    /// The declared (p, h, Z) envelopes (assembly-time sizing/refusals).
    pub fn envelopes(&self) -> [(f64, f64); 3] {
        [self.p_env, self.h_env, self.z_env]
    }

    // --- Raw column queries for the SOLV-4 §3.6 blended occupant ------------
    // `BurnBlendEos` (blend_eos.rs) interrogates *this* surface as one of the
    // two branches; it applies its own `h_offset` when forming the burnt
    // coordinate, so these are RAW (no offset) queries. The blend's mass-
    // weighted specific-volume projection composes the two branches' `ρ`.

    /// Density [kg/m³] at `(p, h, Z)` on this surface (no `h_offset`).
    pub(crate) fn rho_at(&self, p: f64, h: f64, z: f64) -> Result<f64, TableError> {
        self.rho.interpolate(&[p, h, z])
    }
    /// Sound speed [m/s] at `(p, h, Z)` (no `h_offset`).
    pub(crate) fn sound_at(&self, p: f64, h: f64, z: f64) -> Result<f64, TableError> {
        self.sound.interpolate(&[p, h, z])
    }
    /// Temperature [K] at `(p, h, Z)` (no `h_offset`). The blend reads the
    /// unburnt branch here for `T_u` (the SOLV-4 §3.6 rate-law coordinate)
    /// and both branches for its mass-weighted diagnostic temperature (S7).
    pub(crate) fn temp_at(&self, p: f64, h: f64, z: f64) -> Result<f64, TableError> {
        self.temperature.interpolate(&[p, h, z])
    }
    /// The density column's rule-space (log) interpolation bound, if the
    /// producer stamped one — the blend's mid-transition tangency acceptance
    /// composes the two branches' bounds mass-weighted (S7).
    pub(crate) fn rho_bound_log(&self) -> Option<f64> {
        self.rho_err_bound_log
    }
}

/// S13c GPU cross-check for the real (p,h,Z) projection (doc-hidden;
/// additive). Marshals the three surface columns + the (p,h,Z) envelopes for
/// the device `gpu_table_project`, and exposes the oracle projection at one
/// state (the same `project_pressure_hinted`/`project_pressure` +
/// sound/temperature lookups `prim_checked_impl` runs). Changes no production
/// number.
#[doc(hidden)]
pub struct TableEosMarshal {
    pub rho: ColumnMarshal,
    pub sound: ColumnMarshal,
    pub temperature: ColumnMarshal,
    pub p_env: (f64, f64),
    pub h_env: (f64, f64),
    pub z_env: (f64, f64),
}

impl TableEos<'_> {
    #[doc(hidden)]
    pub fn xcheck_marshal(&self) -> TableEosMarshal {
        TableEosMarshal {
            rho: self.rho.marshal(),
            sound: self.sound.marshal(),
            temperature: self.temperature.marshal(),
            p_env: self.p_env,
            h_env: self.h_env,
            z_env: self.z_env,
        }
    }

    /// The oracle projection at one `(ρ, e_q, Z)` state (+ optional pressure
    /// hint): returns `(p, sound, temperature, Γ₁)` — exactly what
    /// `prim_checked_impl` computes, minus the kinematics. `e_q` already
    /// includes any `h_offset`.
    #[doc(hidden)]
    pub fn xcheck_project(
        &self,
        rho: f64,
        e_q: f64,
        z: f64,
        hint: Option<f64>,
    ) -> Result<(f64, f64, f64, f64), &'static str> {
        let inv = 1.0 / rho;
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
        let t = self
            .temperature
            .interpolate(&[p, h, z])
            .map_err(|_| "temperature query failed at the projected state")?;
        Ok((p, a, t, rho * a * a / p))
    }

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
        // Burn progress advects as a passive scalar; the shifting occupant
        // never reads it (a plain shifting run is the inert b ≡ 0 corner —
        // the blended occupant SOLV-4 §3.6 reads it). Carried so the state
        // width is uniform and the projection converts round-trip.
        let b = u[I_RB] * inv;
        if !b.is_finite() {
            return Err("non-finite burn progress");
        }
        Ok([rho, ur, ut, uz, p, z, b, e, g1])
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
        // Session-12 review fix: the gamma-law identity
        // `c² = (Γ₁−1)(h − q²/2)` presumes a ZERO-REFERENCED calorically-
        // perfect enthalpy; on the CEA formation-referenced surface
        // (h ~ −10⁷ J/kg at chamber states) it is hugely negative and the
        // old `.max(0.0)` silently returned c_roe = 0 at EVERY face —
        // degrading the Batten bounds to local-only estimates. The
        // datum-free general-EOS form: Roe-average the squared sound
        // speeds themselves, `c̃² = (√ρ_l·c_l² + √ρ_r·c_r²)/(√ρ_l+√ρ_r)`,
        // with c² = Γ₁·p/ρ from each side's own projected state — always
        // positive, reduces to the classic value for identical states.
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
        // burn = 0: the injected propellant is unburnt (the blended occupant
        // burns it in the chamber); the shifting occupant ignores the slot.
        let mut m = [rho, 0.0, 0.0, 0.0, p_int, c_frac, 0.0, e_true, g1];
        m[normal] = sign * u;
        Ok(m)
    }
}
