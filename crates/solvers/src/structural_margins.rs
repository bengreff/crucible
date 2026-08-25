//! SOLV-6 §3.2 — **analytic quasi-static structural margins** (the v1
//! subset): thin-shell hoop, longitudinal, and through-wall-gradient
//! thermal stresses at a liner shell element, superposed per direction,
//! reduced by the stated von Mises plane-stress criterion, and reported
//! as yield/ultimate margins against T-interpolated allowables.
//!
//! **Ownership (SOLV-6 §1):** this module owns the analytic stress
//! formulas (Roark; META-3 `roark`) and the margin output — nothing
//! else. The p and T operands it consumes come from the solvers
//! (SOLV-1/conduction via the §3.1 structural annotation); the material
//! allowables are FND-7's (here a two-point cited stand-in,
//! [`Allowables`]); the **halt decision is COUP-4's / the engine's** —
//! a margin < 1 is a halt *input*, never acted on here. It is a
//! **constraint check, not a structural simulation** (Failure-Mode
//! Razor, SOLV-6 §0): no FEM, no fatigue, no fracture.
//!
//! **Load convention:** the shell is loaded with whatever pressure the
//! CALLER passes as `p_pa`. The engine (SOLV-6 0.3) passes the DIFFERENCE
//! `|p_gas − p_coolant|` with a declared cited coolant backpressure — a
//! regen liner is loaded by the difference, and an expander jacket runs
//! above chamber pressure; a caller passing absolute internal pressure is
//! making the conservative vacuum-backed choice, declared at ITS site.
//!
//! **Validity (SOLV-6 §5):** thin shell (`2R/t > 20`), elastic,
//! quasi-static. Outside the thin-shell gate the formula is **refused**
//! ([`MarginError::ThinShellViolated`]), never silently applied
//! (META-1 P6). Likewise the allowables refuse interrogation beyond
//! their cited temperature range — extrapolating a handbook strength is
//! guessing.
//!
//! The v1 thermal term is the **through-wall linear-gradient surface
//! stress** — the cooled-liner case (SOLV-6 §3.2, SOLV-6.3). The
//! uniform-constrained biaxial form `E·α·ΔT/(1−ν)` is the annotation's
//! other declared case and arrives with the §3.1 annotation machinery;
//! the uniaxial `E·α·ΔT` is **retired** (doc N8).

/// Yield safety factor (SOLV-6.5, N10): NASA-STD-5012-class factor for
/// liquid-fueled propulsion engine structures [META-3: `nasa-std-5012`].
/// Stated once, here.
pub const FS_YIELD: f64 = 1.1;
/// Ultimate (burst) safety factor (SOLV-6.5, N10): NASA-STD-5012-class
/// factor for liquid-fueled propulsion engine structures
/// [META-3: `nasa-std-5012`]. Stated once, here.
pub const FS_ULT: f64 = 1.4;
/// Thin-shell validity gate (SOLV-6.1, §5): the hoop/longitudinal
/// formulas apply only for `2R/t >` this; below it the check refuses.
pub const THIN_SHELL_MIN_2R_OVER_T: f64 = 20.0;

/// A refused margin evaluation — refuse, never clamp (META-1 P6).
#[derive(Debug, Clone, PartialEq)]
pub enum MarginError {
    /// The thin-shell validity gate `2R/t > 20` failed (SOLV-6.1, §5) —
    /// the formula is not applied outside its idealization.
    ThinShellViolated { two_r_over_t: f64 },
    /// Allowables interrogated at a temperature outside their cited
    /// range `[lo, hi]` — a handbook strength is never extrapolated.
    AllowablesOutOfRange { t_k: f64, lo: f64, hi: f64 },
    /// A non-finite or non-physical (negative where positivity is
    /// physical) operand.
    NonPhysical { what: &'static str },
}

impl std::fmt::Display for MarginError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::ThinShellViolated { two_r_over_t } => write!(
                f,
                "thin-shell validity violated: 2R/t = {two_r_over_t} ≤ \
                 {THIN_SHELL_MIN_2R_OVER_T} — refusing to apply the shell formulas \
                 (SOLV-6 §5, META-1 P6)"
            ),
            Self::AllowablesOutOfRange { t_k, lo, hi } => write!(
                f,
                "allowables interrogated at T = {t_k} K outside their cited range \
                 [{lo}, {hi}] K — refusing to extrapolate a handbook strength \
                 (SOLV-6 §3.2, META-1 P6)"
            ),
            Self::NonPhysical { what } => write!(
                f,
                "non-physical structural operand: {what} — refusing (META-1 P6)"
            ),
        }
    }
}

impl std::error::Error for MarginError {}

/// T-dependent material allowables (SOLV-6 §3.2 last bullet): A/B-basis
/// yield and UTS at two cited temperature points, interpolated linearly
/// between them and **refused** outside them. The FND-7 spine's
/// degradation-state allowable law supersedes this two-point stand-in
/// (Stage-2, SOLV-8); the interface — strength at a temperature, refusal
/// beyond the citation — is the contract that survives.
///
/// Fields are private: an `Allowables` exists only via [`Allowables::new`],
/// which validates, so a held value is always self-consistent.
#[derive(Debug, Clone, Copy)]
pub struct Allowables {
    yield_cold_pa: f64,
    yield_hot_pa: f64,
    uts_cold_pa: f64,
    uts_hot_pa: f64,
    t_cold_k: f64,
    t_hot_k: f64,
    t_solidus_k: f64,
}

impl Allowables {
    /// Validated construction: all strengths positive and finite, the
    /// cited temperature points ordered (`t_cold_k < t_hot_k`), the
    /// solidus positive and finite.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        yield_cold_pa: f64,
        yield_hot_pa: f64,
        uts_cold_pa: f64,
        uts_hot_pa: f64,
        t_cold_k: f64,
        t_hot_k: f64,
        t_solidus_k: f64,
    ) -> Result<Allowables, MarginError> {
        for (what, v) in [
            ("yield_cold_pa", yield_cold_pa),
            ("yield_hot_pa", yield_hot_pa),
            ("uts_cold_pa", uts_cold_pa),
            ("uts_hot_pa", uts_hot_pa),
            ("t_cold_k", t_cold_k),
            ("t_hot_k", t_hot_k),
            ("t_solidus_k", t_solidus_k),
        ] {
            if !(v.is_finite() && v > 0.0) {
                return Err(MarginError::NonPhysical { what });
            }
        }
        if t_cold_k >= t_hot_k {
            return Err(MarginError::NonPhysical {
                what: "allowables temperature points not ordered (t_cold_k < t_hot_k)",
            });
        }
        Ok(Allowables {
            yield_cold_pa,
            yield_hot_pa,
            uts_cold_pa,
            uts_hot_pa,
            t_cold_k,
            t_hot_k,
            t_solidus_k,
        })
    }

    fn interp(&self, cold: f64, hot: f64, t_k: f64) -> Result<f64, MarginError> {
        if !t_k.is_finite() || t_k < self.t_cold_k || t_k > self.t_hot_k {
            return Err(MarginError::AllowablesOutOfRange {
                t_k,
                lo: self.t_cold_k,
                hi: self.t_hot_k,
            });
        }
        let frac = (t_k - self.t_cold_k) / (self.t_hot_k - self.t_cold_k);
        Ok(cold + (hot - cold) * frac)
    }

    /// A/B-basis yield strength [Pa] at `t_k`, linear between the two
    /// cited points; refuses outside `[t_cold_k, t_hot_k]`.
    pub fn yield_at(&self, t_k: f64) -> Result<f64, MarginError> {
        self.interp(self.yield_cold_pa, self.yield_hot_pa, t_k)
    }

    /// A/B-basis ultimate strength [Pa] at `t_k`, linear between the two
    /// cited points; refuses outside `[t_cold_k, t_hot_k]`.
    pub fn uts_at(&self, t_k: f64) -> Result<f64, MarginError> {
        self.interp(self.uts_cold_pa, self.uts_hot_pa, t_k)
    }

    /// Melt predicate: `t_k ≥ t_solidus_k`. The melt **halt** is the
    /// caller's (COUP-4/engine) — this only answers the question.
    pub fn melted(&self, t_k: f64) -> bool {
        t_k >= self.t_solidus_k
    }
}

/// One liner shell element (the SOLV-6 §3.1 annotation's operands, SI):
/// local radius `R`, **real** wall thickness `t`, the shell pressure load
/// `p` (the caller's declared convention — module-doc load note),
/// through-wall
/// `ΔT = T_inner − T_outer` paired along the shell normal, the surface
/// temperature the allowables are interrogated at, and the elastic set
/// `(E, α, ν)`.
#[derive(Debug, Clone, Copy)]
pub struct ShellElement {
    pub r_m: f64,
    pub t_m: f64,
    pub p_pa: f64,
    /// `T_inner − T_outer` [K]; sign carries the load sense (SOLV-6.4) —
    /// a cooled liner (hot inside) has ΔT > 0 and compressive-free
    /// positive surface stress in this convention.
    pub dt_wall_k: f64,
    pub t_surface_k: f64,
    pub e_pa: f64,
    pub alpha_per_k: f64,
    pub nu: f64,
}

/// The margin field one element reports (SOLV-6 §2): the equivalent
/// stress and the two factored margins. Margin < 1 is a COUP-4 halt
/// *input*; this module never halts.
#[derive(Debug, Clone, Copy)]
pub struct MarginReport {
    pub sigma_eq_pa: f64,
    pub yield_margin: f64,
    pub ultimate_margin: f64,
}

/// (SOLV-6.1) thin-shell hoop stress `σ_hoop = p·R/t` [Pa]. Pure
/// algebra; the thin-shell validity gate lives in [`margin_check`].
pub fn sigma_hoop(p_pa: f64, r_m: f64, t_m: f64) -> f64 {
    p_pa * r_m / t_m
}

/// (SOLV-6.2) longitudinal stress `σ_long = p·R/(2t)` [Pa].
pub fn sigma_long(p_pa: f64, r_m: f64, t_m: f64) -> f64 {
    p_pa * r_m / (2.0 * t_m)
}

/// (SOLV-6.3) through-wall **linear-gradient** thermal surface stress
/// `σ_therm = E·α·ΔT/(2(1−ν))` [Pa] (Roark; the cooled-liner case,
/// equal in hoop and longitudinal). The retired uniaxial `E·α·ΔT`
/// mis-states this case by the factor `2(1−ν)` (doc N8).
pub fn sigma_thermal_gradient(e_pa: f64, alpha_per_k: f64, dt_wall_k: f64, nu: f64) -> f64 {
    e_pa * alpha_per_k * dt_wall_k / (2.0 * (1.0 - nu))
}

/// (SOLV-6.4) von Mises plane-stress equivalent
/// `σ_eq = √(σ_θ² − σ_θ·σ_z + σ_z²)` [Pa] — the stated criterion;
/// never a scalar sum of orthogonal components.
pub fn von_mises_plane(sigma_theta: f64, sigma_z: f64) -> f64 {
    (sigma_theta * sigma_theta - sigma_theta * sigma_z + sigma_z * sigma_z).sqrt()
}

/// The SOLV-6 §3.2 chain for one shell element: per-direction
/// superposition (SOLV-6.4) of the pressure (SOLV-6.1/6.2) and
/// linear-gradient thermal (SOLV-6.3) stresses, von Mises reduction,
/// and the factored margins (SOLV-6.5)
/// `margin = allowable(T) / (FS · σ_eq)`.
///
/// Refuses (never clamps): non-finite/non-physical operands, a
/// thin-shell violation (`2R/t ≤ 20`), allowables outside their cited
/// range. An unloaded element (`σ_eq = 0`) reports infinite margins —
/// honest, and the caller's melt/halt logic is unaffected.
pub fn margin_check(inp: &ShellElement, allow: &Allowables) -> Result<MarginReport, MarginError> {
    // Positivity where physical; `p_pa` may be 0 (unpressurized), and
    // `dt_wall_k` is signed by convention (see `ShellElement`).
    for (what, v, strictly_positive) in [
        ("r_m", inp.r_m, true),
        ("t_m", inp.t_m, true),
        ("p_pa", inp.p_pa, false),
        ("dt_wall_k", inp.dt_wall_k.abs(), false),
        ("t_surface_k", inp.t_surface_k, true),
        ("e_pa", inp.e_pa, true),
        ("alpha_per_k", inp.alpha_per_k, false),
    ] {
        if !(v.is_finite() && (v > 0.0 || (!strictly_positive && v >= 0.0))) {
            return Err(MarginError::NonPhysical { what });
        }
    }
    // ν on the physical open interval where the plane form is defined;
    // negative (auxetic) liners are outside this module's citation.
    if !(inp.nu.is_finite() && (0.0..0.5).contains(&inp.nu)) {
        return Err(MarginError::NonPhysical { what: "nu" });
    }
    let two_r_over_t = 2.0 * inp.r_m / inp.t_m;
    if two_r_over_t <= THIN_SHELL_MIN_2R_OVER_T {
        return Err(MarginError::ThinShellViolated { two_r_over_t });
    }
    let s_therm = sigma_thermal_gradient(inp.e_pa, inp.alpha_per_k, inp.dt_wall_k, inp.nu);
    let s_theta = sigma_hoop(inp.p_pa, inp.r_m, inp.t_m) + s_therm;
    let s_z = sigma_long(inp.p_pa, inp.r_m, inp.t_m) + s_therm;
    let sigma_eq_pa = von_mises_plane(s_theta, s_z);
    let yield_allow = allow.yield_at(inp.t_surface_k)?;
    let uts_allow = allow.uts_at(inp.t_surface_k)?;
    Ok(MarginReport {
        sigma_eq_pa,
        yield_margin: yield_allow / (FS_YIELD * sigma_eq_pa),
        ultimate_margin: uts_allow / (FS_ULT * sigma_eq_pa),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A round-number allowables set: yield 400→200 MPa and UTS
    /// 800→400 MPa over 300→900 K, solidus 1560 K.
    fn allow() -> Allowables {
        Allowables::new(4.0e8, 2.0e8, 8.0e8, 4.0e8, 300.0, 900.0, 1560.0)
            .expect("valid allowables fixture")
    }

    fn element() -> ShellElement {
        ShellElement {
            r_m: 1.0,
            t_m: 0.05, // 2R/t = 40 — thin
            p_pa: 2.0e6,
            dt_wall_k: 100.0,
            t_surface_k: 600.0,
            e_pa: 2.0e11,
            alpha_per_k: 1.6e-5,
            nu: 0.3,
        }
    }

    /// (SOLV-6.1)/(SOLV-6.2) — validation plan item 1: the closed-form
    /// pressure-vessel stresses, exactly, for round numbers.
    #[test]
    fn pressure_vessel_closed_forms() {
        // p = 2 MPa, R = 1 m, t = 0.05 m by independent arithmetic:
        // hoop = 2e6·1/0.05 = 4e7 Pa; long = 2e6·1/(2·0.05) = 2e7 Pa.
        assert_eq!(sigma_hoop(2.0e6, 1.0, 0.05), 4.0e7);
        assert_eq!(sigma_long(2.0e6, 1.0, 0.05), 2.0e7);
    }

    /// (SOLV-6.3) — validation plan item 2: the Roark linear-gradient
    /// surface stress, and the retired uniaxial form off by 2(1−ν).
    #[test]
    fn thermal_gradient_closed_form_and_retired_uniaxial_factor() {
        // E = 200 GPa, α = 1.6e-5 /K, ΔT = 100 K, ν = 0.3:
        // E·α·ΔT = 2e11·1.6e-5·100 = 3.2e8; /(2·0.7) = 3.2e8/1.4.
        let s = sigma_thermal_gradient(2.0e11, 1.6e-5, 100.0, 0.3);
        assert_eq!(s, 2.0e11 * 1.6e-5 * 100.0 / (2.0 * (1.0 - 0.3)));
        // The retired uniaxial E·α·ΔT differs by exactly the documented
        // factor 2(1−ν) (doc N8).
        let uniaxial = 2.0e11 * 1.6e-5 * 100.0;
        assert!((uniaxial / s - 2.0 * (1.0 - 0.3)).abs() < 1.0e-12);
    }

    /// (SOLV-6.4) — validation plan item 3: hand-computed von Mises for
    /// a combined hoop+thermal case.
    #[test]
    fn von_mises_combined_by_hand() {
        // σ_θ = 4e7 + 1e7 = 5e7; σ_z = 2e7 + 1e7 = 3e7 (independent
        // arithmetic): σ_eq = √(25e14 − 15e14 + 9e14) = √(1.9e15).
        let a: f64 = 5.0e7;
        let b: f64 = 3.0e7;
        let hand = (a * a - a * b + b * b).sqrt();
        assert_eq!(von_mises_plane(5.0e7, 3.0e7), hand);
        assert!((hand - 1.9e15_f64.sqrt()).abs() < 1.0e-3);
    }

    /// (SOLV-6.5): a case engineered so the yield margin is exactly 2 —
    /// allowable = 2·FS_YIELD·σ_eq at the interrogation temperature.
    #[test]
    fn yield_margin_exactly_two() {
        let inp = ShellElement {
            dt_wall_k: 0.0, // pressure only: σ_θ = 2e7, σ_z = 1e7
            p_pa: 1.0e6,
            ..element()
        };
        // Independent σ_eq: √(4e14 − 2e14 + 1e14) = 1e7·√3.
        let a: f64 = 2.0e7;
        let b: f64 = 1.0e7;
        let sigma_eq = (a * a - a * b + b * b).sqrt();
        let denom = FS_YIELD * sigma_eq;
        // Interrogate exactly at the cold point so yield_at is the cold
        // value bit-for-bit (interpolation fraction is exactly 0).
        let allow = Allowables::new(2.0 * denom, denom, 8.0e8, 4.0e8, 300.0, 900.0, 1560.0)
            .expect("valid allowables");
        let rep = margin_check(
            &ShellElement {
                t_surface_k: 300.0,
                ..inp
            },
            &allow,
        )
        .expect("margin check succeeds");
        assert_eq!(rep.sigma_eq_pa, sigma_eq);
        assert_eq!(rep.yield_margin, 2.0);
        assert!(rep.ultimate_margin > 0.0);
    }

    /// (SOLV-6.1) validity — validation plan item 5's refusal half:
    /// 2R/t = 10 is not a thin shell; the formula is refused, never
    /// silently applied.
    #[test]
    fn thick_shell_refuses() {
        let inp = ShellElement {
            t_m: 0.2, // 2R/t = 10
            ..element()
        };
        match margin_check(&inp, &allow()) {
            Err(MarginError::ThinShellViolated { two_r_over_t }) => {
                assert_eq!(two_r_over_t, 10.0);
            }
            other => panic!("expected ThinShellViolated, got {other:?}"),
        }
    }

    /// Allowables: range refusal, melt predicate, and the linear
    /// midpoint (independent arithmetic: (400+200)/2 = 300 MPa yield,
    /// (800+400)/2 = 600 MPa UTS at 600 K).
    #[test]
    fn allowables_range_melt_and_midpoint() {
        let a = allow();
        assert!(matches!(
            a.yield_at(250.0),
            Err(MarginError::AllowablesOutOfRange {
                t_k, lo, hi
            }) if t_k == 250.0 && lo == 300.0 && hi == 900.0
        ));
        assert!(a.uts_at(901.0).is_err());
        assert!(a.yield_at(f64::NAN).is_err());
        // Endpoints are in range; midpoint is the arithmetic mean.
        assert_eq!(a.yield_at(300.0).expect("cold endpoint"), 4.0e8);
        assert_eq!(a.yield_at(600.0).expect("midpoint"), 3.0e8);
        assert_eq!(a.uts_at(600.0).expect("midpoint"), 6.0e8);
        assert_eq!(a.yield_at(900.0).expect("hot endpoint"), 2.0e8);
        // Melt is the caller's check; the predicate is ≥ solidus.
        assert!(!a.melted(1559.9));
        assert!(a.melted(1560.0));
        assert!(a.melted(2000.0));
        // Construction refusals: unordered T points, non-positive strength.
        assert!(Allowables::new(4.0e8, 2.0e8, 8.0e8, 4.0e8, 900.0, 300.0, 1560.0).is_err());
        assert!(Allowables::new(-1.0, 2.0e8, 8.0e8, 4.0e8, 300.0, 900.0, 1560.0).is_err());
    }

    /// Non-physical operands refuse (META-1 P6), and the margins react
    /// to the load in the right direction.
    #[test]
    fn refusals_and_monotonicity() {
        let good = element();
        for bad in [
            ShellElement { r_m: -1.0, ..good },
            ShellElement {
                p_pa: f64::NAN,
                ..good
            },
            ShellElement { e_pa: 0.0, ..good },
            ShellElement { nu: 0.5, ..good },
            ShellElement {
                dt_wall_k: f64::INFINITY,
                ..good
            },
        ] {
            assert!(matches!(
                margin_check(&bad, &allow()),
                Err(MarginError::NonPhysical { .. })
            ));
        }
        let base = margin_check(&good, &allow()).expect("base case");
        let hotter_wall = margin_check(
            &ShellElement {
                dt_wall_k: 200.0,
                ..good
            },
            &allow(),
        )
        .expect("hotter wall");
        assert!(hotter_wall.sigma_eq_pa > base.sigma_eq_pa);
        assert!(hotter_wall.yield_margin < base.yield_margin);
        assert!(hotter_wall.ultimate_margin < base.ultimate_margin);
    }
}
