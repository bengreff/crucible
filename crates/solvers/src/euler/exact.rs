//! Exact Riemann solver for the gamma-law gas — the Station-1 **oracle**
//! (VAL-2 §3.3 `sod-shock`; META-3 `sod-shock` — Sod, *JCP* **27** (1978);
//! Toro ch. 4 exact solver). Runtime physics never calls this; it exists so
//! the certificate compares the captured waves against the exact
//! shock/contact/rarefaction structure.
//!
//! Determinism (META-1 §2.2): the star-pressure Newton iteration uses a
//! fixed absolute+relative tolerance and a fixed max iteration count; not
//! converging is a panic (an oracle defect, never a tolerated event).

/// One side of the initial discontinuity: (ρ, u, p) with u the velocity
/// toward +x.
#[derive(Debug, Clone, Copy)]
pub struct RiemannSide {
    pub rho: f64,
    pub u: f64,
    pub p: f64,
}

/// The solved wave structure: star-region pressure/velocity and the
/// densities either side of the contact.
#[derive(Debug, Clone, Copy)]
pub struct RiemannSolution {
    pub gamma: f64,
    pub left: RiemannSide,
    pub right: RiemannSide,
    pub p_star: f64,
    pub u_star: f64,
    pub rho_star_l: f64,
    pub rho_star_r: f64,
}

const NEWTON_TOL: f64 = 1e-13;
const NEWTON_MAX_ITERS: usize = 100;

/// Toro's f_K(p): velocity change across the left/right nonlinear wave, and
/// its derivative. Shock branch for p > p_K, rarefaction otherwise.
fn f_and_df(p: f64, s: &RiemannSide, gamma: f64) -> (f64, f64) {
    let c = (gamma * s.p / s.rho).sqrt();
    if p > s.p {
        // Shock (Rankine-Hugoniot).
        let a = 2.0 / ((gamma + 1.0) * s.rho);
        let b = (gamma - 1.0) / (gamma + 1.0) * s.p;
        let q = (a / (p + b)).sqrt();
        ((p - s.p) * q, q * (1.0 - 0.5 * (p - s.p) / (p + b)))
    } else {
        // Rarefaction (isentrope).
        let ex = (gamma - 1.0) / (2.0 * gamma);
        (
            2.0 * c / (gamma - 1.0) * ((p / s.p).powf(ex) - 1.0),
            (p / s.p).powf(-(gamma + 1.0) / (2.0 * gamma)) / (s.rho * c),
        )
    }
}

/// Solve the Riemann problem. Panics on vacuum-generating data or a
/// non-converged iteration (oracle defects).
pub fn solve(left: RiemannSide, right: RiemannSide, gamma: f64) -> RiemannSolution {
    let cl = (gamma * left.p / left.rho).sqrt();
    let cr = (gamma * right.p / right.rho).sqrt();
    let du = right.u - left.u;
    assert!(
        2.0 * (cl + cr) / (gamma - 1.0) > du,
        "vacuum-generating Riemann data — outside the oracle's scope"
    );

    // Two-rarefaction initial guess (Toro eq. 4.46), positive by
    // construction; Newton on f(p) = f_L + f_R + Δu.
    let z = (gamma - 1.0) / (2.0 * gamma);
    let mut p = ((cl + cr - 0.5 * (gamma - 1.0) * du)
        / (cl / left.p.powf(z) + cr / right.p.powf(z)))
    .powf(1.0 / z);
    let mut converged = false;
    for _ in 0..NEWTON_MAX_ITERS {
        let (fl, dfl) = f_and_df(p, &left, gamma);
        let (fr, dfr) = f_and_df(p, &right, gamma);
        let f = fl + fr + du;
        let step = f / (dfl + dfr);
        let p_new = (p - step).max(NEWTON_TOL * p);
        let change = (p_new - p).abs();
        p = p_new;
        if change <= NEWTON_TOL * p.max(1.0) {
            converged = true;
            break;
        }
    }
    assert!(converged, "star-pressure iteration did not converge");

    let (fl, _) = f_and_df(p, &left, gamma);
    let (fr, _) = f_and_df(p, &right, gamma);
    let u_star = 0.5 * (left.u + right.u) + 0.5 * (fr - fl);

    let star_density = |s: &RiemannSide| -> f64 {
        let g = (gamma - 1.0) / (gamma + 1.0);
        if p > s.p {
            // Shock compression.
            s.rho * (p / s.p + g) / (g * p / s.p + 1.0)
        } else {
            // Isentropic expansion.
            s.rho * (p / s.p).powf(1.0 / gamma)
        }
    };
    RiemannSolution {
        gamma,
        left,
        right,
        p_star: p,
        u_star,
        rho_star_l: star_density(&left),
        rho_star_r: star_density(&right),
    }
}

impl RiemannSolution {
    /// Sample (ρ, u, p) at similarity coordinate ξ = x/t (Toro §4.5): the
    /// full shock/contact/rarefaction structure, including fan interiors.
    pub fn sample(&self, xi: f64) -> (f64, f64, f64) {
        let g = self.gamma;
        if xi <= self.u_star {
            // Left of the contact.
            let s = &self.left;
            let c = (g * s.p / s.rho).sqrt();
            if self.p_star > s.p {
                // Left shock.
                let sp = s.u
                    - c * ((g + 1.0) / (2.0 * g) * self.p_star / s.p + (g - 1.0) / (2.0 * g))
                        .sqrt();
                if xi <= sp {
                    (s.rho, s.u, s.p)
                } else {
                    (self.rho_star_l, self.u_star, self.p_star)
                }
            } else {
                // Left rarefaction: head, interior fan, tail.
                let c_star = c * (self.p_star / s.p).powf((g - 1.0) / (2.0 * g));
                let head = s.u - c;
                let tail = self.u_star - c_star;
                if xi <= head {
                    (s.rho, s.u, s.p)
                } else if xi >= tail {
                    (self.rho_star_l, self.u_star, self.p_star)
                } else {
                    let u = 2.0 / (g + 1.0) * (c + 0.5 * (g - 1.0) * s.u + xi);
                    let cf = 2.0 / (g + 1.0) * (c + 0.5 * (g - 1.0) * (s.u - xi));
                    let rho = s.rho * (cf / c).powf(2.0 / (g - 1.0));
                    let p = s.p * (cf / c).powf(2.0 * g / (g - 1.0));
                    (rho, u, p)
                }
            }
        } else {
            // Right of the contact (mirror).
            let s = &self.right;
            let c = (g * s.p / s.rho).sqrt();
            if self.p_star > s.p {
                // Right shock.
                let sp = s.u
                    + c * ((g + 1.0) / (2.0 * g) * self.p_star / s.p + (g - 1.0) / (2.0 * g))
                        .sqrt();
                if xi >= sp {
                    (s.rho, s.u, s.p)
                } else {
                    (self.rho_star_r, self.u_star, self.p_star)
                }
            } else {
                let c_star = c * (self.p_star / s.p).powf((g - 1.0) / (2.0 * g));
                let head = s.u + c;
                let tail = self.u_star + c_star;
                if xi >= head {
                    (s.rho, s.u, s.p)
                } else if xi <= tail {
                    (self.rho_star_r, self.u_star, self.p_star)
                } else {
                    let u = 2.0 / (g + 1.0) * (-c + 0.5 * (g - 1.0) * s.u + xi);
                    let cf = 2.0 / (g + 1.0) * (c - 0.5 * (g - 1.0) * (s.u - xi));
                    let rho = s.rho * (cf / c).powf(2.0 / (g - 1.0));
                    let p = s.p * (cf / c).powf(2.0 * g / (g - 1.0));
                    (rho, u, p)
                }
            }
        }
    }

    /// Right-running shock speed (valid when the right wave is a shock —
    /// the Sod configuration).
    pub fn right_shock_speed(&self) -> f64 {
        let g = self.gamma;
        let s = &self.right;
        let c = (g * s.p / s.rho).sqrt();
        assert!(self.p_star > s.p, "right wave is not a shock");
        s.u + c * ((g + 1.0) / (2.0 * g) * self.p_star / s.p + (g - 1.0) / (2.0 * g)).sqrt()
    }
}

// Covers VAL-2 §3.3 Sod oracle values — named per META-2 §4.
#[cfg(test)]
mod tests {
    use super::*;

    /// The canonical Sod data (META-3 `sod-shock`).
    fn sod() -> RiemannSolution {
        solve(
            RiemannSide {
                rho: 1.0,
                u: 0.0,
                p: 1.0,
            },
            RiemannSide {
                rho: 0.125,
                u: 0.0,
                p: 0.1,
            },
            1.4,
        )
    }

    #[test]
    fn val2_s33_sod_star_state_matches_toro_table_4_2() {
        // Toro (3rd ed.) Table 4.2, test 1: p* = 0.30313, u* = 0.92745,
        // ρ*L = 0.42632, ρ*R = 0.26557 (5 significant figures).
        let s = sod();
        assert!((s.p_star - 0.30313).abs() < 5e-6, "p* = {}", s.p_star);
        assert!((s.u_star - 0.92745).abs() < 5e-6, "u* = {}", s.u_star);
        assert!(
            (s.rho_star_l - 0.42632).abs() < 5e-6,
            "ρ*L = {}",
            s.rho_star_l
        );
        assert!(
            (s.rho_star_r - 0.26557).abs() < 5e-6,
            "ρ*R = {}",
            s.rho_star_r
        );
    }

    #[test]
    fn val2_s33_sampling_is_consistent_at_the_waves() {
        let s = sod();
        // Far left: undisturbed left state; far right: undisturbed right.
        assert_eq!(s.sample(-10.0), (1.0, 0.0, 1.0));
        assert_eq!(s.sample(10.0), (0.125, 0.0, 0.1));
        // Just left of the contact: (ρ*L, u*, p*).
        let (rho, u, p) = s.sample(s.u_star - 1e-9);
        assert!((rho - s.rho_star_l).abs() < 1e-6);
        assert!((u - s.u_star).abs() < 1e-9);
        assert!((p - s.p_star).abs() < 1e-9);
        // Pressure and velocity are continuous across the contact.
        let (_, u2, p2) = s.sample(s.u_star + 1e-9);
        assert!((u2 - u).abs() < 1e-8);
        assert!((p2 - p).abs() < 1e-8);
    }
}
