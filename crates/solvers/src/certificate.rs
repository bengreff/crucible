//! The Goal-A **convergence certificate** studies (VAL-1 ladder, rung
//! "analytic"; VAL-3 §3.2 MMS gate). Each runner returns structured data;
//! `tests/` asserts the pass criteria and `bin/convergence_certificate`
//! renders the committed artifact — both from the SAME named constants
//! below, so the artifact can never assert criteria the tests don't
//! enforce. Everything here is deterministic: fixed grids, fixed step
//! counts, no RNG, no wall-clock.

use crate::conduction::{Bcs, Conduction, Domain, FaceBc, InteriorFaces};
use crucible_grid::{Grid, GridSpec};

pub const FIELDS: &[&str] = &["T", "rate"];

// --- Named certificate constants (META-2 §4: no magic numbers) --------------

/// Observed-order acceptance band around the formal order 2 (VAL-3 §3.2:
/// "observed order within ~5–10% of formal"; ±0.2 is the 10% band).
pub const MMS_ORDER_MIN: f64 = 1.8;
pub const MMS_ORDER_MAX: f64 = 2.2;

/// MMS refinement ladder and step base: level ℓ runs `MMS_BASE_STEPS·4^ℓ`
/// steps to `MMS_T_FINAL`, so dt ∝ h² and the O(dt) Euler error refines at
/// the same 2nd-order rate as space. At every level dt sits at ≈0.57× the
/// spectral stability limit (asserted at run time against `stable_dt`).
pub const MMS_LEVELS: [usize; 3] = [8, 16, 32];
pub const MMS_BASE_STEPS: usize = 25;
pub const MMS_T_FINAL: f64 = 0.05;

/// March safety factor for anchor runs: 0.4× the spectral limit also sits
/// below the ~0.66× discrete-maximum-principle threshold, so anchor data is
/// free of bounded transient overshoot, not just divergence.
pub const CFL_SAFETY: f64 = 0.4;

/// Steady annulus anchor tolerance, relative to ΔT: the measured 2nd-order
/// discretization error at 32 radial cells is ≈4.3e-4 (see the committed
/// certificate); 2e-3 gives ~5× regression headroom while still failing on
/// any loss of 2nd-order boundary treatment.
pub const ANNULUS_TOL_REL: f64 = 2e-3;
pub const ANNULUS_CELLS: usize = 32;

/// Bessel-cylinder anchor tolerance [K] on T₀ = 100 K: spatial error at
/// 48 cells is ≈6e-3 K and series truncation ≈1e-8 K, so 0.5 K (0.5%)
/// passes with ~80× margin yet catches any axis-treatment defect, which
/// produces O(1) K errors at the centerline.
pub const BESSEL_TOL_K: f64 = 0.5;

/// Closed-sweep conservation tolerance: FND-2 §6-2 asks ~1e-10 relative;
/// flux-form telescoping actually delivers round-off (~1e-16), so 1e-12
/// enforces 4 decades tighter than the doc's requirement.
pub const CONSERVATION_TOL_REL: f64 = 1e-12;
pub const CONSERVATION_STEPS: usize = 500;

/// One refinement level of an MMS study: grid spacing and L2 error.
#[derive(Debug, Clone, PartialEq)]
pub struct MmsLevel {
    pub n: usize,
    pub h: f64,
    pub l2_error: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct MmsStudy {
    pub label: &'static str,
    pub levels: Vec<MmsLevel>,
}

impl MmsStudy {
    /// Observed order between consecutive levels: p = ln(e_c/e_f)/ln(2).
    pub fn observed_orders(&self) -> Vec<f64> {
        self.levels
            .windows(2)
            .map(|w| (w[0].l2_error / w[1].l2_error).ln() / std::f64::consts::LN_2)
            .collect()
    }
}

fn annulus_spec(n: usize, n_theta: u32, axisym: bool) -> GridSpec {
    GridSpec {
        r_min: 0.5,
        dr: 1.0 / n as f64,
        n_r: n,
        z_min: 0.0,
        dz: 1.0 / n as f64,
        n_z: n,
        n_theta_max: n_theta,
        axisymmetry_assertion: axisym,
    }
}

/// Volume-weighted L2 norm of (field − exact), in the canonical traversal
/// order (fixed-shape reduction; bit-reproducible).
fn l2_error<F: Fn(f64, f64, f64) -> f64>(g: &Grid, f_id: crucible_grid::FieldId, exact: F) -> f64 {
    let mut num = 0.0f64;
    let mut den = 0.0f64;
    g.for_each_active_cell(|c| {
        let b = g.brick(c.bi);
        let v = g.cell_volume(c.i_r, b.n_theta());
        let d = b.field(f_id)[c.idx] - exact(c.r, c.theta, c.z);
        num += v * d * d;
        den += v;
    });
    (num / den).sqrt()
}

/// 2-D axisymmetric MMS (N_θ = 1 under the recorded assertion): the
/// manufactured solution `T = e^{−λt}·cos(αr)·cos(βz)` implies the source
/// `S = ρc_p ∂T/∂t − k∇²T` with the exact cylindrical Laplacian
/// `∇²T = [−α²cos(αr) − (α/r)sin(αr)]cos(βz) − β²cos(αr)cos(βz)` (× e^{−λt});
/// feeding S back must reproduce T to 2nd order in h.
pub fn mms_axisymmetric() -> MmsStudy {
    mms_run("MMS 2-D axisymmetric (N_θ = 1 asserted)", &MMS_LEVELS, 1, 0)
}

/// 3-D MMS with an m = 2 azimuthal mode on the same annulus — exercises the
/// θ-flux path at matching θ refinement (N_θ = n).
pub fn mms_theta_mode() -> MmsStudy {
    mms_run("MMS 3-D with m = 2 θ-mode", &MMS_LEVELS, 0, 2)
}

fn mms_run(label: &'static str, ns: &[usize], n_theta_fixed: u32, m: u32) -> MmsStudy {
    let (kappa, rho_cp, lambda) = (1.0f64, 1.0f64, 1.0f64);
    let (alpha, beta) = (std::f64::consts::PI, std::f64::consts::PI);
    let mf = f64::from(m);

    let spatial = move |r: f64, theta: f64, z: f64| -> f64 {
        (alpha * r).cos() * (beta * z).cos() * (1.0 + 0.5 * (mf * theta).cos())
    };
    let exact = move |r: f64, theta: f64, z: f64, t: f64| -> f64 {
        (-lambda * t).exp() * spatial(r, theta, z)
    };
    let laplacian = move |r: f64, theta: f64, z: f64, t: f64| -> f64 {
        let ang = 1.0 + 0.5 * (mf * theta).cos();
        let radial = (-alpha * alpha * (alpha * r).cos() - alpha / r * (alpha * r).sin())
            * (beta * z).cos()
            * ang;
        let axial = -beta * beta * (alpha * r).cos() * (beta * z).cos() * ang;
        let azimuthal =
            -(mf * mf) / (r * r) * 0.5 * (mf * theta).cos() * (alpha * r).cos() * (beta * z).cos();
        (-lambda * t).exp() * (radial + axial + azimuthal)
    };
    let source = move |r: f64, theta: f64, z: f64, t: f64| -> f64 {
        -lambda * rho_cp * exact(r, theta, z, t) - kappa * laplacian(r, theta, z, t)
    };
    let dirichlet = move |r: f64, theta: f64, z: f64, t: f64| -> f64 { exact(r, theta, z, t) };

    let mut levels = Vec::new();
    for (lvl, &n) in ns.iter().enumerate() {
        let n_theta = if n_theta_fixed > 0 {
            n_theta_fixed
        } else {
            n as u32
        };
        let spec = annulus_spec(n, n_theta, n_theta == 1);
        let mut g = Grid::build(spec, FIELDS).expect("valid spec");
        let t_id = g.field_id("T").unwrap();
        let rate_id = g.field_id("rate").unwrap();
        g.fill_field(t_id, |r, th, z| exact(r, th, z, 0.0));

        let op = Conduction {
            kappa,
            rho_cp,
            source: &source,
            domain: Domain::FlowActive,
            interior: InteriorFaces::refuse(),
            bcs: Bcs {
                r_inner: FaceBc::Dirichlet(&dirichlet),
                r_outer: FaceBc::Dirichlet(&dirichlet),
                z_lo: FaceBc::Dirichlet(&dirichlet),
                z_hi: FaceBc::Dirichlet(&dirichlet),
            },
        };
        let n_steps = MMS_BASE_STEPS * 4usize.pow(lvl as u32);
        let dt = MMS_T_FINAL / n_steps as f64;
        // Guard against the true spectral limit (factor 1.0 — the bound
        // itself, review-verified sharp), not a safety-scaled one.
        assert!(
            dt < op.stable_dt(&g, 1.0),
            "certificate step must be stable"
        );
        op.advance(&mut g, t_id, rate_id, 0.0, dt, n_steps)
            .expect("advance");

        levels.push(MmsLevel {
            n,
            h: 1.0 / n as f64,
            l2_error: l2_error(&g, t_id, |r, th, z| exact(r, th, z, MMS_T_FINAL)),
        });
    }
    MmsStudy { label, levels }
}

/// Steady annulus anchor: inner wall 500 K, outer wall 300 K, insulated z.
/// Exact steady profile `T(r) = (T1·ln(r2/r) + T2·ln(r/r1))/ln(r2/r1)`.
/// Returns (max relative error, grid). Marches to steady state.
pub fn annulus_anchor() -> (f64, Grid) {
    let (t1, t2) = (500.0f64, 300.0f64);
    let (r1, r2) = (0.05f64, 0.15f64);
    let spec = GridSpec {
        r_min: r1,
        dr: (r2 - r1) / ANNULUS_CELLS as f64,
        n_r: ANNULUS_CELLS,
        z_min: 0.0,
        dz: 0.05,
        n_z: 2,
        n_theta_max: 1,
        axisymmetry_assertion: true,
    };
    let mut g = Grid::build(spec, FIELDS).expect("valid spec");
    let t_id = g.field_id("T").unwrap();
    let rate_id = g.field_id("rate").unwrap();
    g.fill_field(t_id, |_, _, _| 0.5 * (t1 + t2));

    let inner = move |_r: f64, _th: f64, _z: f64, _t: f64| t1;
    let outer = move |_r: f64, _th: f64, _z: f64, _t: f64| t2;
    let zero_src = |_: f64, _: f64, _: f64, _: f64| 0.0;
    let op = Conduction {
        kappa: 20.0,
        rho_cp: 4.0e6,
        source: &zero_src,
        domain: Domain::FlowActive,
        interior: InteriorFaces::refuse(),
        bcs: Bcs {
            r_inner: FaceBc::Dirichlet(&inner),
            r_outer: FaceBc::Dirichlet(&outer),
            z_lo: FaceBc::HeatFlux(0.0),
            z_hi: FaceBc::HeatFlux(0.0),
        },
    };
    let dt = op.stable_dt(&g, CFL_SAFETY);
    // ~8 diffusion times across the gap: transients decay like e^{-t/τ}, so
    // 8τ leaves relative transient content ~e⁻⁸ ≈ 3e-4 of the initial
    // offset — an order below the anchor tolerance.
    let tau = (r2 - r1) * (r2 - r1) * op.rho_cp / op.kappa;
    let n_steps = (8.0 * tau / dt).ceil() as usize;
    op.advance(&mut g, t_id, rate_id, 0.0, dt, n_steps)
        .expect("advance");

    let exact = move |r: f64| (t1 * (r2 / r).ln() + t2 * (r / r1).ln()) / (r2 / r1).ln();
    let mut worst = 0.0f64;
    for i_r in 0..ANNULUS_CELLS {
        let r = g.r_center(i_r);
        let bi = g.brick_index(i_r, 0).unwrap();
        let b = g.brick(bi);
        let t = b.field(t_id)[b.cell_index(0, Grid::local_rz(i_r, 0))];
        worst = worst.max((t - exact(r)).abs() / (t1 - t2));
    }
    (worst, g)
}

/// Transient full-cylinder anchor (exercises the r = 0 axis): uniform
/// initial T₀, surface held at 0, compared against the Bessel series
/// `T = T₀·Σ 2/(λ_n J₁(λ_n))·J₀(λ_n r/R)·e^{−λ_n² t̃}` at t̃ = κ̃t/R² = 0.1.
/// Returns the max absolute error [K] (T₀ = 100 K).
pub fn bessel_cylinder_anchor() -> f64 {
    let t0 = 100.0f64;
    let radius = 0.1f64;
    let n_r = 48usize;
    let spec = GridSpec {
        r_min: 0.0,
        dr: radius / n_r as f64,
        n_r,
        z_min: 0.0,
        dz: 0.05,
        n_z: 2,
        n_theta_max: 4,
        axisymmetry_assertion: false,
    };
    let mut g = Grid::build(spec, FIELDS).expect("valid spec");
    let t_id = g.field_id("T").unwrap();
    let rate_id = g.field_id("rate").unwrap();
    g.fill_field(t_id, |_, _, _| t0);

    let cold = |_: f64, _: f64, _: f64, _: f64| 0.0;
    let zero_src = |_: f64, _: f64, _: f64, _: f64| 0.0;
    let op = Conduction {
        kappa: 1.0,
        rho_cp: 1.0,
        source: &zero_src,
        domain: Domain::FlowActive,
        interior: InteriorFaces::refuse(),
        bcs: Bcs {
            r_inner: FaceBc::HeatFlux(0.0), // ignored: zero-area axis face
            r_outer: FaceBc::Dirichlet(&cold),
            z_lo: FaceBc::HeatFlux(0.0),
            z_hi: FaceBc::HeatFlux(0.0),
        },
    };
    let t_tilde = 0.1f64;
    let t_final = t_tilde * radius * radius; // κ̃ = 1
    let dt = op.stable_dt(&g, CFL_SAFETY);
    let n_steps = (t_final / dt).ceil() as usize;
    let dt = t_final / n_steps as f64;
    op.advance(&mut g, t_id, rate_id, 0.0, dt, n_steps)
        .expect("advance");

    // First five positive zeros of J₀ [META-3: `bessel-j0-zeros` — DLMF
    // §10.21 / Abramowitz & Stegun Table 9.5]. Five terms suffice: at
    // t̃ = 0.1 the n = 5 mode carries e^{−λ₅²·0.1} ≈ 2e-10 of T₀.
    const J0_ZEROS: [f64; 5] = [
        2.404_825_557_695_773,
        5.520_078_110_286_311,
        8.653_727_912_911_013,
        11.791_534_439_014_281,
        14.930_917_708_487_786,
    ];
    let series = |r: f64| -> f64 {
        let mut sum = 0.0f64;
        for &lam in &J0_ZEROS {
            sum += 2.0 / (lam * bessel_j1(lam))
                * bessel_j0(lam * r / radius)
                * (-lam * lam * t_tilde).exp();
        }
        t0 * sum
    };

    let mut worst = 0.0f64;
    for i_r in 0..n_r {
        let r = g.r_center(i_r);
        let bi = g.brick_index(i_r, 0).unwrap();
        let b = g.brick(bi);
        let t = b.field(t_id)[b.cell_index(0, Grid::local_rz(i_r, 0))];
        worst = worst.max((t - series(r)).abs());
    }
    worst
}

/// Power-series J₀ — adequate to x ≈ 15 at f64 (largest argument here is
/// λ₅ ≈ 14.93; worst-case cancellation leaves ~1e-11 absolute error, far
/// below the anchor tolerance).
fn bessel_j0(x: f64) -> f64 {
    let q = 0.25 * x * x;
    let mut term = 1.0f64;
    let mut sum = 1.0f64;
    for m in 1..64 {
        term *= -q / ((m * m) as f64);
        sum += term;
        if term.abs() < 1e-18 {
            break;
        }
    }
    sum
}

fn bessel_j1(x: f64) -> f64 {
    let q = 0.25 * x * x;
    let mut term = 0.5 * x;
    let mut sum = term;
    for m in 1..64 {
        term *= -q / ((m * (m + 1)) as f64);
        sum += term;
        if term.abs() < 1e-18 {
            break;
        }
    }
    sum
}

/// Closed insulated annulus with a smooth bump: relative drift of the total
/// thermal energy after `n_steps` (FND-2 §6-2's closed-sweep half).
pub fn conservation_drift(n_steps: usize) -> f64 {
    let spec = annulus_spec(16, 8, false);
    let mut g = Grid::build(spec, FIELDS).expect("valid spec");
    let t_id = g.field_id("T").unwrap();
    let rate_id = g.field_id("rate").unwrap();
    g.fill_field(t_id, |r, th, z| {
        300.0
            + 50.0
                * ((r - 1.0) * std::f64::consts::PI).cos()
                * (z * 3.0).sin()
                * (1.0 + 0.3 * th.sin())
    });
    let zero_src = |_: f64, _: f64, _: f64, _: f64| 0.0;
    let op = Conduction {
        kappa: 5.0,
        rho_cp: 1.0e3,
        source: &zero_src,
        domain: Domain::FlowActive,
        interior: InteriorFaces::refuse(),
        bcs: Bcs {
            r_inner: FaceBc::HeatFlux(0.0),
            r_outer: FaceBc::HeatFlux(0.0),
            z_lo: FaceBc::HeatFlux(0.0),
            z_hi: FaceBc::HeatFlux(0.0),
        },
    };
    let before = g.reduce_volume_weighted(t_id);
    let dt = op.stable_dt(&g, CFL_SAFETY);
    op.advance(&mut g, t_id, rate_id, 0.0, dt, n_steps)
        .expect("advance");
    let after = g.reduce_volume_weighted(t_id);
    ((after - before) / before).abs()
}
