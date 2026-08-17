//! SOLV-1 §6-2 — **manufactured-solution verification of the whole Euler
//! operator** (META-3 `mms`; VAL-3 §3.3 "whole PDE operator" row): a smooth
//! non-natural field is substituted into the full cylindrical Euler system
//! — every flux direction AND all three geometric source terms (pressure,
//! centrifugal ρu_θ², swirl advection ρu_ru_θ) active — and the analytic
//! residual is fed back through the operator's external-source intake. The
//! computed solution must recover the manufactured one at formal order.
//! This certifies the terms the Sod tube structurally cannot: Sod's radial
//! velocity is identically zero, so its run never exercises radial
//! advection or the swirl sources with nonzero operands.
//!
//! The manufactured field is one shared spatial mode with per-component
//! amplitudes: `W_i = base_i + amp_i·φ`,
//! `φ = sin(a(r−r₀))·cos(b·z)·(1 + ε·cos(mθ))·e^{−λt}` — monotone per
//! pencil in r and z over the domain (quarter waves), so the limiters stay
//! inactive there and the observed order is the scheme's; the θ mode's
//! extrema cost only O(h³) in L2 (isolated-extremum clipping). The source
//! is assembled **exactly** by product rule from the mode's analytic
//! derivatives — no symbolic algebra by hand per term, no numerical
//! differentiation anywhere.

use crate::euler::{
    Cons, Euler, EulerFields, FlowBc, FlowBcs, GammaLaw, NCOMP, Prim, fill_from_prim,
};
use crate::station1_sod::{FIELDS, march_to};
use crucible_grid::{Grid, GridSpec};

// --- Named study constants (META-2 §4) --------------------------------------

pub const MMS_GAMMA: f64 = 1.4;
/// Base state and mode amplitude per primitive [ρ, u_r, u_θ, u_z, p, C]:
/// everything positive and subsonic over the whole run; swirl and radial
/// flow both active.
pub const MMS_BASE: Prim = [1.0, 0.0, 0.0, 0.1, 1.0, 0.5];
pub const MMS_AMP: Prim = [0.2, 0.15, 0.2, 0.15, 0.3, 0.3];
/// Quarter-wave numbers over the unit r/z extents (monotone per pencil).
pub const MMS_A: f64 = std::f64::consts::FRAC_PI_2;
pub const MMS_B: f64 = std::f64::consts::FRAC_PI_2;
pub const MMS_LAMBDA: f64 = 0.5;
/// Azimuthal mode for the 3-D study: m = 2 at relative amplitude ε.
pub const MMS_M: u32 = 2;
pub const MMS_EPS: f64 = 0.3;
pub const MMS_T_FINAL: f64 = 0.1;
/// Annulus off the axis (r_min > 0): the N_θ > 1 study needs no cross-axis
/// parity gather, and the metric sources are O(1) everywhere.
pub const MMS_R_MIN: f64 = 0.5;
/// Refinement ladder; every n is θ-ladder-aligned (4·2^k) for the 3-D study.
/// Starts at 16: the 8² annulus is measurably pre-asymptotic (orders
/// 1.45–1.76 vs 1.82–2.05 on the finer pairs).
pub const MMS_LEVELS: [usize; 3] = [16, 32, 64];

/// Observed-order acceptance band per conserved component (VAL-3 §3.2:
/// observed within ~10% of formal order 2 on the low side; measured
/// 1.88–2.08 (axisym-swirl) and 1.95–2.22 (θ-mode) across all six
/// components — mild pre-asymptotic superconvergence caps the band at 2.3).
pub const MMS_ORDER_MIN: f64 = 1.8;
pub const MMS_ORDER_MAX: f64 = 2.3;

/// The manufactured field and its exact partial derivatives at one point:
/// value plus ∂/∂r, ∂/∂θ, ∂/∂z, ∂/∂t of every primitive.
struct Manufactured {
    w: Prim,
    dr: Prim,
    dth: Prim,
    dz: Prim,
    dt: Prim,
}

fn manufactured(r: f64, theta: f64, z: f64, t: f64, eps: f64) -> Manufactured {
    let m = f64::from(MMS_M);
    let (sr, cr) = (MMS_A * (r - MMS_R_MIN)).sin_cos();
    let (sz, cz) = (MMS_B * z).sin_cos();
    let (st, ct) = (m * theta).sin_cos();
    let decay = (-MMS_LAMBDA * t).exp();
    let ang = 1.0 + eps * ct;

    let phi = sr * cz * ang * decay;
    let phi_r = MMS_A * cr * cz * ang * decay;
    let phi_th = sr * cz * (-eps * m * st) * decay;
    let phi_z = -MMS_B * sr * sz * ang * decay;
    let phi_t = -MMS_LAMBDA * phi;

    let build = |scale: f64| -> Prim { std::array::from_fn(|k| MMS_AMP[k] * scale) };
    let mut w = build(phi);
    for k in 0..NCOMP {
        w[k] += MMS_BASE[k];
    }
    Manufactured {
        w,
        dr: build(phi_r),
        dth: build(phi_th),
        dz: build(phi_z),
        dt: build(phi_t),
    }
}

/// The exact MMS residual `S = ∂U/∂t + ∇·F(U) − S_geom(U)` of the full
/// cylindrical Euler system, assembled by product rule from the analytic
/// primitive derivatives. Feeding this to the operator's source intake
/// makes the manufactured field an exact solution of the forced system.
#[allow(clippy::similar_names)]
fn mms_source(r: f64, theta: f64, z: f64, t: f64, eps: f64, eos: &GammaLaw) -> Cons {
    let mf = manufactured(r, theta, z, t, eps);
    let [rho, ur, ut, uz, p, c] = mf.w;
    let gm1 = eos.gamma - 1.0;

    // Per-direction primitive derivative bundles.
    let d = [mf.dr, mf.dth, mf.dz, mf.dt];
    let (rho_x, ur_x, ut_x, uz_x, p_x, c_x) = (
        [d[0][0], d[1][0], d[2][0], d[3][0]],
        [d[0][1], d[1][1], d[2][1], d[3][1]],
        [d[0][2], d[1][2], d[2][2], d[3][2]],
        [d[0][3], d[1][3], d[2][3], d[3][3]],
        [d[0][4], d[1][4], d[2][4], d[3][4]],
        [d[0][5], d[1][5], d[2][5], d[3][5]],
    );

    // Total energy and its derivatives.
    let q2 = ur * ur + ut * ut + uz * uz;
    let e_tot = p / gm1 + 0.5 * rho * q2;
    let e_x = |x: usize| -> f64 {
        p_x[x] / gm1 + 0.5 * rho_x[x] * q2 + rho * (ur * ur_x[x] + ut * ut_x[x] + uz * uz_x[x])
    };

    // Directional indices into the derivative bundles.
    const R: usize = 0;
    const TH: usize = 1;
    const Z: usize = 2;
    const T: usize = 3;
    let inv_r = 1.0 / r;

    // ∂x(a·b) helpers over the bundles.
    let dx = |ax: &[f64; 4], a: f64, bx: &[f64; 4], b: f64, x: usize| ax[x] * b + a * bx[x];

    // Mass: ρ_t + (ρu_r)/r + ∂r(ρu_r) + (1/r)∂θ(ρu_θ) + ∂z(ρu_z).
    let m_r = rho * ur;
    let s_rho = rho_x[T]
        + m_r * inv_r
        + dx(&rho_x, rho, &ur_x, ur, R)
        + inv_r * dx(&rho_x, rho, &ut_x, ut, TH)
        + dx(&rho_x, rho, &uz_x, uz, Z);

    // r-momentum: ∂t(ρu_r) + (ρu_r²+p)/r + ∂r(ρu_r²+p) + (1/r)∂θ(ρu_θu_r)
    //             + ∂z(ρu_zu_r) − (p + ρu_θ²)/r.
    let dt_mr = rho_x[T] * ur + rho * ur_x[T];
    let dr_flux = rho_x[R] * ur * ur + 2.0 * rho * ur * ur_x[R] + p_x[R];
    let dth_flux = rho_x[TH] * ut * ur + rho * ut_x[TH] * ur + rho * ut * ur_x[TH];
    let dz_flux = rho_x[Z] * uz * ur + rho * uz_x[Z] * ur + rho * uz * ur_x[Z];
    let s_mr = dt_mr + (rho * ur * ur + p) * inv_r + dr_flux + inv_r * dth_flux + dz_flux
        - (p + rho * ut * ut) * inv_r;

    // θ-momentum: ∂t(ρu_θ) + (ρu_ru_θ)/r + ∂r(ρu_ru_θ) + (1/r)∂θ(ρu_θ²+p)
    //             + ∂z(ρu_zu_θ) + ρu_ru_θ/r.
    let dt_mt = rho_x[T] * ut + rho * ut_x[T];
    let dr_flux = rho_x[R] * ur * ut + rho * ur_x[R] * ut + rho * ur * ut_x[R];
    let dth_flux = rho_x[TH] * ut * ut + 2.0 * rho * ut * ut_x[TH] + p_x[TH];
    let dz_flux = rho_x[Z] * uz * ut + rho * uz_x[Z] * ut + rho * uz * ut_x[Z];
    let s_mt = dt_mt
        + rho * ur * ut * inv_r
        + dr_flux
        + inv_r * dth_flux
        + dz_flux
        + rho * ur * ut * inv_r;

    // z-momentum: ∂t(ρu_z) + (ρu_ru_z)/r + ∂r(ρu_ru_z) + (1/r)∂θ(ρu_θu_z)
    //             + ∂z(ρu_z²+p).
    let dt_mz = rho_x[T] * uz + rho * uz_x[T];
    let dr_flux = rho_x[R] * ur * uz + rho * ur_x[R] * uz + rho * ur * uz_x[R];
    let dth_flux = rho_x[TH] * ut * uz + rho * ut_x[TH] * uz + rho * ut * uz_x[TH];
    let dz_flux = rho_x[Z] * uz * uz + 2.0 * rho * uz * uz_x[Z] + p_x[Z];
    let s_mz = dt_mz + rho * ur * uz * inv_r + dr_flux + inv_r * dth_flux + dz_flux;

    // Energy: E_t + u_r(E+p)/r + ∂r(u_r(E+p)) + (1/r)∂θ(u_θ(E+p)) + ∂z(u_z(E+p)).
    let h = e_tot + p;
    let h_x = |x: usize| e_x(x) + p_x[x];
    let s_en = e_x(T)
        + ur * h * inv_r
        + (ur_x[R] * h + ur * h_x(R))
        + inv_r * (ut_x[TH] * h + ut * h_x(TH))
        + (uz_x[Z] * h + uz * h_x(Z));

    // Species: ∂t(ρC) + (ρCu_r)/r + ∂r(ρCu_r) + (1/r)∂θ(ρCu_θ) + ∂z(ρCu_z).
    let g = rho * c;
    let g_x = |x: usize| rho_x[x] * c + rho * c_x[x];
    let s_rc = g_x(T)
        + g * ur * inv_r
        + (g_x(R) * ur + g * ur_x[R])
        + inv_r * (g_x(TH) * ut + g * ut_x[TH])
        + (g_x(Z) * uz + g * uz_x[Z]);

    [s_rho, s_mr, s_mt, s_mz, s_en, s_rc]
}

#[derive(Debug, Clone, PartialEq)]
pub struct MmsEulerLevel {
    pub n: usize,
    pub h: f64,
    /// Volume-weighted L1 error per conserved component.
    pub l1: [f64; NCOMP],
}

#[derive(Debug, Clone, PartialEq)]
pub struct MmsEulerStudy {
    pub label: &'static str,
    pub levels: Vec<MmsEulerLevel>,
}

impl MmsEulerStudy {
    /// Observed order per refinement for component `k`.
    pub fn observed_orders(&self, k: usize) -> Vec<f64> {
        self.levels
            .windows(2)
            .map(|w| (w[0].l1[k] / w[1].l1[k]).ln() / std::f64::consts::LN_2)
            .collect()
    }
}

/// 2-D axisymmetric-with-swirl MMS (N_θ = 1 under the recorded assertion):
/// all radial/axial fluxes and all three geometric sources active.
pub fn mms_axisym_swirl() -> MmsEulerStudy {
    mms_run(
        "Euler MMS 2-D axisymmetric with swirl (N_θ = 1)",
        0.0,
        &MMS_LEVELS,
    )
}

/// 3-D MMS with an m = 2 azimuthal mode at N_θ = n — adds the θ-flux path.
pub fn mms_theta_mode() -> MmsEulerStudy {
    mms_run(
        "Euler MMS 3-D with m = 2 θ-mode (N_θ = n)",
        MMS_EPS,
        &MMS_LEVELS,
    )
}

/// Ladder-parameterized variant (asymptotic-trend probes).
pub fn mms_theta_mode_at(levels: &[usize]) -> MmsEulerStudy {
    mms_run("Euler MMS 3-D with m = 2 θ-mode (N_θ = n)", MMS_EPS, levels)
}

fn mms_run(label: &'static str, eps: f64, ns: &[usize]) -> MmsEulerStudy {
    let eos = GammaLaw { gamma: MMS_GAMMA };
    let mut levels = Vec::new();
    for &n in ns {
        let n_theta = if eps == 0.0 { 1 } else { n as u32 };
        let spec = GridSpec {
            r_min: MMS_R_MIN,
            dr: 1.0 / n as f64,
            n_r: n,
            z_min: 0.0,
            dz: 1.0 / n as f64,
            n_z: n,
            n_theta_max: n_theta,
            axisymmetry_assertion: n_theta == 1,
        };
        let mut g = Grid::build(spec, FIELDS).expect("valid spec");
        let f = EulerFields::resolve(&g).expect("fields");
        fill_from_prim(&mut g, &f, &eos, |r, th, z| {
            manufactured(r, th, z, 0.0, eps).w
        });

        let exact_bc = move |r: f64, th: f64, z: f64, t: f64| manufactured(r, th, z, t, eps).w;
        let source = move |r: f64, th: f64, z: f64, t: f64| mms_source(r, th, z, t, eps, &eos);
        let op = Euler {
            eos,
            source: &source,
            bcs: FlowBcs {
                r_inner: FlowBc::Prescribed(&exact_bc),
                r_outer: FlowBc::Prescribed(&exact_bc),
                z_lo: FlowBc::Prescribed(&exact_bc),
                z_hi: FlowBc::Prescribed(&exact_bc),
            },
            wall_normal: None,
        };
        march_to(&op, &mut g, &f, 0.0, MMS_T_FINAL).expect("march");

        // Volume-weighted L1 per conserved component (canonical traversal).
        // L1 is the shock-capturing verification norm (SOLV-1 §6-1 "L1 at
        // formal order"): the limiter's clipping at the θ-mode's smooth
        // extrema is locally 1st-order over an O(h) measure, which L2
        // amplifies to a measured ~O(h^1.6) tail while L1 keeps the formal
        // order (measured; the artifact records this caveat).
        let ids = f.ids();
        let mut num = [0.0f64; NCOMP];
        let mut den = 0.0f64;
        g.for_each_active_cell(|cell| {
            let b = g.brick(cell.bi);
            let v = g.cell_volume(cell.i_r, b.n_theta());
            let ue =
                eos.prim_to_cons(&manufactured(cell.r, cell.theta, cell.z, MMS_T_FINAL, eps).w);
            for k in 0..NCOMP {
                num[k] += v * (b.field(ids[k])[cell.idx] - ue[k]).abs();
            }
            den += v;
        });
        levels.push(MmsEulerLevel {
            n,
            h: 1.0 / n as f64,
            l1: std::array::from_fn(|k| num[k] / den),
        });
    }
    MmsEulerStudy { label, levels }
}
