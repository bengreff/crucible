//! Goal-B **Station 2: the De Laval nozzle** (SOLV-1 §6-4; VAL-2/SOLV-7
//! §6-1; META-3 `maccormack-nozzle` lineage). A stagnation reservoir feeds
//! an axisymmetric converging–diverging contour; the flow accelerates,
//! **chokes at the throat on its own**, and exits supersonic — nothing
//! about the operating point is imposed. The oracle is quasi-1-D isentropic
//! theory: the area–Mach relation, the choked mass flow, and the vacuum
//! thrust coefficient.
//!
//! Geometry enters through the grid's config-time activity seam
//! (`Grid::build_with_activity`) in its **binary (stair-step) degenerate
//! form**: a cell is fluid iff its center lies inside the revolved contour,
//! and stair wall faces use **slip-ghost mirrors about the true contour
//! normal** (`Euler::wall_normal`), cutting spurious wave generation from
//! O(wall slope) to O(h·curvature) at the cost of a small stair-face
//! transpiration flux (reported as the plane-ṁ spread, shrinking with h).
//! Both are interim by design — superseded by FND-3's partial fractions +
//! cut cells (with State Redistribution) when that wave lands. The
//! quasi-1-D oracle itself carries its own model gap vs the true 2-D
//! axisymmetric flow (centerline ≠ area mean); every gap lands in the
//! measured bands below — stated, not hidden (META-1 P6/S8 spirit).
//!
//! Plane diagnostics implement the **SOLV-7 §3.1/§3.2 subset** this station
//! needs (mass/momentum+pressure flux integrals over declared z-planes;
//! emergent p_c = area-averaged *stagnation* pressure at the inlet-end
//! reference plane, the N11 convention) — promoted to the full SOLV-7
//! performance object in its own wave.

use crate::sdc::{FlowClass, Sdc, SdcError};

use crate::euler::{
    Cons, Euler, EulerFields, FlowBc, FlowBcs, GammaLaw, I_RHO, NCOMP, Prim, fill_from_prim, prim6,
};
use crate::station1_sod::{FIELDS, march_to};
use crucible_grid::{Grid, GridSpec};

// --- Named study constants (META-2 §4) --------------------------------------

pub const NOZZLE_GAMMA: f64 = 1.4;
/// Reservoir (stagnation) state — nondimensional cold gas: c0 = √γ ≈ 1.183.
pub const RESERVOIR_P0: f64 = 1.0;
pub const RESERVOIR_RHO0: f64 = 1.0;
pub const RESERVOIR_C: f64 = 0.5;

/// The contour: a **parabolic radius** `R(z) = R*·(1 + RADIUS_COEFF·((z −
/// THROAT_Z)/THROAT_Z)²)` over z ∈ [0, NOZZLE_LENGTH], throat mid-length.
/// Maximum wall slope 2·RADIUS_COEFF·R*/THROAT_Z ≈ 0.167 (a 9.5° half
/// angle) — the regime where the quasi-1-D oracle is honest for the true
/// 2-D axisymmetric flow. (The first fixture used Anderson's didactic
/// `A/A* = 1 + 2.2(z−1.5)²` profile revolved literally — a 44°-half-angle
/// wall whose 2-D flow legitimately departs ~2× from quasi-1-D theory with
/// steep-wall shocks and unsteadiness; that is a statement about the
/// oracle's validity envelope, not about the solver, so the station uses a
/// contour inside the oracle's envelope. META-3 `maccormack-nozzle` remains
/// the quasi-1-D lineage citation.) Inlet/exit area ratio 1.5625 ⇒ inlet
/// M ≈ 0.43, exit M ≈ 1.93 on the supersonic branch.
pub const THROAT_RADIUS: f64 = 1.0;
pub const NOZZLE_LENGTH: f64 = 6.0;
pub const THROAT_Z: f64 = 3.0;
pub const RADIUS_COEFF: f64 = 0.25;

/// Resolution ladder (n_r, n_z): radial cells span [0, R(0)] so the throat
/// radius holds 12.8 / 19.2 / 25.6 cells across the ladder.
pub const NOZZLE_LEVELS: [(usize, usize); 3] = [(16, 96), (24, 144), (32, 192)];

/// Physical march time to steady state from the quasi-1-D initial guess
/// (several flow-through times of the length-6 nozzle), then a further
/// `STEADY_CHECK_TIME` march whose max relative density change is the
/// reported steadiness residual.
pub const SETTLE_TIME: f64 = 25.0;
pub const STEADY_CHECK_TIME: f64 = 1.0;

#[inline]
pub fn wall_radius(z: f64) -> f64 {
    let s = (z - THROAT_Z) / THROAT_Z;
    THROAT_RADIUS * (1.0 + RADIUS_COEFF * s * s)
}

#[inline]
pub fn area_ratio(z: f64) -> f64 {
    let r = wall_radius(z) / THROAT_RADIUS;
    r * r
}

/// Isentropic area–Mach inversion on the chosen branch (deterministic
/// fixed-count bisection; the function is monotone on each branch).
pub fn mach_from_area_ratio(a_ratio: f64, supersonic: bool) -> f64 {
    let ga = NOZZLE_GAMMA;
    let f = |m: f64| -> f64 {
        let t = (2.0 / (ga + 1.0)) * (1.0 + 0.5 * (ga - 1.0) * m * m);
        t.powf(0.5 * (ga + 1.0) / (ga - 1.0)) / m - a_ratio
    };
    let (mut lo, mut hi) = if supersonic { (1.0, 50.0) } else { (1e-8, 1.0) };
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        // f < 0 between the roots' inner side: on the subsonic branch f
        // decreases with M; on the supersonic branch f increases.
        let positive = f(mid) > 0.0;
        if positive == supersonic {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    0.5 * (lo + hi)
}

/// Static primitive state on the reservoir isentrope at Mach `m`, flowing
/// axially (+z).
pub fn isentropic_prim(m: f64) -> Prim {
    let ga = NOZZLE_GAMMA;
    let fac = 1.0 + 0.5 * (ga - 1.0) * m * m;
    let p = RESERVOIR_P0 * fac.powf(-ga / (ga - 1.0));
    let rho = RESERVOIR_RHO0 * fac.powf(-1.0 / (ga - 1.0));
    let c = (ga * p / rho).sqrt();
    prim6(rho, 0.0, 0.0, m * c, p, RESERVOIR_C)
}

/// Ideal choked mass flow through a throat of area `a_star` (sonic-state
/// identity on the reservoir isentrope).
pub fn ideal_choked_mdot(a_star: f64) -> f64 {
    let ga = NOZZLE_GAMMA;
    let rho_s = RESERVOIR_RHO0 * (2.0 / (ga + 1.0)).powf(1.0 / (ga - 1.0));
    let c_s = (ga * RESERVOIR_P0 / RESERVOIR_RHO0).sqrt() * (2.0 / (ga + 1.0)).sqrt();
    rho_s * c_s * a_star
}

/// Closed-form vacuum thrust coefficient for exit Mach `m_e` and exit area
/// ratio ε (SOLV-7 §6-1's C_F(γ, p_e/p_c, ε), p_a = 0).
pub fn ideal_vacuum_cf(m_e: f64, eps: f64) -> f64 {
    let ga = NOZZLE_GAMMA;
    let pe_pc = (1.0 + 0.5 * (ga - 1.0) * m_e * m_e).powf(-ga / (ga - 1.0));
    let a = 2.0 * ga * ga / (ga - 1.0) * (2.0 / (ga + 1.0)).powf((ga + 1.0) / (ga - 1.0));
    (a * (1.0 - pe_pc.powf((ga - 1.0) / ga))).sqrt() + pe_pc * eps
}

const ZERO_SRC: fn(f64, f64, f64, f64) -> Cons = |_, _, _, _| [0.0; NCOMP];

/// Outward unit normal of the true (smooth) contour at axial position z:
/// the wall is r = R(z), so ∇(r − R) = (1, −R′) normalized.
fn contour_normal(_r: f64, z: f64) -> (f64, f64) {
    let dr_dz = 2.0 * RADIUS_COEFF * THROAT_RADIUS * (z - THROAT_Z) / (THROAT_Z * THROAT_Z);
    let inv = 1.0 / (1.0 + dr_dz * dr_dz).sqrt();
    (inv, -dr_dz * inv)
}

fn nozzle_op(eos: GammaLaw) -> Euler<'static> {
    Euler {
        eos,
        source: &ZERO_SRC,
        bcs: FlowBcs {
            r_inner: FlowBc::Reflecting, // ignored: zero-area axis face
            r_outer: FlowBc::Reflecting, // domain edge = the inlet lip wall
            z_lo: FlowBc::StagnationInflow {
                p0: RESERVOIR_P0,
                rho0: RESERVOIR_RHO0,
                c_frac: RESERVOIR_C,
            },
            z_hi: FlowBc::Transmissive, // exit is supersonic once choked
        },
        wall_normal: Some(&contour_normal),
        slip_wall_z_faces: true, // certified station behavior (slip everywhere)
        combustion: None,
    }
}

pub fn build_nozzle(n_r: usize, n_z: usize) -> (Grid, EulerFields) {
    let r_max = wall_radius(0.0);
    let (dr, dz) = (r_max / n_r as f64, NOZZLE_LENGTH / n_z as f64);
    let spec = GridSpec {
        r_min: 0.0,
        dr,
        n_r,
        z_min: 0.0,
        dz,
        n_z,
        n_theta_max: 1,
        axisymmetry_assertion: true,
    };
    let g = Grid::build_with_activity(spec, FIELDS, |i_r, i_z| {
        let r = (i_r as f64 + 0.5) * dr;
        let z = (i_z as f64 + 0.5) * dz;
        r < wall_radius(z)
    })
    .expect("valid spec");
    let f = EulerFields::resolve(&g).expect("fields registered");
    (g, f)
}

/// Initialize with the quasi-1-D isentropic solution (subsonic branch
/// before the throat, supersonic after) and march to `SETTLE_TIME`, then a
/// further `STEADY_CHECK_TIME` to measure the steadiness residual.
/// Returns (grid, fields, total steps, steadiness residual).
pub fn run_nozzle(n_r: usize, n_z: usize) -> Result<(Grid, EulerFields, usize, f64), SdcError> {
    let (mut g, f) = build_nozzle(n_r, n_z);
    let eos = GammaLaw {
        gamma: NOZZLE_GAMMA,
    };
    fill_from_prim(&mut g, &f, &eos, |_, _, z| {
        isentropic_prim(mach_from_area_ratio(area_ratio(z), z > THROAT_Z))
    });
    let op = nozzle_op(eos);
    let steps = march_to(&op, &mut g, &f, 0.0, SETTLE_TIME)?;

    let rho_id = f.ids()[I_RHO];
    let before: Vec<Vec<f64>> = g
        .bricks()
        .iter()
        .map(|b| b.field(rho_id).to_vec())
        .collect();
    let steps2 = march_to(
        &op,
        &mut g,
        &f,
        SETTLE_TIME,
        SETTLE_TIME + STEADY_CHECK_TIME,
    )?;
    let mut resid = 0.0f64;
    for (bi, b) in g.bricks().iter().enumerate() {
        for (a, r) in before[bi].iter().zip(b.field(rho_id)) {
            if *a > 0.0 {
                resid = resid.max(((r - a) / a).abs());
            }
        }
    }
    Ok((g, f, steps + steps2, resid))
}

// --- SOLV-7 §3.1/§3.2 subset: plane integrals over declared z-planes --------

/// Mass flow through z-plane `i_z`: `ṁ = Σ ρu_z·A_z` over active cells
/// (canonical order, fixed-shape reduction not needed at these sizes —
/// serial accumulation in Morton-consistent i_r order).
pub fn plane_mdot(g: &Grid, f: &EulerFields, eos: &GammaLaw, i_z: usize) -> f64 {
    plane_sum(g, f, eos, i_z, |w, _| w[I_RHO] * w[3])
}

/// Momentum-plus-pressure flux through z-plane `i_z` with p_a = 0
/// (SOLV-7.1): `F = Σ (ρu_z² + p)·A_z`.
pub fn plane_thrust(g: &Grid, f: &EulerFields, eos: &GammaLaw, i_z: usize) -> f64 {
    plane_sum(g, f, eos, i_z, |w, _| w[I_RHO] * w[3] * w[3] + w[4])
}

/// Area-averaged **stagnation** pressure at z-plane `i_z` — the SOLV-7
/// §3.2 (N11) emergent-p_c convention.
pub fn plane_stagnation_p(g: &Grid, f: &EulerFields, eos: &GammaLaw, i_z: usize) -> f64 {
    let ga = eos.gamma;
    let num = plane_sum(g, f, eos, i_z, |w, _| {
        let csq = ga * w[4] / w[I_RHO];
        let msq = (w[1] * w[1] + w[2] * w[2] + w[3] * w[3]) / csq;
        w[4] * (1.0 + 0.5 * (ga - 1.0) * msq).powf(ga / (ga - 1.0))
    });
    num / plane_area(g, i_z)
}

/// Open flow area of z-plane `i_z`.
pub fn plane_area(g: &Grid, i_z: usize) -> f64 {
    let mut a = 0.0f64;
    for i_r in 0..g.spec().n_r {
        if g.is_active(i_r, i_z) {
            a += g.face_area_z(i_r, 1);
        }
    }
    a
}

fn plane_sum(
    g: &Grid,
    f: &EulerFields,
    eos: &GammaLaw,
    i_z: usize,
    integrand: impl Fn(&Prim, f64) -> f64,
) -> f64 {
    // The caller's EOS decodes the conserved state (review finding: a
    // hardcoded fixture γ here would silently mix two gammas the moment a
    // non-1.4 gas reuses these SOLV-7 diagnostics).
    let ids = f.ids();
    let mut acc = 0.0f64;
    for i_r in 0..g.spec().n_r {
        if !g.is_active(i_r, i_z) {
            continue;
        }
        let bi = g.brick_index(i_r, i_z).expect("active cell's brick");
        let b = g.brick(bi);
        let idx = b.cell_index(0, Grid::local_rz(i_r, i_z));
        let u: Cons = std::array::from_fn(|k| b.field(ids[k])[idx]);
        let w = eos.prim_checked(&u).expect("physical state at plane");
        acc += integrand(&w, g.r_center(i_r)) * g.face_area_z(i_r, 1);
    }
    acc
}

/// The narrowest open plane of the discrete (stair-stepped) nozzle: its
/// area and z-index — the throat the flow actually sees.
pub fn discrete_throat(g: &Grid) -> (f64, usize) {
    let mut best = (f64::INFINITY, 0usize);
    for i_z in 0..g.spec().n_z {
        let a = plane_area(g, i_z);
        if a < best.0 {
            best = (a, i_z);
        }
    }
    best
}

/// Centerline Mach number per z-column (innermost ring).
pub fn centerline_mach(g: &Grid, f: &EulerFields) -> Vec<f64> {
    let eos = GammaLaw {
        gamma: NOZZLE_GAMMA,
    };
    let ids = f.ids();
    (0..g.spec().n_z)
        .map(|i_z| {
            let bi = g.brick_index(0, i_z).expect("axis column active");
            let b = g.brick(bi);
            let idx = b.cell_index(0, Grid::local_rz(0, i_z));
            let u: Cons = std::array::from_fn(|k| b.field(ids[k])[idx]);
            let w = eos.prim_checked(&u).expect("physical state on axis");
            let c = eos.sound_speed(w[I_RHO], w[4]);
            (w[1] * w[1] + w[2] * w[2] + w[3] * w[3]).sqrt() / c
        })
        .collect()
}

// --- The study ---------------------------------------------------------------

#[derive(Debug, Clone, PartialEq)]
pub struct NozzleLevel {
    pub n_r: usize,
    pub n_z: usize,
    pub steps: usize,
    /// Steadiness: max relative ρ change over the final unit-time march.
    pub steady_resid: f64,
    /// Discharge coefficient vs the discrete (stair-step) throat area and
    /// vs the analytic πR*² — geometric + flow error separated.
    pub cd_discrete: f64,
    pub cd_analytic: f64,
    /// ṁ uniformity along the nozzle: max |ṁ(z) − ṁ_exit| / ṁ_exit.
    pub mdot_spread: f64,
    /// Emergent chamber pressure (inlet-plane stagnation avg) vs p0.
    pub p_c_over_p0: f64,
    /// Area-averaged... exit-plane diagnostics:
    pub exit_mach_centerline: f64,
    /// max |M_centerline(z) − M_1D(z)|/M_1D(z) excluding only the inlet
    /// band (z < INLET_EXCLUSION). The throat-adjacent band IS measured —
    /// the committed certificate's 6% gate absorbs the sonic-transition
    /// deviation (review fix: an earlier doc claimed a THROAT_EXCLUSION
    /// that never existed in code; the measurement, not the doc, is the
    /// certified behavior).
    pub mach_dev_max: f64,
    /// Thrust coefficient F/(p_c·A*_discrete) vs the ideal vacuum C_F at
    /// the discrete exit conditions.
    pub cf: f64,
    pub cf_ideal: f64,
}

/// The entrance band excluded from the pointwise area–Mach comparison: the
/// first-order stagnation-inflow ghosts and the 2-D entrance turning (the
/// wall's slope is largest at the inlet) are not in quasi-1-D theory; the
/// measured deviation decays from ~+14% at the first cell to ~+1% by
/// z ≈ 1 (recorded in the artifact profile).
pub const INLET_EXCLUSION: f64 = 1.0;

// --- Pass criteria (measured 2026-08-17; gates set with headroom; tests and
// --- the committed artifact both read exactly these constants) --------------

/// Discharge coefficient ṁ_exit/ṁ_ideal(πR*²) — THE choked-flow criterion.
/// Measured 1.0060 / 1.0022 / 1.0020 across the ladder: band ±1.5% at every
/// level, ≤0.8% at the finest (measured 0.20%), and |Cd−1| non-increasing.
pub const CD_BAND: (f64, f64) = (0.985, 1.015);
pub const CD_FINEST_TOL: f64 = 0.008;

/// Emergent chamber pressure (inlet-plane stagnation average) vs the
/// reservoir. Measured 0.99999 — the flow reports the reservoir back.
pub const PC_BAND: (f64, f64) = (0.998, 1.001);

/// Steadiness residual (max relative ρ change over the final unit-time
/// march). Measured 1.7e-2 / 5.3e-3 / 7.6e-4 — converging with h.
pub const STEADY_RESID_MAX: f64 = 3e-2;
pub const STEADY_RESID_FINEST: f64 = 3e-3;

/// Centerline area–Mach deviation past the entrance band. Measured
/// 3.6–4.2%: the compound of stair-wall resolution and the quasi-1-D
/// oracle's own 2-D model floor (centerline ≠ area mean in a 9.5° nozzle)
/// — declared, not hidden.
pub const MACH_DEV_MAX: f64 = 0.06;

/// Exit centerline Mach vs the 1-D exit value (measured −4 to −6%: the 2-D
/// centerline lag), and thrust coefficient vs ideal vacuum C_F (measured
/// ≤4.3%).
pub const EXIT_MACH_DEV_MAX: f64 = 0.08;
pub const CF_DEV_MAX: f64 = 0.06;

/// Plane-ṁ uniformity along the nozzle: the interim slip-ghost wall
/// admits a small stair-face transpiration flux (the declared trade that
/// buys O(slope)→O(h·curvature) spurious-wave reduction; FND-3 cut cells
/// retire it). Measured 7.2% / 5.2% / 3.8% — must stay bounded and shrink.
pub const MDOT_SPREAD_MAX: f64 = 0.10;

/// A uniform gas at rest inside the masked (stair-stepped) nozzle cavity,
/// all boundaries walls, must be a **bitwise** fixed point — with the
/// grid-aligned mirror AND with the slip-ghost wall (reflecting a zero
/// velocity about any normal is the identity). This is the exactness
/// certificate for the run-decomposed masked sweeps themselves.
pub fn masked_uniform_fixed_point(slip_wall: bool) -> bool {
    let (mut g, f) = build_nozzle(12, 48);
    let eos = GammaLaw {
        gamma: NOZZLE_GAMMA,
    };
    fill_from_prim(&mut g, &f, &eos, |_, _, _| {
        prim6(1.3, 0.0, 0.0, 0.0, 2.7, 0.5)
    });
    let op = Euler {
        eos,
        source: &ZERO_SRC,
        bcs: FlowBcs {
            r_inner: FlowBc::Reflecting,
            r_outer: FlowBc::Reflecting,
            z_lo: FlowBc::Reflecting,
            z_hi: FlowBc::Reflecting,
        },
        wall_normal: if slip_wall {
            Some(&contour_normal)
        } else {
            None
        },
        slip_wall_z_faces: true, // certified station behavior (slip everywhere)
        combustion: None,
    };
    let bits = |g: &Grid| -> Vec<u64> {
        g.bricks()
            .iter()
            .flat_map(|b| {
                f.ids()
                    .into_iter()
                    .flat_map(|id| b.field(id).iter().map(|v| v.to_bits()).collect::<Vec<_>>())
            })
            .collect()
    };
    let before = bits(&g);
    let mut sdc = Sdc::new();
    let flow = FlowClass {
        op: &op,
        fields: &f,
    };
    let dt = sdc.stable_dt(&g, &flow, 0.4).expect("dt");
    let mut t = 0.0;
    for _ in 0..50 {
        sdc.step_flow(&mut g, &flow, t, dt).expect("step");
        t += dt;
    }
    bits(&g) == before
}

pub fn nozzle_study() -> Result<Vec<NozzleLevel>, SdcError> {
    let eos = GammaLaw {
        gamma: NOZZLE_GAMMA,
    };
    let mut out = Vec::new();
    for &(n_r, n_z) in &NOZZLE_LEVELS {
        let (g, f, steps, steady_resid) = run_nozzle(n_r, n_z)?;
        let n_zc = g.spec().n_z;
        let exit = n_zc - 1;
        let (a_star_disc, _) = discrete_throat(&g);
        let a_star_analytic = std::f64::consts::PI * THROAT_RADIUS * THROAT_RADIUS;

        let mdot_exit = plane_mdot(&g, &f, &eos, exit);
        let mut mdot_spread = 0.0f64;
        for i_z in 0..n_zc {
            let m = plane_mdot(&g, &f, &eos, i_z);
            mdot_spread = mdot_spread.max(((m - mdot_exit) / mdot_exit).abs());
        }

        let p_c = plane_stagnation_p(&g, &f, &eos, 0);
        let mach = centerline_mach(&g, &f);
        let dz = g.spec().dz;
        let mut mach_dev_max = 0.0f64;
        for (i_z, m_num) in mach.iter().enumerate() {
            let z = (i_z as f64 + 0.5) * dz;
            if z < INLET_EXCLUSION {
                continue;
            }
            let m_1d = mach_from_area_ratio(area_ratio(z), z > THROAT_Z);
            mach_dev_max = mach_dev_max.max(((m_num - m_1d) / m_1d).abs());
        }

        let thrust = plane_thrust(&g, &f, &eos, exit);
        let eps_disc = plane_area(&g, exit) / a_star_disc;
        let m_e_1d = mach_from_area_ratio(eps_disc, true);

        out.push(NozzleLevel {
            n_r,
            n_z,
            steps,
            steady_resid,
            cd_discrete: mdot_exit / ideal_choked_mdot(a_star_disc),
            cd_analytic: mdot_exit / ideal_choked_mdot(a_star_analytic),
            mdot_spread,
            p_c_over_p0: p_c / RESERVOIR_P0,
            exit_mach_centerline: mach[exit],
            mach_dev_max,
            cf: thrust / (p_c * a_star_disc),
            cf_ideal: ideal_vacuum_cf(m_e_1d, eps_disc),
        });
    }
    Ok(out)
}
