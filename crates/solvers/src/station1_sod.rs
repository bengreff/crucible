//! Goal-B **Station 1: the bursting diaphragm** (Sod shock tube) — SOLV-1
//! §6-1 / VAL-2 §3.3 `sod-shock`. A tube of gas along z with a diaphragm at
//! z = 0.5: high pressure left, low right. Burst it and exactly three waves
//! must emerge — a right-running shock, a contact surface carrying the
//! composition jump, and a left-running rarefaction fan — all checked
//! against the exact Riemann solution (the oracle in `euler::exact`).
//!
//! Everything runs on the one cylindrical operator at N_θ = 1 (recorded
//! axisymmetry assertion) with the r = 0 axis *inside* the domain: the tube
//! is a full cylinder, radially uniform, so the run simultaneously certifies
//! that the axis metric and the well-balanced geometric sources leave a
//! radially-uniform state radially uniform — bitwise (asserted).
//!
//! Same contract as the Goal-A certificate: runners return structured data;
//! `tests/solv1_station1_sod.rs` asserts the criteria and
//! `bin/station1_sod_certificate` renders the committed artifact from the
//! SAME named constants, gate 4 of `scripts/check.sh` diffs it. Everything
//! is deterministic: fixed grids, CFL-derived step sequences, no RNG.

use crate::sdc::{FlowClass, Sdc, SdcError};

use crate::euler::{
    Cons, Euler, EulerFields, FlowBc, FlowBcs, GammaLaw, I_RC, I_RHO, NCOMP, Prim, RiemannSide,
    RiemannSolution, fill_from_prim, prim6, solve_riemann,
};
use crucible_grid::{Grid, GridSpec};

pub use crate::euler::EULER_FIELDS as FIELDS;

// --- Named certificate constants (META-2 §4: no magic numbers) --------------

/// The canonical Sod data (META-3 `sod-shock`: Sod 1978; Toro Table 4.2).
pub const SOD_GAMMA: f64 = 1.4;
pub const SOD_LEFT: RiemannSide = RiemannSide {
    rho: 1.0,
    u: 0.0,
    p: 1.0,
};
pub const SOD_RIGHT: RiemannSide = RiemannSide {
    rho: 0.125,
    u: 0.0,
    p: 0.1,
};
pub const SOD_T_FINAL: f64 = 0.2;
pub const SOD_DIAPHRAGM_Z: f64 = 0.5;

/// CFL safety for the SDC-IMEX march (the explicit-class stability
/// polynomial 1 + z + z²/2 + z³/4 is imaginary-axis stable to |z| ≤ 2 and
/// upwind-stable to CFL 1 in 1-D; 0.4 covers the summed multi-direction
/// bound with margin — same safety philosophy as Goal-A).
pub const CFL_FLOW: f64 = 0.4;

/// Sod refinement ladder (axial cells; the tube is 4 radial × n_z cells).
pub const SOD_LEVELS: [usize; 4] = [100, 200, 400, 800];

/// Smooth-region windows at t = 0.2 (fixed physical bounds with margin from
/// every wave): the rarefaction-fan interior [0.30, 0.45] (head at z ≈
/// 0.263, tail at z ≈ 0.486) and the star-region plateaus — left of the
/// contact [0.52, 0.66] (contact at z ≈ 0.686), right of it [0.72, 0.83]
/// (shock at z ≈ 0.850).
pub const FAN_WINDOW: (f64, f64) = (0.30, 0.45);
pub const STAR_L_WINDOW: (f64, f64) = (0.52, 0.66);
pub const STAR_R_WINDOW: (f64, f64) = (0.72, 0.83);

// --- Pass criteria (measured 2026-08-16, tolerances set with headroom; the
// --- tests and the committed artifact both read exactly these constants) ---

/// VAL-2 §3.3: "shock within 1 cell". Measured 0.16 cells at n_z = 800.
pub const SHOCK_POS_TOL_CELLS: f64 = 1.0;

/// The contact is HLLC-Batten's protected wave but smears wider than the
/// shock (no self-steepening). Measured 1.11 cells; 2 cells catches any
/// loss of the Batten contact restoration (HLLE-class smearing is ≫ this).
pub const CONTACT_POS_TOL_CELLS: f64 = 2.0;

/// Star-plateau max errors (ρ*L, ρ*R, u*, p* windows). Measured ≤ 8.9e-4
/// at n_z = 800 on the S2 SDC spine (was ≤ 7.3e-5 under SSP-RK2 — the
/// less-dissipative stage structure rings slightly more behind the
/// captured shock; global L1 and smooth-region orders are unchanged, so
/// this is a dissipation-profile shift, not an accuracy loss). 2.5e-3
/// keeps ~2.8× regression headroom on states of order 0.3–1.
pub const STAR_PLATEAU_TOL: f64 = 2.5e-3;

/// Composition boundedness beyond [0, 1] over the whole tube. Measured
/// 4e-16 (pure round-off — HLLC's upwind species flux creates no new
/// extrema); 1e-12 is 4 decades of headroom yet fails on any real
/// overshoot, which enters at limiter scale (≥1e-6).
pub const C_BOUNDS_TOL: f64 = 1e-12;

/// Global L1(ρ) at the finest level (measured 5.0e-4) and the minimum
/// decrease factor per refinement (measured ≥ 1.84×; discontinuity-limited
/// convergence at observed order ~0.9–1.1).
pub const SOD_FINEST_L1_MAX: f64 = 1e-3;
pub const SOD_L1_DECREASE_MIN: f64 = 1.5;

/// Smooth-region orders. Star-left plateau (truly smooth, all waves gone):
/// measured 2.18–2.57 per refinement — the formal-order-in-smooth-regions
/// criterion. Fan interior: measured 0.98–1.06 — the known first-order
/// behavior of centered rarefactions (the t → 0 startup singularity
/// contaminates the fan; not a smooth-region defect), gated so it cannot
/// silently degrade further.
pub const STAR_WINDOW_ORDER_MIN: f64 = 1.8;
pub const FAN_ORDER_MIN: f64 = 0.8;

/// Smooth advection study: mean observed order across the whole ladder
/// (measured 2.09 for ρ) and the per-refinement band for the passively
/// advected composition (measured 1.99–2.19 — textbook 2nd order; ρ
/// couples to the acoustic families and oscillates 1.6–2.8 per step, so ρ
/// is gated on the mean plus a per-step floor of 1.4).
pub const ADV_MEAN_ORDER_MIN: f64 = 1.8;
pub const ADV_RHO_STEP_ORDER_MIN: f64 = 1.4;
pub const ADV_C_ORDER_MIN: f64 = 1.8;
pub const ADV_C_ORDER_MAX: f64 = 2.3;

/// One refinement level: axial resolution, global L1(ρ) error, and the
/// smooth-window L1(ρ) errors (fan interior; star-left plateau).
#[derive(Debug, Clone, PartialEq)]
pub struct SodLevel {
    pub n_z: usize,
    pub steps: usize,
    pub l1_rho: f64,
    pub l1_fan: f64,
    pub l1_star_l: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SodStudy {
    pub levels: Vec<SodLevel>,
}

impl SodStudy {
    /// Observed order between consecutive levels for a chosen error field.
    pub fn observed_orders(&self, err: impl Fn(&SodLevel) -> f64) -> Vec<f64> {
        self.levels
            .windows(2)
            .map(|w| (err(&w[0]) / err(&w[1])).ln() / std::f64::consts::LN_2)
            .collect()
    }
}

/// Finest-level wave diagnostics against the exact Riemann solution.
#[derive(Debug, Clone, PartialEq)]
pub struct SodWaves {
    pub n_z: usize,
    /// |captured shock front − exact position| in cells.
    pub shock_pos_err_cells: f64,
    /// |composition-midpoint crossing − exact contact position| in cells.
    pub contact_pos_err_cells: f64,
    /// max |ρ − ρ*L| over the star-left window.
    pub rho_star_l_err: f64,
    /// max |ρ − ρ*R| over the star-right window.
    pub rho_star_r_err: f64,
    /// max |u_z − u*| over the FULL star span [STAR_L_WINDOW.0,
    /// STAR_R_WINDOW.1] — deliberately including the smeared contact:
    /// u and p are continuous across a contact discontinuity, so flatness
    /// through it is part of what HLLC must deliver (ρ, which jumps there,
    /// is measured per-window). Review clarification: the certified
    /// 7.3e-5 plateau number has always included the contact span.
    pub u_star_err: f64,
    /// max |p − p*| over the full star span (see `u_star_err`).
    pub p_star_err: f64,
    /// Composition bounds over the whole tube (must stay in [0, 1] up to
    /// round-off — HLLC-Batten's contact treatment is what protects this).
    pub c_min: f64,
    pub c_max: f64,
    /// Bitwise radial uniformity held? (every ring identical to ring 0)
    pub radially_uniform_bitwise: bool,
}

/// Radial extent of the tube: 4 rings crossing the r = 0 axis. The radial
/// spacing is much coarser than dz so the axial CFL governs.
const TUBE_N_R: usize = 4;
const TUBE_DR: f64 = 0.025;

fn tube_spec(n_z: usize) -> GridSpec {
    GridSpec {
        r_min: 0.0,
        dr: TUBE_DR,
        n_r: TUBE_N_R,
        z_min: 0.0,
        dz: 1.0 / n_z as f64,
        n_z,
        n_theta_max: 1,
        axisymmetry_assertion: true,
    }
}

pub fn sod_exact() -> RiemannSolution {
    solve_riemann(SOD_LEFT, SOD_RIGHT, SOD_GAMMA)
}

fn sod_ic(z: f64) -> Prim {
    if z < SOD_DIAPHRAGM_Z {
        prim6(SOD_LEFT.rho, 0.0, 0.0, 0.0, SOD_LEFT.p, 0.0)
    } else {
        prim6(SOD_RIGHT.rho, 0.0, 0.0, 0.0, SOD_RIGHT.p, 1.0)
    }
}

const ZERO_SRC: fn(f64, f64, f64, f64) -> Cons = |_, _, _, _| [0.0; NCOMP];

fn closed_tube_op(eos: GammaLaw) -> Euler<'static> {
    Euler {
        eos,
        source: &ZERO_SRC,
        bcs: FlowBcs {
            r_inner: FlowBc::Reflecting, // ignored: zero-area axis face
            r_outer: FlowBc::Reflecting,
            z_lo: FlowBc::Reflecting,
            z_hi: FlowBc::Reflecting,
        },
        wall_normal: None,
        slip_wall_z_faces: true, // certified station behavior (slip everywhere)
    }
}

/// Deterministic CFL march to `t_final` on the production integrator (the
/// SDC-IMEX step, flow-only schedule; COUP-2's audit armed every step): dt
/// re-derived from the current wave speeds each step (a fixed rule over
/// the data — COUP-3 §2), clipped to land on `t_final` exactly. Returns
/// the step count.
pub fn march_to(
    op: &Euler<'_>,
    g: &mut Grid,
    f: &EulerFields,
    t0: f64,
    t_final: f64,
) -> Result<usize, SdcError> {
    let mut sdc = Sdc::new();
    let flow = FlowClass { op, fields: f };
    let mut t = t0;
    let mut steps = 0usize;
    while t_final - t > 1e-12 * t_final {
        let dt = sdc.stable_dt(g, &flow, CFL_FLOW)?.min(t_final - t);
        sdc.step_flow(g, &flow, t, dt)?;
        t += dt;
        steps += 1;
        assert!(steps < 1_000_000, "runaway march — CFL collapse");
    }
    Ok(steps)
}

/// Run the Sod tube at `n_z` to t = 0.2; returns (grid, fields, steps).
pub fn run_sod(n_z: usize) -> (Grid, EulerFields, usize) {
    let mut g = Grid::build(tube_spec(n_z), FIELDS).expect("valid spec");
    let f = EulerFields::resolve(&g).expect("fields registered");
    let eos = GammaLaw { gamma: SOD_GAMMA };
    fill_from_prim(&mut g, &f, &eos, |_, _, z| sod_ic(z));
    let op = closed_tube_op(eos);
    let steps = march_to(&op, &mut g, &f, 0.0, SOD_T_FINAL).expect("march");
    (g, f, steps)
}

/// ρ of ring 0 along z (the tube profile), plus a bitwise check that every
/// other ring and every component agrees with ring 0 exactly.
fn tube_profile(g: &Grid, f: &EulerFields) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, bool) {
    let n_z = g.spec().n_z;
    let ids = f.ids();
    let eos = GammaLaw { gamma: SOD_GAMMA };
    let mut rho = vec![0.0; n_z];
    let mut u_z = vec![0.0; n_z];
    let mut p = vec![0.0; n_z];
    let mut c_frac = vec![0.0; n_z];
    let mut uniform = true;
    for i_z in 0..n_z {
        let mut ring0: Option<Cons> = None;
        for i_r in 0..g.spec().n_r {
            let bi = g.brick_index(i_r, i_z).expect("active brick");
            let b = g.brick(bi);
            let idx = Grid::local_rz(i_r, i_z);
            let u: Cons = std::array::from_fn(|k| b.field(ids[k])[b.cell_index(0, idx)]);
            match &ring0 {
                None => {
                    let w = eos.prim_checked(&u).expect("physical state");
                    rho[i_z] = w[I_RHO];
                    u_z[i_z] = w[3];
                    p[i_z] = w[4];
                    c_frac[i_z] = w[I_RC];
                    ring0 = Some(u);
                }
                Some(r0) => {
                    if u.iter()
                        .zip(r0.iter())
                        .any(|(a, b)| a.to_bits() != b.to_bits())
                    {
                        uniform = false;
                    }
                }
            }
        }
    }
    (rho, u_z, p, c_frac, uniform)
}

fn l1_window(rho: &[f64], dz: f64, t: f64, exact: &RiemannSolution, lo: f64, hi: f64) -> f64 {
    let mut acc = 0.0f64;
    for (i, r) in rho.iter().enumerate() {
        let z = (i as f64 + 0.5) * dz;
        if z < lo || z > hi {
            continue;
        }
        let (re, _, _) = exact.sample((z - SOD_DIAPHRAGM_Z) / t);
        acc += (r - re).abs() * dz;
    }
    acc
}

/// The Sod refinement study: global + smooth-window L1(ρ) per level.
pub fn sod_convergence() -> SodStudy {
    let exact = sod_exact();
    let mut levels = Vec::new();
    for &n_z in &SOD_LEVELS {
        let (g, f, steps) = run_sod(n_z);
        let dz = g.spec().dz;
        let (rho, _, _, _, _) = tube_profile(&g, &f);
        levels.push(SodLevel {
            n_z,
            steps,
            l1_rho: l1_window(&rho, dz, SOD_T_FINAL, &exact, 0.0, 1.0),
            l1_fan: l1_window(&rho, dz, SOD_T_FINAL, &exact, FAN_WINDOW.0, FAN_WINDOW.1),
            l1_star_l: l1_window(
                &rho,
                dz,
                SOD_T_FINAL,
                &exact,
                STAR_L_WINDOW.0,
                STAR_L_WINDOW.1,
            ),
        });
    }
    SodStudy { levels }
}

/// Wave diagnostics at the finest ladder level.
pub fn sod_waves() -> SodWaves {
    let n_z = SOD_LEVELS[SOD_LEVELS.len() - 1];
    let exact = sod_exact();
    let (g, f, _) = run_sod(n_z);
    let dz = g.spec().dz;
    let (rho, u_z, p, c_frac, uniform) = tube_profile(&g, &f);

    // Shock front: scanning right-to-left, the first cell above the
    // midpoint between the undisturbed and shocked densities.
    let mid_shock = 0.5 * (SOD_RIGHT.rho + exact.rho_star_r);
    let z_shock_num = (0..n_z)
        .rev()
        .find(|&i| rho[i] > mid_shock)
        .map(|i| (i as f64 + 0.5) * dz)
        .expect("shock in domain");
    let z_shock_exact = SOD_DIAPHRAGM_Z + exact.right_shock_speed() * SOD_T_FINAL;

    // Contact: the composition-midpoint crossing (C goes 0 → 1).
    let z_contact_num = (0..n_z)
        .find(|&i| c_frac[i] > 0.5)
        .map(|i| (i as f64 + 0.5) * dz)
        .expect("contact in domain");
    let z_contact_exact = SOD_DIAPHRAGM_Z + exact.u_star * SOD_T_FINAL;

    let window_max = |vals: &[f64], target: f64, lo: f64, hi: f64| -> f64 {
        let mut worst = 0.0f64;
        for (i, v) in vals.iter().enumerate() {
            let z = (i as f64 + 0.5) * dz;
            if z >= lo && z <= hi {
                worst = worst.max((v - target).abs());
            }
        }
        worst
    };

    let mut c_min = f64::INFINITY;
    let mut c_max = f64::NEG_INFINITY;
    for &c in &c_frac {
        c_min = c_min.min(c);
        c_max = c_max.max(c);
    }

    SodWaves {
        n_z,
        shock_pos_err_cells: (z_shock_num - z_shock_exact).abs() / dz,
        contact_pos_err_cells: (z_contact_num - z_contact_exact).abs() / dz,
        rho_star_l_err: window_max(&rho, exact.rho_star_l, STAR_L_WINDOW.0, STAR_L_WINDOW.1),
        rho_star_r_err: window_max(&rho, exact.rho_star_r, STAR_R_WINDOW.0, STAR_R_WINDOW.1),
        u_star_err: window_max(&u_z, exact.u_star, STAR_L_WINDOW.0, STAR_R_WINDOW.1),
        p_star_err: window_max(&p, exact.p_star, STAR_L_WINDOW.0, STAR_R_WINDOW.1),
        c_min,
        c_max,
        radially_uniform_bitwise: uniform,
    }
}

// --- Smooth formal-order study (advected composition ramp) ------------------

/// The smooth study: a resolved tanh density/composition ramp advecting at
/// uniform (u, p) — the contact-family solution, exact profile known, no
/// interior extrema so the limiters stay inactive and the observed order is
/// the scheme's, not the limiter's.
pub const ADV_LEVELS: [usize; 4] = [50, 100, 200, 400];
pub const ADV_T_FINAL: f64 = 0.25;
pub const ADV_U: f64 = 1.0;
pub const ADV_P: f64 = 1.0;
pub const ADV_RAMP_Z0: f64 = 0.3;
pub const ADV_RAMP_WIDTH: f64 = 0.1;

fn adv_prim(z: f64, t: f64) -> Prim {
    let s = 0.5 * (1.0 + ((z - ADV_RAMP_Z0 - ADV_U * t) / ADV_RAMP_WIDTH).tanh());
    prim6(1.0 + 0.5 * s, 0.0, 0.0, ADV_U, ADV_P, s)
}

#[derive(Debug, Clone, PartialEq)]
pub struct AdvLevel {
    pub n_z: usize,
    pub l1_rho: f64,
    pub l1_c: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AdvStudy {
    pub levels: Vec<AdvLevel>,
}

impl AdvStudy {
    pub fn observed_orders(&self, err: impl Fn(&AdvLevel) -> f64) -> Vec<f64> {
        self.levels
            .windows(2)
            .map(|w| (err(&w[0]) / err(&w[1])).ln() / std::f64::consts::LN_2)
            .collect()
    }

    /// Mean observed order across the whole ladder (coarsest vs finest).
    pub fn mean_order(&self, err: impl Fn(&AdvLevel) -> f64) -> f64 {
        let first = &self.levels[0];
        let last = &self.levels[self.levels.len() - 1];
        (err(first) / err(last)).ln() / (last.n_z as f64 / first.n_z as f64).ln()
    }
}

pub fn advection_order() -> AdvStudy {
    let eos = GammaLaw { gamma: SOD_GAMMA };
    let exact_bc = |_r: f64, _th: f64, z: f64, t: f64| adv_prim(z, t);
    let mut levels = Vec::new();
    for &n_z in &ADV_LEVELS {
        let spec = GridSpec {
            r_min: 0.0,
            dr: 0.05,
            n_r: 2,
            z_min: 0.0,
            dz: 1.0 / n_z as f64,
            n_z,
            n_theta_max: 1,
            axisymmetry_assertion: true,
        };
        let mut g = Grid::build(spec, FIELDS).expect("valid spec");
        let f = EulerFields::resolve(&g).expect("fields");
        fill_from_prim(&mut g, &f, &eos, |_, _, z| adv_prim(z, 0.0));
        let op = Euler {
            eos,
            source: &ZERO_SRC,
            bcs: FlowBcs {
                r_inner: FlowBc::Reflecting,
                r_outer: FlowBc::Reflecting,
                z_lo: FlowBc::Prescribed(&exact_bc),
                z_hi: FlowBc::Prescribed(&exact_bc),
            },
            wall_normal: None,
            slip_wall_z_faces: true, // certified station behavior (slip everywhere)
        };
        march_to(&op, &mut g, &f, 0.0, ADV_T_FINAL).expect("march");
        let (rho, _, _, c_frac, _) = tube_profile_any(&g, &f);
        let dz = g.spec().dz;
        let mut l1_rho = 0.0;
        let mut l1_c = 0.0;
        for i in 0..n_z {
            let z = (i as f64 + 0.5) * dz;
            let we = adv_prim(z, ADV_T_FINAL);
            l1_rho += (rho[i] - we[I_RHO]).abs() * dz;
            l1_c += (c_frac[i] - we[I_RC]).abs() * dz;
        }
        levels.push(AdvLevel { n_z, l1_rho, l1_c });
    }
    AdvStudy { levels }
}

/// Ring-0 profile without the Sod-specific EOS assumption baked in
/// (same γ here, but kept separate for clarity).
fn tube_profile_any(g: &Grid, f: &EulerFields) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, bool) {
    tube_profile(g, f)
}

// --- Well-balance & conservation --------------------------------------------

/// Steps for the uniform-fixed-point and conservation studies.
pub const UNIFORM_STEPS: usize = 50;
pub const CONSERVATION_T_FINAL: f64 = 0.6;
pub const CONSERVATION_TOL_REL: f64 = 1e-12;

/// A uniform gas at rest must be a **bitwise** fixed point of the operator
/// — the well-balanced geometric sources against the exact cylindrical
/// metric, axis included (`axis = true`) or annulus with N_θ = 8
/// (`axis = false`, exercising the θ sweep). Returns true iff every field
/// is bit-identical after `UNIFORM_STEPS`.
pub fn uniform_state_is_bitwise_fixed_point(axis: bool) -> bool {
    let spec = if axis {
        GridSpec {
            r_min: 0.0,
            dr: 0.01,
            n_r: 16,
            z_min: 0.0,
            dz: 0.02,
            n_z: 8,
            n_theta_max: 1,
            axisymmetry_assertion: true,
        }
    } else {
        GridSpec {
            r_min: 0.3,
            dr: 0.01,
            n_r: 16,
            z_min: 0.0,
            dz: 0.02,
            n_z: 8,
            n_theta_max: 8,
            axisymmetry_assertion: false,
        }
    };
    let mut g = Grid::build(spec, FIELDS).expect("valid spec");
    let f = EulerFields::resolve(&g).expect("fields");
    let eos = GammaLaw { gamma: SOD_GAMMA };
    fill_from_prim(&mut g, &f, &eos, |_, _, _| {
        prim6(1.3, 0.0, 0.0, 0.0, 2.7, 0.5)
    });
    let before: Vec<Vec<u64>> = snapshot_bits(&g, &f);
    let op = closed_tube_op(eos);
    let mut sdc = Sdc::new();
    let flow = FlowClass {
        op: &op,
        fields: &f,
    };
    let dt = sdc.stable_dt(&g, &flow, CFL_FLOW).expect("dt");
    let mut t = 0.0;
    for _ in 0..UNIFORM_STEPS {
        sdc.step_flow(&mut g, &flow, t, dt).expect("step");
        t += dt;
    }
    snapshot_bits(&g, &f) == before
}

fn snapshot_bits(g: &Grid, f: &EulerFields) -> Vec<Vec<u64>> {
    let ids = f.ids();
    g.bricks()
        .iter()
        .map(|b| {
            ids.iter()
                .flat_map(|&id| b.field(id).iter().map(|v| v.to_bits()))
                .collect()
        })
        .collect()
}

/// Closed reflecting tube, Sod burst, marched to t = 0.6 — through the
/// first wall reflections. Returns (mass drift, energy drift) relative to
/// the initial totals; flux-form telescoping + wall fluxes from mirrored
/// states must hold both at round-off. (Momentum is NOT conserved in a
/// closed tube — the walls push back; that ledger is COUP-2 §3.1.2's
/// mount-reaction term, a later wave.)
pub fn closed_tube_conservation() -> Result<(f64, f64), SdcError> {
    let n_z = 200;
    let mut g = Grid::build(tube_spec(n_z), FIELDS).expect("valid spec");
    let f = EulerFields::resolve(&g).expect("fields");
    let eos = GammaLaw { gamma: SOD_GAMMA };
    fill_from_prim(&mut g, &f, &eos, |_, _, z| sod_ic(z));
    let ids = f.ids();
    let mass0 = g.reduce_volume_weighted(ids[I_RHO]);
    let energy0 = g.reduce_volume_weighted(ids[4]);
    let op = closed_tube_op(eos);
    march_to(&op, &mut g, &f, 0.0, CONSERVATION_T_FINAL)?;
    let mass1 = g.reduce_volume_weighted(ids[I_RHO]);
    let energy1 = g.reduce_volume_weighted(ids[4]);
    Ok((
        ((mass1 - mass0) / mass0).abs(),
        ((energy1 - energy0) / energy0).abs(),
    ))
}
