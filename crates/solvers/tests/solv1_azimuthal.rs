//! SOLV-1 §3.3 / FND-2 §3.2–3.4 (plan S8) — **the azimuthal capability
//! gates**: r = 0 axis crossing at N_θ > 1 (the θ↔θ+π parity-pair gather)
//! and mixed-N_θ brick refluxing (the ring-interface exchange). Plan §7
//! names this area "the historically buggiest class" — these are GATE
//! tests, not demos:
//!
//! * `axis_pulse_3d` — an acoustic pulse crosses the axis cleanly: exact
//!   conservation (the armed COUP-2 audit every step), the signal genuinely
//!   reaches the far side, mirror symmetry holds, and an axisymmetric pulse
//!   at N_θ = 8 reproduces the N_θ = 1 march **bitwise** per θ-plane (the
//!   metric factors differ by exact powers of two, so the arithmetic
//!   cancels exactly — any deviation is an axis-machinery defect).
//! * `theta_reflux` — conservation at N_θ interfaces to round-off with the
//!   audit armed; a uniform state is a bitwise fixed point on a mixed-N_θ
//!   world; an axisymmetric transient on a mixed-N_θ world reproduces the
//!   N_θ = 1 march bitwise (the fine-side-owns-the-flux aggregate is exact
//!   on axisymmetric data).
//! * the plan §3 **Δt-vs-N_θ-profile** gate: a coarse-inner N_θ(r) profile
//!   holds Δt against the uniform-fine θ-CFL crush at small radius.
//! * the **symmetry controller in anger**: collapse and expand decisions
//!   fire mid-march, the Euler workspace re-keys (staleness rebuild), and
//!   the march stays audited.
//! * bit-determinism at any thread count on a mixed-N_θ axis world.

use crucible_grid::{Grid, GridSpec, MomentumFields, N_THETA_GUARD, ThetaAction, ThetaController};
use crucible_solvers::euler::{
    Cons, Euler, EulerFields, FlowBc, FlowBcs, GammaLaw, I_RHO, NCOMP, Prim, fill_from_prim, prim6,
};
use crucible_solvers::sdc::{FlowClass, Sdc};
use crucible_solvers::station1_sod::FIELDS;

const GAMMA: f64 = 1.4;
const ZERO_SRC: fn(f64, f64, f64, f64) -> Cons = |_, _, _, _| [0.0; NCOMP];

fn closed_box_op(eos: GammaLaw) -> Euler<'static, GammaLaw> {
    Euler {
        eos,
        source: &ZERO_SRC,
        bcs: FlowBcs {
            r_inner: FlowBc::Reflecting,
            r_outer: FlowBc::Reflecting,
            z_lo: FlowBc::Reflecting,
            z_hi: FlowBc::Reflecting,
        },
        wall_normal: None,
        slip_wall_z_faces: true,
        combustion: None,
    }
}

fn axis_spec(n_theta_max: u32) -> GridSpec {
    GridSpec {
        r_min: 0.0,
        dr: 0.01,
        n_r: 16,
        z_min: 0.0,
        dz: 0.01,
        n_z: 8,
        n_theta_max,
        axisymmetry_assertion: n_theta_max == 1,
    }
}

/// Axisymmetric Gaussian pressure pulse centered on the axis.
fn axisym_pulse(r: f64, _th: f64, z: f64) -> Prim {
    let dp = 0.4 * (-((r / 0.04).powi(2) + ((z - 0.04) / 0.04).powi(2))).exp();
    prim6(1.0, 0.0, 0.0, 0.0, 1.0 + dp, 0.3)
}

/// Off-axis Gaussian pressure pulse centered at (r0, θ = 0) — genuinely 3-D;
/// symmetric under θ → −θ.
fn offaxis_pulse(r: f64, th: f64, z: f64) -> Prim {
    let (x, y) = (r * th.cos(), r * th.sin());
    let (x0, w) = (0.06, 0.03);
    let d2 = (x - x0).powi(2) + y * y + (z - 0.04).powi(2);
    let dp = 0.5 * (-(d2 / (w * w))).exp();
    prim6(1.0, 0.0, 0.0, 0.0, 1.0 + dp, 0.3)
}

/// Per-θ-plane bitwise snapshot: `out[bi][j][k·64 + local]` bits.
fn plane_bits(g: &Grid, f: &EulerFields) -> Vec<Vec<Vec<u64>>> {
    let ids = f.ids();
    g.bricks()
        .iter()
        .map(|b| {
            (0..b.n_theta())
                .map(|j| {
                    ids.iter()
                        .flat_map(|&id| {
                            (0..64).map(move |local| b.field(id)[b.cell_index(j, local)].to_bits())
                        })
                        .collect()
                })
                .collect()
        })
        .collect()
}

/// March `n` steps with the fine grid's own CFL Δt, feeding the SAME Δt to
/// a shadow N_θ = 1 grid; states must stay plane-wise bitwise equal.
fn march_pair_bitwise(
    fill: impl Fn(f64, f64, f64) -> Prim + Copy,
    mixed: impl Fn(&mut Grid),
    n_steps: usize,
) {
    let eos = GammaLaw { gamma: GAMMA };
    let mut g1 = Grid::build(axis_spec(1), FIELDS).expect("N_θ=1 grid");
    let mut g8 = Grid::build(axis_spec(8), FIELDS).expect("N_θ=8 grid");
    mixed(&mut g8);
    let f1 = EulerFields::resolve(&g1).expect("fields");
    let f8 = EulerFields::resolve(&g8).expect("fields");
    fill_from_prim(&mut g1, &f1, &eos, fill);
    fill_from_prim(&mut g8, &f8, &eos, fill);
    let op = closed_box_op(eos);
    let (mut sdc1, mut sdc8) = (Sdc::new(), Sdc::new());
    let flow1 = FlowClass {
        op: &op,
        fields: &f1,
    };
    let flow8 = FlowClass {
        op: &op,
        fields: &f8,
    };
    let mut t = 0.0;
    for step in 0..n_steps {
        let dt = sdc8.stable_dt(&g8, &flow8, 0.4).expect("dt");
        sdc8.step_flow(&mut g8, &flow8, t, dt).expect("3-D step");
        sdc1.step_flow(&mut g1, &flow1, t, dt).expect("2-D step");
        t += dt;
        let b1 = plane_bits(&g1, &f1);
        let b8 = plane_bits(&g8, &f8);
        for (bi, planes) in b8.iter().enumerate() {
            for (j, plane) in planes.iter().enumerate() {
                assert_eq!(
                    plane, &b1[bi][0],
                    "θ-plane {j} of brick {bi} diverged from the N_θ=1 march at step {step}"
                );
            }
        }
    }
}

// Covers FND-2 §3.2 (axis parity pairing at N_θ > 1) — named per META-2 §4.
#[test]
fn fnd2_s32_axis_pulse_axisymmetric_ntheta8_matches_ntheta1_bitwise() {
    // The exact powers-of-two metric cancellation makes the N_θ = 8 march of
    // axisymmetric data the N_θ = 1 march replicated per plane, THROUGH the
    // axis machinery (parity-pair gather active every step). 60 steps carries
    // the pulse through the axis and back off the outer wall.
    march_pair_bitwise(axisym_pulse, |_| {}, 60);
}

// Covers FND-2 §3.4 (ring-interface exchange) — the same bitwise argument
// with the refluxing active: brick (0,0) coarsened to N_θ = 4 creates one
// r-interface and one z-interface against its N_θ = 8 neighbors.
#[test]
fn fnd2_s34_theta_reflux_axisymmetric_transient_matches_ntheta1_bitwise() {
    march_pair_bitwise(
        axisym_pulse,
        |g| {
            let bi = g.brick_index(0, 0).expect("brick (0,0)");
            g.coarsen_theta(bi, None).expect("8 → 4");
        },
        60,
    );
}

// Covers FND-2 §3.4 / COUP-2 §3.1 (mixed-N_θ uniform fixed point).
#[test]
fn coup2_s31_uniform_state_is_bitwise_fixed_point_on_mixed_ntheta_world() {
    let eos = GammaLaw { gamma: GAMMA };
    let mut g = Grid::build(axis_spec(8), FIELDS).expect("grid");
    let bi = g.brick_index(0, 0).expect("brick (0,0)");
    g.coarsen_theta(bi, None).expect("8 → 4");
    let f = EulerFields::resolve(&g).expect("fields");
    fill_from_prim(&mut g, &f, &eos, |_, _, _| {
        prim6(1.3, 0.0, 0.0, 0.0, 2.7, 0.5)
    });
    let before = plane_bits(&g, &f);
    let op = closed_box_op(eos);
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
    assert_eq!(
        plane_bits(&g, &f),
        before,
        "uniform state drifted on the mixed-N_θ axis world"
    );
}

// Covers SOLV-1 §3.3 / FND-2 §3.2 (the S8 axis crossing, genuinely 3-D).
#[test]
fn solv1_s33_axis_pulse_3d_crosses_cleanly() {
    let eos = GammaLaw { gamma: GAMMA };
    let mut g = Grid::build(axis_spec(8), FIELDS).expect("grid");
    let f = EulerFields::resolve(&g).expect("fields");
    fill_from_prim(&mut g, &f, &eos, offaxis_pulse);
    let ids = f.ids();
    let mass0 = g.reduce_volume_weighted(ids[I_RHO]);
    let en0 = g.reduce_volume_weighted(ids[4]);

    // The far-side probe: max density excursion over the θ ≈ π half before
    // the march (the pulse is at θ = 0; the far side starts quiescent).
    let far_excursion = |g: &Grid| -> f64 {
        let mut worst = 0.0f64;
        g.for_each_active_cell(|cell| {
            let b = g.brick(cell.bi);
            let nt = b.n_theta();
            let j_quarter = nt / 4;
            if cell.i_theta >= nt / 2 - j_quarter && cell.i_theta < nt / 2 + j_quarter {
                let rho = b.field(ids[I_RHO])[cell.idx];
                worst = worst.max((rho - 1.0).abs());
            }
        });
        worst
    };
    assert!(
        far_excursion(&g) < 1e-6,
        "IC leaked to the far side: the probe is meaningless"
    );

    let op = closed_box_op(eos);
    let mut sdc = Sdc::new();
    let flow = FlowClass {
        op: &op,
        fields: &f,
    };
    // ~1.5 crossing times of the ~0.16 m domain at c ≈ √1.4.
    let t_final = 0.2;
    let mut t = 0.0;
    let mut steps = 0usize;
    while t < t_final {
        let dt = sdc.stable_dt(&g, &flow, 0.4).expect("dt").min(t_final - t);
        // The armed COUP-2 audit IS the conservation gate, every step.
        sdc.step_flow(&mut g, &flow, t, dt).expect("audited step");
        t += dt;
        steps += 1;
        assert!(steps < 20_000, "runaway march");
    }

    // (1) Exact conservation (closed box; reductions to round-off).
    let mass1 = g.reduce_volume_weighted(ids[I_RHO]);
    let en1 = g.reduce_volume_weighted(ids[4]);
    assert!(
        ((mass1 - mass0) / mass0).abs() < 1e-12,
        "mass drifted: {:.3e}",
        (mass1 - mass0) / mass0
    );
    assert!(
        ((en1 - en0) / en0).abs() < 1e-12,
        "energy drifted: {:.3e}",
        (en1 - en0) / en0
    );

    // (2) The pulse genuinely crossed the axis: the far side saw a real
    // signal (the crossing path through r = 0 is shorter than around).
    let far = far_excursion(&g);
    assert!(
        far > 1e-3,
        "no signal reached the far side (θ ≈ π): far excursion {far:.3e} — \
         the axis is acting as a wall"
    );

    // (3) Mirror symmetry: the IC is symmetric under θ → −θ, and the
    // operator must preserve it (ρ, u_r, u_z, p, C even; u_θ odd). The
    // sweeps' pencil arithmetic is not bit-symmetric under reflection, so
    // this is a tolerance gate, not a bitwise one.
    let mut worst_sym = 0.0f64;
    for b in g.bricks() {
        let nt = b.n_theta();
        for j in 0..nt / 2 {
            let jm = nt - 1 - j; // mirror partner of θ_j under θ → −θ
            for local in 0..64 {
                if b.mask() & (1u64 << local) == 0 {
                    continue;
                }
                for (k, &id) in ids.iter().enumerate() {
                    let a = b.field(id)[b.cell_index(j, local)];
                    let c = b.field(id)[b.cell_index(jm, local)];
                    let d = if k == 2 { (a + c).abs() } else { (a - c).abs() };
                    worst_sym = worst_sym.max(d);
                }
            }
        }
    }
    assert!(
        worst_sym < 1e-10,
        "θ → −θ mirror symmetry broken: worst {worst_sym:.3e}"
    );
}

// Covers FND-2 §3.4 (mixed-N_θ conservation on genuinely 3-D data): the
// audit is armed every step; the explicit reductions close the loop.
#[test]
fn fnd2_s34_theta_reflux_conserves_on_3d_data() {
    let eos = GammaLaw { gamma: GAMMA };
    // 16 z-cells = two brick rows, so coarsening (0,0) makes BOTH an
    // r-interface (against (1,0)) and a z-interface (against (0,1)).
    let spec = GridSpec {
        n_z: 16,
        ..axis_spec(8)
    };
    let mut g = Grid::build(spec, FIELDS).expect("grid");
    let bi = g.brick_index_by_coords(0, 0).expect("inner brick present");
    g.coarsen_theta(bi, None).expect("8 → 4");
    let f = EulerFields::resolve(&g).expect("fields");
    fill_from_prim(&mut g, &f, &eos, offaxis_pulse);
    let ids = f.ids();
    let mass0 = g.reduce_volume_weighted(ids[I_RHO]);
    let en0 = g.reduce_volume_weighted(ids[4]);
    let op = closed_box_op(eos);
    let mut sdc = Sdc::new();
    let flow = FlowClass {
        op: &op,
        fields: &f,
    };
    let mut t = 0.0;
    for _ in 0..80 {
        let dt = sdc.stable_dt(&g, &flow, 0.4).expect("dt");
        sdc.step_flow(&mut g, &flow, t, dt)
            .expect("audited step across N_θ interfaces");
        t += dt;
    }
    let mass1 = g.reduce_volume_weighted(ids[I_RHO]);
    let en1 = g.reduce_volume_weighted(ids[4]);
    assert!(
        ((mass1 - mass0) / mass0).abs() < 1e-12,
        "mass drifted across the N_θ interfaces: {:.3e}",
        (mass1 - mass0) / mass0
    );
    assert!(
        ((en1 - en0) / en0).abs() < 1e-12,
        "energy drifted across the N_θ interfaces: {:.3e}",
        (en1 - en0) / en0
    );
}

// Covers plan §3 (θ-CFL at inner rings) / COUP-3 §3.4 (the azimuthal Δt
// member per brick): the coarse-inner N_θ(r) profile holds Δt.
#[test]
fn coup3_s34_ntheta_profile_holds_dt_against_the_axis_theta_cfl() {
    let eos = GammaLaw { gamma: GAMMA };
    let fill = |_r: f64, _th: f64, _z: f64| prim6(1.0, 0.0, 0.0, 0.0, 1.0, 0.3);
    let dt_at = |profile: bool| -> f64 {
        let mut g = Grid::build(axis_spec(16), FIELDS).expect("grid");
        if profile {
            for bz in 0..1u32 {
                let bi = g.brick_index_by_coords(0, bz).expect("inner brick");
                g.coarsen_theta(bi, None).expect("16 → 8");
            }
        }
        let f = EulerFields::resolve(&g).expect("fields");
        fill_from_prim(&mut g, &f, &eos, fill);
        let op = closed_box_op(eos);
        let mut sdc = Sdc::new();
        let flow = FlowClass {
            op: &op,
            fields: &f,
        };
        sdc.stable_dt(&g, &flow, 0.4).expect("dt")
    };
    let dt_uniform = dt_at(false);
    let dt_profile = dt_at(true);
    println!("θ-CFL: uniform-16 dt = {dt_uniform:.3e}, profiled dt = {dt_profile:.3e}");
    // The innermost ring's arc doubles under the profile; the binding signal
    // is azimuthal there, so Δt must grow substantially (not exactly 2× —
    // the r/z members dilute it).
    assert!(
        dt_profile > 1.5 * dt_uniform,
        "N_θ(r) profile failed to hold Δt: {dt_profile:.3e} vs uniform {dt_uniform:.3e}"
    );
}

// Covers FND-2 §3.4 (the adaptive symmetry controller, exercised in anger:
// decisions fire mid-march, the flow workspace re-keys, the audit stays
// armed).
#[test]
fn fnd2_s34_symmetry_controller_collapses_and_expands_mid_march() {
    let eos = GammaLaw { gamma: GAMMA };
    let mut g = Grid::build(axis_spec(8), FIELDS).expect("grid");
    let f = EulerFields::resolve(&g).expect("fields");
    // Phase A: axisymmetric flow — the indicator reads ~0, the controller
    // must decide Collapse after N_DWELL evaluations, and the march must
    // continue clean on the collapsed world.
    fill_from_prim(&mut g, &f, &eos, axisym_pulse);
    let ids = f.ids();
    let mom = MomentumFields {
        rho: ids[0],
        mom_r: ids[1],
        mom_theta: ids[2],
        mom_z: ids[3],
    };
    let fields_floors: Vec<_> = ids.iter().map(|&id| (id, 1e-6)).collect();
    let op = closed_box_op(eos);
    let mut sdc = Sdc::new();
    let flow = FlowClass {
        op: &op,
        fields: &f,
    };
    let mut controllers = vec![ThetaController::default(); g.n_bricks()];
    let mut t = 0.0;
    let mut collapsed = 0usize;
    for step in 0..40 {
        let dt = sdc.stable_dt(&g, &flow, 0.4).expect("dt");
        sdc.step_flow(&mut g, &flow, t, dt).expect("audited step");
        t += dt;
        if step % 4 == 3 {
            // Accelerated cadence for the test (production: N_SYM_CADENCE).
            for (bi, ctl) in controllers.iter_mut().enumerate() {
                let a = g.symmetry_indicator(bi, &fields_floors).expect("indicator");
                if ctl.observe(a) == ThetaAction::Collapse && g.brick(bi).n_theta() > N_THETA_GUARD
                {
                    g.coarsen_theta(bi, Some(mom)).expect("collapse");
                    collapsed += 1;
                }
            }
        }
    }
    assert!(
        collapsed >= g.n_bricks(),
        "the controller failed to collapse the axisymmetric flow \
         (collapsed {collapsed} of {} bricks)",
        g.n_bricks()
    );
    assert!(
        g.bricks().iter().all(|b| b.n_theta() == N_THETA_GUARD),
        "collapse did not reach the guard resolution"
    );

    // Phase B: strong m = 1 content — the controller must decide Expand,
    // and the march must continue clean on the re-refined world.
    fill_from_prim(&mut g, &f, &eos, |r, th, z| {
        let mut w = axisym_pulse(r, th, z);
        w[0] *= 1.0 + 0.3 * th.cos();
        w
    });
    let mut expanded = 0usize;
    for _ in 0..8 {
        for (bi, ctl) in controllers.iter_mut().enumerate() {
            let a = g.symmetry_indicator(bi, &fields_floors).expect("indicator");
            if ctl.observe(a) == ThetaAction::Expand && g.brick(bi).n_theta() < g.spec().n_theta_max
            {
                g.refine_theta(bi).expect("expand");
                expanded += 1;
            }
        }
        let dt = sdc.stable_dt(&g, &flow, 0.4).expect("dt");
        sdc.step_flow(&mut g, &flow, t, dt)
            .expect("audited step after re-expansion");
        t += dt;
    }
    assert!(
        expanded >= g.n_bricks(),
        "the controller failed to re-expand under m = 1 content \
         (expanded {expanded} of {} bricks)",
        g.n_bricks()
    );
}

// Covers FND-2 §3.7 / SOLV-1 §3.6 (bit-determinism at any thread count on
// the S8 machinery: axis parity gather + mixed-N_θ refluxing active).
#[test]
fn solv1_s36_mixed_ntheta_axis_march_is_bitwise_thread_invariant() {
    let run = |threads: usize| -> Vec<Vec<Vec<u64>>> {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .expect("pool");
        pool.install(|| {
            let eos = GammaLaw { gamma: GAMMA };
            let mut g = Grid::build(axis_spec(8), FIELDS).expect("grid");
            let bi = g.brick_index(0, 0).expect("brick (0,0)");
            g.coarsen_theta(bi, None).expect("8 → 4");
            let f = EulerFields::resolve(&g).expect("fields");
            fill_from_prim(&mut g, &f, &eos, offaxis_pulse);
            let op = closed_box_op(eos);
            let mut sdc = Sdc::new();
            let flow = FlowClass {
                op: &op,
                fields: &f,
            };
            let mut t = 0.0;
            for _ in 0..25 {
                let dt = sdc.stable_dt(&g, &flow, 0.4).expect("dt");
                sdc.step_flow(&mut g, &flow, t, dt).expect("step");
                t += dt;
            }
            plane_bits(&g, &f)
        })
    };
    assert_eq!(
        run(1),
        run(4),
        "mixed-N_θ axis march differs between 1 and 4 threads"
    );
}

// Covers FND-2 §3.4's OTHER fix-up branch (S8 review: every earlier gate
// coarsened brick (0,0), so only the "right finer" aggregate ran): with the
// OUTER brick coarsened, the fine segment sits LEFT of the jump in every
// r-pencil — the "left finer" indices carry the flux.
#[test]
fn fnd2_s34_theta_reflux_left_finer_r_interface_matches_ntheta1_bitwise() {
    march_pair_bitwise(
        axisym_pulse,
        |g| {
            let bi = g.brick_index(8, 0).expect("outer brick (1,0)");
            g.coarsen_theta(bi, None).expect("8 → 4");
        },
        60,
    );
}

// The z-direction left-finer branch: coarsening the UPPER-z brick puts the
// fine segment first along every z-pencil of the inner column.
#[test]
fn fnd2_s34_theta_reflux_left_finer_z_interface_conserves() {
    let eos = GammaLaw { gamma: GAMMA };
    let spec = GridSpec {
        n_z: 16,
        ..axis_spec(8)
    };
    let mut g = Grid::build(spec, FIELDS).expect("grid");
    for (br, bz) in [(0u32, 1u32), (1, 1)] {
        let bi = g.brick_index_by_coords(br, bz).expect("upper brick");
        g.coarsen_theta(bi, None).expect("8 → 4");
    }
    let f = EulerFields::resolve(&g).expect("fields");
    fill_from_prim(&mut g, &f, &eos, offaxis_pulse);
    let ids = f.ids();
    let mass0 = g.reduce_volume_weighted(ids[I_RHO]);
    let en0 = g.reduce_volume_weighted(ids[4]);
    let op = closed_box_op(eos);
    let mut sdc = Sdc::new();
    let flow = FlowClass {
        op: &op,
        fields: &f,
    };
    let mut t = 0.0;
    for _ in 0..60 {
        let dt = sdc.stable_dt(&g, &flow, 0.4).expect("dt");
        sdc.step_flow(&mut g, &flow, t, dt)
            .expect("audited step across the left-finer z-interface");
        t += dt;
    }
    let mass1 = g.reduce_volume_weighted(ids[I_RHO]);
    let en1 = g.reduce_volume_weighted(ids[4]);
    assert!(((mass1 - mass0) / mass0).abs() < 1e-12);
    assert!(((en1 - en0) / en0).abs() < 1e-12);
}

// Covers FND-2 §3.2's u_θ SIGN at the axis (S8 review: the earlier gates
// cannot discriminate it — axisymmetric data makes the parity gather
// value-identical to the self-gather, and a wrong u_θ sign is itself
// θ→−θ mirror-symmetric AND exactly conservative). A uniform TRANSVERSE
// Cartesian flow (u_x = U ⇒ u_r = U·cosθ, u_θ = −U·sinθ) is an exact
// steady solution whose axis-crossing stencils read O(1)-wrong ghosts if
// either sign of the basis flip is wrong; with the right signs the state
// drifts only at truncation level.
#[test]
fn fnd2_s32_transverse_flow_through_axis_holds_steady() {
    let eos = GammaLaw { gamma: GAMMA };
    let mut g = Grid::build(axis_spec(16), FIELDS).expect("grid");
    let f = EulerFields::resolve(&g).expect("fields");
    const U: f64 = 0.2;
    let transverse = |_r: f64, th: f64, _z: f64| -> Prim {
        prim6(1.0, U * th.cos(), -U * th.sin(), 0.0, 1.0, 0.3)
    };
    fill_from_prim(&mut g, &f, &eos, transverse);
    let exact_bc = move |r: f64, th: f64, z: f64, _t: f64| transverse(r, th, z);
    let op = Euler {
        eos,
        source: &ZERO_SRC,
        bcs: FlowBcs {
            r_inner: FlowBc::Reflecting, // unused: the axis is geometric
            r_outer: FlowBc::Prescribed(&exact_bc),
            z_lo: FlowBc::Prescribed(&exact_bc),
            z_hi: FlowBc::Prescribed(&exact_bc),
        },
        wall_normal: None,
        slip_wall_z_faces: true,
        combustion: None,
    };
    let mut sdc = Sdc::new();
    let flow = FlowClass {
        op: &op,
        fields: &f,
    };
    let mut t = 0.0;
    for _ in 0..30 {
        let dt = sdc.stable_dt(&g, &flow, 0.4).expect("dt");
        sdc.step_flow(&mut g, &flow, t, dt).expect("audited step");
        t += dt;
    }
    let ids = f.ids();
    let mut worst = 0.0f64;
    g.for_each_active_cell(|cell| {
        let b = g.brick(cell.bi);
        let rho = b.field(ids[0])[cell.idx];
        let ut_num = b.field(ids[2])[cell.idx] / rho;
        let ut_exact = -U * cell.theta.sin();
        worst = worst.max((ut_num - ut_exact).abs());
    });
    println!("transverse-through-axis worst u_θ drift: {worst:.3e} (U = {U})");
    assert!(
        worst < 0.05 * U,
        "transverse flow drifted at the axis: worst u_θ error {worst:.3e} \
         — the parity gather's basis flip is suspect"
    );
}
