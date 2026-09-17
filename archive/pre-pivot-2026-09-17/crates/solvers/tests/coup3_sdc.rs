//! COUP-3 §6 — the S2 SDC-IMEX battery: temporal order (items 1/6),
//! class-`D` stiffness far beyond the explicit limit (item 6), the COUP-2
//! audit's halt mechanism (COUP-2 §6-1/2/5 subset), and 1-vs-N-thread bit
//! identity through the FULL coupled step (item 5; META-1 §2).
//!
//! (The other §6 items live where their subjects live: manufactured
//! diffusion orders + analytic anchors = the Goal-A battery, re-earned on
//! the implicit path; the cooled-duct fixture = station 4; conservation
//! through cut geometry = the cut-cell battery; the expander fixed point's
//! numerical-vs-physical diagnosis = the engine tests. All of those now
//! march with the audit armed every step — the closure check is not one
//! test, it is the standing condition of every march in the battery.)

use crucible_grid::{Grid, GridSpec};
use crucible_solvers::euler::{
    Cons, Euler, EulerFields, FlowBc, FlowBcs, GammaLaw, NCOMP, fill_from_prim, prim6,
};
use crucible_solvers::sdc::{AuditSpec, DiffusionClass, FlowClass, Sdc, SdcError};
use crucible_solvers::station4_cooled_wall as s4;
use crucible_solvers::{Bcs, Conduction, Domain, FaceBc, InteriorFaces};

const ZERO_SRC: fn(f64, f64, f64, f64) -> Cons = |_, _, _, _| [0.0; NCOMP];

fn smooth_tube() -> (Grid, EulerFields, Euler<'static>) {
    let spec = GridSpec {
        r_min: 0.0,
        dr: 0.025,
        n_r: 4,
        z_min: 0.0,
        dz: 1.0 / 48.0,
        n_z: 48,
        n_theta_max: 1,
        axisymmetry_assertion: true,
    };
    let mut g = Grid::build(spec, crucible_solvers::station1_sod::FIELDS).expect("grid");
    let f = EulerFields::resolve(&g).expect("fields");
    let eos = GammaLaw { gamma: 1.4 };
    // A smooth acoustic pulse in a closed tube: every family runs, nothing
    // shocks over the short window.
    fill_from_prim(&mut g, &f, &eos, |_, _, z| {
        let s = (std::f64::consts::TAU * z).sin();
        prim6(1.0 + 0.05 * s, 0.0, 0.0, 0.02 * s, 1.0 + 0.05 * s, 0.5)
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
        wall_normal: None,
        slip_wall_z_faces: true,
        combustion: None,
    };
    (g, f, op)
}

/// L2 of the field difference between two grids (same layout).
fn field_l2_diff(a: &Grid, b: &Grid, f: &EulerFields) -> f64 {
    let mut acc = 0.0f64;
    for (ba, bb) in a.bricks().iter().zip(b.bricks()) {
        for id in f.ids() {
            for (x, y) in ba.field(id).iter().zip(bb.field(id)) {
                acc += (x - y) * (x - y);
            }
        }
    }
    acc.sqrt()
}

/// COUP-3 §6-1: the SDC step recovers 2nd temporal order. Fixed grid, one
/// smooth transient, three step counts (N, 2N, 4N) to the same T —
/// Richardson on the run-to-run differences isolates the TEMPORAL error
/// (the identical spatial part cancels): order = log2(‖u_N − u_2N‖ /
/// ‖u_2N − u_4N‖) → 2.
#[test]
fn sdc_flow_temporal_order_is_two() {
    let t_final = 0.05;
    let run = |n_steps: usize| -> Grid {
        let (mut g, f, op) = smooth_tube();
        let mut sdc = Sdc::new();
        let flow = FlowClass {
            op: &op,
            fields: &f,
        };
        let dt = t_final / n_steps as f64;
        let mut t = 0.0;
        for _ in 0..n_steps {
            sdc.step_flow(&mut g, &flow, t, dt).expect("step");
            t += dt;
        }
        g
    };
    let (_, f, _) = smooth_tube();
    let (a, b, c) = (run(20), run(40), run(80));
    let e1 = field_l2_diff(&a, &b, &f);
    let e2 = field_l2_diff(&b, &c, &f);
    let order = (e1 / e2).log2();
    println!("SDC flow temporal order: {order:.3} (diffs {e1:.3e}, {e2:.3e})");
    assert!(
        (1.8..=2.3).contains(&order),
        "temporal order {order:.3} outside [1.8, 2.3]"
    );
}

fn annulus_op(source: &dyn Fn(f64, f64, f64, f64) -> f64) -> Conduction<'_> {
    static HOT: fn(f64, f64, f64, f64) -> f64 = |_, _, _, _| 500.0;
    static COLD: fn(f64, f64, f64, f64) -> f64 = |_, _, _, _| 300.0;
    Conduction {
        kappa: 20.0,
        rho_cp: 4.0e6,
        source,
        domain: Domain::FlowActive,
        interior: InteriorFaces::refuse(),
        bcs: Bcs {
            r_inner: FaceBc::Dirichlet(&HOT),
            r_outer: FaceBc::Dirichlet(&COLD),
            z_lo: FaceBc::HeatFlux(0.0),
            z_hi: FaceBc::HeatFlux(0.0),
        },
    }
}

fn annulus_grid() -> (Grid, crucible_grid::FieldId, crucible_grid::FieldId) {
    let spec = GridSpec {
        r_min: 0.05,
        dr: 0.1 / 32.0,
        n_r: 32,
        z_min: 0.0,
        dz: 0.05,
        n_z: 2,
        n_theta_max: 1,
        axisymmetry_assertion: true,
    };
    let mut g = Grid::build(spec, &["T", "rate"]).expect("grid");
    let t_id = g.field_id("T").unwrap();
    let rate_id = g.field_id("rate").unwrap();
    // A steep initial profile: the transient carries real temporal content.
    g.fill_field(t_id, |r, _, _| {
        400.0 + 80.0 * ((r - 0.05) / 0.1 * std::f64::consts::PI * 3.0).sin()
    });
    (g, t_id, rate_id)
}

/// COUP-3 §6-1/§6-6 (temporal half): the class-`D` implicit path is 2nd
/// order in time — Richardson at fixed grid, with every dt far above the
/// explicit stability bound (the regime the class exists for).
#[test]
fn sdc_diffusion_temporal_order_is_two() {
    let zero = |_: f64, _: f64, _: f64, _: f64| 0.0;
    let op = annulus_op(&zero);
    let (probe, _, _) = annulus_grid();
    let dt_expl = op.stable_dt(&probe, 1.0);
    let t_final = 128.0 * dt_expl; // a real transient window
    let run = |n_steps: usize| -> Grid {
        let (mut g, t_id, rate_id) = annulus_grid();
        let dc = DiffusionClass {
            op: &op,
            t_field: t_id,
            scratch_field: rate_id,
        };
        Sdc::new()
            .advance_diffusion(&mut g, &dc, 0.0, t_final / n_steps as f64, n_steps)
            .expect("march");
        g
    };
    // dt = 32×, 16×, 8× the explicit bound — all implicit-only territory.
    let (a, b, c) = (run(4), run(8), run(16));
    let (_, t_id, _) = annulus_grid();
    let diff = |x: &Grid, y: &Grid| -> f64 {
        let mut acc = 0.0f64;
        for (bx, by) in x.bricks().iter().zip(y.bricks()) {
            for (u, v) in bx.field(t_id).iter().zip(by.field(t_id)) {
                acc += (u - v) * (u - v);
            }
        }
        acc.sqrt()
    };
    let (e1, e2) = (diff(&a, &b), diff(&b, &c));
    let order = (e1 / e2).log2();
    println!("class-D temporal order: {order:.3} (diffs {e1:.3e}, {e2:.3e})");
    // ≥ 2 is the claim (COUP-3 §6-1). On this LINEAR fixture the truncated
    // sweep composition matches the trapezoidal amplification through z³,
    // so the measured order legitimately EXCEEDS 2 (measured ~3.0); the
    // nonlinear flow test above is the strict 2nd-order gate.
    assert!(
        (1.7..=3.5).contains(&order),
        "class-D temporal order {order:.3} outside [1.7, 3.5]"
    );
}

/// COUP-3 §6-6 (stiffness half): the fixed-cycle CG solve stays stable and
/// steady-state-accurate at 512× the explicit bound — where any explicit
/// or Jacobi/cell-local treatment diverges outright. (The committed Goal-A
/// anchors march at 4–32×; this probes deep into the stiff regime.)
#[test]
fn class_d_is_stable_and_accurate_far_beyond_the_explicit_limit() {
    let zero = |_: f64, _: f64, _: f64, _: f64| 0.0;
    let op = annulus_op(&zero);
    let (mut g, t_id, rate_id) = annulus_grid();
    let dt = 512.0 * op.stable_dt(&g, 1.0);
    let tau = 0.1 * 0.1 * op.rho_cp / op.kappa;
    let n_steps = (8.0 * tau / dt).ceil() as usize;
    let dc = DiffusionClass {
        op: &op,
        t_field: t_id,
        scratch_field: rate_id,
    };
    Sdc::new()
        .advance_diffusion(&mut g, &dc, 0.0, dt, n_steps)
        .expect("march");
    // Steady log profile between the Dirichlet walls (r1 = 0.05, r2 = 0.15).
    let (t1, t2, r1, r2) = (500.0, 300.0, 0.05, 0.15);
    let exact = |r: f64| (t1 * (r2 / r).ln() + t2 * (r / r1).ln()) / (r2 / r1).ln();
    let mut worst = 0.0f64;
    for i_r in 0..32 {
        let r = g.r_center(i_r);
        let bi = g.brick_index(i_r, 0).unwrap();
        let b = g.brick(bi);
        let t = b.field(t_id)[b.cell_index(0, Grid::local_rz(i_r, 0))];
        worst = worst.max((t - exact(r)).abs() / (t1 - t2));
    }
    println!("512× explicit limit: {n_steps} steps, worst rel {worst:.2e}");
    assert!(worst < 2e-3, "steady profile off by {worst:.3e}");
}

/// COUP-2 §6-1/§6-2 made explicit: a marching step's audit rows exist,
/// carry real port/source activity, and close within the derived
/// tolerance (the assertion every march makes implicitly — here read out).
#[test]
fn audit_rows_close_on_a_real_transient() {
    let (mut g, f, op) = smooth_tube();
    let mut sdc = Sdc::new();
    let flow = FlowClass {
        op: &op,
        fields: &f,
    };
    let dt = sdc.stable_dt(&g, &flow, 0.4).expect("dt");
    let report = sdc.step_flow(&mut g, &flow, 0.0, dt).expect("step");
    assert_eq!(report.audit.len(), NCOMP, "one row per conserved component");
    for row in &report.audit {
        let gap = (row.delta - row.applied).abs();
        assert!(
            gap <= row.tol,
            "{}: gap {gap:.3e} vs tol {:.3e}",
            row.quantity,
            row.tol
        );
        // A zero tolerance is legal ONLY for an identically-zero quantity
        // (θ-momentum on this no-swirl fixture: stored, ports, and sources
        // all exactly 0 — the identity is 0 = 0).
        assert!(
            row.tol > 0.0 || (row.delta == 0.0 && row.applied == 0.0),
            "{}: vanishing tolerance on a live quantity",
            row.quantity
        );
    }
    // The momentum rows must show REAL wall/geometric activity (walls push
    // back; the ledger is not vacuously zero).
    let mz = &report.audit[3];
    assert!(
        mz.applied.abs() > 0.0 || mz.delta.abs() > 0.0,
        "z-momentum ledger vacuous"
    );
}

/// COUP-2 §2: a violation beyond TOL_AUDIT halts with a diagnosis. Planted
/// by shrinking the safety factor below the genuine rounding floor — the
/// halt path itself is what's under test (the physics is untouched).
#[test]
fn audit_violation_halts_with_diagnosis() {
    let (mut g, f, op) = smooth_tube();
    let mut sdc = Sdc::with_audit(AuditSpec {
        k_audit: 1e-6, // far below the Higham √N floor — must trip
        ref_scale: [0.0; 4],
    });
    let flow = FlowClass {
        op: &op,
        fields: &f,
    };
    let dt = sdc.stable_dt(&g, &flow, 0.4).expect("dt");
    let mut tripped = false;
    let mut t = 0.0;
    for _ in 0..50 {
        match sdc.step_flow(&mut g, &flow, t, dt) {
            Ok(_) => t += dt,
            Err(SdcError::AuditViolation { quantity, tol, .. }) => {
                println!("tripped on {quantity} at tol {tol:.3e}");
                tripped = true;
                break;
            }
            Err(e) => panic!("wrong error class: {e}"),
        }
    }
    assert!(
        tripped,
        "a 1e-6 safety factor must sit below real rounding on a 50-step march"
    );
}

/// META-1 §2 / COUP-3 §6-5: the FULL coupled step (explicit sweeps +
/// class-D CG + Robin-Robin exchange + audit reductions) is bit-identical
/// at any thread count — the station-4 duct marched on 1 vs 4 rayon
/// threads.
#[test]
fn coupled_step_is_bit_identical_at_any_thread_count() {
    let run = |threads: usize| -> Vec<u64> {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .expect("pool");
        pool.install(|| {
            let (duct, _, _, _) = s4::run_duct_for(5.0e-4, 1.0e-4).expect("march");
            let mut bits = Vec::new();
            for b in duct.grid.bricks() {
                for id in duct.flow.ids() {
                    bits.extend(b.field(id).iter().map(|v| v.to_bits()));
                }
                bits.extend(b.field(duct.t_solid).iter().map(|v| v.to_bits()));
            }
            bits
        })
    };
    let one = run(1);
    let four = run(4);
    assert_eq!(one.len(), four.len());
    let diff = one.iter().zip(&four).filter(|(a, b)| a != b).count();
    assert_eq!(
        diff, 0,
        "{diff} values differ between 1 and 4 threads through the coupled step"
    );
}
