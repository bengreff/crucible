//! GOAL-B station 4 — the cooled wall (SOLV-1 §3.5, COUP-2 §3.5 subset).
//! Criteria are the named constants in `station4_cooled_wall.rs`; the
//! certificate binary prints the same numbers this battery asserts.

use crucible_solvers::station4_cooled_wall::*;
use crucible_solvers::{Bcs, Conduction, Domain, FaceBc, InteriorFaces, SolverError};

// Acceptance gates are the named shared constants in station4_cooled_wall.rs
// (session-6 convention: one owner for test + certificate binary).

#[test]
fn station4_duct_steady_ledger_and_series_resistance_oracle() {
    let (duct, rec, steps, resid) = run_duct().expect("coupled march");
    let r = duct_report(&duct, &rec, steps, resid);
    println!("{r:#?}");
    assert!(r.resid < RESID_MAX, "not steady: {:.3e}", r.resid);

    // Three independent steady rates must agree: gas enthalpy deficit,
    // wall-face exchange, coolant extraction.
    assert!(((r.wall_watts - r.gas_watts) / r.wall_watts).abs() < LEDGER_REL_TOL);
    assert!(((r.wall_watts - r.coolant_watts) / r.wall_watts).abs() < LEDGER_REL_TOL);

    // Series-resistance oracle, pointwise past the entrance band (one face
    // per i_z on the straight liner).
    assert_eq!(duct.faces.len(), N_Z);
    assert!(
        r.oracle_worst < ORACLE_REL_TOL,
        "oracle deviation {:.4}",
        r.oracle_worst
    );

    // Physics rows the certificate quotes: near-wall recovery sits inside
    // the declared band and the flux is rocket-scale (~1 MW/m²).
    assert!(
        r.t_aw_ratio_mid > T_AW_RATIO_BAND.0 && r.t_aw_ratio_mid < T_AW_RATIO_BAND.1,
        "near-wall recovery ratio {:.3} outside the declared band",
        r.t_aw_ratio_mid
    );
    assert!(
        r.q_mid > Q_MID_MIN,
        "rocket-scale flux expected, got {}",
        r.q_mid
    );

    // Review regression (session-10 review): the solid operator must never
    // touch T_solid outside its own domain — those cells hold their exact
    // initial value (0.0; `fill_solid` writes solid cells only) after the
    // whole march. The old unmasked pass-2 integrated stale rates there.
    for i_z in 0..N_Z {
        for i_r in 0..N_R_GAS {
            let v = cell_value(&duct.grid, duct.t_solid, i_r, i_z, 0);
            assert_eq!(v, 0.0, "T_solid drifted on gas cell ({i_r}, {i_z})");
        }
    }
}

#[test]
fn station4_uniform_rest_is_a_bitwise_fixed_point() {
    // Gas at rest at T_COOL, solid at T_COOL, coolant at T_COOL: every
    // exchange is exactly zero and the coupled step must not move one bit.
    let mut duct = build_duct();
    let rest: crucible_solvers::euler::Prim = [
        T_COOL * 0.0 + P_IN / (R_SPECIFIC * T_COOL),
        0.0,
        0.0,
        0.0,
        P_IN,
        0.0,
    ];
    let eos = duct.eos;
    crucible_solvers::euler::fill_from_prim(&mut duct.grid, &duct.flow, &eos, |_, _, _| rest);
    fill_solid(&mut duct.grid, duct.t_solid, T_COOL);

    // Sealed box: reflecting everywhere (no inflow/outflow), same wall law.
    static ZERO_SRC: fn(f64, f64, f64, f64) -> crucible_solvers::euler::Cons =
        |_, _, _, _| [0.0; crucible_solvers::euler::NCOMP];
    let op = crucible_solvers::euler::Euler {
        eos,
        source: &ZERO_SRC,
        bcs: crucible_solvers::euler::FlowBcs {
            r_inner: crucible_solvers::euler::FlowBc::Reflecting,
            r_outer: crucible_solvers::euler::FlowBc::Reflecting,
            z_lo: crucible_solvers::euler::FlowBc::Reflecting,
            z_hi: crucible_solvers::euler::FlowBc::Reflecting,
        },
        wall_normal: None,
    };
    // The fixture coolant is T_COOL already; the exchange at equal
    // temperatures is exactly 0 (unit-tested in wall_heat) — so the march
    // must be bit-frozen.
    let before: Vec<Vec<Vec<f64>>> = duct
        .grid
        .bricks()
        .iter()
        .map(|b| {
            let mut fields: Vec<Vec<f64>> = duct
                .flow
                .ids()
                .iter()
                .map(|id| b.field(*id).to_vec())
                .collect();
            fields.push(b.field(duct.t_solid).to_vec());
            fields
        })
        .collect();
    let mut rec = ExchangeRecord::default();
    let mut t = 0.0;
    for _ in 0..200 {
        t += coupled_step(&mut duct, &op, t, f64::INFINITY, &mut rec).expect("step");
    }
    for (bi, b) in duct.grid.bricks().iter().enumerate() {
        for (k, id) in duct.flow.ids().iter().enumerate() {
            assert_eq!(b.field(*id), &before[bi][k][..], "gas field {k} drifted");
        }
        assert_eq!(
            b.field(duct.t_solid),
            &before[bi][crucible_solvers::euler::NCOMP][..],
            "solid T drifted"
        );
    }
    assert_eq!(rec.wall_joules, 0.0, "exchange must be exactly zero");
}

#[test]
fn station4_rerun_is_bit_identical() {
    // Short window: bit-identity needs a real march, not a steady one.
    let (a, _, _, _) = run_duct_for(5.0e-4, 1.0e-4).expect("first run");
    let (b, _, _, _) = run_duct_for(5.0e-4, 1.0e-4).expect("second run");
    for (ba, bb) in a.grid.bricks().iter().zip(b.grid.bricks()) {
        for id in a.flow.ids() {
            assert_eq!(ba.field(id), bb.field(id));
        }
        assert_eq!(ba.field(a.t_solid), bb.field(b.t_solid));
    }
}

#[test]
fn station4_stepped_cavity_stair_interface_conserves_to_round_off() {
    // Hot gas at rest in a closed stepped cavity, cold stair liner,
    // insulated exterior: the interface presents r-faces, z-faces at the
    // steps, and solid↔exterior faces — the general machinery. Whatever
    // the gas loses the liner must hold, to accumulation round-off (the
    // same q·A applied to both sides; face areas identical by the
    // session-7 `face_radius` single-owner guarantee).
    let (mut duct, op) = build_stepped_cavity();
    let dirs: std::collections::BTreeSet<_> =
        duct.faces.iter().map(|f| format!("{:?}", f.dir)).collect();
    assert!(
        dirs.len() >= 2,
        "the stepped contour must exercise multiple face directions, got {dirs:?}"
    );
    let (e_gas0, e_solid0) = cavity_energies(&duct);
    let mut rec = ExchangeRecord::default();
    let mut t = 0.0;
    for _ in 0..500 {
        t += stepped_cavity_step(&mut duct, &op, t, &mut rec).expect("step");
    }
    let (e_gas1, e_solid1) = cavity_energies(&duct);
    let lost = e_gas0 - e_gas1;
    let gained = e_solid1 - e_solid0;
    println!(
        "cavity: gas lost {lost:.6e} J, solid gained {gained:.6e} J, \
         rel mismatch {:.2e}",
        ((lost - gained) / gained).abs()
    );
    assert!(gained > 0.0, "the cold liner must heat up");
    assert!(
        ((lost - gained) / gained).abs() < CAVITY_CONSERVATION_TOL,
        "stair-interface exchange must conserve to round-off"
    );
}

#[test]
fn station4_solid_step_without_gas_closure_refuses() {
    let mut duct = build_duct();
    let zero = |_: f64, _: f64, _: f64, _: f64| 0.0;
    let hold = |_: f64, _: f64, _: f64, _: f64| T_SOLID_INIT;
    let op = Conduction {
        kappa: KAPPA_S,
        rho_cp: RHO_CP_S,
        source: &zero,
        domain: Domain::Solid,
        interior: InteriorFaces {
            gas: None, // the refusal under test
            exterior: Some(FaceBc::HeatFlux(0.0)),
        },
        bcs: Bcs {
            r_inner: FaceBc::Dirichlet(&hold),
            r_outer: FaceBc::Dirichlet(&hold),
            z_lo: FaceBc::HeatFlux(0.0),
            z_hi: FaceBc::HeatFlux(0.0),
        },
    };
    let (tf, rf) = (duct.t_solid, duct.rate_solid);
    match op.step(&mut duct.grid, tf, rf, 0.0, 1e-9) {
        Err(SolverError::UnhandledInteriorFace { .. }) => {}
        other => panic!("expected UnhandledInteriorFace, got {other:?}"),
    }
}

#[test]
fn station4_robin_annulus_matches_the_exact_steady_profile() {
    // Conduction-only anchor for the new Robin face: annulus r ∈ [r1, r2],
    // Dirichlet T1 inside, Robin (h, T∞) outside. Exact steady solution:
    // T(r) = T1 + (T∞ − T1)·ln(r/r1) / (ln(r2/r1) + κ/(h·r2)).
    use crucible_grid::{Grid, GridSpec};
    let (r1, r2) = (0.05, 0.09);
    let (t1, t_inf, h, kappa) = (500.0, 300.0, 800.0, 15.0);
    let n_r = 32;
    let dr = (r2 - r1) / n_r as f64;
    let spec = GridSpec {
        r_min: r1,
        dr,
        n_r,
        z_min: 0.0,
        dz: 0.01,
        n_z: 4,
        n_theta_max: 1,
        axisymmetry_assertion: true,
    };
    let mut g = Grid::build(spec, &["T", "rate"]).expect("grid");
    let t_id = g.field_id("T").unwrap();
    let rate_id = g.field_id("rate").unwrap();
    g.fill_field(t_id, |_, _, _| t_inf);
    let zero = |_: f64, _: f64, _: f64, _: f64| 0.0;
    let hot = move |_: f64, _: f64, _: f64, _: f64| t1;
    let op = Conduction {
        kappa,
        rho_cp: 1000.0,
        source: &zero,
        domain: Domain::FlowActive,
        interior: InteriorFaces::refuse(),
        bcs: Bcs {
            r_inner: FaceBc::Dirichlet(&hot),
            r_outer: FaceBc::Robin { h, t_inf },
            z_lo: FaceBc::HeatFlux(0.0),
            z_hi: FaceBc::HeatFlux(0.0),
        },
    };
    let dt = op.stable_dt(&g, 0.9);
    // Diffusive settling: a few τ = (r2−r1)²/α.
    let tau = (r2 - r1) * (r2 - r1) * 1000.0 / kappa;
    let steps = (6.0 * tau / dt).ceil() as usize;
    op.advance(&mut g, t_id, rate_id, 0.0, dt, steps)
        .expect("march");

    let denom = (r2 / r1).ln() + kappa / (h * r2);
    let mut worst = 0.0f64;
    for i_r in 0..n_r {
        let r = r1 + (i_r as f64 + 0.5) * dr;
        let exact = t1 + (t_inf - t1) * (r / r1).ln() / denom;
        let got = cell_value(&g, t_id, i_r, 2, 0);
        worst = worst.max((got - exact).abs() / (t1 - t_inf).abs());
    }
    println!("robin annulus worst rel error = {worst:.2e} ({steps} steps)");
    assert!(worst < ANNULUS_ROBIN_TOL, "worst {worst:.3e}");
}
