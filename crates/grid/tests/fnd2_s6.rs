//! FND-2 v0.5 §6 validation plan — the items testable on the structural
//! core (no solvers yet). Numbering follows the doc:
//! - §6-3 determinism (repeat + partial-order independence);
//! - §6-7 conservative θ-coarsen/refine round trip;
//! - §6-8 collapse thermalization accounting;
//! - §6-9 guard-band re-expansion fires (+ axis metric treatment);
//! - plus metric exactness, the guard/assertion rules, and the post-review
//!   fail-loud regressions (NaN spec, NaN indicator).
//!
//! Deferred with their owners: §6-1 voxelization (FND-3), §6-2 conservation
//! sweep + flux closure across an N_θ jump (conduction session), §6-4
//! recession (SOLV-8), §6-5 tiles (grouping wave), §6-6 multi-material
//! (FND-7/SOLV-1 wave).

use crucible_grid::{
    BRICK_CELLS, Grid, GridError, GridSpec, MomentumFields, N_THETA_GUARD, TAU_EXPAND, ThetaAction,
    ThetaController, morton2, tree_combine,
};

fn spec(n_r: usize, n_z: usize, n_theta_max: u32) -> GridSpec {
    GridSpec {
        r_min: 0.0,
        dr: 0.01,
        n_r,
        z_min: 0.0,
        dz: 0.02,
        n_z,
        n_theta_max,
        axisymmetry_assertion: false,
    }
}

const FIELDS: &[&str] = &["rho", "mom_r", "mom_theta", "mom_z", "rho_e"];

fn momentum(g: &Grid) -> MomentumFields {
    MomentumFields {
        rho: g.field_id("rho").unwrap(),
        mom_r: g.field_id("mom_r").unwrap(),
        mom_theta: g.field_id("mom_theta").unwrap(),
        mom_z: g.field_id("mom_z").unwrap(),
    }
}

/// Deterministic pseudo-data: a fixed function of the cell position (no RNG
/// — counter-style, reproducible by construction).
fn fill_deterministic(g: &mut Grid) {
    for (k, name) in FIELDS.iter().enumerate() {
        let id = g.field_id(name).unwrap();
        g.fill_field(id, move |r, theta, z| {
            let x = r * 37.0 + theta * 1.1 + z * 53.0 + k as f64;
            1.5 + x.sin() * 0.25
        });
    }
}

// --- Metric exactness (§3.2) -------------------------------------------------

#[test]
fn fnd2_s32_ring_metric_sums_to_exact_cylinder_volume() {
    let g = Grid::build(spec(16, 4, 8), FIELDS).unwrap();
    // Σ over all cells of one z-layer & all θ = π R² Δz per layer.
    let mut v_sum = 0.0f64;
    for i_r in 0..16 {
        v_sum += g.cell_volume(i_r, 8) * 8.0;
    }
    let r_outer = 16.0 * 0.01;
    let exact = std::f64::consts::PI * r_outer * r_outer * 0.02;
    assert!(
        ((v_sum - exact) / exact).abs() < 1e-14,
        "ring volumes must telescope to the exact cylinder: {v_sum} vs {exact}"
    );

    // Axis face has exactly zero area — drops out of any flux stencil.
    assert_eq!(g.face_area_r(0, false, 8), 0.0);
    // Parity pairing: cross-axis neighbor is θ+π.
    assert_eq!(Grid::axis_pair(0, 8), 4);
    assert_eq!(Grid::axis_pair(5, 8), 1);
}

// --- §6-3: determinism -------------------------------------------------------

#[test]
fn fnd2_s6_3_reductions_are_bit_identical_and_partial_order_independent() {
    let mut g = Grid::build(spec(20, 12, 8), FIELDS).unwrap();
    fill_deterministic(&mut g);
    let f = g.field_id("rho_e").unwrap();

    // Repeat runs: bit-identical.
    let a = g.reduce_volume_weighted(f);
    let b = g.reduce_volume_weighted(f);
    assert_eq!(a.to_bits(), b.to_bits());

    // The tree combine is a pure function of the partial slots — identical
    // no matter which "thread" computed which partial (§3.7).
    let partials_fwd: Vec<f64> = (0..8).map(|i| (i as f64 + 0.1).sin()).collect();
    let mut partials_rev = vec![0.0f64; 8];
    for i in (0..8).rev() {
        partials_rev[i] = (i as f64 + 0.1).sin();
    }
    assert_eq!(
        tree_combine(&partials_fwd).to_bits(),
        tree_combine(&partials_rev).to_bits()
    );

    // Morton order is canonical: bricks are strictly sorted.
    assert!(g.bricks().windows(2).all(|w| w[0].morton() < w[1].morton()));
    assert_eq!(morton2(3, 5), 0b100111);
}

// --- §6-7: conservative θ-coarsen/refine round trip --------------------------

#[test]
fn fnd2_s6_7_theta_round_trip_preserves_ring_integrals() {
    let mut g = Grid::build(spec(10, 6, 16), FIELDS).unwrap();
    fill_deterministic(&mut g);
    let ids: Vec<_> = FIELDS.iter().map(|f| g.field_id(f).unwrap()).collect();

    let before: Vec<f64> = ids.iter().map(|f| g.reduce_volume_weighted(*f)).collect();
    let mom = momentum(&g);
    for bi in 0..g.n_bricks() {
        g.coarsen_theta(bi, Some(mom)).expect("16 → 8");
        g.coarsen_theta(bi, Some(mom)).expect("8 → 4");
    }
    let coarse: Vec<f64> = ids.iter().map(|f| g.reduce_volume_weighted(*f)).collect();
    for bi in 0..g.n_bricks() {
        g.refine_theta(bi).expect("4 → 8");
        g.refine_theta(bi).expect("8 → 16");
    }
    let after: Vec<f64> = ids.iter().map(|f| g.reduce_volume_weighted(*f)).collect();

    for ((name, b), (c, a)) in FIELDS.iter().zip(&before).zip(coarse.iter().zip(&after)) {
        let rel = |x: f64, y: f64| ((x - y) / y).abs();
        assert!(rel(*b, *c) < 1e-13, "{name}: coarsen lost {b} → {c}");
        assert!(rel(*b, *a) < 1e-13, "{name}: round trip lost {b} → {a}");
    }
}

// --- §6-8: collapse thermalization accounting --------------------------------

#[test]
fn fnd2_s6_8_collapse_logs_thermalized_ke_exactly() {
    let mut g = Grid::build(spec(8, 8, 8), FIELDS).unwrap();
    let mom = momentum(&g);
    let rho_id = g.field_id("rho").unwrap();
    let mt_id = g.field_id("mom_theta").unwrap();

    // Uniform ρ = 2, seeded m = 1 azimuthal velocity perturbation in ρu_θ.
    g.fill_field(rho_id, |_, _, _| 2.0);
    g.fill_field(mt_id, |_, theta, _| 0.3 * theta.sin());

    // Independent KE bookkeeping: resolved KE before and after must differ
    // by exactly the logged ΔKE.
    let resolved_ke = |g: &Grid| -> f64 {
        let mut ke = 0.0f64;
        g.for_each_active_cell(|c| {
            let b = g.brick(c.bi);
            let v = g.cell_volume(c.i_r, b.n_theta());
            let m = b.field(mt_id)[c.idx];
            ke += v * m * m / (2.0 * b.field(rho_id)[c.idx]);
        });
        ke
    };

    let ke_before = resolved_ke(&g);
    let mass_before = g.reduce_volume_weighted(rho_id);
    let mut logged = 0.0f64;
    for bi in 0..g.n_bricks() {
        logged += g
            .coarsen_theta(bi, Some(mom))
            .expect("8 → 4")
            .thermalized_ke;
    }
    let ke_after = resolved_ke(&g);
    let mass_after = g.reduce_volume_weighted(rho_id);

    assert!(
        logged > 0.0,
        "an m=1 perturbation must thermalize KE on merge"
    );
    let dke = ke_before - ke_after;
    assert!(
        ((logged - dke) / dke).abs() < 1e-12,
        "logged ΔKE {logged} must equal the resolved-KE difference {dke}"
    );
    assert!(((mass_before - mass_after) / mass_before).abs() < 1e-14);
}

// --- §6-9: guard band + re-expansion machinery -------------------------------

#[test]
fn fnd2_s6_9_indicator_sees_m1_at_guard_and_controller_fires_after_dwell() {
    let mut g = Grid::build(spec(8, 8, 8), FIELDS).unwrap();
    let mom = momentum(&g);
    let rho_id = g.field_id("rho").unwrap();
    let mt_id = g.field_id("mom_theta").unwrap();
    g.fill_field(rho_id, |_, _, _| 2.0);
    g.fill_field(mt_id, |_, _, _| 0.0);
    for bi in 0..g.n_bricks() {
        g.coarsen_theta(bi, Some(mom)).expect("to guard");
        assert_eq!(g.brick(bi).n_theta(), N_THETA_GUARD);
        // Adaptive collapse below the guard is refused (S4).
        assert!(matches!(
            g.coarsen_theta(bi, Some(mom)),
            Err(GridError::BadThetaResolution { .. })
        ));
    }

    // Symmetric state: indicator ~ 0 at guard resolution.
    let fields = [(mt_id, 1e-6), (rho_id, 2e-6)];
    assert!(g.symmetry_indicator(0, &fields).unwrap() < 1e-30);

    // A seeded m=1 perturbation at N_θ^guard = 4 is visible (both phases
    // carried — the S4 rationale) and can drive the indicator past τ_expand.
    let nt = g.brick(0).n_theta();
    let (data_start, cells) = (0usize, BRICK_CELLS);
    let mt = g.brick_field_mut(0, mt_id);
    for j in 0..nt {
        let theta = std::f64::consts::TAU * (f64::from(j) + 0.5) / f64::from(nt);
        for local in 0..cells {
            mt[data_start + j as usize * BRICK_CELLS + local] = 0.05 * theta.cos();
        }
    }
    let a = g.symmetry_indicator(0, &fields).unwrap();
    assert!(a > TAU_EXPAND, "m=1 at guard must be detectable: A = {a}");

    // Dwell: N_DWELL consecutive over-threshold observations before Expand.
    let mut ctl = ThetaController::default();
    assert_eq!(ctl.observe(a), ThetaAction::Hold);
    assert_eq!(ctl.observe(a), ThetaAction::Hold);
    assert_eq!(ctl.observe(a), ThetaAction::Hold);
    assert_eq!(ctl.observe(a), ThetaAction::Expand);
    // Hysteresis: mid-band observations reset the dwell counters.
    let mut ctl2 = ThetaController::default();
    ctl2.observe(a);
    ctl2.observe(1e-5); // between τ_collapse and τ_expand
    assert_eq!(ctl2.observe(a), ThetaAction::Hold);

    g.refine_theta(0).expect("expand 4 → 8 after the trigger");
    assert_eq!(g.brick(0).n_theta(), 8);
    assert!(
        matches!(g.refine_theta(0), Err(GridError::BadThetaResolution { .. })),
        "cannot refine past N_θ^max"
    );
}

// --- θ-ladder & assertion rules (§3.4 / FND-4 §3.4-5b) -----------------------

#[test]
fn fnd2_s34_ladder_and_axisymmetry_assertion_rules() {
    // Off-ladder N_θ^max refused at construction.
    assert!(matches!(
        Grid::build(spec(4, 4, 24), FIELDS),
        Err(GridError::BadThetaResolution { .. })
    ));
    // N_θ = 1 without the assertion refused.
    assert!(matches!(
        Grid::build(spec(4, 4, 1), FIELDS),
        Err(GridError::BadThetaResolution { .. })
    ));

    // With the recorded assertion: build at N_θ^max, then collapse to 1.
    let mut s = spec(4, 4, 8);
    s.axisymmetry_assertion = true;
    let mut g = Grid::build(s, FIELDS).unwrap();
    fill_deterministic(&mut g);
    let rho_e = g.field_id("rho_e").unwrap();
    let before = g.reduce_volume_weighted(rho_e);
    let mom = momentum(&g);
    let log = g
        .assert_axisymmetric(0, Some(mom))
        .expect("asserted collapse");
    assert!(log.by_assertion);
    assert_eq!(g.brick(0).n_theta(), 1);
    let after = g.reduce_volume_weighted(rho_e);
    assert!(
        ((before - after) / before).abs() < 1e-13,
        "assertion collapse conserves"
    );
    // An asserted region never re-expands adaptively.
    assert!(g.refine_theta(0).is_err());

    // Without the assertion, assert_axisymmetric is refused.
    let mut g2 = Grid::build(spec(4, 4, 8), FIELDS).unwrap();
    assert!(g2.assert_axisymmetric(0, None).is_err());
}

// --- Post-review fail-loud regressions ---------------------------------------

#[test]
fn review_nan_and_inf_spec_values_are_refused_at_build() {
    for bad in [f64::NAN, f64::INFINITY] {
        let mut s = spec(4, 4, 8);
        s.dr = bad;
        assert!(matches!(Grid::build(s, FIELDS), Err(GridError::BadSpec(_))));
        let mut s = spec(4, 4, 8);
        s.dz = bad;
        assert!(matches!(Grid::build(s, FIELDS), Err(GridError::BadSpec(_))));
    }
}

#[test]
fn review_symmetry_indicator_fails_loud_on_nan_and_bad_floor() {
    let mut g = Grid::build(spec(8, 8, 8), FIELDS).unwrap();
    let rho_id = g.field_id("rho").unwrap();
    g.fill_field(rho_id, |_, _, _| 2.0);

    // Zero/negative/NaN floors are refused (a zero floor makes an all-zero
    // field divide 0/0 — the exact case f64::max used to swallow).
    for bad_floor in [0.0, -1.0, f64::NAN] {
        assert!(g.symmetry_indicator(0, &[(rho_id, bad_floor)]).is_err());
    }

    // NaN in a monitored field is a halt, never "perfectly symmetric".
    g.brick_field_mut(0, rho_id)[3] = f64::NAN;
    match g.symmetry_indicator(0, &[(rho_id, 1e-6)]) {
        Err(GridError::NonFinite { .. }) => {}
        other => panic!("expected NonFinite halt, got {other:?}"),
    }

    // And the collapse ΔKE path refuses NaN density rather than merging it.
    let mom = momentum(&g);
    assert!(matches!(
        g.coarsen_theta(0, Some(mom)),
        Err(GridError::NonPositiveDensity { .. })
    ));
}
