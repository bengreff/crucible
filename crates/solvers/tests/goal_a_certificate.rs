//! The Goal-A convergence-certificate pass criteria (VAL-1 rung "analytic";
//! VAL-3 §3.2 MMS + envelope + determinism gates). The committed artifact in
//! `certificates/` is rendered by `bin/convergence_certificate` from the
//! same study runners AND the same named criteria constants asserted here
//! (gate 4 of `scripts/check.sh` regenerates the artifact and diffs it, so
//! the committed record can never silently diverge from the code).

use crucible_solvers::certificate::{
    self, ANNULUS_TOL_REL, BESSEL_TOL_K, CFL_SAFETY, CONSERVATION_STEPS, CONSERVATION_TOL_REL,
    FIELDS, MMS_ORDER_MAX, MMS_ORDER_MIN, annulus_anchor, bessel_cylinder_anchor,
    conservation_drift, mms_axisymmetric, mms_theta_mode,
};
use crucible_solvers::{Bcs, Conduction, FaceBc, SetupError, from_loaded, registry};

// --- MMS order-of-accuracy (the heart of the certificate) --------------------

#[test]
fn cert_mms_axisymmetric_observed_order_is_two() {
    let study = mms_axisymmetric();
    let orders = study.observed_orders();
    assert_eq!(orders.len(), 2);
    for (i, p) in orders.iter().enumerate() {
        assert!(
            (MMS_ORDER_MIN..=MMS_ORDER_MAX).contains(p),
            "observed order {p} at refinement {i} outside \
             [{MMS_ORDER_MIN}, {MMS_ORDER_MAX}]: {study:?}"
        );
    }
    // Errors must actually decrease ~4× per level (2^MMS_ORDER_MIN ≈ 3.5;
    // 3× is that bound with rounding headroom — order alone could be met
    // by two tiny errors at the noise floor).
    assert!(study.levels[0].l2_error > 3.0 * study.levels[1].l2_error);
    assert!(study.levels[1].l2_error > 3.0 * study.levels[2].l2_error);
}

#[test]
fn cert_mms_theta_mode_observed_order_is_two() {
    let study = mms_theta_mode();
    for (i, p) in study.observed_orders().iter().enumerate() {
        assert!(
            (MMS_ORDER_MIN..=MMS_ORDER_MAX).contains(p),
            "θ-mode observed order {p} at refinement {i} outside band: {study:?}"
        );
    }
}

// --- Analytic anchors --------------------------------------------------------

#[test]
fn cert_annulus_log_profile_anchor() {
    let (worst_rel, _g) = annulus_anchor();
    assert!(
        worst_rel < ANNULUS_TOL_REL,
        "steady annulus profile off by {worst_rel} (rel to ΔT = 200 K)"
    );
}

#[test]
fn cert_bessel_cylinder_anchor_exercises_the_axis() {
    let worst_abs = bessel_cylinder_anchor();
    assert!(
        worst_abs < BESSEL_TOL_K,
        "transient cylinder vs Bessel series: max error {worst_abs} K (T0 = 100 K)"
    );
}

// --- Conservation (FND-2 §6-2, closed-sweep half) ----------------------------

#[test]
fn cert_closed_sweep_conserves_energy() {
    let drift = conservation_drift(CONSERVATION_STEPS);
    assert!(
        drift < CONSERVATION_TOL_REL,
        "closed insulated sweep drifted {drift} relative"
    );
}

// --- End-to-end: config → loader → registry → grid → operator ---------------

const E2E_CONFIG: &str = r#"
schema_version = 1

[meta]
name = "goal-a-certificate"
description = "conduction end-to-end: the convergence-certificate config"

[geometry]
n_theta_max = 1
axisymmetric = true
r_min = 0.05
dr = 0.003125
n_r = 32
z_min = 0.0
dz = 0.05
n_z = 2

[mechanisms.wall]
type = "conduction"
kappa_w_per_m_k = 20.0
rho_cp_j_per_m3_k = 4.0e6

[determinism]
mode = "fixed-order"

[rng]
master_seed = 1
"#;

#[test]
fn cert_config_drives_a_run_and_reruns_bit_identically() {
    let run = |steps: usize| -> (Vec<u64>, String) {
        let loaded = crucible_config::load_str(E2E_CONFIG, &registry()).expect("config loads");
        let setup = from_loaded(&loaded, FIELDS).expect("setup");
        let mut g = setup.grid;
        assert_eq!(setup.instance, "wall");
        // The setup boundary is unit-typed (META-2 §4 ★); the kernel
        // receives documented-SI f64 via the one extraction point.
        assert_eq!(crucible_units::si(setup.kappa), 20.0);
        let t_id = g.field_id("T").unwrap();
        let rate_id = g.field_id("rate").unwrap();
        g.fill_field(t_id, |_, _, _| 400.0);
        let hot = |_: f64, _: f64, _: f64, _: f64| 500.0;
        let cold = |_: f64, _: f64, _: f64, _: f64| 300.0;
        let zero_src = |_: f64, _: f64, _: f64, _: f64| 0.0;
        let op = Conduction {
            kappa: crucible_units::si(setup.kappa),
            rho_cp: crucible_units::si(setup.rho_cp),
            source: &zero_src,
            domain: crucible_solvers::Domain::FlowActive,
            interior: crucible_solvers::InteriorFaces::refuse(),
            bcs: Bcs {
                r_inner: FaceBc::Dirichlet(&hot),
                r_outer: FaceBc::Dirichlet(&cold),
                z_lo: FaceBc::HeatFlux(0.0),
                z_hi: FaceBc::HeatFlux(0.0),
            },
        };
        let dt = op.stable_dt(&g, CFL_SAFETY);
        op.advance(&mut g, t_id, rate_id, 0.0, dt, steps)
            .expect("advance");
        let bits: Vec<u64> = g
            .bricks()
            .iter()
            .flat_map(|b| b.field(t_id).iter().map(|v| v.to_bits()))
            .collect();
        (bits, loaded.manifest.config_content_hash)
    };

    let (bits_a, hash_a) = run(200);
    let (bits_b, hash_b) = run(200);
    assert_eq!(bits_a, bits_b, "rerun must be byte-identical (S6)");
    assert_eq!(hash_a, hash_b);

    // The same config without the assertion is refused at the front door
    // (regression tie to FND-4 §6-8).
    let bad = E2E_CONFIG.replace("axisymmetric = true", "axisymmetric = false");
    assert!(crucible_config::load_str(&bad, &registry()).is_err());
}

#[test]
fn cert_resolved_config_with_extents_is_a_fixed_point() {
    // Review regression: the resolved form used to serialize extents as a
    // nested table the author grammar refused — breaking manifest replay
    // for every config that declares a world.
    let loaded = crucible_config::load_str(E2E_CONFIG, &registry()).expect("loads");
    let replay = crucible_config::load_str(&loaded.resolved.to_toml(), &registry())
        .expect("resolved config with extents must reload (S6 replay)");
    assert_eq!(replay.resolved, loaded.resolved);
}

#[test]
fn cert_two_conduction_instances_are_refused_not_silently_picked() {
    // Review regression: from_loaded used to take the alphabetically-first
    // instance and silently ignore the rest.
    let two = format!(
        "{E2E_CONFIG}\n[mechanisms.liner]\ntype = \"conduction\"\n\
         kappa_w_per_m_k = 400.0\nrho_cp_j_per_m3_k = 3.0e6\n"
    );
    let loaded = crucible_config::load_str(&two, &registry()).expect("loader accepts both");
    match from_loaded(&loaded, FIELDS) {
        Err(SetupError::Ambiguous { instances, .. }) => {
            assert_eq!(instances, ["liner", "wall"], "both named, canonical order");
        }
        other => panic!("expected Ambiguous refusal, got {other:?}"),
    }
}

// --- The certificate studies are themselves deterministic --------------------

#[test]
fn cert_studies_are_reproducible() {
    let a = certificate::mms_axisymmetric();
    let b = certificate::mms_axisymmetric();
    for (x, y) in a.levels.iter().zip(&b.levels) {
        assert_eq!(x.l2_error.to_bits(), y.l2_error.to_bits());
    }
}
