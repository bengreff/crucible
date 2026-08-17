//! Goal-B **Station 1** pass criteria (SOLV-1 §6-1; VAL-2 §3.3 `sod-shock`).
//! The committed artifact in `certificates/` is rendered by
//! `bin/station1_sod_certificate` from the same study runners AND the same
//! named criteria constants asserted here; gate 4 of `scripts/check.sh`
//! regenerates it and diffs, so the record cannot silently diverge.

use crucible_solvers::station1_sod::{
    ADV_C_ORDER_MAX, ADV_C_ORDER_MIN, ADV_MEAN_ORDER_MIN, ADV_RHO_STEP_ORDER_MIN, C_BOUNDS_TOL,
    CONSERVATION_TOL_REL, CONTACT_POS_TOL_CELLS, FAN_ORDER_MIN, FIELDS, SHOCK_POS_TOL_CELLS,
    SOD_FINEST_L1_MAX, SOD_L1_DECREASE_MIN, SOD_LEVELS, STAR_PLATEAU_TOL, STAR_WINDOW_ORDER_MIN,
    advection_order, closed_tube_conservation, run_sod, sod_convergence, sod_waves,
    uniform_state_is_bitwise_fixed_point,
};
use crucible_solvers::{SetupError, flow_from_loaded, registry};

// --- The waves against the exact Riemann oracle ------------------------------

#[test]
fn station1_sod_convergence_ladder() {
    let study = sod_convergence();
    assert_eq!(study.levels.len(), SOD_LEVELS.len());

    // Global L1(ρ): monotone decrease at the discontinuity-limited rate,
    // absolute accuracy at the finest level.
    for w in study.levels.windows(2) {
        assert!(
            w[0].l1_rho > SOD_L1_DECREASE_MIN * w[1].l1_rho,
            "global L1 must fall ≥{SOD_L1_DECREASE_MIN}× per refinement: {study:?}"
        );
    }
    let finest = &study.levels[study.levels.len() - 1];
    assert!(
        finest.l1_rho < SOD_FINEST_L1_MAX,
        "finest global L1 {} exceeds {SOD_FINEST_L1_MAX}",
        finest.l1_rho
    );

    // Formal order in smooth regions: the star-left plateau window.
    for (i, p) in study.observed_orders(|l| l.l1_star_l).iter().enumerate() {
        assert!(
            *p >= STAR_WINDOW_ORDER_MIN,
            "star-window order {p} at refinement {i} below {STAR_WINDOW_ORDER_MIN}: {study:?}"
        );
    }
    // Fan interior: first-order (startup singularity), gated against
    // silent further degradation.
    for (i, p) in study.observed_orders(|l| l.l1_fan).iter().enumerate() {
        assert!(
            *p >= FAN_ORDER_MIN,
            "fan-interior order {p} at refinement {i} below {FAN_ORDER_MIN}: {study:?}"
        );
    }
}

#[test]
fn station1_sod_wave_structure() {
    let w = sod_waves();
    assert!(
        w.shock_pos_err_cells <= SHOCK_POS_TOL_CELLS,
        "shock position off by {} cells (VAL-2: within {SHOCK_POS_TOL_CELLS})",
        w.shock_pos_err_cells
    );
    assert!(
        w.contact_pos_err_cells <= CONTACT_POS_TOL_CELLS,
        "contact position off by {} cells",
        w.contact_pos_err_cells
    );
    for (name, err) in [
        ("ρ*L", w.rho_star_l_err),
        ("ρ*R", w.rho_star_r_err),
        ("u*", w.u_star_err),
        ("p*", w.p_star_err),
    ] {
        assert!(
            err <= STAR_PLATEAU_TOL,
            "star plateau {name} off by {err} (tol {STAR_PLATEAU_TOL})"
        );
    }
    // The composition rides the contact and stays in [0, 1] to round-off —
    // HLLC-Batten's species-advection guarantee.
    assert!(w.c_min >= -C_BOUNDS_TOL, "C undershoot: {}", w.c_min);
    assert!(w.c_max <= 1.0 + C_BOUNDS_TOL, "C overshoot: {}", w.c_max);
    // One law on the one grid: a radially-uniform tube stays radially
    // uniform BITWISE through the whole shock evolution.
    assert!(w.radially_uniform_bitwise, "radial uniformity broke");
}

// --- Formal order on the smooth advection ladder -----------------------------

#[test]
fn station1_smooth_advection_order() {
    let adv = advection_order();
    let mean_rho = adv.mean_order(|l| l.l1_rho);
    assert!(
        mean_rho >= ADV_MEAN_ORDER_MIN,
        "mean ρ order {mean_rho} below {ADV_MEAN_ORDER_MIN}: {adv:?}"
    );
    for (i, p) in adv.observed_orders(|l| l.l1_rho).iter().enumerate() {
        assert!(
            *p >= ADV_RHO_STEP_ORDER_MIN,
            "ρ order {p} at refinement {i} below the {ADV_RHO_STEP_ORDER_MIN} floor: {adv:?}"
        );
    }
    for (i, p) in adv.observed_orders(|l| l.l1_c).iter().enumerate() {
        assert!(
            (ADV_C_ORDER_MIN..=ADV_C_ORDER_MAX).contains(p),
            "composition order {p} at refinement {i} outside \
             [{ADV_C_ORDER_MIN}, {ADV_C_ORDER_MAX}]: {adv:?}"
        );
    }
}

// --- Well-balance, conservation, determinism ---------------------------------

#[test]
fn station1_uniform_state_is_a_bitwise_fixed_point() {
    // The well-balanced geometric sources against the exact cylindrical
    // metric: a uniform gas at rest does not move, bit for bit — with the
    // r = 0 axis in the domain (N_θ = 1) and on an annulus at N_θ = 8
    // (exercising the θ sweep).
    assert!(uniform_state_is_bitwise_fixed_point(true), "axis config");
    assert!(
        uniform_state_is_bitwise_fixed_point(false),
        "annulus config"
    );
}

#[test]
fn station1_closed_tube_conserves_mass_and_energy() {
    let (dm, de) = closed_tube_conservation().expect("run");
    assert!(
        dm < CONSERVATION_TOL_REL,
        "mass drift {dm} (tol {CONSERVATION_TOL_REL})"
    );
    assert!(
        de < CONSERVATION_TOL_REL,
        "energy drift {de} (tol {CONSERVATION_TOL_REL})"
    );
}

#[test]
fn station1_sod_rerun_is_byte_identical() {
    let bits = |n: usize| -> Vec<u64> {
        let (g, f, _) = run_sod(n);
        g.bricks()
            .iter()
            .flat_map(|b| {
                f.ids()
                    .into_iter()
                    .flat_map(|id| b.field(id).iter().map(|v| v.to_bits()).collect::<Vec<_>>())
            })
            .collect()
    };
    assert_eq!(bits(200), bits(200), "rerun must be byte-identical (S6)");
}

// --- End-to-end: config → loader → registry → grid → operator ---------------

const E2E_CONFIG: &str = r#"
schema_version = 1

[meta]
name = "station1-sod"
description = "bursting diaphragm: the station-1 certificate config"

[geometry]
n_theta_max = 1
axisymmetric = true
r_min = 0.0
dr = 0.025
n_r = 4
z_min = 0.0
dz = 0.01
n_z = 100

[mechanisms.tube]
type = "flow"
gamma = 1.4

[determinism]
mode = "fixed-order"

[rng]
master_seed = 1
"#;

#[test]
fn station1_config_drives_a_flow_run() {
    let loaded = crucible_config::load_str(E2E_CONFIG, &registry()).expect("config loads");
    let setup = flow_from_loaded(&loaded, FIELDS).expect("setup");
    assert_eq!(setup.instance, "tube");
    assert_eq!(setup.gamma, 1.4);
    assert_eq!(setup.grid.spec().n_z, 100);
}

#[test]
fn station1_two_flow_instances_are_refused_not_silently_picked() {
    let two = format!("{E2E_CONFIG}\n[mechanisms.tube2]\ntype = \"flow\"\ngamma = 1.4\n");
    let loaded = crucible_config::load_str(&two, &registry()).expect("loader accepts both");
    match flow_from_loaded(&loaded, FIELDS) {
        Err(SetupError::Ambiguous { instances, .. }) => {
            assert_eq!(instances, ["tube", "tube2"]);
        }
        other => panic!("expected Ambiguous refusal, got {other:?}"),
    }
}

#[test]
fn station1_gamma_outside_validity_range_is_refused_at_load() {
    // META-1 P8: the γ validity range in the flow Manifest is machinery,
    // not discipline — an isothermal-limit γ refuses at the front door.
    let bad = E2E_CONFIG.replace("gamma = 1.4", "gamma = 1.0");
    assert!(
        crucible_config::load_str(&bad, &registry()).is_err(),
        "γ = 1.0 must refuse at load (validity range [1.001, 1.667])"
    );
}
