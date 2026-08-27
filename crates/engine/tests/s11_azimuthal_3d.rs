//! S11 / ◆C3: the engine assembles and marches a genuinely 3-D (N_θ > 1)
//! world. The flow operator's azimuthal machinery (axis crossing, per-sector
//! apertures, reflux) was built and gated in the SOLVER at S8/S9; this proves
//! it now reaches the ENGINE — config → GridSpec with n_theta_max > 1 →
//! build_with_geometry_theta on the revolved contour → the SDC step with the
//! axis machinery live → θ-summed SOLV-7 readout — and that the deferred
//! couplings (cooled wall, F_visc on cut θ-faces) refuse loudly at assembly.

use crucible_config::load_str_with_sidecars;

const ROOT: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../..");

fn read_rel(path: &str) -> Result<String, String> {
    std::fs::read_to_string(format!("{ROOT}/{path}")).map_err(|e| format!("{path}: {e}"))
}

/// A coarse ADIABATIC 3-D preset: the RL10 contour at a tiny dial, N_θ = 8,
/// no cooled wall (liner_thickness_m = 0) — the ◆C3 flow core. Short march.
fn author_3d(extra: &str) -> String {
    format!(
        r#"
        schema_version = 1
        [meta]
        name = "s11-3d-smoke"
        [geometry]
        n_theta_max = 8
        axisymmetric = false
        contour = "data/anchors/rl10_contour.csv"
        contour_units = "in"
        liner_thickness_m = 0.0
        cells_across_throat = 2.0
        [mechanisms.flow]
        type = "flow_shifting"
        [mechanisms.injector]
        type = "injector_prior"
        mdot_kg_per_s = 16.9465
        mixture_ratio = 5.0
        h_inj_j_per_kg = -1083091.499726406
        [mechanisms.transport]
        type = "transport_constant"
        cp_j_per_kg_k = 5000.0
        mu_pa_s = 1.0e-4
        pr = 0.6
        gamma = 1.2
        schmidt = 0.5
        {extra}
        [engine.bindings]
        "flow.inflow" = "injector.combustion_inlet"
        [operating_profile]
        mode = "steady_march"
        flowthroughs = 0.4
        fill_p_pa = 3.0e6
        cfl = 0.4
        [tables.chem_equilibrium]
        file = "tables/chem/lox_lh2_v0.4.0.h5"
        group = "/chem/lox_lh2/equilibrium"
        pins = "tables/chem/lox_lh2_v0.4.0.pins.toml"
        [determinism]
        mode = "fixed-order"
        [rng]
        master_seed = 11
    "#
    )
}

#[test]
fn coarse_3d_rl10_assembles_marches_audited_and_reads_out_theta_summed() {
    std::env::set_current_dir(ROOT).expect("repo root exists");
    let registry = crucible_engine::registry();
    let loaded =
        load_str_with_sidecars(&author_3d(""), &registry, &read_rel).expect("3-D config loads");
    let mut spec = crucible_engine::assembly::assemble(&loaded, &read_rel).expect("assembles");

    // THE 3-D dimension reached the engine grid.
    assert_eq!(spec.grid.spec().n_theta_max, 8, "N_θ = 8 world");
    assert!(
        !spec.grid.spec().axisymmetry_assertion,
        "no axisym assertion"
    );
    assert!(
        spec.jacket.is_none() && spec.wall_law.is_none(),
        "adiabatic"
    );

    let table = crucible_engine::run::open_pinned_table(&spec).expect("pinned table opens");
    // A few flow-throughs on the real SDC step with the axis/reflux machinery
    // live; the COUP-2 audit is armed EVERY step, so a clean return IS the
    // audit-closes gate.
    let report = crucible_engine::run::run(&mut spec, &table, None, None, &mut |_| {})
        .expect("3-D march completes, audit closes every step");

    assert!(
        report.active_gas_cells > 200,
        "contour interior is gas (3-D)"
    );
    assert!(
        report.thrust_n > 0.0,
        "θ-summed exit momentum+pressure flux > 0"
    );
    assert!(report.mdot_exit_kg_per_s > 0.0, "gas leaves the nozzle");
    assert!(report.p_c_pa > 0.0, "θ-averaged chamber pressure > 0");
    // The 3-D viz feed carries the θ column + swirl velocity.
    assert!(
        report
            .fields_csv
            .starts_with("r_m,theta_rad,z_m,region,rho,u_r,u_theta,u_z,"),
        "3-D fields CSV must carry the θ column: {:?}",
        report.fields_csv.lines().next()
    );
}

#[test]
fn n_theta_gt_1_with_a_cooled_wall_refuses_the_s11_deferral() {
    std::env::set_current_dir(ROOT).expect("repo root exists");
    let registry = crucible_engine::registry();
    // A liner + jacket + wall + conduction on a 3-D world — the per-θ wall
    // patch / coupled conduction wave (typed S11 refusal).
    let cooled = r#"
        [mechanisms.liner]
        type = "conduction"
        kappa_w_per_m_k = 2000.0
        rho_cp_j_per_m3_k = 1000.0
        [mechanisms.wall]
        type = "wall_heat"
        [mechanisms.jacket]
        type = "jacket_coolant"
        h_w_per_m2_k = 30000.0
        t_coolant_k = 120.0"#;
    // liner_thickness_m must be > 0 for a cooled build, so patch it in too.
    let author = author_3d(cooled).replace("liner_thickness_m = 0.0", "liner_thickness_m = 0.02");
    let loaded = load_str_with_sidecars(&author, &registry, &read_rel).expect("loads");
    let Err(err) = crucible_engine::assembly::assemble(&loaded, &read_rel) else {
        panic!("cooled 3-D must refuse at assembly");
    };
    assert!(
        err.contains("N_θ > 1") && err.contains("cooled"),
        "the refusal must name the cooled-wall S11 deferral: {err}"
    );
}

#[test]
fn n_theta_gt_1_with_gas_diffusion_refuses_the_s11_deferral() {
    std::env::set_current_dir(ROOT).expect("repo root exists");
    let registry = crucible_engine::registry();
    let author = author_3d(
        r#"[mechanisms.gasvisc]
        type = "gas_diffusion""#,
    );
    let loaded = load_str_with_sidecars(&author, &registry, &read_rel).expect("loads");
    let Err(err) = crucible_engine::assembly::assemble(&loaded, &read_rel) else {
        panic!("F_visc on cut θ > 1 must refuse at assembly");
    };
    assert!(
        err.contains("N_θ > 1") && err.contains("gas_diffusion"),
        "the refusal must name the F_visc S11 deferral: {err}"
    );
}
