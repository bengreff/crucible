//! The sandbox seam end-to-end on the REAL geometry-of-record and the REAL
//! pinned table: load a config-driven engine (contour CSV + TOML + pins),
//! assemble the ternary grid, march a short budget, and read out SOLV-7.
//! A tiny dial + short settle keep this in the fast battery (~seconds);
//! steadiness and anchor comparison belong to the run presets and the
//! cycle-wave certificate, not here.

use crucible_config::load_str_with_sidecars;

const ROOT: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../..");

fn read_rel(path: &str) -> Result<String, String> {
    std::fs::read_to_string(format!("{ROOT}/{path}")).map_err(|e| format!("{path}: {e}"))
}

/// A coarse smoke preset: the RL10 contour at dial 2.5, cold-start march of
/// 2 flow-throughs — enough for the nozzle to choke, not enough to settle.
fn smoke_author() -> String {
    r#"
        schema_version = 1
        [meta]
        name = "engine-smoke"
        [geometry]
        n_theta_max = 1
        axisymmetric = true
        contour = "data/anchors/rl10_contour.csv"
        contour_units = "in"
        liner_thickness_m = 0.05
        cells_across_throat = 2.5
        [mechanisms.flow]
        type = "flow_shifting"
        [mechanisms.injector]
        type = "injector_prior"
        mdot_kg_per_s = 16.9465
        mixture_ratio = 5.0
        h_inj_j_per_kg = -1083091.499726406
        [mechanisms.liner]
        type = "conduction"
        kappa_w_per_m_k = 2000.0
        rho_cp_j_per_m3_k = 1000.0
        [mechanisms.wall]
        type = "wall_heat"
        [mechanisms.transport]
        type = "transport_constant"
        cp_j_per_kg_k = 5000.0
        mu_pa_s = 1.0e-4
        pr = 0.6
        gamma = 1.2
        schmidt = 0.5
        [mechanisms.jacket]
        type = "jacket_coolant"
        h_w_per_m2_k = 30000.0
        t_coolant_k = 120.0
        [engine.bindings]
        "flow.inflow" = "injector.combustion_inlet"
        [operating_profile]
        mode = "steady_march"
        flowthroughs = 2.0
        fill_p_pa = 3.0e6
        cfl = 0.4
        [tables.chem_equilibrium]
        file = "tables/chem/lox_lh2_v0.3.2.h5"
        group = "/chem/lox_lh2/equilibrium"
        pins = "tables/chem/lox_lh2_v0.3.2.pins.toml"
        [determinism]
        mode = "fixed-order"
        [rng]
        master_seed = 11
    "#
    .to_string()
}

#[test]
fn config_driven_engine_assembles_marches_and_reads_out() {
    // Table paths inside the config are CWD-relative; tests run from the
    // crate dir, so pin the process CWD to the repo root for the duration.
    std::env::set_current_dir(ROOT).expect("repo root exists");
    let registry = crucible_engine::registry();
    let loaded =
        load_str_with_sidecars(&smoke_author(), &registry, &read_rel).expect("smoke config loads");

    let mut spec = crucible_engine::assembly::assemble(&loaded, &read_rel).expect("assembles");
    assert!(
        spec.jacket.is_some() && spec.wall_law.is_some(),
        "cooled build"
    );
    // The geometry-of-record: throat radius 2.47 in (CSV-header erratum).
    assert!(
        (spec.contour.r_throat_m - 2.47 * crucible_config::INCH_M).abs() < 1e-12,
        "throat radius of record"
    );

    let table = crucible_engine::run::open_pinned_table(&spec).expect("pinned table opens");
    let report = crucible_engine::run::run(&mut spec, &table, None, &mut |_| {})
        .expect("short march completes");

    // Machinery truths, not steadiness claims (2 flow-throughs only):
    assert!(report.mdot_exit_kg_per_s > 0.0, "gas leaves the nozzle");
    assert!(
        report.thrust_n > 0.0,
        "positive exit momentum+pressure flux"
    );
    assert!(
        report.p_c_pa > 1.0e6,
        "chamber holds MPa-class pressure (emergent): {}",
        report.p_c_pa
    );
    assert!(report.jacket_watts > 0.0, "hot gas heats the liner");
    assert!(
        report.liner_t_max_k > 120.0 && report.liner_t_max_k < 4000.0,
        "liner temperature sane: {}",
        report.liner_t_max_k
    );
    assert!(report.active_gas_cells > 300, "contour interior is gas");
    assert!(!report.fields_csv.is_empty(), "viz feed emitted");
}

/// The **S4 configuration end to end**: the same seam with the two things
/// plan S4 added — the tabulated FND-7 spine (real μ/k/c_p/c_v/∂h∂Z over
/// the local state) and SOLV-1 §3.1's `F_visc` scheduled. This is ◆C1 in
/// miniature and it is a gate, not a demo: what it proves is that the
/// engine composes them **from config alone**, that the two pinned surfaces
/// agree on their envelope, and that a march with the missing forces on
/// still chokes the nozzle and heats the liner — with no schedule tuning
/// separating it from the constant-transport preset above.
fn s4_author() -> String {
    smoke_author()
        .replace(
            r#"[mechanisms.transport]
        type = "transport_constant"
        cp_j_per_kg_k = 5000.0
        mu_pa_s = 1.0e-4
        pr = 0.6
        gamma = 1.2
        schmidt = 0.5"#,
            r#"[mechanisms.transport]
        type = "transport_table"
        schmidt = 0.5
        [mechanisms.gasvisc]
        type = "gas_diffusion""#,
        )
        .replace(
            r#"pins = "tables/chem/lox_lh2_v0.3.2.pins.toml""#,
            r#"pins = "tables/chem/lox_lh2_v0.3.2.pins.toml"
        [tables.spine_transport]
        file = "tables/spine/lox_lh2_transport_v0.1.0.h5"
        group = "/spine/lox_lh2/transport"
        pins = "tables/spine/lox_lh2_transport_v0.1.0.pins.toml""#,
        )
}

#[test]
fn s4_config_composes_the_spine_and_the_missing_forces() {
    std::env::set_current_dir(ROOT).expect("repo root exists");
    let registry = crucible_engine::registry();
    let loaded =
        load_str_with_sidecars(&s4_author(), &registry, &read_rel).expect("S4 config loads");
    let mut spec = crucible_engine::assembly::assemble(&loaded, &read_rel).expect("assembles");
    assert!(spec.gas_diffusion, "F_visc is scheduled");
    assert!(
        matches!(
            spec.transport,
            crucible_engine::assembly::TransportSpec::Table { .. }
        ),
        "the tabulated spine occupant is selected"
    );

    let table = crucible_engine::run::open_pinned_table(&spec).expect("EOS table opens");
    let transport = crucible_engine::run::open_transport_table(&spec)
        .expect("spine transport table opens")
        .expect("the tabulated occupant needs one");
    let report = crucible_engine::run::run(&mut spec, &table, Some(&transport), &mut |_| {})
        .expect("short march completes with the missing forces on");

    assert!(report.mdot_exit_kg_per_s > 0.0, "gas leaves the nozzle");
    assert!(
        report.thrust_n > 0.0,
        "positive exit momentum+pressure flux"
    );
    assert!(
        report.p_c_pa > 1.0e6,
        "chamber holds MPa-class pressure (emergent): {}",
        report.p_c_pa
    );
    assert!(report.jacket_watts > 0.0, "hot gas heats the liner");
}

/// The envelope-identity guard fires before any march (OFFL-5 §3.1a). A
/// transport surface that did not cover the EOS surface would let a march
/// walk off one while the other still answered — so the engine refuses at
/// assembly. Proven by pointing the transport pin at the EOS table's own
/// group, whose axes are the same but whose columns are not.
#[test]
fn a_transport_surface_that_is_not_the_spine_surface_refuses_at_bind() {
    std::env::set_current_dir(ROOT).expect("repo root exists");
    let registry = crucible_engine::registry();
    let author = s4_author().replace(
        r#"file = "tables/spine/lox_lh2_transport_v0.1.0.h5"
        group = "/spine/lox_lh2/transport"
        pins = "tables/spine/lox_lh2_transport_v0.1.0.pins.toml""#,
        r#"file = "tables/chem/lox_lh2_v0.3.2.h5"
        group = "/chem/lox_lh2/equilibrium"
        pins = "tables/chem/lox_lh2_v0.3.2.pins.toml""#,
    );
    let loaded = load_str_with_sidecars(&author, &registry, &read_rel).expect("loads");
    let mut spec = crucible_engine::assembly::assemble(&loaded, &read_rel).expect("assembles");
    let table = crucible_engine::run::open_pinned_table(&spec).expect("EOS table opens");
    let wrong = crucible_engine::run::open_transport_table(&spec)
        .expect("the equilibrium table opens under its own pin")
        .expect("present");
    let halt = crucible_engine::run::run(&mut spec, &table, Some(&wrong), &mut |_| {})
        .expect_err("a surface without the transport columns must refuse");
    assert!(
        halt.message.contains("spine transport bind"),
        "the refusal must name the bind, not fail somewhere downstream: {}",
        halt.message
    );
}
