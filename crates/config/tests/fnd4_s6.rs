//! FND-4 v0.2 §6 validation plan, test-first (VAL-3 §3.1). Item numbering
//! follows the doc. **Item 4 (table-pin content-hash mismatch) is deferred
//! to the FND-5 loader session** — until then the loader refuses non-empty
//! `[tables]`, which is asserted here instead.
//!
//! The registry is a test double: COUP-8-shaped `Manifest`s exercising the
//! dispatch, O20/O21, and coupler-balance paths without any real solver.

use crucible_config::{DEFAULT_DETERMINISM_MODE, DEFAULT_RNG_ALGORITHM, Diagnostics, load_str};
use crucible_registry::{
    ChaoticClass, CouplerKind, CouplerSpec, Direction, InterfaceVersion, Manifest, ParamSpec,
    ParamType, ParamValue, PortKind, PortRole, PortSpec, Quantity, Regime, Registry,
};

const IV: InterfaceVersion = InterfaceVersion { major: 1, minor: 0 };
const NON_CHAOTIC: &[(Regime, ChaoticClass)] = &[
    (Regime::Steady, ChaoticClass::NonChaotic),
    (Regime::Transient, ChaoticClass::NonChaotic),
];

/// Diffusion-like mechanism: one required + one defaulted param.
static HEAT: Manifest = Manifest {
    id: "heat-diffusion",
    interface_version: IV,
    tables: &[],
    couplers: &[],
    ports: &[],
    chaotic_class: NON_CHAOTIC,
    params: &[
        ParamSpec {
            name: "kappa_w_per_m_k",
            ty: ParamType::Float,
            default: None, // required
            range: Some((0.0, 1.0e5)),
        },
        ParamSpec {
            name: "limiter",
            ty: ParamType::Str,
            default: Some(ParamValue::Str("minmod")),
            range: None,
        },
    ],
};

/// Turbulence-resolving mechanism: chaotic in the Transient regime (O21 prey).
static TURB: Manifest = Manifest {
    id: "turb-combustion",
    interface_version: IV,
    tables: &[],
    couplers: &[CouplerSpec {
        kind: CouplerKind::WallExchange,
        quantity: Quantity::Heat,
        direction: Direction::Source,
    }],
    ports: &[],
    chaotic_class: &[
        (Regime::Steady, ChaoticClass::NonChaotic),
        (Regime::Transient, ChaoticClass::Chaotic),
    ],
    params: &[],
};

/// Wall counterpart so TURB's coupler edge balances.
static WALL: Manifest = Manifest {
    id: "wall-conduction",
    interface_version: IV,
    tables: &[],
    couplers: &[CouplerSpec {
        kind: CouplerKind::WallExchange,
        quantity: Quantity::Heat,
        direction: Direction::Sink,
    }],
    ports: &[],
    chaotic_class: NON_CHAOTIC,
    params: &[],
};

/// Consumer with a Require port of kind Pump (O20 prey).
static FEED: Manifest = Manifest {
    id: "feed-consumer",
    interface_version: IV,
    tables: &[],
    couplers: &[],
    ports: &[PortSpec {
        name: "inflow",
        role: PortRole::Require,
        kind: PortKind::Pump,
    }],
    chaotic_class: NON_CHAOTIC,
    params: &[],
};

/// Pump provider (two instances of this = two same-kind candidates).
static PUMP: Manifest = Manifest {
    id: "pump",
    interface_version: IV,
    tables: &[],
    couplers: &[],
    ports: &[PortSpec {
        name: "outflow",
        role: PortRole::Provide,
        kind: PortKind::Pump,
    }],
    chaotic_class: NON_CHAOTIC,
    params: &[],
};

// Sorted by id (COUP-8 §3.2).
static MECHS: [&Manifest; 5] = [&FEED, &HEAT, &PUMP, &TURB, &WALL];

fn registry() -> Registry {
    Registry::new(&MECHS, &[])
}

fn messages(d: &Diagnostics) -> Vec<String> {
    d.iter()
        .map(|x| format!("[{}] {}", x.path, x.message))
        .collect()
}

// --- §6-1: config↔manifest round-trip, zero implicit values -----------------

#[test]
fn fnd4_s6_1_resolved_config_is_a_fixed_point_with_zero_implicit_values() {
    // Omits: limiter (defaulted), [determinism], rng.algorithm, [meta].
    let author = r#"
        schema_version = 1
        [mechanisms.core]
        type = "heat-diffusion"
        kappa_w_per_m_k = 45.0
        [rng]
        master_seed = 7
    "#;
    let loaded = load_str(author, &registry()).expect("valid config must load");

    // Every effective default materialized.
    let core = &loaded.resolved.mechanisms["core"];
    assert_eq!(core["limiter"].as_str(), Some("minmod"));
    assert_eq!(loaded.resolved.determinism.mode, DEFAULT_DETERMINISM_MODE);
    assert_eq!(loaded.resolved.rng.algorithm, DEFAULT_RNG_ALGORITHM);
    assert_eq!(loaded.resolved.rng.master_seed, 7);

    // Fixed point: reloading the serialized resolved config reproduces it
    // exactly — nothing was left implicit.
    let round = load_str(&loaded.resolved.to_toml(), &registry())
        .expect("resolved config must itself be a valid author config");
    assert_eq!(round.resolved, loaded.resolved);
}

// --- §6-2: unknown key / unknown mechanism fail loud ------------------------

#[test]
fn fnd4_s6_2a_unknown_top_level_key_is_a_hard_error() {
    let err = load_str("schema_version = 1\nturbo_mode = true\n", &registry()).unwrap_err();
    assert!(
        format!("{err}").contains("turbo_mode"),
        "error must name the unknown key: {err}"
    );
}

#[test]
fn fnd4_s6_2b_unknown_mechanism_lists_valid_ids() {
    let author = r#"
        schema_version = 1
        [mechanisms.x]
        type = "warp-drive"
    "#;
    let err = load_str(author, &registry()).unwrap_err();
    let all = messages(&err).join("\n");
    assert!(all.contains("warp-drive"), "{all}");
    for id in ["feed-consumer", "heat-diffusion", "pump", "turb-combustion"] {
        assert!(all.contains(id), "valid ids must be listed: {all}");
    }
}

#[test]
fn fnd4_s6_2c_unknown_and_missing_params_fail_loud() {
    let author = r#"
        schema_version = 1
        [mechanisms.core]
        type = "heat-diffusion"
        kappa_typo = 45.0
    "#;
    let err = load_str(author, &registry()).unwrap_err();
    let all = messages(&err).join("\n");
    assert!(all.contains("kappa_typo"), "unknown param named: {all}");
    assert!(
        all.contains("kappa_w_per_m_k") && all.contains("required"),
        "missing required param named: {all}"
    );
}

// --- §6-3: all errors surfaced together, deterministic order ----------------

#[test]
fn fnd4_s6_3_multi_fault_config_reports_every_fault_in_deterministic_order() {
    // Four independent faults: unknown mechanism, bad θ-ladder, relaxed+chaotic,
    // out-of-range param.
    let author = r#"
        schema_version = 1
        [geometry]
        n_theta_max = 24
        [determinism]
        mode = "relaxed"
        [mechanisms.bad]
        type = "no-such-thing"
        [mechanisms.core]
        type = "heat-diffusion"
        kappa_w_per_m_k = -3.0
        [mechanisms.flame]
        type = "turb-combustion"
        [mechanisms.wall]
        type = "wall-conduction"
    "#;
    let err = load_str(author, &registry()).unwrap_err();
    assert!(err.len() >= 4, "all faults reported at once, got: {err}");
    let all = messages(&err).join("\n");
    assert!(all.contains("no-such-thing"), "{all}");
    assert!(all.contains("n_theta_max"), "{all}");
    assert!(all.contains("relaxed"), "{all}");
    assert!(all.contains("outside validity range"), "{all}");

    // Deterministic: two loads give byte-identical reports.
    let err2 = load_str(author, &registry()).unwrap_err();
    assert_eq!(format!("{err}"), format!("{err2}"));
}

// --- §6-4 (deferred to FND-5): [tables] currently refuses -------------------

#[test]
fn fnd4_s6_4_placeholder_nonempty_tables_block_refuses() {
    let author = r#"
        schema_version = 1
        [tables]
        eos = "eos@3"
    "#;
    let err = load_str(author, &registry()).unwrap_err();
    assert!(format!("{err}").contains("FND-5"), "{err}");
}

// --- §6-5: schema versioning ------------------------------------------------

#[test]
fn fnd4_s6_5_version_gates_refuse_loudly() {
    let newer = load_str("schema_version = 2\n", &registry()).unwrap_err();
    assert!(format!("{newer}").contains("newer"), "{newer}");
    let older = load_str("schema_version = 0\n", &registry()).unwrap_err();
    assert!(format!("{older}").contains("older"), "{older}");
    let missing = load_str("[meta]\nname = \"x\"\n", &registry()).unwrap_err();
    assert!(format!("{missing}").contains("schema_version"), "{missing}");
}

// --- §6-6: single source of intent (no env influence) -----------------------

// `unsafe` justification (META-2 §4): `std::env::set_var` is unsafe in
// edition 2024 (thread-safety of the process environment). This test is the
// one place the project *deliberately* mutates the environment — to prove the
// loader ignores it. Test-only, no physics path, single-threaded use.
#[allow(unsafe_code)]
#[test]
fn fnd4_s6_6_environment_variables_do_not_change_the_result() {
    let author = r#"
        schema_version = 1
        [mechanisms.core]
        type = "heat-diffusion"
        kappa_w_per_m_k = 45.0
        [rng]
        master_seed = 42
    "#;
    let before = load_str(author, &registry()).expect("loads");
    // SAFETY-adjacent note: set_var is fine in a single-threaded test binary
    // section; the point is the loader must not read it.
    unsafe { std::env::set_var("CRUCIBLE_KAPPA_OVERRIDE", "9999") };
    let after = load_str(author, &registry()).expect("loads");
    unsafe { std::env::remove_var("CRUCIBLE_KAPPA_OVERRIDE") };
    assert_eq!(before.resolved, after.resolved);
    assert_eq!(
        before.manifest.config_content_hash,
        after.manifest.config_content_hash
    );
}

// --- §6-7: O21 determinism refusal ------------------------------------------

#[test]
fn fnd4_s6_7_relaxed_plus_chaotic_is_refused_naming_the_mechanism() {
    let base = r#"
        schema_version = 1
        [determinism]
        mode = "MODE"
        [mechanisms.flame]
        type = "turb-combustion"
        [mechanisms.wall]
        type = "wall-conduction"
    "#;
    let relaxed = base.replace("MODE", "relaxed");
    let err = load_str(&relaxed, &registry()).unwrap_err();
    let all = format!("{err}");
    assert!(
        all.contains("flame") && all.contains("turb-combustion"),
        "refusal must name the chaotic mechanism: {all}"
    );
    assert!(
        !all.contains("wall-conduction ("),
        "non-chaotic not blamed: {all}"
    );

    let fixed = base.replace("MODE", "fixed-order");
    let loaded = load_str(&fixed, &registry()).expect("fixed-order must load");
    assert_eq!(loaded.manifest.determinism_mode, "fixed-order");
    // The manifest records the classification in force (§3.6).
    let flame = loaded
        .manifest
        .chaotic_class_in_force
        .iter()
        .find(|r| r.instance == "flame")
        .expect("flame recorded");
    assert_eq!(flame.classes["Transient"], "chaotic");
}

#[test]
fn fnd4_s6_7b_relaxed_with_only_nonchaotic_mechanisms_loads() {
    let author = r#"
        schema_version = 1
        [determinism]
        mode = "relaxed"
        [mechanisms.core]
        type = "heat-diffusion"
        kappa_w_per_m_k = 45.0
    "#;
    let loaded = load_str(author, &registry()).expect("relaxed+non-chaotic is admissible");
    assert_eq!(loaded.manifest.determinism_mode, "relaxed");
}

// --- §6-8: θ-ladder alignment ------------------------------------------------

#[test]
fn fnd4_s6_8_theta_ladder_24_refused_naming_16_and_32() {
    let author = "schema_version = 1\n[geometry]\nn_theta_max = 24\n";
    let err = load_str(author, &registry()).unwrap_err();
    let all = format!("{err}");
    assert!(all.contains("16") && all.contains("32"), "{all}");

    for good in [4, 8, 32, 256] {
        let cfg = format!("schema_version = 1\n[geometry]\nn_theta_max = {good}\n");
        let loaded = load_str(&cfg, &registry()).expect("aligned N_theta_max loads");
        assert_eq!(loaded.resolved.geometry.as_ref().unwrap().n_theta_max, good);
    }
}

// --- §6-9: O20 explicit bindings ---------------------------------------------

#[test]
fn fnd4_s6_9_unbound_require_with_two_candidates_lists_both() {
    let author = r#"
        schema_version = 1
        [mechanisms.chamber]
        type = "feed-consumer"
        [mechanisms.pump_a]
        type = "pump"
        [mechanisms.pump_b]
        type = "pump"
    "#;
    let err = load_str(author, &registry()).unwrap_err();
    let all = format!("{err}");
    assert!(
        all.contains("pump_a.outflow") && all.contains("pump_b.outflow"),
        "both same-kind candidates must be listed: {all}"
    );

    // With an explicit binding it loads; nothing auto-matched.
    let bound = format!("{author}\n[engine.bindings]\n\"chamber.inflow\" = \"pump_a.outflow\"\n");
    let loaded = load_str(&bound, &registry()).expect("explicitly bound config loads");
    assert_eq!(
        loaded.resolved.engine.bindings["chamber.inflow"],
        "pump_a.outflow"
    );
}

#[test]
fn fnd4_s6_9b_dangling_binding_and_kind_mismatch_fail_loud() {
    let author = r#"
        schema_version = 1
        [mechanisms.chamber]
        type = "feed-consumer"
        [engine.bindings]
        "chamber.inflow" = "ghost.outflow"
        "chamber.typo" = "chamber.inflow"
    "#;
    let err = load_str(author, &registry()).unwrap_err();
    let all = format!("{err}");
    assert!(all.contains("ghost"), "dangling target named: {all}");
    assert!(
        all.contains("chamber.typo"),
        "binding key naming no Require port refused: {all}"
    );
}

// --- Coupler balance (COUP-8 §3.3-3, exercised through the loader) ----------

#[test]
fn coup8_s33_3_one_sided_coupler_edge_fails_loud() {
    let author = r#"
        schema_version = 1
        [mechanisms.flame]
        type = "turb-combustion"
    "#;
    let err = load_str(author, &registry()).unwrap_err();
    let all = format!("{err}");
    assert!(
        all.contains("WallExchange/Heat") && all.contains("no Sink"),
        "{all}"
    );
}

// --- §3.8: canonical hashing -------------------------------------------------

#[test]
fn fnd4_s38_semantically_identical_configs_hash_identically() {
    let a = "schema_version = 1\n[rng]\nmaster_seed = 1\nalgorithm = \"philox4x32-10\"\n";
    let b = "schema_version = 1\n[rng]\nalgorithm = \"philox4x32-10\"\nmaster_seed = 1\n";
    let ha = load_str(a, &registry())
        .unwrap()
        .manifest
        .config_content_hash;
    let hb = load_str(b, &registry())
        .unwrap()
        .manifest
        .config_content_hash;
    assert_eq!(ha, hb, "key order must not affect the content hash");

    let c = "schema_version = 1\n[rng]\nmaster_seed = 2\nalgorithm = \"philox4x32-10\"\n";
    let hc = load_str(c, &registry())
        .unwrap()
        .manifest
        .config_content_hash;
    assert_ne!(ha, hc, "different intent must hash differently");
}

// --- Auto-drawn seed is recorded, not wall-clock -----------------------------

#[test]
fn fnd4_s31_omitted_seed_is_auto_drawn_and_recorded() {
    let author = "schema_version = 1\n";
    let loaded = load_str(author, &registry()).unwrap();
    assert!(loaded.manifest.master_seed >= 0);
    assert_eq!(loaded.manifest.master_seed, loaded.resolved.rng.master_seed);
    // Replaying the *resolved* config reproduces the drawn seed exactly.
    let replay = load_str(&loaded.resolved.to_toml(), &registry()).unwrap();
    assert_eq!(
        replay.resolved.rng.master_seed,
        loaded.resolved.rng.master_seed
    );
}
