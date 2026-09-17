//! FND-5 §3.3 — `BoundColumn` contract on the *production* equilibrium
//! surface: a bound query must be **bit-identical** to the general
//! `Table::interpolate` path (same fixed-order kernel, the parse/allocation
//! merely hoisted to bind time), and every refusal class must survive the
//! hoisting unchanged (units gate at bind; envelope/domain/arity per query).

use crucible_tables::{Pin, Table, TableError};

const FILE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../tables/chem/lox_lh2_v0.1.0.h5"
);
const PINS_TOML: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../tables/chem/lox_lh2_v0.1.0.pins.toml"
));

fn open_equilibrium() -> Table {
    let doc: toml::Table = PINS_TOML.parse().expect("pins sidecar parses as TOML");
    let entry = doc["/chem/lox_lh2/equilibrium"]
        .as_table()
        .expect("sidecar entry");
    let pin = Pin {
        data_version: entry["data_version"].as_str().expect("string").to_string(),
        content_digest: Some(
            entry["content_digest"]
                .as_str()
                .expect("string")
                .to_string(),
        ),
    };
    Table::open(FILE, "/chem/lox_lh2/equilibrium", &pin)
        .expect("production equilibrium surface must load under its pin")
}

/// An irrationally-spaced sweep of the declared envelope (never on grid
/// nodes), covering the log axis (p), the lin axes (h, Z), a log-valued
/// column (`density`) and lin-valued columns.
fn sweep_points() -> Vec<[f64; 3]> {
    // Envelope: p ∈ [2e3, 7e6], h ∈ [-11.5e6, -2e5], Z ∈ [1/9, 1/4].
    let mut pts = Vec::new();
    for ip in 0..7 {
        // log-spaced with an irrational offset so no query lands on a node
        let fp = (ip as f64 + std::f64::consts::FRAC_1_PI) / 7.0;
        let p = 2.0e3 * (7.0e6f64 / 2.0e3).powf(fp);
        for ih in 0..5 {
            let fh = (ih as f64 + 0.577_215_664) / 5.0; // Euler–Mascheroni
            let h = -11.5e6 + fh * (11.5e6 - 2.0e5);
            for iz in 0..5 {
                let fz = (iz as f64 + std::f64::consts::FRAC_1_SQRT_2) / 5.0;
                let z = 1.0 / 9.0 + fz * (0.25 - 1.0 / 9.0);
                pts.push([p, h, z]);
            }
        }
    }
    pts
}

#[test]
fn bound_query_is_bit_identical_to_the_general_path() {
    let t = open_equilibrium();
    let columns = [
        ("temperature", "K"),
        ("density", "kg/m^3"), // log-valued: exercises the ln-hoisted data
        ("gamma_eff", "1"),
        ("sound_speed", "m/s"),
    ];
    let mut checked = 0usize;
    for (name, units) in columns {
        let b = t.bind(name, units).expect("bind must succeed");
        assert_eq!(b.name(), name);
        assert!(b.interp_error_bound() > 0.0);
        for q in sweep_points() {
            let general = t.interpolate(name, &q).expect("general path in-envelope");
            let bound = b.interpolate(&q).expect("bound path in-envelope");
            assert_eq!(
                general.to_bits(),
                bound.to_bits(),
                "{name} at {q:?}: general {general} vs bound {bound} — the bound kernel \
                 must be the identical fixed-order arithmetic"
            );
            checked += 1;
        }
    }
    assert_eq!(checked, 4 * 7 * 5 * 5, "sweep coverage");
}

#[test]
fn bind_runs_the_units_gate() {
    let t = open_equilibrium();
    match t.bind("temperature", "degC") {
        Err(TableError::UnitsMismatch {
            value, expected, ..
        }) => {
            assert_eq!(value, "temperature");
            assert_eq!(expected, "degC");
        }
        other => panic!("relabel bind must refuse: {other:?}"),
    }
    match t.bind("no_such_column", "K") {
        Err(TableError::UnknownValue { name }) => assert_eq!(name, "no_such_column"),
        other => panic!("unknown column bind must refuse: {other:?}"),
    }
}

#[test]
fn bound_refusals_match_the_general_path() {
    let t = open_equilibrium();
    let b = t.bind("temperature", "K").expect("bind");
    // In grid domain but outside the declared envelope (p grid reaches
    // 1.5e3 but the envelope floor is 2e3) — Refuse policy hard-errors.
    let q_env = [1.6e3, -1.0e6, 1.0 / 6.0];
    match (t.interpolate("temperature", &q_env), b.interpolate(&q_env)) {
        (
            Err(TableError::OutOfEnvelope { axis: a1, .. }),
            Err(TableError::OutOfEnvelope { axis: a2, .. }),
        ) => {
            assert_eq!(a1, a2);
        }
        other => panic!("both paths must refuse out-of-envelope identically: {other:?}"),
    }
    // Outside the grid domain entirely.
    let q_dom = [9.0e6, -1.0e6, 1.0 / 6.0];
    match b.interpolate(&q_dom) {
        Err(TableError::OutOfDomain { axis, .. }) => assert_eq!(axis, "p"),
        other => panic!("out-of-domain must refuse: {other:?}"),
    }
    // NaN refuses (never a silent poison value).
    match b.interpolate(&[f64::NAN, -1.0e6, 1.0 / 6.0]) {
        Err(TableError::OutOfDomain { .. }) => {}
        other => panic!("NaN query must refuse: {other:?}"),
    }
    // Arity.
    match b.interpolate(&[1.0e6, -1.0e6]) {
        Err(TableError::QueryArity { expected, found }) => {
            assert_eq!((expected, found), (3, 2));
        }
        other => panic!("arity must refuse: {other:?}"),
    }
}

#[test]
fn axis_domain_reports_the_grid_span() {
    let t = open_equilibrium();
    let b = t.bind("temperature", "K").expect("bind");
    let (p_lo, p_hi) = b.axis_domain(0);
    assert!(p_lo < 2.0e3 && p_hi >= 7.0e6, "p grid spans the envelope");
}
