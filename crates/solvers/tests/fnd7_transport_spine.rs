//! FND-7 §3.3 / OFFL-5 §3.1a — the spine's transport slot on the
//! **production** artifact `tables/spine/lox_lh2_transport_v0.2.0.h5`
//! (regenerable via `offline/scripts/make_spine_transport_tables.py`),
//! loaded through the full FND-5 §3.2 pin gate.
//!
//! What this file is for: the tabulated occupant is where "real properties"
//! stops being a claim and becomes a number. The tests below check the
//! things that would let a wrong number through quietly —
//!
//! - the **units gate** at bind (a producer-side relabel must refuse, not
//!   misread);
//! - the **derived-group identities** the module promises (`Pr = μc_p/k`
//!   both directions; `ρD = μ/Sc`);
//! - the **equilibrium-conductivity composition** actually firing, and in
//!   the right direction and magnitude — this is the piece with no analytic
//!   fixture behind it, so it is checked against the physics it encodes;
//! - **envelope refusal** rather than extrapolation;
//! - and the property that motivated the whole session: the tabulated
//!   occupant must genuinely **differ from the constant one** across the
//!   engine's state range. A spine that returned chamber-like numbers
//!   everywhere would pass every other test in this repo.
//!
//! Pins are READ from the machine-written sidecar — one owner, no
//! hand-copied digests (the station-3 review finding).

use crucible_solvers::transport::{
    COLUMNS, ConstantTransport, MediumState, TabulatedTransport, TransportSpine,
};
use crucible_tables::{Pin, Table};
use crucible_units::{dynamic_viscosity_pa_s, specific_heat_capacity_j_per_kg_k};

const FILE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../tables/spine/lox_lh2_transport_v0.2.0.h5"
);
const GROUP: &str = "/spine/lox_lh2/transport";
const PINS_TOML: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../tables/spine/lox_lh2_transport_v0.2.0.pins.toml"
));
/// The Schmidt number the RL10 configs declare (META-3
/// `schmidt-combustion-gas`).
const SC: f64 = 0.5;

fn table() -> Table {
    let doc: toml::Table = PINS_TOML.parse().expect("pins sidecar parses as TOML");
    let entry = doc
        .get(GROUP)
        .and_then(|v| v.as_table())
        .expect("pins sidecar has the transport group");
    let pin = Pin {
        data_version: entry["data_version"].as_str().expect("string").to_string(),
        content_digest: Some(
            entry["content_digest"]
                .as_str()
                .expect("string")
                .to_string(),
        ),
    };
    Table::open(FILE, GROUP, &pin).expect("production transport surface opens under its pin")
}

/// Chamber-class state on the RL10 design line (p ≈ 3 MPa, injection
/// enthalpy, MR 5) and a deep plume-fringe state — the two ends of what a
/// march actually visits.
fn chamber() -> MediumState {
    MediumState {
        p: 3.0e6,
        h: -1.083e6,
        z: 1.0 / 6.0,
    }
}

fn fringe() -> MediumState {
    MediumState {
        p: 1.0e3,
        h: -1.0e7,
        z: 1.0 / 6.0,
    }
}

#[test]
fn binds_the_declared_columns_and_refuses_a_units_drift() {
    let t = table();
    let spine = TransportSpine::Tabulated(
        TabulatedTransport::bind(&t, SC).expect("the six declared columns bind"),
    );
    assert!(spine.is_tabulated());
    // Every declared column must exist under exactly its declared units —
    // the string is the cross-language contract with
    // `offline/crucible_offl/transport.py`'s `_TR_COLUMNS`.
    for (name, units) in COLUMNS {
        t.bind(name, units).unwrap_or_else(|e| {
            panic!("column {name} must bind as {units}: {e}");
        });
        assert!(
            t.bind(name, "totally-not-the-units").is_err(),
            "{name} must refuse a units drift at bind, not misread at kernel rate"
        );
    }
}

#[test]
fn refuses_a_schmidt_that_is_not_a_number() {
    let t = table();
    assert!(TabulatedTransport::bind(&t, 0.0).is_err());
    assert!(TabulatedTransport::bind(&t, f64::NAN).is_err());
    assert!(TabulatedTransport::bind(&t, -0.5).is_err());
}

#[test]
fn derived_groups_hold_at_every_probed_state() {
    let t = table();
    let spine = TransportSpine::Tabulated(TabulatedTransport::bind(&t, SC).expect("bind"));
    for st in [chamber(), fringe()] {
        let p = spine.eval(st).expect("in-envelope state");
        assert!(p.mu > 0.0 && p.k > 0.0 && p.cp > 0.0 && p.cv > 0.0);
        // ρD = μ/Sc, exactly.
        assert_eq!(p.rho_d, p.mu / SC);
        // Pr = μc_p/k, to round-off (the tabulated occupant derives it in
        // exactly this direction; the constant occupant derives k from a
        // declared Pr and the same identity is asserted there).
        let pr = p.mu * p.cp / p.k;
        assert!((pr - p.pr).abs() <= 8.0 * f64::EPSILON * p.pr);
        // c_p > c_v for any gas.
        assert!(p.cp > p.cv, "c_p must exceed c_v");
    }
}

/// The effective-conductivity composition (`crate::transport` module doc)
/// has no analytic fixture behind it, so check the physics it encodes:
/// where the gas is barely dissociated the equilibrium and molecular
/// conductivities nearly coincide, and where it is strongly dissociated the
/// recombination limb dominates — by a large, one-signed factor.
#[test]
fn equilibrium_conductivity_fires_where_dissociation_runs_and_not_elsewhere() {
    let t = table();
    let spine = TransportSpine::Tabulated(TabulatedTransport::bind(&t, SC).expect("bind"));
    let k_frozen = t.bind("conductivity_frozen", "W/(m*K)").expect("bind");
    let ratio = |st: MediumState| {
        let p = spine.eval(st).expect("in-envelope");
        let k_fr = k_frozen
            .interpolate(&[st.p, st.h, st.z])
            .expect("in-envelope");
        (p.k / k_fr, p.pr)
    };
    // Cold fringe: recombination is done, so k_eff ≈ k_frozen.
    let (cold, pr_cold) = ratio(fringe());
    assert!(
        (1.0..1.15).contains(&cold),
        "a cold, recombined gas has no reaction conductivity to add (got {cold:.3}×)"
    );
    // Hot, low-pressure: dissociation is strong and the recombination limb
    // must dominate the molecular one.
    let (hot, pr_hot) = ratio(MediumState {
        p: 1.0e3,
        h: 0.0,
        z: 1.0 / 6.0,
    });
    assert!(
        hot > 3.0,
        "a strongly dissociated gas must conduct far more than its molecular \
         conductivity alone (got {hot:.3}×) — this is the classical \
         equilibrium-conductivity peak, and the wall-heat flux depends on it"
    );
    // And the whole point of composing k and c_p consistently: the Prandtl
    // number stays a gas Prandtl number at both ends. If k had been left
    // frozen while c_p went equilibrium, Pr would blow up with the
    // dissociation and the Colburn analogy would silently misfire.
    for pr in [pr_cold, pr_hot] {
        assert!(
            (0.2..1.5).contains(&pr),
            "Pr = {pr:.3} is not a gas Prandtl number — the k/c_p pairing is inconsistent"
        );
    }
}

/// The session's motivating claim, as a test: real properties are not the
/// constants they replace. The RL10 configs declared c_p = 5000, μ = 1e-4,
/// Pr = 0.6 — good at the chamber, and that is the point; a spine that
/// merely reproduced them everywhere would have been wasted work.
#[test]
fn tabulated_transport_departs_from_the_retired_constants_where_it_should() {
    let t = table();
    let tab = TransportSpine::Tabulated(TabulatedTransport::bind(&t, SC).expect("bind"));
    let konst = TransportSpine::Constant(
        ConstantTransport::new(
            specific_heat_capacity_j_per_kg_k(5000.0),
            dynamic_viscosity_pa_s(1.0e-4),
            0.6,
            1.2,
            SC,
        )
        .expect("the retired RL10 constants"),
    );
    let k_ref = konst.eval(chamber()).expect("state-independent");
    let at_chamber = tab.eval(chamber()).expect("in-envelope");
    let at_fringe = tab.eval(fringe()).expect("in-envelope");

    // At the chamber the declared constants were a defensible fit — within
    // a factor of ~2 of the truth on every quantity.
    for (name, got, want) in [
        ("mu", at_chamber.mu, k_ref.mu),
        ("cp", at_chamber.cp, k_ref.cp),
        ("k", at_chamber.k, k_ref.k),
    ] {
        let r = got / want;
        assert!(
            (0.4..2.5).contains(&r),
            "{name} at the chamber: tabulated {got:.4e} vs the retired constant \
             {want:.4e} ({r:.2}×) — the constants were sized for THIS state, so a \
             large gap here means the surface, not the constants, is wrong"
        );
    }
    // At the plume fringe they were not: the cold, recombined gas is far
    // less viscous and far less conductive than the chamber-sized numbers.
    assert!(
        at_fringe.mu < 0.5 * k_ref.mu,
        "fringe μ {:.3e} should be far below the chamber constant {:.3e}",
        at_fringe.mu,
        k_ref.mu
    );
    assert!(
        at_fringe.k < 0.5 * k_ref.k,
        "fringe k {:.3e} should be far below the chamber constant {:.3e}",
        at_fringe.k,
        k_ref.k
    );
    // And c_v — the class-D temperature solve's linearization slope —
    // varies by a large factor across the same span. This is the deferral
    // the S3 module header recorded: a constant c_v is exact only for the
    // gamma-law class, and off by this much for a real dissociating gas.
    let hot = tab
        .eval(MediumState {
            p: 1.0e3,
            h: 0.0,
            z: 1.0 / 6.0,
        })
        .expect("in-envelope");
    assert!(
        hot.cv / at_fringe.cv > 3.0,
        "c_v spans only {:.2}× between the recombined fringe and the dissociated \
         core — a constant slope would have been fine and S4's T-solve work \
         unnecessary",
        hot.cv / at_fringe.cv
    );
}

#[test]
fn out_of_envelope_refuses_rather_than_extrapolating() {
    let t = table();
    let spine = TransportSpine::Tabulated(TabulatedTransport::bind(&t, SC).expect("bind"));
    let env = match &spine {
        TransportSpine::Tabulated(x) => x.envelopes(),
        TransportSpine::Constant(_) => unreachable!(),
    };
    let good = chamber();
    assert!(spine.eval(good).is_ok());
    for bad in [
        MediumState {
            p: env[0].0 * 0.5,
            ..good
        },
        MediumState {
            p: env[0].1 * 2.0,
            ..good
        },
        MediumState {
            h: env[1].0 - 1.0e6,
            ..good
        },
        MediumState {
            h: env[1].1 + 1.0e6,
            ..good
        },
        MediumState {
            z: env[2].0 - 0.05,
            ..good
        },
        MediumState {
            z: env[2].1 + 0.05,
            ..good
        },
        MediumState {
            p: f64::NAN,
            ..good
        },
    ] {
        assert!(
            spine.eval(bad).is_err(),
            "off-surface query {bad:?} must refuse, never extrapolate (FND-5 §3.5)"
        );
    }
}

/// The measured interpolation bounds must be present and comfortably below
/// the declared 10–20% physics band (FND-5 §3.4's grid-sizing rule). This
/// is the Rust half of the Python gate — it reads the shipped metadata the
/// runtime actually loads.
#[test]
fn interpolation_bounds_are_declared_and_sit_inside_the_physics_band() {
    let t = table();
    let bound = TabulatedTransport::bind(&t, SC)
        .expect("bind")
        .interp_error_bounds();
    assert!(bound.iter().all(|b| b.is_finite() && *b > 0.0));
    for (i, (name, units)) in COLUMNS.iter().enumerate() {
        let c = t.bind(name, units).expect("bind");
        if let Some(rel) = c.interp_error_bound_log() {
            assert!(
                rel < 0.10,
                "{name}: rule-space bound {:.3} rivals the declared 10% band",
                rel
            );
        }
        assert_eq!(bound[i], c.interp_error_bound());
    }
}
