//! The Rust half of the cross-language seam proof (FND-5 §3.2, VISION_SCOPE
//! §6): `tests/fixtures/python_seam.h5` was written by the *Python* side
//! (`offline/crucible_offl/seam_fixture.py`, regenerable via
//! `offline/scripts/make_seam_fixture.py`) and is opened here with the same
//! pin the Python tests assert. If h5py and this reader ever disagree about
//! string encodings, attribute layout, dataset bytes, or the digest stream,
//! this test halts the battery — a table the Rust loader can't verify is
//! worthless, so this is the seam's load-bearing test.

use crucible_tables::{Pin, Table, TableError};

const FIXTURE: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/python_seam.h5");
const GROUP: &str = "/chem/seam_probe";
// Must equal crucible_offl.seam_fixture.{DATA_VERSION, PIN_DIGEST}.
const DATA_VERSION: &str = "0.1.0";
const PIN_DIGEST: &str = "sha256:5f91b7ab6cbd34c836fa214175553de5f18cee3c6842837443e0f4ace49a78b4";

fn pin() -> Pin {
    Pin {
        data_version: DATA_VERSION.into(),
        content_digest: Some(PIN_DIGEST.into()),
    }
}

#[test]
fn python_written_table_opens_and_digest_verifies() {
    let t = Table::open(FIXTURE, GROUP, &pin()).expect("python-written table must load");
    assert_eq!(t.content_digest, PIN_DIGEST);

    // The features the fixture deliberately exercises, read back exactly.
    assert_eq!(t.kind, "regular");
    assert_eq!(t.provenance.producer, "crucible-offl/seam-fixture");
    assert_eq!(
        t.provenance.source_library, "Glenn-α (UTF-8 probe)",
        "UTF-8 length-prefixing must count bytes, not chars"
    );
    assert_eq!(t.provenance.rng_seed, Some(20260817));
    assert_eq!(t.axes.len(), 2);
    assert_eq!(t.axes[0].name, "p");
    assert_eq!(t.axes[0].envelope_min, 2.0e4);
    assert_eq!(t.values.len(), 3);
    let temp = t
        .values
        .iter()
        .find(|v| v.name == "temperature")
        .expect("temperature");
    assert_eq!(temp.units, "K");
    assert!(temp.sigma.is_some(), "per-point sigma dataset must survive");
    let gamma = t
        .values
        .iter()
        .find(|v| v.name == "gamma_eff")
        .expect("gamma_eff");
    assert_eq!(gamma.sigma_scalar, Some(0.005));

    // And it interpolates: a grid point returns its tabulated value exactly
    // (multilinear at a node), through the log-p rule.
    let v = t
        .interpolate("temperature", &[1.0e5, 1.0e6])
        .expect("in-envelope query");
    assert_eq!(v, 300.0 + 10.0 * 1.0 + 1.5 * 1.0);
}

#[test]
fn python_written_table_wrong_digest_refuses() {
    let bad = Pin {
        data_version: DATA_VERSION.into(),
        content_digest: Some(
            "sha256:0000000000000000000000000000000000000000000000000000000000000000".into(),
        ),
    };
    match Table::open(FIXTURE, GROUP, &bad) {
        Err(TableError::PinDigestMismatch { found, .. }) => assert_eq!(found, PIN_DIGEST),
        other => panic!("expected PinDigestMismatch, got {other:?}"),
    }
}
