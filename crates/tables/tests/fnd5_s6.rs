//! FND-5 v0.3 §6 validation plan, test-first (VAL-3 §3.1). Item numbering
//! follows the doc. Item 7 (thermo-potential Maxwell consistency) is
//! deferred with the `thermo_potential` kind itself — asserted here as a
//! loud load refusal instead; the full test lands with FND-7/OFFL-5.

use crucible_tables::writer::{WriteSpec, write_table};
use crucible_tables::{Axis, Pin, Provenance, Table, TableError, TableValue};
use std::path::PathBuf;

fn scratch(name: &str) -> String {
    let dir = std::env::temp_dir().join(format!("crucible-fnd5-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("scratch dir");
    let p: PathBuf = dir.join(name);
    let _ = std::fs::remove_file(&p);
    p.to_str().expect("utf-8 temp path").to_string()
}

fn provenance() -> Provenance {
    Provenance {
        producer: "crucible-test".into(),
        producer_version: "0.1.0".into(),
        input_deck_hash: "sha256:0000".into(),
        source_library: "analytic".into(),
        generator_commit: "deadbeef".into(),
        rng_seed: None,
    }
}

fn axis(name: &str, points: Vec<f64>) -> Axis {
    let (lo, hi) = (points[0], points[points.len() - 1]);
    Axis {
        name: name.into(),
        points,
        envelope_min: lo,
        envelope_max: hi,
    }
}

fn value(name: &str, data: Vec<f64>, rule: &str, bound: f64) -> TableValue {
    TableValue {
        name: name.into(),
        data,
        units: "K".into(),
        interp_rule: rule.into(),
        interp_error_bound: bound,
        interp_error_bound_log: None,
        sigma: None,
        sigma_scalar: None,
    }
}

fn spec(axes: Vec<Axis>, values: Vec<TableValue>) -> WriteSpec {
    WriteSpec {
        kind: "regular".into(),
        schema_version: "1.0".into(),
        data_version: "1.0.0".into(),
        interp_method: "multilinear".into(),
        provenance: provenance(),
        axes,
        values,
    }
}

fn pin(digest: &str) -> Pin {
    Pin {
        data_version: "1.0.0".into(),
        content_digest: Some(digest.to_string()),
    }
}

/// 2-D test surface: smooth, positive, non-separable.
fn f(x: f64, y: f64) -> f64 {
    (x / 2.0).exp() * (1.0 + 0.5 * y.cos())
}

fn grid2d(nx: usize, ny: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let xs: Vec<f64> = (0..nx).map(|i| i as f64 / (nx - 1) as f64).collect();
    let ys: Vec<f64> = (0..ny)
        .map(|j| 2.0 * (j as f64) / (ny - 1) as f64)
        .collect();
    let mut data = Vec::with_capacity(nx * ny);
    for &x in &xs {
        for &y in &ys {
            data.push(f(x, y));
        }
    }
    (xs, ys, data)
}

// --- §6-1: round-trip --------------------------------------------------------

#[test]
fn fnd5_s6_1_write_load_read_reproduces_values_and_metadata_exactly() {
    let path = scratch("roundtrip.h5");
    let (xs, ys, data) = grid2d(5, 7);
    let mut v = value("temperature", data.clone(), "lin-lin-lin", 1e-2);
    v.sigma = Some(vec![0.1; data.len()]);
    let s = spec(vec![axis("p", xs.clone()), axis("h", ys.clone())], vec![v]);
    let digest = write_table(&path, "/eos/test", &s).expect("write");

    let t = Table::open(&path, "/eos/test", &pin(&digest)).expect("load");
    assert_eq!(t.content_digest, digest);
    assert_eq!(t.kind, "regular");
    assert_eq!(t.data_version, "1.0.0");
    assert_eq!(t.provenance, provenance());
    assert_eq!(t.axes.len(), 2);
    assert_eq!(t.axes[0].points, xs);
    assert_eq!(t.axes[1].points, ys);
    assert_eq!(t.values[0].data, data);
    assert_eq!(t.values[0].units, "K");
    assert_eq!(t.values[0].interp_error_bound, 1e-2);
    assert_eq!(
        t.values[0].sigma.as_deref(),
        Some(&vec![0.1; data.len()][..])
    );

    // Grid nodes reproduce exactly (multilinear is interpolatory).
    let q = t.interpolate("temperature", &[xs[2], ys[3]]).unwrap();
    assert_eq!(q, f(xs[2], ys[3]));
}

// --- §6-2: interpolation accuracy meets the stored bound ---------------------

#[test]
fn fnd5_s6_2_interpolated_values_meet_the_stored_error_bound() {
    // The producer's obligation (§3.4): measure the actual error and write
    // the bound. We do the same here — measure on a dense sample, then
    // assert the loaded table honors it everywhere on a *different* sample.
    let path = scratch("bound.h5");
    let (xs, ys, data) = grid2d(21, 41);

    // Offline-style measurement pass.
    let probe = |t: &Table, n: usize, seed_off: f64| -> f64 {
        let mut worst = 0.0f64;
        for i in 0..n {
            for j in 0..n {
                let x = (i as f64 + seed_off) / n as f64 * 0.98;
                let y = (j as f64 + seed_off) / n as f64 * 1.96;
                let e = (t.interpolate("temperature", &[x, y]).unwrap() - f(x, y)).abs();
                worst = worst.max(e);
            }
        }
        worst
    };
    let s0 = spec(
        vec![axis("p", xs.clone()), axis("h", ys.clone())],
        vec![value("temperature", data.clone(), "lin-lin-lin", f64::MAX)],
    );
    let d0 = write_table(&path, "/pass0", &s0).expect("write");
    let t0 = Table::open(&path, "/pass0", &pin(&d0)).expect("load");
    let measured = probe(&t0, 100, 0.0);
    assert!(measured > 0.0, "probe must see real interpolation error");

    // Ship with the measured bound (+5% producer margin), verify on an
    // offset sample.
    let bound = measured * 1.05;
    let s1 = spec(
        vec![axis("p", xs), axis("h", ys)],
        vec![value("temperature", data, "lin-lin-lin", bound)],
    );
    let d1 = write_table(&path, "/pass1", &s1).expect("write");
    let t1 = Table::open(&path, "/pass1", &pin(&d1)).expect("load");
    assert_eq!(t1.values[0].interp_error_bound, bound);
    assert!(
        probe(&t1, 137, 0.37) <= bound,
        "declared interp_error_bound must hold on unseen queries"
    );
}

// --- §6-3: monotonicity / positivity (no overshoot) --------------------------

#[test]
fn fnd5_s6_3_multilinear_never_overshoots_on_steep_gradients() {
    // Steep near-step data (tanh with width 0.02 on a coarse grid).
    let path = scratch("steep.h5");
    let xs: Vec<f64> = (0..11).map(|i| i as f64 / 10.0).collect();
    let data: Vec<f64> = xs
        .iter()
        .map(|&x| 1.0 + ((x - 0.5) / 0.02).tanh())
        .collect();
    let (dmin, dmax) = (
        data.iter().cloned().fold(f64::INFINITY, f64::min),
        data.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
    );
    let s = spec(vec![axis("x", xs)], vec![value("q", data, "lin-lin", 1.0)]);
    let d = write_table(&path, "/steep", &s).expect("write");
    let t = Table::open(&path, "/steep", &pin(&d)).expect("load");

    let mut prev = f64::NEG_INFINITY;
    for k in 0..=1000 {
        let x = k as f64 / 1000.0;
        let q = t.interpolate("q", &[x]).unwrap();
        assert!(q >= dmin && q <= dmax, "overshoot at x={x}: {q}");
        assert!(q >= prev, "monotone data must interpolate monotonically");
        prev = q;
    }
}

#[test]
fn fnd5_s6_3b_log_log_interpolation_preserves_positivity() {
    // Power-law data spanning 6 decades: log-log is exact for pure power
    // laws, and positivity must hold everywhere.
    let path = scratch("loglog.h5");
    let xs: Vec<f64> = (0..7).map(|i| 10f64.powi(i - 3)).collect();
    let data: Vec<f64> = xs.iter().map(|&x| 2.5 * x.powf(-1.7)).collect();
    let s = spec(
        vec![axis("e", xs)],
        vec![value("xs", data, "log-log", 1e-12)],
    );
    let d = write_table(&path, "/xs", &s).expect("write");
    let t = Table::open(&path, "/xs", &pin(&d)).expect("load");
    for k in 1..300 {
        let x = 1e-3 * (1e6f64).powf(k as f64 / 300.0);
        let q = t.interpolate("xs", &[x]).unwrap();
        assert!(q > 0.0);
        let exact = 2.5 * x.powf(-1.7);
        assert!(
            (q - exact).abs() / exact < 1e-12,
            "log-log must reproduce a power law to round-off: x={x} q={q} exact={exact}"
        );
    }
}

// --- §6-4: envelope refusal ---------------------------------------------------

#[test]
fn fnd5_s6_4_out_of_envelope_refuses_flag_policy_reports() {
    let path = scratch("envelope.h5");
    let xs: Vec<f64> = (0..11).map(|i| i as f64 / 10.0).collect();
    let data: Vec<f64> = xs.iter().map(|&x| 300.0 + 100.0 * x).collect();
    let mut ax = axis("x", xs);
    // Envelope deliberately tighter than the grid.
    ax.envelope_min = 0.2;
    ax.envelope_max = 0.8;
    let s = spec(vec![ax], vec![value("t", data, "lin-lin", 1e-9)]);
    let d = write_table(&path, "/env", &s).expect("write");
    let t = Table::open(&path, "/env", &pin(&d)).expect("load");

    // Inside envelope: fine under both policies.
    assert!(t.interpolate("t", &[0.5]).is_ok());

    // Outside envelope, inside grid: Refuse errors; Flag computes + reports.
    match t.interpolate("t", &[0.1]) {
        Err(TableError::OutOfEnvelope { axis, .. }) => assert_eq!(axis, "x"),
        other => panic!("expected OutOfEnvelope, got {other:?}"),
    }
    let (v, hits) = t
        .interpolate_flagged("t", &[0.1])
        .expect("flag policy computes");
    assert_eq!(v, 310.0);
    assert_eq!(hits.len(), 1);
    assert_eq!(hits[0].axis, "x");

    // Outside the grid: hard error under BOTH policies — never extrapolate.
    assert!(matches!(
        t.interpolate("t", &[1.5]),
        Err(TableError::OutOfDomain { .. })
    ));
    assert!(matches!(
        t.interpolate_flagged("t", &[1.5]),
        Err(TableError::OutOfDomain { .. })
    ));
}

// --- §6-5: version / digest / provenance gates -------------------------------

#[test]
fn fnd5_s6_5_version_and_digest_mismatches_halt() {
    let path = scratch("pins.h5");
    let xs: Vec<f64> = vec![0.0, 1.0];
    let s = spec(
        vec![axis("x", xs)],
        vec![value("t", vec![1.0, 2.0], "lin-lin", 1e-9)],
    );
    let digest = write_table(&path, "/t", &s).expect("write");

    // Wrong pinned data_version.
    let bad_version = Pin {
        data_version: "2.0.0".into(),
        content_digest: Some(digest.clone()),
    };
    assert!(matches!(
        Table::open(&path, "/t", &bad_version),
        Err(TableError::PinVersionMismatch { .. })
    ));

    // Same label, different bytes: digest catches it (the DVC lesson).
    let tampered = Pin {
        data_version: "1.0.0".into(),
        content_digest: Some(
            "sha256:0000000000000000000000000000000000000000000000000000000000000000".into(),
        ),
    };
    assert!(matches!(
        Table::open(&path, "/t", &tampered),
        Err(TableError::PinDigestMismatch { .. })
    ));

    // First-load pin (no digest) loads and reports the digest for pinning.
    let first = Pin {
        data_version: "1.0.0".into(),
        content_digest: None,
    };
    let t = Table::open(&path, "/t", &first).expect("first load");
    assert_eq!(t.content_digest, digest);
}

#[test]
fn fnd5_s6_5b_schema_major_mismatch_and_missing_metadata_halt() {
    let path = scratch("schema.h5");
    let xs: Vec<f64> = vec![0.0, 1.0];
    let mut s = spec(
        vec![axis("x", xs)],
        vec![value("t", vec![1.0, 2.0], "lin-lin", 1e-9)],
    );
    s.schema_version = "2.0".into();
    let digest = write_table(&path, "/t2", &s).expect("write");
    assert!(matches!(
        Table::open(&path, "/t2", &pin(&digest)),
        Err(TableError::SchemaMajorMismatch { .. })
    ));

    // A table missing a mandatory attribute cannot load: strip `units`.
    let path2 = scratch("missing.h5");
    let mut s2 = spec(
        vec![axis("x", vec![0.0, 1.0])],
        vec![value("t", vec![1.0, 2.0], "lin-lin", 1e-9)],
    );
    s2.schema_version = "1.0".into();
    let d2 = write_table(&path2, "/t", &s2).expect("write");
    {
        let f = hdf5::File::append(&path2).expect("reopen");
        f.dataset("/t/values/t")
            .expect("ds")
            .delete_attr("units")
            .expect("strip");
    }
    match Table::open(&path2, "/t", &pin(&d2)) {
        Err(TableError::MissingAttribute { attribute, .. }) => assert_eq!(attribute, "units"),
        other => panic!("expected MissingAttribute, got {other:?}"),
    }
}

// --- §6-6: determinism --------------------------------------------------------

#[test]
fn fnd5_s6_6_identical_results_across_threads_and_repeats() {
    let path = scratch("det.h5");
    let (xs, ys, data) = grid2d(9, 13);
    let s = spec(
        vec![axis("p", xs), axis("h", ys)],
        vec![value("temperature", data, "lin-lin-lin", 1e-1)],
    );
    let d = write_table(&path, "/det", &s).expect("write");
    let t = std::sync::Arc::new(Table::open(&path, "/det", &pin(&d)).expect("load"));

    let queries: Vec<[f64; 2]> = (0..500)
        .map(|k| {
            let u = k as f64 / 500.0;
            [0.97 * u, 1.9 * (1.0 - u * u)]
        })
        .collect();
    let reference: Vec<u64> = queries
        .iter()
        .map(|q| t.interpolate("temperature", q).unwrap().to_bits())
        .collect();

    let handles: Vec<_> = (0..8)
        .map(|_| {
            let t = t.clone();
            let queries = queries.clone();
            std::thread::spawn(move || {
                queries
                    .iter()
                    .map(|q| t.interpolate("temperature", q).unwrap().to_bits())
                    .collect::<Vec<u64>>()
            })
        })
        .collect();
    for h in handles {
        assert_eq!(
            h.join().expect("thread"),
            reference,
            "bit-identical on every thread"
        );
    }
}

// --- §6-7 (deferred): thermo_potential refuses loudly ------------------------

#[test]
fn fnd5_s6_7_placeholder_thermo_potential_kind_refuses_to_load() {
    let path = scratch("thermo.h5");
    let mut s = spec(
        vec![axis("rho", vec![1.0, 2.0])],
        vec![value("free_energy", vec![1.0, 2.0], "lin-lin", 1e-9)],
    );
    s.kind = "thermo_potential".into();
    let d = write_table(&path, "/F", &s).expect("write");
    match Table::open(&path, "/F", &pin(&d)) {
        Err(TableError::NotYetImplemented { what, .. }) => {
            assert!(what.contains("thermo_potential"), "{what}")
        }
        other => panic!("expected loud refusal, got {other:?}"),
    }
}

// --- Digest is the cross-language contract: golden vector --------------------

#[test]
fn fnd5_s32_digest_golden_vector_guards_the_python_contract() {
    // If this test breaks, the digest encoding changed — which breaks every
    // Python-side stamp. Change requires bumping the digest tag string AND
    // the OFFL implementations together.
    let path = scratch("golden.h5");
    let s = spec(
        vec![axis("x", vec![0.0, 1.0])],
        vec![value("t", vec![10.0, 20.0], "lin-lin", 0.5)],
    );
    let digest = write_table(&path, "/g", &s).expect("write");
    assert_eq!(
        digest,
        "sha256:a3b0bc5987ad11ba166d006ab7103c16ff1d3a8f17f64a248d4b11d8091023f2"
    );
}

// --- Post-review seam regressions --------------------------------------------

#[test]
fn review_newline_in_provenance_changes_the_digest() {
    // v1 digest terminated strings with \n, so a newline INSIDE a field
    // shifted bytes into the next field and two different tables hashed
    // identically. v2 length-prefixes every string.
    let path = scratch("nl.h5");
    let base = spec(
        vec![axis("x", vec![0.0, 1.0])],
        vec![value("t", vec![1.0, 2.0], "lin-lin", 1e-9)],
    );
    let mut a = base.clone();
    a.provenance.producer = "cea\nv2".into();
    a.provenance.producer_version = "1.0".into();
    let mut b = base;
    b.provenance.producer = "cea".into();
    b.provenance.producer_version = "v2\n1.0".into();
    let da = write_table(&path, "/a", &a).expect("write a");
    let db = write_table(&path, "/b", &b).expect("write b");
    assert_ne!(da, db, "different provenance must hash differently");
}

#[test]
fn review_sigma_prefix_is_reserved_and_orphans_refuse() {
    // Writer side: a physics value named sigma_* is refused up front.
    let path = scratch("sigma.h5");
    let s = spec(
        vec![axis("e", vec![0.0, 1.0])],
        vec![value("sigma_t", vec![1.0, 2.0], "lin-lin", 1e-9)],
    );
    assert!(matches!(
        write_table(&path, "/xs", &s),
        Err(TableError::UnexpectedMember { .. })
    ));

    // Reader side: an orphan sigma_ dataset smuggled in refuses at load
    // (v1 silently ignored it AND left it out of the digest).
    let ok = spec(
        vec![axis("e", vec![0.0, 1.0])],
        vec![value("t", vec![1.0, 2.0], "lin-lin", 1e-9)],
    );
    let d = write_table(&path, "/t", &ok).expect("write");
    {
        let f = hdf5::File::append(&path).expect("reopen");
        f.group("/t/values")
            .expect("group")
            .new_dataset_builder()
            .with_data(&[9.0f64])
            .create("sigma_ghost")
            .expect("orphan");
    }
    assert!(matches!(
        Table::open(&path, "/t", &pin(&d)),
        Err(TableError::UnexpectedMember { .. })
    ));
}

#[test]
fn review_unaccounted_datasets_refuse_instead_of_evading_the_digest() {
    let path = scratch("orphan.h5");
    let ok = spec(
        vec![axis("e", vec![0.0, 1.0])],
        vec![value("t", vec![1.0, 2.0], "lin-lin", 1e-9)],
    );
    let d = write_table(&path, "/t", &ok).expect("write");
    // A dataset under axes/ not listed in axis_order.
    {
        let f = hdf5::File::append(&path).expect("reopen");
        f.group("/t/axes")
            .expect("group")
            .new_dataset_builder()
            .with_data(&[0.0f64, 1.0])
            .create("stowaway")
            .expect("extra axis ds");
    }
    assert!(matches!(
        Table::open(&path, "/t", &pin(&d)),
        Err(TableError::UnexpectedMember { .. })
    ));

    // An unknown top-level member of the table group.
    let d2 = write_table(
        &path,
        "/u",
        &spec(
            vec![axis("e", vec![0.0, 1.0])],
            vec![value("t", vec![1.0, 2.0], "lin-lin", 1e-9)],
        ),
    )
    .expect("write");
    {
        let f = hdf5::File::append(&path).expect("reopen");
        f.group("/u")
            .expect("group")
            .create_group("extra")
            .expect("subgroup");
    }
    assert!(matches!(
        Table::open(&path, "/u", &pin(&d2)),
        Err(TableError::UnexpectedMember { .. })
    ));
}

#[test]
fn review_both_sigma_forms_refuse_at_the_writer() {
    // The digest hashes exactly one sigma marker (digest.rs doc); a value
    // carrying both would put unpinned bytes in the file. The reader holds
    // the mirror check for foreign files.
    let path = scratch("both_sigma.h5");
    let mut v = value("t", vec![1.0, 2.0], "lin-lin", 1e-9);
    v.sigma = Some(vec![0.1, 0.2]);
    v.sigma_scalar = Some(0.1);
    let s = spec(vec![axis("x", vec![0.0, 1.0])], vec![v]);
    match write_table(&path, "/b", &s) {
        Err(TableError::UnexpectedMember { reason, .. }) => {
            assert!(reason.contains("exactly one form"), "{reason}")
        }
        other => panic!("expected refusal, got {other:?}"),
    }
}

#[test]
fn review_sigma_length_mismatch_refuses_at_the_writer() {
    let path = scratch("short_sigma.h5");
    let mut v = value("t", vec![1.0, 2.0], "lin-lin", 1e-9);
    v.sigma = Some(vec![0.1]);
    let s = spec(vec![axis("x", vec![0.0, 1.0])], vec![v]);
    match write_table(&path, "/b", &s) {
        Err(TableError::ShapeMismatch {
            expected, found, ..
        }) => {
            assert_eq!((expected, found), (2, 1));
        }
        other => panic!("expected ShapeMismatch, got {other:?}"),
    }
}
