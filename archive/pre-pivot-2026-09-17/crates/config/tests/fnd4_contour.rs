//! FND-3/FND-4 contour-of-revolution grammar: the fidelity dial derives the
//! grid from the cited r(z) CSV (content-addressed like a table pin), the
//! resolved form replays purely, and every misdeclaration refuses loudly.
//! Exercised against the REAL geometry-of-record
//! (`data/anchors/rl10_contour.csv`).

use crucible_config::{INCH_M, load_str, load_str_with_sidecars};
use crucible_registry::Registry;

static EMPTY: [&crucible_registry::Manifest; 0] = [];

fn registry() -> Registry {
    Registry::new(&EMPTY, &EMPTY)
}

const ROOT: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../..");

fn fs_reader(path: &str) -> Result<String, String> {
    std::fs::read_to_string(format!("{ROOT}/{path}")).map_err(|e| e.to_string())
}

fn rl10_geometry_author(dial: f64) -> String {
    format!(
        r#"
        schema_version = 1
        [geometry]
        n_theta_max = 1
        axisymmetric = true
        contour = "data/anchors/rl10_contour.csv"
        contour_units = "in"
        liner_thickness_m = 0.003
        cells_across_throat = {dial}
        [rng]
        master_seed = 7
    "#
    )
}

#[test]
fn contour_derives_the_grid_from_the_dial() {
    let loaded = load_str_with_sidecars(&rl10_geometry_author(5.0), &registry(), &fs_reader)
        .expect("RL10 contour config loads");
    let geo = loaded.resolved.geometry.as_ref().expect("geometry");
    let ext = geo.extents.as_ref().expect("derived extents");
    let contour = geo.contour.as_ref().expect("resolved contour");

    // Dial semantics: cell size = r_throat / dial, exactly.
    let r_throat = 2.47 * INCH_M; // RADIUS of record (CSV-header erratum)
    assert!((ext.dr - r_throat / 5.0).abs() < 1e-15, "dr = {}", ext.dr);
    // Radial cells cover exit radius + liner.
    let r_exit = 19.2914 * INCH_M;
    assert_eq!(ext.n_r, ((r_exit + 0.003) / ext.dr).ceil() as i64);
    // Axial span: injector plane (−13 in) to exit closure (+45.3 in).
    assert!((ext.z_min - (-13.0 * INCH_M)).abs() < 1e-12);
    let span = (45.3 + 13.0) * INCH_M;
    assert_eq!(ext.n_z, (span / ext.dr).ceil() as i64);
    assert!(
        (ext.dz * ext.n_z as f64 - span).abs() < 1e-12,
        "dz covers the span exactly"
    );
    assert!(
        ext.n_z > 100 && ext.n_z < 140,
        "coarse-tier axial count: {}",
        ext.n_z
    );
    assert!(
        contour.contour_digest.starts_with("sha256:"),
        "geometry is content-addressed"
    );

    // Replay: the resolved form (extents + digest materialized) reloads
    // through the PURE entry point and reproduces itself.
    let round = load_str(&loaded.resolved.to_toml(), &registry())
        .expect("resolved contour config replays without file access");
    assert_eq!(round.resolved, loaded.resolved);
}

#[test]
fn the_dial_scales_compute_not_geometry() {
    let coarse = load_str_with_sidecars(&rl10_geometry_author(5.0), &registry(), &fs_reader)
        .expect("coarse loads");
    let full = load_str_with_sidecars(&rl10_geometry_author(16.0), &registry(), &fs_reader)
        .expect("full loads");
    let (ce, fe) = (
        coarse
            .resolved
            .geometry
            .as_ref()
            .unwrap()
            .extents
            .as_ref()
            .unwrap(),
        full.resolved
            .geometry
            .as_ref()
            .unwrap()
            .extents
            .as_ref()
            .unwrap(),
    );
    // Same world, finer cells: spans agree, counts scale ~16/5.
    assert!((ce.z_min - fe.z_min).abs() < 1e-12);
    let ratio = fe.n_z as f64 / ce.n_z as f64;
    assert!((ratio - 3.2).abs() < 0.05, "n_z ratio {ratio}");
}

#[test]
fn contour_refusals_are_diagnosed() {
    // Missing dial.
    let author = r#"
        schema_version = 1
        [geometry]
        n_theta_max = 1
        axisymmetric = true
        contour = "data/anchors/rl10_contour.csv"
        contour_units = "in"
        liner_thickness_m = 0.003
    "#;
    let err = load_str_with_sidecars(author, &registry(), &fs_reader).unwrap_err();
    assert!(format!("{err}").contains("cells_across_throat"), "{err}");

    // Unknown units.
    let author = rl10_geometry_author(5.0).replace("\"in\"", "\"furlong\"");
    let err = load_str_with_sidecars(&author, &registry(), &fs_reader).unwrap_err();
    assert!(format!("{err}").contains("furlong"), "{err}");

    // Stale author digest (same-label/different-bytes doctrine).
    let author = rl10_geometry_author(5.0).replace(
        "contour_units = \"in\"",
        "contour_units = \"in\"\n        contour_digest = \"sha256:deadbeef\"",
    );
    let err = load_str_with_sidecars(&author, &registry(), &fs_reader).unwrap_err();
    assert!(format!("{err}").contains("stale"), "{err}");

    // Explicit extents + contour without a digest: not an author form.
    let author = rl10_geometry_author(5.0).replace(
        "[rng]",
        "r_min = 0.0\n        dr = 0.01\n        n_r = 4\n        z_min = 0.0\n        dz = 0.01\n        n_z = 4\n        [rng]",
    );
    let err = load_str_with_sidecars(&author, &registry(), &fs_reader).unwrap_err();
    assert!(format!("{err}").contains("replay form"), "{err}");

    // Contour-companion keys without a contour.
    let author = r#"
        schema_version = 1
        [geometry]
        n_theta_max = 1
        axisymmetric = true
        cells_across_throat = 5.0
    "#;
    let err = load_str(author, &registry()).unwrap_err();
    assert!(format!("{err}").contains("require a `contour`"), "{err}");

    // The pure entry point cannot read the CSV.
    let err = load_str(&rl10_geometry_author(5.0), &registry()).unwrap_err();
    assert!(format!("{err}").contains("load_str_with_sidecars"), "{err}");
}

#[test]
fn contour_dial_below_floor_refuses() {
    let err =
        load_str_with_sidecars(&rl10_geometry_author(1.0), &registry(), &fs_reader).unwrap_err();
    assert!(format!("{err}").contains("cells_across_throat"), "{err}");
}
