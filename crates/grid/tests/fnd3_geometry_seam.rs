//! FND-3 §3.3 → FND-2 §3.6 ingest-seam tests: geometry-mode build
//! validation (Gas ⇔ κ > 0, shared-face bitwise coherence, covered-face
//! rule), accessor defaults on full-box worlds, and the wall-closure
//! identity (zero for uncut interior cells; the exact covered-face vector
//! in the stair-degenerate case).

use crucible_grid::{CellGeom, FaceDir, Grid, GridSpec, Region};

fn spec(n_r: usize, n_z: usize) -> GridSpec {
    GridSpec {
        r_min: 0.0,
        dr: 0.1,
        n_r,
        z_min: 0.0,
        dz: 0.1,
        n_z,
        n_theta_max: 1,
        axisymmetry_assertion: true,
    }
}

/// A cylinder wall at r = 0.55 across a 8×8 world: rings 0–4 gas (ring 5
/// cut with κ = 0.5-ish), ring 5+ exterior. Face-coordinate arithmetic
/// shared by both sides.
fn cylinder_geom(i_r: usize, _i_z: usize) -> CellGeom {
    let (r0, r1) = (i_r as f64 * 0.1, (i_r + 1) as f64 * 0.1);
    let w = 0.55f64;
    let x = w.clamp(r0, r1);
    let kappa = ((x * x - r0 * r0) / (r1 * r1 - r0 * r0)).clamp(0.0, 1.0);
    let ap_r = |rf: f64| if w > rf { 1.0 } else { 0.0 };
    if kappa > 0.0 {
        CellGeom {
            region: Region::Gas,
            kappa,
            aperture: [ap_r(r0), ap_r(r1), kappa, kappa],
        }
    } else {
        CellGeom {
            region: Region::Exterior,
            kappa: 0.0,
            aperture: [0.0; 4],
        }
    }
}

#[test]
fn geometry_build_masks_follow_kappa() {
    let g = Grid::build_with_geometry(spec(8, 8), &["q"], cylinder_geom).expect("builds");
    assert!(g.has_cut_geometry());
    for i_r in 0..8 {
        let active = g.is_active(i_r, 3);
        assert_eq!(active, i_r <= 5, "ring {i_r} activity");
        if active {
            assert!(g.kappa(i_r, 3) > 0.0);
        } else {
            assert_eq!(g.kappa(i_r, 3), 0.0);
        }
    }
    // The cut ring's κ: (0.55² − 0.5²)/(0.6² − 0.5²) = 0.0525/0.11.
    let exact = (0.55f64 * 0.55 - 0.25) / (0.36 - 0.25);
    assert!((g.kappa(5, 0) - exact).abs() < 1e-15);
}

#[test]
fn full_box_worlds_default_to_unit_kappa_and_apertures() {
    let g = Grid::build(spec(8, 8), &["q"]).expect("builds");
    assert!(!g.has_cut_geometry());
    assert_eq!(g.kappa(3, 3), 1.0);
    assert_eq!(g.aperture(3, 3, FaceDir::RPlus), 1.0);
    let (w_r, w_z) = g.wall_closure(3, 3, 1);
    assert_eq!((w_r, w_z), (0.0, 0.0), "no closure on full boxes");
}

#[test]
fn wall_closure_zero_interior_exact_stair_on_covered_faces() {
    let g = Grid::build_with_geometry(spec(8, 8), &["q"], cylinder_geom).expect("builds");
    // Interior full cell: exactly zero (bitwise cancellation).
    assert_eq!(g.wall_closure(2, 3, 1), (0.0, 0.0));
    // The cut ring (5): W_r = 0·A_out − 1·A_in − κ(A_out − A_in) < 0
    // (outward wall), W_z = 0 (cylinder wall parallel to z).
    let (w_r, w_z) = g.wall_closure(5, 3, 1);
    assert_eq!(w_z, 0.0);
    let a_in = g.face_area_r(5, false, 1);
    let a_out = g.face_area_r(5, true, 1);
    let expected = -a_in - g.kappa(5, 3) * (a_out - a_in);
    assert!((w_r - expected).abs() < 1e-15 * a_in.abs());
    // The closure magnitude ≈ the cylinder lateral area τ·w·dz.
    let lateral = std::f64::consts::TAU * 0.55 * 0.1;
    assert!(
        (w_r.abs() - lateral).abs() < 0.01 * lateral,
        "|W_r| {} vs lateral {lateral}",
        w_r.abs()
    );
}

#[test]
fn incoherent_suppliers_refuse() {
    // Gas with κ = 0: refused.
    let bad_region = |_: usize, _: usize| CellGeom {
        region: Region::Gas,
        kappa: 0.0,
        aperture: [0.0; 4],
    };
    assert!(Grid::build_with_geometry(spec(4, 4), &["q"], bad_region).is_err());

    // Shared-face mismatch: cell (1, z) says its r+ is 0.7, cell (2, z)
    // says its r− is 0.6 — refused.
    let mismatch = |i_r: usize, _: usize| CellGeom {
        region: Region::Gas,
        kappa: 1.0,
        aperture: if i_r == 1 {
            [1.0, 0.7, 1.0, 1.0]
        } else if i_r == 2 {
            [0.6, 1.0, 1.0, 1.0]
        } else {
            [1.0; 4]
        },
    };
    assert!(Grid::build_with_geometry(spec(4, 4), &["q"], mismatch).is_err());

    // Open aperture against a κ = 0 cell: refused.
    let open_into_wall = |i_r: usize, _: usize| {
        if i_r < 2 {
            CellGeom {
                region: Region::Gas,
                kappa: 1.0,
                aperture: [1.0; 4],
            }
        } else {
            CellGeom {
                region: Region::Exterior,
                kappa: 0.0,
                aperture: [0.0; 4],
            }
        }
    };
    // Cell 1's r+ = 1.0 faces cell 2's κ = 0 — and cell 2 reports 0.0, so
    // this trips BOTH the mismatch and covered-face rules; either way loud.
    assert!(Grid::build_with_geometry(spec(4, 4), &["q"], open_into_wall).is_err());

    // Geometry at N_θ > 1: deferred loudly.
    let mut s = spec(4, 4);
    s.n_theta_max = 4;
    s.axisymmetry_assertion = false;
    let full = |_: usize, _: usize| CellGeom {
        region: Region::Gas,
        kappa: 1.0,
        aperture: [1.0; 4],
    };
    assert!(Grid::build_with_geometry(s, &["q"], full).is_err());
}
