//! FND-3 §3.3 → FND-2 §3.6 ingest-seam tests: geometry-mode build
//! validation (Gas ⇔ κ > 0, shared-face bitwise coherence, covered-face
//! rule), accessor defaults on full-box worlds, and the wall-closure
//! identity (zero for uncut interior cells; the exact covered-face vector
//! in the stair-degenerate case). S9 (FND-2 §3.4 3-D aperture wave):
//! `build_with_geometry_theta` validation — the per-sector scope rule
//! (θ-sector coverage inside a gas ring refuses), θ-face bitwise
//! coherence, the geometry-floor coarsen/refine refusal — and the
//! `wall_closure_cell` θ-limb identities.

use crucible_grid::{CellGeom, CellGeomTheta, FaceDir, Grid, GridError, GridSpec, Region};

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

    // The revolved per-(r,z) path at N_θ > 1: refused, pointing at the S9
    // 3-D entry point.
    let mut s = spec(4, 4);
    s.n_theta_max = 4;
    s.axisymmetry_assertion = false;
    let full = |_: usize, _: usize| CellGeom {
        region: Region::Gas,
        kappa: 1.0,
        aperture: [1.0; 4],
    };
    let err = Grid::build_with_geometry(s, &["q"], full).expect_err("must refuse");
    assert!(
        format!("{err}").contains("build_with_geometry_theta"),
        "refusal must name the 3-D entry point: {err}"
    );
}

// --- S9: the θ-sector builder (FND-3 §3.3 3-D apertures) --------------------

fn spec_theta(n_r: usize, n_z: usize, n_theta: u32) -> GridSpec {
    GridSpec {
        r_min: 0.0,
        dr: 0.1,
        n_r,
        z_min: 0.0,
        dz: 0.1,
        n_z,
        n_theta_max: n_theta,
        axisymmetry_assertion: false,
    }
}

/// A θ-varying legal world on 4×4×N_θ4: every cell gas, per-sector κ in
/// (0, 1], θ-face apertures computed from the FACE index (coherent by
/// construction), r/z faces fully open in the interior.
fn bumpy_supplier(i_r: usize, j: u32, _i_z: usize) -> CellGeomTheta {
    let nt = 4u32;
    // θ-face aperture as a function of the face (between sector fi−1 and
    // fi): one value per face index, shared by both adjacent sectors.
    let theta_face = |fi: u32| 0.7 + 0.05 * f64::from((fi + i_r as u32) % nt);
    CellGeomTheta {
        kappa: 0.6 + 0.08 * f64::from(j),
        aperture: [
            1.0,
            1.0,
            1.0,
            1.0,
            theta_face(j),            // θ−: the face below sector j
            theta_face((j + 1) % nt), // θ+: the face above sector j
        ],
    }
}

#[test]
fn fnd3_s33_theta_builder_accepts_a_theta_varying_world() {
    let g = Grid::build_with_geometry_theta(spec_theta(4, 4, 4), &["q"], bumpy_supplier, |_, _| {
        Region::Exterior
    })
    .expect("legal θ-varying world builds");
    assert!(g.has_cut_geometry());
    // Per-sector accessors read the sector's own values; plane-0 forms read
    // sector 0 (the documented contract).
    assert_eq!(g.kappa_at(1, 3, 2), 0.6 + 0.08 * 3.0);
    assert_eq!(g.kappa(1, 2), 0.6);
    assert_eq!(
        g.aperture_at(1, 2, 2, FaceDir::ThetaMinus).to_bits(),
        g.aperture_at(1, 1, 2, FaceDir::ThetaPlus).to_bits(),
        "one θ-face, one bit pattern on both sides"
    );
    // The geometry floor pins every brick at the built N_θ.
    for b in g.bricks() {
        assert_eq!(b.n_theta_geom_floor(), 4);
    }
}

#[test]
fn fnd3_s33_theta_sector_coverage_inside_a_gas_ring_refuses_typed() {
    // Sector 2 of one ring cell fully covered (κ = 0) while other sectors
    // hold gas — the S9 scope rule's typed refusal (S10/S11 wave).
    let pillar = |i_r: usize, j: u32, i_z: usize| {
        let mut c = bumpy_supplier(i_r, j, i_z);
        if (i_r, i_z) == (2, 2) && j == 2 {
            c.kappa = 0.0;
        }
        c
    };
    let err = Grid::build_with_geometry_theta(spec_theta(4, 4, 4), &["q"], pillar, |_, _| {
        Region::Exterior
    })
    .expect_err("sector coverage must refuse");
    match err {
        GridError::ThetaSectorCovered { i_r, i_z, i_theta } => {
            assert_eq!((i_r, i_z, i_theta), (2, 2, 2));
        }
        other => panic!("wrong refusal type: {other}"),
    }
    assert!(
        format!(
            "{}",
            GridError::ThetaSectorCovered {
                i_r: 2,
                i_z: 2,
                i_theta: 2
            }
        )
        .contains("per-sector activity masks"),
        "the refusal must name the missing capability"
    );
}

#[test]
fn fnd3_s33_theta_face_bit_mismatch_refuses() {
    // Sector 1's θ+ disagrees with sector 2's θ− by one ulp-scale amount:
    // the bitwise shared-face rule refuses.
    let mismatch = |i_r: usize, j: u32, i_z: usize| {
        let mut c = bumpy_supplier(i_r, j, i_z);
        if (i_r, i_z) == (1, 1) && j == 1 {
            c.aperture[FaceDir::ThetaPlus.index()] += 1e-12;
        }
        c
    };
    assert!(
        Grid::build_with_geometry_theta(spec_theta(4, 4, 4), &["q"], mismatch, |_, _| {
            Region::Exterior
        })
        .is_err(),
        "θ-face bit mismatch must refuse"
    );
}

#[test]
fn fnd3_s34_theta_regrid_refuses_on_geometry_bricks() {
    // FND-3 §3.4 geometry floor, S9: the controller can never re-grid cut
    // geometry — coarsen AND refine refuse, naming plan S11.
    let mut g = Grid::build_with_geometry_theta(
        spec_theta(4, 4, 8),
        &["q"],
        |i_r, j, i_z| bumpy_supplier(i_r, j % 4, i_z),
        |_, _| Region::Exterior,
    )
    .expect("builds at N_θ = 8");
    let coarsen = g.coarsen_theta(0, None).expect_err("coarsen must refuse");
    assert!(
        format!("{coarsen}").contains("S11"),
        "coarsen refusal names S11: {coarsen}"
    );
    let refine = g.refine_theta(0).expect_err("refine must refuse");
    assert!(
        format!("{refine}").contains("S11"),
        "refine refusal names S11: {refine}"
    );
}

#[test]
fn fnd3_s33_wall_closure_cell_theta_limb_exact_stair_identity() {
    // A hand-built θ-varying world: cell (1, 1)'s sector 1 has θ− = 0.25
    // and θ+ = 1.0 ⇒ W_θ = (1.0 − 0.25)·A_θ exactly; its r/z faces are
    // symmetric so W_z = 0 and W_r keeps the revolved identity. Uncut
    // interior cells return the exact (0, 0, 0).
    let nt = 4u32;
    let supplier = |i_r: usize, j: u32, i_z: usize| {
        if (i_r, i_z) == (1, 1) {
            // θ-face apertures per face index: face 1 (between sectors 0
            // and 1) is 0.25, every other face 1.0.
            let face = |fi: u32| if fi == 1 { 0.25 } else { 1.0 };
            CellGeomTheta {
                kappa: 0.8,
                aperture: [1.0, 1.0, 1.0, 1.0, face(j), face((j + 1) % nt)],
            }
        } else {
            CellGeomTheta {
                kappa: 1.0,
                aperture: [1.0; 6],
            }
        }
    };
    let g = Grid::build_with_geometry_theta(spec_theta(4, 4, nt), &["q"], supplier, |_, _| {
        Region::Exterior
    })
    .expect("builds");
    let a_th = g.face_area_theta();
    // Sector 1: θ− is face 1 (0.25), θ+ is face 2 (1.0).
    let (_, w_theta, w_z) = g.wall_closure_cell(1, 1, 1, nt);
    assert_eq!(w_z, 0.0);
    assert_eq!(
        w_theta.to_bits(),
        ((1.0 - 0.25) * a_th).to_bits(),
        "the θ-limb is the exact stair identity (a_θ₊ − a_θ₋)·A_θ"
    );
    // Sector 0: θ− is face 0 (1.0), θ+ is face 1 (0.25) — the mirror sign.
    let (_, w_theta0, _) = g.wall_closure_cell(1, 0, 1, nt);
    assert_eq!(w_theta0.to_bits(), ((0.25 - 1.0) * a_th).to_bits());
    // Uncut interior: exact (0, 0, 0) by bitwise cancellation.
    assert_eq!(g.wall_closure_cell(2, 2, 2, nt), (0.0, 0.0, 0.0));
}
