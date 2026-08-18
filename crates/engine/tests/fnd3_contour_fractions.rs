//! FND-3 §6.1-class validation of the analytic contour clip: recovered
//! fractions of a revolved linear profile (cone frustum) reproduce the
//! exact volume to round-off — the analytic path's version of the
//! "analytic volumes" rung (no sampling, so no convergence study: the
//! answer is closed-form and must be exact).

use crucible_engine::geometry::Contour;

fn cone() -> Contour {
    // r_wall: 0.5 → 1.0 over z ∈ [0, 1] (linear), no liner.
    Contour::new(vec![(0.0, 0.5), (1.0, 1.0)], 0.0).expect("valid contour")
}

#[test]
fn frustum_volume_recovered_exactly() {
    let c = cone();
    // Incommensurate cell sizes so faces never align with the wall.
    let (dr, dz) = (0.11, 0.13);
    let (n_r, n_z) = (11, 8); // covers r ∈ [0, 1.21], z ∈ [0, 1.04]
    let mut vol = 0.0f64;
    for i_r in 0..n_r {
        for i_z in 0..n_z {
            let (r0, r1) = (i_r as f64 * dr, (i_r + 1) as f64 * dr);
            let (z0, z1) = ((i_z as f64 * dz).min(1.0), ((i_z + 1) as f64 * dz).min(1.0));
            if z1 <= z0 {
                continue;
            }
            let kappa = c.gas_volume_fraction(r0, r1, z0, z1);
            vol += kappa * 0.5 * (r1 * r1 - r0 * r0) * (z1 - z0);
        }
    }
    // ∫₀¹ ½ w(z)² dz, w linear 0.5 → 1.0: ½·(wa² + wa·wb + wb²)/3.
    let exact = 0.5 * (0.25 + 0.5 + 1.0) / 3.0;
    assert!(
        ((vol - exact) / exact).abs() < 1e-13,
        "frustum volume {vol} vs exact {exact}"
    );
}

#[test]
fn kappa_is_additive_across_a_station_kink() {
    // A kink at z = 0.3 inside the cell [0.55, 0.75] × [0.0, 1.0]:
    // the z-subdivided fractions must recombine exactly.
    let c = Contour::new(vec![(0.0, 0.6), (0.3, 0.72), (1.0, 0.58)], 0.0).expect("valid");
    let (r0, r1) = (0.55, 0.75);
    let full = c.gas_volume_fraction(r0, r1, 0.0, 1.0);
    let a = c.gas_volume_fraction(r0, r1, 0.0, 0.3);
    let b = c.gas_volume_fraction(r0, r1, 0.3, 1.0);
    let recombined = a * 0.3 + b * 0.7;
    assert!(
        (full - recombined).abs() < 1e-14,
        "κ additivity: {full} vs {recombined}"
    );
}

#[test]
fn face_apertures_match_closed_forms() {
    let c = cone();
    // w crosses r_f = 0.75 at z = 0.5 exactly.
    let a = c.r_face_aperture(0.75, 0.0, 1.0);
    assert!((a - 0.5).abs() < 1e-14, "r-face aperture {a}");
    // Below the wall everywhere → fully open; above → fully covered.
    assert_eq!(c.r_face_aperture(0.4, 0.0, 1.0), 1.0);
    assert_eq!(c.r_face_aperture(1.1, 0.0, 1.0), 0.0);
    // Sub-span where w ∈ [0.6, 0.7] stays below 0.75 → covered face.
    assert_eq!(c.r_face_aperture(0.75, 0.2, 0.4), 0.0);
    // z-face at z = 0.5 (w = 0.75), ring [0.7, 0.8]:
    // (0.75² − 0.7²)/(0.8² − 0.7²).
    let az = c.z_face_aperture(0.7, 0.8, 0.5);
    let exact = (0.75f64 * 0.75 - 0.49) / (0.64 - 0.49);
    assert!(
        (az - exact).abs() < 1e-14,
        "z-face aperture {az} vs {exact}"
    );
    // Fully-gas and fully-covered rings.
    assert_eq!(c.z_face_aperture(0.1, 0.2, 0.5), 1.0);
    assert_eq!(c.z_face_aperture(0.9, 1.0, 0.5), 0.0);
}

#[test]
fn coherence_kappa_zero_iff_faces_covered() {
    // The consistency the grid seam validates: a cell with κ = 0 must see
    // covered faces from its gas neighbor's side. Probe the boundary ring
    // along the cone.
    let c = cone();
    let (dr, dz) = (0.11, 0.13);
    for i_r in 0..11usize {
        for i_z in 0..7usize {
            let (r0, r1) = (i_r as f64 * dr, (i_r + 1) as f64 * dr);
            let (z0, z1) = (i_z as f64 * dz, (i_z + 1) as f64 * dz);
            let kappa = c.gas_volume_fraction(r0, r1, z0, z1);
            if kappa == 0.0 {
                // Its inner face (the gas neighbor's outer face) is covered.
                assert_eq!(
                    c.r_face_aperture(r0, z0, z1),
                    0.0,
                    "cell ({i_r},{i_z}) has κ = 0 but an open inner face"
                );
            }
        }
    }
}
