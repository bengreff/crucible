//! FND-3 §3/§6 gate battery for the S9 `geom3d` kernel: analytic-volume
//! convergence at the declared jittered rate (§6.1 — the C_JITTER
//! calibration), imperfect-STL winding robustness where ray-parity fails
//! (§6.2), PLIC area/volume consistency (§6.3), CSG↔STL cross-check
//! (§6.4), determinism (§6.5), and the cylindrical metric + azimuthal
//! geometry floors (§6.6).

use crucible_grid::geom3d::{
    Sdf, Solid, TriMesh, VoxelSpec, VoxelWorld, theta_geom_floor, voxelize, voxelize_with_samples,
};
use std::f64::consts::{PI, TAU};

// ---------------------------------------------------------------- helpers

fn cell_volume(spec: &VoxelSpec, i_r: usize, i_z: usize) -> f64 {
    let r0 = spec.r_min + i_r as f64 * spec.dr;
    let r1 = r0 + spec.dr;
    let _ = i_z;
    0.5 * (r1 * r1 - r0 * r0) * (TAU / f64::from(spec.n_theta)) * spec.dz
}

/// Exact SOLID fraction of the ring cell [r0,r1]×[z0,z1] (any θ) cut by a
/// sphere of radius `rad` centered on the axis at height `zc`:
/// frac = ∫ (clamp(ρ²(z), r0², r1²) − r0²) dz / ((r1²−r0²)(z1−z0)) with
/// ρ²(z) = rad² − (z−zc)² — piecewise quadratic, closed form per regime.
fn sphere_ring_solid_frac(r0: f64, r1: f64, z0: f64, z1: f64, zc: f64, rad: f64) -> f64 {
    let mut cuts = vec![z0, z1];
    for rc in [r0, r1] {
        if rad * rad > rc * rc {
            let dz = (rad * rad - rc * rc).sqrt();
            for z in [zc - dz, zc + dz] {
                if z > z0 && z < z1 {
                    cuts.push(z);
                }
            }
        }
    }
    // Sphere edge itself (ρ² = 0 crossing).
    for z in [zc - rad, zc + rad] {
        if z > z0 && z < z1 {
            cuts.push(z);
        }
    }
    cuts.sort_by(|a, b| a.partial_cmp(b).expect("finite z"));
    let mut acc = 0.0f64;
    for w in cuts.windows(2) {
        let (za, zb) = (w[0], w[1]);
        let len = zb - za;
        if len <= 0.0 {
            continue;
        }
        let zm = 0.5 * (za + zb);
        let rho2 = rad * rad - (zm - zc) * (zm - zc);
        if rho2 <= r0 * r0 {
            // no overlap on this span
        } else if rho2 >= r1 * r1 {
            acc += (r1 * r1 - r0 * r0) * len;
        } else {
            // ∫ (rad² − (z−zc)² − r0²) dz over [za, zb]
            let a3 = (zb - zc).powi(3) - (za - zc).powi(3);
            acc += (rad * rad - r0 * r0) * len - a3 / 3.0;
        }
    }
    acc / ((r1 * r1 - r0 * r0) * (z1 - z0))
}

/// Least-squares slope of ln(y) vs ln(x).
fn log_slope(xs: &[f64], ys: &[f64]) -> f64 {
    let n = xs.len() as f64;
    let lx: Vec<f64> = xs.iter().map(|v| v.ln()).collect();
    let ly: Vec<f64> = ys.iter().map(|v| v.ln()).collect();
    let mx = lx.iter().sum::<f64>() / n;
    let my = ly.iter().sum::<f64>() / n;
    let num: f64 = lx.iter().zip(&ly).map(|(x, y)| (x - mx) * (y - my)).sum();
    let den: f64 = lx.iter().map(|x| (x - mx) * (x - mx)).sum();
    num / den
}

/// The outward-oriented unit tetrahedron (0,0,0)-(1,0,0)-(0,1,0)-(0,0,1)
/// as four triangles; `holed` drops the slanted face.
fn tet_tris(holed: bool) -> Vec<[[f64; 3]; 3]> {
    let v0 = [0.0, 0.0, 0.0];
    let v1 = [1.0, 0.0, 0.0];
    let v2 = [0.0, 1.0, 0.0];
    let v3 = [0.0, 0.0, 1.0];
    let mut tris = vec![[v0, v2, v1], [v0, v1, v3], [v0, v3, v2]];
    if !holed {
        tris.push([v1, v2, v3]);
    }
    tris
}

/// Ray-parity along +z from `p` (the classifier FND-3 §3.2 rejects as
/// primary): count strict upward crossings by 2-D barycentric test.
fn z_ray_parity_inside(tris: &[[[f64; 3]; 3]], p: [f64; 3]) -> bool {
    let mut crossings = 0usize;
    for t in tris {
        let (a, b, c) = (t[0], t[1], t[2]);
        let det = (b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]);
        if det == 0.0 {
            continue; // vertical triangle: degenerate in projection
        }
        let l1 = ((p[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (p[1] - a[1])) / det;
        let l2 = ((b[0] - a[0]) * (p[1] - a[1]) - (p[0] - a[0]) * (b[1] - a[1])) / det;
        let l0 = 1.0 - l1 - l2;
        if l0 > 0.0 && l1 > 0.0 && l2 > 0.0 {
            let z_hit = l0 * a[2] + l1 * b[2] + l2 * c[2];
            if z_hit > p[2] {
                crossings += 1;
            }
        }
    }
    crossings % 2 == 1
}

/// 12-triangle outward-oriented mesh of an axis-aligned box.
fn box_tri_mesh(min: [f64; 3], max: [f64; 3]) -> TriMesh {
    let v = |ix: usize, iy: usize, iz: usize| -> [f64; 3] {
        [
            if ix == 0 { min[0] } else { max[0] },
            if iy == 0 { min[1] } else { max[1] },
            if iz == 0 { min[2] } else { max[2] },
        ]
    };
    // Each face as two triangles, outward-oriented (CCW from outside).
    let quads: [[[usize; 3]; 4]; 6] = [
        [[0, 0, 0], [0, 0, 1], [0, 1, 1], [0, 1, 0]], // x−
        [[1, 0, 0], [1, 1, 0], [1, 1, 1], [1, 0, 1]], // x+
        [[0, 0, 0], [1, 0, 0], [1, 0, 1], [0, 0, 1]], // y−
        [[0, 1, 0], [0, 1, 1], [1, 1, 1], [1, 1, 0]], // y+
        [[0, 0, 0], [0, 1, 0], [1, 1, 0], [1, 0, 0]], // z−
        [[0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], // z+
    ];
    let mut tris = Vec::with_capacity(12);
    for [a, b, c, d] in quads {
        let (pa, pb, pc, pd) = (
            v(a[0], a[1], a[2]),
            v(b[0], b[1], b[2]),
            v(c[0], c[1], c[2]),
            v(d[0], d[1], d[2]),
        );
        tris.push([pa, pb, pc]);
        tris.push([pa, pc, pd]);
    }
    TriMesh {
        tris,
        dropped_degenerate: 0,
    }
}

/// Bit-serialize a world for exact comparison.
fn world_bits(w: &VoxelWorld) -> Vec<u64> {
    let mut out = vec![w.eps_alpha.to_bits()];
    for c in &w.cells {
        out.push(c.kappa.to_bits());
        for a in &c.aperture {
            out.push(a.to_bits());
        }
        match &c.plic {
            None => out.push(0),
            Some(p) => {
                out.push(1);
                for v in p.normal {
                    out.push(v.to_bits());
                }
                out.push(p.d.to_bits());
                for v in p.centroid {
                    out.push(v.to_bits());
                }
                out.push(p.area.to_bits());
            }
        }
    }
    out
}

/// The annular-cylinder fixture shared by gates 2/3/7/8: the revolved
/// rectangle r ∈ [0.35, 0.55], z ∈ [0.25, 0.75].
const RA: f64 = 0.35;
const RB: f64 = 0.55;
const ZA: f64 = 0.25;
const ZB: f64 = 0.75;

fn annular_cylinder() -> Solid {
    Solid::Csg(Sdf::Revolved {
        profile: vec![(RA, ZA), (RB, ZA), (RB, ZB), (RA, ZB)],
    })
}

// ------------------------------------------------------------------ gates

/// §6.1 + §3.1: per-cell fraction errors of an axis-centered sphere decay
/// at the jittered-stratified rate N^(−2/3), and the DECLARED ε_α =
/// C_JITTER·N^(−2/3) bounds the measured per-cell max error at N = 512 —
/// this is the C_JITTER calibration of record (the measured constant is
/// printed; the shipped constant is measured × 1.5, rounded up).
#[test]
fn sphere_fractions_converge_at_the_jittered_rate() {
    let (zc, rad) = (0.8, 0.5);
    let solid = Solid::Csg(Sdf::Sphere {
        center: [0.0, 0.0, zc],
        radius: rad,
    });
    let spec = VoxelSpec {
        r_min: 0.0,
        dr: 0.05,
        n_r: 16,
        z_min: 0.0,
        dz: 0.05,
        n_z: 32,
        n_theta: 1,
        seed: 0xC0FFEE,
    };
    let ns = [64usize, 512, 4096];
    let mut rms = Vec::new();
    let mut max_err_512 = 0.0f64;
    let mut eps_512 = 0.0f64;
    let mut w512: Option<VoxelWorld> = None;
    for &n in &ns {
        let w = voxelize_with_samples(&solid, &spec, n).expect("voxelizes");
        let mut sum2 = 0.0f64;
        let mut cnt = 0usize;
        let mut maxe = 0.0f64;
        for i_r in 0..spec.n_r {
            for i_z in 0..spec.n_z {
                let r0 = i_r as f64 * spec.dr;
                let z0 = i_z as f64 * spec.dz;
                let exact_solid =
                    sphere_ring_solid_frac(r0, r0 + spec.dr, z0, z0 + spec.dz, zc, rad);
                if exact_solid <= 0.0 || exact_solid >= 1.0 {
                    continue; // pure cells: decided by bound tests, exact
                }
                let idx = w.cell_index(i_r, 0, i_z);
                let err = (w.cells[idx].kappa - (1.0 - exact_solid)).abs();
                sum2 += err * err;
                cnt += 1;
                maxe = maxe.max(err);
            }
        }
        rms.push((sum2 / cnt as f64).sqrt());
        if n == 512 {
            max_err_512 = maxe;
            eps_512 = w.eps_alpha;
            w512 = Some(w);
        }
    }
    // Rate: jittered-stratified indicator sampling in 3-D → N^(−2/3).
    let slope = log_slope(&[64.0, 512.0, 4096.0], &rms);
    println!("jitter rate exponent = {slope:.3} (rms = {rms:?})");
    assert!(
        (-0.85..=-0.5).contains(&slope),
        "jitter rate exponent {slope} outside [-0.85, -0.5]; rms = {rms:?}"
    );
    // Calibration: the declared bound covers the measured max error.
    let c_measured = max_err_512 * 512f64.powf(2.0 / 3.0);
    println!("measured C_jitter = {c_measured:.4} (max per-cell err at 512 = {max_err_512:.3e})");
    assert!(
        max_err_512 <= eps_512,
        "declared eps_alpha {eps_512:.3e} does not bound measured max {max_err_512:.3e}"
    );
    // Total gas volume vs exact, within the coherent worst-case bound.
    let w = w512.expect("512-sample world");
    let mut v_gas = 0.0f64;
    let mut bound = 1e-12;
    for i_r in 0..spec.n_r {
        for i_z in 0..spec.n_z {
            let idx = w.cell_index(i_r, 0, i_z);
            let vc = cell_volume(&spec, i_r, i_z);
            v_gas += w.cells[idx].kappa * vc;
            if w.cells[idx].plic.is_some() {
                bound += w.eps_alpha * vc;
            }
        }
    }
    let r_dom = spec.dr * spec.n_r as f64;
    let v_exact = PI * r_dom * r_dom * (spec.dz * spec.n_z as f64) - 4.0 / 3.0 * PI * rad.powi(3);
    assert!(
        (v_gas - v_exact).abs() <= bound,
        "total gas volume {v_gas} vs exact {v_exact}, bound {bound:.3e}"
    );
}

/// §6.1: exact cone volume + a revolved-rectangle (annular-cylinder)
/// volume recovered within the ε_α-propagated (coherent worst-case) bound.
#[test]
fn cone_and_revolved_profile_volumes_within_bound() {
    let spec = VoxelSpec {
        r_min: 0.0,
        dr: 0.05,
        n_r: 16,
        z_min: 0.0,
        dz: 0.05,
        n_z: 24,
        n_theta: 1,
        seed: 7,
    };
    let h = 0.6f64;
    let r_base = 0.4f64;
    let half_angle = (r_base / h).atan();
    let cases: [(Solid, f64); 2] = [
        (
            Solid::Csg(Sdf::ConeZ {
                apex: [0.0, 0.0, 0.9],
                half_angle,
                height: h,
            }),
            PI * r_base * r_base * h / 3.0,
        ),
        (annular_cylinder(), PI * (RB * RB - RA * RA) * (ZB - ZA)),
    ];
    for (solid, v_exact) in cases {
        let w = voxelize(&solid, &spec).expect("voxelizes");
        let mut v_solid = 0.0f64;
        let mut bound = 1e-12;
        for i_r in 0..spec.n_r {
            for i_z in 0..spec.n_z {
                let idx = w.cell_index(i_r, 0, i_z);
                let vc = cell_volume(&spec, i_r, i_z);
                v_solid += (1.0 - w.cells[idx].kappa) * vc;
                if w.cells[idx].plic.is_some() {
                    bound += w.eps_alpha * vc;
                }
            }
        }
        assert!(
            (v_solid - v_exact).abs() <= bound,
            "solid volume {v_solid} vs exact {v_exact}, bound {bound:.3e}"
        );
    }
}

/// §6.6 (the fnd3_contour_fractions style): the sampled path cross-checks
/// the exact closed forms of an annular cylinder — per-cell κ within the
/// declared ε_α, apertures within the 8×8 boundary-stratum tolerance
/// (one axis-aligned interface crosses at most one stratum row; its 8
/// jittered Bernoulli samples give σ ≈ √(p(1−p)/8)/8 ≈ 0.02 — 0.06 ≈ 3σ,
/// and the whole check is deterministic at the fixed seed anyway).
#[test]
fn revolved_profile_matches_exact_ring_forms() {
    let spec = VoxelSpec {
        r_min: 0.0,
        dr: 0.1,
        n_r: 8,
        z_min: 0.0,
        dz: 0.1,
        n_z: 10,
        n_theta: 1,
        seed: 42,
    };
    let w = voxelize(&annular_cylinder(), &spec).expect("voxelizes");
    let tol_ap = 0.06;
    let ov = |a0: f64, a1: f64, b0: f64, b1: f64| (a1.min(b1) - a0.max(b0)).max(0.0);
    for i_r in 0..spec.n_r {
        for i_z in 0..spec.n_z {
            let r0 = i_r as f64 * spec.dr;
            let r1 = r0 + spec.dr;
            let z0 = i_z as f64 * spec.dz;
            let z1 = z0 + spec.dz;
            let frac_z = ov(z0, z1, ZA, ZB) / spec.dz;
            // κ: solid occupies the r² band × z band of the cell.
            let r_lo = RA.max(r0);
            let r_hi = RB.min(r1);
            let frac_r2 = if r_hi > r_lo {
                (r_hi * r_hi - r_lo * r_lo) / (r1 * r1 - r0 * r0)
            } else {
                0.0
            };
            let kappa_exact = 1.0 - frac_r2 * frac_z;
            let cell = &w.cells[w.cell_index(i_r, 0, i_z)];
            assert!(
                (cell.kappa - kappa_exact).abs() <= w.eps_alpha.max(1e-15),
                "cell ({i_r},{i_z}) kappa {} vs exact {kappa_exact}",
                cell.kappa
            );
            // r-faces: fraction of the z span outside the solid band
            // (1 everywhere if the face radius is outside [RA, RB]).
            for (slot, rf) in [(0usize, r0), (1usize, r1)] {
                let exact = if rf > RA && rf < RB {
                    1.0 - frac_z
                } else {
                    1.0
                };
                assert!(
                    (cell.aperture[slot] - exact).abs() <= tol_ap,
                    "cell ({i_r},{i_z}) r-face {slot} {} vs exact {exact}",
                    cell.aperture[slot]
                );
            }
            // θ-faces: PLANAR (r, z) measure — linear r overlap. FaceDir
            // order (S10): θ± are slots 4, 5.
            let frac_r_lin = (r_hi - r_lo).max(0.0) / spec.dr;
            let th_exact = 1.0 - frac_r_lin * frac_z;
            for slot in [4usize, 5] {
                assert!(
                    (cell.aperture[slot] - th_exact).abs() <= tol_ap,
                    "cell ({i_r},{i_z}) θ-face {slot} {} vs exact {th_exact}",
                    cell.aperture[slot]
                );
            }
            // z-faces: annular (r²) measure at the face height. FaceDir order
            // (S10): z± are slots 2, 3.
            for (slot, zf) in [(2usize, z0), (3usize, z1)] {
                let exact = if zf > ZA && zf < ZB {
                    1.0 - frac_r2
                } else {
                    1.0
                };
                assert!(
                    (cell.aperture[slot] - exact).abs() <= tol_ap,
                    "cell ({i_r},{i_z}) z-face {slot} {} vs exact {exact}",
                    cell.aperture[slot]
                );
            }
        }
    }
}

/// §6.2: the winding number classifies where ray-parity fails. A closed
/// tetrahedron classifies exactly; removing one face (a hole) leaves the
/// winding classification of far-side interior points intact (w > 0.5)
/// while a +z ray through the hole counts zero crossings and flips the
/// parity verdict — the silent-catastrophe FND-3 §3.2 bans.
#[test]
fn winding_number_classifies_where_ray_parity_fails() {
    let closed = TriMesh {
        tris: tet_tris(false),
        dropped_degenerate: 0,
    };
    let holed = TriMesh {
        tris: tet_tris(true),
        dropped_degenerate: 0,
    };
    let p_in = [0.1, 0.1, 0.05]; // interior, near the corner far from the hole
    let p_out = [0.5, 0.5, 0.5]; // outside (x+y+z > 1)

    // Closed mesh: exact integers, both classifiers agree.
    assert!((closed.winding_number(p_in) - 1.0).abs() < 1e-9);
    assert!(closed.winding_number(p_out).abs() < 1e-9);
    assert!(z_ray_parity_inside(&closed.tris, p_in));

    // Holed mesh: the +z ray exits through the missing slant face —
    // parity misclassifies the interior point as outside…
    assert!(
        !z_ray_parity_inside(&holed.tris, p_in),
        "expected the demonstration ray to escape through the hole"
    );
    // …the winding number does not.
    let w = holed.winding_number(p_in);
    assert!(
        w > 0.5,
        "holed-mesh winding {w} should still classify inside"
    );
}

/// §6.4: the CSG↔STL cross-check — one part authored both ways agrees per
/// cell within sampling tolerance (identical jitter keys ⇒ identical
/// sample points; the two classifiers may only disagree on measure-zero
/// boundary hits, and the pure-cell provers differ — hence the 2ε_α + abs
/// floor rather than bitwise).
#[test]
fn csg_and_stl_agree_on_the_same_part() {
    let bmin = [-0.62, -0.57, 0.17];
    let bmax = [0.41, 0.33, 0.58];
    let spec = VoxelSpec {
        r_min: 0.0,
        dr: 0.1,
        n_r: 8,
        z_min: 0.0,
        dz: 0.1,
        n_z: 8,
        n_theta: 8,
        seed: 99,
    };
    let csg = Solid::Csg(Sdf::BoxAxis {
        min: bmin,
        max: bmax,
    });
    let stl = Solid::Mesh(box_tri_mesh(bmin, bmax));
    let wc = voxelize(&csg, &spec).expect("csg voxelizes");
    let ws = voxelize(&stl, &spec).expect("stl voxelizes");
    for (i, (a, b)) in wc.cells.iter().zip(&ws.cells).enumerate() {
        assert!(
            (a.kappa - b.kappa).abs() <= 2.0 * wc.eps_alpha + 1e-3,
            "cell {i}: csg kappa {} vs stl kappa {}",
            a.kappa,
            b.kappa
        );
    }
}

/// §6.3: a half-space cut (the huge-box face as the plane z = 0.437) —
/// per cut cell the PLIC normal is the plane normal, the volume-matched d
/// reproduces κ to 1e-12 (for n = ẑ, V(d)/V_box = (d + Δz/2)/Δz), and the
/// area equals the exact ring-sector slice ½(r_o²−r_i²)·Δθ.
#[test]
fn plic_recovers_a_planar_cut_exactly() {
    let z_cut = 0.437;
    let solid = Solid::Csg(Sdf::BoxAxis {
        min: [-5.0, -5.0, -5.0],
        max: [5.0, 5.0, z_cut],
    });
    let spec = VoxelSpec {
        r_min: 0.0,
        dr: 0.1,
        n_r: 4,
        z_min: 0.0,
        dz: 0.1,
        n_z: 8,
        n_theta: 4,
        seed: 3,
    };
    let w = voxelize(&solid, &spec).expect("voxelizes");
    let dth = TAU / f64::from(spec.n_theta);
    let mut cut_cells = 0usize;
    for i_r in 0..spec.n_r {
        for i_t in 0..spec.n_theta as usize {
            for i_z in 0..spec.n_z {
                let cell = &w.cells[w.cell_index(i_r, i_t, i_z)];
                let Some(p) = &cell.plic else { continue };
                cut_cells += 1;
                // Normal = +ẑ in the local frame (field increases upward).
                assert!(
                    (p.normal[0]).abs() < 1e-6
                        && (p.normal[1]).abs() < 1e-6
                        && (p.normal[2] - 1.0).abs() < 1e-6,
                    "normal {:?}",
                    p.normal
                );
                // d ⇒ κ consistency to the bisection floor.
                let solid_frac_from_d = (p.d + 0.5 * spec.dz) / spec.dz;
                assert!(
                    (solid_frac_from_d - (1.0 - cell.kappa)).abs() < 1e-12,
                    "d {} vs kappa {}",
                    p.d,
                    cell.kappa
                );
                // Exact slice area of the ring sector.
                let r0 = i_r as f64 * spec.dr;
                let r1 = r0 + spec.dr;
                let a_exact = 0.5 * (r1 * r1 - r0 * r0) * dth;
                assert!(
                    (p.area - a_exact).abs() < 1e-10 * a_exact,
                    "area {} vs exact {a_exact}",
                    p.area
                );
                // And d recovers the plane height within the sampled-κ error.
                let zc = (i_z as f64 + 0.5) * spec.dz;
                assert!(
                    (p.d - (z_cut - zc)).abs() <= w.eps_alpha * spec.dz + 1e-12,
                    "d {} vs plane offset {}",
                    p.d,
                    z_cut - zc
                );
            }
        }
    }
    assert!(cut_cells > 0, "the fixture must actually cut cells");
}

/// §6.3: Σ PLIC areas over the cut cells of an axis-centered sphere ≈
/// 4πR². Tolerance rationale: PLIC is first-order in interface curvature
/// (each plane is the tangent-plane surrogate over one cell, O(h/R) area
/// defect per cell) and the local box models the sector arc as flat
/// (O(Δθ²/24) ≈ 0.9% at N_θ = 8); a few % is the honest S9 gate.
#[test]
fn sphere_plic_area_sums_to_the_sphere() {
    let (zc, rad) = (0.8, 0.5);
    let solid = Solid::Csg(Sdf::Sphere {
        center: [0.0, 0.0, zc],
        radius: rad,
    });
    let spec = VoxelSpec {
        r_min: 0.0,
        dr: 0.025,
        n_r: 32,
        z_min: 0.0,
        dz: 0.025,
        n_z: 64,
        n_theta: 8,
        seed: 11,
    };
    let w = voxelize(&solid, &spec).expect("voxelizes");
    let a_sum: f64 = w
        .cells
        .iter()
        .filter_map(|c| c.plic.as_ref().map(|p| p.area))
        .sum();
    let a_exact = 4.0 * PI * rad * rad;
    let rel = (a_sum - a_exact).abs() / a_exact;
    println!("PLIC area sum {a_sum:.6} vs sphere {a_exact:.6} (rel {rel:.4})");
    assert!(
        rel < 0.05,
        "PLIC area sum {a_sum} vs 4πR² {a_exact} (rel {rel:.4})"
    );
}

/// §6.5 + §3.5: voxelization is a pure function of {geometry, resolution,
/// seed} — two identical calls are bit-identical over every emitted f64;
/// a different seed produces different bits (the seed is live).
#[test]
fn voxelization_is_deterministic() {
    let solid = Solid::Csg(Sdf::Sphere {
        center: [0.3, 0.0, 0.4],
        radius: 0.15,
    });
    let spec = VoxelSpec {
        r_min: 0.0,
        dr: 0.1,
        n_r: 8,
        z_min: 0.0,
        dz: 0.1,
        n_z: 8,
        n_theta: 8,
        seed: 1234,
    };
    let w1 = voxelize(&solid, &spec).expect("voxelizes");
    let w2 = voxelize(&solid, &spec).expect("voxelizes");
    assert_eq!(
        world_bits(&w1),
        world_bits(&w2),
        "identical inputs must be bit-identical"
    );
    let spec_b = VoxelSpec { seed: 1235, ..spec };
    let w3 = voxelize(&solid, &spec_b).expect("voxelizes");
    assert_ne!(
        world_bits(&w1),
        world_bits(&w3),
        "a different seed must move some bits"
    );
}

/// §6.6 + §3.4: an axisymmetric solid at N_θ = 8 yields per-θ-identical
/// fractions and apertures BITWISE (the θ-independent jitter entity keys
/// make θ-congruence exact, not statistical) and no azimuthal floor;
/// PLIC diagnostics, whose gradient probes are evaluated at rotated
/// Cartesian points, agree to round-off (1e-9) in the local frame. An
/// off-axis sphere earns the floor N_θ = 8 on the (r,z) rows it cuts.
#[test]
fn axisymmetric_solid_has_zero_theta_variance_and_no_floor() {
    let spec = VoxelSpec {
        r_min: 0.0,
        dr: 0.1,
        n_r: 8,
        z_min: 0.0,
        dz: 0.1,
        n_z: 10,
        n_theta: 8,
        seed: 42,
    };
    let w = voxelize(&annular_cylinder(), &spec).expect("voxelizes");
    let nt = spec.n_theta as usize;
    for i_r in 0..spec.n_r {
        for i_z in 0..spec.n_z {
            let base = &w.cells[w.cell_index(i_r, 0, i_z)];
            for i_t in 1..nt {
                let c = &w.cells[w.cell_index(i_r, i_t, i_z)];
                assert_eq!(
                    c.kappa.to_bits(),
                    base.kappa.to_bits(),
                    "κ bits differ at ({i_r},{i_t},{i_z})"
                );
                for f in 0..6 {
                    assert_eq!(
                        c.aperture[f].to_bits(),
                        base.aperture[f].to_bits(),
                        "aperture {f} bits differ at ({i_r},{i_t},{i_z})"
                    );
                }
                match (&base.plic, &c.plic) {
                    (None, None) => {}
                    (Some(a), Some(b)) => {
                        for k in 0..3 {
                            assert!(
                                (a.normal[k] - b.normal[k]).abs() < 1e-9,
                                "normal[{k}] {} vs {} at ({i_r},{i_t},{i_z})",
                                a.normal[k],
                                b.normal[k]
                            );
                        }
                        assert!(
                            (a.d - b.d).abs() < 1e-9,
                            "d {} vs {} at ({i_r},{i_t},{i_z}); n {:?} vs {:?}",
                            a.d,
                            b.d,
                            a.normal,
                            b.normal
                        );
                        assert!((a.area - b.area).abs() < 1e-9 * a.area.abs().max(1.0));
                    }
                    _ => panic!("PLIC presence differs across θ at ({i_r},{i_t},{i_z})"),
                }
            }
        }
    }
    let floors = theta_geom_floor(&w);
    assert!(
        floors.iter().all(|&f| f == 1),
        "axisymmetric solid must impose no floor"
    );

    // Off-axis sphere (S10 coarsest-reproducing form): the rows it cuts earn
    // a floor equal to the COARSEST N_θ whose θ-coarsening still reproduces
    // their fractions within ε_α — not the blanket N_θ the S9 binary form
    // pinned on any variation. A near-axis row where the sphere localizes
    // into ~one sector needs the full 8; a row it grazes may reproduce at a
    // coarser ladder rung. Every floor is a valid ladder value ≤ N_θ, a
    // θ-uniform row stays floor-free, and a strongly-cut row demands θ.
    let off = Solid::Csg(Sdf::Sphere {
        center: [0.3, 0.0, 0.4],
        radius: 0.15,
    });
    let w2 = voxelize(&off, &spec).expect("voxelizes");
    let floors2 = theta_geom_floor(&w2);
    let cut_row = floors2[2 * spec.n_z + 3];
    println!(
        "coarsest-reproducing floors (off-axis sphere): cut row (2,3) = {cut_row}, \
         distinct = {:?}",
        {
            let mut v: Vec<u32> = floors2.to_vec();
            v.sort_unstable();
            v.dedup();
            v
        }
    );
    assert!(
        cut_row >= 2,
        "a strongly-cut (r,z) row must floor above 1 (got {cut_row})"
    );
    assert_eq!(
        floors2[7 * spec.n_z + 9],
        1,
        "a far row must stay floor-free"
    );
    // Every floor is a power-of-two ladder value in [1, N_θ] (never 3, 5, …,
    // never above the sampled resolution).
    assert!(
        floors2
            .iter()
            .all(|&f| f >= 1 && f <= spec.n_theta && f.is_power_of_two()),
        "a floor escaped the θ-ladder: {floors2:?}"
    );
    // The coarsest-reproducing search is a refinement of (never coarser than
    // "no floor", never finer than) the S9 binary form: 1 ≤ f ≤ 8.
    assert!(
        floors2.iter().any(|&f| f > 1),
        "the sphere must floor SOMETHING"
    );
}

/// §3.2 import: one tetrahedron authored as binary STL bytes (with a
/// header that deliberately begins with "solid" — the record-count size
/// identity must win the disambiguation) and as ASCII text; both parse to
/// the same triangles, and the deliberately-degenerate triangle in each is
/// dropped WITH ITS COUNT RECORDED, never silently.
#[test]
fn stl_round_trip_parses_binary_and_ascii() {
    let tris = tet_tris(false);
    let degenerate = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];

    // Binary: 80-byte header + u32 count + 50-byte records.
    let mut bytes = Vec::new();
    let mut header = [b' '; 80];
    header[..17].copy_from_slice(b"solid binary-trap");
    bytes.extend_from_slice(&header);
    let all: Vec<[[f64; 3]; 3]> = tris.iter().copied().chain([degenerate]).collect();
    bytes.extend_from_slice(&(all.len() as u32).to_le_bytes());
    for t in &all {
        for _ in 0..3 {
            bytes.extend_from_slice(&0f32.to_le_bytes()); // normal (ignored)
        }
        for v in t {
            for c in v {
                bytes.extend_from_slice(&(*c as f32).to_le_bytes());
            }
        }
        bytes.extend_from_slice(&0u16.to_le_bytes());
    }
    let bin = TriMesh::from_stl_bytes(&bytes).expect("binary parses");
    assert_eq!(bin.tris.len(), 4);
    assert_eq!(bin.dropped_degenerate, 1);

    // ASCII.
    let mut text = String::from("solid tet\n");
    for t in &all {
        text.push_str("  facet normal 0 0 0\n    outer loop\n");
        for v in t {
            text.push_str(&format!("      vertex {} {} {}\n", v[0], v[1], v[2]));
        }
        text.push_str("    endloop\n  endfacet\n");
    }
    text.push_str("endsolid tet\n");
    let asc = TriMesh::from_stl_bytes(text.as_bytes()).expect("ascii parses");
    assert_eq!(asc.tris.len(), 4);
    assert_eq!(asc.dropped_degenerate, 1);
    assert_eq!(
        bin.tris, asc.tris,
        "binary and ASCII must yield identical triangles"
    );
}
