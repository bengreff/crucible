//! SOLV-1 §3.6 cut-cell tests: the aperture-weighted operator + embedded
//! wall closure + State Redistribution (META-3 `state-redistribution`) on
//! a revolved LINEAR wall (cone) — exact test-local fractions.
//!
//! 1. Well-balance: a uniform state at rest in the cut cone is preserved
//!    to round-off (the closure vector is DEFINED by this identity; the
//!    test guards the implementation).
//! 2. Conservation: a closed cut domain under a strong internal transient
//!    conserves mass/energy/composition to ~1e-12 with SRD active — SRD
//!    redistributes, never creates.
//! 3. Stability: cells with κ ~ 1e-3 march at the UNCUT CFL without
//!    density runaway (the session-11 exit-lip pathology class).
//! 4. An enclosed sliver with no flow-connected neighborhood refuses.

use crucible_grid::{CellGeom, Grid, GridSpec, Region};
use crucible_solvers::euler::{
    Cons, Euler, EulerFields, FlowBc, FlowBcs, GammaLaw, NCOMP, fill_from_prim, prim6,
};

const GAMMA: f64 = 1.4;

/// Exact fractions for a linear wall w(z) = C0 + C1·z over one cell
/// (single segment — the test-local oracle; the production contour path
/// has its own analytic tests in the engine crate).
fn linear_wall_geom(c0: f64, c1: f64, r0: f64, r1: f64, z0: f64, z1: f64) -> (f64, [f64; 4]) {
    let w = |z: f64| c0 + c1 * z;
    // κ: integrate ½(clamp(w,r0,r1)² − r0²) dz piecewise.
    let mut pts = vec![z0];
    for rc in [r0, r1] {
        if (w(z0) - rc) * (w(z1) - rc) < 0.0 {
            pts.push(z0 + (z1 - z0) * (rc - w(z0)) / (w(z1) - w(z0)));
        }
    }
    pts.push(z1);
    pts.sort_by(|a, b| a.partial_cmp(b).expect("finite"));
    let mut acc = 0.0f64;
    for p in pts.windows(2) {
        let (pa, pb) = (p[0], p[1]);
        let len = pb - pa;
        if len <= 0.0 {
            continue;
        }
        let wm = w(0.5 * (pa + pb));
        if wm <= r0 {
        } else if wm >= r1 {
            acc += 0.5 * (r1 * r1 - r0 * r0) * len;
        } else {
            let (va, vb) = (w(pa), w(pb));
            acc += 0.5 * (len * (va * va + va * vb + vb * vb) / 3.0 - r0 * r0 * len);
        }
    }
    let kappa = (acc / (0.5 * (r1 * r1 - r0 * r0) * (z1 - z0))).clamp(0.0, 1.0);
    let ap_r = |rf: f64| -> f64 {
        let (wa, wb) = (w(z0), w(z1));
        if wa > rf && wb > rf {
            1.0
        } else if wa <= rf && wb <= rf {
            0.0
        } else {
            let t = (rf - wa) / (wb - wa);
            if wa > rf { t } else { 1.0 - t }
        }
    };
    let ap_z = |zf: f64| -> f64 {
        let x = w(zf).clamp(r0, r1);
        ((x * x - r0 * r0) / (r1 * r1 - r0 * r0)).clamp(0.0, 1.0)
    };
    (kappa, [ap_r(r0), ap_r(r1), ap_z(z0), ap_z(z1)])
}

/// Cut-cone world: wall w(z) = c0 + c1·z on a (n_r × n_z) grid with
/// dr = dz = h, r_min = 0 (axis).
fn cone_grid(c0: f64, c1: f64, h: f64, n_r: usize, n_z: usize) -> Grid {
    let spec = GridSpec {
        r_min: 0.0,
        dr: h,
        n_r,
        z_min: 0.0,
        dz: h,
        n_z,
        n_theta_max: 1,
        axisymmetry_assertion: true,
    };
    let names = ["rho", "mom_r", "mom_theta", "mom_z", "rho_e", "rho_c"];
    Grid::build_with_geometry(spec, &names, |i_r, i_z| {
        let (r0, r1) = (i_r as f64 * h, (i_r + 1) as f64 * h);
        let (z0, z1) = (i_z as f64 * h, (i_z + 1) as f64 * h);
        let (kappa, aperture) = linear_wall_geom(c0, c1, r0, r1, z0, z1);
        if kappa > 0.0 {
            CellGeom {
                region: Region::Gas,
                kappa,
                aperture,
            }
        } else {
            CellGeom {
                region: Region::Exterior,
                kappa: 0.0,
                aperture: [0.0; 4],
            }
        }
    })
    .expect("cone world builds")
}

fn closed_cone_op<'a>(src: &'a (dyn Fn(f64, f64, f64, f64) -> Cons + Sync)) -> Euler<'a, GammaLaw> {
    Euler {
        eos: GammaLaw { gamma: GAMMA },
        source: src,
        bcs: FlowBcs {
            r_inner: FlowBc::Reflecting,
            r_outer: FlowBc::Reflecting,
            z_lo: FlowBc::Reflecting,
            z_hi: FlowBc::Reflecting,
        },
        wall_normal: None,
        slip_wall_z_faces: false,
    }
}

fn gas_mass_energy(g: &Grid, f: &EulerFields) -> (f64, f64, f64) {
    let ids = f.ids();
    let (mut m, mut e, mut c) = (0.0f64, 0.0f64, 0.0f64);
    g.for_each_active_cell(|cell| {
        let kv = g.kappa(cell.i_r, cell.i_z) * g.cell_volume(cell.i_r, 1);
        let b = g.brick(cell.bi);
        m += kv * b.field(ids[0])[cell.idx];
        e += kv * b.field(ids[4])[cell.idx];
        c += kv * b.field(ids[5])[cell.idx];
    });
    (m, e, c)
}

#[test]
fn uniform_state_at_rest_is_preserved_in_the_cut_cone() {
    // Wall from 0.62·h·n to steeper — cuts cells at varying κ, both r- and
    // z-faces partially open (c1 = 0.35 slope).
    let h = 0.1;
    let g0 = cone_grid(0.62, 0.35, h, 12, 16);
    let zero: fn(f64, f64, f64, f64) -> Cons = |_, _, _, _| [0.0; NCOMP];
    let op = closed_cone_op(&zero);
    let mut g = g0;
    let f = EulerFields::resolve(&g).expect("fields");
    let eos = GammaLaw { gamma: GAMMA };
    fill_from_prim(&mut g, &f, &eos, |_, _, _| {
        prim6(1.3, 0.0, 0.0, 0.0, 2.7e5, 0.4)
    });
    let dt = op.stable_dt(&g, &f, 0.4).expect("dt");
    op.advance(&mut g, &f, 0.0, dt, 25).expect("marches");
    let ids = f.ids();
    let mut worst = 0.0f64;
    g.for_each_active_cell(|cell| {
        let b = g.brick(cell.bi);
        let rho = b.field(ids[0])[cell.idx];
        let mr = b.field(ids[1])[cell.idx];
        let mz = b.field(ids[3])[cell.idx];
        worst = worst
            .max(((rho - 1.3) / 1.3).abs())
            .max(mr.abs() / 1.3)
            .max(mz.abs() / 1.3);
    });
    assert!(
        worst < 1e-11,
        "uniform rest state drifted by {worst:.3e} in the cut cone"
    );
}

#[test]
fn closed_cut_domain_conserves_through_a_shock_transient() {
    let h = 0.1;
    let g0 = cone_grid(0.62, 0.35, h, 12, 16);
    let zero: fn(f64, f64, f64, f64) -> Cons = |_, _, _, _| [0.0; NCOMP];
    let op = closed_cone_op(&zero);
    let mut g = g0;
    let f = EulerFields::resolve(&g).expect("fields");
    let eos = GammaLaw { gamma: GAMMA };
    // A strong off-center pressure bump drives shocks into the cut wall.
    fill_from_prim(&mut g, &f, &eos, |r, _, z| {
        let d2 = (r - 0.2) * (r - 0.2) + (z - 0.5) * (z - 0.5);
        let p = 1.0e5 * (1.0 + 9.0 * (-d2 / 0.01).exp());
        prim6(1.0, 0.0, 0.0, 0.0, p, 0.3)
    });
    let (m0, e0, c0) = gas_mass_energy(&g, &f);
    for _ in 0..60 {
        let dt = op.stable_dt(&g, &f, 0.4).expect("dt");
        op.step(&mut g, &f, 0.0, dt).expect("steps");
    }
    let (m1, e1, c1) = gas_mass_energy(&g, &f);
    assert!(
        ((m1 - m0) / m0).abs() < 1e-12,
        "mass drift {:.3e}",
        (m1 - m0) / m0
    );
    assert!(
        ((e1 - e0) / e0).abs() < 1e-12,
        "energy drift {:.3e}",
        (e1 - e0) / e0
    );
    assert!(
        ((c1 - c0) / c0).abs() < 1e-12,
        "composition drift {:.3e}",
        (c1 - c0) / c0
    );
}

#[test]
fn sliver_cells_march_at_the_uncut_cfl_without_runaway() {
    // Wall grazing a face: w = 0.60002 + 0.01·z over a [0.6, 0.7] ring —
    // κ down to ~4e-3 at z = 0. The session-11 pathology died here; SRD
    // must march it at the uncut CFL.
    let h = 0.1;
    let g0 = cone_grid(0.60002, 0.01, h, 8, 16);
    let mut min_kappa = f64::INFINITY;
    for i_z in 0..16 {
        let k = g0.kappa(6, i_z);
        if k > 0.0 {
            min_kappa = min_kappa.min(k);
        }
    }
    assert!(
        min_kappa < 5e-3,
        "fixture must contain sliver cells (min κ {min_kappa:.2e})"
    );
    let zero: fn(f64, f64, f64, f64) -> Cons = |_, _, _, _| [0.0; NCOMP];
    let op = closed_cone_op(&zero);
    let mut g = g0;
    let f = EulerFields::resolve(&g).expect("fields");
    let eos = GammaLaw { gamma: GAMMA };
    fill_from_prim(&mut g, &f, &eos, |_, _, z| {
        // An axial pressure gradient drives flow along (and into) the wall.
        prim6(1.0, 0.0, 0.0, 0.0, 1.0e5 * (1.0 + 0.5 * (1.6 - z)), 0.3)
    });
    let dt0 = op.stable_dt(&g, &f, 0.4).expect("dt");
    for _ in 0..80 {
        let dt = op.stable_dt(&g, &f, 0.4).expect("dt stays evaluable");
        assert!(
            dt > 0.25 * dt0,
            "CFL collapsed: dt {dt:.3e} vs initial {dt0:.3e} — small-cell runaway"
        );
        op.step(&mut g, &f, 0.0, dt).expect("marches");
    }
    // Densities stay physical everywhere, slivers included.
    let ids = f.ids();
    g.for_each_active_cell(|cell| {
        let rho = g.brick(cell.bi).field(ids[0])[cell.idx];
        assert!(
            rho.is_finite() && rho > 1e-3,
            "cell ({}, {}) density {rho:.3e} ran away",
            cell.i_r,
            cell.i_z
        );
    });
}

#[test]
fn an_enclosed_sliver_refuses_loudly() {
    // A κ = 0.1 pocket with all faces sealed (legal geometry input — a
    // feature thinner than a cell) has no merge neighborhood: refuse.
    let spec = GridSpec {
        r_min: 0.0,
        dr: 0.1,
        n_r: 4,
        z_min: 0.0,
        dz: 0.1,
        n_z: 4,
        n_theta_max: 1,
        axisymmetry_assertion: true,
    };
    let names = ["rho", "mom_r", "mom_theta", "mom_z", "rho_e", "rho_c"];
    let g = Grid::build_with_geometry(spec, &names, |i_r, i_z| {
        if (i_r, i_z) == (2, 2) {
            CellGeom {
                region: Region::Gas,
                kappa: 0.1,
                aperture: [0.0; 4],
            }
        } else {
            CellGeom {
                region: Region::Exterior,
                kappa: 0.0,
                aperture: [0.0; 4],
            }
        }
    })
    .expect("a sealed pocket is legal geometry");
    let mut g = g;
    let f = EulerFields::resolve(&g).expect("fields");
    let eos = GammaLaw { gamma: GAMMA };
    fill_from_prim(&mut g, &f, &eos, |_, _, _| {
        prim6(1.0, 0.0, 0.0, 0.0, 1.0e5, 0.0)
    });
    let zero: fn(f64, f64, f64, f64) -> Cons = |_, _, _, _| [0.0; NCOMP];
    let op = closed_cone_op(&zero);
    let err = op.step(&mut g, &f, 0.0, 1e-6).expect_err("must refuse");
    assert!(
        format!("{err}").contains("no ") || format!("{err}").contains("neighborhood"),
        "wrong refusal: {err}"
    );
}

#[test]
fn parallel_march_is_bit_identical_at_any_thread_count() {
    // FND-2 §3.7: the parallel path partitions work by data ownership
    // (brick rows/columns), so every cell's accumulation order — and hence
    // every bit — is thread-count-independent. March the cut cone on 1 vs
    // 4 rayon threads and compare every field bitwise.
    let run = |threads: usize| -> Vec<u64> {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .expect("pool");
        pool.install(|| {
            let h = 0.1;
            let mut g = cone_grid(0.62, 0.35, h, 12, 16);
            let f = EulerFields::resolve(&g).expect("fields");
            let eos = GammaLaw { gamma: GAMMA };
            fill_from_prim(&mut g, &f, &eos, |r, _, z| {
                let d2 = (r - 0.2) * (r - 0.2) + (z - 0.5) * (z - 0.5);
                prim6(
                    1.0,
                    0.0,
                    0.0,
                    0.0,
                    1.0e5 * (1.0 + 4.0 * (-d2 / 0.02).exp()),
                    0.3,
                )
            });
            let zero: fn(f64, f64, f64, f64) -> Cons = |_, _, _, _| [0.0; NCOMP];
            let op = closed_cone_op(&zero);
            let mut ws = op.workspace(&g).expect("workspace");
            for _ in 0..25 {
                let dt = op.stable_dt_ws(&g, &f, &ws, 0.4).expect("dt");
                op.step_ws(&mut g, &f, &mut ws, 0.0, dt).expect("steps");
            }
            let ids = f.ids();
            let mut bits = Vec::new();
            g.for_each_active_cell(|c| {
                for id in ids {
                    bits.push(g.brick(c.bi).field(id)[c.idx].to_bits());
                }
            });
            bits
        })
    };
    let one = run(1);
    let four = run(4);
    assert_eq!(one.len(), four.len());
    let diff = one.iter().zip(&four).filter(|(a, b)| a != b).count();
    assert_eq!(
        diff, 0,
        "{diff} field values differ between 1 and 4 threads"
    );
}
