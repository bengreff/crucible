//! S13c residency cross-check: the class-A step on the REAL 3-D CUT GEOMETRY
//! (`cuda/residency_geometry.cu`) scored against the bit-exact CPU reference —
//! `Euler::eval_rhs`, `Euler::apply_srd`, `Sdc::step_flow` — on the ◆C3 world:
//! the geometry-of-record RL10 contour (`data/anchors/rl10_contour.csv`)
//! revolved at uniform N_θ = 8 through `Grid::build_with_geometry_theta` at the
//! ◆C3 desktop dial (5 cells across the throat), r_min = 0 (the axis parity
//! pair live), the contour's true wall normal for the slip ghosts, and the
//! cut cells' State Redistribution.
//!
//! PHASE 1: one class-A RHS on ALL active cells (every ghost rule is on the
//! device now, so nothing is excluded). PHASE 2: one SRD pass on a perturbed
//! state. PHASE 3: the resident marched SDC step with stagewise SRD vs
//! `step_flow`, all cells, + the FND-6 checkpoint split + rerun bit-identity.
//! PHASE 4: throughput of the 3-D rate on the resident world.
//!
//! Scope (Ben, session 28): interior + reflective cut walls + the Reflecting /
//! Transmissive domain BCs; GammaLaw EOS; no ledger reductions.
use crucible_engine::geometry::Contour;
use crucible_grid::{CellGeomTheta, FaceDir, Grid, GridSpec, Region};
use crucible_solvers::euler::{
    EULER_FIELDS, Euler, EulerFields, FlowBc, FlowBcs, GammaLaw, I_MZ, NCOMP, NPRIM,
    xcheck_srd_neighborhoods,
};
use crucible_solvers::sdc::{FlowClass, Sdc};

const GAMMA: f64 = 1.4;
const ECT: f64 = 1e-9;
/// Bound on the per-cell relative error, which the well-balanced cancellation
/// inflates (an FMA-order absolute difference over a near-zero net rate): a
/// sanity rail, not the ECT — the ECT is the component-scaled measure.
const CANCEL_BOUND: f64 = 1e-6;
/// ◆C3 desktop-spec dial (`configs/rl10_startup_3d.toml`).
const DIAL: f64 = 5.0;
const N_THETA: u32 = 8;

const GK_DOMAIN_REFLECT: i32 = 0;
const GK_DOMAIN_TRANSMISSIVE: i32 = 1;
const GK_WALL_SLIP: i32 = 3;
const GK_AXIS: i32 = 4;

#[repr(C)]
struct HostWorld {
    n_r: i32,
    n_z: i32,
    nt: i32,
    r_min: f64,
    dr: f64,
    z_min: f64,
    dz: f64,
    gamma: f64,
    has_geom: i32,
    act: *const i32,
    kappa: *const f64,
    ap: *const f64,
    rs_r: *const i32,
    rl_r: *const i32,
    klo_r: *const i32,
    khi_r: *const i32,
    nlo_r: *const f64,
    nhi_r: *const f64,
    rs_z: *const i32,
    rl_z: *const i32,
    klo_z: *const i32,
    khi_z: *const i32,
    nlo_z: *const f64,
    nhi_z: *const f64,
}
#[repr(C)]
struct HostSrd {
    ns: i32,
    na: i32,
    mem_off: *const i32,
    mem: *const i32,
    mem_kv: *const f64,
    cnt: *const i32,
    aff: *const i32,
    aff_owner: *const i32,
    inv_off: *const i32,
    inv: *const i32,
}

unsafe extern "C" {
    fn gpu_class_a_rhs_3d(cons: *const f64, hw: *const HostWorld, rate: *mut f64) -> i32;
    fn gpu_srd_3d(cons: *mut f64, hw: *const HostWorld, hs: *const HostSrd);
    fn gpu_class_a_march_3d(
        cons: *mut f64, hw: *const HostWorld, hs: *const HostSrd, dt: f64, nsteps: i32,
    ) -> i32;
    fn gpu_class_a_bench_3d(cons: *const f64, hw: *const HostWorld, iters: i32) -> f64;
    fn gpu_stable_dt_3d(cons: *const f64, hw: *const HostWorld, cfl: f64, bad: *mut i32) -> f64;
}

/// Smooth, subsonic, strictly-positive fixture with genuine θ-variation
/// (the m = 1 content the axis parity pair and the θ-sweep act on).
fn fixture(r: f64, th: f64, z: f64) -> [f64; NPRIM] {
    let rho = 1.0 + 0.20 * (5.7 * r + 2.9 * z).sin() * th.cos();
    let ur = 0.30 * (3.1 * z - 4.4 * r).sin() * (0.5 * th).cos();
    let ut = 0.25 * (2.8 * r + 1.5 * z).cos() + 0.1 * th.sin();
    let uz = 0.35 * (2.9 * r - 3.3 * z).sin();
    let p = 2.0 + 0.30 * (3.3 * z - 1.7 * r).cos() * (th + 0.3).cos();
    let c = 0.5 + 0.1 * (r * 3.0).sin();
    crucible_solvers::euler::prim6(rho, ur, ut, uz, p, c)
}
fn fixture_b(r: f64, th: f64, z: f64) -> [f64; NPRIM] {
    let a = fixture(r, th, z);
    let mut b = a;
    b[0] = a[0] * (1.0 + 0.15 * (7.0 * z + th).sin());
    b[4] = a[4] * (1.0 + 0.10 * (6.0 * r - th).cos());
    b
}

/// Host mirror of the host-side pieces of `HostWorld`/`HostSrd` (owning buffers).
struct Tables {
    act: Vec<i32>,
    kappa: Vec<f64>,
    ap: Vec<f64>,
    rs_r: Vec<i32>,
    rl_r: Vec<i32>,
    klo_r: Vec<i32>,
    khi_r: Vec<i32>,
    nlo_r: Vec<f64>,
    nhi_r: Vec<f64>,
    rs_z: Vec<i32>,
    rl_z: Vec<i32>,
    klo_z: Vec<i32>,
    khi_z: Vec<i32>,
    nlo_z: Vec<f64>,
    nhi_z: Vec<f64>,
    // SRD
    mem_off: Vec<i32>,
    mem: Vec<i32>,
    mem_kv: Vec<f64>,
    cnt: Vec<i32>,
    aff: Vec<i32>,
    aff_owner: Vec<i32>,
    inv_off: Vec<i32>,
    inv: Vec<i32>,
}

fn main() {
    // --- The geometry of record, revolved at N_θ = 8 (engine assembly's clip).
    let path = format!("{}/../../data/anchors/rl10_contour.csv", env!("CARGO_MANIFEST_DIR"));
    let content = std::fs::read_to_string(&path).expect("contour of record");
    let stations = crucible_config::parse_contour_csv(&content, crucible_config::INCH_M)
        .expect("contour parses");
    let contour = Contour::new(stations.clone(), 0.0).expect("contour");
    let r_throat = stations.iter().map(|&(_, r)| r).fold(f64::INFINITY, f64::min);
    let r_max = stations.iter().map(|&(_, r)| r).fold(0.0f64, f64::max);
    let z_min = stations[0].0;
    let span = stations[stations.len() - 1].0 - z_min;
    let dr = r_throat / DIAL;
    let n_r = (r_max / dr).ceil() as usize;
    let n_z = (span / dr).ceil() as usize;
    let dz = span / n_z as f64;
    let spec = GridSpec {
        r_min: 0.0,
        dr,
        n_r,
        z_min,
        dz,
        n_z,
        n_theta_max: N_THETA,
        axisymmetry_assertion: false,
    };
    let (r_min, z0) = (spec.r_min, spec.z_min);
    let c = contour.clone();
    let geom_rz = move |i_r: usize, i_z: usize| -> (Region, f64, [f64; 4]) {
        let r0 = r_min + i_r as f64 * dr;
        let r1 = r_min + (i_r + 1) as f64 * dr;
        let za = z0 + i_z as f64 * dz;
        let zb = z0 + (i_z + 1) as f64 * dz;
        let kappa = c.gas_volume_fraction(r0, r1, za, zb);
        if kappa > 0.0 {
            (
                Region::Gas,
                kappa,
                [
                    c.r_face_aperture(r0, za, zb),
                    c.r_face_aperture(r1, za, zb),
                    c.z_face_aperture(r0, r1, za),
                    c.z_face_aperture(r0, r1, zb),
                ],
            )
        } else {
            (Region::Exterior, 0.0, [0.0; 4])
        }
    };
    let mut g = Grid::build_with_geometry_theta(
        spec,
        EULER_FIELDS,
        |i_r, _j, i_z| {
            let (_, kappa, a) = geom_rz(i_r, i_z);
            CellGeomTheta { kappa, aperture: [a[0], a[1], a[2], a[3], kappa, kappa] }
        },
        |i_r, i_z| geom_rz(i_r, i_z).0,
    )
    .expect("◆C3 grid builds");
    let f = EulerFields::resolve(&g).expect("fields");
    let ids = f.ids();
    let eos = GammaLaw { gamma: GAMMA };
    crucible_solvers::euler::fill_from_prim(&mut g, &f, &eos, fixture);
    let nt = N_THETA as usize;
    let nrz = n_r * n_z;
    let n = nt * nrz;
    let n_active: usize = (0..nrz).filter(|&rz| g.is_active(rz / n_z, rz % n_z)).count();
    println!("S13c — class-A step on the REAL 3-D CUT GEOMETRY (◆C3 world) on-device");
    println!(
        "  RL10 contour of record, dial {DIAL}: {n_r}×{n_z} (r,z) × N_θ={nt} = {n} cells, \
         {} active gas cells, r_min = 0 (axis live)",
        n_active * nt
    );

    // --- The CPU operator (the ◆C3 run's own BC/wall settings, no inflow).
    let zero = |_: f64, _: f64, _: f64, _: f64| [0.0; NCOMP];
    let normal_fn = move |r: f64, z: f64| contour.wall_normal(r, z);
    let op = Euler {
        eos,
        source: &zero,
        bcs: FlowBcs {
            r_inner: FlowBc::Reflecting,
            r_outer: FlowBc::Reflecting,
            z_lo: FlowBc::Reflecting,
            z_hi: FlowBc::Transmissive,
        },
        wall_normal: Some(&normal_fn),
        slip_wall_z_faces: true,
        combustion: None,
    };
    let mut ws = op.workspace(&g).expect("workspace");

    // --- Host geometry-time tables (a pure function of the grid).
    let dense = |j: usize, i_r: usize, i_z: usize| j * nrz + i_r * n_z + i_z;
    let mut t = Tables {
        act: vec![0; nrz],
        kappa: vec![0.0; n],
        ap: vec![0.0; 6 * n],
        rs_r: vec![0; nrz],
        rl_r: vec![0; nrz],
        klo_r: vec![0; nrz],
        khi_r: vec![0; nrz],
        nlo_r: vec![0.0; 2 * nrz],
        nhi_r: vec![0.0; 2 * nrz],
        rs_z: vec![0; nrz],
        rl_z: vec![0; nrz],
        klo_z: vec![0; nrz],
        khi_z: vec![0; nrz],
        nlo_z: vec![0.0; 2 * nrz],
        nhi_z: vec![0.0; 2 * nrz],
        mem_off: vec![0],
        mem: vec![],
        mem_kv: vec![],
        cnt: vec![1; n],
        aff: vec![],
        aff_owner: vec![],
        inv_off: vec![0],
        inv: vec![],
    };
    let dirs = [
        FaceDir::RMinus,
        FaceDir::RPlus,
        FaceDir::ZMinus,
        FaceDir::ZPlus,
        FaceDir::ThetaMinus,
        FaceDir::ThetaPlus,
    ];
    for i_r in 0..n_r {
        for i_z in 0..n_z {
            let rz = i_r * n_z + i_z;
            if !g.is_active(i_r, i_z) {
                continue;
            }
            t.act[rz] = 1;
            for j in 0..nt {
                let cc = dense(j, i_r, i_z);
                t.kappa[cc] = g.kappa_at(i_r, j as u32, i_z);
                for (d, &dir) in dirs.iter().enumerate() {
                    t.ap[d * n + cc] = g.aperture_at(i_r, j as u32, i_z, dir);
                }
            }
        }
    }
    // Runs along r (per i_z line) — sweep_r's decomposition + ghost rules.
    let on_axis = r_min == 0.0;
    let bc_kind = |bc: &FlowBc<'_>| match bc {
        FlowBc::Reflecting => GK_DOMAIN_REFLECT,
        FlowBc::Transmissive => GK_DOMAIN_TRANSMISSIVE,
        _ => panic!("fixture BCs are Reflecting/Transmissive only"),
    };
    for i_z in 0..n_z {
        let z = g.z_center(i_z);
        let mut i = 0usize;
        while i < n_r {
            if t.act[i * n_z + i_z] == 0 {
                i += 1;
                continue;
            }
            let start = i;
            while i < n_r && t.act[i * n_z + i_z] == 1 {
                i += 1;
            }
            let len = i - start;
            let (klo, nlo) = if start == 0 && on_axis {
                (GK_AXIS, (0.0, 0.0))
            } else if start == 0 {
                (bc_kind(&op.bcs.r_inner), (0.0, 0.0))
            } else {
                // wall_ghosts_low(w, n, I_MR, r0 + start·dr, z): slip (r-faces always).
                (GK_WALL_SLIP, normal_fn(r_min + start as f64 * dr, z))
            };
            let (khi, nhi) = if start + len == n_r {
                (bc_kind(&op.bcs.r_outer), (0.0, 0.0))
            } else {
                (GK_WALL_SLIP, normal_fn(r_min + (start + len) as f64 * dr, z))
            };
            for ii in start..start + len {
                let rz = ii * n_z + i_z;
                t.rs_r[rz] = start as i32;
                t.rl_r[rz] = len as i32;
                t.klo_r[rz] = klo;
                t.khi_r[rz] = khi;
                t.nlo_r[2 * rz] = nlo.0;
                t.nlo_r[2 * rz + 1] = nlo.1;
                t.nhi_r[2 * rz] = nhi.0;
                t.nhi_r[2 * rz + 1] = nhi.1;
            }
        }
    }
    // Runs along z (per i_r line) — sweep_z's decomposition + ghost rules.
    let slip_z = op.slip_wall_z_faces; // wall_ghosts_*: slip iff normal != I_MZ || flag
    let _ = I_MZ;
    for i_r in 0..n_r {
        let r = g.r_center(i_r);
        let mut i = 0usize;
        while i < n_z {
            if t.act[i_r * n_z + i] == 0 {
                i += 1;
                continue;
            }
            let start = i;
            while i < n_z && t.act[i_r * n_z + i] == 1 {
                i += 1;
            }
            let len = i - start;
            let wall = |zf: f64| -> (i32, (f64, f64)) {
                if slip_z {
                    (GK_WALL_SLIP, normal_fn(r, zf))
                } else {
                    (2, (0.0, 0.0)) // GK_WALL_MIRROR
                }
            };
            let (klo, nlo) = if start == 0 {
                (bc_kind(&op.bcs.z_lo), (0.0, 0.0))
            } else {
                wall(z0 + start as f64 * dz)
            };
            let (khi, nhi) = if start + len == n_z {
                (bc_kind(&op.bcs.z_hi), (0.0, 0.0))
            } else {
                wall(z0 + (start + len) as f64 * dz)
            };
            for ii in start..start + len {
                let rz = i_r * n_z + ii;
                t.rs_z[rz] = start as i32;
                t.rl_z[rz] = len as i32;
                t.klo_z[rz] = klo;
                t.khi_z[rz] = khi;
                t.nlo_z[2 * rz] = nlo.0;
                t.nlo_z[2 * rz + 1] = nlo.1;
                t.nhi_z[2 * rz] = nhi.0;
                t.nhi_z[2 * rz + 1] = nhi.1;
            }
        }
    }
    // SRD tables in the CPU's own build order.
    let small = xcheck_srd_neighborhoods(&g).expect("SRD neighborhoods");
    let mut inv_lists: std::collections::BTreeMap<usize, (bool, Vec<i32>)> =
        std::collections::BTreeMap::new();
    for (s, (cells, kv)) in small.iter().enumerate() {
        for (m, &(r, z, j)) in cells.iter().enumerate() {
            let cc = dense(j as usize, r, z);
            t.mem.push(cc as i32);
            t.mem_kv.push(kv[m]);
            if m > 0 {
                t.cnt[cc] += 1;
            }
            let e = inv_lists.entry(cc).or_insert((false, vec![]));
            if m == 0 {
                e.0 = true;
            }
            e.1.push(s as i32);
        }
        t.mem_off.push(t.mem.len() as i32);
    }
    for (&cc, (owner, list)) in &inv_lists {
        t.aff.push(cc as i32);
        t.aff_owner.push(i32::from(*owner));
        t.inv.extend_from_slice(list);
        t.inv_off.push(t.inv.len() as i32);
    }
    let n_small = small.len();
    let hw = HostWorld {
        n_r: n_r as i32,
        n_z: n_z as i32,
        nt: nt as i32,
        r_min,
        dr,
        z_min: z0,
        dz,
        gamma: GAMMA,
        has_geom: i32::from(g.has_cut_geometry()),
        act: t.act.as_ptr(),
        kappa: t.kappa.as_ptr(),
        ap: t.ap.as_ptr(),
        rs_r: t.rs_r.as_ptr(),
        rl_r: t.rl_r.as_ptr(),
        klo_r: t.klo_r.as_ptr(),
        khi_r: t.khi_r.as_ptr(),
        nlo_r: t.nlo_r.as_ptr(),
        nhi_r: t.nhi_r.as_ptr(),
        rs_z: t.rs_z.as_ptr(),
        rl_z: t.rl_z.as_ptr(),
        klo_z: t.klo_z.as_ptr(),
        khi_z: t.khi_z.as_ptr(),
        nlo_z: t.nlo_z.as_ptr(),
        nhi_z: t.nhi_z.as_ptr(),
    };
    let hs = HostSrd {
        ns: n_small as i32,
        na: t.aff.len() as i32,
        mem_off: t.mem_off.as_ptr(),
        mem: t.mem.as_ptr(),
        mem_kv: t.mem_kv.as_ptr(),
        cnt: t.cnt.as_ptr(),
        aff: t.aff.as_ptr(),
        aff_owner: t.aff_owner.as_ptr(),
        inv_off: t.inv_off.as_ptr(),
        inv: t.inv.as_ptr(),
    };
    println!(
        "  SRD: {n_small} small sector-cells, {} affected cells, {} memberships",
        t.aff.len(),
        t.mem.len()
    );

    // Dense state <-> grid helpers.
    let read_dense = |g: &Grid| -> Vec<f64> {
        let mut v = vec![0.0f64; n * NCOMP];
        g.for_each_active_cell(|cell| {
            let cc = dense(cell.i_theta as usize, cell.i_r, cell.i_z);
            let b = g.brick(cell.bi);
            for k in 0..NCOMP {
                v[cc * NCOMP + k] = b.field(ids[k])[cell.idx];
            }
        });
        v
    };
    // Two measures per component k: the per-cell relative error (inflated
    // wherever the well-balanced cancellation — pressure flux vs the
    // geometric/wall-closure source, each O(10²) summing to ~0 — divides an
    // FMA-order absolute difference by a near-zero result), and the
    // COMPONENT-SCALED error |Δ| / max_cells |rate_k| — the ECT measure that
    // does not depend on where a cancellation lands. Returns (scaled, per-cell).
    let compare = |label: &str, cpu: &[f64], gpu: &[f64]| -> (f64, f64) {
        let mut comp_max = [0.0f64; NCOMP];
        for j in 0..nt {
            for rz in 0..nrz {
                if t.act[rz] == 0 {
                    continue;
                }
                let cc = j * nrz + rz;
                for k in 0..NCOMP {
                    comp_max[k] = comp_max[k].max(cpu[cc * NCOMP + k].abs());
                }
            }
        }
        let mut worst_abs = 0.0f64;
        let mut worst_rel = 0.0f64;
        let mut worst_scaled = 0.0f64;
        let mut n_cmp = 0usize;
        let mut worst_at = (0usize, 0usize, 0.0f64);
        for j in 0..nt {
            for rz in 0..nrz {
                if t.act[rz] == 0 {
                    continue;
                }
                let cc = j * nrz + rz;
                for k in 0..NCOMP {
                    let a = cpu[cc * NCOMP + k];
                    let b = gpu[cc * NCOMP + k];
                    let d = (a - b).abs();
                    worst_abs = worst_abs.max(d);
                    if comp_max[k] > 0.0 {
                        worst_scaled = worst_scaled.max(d / comp_max[k]);
                    }
                    let s = a.abs().max(b.abs());
                    if s > 1e-6 && d / s > worst_rel {
                        worst_rel = d / s;
                        worst_at = (cc, k, a);
                    }
                    n_cmp += 1;
                }
            }
        }
        let (cc, k, a) = worst_at;
        println!(
            "  {label}: {n_cmp} scalar comparisons (all active cells); worst abs {worst_abs:.3e}; \
             worst component-scaled rel {worst_scaled:.3e}; worst per-cell rel {worst_rel:.3e} \
             at (j={}, i_r={}, i_z={}, k={k}: value {a:.3e} vs component max {:.3e})",
            cc / nrz,
            (cc % nrz) / n_z,
            cc % n_z,
            comp_max[k]
        );
        (worst_scaled, worst_rel)
    };

    // ================================================================
    // PHASE 1 — one class-A RHS on the cut 3-D world, ALL active cells.
    // ================================================================
    op.eval_rhs(&g, &f, &mut ws, 0.0).expect("cpu eval_rhs");
    let cons0 = read_dense(&g);
    let mut cpu_rate = vec![0.0f64; n * NCOMP];
    let rates = ws.rates();
    g.for_each_active_cell(|cell| {
        let cc = dense(cell.i_theta as usize, cell.i_r, cell.i_z);
        for k in 0..NCOMP {
            cpu_rate[cc * NCOMP + k] = rates[cell.bi][cell.idx][k];
        }
    });
    let mut gpu_rate = vec![0.0f64; n * NCOMP];
    let bad = unsafe { gpu_class_a_rhs_3d(cons0.as_ptr(), &hw, gpu_rate.as_mut_ptr()) };
    assert_eq!(bad, 0, "GPU flagged a non-physical cell");
    let mut gpu_rate2 = vec![0.0f64; n * NCOMP];
    unsafe { gpu_class_a_rhs_3d(cons0.as_ptr(), &hw, gpu_rate2.as_mut_ptr()) };
    println!("PHASE 1 — class-A RHS: r/θ/z sweeps + axis parity + wall ghosts + closure");
    let (w1, c1) = compare("rhs", &cpu_rate, &gpu_rate);
    let bit1 = gpu_rate == gpu_rate2;
    println!("  GPU same-build rerun: {}", if bit1 { "BIT-IDENTICAL" } else { "*** DIFFERS ***" });
    assert!(w1 < ECT, "CPU↔GPU 3-D class-A RHS diverged beyond ECT (component-scaled): {w1:.3e}");
    assert!(c1 < CANCEL_BOUND, "per-cell rel {c1:.3e} beyond the cancellation-inflated bound");
    assert!(bit1, "GPU rerun not deterministic");
    println!("  PHASE 1 PASS");

    // ================================================================
    // PHASE 1b — stable_dt with the θ-arc member (the 3-D CFL clock).
    // ================================================================
    let cfl = 0.4;
    let cpu_dt = op.stable_dt(&g, &f, cfl).expect("cpu stable_dt (θ-CFL live)");
    let mut sb = 0i32;
    let gpu_dt = unsafe { gpu_stable_dt_3d(cons0.as_ptr(), &hw, cfl, &mut sb) };
    assert_eq!(sb, 0, "GPU stable_dt flagged a non-physical cell");
    let mut sb2 = 0i32;
    let gpu_dt2 = unsafe { gpu_stable_dt_3d(cons0.as_ptr(), &hw, cfl, &mut sb2) };
    let sdt_rel = (cpu_dt - gpu_dt).abs() / cpu_dt.abs().max(gpu_dt.abs());
    println!("PHASE 1b — stable_dt with the θ-arc member, cfl = {cfl}");
    println!("  cpu Δt = {cpu_dt:.9e}   gpu Δt = {gpu_dt:.9e}   rel = {sdt_rel:.3e}");
    println!(
        "  GPU same-build rerun: {}",
        if gpu_dt == gpu_dt2 { "BIT-IDENTICAL" } else { "*** DIFFERS ***" }
    );
    assert!(sdt_rel < ECT, "CPU↔GPU 3-D stable_dt diverged beyond ECT");
    assert!(gpu_dt == gpu_dt2, "GPU stable_dt rerun not deterministic");
    println!("  PHASE 1b PASS");

    // ================================================================
    // PHASE 2 — one State-Redistribution pass on a perturbed state.
    // ================================================================
    let mut g2 = g.clone();
    crucible_solvers::euler::fill_from_prim(&mut g2, &f, &eos, fixture_b);
    let pre = read_dense(&g2);
    op.apply_srd(&mut g2, &f, &ws);
    let cpu_srd = read_dense(&g2);
    let mut gpu_srd = pre.clone();
    unsafe { gpu_srd_3d(gpu_srd.as_mut_ptr(), &hw, &hs) };
    let mut gpu_srd2 = pre.clone();
    unsafe { gpu_srd_3d(gpu_srd2.as_mut_ptr(), &hw, &hs) };
    // The pass must actually move something (the fixture is non-uniform).
    let moved = pre.iter().zip(&cpu_srd).filter(|(a, b)| a != b).count();
    println!("PHASE 2 — State Redistribution pass ({moved} scalars moved on the CPU)");
    let (w2, c2) = compare("srd", &cpu_srd, &gpu_srd);
    let bit2 = gpu_srd == gpu_srd2;
    println!("  GPU same-build rerun: {}", if bit2 { "BIT-IDENTICAL" } else { "*** DIFFERS ***" });
    assert!(moved > 0, "SRD pass moved nothing — the fixture is not exercising small cells");
    assert!(w2 < ECT, "CPU↔GPU SRD diverged beyond ECT: {w2:.3e}");
    assert!(c2 < CANCEL_BOUND, "SRD per-cell rel {c2:.3e} beyond the cancellation-inflated bound");
    assert!(bit2, "GPU SRD rerun not deterministic");
    println!("  PHASE 2 PASS");

    // ================================================================
    // PHASE 3 — the resident marched SDC step with stagewise SRD.
    // ================================================================
    let m = 5i32;
    let dt = op.stable_dt(&g, &f, 0.4).expect("cpu stable_dt (θ-CFL live)");
    let mut sdc = Sdc::new();
    let flow = FlowClass { op: &op, fields: &f };
    let mut tt = 0.0;
    for _ in 0..m {
        sdc.step_flow(&mut g, &flow, tt, dt).expect("cpu step_flow");
        tt += dt;
    }
    let cpu_final = read_dense(&g);
    let mut gpu_final = cons0.clone();
    let bad = unsafe { gpu_class_a_march_3d(gpu_final.as_mut_ptr(), &hw, &hs, dt, m) };
    assert_eq!(bad, 0, "GPU march flagged a non-physical cell");
    let mut split = cons0.clone();
    unsafe {
        gpu_class_a_march_3d(split.as_mut_ptr(), &hw, &hs, dt, 2);
        gpu_class_a_march_3d(split.as_mut_ptr(), &hw, &hs, dt, 3);
    }
    let mut gpu_final2 = cons0.clone();
    unsafe { gpu_class_a_march_3d(gpu_final2.as_mut_ptr(), &hw, &hs, dt, m) };
    println!("PHASE 3 — resident marched SDC step + stagewise SRD ({m} steps, dt = {dt:.4e} s)");
    let (w3, c3) = compare("march", &cpu_final, &gpu_final);
    let bit3 = gpu_final == gpu_final2;
    let ck3 = split == gpu_final;
    println!("  GPU resident-march rerun: {}", if bit3 { "BIT-IDENTICAL" } else { "*** DIFFERS ***" });
    println!(
        "  checkpoint march(2)+march(3) == march(5): {}",
        if ck3 { "BIT-IDENTICAL" } else { "*** DIFFERS ***" }
    );
    assert!(w3 < 1e-8, "marched CPU↔GPU diverged beyond the FMA-order band: {w3:.3e}");
    assert!(c3 < CANCEL_BOUND, "march per-cell rel {c3:.3e} beyond the cancellation-inflated bound");
    assert!(bit3, "GPU resident march not deterministic");
    assert!(ck3, "checkpoint/restart not bit-faithful");
    println!("  PHASE 3 PASS");

    // ================================================================
    // PHASE 4 — throughput of the 3-D rate on the resident ◆C3 world.
    // ================================================================
    let iters = 200i32;
    let ms = unsafe { gpu_class_a_bench_3d(cons0.as_ptr(), &hw, iters) };
    let evals = iters as f64 * (n_active * nt) as f64;
    println!("PHASE 4 — 3-D class-A rate throughput ({} active cells, {iters} evals)", n_active * nt);
    println!(
        "  {ms:.2} ms total → {:.3e} cell-RHS/s (3 sweeps + sources per eval)",
        evals / (ms / 1000.0)
    );
    println!("ALL PASS");
}
