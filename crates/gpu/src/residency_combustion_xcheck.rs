//! S13c residency cross-check: the COMBUSTION SOURCE (SOLV-4.4 bistable-Nagumo
//! propagation + matched front-thickening diffusion) ported to device
//! (`cuda/residency_combustion.cu`: `k_comb_source`) scored against the
//! bit-exact CPU reference `Combustion::accumulate_inner` (via the additive
//! doc-hidden `xcheck_source_dense` accessor, which runs the IDENTICAL code the
//! production `accumulate` runs) on the PRODUCTION LOX/LH₂ combustion surfaces.
//!
//! This is the flame-propagation physics — the ◆C3 ignition headline: the
//! pushed front ρ_u·K·b(1−b)(b−a) + ∇·(ρD_c∇b), reading the unburnt (p,h,Z)
//! T_u/ρ_u surfaces + the ignition (p,T_u,Z) S_L surface. A per-cell gather.
//!
//! SCOPE: the class-A propagation source; the class-R implicit auto-ignition
//! node solve (spontaneous light) is the remaining combustion leg.
use crucible_grid::{Grid, GridSpec};
use crucible_solvers::euler::{
    BurnBlendEos, CombMarshal, Combustion, EULER_FIELDS, EulerFields, I_RB, I_RHO, IgnitionColumns,
    NPRIM, TableEos, THETA_CELLS,
};
use crucible_tables::{ColumnMarshal, Pin, Table};

const NCOMP: usize = 7;
const ECT: f64 = 1e-9;

// Design flame state (solv4_combustion.rs): flammable unburnt at 1 atm.
const P0: f64 = 1.0e5;
const H0: f64 = 4.364e4;
const Z0: f64 = 0.167;

unsafe extern "C" {
    #[allow(clippy::too_many_arguments)]
    fn gpu_combustion_source(
        rho: *const f64, p: *const f64, z: *const f64, b: *const f64, e: *const f64,
        gas: *const f64, n_r: i32, n_z: i32, r_min: f64, dr: f64, dz: f64,
        up: *const f64, unp: i32, ulp: i32, uh: *const f64, unh: i32, ulh: i32,
        uz: *const f64, unz: i32, ulz: i32, ustr: *const i32,
        temp_data: *const f64, temp_vlog: i32, rhou_data: *const f64, rhou_vlog: i32,
        ip: *const f64, inp: i32, ilp: i32, it: *const f64, inh: i32, ilh: i32,
        iz: *const f64, inz: i32, ilz: i32, istr: *const i32, flame_data: *const f64, flame_vlog: i32,
        hu_floor: f64, hu_ceil: f64, p_floor: f64, tu_floor: f64, wrinkling: f64, theta: f64,
        rate: *mut f64,
    );
}

fn open(file: &str, group: &str, pins_toml: &str) -> Table {
    let path = format!("{}/../../tables/chem/{file}", env!("CARGO_MANIFEST_DIR"));
    let doc: toml::Table = pins_toml.parse().expect("pins parse");
    let entry = doc[group].as_table().expect("group entry");
    let pin = Pin {
        data_version: entry["data_version"].as_str().expect("ver").to_string(),
        content_digest: entry.get("content_digest").and_then(|v| v.as_str()).map(str::to_string),
    };
    Table::open(&path, group, &pin).expect("table loads under its pin")
}

fn main() {
    let ut = open(
        "lox_lh2_unburnt_v0.3.0.h5", "/chem/lox_lh2/unburnt",
        include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tables/chem/lox_lh2_unburnt_v0.3.0.pins.toml")),
    );
    let bt = open(
        "lox_lh2_v0.4.0.h5", "/chem/lox_lh2/equilibrium",
        include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tables/chem/lox_lh2_v0.4.0.pins.toml")),
    );
    let it = open(
        "lox_lh2_ignition_v0.3.0.h5", "/chem/lox_lh2/ignition",
        include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tables/chem/lox_lh2_ignition_v0.3.0.pins.toml")),
    );
    let blend = BurnBlendEos::new(
        TableEos::bind(&ut).expect("unburnt binds"),
        TableEos::bind(&bt).expect("burnt binds"),
    );
    let ign = IgnitionColumns::bind(&it).expect("ignition binds");
    let comb = Combustion { blend: &blend, ignition: ign, wrinkling: 1.0, theta: THETA_CELLS };
    let m: CombMarshal = comb.xcheck_marshal();

    // N_θ=1 box; b varies in both r and z (exercise both diffusion directions);
    // p, z, h uniform at the design flame state (T_u flammable, reaction active
    // where b ∈ (a, 1−BURN_COMPLETE)).
    let (n_r, n_z) = (32usize, 48usize);
    let (r_min, dr, dz) = (0.5, 1.0 / n_r as f64, 1.0 / n_z as f64);
    let spec = GridSpec { r_min, dr, n_r, z_min: 0.0, dz, n_z, n_theta_max: 1, axisymmetry_assertion: true };
    let mut g = Grid::build(spec, EULER_FIELDS).expect("grid");
    let f = EulerFields::resolve(&g).expect("fields");
    let ids = f.ids();

    let rho0 = 0.08f64; // ~H2/O2 at 1 atm, 312 K, MR5 (any positive value is valid)
    let b_of = |i_r: usize, i_z: usize| {
        let s = 0.5
            + 0.42 * (2.3 * i_r as f64 / n_r as f64 - 0.5).sin()
            * (3.1 * i_z as f64 / n_z as f64).cos();
        s.clamp(0.05, 0.95)
    };
    let e0 = H0 - P0 / rho0; // so h = e + p/ρ = H0

    // Fill the grid ρ and ρb fields (b_at reads them for the neighbour gradient).
    let cells: Vec<(usize, usize, usize, usize)> = {
        let mut v = vec![];
        g.for_each_active_cell(|c| v.push((c.bi, c.idx, c.i_r, c.i_z)));
        v
    };
    for &(bi, idx, i_r, i_z) in &cells {
        let b = b_of(i_r, i_z);
        g.brick_field_mut(bi, ids[I_RHO])[idx] = rho0;
        g.brick_field_mut(bi, ids[I_RB])[idx] = rho0 * b;
    }

    // CPU oracle: the REAL combustion source via accumulate_inner.
    let prim_of = |i_r: usize, i_z: usize| -> [f64; NPRIM] {
        let b = b_of(i_r, i_z);
        // [rho, ur, ut, uz, p, z, b, e, g1]
        [rho0, 0.0, 0.0, 0.0, P0, Z0, b, e0, 0.0]
    };
    let cpu = comb.xcheck_source_dense(&g, &ids, prim_of).expect("cpu combustion source");

    // Dense prim arrays for the GPU.
    let ncell = n_r * n_z;
    let rho = vec![rho0; ncell];
    let p = vec![P0; ncell];
    let z = vec![Z0; ncell];
    let e = vec![e0; ncell];
    let gas = vec![1.0f64; ncell];
    let mut b = vec![0.0f64; ncell];
    for i_r in 0..n_r {
        for i_z in 0..n_z {
            b[i_r * n_z + i_z] = b_of(i_r, i_z);
        }
    }

    let ax = |c: &ColumnMarshal, d: usize| (c.axis_points[d].clone(), c.axis_is_log[d] as i32);
    let (up, ulp) = ax(&m.unburnt_temp, 0);
    let (uh, ulh) = ax(&m.unburnt_temp, 1);
    let (uz, ulz) = ax(&m.unburnt_temp, 2);
    let ustr: Vec<i32> = m.unburnt_temp.strides.iter().map(|&s| s as i32).collect();
    let (ip, ilp) = ax(&m.flame, 0);
    let (itax, ilh) = ax(&m.flame, 1);
    let (iz, ilz) = ax(&m.flame, 2);
    let istr: Vec<i32> = m.flame.strides.iter().map(|&s| s as i32).collect();

    let run = |rate: &mut [f64]| unsafe {
        gpu_combustion_source(
            rho.as_ptr(), p.as_ptr(), z.as_ptr(), b.as_ptr(), e.as_ptr(), gas.as_ptr(),
            n_r as i32, n_z as i32, r_min, dr, dz,
            up.as_ptr(), up.len() as i32, ulp, uh.as_ptr(), uh.len() as i32, ulh,
            uz.as_ptr(), uz.len() as i32, ulz, ustr.as_ptr(),
            m.unburnt_temp.data.as_ptr(), m.unburnt_temp.value_is_log as i32,
            m.unburnt_rho.data.as_ptr(), m.unburnt_rho.value_is_log as i32,
            ip.as_ptr(), ip.len() as i32, ilp, itax.as_ptr(), itax.len() as i32, ilh,
            iz.as_ptr(), iz.len() as i32, ilz, istr.as_ptr(),
            m.flame.data.as_ptr(), m.flame.value_is_log as i32,
            m.hu_floor, m.hu_ceil, m.p_floor, m.tu_floor, m.wrinkling, m.theta,
            rate.as_mut_ptr(),
        )
    };
    let mut grate = vec![0.0f64; ncell * NCOMP];
    let mut grate2 = vec![0.0f64; ncell * NCOMP];
    run(&mut grate);
    run(&mut grate2);
    let bit_identical = grate == grate2;

    // Compare the I_RB source over interior cells (the compact diffusion stencil
    // is all-interior there — domain-edge cells' zero-flux BC is the CPU's too,
    // but stay clean of it). Also confirm the source is actually non-trivial.
    let mut worst_abs = 0.0f64;
    let mut worst_rel = 0.0f64;
    let mut max_src = 0.0f64;
    let mut n_active = 0usize;
    for i_r in 1..n_r - 1 {
        for i_z in 1..n_z - 1 {
            let c = i_r * n_z + i_z;
            let a = cpu[c];
            let bb = grate[c * NCOMP + I_RB];
            worst_abs = worst_abs.max((a - bb).abs());
            let s = a.abs().max(bb.abs());
            if s > 1e-6 {
                worst_rel = worst_rel.max((a - bb).abs() / s);
            }
            max_src = max_src.max(a.abs());
            if a.abs() > 1e-3 * max_src.max(1e-30) {
                n_active += 1;
            }
        }
    }
    println!("S13c — COMBUSTION SOURCE (Nagumo + front diffusion) on-device");
    println!("  grid {n_r}×{n_z} (N_θ=1 box), design flame state P0={P0:.0e} H0={H0:.3e} Z0={Z0}");
    println!("  max |source| = {max_src:.3e} kg/m³/s, active interior cells = {n_active}");
    println!("  worst abs diff = {worst_abs:.3e}");
    println!("  worst rel diff = {worst_rel:.3e}   (declared ECT {ECT:.0e})");
    println!(
        "  GPU same-build rerun: {}",
        if bit_identical { "BIT-IDENTICAL" } else { "*** DIFFERS ***" }
    );
    assert!(max_src > 0.0, "combustion source identically zero — fixture not exercising the reaction");
    assert!(worst_rel < ECT, "CPU↔GPU combustion source diverged beyond ECT: {worst_rel:.3e}");
    assert!(bit_identical, "GPU combustion source rerun not deterministic");
    println!("ALL PASS");
}
