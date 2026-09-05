//! S13c residency cross-check: the CLASS-R implicit AUTO-IGNITION node solve
//! ported to device (`cuda/residency_class_r.cu`: `k_class_r`) scored against
//! the bit-exact CPU reference `Combustion::implicit_auto_update` (pub — called
//! directly) on the PRODUCTION LOX/LH₂ surfaces. The BE-at-frozen-τ loop with
//! `N_TAU_REFREEZE = 2`, the exact cap parking, the symmetric base projection,
//! the cold floors — and THE new device surface: the blend's full `prim_checked`
//! three-way branch select (pure-unburnt / mid-b / pure-burnt) that a stiff node
//! walks within one solve as x: base → mid → cap.
//!
//! Fixture: the hot unburnt design state (T_u ≈ 1050 K, real finite τ_ign — the
//! `implicit_node_solve_parks_exactly_at_any_stiffness` state) with the initial
//! burn fraction b₀ ∈ {0, 0.3, 0.6}, p and h spread, and the node weight swept
//! over w/τ ∈ {1e-3 (mild: explicit limit), 0.3, 1, 3 (the BE knee band), 30,
//! 1e6 (parks exactly)}; plus the three guard corners (advected ρb past the cap,
//! a negative quadrature base, base == cap). Every regime must be represented
//! in the compared set or the harness refuses.
use crucible_solvers::euler::{
    BURN_COMPLETE, BurnBlendEos, Combustion, Cons, EosLaw, I_RB, IgnitionColumns, THETA_CELLS,
    TableEos,
};
use crucible_tables::{ColumnMarshal, Pin, Table};

const NC: usize = 7;
const ECT: f64 = 1e-9;

// Design hot unburnt state (solv4_combustion.rs): P0, Z0 + H_HOT (T_u ≈ 1050 K).
const P0: f64 = 1.0e5;
const Z0: f64 = 0.167;
const H_HOT: f64 = 2.5e6;

unsafe extern "C" {
    #[allow(clippy::too_many_arguments)]
    fn gpu_class_r_update(
        u: *const f64,
        base: *const f64,
        n: i32,
        w_new: f64,
        up: *const f64,
        unp: i32,
        ulp: i32,
        uh: *const f64,
        unh: i32,
        ulh: i32,
        uz: *const f64,
        unz: i32,
        ulz: i32,
        ustr: *const i32,
        urho_data: *const f64,
        urho_vlog: i32,
        utemp_data: *const f64,
        utemp_vlog: i32,
        bp: *const f64,
        bnp: i32,
        blp: i32,
        bh: *const f64,
        bnh: i32,
        blh: i32,
        bz: *const f64,
        bnz: i32,
        blz: i32,
        bstr: *const i32,
        brho_data: *const f64,
        brho_vlog: i32,
        ip: *const f64,
        inp: i32,
        ilp: i32,
        it: *const f64,
        inh: i32,
        ilh: i32,
        iz: *const f64,
        inz: i32,
        ilz: i32,
        istr: *const i32,
        dly_data: *const f64,
        dly_vlog: i32,
        pu_lo: f64,
        pu_hi: f64,
        hu_floor: f64,
        hu_ceil: f64,
        pb_lo: f64,
        pb_hi: f64,
        hb_lo: f64,
        hb_hi: f64,
        h_off: f64,
        p_floor: f64,
        tu_floor: f64,
        x_out: *mut f64,
        rate_out: *mut f64,
    );
}

fn open(file: &str, group: &str, pins_toml: &str) -> Table {
    let path = format!("{}/../../tables/chem/{file}", env!("CARGO_MANIFEST_DIR"));
    let doc: toml::Table = pins_toml.parse().expect("pins parse");
    let entry = doc[group].as_table().expect("group entry");
    let pin = Pin {
        data_version: entry["data_version"].as_str().expect("ver").to_string(),
        content_digest: entry
            .get("content_digest")
            .and_then(|v| v.as_str())
            .map(str::to_string),
    };
    Table::open(&path, group, &pin).expect("table loads under its pin")
}

/// Which `prim_checked` branch a burn fraction selects (the CPU's need_u/need_b).
fn branch_of(b: f64) -> &'static str {
    if b <= 1.0e-9 {
        "pure-unburnt"
    } else if b >= 1.0 - 2.0 * BURN_COMPLETE {
        "pure-burnt"
    } else {
        "mid-b"
    }
}

fn main() {
    crucible_gpu::keep_kernels(); // keeps this bin's CUDA entry points in the link (lib.rs)
    let ut = open(
        "lox_lh2_unburnt_v0.3.0.h5",
        "/chem/lox_lh2/unburnt",
        include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../tables/chem/lox_lh2_unburnt_v0.3.0.pins.toml"
        )),
    );
    let bt = open(
        "lox_lh2_v0.4.0.h5",
        "/chem/lox_lh2/equilibrium",
        include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../tables/chem/lox_lh2_v0.4.0.pins.toml"
        )),
    );
    let it = open(
        "lox_lh2_ignition_v0.3.0.h5",
        "/chem/lox_lh2/ignition",
        include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../tables/chem/lox_lh2_ignition_v0.3.0.pins.toml"
        )),
    );
    let blend = BurnBlendEos::new(
        TableEos::bind(&ut).expect("unburnt binds"),
        TableEos::bind(&bt).expect("burnt binds"),
    );
    let ign = IgnitionColumns::bind(&it).expect("ignition binds");
    let comb = Combustion {
        blend: &blend,
        ignition: ign,
        wrinkling: 1.0,
        theta: THETA_CELLS,
    };

    // Marshal: both ρ columns + envelopes + h_off; the unburnt T column; the τ_ign column.
    let (urho, brho, uenv, benv, h_off) = blend.xcheck_blend_marshal();
    let (utemp, _, _) = blend.xcheck_unburnt_marshal();
    let (_, p_floor, tu_floor) = comb.ignition.xcheck_marshal_flame();
    let dly = comb.ignition.xcheck_marshal_delay();

    // τ at the design state (the same query path the solve uses).
    let u0 = blend
        .cons_from_phzb(P0, H_HOT, Z0, 0.0, [0.0, 0.0, 0.0])
        .expect("design state");
    let w0 = blend.prim_checked(&u0).expect("design prim");
    let t_u0 = blend.unburnt_temperature(&w0).expect("T_u");
    let tau0 = comb.ignition.induction_time(P0, t_u0, Z0).expect("τ_ign");

    // Cells: (b₀, p, h) spread + the guard corners. The conserved state comes
    // from the blend's own constructor (on-surface by construction).
    let mut cells: Vec<Cons> = vec![];
    let mut bases: Vec<f64> = vec![];
    let mut tags: Vec<&'static str> = vec![];
    for &b0 in &[0.0f64, 0.3, 0.6] {
        for ip in 0..3 {
            for ih in 0..3 {
                let p = P0 * (1.0 + 0.5 * ip as f64);
                let h = H_HOT * (0.94 + 0.06 * ih as f64);
                let Ok(u) = blend.cons_from_phzb(p, h, Z0, b0, [30.0, 0.0, -12.0]) else {
                    continue;
                };
                cells.push(u);
                bases.push(u[I_RB]); // base = the advected ρb (no quadrature terms)
                tags.push("regular");
            }
        }
    }
    // Guard corners on the design cell.
    let rho0 = u0[0];
    let cap0 = rho0 * (1.0 - BURN_COMPLETE);
    {
        let mut u = u0;
        u[I_RB] = cap0 * 1.0000001; // advected ρb past the cap: state unchanged
        cells.push(u);
        bases.push(u[I_RB]);
        tags.push("adv-past-cap");
    }
    {
        cells.push(u0);
        bases.push(-0.08 * rho0); // the measured ◆C2 negative quadrature base
        tags.push("negative-base");
    }
    {
        cells.push(u0);
        bases.push(cap0); // base exactly at the cap: parked
        tags.push("base-at-cap");
    }
    {
        let u = blend
            .cons_from_phzb(P0, H_HOT, Z0, 0.3, [0.0, 0.0, 0.0])
            .expect("mid cell");
        cells.push(u);
        bases.push(u[I_RB] + 0.4 * rho0); // quadrature-lifted base (mid-b start)
        tags.push("lifted-base");
    }
    let n = cells.len();
    assert!(n >= 16, "too few fixture cells: {n}");
    let mut u_flat = vec![0.0f64; n * NC];
    for (i, u) in cells.iter().enumerate() {
        u_flat[i * NC..(i + 1) * NC].copy_from_slice(u);
    }

    let ax = |c: &ColumnMarshal, d: usize| (c.axis_points[d].clone(), c.axis_is_log[d] as i32);
    let (up, ulp) = ax(&urho, 0);
    let (uh, ulh) = ax(&urho, 1);
    let (uz, ulz) = ax(&urho, 2);
    let ustr: Vec<i32> = urho.strides.iter().map(|&s| s as i32).collect();
    let (bp, blp) = ax(&brho, 0);
    let (bh, blh) = ax(&brho, 1);
    let (bz, blz) = ax(&brho, 2);
    let bstr: Vec<i32> = brho.strides.iter().map(|&s| s as i32).collect();
    let (ip, ilp) = ax(&dly, 0);
    let (itx, ilh) = ax(&dly, 1);
    let (iz, ilz) = ax(&dly, 2);
    let istr: Vec<i32> = dly.strides.iter().map(|&s| s as i32).collect();
    assert_eq!(
        utemp.strides, urho.strides,
        "unburnt T and ρ columns must share axes"
    );

    let run = |w_new: f64, x: &mut [f64], r: &mut [f64]| unsafe {
        gpu_class_r_update(
            u_flat.as_ptr(),
            bases.as_ptr(),
            n as i32,
            w_new,
            up.as_ptr(),
            up.len() as i32,
            ulp,
            uh.as_ptr(),
            uh.len() as i32,
            ulh,
            uz.as_ptr(),
            uz.len() as i32,
            ulz,
            ustr.as_ptr(),
            urho.data.as_ptr(),
            urho.value_is_log as i32,
            utemp.data.as_ptr(),
            utemp.value_is_log as i32,
            bp.as_ptr(),
            bp.len() as i32,
            blp,
            bh.as_ptr(),
            bh.len() as i32,
            blh,
            bz.as_ptr(),
            bz.len() as i32,
            blz,
            bstr.as_ptr(),
            brho.data.as_ptr(),
            brho.value_is_log as i32,
            ip.as_ptr(),
            ip.len() as i32,
            ilp,
            itx.as_ptr(),
            itx.len() as i32,
            ilh,
            iz.as_ptr(),
            iz.len() as i32,
            ilz,
            istr.as_ptr(),
            dly.data.as_ptr(),
            dly.value_is_log as i32,
            uenv[0].0,
            uenv[0].1,
            uenv[1].0,
            uenv[1].1,
            benv[0].0,
            benv[0].1,
            benv[1].0,
            benv[1].1,
            h_off,
            p_floor,
            tu_floor,
            x.as_mut_ptr(),
            r.as_mut_ptr(),
        )
    };

    const W_OVER_TAU: [f64; 6] = [1.0e-3, 0.3, 1.0, 3.0, 30.0, 1.0e6];
    let mut worst_rel_x = 0.0f64;
    let mut worst_rel_r = 0.0f64;
    let mut n_cmp = 0usize;
    let mut n_skip = 0usize;
    let mut bit_identical = true;
    let mut n_parked = 0usize;
    let mut seen_branch = [0usize; 3]; // final-state branch: unburnt / mid / burnt
    println!("S13c — CLASS-R implicit auto-ignition node solve on-device");
    println!(
        "  surfaces: unburnt v0.3.0 ⊕ burnt v0.4.0 + ignition v0.3.0; {n} cells; τ(design) = {tau0:.3e} s"
    );
    for &ratio in &W_OVER_TAU {
        let w_new = ratio * tau0;
        let mut gx = vec![0.0f64; n];
        let mut gr = vec![0.0f64; n];
        let mut gx2 = vec![0.0f64; n];
        let mut gr2 = vec![0.0f64; n];
        run(w_new, &mut gx, &mut gr);
        run(w_new, &mut gx2, &mut gr2);
        bit_identical &= gx == gx2 && gr == gr2;
        let mut wr_x = 0.0f64;
        let mut wr_r = 0.0f64;
        let mut parked = 0usize;
        for i in 0..n {
            let Ok((cx, cr)) = comb.implicit_auto_update(&cells[i], bases[i], w_new) else {
                n_skip += 1;
                continue;
            };
            if !gx[i].is_finite() {
                panic!(
                    "device returned NaN at cell {i} ({}) w/τ={ratio:.1e}: tangency-corner sentinel \
                     reached inside the fixture — the fixture must stay interior",
                    tags[i]
                );
            }
            let rho = cells[i][0];
            let cap = rho * (1.0 - BURN_COMPLETE);
            if cx == cap {
                parked += 1;
                // The exact parking must be reproduced EXACTLY (the fixed point
                // is the declared law's integral, not a tolerance).
                assert!(
                    gx[i] == cap,
                    "cell {i} ({}) w/τ={ratio:.1e}: CPU parked at cap, GPU x = {:.17e} vs cap {:.17e}",
                    tags[i],
                    gx[i],
                    cap
                );
            }
            let sx = cx.abs().max(gx[i].abs());
            if sx > 1e-300 {
                wr_x = wr_x.max((cx - gx[i]).abs() / sx);
            }
            let sr = cr.abs().max(gr[i].abs());
            if sr > 1e-300 {
                wr_r = wr_r.max((cr - gr[i]).abs() / sr);
            }
            let b_final = cx / rho;
            let k = match branch_of(b_final) {
                "pure-unburnt" => 0,
                "mid-b" => 1,
                _ => 2,
            };
            seen_branch[k] += 1;
            n_cmp += 1;
        }
        n_parked += parked;
        worst_rel_x = worst_rel_x.max(wr_x);
        worst_rel_r = worst_rel_r.max(wr_r);
        println!(
            "  w/τ = {ratio:>7.1e}: worst rel Δρb = {wr_x:.3e}, worst rel Δrate = {wr_r:.3e}, parked {parked}/{n}"
        );
    }
    println!("  compared {n_cmp} (skipped {n_skip} CPU-refused), parked total {n_parked}");
    println!(
        "  final-state branch coverage: pure-unburnt {}, mid-b {}, pure-burnt {}",
        seen_branch[0], seen_branch[1], seen_branch[2]
    );
    println!(
        "  worst rel diff ρb = {worst_rel_x:.3e}, rate = {worst_rel_r:.3e}   (declared ECT {ECT:.0e})"
    );
    println!(
        "  GPU same-build rerun: {}",
        if bit_identical {
            "BIT-IDENTICAL"
        } else {
            "*** DIFFERS ***"
        }
    );
    assert!(n_cmp >= 6 * 16, "too few valid comparisons: {n_cmp}");
    assert!(
        seen_branch[1] > 0 && seen_branch[2] > 0,
        "fixture must exercise BOTH the mid-b and the parked pure-burnt regimes"
    );
    assert!(
        n_parked > 0,
        "fixture never parks — the stiff regime is not exercised"
    );
    assert!(
        worst_rel_x < ECT,
        "CPU↔GPU class-R ρb diverged beyond ECT: {worst_rel_x:.3e}"
    );
    assert!(
        worst_rel_r < ECT,
        "CPU↔GPU class-R rate diverged beyond ECT: {worst_rel_r:.3e}"
    );
    assert!(bit_identical, "GPU class-R rerun not deterministic");
    println!("ALL PASS");
}
