//! S13c residency cross-check: the BLEND EOS (p,h,Z) mid-b projection ported to
//! device (`cuda/residency_blend_eos.cu`: `k_blend_project`) scored against the
//! bit-exact CPU reference `BurnBlendEos::prim_checked` (the two-branch
//! shifting-equilibrium projection) on the PRODUCTION LOX/LH₂ surfaces. This is
//! the class-R auto-ignition prerequisite: the per-refreeze state query is this
//! projection. The mid-b path (both branches live) is the genuinely new blend
//! physics; the pure limits delegate to the resident TableEos projection.
//!
//! Two-root resolves by CONTINUITY (first scan crossing) — Ben ruling
//! 2026-08-31; the fixture stays unique-root so it matches the CPU
//! bit-for-formula. Machine-precision ECT expected (interp reduction order
//! identical CPU↔GPU).
use crucible_solvers::euler::{BurnBlendEos, Cons, EosLaw, TableEos};
use crucible_tables::{ColumnMarshal, Pin, Table};

const ECT: f64 = 1e-9;

unsafe extern "C" {
    #[allow(clippy::too_many_arguments)]
    fn gpu_blend_project(
        rho: *const f64,
        e: *const f64,
        z: *const f64,
        b: *const f64,
        n: i32,
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
        pu_lo: f64,
        pu_hi: f64,
        huf: f64,
        huc: f64,
        pb_lo: f64,
        pb_hi: f64,
        hbf: f64,
        hbc: f64,
        h_off: f64,
        p_out: *mut f64,
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

// Host mirror of BoundColumn::interpolate (for building on-surface fixtures).
fn interp(col: &ColumnMarshal, qp: f64, qh: f64, qz: f64) -> f64 {
    let find = |pts: &[f64], q: f64| {
        pts.iter()
            .filter(|&&p| p <= q)
            .count()
            .min(pts.len() - 1)
            .max(1)
            - 1
    };
    let frac = |pts: &[f64], i: usize, lg: bool, q: f64| {
        if lg {
            (q.ln() - pts[i].ln()) / (pts[i + 1].ln() - pts[i].ln())
        } else {
            (q - pts[i]) / (pts[i + 1] - pts[i])
        }
    };
    let (pp, hp, zp) = (
        &col.axis_points[0],
        &col.axis_points[1],
        &col.axis_points[2],
    );
    let (ip, ih, iz) = (find(pp, qp), find(hp, qh), find(zp, qz));
    let cell = [ip, ih, iz];
    let tt = [
        frac(pp, ip, col.axis_is_log[0], qp),
        frac(hp, ih, col.axis_is_log[1], qh),
        frac(zp, iz, col.axis_is_log[2], qz),
    ];
    let mut acc = 0.0;
    for corner in 0..8usize {
        let mut w = 1.0;
        let mut idx = 0usize;
        for d in 0..3 {
            let up = (corner >> d) & 1 == 1;
            w *= if up { tt[d] } else { 1.0 - tt[d] };
            idx += (cell[d] + usize::from(up)) * col.strides[d];
        }
        let v = if col.value_is_log {
            col.data[idx].ln()
        } else {
            col.data[idx]
        };
        acc += w * v;
    }
    if col.value_is_log { acc.exp() } else { acc }
}

const EPS_B_PURE_UNBURNT: f64 = 1.0e-9;
const EPS_B_PURE_BURNT: f64 = 2.0e-3;

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
    let blend = BurnBlendEos::new(
        TableEos::bind(&ut).expect("unburnt binds"),
        TableEos::bind(&bt).expect("burnt binds"),
    );
    let (urho, brho, uenv, benv, h_off) = blend.xcheck_blend_marshal();
    let (pu, hu, zu) = (uenv[0], uenv[1], uenv[2]);
    let (pb, hb, zb) = (benv[0], benv[1], benv[2]);

    // Host partition_h → h_u (blend_eos.rs).
    let partition_h = |h: f64, b: f64| -> (f64, f64) {
        let (huf, huc) = (hu.0, hu.1);
        if (h >= huf && h <= huc) || b <= EPS_B_PURE_UNBURNT {
            return (h, h);
        }
        let hu_pin = h.clamp(huf, huc);
        let (hbl, hbh) = (hb.0 - h_off, hb.1 - h_off);
        let balance = (h - (1.0 - b) * hu_pin) / b;
        (hu_pin, balance.clamp(hbl, hbh))
    };
    let blend_inv_rho = |p: f64, h: f64, z: f64, b: f64| -> f64 {
        let (h_u, h_b) = partition_h(h, b);
        let mut v = 0.0;
        if b < 1.0 - EPS_B_PURE_BURNT {
            v += (1.0 - b) / interp(&urho, p, h_u, z);
        }
        if b > EPS_B_PURE_UNBURNT {
            v += b / interp(&brho, p, h_b + h_off, z);
        }
        v
    };

    // Mid-b ON-SURFACE states: p,z interior to both envelopes; h ABOVE the
    // unburnt ceiling (so the partition genuinely splits — the two-sub-state
    // physics) but inside the products window; b mid. rho = 1/blend_inv_rho.
    let p_lo = pu.0.max(pb.0);
    let p_hi = pu.1.min(pb.1);
    let z_lo = zu.0.max(zb.0);
    let z_hi = zu.1.min(zb.1);
    let inset = |a: f64, b: f64, f: f64| a + (b - a) * (0.2 + 0.6 * f);
    let mut rho = vec![];
    let mut e = vec![];
    let mut z = vec![];
    let mut bb = vec![];
    for ip in 0..4 {
        for ib in 0..4 {
            for ih in 0..4 {
                let p = inset(p_lo, p_hi, ip as f64 / 3.0);
                let zz = inset(z_lo, z_hi, 0.5);
                let b = 0.2 + 0.6 * ib as f64 / 3.0; // 0.2..0.8
                // h between just above the unburnt ceiling and mid products window
                let h_top = (hb.1 - h_off).min(hu.1 * 4.0);
                let h = inset(hu.1 * 1.02, h_top, ih as f64 / 3.0);
                let inv = blend_inv_rho(p, h, zz, b);
                if !(inv.is_finite() && inv > 0.0) {
                    continue;
                }
                let r = 1.0 / inv;
                let ee = h - p / r;
                if !(r.is_finite() && ee.is_finite()) {
                    continue;
                }
                rho.push(r);
                e.push(ee);
                z.push(zz);
                bb.push(b);
                let _ = ih;
            }
        }
    }
    let n = rho.len();
    assert!(n >= 8, "too few valid mid-b fixture states: {n}");

    // CPU oracle via the real blend projection (prim_checked).
    let mut cp = vec![0.0f64; n];
    let mut ok = vec![true; n];
    for i in 0..n {
        let u: Cons = [
            rho[i],
            0.0,
            0.0,
            0.0,
            rho[i] * e[i],
            rho[i] * z[i],
            rho[i] * bb[i],
        ];
        match blend.prim_checked(&u) {
            Ok(w) => cp[i] = w[4],
            Err(_) => ok[i] = false,
        }
    }

    // GPU.
    let ax = |c: &ColumnMarshal, d: usize| (c.axis_points[d].clone(), c.axis_is_log[d] as i32);
    let (up, ulp) = ax(&urho, 0);
    let (uh, ulh) = ax(&urho, 1);
    let (uz, ulz) = ax(&urho, 2);
    let ustr: Vec<i32> = urho.strides.iter().map(|&s| s as i32).collect();
    let (bp, blp) = ax(&brho, 0);
    let (bh, blh) = ax(&brho, 1);
    let (bz, blz) = ax(&brho, 2);
    let bstr: Vec<i32> = brho.strides.iter().map(|&s| s as i32).collect();
    let mut gp = vec![0.0f64; n];
    let mut gp2 = vec![0.0f64; n];
    let run = |out: &mut [f64]| unsafe {
        gpu_blend_project(
            rho.as_ptr(),
            e.as_ptr(),
            z.as_ptr(),
            bb.as_ptr(),
            n as i32,
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
            pu.0,
            pu.1,
            hu.0,
            hu.1,
            pb.0,
            pb.1,
            hb.0,
            hb.1,
            h_off,
            out.as_mut_ptr(),
        )
    };
    run(&mut gp);
    run(&mut gp2);
    let bit_identical = gp == gp2;

    let mut worst_rel = 0.0f64;
    let mut worst_abs = 0.0f64;
    let mut n_cmp = 0usize;
    for i in 0..n {
        if !ok[i] || !gp[i].is_finite() {
            continue;
        }
        worst_abs = worst_abs.max((cp[i] - gp[i]).abs());
        let s = cp[i].abs().max(gp[i].abs());
        if s > 1e-6 {
            worst_rel = worst_rel.max((cp[i] - gp[i]).abs() / s);
        }
        n_cmp += 1;
    }
    println!("S13c — BLEND EOS mid-b (p,h,Z) projection on-device");
    println!(
        "  surfaces: unburnt v0.3.0 ⊕ burnt v0.4.0; {n} states, {n_cmp} compared (mid-b, both branches)"
    );
    println!("  worst abs diff = {worst_abs:.3e}");
    println!("  worst rel diff = {worst_rel:.3e}   (declared ECT {ECT:.0e})");
    println!(
        "  GPU same-build rerun: {}",
        if bit_identical {
            "BIT-IDENTICAL"
        } else {
            "*** DIFFERS ***"
        }
    );
    assert!(n_cmp >= 8, "too few valid comparisons: {n_cmp}");
    assert!(
        worst_rel < ECT,
        "CPU↔GPU blend projection diverged beyond ECT: {worst_rel:.3e}"
    );
    assert!(
        bit_identical,
        "GPU blend projection rerun not deterministic"
    );
    println!("ALL PASS");
}
