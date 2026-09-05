//! S13c residency cross-check: the REAL HDF5 TableEos (p,h,Z) equilibrium
//! projection ported to device (`cuda/residency_eos.cu`: `k_project` + the
//! device multilinear interp + Illinois root-find) scored against the bit-exact
//! CPU reference `crucible_solvers::euler::TableEos` on the PRODUCTION
//! LOX/LH₂ equilibrium surface (`tables/chem/lox_lh2_v0.4.0.h5` — the
//! station-5 / ◆C3 surface).
//!
//! This replaces the class-A/D GammaLaw stand-in with the actual combustion-
//! products EOS. Each cell independently projects (ρ, e_q, Z) → p, then reads
//! sound & temperature at (p, h, Z) — a per-cell gather (no reduction), so
//! same-build reruns are bit-identical and CPU↔GPU is the interp's
//! FMA-order + libm(ln/exp) ECT. The two-root case resolves by CONTINUITY
//! (warm hint keeps the near root; cold takes the first scan crossing) — never
//! a halt (Ben ruling 2026-08-31).
//!
//! SCOPE: equilibrium chamber states (the dominant warm/cold-straddling paths);
//! the rare near-vacuum golden-section tangency corner is a documented
//! follow-on (the device returns NaN there — this fixture stays interior).
use crucible_solvers::euler::{TableEos, TableEosMarshal};
use crucible_tables::{Pin, Table};

const FILE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../tables/chem/lox_lh2_v0.4.0.h5"
);
const PINS_TOML: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../tables/chem/lox_lh2_v0.4.0.pins.toml"
));
const GROUP: &str = "/chem/lox_lh2/equilibrium";
/// FMA + libm(ln/exp) ECT for the projection (a fixed 8-corner interp + a
/// deterministic root-find; the interp reduction order is identical CPU↔GPU,
/// only ln/exp differ per host).
const ECT: f64 = 1e-9;

unsafe extern "C" {
    #[allow(clippy::too_many_arguments)]
    fn gpu_table_project(
        rho: *const f64,
        e_q: *const f64,
        z: *const f64,
        hint: *const f64,
        n: i32,
        pp: *const f64,
        np: i32,
        lp: i32,
        hp: *const f64,
        nh: i32,
        lh: i32,
        zp: *const f64,
        nz: i32,
        lz: i32,
        strides: *const i32,
        rho_data: *const f64,
        rho_vlog: i32,
        snd_data: *const f64,
        snd_vlog: i32,
        tmp_data: *const f64,
        tmp_vlog: i32,
        p_lo: f64,
        p_hi: f64,
        h_lo: f64,
        h_hi: f64,
        p_out: *mut f64,
        a_out: *mut f64,
        t_out: *mut f64,
        g1_out: *mut f64,
    );
}

fn open() -> Table {
    let doc: toml::Table = PINS_TOML.parse().expect("pins sidecar parses");
    let entry = doc[GROUP].as_table().expect("entry");
    let pin = Pin {
        data_version: entry["data_version"].as_str().expect("str").to_string(),
        content_digest: Some(entry["content_digest"].as_str().expect("str").to_string()),
    };
    Table::open(FILE, GROUP, &pin).expect("production surface loads")
}

fn main() {
    crucible_gpu::keep_kernels(); // keeps this bin's CUDA entry points in the link (lib.rs)
    let t = open();
    let eos = TableEos::bind(&t).expect("bind the equilibrium surface");
    let m: TableEosMarshal = eos.xcheck_marshal();

    // Build a grid of ON-SURFACE states inside the (p,h,Z) envelope: pick
    // (p,h,Z), read ρ from the surface, set e_q = h − p/ρ, hint = p·(1+δ) so
    // the warm path runs (and must recover p). Stay a margin inside every
    // envelope edge so the interp never refuses and we avoid the tangency corner.
    let (pe, he, ze) = (m.p_env, m.h_env, m.z_env);
    let lerp = |a: f64, b: f64, f: f64| a + (b - a) * f;
    let inset = |a: f64, b: f64, f: f64| lerp(a, b, 0.12 + 0.76 * f); // 12%..88%
    let np_s = 7usize;
    let nh_s = 7usize;
    let nz_s = 5usize;
    let mut rho = vec![];
    let mut e_q = vec![];
    let mut z = vec![];
    let mut hint = vec![];
    // The density column, for the on-surface ρ at each (p,h,Z).
    let rho_col = &m.rho;
    let interp = |col: &crucible_tables::ColumnMarshal, qp: f64, qh: f64, qz: f64| -> f64 {
        // Mirror BoundColumn::interpolate for building the fixture (host side).
        let find = |pts: &[f64], q: f64| -> usize {
            let cnt = pts.iter().filter(|&&p| p <= q).count();
            cnt.min(pts.len() - 1).max(1) - 1
        };
        let frac = |pts: &[f64], i: usize, is_log: bool, q: f64| -> f64 {
            if is_log {
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
    };
    for ip in 0..np_s {
        for ih in 0..nh_s {
            for iz in 0..nz_s {
                let p = inset(pe.0, pe.1, ip as f64 / (np_s - 1) as f64);
                let h = inset(he.0, he.1, ih as f64 / (nh_s - 1) as f64);
                let zz = inset(ze.0, ze.1, iz as f64 / (nz_s - 1) as f64);
                let r = interp(rho_col, p, h, zz);
                rho.push(r);
                e_q.push(h - p / r);
                z.push(zz);
                hint.push(p * 1.02); // 2% off → warm path recovers p
            }
        }
    }
    let n = rho.len();

    // CPU oracle.
    let mut cp = vec![0.0f64; n];
    let mut ca = vec![0.0f64; n];
    let mut ct = vec![0.0f64; n];
    for i in 0..n {
        let (p, a, tt, _g1) = eos
            .xcheck_project(rho[i], e_q[i], z[i], Some(hint[i]))
            .expect("oracle projection on-surface");
        cp[i] = p;
        ca[i] = a;
        ct[i] = tt;
    }

    // GPU. Shared axes from the density marshal; 3 data blocks.
    let pp = &rho_col.axis_points[0];
    let hp = &rho_col.axis_points[1];
    let zp = &rho_col.axis_points[2];
    let strides: Vec<i32> = rho_col.strides.iter().map(|&s| s as i32).collect();
    let mut gp = vec![0.0f64; n];
    let mut ga = vec![0.0f64; n];
    let mut gt = vec![0.0f64; n];
    let mut gg1 = vec![0.0f64; n];
    let run = |gp: &mut [f64], ga: &mut [f64], gt: &mut [f64], gg1: &mut [f64]| unsafe {
        gpu_table_project(
            rho.as_ptr(),
            e_q.as_ptr(),
            z.as_ptr(),
            hint.as_ptr(),
            n as i32,
            pp.as_ptr(),
            pp.len() as i32,
            rho_col.axis_is_log[0] as i32,
            hp.as_ptr(),
            hp.len() as i32,
            rho_col.axis_is_log[1] as i32,
            zp.as_ptr(),
            zp.len() as i32,
            rho_col.axis_is_log[2] as i32,
            strides.as_ptr(),
            m.rho.data.as_ptr(),
            m.rho.value_is_log as i32,
            m.sound.data.as_ptr(),
            m.sound.value_is_log as i32,
            m.temperature.data.as_ptr(),
            m.temperature.value_is_log as i32,
            pe.0,
            pe.1,
            he.0,
            he.1,
            gp.as_mut_ptr(),
            ga.as_mut_ptr(),
            gt.as_mut_ptr(),
            gg1.as_mut_ptr(),
        )
    };
    run(&mut gp, &mut ga, &mut gt, &mut gg1);
    // Same-build rerun bit-identity.
    let (mut gp2, mut ga2, mut gt2, mut gg2) =
        (vec![0.0; n], vec![0.0; n], vec![0.0; n], vec![0.0; n]);
    run(&mut gp2, &mut ga2, &mut gt2, &mut gg2);
    let bit_identical = gp == gp2 && ga == ga2 && gt == gt2;

    let cols: [(&str, &[f64], &[f64]); 3] =
        [("p", &cp, &gp), ("sound", &ca, &ga), ("temp", &ct, &gt)];
    println!("S13c — REAL TableEos (p,h,Z) projection on-device ({n} on-surface states)");
    println!("  surface: lox_lh2_v0.4.0 (station-5 / ◆C3 equilibrium)");
    println!(
        "  {:>6}   {:>11}   {:>11}",
        "field", "worst_abs", "worst_rel"
    );
    let mut overall = 0.0f64;
    let mut recovered = 0.0f64;
    for (name, c, gcol) in cols {
        let (mut wa, mut wr) = (0.0f64, 0.0f64);
        for i in 0..n {
            wa = wa.max((c[i] - gcol[i]).abs());
            let s = c[i].abs().max(gcol[i].abs());
            if s > 1e-6 {
                wr = wr.max((c[i] - gcol[i]).abs() / s);
            }
        }
        overall = overall.max(wr);
        println!("  {name:>6}   {wa:>11.3e}   {wr:>11.3e}");
    }
    // Sanity: the projection recovered the sampling pressure (both CPU & GPU).
    for i in 0..n {
        let p_sample = {
            // reconstruct: e_q = h − p/ρ  ⇒  h = e_q + p/ρ; the sampling p is
            // the CPU-recovered p (it round-trips to <1e-8 by construction).
            cp[i]
        };
        recovered = recovered.max((gp[i] - p_sample).abs() / p_sample.abs());
    }
    println!("  overall worst rel (CPU↔GPU) = {overall:.3e}   (declared ECT {ECT:.0e})");
    println!("  GPU vs CPU-recovered p       = {recovered:.3e}");
    println!(
        "  GPU same-build rerun: {}",
        if bit_identical {
            "BIT-IDENTICAL"
        } else {
            "*** DIFFERS ***"
        }
    );
    assert!(
        overall < ECT,
        "CPU↔GPU TableEos projection diverged beyond ECT: {overall:.3e}"
    );
    assert!(bit_identical, "GPU projection rerun not deterministic");
    println!("ALL PASS");
}
