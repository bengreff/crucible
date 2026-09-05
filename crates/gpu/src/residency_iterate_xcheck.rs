//! S13c residency cross-check: the FULL RESIDENT class-D gas Picard iterate —
//! the composition the class-D residency exists for. It marches one SDC-inner
//! implicit-diffusion sweep entirely on-device (fields resident across the
//! whole iterate): fill_lag_gradients → assemble(sol,lag) → {Ur,Uz,Om} CG
//! solves → re-assemble → {T,C} CG solves — the FORCING (S13c assemble) + the
//! SOLVER (S13b CG) + the RHS builder (fill_gas_rhs) composed. Scored against
//! the bit-exact CPU reference `crucible_solvers::sdc::xcheck_class_d_iterate_dense`
//! (which runs the real production `assemble_rates`/`fill_gas_rhs`/`cg_solve`
//! in the identical order).
//!
//! Determinism (META-1 §2.5): the iterate chains gather kernels + the
//! fixed-topology CG reduction; same-build reruns bit-identical. CPU↔GPU is the
//! accumulated FMA-order/reduction-shape tolerance over 5 converged CG solves
//! (declared looser than a single op — the same ECT class as the S13b CG).
//!
//! SCOPE: N_θ=1 box, free BCs, GammaLaw-class operands seeded directly (the
//! EOS-seam-agnostic resident step — the real TableEos just swaps the operand
//! decode). Real walls / θ-tensor / combustion / cut geometry are S13c's
//! remaining legs.
use crucible_grid::{Grid, GridSpec};
use crucible_solvers::gas_diffusion::{FaceGasBc, GasDiffBcs, GasDiffusion, GasTransportField};
use crucible_solvers::sdc::xcheck_class_d_iterate_dense;
use crucible_solvers::transport::TransportProps;

const NCOMP: usize = 7;
/// Converged-solve ECT (META-1 §2.5): 5 CG solves chained, so the accumulated
/// reduction-shape/FMA difference matches the S13b CG's 1e-8 band.
const ECT: f64 = 1e-8;

unsafe extern "C" {
    fn gpu_class_d_iterate(
        sol_rho: *mut f64,
        sol_ur: *mut f64,
        sol_om: *mut f64,
        sol_uz: *mut f64,
        sol_tt: *mut f64,
        sol_cc: *mut f64,
        lag_ur: *const f64,
        lag_uz: *const f64,
        mu: *const f64,
        kk: *const f64,
        rhod: *const f64,
        dhdz: *const f64,
        cv: *const f64,
        dlag: *const f64,
        gas: *const f64,
        n_r: i32,
        n_z: i32,
        r_min: f64,
        dr: f64,
        dz: f64,
        wqnew: f64,
    ) -> i32;
}

fn main() {
    crucible_gpu::keep_kernels(); // keeps this bin's CUDA entry points in the link (lib.rs)
    let (n_r, n_z) = (48usize, 96usize);
    let (r_min, dr, dz) = (0.5, 1.0 / n_r as f64, 1.0 / n_z as f64);
    let spec = GridSpec {
        r_min,
        dr,
        n_r,
        z_min: 0.0,
        dz,
        n_z,
        n_theta_max: 1,
        axisymmetry_assertion: true,
    };
    let g = Grid::build(spec, &["dummy"]).expect("valid spec");

    const MU: f64 = 0.01;
    const K: f64 = 0.02;
    const RHOD: f64 = 0.0125;
    const CV: f64 = 720.0;
    let mut tr = GasTransportField::alloc(&g);
    g.for_each_active_cell(|cell| {
        let f = 1.0 + 0.31 * (cell.i_r as f64) + 0.17 * (cell.i_z as f64);
        let props = TransportProps {
            mu: MU * f,
            k: K * f,
            cp: 1005.0,
            cp_film: 1005.0,
            cv: CV,
            rho_d: RHOD * f,
            dh_dz: 3.0e5 * (1.0 + 0.05 * (cell.i_r as f64 - cell.i_z as f64)),
            pr: 0.71,
        };
        tr.set(cell.bi, cell.idx, cell.i_r, cell.i_z, &props)
            .expect("valid transport");
    });

    let op = GasDiffusion::new(GasDiffBcs {
        r_inner: FaceGasBc::free(),
        r_outer: FaceGasBc::free(),
        z_lo: FaceGasBc::free(),
        z_hi: FaceGasBc::free(),
    });

    let sol_rho = |r: usize, z: usize| 1.0 + 0.2 * (1.7 * r as f64 + 0.9 * z as f64).sin();
    let sol_ur = |r: usize, z: usize| 0.30 * (1.1 * z as f64 - 0.4 * r as f64).sin();
    let sol_om = |r: usize, z: usize| 0.25 * (0.8 * r as f64 + 0.5 * z as f64).cos();
    let sol_uz = |r: usize, z: usize| 0.35 * (0.9 * r as f64 - 1.3 * z as f64).sin();
    let sol_tt = |r: usize, z: usize| 300.0 + 20.0 * (0.3 * r as f64 + 0.2 * z as f64).sin();
    let sol_cc = |r: usize, z: usize| 0.5 + 0.1 * (0.25 * r as f64).sin() * (0.2 * z as f64).cos();
    let lag_ur = |r: usize, z: usize| sol_ur(r, z) + 0.05 * (0.7 * r as f64 + 0.3 * z as f64).cos();
    let lag_uz = |r: usize, z: usize| sol_uz(r, z) - 0.04 * (0.5 * r as f64 - 0.6 * z as f64).sin();
    // A nonzero lagged rate so the (dstage − dlag) subtraction is exercised.
    let dlag_of = |r: usize, z: usize| {
        let mut a = [0.0f64; NCOMP];
        a[1] = 1e-3 * (0.4 * r as f64).sin();
        a[2] = 5e-4 * (0.3 * z as f64).cos();
        a[3] = 8e-4 * (0.2 * r as f64 - 0.1 * z as f64).sin();
        a[4] = 2e-2 * (0.15 * r as f64).cos();
        a[5] = 3e-4 * (0.25 * z as f64).sin();
        a
    };
    let wqnew = 5.0e-2;

    let d = xcheck_class_d_iterate_dense(
        &g, &op, &tr, wqnew, sol_rho, sol_ur, sol_om, sol_uz, sol_tt, sol_cc, lag_ur, lag_uz,
        dlag_of,
    );
    let ncell = n_r * n_z;

    // Flatten dlag for the FFI (stride NCOMP).
    let mut dlag_flat = vec![0.0f64; ncell * NCOMP];
    for c in 0..ncell {
        dlag_flat[c * NCOMP..c * NCOMP + NCOMP].copy_from_slice(&d.dlag[c]);
    }

    // GPU: seed the sol buffers with the INITIAL operands, run the iterate.
    let run = |rho: &Vec<f64>| -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
        let mut gr = rho.clone();
        let mut gur = d.init_ur.clone();
        let mut gom = d.init_om.clone();
        let mut guz = d.init_uz.clone();
        let mut gtt = d.init_tt.clone();
        let mut gcc = d.init_cc.clone();
        unsafe {
            gpu_class_d_iterate(
                gr.as_mut_ptr(),
                gur.as_mut_ptr(),
                gom.as_mut_ptr(),
                guz.as_mut_ptr(),
                gtt.as_mut_ptr(),
                gcc.as_mut_ptr(),
                d.lag_ur.as_ptr(),
                d.lag_uz.as_ptr(),
                d.mu.as_ptr(),
                d.k.as_ptr(),
                d.rhod.as_ptr(),
                d.dhdz.as_ptr(),
                d.cv.as_ptr(),
                dlag_flat.as_ptr(),
                d.gas.as_ptr(),
                n_r as i32,
                n_z as i32,
                r_min,
                dr,
                dz,
                wqnew,
            );
        }
        (gur, gom, guz, gtt, gcc)
    };
    let (gur, gom, guz, gtt, gcc) = run(&d.rho);
    let (gur2, gom2, guz2, gtt2, gcc2) = run(&d.rho);
    let bit_identical = gur == gur2 && gom == gom2 && guz == guz2 && gtt == gtt2 && gcc == gcc2;

    // Compare the updated sol over all gas cells, per component.
    let cpu = [&d.ur, &d.om, &d.uz, &d.tt, &d.cc];
    let gpu = [&gur, &gom, &guz, &gtt, &gcc];
    let names = ["u_r", "ω", "u_z", "T", "C"];
    println!("S13c — FULL RESIDENT class-D iterate (grid {n_r}×{n_z}, N_θ=1 box, wqnew={wqnew})");
    println!(
        "  {:>4}   {:>11}   {:>11}",
        "comp", "worst_abs", "worst_rel"
    );
    let mut overall = 0.0f64;
    for j in 0..5 {
        let (mut wa, mut wr) = (0.0f64, 0.0f64);
        for c in 0..ncell {
            if d.gas[c] == 0.0 {
                continue;
            }
            let (a, b) = (cpu[j][c], gpu[j][c]);
            wa = wa.max((a - b).abs());
            let s = a.abs().max(b.abs());
            if s > 1e-6 {
                wr = wr.max((a - b).abs() / s);
            }
        }
        overall = overall.max(wr);
        println!("  {:>4}   {wa:>11.3e}   {wr:>11.3e}", names[j]);
    }
    println!("  overall worst rel = {overall:.3e}   (declared ECT {ECT:.0e})");
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
        "CPU↔GPU resident class-D iterate diverged beyond ECT: {overall:.3e}"
    );
    assert!(bit_identical, "GPU iterate rerun not deterministic");
    println!("ALL PASS");
}
