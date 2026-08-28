//! S12 GPU-spike harness: the CPU↔GPU tolerance cross-check + same-build
//! rerun bit-identity for the HLLC-Batten flux (hot kernel a's Riemann core).
//! The CPU side is the bit-exact reference `crucible_solvers::euler::hllc_flux`
//! (GammaLaw); the GPU side is the CUDA port in `cuda/hllc_kernel.cu`. We do NOT
//! promise cross-device bit-identity (FMA contraction differs) — this measures
//! the tolerance to declare (META-1 §2.5 ECT), and proves the GPU path is
//! deterministic run-to-run (same build, same device).
use crucible_solvers::euler::{GammaLaw, NCOMP, NPRIM, Prim, hllc_flux};

unsafe extern "C" {
    fn gpu_hllc(wl: *const f64, wr: *const f64, n: i32, dir: i32, out: *mut f64);
}

fn main() {
    let n: usize = 1 << 20; // 1,048,576 face pairs
    let eos = GammaLaw { gamma: 1.4 };

    // Deterministic LCG fixture of admissible (rho>0, p>0) primitive pairs.
    let mut s: u64 = 0x0123_4567_89ab_cdef;
    let mut u01 = || {
        s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((s >> 33) as f64) / ((1u64 << 31) as f64)
    };
    let mut mk = || -> Prim {
        let rho = 0.1 + 3.0 * u01();
        let ur = -2.0 + 4.0 * u01();
        let ut = -2.0 + 4.0 * u01();
        let uz = -2.0 + 4.0 * u01();
        let p = 0.1 + 5.0 * u01();
        let c = u01();
        let b = u01();
        [rho, ur, ut, uz, p, c, b, 0.0, 0.0]
    };
    let wl: Vec<Prim> = (0..n).map(|_| mk()).collect();
    let wr: Vec<Prim> = (0..n).map(|_| mk()).collect();

    // flatten for the FFI (stride NPRIM)
    let mut wl_flat = vec![0.0f64; n * NPRIM];
    let mut wr_flat = vec![0.0f64; n * NPRIM];
    for i in 0..n {
        wl_flat[i * NPRIM..i * NPRIM + NPRIM].copy_from_slice(&wl[i]);
        wr_flat[i * NPRIM..i * NPRIM + NPRIM].copy_from_slice(&wr[i]);
    }

    let mut worst_abs = 0.0f64;
    let mut worst_rel = 0.0f64;
    let mut gpu_first: Vec<f64> = Vec::new();
    for dir in 1..=3usize {
        // CPU reference
        let mut cpu = vec![0.0f64; n * NCOMP];
        for i in 0..n {
            let f = hllc_flux(&wl[i], &wr[i], dir, &eos);
            cpu[i * NCOMP..i * NCOMP + NCOMP].copy_from_slice(&f);
        }
        // GPU
        let mut gpu = vec![0.0f64; n * NCOMP];
        unsafe {
            gpu_hllc(
                wl_flat.as_ptr(),
                wr_flat.as_ptr(),
                n as i32,
                dir as i32,
                gpu.as_mut_ptr(),
            );
        }
        for k in 0..n * NCOMP {
            let d = (cpu[k] - gpu[k]).abs();
            if d > worst_abs {
                worst_abs = d;
            }
            let denom = cpu[k].abs().max(1e-12);
            let rel = d / denom;
            if rel > worst_rel {
                worst_rel = rel;
            }
        }
        if dir == 3 {
            gpu_first = gpu;
        }
    }

    // same-build rerun bit-identity (GPU, dir=3)
    let mut gpu_second = vec![0.0f64; n * NCOMP];
    unsafe {
        gpu_hllc(
            wl_flat.as_ptr(),
            wr_flat.as_ptr(),
            n as i32,
            3,
            gpu_second.as_mut_ptr(),
        );
    }
    let bit_identical = gpu_first == gpu_second;

    println!("HLLC CPU↔GPU cross-check over {n} face pairs × 3 directions:");
    println!("  worst abs diff = {worst_abs:.3e}");
    println!("  worst rel diff = {worst_rel:.3e}   (declared cross-device ECT tolerance)");
    println!(
        "  GPU same-build rerun: {}",
        if bit_identical {
            "BIT-IDENTICAL"
        } else {
            "*** DIFFERS ***"
        }
    );
    assert!(
        worst_rel < 1e-9,
        "CPU↔GPU HLLC diverged beyond 1e-9 — not an FMA-order effect"
    );
    assert!(bit_identical, "GPU rerun not deterministic");
    println!("PASS");
}
