//! S13b residency cross-check: the class-D IMPLICIT DIFFUSION per-component
//! symmetric CG ported to a device-resident kernel set —
//! `cuda/residency_diffusion.cu` — scored against the bit-exact CPU reference
//! `crucible_solvers::gas_diffusion::GasDiffusion::cg_solve` (via the additive
//! `xcheck_cg_dense` accessor) on a real (r,z) fixture.
//!
//! This is the "single biggest piece" of the S13b split (the brief): the
//! fixed-structure Jacobi-preconditioned CG that solves `(mass − wq·L)·δ = b`
//! for one gas component. The genuinely NEW on-device surface vs the S13
//! class-A residency is the REDUCTION — the CG dot products — realized as a
//! fixed-topology tree reduction (META-1 §2.5). The five gas components
//! (Ur,Uz,Om,T,C) all use the ONE component-generic solve; only `face_coef`
//! and `mass` differ, so validating all five validates the whole class-D CG.
//!
//! Determinism (META-1 §2.5): same-build GPU rerun is bit-identical; CPU↔GPU
//! agree to the declared ECT — but this is a CONVERGED SOLVE, not a single op,
//! so the tolerance is the accumulation over the CG iterations of the
//! reduction-shape/FMA-order difference, declared looser than the per-op ECT.
//!
//! SCOPE: N_θ = 1 box, free (zero-flux Neumann) BCs, constant/varying transport
//! read as the two-cell face average. The RHS `b`-assembly, Robin-Robin wall
//! coupling, solid-conduction CG, combustion, real HDF5 TableEos, and
//! cut/mixed-N_θ geometry are the S13c split (recorded in the SESSION_LOG).
use crucible_grid::{Grid, GridSpec};
use crucible_solvers::gas_diffusion::{FaceGasBc, GasDiffBcs, GasDiffusion, GasTransportField};
use crucible_solvers::transport::TransportProps;

/// Declared cross-device ECT for the converged CG solve (META-1 §2.5). The
/// per-op class-A RHS held 1e-9; a CG solve accumulates the reduction-shape/
/// FMA-order difference over its iterations to a tight residual floor, so the
/// solve tolerance is declared one decade looser.
const ECT_SOLVE: f64 = 1e-8;

unsafe extern "C" {
    fn gpu_class_d_cg(
        x: *mut f64,
        b: *const f64,
        rho: *const f64,
        cv: *const f64,
        mu: *const f64,
        kk: *const f64,
        rhod: *const f64,
        gas: *const f64,
        comp: i32,
        n_r: i32,
        n_z: i32,
        r_min: f64,
        dr: f64,
        dz: f64,
        wq: f64,
        fixed_iters: i32,
        iters_out: *mut i32,
        resid_out: *mut f64,
    ) -> i32;
}

/// Spatially-varying, nowhere-symmetric transport scaling — a face coefficient
/// that used one cell's value instead of the two-cell average would show up
/// (the S4 blindness test's logic, now cross-device).
fn scale(i_r: usize, i_z: usize) -> f64 {
    1.0 + 0.31 * (i_r as f64) + 0.17 * (i_z as f64)
}

fn main() {
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

    // Varying transport over every gas cell (base μ,k,ρD × scale; cv constant).
    const MU: f64 = 0.01;
    const K: f64 = 0.02;
    const RHOD: f64 = 0.0125;
    const CV: f64 = 720.0;
    let mut tr = GasTransportField::alloc(&g);
    g.for_each_active_cell(|cell| {
        let f = scale(cell.i_r, cell.i_z);
        let props = TransportProps {
            mu: MU * f,
            k: K * f,
            cp: 1005.0,
            cp_film: 1005.0,
            cv: CV,
            rho_d: RHOD * f,
            dh_dz: 0.0,
            pr: 0.71,
        };
        // N_θ = 1 ⇒ the plane index is the brick-local index.
        tr.set(cell.bi, cell.idx, cell.i_r, cell.i_z, &props)
            .expect("valid transport");
    });

    let op = GasDiffusion::new(GasDiffBcs {
        r_inner: FaceGasBc::free(),
        r_outer: FaceGasBc::free(),
        z_lo: FaceGasBc::free(),
        z_hi: FaceGasBc::free(),
    });

    // The CG problem: a positive varying density, a nonzero smooth RHS b, a
    // smooth nonzero initial x0 (x0 only shifts the result: δ solves the
    // system, x_final = x0 + δ — same on CPU and GPU). wq = a representative
    // implicit weight (backward-Euler-class); large enough to make the solve
    // take real iterations across components (velocity/species masses are
    // lighter than the c_v-weighted T mass, so they are stiffer).
    let rho_of = |i_r: usize, i_z: usize| 1.0 + 0.2 * (1.7 * i_r as f64 + 0.9 * i_z as f64).sin();
    let b_of = |i_r: usize, i_z: usize| {
        0.7 * (0.30 * i_r as f64).sin() * (0.21 * i_z as f64).cos()
            + 0.15 * (0.9 * i_r as f64 - 0.5 * i_z as f64).sin()
    };
    let x0_of = |i_r: usize, i_z: usize| 0.4 + 0.1 * (0.5 * i_r as f64 + 0.3 * i_z as f64).cos();
    let wq = 5.0e-2;

    let comp_names = ["Ur", "Uz", "Om", "T", "C"];
    let mut overall_worst = 0.0f64;
    let mut all_bit_identical = true;

    println!("S13b — class-D per-component CG residency (grid {n_r}×{n_z}, N_θ=1 box, wq={wq})");
    println!(
        "  {:>3}  {:>6} {:>6}   {:>11}   {:>11}   {:>10}/{:>10}   {:>9}",
        "cmp", "cpu_it", "gpu_it", "cpu_resid", "gpu_resid", "x_rel", "d_rel", "rerun"
    );

    for comp in 0..5usize {
        let d = op.xcheck_cg_dense(&g, comp, &tr, rho_of, b_of, x0_of, wq);
        let ncell = d.n_r * d.n_z;

        // GPU with the SAME iteration count as the CPU (the clean ECT compare:
        // identical work, only arithmetic order differs).
        let mut xg = d.x0.clone();
        let (mut gi, mut gr) = (0i32, 0.0f64);
        let fail = unsafe {
            gpu_class_d_cg(
                xg.as_mut_ptr(),
                d.b.as_ptr(),
                d.rho.as_ptr(),
                d.cv.as_ptr(),
                d.mu.as_ptr(),
                d.k.as_ptr(),
                d.rhod.as_ptr(),
                d.gas.as_ptr(),
                comp as i32,
                d.n_r as i32,
                d.n_z as i32,
                r_min,
                dr,
                dz,
                wq,
                d.iters as i32,
                &mut gi,
                &mut gr,
            )
        };
        assert_eq!(
            fail, 0,
            "GPU CG ({}) missed EPS acceptance",
            comp_names[comp]
        );

        // GPU with its OWN data-dependent termination (sanity: should land at
        // the same iteration count / residual floor as the CPU, ±1).
        let mut xg_own = d.x0.clone();
        let (mut gi_own, mut gr_own) = (0i32, 0.0f64);
        unsafe {
            gpu_class_d_cg(
                xg_own.as_mut_ptr(),
                d.b.as_ptr(),
                d.rho.as_ptr(),
                d.cv.as_ptr(),
                d.mu.as_ptr(),
                d.k.as_ptr(),
                d.rhod.as_ptr(),
                d.gas.as_ptr(),
                comp as i32,
                d.n_r as i32,
                d.n_z as i32,
                r_min,
                dr,
                dz,
                wq,
                -1,
                &mut gi_own,
                &mut gr_own,
            );
        }

        // Same-build rerun bit-identity (fixed-iters path).
        let mut xg2 = d.x0.clone();
        let (mut gi2, mut gr2) = (0i32, 0.0f64);
        unsafe {
            gpu_class_d_cg(
                xg2.as_mut_ptr(),
                d.b.as_ptr(),
                d.rho.as_ptr(),
                d.cv.as_ptr(),
                d.mu.as_ptr(),
                d.k.as_ptr(),
                d.rhod.as_ptr(),
                d.gas.as_ptr(),
                comp as i32,
                d.n_r as i32,
                d.n_z as i32,
                r_min,
                dr,
                dz,
                wq,
                d.iters as i32,
                &mut gi2,
                &mut gr2,
            );
        }
        let bit_identical = xg == xg2;
        all_bit_identical &= bit_identical;

        // Compare over gas cells, both on the deliverable field x_final and on
        // the solver output δ = x_final − x0 (sharper — x0 is identical on both
        // sides, so the x_final scale would otherwise dilute a δ-level error).
        let mut worst_rel = 0.0f64;
        let mut worst_delta_rel = 0.0f64;
        for c in 0..ncell {
            if d.gas[c] == 0.0 {
                continue;
            }
            let (a, b) = (d.x_final[c], xg[c]);
            let scale = a.abs().max(b.abs());
            if scale > 1e-6 {
                worst_rel = worst_rel.max((a - b).abs() / scale);
            }
            let (da, db) = (a - d.x0[c], b - d.x0[c]);
            let dscale = da.abs().max(db.abs());
            if dscale > 1e-9 {
                worst_delta_rel = worst_delta_rel.max((da - db).abs() / dscale);
            }
        }
        overall_worst = overall_worst.max(worst_rel).max(worst_delta_rel);

        println!(
            "  {:>3}  {:>6} {:>6}   {:>11.3e}   {:>11.3e}   {:>10.3e}/{:>10.3e}   {:>9}",
            comp_names[comp],
            d.iters,
            gi_own,
            d.resid,
            gr_own,
            worst_rel,
            worst_delta_rel,
            if bit_identical {
                "BIT-IDENT"
            } else {
                "*DIFFERS*"
            }
        );

        assert!(
            worst_rel < ECT_SOLVE && worst_delta_rel < ECT_SOLVE,
            "CPU↔GPU class-D CG ({}) diverged beyond ECT {ECT_SOLVE:.0e}: x {worst_rel:.3e}, δ {worst_delta_rel:.3e}",
            comp_names[comp]
        );
        assert!(
            bit_identical,
            "GPU CG ({}) rerun not deterministic",
            comp_names[comp]
        );
    }

    println!("  overall worst rel = {overall_worst:.3e}   (declared ECT {ECT_SOLVE:.0e})");
    println!(
        "  GPU same-build reruns: {}",
        if all_bit_identical {
            "ALL BIT-IDENTICAL"
        } else {
            "*** SOME DIFFER ***"
        }
    );
    assert!(overall_worst < ECT_SOLVE);
    assert!(all_bit_identical);
    println!("ALL PASS");
}
