//! S13c residency cross-check: the class-D diffusion FORCING (the affine
//! viscous-stress physics — `fill_lag_gradients` + `assemble_rates`) ported to
//! a device kernel set (`cuda/residency_diffusion.cu`: `kd_lag_grads` +
//! `kd_assemble_rates`, FFI `gpu_class_d_assemble`) scored against the bit-exact
//! CPU reference `GasDiffusion::{fill_lag_gradients, assemble_rates}` (via the
//! additive doc-hidden `xcheck_assemble_dense` accessor) on a real (r,z) fixture.
//!
//! This is the RHS-builder that feeds the S13b class-D CG: the affine rate
//! `R_affine(x₀)` whose `wqnew·(unit·(R_affine − d_lag))` is the CG's `b`. It is
//! the intricate cross-term core of F_visc — compressible viscous stress
//! (τ_rr/τ_zz/τ_rz/τ_θθ), Fourier conduction, species diffusion + its enthalpy
//! flux. A pure GATHER kernel (no reduction, no Picard), so — like the S13
//! class-A RHS — CPU↔GPU is a per-cell FMA-order tolerance.
//!
//! SCOPE: N_θ = 1 box, free (FreeSlip/Adiabatic/ZeroFlux) BCs; INTERIOR cells
//! only (≥ 3 from every edge so the lag-gradient + face stencil is all real
//! interior data). Real NoSlip/Robin walls, the θ-stress tensor, cut apertures,
//! and the full Picard→CG wiring are the remaining S13c legs.
use crucible_grid::{Grid, GridSpec};
use crucible_solvers::gas_diffusion::{FaceGasBc, GasDiffBcs, GasDiffusion, GasTransportField};
use crucible_solvers::transport::TransportProps;

const NGH: usize = 3;
const NCOMP: usize = 7;
/// Declared cross-device ECT (META-1 §2.5): the affine rate is a single gather
/// evaluation (like the class-A RHS, which held 1.2e-10), so the same FMA-order
/// band applies.
const ECT: f64 = 1e-9;

unsafe extern "C" {
    fn gpu_class_d_assemble(
        s_ur: *const f64,
        s_om: *const f64,
        s_uz: *const f64,
        s_tt: *const f64,
        s_cc: *const f64,
        lag_ur: *const f64,
        lag_uz: *const f64,
        mu: *const f64,
        kk: *const f64,
        rhod: *const f64,
        dhdz: *const f64,
        gas: *const f64,
        n_r: i32,
        n_z: i32,
        r_min: f64,
        dr: f64,
        dz: f64,
        rate: *mut f64,
    );
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

    // Varying, nowhere-symmetric transport (a one-sided face coefficient would
    // show), plus a nonzero ∂h/∂Z so the species-enthalpy flux limb is live
    // (the constant occupant zeroes it — this fixture exercises it directly).
    const MU: f64 = 0.01;
    const K: f64 = 0.02;
    const RHOD: f64 = 0.0125;
    let mut tr = GasTransportField::alloc(&g);
    g.for_each_active_cell(|cell| {
        let f = 1.0 + 0.31 * (cell.i_r as f64) + 0.17 * (cell.i_z as f64);
        let props = TransportProps {
            mu: MU * f,
            k: K * f,
            cp: 1005.0,
            cp_film: 1005.0,
            cv: 720.0,
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

    // Smooth, nonzero sol (current iterate) and a DISTINCT lag state (the
    // cross-term lag ≠ sol — so the lagged pieces are genuinely exercised, not
    // silently equal to the implicit ones). Positive T (assemble refuses T≤0).
    let sol_ur = |r: usize, z: usize| 0.30 * (1.1 * z as f64 - 0.4 * r as f64).sin();
    let sol_om = |r: usize, z: usize| 0.25 * (0.8 * r as f64 + 0.5 * z as f64).cos();
    let sol_uz = |r: usize, z: usize| 0.35 * (0.9 * r as f64 - 1.3 * z as f64).sin();
    let sol_tt = |r: usize, z: usize| 300.0 + 20.0 * (0.3 * r as f64 + 0.2 * z as f64).sin();
    let sol_cc = |r: usize, z: usize| 0.5 + 0.1 * (0.25 * r as f64).sin() * (0.2 * z as f64).cos();
    // Lag = sol perturbed (a different Picard iterate).
    let lag_ur = |r: usize, z: usize| sol_ur(r, z) + 0.05 * (0.7 * r as f64 + 0.3 * z as f64).cos();
    let lag_uz = |r: usize, z: usize| sol_uz(r, z) - 0.04 * (0.5 * r as f64 - 0.6 * z as f64).sin();

    let d = op.xcheck_assemble_dense(
        &g, &tr, sol_ur, sol_om, sol_uz, sol_tt, sol_cc, lag_ur, lag_uz,
    );
    let ncell = n_r * n_z;

    let mut gpu_rate = vec![0.0f64; ncell * NCOMP];
    unsafe {
        gpu_class_d_assemble(
            d.ur.as_ptr(),
            d.om.as_ptr(),
            d.uz.as_ptr(),
            d.tt.as_ptr(),
            d.cc.as_ptr(),
            d.lag_ur.as_ptr(),
            d.lag_uz.as_ptr(),
            d.mu.as_ptr(),
            d.k.as_ptr(),
            d.rhod.as_ptr(),
            d.dhdz.as_ptr(),
            d.gas.as_ptr(),
            n_r as i32,
            n_z as i32,
            r_min,
            dr,
            dz,
            gpu_rate.as_mut_ptr(),
        );
    }

    // Same-build rerun bit-identity.
    let mut gpu_rate2 = vec![0.0f64; ncell * NCOMP];
    unsafe {
        gpu_class_d_assemble(
            d.ur.as_ptr(),
            d.om.as_ptr(),
            d.uz.as_ptr(),
            d.tt.as_ptr(),
            d.cc.as_ptr(),
            d.lag_ur.as_ptr(),
            d.lag_uz.as_ptr(),
            d.mu.as_ptr(),
            d.k.as_ptr(),
            d.rhod.as_ptr(),
            d.dhdz.as_ptr(),
            d.gas.as_ptr(),
            n_r as i32,
            n_z as i32,
            r_min,
            dr,
            dz,
            gpu_rate2.as_mut_ptr(),
        );
    }
    let bit_identical = gpu_rate == gpu_rate2;

    // Compare interior cells, per component, so a single wrong term is visible.
    let names = ["rho", "m_r", "m_th", "m_z", "en", "rho_c", "rho_b"];
    let mut worst_rel = [0.0f64; NCOMP];
    let mut worst_abs = [0.0f64; NCOMP];
    let mut n_cmp = 0usize;
    for i_r in NGH..n_r - NGH {
        for i_z in NGH..n_z - NGH {
            let c = i_r * n_z + i_z;
            for k in 0..NCOMP {
                let a = d.rate[c][k];
                let b = gpu_rate[c * NCOMP + k];
                let ad = (a - b).abs();
                if ad > worst_abs[k] {
                    worst_abs[k] = ad;
                }
                let scale = a.abs().max(b.abs());
                if scale > 1e-6 {
                    worst_rel[k] = worst_rel[k].max(ad / scale);
                }
            }
            n_cmp += 1;
        }
    }

    println!("S13c — class-D affine FORCING (fill_lag_gradients + assemble_rates)");
    println!("  grid {n_r}×{n_z} (N_θ=1 box), {n_cmp} interior cells compared");
    println!(
        "  {:>5}   {:>11}   {:>11}",
        "comp", "worst_abs", "worst_rel"
    );
    let mut overall = 0.0f64;
    for k in 0..NCOMP {
        println!(
            "  {:>5}   {:>11.3e}   {:>11.3e}",
            names[k], worst_abs[k], worst_rel[k]
        );
        overall = overall.max(worst_rel[k]);
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
        "CPU↔GPU class-D affine rate diverged beyond ECT: {overall:.3e}"
    );
    assert!(bit_identical, "GPU assemble rerun not deterministic");
    println!("ALL PASS");
}
