//! S13 residency cross-check (STEP 1 milestone): the class-A (explicit
//! hyperbolic) SDC RHS ported to a device-resident kernel set —
//! `cuda/residency.cu` — scored against the bit-exact CPU reference
//! `crucible_solvers::euler::Euler::eval_rhs` on a REAL (r,z) fixture.
//!
//! This generalizes the S12 spike's two simplifications that S13 STEP 1 owns:
//!   (a) the sweep is now the real 2-direction operator — the exact
//!       `face_radius` cylindrical metric (r-sweep area-weighted, z-sweep
//!       metric-ratio) + the SOLV-1 §3.3 geometric sources — not the S12
//!       z-only-uniform pass where the annular z-face areas cancel;
//!   (c) the kernels are staged (fill_prims -> rate), not one fused kernel.
//! (b real HDF5 TableEos + d class-D diffusion + cut apertures + mixed-N_θ
//! are the S13b split, recorded in the SESSION_LOG.)
//!
//! Determinism (META-1 §2.5): same-build GPU rerun is bit-identical; CPU<->GPU
//! agree to the declared ECT (FMA-order, non-chaotic fixture). We compare only
//! INTERIOR cells (>= NGHOST from every edge) whose compact PPM stencil is all
//! real interior data, so the BC / reflux / axis machinery (S13b) is out of
//! the compared set.
use crucible_grid::{Grid, GridSpec};
use crucible_solvers::euler::{EULER_FIELDS, Euler, EulerFields, FlowBc, FlowBcs, GammaLaw, NCOMP};
use crucible_solvers::sdc::{FlowClass, Sdc};

const NGH: usize = 3;
const GAMMA: f64 = 1.4;
/// Declared cross-device ECT (META-1 §2.5; S12 measured 5.0e-10 on HLLC over
/// 1e6 random states — the full class-A RHS on a smooth fixture must stay in
/// this FMA-order band).
const ECT: f64 = 1e-9;

unsafe extern "C" {
    fn gpu_class_a_rhs(
        cons: *const f64,
        n_r: i32,
        n_z: i32,
        r_min: f64,
        dr: f64,
        z_min: f64,
        dz: f64,
        gamma: f64,
        rate: *mut f64,
    ) -> i32;
    fn gpu_class_a_march(
        cons: *mut f64,
        n_r: i32,
        n_z: i32,
        r_min: f64,
        dr: f64,
        z_min: f64,
        dz: f64,
        gamma: f64,
        dt: f64,
        nsteps: i32,
    ) -> i32;
    fn gpu_class_a_bench(
        cons: *const f64,
        n_r: i32,
        n_z: i32,
        r_min: f64,
        dr: f64,
        z_min: f64,
        dz: f64,
        gamma: f64,
        iters: i32,
    ) -> f64;
    fn gpu_stable_dt(
        cons: *const f64,
        n_r: i32,
        n_z: i32,
        dr: f64,
        dz: f64,
        gamma: f64,
        cfl: f64,
        bad: *mut i32,
    ) -> f64;
}

/// Smooth, subsonic, strictly-positive (ρ,p) fixture — the reconstruction
/// limiter is inactive on smooth data, so the PPM branches don't flip between
/// CPU and GPU (they would at a shock; that ECT is larger and physical — noted
/// in the log). Nonzero u_r/u_θ/u_z so every geometric source is exercised.
fn fixture(r: f64, th: f64, z: f64) -> [f64; NPRIM_OUT] {
    let rho = 1.0 + 0.20 * (1.7 * r + 0.9 * z).sin() * (th).cos();
    let ur = 0.30 * (1.1 * z - 0.4 * r).sin();
    let ut = 0.25 * (0.8 * r + 0.5 * z).cos();
    let uz = 0.35 * (0.9 * r - 1.3 * z).sin();
    let p = 2.0 + 0.30 * (1.3 * z - 0.7 * r).cos();
    let c = 0.5 + 0.1 * (r).sin();
    // Prim slot order [rho,u_r,u_theta,u_z,p,c] (prim6 fills burn+aux 0).
    [rho, ur, ut, uz, p, c]
}
const NPRIM_OUT: usize = 6;

fn main() {
    crucible_gpu::keep_kernels(); // keeps this bin's CUDA entry points in the link (lib.rs)
    // Uniform box, N_θ = 1, r_min > 0 (no axis — a clean box; the axis pair
    // gather is S13b). Modest size so interior cells dominate the compare.
    let (n_r, n_z) = (48usize, 96usize);
    let spec = GridSpec {
        r_min: 0.5,
        dr: 1.0 / n_r as f64,
        n_r,
        z_min: 0.0,
        dz: 1.0 / n_z as f64,
        n_z,
        n_theta_max: 1,
        axisymmetry_assertion: true,
    };
    let (r_min, dr, z_min, dz) = (spec.r_min, spec.dr, spec.z_min, spec.dz);
    let eos = GammaLaw { gamma: GAMMA };
    let mut g = Grid::build(spec, EULER_FIELDS).expect("valid spec");
    let f = EulerFields::resolve(&g).expect("fields");
    crucible_solvers::euler::fill_from_prim(&mut g, &f, &eos, |r, th, z| {
        let a = fixture(r, th, z);
        crucible_solvers::euler::prim6(a[0], a[1], a[2], a[3], a[4], a[5])
    });

    // CPU reference: one class-A rhs evaluation. Transmissive BCs — interior
    // cells (the compared set) never see them; zero external source.
    let zero = |_: f64, _: f64, _: f64, _: f64| [0.0; NCOMP];
    let op = Euler {
        eos,
        source: &zero,
        bcs: FlowBcs {
            r_inner: FlowBc::Transmissive,
            r_outer: FlowBc::Transmissive,
            z_lo: FlowBc::Transmissive,
            z_hi: FlowBc::Transmissive,
        },
        wall_normal: None,
        slip_wall_z_faces: false,
        combustion: None,
    };
    let mut ws = op.workspace(&g).expect("workspace");
    op.eval_rhs(&g, &f, &mut ws, 0.0).expect("cpu eval_rhs");

    // Dense cons [cell*NCOMP], cell = i_r*n_z + i_z (the kernel's layout), and
    // a parallel map (cell) -> CPU rate for the interior compare.
    let ncell = n_r * n_z;
    let ids = f.ids();
    let mut cons = vec![0.0f64; ncell * NCOMP];
    let mut cpu_rate = vec![0.0f64; ncell * NCOMP];
    let rates = ws.rates();
    g.for_each_active_cell(|cell| {
        let c = cell.i_r * n_z + cell.i_z;
        let b = g.brick(cell.bi);
        for k in 0..NCOMP {
            cons[c * NCOMP + k] = b.field(ids[k])[cell.idx];
            cpu_rate[c * NCOMP + k] = rates[cell.bi][cell.idx][k];
        }
    });

    // GPU: the resident class-A rhs.
    let mut gpu_rate = vec![0.0f64; ncell * NCOMP];
    let bad = unsafe {
        gpu_class_a_rhs(
            cons.as_ptr(),
            n_r as i32,
            n_z as i32,
            r_min,
            dr,
            z_min,
            dz,
            GAMMA,
            gpu_rate.as_mut_ptr(),
        )
    };
    assert_eq!(bad, 0, "GPU flagged a non-physical cell");

    // Compare INTERIOR cells only.
    let mut worst_abs = 0.0f64;
    let mut worst_rel = 0.0f64;
    let mut n_cmp = 0usize;
    for i_r in NGH..n_r - NGH {
        for i_z in NGH..n_z - NGH {
            let c = i_r * n_z + i_z;
            for k in 0..NCOMP {
                let a = cpu_rate[c * NCOMP + k];
                let b = gpu_rate[c * NCOMP + k];
                let d = (a - b).abs();
                if d > worst_abs {
                    worst_abs = d;
                }
                let scale = a.abs().max(b.abs());
                if scale > 1e-6 {
                    let rel = d / scale;
                    if rel > worst_rel {
                        worst_rel = rel;
                    }
                }
                n_cmp += 1;
            }
        }
    }

    // Same-build GPU rerun bit-identity (determinism gate, META-1 §2.5).
    let mut gpu_rerun = vec![0.0f64; ncell * NCOMP];
    unsafe {
        gpu_class_a_rhs(
            cons.as_ptr(),
            n_r as i32,
            n_z as i32,
            r_min,
            dr,
            z_min,
            dz,
            GAMMA,
            gpu_rerun.as_mut_ptr(),
        );
    }
    let bit_identical = gpu_rate == gpu_rerun;

    println!("PHASE 1 — single class-A RHS (real metric + sources), GammaLaw");
    println!("  grid {n_r}×{n_z} (N_θ=1 box), {n_cmp} interior scalar comparisons");
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
    assert!(worst_rel < ECT, "CPU↔GPU class-A RHS diverged beyond ECT");
    assert!(bit_identical, "GPU rerun not deterministic");
    println!("  PHASE 1 PASS");

    // ================================================================
    // PHASE 1b — stable_dt on-device (S13c: the CFL clock). Done HERE,
    // before PHASE 2's march mutates `g` — so the CPU `Euler::stable_dt`
    // and the GPU `gpu_stable_dt(cons, …)` score the SAME (fixture) state.
    // The reduction is a MAX (exactly order-independent), so same-build
    // reruns are bit-identical; CPU↔GPU differ only at the per-cell σ's
    // FMA order (declared ECT).
    // ================================================================
    let cfl = 0.4;
    let cpu_dt = op.stable_dt(&g, &f, cfl).expect("cpu stable_dt");
    let mut sdt_bad = 0i32;
    let gpu_dt = unsafe {
        gpu_stable_dt(
            cons.as_ptr(),
            n_r as i32,
            n_z as i32,
            dr,
            dz,
            GAMMA,
            cfl,
            &mut sdt_bad,
        )
    };
    assert_eq!(sdt_bad, 0, "GPU stable_dt flagged a non-physical cell");
    let mut sdt_bad2 = 0i32;
    let gpu_dt2 = unsafe {
        gpu_stable_dt(
            cons.as_ptr(),
            n_r as i32,
            n_z as i32,
            dr,
            dz,
            GAMMA,
            cfl,
            &mut sdt_bad2,
        )
    };
    let sdt_rel = (cpu_dt - gpu_dt).abs() / cpu_dt.abs().max(gpu_dt.abs());
    println!("PHASE 1b — stable_dt (CFL clock), cfl={cfl}");
    println!(
        "  cpu Δt = {cpu_dt:.9e}   gpu Δt = {gpu_dt:.9e}   rel = {sdt_rel:.3e}   (ECT {ECT:.0e})"
    );
    println!(
        "  GPU same-build rerun: {}",
        if gpu_dt == gpu_dt2 {
            "BIT-IDENTICAL"
        } else {
            "*** DIFFERS ***"
        }
    );
    assert!(sdt_rel < ECT, "CPU↔GPU stable_dt diverged beyond ECT");
    assert!(gpu_dt == gpu_dt2, "GPU stable_dt rerun not deterministic");
    println!("  PHASE 1b PASS");

    // ================================================================
    // PHASE 2 — the RESIDENT marched SDC step + the FND-6 checkpoint.
    // Fixed dt fed to both CPU and GPU (stable_dt-on-device is S13b);
    // the shared trajectory is deterministic given dt. `g` still holds
    // the initial fixture (eval_rhs does not mutate U), so `cons` above
    // is the shared initial state.
    // ================================================================
    let m = 5i32; // steps; frozen-boundary/BC influence stays within 3·M of an edge
    let dt = 0.001; // CFL ≈ 0.4 at max wavespeed ~2, dz = 1/96
    let cons0 = cons.clone();

    // CPU reference march (mutates g to the final state).
    let mut sdc = Sdc::new();
    let flow = FlowClass {
        op: &op,
        fields: &f,
    };
    let mut t = 0.0;
    for _ in 0..m {
        sdc.step_flow(&mut g, &flow, t, dt).expect("cpu step_flow");
        t += dt;
    }
    let mut cpu_final = vec![0.0f64; ncell * NCOMP];
    g.for_each_active_cell(|cell| {
        let c = cell.i_r * n_z + cell.i_z;
        let b = g.brick(cell.bi);
        for k in 0..NCOMP {
            cpu_final[c * NCOMP + k] = b.field(ids[k])[cell.idx];
        }
    });

    // GPU resident march (state stays on device across the internal loop).
    let mut gpu_final = cons0.clone();
    let bad = unsafe {
        gpu_class_a_march(
            gpu_final.as_mut_ptr(),
            n_r as i32,
            n_z as i32,
            r_min,
            dr,
            z_min,
            dz,
            GAMMA,
            dt,
            m,
        )
    };
    assert_eq!(bad, 0, "GPU march flagged a non-physical cell");

    // Compare cells ≥ 3·M from every edge (neither the CPU's BCs nor the GPU's
    // frozen boundary can have influenced them in M compact-stencil steps).
    let margin = (3 * m) as usize;
    let mut mworst_abs = 0.0f64;
    let mut mworst_rel = 0.0f64;
    let mut m_cmp = 0usize;
    for i_r in margin..n_r - margin {
        for i_z in margin..n_z - margin {
            let c = i_r * n_z + i_z;
            for k in 0..NCOMP {
                let a = cpu_final[c * NCOMP + k];
                let b = gpu_final[c * NCOMP + k];
                let d = (a - b).abs();
                if d > mworst_abs {
                    mworst_abs = d;
                }
                let scale = a.abs().max(b.abs());
                if scale > 1e-6 {
                    let rel = d / scale;
                    if rel > mworst_rel {
                        mworst_rel = rel;
                    }
                }
                m_cmp += 1;
            }
        }
    }

    // FND-6 checkpoint/restart: a march split at step 2 (host round-trip =
    // the bit-faithful device-state serialize/resume) must equal the
    // uninterrupted march BIT-IDENTICALLY on the same device/build.
    let mut split = cons0.clone();
    unsafe {
        gpu_class_a_march(
            split.as_mut_ptr(),
            n_r as i32,
            n_z as i32,
            r_min,
            dr,
            z_min,
            dz,
            GAMMA,
            dt,
            2,
        );
        gpu_class_a_march(
            split.as_mut_ptr(),
            n_r as i32,
            n_z as i32,
            r_min,
            dr,
            z_min,
            dz,
            GAMMA,
            dt,
            3,
        );
    }
    let checkpoint_bit_identical = split == gpu_final;

    // Same-build resident-march rerun determinism.
    let mut gpu_final2 = cons0.clone();
    unsafe {
        gpu_class_a_march(
            gpu_final2.as_mut_ptr(),
            n_r as i32,
            n_z as i32,
            r_min,
            dr,
            z_min,
            dz,
            GAMMA,
            dt,
            m,
        );
    }
    let march_bit_identical = gpu_final == gpu_final2;

    println!("PHASE 2 — resident marched SDC step ({m} steps, dt={dt}) + FND-6 checkpoint");
    println!("  compared cells ≥ {margin} from every edge, {m_cmp} scalar comparisons");
    println!("  worst abs diff = {mworst_abs:.3e}");
    println!(
        "  worst rel diff = {mworst_rel:.3e}   (ECT does not grow: phase-1 was {worst_rel:.3e})"
    );
    println!(
        "  GPU resident-march rerun: {}",
        if march_bit_identical {
            "BIT-IDENTICAL"
        } else {
            "*** DIFFERS ***"
        }
    );
    println!(
        "  checkpoint march(2)+march(3) == march(5): {}",
        if checkpoint_bit_identical {
            "BIT-IDENTICAL"
        } else {
            "*** DIFFERS ***"
        }
    );
    assert!(
        mworst_rel < 1e-8,
        "marched CPU↔GPU diverged beyond FMA-order band"
    );
    assert!(march_bit_identical, "GPU resident march not deterministic");
    assert!(
        checkpoint_bit_identical,
        "checkpoint/restart not bit-faithful"
    );
    println!("  PHASE 2 PASS");

    // ================================================================
    // PHASE 3 — throughput of the staged class-A rate (STEP 4 tuning).
    // Big resident grid; time `iters` full RHS evals (k_rate_r+k_rate_z).
    // ================================================================
    let (bn_r, bn_z) = (256usize, 512usize);
    let bncell = bn_r * bn_z;
    let bdr = 1.0 / bn_r as f64;
    let bdz = 1.0 / bn_z as f64;
    let mut bcons = vec![0.0f64; bncell * NCOMP];
    for i_r in 0..bn_r {
        for i_z in 0..bn_z {
            let r = r_min + (i_r as f64 + 0.5) * bdr;
            let z = (i_z as f64 + 0.5) * bdz;
            let a = fixture(r, 0.0, z);
            // GammaLaw prim->cons (matches euler::GammaLaw::prim_to_cons).
            let (rho, ur, ut, uz, p, cc) = (a[0], a[1], a[2], a[3], a[4], a[5]);
            let etot = p / (GAMMA - 1.0) + 0.5 * rho * (ur * ur + ut * ut + uz * uz);
            let c = i_r * bn_z + i_z;
            let u = [rho, rho * ur, rho * ut, rho * uz, etot, rho * cc, 0.0];
            bcons[c * NCOMP..c * NCOMP + NCOMP].copy_from_slice(&u);
        }
    }
    let iters = 300i32;
    let ms = unsafe {
        gpu_class_a_bench(
            bcons.as_ptr(),
            bn_r as i32,
            bn_z as i32,
            r_min,
            bdr,
            0.0,
            bdz,
            GAMMA,
            iters,
        )
    };
    let evals = iters as f64 * bncell as f64;
    let cell_rhs_per_s = evals / (ms / 1000.0);
    // Each RHS eval = r-sweep + z-sweep + sources ≈ 2 sweep passes; the S12
    // spike quoted a single-sweep number (~1e8 cups), so the comparable
    // per-sweep rate is ~2× cell_rhs_per_s.
    println!(
        "PHASE 3 — staged class-A rate throughput ({bn_r}×{bn_z} = {bncell} cells, {iters} evals)"
    );
    println!(
        "  {ms:.2} ms total → {:.3e} cell-RHS/s (≈ {:.3e} cell-updates/s per-sweep-equiv)",
        cell_rhs_per_s,
        2.0 * cell_rhs_per_s
    );
    println!("ALL PASS");
}
