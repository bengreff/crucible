//! S13c CLOSE cross-check: the COMPOSED full-physics 3-D resident engine step
//! (`cuda/residency_engine.cu`, host bridge `engine_host`) scored against the
//! CPU production path — `Euler::eval_rhs` (blend EOS + combustion + igniter +
//! inflow/outflow) and `Sdc::step` (flow + reaction classes, stagewise SRD, the
//! COUP-2 audit) — on the REAL `configs/rl10_startup_3d.toml` assembly (the
//! ◆C3 world), tables of record, the run's own operator construction
//! (`engine_host::Schedule` = `crucible_engine::run`'s schedule, one owner).
//!
//! The CPU pre-marches the real start (spark at the config's own window) to a
//! target time (XCHECK_T1, default 3 ms) so the cross-check sits in the
//! reacting regime. PHASE 1: one full RHS (+ ledger). PHASE 1b: `stable_dt`
//! (θ-arc + front carrier). PHASE 2: M2 coupled steps — per-step audit rows
//! (the CPU's own `Sdc` audit arithmetic on the device's reduced operands vs
//! the CPU report), the final state, rerun bit-identity, the FND-6 checkpoint
//! split (cons + the warm-start hints). PHASE 3: throughput.
#![allow(clippy::needless_range_loop)] // dense (j, rz, k) index arithmetic mirrors the device layout
use crucible_engine::eos_sel::ChemEos;
use crucible_gpu::engine_host::{DeviceEngine, Schedule, load_run};
use crucible_solvers::euler::{
    Combustion, Cons, EosLaw, Euler, FlowBc, FlowBcs, I_RB, IgnitionColumns, NCOMP, NPRIM,
    THETA_CELLS, reacting_measure,
};
use crucible_solvers::sdc::{FlowClass, ReactionClass, Sdc};

/// The composed step's declared cross-device ECT (component-scaled): the
/// (p,h,Z) projection converges to EPS_P_PROJECTION = 1e-11 relative on BOTH
/// sides but by FMA-distinct Illinois paths, and the well-balanced θ-momentum
/// at the axis ring divides that by a ~10⁴ cancellation (the θ-flux difference
/// × 1/(r̄Δθ)) — measured 9.5e-9 on the reacting ◆C3 state. Declared 1e-7.
const ECT: f64 = 1e-7;
/// The prims themselves (no cancellation): the projection-tolerance class.
const ECT_PRIM: f64 = 1e-9;
/// Per-cell rail, applied only where the value is ≥ 1e-6 of the component max
/// (below that a well-balance residual's sign is rounding, not physics).
const CANCEL_BOUND: f64 = 1e-6;
const CONFIG: &str = "configs/rl10_startup_3d.toml";

fn main() {
    let root = format!("{}/../..", env!("CARGO_MANIFEST_DIR"));
    std::env::set_current_dir(&root).expect("repo root");
    let m1: usize = std::env::var("XCHECK_M1")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(20000);
    let t1_target: f64 = std::env::var("XCHECK_T1")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(3.0e-3);
    let m2: usize = std::env::var("XCHECK_M2")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(10);

    // --- The real assembly + occupants (engine_host = run.rs's construction).
    let lr = load_run(CONFIG).expect("load");
    let mut spec = lr.spec;
    let eos = ChemEos::bind_blend(&lr.unburnt, &lr.table, spec.injector.h_offset_j_per_kg)
        .expect("blend binds");
    let blend = eos.blend().expect("blend");
    let bl = spec.blend.as_ref().expect("blend spec");
    let comb = Combustion {
        blend,
        ignition: IgnitionColumns::bind(&lr.ignition).expect("ignition binds"),
        wrinkling: bl.wrinkling,
        theta: THETA_CELLS,
    };
    let ids = spec.fields.ids();
    let u_fill: Cons = eos
        .fill_cons(
            spec.fill_p_pa,
            spec.injector.h_inj_j_per_kg,
            spec.injector.z_frac,
        )
        .expect("fill state");
    for (k, &v) in u_fill.iter().enumerate() {
        spec.grid.fill_field(ids[k], move |_, _, _| v);
    }
    let sched = Schedule::new(&spec, &eos, &u_fill).expect("schedule");
    let dev = DeviceEngine::new(&spec, blend, &comb, &sched).expect("device engine");
    let lay = dev.layout;
    let (n, nrz, nt, n_z) = (lay.n, lay.nrz, lay.nt, lay.n_z);
    let act = dev.act.clone();

    // The CPU operator, exactly as run.rs builds it (the schedule's own closures).
    let source_fn = sched.source_fn();
    let contour = spec.contour.clone();
    let normal_fn = move |r: f64, z: f64| contour.wall_normal(r, z);
    let sched_p = sched.clone();
    let pump = move |t: f64| sched_p.p_amb(t);
    let (h_total, c_frac) = (sched.h_total, sched.c_frac);
    let mut op = Euler {
        eos: eos.clone(),
        source: &source_fn,
        bcs: FlowBcs {
            r_inner: FlowBc::Reflecting,
            r_outer: FlowBc::Reflecting,
            z_lo: FlowBc::MassFlowInflow {
                mdot_per_area: 0.0,
                h_total,
                c_frac,
            },
            z_hi: FlowBc::PressureOutflow(&pump),
        },
        wall_normal: Some(&normal_fn),
        slip_wall_z_faces: true,
        combustion: Some(&comb),
    };
    let reaction = ReactionClass { op: &comb };
    let mut sdc = Sdc::with_audit(sched.audit);
    let sdc_gpu = Sdc::with_audit(sched.audit);
    println!("S13c CLOSE — the COMPOSED full-physics 3-D resident engine step on `{CONFIG}`");
    println!(
        "  {}×{}×{} = {n} cells, {} active; t_final {:.4e} s, ramp {:.4e} s, spark at {:.3e} s \
         (q0 = {:.3e} W/m³); {} SRD small cells; pre-march to {t1_target:.3e} s",
        lay.n_r,
        lay.n_z,
        nt,
        sched.n_cells,
        sched.t_final,
        sched.t_ramp_window,
        sched.igniter.t_on,
        sched.igniter.q0,
        dev.n_small()
    );

    let compare = |label: &str, cpu: &[f64], gpu: &[f64], stride: usize| -> (f64, f64) {
        let mut comp_max = vec![0.0f64; stride];
        for j in 0..nt {
            for rz in 0..nrz {
                if act[rz] == 0 {
                    continue;
                }
                let cc = j * nrz + rz;
                for k in 0..stride {
                    comp_max[k] = comp_max[k].max(cpu[cc * stride + k].abs());
                }
            }
        }
        let (mut worst_abs, mut worst_rel, mut worst_scaled, mut n_cmp) =
            (0.0f64, 0.0f64, 0.0f64, 0usize);
        let mut at = (0usize, 0usize, 0.0f64, 0.0f64);
        let mut at_rel = (0usize, 0usize, 0.0f64, 0.0f64);
        for j in 0..nt {
            for rz in 0..nrz {
                if act[rz] == 0 {
                    continue;
                }
                let cc = j * nrz + rz;
                for k in 0..stride {
                    let a = cpu[cc * stride + k];
                    let b = gpu[cc * stride + k];
                    let d = (a - b).abs();
                    worst_abs = worst_abs.max(d);
                    if comp_max[k] > 0.0 && d / comp_max[k] > worst_scaled {
                        worst_scaled = d / comp_max[k];
                        at = (cc, k, a, b);
                    }
                    let s = a.abs().max(b.abs());
                    if s > 1e-6 * comp_max[k] && d / s > worst_rel {
                        worst_rel = d / s;
                        at_rel = (cc, k, a, b);
                    }
                    n_cmp += 1;
                }
            }
        }
        let (cc, k, a, b) = at;
        println!(
            "  {label}: {n_cmp} scalars; worst abs {worst_abs:.3e}; worst component-scaled rel {worst_scaled:.3e} \
             at (j={}, i_r={}, i_z={}, k={k}: cpu {a:.9e} gpu {b:.9e}, comp max {:.3e})",
            cc / nrz,
            (cc % nrz) / n_z,
            cc % n_z,
            comp_max[k]
        );
        let (cc, k, a, b) = at_rel;
        println!(
            "    worst per-cell rel {worst_rel:.3e} at (j={}, i_r={}, i_z={}, k={k}: cpu {a:.9e} gpu {b:.9e})",
            cc / nrz,
            (cc % nrz) / n_z,
            cc % n_z
        );
        (worst_scaled, worst_rel)
    };

    // --- CPU pre-march to t1_target (the real start) --------------------------
    let mut t = 0.0f64;
    let started = std::time::Instant::now();
    let mut m1_done = 0usize;
    if let Ok(path) = std::env::var("XCHECK_LOAD") {
        let bytes = std::fs::read(&path).expect("load state");
        let vals: Vec<f64> = bytes
            .chunks_exact(8)
            .map(|c| f64::from_le_bytes(c.try_into().unwrap()))
            .collect();
        t = vals[0];
        lay.write_dense(&mut spec.grid, &ids, &vals[1..]);
        println!(
            "  loaded the pre-marched state from {path} (t = {t:.4e} s) — the CPU workspace is cold"
        );
    }
    while t < t1_target && m1_done < m1 {
        op.bcs.z_lo = FlowBc::MassFlowInflow {
            mdot_per_area: sched.mdot_per_area(t),
            h_total,
            c_frac,
        };
        let flow = FlowClass {
            op: &op,
            fields: &spec.fields,
        };
        let dt = sdc
            .stable_dt(&spec.grid, &flow, sched.cfl)
            .expect("cpu stable_dt")
            .min(sched.t_final - t);
        sdc.step(
            &mut spec.grid,
            Some(&flow),
            None,
            None,
            None,
            Some(&reaction),
            t,
            dt,
        )
        .expect("cpu step");
        t += dt;
        m1_done += 1;
        if m1_done.is_multiple_of(500) {
            let r = reacting_measure(
                &spec.grid,
                &ids,
                blend,
                &comb.ignition,
                comb.wrinkling,
                comb.theta,
            )
            .unwrap_or(f64::NAN);
            println!(
                "    pre-march step {m1_done}: t = {t:.4e} s, dt = {dt:.3e}, R = {r:.3e} kg/s"
            );
        }
    }
    let t1 = t;
    if let Ok(path) = std::env::var("XCHECK_SAVE") {
        let cons_now = lay.read_dense(&spec.grid, &ids);
        let mut bytes = Vec::with_capacity(8 * (cons_now.len() + 1));
        bytes.extend_from_slice(&t1.to_le_bytes());
        for v in &cons_now {
            bytes.extend_from_slice(&v.to_le_bytes());
        }
        std::fs::write(&path, bytes).expect("save state");
        println!("  saved the pre-marched state to {path}");
    }
    let r_meas = reacting_measure(
        &spec.grid,
        &ids,
        blend,
        &comb.ignition,
        comb.wrinkling,
        comb.theta,
    )
    .expect("R");
    let b_max = {
        let mut m = 0.0f64;
        spec.grid.for_each_active_cell(|c| {
            let b = spec.grid.brick(c.bi);
            m = m.max(b.field(ids[I_RB])[c.idx] / b.field(ids[0])[c.idx]);
        });
        m
    };
    println!(
        "  CPU pre-march: {m1_done} steps to t = {t1:.4e} s in {:.1} s ({:.3} s/step); reacting measure {r_meas:.3e} kg/s, max b {b_max:.3e}",
        started.elapsed().as_secs_f64(),
        started.elapsed().as_secs_f64() / m1_done.max(1) as f64
    );

    // ================================================================
    // PHASE 1 — one full RHS (+ ledger) on the pre-marched state.
    // ================================================================
    let cons1 = lay.read_dense(&spec.grid, &ids);
    let mut ws = op.workspace(&spec.grid).expect("workspace");
    op.bcs.z_lo = FlowBc::MassFlowInflow {
        mdot_per_area: sched.mdot_per_area(t1),
        h_total,
        c_frac,
    };
    op.eval_rhs(&spec.grid, &spec.fields, &mut ws, t1)
        .expect("cpu eval_rhs");
    let cpu_led = *ws.ledger();
    let mut cpu_rate = vec![0.0f64; n * NCOMP];
    let rates = ws.rates();
    spec.grid.for_each_active_cell(|cell| {
        let cc = lay.dense(cell.i_theta as usize, cell.i_r, cell.i_z);
        for k in 0..NCOMP {
            cpu_rate[cc * NCOMP + k] = rates[cell.bi][cell.idx][k];
        }
    });
    dev.upload(&cons1, None);
    let (gpu_rate, gl) = dev.rhs(&sched, t1).expect("gpu rhs");
    // PHASE 0 — the prims themselves (cold projections on both sides): the
    // device prim cache after the RHS vs the CPU `prim_checked` per cell.
    {
        let mut gcons = vec![0.0f64; n * NCOMP];
        let mut gprim = vec![0.0f64; n * NPRIM];
        dev.download(&mut gcons, &mut gprim);
        let mut cprim = vec![0.0f64; n * NPRIM];
        for j in 0..nt {
            for rz in 0..nrz {
                if act[rz] == 0 {
                    continue;
                }
                let cc = j * nrz + rz;
                let u: Cons = std::array::from_fn(|k| cons1[cc * NCOMP + k]);
                let w = eos.prim_checked(&u).expect("cpu prim");
                cprim[cc * NPRIM..(cc + 1) * NPRIM].copy_from_slice(&w);
            }
        }
        println!("PHASE 0 — prims (cold) at t = {t1:.4e}");
        let (w0, c0) = compare("prim", &cprim, &gprim, NPRIM);
        assert!(
            w0 < ECT_PRIM,
            "prims diverged beyond the projection-tolerance class: {w0:.3e}"
        );
        assert!(
            c0 < CANCEL_BOUND,
            "prim per-cell rel {c0:.3e} (prims carry no cancellation)"
        );
        // The worst RHS cell's neighbourhood prims (θ-momentum diagnosis).
        let (j, i_r, i_z) = (2usize, 0usize, 4usize);
        for dj in [nt - 1, 0, 1, nt / 2] {
            let jj = (j + dj) % nt;
            let cc = lay.dense(jj, i_r, i_z);
            println!(
                "    (j={jj}, i_r={i_r}, i_z={i_z}) cpu p {:.12e} e {:.12e} g1 {:.12e} b {:.3e} | gpu p {:.12e} e {:.12e} g1 {:.12e}",
                cprim[cc * NPRIM + 4],
                cprim[cc * NPRIM + 7],
                cprim[cc * NPRIM + 8],
                cprim[cc * NPRIM + 6],
                gprim[cc * NPRIM + 4],
                gprim[cc * NPRIM + 7],
                gprim[cc * NPRIM + 8]
            );
        }
    }
    dev.upload(&cons1, None);
    let (gpu_rate2, gl2) = dev.rhs(&sched, t1).expect("gpu rhs rerun");
    println!(
        "PHASE 1 — full RHS (blend EOS + general HLLC + inflow/outflow + igniter + combustion) at t = {t1:.4e}"
    );
    let (w1, c1) = compare("rhs", &cpu_rate, &gpu_rate, NCOMP);
    let mut led_worst = 0.0f64;
    for k in 0..NCOMP {
        // Net sums are scaled by their own GROSS (a cancelled net of 1e-13 over a
        // gross of 0.1 is rounding, exactly the COUP-2 §3.1.1 scale rule).
        for (a, b, scale) in [
            (cpu_led.port_net[k], gl.port_net[k], cpu_led.port_abs[k]),
            (cpu_led.port_abs[k], gl.port_abs[k], cpu_led.port_abs[k]),
            (cpu_led.src_net[k], gl.src_net[k], cpu_led.src_abs[k]),
            (cpu_led.src_abs[k], gl.src_abs[k], cpu_led.src_abs[k]),
        ] {
            let s = a.abs().max(b.abs()).max(scale);
            if s > 1e-300 {
                led_worst = led_worst.max((a - b).abs() / s);
            }
        }
        println!(
            "    k={k}: port_net {:+.6e} {:+.6e} | port_abs {:.6e} {:.6e} | src_net {:+.6e} {:+.6e} | src_abs {:.6e} {:.6e}",
            cpu_led.port_net[k],
            gl.port_net[k],
            cpu_led.port_abs[k],
            gl.port_abs[k],
            cpu_led.src_net[k],
            gl.src_net[k],
            cpu_led.src_abs[k],
            gl.src_abs[k]
        );
    }
    println!("  ledger (port/src net+gross × 7, gross-scaled): worst rel {led_worst:.3e}");
    let bit1 = gpu_rate == gpu_rate2 && gl.port_net == gl2.port_net && gl.src_abs == gl2.src_abs;
    println!(
        "  GPU same-build rerun: {}",
        if bit1 {
            "BIT-IDENTICAL"
        } else {
            "*** DIFFERS ***"
        }
    );
    assert!(
        w1 < ECT,
        "RHS diverged beyond ECT (component-scaled): {w1:.3e}"
    );
    let _ = c1; // per-cell rel is diagnostic only: ECT × comp_max / |value| at cancellation cells
    assert!(led_worst < 1e-8, "ledger diverged: {led_worst:.3e}");
    assert!(bit1, "GPU RHS rerun not deterministic");
    println!("  PHASE 1 PASS");

    // ================================================================
    // PHASE 1c — the S14 flux-buffer sweep path vs the fused kernels: the
    // same formulas on the same device ⇒ BYTE-IDENTICAL rate + ledger.
    // ================================================================
    dev.set_fluxbuf(true);
    dev.upload(&cons1, None);
    let (gpu_rate_fb, gl_fb) = dev.rhs(&sched, t1).expect("gpu rhs (flux buffer)");
    // The two paths are the same formulas, but nvcc's FMA contraction is
    // decided per expression context (the per-component scalar PPM vs the
    // array form) — measured NOT bit-identical (4e-11 per-cell), so the
    // contract is: bit-identical WITHIN a path (same-build reruns), ECT-class
    // ACROSS paths (a build-variant, like CPU↔GPU). Both compared to the CPU.
    let fb_same = gpu_rate_fb == gpu_rate;
    println!(
        "PHASE 1c — S14 flux-buffer sweeps vs fused kernels (same device): {}",
        if fb_same {
            "BIT-IDENTICAL"
        } else {
            "not bit-identical (FMA contraction differs per path — ECT-class)"
        }
    );
    let (w1c, _) = compare("flux-buffer vs fused", &gpu_rate, &gpu_rate_fb, NCOMP);
    let (w1d, _) = compare("flux-buffer vs CPU", &cpu_rate, &gpu_rate_fb, NCOMP);
    let mut fb_led = 0.0f64;
    for k in 0..NCOMP {
        for (a, b, sc) in [
            (gl.port_net[k], gl_fb.port_net[k], gl.port_abs[k]),
            (gl.src_net[k], gl_fb.src_net[k], gl.src_abs[k]),
        ] {
            let s = a.abs().max(b.abs()).max(sc);
            if s > 1e-300 {
                fb_led = fb_led.max((a - b).abs() / s);
            }
        }
    }
    println!("  flux-buffer ledger vs fused (gross-scaled): worst rel {fb_led:.3e}");
    assert!(
        w1c < ECT,
        "flux-buffer path diverged from the fused path beyond ECT: {w1c:.3e}"
    );
    assert!(
        w1d < ECT,
        "flux-buffer path diverged from the CPU beyond ECT: {w1d:.3e}"
    );
    dev.set_fluxbuf(false);
    println!("  PHASE 1c PASS");

    // ================================================================
    // PHASE 1b — stable_dt (θ-arc + front carrier) at the pre-marched state.
    // ================================================================
    let cpu_dt = {
        let flow = FlowClass {
            op: &op,
            fields: &spec.fields,
        };
        sdc.stable_dt(&spec.grid, &flow, sched.cfl)
            .expect("cpu stable_dt")
    };
    let gpu_dt = dev.stable_dt(t1, sched.cfl).expect("gpu stable_dt");
    let sdt_rel = (cpu_dt - gpu_dt).abs() / cpu_dt.abs().max(gpu_dt.abs());
    println!("PHASE 1b — stable_dt: cpu {cpu_dt:.9e} gpu {gpu_dt:.9e} rel {sdt_rel:.3e}");
    assert!(sdt_rel < ECT, "stable_dt diverged");
    println!("  PHASE 1b PASS");

    // ================================================================
    // PHASE 2 — M2 coupled steps: audit rows per step, final state, rerun +
    // checkpoint bit-identity.
    // ================================================================
    let gpu_march = |from: &[f64],
                     t_start: f64,
                     dts: &[f64]|
     -> (Vec<f64>, Vec<f64>, Vec<crucible_gpu::engine_host::StepOut>) {
        dev.upload(from, None);
        let mut tt = t_start;
        let mut outs = Vec::new();
        for &dt in dts {
            outs.push(
                dev.step(&sched, tt, dt)
                    .unwrap_or_else(|e| panic!("GPU step at t = {tt:.4e}: {e}")),
            );
            tt += dt;
        }
        let mut cons = vec![0.0f64; n * NCOMP];
        let mut hint = vec![0.0f64; n * NPRIM];
        dev.download(&mut cons, &mut hint);
        (cons, hint, outs)
    };
    let mut dts2 = Vec::new();
    let mut cpu_rows = Vec::new();
    let mut t2 = t1;
    for _ in 0..m2 {
        op.bcs.z_lo = FlowBc::MassFlowInflow {
            mdot_per_area: sched.mdot_per_area(t2),
            h_total,
            c_frac,
        };
        let flow = FlowClass {
            op: &op,
            fields: &spec.fields,
        };
        let dt = sdc
            .stable_dt(&spec.grid, &flow, sched.cfl)
            .expect("cpu stable_dt")
            .min(sched.t_final - t2);
        let rep = sdc
            .step(
                &mut spec.grid,
                Some(&flow),
                None,
                None,
                None,
                Some(&reaction),
                t2,
                dt,
            )
            .expect("cpu step");
        cpu_rows.push(rep.audit.clone());
        dts2.push(dt);
        t2 += dt;
    }
    let cpu_final = lay.read_dense(&spec.grid, &ids);
    let (gpu_final, _hint_final, outs) = gpu_march(&cons1, t1, &dts2);
    println!(
        "PHASE 2 — {m2} coupled SDC steps (flow + class-R, stagewise SRD, audited) from t = {t1:.4e}"
    );
    let mut worst_gap_tol = 0.0f64;
    let mut worst_row_rel = 0.0f64;
    for (s, o) in outs.iter().enumerate() {
        let rows = sdc_gpu
            .xcheck_audit_flow(
                sched.n_cells,
                dts2[s],
                o.before,
                o.after,
                &o.l0,
                &o.l_last,
                o.burn_applied,
                o.burn_gross,
            )
            .unwrap_or_else(|e| panic!("GPU-operand audit VIOLATION at step {s}: {e}"));
        for (rc, rg) in cpu_rows[s].iter().zip(&rows) {
            worst_gap_tol = worst_gap_tol.max((rg.delta - rg.applied).abs() / rg.tol);
            for (a, b) in [
                (rc.delta, rg.delta),
                (rc.applied, rg.applied),
                (rc.tol, rg.tol),
            ] {
                // In units of the row's own tolerance: a cancelled delta of
                // 1e-13 differing by 30 % is rounding, not drift.
                worst_row_rel = worst_row_rel.max((a - b).abs() / rc.tol);
            }
        }
    }
    println!(
        "  audit: every GPU step passes the CPU's own identity; worst |Δ−applied|/tol = {worst_gap_tol:.3e}; \
         worst row-value |Δ|/tol vs the CPU report {worst_row_rel:.3e}"
    );
    let (w2, c2) = compare("state", &cpu_final, &gpu_final, NCOMP);
    let (gpu_final2, _, _) = gpu_march(&cons1, t1, &dts2);
    let bit2 = gpu_final == gpu_final2;
    let h1 = m2 / 2;
    let (mid, mid_hint, _) = gpu_march(&cons1, t1, &dts2[..h1]);
    dev.upload(&mid, Some(&mid_hint));
    let mut tt = t1 + dts2[..h1].iter().sum::<f64>();
    for &dt in &dts2[h1..] {
        dev.step(&sched, tt, dt).expect("resumed step");
        tt += dt;
    }
    let mut split = vec![0.0f64; n * NCOMP];
    let mut split_hint = vec![0.0f64; n * NPRIM];
    dev.download(&mut split, &mut split_hint);
    let ck2 = split == gpu_final;
    println!(
        "  GPU resident-march rerun: {}",
        if bit2 {
            "BIT-IDENTICAL"
        } else {
            "*** DIFFERS ***"
        }
    );
    println!(
        "  checkpoint (cons + hints) split at step {h1}: {}",
        if ck2 {
            "BIT-IDENTICAL"
        } else {
            "*** DIFFERS ***"
        }
    );
    assert!(
        w2 < 1e-8,
        "marched state diverged beyond the FMA-order band: {w2:.3e}"
    );
    let _ = c2;
    assert!(bit2, "GPU engine march not deterministic");
    assert!(ck2, "checkpoint/restart not bit-faithful");
    println!("  PHASE 2 PASS");

    // ================================================================
    // PHASE 3 — throughput of the composed step (resident, own dt).
    // ================================================================
    let iters = 50usize;
    for fluxbuf in [false, true] {
        dev.set_fluxbuf(fluxbuf);
        dev.upload(&cons1, None);
        let mut tt = t1;
        let t0 = std::time::Instant::now();
        for _ in 0..iters {
            let dt = dev.stable_dt(tt, sched.cfl).expect("dt");
            dev.step(&sched, tt, dt).expect("step");
            tt += dt;
        }
        let el = t0.elapsed().as_secs_f64();
        println!(
            "PHASE 3 — {iters} full engine steps ({}) in {el:.3} s: {:.4} s/step, {:.3e} cell-steps/s",
            if fluxbuf {
                "S14 flux-buffer sweeps"
            } else {
                "fused sweep kernels"
            },
            el / iters as f64,
            iters as f64 * sched.n_cells as f64 / el
        );
    }
    println!("ALL PASS");
}
