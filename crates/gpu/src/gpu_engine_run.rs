//! S14 — the OVERNIGHT HARNESS: the resident full-physics engine step
//! (`cuda/residency_engine.cu`, validated by `residency_engine_xcheck`)
//! marched under the run's own schedule (`crucible_engine::run` replicated in
//! `engine_host::Schedule`), with:
//!   • the COUP-2 audit every step — the CPU's own `Sdc` identity on the
//!     device's reduced operands (a violation halts, never continues);
//!   • the COUP-4 verdict machinery at the run's `PROBE_EVERY` cadence on the
//!     downloaded state, through the engine's own probe functions
//!     (`plane_mdot`, `plane_thrust`, `injector_end_stagnation_p`,
//!     `reacting_measure`) — ignition floor, NEVER_IGNITED, FLAMEOUT, the dwell
//!     → WORKS, FAILED_TO_REACH at the horizon;
//!   • FND-6 §3.8 checkpoints to disk (the dense state + the warm-start prim
//!     cache + the clock + the verdict trackers + the pinned-input digests),
//!     written atomically; `--resume` verifies the digests and continues the
//!     SAME trajectory (the device round-trip is bit-faithful);
//!   • halt artifacts: the crash/final fields CSV (`run::fields_csv`) + a
//!     verdict file; a progress CSV at every probe.
//!
//! Usage: gpu_engine_run [config.toml] [--out DIR] [--probe-every N]
//!        [--checkpoint-every N] [--resume FILE] [--max-steps N] [--max-wall-s S]
//! Paths are repo-root relative (the harness chdirs to the repo root).
use crucible_engine::eos_sel::ChemEos;
use crucible_engine::run::{
    EPS_IGNITED_FRAC, NEVER_IGNITED_GRACE_S, PROBE_EVERY, exit_plane, fields_csv,
    injector_end_stagnation_p, plane_mdot, plane_thrust,
};
use crucible_gpu::engine_host::{DeviceEngine, Schedule, load_run};
use crucible_grid::FaceDir;
use crucible_solvers::euler::{
    Combustion, Cons, EPS_IGNITED, I_RB, IgnitionColumns, NCOMP, NPRIM, THETA_CELLS,
    reacting_measure,
};
use crucible_solvers::sdc::Sdc;
use std::io::Write;

const CKPT_MAGIC: &str = "CRUCIBLE-GPU-CKPT v1";

struct Args {
    config: String,
    out: Option<String>,
    probe_every: usize,
    checkpoint_every: usize,
    resume: Option<String>,
    max_steps: usize,
    max_wall_s: f64,
}
fn parse_args() -> Args {
    let mut a = Args {
        config: "configs/rl10_startup_3d.toml".into(),
        out: None,
        probe_every: PROBE_EVERY,
        checkpoint_every: 2000,
        resume: None,
        max_steps: usize::MAX,
        max_wall_s: f64::INFINITY,
    };
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let mut i = 0;
    while i < argv.len() {
        let next = |i: &mut usize| -> String {
            *i += 1;
            argv.get(*i)
                .cloned()
                .unwrap_or_else(|| panic!("missing value after {}", argv[*i - 1]))
        };
        match argv[i].as_str() {
            "--out" => a.out = Some(next(&mut i)),
            "--probe-every" => a.probe_every = next(&mut i).parse().expect("probe-every"),
            "--checkpoint-every" => {
                a.checkpoint_every = next(&mut i).parse().expect("checkpoint-every")
            }
            "--resume" => a.resume = Some(next(&mut i)),
            "--max-steps" => a.max_steps = next(&mut i).parse().expect("max-steps"),
            "--max-wall-s" => a.max_wall_s = next(&mut i).parse().expect("max-wall-s"),
            s if !s.starts_with("--") => a.config = s.to_string(),
            s => panic!("unknown flag {s}"),
        }
        i += 1;
    }
    a
}

/// The march clock + the COUP-4 trackers — everything a resume must restore.
#[derive(Debug, Clone)]
struct Clock {
    t: f64,
    steps: usize,
    dwell_start: Option<f64>,
    ignited: bool,
    peak_r: f64,
}

/// Checkpoint: text header + raw little-endian f64 blocks (cons, prim hints).
fn write_checkpoint(
    path: &str,
    digests: &str,
    clk: &Clock,
    primed: bool,
    cons: &[f64],
    hint: &[f64],
) -> Result<(), String> {
    let tmp = format!("{path}.tmp");
    let mut f =
        std::io::BufWriter::new(std::fs::File::create(&tmp).map_err(|e| format!("{tmp}: {e}"))?);
    let hdr = format!(
        "{CKPT_MAGIC}\ndigests={digests}\nt={:e}\nsteps={}\ndwell_start={}\nignited={}\npeak_r={:e}\nprimed={}\nn_cons={}\nn_hint={}\n---\n",
        clk.t,
        clk.steps,
        clk.dwell_start
            .map_or("none".to_string(), |v| format!("{v:e}")),
        clk.ignited,
        clk.peak_r,
        primed,
        cons.len(),
        hint.len()
    );
    f.write_all(hdr.as_bytes()).map_err(|e| e.to_string())?;
    for v in cons.iter().chain(hint) {
        f.write_all(&v.to_le_bytes()).map_err(|e| e.to_string())?;
    }
    f.flush().map_err(|e| e.to_string())?;
    drop(f);
    std::fs::rename(&tmp, path).map_err(|e| format!("rename {tmp}: {e}"))
}
fn read_checkpoint(
    path: &str,
    digests: &str,
    n_cons: usize,
    n_hint: usize,
) -> Result<(Clock, bool, Vec<f64>, Vec<f64>), String> {
    let bytes = std::fs::read(path).map_err(|e| format!("{path}: {e}"))?;
    let sep = b"\n---\n";
    let pos = bytes
        .windows(sep.len())
        .position(|w| w == sep)
        .ok_or("checkpoint: no header separator")?;
    let hdr = std::str::from_utf8(&bytes[..pos]).map_err(|e| e.to_string())?;
    let mut lines = hdr.lines();
    if lines.next() != Some(CKPT_MAGIC) {
        return Err("checkpoint: bad magic".into());
    }
    let mut kv = std::collections::BTreeMap::new();
    for l in lines {
        if let Some((k, v)) = l.split_once('=') {
            kv.insert(k.to_string(), v.to_string());
        }
    }
    if kv.get("digests").map(String::as_str) != Some(digests) {
        return Err(format!(
            "checkpoint REFUSED: pinned-input digests differ (checkpoint {:?} vs this build/config {digests:?}) — FND-6 §3.8",
            kv.get("digests")
        ));
    }
    let g = |k: &str| kv.get(k).ok_or_else(|| format!("checkpoint: missing {k}"));
    let clk = Clock {
        t: g("t")?.parse().map_err(|e| format!("t: {e}"))?,
        steps: g("steps")?.parse().map_err(|e| format!("steps: {e}"))?,
        dwell_start: match g("dwell_start")?.as_str() {
            "none" => None,
            s => Some(s.parse().map_err(|e| format!("dwell_start: {e}"))?),
        },
        ignited: g("ignited")? == "true",
        peak_r: g("peak_r")?.parse().map_err(|e| format!("peak_r: {e}"))?,
    };
    let primed = g("primed")? == "true";
    let nc: usize = g("n_cons")?.parse().map_err(|e| format!("n_cons: {e}"))?;
    let nh: usize = g("n_hint")?.parse().map_err(|e| format!("n_hint: {e}"))?;
    if nc != n_cons || nh != n_hint {
        return Err(format!(
            "checkpoint REFUSED: state size {nc}/{nh} vs the assembled world {n_cons}/{n_hint}"
        ));
    }
    let body = &bytes[pos + sep.len()..];
    if body.len() != 8 * (nc + nh) {
        return Err("checkpoint: truncated body".into());
    }
    let mut vals = Vec::with_capacity(nc + nh);
    for ch in body.chunks_exact(8) {
        vals.push(f64::from_le_bytes(ch.try_into().unwrap()));
    }
    let hint = vals.split_off(nc);
    Ok((clk, primed, vals, hint))
}

fn main() {
    let root = format!("{}/../..", env!("CARGO_MANIFEST_DIR"));
    std::env::set_current_dir(&root).expect("repo root");
    let args = parse_args();
    let started = std::time::Instant::now();

    // --- Load + assemble (the CLI's path) -------------------------------------
    let lr = load_run(&args.config).unwrap_or_else(|e| {
        eprintln!("error: {e}");
        std::process::exit(1)
    });
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
    let sdc = Sdc::with_audit(sched.audit);
    // Pinned-input digests the checkpoint is bound to (FND-6 §3.8): the config
    // text + the table pins (file, version, digest) + the world shape.
    let digests = {
        let cfg = crucible_config::sha256_hex(
            crucible_gpu::engine_host::read_rel(&args.config)
                .unwrap()
                .as_bytes(),
        );
        let mut d = format!("cfg:{}", &cfg[..16]);
        d.push_str(&format!(
            ";eq:{}",
            &spec.table_pin.3.chars().take(24).collect::<String>()
        ));
        d.push_str(&format!(
            ";unb:{}",
            &bl.unburnt_pin.3.chars().take(24).collect::<String>()
        ));
        d.push_str(&format!(
            ";ign:{}",
            &bl.ignition_pin.3.chars().take(24).collect::<String>()
        ));
        d.push_str(&format!(";world:{}x{}x{}", lay.n_r, lay.n_z, lay.nt));
        d
    };
    let out_dir = args
        .out
        .clone()
        .unwrap_or_else(|| format!("runs/{}-gpu", lr.name));
    std::fs::create_dir_all(&out_dir).expect("out dir");
    let ckpt_path = format!("{out_dir}/checkpoint.bin");
    let n_cons = lay.n * NCOMP;
    let n_hint = lay.n * NPRIM;

    // --- Initial state: the fill, or the resumed checkpoint -------------------
    let mut cons = lay.read_dense(&spec.grid, &ids);
    let mut hint = vec![0.0f64; n_hint];
    let mut clk = Clock {
        t: 0.0,
        steps: 0,
        dwell_start: None,
        ignited: false,
        peak_r: f64::NAN,
    };
    let mut primed = false;
    if let Some(r) = &args.resume {
        let (c, p, cv, hv) = read_checkpoint(r, &digests, n_cons, n_hint).unwrap_or_else(|e| {
            eprintln!("error: {e}");
            std::process::exit(1)
        });
        clk = c;
        primed = p;
        cons = cv;
        hint = hv;
        println!(
            "RESUMED from {r}: step {} t = {:.6e} s (primed {primed})",
            clk.steps, clk.t
        );
    }
    dev.upload(&cons, if primed { Some(&hint) } else { None });
    println!(
        "CRUCIBLE GPU run `{}` on `{}`: {}×{}×{} = {} cells ({} active), {} SRD small cells; t_final {:.4e} s; \
         probe every {} steps, checkpoint every {} steps → {out_dir}",
        lr.name,
        args.config,
        lay.n_r,
        lay.n_z,
        lay.nt,
        lay.n,
        sched.n_cells,
        dev.n_small(),
        sched.t_final,
        args.probe_every,
        args.checkpoint_every
    );
    println!("  digests: {digests}");
    let mut progress = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(format!("{out_dir}/progress.csv"))
        .expect("progress.csv");
    if args.resume.is_none() {
        writeln!(
            progress,
            "step,t_s,dt_s,resid,mdot_exit_kg_s,thrust_n,p_c_pa,reacting_r_kg_s,wall_s"
        )
        .unwrap();
    }

    let exit_i = exit_plane(&spec.grid);
    let write_state_to_grid = |spec: &mut crucible_engine::assembly::EngineSpec, cons: &[f64]| {
        lay.write_dense(&mut spec.grid, &ids, cons);
    };
    let mut rho_probe: Vec<f64> = cons.iter().step_by(NCOMP).copied().collect();
    let mut halt: Option<(String, Option<(usize, usize)>)> = None;
    let mut works: Option<String> = None;
    let mut last_dt = 0.0f64;

    // --- The march --------------------------------------------------------------
    while sched.t_final - clk.t > 1e-12 * sched.t_final {
        if clk.steps >= args.max_steps || started.elapsed().as_secs_f64() > args.max_wall_s {
            println!(
                "STOP: step/wall cap reached at step {} t = {:.6e} s — checkpointing",
                clk.steps, clk.t
            );
            break;
        }
        let dt = match dev.stable_dt(clk.t, sched.cfl) {
            Ok(d) => d.min(sched.t_final - clk.t),
            Err(e) => {
                halt = Some((format!("stable_dt: {e}"), None));
                break;
            }
        };
        let out = match dev.step(&sched, clk.t, dt) {
            Ok(o) => o,
            Err(e) => {
                halt = Some((format!("step {}: {e}", clk.steps), None));
                break;
            }
        };
        // The COUP-2 audit — the CPU's own identity on the device operands.
        if let Err(e) = sdc.xcheck_audit_flow(
            sched.n_cells,
            dt,
            out.before,
            out.after,
            &out.l0,
            &out.l_last,
            out.burn_applied,
            out.burn_gross,
        ) {
            halt = Some((
                format!("COUP-2 audit VIOLATION at step {}: {e}", clk.steps),
                None,
            ));
            break;
        }
        clk.t += dt;
        clk.steps += 1;
        last_dt = dt;

        let probe_now = clk.steps.is_multiple_of(args.probe_every);
        let ckpt_now = clk.steps.is_multiple_of(args.checkpoint_every);
        if probe_now || ckpt_now {
            primed = dev.download(&mut cons, &mut hint);
        }
        if probe_now {
            write_state_to_grid(&mut spec, &cons);
            let now: Vec<f64> = cons.iter().step_by(NCOMP).copied().collect();
            let resid = rho_probe
                .iter()
                .zip(&now)
                .map(|(a, b)| ((b - a) / a.abs().max(1e-300)).abs())
                .fold(0.0f64, f64::max);
            rho_probe = now;
            let mdot_exit = plane_mdot(&spec.grid, &spec.fields, exit_i, FaceDir::ZPlus);
            let thrust = plane_thrust(&spec.grid, &spec.fields, &eos, exit_i, FaceDir::ZPlus)
                .unwrap_or(f64::NAN);
            let p_c =
                injector_end_stagnation_p(&spec.grid, &spec.fields, &eos, 0).unwrap_or(f64::NAN);
            let r_meas = match reacting_measure(
                &spec.grid,
                &ids,
                blend,
                &comb.ignition,
                comb.wrinkling,
                comb.theta,
            ) {
                Ok(r) => r,
                Err(e) => {
                    halt = Some((format!("reacting measure: {e}"), None));
                    break;
                }
            };
            clk.peak_r = if clk.peak_r.is_nan() {
                r_meas
            } else {
                clk.peak_r.max(r_meas)
            };
            let delivered = sched.mdot_delivered(clk.t);
            let floor = EPS_IGNITED.max(EPS_IGNITED_FRAC * delivered);
            if r_meas >= floor {
                clk.ignited = true;
            } else if clk.ignited {
                halt = Some((
                    format!(
                        "FLAMEOUT: the reacting measure collapsed to {r_meas:.3e} kg/s below its floor {floor:.3e} \
                         after establishment (peak {:.3e}) — quench/flammability physics (SOLV-4 §3.6)",
                        clk.peak_r
                    ),
                    max_flame_cell(&spec.grid, &ids),
                ));
                break;
            }
            let w_end = sched.igniter.t_off;
            if !clk.ignited && clk.t > w_end + NEVER_IGNITED_GRACE_S {
                halt = Some((
                    format!(
                        "NEVER_IGNITED: the igniter schedule is exhausted (window ended at {w_end:.3e} s + the \
                         establishment grace) and the reacting measure never exceeded its floor (peak {:.3e} kg/s vs \
                         floor {floor:.3e}) — no self-sustaining burn front exists (SOLV-4 §3.6 / COUP-4 §3.2)",
                        clk.peak_r
                    ),
                    None,
                ));
                break;
            }
            let wall = started.elapsed().as_secs_f64();
            println!(
                "  step {:>8}  t = {:>9.4} ms / {:.4} ms  dt {:.3e}  resid {:.2e}  mdot_exit {:>7.3} kg/s  F {:>9.1} N  \
                 p_c {:>9.4} MPa  R {:.2e}  [{:.0} s wall, {:.4} s/step]",
                clk.steps,
                clk.t * 1e3,
                sched.t_final * 1e3,
                dt,
                resid,
                mdot_exit,
                thrust,
                p_c / 1e6,
                r_meas,
                wall,
                wall / clk.steps as f64
            );
            writeln!(
                progress,
                "{},{:e},{:e},{:e},{:e},{:e},{:e},{:e},{:.1}",
                clk.steps, clk.t, dt, resid, mdot_exit, thrust, p_c, r_meas, wall
            )
            .unwrap();
            let _ = std::io::stdout().flush();
            // COUP-4 dwell tracker → WORKS.
            if let (Some(v), Some(t_dwell)) = (&spec.verdict, sched.t_dwell_s) {
                let mut all_in = true;
                for (cmd, meas) in [(v.commanded_p_c_pa, p_c), (v.commanded_thrust_n, thrust)] {
                    if let Some(cmd) = cmd {
                        let rel = ((meas - cmd) / cmd).abs();
                        // NaN-safe: a NaN readout is "not holding" (run.rs).
                        if rel.is_nan() || rel > v.eps_works {
                            all_in = false;
                        }
                    }
                }
                if all_in {
                    let since = *clk.dwell_start.get_or_insert(clk.t);
                    if clk.t - since >= t_dwell {
                        works = Some(format!(
                            "WORKS: dwell held — every commanded quantity inside eps_works = {} for {:.4e} s (t = {:.6e} s)",
                            v.eps_works, t_dwell, clk.t
                        ));
                        break;
                    }
                } else {
                    clk.dwell_start = None;
                }
            }
        }
        if ckpt_now {
            write_checkpoint(&ckpt_path, &digests, &clk, primed, &cons, &hint)
                .unwrap_or_else(|e| eprintln!("checkpoint: {e}"));
            println!(
                "  checkpoint written at step {} (t = {:.6e} s)",
                clk.steps, clk.t
            );
        }
    }

    // --- Close-out: final state, checkpoint, artifacts --------------------------
    primed = dev.download(&mut cons, &mut hint);
    write_checkpoint(&ckpt_path, &digests, &clk, primed, &cons, &hint)
        .unwrap_or_else(|e| eprintln!("checkpoint: {e}"));
    write_state_to_grid(&mut spec, &cons);
    let csv = fields_csv(&spec.grid, &spec.fields, &eos, spec.t_solid)
        .unwrap_or_else(|e| format!("# fields csv failed: {e}\n"));
    let horizon_done = sched.t_final - clk.t <= 1e-12 * sched.t_final;
    let verdict = if let Some((msg, loc)) = &halt {
        format!(
            "DOESN'T WORK — {msg}{}",
            loc.map_or(String::new(), |(r, z)| format!(" at cell ({r}, {z})"))
        )
    } else if let Some(w) = &works {
        w.clone()
    } else if horizon_done {
        let mut offenders = Vec::new();
        if let Some(v) = &spec.verdict {
            let thrust = plane_thrust(&spec.grid, &spec.fields, &eos, exit_i, FaceDir::ZPlus)
                .unwrap_or(f64::NAN);
            let p_c =
                injector_end_stagnation_p(&spec.grid, &spec.fields, &eos, 0).unwrap_or(f64::NAN);
            for (name, cmd, meas) in [
                ("p_c_pa", v.commanded_p_c_pa, p_c),
                ("thrust_n", v.commanded_thrust_n, thrust),
            ] {
                if let Some(cmd) = cmd {
                    let rel = (meas - cmd) / cmd;
                    if rel.is_nan() || rel.abs() > v.eps_works {
                        offenders.push(format!(
                            "{name} {meas:.6e} vs commanded {cmd:.6e} ({:+.2}%)",
                            rel * 100.0
                        ));
                    }
                }
            }
        }
        format!(
            "DOESN'T WORK — FAILED_TO_REACH: the Stage-1 horizon ({:.4e} s) expired without a completed dwell — {} (COUP-4 §3.2)",
            sched.t_final,
            if offenders.is_empty() {
                "inside tolerance at the horizon but the dwell never ran its full span".to_string()
            } else {
                offenders.join("; ")
            }
        )
    } else {
        format!(
            "PAUSED at step {} t = {:.6e} s (step/wall cap) — resume with --resume {ckpt_path}",
            clk.steps, clk.t
        )
    };
    let name = if halt.is_some() {
        "crash_fields.csv"
    } else {
        "fields.csv"
    };
    std::fs::write(format!("{out_dir}/{name}"), csv).expect("fields csv");
    std::fs::write(
        format!("{out_dir}/verdict.txt"),
        format!(
            "{verdict}\nsteps={} t={:.6e} last_dt={:.3e} wall_s={:.1}\n",
            clk.steps,
            clk.t,
            last_dt,
            started.elapsed().as_secs_f64()
        ),
    )
    .expect("verdict");
    println!("{verdict}");
    println!(
        "  {} steps, t = {:.6e} s, {:.1} s wall ({:.4} s/step); artifacts in {out_dir}",
        clk.steps,
        clk.t,
        started.elapsed().as_secs_f64(),
        started.elapsed().as_secs_f64() / clk.steps.max(1) as f64
    );
    if halt.is_some() {
        std::process::exit(2);
    }
}

/// The strongest front remnant (run.rs `max_flame_cell`): the cell with the
/// largest b(1−b).
fn max_flame_cell(
    g: &crucible_grid::Grid,
    ids: &[crucible_grid::FieldId; NCOMP],
) -> Option<(usize, usize)> {
    let mut best: Option<((usize, usize), f64)> = None;
    g.for_each_active_cell(|c| {
        let b = g.brick(c.bi);
        let bb = b.field(ids[I_RB])[c.idx] / b.field(ids[0])[c.idx];
        let s = bb * (1.0 - bb);
        if best.is_none_or(|(_, v)| s > v) {
            best = Some(((c.i_r, c.i_z), s));
        }
    });
    best.map(|(loc, _)| loc)
}
