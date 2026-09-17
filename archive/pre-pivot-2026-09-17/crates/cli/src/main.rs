//! CRUCIBLE CLI — the VISION_SCOPE §4.3 workflow contract: *define geometry
//! → pick mechanisms → set an operating profile → run → receive
//! distributions*. `crucible run <config.toml>` executes one config-driven
//! engine assembly and prints the SOLV-7 performance readout; the fields
//! CSV + run manifest land in `runs/<name>/`.
//!
//! Paths inside a config (contour CSV, table artifacts, pin sidecars) are
//! resolved against the process working directory — run from the repository
//! root. Scalar report only this wave; the p-box wrapper (COUP-5 brackets)
//! and the FND-6 results bundle arrive with the cycle wave.

use std::io::Write as _;

fn main() {
    std::process::exit(real_main());
}

fn real_main() -> i32 {
    let args: Vec<String> = std::env::args().collect();
    match args.get(1).map(String::as_str) {
        Some("run") if args.len() == 3 => run_config(&args[2]),
        _ => {
            eprintln!(
                "CRUCIBLE {} — usage: crucible run <config.toml>\nConstants: {}",
                env!("CARGO_PKG_VERSION"),
                crucible_constants::SOURCE,
            );
            2
        }
    }
}

fn read_rel(path: &str) -> Result<String, String> {
    std::fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))
}

fn run_config(path: &str) -> i32 {
    let author = match read_rel(path) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("error: {e}");
            return 1;
        }
    };
    let registry = crucible_engine::registry();
    let loaded = match crucible_config::load_str_with_sidecars(&author, &registry, &read_rel) {
        Ok(l) => l,
        Err(diags) => {
            eprintln!("{diags}");
            return 1;
        }
    };
    let name = if loaded.resolved.meta.name.is_empty() {
        "run".to_string()
    } else {
        loaded.resolved.meta.name.clone()
    };

    let mut spec = match crucible_engine::assembly::assemble(&loaded, &read_rel) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("assembly error: {e}");
            return 1;
        }
    };
    let ext = loaded
        .resolved
        .geometry
        .as_ref()
        .and_then(|g| g.extents.as_ref())
        .expect("assembled config has extents");
    println!(
        "CRUCIBLE run `{name}`: grid {} x {} (dr = {:.4} mm), {} steady-march flow-throughs, CFL {}",
        ext.n_r,
        ext.n_z,
        ext.dr * 1e3,
        spec.flowthroughs,
        spec.cfl,
    );

    let table = match crucible_engine::run::open_pinned_table(&spec) {
        Ok(t) => t,
        Err(e) => {
            eprintln!("table error: {e}");
            return 1;
        }
    };
    let transport_table = match crucible_engine::run::open_transport_table(&spec) {
        Ok(t) => t,
        Err(e) => {
            eprintln!("table error: {e}");
            return 1;
        }
    };
    let blend_tables = match crucible_engine::run::open_blend_tables(&spec) {
        Ok(t) => t,
        Err(e) => {
            eprintln!("table error: {e}");
            return 1;
        }
    };
    let started = std::time::Instant::now();
    let mut on_progress = |p: &crucible_engine::run::Progress| {
        let r = if p.reacting_r.is_finite() {
            format!("  R {:.2e}", p.reacting_r)
        } else {
            String::new()
        };
        println!(
            "  step {:>7}  t = {:>8.4} ms / {:.4} ms  resid {:.2e}  mdot_exit {:>7.3} kg/s  F {:>9.1} N  p_c {:>9.4} MPa{r}",
            p.step,
            p.t * 1e3,
            p.t_final * 1e3,
            p.resid,
            p.mdot_exit,
            p.thrust_n,
            p.p_c_pa / 1e6,
        );
        let _ = std::io::stdout().flush();
    };
    let report = match crucible_engine::run::run(
        &mut spec,
        &table,
        transport_table.as_ref(),
        blend_tables.as_ref(),
        &mut on_progress,
    ) {
        Ok(r) => r,
        Err(halt) => {
            eprintln!("run halted: {halt}");
            if let Some(v) = &halt.verdict {
                eprintln!("\n== COUP-4 VERDICT ==\n  {v}");
                let dir = format!("runs/{name}");
                let path = format!("{dir}/verdict.txt");
                if std::fs::create_dir_all(&dir).is_ok()
                    && std::fs::write(&path, format!("{v}\n")).is_ok()
                {
                    eprintln!("  verdict artifact: {path}");
                }
            }
            if !halt.crash_csv.is_empty() {
                let dir = format!("runs/{name}");
                let path = format!("{dir}/crash_fields.csv");
                match std::fs::create_dir_all(&dir)
                    .map_err(|e| e.to_string())
                    .and_then(|()| {
                        std::fs::write(&path, &halt.crash_csv).map_err(|e| e.to_string())
                    }) {
                    Ok(()) => eprintln!("crash artifact: {path}"),
                    Err(e) => eprintln!("warning: crash artifact {path}: {e}"),
                }
            }
            return 1;
        }
    };
    let wall = started.elapsed().as_secs_f64();

    println!("\n== SOLV-7 performance readout (emergent, never imposed) ==");
    println!("  thrust           {:>12.1} N", report.thrust_n);
    println!("  Isp              {:>12.2} s", report.isp_s);
    println!("  v_e              {:>12.1} m/s", report.v_e_m_per_s);
    println!("  c*               {:>12.1} m/s", report.c_star_m_per_s);
    println!("  C_F              {:>12.4}", report.c_f);
    println!(
        "  p_c (N11 stag.)  {:>12.4} MPa  ({:.1} psia)",
        report.p_c_pa / 1e6,
        report.p_c_pa / 6894.757,
    );
    println!(
        "  mdot exit/inj    {:>12.3} / {:.3} kg/s  (inflow-plane measured {:.3})",
        report.mdot_exit_kg_per_s, report.mdot_injected_kg_per_s, report.mdot_inflow_plane_kg_per_s,
    );
    if ((report.mdot_inflow_plane_kg_per_s - report.mdot_injected_kg_per_s)
        / report.mdot_injected_kg_per_s)
        .abs()
        > 0.01
    {
        println!(
            "  WARNING: inflow plane delivers {:.3} kg/s vs the declared/solved {:.3} — \
             the injector face is still choked (sonic startup cap active at readout); \
             the run is NOT at its declared operating point",
            report.mdot_inflow_plane_kg_per_s, report.mdot_injected_kg_per_s,
        );
    }
    println!(
        "  jacket heat      {:>12.3} MW   liner T_max {:.1} K",
        report.jacket_watts / 1e6,
        report.liner_t_max_k,
    );
    if let Some(ex) = &report.expander {
        println!(
            "  expander (CLOSED mode): delivered mdot {:.4} kg/s, turbine {:.1} kW, \
             T_turb_in {:.1} K, resid {:.2e}",
            ex.mdot_kg_per_s,
            ex.turbine_power_w / 1e3,
            ex.t_turbine_in_k,
            ex.resid,
        );
    }
    println!(
        "  {} steps to t = {:.3} ms, steadiness resid {:.2e}, {} active gas cells, {:.1} s wall clock",
        report.steps,
        report.t_end * 1e3,
        report.steady_resid,
        report.active_gas_cells,
        wall,
    );
    if report.peak_reacting_measure.is_finite() {
        println!(
            "  peak reacting measure R {:.3e} kg/s (COUP-4 ignition witness)",
            report.peak_reacting_measure
        );
    }
    if let Some(v) = &report.verdict {
        println!("\n== COUP-4 VERDICT ==\n  {v}");
    }

    // Run artifacts: fields CSV (the viz feed) + the FND-4 manifest.
    let dir = format!("runs/{name}");
    if let Err(e) = std::fs::create_dir_all(&dir) {
        eprintln!("warning: cannot create {dir}: {e}");
        return 0;
    }
    let fields_path = format!("{dir}/fields.csv");
    let manifest_path = format!("{dir}/manifest.toml");
    if let Err(e) = std::fs::write(&fields_path, &report.fields_csv) {
        eprintln!("warning: {fields_path}: {e}");
    }
    if let Err(e) = std::fs::write(&manifest_path, loaded.manifest.to_toml()) {
        eprintln!("warning: {manifest_path}: {e}");
    }
    if let Some(v) = &report.verdict {
        let vpath = format!("{dir}/verdict.txt");
        if let Err(e) = std::fs::write(&vpath, format!("{v}\n")) {
            eprintln!("warning: {vpath}: {e}");
        }
    }
    println!("  artifacts: {fields_path}, {manifest_path}");
    0
}
