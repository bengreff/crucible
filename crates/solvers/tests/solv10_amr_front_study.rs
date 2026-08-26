//! Plan S10 — **the AMR gate**: the 2-D smeared-vs-sharp front study, a
//! MEASURED go/no-go on dynamic front-tracking refinement (plan §1 ruling #7,
//! §5, §7 "dynamic AMR complexity spiral"). Not an assumption — a number.
//!
//! The closure-set pushed front (SOLV-4 §3.6) travels at `S_T` with a width of
//! a fixed `Θ` cells at **every** resolution. So on a coarse grid the front is
//! physically WIDE (Θ·h_c — "smeared"); on a fine grid it is physically SHARP
//! (Θ·h_f). Dynamic front-tracking would keep the front band at the fine
//! resolution while it moves — the question this study answers is: **what does
//! that buy the ANSWER?** The startup mission's verdict (COUP-4) is timeline-
//! driven — *when* does the chamber fill / the pressure rise / the flame reach a
//! plane. So the operative measurement is not the front's sharpness but whether
//! sharpening it moves the **timeline**.
//!
//! Method (real SDC combustion marches, audit armed — the `flame_1d` fixture's
//! machinery, extended to track the front over time):
//!   * march the identical physical flame tube at coarse `h_c` and fine
//!     `h_f = h_c/2`, same IC, same settle time;
//!   * sample the front position `z*(t)` (where `b` crosses ½) at a schedule of
//!     times → the **timeline delta** (does the sharp front arrive earlier?);
//!   * measure the front thickness in physical units on each → the **sharpness**
//!     refinement actually buys;
//!   * count `steps × cells` for the equal physical march → the **cost** a
//!     front-tracking scheme would pay to buy that sharpness.
//!
//! The uniform-coarse vs uniform-fine pair **brackets** static front refinement
//! (a static fine zone at the front is the fine resolution locally) and bounds
//! dynamic front-tracking from above (it pays fine cost only in the band, and
//! can buy at most what uniform-fine buys). So: if the fine timeline ≈ the
//! coarse timeline, dynamic front-tracking cannot move the verdict — NO-GO, and
//! the plan's declared fallback (static refinement + closure-set speed = blurry
//! front, correct timeline; §7) stands measured, not assumed.

use crucible_grid::{Grid, GridSpec};
use crucible_solvers::euler::{
    BurnBlendEos, Combustion, Cons, Euler, EulerFields, FlowBc, FlowBcs, I_RB, I_RHO, NCOMP,
    THETA_CELLS, TableEos,
};
use crucible_solvers::sdc::{FlowClass, ReactionClass, Sdc};
use crucible_tables::{Pin, Table};

const P0: f64 = 1.0e5;
const H0: f64 = 4.364e4;
const Z0: f64 = 0.167;
const ZERO_SRC: fn(f64, f64, f64, f64) -> Cons = |_, _, _, _| [0.0; NCOMP];

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
fn unburnt() -> Table {
    open(
        "lox_lh2_unburnt_v0.3.0.h5",
        "/chem/lox_lh2/unburnt",
        include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../tables/chem/lox_lh2_unburnt_v0.3.0.pins.toml"
        )),
    )
}
fn burnt() -> Table {
    open(
        "lox_lh2_v0.4.0.h5",
        "/chem/lox_lh2/equilibrium",
        include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../tables/chem/lox_lh2_v0.4.0.pins.toml"
        )),
    )
}
fn ignition() -> Table {
    open(
        "lox_lh2_ignition_v0.3.0.h5",
        "/chem/lox_lh2/ignition",
        include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../tables/chem/lox_lh2_ignition_v0.3.0.pins.toml"
        )),
    )
}

/// A closed off-axis tube (flat metric, r_min = 1.0), uniform in r so the flame
/// is purely axial — `n_z` cells of size `h` (square dr = dz = h). Same fixture
/// as the `flame_1d` gate.
fn tube(n_z: usize, h: f64) -> (Grid, EulerFields) {
    let spec = GridSpec {
        r_min: 1.0,
        dr: h,
        n_r: 4,
        z_min: 0.0,
        dz: h,
        n_z,
        n_theta_max: 1,
        axisymmetry_assertion: true,
    };
    let g = Grid::build(spec, crucible_solvers::euler::EULER_FIELDS).expect("grid");
    let f = EulerFields::resolve(&g).expect("fields");
    (g, f)
}

/// Fill the exact Nagumo travelling-wave front `b(z) = 1/(1+exp((z−z₀)/w))`,
/// burnt below, unburnt above, at rest — no relaxation transient.
fn fill_flame(g: &mut Grid, f: &EulerFields, blend: &BurnBlendEos, z_flame: f64, w: f64) {
    let ids = f.ids();
    for (k, &id) in ids.iter().enumerate() {
        g.fill_field(id, |_r, _th, z| {
            let b = 1.0 / (1.0 + ((z - z_flame) / w).exp());
            blend
                .cons_from_phzb(P0, H0, Z0, b, [0.0, 0.0, 0.0])
                .expect("valid flame IC")[k]
        });
    }
}

/// The axial `b(z)` profile (burnt-mass fraction) along the tube centre column.
fn b_profile(g: &Grid, f: &EulerFields) -> Vec<f64> {
    let ids = f.ids();
    (0..g.spec().n_z)
        .map(|i_z| {
            let bi = g.brick_index(0, i_z).expect("active");
            let br = g.brick(bi);
            let cell = br.cell_index(0, Grid::local_rz(0, i_z));
            br.field(ids[I_RB])[cell] / br.field(ids[I_RHO])[cell]
        })
        .collect()
}

/// The front position `z*` = the (interpolated) `z` where `b` crosses ½,
/// scanning from the burnt end. `NaN` if no crossing (front left the domain).
fn front_z(b: &[f64], h: f64) -> f64 {
    for i in 0..b.len() - 1 {
        if b[i] >= 0.5 && b[i + 1] < 0.5 {
            let frac = (b[i] - 0.5) / (b[i] - b[i + 1]);
            return (i as f64 + 0.5 + frac) * h;
        }
    }
    f64::NAN
}

/// The front thickness in physical metres: `∫ b(1−b) dz` (a resolution-
/// independent-in-Θ-cells witness scaled back to metres — the "sharpness").
fn front_thickness_m(b: &[f64], h: f64) -> f64 {
    b.iter().map(|&x| x * (1.0 - x) * h).sum()
}

/// One tracked march: returns (front-position track at the sample times, final
/// physical thickness, total step count). The front starts at `z_flame` and
/// travels up-tube; sampled at each `t` in `sample_times`.
fn track_front(n_z: usize, h: f64, z_flame: f64, sample_times: &[f64]) -> (Vec<f64>, f64, usize) {
    let ut = unburnt();
    let bt = burnt();
    let it = ignition();
    let blend = BurnBlendEos::new(
        TableEos::bind(&ut).expect("unburnt binds"),
        TableEos::bind(&bt).expect("burnt binds"),
    );
    let ign = crucible_solvers::euler::IgnitionColumns::bind(&it).expect("ignition binds");
    let comb = Combustion {
        blend: &blend,
        ignition: ign,
        wrinkling: 1.0,
        theta: THETA_CELLS,
    };
    let (mut g, f) = tube(n_z, h);
    fill_flame(&mut g, &f, &blend, z_flame, THETA_CELLS * h);
    let op = Euler {
        eos: blend.clone(),
        source: &ZERO_SRC,
        bcs: FlowBcs {
            r_inner: FlowBc::Reflecting,
            r_outer: FlowBc::Reflecting,
            z_lo: FlowBc::Reflecting,
            z_hi: FlowBc::Reflecting,
        },
        wall_normal: None,
        slip_wall_z_faces: true,
        combustion: Some(&comb),
    };
    let flow = FlowClass {
        op: &op,
        fields: &f,
    };
    let reaction = ReactionClass { op: &comb };
    let mut sdc = Sdc::new();
    let mut t = 0.0f64;
    let mut steps = 0usize;
    let mut track = Vec::with_capacity(sample_times.len());
    let mut si = 0usize;
    let t_end = *sample_times.last().unwrap();
    while t < t_end {
        // Stop exactly on each sample time so the tracks are compared at the
        // identical physical instants on both grids.
        let next_sample = sample_times[si];
        let dt = sdc
            .stable_dt(&g, &flow, 0.4)
            .expect("dt")
            .min(next_sample - t);
        sdc.step(
            &mut g,
            Some(&flow),
            None,
            None,
            None,
            Some(&reaction),
            t,
            dt,
        )
        .expect("combustion step (audit armed)");
        t += dt;
        steps += 1;
        assert!(steps < 2_000_000, "runaway flame march");
        if t + 1e-15 >= next_sample {
            track.push(front_z(&b_profile(&g, &f), h));
            si += 1;
            if si == sample_times.len() {
                break;
            }
        }
    }
    let thick = front_thickness_m(&b_profile(&g, &f), h);
    (track, thick, steps)
}

/// THE AMR gate (plan S10): the smeared-vs-sharp front study. Marches the
/// identical physical flame the same times at coarse `h_c` and fine `h_c/2`,
/// and measures the timeline delta, the sharpness gained, and the cost paid.
#[test]
fn amr_gate_smeared_vs_sharp_front_timeline() {
    // A tube long enough for the front to travel many coarse cells. The front
    // starts a quarter up so it has clean run ahead of it and never reaches the
    // closed far end within the window (which would compress the reactant and
    // confound the arrival measurement).
    // Kept deliberately small so it fits the per-commit battery (~1 min, the
    // heavy-item budget) while the front still travels several coarse cells:
    // the closed-tube front runs at ~19 m/s (consumption speed + expansion
    // push), so ~5×10⁻⁵ s carries it ~5 coarse cells, clear of the far wall.
    let h_c = 2.0e-4;
    let n_c = 32usize; // 6.4 mm tube
    let z_flame = 0.25 * (n_c as f64) * h_c;
    let samples = [2.0e-5, 3.5e-5, 5.0e-5];

    let (track_c, thick_c, steps_c) = track_front(n_c, h_c, z_flame, &samples);
    // Fine: half the cell size, twice the cells, identical physical domain.
    let (track_f, thick_f, steps_f) = track_front(2 * n_c, 0.5 * h_c, z_flame, &samples);

    println!("=== AMR gate: smeared-vs-sharp front study (plan S10) ===");
    println!("coarse: h = {h_c:.2e} m, {n_c} cells, {steps_c} steps");
    println!(
        "fine:   h = {:.2e} m, {} cells, {steps_f} steps",
        0.5 * h_c,
        2 * n_c
    );
    println!("  t [s]      z*_coarse   z*_fine     Δz [m]      Δz / h_c");
    let mut max_dz_cells = 0.0f64;
    for (i, &t) in samples.iter().enumerate() {
        let (zc, zf) = (track_c[i], track_f[i]);
        let dz = (zc - zf).abs();
        let dz_cells = dz / h_c;
        max_dz_cells = max_dz_cells.max(dz_cells);
        println!("  {t:.2e}   {zc:.6}   {zf:.6}   {dz:.3e}   {dz_cells:.3}");
    }
    println!(
        "front thickness: coarse {:.3e} m, fine {:.3e} m (ratio {:.2})",
        thick_c,
        thick_f,
        thick_c / thick_f
    );
    println!(
        "cost (steps × cells): coarse {}, fine {} ({:.1}× more for the sharp front)",
        steps_c * n_c,
        steps_f * 2 * n_c,
        (steps_f * 2 * n_c) as f64 / (steps_c * n_c) as f64
    );

    // Every sample must have a live front on both grids (the fixture is sized
    // so the front never leaves the domain).
    for (i, &t) in samples.iter().enumerate() {
        assert!(
            track_c[i].is_finite() && track_f[i].is_finite(),
            "front left the domain at t = {t:.2e} (coarse {}, fine {})",
            track_c[i],
            track_f[i]
        );
    }

    // ---- THE MEASURED VERDICT ------------------------------------------------
    // (1) TIMELINE: the front position agrees between the smeared (coarse) and
    //     sharp (fine) grids to within ~1 coarse cell over the whole march. The
    //     timeline — when the front reaches any plane — is closure-set, not
    //     grid-set. This is the go/no-go number: refinement does NOT move the
    //     arrival timeline that drives the COUP-4 verdict.
    assert!(
        max_dz_cells < 1.5,
        "front timeline is grid-DEPENDENT: max |Δz| = {max_dz_cells:.3} coarse cells \
         (a GO signal for dynamic front-tracking — escalate to Ben, §5)"
    );
    // (2) SHARPNESS: the fine front is ~2× thinner in physical units (fixed Θ
    //     cells × half h). This is ALL that refinement buys — localization, not
    //     timeline. Bracketed [1.6, 2.4] around the ideal 2.0 (discrete-front
    //     ∫b(1−b) sampling + the finite settle).
    let ratio = thick_c / thick_f;
    assert!(
        (1.6..=2.4).contains(&ratio),
        "front-thickness ratio {ratio:.2} not ≈ 2 — the front is not a fixed cell count wide"
    );
    // The verdict, in the log: sharpness gained, timeline unchanged ⇒ NO-GO on
    // dynamic front-tracking for the startup-verification mission; the plan's
    // declared static-refinement fallback stands measured. (Full numbers +
    // disposition in SESSION_LOG / PLAN §8.)
}
