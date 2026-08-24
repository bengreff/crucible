//! SOLV-4 §3.6 (plan S6) — the **burn-progress ignition battery** on the
//! production artifacts: the blended EOS (unburnt ⊕ burnt) + the SOLV-4.4
//! reaction-diffusion source + the ignition closures.
//!
//! * `flame_1d` — **THE acceptance gate**: the front's propagation speed is
//!   closure-set (`S_L` from the surface), **grid-independent** — a coarse
//!   and a fine grid burn at the same rate. Refinement sharpens *where*, never
//!   *what* (SOLV-4 §3.6).
//! * `spark_box` — an igniter energy deposit lights a cold flammable mixture
//!   (auto-ignition seeds `c` past the bistable threshold `a`, the pushed
//!   front carries it) → the reacting measure `R` exceeds the ignition floor.
//! * `lean_no_light` — a deposit too weak to reach the induction threshold
//!   never lifts `c` out of the metastable well, so `R` never exceeds the
//!   floor (the `NEVER_IGNITED` outcome, COUP-4).
//! * `adiabatic_box_cannot_flame_out` — the S6-close finding (SOLV-4 0.4.3):
//!   quenching is a HEAT-LOSS phenomenon, so a lit **closed adiabatic** box
//!   must burn to completion — `FLAMEOUT` is unreachable without a loss
//!   channel. Pinned as a converse test; the real `quench_box` (conductive
//!   loss to cold walls, class-D coupling) rides S7 with its consumer, the
//!   COUP-4 `FLAMEOUT` verdict object.
//! * `ignition_delay` — a 0-D cell auto-ignites on the `τ_ign` timescale.
//!
//! All are marched through the REAL SDC step (the combustion source rides the
//! class-A rate inside `Euler::eval_rhs`), with the COUP-2 conservation audit
//! armed — closure is the standing condition.

use crucible_grid::{Grid, GridSpec};
use crucible_solvers::euler::{
    BurnBlendEos, Combustion, Cons, EPS_IGNITED, EosLaw, Euler, EulerFields, FlowBc, FlowBcs, I_RB,
    I_RHO, IgnitionColumns, NCOMP, THETA_CELLS, TableEos, consumption_rate, reacting_measure,
};
use crucible_solvers::sdc::{FlowClass, Sdc};
use crucible_tables::{Pin, Table};

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
        "lox_lh2_unburnt_v0.2.0.h5",
        "/chem/lox_lh2/unburnt",
        include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../tables/chem/lox_lh2_unburnt_v0.2.0.pins.toml"
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
        "lox_lh2_ignition_v0.2.0.h5",
        "/chem/lox_lh2/ignition",
        include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../tables/chem/lox_lh2_ignition_v0.2.0.pins.toml"
        )),
    )
}

// --- The design flame state (measured from the artifacts) --------------------
// Z = 0.167 (MR 5, in both surfaces' Z-envelopes); h picked so the unburnt
// branch reads ~312 K (flammable) and the burnt branch ~2971 K at 1 atm.
const P0: f64 = 1.0e5;
const H0: f64 = 4.364e4;
const Z0: f64 = 0.167;

const ZERO_SRC: fn(f64, f64, f64, f64) -> Cons = |_, _, _, _| [0.0; NCOMP];

/// A closed off-axis tube (flat metric, r_min = 1.0): uniform in r, so the
/// flame is purely axial and the radial combustion diffusion is identically
/// zero. `n_z` cells of size `h` (square: dr = dz = h).
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

/// March a closed flame tube of `n_z × h` for `t_settle`, then return the
/// **consumption speed** `S_c = ∫ω dV / (ρ_u·A)` — the frame- and
/// volume-independent burning speed. With the exact-Nagumo IC (no relaxation)
/// and the bistable pushed front, both coarse and fine evolve identically, so
/// `S_c` is closure-set (grid-independent). Closed ends keep `p` inside the
/// ignition-surface envelope (compression only, never rarefaction below its
/// floor). Returns `(S_c, ∫b(1−b)dz / (Θh))` — the second is the front
/// thickness in units of the declared width (a resolution witness).
fn consumption_speed(n_z: usize, h: f64, t_settle: f64) -> (f64, f64) {
    let ut = unburnt();
    let bt = burnt();
    let it = ignition();
    let blend = BurnBlendEos::new(
        TableEos::bind(&ut).expect("unburnt binds"),
        TableEos::bind(&bt).expect("burnt binds"),
    );
    let ign = IgnitionColumns::bind(&it).expect("ignition binds");
    let comb = Combustion {
        blend: &blend,
        ignition: ign,
        wrinkling: 1.0,
        theta: THETA_CELLS,
    };
    let (mut g, f) = tube(n_z, h);
    fill_flame(&mut g, &f, &blend, 0.5 * (n_z as f64) * h, THETA_CELLS * h);
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
    let mut sdc = Sdc::new();
    let mut t = 0.0;
    let mut steps = 0usize;
    while t < t_settle {
        let dt = sdc.stable_dt(&g, &flow, 0.4).expect("dt").min(t_settle - t);
        sdc.step(&mut g, Some(&flow), None, None, None, t, dt)
            .expect("combustion step (audit armed)");
        t += dt;
        steps += 1;
        assert!(steps < 500_000, "runaway flame march");
    }
    let ids = f.ids();
    let nt = g.brick(0).n_theta();
    let mut rho_u = 0.0f64;
    let mut thickness = 0.0;
    for i_z in 0..g.spec().n_z {
        let bi = g.brick_index(0, i_z).expect("active");
        let br = g.brick(bi);
        let cell = br.cell_index(0, Grid::local_rz(0, i_z));
        let rho = br.field(ids[I_RHO])[cell];
        let b = br.field(ids[I_RB])[cell] / rho;
        rho_u = rho_u.max(rho);
        thickness += b * (1.0 - b) * h;
    }
    let area: f64 = (0..g.spec().n_r).map(|i_r| g.face_area_z(i_r, nt)).sum();
    let cr = consumption_rate(&g, &ids, &blend, &comb.ignition, 1.0, THETA_CELLS)
        .expect("consumption rate");
    (cr / (rho_u * area), thickness / (THETA_CELLS * h))
}

/// Fill the tube with the **exact Nagumo travelling-wave** burn-progress front
/// `b(z) = 1/(1 + exp((z − z_flame)/w))` (width scale `w`), burnt (b→1) below
/// it and unburnt (b→0) above, both at (P0, H0, Z0), at rest. Initializing at
/// the settled profile means no pulled-/pushed-front relaxation transient.
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

#[test]
fn flame_1d_front_speed_is_grid_independent() {
    // THE acceptance gate (SOLV-4 §3.6): the front's consumption speed is
    // closure-set, not grid-set. Coarse and fine grids over the same physical
    // domain, marched the same settle time, must burn at the same speed —
    // refinement sharpens *where* (the front thins with h), never *what*.
    let t_settle = 1.2e-4;
    let (coarse, thk_c) = consumption_speed(48, 2.0e-4, t_settle);
    let (fine, thk_f) = consumption_speed(96, 1.0e-4, t_settle);
    println!(
        "flame_1d S_c: coarse = {coarse:.3} m/s (front {thk_c:.2} Θh), \
         fine = {fine:.3} m/s (front {thk_f:.2} Θh)"
    );
    assert!(
        coarse > 1.0 && fine > 1.0,
        "the flame must propagate (coarse {coarse}, fine {fine})"
    );
    // The front thickness scales with h (fixed Θ cells): its ∫b(1−b)/Θh
    // witness is the same on both grids.
    assert!(
        (thk_c - thk_f).abs() < 0.15,
        "front is not a fixed cell-count wide: {thk_c:.3} vs {thk_f:.3}"
    );
    // THE gate: the consumption speed is grid-independent.
    let rel = (coarse - fine).abs() / (0.5 * (coarse + fine));
    assert!(
        rel < 0.10,
        "front speed is grid-DEPENDENT: coarse {coarse:.3} vs fine {fine:.3} (rel {rel:.3})"
    );
}

// --- Igniter + 0-D auto-ignition scenarios -----------------------------------

/// Fill the whole tube uniformly at (p, h, Z, b), at rest.
fn fill_uniform(g: &mut Grid, f: &EulerFields, blend: &BurnBlendEos, h: f64, b: f64) {
    let ids = f.ids();
    for (k, &id) in ids.iter().enumerate() {
        g.fill_field(id, |_r, _th, _z| {
            blend
                .cons_from_phzb(P0, h, Z0, b, [0.0, 0.0, 0.0])
                .expect("valid uniform IC")[k]
        });
    }
}

/// Mean burnt mass fraction `∫ρb dV / ∫ρ dV`.
fn mean_burn(g: &Grid, f: &EulerFields) -> f64 {
    let ids = f.ids();
    let nt = g.brick(0).n_theta();
    let (mut mb, mut m) = (0.0, 0.0);
    for bi in 0..g.n_bricks() {
        let br = g.brick(bi);
        let mask = br.mask();
        let (rho, rhob) = (br.field(ids[I_RHO]), br.field(ids[I_RB]));
        for local in 0..64 {
            if mask & (1u64 << local) == 0 {
                continue;
            }
            let (i_r, _) = br.global_rz(local);
            let cell = br.cell_index(0, local);
            let v = g.cell_volume(i_r, nt);
            mb += rhob[cell] * v;
            m += rho[cell] * v;
        }
    }
    mb / m
}

#[test]
fn ignition_delay_reproduces_the_tau_surface() {
    // A 0-D-in-effect uniform hot reactant cell (T_u ≈ 1150 K, above the flame
    // crossover so S_L = 0 and only the auto-ignition term acts) auto-ignites
    // on the τ_ign timescale: dc/dt = (1−c)/τ ⇒ c crosses 0.5 at ≈ 0.69 τ.
    const H_HOT: f64 = 2.5e6;
    let ut = unburnt();
    let bt = burnt();
    let it = ignition();
    let blend = BurnBlendEos::new(TableEos::bind(&ut).unwrap(), TableEos::bind(&bt).unwrap());
    let ign = IgnitionColumns::bind(&it).unwrap();
    // The surface τ_ign at the uniform state (T_u from the unburnt branch).
    let w0 = blend
        .cons_from_phzb(P0, H_HOT, Z0, 0.0, [0.0, 0.0, 0.0])
        .and_then(|u| blend.prim_checked(&u))
        .unwrap();
    let t_u = blend.unburnt_temperature(&w0).unwrap();
    let comb = Combustion {
        blend: &blend,
        ignition: ign,
        wrinkling: 1.0,
        theta: THETA_CELLS,
    };
    let (mut g, f) = tube(8, 5.0e-4);
    fill_uniform(&mut g, &f, &blend, H_HOT, 0.0);
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
    let mut sdc = Sdc::new();
    let (mut t, mut t_ign) = (0.0f64, f64::NAN);
    while t < 5.0e-3 {
        let dt = sdc.stable_dt(&g, &flow, 0.4).unwrap();
        sdc.step(&mut g, Some(&flow), None, None, None, t, dt)
            .unwrap();
        t += dt;
        if t_ign.is_nan() && mean_burn(&g, &f) >= 0.5 {
            t_ign = t;
            break;
        }
    }
    // The surface τ_ign at this state (interrogated the same way the operator
    // does): the model must ignite on that timescale.
    let tau = comb.ignition.induction_time(P0, t_u, Z0).unwrap();
    println!(
        "ignition_delay: T_u={t_u:.0} K  τ_surface={tau:.3e} s  t_ign={t_ign:.3e} s  (t/τ={:.2})",
        t_ign / tau
    );
    assert!(t_ign.is_finite(), "the hot cell must auto-ignite");
    // dc/dt = (1−c)/τ gives t(c=0.5) ≈ 0.69 τ, but the induction is a thermal
    // RUNAWAY — as c rises the blend reads hotter, T_u climbs, τ_ign collapses
    // — so the cell ignites in a fraction of the initial τ (here ~0.13). The
    // model must ignite ON the τ_ign timescale (within ~an order), not 100×
    // off: that is the reproduction claim (the surface's own band is ~4×).
    let ratio = t_ign / tau;
    assert!(
        (0.05..=3.0).contains(&ratio),
        "ignition time {t_ign:.3e} is not on the τ_ign={tau:.3e} timescale (ratio {ratio:.2})"
    );
}

/// Run a cold flammable box (b = 0, T_u ≈ 312 K) with a **scheduled, ramped
/// igniter** energy deposit (COUP-7 §3.3): a kernel of half-width `half_cells`
/// centred in the tube, firing for `t_fire` at `power` [W/m³] (ramped over the
/// first fifth). March to `t_end`; return `(peak R, final mean burnt
/// fraction)` — the reacting measure's peak (COUP-4) and how much burned.
fn igniter_run(
    n_z: usize,
    h: f64,
    power: f64,
    half_cells: f64,
    t_fire: f64,
    t_end: f64,
) -> (f64, f64) {
    let ut = unburnt();
    let bt = burnt();
    let it = ignition();
    let blend = BurnBlendEos::new(TableEos::bind(&ut).unwrap(), TableEos::bind(&bt).unwrap());
    let ign = IgnitionColumns::bind(&it).unwrap();
    let comb = Combustion {
        blend: &blend,
        ignition: ign,
        wrinkling: 1.0,
        theta: THETA_CELLS,
    };
    let (mut g, f) = tube(n_z, h);
    fill_uniform(&mut g, &f, &blend, H0, 0.0);
    let z_c = 0.5 * (n_z as f64) * h;
    let half = half_cells * h;
    let t_ramp = 0.3 * t_fire;
    let igniter = move |_r: f64, _th: f64, z: f64, t: f64| -> Cons {
        let on = t < t_fire && (z - z_c).abs() < half;
        let q = if on {
            power * (t / t_ramp).min(1.0)
        } else {
            0.0
        };
        [0.0, 0.0, 0.0, 0.0, q, 0.0, 0.0]
    };
    let op = Euler {
        eos: blend.clone(),
        source: &igniter,
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
    let mut sdc = Sdc::new();
    let (mut t, mut peak_r) = (0.0f64, 0.0f64);
    let ids = f.ids();
    while t < t_end {
        let dt = sdc.stable_dt(&g, &flow, 0.4).unwrap().min(t_end - t);
        sdc.step(&mut g, Some(&flow), None, None, None, t, dt)
            .unwrap();
        t += dt;
        let r = reacting_measure(&g, &ids, &blend, &comb.ignition, 1.0, THETA_CELLS).unwrap();
        peak_r = peak_r.max(r);
    }
    (peak_r, mean_burn(&g, &f))
}

#[test]
fn spark_box_lights() {
    // A spark-class pulse lights the flammable box: R exceeds the ignition
    // floor and a substantial fraction burns. Sizing (`spark-igniter-class`,
    // META-3 §6.9): ~7.9 J delivered into the 12-cell kernel in 20 µs
    // (E/V ≈ 6.5e5 J/m³ ⇒ kernel h ≈ 2.5e6 J/kg, T_u ≈ 1050 K, τ_ign ~
    // 1e-4 s) — inside the aerospace exciter class (~0.1–20 J/discharge)
    // and inside every declared surface envelope. The pulse ENDS at about
    // the ignition time: a spark is a bounded electrical deposit (Ben
    // ruling, S6), not a sustained drive — continued firing into the
    // already-burnt kernel superheats it past the unburnt/ignition
    // envelopes (the metastable-reactant validity edge; the hot ASI-torch
    // regime is the declared S7 hardening).
    let (peak_r, burned) = igniter_run(40, 2.0e-4, 3.8e10, 6.0, 2.0e-5, 3.0e-4);
    println!("spark_box: peak_R={peak_r:.3e}  burned={burned:.3}");
    assert!(
        peak_r > EPS_IGNITED,
        "the spark must light the box (R {peak_r:.2e})"
    );
    assert!(
        burned > 0.3,
        "a lit box must burn substantially (burned {burned:.3})"
    );
}

#[test]
fn lean_no_light_never_ignites() {
    // Too weak a deposit never lifts the kernel to the induction threshold:
    // c stays below the metastable well `a`, R never exceeds the floor —
    // the NEVER_IGNITED outcome (COUP-4).
    let (peak_r, burned) = igniter_run(40, 2.0e-4, 6.0e9, 6.0, 5.0e-5, 4.0e-4);
    println!("lean_no_light: peak_R={peak_r:.3e}  burned={burned:.3}");
    assert!(
        peak_r < EPS_IGNITED,
        "a sub-threshold deposit must NOT ignite (R {peak_r:.2e})"
    );
    assert!(burned < 0.05, "nothing should burn (burned {burned:.3})");
}

#[test]
fn adiabatic_box_cannot_flame_out() {
    // The S6-close finding (SOLV-4 0.4.3 §6.5): quenching is a HEAT-LOSS
    // phenomenon, and this battery's box is closed and adiabatic — so a lit
    // kernel, however small (here ONE cell, ~0.7 J; half = 0.6 cells keeps the one-cell kernel off the exact cell-face ulp edge), must burn the box to
    // completion: the deposited energy has nowhere to go, the pressure rise
    // compression-heats the reactants, and the auto-ignition term finishes
    // whatever the front does not. FLAMEOUT (R collapsing after exceeding
    // the floor) is UNREACHABLE without a loss channel — this test pins that
    // converse so the physics cannot silently drift. The real `quench_box`
    // (conductive loss to cold isothermal walls through the class-D
    // gas-diffusion coupling; the fixture's 0.8 mm gap is already at the
    // H₂/O₂ quenching-distance scale) rides S7 with the blend↔diffusion
    // march its startup builds regardless — landing with its consumer, the
    // COUP-4 FLAMEOUT verdict object.
    let (peak_r, burned) = igniter_run(40, 2.0e-4, 3.8e10, 0.6, 2.0e-5, 4.0e-4);
    println!("adiabatic_box: peak_R={peak_r:.3e}  burned={burned:.3}");
    assert!(
        peak_r > EPS_IGNITED,
        "the kernel must establish a burn (R {peak_r:.2e})"
    );
    assert!(
        burned > 0.9,
        "a lit closed ADIABATIC box must burn to completion — a partial burn \
         means a spurious loss channel has appeared (burned {burned:.3})"
    );
}
