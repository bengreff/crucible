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
use crucible_solvers::gas_diffusion::{
    FaceGasBc, GasDiffBcs, GasDiffusion, SpeciesBc, ThermalBc, VelocityBc,
};
use crucible_solvers::sdc::{FlowClass, GasDiffusionClass, ReactionClass, Sdc};
use crucible_solvers::transport::{ConstantTransport, TransportProps};
use crucible_tables::{Pin, Table};
use crucible_units::{dynamic_viscosity_pa_s, specific_heat_capacity_j_per_kg_k};

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
    let reaction = ReactionClass { op: &comb };
    let mut sdc = Sdc::new();
    let mut t = 0.0;
    let mut steps = 0usize;
    while t < t_settle {
        let dt = sdc.stable_dt(&g, &flow, 0.4).expect("dt").min(t_settle - t);
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
    let reaction = ReactionClass { op: &comb };
    let mut sdc = Sdc::new();
    let (mut t, mut t_ign) = (0.0f64, f64::NAN);
    while t < 5.0e-3 {
        let dt = sdc.stable_dt(&g, &flow, 0.4).unwrap();
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
    let reaction = ReactionClass { op: &comb };
    let mut sdc = Sdc::new();
    let (mut t, mut peak_r) = (0.0f64, 0.0f64);
    let ids = f.ids();
    while t < t_end {
        let dt = sdc.stable_dt(&g, &flow, 0.4).unwrap().min(t_end - t);
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

/// March a uniform box at `(p, h, Z, b = 0)` with the reaction class
/// scheduled, no igniter; return `(final mean burnt fraction, max b)`.
fn uniform_reactive_march(p: f64, h_spec: f64, n_steps: usize) -> (f64, f64) {
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
    let (mut g, f) = tube(8, 5.0e-4);
    let ids = f.ids();
    for (k, &id) in ids.iter().enumerate() {
        g.fill_field(id, |_r, _th, _z| {
            blend
                .cons_from_phzb(p, h_spec, Z0, 0.0, [0.0, 0.0, 0.0])
                .expect("valid uniform IC")[k]
        });
    }
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
    let mut t = 0.0;
    for _ in 0..n_steps {
        let dt = sdc.stable_dt(&g, &flow, 0.4).unwrap();
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
        .unwrap();
        t += dt;
    }
    let mut b_max = 0.0f64;
    for bi in 0..g.n_bricks() {
        let br = g.brick(bi);
        let mask = br.mask();
        let (rho, rhob) = (br.field(ids[I_RHO]), br.field(ids[I_RB]));
        for local in 0..64 {
            if mask & (1u64 << local) != 0 {
                let cell = br.cell_index(0, local);
                b_max = b_max.max(rhob[cell] / rho[cell]);
            }
        }
    }
    (mean_burn(&g, &f), b_max)
}

#[test]
fn stiff_auto_ignition_is_stable_past_the_explicit_bound() {
    // THE class-R acceptance (COUP-3 §3.3, S7): superheated reactants at
    // elevated pressure have τ_ign at/below the acoustic Δt — dt/τ ≳ 1,
    // exactly where the S6 explicit advance overshot the positivity guard
    // and halted (Δb = dt·(1−b)/τ > 1 in one step). The implicit node solve
    // must march it: unconditionally stable, the burn completing on the
    // guarded law's own fixed point (b parks at 1 − BURN_COMPLETE from
    // below, never past). The tube is OPEN (constant-pressure outflow at the
    // fill pressure) — an engine is open, so the fixture must not compress
    // itself into confined-blast corners no chamber reaches.
    const P_CHAMBER: f64 = 1.5e6; // 15 bar — τ_ign shrinks with p
    const H_SUPER: f64 = 5.5e6; // superheated reactants (T_u ≈ 1890 K)
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
    let (mut g, f) = tube(8, 4.0e-3);
    let ids = f.ids();
    for (k, &id) in ids.iter().enumerate() {
        g.fill_field(id, |_r, _th, _z| {
            blend
                .cons_from_phzb(P_CHAMBER, H_SUPER, Z0, 0.0, [0.0, 0.0, 0.0])
                .expect("valid superheated IC")[k]
        });
    }
    let p_amb = |_t: f64| P_CHAMBER;
    let op = Euler {
        eos: blend.clone(),
        source: &ZERO_SRC,
        bcs: FlowBcs {
            r_inner: FlowBc::Reflecting,
            r_outer: FlowBc::Reflecting,
            z_lo: FlowBc::Reflecting,
            z_hi: FlowBc::PressureOutflow(&p_amb),
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
    // Self-verify the stiffness claim: this fixture must genuinely sit past
    // the explicit bound (dt/τ > 1 ⇒ the S6 explicit advance would halt).
    let w0 = blend
        .cons_from_phzb(P_CHAMBER, H_SUPER, Z0, 0.0, [0.0, 0.0, 0.0])
        .and_then(|u| blend.prim_checked(&u))
        .unwrap();
    let t_u = blend.unburnt_temperature(&w0).unwrap();
    let tau = comb.ignition.induction_time(P_CHAMBER, t_u, Z0).unwrap();
    let dt0 = sdc.stable_dt(&g, &flow, 0.4).unwrap();
    println!(
        "stiff_auto_ignition: T_u={t_u:.0} K  τ={tau:.3e} s  dt={dt0:.3e} s  dt/τ={:.2}",
        dt0 / tau
    );
    assert!(
        dt0 / tau > 1.0,
        "fixture not stiff enough to exercise the class-R claim (dt/τ = {:.2})",
        dt0 / tau
    );
    let mut t = 0.0;
    for _ in 0..12 {
        let dt = sdc.stable_dt(&g, &flow, 0.4).unwrap();
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
        .expect("stiff auto-ignition step (the S6 explicit tier halted here)");
        t += dt;
    }
    let burned = mean_burn(&g, &f);
    let mut b_max = 0.0f64;
    for bi in 0..g.n_bricks() {
        let br = g.brick(bi);
        let mask = br.mask();
        let (rho, rhob) = (br.field(ids[I_RHO]), br.field(ids[I_RB]));
        for local in 0..64 {
            if mask & (1u64 << local) != 0 {
                let cell = br.cell_index(0, local);
                b_max = b_max.max(rhob[cell] / rho[cell]);
            }
        }
    }
    println!("stiff_auto_ignition: burned={burned:.6}  b_max={b_max:.7}");
    assert!(
        burned > 0.95,
        "a superheated tube must auto-ignite essentially instantly (burned {burned:.4})"
    );
    assert!(
        b_max <= 1.0 - crucible_solvers::euler::BURN_COMPLETE + 1.0e-12,
        "the implicit solve must park AT the guarded fixed point, never past it \
         (b_max {b_max:.9})"
    );
}

#[test]
fn implicit_node_solve_parks_exactly_at_any_stiffness() {
    // The class-R node solve as a pure function, at stiffness no march can
    // reach (w/τ ~ 10⁶): the update must land EXACTLY on the guarded law's
    // parking point ρ(1 − BURN_COMPLETE) — the exact integral of the
    // declared discontinuous rate law — and the mild limit must reduce to
    // the explicit rate to O((w/τ)²).
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
    const H_HOT: f64 = 2.5e6; // T_u ≈ 1050 K: real finite τ_ign
    let u = blend
        .cons_from_phzb(P0, H_HOT, Z0, 0.0, [0.0, 0.0, 0.0])
        .unwrap();
    let w = blend.prim_checked(&u).unwrap();
    let t_u = blend.unburnt_temperature(&w).unwrap();
    let tau = comb.ignition.induction_time(P0, t_u, Z0).unwrap();
    let rho = u[0];
    let cap = rho * (1.0 - crucible_solvers::euler::BURN_COMPLETE);
    // Extreme stiffness: parks exactly at the fixed point.
    let (x, r) = comb.implicit_auto_update(&u, 0.0, 1.0e6 * tau).unwrap();
    assert!(
        x == cap,
        "extreme-stiffness update must park exactly at ρ(1−BURN_COMPLETE): \
         x = {x:.12e} vs cap = {cap:.12e}"
    );
    assert!(r > 0.0 && r.is_finite());
    // Mild limit: reduces to the explicit rate. The residual deviation is
    // NOT the BE factor (that is O(w/τ) ~ 1e-3 here) but the τ-refreeze
    // evaluating τ at the advanced state — τ is Arrhenius-steep in the
    // enthalpy rise the tiny burn causes (measured d(ln τ)/db ≈ 12 at this
    // state), which is the backward-Euler semantics working as declared.
    // Assert same-value-to-a-few-percent, which pins the reduction without
    // faking a tighter identity than the scheme claims.
    let w_mild = 1.0e-3 * tau;
    let (x_mild, _) = comb.implicit_auto_update(&u, 0.0, w_mild).unwrap();
    let explicit = w_mild * rho / tau;
    let rel = ((x_mild - explicit) / explicit).abs();
    assert!(
        rel < 5.0e-2,
        "mild-limit update must reduce to the explicit rate (rel dev {rel:.2e})"
    );
    // A base at/past the cap is a zero-source fixed point.
    let (x_at, r_at) = comb.implicit_auto_update(&u, cap, tau).unwrap();
    assert!(x_at == cap && r_at == 0.0);
}

/// March a lit flame tube WITH the class-D gas-diffusion coupling — the
/// blend↔diffusion march S7's startup runs — under the given wall thermal
/// condition; return `(R_initial, R_final, final mean burnt fraction)`.
/// The tube's r-gap is the quench-scale dimension (n_r cells of `h`).
fn quench_march(p_fill: f64, wall_t_k: Option<f64>, n_steps: usize) -> (f64, f64, f64) {
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
    let (n_z, h) = (48, 2.0e-4); // r-gap = 4·h = 0.8 mm (the quench scale)
    let (mut g, f) = tube(n_z, h);
    let ids = f.ids();
    for (k, &id) in ids.iter().enumerate() {
        g.fill_field(id, |_r, _th, z| {
            let b = 1.0 / (1.0 + ((z - 0.5 * (n_z as f64) * h) / (THETA_CELLS * h)).exp());
            blend
                .cons_from_phzb(p_fill, H0, Z0, b, [0.0, 0.0, 0.0])
                .expect("valid flame IC")[k]
        });
    }
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
    // The class-D gas occupant: viscous walls at the r-gap faces carrying
    // the declared wall thermal condition; closed adiabatic z ends.
    let no_slip = |_r: f64, _th: f64, _z: f64, _t: f64| (0.0, 0.0, 0.0);
    let t_wall = move |_r: f64, _th: f64, _z: f64, _t: f64| wall_t_k.unwrap_or(f64::NAN);
    let thermal = |on: bool| -> ThermalBc<'_> {
        if on {
            ThermalBc::Isothermal(&t_wall)
        } else {
            ThermalBc::Adiabatic
        }
    };
    let gas_op = GasDiffusion::new(GasDiffBcs {
        r_inner: FaceGasBc {
            velocity: VelocityBc::NoSlip(&no_slip),
            thermal: thermal(wall_t_k.is_some()),
            species: SpeciesBc::ZeroFlux,
        },
        r_outer: FaceGasBc {
            velocity: VelocityBc::NoSlip(&no_slip),
            thermal: thermal(wall_t_k.is_some()),
            species: SpeciesBc::ZeroFlux,
        },
        z_lo: FaceGasBc {
            velocity: VelocityBc::FreeSlip,
            thermal: ThermalBc::Adiabatic,
            species: SpeciesBc::ZeroFlux,
        },
        z_hi: FaceGasBc {
            velocity: VelocityBc::FreeSlip,
            thermal: ThermalBc::Adiabatic,
            species: SpeciesBc::ZeroFlux,
        },
    });
    // The blend↔class-D seam (S7): the gas operator reads the BLEND's
    // mass-weighted temperature and the spine at the blend's own
    // interrogation coordinate — here the declared-constant occupant with
    // H₂/O₂-mixture-class values (a state-varying spine would not change
    // what this fixture demonstrates: the loss channel).
    let props: TransportProps = ConstantTransport::new(
        specific_heat_capacity_j_per_kg_k(3.8e3),
        dynamic_viscosity_pa_s(2.0e-5),
        0.7, // Pr
        1.4,
        0.5, // Sc
    )
    .expect("transport set")
    .into_props();
    let temperature = |w: &crucible_solvers::euler::Prim| blend.temperature_w(w);
    let transport =
        move |_w: &crucible_solvers::euler::Prim| -> Result<TransportProps, &'static str> {
            Ok(props)
        };
    let gas = GasDiffusionClass {
        op: &gas_op,
        temperature: &temperature,
        transport: &transport,
    };
    let mut sdc = Sdc::new();
    let r_initial = reacting_measure(&g, &ids, &blend, &comb.ignition, 1.0, THETA_CELLS).unwrap();
    let mut t = 0.0;
    for _step in 0..n_steps {
        let dt = sdc.stable_dt(&g, &flow, 0.4).unwrap();
        sdc.step(
            &mut g,
            Some(&flow),
            None,
            Some(&gas),
            None,
            Some(&reaction),
            t,
            dt,
        )
        .expect("blend + class-D quench march (audit armed)");
        t += dt;
    }
    let r_final = reacting_measure(&g, &ids, &blend, &comb.ignition, 1.0, THETA_CELLS).unwrap();
    (r_initial, r_final, mean_burn(&g, &f))
}

#[test]
fn quench_box_flames_out_on_cold_walls() {
    // THE FLAMEOUT demonstration (SOLV-4 §6.5, re-scoped S6→S7): quenching
    // is a heat-loss phenomenon, so it needs the loss channel — conductive
    // loss to cold isothermal walls through the class-D gas-diffusion
    // coupling (the blend↔diffusion march the startup runs anyway). At
    // 0.8 mm gap and reduced pressure (quenching distance ∝ 1/p: ~0.2 mm at
    // 1 atm for H₂/O₂, ~the gap at ~0.25 atm), a lit front between cold
    // walls must DIE: the reacting measure collapses after having exceeded
    // the floor — the exact R-trajectory the COUP-4 FLAMEOUT verdict
    // consumes. The adiabatic control on the same fixture keeps burning
    // (the S6 converse), so the loss channel — not the fixture — is what
    // kills the flame.
    // 420 K walls at 0.1 atm: the coldest wall the GAS-PHASE model can
    // honestly quench against — the products surface's declared envelope
    // floor reads ~407 K at the fixture state (grid floor ~332 K; H2O
    // condensation onset ~310 K — measured from the v0.4.0 artifact), so
    // colder walls chase burnt gas off its own surface (ice — the S15
    // two-phase wave's territory). Quenching needs no cold wall, only a
    // loss rate the front cannot outrun: quenching distance ∝ 1/p, so
    // 0.1 atm puts the H2/O2 quench scale well above the 0.8 mm gap.
    let (r0, r_end, burned) = quench_march(1.0e4, Some(420.0), 2600);
    println!("quench_box (cold walls): R {r0:.3e} -> {r_end:.3e} kg/s, burned {burned:.3}");
    assert!(
        r0 > EPS_IGNITED,
        "the initialized front must start alive (R {r0:.2e})"
    );
    assert!(
        r_end < EPS_IGNITED,
        "cold walls at the quench gap must extinguish the front — R stayed at \
         {r_end:.2e} (FLAMEOUT unreached)"
    );
    assert!(
        burned < 0.8,
        "the front must die before consuming the tube (burned {burned:.3})"
    );

    let (c0, c_end, c_burned) = quench_march(1.0e4, None, 2600);
    println!(
        "quench_box (adiabatic control): R {c0:.3e} -> {c_end:.3e} kg/s, burned {c_burned:.3}"
    );
    // The control discriminator is the reacting measure, not the burned
    // fraction: at 0.1 atm the laminar front crosses ~half a cell in this
    // march window, so consumption cannot separate the cases — but a LIVE
    // front holds R at its initial scale while the quenched one collapsed
    // to zero. Same fixture, same march, only the loss channel differs.
    assert!(
        c0 == r0,
        "the two fixtures must start from the same state (R {c0:.3e} vs {r0:.3e})"
    );
    assert!(
        c_end > 0.5 * c0,
        "the adiabatic control's front must stay alive (R {c0:.3e} -> {c_end:.3e}) — \
         otherwise the loss channel is not what killed the quenched flame"
    );
    assert!(
        c_burned >= burned,
        "the adiabatic control must burn at least as far as the quenched case \
         (control {c_burned:.3} vs quenched {burned:.3})"
    );
}

#[test]
fn cold_floor_cell_is_declared_non_reactive() {
    // The cold-side non-reactive floor (SOLV-4 §3.6 / OFFL-3 §3.3, S7): a
    // cell below the ignition surface's declared p-envelope floor (here
    // ~1 kPa — near-vacuum fill territory) is declared no-burn: the march
    // neither refuses on the closure query (the S6 behavior this cures) nor
    // burns anything. The real S_L down to the unburnt branch's own cold
    // validity floor is carried by the surface's cold rows (ignition
    // v0.3.0), not by this guard.
    const P_FILL: f64 = 1.0e3; // below the ignition surface's p floor
    const H_COLD: f64 = -1.0e5; // cold reactants (~250 K class)
    let (burned, b_max) = uniform_reactive_march(P_FILL, H_COLD, 5);
    println!("cold_floor: burned={burned:.3e}  b_max={b_max:.3e}");
    assert!(
        burned == 0.0 && b_max == 0.0,
        "a below-floor cell must be exactly non-reactive (burned {burned:.2e}, \
         b_max {b_max:.2e})"
    );
}

// --- S8: the azimuthal re-key + the S_T-CFL --------------------------------

/// A thin annulus whose θ-arcs are at the cell scale (r̄·Δθ ≈ h), so a
/// point-in-θ kernel's azimuthal spread is resolved front physics, not a
/// geometric artifact. Off-axis (r_min > 0): the axis machinery has its own
/// gates (`solv1_azimuthal`).
fn theta_tube(n_z: usize, h: f64, n_theta: u32) -> (Grid, EulerFields) {
    let spec = GridSpec {
        r_min: 0.5 * h,
        dr: h,
        n_r: 4,
        z_min: 0.0,
        dz: h,
        n_z,
        n_theta_max: n_theta,
        axisymmetry_assertion: false,
    };
    let g = Grid::build(spec, crucible_solvers::euler::EULER_FIELDS).expect("grid");
    let f = EulerFields::resolve(&g).expect("fields");
    (g, f)
}

/// Mean burnt fraction per θ-sector: `∫ρb dV / ∫ρ dV` over each j.
fn sector_burn(g: &Grid, f: &EulerFields) -> Vec<f64> {
    let ids = f.ids();
    let nt = g.brick(0).n_theta();
    let (mut mb, mut m) = (vec![0.0f64; nt as usize], vec![0.0f64; nt as usize]);
    for bi in 0..g.n_bricks() {
        let br = g.brick(bi);
        let mask = br.mask();
        let (rho, rhob) = (br.field(ids[I_RHO]), br.field(ids[I_RB]));
        for local in 0..64 {
            if mask & (1u64 << local) == 0 {
                continue;
            }
            let (i_r, _) = br.global_rz(local);
            let v = g.cell_volume(i_r, nt);
            for j in 0..nt {
                let cell = br.cell_index(j, local);
                mb[j as usize] += rhob[cell] * v;
                m[j as usize] += rho[cell] * v;
            }
        }
    }
    mb.iter().zip(&m).map(|(a, b)| a / b).collect()
}

#[test]
fn spark_point_in_theta_spreads_azimuthally() {
    // SOLV-4 0.4.6 (plan S8): a spark kernel that is a POINT in θ — not the
    // S7 ring — lights its sector, and the front's θ-direction D_c faces
    // carry it to the ADJACENT sectors before the opposite one: azimuthal
    // propagation is resolved front physics (◆C3's prerequisite). Marched
    // through the real SDC step (class A + class R), audit armed.
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
    // WARM reactants (T_u ~ 600 K class): S_L is tens of m/s there, so the
    // sector-crossing time is ~50 acoustic steps, not ~800 — the same
    // physics gate at fast-battery cost (VAL-3 §3.2 budget).
    const H_WARM: f64 = 1.3e6;
    let h = 1.0e-3;
    let n_z = 6usize;
    let nt = 8u32;
    let (mut g, f) = theta_tube(n_z, h, nt);
    let ids = f.ids();
    for (k, &id) in ids.iter().enumerate() {
        g.fill_field(id, |_r, _th, _z| {
            blend
                .cons_from_phzb(P0, H_WARM, Z0, 0.0, [0.0, 0.0, 0.0])
                .expect("valid uniform IC")[k]
        });
    }
    // The kernel: θ-sector 2 only, z middle ±2 cells, all r — a point in θ.
    let z_c = 0.5 * (n_z as f64) * h;
    let sector = std::f64::consts::TAU / f64::from(nt);
    let igniter = move |_r: f64, th: f64, z: f64, t: f64| -> Cons {
        let in_theta = (th / sector).floor() as i32 == 2;
        let on = t < 1.2e-5 && in_theta && (z - z_c).abs() < 2.0 * h;
        let q = if on {
            5.0e10 * (t / 4.0e-6).min(1.0)
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
    let reaction = ReactionClass { op: &comb };
    let mut sdc = Sdc::new();
    let (mut t, mut peak_r) = (0.0f64, 0.0f64);
    let mut adjacency_checked = false;
    let mut steps = 0usize;
    while t < 1.5e-4 {
        let dt = sdc
            .stable_dt(&g, &flow, 0.4)
            .expect("dt (σ_front member live)");
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
        .expect("audited 3-D combustion step");
        t += dt;
        steps += 1;
        assert!(steps < 60_000, "runaway spark-point march");
        if !steps.is_multiple_of(4) {
            continue; // probe cadence (slow-clock members — the S7 pattern)
        }
        let r = reacting_measure(&g, &ids, &blend, &comb.ignition, 1.0, THETA_CELLS).unwrap();
        peak_r = peak_r.max(r);
        if !adjacency_checked {
            let s = sector_burn(&g, &f);
            let adj = s[1].min(s[3]);
            if adj >= 0.05 {
                // THE ordering gate, taken at first adjacent light-off: the
                // opposite sector must still be (relatively) dark — the
                // spread is a propagating front, not a uniform volumetric
                // ignition.
                assert!(
                    s[6] < 0.5 * adj,
                    "opposite sector lit with the adjacent ones \
                     (adj {adj:.3} vs opposite {:.3}) — no propagating θ-front",
                    s[6]
                );
                adjacency_checked = true;
            }
        }
    }
    let s = sector_burn(&g, &f);
    println!("spark_point sectors: {s:?}  peak_R={peak_r:.3e}");
    assert!(
        peak_r > EPS_IGNITED,
        "the point spark must light (R {peak_r:.2e})"
    );
    assert!(
        s[2] > 0.2,
        "the sparked sector must burn substantially ({:.3})",
        s[2]
    );
    assert!(
        adjacency_checked,
        "adjacent sectors never reached the light-off threshold — no azimuthal \
         spread (sectors {s:?})"
    );
}

#[test]
fn front_carrier_joins_the_dt_rule_and_separation_guard_refuses() {
    // SOLV-4 0.4.6 (plan S8), the S_T-CFL: (a) with the rate law live, a
    // larger declared wrinkling shortens stable_dt (σ_front is a Δt-rule
    // member — wrinkling > 1 marches honestly instead of tripping the
    // positivity guard); (b) a wrinkling that drives S_T toward the sound
    // speed trips the K_FRONT_SEP scale-separation refusal, typed and
    // cell-named — the fast-deflagration/DDT class is out of the declared
    // model form, never marched through silently.
    let ut = unburnt();
    let bt = burnt();
    let it = ignition();
    let blend = BurnBlendEos::new(TableEos::bind(&ut).unwrap(), TableEos::bind(&bt).unwrap());
    let dt_at = |wrinkling: f64| -> Result<f64, crucible_solvers::sdc::SdcError> {
        let ign = IgnitionColumns::bind(&it).unwrap();
        let comb = Combustion {
            blend: &blend,
            ignition: ign,
            wrinkling,
            theta: THETA_CELLS,
        };
        let (mut g, f) = tube(8, 2.0e-4);
        // A half-burnt mid-flame state: the rate law is live (b well inside
        // (0,1), T_u flammable), so σ_front is nonzero.
        let ids = f.ids();
        for (k, &id) in ids.iter().enumerate() {
            g.fill_field(id, |_r, _th, _z| {
                blend
                    .cons_from_phzb(P0, H0, Z0, 0.5, [0.0, 0.0, 0.0])
                    .expect("valid IC")[k]
            });
        }
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
        Sdc::new().stable_dt(&g, &flow, 0.4)
    };
    let dt_laminar = dt_at(1.0).expect("laminar dt");
    let dt_wrinkled = dt_at(8.0).expect("wrinkling 8 must march (the Δt member absorbs it)");
    println!("S_T-CFL: dt(w=1)={dt_laminar:.3e}  dt(w=8)={dt_wrinkled:.3e}");
    assert!(
        dt_wrinkled < dt_laminar,
        "the front-carrier member must shorten Δt at wrinkling 8 \
         ({dt_wrinkled:.3e} vs {dt_laminar:.3e})"
    );
    // (b) the sonic-end refusal: S_L ≈ 2.6 m/s here, so wrinkling 2e4 puts
    // S_T in the 5e4 m/s class — far past any sound speed.
    match dt_at(2.0e4) {
        Err(e) => {
            let msg = format!("{e}");
            assert!(
                msg.contains("scale separation") || msg.contains("front-carrier"),
                "wrong refusal: {msg}"
            );
        }
        Ok(dt) => panic!("sonic-class S_T must refuse, got dt = {dt:.3e}"),
    }
}

#[test]
fn near_axis_combustion_dt_does_not_spuriously_refuse() {
    // S8 review finding (cured): the first guard compared MESH rates, and the
    // θ-arc's 1/arc² carrier-diffusion rate outgrows the 1/arc acoustic rate —
    // at fine N_θ near the axis it refused mild flames the certified tier
    // marches. The model-form guard (S_T vs the LOCAL sound speed) is
    // resolution-independent: a mid-flame state on an AXIS world at N_θ = 16
    // must produce a finite Δt with no refusal (S_T/c ~ 1e-3 here), while the
    // σ_front member still pays the honest near-axis θ-diffusion cost.
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
    let h = 2.0e-4;
    let spec = GridSpec {
        r_min: 0.0,
        dr: h,
        n_r: 4,
        z_min: 0.0,
        dz: h,
        n_z: 4,
        n_theta_max: 16,
        axisymmetry_assertion: false,
    };
    let mut g = Grid::build(spec, crucible_solvers::euler::EULER_FIELDS).expect("grid");
    let f = EulerFields::resolve(&g).expect("fields");
    let ids = f.ids();
    for (k, &id) in ids.iter().enumerate() {
        g.fill_field(id, |_r, _th, _z| {
            blend
                .cons_from_phzb(P0, H0, Z0, 0.5, [0.0, 0.0, 0.0])
                .expect("valid IC")[k]
        });
    }
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
    let dt = Sdc::new()
        .stable_dt(&g, &flow, 0.4)
        .expect("a mild flame near the axis at N_θ = 16 must not refuse");
    assert!(dt.is_finite() && dt > 0.0);
}
