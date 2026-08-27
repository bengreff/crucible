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

use crucible_grid::{CellGeom, CellGeomTheta, Grid, GridSpec, Region};
use crucible_solvers::euler::{
    BurnBlendEos, Combustion, Cons, EPS_IGNITED, EosLaw, Euler, EulerFields, FlowBc, FlowBcs, I_RB,
    I_RC, I_RHO, IgnitionColumns, NCOMP, THETA_CELLS, TableEos, consumption_rate, reacting_measure,
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

/// Run-to-run L2 differencer over all conserved fields (the COUP-3 §6.1
/// dt-Richardson pattern, copied in from coup3_sdc.rs — same-layout grids).
fn field_l2_diff(a: &Grid, b: &Grid, f: &EulerFields) -> f64 {
    let mut acc = 0.0f64;
    for (ba, bb) in a.bricks().iter().zip(b.bricks()) {
        for id in f.ids() {
            for (x, y) in ba.field(id).iter().zip(bb.field(id)) {
                acc += (x - y) * (x - y);
            }
        }
    }
    acc.sqrt()
}

#[test]
fn class_r_temporal_order_holds_at_mid_stiffness() {
    // SOLV-4 §3.6 v0.4.8 (the recorded S8 carry) + COUP-3 §6.1: the class-R
    // dt-Richardson temporal-order gate at MID stiffness (dt/τ ~ 1) — the
    // band where the backward-Euler truncation is the OBSERVABLE error. The
    // parking test (extreme stiffness, exact fixed point) and the mild-limit
    // test (explicit reduction) bracket the ends; both are blind to the
    // middle, where the composed SDC step (predictor + 2 trapezoid sweeps)
    // must lift the BE node solve to 2nd order. Fixed grid, one fixed
    // t_final, three step counts (N, 2N, 4N): the run-to-run differences
    // cancel the identical spatial error and isolate the temporal one —
    // order = log2(e1/e2).
    //
    // Fixture: the stiff test's open superheated tube held in-band by dt
    // sizing (fixed dt = t_final/n, NOT the CFL clock — far below the
    // acoustic bound, so no guard is near tripping) and stopped well before
    // parking (the order genuinely collapses in the saturated regime where
    // every trajectory lands on the same cap).
    const P_CHAMBER: f64 = 1.5e6; // 15 bar — τ_ign in the acoustic-dt decade
    const H_MID: f64 = 5.5e6; // T_u ≈ 1888 K: τ = 1.758e-7 s (printed below)
    /// Coarsest run's step count; 2N and 4N refine it.
    const N_BASE: usize = 3;
    /// March horizon in units of the ENTRY τ: long enough for real reactive
    /// content (measured mean burn 0.916), short enough that no cell parks
    /// at the cap (measured: a 3.2τ horizon saturates at b ≈ 0.996 and the
    /// order degrades to 1.72 — the parked regime the doc excludes). τ
    /// itself collapses through the run (reaction-driven compression
    /// heating), so the marched nodes sweep the LOCAL dt/τ from the printed
    /// entry value up through ≫ 1 — the mid band is genuinely traversed.
    const T_FINAL_OVER_TAU: f64 = 1.2;

    // τ at the fixture state, via the operator's own query path.
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
    let w0 = blend
        .cons_from_phzb(P_CHAMBER, H_MID, Z0, 0.0, [0.0, 0.0, 0.0])
        .and_then(|u| blend.prim_checked(&u))
        .unwrap();
    let t_u = blend.unburnt_temperature(&w0).unwrap();
    let tau = comb.ignition.induction_time(P_CHAMBER, t_u, Z0).unwrap();
    let t_final = T_FINAL_OVER_TAU * tau;
    let dt_coarse = t_final / N_BASE as f64;

    let run = |n_steps: usize| -> Grid {
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
        for (k, &id) in f.ids().iter().enumerate() {
            g.fill_field(id, |_r, _th, _z| {
                blend
                    .cons_from_phzb(P_CHAMBER, H_MID, Z0, 0.0, [0.0, 0.0, 0.0])
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
        // Fixture self-check: the fixed Richardson dt must sit far below the
        // acoustic CFL bound (no guard near tripping — the doc's constraint).
        let dt_cfl = sdc.stable_dt(&g, &flow, 0.4).expect("dt");
        let dt = t_final / n_steps as f64;
        assert!(
            dt < 0.75 * dt_cfl,
            "Richardson dt {dt:.3e} too close to the CFL bound {dt_cfl:.3e}"
        );
        let mut t = 0.0;
        for _ in 0..n_steps {
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
            .expect("mid-stiffness class-R step");
            t += dt;
        }
        g
    };

    let (a, b, c) = (run(N_BASE), run(2 * N_BASE), run(4 * N_BASE));
    let (_, f) = tube(8, 4.0e-3);
    // The band + no-parking preconditions, measured on the coarsest run.
    let burned = mean_burn(&a, &f);
    let ids = f.ids();
    let mut b_max = 0.0f64;
    for bi in 0..a.n_bricks() {
        let br = a.brick(bi);
        let mask = br.mask();
        let (rho, rhob) = (br.field(ids[I_RHO]), br.field(ids[I_RB]));
        for local in 0..64 {
            if mask & (1u64 << local) != 0 {
                let cell = br.cell_index(0, local);
                b_max = b_max.max(rhob[cell] / rho[cell]);
            }
        }
    }
    let e1 = field_l2_diff(&a, &b, &f);
    let e2 = field_l2_diff(&b, &c, &f);
    let order = (e1 / e2).log2();
    println!(
        "class_r_order: T_u={t_u:.0} K  τ={tau:.3e} s  dt_coarse/τ={:.2}  \
         burned={burned:.3}  b_max={b_max:.4}  e1={e1:.3e}  e2={e2:.3e}  order={order:.3}",
        dt_coarse / tau
    );
    // Mid-stiffness fixture self-check: the coarsest run's dt/τ must sit in
    // the ~1 band (the mid-band the extreme-end tests were blind to).
    let stiff = dt_coarse / tau;
    assert!(
        (0.2..=3.0).contains(&stiff),
        "fixture not in the mid-stiffness band (dt/τ = {stiff:.2})"
    );
    // Stopped before parking: order measurement is meaningless once cells
    // saturate at the cap (every trajectory lands on the same fixed point).
    assert!(
        b_max < 1.0 - 10.0 * crucible_solvers::euler::BURN_COMPLETE,
        "fixture parked (b_max {b_max:.6}) — the order gate needs a live transient"
    );
    assert!(
        burned > 0.1,
        "fixture carries no reactive content (burned {burned:.3})"
    );
    // Measured: 1.80 at the committed fixture (2.09 at a 1.0τ horizon, 1.84
    // at N_BASE = 4 — consistently ~2, degrading only toward saturation), so
    // the band pins genuine 2nd order without faking tightness the composed
    // nonlinear step does not claim.
    assert!(
        (1.5..=2.8).contains(&order),
        "class-R temporal order {order:.3} outside [1.5, 2.8] at dt/τ = {stiff:.2} \
         (SOLV-4 §3.6 v0.4.8 + COUP-3 §6.1: the composed step must hold ~2nd order \
         where the class-R truncation dominates)"
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

#[test]
fn n_tau_refreeze_node_lag_stays_within_the_recorded_envelope() {
    // SOLV-4 §3.6 v0.4.8 (the recorded S8 carry): the measured N_TAU_REFREEZE
    // witness in the Δt/τ ~ 1 mid-band — and the S9 FINDING it produced.
    // The 0.4.4 sizing premise ("τ depends on the unknown only weakly —
    // Δp/p per node is CFL-bounded") is FALSE in the reaction-driven-
    // compression band: within ONE node solve the burn's constant-volume
    // compression heating drives T_u 1057 → 1512 K and τ down ×1/93 at this
    // fixture state, so the 2-refreeze solve lags the converged frozen-τ
    // fixed point by rel Δ ≈ 3.0e-1 at w/τ = 0.3 (2.2e-2 at 1, 2.1e-3 at 3;
    // the reference converges monotonically to machine precision by pass
    // ~13, so the witness is well-posed). The contract is therefore
    // RESTATED (SOLV-4 0.4.8 as amended): the fixed refreeze count is
    // STRUCTURE, not a convergence claim — the per-node lag is a
    // temporal-truncation-class term absorbed by the SDC sweeps' own
    // re-evaluations, and ACCURACY is owned by the composed-order gate
    // (`class_r_temporal_order_holds_at_mid_stiffness`, measured 1.8–2.1
    // through this same band) — exactly the S3 truncated-Picard idiom.
    // This test pins the MEASURED lag envelope so a regression that grows
    // the lag (a τ-surface change, a projection change) fails loudly.
    //
    // The production 2-refreeze node solve is compared against the SAME
    // frozen-τ BE fixed-point iteration run to convergence (a faithful
    // replica of the implicit_auto_update loop — same prim_checked
    // projection, same floor guards, same induction-time query path, same
    // closed form and cap).
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
    const H_HOT: f64 = 2.5e6; // T_u ≈ 1050 K: real finite τ_ign (the parks-test state)
    /// Reference pass count: far past N_TAU_REFREEZE = 2. Measured: the
    /// frozen-τ map converges monotonically to machine precision by pass
    /// ~13 at every probed band point, so 20 passes ARE the fixed point.
    const N_REF_ITERS: usize = 20;
    /// The mid-stiffness band (w/τ): below, at, and above the BE knee.
    const W_OVER_TAU: [f64; 3] = [0.3, 1.0, 3.0];
    /// The recorded lag envelope per band point (measured 2026-08-25:
    /// {3.022e-1, 2.242e-2, 2.077e-3} at w/τ = {0.3, 1, 3}, ×~1.3 headroom
    /// for fixture drift). NOT an accuracy claim — the composed-order gate
    /// owns accuracy; this pins the measurement so lag GROWTH is caught.
    const LAG_ENVELOPE: [f64; 3] = [4.0e-1, 3.0e-2, 3.0e-3];

    let u = blend
        .cons_from_phzb(P0, H_HOT, Z0, 0.0, [0.0, 0.0, 0.0])
        .unwrap();
    let w0 = blend.prim_checked(&u).unwrap();
    let t_u = blend.unburnt_temperature(&w0).unwrap();
    // τ at the state, via the same induction-time query path the solve uses.
    let tau = comb.ignition.induction_time(P0, t_u, Z0).unwrap();
    let rho = u[I_RHO];
    let cap = rho * (1.0 - crucible_solvers::euler::BURN_COMPLETE);

    for (bi, r) in W_OVER_TAU.into_iter().enumerate() {
        let w_new = r * tau;
        // The production solve: N_TAU_REFREEZE = 2.
        let (x_2, _) = comb.implicit_auto_update(&u, 0.0, w_new).unwrap();
        // The reference: the identical frozen-τ BE fixed-point map, run to
        // N_REF_ITERS passes (combustion.rs implicit_auto_update, replicated
        // faithfully; base = 0 is inside the invariant set, so the entry
        // projection is the identity here).
        let base = 0.0f64;
        let mut x_ref = base;
        for _ in 0..N_REF_ITERS {
            let mut ut_ref = u;
            ut_ref[I_RB] = x_ref;
            let wp = blend.prim_checked(&ut_ref).unwrap();
            // The solve's floor guards, replicated: neither may trip at this
            // hot fixture state (they would return `base` and the witness
            // would compare guards, not refreeze convergence).
            assert!(
                !blend.below_unburnt_floor(&wp),
                "fixture fell below the unburnt h-floor mid-iteration"
            );
            let (p, z) = (wp[4], wp[I_RC]);
            let t_u_i = blend.unburnt_temperature(&wp).unwrap();
            assert!(
                !comb.ignition.non_reactive_floor(p, t_u_i),
                "fixture fell below the non-reactive floor mid-iteration"
            );
            let tau_i = comb.ignition.induction_time(p, t_u_i, z).unwrap();
            // BE at frozen τ: x = (base + w·ρ/τ)/(1 + w/τ), parked at the cap.
            x_ref = ((base + w_new * rho / tau_i) / (1.0 + w_new / tau_i)).min(cap);
            if x_ref == cap {
                break; // parked: further τ refreshes cannot move it
            }
        }
        let delta = (x_2 - x_ref).abs() / x_ref.abs().max(f64::MIN_POSITIVE);
        println!(
            "refreeze_witness: w/τ={r}  x_2/ρ={:.15}  x_ref/ρ={:.15}  rel Δ={delta:.3e}",
            x_2 / rho,
            x_ref / rho
        );
        assert!(
            delta < LAG_ENVELOPE[bi],
            "N_TAU_REFREEZE = 2 node lag at w/τ = {r} grew past the recorded \
             envelope (rel Δ = {delta:.3e} vs pinned {:.1e}; SOLV-4 §3.6 v0.4.8 \
             as amended: the lag is a truncation-class term absorbed by the SDC \
             sweeps — accuracy is owned by the composed-order gate — but GROWTH \
             here means the τ surface or the projection changed character)",
            LAG_ENVELOPE[bi]
        );
    }
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

// --- S11 / ◆C3: combustion on a REVOLVED CUT world at N_θ > 1 ----------------
// The S9-standing refusal ("combustion on cut θ > 1 worlds") is retired here
// for REVOLVED (θ-uniform) cut walls — every world the ◆C3 engine builds:
// the operator's D_c stencil now reads per-sector κ + apertures (kv, the
// meridional faces, and the θ-face aperture that weights the ring flux). A
// genuinely θ-VARYING wall (CSG/STL) stays refused (test below).

const CUT_H: f64 = 1.0e-3;
/// Wall at 3.7 cells: ring 3 is a κ ≈ 0.67 gas sliver (cut, but above the SRD
/// threshold so this gate isolates the combustion θ-faces), ring 4+ is wall.
const CUT_WALL_R: f64 = 3.7 * CUT_H;

fn cut_cyl_rz(i_r: usize) -> (f64, [f64; 4]) {
    let r0 = i_r as f64 * CUT_H;
    let r1 = (i_r as f64 + 1.0) * CUT_H;
    let x = CUT_WALL_R.clamp(r0, r1);
    let kappa = ((x * x - r0 * r0) / (r1 * r1 - r0 * r0)).clamp(0.0, 1.0);
    let ap_r = |rf: f64| if CUT_WALL_R > rf { 1.0 } else { 0.0 };
    (kappa, [ap_r(r0), ap_r(r1), kappa, kappa])
}

fn cut_cyl_spec(n_z: usize, n_theta: u32) -> GridSpec {
    GridSpec {
        r_min: 0.0,
        dr: CUT_H,
        n_r: 5,
        z_min: 0.0,
        dz: CUT_H,
        n_z,
        n_theta_max: n_theta,
        axisymmetry_assertion: n_theta == 1,
    }
}

fn cut_cyl_grid_1(n_z: usize) -> Grid {
    Grid::build_with_geometry(
        cut_cyl_spec(n_z, 1),
        crucible_solvers::euler::EULER_FIELDS,
        |i_r, _| {
            let (kappa, aperture) = cut_cyl_rz(i_r);
            let region = if kappa > 0.0 {
                Region::Gas
            } else {
                Region::Exterior
            };
            CellGeom {
                region,
                kappa,
                aperture,
            }
        },
    )
    .expect("revolved cut cylinder at N_θ = 1")
}

fn cut_cyl_grid_theta(n_z: usize, n_theta: u32) -> Grid {
    Grid::build_with_geometry_theta(
        cut_cyl_spec(n_z, n_theta),
        crucible_solvers::euler::EULER_FIELDS,
        |i_r, _j, _i_z| {
            let (kappa, ap) = cut_cyl_rz(i_r);
            CellGeomTheta {
                kappa,
                aperture: [ap[0], ap[1], ap[2], ap[3], kappa, kappa],
            }
        },
        |_, _| Region::Exterior,
    )
    .expect("revolved cut cylinder through the θ-builder")
}

const H_CUT_WARM: f64 = 1.3e6;

/// An axisymmetric burning IC: a z-front (burnt below, unburnt above) in warm
/// reactants — θ-uniform, so the reaction + meridional front diffusion are
/// live but the θ-direction flux is identically zero.
fn fill_axisym_front(g: &mut Grid, f: &EulerFields, blend: &BurnBlendEos<'_>) {
    let ids = f.ids();
    for (k, &id) in ids.iter().enumerate() {
        g.fill_field(id, |_r, _th, z| {
            let b0 = if z < 3.0 * CUT_H { 1.0 } else { 0.0 };
            blend
                .cons_from_phzb(P0, H_CUT_WARM, Z0, b0, [0.0, 0.0, 0.0])
                .expect("valid axisym IC")[k]
        });
    }
}

// THE gate: an axisymmetric burn on a revolved CUT world at N_θ = 8 stays
// exactly θ-symmetric (every sector bit-identical to sector 0 — the per-sector
// D_c stencil and the zero θ-flux introduce no spurious asymmetry) and reduces
// to the certified N_θ = 1 cut march (tight tolerance early; unlike the flow
// sweep, the per-θ face-area accumulation is not power-of-two bit-reducible, so
// the sub-ULP round-off amplifies chaotically — different N_θ discretizations
// need only converge, not bit-match). Marched through the real SDC step
// (class A + class R), audit armed.
#[test]
fn s11_revolved_cut_combustion_n_theta_8_is_theta_symmetric_and_reduces() {
    let ut = unburnt();
    let bt = burnt();
    let it = ignition();
    let blend = BurnBlendEos::new(TableEos::bind(&ut).unwrap(), TableEos::bind(&bt).unwrap());
    let comb = Combustion {
        blend: &blend,
        ignition: IgnitionColumns::bind(&it).unwrap(),
        wrinkling: 1.0,
        theta: THETA_CELLS,
    };
    let n_z = 6usize;
    let mut g1 = cut_cyl_grid_1(n_z);
    let mut g8 = cut_cyl_grid_theta(n_z, 8);
    // The fixture must actually be cut (a partial ring) yet θ-uniform.
    assert!(g8.has_cut_geometry(), "fixture must carry cut geometry");
    assert!(
        g8.geometry_is_theta_uniform(),
        "a revolved cut world must read θ-uniform"
    );
    assert!(
        (0.6..0.75).contains(&g8.kappa_at(3, 0, 0)),
        "ring 3 must be the declared κ ≈ 0.67 sliver"
    );
    let f1 = EulerFields::resolve(&g1).expect("fields");
    let f8 = EulerFields::resolve(&g8).expect("fields");
    fill_axisym_front(&mut g1, &f1, &blend);
    fill_axisym_front(&mut g8, &f8, &blend);
    let op1 = Euler {
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
    let op8 = Euler {
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
    let flow1 = FlowClass {
        op: &op1,
        fields: &f1,
    };
    let flow8 = FlowClass {
        op: &op8,
        fields: &f8,
    };
    let react1 = ReactionClass { op: &comb };
    let react8 = ReactionClass { op: &comb };
    let (mut sdc1, mut sdc8) = (Sdc::new(), Sdc::new());
    let mut t = 0.0;
    for step in 0..24 {
        let dt = sdc8.stable_dt(&g8, &flow8, 0.4).expect("dt");
        sdc8.step(
            &mut g8,
            Some(&flow8),
            None,
            None,
            None,
            Some(&react8),
            t,
            dt,
        )
        .expect("audited 3-D cut combustion step");
        sdc1.step(
            &mut g1,
            Some(&flow1),
            None,
            None,
            None,
            Some(&react1),
            t,
            dt,
        )
        .expect("audited 2-D cut combustion step");
        t += dt;
        let ids = f8.ids();
        for bi in 0..g8.n_bricks() {
            let b1 = g1.brick(bi);
            let b8 = g8.brick(bi);
            let nt = b8.n_theta();
            for &id in ids.iter() {
                let d1 = b1.field(id);
                let d8 = b8.field(id);
                for local in 0..crucible_grid::BRICK_CELLS {
                    // (A) THE θ-symmetry gate (bitwise, every step): an
                    // axisymmetric burn must keep every sector bit-identical
                    // to sector 0 — the per-sector D_c stencil and the zero
                    // θ-flux introduce NO spurious azimuthal asymmetry. A bug
                    // in the θ-faces would break this immediately.
                    let v0 = d8[b8.cell_index(0, local)];
                    for j in 1..nt {
                        assert_eq!(
                            v0.to_bits(),
                            d8[b8.cell_index(j, local)].to_bits(),
                            "sector {j} of brick {bi} field {id:?} cell {local} broke θ-symmetry \
                             at step {step} — a spurious azimuthal asymmetry"
                        );
                    }
                    // (B) Reduction to the certified N_θ = 1 cut march, checked
                    // in the first steps before chaotic amplification of the
                    // sub-ULP round-off differences (the per-θ face-area
                    // accumulation is not power-of-two bit-reducible the way
                    // the flow sweep is; the physics is identical).
                    if step < 3 {
                        let v1 = d1[b1.cell_index(0, local)];
                        let denom = v1.abs().max(1e-6);
                        assert!(
                            (v0 - v1).abs() / denom < 1e-9,
                            "brick {bi} field {id:?} cell {local}: N_θ = 8 {v0:.9e} vs \
                             N_θ = 1 {v1:.9e} at step {step} — combustion did not reduce"
                        );
                    }
                }
            }
        }
    }
}

// A genuinely θ-VARYING cut wall stays refused (the CSG/STL wave, plan S11):
// per-sector apertures differ, so the D_c stencil's shared-face assumptions
// are not yet the whole story — refuse loudly, cell/owner named.
#[test]
fn s11_theta_varying_cut_combustion_still_refuses() {
    let ut = unburnt();
    let bt = burnt();
    let it = ignition();
    let blend = BurnBlendEos::new(TableEos::bind(&ut).unwrap(), TableEos::bind(&bt).unwrap());
    let comb = Combustion {
        blend: &blend,
        ignition: IgnitionColumns::bind(&it).unwrap(),
        wrinkling: 1.0,
        theta: THETA_CELLS,
    };
    // All-gas, but κ + apertures VARY by sector — θ-varying cut geometry.
    let nt = 8u32;
    let spec = GridSpec {
        r_min: 0.5 * CUT_H,
        dr: CUT_H,
        n_r: 4,
        z_min: 0.0,
        dz: CUT_H,
        n_z: 4,
        n_theta_max: nt,
        axisymmetry_assertion: false,
    };
    let mut g = Grid::build_with_geometry_theta(
        spec,
        crucible_solvers::euler::EULER_FIELDS,
        |_i_r, j, _i_z| {
            // κ dips in one sector (all gas) — θ-VARYING — while every face
            // aperture stays 1.0 so the θ-pair + shared-face coherence rules
            // are satisfied and the world builds. Only κ varies by sector.
            let kappa = if j == 3 { 0.6 } else { 0.8 };
            CellGeomTheta {
                kappa,
                aperture: [1.0; 6],
            }
        },
        |_, _| Region::Exterior,
    )
    .expect("θ-varying cut world builds");
    assert!(g.has_cut_geometry());
    assert!(
        !g.geometry_is_theta_uniform(),
        "the fixture must read θ-VARYING"
    );
    let f = EulerFields::resolve(&g).expect("fields");
    for (k, &id) in f.ids().iter().enumerate() {
        g.fill_field(id, |_r, _th, _z| {
            blend
                .cons_from_phzb(P0, H_CUT_WARM, Z0, 0.5, [0.0, 0.0, 0.0])
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
    let react = ReactionClass { op: &comb };
    let err = Sdc::new()
        .step(
            &mut g,
            Some(&flow),
            None,
            None,
            None,
            Some(&react),
            0.0,
            1.0e-9,
        )
        .expect_err("θ-varying cut combustion must refuse");
    let msg = err.to_string();
    assert!(
        msg.contains("VARYING") || msg.contains("θ-VARYING") || msg.contains("combustion"),
        "the refusal must name the θ-varying-combustion owner, got: {msg}"
    );
}
