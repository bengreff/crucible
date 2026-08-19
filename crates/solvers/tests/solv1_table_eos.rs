//! SOLV-1 §3.4 — shifting-equilibrium mode on the *production* LOX/LH₂
//! surface: the per-cell equilibrium projection (fixed-count p-iteration at
//! (p, h = e + p/ρ, Z)) as the [`EosLaw`] occupant of the one field
//! operator. Combustion lives in the EOS: the RL10 chamber state emerges
//! from (ρ, e, Z) alone — no burn branch anywhere.
//!
//! Covers: projection round-trip on and off the design line; consistency of
//! the aux (e, Γ₁) slots; refusal classes (envelope, no-bracket, unsupported
//! BC); the whole-operator seam (uniform equilibrium gas at rest is a
//! bitwise fixed point of the full step; a closed hot/cold tube conserves
//! mass, energy, and elemental Z through real wave dynamics); determinism.

use crucible_grid::{Grid, GridSpec};
use crucible_solvers::euler::{
    Cons, EosLaw, Euler, EulerFields, FlowBc, FlowBcs, FlowError, I_EN, I_G1, I_RC, I_RHO, NCOMP,
    TableEos,
};
use crucible_solvers::sdc::{FlowClass, Sdc};
use crucible_tables::{Pin, Table};

const FILE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../tables/chem/lox_lh2_v0.1.0.h5"
);
const PINS_TOML: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../tables/chem/lox_lh2_v0.1.0.pins.toml"
));

// The RL10-class chamber point (station-3 pinned values, same source):
// p_c = 32.75 bar, MR = 5 (Z = 1/6), h at the liquid-injection enthalpy.
const P_C: f64 = 32.75e5;
const H_INJ: f64 = -1083091.499726406;
const Z_MR5: f64 = 1.0 / 6.0;
const T_CHAMBER_CEA: f64 = 3225.4141325483242;

fn open_equilibrium() -> Table {
    let doc: toml::Table = PINS_TOML.parse().expect("pins sidecar parses");
    let entry = doc["/chem/lox_lh2/equilibrium"].as_table().expect("entry");
    let pin = Pin {
        data_version: entry["data_version"].as_str().expect("str").to_string(),
        content_digest: Some(entry["content_digest"].as_str().expect("str").to_string()),
    };
    Table::open(FILE, "/chem/lox_lh2/equilibrium", &pin).expect("production surface loads")
}

const ZERO_SRC: fn(f64, f64, f64, f64) -> Cons = |_, _, _, _| [0.0; NCOMP];

fn closed_tube_op(eos: TableEos<'_>) -> Euler<'_, TableEos<'_>> {
    Euler {
        eos,
        source: &ZERO_SRC,
        bcs: FlowBcs {
            r_inner: FlowBc::Reflecting,
            r_outer: FlowBc::Reflecting,
            z_lo: FlowBc::Reflecting,
            z_hi: FlowBc::Reflecting,
        },
        wall_normal: None,
        slip_wall_z_faces: true, // certified station behavior (slip everywhere)
    }
}

fn small_tube() -> (Grid, EulerFields) {
    let spec = GridSpec {
        r_min: 0.0,
        dr: 0.01,
        n_r: 8,
        z_min: 0.0,
        dz: 0.02,
        n_z: 24,
        n_theta_max: 1,
        axisymmetry_assertion: true,
    };
    let mut g = Grid::build(spec, crucible_solvers::euler::EULER_FIELDS).expect("valid spec");
    let f = EulerFields::resolve(&g).expect("fields");
    // Placeholder fill; tests overwrite via fill_from_prim/cons.
    for k in 0..NCOMP {
        g.fill_field(f.ids()[k], |_, _, _| 0.0);
    }
    (g, f)
}

#[test]
fn solv1_s34_projection_recovers_the_rl10_chamber() {
    let t = open_equilibrium();
    let eos = TableEos::bind(&t).expect("bind");
    let u = eos
        .cons_from_phz(P_C, H_INJ, Z_MR5, [0.0, 0.0, 0.0])
        .expect("chamber state on the surface");
    let w = eos.prim_checked(&u).expect("projection succeeds");
    assert!(
        ((w[4] - P_C) / P_C).abs() < 1e-8,
        "projected p {} vs chamber {}",
        w[4],
        P_C
    );
    // The flame is in the EOS: temperature at the projected state is the
    // CEA chamber solution (station-3 seam agreement was 0.42 K).
    let temp = eos.temperature_w(&w).expect("T at projected state");
    assert!(
        (temp - T_CHAMBER_CEA).abs() < 5.0,
        "T {temp} vs CEA {T_CHAMBER_CEA}"
    );
    // Aux-slot consistency: Γ₁ reproduces the tabulated sound speed.
    let a_slot = eos.sound_speed_w(&w);
    assert!(w[I_G1] > 1.0 && w[I_G1] < 2.0, "Γ₁ sane: {}", w[I_G1]);
    assert!(a_slot > 800.0, "chamber sound speed {a_slot} m/s");
    // Conserved round trip through the aux slots.
    let back = eos.prim_to_cons(&w);
    for k in 0..NCOMP {
        let scale = u[k].abs().max(1e-30);
        assert!(
            ((back[k] - u[k]) / scale).abs() < 1e-12,
            "component {k}: {} vs {}",
            back[k],
            u[k]
        );
    }
}

#[test]
fn solv1_s34_projection_round_trips_across_the_envelope() {
    let t = open_equilibrium();
    let eos = TableEos::bind(&t).expect("bind");
    let mut checked = 0usize;
    for p in [3.0e3, 5.0e3, 1.0e5, 1.0e6, 3.275e6, 6.5e6] {
        for h in [-11.0e6, -5.0e6, H_INJ, -1.0e6, -4.0e5] {
            for z in [0.115, Z_MR5, 0.24] {
                let Ok(u) = eos.cons_from_phz(p, h, z, [120.0, -40.0, 900.0]) else {
                    continue; // corner outside the (p,h) envelope rectangle
                };
                let w = eos
                    .prim_checked(&u)
                    .unwrap_or_else(|e| panic!("projection failed at (p={p}, h={h}, z={z}): {e}"));
                assert!(
                    ((w[4] - p) / p).abs() < 1e-8,
                    "(p={p}, h={h}, z={z}): projected {}",
                    w[4]
                );
                // Velocity decomposition survived the round trip.
                assert!((w[1] - 120.0).abs() < 1e-9 && (w[3] - 900.0).abs() < 1e-9);
                checked += 1;
            }
        }
    }
    assert!(checked > 30, "sweep coverage: {checked}");
}

#[test]
fn solv1_s34_refusals_are_loud_and_diagnosed() {
    let t = open_equilibrium();
    let eos = TableEos::bind(&t).expect("bind");
    let u_ok = eos
        .cons_from_phz(P_C, H_INJ, Z_MR5, [0.0, 0.0, 0.0])
        .expect("chamber");

    // Z outside the declared envelope refuses at the projection gate.
    let mut u = u_ok;
    u[I_RC] = u[I_RHO] * 0.5; // Z = 0.5 ≫ envelope max 0.25
    match eos.prim_checked(&u) {
        Err(what) => assert!(what.contains("mixture fraction"), "{what}"),
        Ok(_) => panic!("off-envelope Z must refuse"),
    }

    // An energy far above every tabulated enthalpy has no admissible
    // pressure bracket — refused, never clamped.
    let mut u = u_ok;
    u[I_EN] = u[I_RHO] * 5.0e7;
    match eos.prim_checked(&u) {
        Err(what) => assert!(
            what.contains("bracket") || what.contains("envelope"),
            "{what}"
        ),
        Ok(_) => panic!("off-surface state must refuse"),
    }

    // StagnationInflow has no closed form on a tabulated surface.
    match eos.stagnation_ghost(1.0, 1.0, 0.0, 0.0, 3) {
        Err(FlowError::BcUnsupportedByEos { bc }) => assert_eq!(bc, "StagnationInflow"),
        other => panic!("expected BcUnsupportedByEos, got {other:?}"),
    }
}

#[test]
fn solv1_s34_uniform_equilibrium_rest_is_a_bitwise_fixed_point() {
    let t = open_equilibrium();
    let eos = TableEos::bind(&t).expect("bind");
    let u0 = eos
        .cons_from_phz(P_C, H_INJ, Z_MR5, [0.0, 0.0, 0.0])
        .expect("chamber state");
    let (mut g, f) = small_tube();
    for (k, &v) in u0.iter().enumerate() {
        g.fill_field(f.ids()[k], move |_, _, _| v);
    }
    let snapshot = |g: &Grid| -> Vec<u64> {
        let mut bits = Vec::new();
        for b in g.bricks() {
            for &id in &f.ids() {
                bits.extend(b.field(id).iter().map(|v| v.to_bits()));
            }
        }
        bits
    };
    let before = snapshot(&g);
    let op = closed_tube_op(eos);
    let mut sdc = Sdc::new();
    let flow = FlowClass {
        op: &op,
        fields: &f,
    };
    let dt = sdc.stable_dt(&g, &flow, 0.4).expect("dt");
    let mut t_now = 0.0;
    for _ in 0..5 {
        sdc.step_flow(&mut g, &flow, t_now, dt).expect("march");
        t_now += dt;
    }
    assert_eq!(
        snapshot(&g),
        before,
        "uniform equilibrium gas at rest must be a bitwise fixed point of \
         the full step under the table EOS (well-balance is geometric, not \
         EOS-specific)"
    );
}

#[test]
fn solv1_s34_closed_hot_cold_tube_conserves_through_real_dynamics() {
    let t = open_equilibrium();
    let eos = TableEos::bind(&t).expect("bind");
    // Hot chamber-like gas left, colder same-pressure gas right: a contact
    // + acoustic transient inside closed walls — every wave runs through
    // the projection every stage.
    let hot = eos
        .cons_from_phz(P_C, H_INJ, Z_MR5, [0.0, 0.0, 0.0])
        .expect("hot");
    let cold = eos
        .cons_from_phz(P_C, -4.0e6, 0.14, [0.0, 0.0, 0.0])
        .expect("cold");
    let (mut g, f) = small_tube();
    let mid = 0.02 * 12.0;
    for (k, (&h_k, &c_k)) in hot.iter().zip(cold.iter()).enumerate() {
        g.fill_field(f.ids()[k], move |_, _, z| if z < mid { h_k } else { c_k });
    }
    let totals = |g: &Grid| -> [f64; 3] {
        let mut acc = [0.0f64; 3];
        g.for_each_active_cell(|c| {
            let b = g.brick(c.bi);
            let vol = g.cell_volume(c.i_r, 1);
            acc[0] += b.field(f.ids()[I_RHO])[c.idx] * vol;
            acc[1] += b.field(f.ids()[I_EN])[c.idx] * vol;
            acc[2] += b.field(f.ids()[I_RC])[c.idx] * vol;
        });
        acc
    };
    let before = totals(&g);
    let op = closed_tube_op(eos);
    let mut sdc = Sdc::new();
    let flow = FlowClass {
        op: &op,
        fields: &f,
    };
    let mut t_now = 0.0;
    for _ in 0..40 {
        let dt = sdc.stable_dt(&g, &flow, 0.4).expect("dt");
        sdc.step_flow(&mut g, &flow, t_now, dt).expect("step");
        t_now += dt;
    }
    let after = totals(&g);
    for (q, (b, a)) in before.iter().zip(after.iter()).enumerate() {
        let rel = ((a - b) / b).abs();
        assert!(
            rel < 1e-12,
            "quantity {q} drifted {rel:.2e} in a closed tube (flux-form \
             telescoping must hold under the table EOS)"
        );
    }
    // And the dynamics were real: the axial-momentum field is no longer
    // identically zero (waves ran; this is not a trivial fixed point).
    let mut max_mom = 0.0f64;
    g.for_each_active_cell(|c| {
        max_mom = max_mom.max(g.brick(c.bi).field(f.ids()[3])[c.idx].abs());
    });
    assert!(max_mom > 0.0, "the hot/cold transient must actually move");
}

#[test]
fn solv1_s34_projection_is_deterministic() {
    let t = open_equilibrium();
    let eos = TableEos::bind(&t).expect("bind");
    let u = eos
        .cons_from_phz(1.7e6, -3.3e6, 0.19, [55.0, 5.0, -220.0])
        .expect("state");
    let w1 = eos.prim_checked(&u).expect("first");
    let w2 = eos.prim_checked(&u).expect("second");
    for k in 0..w1.len() {
        assert_eq!(w1[k].to_bits(), w2[k].to_bits(), "slot {k}");
    }
}

#[test]
fn solv1_s18_knockdown_actually_bites() {
    // The S18 η_c* hook must CHANGE the realized equilibrium (session-12
    // regression: the session-11 wiring compensated the offset on both the
    // store and query sides — a pure gauge relabeling that left every
    // readout bit-identical; found when the calibration trial returned the
    // baseline to the last digit). Correct semantics: stored energy true,
    // interrogation shifted — so a knocked chamber state must sit COLDER
    // than full equilibrium at the same (p, h_inj, Z), the stored energy
    // must NOT carry the shift, and the projection must still round-trip.
    let t = open_equilibrium();
    let eos0 = TableEos::bind(&t).expect("bind");
    let mut eosk = TableEos::bind(&t).expect("bind");
    eosk.h_offset = -3.0e5;

    let u0 = eos0
        .cons_from_phz(P_C, H_INJ, Z_MR5, [0.0, 0.0, 0.0])
        .expect("full-equilibrium chamber state");
    let uk = eosk
        .cons_from_phz(P_C, H_INJ, Z_MR5, [0.0, 0.0, 0.0])
        .expect("knocked chamber state");

    // Stored specific energy is TRUE in both (deficit sequestered, not
    // subtracted): e = h_inj − p/ρ at each realization's own density.
    let e0 = u0[I_EN] / u0[I_RHO];
    let ek = uk[I_EN] / uk[I_RHO];
    assert!(
        ((e0 - (H_INJ - P_C / u0[I_RHO])) / e0).abs() < 1e-12,
        "full-eq stored energy is h − p/ρ"
    );
    assert!(
        ((ek - (H_INJ - P_C / uk[I_RHO])) / ek).abs() < 1e-12,
        "knocked stored energy is h − p/ρ (no gauge shift)"
    );

    let w0 = eos0.prim_checked(&u0).expect("full-eq round trip");
    let wk = eosk.prim_checked(&uk).expect("knocked round trip");
    assert!(
        ((w0[4] - P_C) / P_C).abs() < 1e-6 && ((wk[4] - P_C) / P_C).abs() < 1e-6,
        "projection round-trips both states to p_c"
    );
    let t0 = eos0.temperature_w(&w0).expect("T full");
    let tk = eosk.temperature_w(&wk).expect("T knocked");
    assert!(
        (t0 - T_CHAMBER_CEA).abs() < 2.0,
        "full equilibrium reproduces the pinned CEA chamber T: {t0}"
    );
    // δh = −3e5 J/kg at cp_eq ~ 8–9 kJ/kg/K near 3200 K ⇒ ~30–40 K colder.
    assert!(
        tk < t0 - 15.0,
        "the knockdown must bite: T_knocked {tk} vs T_full {t0}"
    );
}
