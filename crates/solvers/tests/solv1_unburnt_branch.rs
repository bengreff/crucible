//! SOLV-1 §3.4 / SOLV-4 §3.6 (plan S5) — the **cold/unburnt branch** on the
//! *production* artifact `tables/chem/lox_lh2_unburnt_v0.2.0.h5` (regenerable
//! via `offline/scripts/make_unburnt_tables.py`), the burn-progress `c = 0`
//! branch.
//!
//! The design claim these tests verify: the unburnt-reactant surface is on
//! the **same (p, h, Z) FND-5 schema** as the equilibrium surface, so the
//! **existing [`TableEos`] occupant binds it with no new code** (OFFL-3
//! §3.3, SOLV-1 §3.4). The frozen-mode `{ρX_k}` *field* widening rides plan
//! S5b; the cold branch needs none of it, because a two-stream reactant
//! mixture is set by Z alone.
//!
//! And the establishment-cure core: **cold reactant gas marches** — the
//! projection, the flux telescoping, and the well-balance all hold on the
//! cold branch exactly as they do on the burnt one, so a startup cell has an
//! honest home instead of being forced onto the burnt surface. (S6's c-blend
//! routes cells to this branch; S5 proves the branch is real and marchable.)

use crucible_grid::{Grid, GridSpec};
use crucible_solvers::euler::{
    Cons, EosLaw, Euler, EulerFields, FlowBc, FlowBcs, I_EN, I_G1, I_RC, I_RHO, NCOMP, TableEos,
};
use crucible_solvers::sdc::{FlowClass, Sdc};
use crucible_tables::{Pin, Table};

const FILE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../tables/chem/lox_lh2_unburnt_v0.2.0.h5"
);
const PINS_TOML: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../tables/chem/lox_lh2_unburnt_v0.2.0.pins.toml"
));
const GROUP: &str = "/chem/lox_lh2/unburnt";

// Cold unburnt states on the surface (measured from the artifact):
// Z = 1/6 (design MR 5); h picked so T lands in the cold-gas band.
const Z_MR5: f64 = 1.0 / 6.0;
const H_COLD: f64 = -4.0e5; // T ≈ 168 K
const H_WARM: f64 = 5.0e5; // T ≈ 455 K
const P_MID: f64 = 1.0e6;

fn open_unburnt() -> Table {
    let doc: toml::Table = PINS_TOML.parse().expect("pins sidecar parses");
    let entry = doc[GROUP].as_table().expect("entry");
    let pin = Pin {
        data_version: entry["data_version"].as_str().expect("str").to_string(),
        content_digest: Some(entry["content_digest"].as_str().expect("str").to_string()),
    };
    Table::open(FILE, GROUP, &pin).expect("production unburnt surface loads under its pin")
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
        slip_wall_z_faces: true,
        combustion: None,
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
    for k in 0..NCOMP {
        g.fill_field(f.ids()[k], |_, _, _| 0.0);
    }
    (g, f)
}

#[test]
fn solv1_s5_table_eos_binds_the_unburnt_surface_with_no_new_occupant() {
    // The whole runtime cost of the cold branch: bind the SAME occupant to a
    // different table. No `if(unburnt)`, no new state — this is the design.
    let t = open_unburnt();
    let eos = TableEos::bind(&t).expect("the existing TableEos binds the unburnt surface");
    let u = eos
        .cons_from_phz(P_MID, H_COLD, Z_MR5, [0.0, 0.0, 0.0])
        .expect("cold unburnt state on the surface");
    let w = eos
        .prim_checked(&u)
        .expect("projection succeeds on the cold branch");
    // The projection recovers the imposed pressure (EOS inversion holds on a
    // frozen surface exactly as on the equilibrium one).
    assert!(
        ((w[4] - P_MID) / P_MID).abs() < 1e-8,
        "projected p {}",
        w[4]
    );
    // The state is genuinely COLD reactant gas, not a burnt-chamber state:
    // the equilibrium surface would read ~3000 K here; the unburnt reads ~168 K.
    let temp = eos
        .temperature_w(&w)
        .expect("T at the projected cold state");
    assert!(
        (100.0..350.0).contains(&temp),
        "cold unburnt gas temperature {temp} K (not a burnt-chamber ~3000 K)"
    );
    // A frozen H/O gas mixture is diatomic-dominated: Γ₁ near 1.4, not the
    // ~1.14 of the dissociated burnt chamber.
    assert!(w[I_G1] > 1.3 && w[I_G1] < 1.5, "frozen Γ₁ {}", w[I_G1]);
    // Conserved round trip through the aux slots (same contract as shifting).
    let back = eos.prim_to_cons(&w);
    for k in 0..NCOMP {
        let scale = u[k].abs().max(1e-30);
        assert!(((back[k] - u[k]) / scale).abs() < 1e-12, "component {k}");
    }
}

#[test]
fn solv1_s5_uniform_cold_unburnt_rest_is_a_bitwise_fixed_point() {
    // Well-balance is geometric, not EOS-specific — so it must hold on the
    // cold branch too. A uniform cold reactant gas at rest cannot drift.
    let t = open_unburnt();
    let eos = TableEos::bind(&t).expect("bind");
    let u0 = eos
        .cons_from_phz(P_MID, H_COLD, Z_MR5, [0.0, 0.0, 0.0])
        .expect("cold state");
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
        "uniform cold unburnt gas at rest must be a bitwise fixed point"
    );
}

#[test]
fn solv1_s5_cold_unburnt_transient_marches_and_conserves() {
    // The establishment-cure core: a cold reactant gas with an internal
    // gradient runs real wave dynamics through the projection every stage,
    // WITHOUT refusing, and telescopes to round-off. This is the behavior a
    // startup cell needs and could not get on the burnt surface.
    let t = open_unburnt();
    let eos = TableEos::bind(&t).expect("bind");
    let cold = eos
        .cons_from_phz(P_MID, H_COLD, Z_MR5, [0.0, 0.0, 0.0])
        .expect("cold");
    let warm = eos
        .cons_from_phz(P_MID, H_WARM, Z_MR5, [0.0, 0.0, 0.0])
        .expect("warm");
    let (mut g, f) = small_tube();
    let mid = 0.02 * 12.0;
    for (k, (&c_k, &w_k)) in cold.iter().zip(warm.iter()).enumerate() {
        g.fill_field(f.ids()[k], move |_, _, z| if z < mid { c_k } else { w_k });
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
        sdc.step_flow(&mut g, &flow, t_now, dt)
            .expect("cold branch marches");
        t_now += dt;
    }
    let after = totals(&g);
    for (q, (b, a)) in before.iter().zip(after.iter()).enumerate() {
        let rel = ((a - b) / b).abs();
        assert!(
            rel < 1e-12,
            "quantity {q} drifted {rel:.2e} on the cold branch"
        );
    }
    // Elemental Z is untouched (source-free advection) and the dynamics were
    // real (waves ran — momentum is no longer identically zero).
    let mut max_mom = 0.0f64;
    let mut z_ok = true;
    g.for_each_active_cell(|c| {
        let b = g.brick(c.bi);
        max_mom = max_mom.max(b.field(f.ids()[3])[c.idx].abs());
        let z = b.field(f.ids()[I_RC])[c.idx] / b.field(f.ids()[I_RHO])[c.idx];
        if (z - Z_MR5).abs() > 1e-9 {
            z_ok = false;
        }
    });
    assert!(max_mom > 0.0, "the cold transient must actually move");
    assert!(
        z_ok,
        "elemental Z is source-free and must stay at the injected value"
    );
}

#[test]
fn solv1_s5_below_the_gas_floor_refuses_rather_than_extrapolates() {
    // The gas-phase branch is a declared model down to its cold floor
    // (~100 K); a state below the envelope refuses (META-1 P6) — the
    // liquid/vapor regime is plan S15, not a silent extrapolation here.
    let t = open_unburnt();
    let eos = TableEos::bind(&t).expect("bind");
    // An enthalpy far below the surface floor has no admissible bracket.
    let u_ok = eos
        .cons_from_phz(P_MID, H_COLD, Z_MR5, [0.0, 0.0, 0.0])
        .expect("cold ok");
    let mut u = u_ok;
    u[I_EN] = u[I_RHO] * (-5.0e6); // e ≈ −5 MJ/kg, far below the cold floor
    match eos.prim_checked(&u) {
        Err(what) => assert!(
            what.contains("bracket") || what.contains("envelope"),
            "cold-floor refusal is diagnosed: {what}"
        ),
        Ok(_) => panic!("a sub-floor state must refuse, never extrapolate"),
    }
}
