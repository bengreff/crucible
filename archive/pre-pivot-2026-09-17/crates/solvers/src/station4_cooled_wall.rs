//! GOAL-B station 4 — the cooled wall (SOLV-1 §3.5 + COUP-2 §3.5 subset).
//!
//! The physical system: supersonic hot gas (M = 2, 800 K static — recovery
//! temperature ≈ 1373 K) flows through a straight cylindrical duct whose
//! metal liner is regeneratively cooled on the outside. The gas-side flux
//! comes from the **one wall-function law** (`wall_heat`, ±20–30% declared
//! band); the liner conducts it radially (the Goal-A-certified operator on
//! the `Solid` region); the coolant side is a Robin film at the outer
//! surface (COUP-7's cooling-jacket object supersedes that closure). One
//! wall exchange is evaluated **once per face per step** and applied with
//! opposite signs to both sides — interface conservation by construction
//! (COUP-2's ledger discipline).
//!
//! S2: the fixture marches on the ONE production integrator — the SDC-IMEX
//! step with the Robin-Robin exchange inside each sweep's class-`D` solve
//! (COUP-2 §3.5; the explicit flux-matched splitting is retired). Δt is
//! the gas CFL alone — the solid/exchange stability limits are gone by
//! construction (implicit), which is what frees S4 to give the liner its
//! physical ρc_p. Here ρc_p stays the declared steady-state continuation
//! value (the steady conjugate solution is independent of it; the physical
//! value changes only how fast the same state is reached — retired at S4
//! per the plan). The COUP-2 audit is armed every step of every march.
//!
//! The oracle: at steady state the coupled system must reproduce the
//! cylindrical series-resistance solution — the same film + ln-annulus +
//! coolant-film network every heat-transfer text builds (VAL-1 rung i,
//! analytic) — using only the simulated *gas* state and the declared
//! coolant data, never the simulated solid temperatures. Axial conduction
//! and discreteness set the tolerance, measured and pinned below.

use crate::conduction::{Bcs, Conduction, Domain, FaceBc, InteriorFaces};
use crate::euler::{
    Cons, Euler, EulerFields, FlowBc, FlowBcs, GammaLaw, I_EN, I_MZ, I_RHO, NCOMP, Prim,
    fill_from_prim, prim6,
};
use crate::sdc::{
    DiffusionClass, ExchangeClass, FlowClass, Sdc, SdcError, WallPatch, build_wall_patches,
};
use crate::station1_sod::FIELDS as EULER_FIELD_NAMES;
use crate::transport::{ConstantTransport, TransportProps};
use crate::wall_heat::WallLaw;
use crucible_grid::{BRICK, FieldId, Grid, GridSpec, InterfaceFace, Region};
use crucible_units::{dynamic_viscosity_pa_s, specific_heat_capacity_j_per_kg_k};

// --- Fixture geometry (SI) ---------------------------------------------------

/// Gas core radius [m] (20 cells) and liner thickness [m] (4 cells).
pub const R_GAS: f64 = 0.05;
pub const LINER: f64 = 0.01;
pub const DUCT_LEN: f64 = 0.4;
pub const N_R_GAS: usize = 20;
pub const N_R_SOLID: usize = 4;
pub const N_Z: usize = 160;
pub const DR: f64 = R_GAS / N_R_GAS as f64; // = DUCT_LEN / N_Z: uniform 2.5 mm
pub const DZ: f64 = DUCT_LEN / N_Z as f64;

// --- Gas (hot combustion-product stand-in; degenerate spine constants) -------

pub const GAMMA: f64 = 1.4;
pub const CP: f64 = 1004.0; // J/(kg·K)
pub const R_SPECIFIC: f64 = CP * (GAMMA - 1.0) / GAMMA; // ≈ 286.9
pub const MU: f64 = 4.0e-5; // Pa·s (hot-gas scale)
pub const PR: f64 = 0.72;
/// Station 4's Schmidt number. The fixture carries no composition gradient
/// (C ≡ 0 everywhere), so the species flux and its enthalpy limb are
/// identically zero here — the value only has to be a legal gas one.
pub const SCHMIDT: f64 = 0.7;

/// Station 4's FND-7 §3.3 spine query: the declared-constant occupant
/// (`CP`, `MU`, `PR`, `GAMMA`, `SCHMIDT`). A single-composition gamma-law
/// gas has state-independent properties, which is what makes this
/// fixture's analytic comparison meaningful — so the bundle is built and
/// validated ONCE and handed out, not reconstructed per wall-patch query
/// (S4 review: this is a coupler-rate call inside every Picard sweep).
static STATION4_PROPS: std::sync::LazyLock<TransportProps> = std::sync::LazyLock::new(|| {
    ConstantTransport::new(
        specific_heat_capacity_j_per_kg_k(CP),
        dynamic_viscosity_pa_s(MU),
        PR,
        GAMMA,
        SCHMIDT,
    )
    .expect("station-4 transport constants are in range")
    .into_props()
});

pub fn station4_transport(_w: &Prim) -> Result<TransportProps, &'static str> {
    Ok(*STATION4_PROPS)
}
pub const MACH_IN: f64 = 2.0;
pub const T_IN: f64 = 800.0; // K static ⇒ T_aw ≈ 1373 K
pub const P_IN: f64 = 2.0e5; // Pa

// --- Liner + coolant ---------------------------------------------------------

pub const KAPPA_S: f64 = 20.0; // W/(m·K), steel-class
/// Steady-state continuation ρc_p [J/(m³·K)] (module doc); the physical
/// liner value (~3.6e6) only changes how fast the same steady state is
/// reached.
pub const RHO_CP_S: f64 = 200.0;
pub const H_COOL: f64 = 2.0e4; // W/(m²·K), regen-channel scale
pub const T_COOL: f64 = 300.0; // K
pub const T_SOLID_INIT: f64 = 400.0; // K

// --- March control -----------------------------------------------------------

pub const CFL_S4: f64 = 0.4;
pub const SETTLE_TIME_S4: f64 = 5.0e-3; // s (≈ 14 transits + liner settling)
pub const STEADY_CHECK_TIME_S4: f64 = 5.0e-4;

/// Entrance cells excluded from the pointwise profile gate (inflow-BC
/// adjacency — the station-2 entrance-band pattern).
pub const ENTRANCE_BAND: usize = 8;

pub const S4_FIELD_T_SOLID: &str = "T_solid";
pub const S4_FIELD_RATE_SOLID: &str = "rate_solid";

// --- Certificate criteria (named, shared by the test battery AND the
// certificate binary — the session-6 convention; META-2 §4 no-magic-numbers).
// Measured values on the committed fixture are recorded next to each gate.

/// Steadiness: max relative change of gas ρ and solid T over the check
/// window (measured 7.3e-4).
pub const RESID_MAX: f64 = 2.0e-3;
/// Energy-ledger closure between the three independent steady rates
/// (measured 3.2e-3 gas-vs-wall — steadiness-limited; 1.4e-7 wall-vs-coolant).
pub const LEDGER_REL_TOL: f64 = 0.02;
/// Pointwise series-resistance-oracle agreement past the entrance band
/// (measured 9.4e-4; margin ×5).
pub const ORACLE_REL_TOL: f64 = 5.0e-3;
/// Near-wall recovery-temperature ratio band: the wall-adjacent cell is
/// itself cooled, reading T_aw below free-stream by O(Δr) — measured 0.94
/// at 2.5 mm; the ±20–30% closure band owns this cell-size dependence.
pub const T_AW_RATIO_BAND: (f64, f64) = (0.85, 1.0);
/// The mid-duct flux must be rocket-scale (measured 1.06 MW/m²): guards a
/// silently de-energized fixture.
pub const Q_MID_MIN: f64 = 0.5e6;
/// Stair-interface closed-ledger conservation (measured 4.2e-12 —
/// accumulation round-off).
pub const CAVITY_CONSERVATION_TOL: f64 = 1.0e-11;
/// Robin-annulus analytic anchor, worst relative error at 32 radial cells
/// (measured well inside; discretization-limited).
pub const ANNULUS_ROBIN_TOL: f64 = 2.0e-3;
/// Runaway-march backstop (a dt collapse is a config/physics problem the
/// typed error reports — far above any legitimate fixture march).
pub const MARCH_STEP_CAP: usize = 2_000_000;

#[derive(Debug)]
pub enum CoupledError {
    /// Any SDC-step failure (flow, class-D, exchange residual, audit).
    Step(SdcError),
    /// Step-count backstop tripped (see [`MARCH_STEP_CAP`]).
    RunawayMarch { steps: usize },
}

impl std::fmt::Display for CoupledError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Step(e) => write!(f, "coupled step: {e}"),
            Self::RunawayMarch { steps } => {
                write!(
                    f,
                    "coupled march exceeded MARCH_STEP_CAP ({steps} steps) — dt collapse"
                )
            }
        }
    }
}

impl std::error::Error for CoupledError {}

impl From<SdcError> for CoupledError {
    fn from(e: SdcError) -> Self {
        Self::Step(e)
    }
}

/// Direct field read at one cell (θ-plane j) — random access for the
/// face-local coupler; sweeps stay on the fast paths.
pub fn cell_value(g: &Grid, f: FieldId, i_r: usize, i_z: usize, j: u32) -> f64 {
    let bi = g
        .brick_index_by_coords((i_r / BRICK) as u32, (i_z / BRICK) as u32)
        .expect("cell in an allocated brick");
    let b = g.brick(bi);
    b.field(f)[b.cell_index(j, (i_r % BRICK) * BRICK + i_z % BRICK)]
}

/// Fill a field over SOLID cells (the gas-mask `fill_field` twin).
pub fn fill_solid(g: &mut Grid, f: FieldId, value: f64) {
    for bi in 0..g.n_bricks() {
        let (solid, nt) = {
            let b = g.brick(bi);
            (b.solid_mask(), b.n_theta())
        };
        if solid == 0 {
            continue;
        }
        // Recompute per θ-plane: mask bits are (r,z)-shaped.
        let data = g.brick_field_mut(bi, f);
        for j in 0..nt {
            for local in 0..crucible_grid::BRICK_CELLS {
                if solid & (1u64 << local) != 0 {
                    data[j as usize * crucible_grid::BRICK_CELLS + local] = value;
                }
            }
        }
    }
}

/// The assembled fixture: grid + fields + the config-time wall-face list
/// and its SOLV-1 §3.5 patch form (box world ⇒ one patch per face, same
/// order — `patches[k]` IS `faces[k]`).
pub struct Duct {
    pub grid: Grid,
    pub flow: EulerFields,
    pub t_solid: FieldId,
    pub rate_solid: FieldId,
    pub faces: Vec<InterfaceFace>,
    pub patches: Vec<WallPatch>,
    pub law: WallLaw,
    pub eos: GammaLaw,
}

/// The FND-7 spine temperature query of this fixture's degenerate
/// constant-transport gas: T = p/(ρ·R_specific).
pub fn gas_temperature(w: &Prim) -> Result<f64, &'static str> {
    Ok(w[4] / (w[0] * R_SPECIFIC))
}

/// Ledger + last-step exchange record for diagnostics and audits. The
/// joules lines are the SDC step's OWN applied-increment records (COUP-2
/// §3.1 ledger discipline — exactly what was integrated, never rate×Δt).
#[derive(Debug, Clone, Default)]
pub struct ExchangeRecord {
    /// Per-face flux q [W/m²] at the LAST step's accepted solve, face
    /// order = `Duct::faces`.
    pub q: Vec<f64>,
    /// Per-face film h and recovery T_aw of the last step.
    pub h: Vec<f64>,
    pub t_aw: Vec<f64>,
    /// Σ applied exchange energy — total heat leaving the gas through the
    /// wall [J].
    pub wall_joules: f64,
    /// Σ applied coolant extraction [J] (recorded positive out of the
    /// solid).
    pub coolant_joules: f64,
}

pub fn inflow_prim() -> Prim {
    let a = (GAMMA * R_SPECIFIC * T_IN).sqrt();
    let rho = P_IN / (R_SPECIFIC * T_IN);
    prim6(rho, 0.0, 0.0, MACH_IN * a, P_IN, 0.0)
}

pub fn build_duct() -> Duct {
    let spec = GridSpec {
        r_min: 0.0,
        dr: DR,
        n_r: N_R_GAS + N_R_SOLID,
        z_min: 0.0,
        dz: DZ,
        n_z: N_Z,
        n_theta_max: 1,
        axisymmetry_assertion: true,
    };
    let mut names: Vec<&str> = EULER_FIELD_NAMES.to_vec();
    names.push(S4_FIELD_T_SOLID);
    names.push(S4_FIELD_RATE_SOLID);
    let mut g = Grid::build_with_regions(spec, &names, |i_r, _| {
        if i_r < N_R_GAS {
            Region::Gas
        } else {
            Region::Solid
        }
    })
    .expect("valid spec");
    let flow = EulerFields::resolve(&g).expect("fields registered");
    let t_solid = g.field_id(S4_FIELD_T_SOLID).expect("registered");
    let rate_solid = g.field_id(S4_FIELD_RATE_SOLID).expect("registered");
    let eos = GammaLaw { gamma: GAMMA };
    let w_in = inflow_prim();
    fill_from_prim(&mut g, &flow, &eos, |_, _, _| w_in);
    fill_solid(&mut g, t_solid, T_SOLID_INIT);
    let faces = g.gas_solid_faces();
    let patches = build_wall_patches(&g).expect("box-world patches");
    let law = WallLaw::new();
    Duct {
        grid: g,
        flow,
        t_solid,
        rate_solid,
        faces,
        patches,
        law,
        eos,
    }
}

fn duct_flow_op(
    eos: GammaLaw,
    inflow: &'static (dyn Fn(f64, f64, f64, f64) -> Prim + Sync),
) -> Euler<'static> {
    const ZERO_SRC: fn(f64, f64, f64, f64) -> Cons = |_, _, _, _| [0.0; NCOMP];
    Euler {
        eos,
        source: &ZERO_SRC,
        bcs: FlowBcs {
            r_inner: FlowBc::Reflecting, // zero-area axis face
            r_outer: FlowBc::Reflecting, // unreached: liner rings interpose
            z_lo: FlowBc::Prescribed(inflow),
            z_hi: FlowBc::Transmissive, // supersonic exit
        },
        wall_normal: None,       // straight duct: grid-aligned mirror is exact
        slip_wall_z_faces: true, // certified station behavior (slip everywhere)
        combustion: None,
    }
}

/// One coupled step at the current state: evaluate every wall exchange
/// once, advance the solid (with the exchanges as gas-face fluxes and the
/// coolant Robin at the outer edge), debit the gas energy, advance the gas.
/// Returns the dt taken.
/// One coupled SDC-IMEX step (COUP-3 §3.1) at the current state: class A =
/// the gas operator, class D = the liner conduction with the coolant Robin
/// on its exterior faces, and the Robin-Robin wall exchange inside the
/// class-D solve (COUP-2 §3.5). Δt = the gas CFL alone (dt_cap-clipped);
/// the audit is armed. Returns the dt taken and updates the record from
/// the step's own applied-increment ledger.
pub fn coupled_step(
    duct: &mut Duct,
    op: &Euler<'_>,
    sdc: &mut Sdc,
    t: f64,
    dt_cap: f64,
    rec: &mut ExchangeRecord,
) -> Result<f64, CoupledError> {
    let zero_src = |_: f64, _: f64, _: f64, _: f64| 0.0;
    let solid_op = Conduction {
        kappa: KAPPA_S,
        rho_cp: RHO_CP_S,
        source: &zero_src,
        domain: Domain::Solid,
        interior: InteriorFaces {
            exterior: Some(FaceBc::Robin {
                h: H_COOL,
                t_inf: T_COOL,
            }),
        },
        bcs: Bcs {
            r_inner: FaceBc::HeatFlux(0.0), // unreached: gas rings interpose
            r_outer: FaceBc::Robin {
                h: H_COOL,
                t_inf: T_COOL,
            },
            z_lo: FaceBc::HeatFlux(0.0),
            z_hi: FaceBc::HeatFlux(0.0),
        },
    };
    let flow = FlowClass {
        op,
        fields: &duct.flow,
    };
    let diffusion = DiffusionClass {
        op: &solid_op,
        t_field: duct.t_solid,
        scratch_field: duct.rate_solid,
    };
    let exchange = ExchangeClass {
        patches: &duct.patches,
        law: &duct.law,
        temperature: &gas_temperature,
        transport: &station4_transport,
    };
    let dt = sdc.stable_dt(&duct.grid, &flow, CFL_S4)?.min(dt_cap);
    let report = sdc.step(
        &mut duct.grid,
        Some(&flow),
        Some(&diffusion),
        None,
        Some(&exchange),
        None,
        t,
        dt,
    )?;
    let ex = report.exchange.expect("exchange scheduled");
    rec.q = duct
        .patches
        .iter()
        .zip(&ex.q_w)
        .map(|(p, q)| q / p.area)
        .collect();
    rec.h = ex.h;
    rec.t_aw = ex.t_aw;
    rec.wall_joules += ex.applied_exchange_j;
    // Coolant extraction = every non-exchange heat OUT of the liner: this
    // duct's liner reaches the domain edge, so its coolant Robin fires on
    // the r_outer BC line (applied_bc_j), not the exterior line — read
    // BOTH (S2 review finding: the exterior line alone records zero here;
    // the z-edge/r-inner HeatFlux(0) faces contribute exact zeros).
    rec.coolant_joules += -(ex.applied_exterior_j + ex.applied_bc_j);
    Ok(dt)
}

/// March the coupled system to `t_final`. Returns (steps, final time).
pub fn march_coupled(
    duct: &mut Duct,
    op: &Euler<'_>,
    t0: f64,
    t_final: f64,
    rec: &mut ExchangeRecord,
) -> Result<(usize, f64), CoupledError> {
    let mut sdc = Sdc::new();
    let mut t = t0;
    let mut steps = 0usize;
    while t_final - t > 1e-12 * t_final {
        let dt = coupled_step(duct, op, &mut sdc, t, t_final - t, rec)?;
        t += dt;
        steps += 1;
        if steps >= MARCH_STEP_CAP {
            return Err(CoupledError::RunawayMarch { steps });
        }
    }
    Ok((steps, t))
}

/// Run the duct to steadiness; returns (duct, record, steps, residual)
/// where the residual is the station-2-style max relative change of gas
/// density AND solid temperature over the check window.
pub fn run_duct() -> Result<(Duct, ExchangeRecord, usize, f64), CoupledError> {
    run_duct_for(SETTLE_TIME_S4, STEADY_CHECK_TIME_S4)
}

/// The duct march with explicit settle/check windows (the determinism
/// rerun uses a short window — bit-identity does not need steadiness).
pub fn run_duct_for(
    settle: f64,
    check: f64,
) -> Result<(Duct, ExchangeRecord, usize, f64), CoupledError> {
    static INFLOW: fn(f64, f64, f64, f64) -> Prim = |_, _, _, _| inflow_prim();
    let mut duct = build_duct();
    let op = duct_flow_op(duct.eos, &INFLOW);
    let mut rec = ExchangeRecord::default();
    let (steps, t) = march_coupled(&mut duct, &op, 0.0, settle, &mut rec)?;

    let rho_id = duct.flow.ids()[I_RHO];
    let before: Vec<(Vec<f64>, Vec<f64>)> = duct
        .grid
        .bricks()
        .iter()
        .map(|b| (b.field(rho_id).to_vec(), b.field(duct.t_solid).to_vec()))
        .collect();
    let (steps2, _) = march_coupled(&mut duct, &op, t, settle + check, &mut rec)?;
    let mut resid = 0.0f64;
    for (bi, b) in duct.grid.bricks().iter().enumerate() {
        let solid = b.solid_mask();
        let gas = b.mask();
        for j in 0..b.n_theta() {
            for local in 0..crucible_grid::BRICK_CELLS {
                let idx = j as usize * crucible_grid::BRICK_CELLS + local;
                if gas & (1u64 << local) != 0 {
                    let (a, r) = (before[bi].0[idx], b.field(rho_id)[idx]);
                    resid = resid.max(((r - a) / a).abs());
                }
                if solid & (1u64 << local) != 0 {
                    let (a, s) = (before[bi].1[idx], b.field(duct.t_solid)[idx]);
                    resid = resid.max(((s - a) / a).abs());
                }
            }
        }
    }
    Ok((duct, rec, steps + steps2, resid))
}

// --- The stepped cavity: stair-interface machinery under a closed ledger -----

/// Gas radius (in cells) of the stepped cavity at axial index `i_z`: three
/// steps down along z, so the liner band presents r-faces AND z-faces
/// (the stair-interface path the straight duct cannot exercise), plus
/// solid↔exterior faces beyond the band.
pub fn stepped_gas_rings(i_z: usize) -> usize {
    match i_z {
        0..=15 => 18,
        16..=31 => 14,
        _ => 10,
    }
}

pub const STEPPED_N_Z: usize = 48;
pub const STEPPED_N_R: usize = 24;
pub const STEPPED_LINER_CELLS: usize = 3;
pub const T_HOT_CAVITY: f64 = 800.0;

/// Build the closed stepped cavity: hot gas at rest against a cold stair
/// liner, exterior beyond, every domain edge reflecting/insulated. The
/// exterior faces are INSULATED here so the two-domain energy ledger is
/// closed: whatever the gas loses the liner must hold, to round-off.
pub fn build_stepped_cavity() -> (Duct, Euler<'static>) {
    let spec = GridSpec {
        r_min: 0.0,
        dr: DR,
        n_r: STEPPED_N_R,
        z_min: 0.0,
        dz: DZ,
        n_z: STEPPED_N_Z,
        n_theta_max: 1,
        axisymmetry_assertion: true,
    };
    let mut names: Vec<&str> = EULER_FIELD_NAMES.to_vec();
    names.push(S4_FIELD_T_SOLID);
    names.push(S4_FIELD_RATE_SOLID);
    let mut g = Grid::build_with_regions(spec, &names, |i_r, i_z| {
        let gas_rings = stepped_gas_rings(i_z);
        if i_r < gas_rings {
            Region::Gas
        } else if i_r < gas_rings + STEPPED_LINER_CELLS {
            Region::Solid
        } else {
            Region::Exterior
        }
    })
    .expect("valid spec");
    let flow = EulerFields::resolve(&g).expect("fields registered");
    let t_solid = g.field_id(S4_FIELD_T_SOLID).expect("registered");
    let rate_solid = g.field_id(S4_FIELD_RATE_SOLID).expect("registered");
    let eos = GammaLaw { gamma: GAMMA };
    let rho_hot = P_IN / (R_SPECIFIC * T_HOT_CAVITY);
    fill_from_prim(&mut g, &flow, &eos, |_, _, _| {
        prim6(rho_hot, 0.0, 0.0, 0.0, P_IN, 0.0)
    });
    fill_solid(&mut g, t_solid, T_SOLID_INIT);
    let faces = g.gas_solid_faces();
    let law = WallLaw::new();
    const ZERO_SRC: fn(f64, f64, f64, f64) -> Cons = |_, _, _, _| [0.0; NCOMP];
    let op = Euler {
        eos,
        source: &ZERO_SRC,
        bcs: FlowBcs {
            r_inner: FlowBc::Reflecting,
            r_outer: FlowBc::Reflecting,
            z_lo: FlowBc::Reflecting,
            z_hi: FlowBc::Reflecting,
        },
        wall_normal: None,       // grid-aligned stair mirror (exact for this test)
        slip_wall_z_faces: true, // certified station behavior (slip everywhere)
        combustion: None,
    };
    let patches = build_wall_patches(&g).expect("box-world patches");
    (
        Duct {
            grid: g,
            flow,
            t_solid,
            rate_solid,
            faces,
            patches,
            law,
            eos,
        },
        op,
    )
}

/// One coupled SDC step of the stepped cavity with INSULATED exterior/
/// edges (the closed-ledger variant of [`coupled_step`] — same exchange,
/// same schedule, `HeatFlux(0)` instead of the coolant Robin). The COUP-2
/// energy audit inside the step IS the closed-ledger claim, asserted at
/// `TOL_AUDIT` every step; the fixture's own two-domain ledger re-checks
/// it independently over the whole march.
pub fn stepped_cavity_step(
    duct: &mut Duct,
    op: &Euler<'_>,
    sdc: &mut Sdc,
    t: f64,
    rec: &mut ExchangeRecord,
) -> Result<f64, CoupledError> {
    let zero_src = |_: f64, _: f64, _: f64, _: f64| 0.0;
    let solid_op = Conduction {
        kappa: KAPPA_S,
        rho_cp: RHO_CP_S,
        source: &zero_src,
        domain: Domain::Solid,
        interior: InteriorFaces {
            exterior: Some(FaceBc::HeatFlux(0.0)), // closed ledger
        },
        bcs: Bcs {
            r_inner: FaceBc::HeatFlux(0.0),
            r_outer: FaceBc::HeatFlux(0.0),
            z_lo: FaceBc::HeatFlux(0.0),
            z_hi: FaceBc::HeatFlux(0.0),
        },
    };
    let flow = FlowClass {
        op,
        fields: &duct.flow,
    };
    let diffusion = DiffusionClass {
        op: &solid_op,
        t_field: duct.t_solid,
        scratch_field: duct.rate_solid,
    };
    let exchange = ExchangeClass {
        patches: &duct.patches,
        law: &duct.law,
        temperature: &gas_temperature,
        transport: &station4_transport,
    };
    let dt = sdc.stable_dt(&duct.grid, &flow, CFL_S4)?;
    let report = sdc.step(
        &mut duct.grid,
        Some(&flow),
        Some(&diffusion),
        None,
        Some(&exchange),
        None,
        t,
        dt,
    )?;
    let ex = report.exchange.expect("exchange scheduled");
    rec.q = duct
        .patches
        .iter()
        .zip(&ex.q_w)
        .map(|(p, q)| q / p.area)
        .collect();
    rec.wall_joules += ex.applied_exchange_j;
    Ok(dt)
}

/// Total gas energy [J] (Σ ρE·V over gas cells) and total solid thermal
/// content ρc_p·Σ T·V — the two sides of the closed cavity's ledger.
pub fn cavity_energies(duct: &Duct) -> (f64, f64) {
    let en_id = duct.flow.ids()[I_EN];
    let (mut e_gas, mut e_solid) = (0.0f64, 0.0f64);
    for i_z in 0..duct.grid.spec().n_z {
        for i_r in 0..duct.grid.spec().n_r {
            match duct.grid.region(i_r, i_z) {
                Region::Gas => {
                    e_gas +=
                        cell_value(&duct.grid, en_id, i_r, i_z, 0) * duct.grid.cell_volume(i_r, 1);
                }
                Region::Solid => {
                    e_solid += RHO_CP_S
                        * cell_value(&duct.grid, duct.t_solid, i_r, i_z, 0)
                        * duct.grid.cell_volume(i_r, 1);
                }
                Region::Exterior => {}
            }
        }
    }
    (e_gas, e_solid)
}

// --- Steady-state diagnostics for the oracle + certificate -------------------

/// Mass-flux-weighted stagnation temperature at z-plane `i_z`:
/// `T0 = Σ ρu_z·A·(T + u²/2c_p) / (Σ ρu_z·A)`.
pub fn plane_t0(duct: &Duct, i_z: usize) -> f64 {
    let ids = duct.flow.ids();
    let (mut num, mut den) = (0.0f64, 0.0f64);
    for i_r in 0..N_R_GAS {
        let mut u = [0.0f64; NCOMP];
        for (k, id) in ids.iter().enumerate() {
            u[k] = cell_value(&duct.grid, *id, i_r, i_z, 0);
        }
        let w = duct.eos.prim_checked(&u).expect("physical state");
        let a_z = duct.grid.face_area_z(i_r, 1);
        let flux = w[0] * w[3] * a_z;
        let t = w[4] / (w[0] * R_SPECIFIC);
        let speed2 = w[1] * w[1] + w[2] * w[2] + w[3] * w[3];
        num += flux * (t + 0.5 * speed2 / CP);
        den += flux;
    }
    num / den
}

/// Mass flow through z-plane `i_z` [kg/s].
pub fn plane_mdot(duct: &Duct, i_z: usize) -> f64 {
    let ids = duct.flow.ids();
    (0..N_R_GAS)
        .map(|i_r| cell_value(&duct.grid, ids[I_MZ], i_r, i_z, 0) * duct.grid.face_area_z(i_r, 1))
        .sum()
}

/// The series-resistance oracle flux at axial index `i_z` [W/m²] on the
/// inner wall: `(T_aw − T_cool) / (1/h + R1·ln(R2/R1)/κ_s + R1/(R2·h_cool))`
/// built from the simulated gas state (h, T_aw of the wall face at this z)
/// and the declared coolant data — the simulated SOLID field never enters.
pub fn oracle_flux(rec: &ExchangeRecord, face_index: usize) -> f64 {
    let r1 = R_GAS;
    let r2 = R_GAS + LINER;
    let resistance = 1.0 / rec.h[face_index] + r1 * (r2 / r1).ln() / KAPPA_S + r1 / (r2 * H_COOL);
    (rec.t_aw[face_index] - T_COOL) / resistance
}

/// The steady-state report the test battery asserts and the certificate
/// prints — one computation, two consumers (the station-2 pattern).
#[derive(Debug, Clone)]
pub struct DuctReport {
    pub steps: usize,
    pub resid: f64,
    pub mdot: f64,
    pub t0_in: f64,
    pub t0_out: f64,
    /// The three independent steady rates [W]: gas enthalpy deficit,
    /// wall-face exchange, coolant extraction.
    pub gas_watts: f64,
    pub wall_watts: f64,
    pub coolant_watts: f64,
    /// Worst |q_sim − q_1D|/q_1D past the entrance band.
    pub oracle_worst: f64,
    /// Mid-duct film h, flux q, and T_aw / free-stream-recovery ratio.
    pub h_mid: f64,
    pub q_mid: f64,
    pub t_aw_ratio_mid: f64,
    /// Liner surface temperatures at mid-duct (inner/outer rings).
    pub t_liner_inner_mid: f64,
    pub t_liner_outer_mid: f64,
}

pub fn duct_report(duct: &Duct, rec: &ExchangeRecord, steps: usize, resid: f64) -> DuctReport {
    let mdot = plane_mdot(duct, N_Z - 1);
    let t0_in = plane_t0(duct, 0);
    let t0_out = plane_t0(duct, N_Z - 1);
    let wall_watts: f64 = duct
        .faces
        .iter()
        .zip(&rec.q)
        .map(|(f, q)| q * duct.grid.interface_area_per_theta(f, 1))
        .sum();
    let outer = N_R_GAS + N_R_SOLID - 1;
    let a_outer = duct.grid.face_area_r(outer, true, 1);
    let coolant_watts: f64 = (0..N_Z)
        .map(|i_z| {
            let t_s = cell_value(&duct.grid, duct.t_solid, outer, i_z, 0);
            a_outer * (t_s - T_COOL) / (1.0 / H_COOL + 0.5 * DR / KAPPA_S)
        })
        .sum();
    let mut oracle_worst = 0.0f64;
    for (k, face) in duct.faces.iter().enumerate() {
        if face.gas.1 < ENTRANCE_BAND {
            continue;
        }
        let q_1d = oracle_flux(rec, k);
        oracle_worst = oracle_worst.max(((rec.q[k] - q_1d) / q_1d).abs());
    }
    let k_mid = N_Z / 2;
    let t_aw_freestream = T_IN
        * (1.0
            + PR.powf(crate::wall_heat::RECOVERY_PR_EXP) * 0.5 * (GAMMA - 1.0) * MACH_IN * MACH_IN);
    DuctReport {
        steps,
        resid,
        mdot,
        t0_in,
        t0_out,
        gas_watts: mdot * CP * (t0_in - t0_out),
        wall_watts,
        coolant_watts,
        oracle_worst,
        h_mid: rec.h[k_mid],
        q_mid: rec.q[k_mid],
        t_aw_ratio_mid: rec.t_aw[k_mid] / t_aw_freestream,
        t_liner_inner_mid: cell_value(&duct.grid, duct.t_solid, N_R_GAS, k_mid, 0),
        t_liner_outer_mid: cell_value(&duct.grid, duct.t_solid, outer, k_mid, 0),
    }
}
