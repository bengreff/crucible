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
//! Honest scaffolding (the session-5/7 pattern, superseded by COUP-3):
//! the coupling is explicit flux-matched operator splitting at the gas CFL
//! dt, not the Robin-Robin Picard sweeps inside a class-`D` implicit solve
//! — valid here because the gas dt is far below every thermal stability
//! limit (guarded each step, fail-loud). The liner's ρc_p is set small as
//! a **steady-state continuation device** (declared): the steady conjugate
//! solution is independent of ρc_p, and the certificate's claims are
//! steady-state claims; transient wall fidelity arrives with class-`D`.
//!
//! The oracle: at steady state the coupled system must reproduce the
//! cylindrical series-resistance solution — the same film + ln-annulus +
//! coolant-film network every heat-transfer text builds (VAL-1 rung i,
//! analytic) — using only the simulated *gas* state and the declared
//! coolant data, never the simulated solid temperatures. Axial conduction
//! and discreteness set the tolerance, measured and pinned below.

use std::collections::BTreeMap;

use crate::conduction::{Bcs, Conduction, Domain, FaceBc, InteriorFaces, SolverError};
use crate::euler::{
    Cons, Euler, EulerFields, FlowBc, FlowBcs, FlowError, GammaLaw, I_EN, I_MZ, I_RHO, NCOMP, Prim,
    fill_from_prim, prim6,
};
use crate::station1_sod::FIELDS as EULER_FIELD_NAMES;
use crate::wall_heat::{NearWallGas, WallHeatError, WallLaw};
use crucible_grid::{BRICK, FaceDir, FieldId, Grid, GridSpec, InterfaceFace, Region};
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
/// Fractions of the two thermal stability limits dt may use (guarded).
pub const SOLID_DT_FRAC: f64 = 0.5;
pub const EXCHANGE_DT_FRAC: f64 = 0.2;
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
    Flow(FlowError),
    Solid(SolverError),
    Wall(WallHeatError),
    /// The explicit-coupling premise (gas CFL dt below every thermal
    /// stability limit) failed — a structured refusal (META-2 §4 halt
    /// condition), because running degraded would silently violate the
    /// scaffolding's declared validity; the class-D implicit solve
    /// (COUP-3) is the correct tool for such parameters.
    ThermalLimitUnderCfl {
        dt_gas: f64,
        dt_solid: f64,
        dt_exchange: f64,
    },
    /// Step-count backstop tripped (see [`MARCH_STEP_CAP`]).
    RunawayMarch {
        steps: usize,
    },
}

impl std::fmt::Display for CoupledError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Flow(e) => write!(f, "gas operator: {e}"),
            Self::Solid(e) => write!(f, "solid operator: {e}"),
            Self::Wall(e) => write!(f, "wall law: {e}"),
            Self::ThermalLimitUnderCfl {
                dt_gas,
                dt_solid,
                dt_exchange,
            } => write!(
                f,
                "explicit-coupling premise violated: a thermal stability limit (solid \
                 {dt_solid:.3e} s, exchange {dt_exchange:.3e} s) undercuts the gas CFL dt \
                 ({dt_gas:.3e} s) — refusing to run degraded; COUP-3's class-D implicit \
                 solve is required for these parameters"
            ),
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

impl From<FlowError> for CoupledError {
    fn from(e: FlowError) -> Self {
        Self::Flow(e)
    }
}
impl From<SolverError> for CoupledError {
    fn from(e: SolverError) -> Self {
        Self::Solid(e)
    }
}
impl From<WallHeatError> for CoupledError {
    fn from(e: WallHeatError) -> Self {
        Self::Wall(e)
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

fn cell_add(g: &mut Grid, f: FieldId, i_r: usize, i_z: usize, j: u32, dv: f64) {
    let bi = g
        .brick_index_by_coords((i_r / BRICK) as u32, (i_z / BRICK) as u32)
        .expect("cell in an allocated brick");
    let idx = {
        let b = g.brick(bi);
        b.cell_index(j, (i_r % BRICK) * BRICK + i_z % BRICK)
    };
    g.brick_field_mut(bi, f)[idx] += dv;
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

/// The assembled fixture: grid + fields + the config-time wall-face list.
pub struct Duct {
    pub grid: Grid,
    pub flow: EulerFields,
    pub t_solid: FieldId,
    pub rate_solid: FieldId,
    pub faces: Vec<InterfaceFace>,
    pub law: WallLaw,
    pub eos: GammaLaw,
}

/// Ledger + last-step exchange record for diagnostics and audits.
#[derive(Debug, Clone, Default)]
pub struct ExchangeRecord {
    /// Per-face flux q [W/m²] of the LAST step, face order = `Duct::faces`.
    pub q: Vec<f64>,
    /// Per-face film h and recovery T_aw of the last step.
    pub h: Vec<f64>,
    pub t_aw: Vec<f64>,
    /// ∫ Σ q·A dt — total heat leaving the gas through the wall [J].
    pub wall_joules: f64,
    /// ∫ coolant extraction dt [J] (negative of heat into coolant would be
    /// symmetric; recorded positive out of the solid).
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
    let law = WallLaw::new(
        specific_heat_capacity_j_per_kg_k(CP),
        dynamic_viscosity_pa_s(MU),
        PR,
    )
    .expect("valid transport set");
    Duct {
        grid: g,
        flow,
        t_solid,
        rate_solid,
        faces,
        law,
        eos,
    }
}

fn duct_flow_op(
    eos: GammaLaw,
    inflow: &'static dyn Fn(f64, f64, f64, f64) -> Prim,
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
    }
}

fn dir_code(d: FaceDir) -> u8 {
    match d {
        FaceDir::RMinus => 0,
        FaceDir::RPlus => 1,
        FaceDir::ZMinus => 2,
        FaceDir::ZPlus => 3,
    }
}

fn opposite(d: FaceDir) -> FaceDir {
    match d {
        FaceDir::RMinus => FaceDir::RPlus,
        FaceDir::RPlus => FaceDir::RMinus,
        FaceDir::ZMinus => FaceDir::ZPlus,
        FaceDir::ZPlus => FaceDir::ZMinus,
    }
}

/// Per-face wall-normal geometry: (u_t from the prim, wall distance y,
/// solid center-to-face distance d_s), sized by the grid's own spacings —
/// never fixture constants (review finding).
fn face_geometry(w: &Prim, dir: FaceDir, dr: f64, dz: f64) -> (f64, f64, f64) {
    match dir {
        FaceDir::RMinus | FaceDir::RPlus => {
            let u_t = (w[2] * w[2] + w[3] * w[3]).sqrt();
            (u_t, 0.5 * dr, 0.5 * dr)
        }
        FaceDir::ZMinus | FaceDir::ZPlus => {
            let u_t = (w[1] * w[1] + w[2] * w[2]).sqrt();
            (u_t, 0.5 * dz, 0.5 * dz)
        }
    }
}

/// One coupled step at the current state: evaluate every wall exchange
/// once, advance the solid (with the exchanges as gas-face fluxes and the
/// coolant Robin at the outer edge), debit the gas energy, advance the gas.
/// Returns the dt taken.
pub fn coupled_step(
    duct: &mut Duct,
    op: &Euler<'_>,
    t: f64,
    dt_cap: f64,
    rec: &mut ExchangeRecord,
) -> Result<f64, CoupledError> {
    let nt = 1u32; // axisymmetric fixture (asserted by the grid build)

    // --- Wall exchanges from the pre-step state (single evaluation) ------
    let ids = duct.flow.ids();
    let mut q_map: BTreeMap<(usize, usize, u8), f64> = BTreeMap::new();
    rec.q.clear();
    rec.h.clear();
    rec.t_aw.clear();
    let mut u_series_max = 0.0f64;
    for face in &duct.faces {
        let mut u = [0.0f64; NCOMP];
        for (k, id) in ids.iter().enumerate() {
            u[k] = cell_value(&duct.grid, *id, face.gas.0, face.gas.1, 0);
        }
        let w = duct
            .eos
            .prim_checked(&u)
            .map_err(|what| FlowError::NonPhysicalState {
                i_r: face.gas.0,
                i_z: face.gas.1,
                i_theta: 0,
                what,
            })?;
        // Wall-normal geometry from the grid's own spec (review finding:
        // fixture constants here would silently mis-size Re_y and the
        // series resistance on any differently-spaced grid).
        let (dr, dz) = (duct.grid.spec().dr, duct.grid.spec().dz);
        let (u_t, y, d_s) = face_geometry(&w, face.dir, dr, dz);
        let temperature = w[4] / (w[0] * R_SPECIFIC);
        let gas = NearWallGas {
            rho: w[0],
            u_t,
            temperature,
            y,
        };
        let t_s = cell_value(&duct.grid, duct.t_solid, face.solid.0, face.solid.1, 0);
        let ex = duct.law.wall_exchange(&gas, t_s, d_s, KAPPA_S)?;
        rec.q.push(ex.q);
        rec.h.push(ex.h);
        rec.t_aw.push(ex.t_aw);
        u_series_max = u_series_max.max(ex.u_series);
        q_map.insert(
            (face.solid.0, face.solid.1, dir_code(opposite(face.dir))),
            ex.q,
        );
    }

    // --- dt: gas CFL, capped by both thermal stability limits ------------
    let zero_src = |_: f64, _: f64, _: f64, _: f64| 0.0;
    let gas_face_q = |i_r: usize, i_z: usize, _j: u32, dir: FaceDir| -> f64 {
        *q_map
            .get(&(i_r, i_z, dir_code(dir)))
            .expect("every gas-facing solid face has an exchange")
    };
    let solid_op = Conduction {
        kappa: KAPPA_S,
        rho_cp: RHO_CP_S,
        source: &zero_src,
        domain: Domain::Solid,
        interior: InteriorFaces {
            gas: Some(&gas_face_q),
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
    let dt_gas = op.stable_dt(&duct.grid, &duct.flow, CFL_S4)?;
    let dt_solid = SOLID_DT_FRAC * solid_op.stable_dt(&duct.grid, 1.0);
    let spec_dr = duct.grid.spec().dr;
    let dt_exchange = if u_series_max > 0.0 {
        EXCHANGE_DT_FRAC * RHO_CP_S * spec_dr.min(duct.grid.spec().dz) / u_series_max
    } else {
        f64::INFINITY
    };
    let dt = dt_gas.min(dt_solid).min(dt_exchange).min(dt_cap);
    if dt != dt_gas.min(dt_cap) {
        return Err(CoupledError::ThermalLimitUnderCfl {
            dt_gas,
            dt_solid,
            dt_exchange,
        });
    }

    // --- Coolant ledger BEFORE the solid advance: the Robin faces inside
    // the step extract heat from the PRE-step temperatures (rate-then-
    // apply), so the ledger must read the same state (review finding: a
    // post-step read biased the transient closure).
    let mut coolant_watts = 0.0f64;
    let outer = N_R_GAS + N_R_SOLID - 1;
    let a_outer = duct.grid.face_area_r(outer, true, nt);
    for i_z in 0..N_Z {
        let t_s = cell_value(&duct.grid, duct.t_solid, outer, i_z, 0);
        coolant_watts += a_outer * (t_s - T_COOL) / (1.0 / H_COOL + 0.5 * spec_dr / KAPPA_S);
    }
    rec.coolant_joules += coolant_watts * dt;

    // --- Solid advance (reads the exchanges through the interface seam) --
    solid_op.step(&mut duct.grid, duct.t_solid, duct.rate_solid, t, dt)?;

    // --- Gas energy debit: the SAME q, opposite sign, per face -----------
    let mut wall_watts = 0.0f64;
    for (face, q) in duct.faces.iter().zip(&rec.q) {
        let area = duct.grid.interface_area_per_theta(face, nt);
        let vol = duct.grid.cell_volume(face.gas.0, nt);
        cell_add(
            &mut duct.grid,
            ids[I_EN],
            face.gas.0,
            face.gas.1,
            0,
            -q * area * dt / vol,
        );
        wall_watts += q * area;
    }
    rec.wall_joules += wall_watts * dt;

    // --- Gas advance ------------------------------------------------------
    op.step(&mut duct.grid, &duct.flow, t, dt)?;
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
    let mut t = t0;
    let mut steps = 0usize;
    while t_final - t > 1e-12 * t_final {
        let dt = coupled_step(duct, op, t, t_final - t, rec)?;
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
    let law = WallLaw::new(
        specific_heat_capacity_j_per_kg_k(CP),
        dynamic_viscosity_pa_s(MU),
        PR,
    )
    .expect("valid transport set");
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
    };
    (
        Duct {
            grid: g,
            flow,
            t_solid,
            rate_solid,
            faces,
            law,
            eos,
        },
        op,
    )
}

/// One coupled step of the stepped cavity with INSULATED exterior/edges
/// (the closed-ledger variant of `coupled_step` — same exchange, same
/// order, HeatFlux(0) instead of the coolant Robin).
pub fn stepped_cavity_step(
    duct: &mut Duct,
    op: &Euler<'_>,
    t: f64,
    rec: &mut ExchangeRecord,
) -> Result<f64, CoupledError> {
    let ids = duct.flow.ids();
    let mut q_map: BTreeMap<(usize, usize, u8), f64> = BTreeMap::new();
    rec.q.clear();
    let mut u_series_max = 0.0f64;
    for face in &duct.faces {
        let mut u = [0.0f64; NCOMP];
        for (k, id) in ids.iter().enumerate() {
            u[k] = cell_value(&duct.grid, *id, face.gas.0, face.gas.1, 0);
        }
        let w = duct
            .eos
            .prim_checked(&u)
            .map_err(|what| FlowError::NonPhysicalState {
                i_r: face.gas.0,
                i_z: face.gas.1,
                i_theta: 0,
                what,
            })?;
        let (dr, dz) = (duct.grid.spec().dr, duct.grid.spec().dz);
        let (u_t, y, d_s) = face_geometry(&w, face.dir, dr, dz);
        let temperature = w[4] / (w[0] * R_SPECIFIC);
        let gas = NearWallGas {
            rho: w[0],
            u_t,
            temperature,
            y,
        };
        let t_s = cell_value(&duct.grid, duct.t_solid, face.solid.0, face.solid.1, 0);
        let ex = duct.law.wall_exchange(&gas, t_s, d_s, KAPPA_S)?;
        rec.q.push(ex.q);
        u_series_max = u_series_max.max(ex.u_series);
        q_map.insert(
            (face.solid.0, face.solid.1, dir_code(opposite(face.dir))),
            ex.q,
        );
    }
    let zero_src = |_: f64, _: f64, _: f64, _: f64| 0.0;
    let gas_face_q = |i_r: usize, i_z: usize, _j: u32, dir: FaceDir| -> f64 {
        *q_map
            .get(&(i_r, i_z, dir_code(dir)))
            .expect("every gas-facing solid face has an exchange")
    };
    let solid_op = Conduction {
        kappa: KAPPA_S,
        rho_cp: RHO_CP_S,
        source: &zero_src,
        domain: Domain::Solid,
        interior: InteriorFaces {
            gas: Some(&gas_face_q),
            exterior: Some(FaceBc::HeatFlux(0.0)), // closed ledger
        },
        bcs: Bcs {
            r_inner: FaceBc::HeatFlux(0.0),
            r_outer: FaceBc::HeatFlux(0.0),
            z_lo: FaceBc::HeatFlux(0.0),
            z_hi: FaceBc::HeatFlux(0.0),
        },
    };
    let dt_gas = op.stable_dt(&duct.grid, &duct.flow, CFL_S4)?;
    let dt_solid = SOLID_DT_FRAC * solid_op.stable_dt(&duct.grid, 1.0);
    let dt_exchange = if u_series_max > 0.0 {
        EXCHANGE_DT_FRAC * RHO_CP_S * duct.grid.spec().dr.min(duct.grid.spec().dz) / u_series_max
    } else {
        f64::INFINITY
    };
    let dt = dt_gas.min(dt_solid).min(dt_exchange);
    solid_op.step(&mut duct.grid, duct.t_solid, duct.rate_solid, t, dt)?;
    let mut wall_watts = 0.0f64;
    for (face, q) in duct.faces.iter().zip(&rec.q) {
        let area = duct.grid.interface_area_per_theta(face, 1);
        let vol = duct.grid.cell_volume(face.gas.0, 1);
        cell_add(
            &mut duct.grid,
            ids[I_EN],
            face.gas.0,
            face.gas.1,
            0,
            -q * area * dt / vol,
        );
        wall_watts += q * area;
    }
    rec.wall_joules += wall_watts * dt;
    op.step(&mut duct.grid, &duct.flow, t, dt)?;
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
