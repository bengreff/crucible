//! SOLV-family operators on the world-state grid: the **conduction
//! operator** (SOLV-1 §3.5's diffusion class; the Goal-A convergence
//! certificate) and the **Euler flux operator** (SOLV-1 §3.1–3.3 reacting-gas
//! subset — PPM/HLLC-Batten on the exact cylindrical metric; the Goal-B
//! Station-1 Sod certificate + the §6-2 MMS), with their COUP-8 registry
//! entries (`conduction`, `flow`).
//!
//! Time integration (S2): the production advance is **COUP-3's SDC-coupled
//! IMEX step** (`sdc` — explicit hyperbolic class `A` + the fixed-cycle
//! implicit diffusion class `D` with the Robin-Robin wall exchange inside,
//! COUP-2's conservation audit armed every step). The explicit scaffolding
//! integrators of sessions 5–12 are retired; the flux-form spatial
//! operators carried over unchanged.
//!
//! The missing forces (S3): `gas_diffusion` — SOLV-1 §3.1's `F_visc`
//! (compressible viscous stress + Fourier conduction + species diffusion
//! on the exact cylindrical metric, swirl included), the gas occupant of
//! class `D`, suppressed at wall-law faces (SOLV-1 §3.5 — the wall
//! function replaces, never adds).
//!
//! Real properties (S4): `transport` — the FND-7 §3.3 spine seam, the ONE
//! provider of the diffusive-flux closure (μ, k, c_p, c_v, ρD, ∂h/∂Z, Pr)
//! over the medium state, with the declared-constant and the OFFL-5 §3.1a
//! tabulated occupants. `gas_diffusion` and `wall_heat` both read it; after
//! S4 neither states a transport constant of its own.
//!
//! Session scope: uniform N_θ across bricks per sweep (asserted); the
//! conservative flux aggregation across an N_θ jump (AMR-refluxing style,
//! FND-2 §3.4) arrives with the plan's 3-D wave (S8).

pub mod certificate;
mod conduction;
pub mod euler;
pub mod euler_mms;
pub mod gas_diffusion;
mod mechanism;
pub mod sdc;
pub mod station1_sod;
pub mod station2_nozzle;
pub mod station4_cooled_wall;
pub mod structural_margins;
pub mod transport;
pub mod wall_heat;

pub use conduction::{
    Bcs, Conduction, Domain, ExchangeKey, FaceBc, GasFaceRobin, HeatLedger, InteriorFaces,
    SolverError,
};
pub use mechanism::{
    CONDUCTION_MANIFEST, ConductionSetup, FLOW_MANIFEST, FlowSetup, GAS_DIFFUSION_MANIFEST,
    SetupError, TRANSPORT_CONSTANT_MANIFEST, TRANSPORT_TABLE_MANIFEST, WALL_HEAT_MANIFEST,
    constant_transport_from_loaded, flow_from_loaded, from_loaded, registry,
    table_transport_schmidt_from_loaded, wall_law_from_loaded,
};
