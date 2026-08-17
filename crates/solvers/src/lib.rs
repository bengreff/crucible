//! SOLV-family operators on the world-state grid: the **conduction
//! operator** (SOLV-1 §3.5's diffusion class; the Goal-A convergence
//! certificate) and the **Euler flux operator** (SOLV-1 §3.1–3.3 reacting-gas
//! subset — PPM/HLLC-Batten on the exact cylindrical metric; the Goal-B
//! Station-1 Sod certificate + the §6-2 MMS), with their COUP-8 registry
//! entries (`conduction`, `flow`).
//!
//! Time integration note (honest scaffolding): the production advance is
//! COUP-3's SDC-coupled IMEX (explicit hyperbolic + fixed-cycle implicit
//! diffusion class `D`). This crate's explicit fixed-order reference
//! integrators exist to drive the certificates' *spatial-operator*
//! verification (MMS + analytic anchors) and are superseded — not extended —
//! when COUP-3 lands. The flux-form spatial operators carry over unchanged.
//!
//! Session scope: uniform N_θ across bricks per sweep (asserted); the
//! conservative flux aggregation across an N_θ jump (AMR-refluxing style,
//! FND-2 §3.4) arrives with the COUP-2/COUP-3 wave.

pub mod certificate;
mod conduction;
pub mod euler;
pub mod euler_mms;
mod mechanism;
pub mod station1_sod;
pub mod station2_nozzle;
pub mod station4_cooled_wall;
pub mod wall_heat;

pub use conduction::{Bcs, Conduction, Domain, FaceBc, InteriorFaces, SolverError};
pub use mechanism::{
    CONDUCTION_MANIFEST, ConductionSetup, FLOW_MANIFEST, FlowSetup, SetupError, WALL_HEAT_MANIFEST,
    flow_from_loaded, from_loaded, registry, wall_law_from_loaded,
};
