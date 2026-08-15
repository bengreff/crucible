//! SOLV-family operators on the world-state grid. First inhabitant: the
//! **conduction operator** (SOLV-1 §3.5's diffusion class on the FND-2
//! cylindrical metric) plus the project's first real COUP-8 registry entry,
//! and the **convergence-certificate studies** (the Goal-A deliverable).
//!
//! Time integration note (honest scaffolding): the production advance is
//! COUP-3's SDC-coupled IMEX with the deterministic fixed-cycle implicit
//! diffusion solve (class `D`). This crate's explicit fixed-order reference
//! integrator exists to drive the certificate's *spatial-operator*
//! verification (MMS + analytic anchors) and is superseded — not extended —
//! when COUP-3 lands. The flux-form spatial operator carries over unchanged.
//!
//! Session scope: uniform N_θ across bricks per sweep (asserted); the
//! conservative flux aggregation across an N_θ jump (AMR-refluxing style,
//! FND-2 §3.4) arrives with the COUP-2/COUP-3 wave.

pub mod certificate;
mod conduction;
mod mechanism;

pub use conduction::{Bcs, Conduction, FaceBc, SolverError};
pub use mechanism::{CONDUCTION_MANIFEST, ConductionSetup, SetupError, from_loaded, registry};
