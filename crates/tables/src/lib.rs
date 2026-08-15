//! Implements **FND-5 v0.3** — the versioned-HDF5 table contract: the only
//! cross-language seam in the project (VISION_SCOPE §6). Python OFFL
//! pipelines write tables; this crate is everything the runtime knows about
//! reading them: the schema (§3.1), the load-time contract — version pins,
//! content digest, provenance, validity envelope (§3.2/§3.5) — and
//! deterministic N-D interpolation returning **scalars only** (§3.3; the
//! table's physical uncertainty rides the outer loop, META-1 §4).
//!
//! Session-scoped deferrals (each refuses loudly at load, never silently):
//! - `kind = "thermo_potential"` (Helmholtz/biquintic, S11) and
//!   `kind = "sample_set"` (SANDY realizations, FND-1 §3.9) — arrive with
//!   FND-7/OFFL-5 and OFFL-2 respectively (§6 test 7 lands then).
//! - `interp_method = "pchip"` (monotone-cubic shipped coefficients) —
//!   arrives with the first C¹ consumer (EOS/reactivity differencing).
//! - Config-side pin wiring (FND-4 `[tables]` → this loader) — lands when
//!   the first mechanism consumes a table (the conduction session).
//!
//! Determinism (§3.6): after load a `Table` is plain `f64` data; queries are
//! fixed-order arithmetic — no hash iteration, no stochastic interpolation,
//! no wall-clock. Same query, same bits, any thread count.

mod digest;
mod interp;
mod model;
mod reader;
pub mod writer;

pub use model::{
    Axis, EnvelopeHit, EnvelopePolicy, Provenance, TABLE_SCHEMA_MAJOR, Table, TableError,
    TableValue,
};
pub use reader::Pin;
