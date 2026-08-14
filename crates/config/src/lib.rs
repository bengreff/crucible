//! Implements **FND-4 v0.2** — the author-config schema (§3.2), the loader
//! pipeline (§3.4), no-hidden-defaults resolution (§3.5), and the
//! config-derived run manifest (§3.6). The single source of a run's intent:
//! one TOML file; no environment variable, CLI overlay, or ambient default
//! may change a result.
//!
//! Session-scoped deferrals (each lands with its owner):
//! - Table pin → version/content-hash resolution (§3.6): needs the FND-5
//!   loader; until then a non-empty `[tables]` block **refuses to load**
//!   (fail loud, never a silently unresolved pin). §6 test 4 lands then.
//! - `[couplers]`/`[operating_profile]`/`[uq]` grammars (deferred in the doc
//!   itself): non-empty blocks refuse to load.
//! - Envelope-coverage check (COUP-8 §3.3-2): needs FND-5 table metadata.
//! - Diagnostics carry deterministic dotted TOML paths (`mechanisms.hydro.type`),
//!   not yet byte spans — the §3.4 miette span annotation is a pending polish
//!   pass (structural serde errors do carry toml's native line/column).

mod diag;
mod loader;
mod manifest;
mod schema;

pub use diag::{Diagnostic, Diagnostics};
pub use loader::{
    CURRENT_SCHEMA_VERSION, DEFAULT_DETERMINISM_MODE, DEFAULT_RNG_ALGORITHM, Loaded,
    OLDEST_SUPPORTED_SCHEMA_VERSION, load_str,
};
pub use manifest::{ChaoticRecord, RunManifest};
pub use schema::ResolvedConfig;
