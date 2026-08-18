//! Implements **FND-4 v0.2** — the author-config schema (§3.2), the loader
//! pipeline (§3.4), no-hidden-defaults resolution (§3.5), and the
//! config-derived run manifest (§3.6). The single source of a run's intent:
//! one TOML file; no environment variable, CLI overlay, or ambient default
//! may change a result.
//!
//! `[tables]` (§6-4) is live: a pin entry states the explicit
//! `data_version` + `content_digest` pair, or delegates to the
//! machine-written `…pins.toml` sidecar (the single pin owner) via
//! [`load_str_with_sidecars`]; resolved configs always carry the explicit
//! pair, so replay is a pure string load. FND-5's open gate re-verifies the
//! digest against the artifact bytes at bind.
//!
//! Session-scoped deferrals (each lands with its owner):
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
    OLDEST_SUPPORTED_SCHEMA_VERSION, load_str, load_str_with_sidecars,
};
pub use manifest::{ChaoticRecord, RunManifest, TablePin};
pub use schema::{ResolvedConfig, ResolvedExtents, ResolvedGeometry, ResolvedTablePin};
/// Resolved block bodies are `toml::Table`s; re-export the crate so
/// downstream consumers name those types without a version-skew risk.
pub use toml;
