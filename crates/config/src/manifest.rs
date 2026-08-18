//! FND-4 §3.6 — the config-derived slice of the run manifest: the
//! regeneration key. FND-6 later adds the build/environment portion
//! (rustc/target/FP-policy/lockfile hash) and bundles it with results.

use crate::schema::ResolvedConfig;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

pub const MANIFEST_SCHEMA_VERSION: i64 = 1;

/// Per-instance record of the chaotic classification in force at load —
/// the record behind the O21 refusal (§3.6).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ChaoticRecord {
    pub instance: String,
    pub mechanism_id: String,
    /// regime name → `"chaotic"` / `"non_chaotic"`.
    pub classes: BTreeMap<String, String>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RunManifest {
    pub manifest_schema_version: i64,
    /// SHA-256 (hex) of the canonicalized author TOML (§3.8: stable key
    /// order — a semantically identical config always hashes identically).
    pub config_content_hash: String,
    /// Author-document `schema_version` and the migration chain applied.
    pub schema_version: i64,
    pub migrations_applied: Vec<String>,
    pub determinism_mode: String,
    pub chaotic_class_in_force: Vec<ChaoticRecord>,
    pub master_seed: i64,
    pub rng_algorithm: String,
    /// Table pins `{logical_name, resolved_version, content_hash}` (§3.6),
    /// resolved through the FND-4 §6-4 `[tables]` grammar (explicit pair or
    /// the machine-written pins sidecar).
    pub table_pins: Vec<TablePin>,
    /// The full resolved config — the actual replay input.
    pub resolved_config: ResolvedConfig,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct TablePin {
    pub logical_name: String,
    pub resolved_version: String,
    pub content_hash: String,
    /// Provenance read from the table itself when the run binds it (FND-5
    /// `Provenance.producer_version`) — `None` at load time; FND-6 records
    /// the bound value into the results bundle. Never fabricated here.
    #[serde(skip_serializing_if = "Option::is_none")]
    #[serde(default)]
    pub producing_generator_version: Option<String>,
}

impl RunManifest {
    pub fn to_toml(&self) -> String {
        toml::to_string(self).expect("manifest is TOML-representable by construction")
    }
}
