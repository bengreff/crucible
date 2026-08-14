//! FND-4 §3.2 — the two artifacts of §3.1: the **author config** (intent,
//! may omit defaulted values) and the **resolved config** (every effective
//! value materialized, §3.5). The resolved form serializes back to the same
//! top-level grammar, so `resolve` is a checkable fixed point:
//! `load(serialize(resolved))` must reproduce `resolved` exactly (§6 test 1).
//!
//! The top level is typed with `deny_unknown_fields` (unknown key = hard
//! error, §3.2); `type`-keyed block *bodies* are open here and validated
//! against the COUP-8 registry `Manifest.params` in the semantic pass (§3.3).

use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// The author document (§3.1: intent — intentionally underspecifiable).
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct AuthorConfig {
    pub schema_version: i64,
    #[serde(default)]
    pub meta: Option<MetaBlock>,
    #[serde(default)]
    pub geometry: Option<GeometryBlock>,
    #[serde(default)]
    pub materials: BTreeMap<String, TypedBlock>,
    #[serde(default)]
    pub mechanisms: BTreeMap<String, TypedBlock>,
    /// Grammar deferred (COUP-2 coupler selections) — must be empty to load.
    #[serde(default)]
    pub couplers: toml::Table,
    #[serde(default)]
    pub engine: Option<EngineBlock>,
    /// Grammar deferred — must be empty to load.
    #[serde(default)]
    pub operating_profile: toml::Table,
    #[serde(default)]
    pub determinism: Option<DeterminismBlock>,
    /// Pins are FND-5-resolved; until that lands, must be empty to load.
    #[serde(default)]
    pub tables: BTreeMap<String, String>,
    /// Grammar deferred (COUP-5) — must be empty to load.
    #[serde(default)]
    pub uq: toml::Table,
    #[serde(default)]
    pub rng: Option<RngBlock>,
}

#[derive(Debug, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct MetaBlock {
    #[serde(default)]
    pub name: Option<String>,
    #[serde(default)]
    pub description: Option<String>,
}

/// `[geometry]` — grammar owned by FND-3; only the pieces the loader checks
/// now are typed (the CSG/STL grammar lands with the FND-3 session).
/// `n_theta_max` is the config-declared finest azimuthal resolution
/// (FND-3 §3.3), gated by the θ-ladder rule (FND-4 §3.4-5b / FND-2 §3.4).
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct GeometryBlock {
    #[serde(default)]
    pub n_theta_max: Option<i64>,
}

/// A `type`-keyed registry-dispatched block (§3.3). Body deliberately open:
/// params are validated against the mechanism's `Manifest.params`, not serde.
#[derive(Debug, Deserialize)]
pub(crate) struct TypedBlock {
    #[serde(rename = "type")]
    pub type_id: String,
    #[serde(flatten)]
    pub params: toml::Table,
}

/// `[engine]` — full composition grammar deferred; the O20 explicit-binding
/// semantics are fixed now: `[engine.bindings]` maps
/// `"<instance>.<require_port>"` → `"<provider instance>.<provide_port>"`.
#[derive(Debug, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct EngineBlock {
    #[serde(default)]
    pub bindings: BTreeMap<String, String>,
}

/// `[determinism]` (O21) — the S6 regime→guarantee declaration surface.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct DeterminismBlock {
    pub mode: String,
}

#[derive(Debug, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct RngBlock {
    /// TOML integers are i64, so the seed space is 63-bit (drawn seeds are
    /// masked to keep the resolved config TOML-representable).
    #[serde(default)]
    pub master_seed: Option<i64>,
    #[serde(default)]
    pub algorithm: Option<String>,
}

// ---------------------------------------------------------------------------
// Resolved form (§3.5): zero implicit values, serializable, fixed point.
// ---------------------------------------------------------------------------

/// The fully-resolved config: every effective default materialized. This is
/// what the manifest embeds (§3.6) and what regeneration replays.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ResolvedConfig {
    pub schema_version: i64,
    pub meta: ResolvedMeta,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub geometry: Option<ResolvedGeometry>,
    /// Block bodies as full TOML tables (`type` key + every param explicit).
    pub materials: BTreeMap<String, toml::Table>,
    pub mechanisms: BTreeMap<String, toml::Table>,
    pub engine: ResolvedEngine,
    pub determinism: ResolvedDeterminism,
    pub rng: ResolvedRng,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ResolvedMeta {
    pub name: String,
    pub description: String,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ResolvedGeometry {
    pub n_theta_max: i64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ResolvedEngine {
    pub bindings: BTreeMap<String, String>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ResolvedDeterminism {
    pub mode: String,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ResolvedRng {
    pub master_seed: i64,
    pub algorithm: String,
}

impl ResolvedConfig {
    /// Serialize to author-grammar TOML (the §6-1 round-trip surface).
    /// All maps are `BTreeMap` and struct fields have fixed order, so this
    /// serialization is canonical by construction (§3.8).
    pub fn to_toml(&self) -> String {
        toml::to_string(self).expect("resolved config is TOML-representable by construction")
    }
}
