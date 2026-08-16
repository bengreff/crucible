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
/// (FND-3 §3.3), gated by the θ-ladder rule (FND-4 §3.4-5b / FND-2 §3.4);
/// `axisymmetric = true` is the FND-2 §3.4 **recorded axisymmetry
/// assertion** (pedigree-visible) that alone admits `n_theta_max = 1`.
/// The explicit grid extents (uniform spacings, FND-2 §3.2) are the minimal
/// pre-CSG world declaration: all six present, or none.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct GeometryBlock {
    #[serde(default)]
    pub n_theta_max: Option<i64>,
    #[serde(default)]
    pub axisymmetric: Option<bool>,
    #[serde(default)]
    pub r_min: Option<f64>,
    #[serde(default)]
    pub dr: Option<f64>,
    #[serde(default)]
    pub n_r: Option<i64>,
    #[serde(default)]
    pub z_min: Option<f64>,
    #[serde(default)]
    pub dz: Option<f64>,
    #[serde(default)]
    pub n_z: Option<i64>,
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
    /// The recorded, pedigree-visible axisymmetry assertion (FND-2 §3.4).
    pub axisymmetric: bool,
    /// Flattened so the resolved form serializes to the SAME flat grammar
    /// the author `GeometryBlock` parses — the §6-1 fixed-point/replay
    /// contract requires load(serialize(resolved)) to succeed (review
    /// finding: a nested `[geometry.extents]` table broke replay for every
    /// extents-declaring config).
    #[serde(flatten)]
    pub extents: Option<ResolvedExtents>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ResolvedExtents {
    pub r_min: f64,
    pub dr: f64,
    pub n_r: i64,
    pub z_min: f64,
    pub dz: f64,
    pub n_z: i64,
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

    /// All mechanism instances of one registry type, in canonical
    /// (`BTreeMap`) order — the seam every solver's setup goes through, so
    /// no consumer hand-scans resolved TOML (and none can silently pick
    /// "the first" of several instances; count before you choose).
    pub fn mechanisms_of_type<'a>(&'a self, type_id: &str) -> Vec<(&'a str, &'a toml::Table)> {
        self.mechanisms
            .iter()
            .filter(|(_, b)| b.get("type").and_then(|v| v.as_str()) == Some(type_id))
            .map(|(name, b)| (name.as_str(), b))
            .collect()
    }

    /// Loader-validated parameter access on a resolved block: every param the
    /// registry `Manifest` declares is present (defaults materialized, §3.5),
    /// so `None` here means the caller asked for an undeclared name.
    pub fn param_f64(block: &toml::Table, name: &str) -> Option<f64> {
        block
            .get(name)
            .and_then(|v| v.as_float().or_else(|| v.as_integer().map(|i| i as f64)))
    }
}
