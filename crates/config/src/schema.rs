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
    /// Minimal steady-march subset of the operating-profile grammar (the
    /// complete grammar remains deferred per the FND-4 skeleton split).
    #[serde(default)]
    pub operating_profile: Option<ProfileBlock>,
    #[serde(default)]
    pub determinism: Option<DeterminismBlock>,
    /// FND-4 §6-4 `[tables]`: logical name → pinned artifact reference.
    #[serde(default)]
    pub tables: BTreeMap<String, TableRefBlock>,
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
    /// FND-3 contour-of-revolution grammar (stair degenerate form): a cited
    /// r(z) CSV (columns `z_*,area_ratio,r_*,source`, `#` comments), from
    /// which the grid extents are DERIVED via the fidelity dial and
    /// materialized into the resolved form (manifest-recorded). The CSV is
    /// content-addressed: its digest is resolved at load (like a table pin)
    /// and re-verified by the engine assembly at build.
    #[serde(default)]
    pub contour: Option<String>,
    /// Units of the CSV's length columns: "in" (converted ×0.0254 exactly,
    /// META-1 §3 edge conversion) or "m".
    #[serde(default)]
    pub contour_units: Option<String>,
    /// sha256 of the CSV content; author-optional (resolved always carries
    /// it — its presence alongside explicit extents marks the replay form,
    /// which needs no file access).
    #[serde(default)]
    pub contour_digest: Option<String>,
    /// Liner ring thickness [m] added outside the contour (0 = no solid
    /// region: cold-flow geometry).
    #[serde(default)]
    pub liner_thickness_m: Option<f64>,
    /// THE compute-fidelity dial: cell size = r_throat / this, isotropic in
    /// (r, z). Physically anchored to the throat radius, so it transfers to
    /// any engine contour; derived n_r/n_z land in the manifest.
    #[serde(default)]
    pub cells_across_throat: Option<f64>,
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

/// One `[tables.<name>]` entry (FND-4 §6-4). Author intent may state the pin
/// two ways: **explicitly** (`data_version` + `content_digest`, both — the
/// grammar the resolved form always emits, so replay needs no file access) or
/// **via the sidecar** (`pins` = path of the machine-written `…pins.toml`,
/// the single pin owner; resolution reads the entry for `group`). One of the
/// two forms is mandatory — a pin-less table reference is a load error, and
/// FND-5's open gate re-verifies the digest against the actual bytes at bind.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct TableRefBlock {
    pub file: String,
    pub group: String,
    #[serde(default)]
    pub pins: Option<String>,
    #[serde(default)]
    pub data_version: Option<String>,
    #[serde(default)]
    pub content_digest: Option<String>,
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

/// `[operating_profile]` — the steady-march subset: march the transient to
/// a fixed settle budget (deterministic, never wall-clock), budget stated in
/// injector-state flow-through times.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct ProfileBlock {
    pub mode: String,
    #[serde(default)]
    pub flowthroughs: Option<f64>,
    #[serde(default)]
    pub cfl: Option<f64>,
    /// Startup device: quiescent fill pressure [Pa] the march starts from
    /// (drains smoothly to the emergent operating point; the steady result
    /// is fill-independent).
    #[serde(default)]
    pub fill_p_pa: Option<f64>,
    /// Startup device: the ambient pump-down window, in flow-through times
    /// (the altitude-cell schedule — back-pressure falls log-linearly from
    /// the fill to the vacuum floor over this window, so the nozzle
    /// establishes quasi-statically).
    #[serde(default)]
    pub pumpdown_flowthroughs: Option<f64>,
    /// The declared altitude-cell ambient floor [Pa] the pump-down bottoms
    /// out at (session 12): a real test cell holds finite pressure (~1 mbar
    /// class), and declaring it keeps the settled plume fringe inside the
    /// equilibrium surface's cold envelope instead of pinning its edge.
    #[serde(default)]
    pub p_amb_floor_pa: Option<f64>,
    /// Startup device (session 12): the injector mass-flow ramp window, in
    /// flow-through times — ṁ(t) = ṁ·min(1, t/t_ramp), the declared start
    /// schedule (a real engine's valve sequence). Kills the injector-piston
    /// shock a step start drives into the fill gas; with a modest fill
    /// pressure the establishment is quasi-static end to end. 0 = step
    /// start (the pre-session-12 behavior).
    #[serde(default)]
    pub injector_ramp_flowthroughs: Option<f64>,
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
    /// Steady-march operating profile, when declared (subset grammar).
    #[serde(skip_serializing_if = "Option::is_none")]
    #[serde(default)]
    pub operating_profile: Option<ResolvedProfile>,
    pub determinism: ResolvedDeterminism,
    /// FND-4 §6-4: every pin fully resolved (explicit version + digest), so
    /// replay of the resolved form needs no sidecar access.
    #[serde(default)]
    pub tables: BTreeMap<String, ResolvedTablePin>,
    pub rng: ResolvedRng,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ResolvedProfile {
    pub mode: String,
    pub flowthroughs: f64,
    pub cfl: f64,
    pub fill_p_pa: f64,
    pub pumpdown_flowthroughs: f64,
    pub p_amb_floor_pa: f64,
    pub injector_ramp_flowthroughs: f64,
}

/// A fully-resolved table pin — the §3.6 regeneration-key row. `pins`
/// records which sidecar resolved it (provenance; `None` = the author stated
/// the pair explicitly).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ResolvedTablePin {
    pub file: String,
    pub group: String,
    /// PROVENANCE ONLY (session-12 review): the sidecar this pin was
    /// originally resolved from, when the author named one. Whenever an
    /// explicit `data_version`/`content_digest` pair is present it is
    /// AUTHORITATIVE and the sidecar is NOT consulted — the resolved form
    /// carries pair + this field and must replay with zero sidecar I/O
    /// (§3.5 fixed point). An author who states a pair alongside a
    /// sidecar is therefore declaring the pair, not the sidecar.
    #[serde(skip_serializing_if = "Option::is_none")]
    #[serde(default)]
    pub pins: Option<String>,
    pub data_version: String,
    pub content_digest: String,
}

impl ResolvedTablePin {
    /// The FND-5 open-gate pin for this entry (the loader-side verification
    /// key: version + digest, both enforced at `Table::open`).
    pub fn fnd5_pin(&self) -> (String, Option<String>) {
        (self.data_version.clone(), Some(self.content_digest.clone()))
    }
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
    /// Contour-of-revolution declaration, digest resolved; flattened for the
    /// same replay reason. When present, `extents` holds the DERIVED grid
    /// (dial → cell size → counts), so replay never re-reads the CSV.
    #[serde(flatten)]
    #[serde(default)]
    pub contour: Option<ResolvedContour>,
}

/// The resolved contour declaration — author-shaped keys, all materialized.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ResolvedContour {
    pub contour: String,
    pub contour_units: String,
    pub contour_digest: String,
    pub liner_thickness_m: f64,
    pub cells_across_throat: f64,
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

/// META-2 §4 ★ — dimensioned parameter access at the parse boundary. The
/// FND-4 param-name convention already spells the unit in the name
/// (`kappa_w_per_m_k`); these accessors are the enforcement: the value
/// leaves the config layer as a typed quantity, so a consumer wanting the
/// wrong dimension fails to compile instead of misreading a number.
/// Dimensionless params (γ, Pr, Z) stay on [`ResolvedConfig::param_f64`].
/// Each accessor is named for its unit — it must agree with the param-name
/// suffix it reads, and the pairing is the reviewable surface.
impl ResolvedConfig {
    pub fn param_length_m(block: &toml::Table, name: &str) -> Option<crucible_units::Length> {
        Self::param_f64(block, name).map(crucible_units::length_m)
    }

    pub fn param_temperature_k(
        block: &toml::Table,
        name: &str,
    ) -> Option<crucible_units::ThermodynamicTemperature> {
        Self::param_f64(block, name).map(crucible_units::temperature_k)
    }

    pub fn param_pressure_pa(block: &toml::Table, name: &str) -> Option<crucible_units::Pressure> {
        Self::param_f64(block, name).map(crucible_units::pressure_pa)
    }

    pub fn param_thermal_conductivity_w_per_m_k(
        block: &toml::Table,
        name: &str,
    ) -> Option<crucible_units::ThermalConductivity> {
        Self::param_f64(block, name).map(crucible_units::thermal_conductivity_w_per_m_k)
    }

    pub fn param_volumetric_heat_capacity_j_per_m3_k(
        block: &toml::Table,
        name: &str,
    ) -> Option<crucible_units::VolumetricHeatCapacity> {
        Self::param_f64(block, name).map(crucible_units::volumetric_heat_capacity_j_per_m3_k)
    }

    pub fn param_heat_transfer_w_per_m2_k(
        block: &toml::Table,
        name: &str,
    ) -> Option<crucible_units::HeatTransfer> {
        Self::param_f64(block, name).map(crucible_units::heat_transfer_w_per_m2_k)
    }

    pub fn param_dynamic_viscosity_pa_s(
        block: &toml::Table,
        name: &str,
    ) -> Option<crucible_units::DynamicViscosity> {
        Self::param_f64(block, name).map(crucible_units::dynamic_viscosity_pa_s)
    }

    pub fn param_specific_heat_capacity_j_per_kg_k(
        block: &toml::Table,
        name: &str,
    ) -> Option<crucible_units::SpecificHeatCapacity> {
        Self::param_f64(block, name).map(crucible_units::specific_heat_capacity_j_per_kg_k)
    }
}

impl ResolvedExtents {
    /// Typed views of the geometry extents (the config→grid seam). The
    /// serialized struct itself stays plain SI `f64` — the §6-1 replay
    /// grammar must not change shape under a units library.
    pub fn r_min_length(&self) -> crucible_units::Length {
        crucible_units::length_m(self.r_min)
    }
    pub fn dr_length(&self) -> crucible_units::Length {
        crucible_units::length_m(self.dr)
    }
    pub fn z_min_length(&self) -> crucible_units::Length {
        crucible_units::length_m(self.z_min)
    }
    pub fn dz_length(&self) -> crucible_units::Length {
        crucible_units::length_m(self.dz)
    }
}
