//! FND-4 §3.4 — the loader pipeline. Ordered; either produces a
//! fully-specified config + manifest or **halts at load with all diagnoses**
//! (never a partial run):
//!
//! parse, then migrate, then structural deserialize (`deny_unknown_fields`,
//! first-error serde layer), then registry dispatch (COUP-8 `Manifest.params`),
//! then wiring validation + FND-4-owned value checks (error-accumulating),
//! then default resolution (§3.5, zero implicit values), then manifest (§3.6).

use crate::diag::Diagnostics;
use crate::manifest::{ChaoticRecord, MANIFEST_SCHEMA_VERSION, RunManifest};
use crate::schema::{
    AuthorConfig, ResolvedConfig, ResolvedDeterminism, ResolvedEngine, ResolvedGeometry,
    ResolvedMeta, ResolvedRng, TypedBlock,
};
use crucible_registry::{
    ChaoticClass, Direction, Manifest, ParamSpec, ParamType, ParamValue, PortRole, Registry,
};
use sha2::{Digest, Sha256};
use std::collections::BTreeMap;

pub const CURRENT_SCHEMA_VERSION: i64 = 1;
pub const OLDEST_SUPPORTED_SCHEMA_VERSION: i64 = 1;
/// FND-4 §3.2 (O21): fixed-order is the default determinism mode.
pub const DEFAULT_DETERMINISM_MODE: &str = "fixed-order";
/// The FND-1 §3.5 counter-based PRNG contract; Philox4x32-10 (Random123) is
/// the pinned engineering choice, implemented when COUP-5 first draws.
pub const DEFAULT_RNG_ALGORITHM: &str = "philox4x32-10";

#[derive(Debug, Clone, PartialEq)]
pub struct Loaded {
    pub resolved: ResolvedConfig,
    pub manifest: RunManifest,
}

/// Load an author config against the registry. The **only** inputs are the
/// two arguments — no environment variable, CLI overlay, wall-clock, or
/// ambient default may influence the result (§2 invariant; the auto-drawn
/// seed, when the author omits one, is OS entropy *recorded into the
/// manifest*, never an input to any other resolution).
pub fn load_str(author_toml: &str, registry: &Registry) -> Result<Loaded, Diagnostics> {
    let mut diags = Diagnostics::new();

    // 1. Parse.
    let doc: toml::Value = match author_toml.parse() {
        Ok(v) => v,
        Err(e) => {
            diags.push("", format!("TOML parse error: {e}"));
            return Err(diags.into_sorted());
        }
    };

    // 2. Schema version gate + migration chain (chain is empty at v1).
    let schema_version = match doc.get("schema_version").and_then(toml::Value::as_integer) {
        Some(v) => v,
        None => {
            diags.push(
                "schema_version",
                "missing or non-integer `schema_version` (mandatory, §3.7)".to_string(),
            );
            return Err(diags.into_sorted());
        }
    };
    if schema_version > CURRENT_SCHEMA_VERSION {
        diags.push(
            "schema_version",
            format!(
                "config schema_version {schema_version} is newer than this build understands \
                 (current: {CURRENT_SCHEMA_VERSION}) — refusing loudly (§3.7)"
            ),
        );
        return Err(diags.into_sorted());
    }
    if schema_version < OLDEST_SUPPORTED_SCHEMA_VERSION {
        diags.push(
            "schema_version",
            format!(
                "config schema_version {schema_version} is older than the oldest supported \
                 ({OLDEST_SUPPORTED_SCHEMA_VERSION}) — refusing loudly (§3.7)"
            ),
        );
        return Err(diags.into_sorted());
    }
    let migrations_applied: Vec<String> = Vec::new(); // ordered chain; empty until a v2 exists

    // 3. Structural deserialize (serde stops at the first error — deliberate
    // two-layer design, §3.4; the semantic pass below accumulates).
    let author: AuthorConfig = match toml::from_str(author_toml) {
        Ok(a) => a,
        Err(e) => {
            diags.push("", format!("structural error: {e}"));
            return Err(diags.into_sorted());
        }
    };
    debug_assert_eq!(
        author.schema_version, schema_version,
        "step-2 pre-parse and step-3 structural view must agree"
    );

    // 4+5. Registry dispatch + wiring/value checks, all errors collected.
    let mut resolved_materials: BTreeMap<String, toml::Table> = BTreeMap::new();
    let mut resolved_mechanisms: BTreeMap<String, toml::Table> = BTreeMap::new();
    let mut selected: BTreeMap<&str, &'static Manifest> = BTreeMap::new();

    for (name, block) in &author.materials {
        let path = format!("materials.{name}");
        check_instance_name(&mut diags, &path, name);
        match registry.material(&block.type_id) {
            None => diags.push(
                format!("{path}.type"),
                unknown_id_message("material", &block.type_id, registry.material_ids()),
            ),
            Some(m) => {
                check_interface_version(&mut diags, &path, m);
                if let Some(t) = validate_and_resolve_params(&mut diags, &path, block, m) {
                    resolved_materials.insert(name.clone(), t);
                }
            }
        }
    }
    for (name, block) in &author.mechanisms {
        let path = format!("mechanisms.{name}");
        check_instance_name(&mut diags, &path, name);
        match registry.mechanism(&block.type_id) {
            None => diags.push(
                format!("{path}.type"),
                unknown_id_message("mechanism", &block.type_id, registry.mechanism_ids()),
            ),
            Some(m) => {
                check_interface_version(&mut diags, &path, m);
                selected.insert(name.as_str(), m);
                if let Some(t) = validate_and_resolve_params(&mut diags, &path, block, m) {
                    resolved_mechanisms.insert(name.clone(), t);
                }
            }
        }
    }

    // Deferred-grammar blocks refuse when non-empty (fail loud, lib.rs note).
    for (path, table) in [
        ("couplers", &author.couplers),
        ("operating_profile", &author.operating_profile),
        ("uq", &author.uq),
    ] {
        if !table.is_empty() {
            diags.push(
                path,
                format!("`[{path}]` grammar is deferred (FND-4 skeleton) — not yet implemented"),
            );
        }
    }
    if !author.tables.is_empty() {
        diags.push(
            "tables",
            "table pin resolution lands with the FND-5 loader — until then `[tables]` must be empty",
        );
    }

    // O21 — determinism mode + chaotic refusal.
    let mode = author
        .determinism
        .as_ref()
        .map_or(DEFAULT_DETERMINISM_MODE, |d| d.mode.as_str());
    match mode {
        "fixed-order" => {}
        "relaxed" => {
            let chaotic: Vec<String> = selected
                .iter()
                .filter(|(_, m)| m.any_chaotic())
                .map(|(name, m)| format!("{name} ({})", m.id))
                .collect();
            if !chaotic.is_empty() {
                diags.push(
                    "determinism.mode",
                    format!(
                        "`relaxed` refused: chaotic-classified mechanism(s) selected — {} — \
                         fixed-order reductions are mandatory in chaotic regimes (O21, META-1 §2.1)",
                        chaotic.join(", ")
                    ),
                );
            }
        }
        other => diags.push(
            "determinism.mode",
            format!("unknown mode {other:?}; expected \"fixed-order\" or \"relaxed\""),
        ),
    }

    // θ-ladder alignment (§3.4-5b): N_θ^max = 4·2^n so factor-2 coarsening
    // lands exactly on the guard resolution N_θ^guard = 4 (FND-2 §3.4).
    let resolved_geometry = match &author.geometry {
        None => None,
        Some(g) => {
            let axisymmetric = g.axisymmetric.unwrap_or(false);
            let n_theta = match g.n_theta_max {
                None => {
                    diags.push(
                        "geometry.n_theta_max",
                        "required when [geometry] is declared (no hidden defaults, §3.5)",
                    );
                    0
                }
                // N_θ = 1 is admissible only under the recorded axisymmetry
                // assertion (FND-2 §3.4, S4) — never silently.
                Some(1) if axisymmetric => 1,
                Some(1) => {
                    diags.push(
                        "geometry.n_theta_max",
                        "N_θ^max = 1 requires `axisymmetric = true` — the recorded, \
                         pedigree-visible axisymmetry assertion (FND-2 §3.4)",
                    );
                    1
                }
                Some(v) => {
                    if !is_theta_ladder_aligned(v) {
                        diags.push("geometry.n_theta_max", theta_ladder_message(v));
                    }
                    v
                }
            };
            // Extents: all six or none (a half-declared world is a fault).
            let ext = [
                g.r_min.is_some(),
                g.dr.is_some(),
                g.n_r.is_some(),
                g.z_min.is_some(),
                g.dz.is_some(),
                g.n_z.is_some(),
            ];
            let extents = if ext.iter().all(|&p| p) {
                let (r_min, dr, n_r) = (g.r_min.unwrap(), g.dr.unwrap(), g.n_r.unwrap());
                let (z_min, dz, n_z) = (g.z_min.unwrap(), g.dz.unwrap(), g.n_z.unwrap());
                if r_min.is_nan() || r_min < 0.0 || dr <= 0.0 || dz <= 0.0 || n_r < 1 || n_z < 1 {
                    diags.push(
                        "geometry",
                        "extents need r_min ≥ 0, dr > 0, dz > 0, n_r ≥ 1, n_z ≥ 1",
                    );
                }
                Some(crate::schema::ResolvedExtents {
                    r_min,
                    dr,
                    n_r,
                    z_min,
                    dz,
                    n_z,
                })
            } else {
                if ext.iter().any(|&p| p) {
                    diags.push(
                        "geometry",
                        "grid extents are all-or-none: declare r_min, dr, n_r, z_min, dz, n_z \
                         together (FND-2 §3.2), or omit all six",
                    );
                }
                None
            };
            Some(ResolvedGeometry {
                n_theta_max: n_theta,
                axisymmetric,
                extents,
            })
        }
    };

    // O20 — explicit Require→provider bindings; nothing auto-matches.
    let bindings = author
        .engine
        .as_ref()
        .map(|e| e.bindings.clone())
        .unwrap_or_default();
    check_bindings(&mut diags, &bindings, &selected);

    // Coupler source/sink balance (COUP-8 §3.3-3) over the selected set.
    check_coupler_balance(&mut diags, &selected);

    // RNG block value checks.
    if let Some(rng) = &author.rng
        && let Some(seed) = rng.master_seed
        && seed < 0
    {
        diags.push("rng.master_seed", "must be non-negative");
    }

    if !diags.is_empty() {
        return Err(diags.into_sorted());
    }

    // 6. Resolve defaults — zero implicit values (§3.5).
    let master_seed = match author.rng.as_ref().and_then(|r| r.master_seed) {
        Some(s) => s,
        None => match draw_seed() {
            Ok(s) => s,
            Err(e) => {
                diags.push("rng.master_seed", format!("auto-draw failed: {e}"));
                return Err(diags.into_sorted());
            }
        },
    };
    let resolved = ResolvedConfig {
        schema_version: CURRENT_SCHEMA_VERSION,
        meta: ResolvedMeta {
            name: author
                .meta
                .as_ref()
                .and_then(|m| m.name.clone())
                .unwrap_or_default(),
            description: author
                .meta
                .as_ref()
                .and_then(|m| m.description.clone())
                .unwrap_or_default(),
        },
        geometry: resolved_geometry,
        materials: resolved_materials,
        mechanisms: resolved_mechanisms,
        engine: ResolvedEngine { bindings },
        determinism: ResolvedDeterminism {
            mode: mode.to_string(),
        },
        rng: ResolvedRng {
            master_seed,
            algorithm: author
                .rng
                .as_ref()
                .and_then(|r| r.algorithm.clone())
                .unwrap_or_else(|| DEFAULT_RNG_ALGORITHM.to_string()),
        },
    };

    // 7. Manifest (§3.6).
    let chaotic_class_in_force = selected
        .iter()
        .map(|(name, m)| ChaoticRecord {
            instance: (*name).to_string(),
            mechanism_id: m.id.to_string(),
            classes: m
                .chaotic_class
                .iter()
                .map(|(r, c)| {
                    (
                        format!("{r:?}"),
                        match c {
                            ChaoticClass::Chaotic => "chaotic".to_string(),
                            ChaoticClass::NonChaotic => "non_chaotic".to_string(),
                        },
                    )
                })
                .collect(),
        })
        .collect();
    let manifest = RunManifest {
        manifest_schema_version: MANIFEST_SCHEMA_VERSION,
        config_content_hash: content_hash(&doc),
        schema_version,
        migrations_applied,
        determinism_mode: resolved.determinism.mode.clone(),
        chaotic_class_in_force,
        master_seed,
        rng_algorithm: resolved.rng.algorithm.clone(),
        table_pins: Vec::new(),
        resolved_config: resolved.clone(),
    };

    Ok(Loaded { resolved, manifest })
}

// ---------------------------------------------------------------------------
// Checks
// ---------------------------------------------------------------------------

fn check_instance_name(diags: &mut Diagnostics, path: &str, name: &str) {
    let ok = !name.is_empty()
        && name
            .chars()
            .all(|c| c.is_ascii_alphanumeric() || c == '_' || c == '-');
    if !ok {
        diags.push(
            path.to_string(),
            "instance names must match [A-Za-z0-9_-]+ (`.` is the O20 binding separator)",
        );
    }
}

fn check_interface_version(diags: &mut Diagnostics, path: &str, m: &Manifest) {
    if !m.interface_version.supported_by_core() {
        diags.push(
            path.to_string(),
            format!(
                "mechanism {} built against interface {}.{}; core requires {}.{}–{}.{} (COUP-8 §3.3-4)",
                m.id,
                m.interface_version.major,
                m.interface_version.minor,
                crucible_registry::CORE_INTERFACE_MAJOR,
                crucible_registry::CORE_SUPPORTED_MINOR_MIN,
                crucible_registry::CORE_INTERFACE_MAJOR,
                crucible_registry::CORE_SUPPORTED_MINOR_MAX,
            ),
        );
    }
}

fn unknown_id_message<'a>(kind: &str, got: &str, valid: impl Iterator<Item = &'a str>) -> String {
    let ids: Vec<&str> = valid.collect();
    if ids.is_empty() {
        format!("unknown {kind} type {got:?}; the registry is empty")
    } else {
        format!("unknown {kind} type {got:?}; valid ids: {}", ids.join(", "))
    }
}

/// Validate a `type`-keyed block body against `Manifest.params` and return
/// the fully-materialized table (`type` + every param explicit, defaults from
/// the registry descriptor — §3.3/§3.5). Numeric params declared Float accept
/// TOML integers and are normalized to floats (canonical resolved form).
fn validate_and_resolve_params(
    diags: &mut Diagnostics,
    path: &str,
    block: &TypedBlock,
    m: &Manifest,
) -> Option<toml::Table> {
    let mut ok = true;
    // Unknown params (author keys not in the schema), in sorted order.
    for key in block.params.keys() {
        if !m.params.iter().any(|p| p.name == key) {
            let valid: Vec<&str> = m.params.iter().map(|p| p.name).collect();
            diags.push(
                format!("{path}.{key}"),
                format!(
                    "unknown parameter for {}; valid: [{}]",
                    m.id,
                    valid.join(", ")
                ),
            );
            ok = false;
        }
    }
    // Declared params: presence, type, range; materialize value-or-default.
    let mut out: BTreeMap<String, toml::Value> = BTreeMap::new();
    out.insert("type".into(), toml::Value::String(m.id.to_string()));
    for spec in m.params {
        let ppath = format!("{path}.{}", spec.name);
        match (block.params.get(spec.name), spec.default) {
            (None, None) => {
                diags.push(ppath, format!("required parameter of {} is missing", m.id));
                ok = false;
            }
            (None, Some(d)) => {
                out.insert(spec.name.into(), default_to_value(d));
            }
            (Some(v), _) => match coerce(v, spec) {
                Ok(cv) => {
                    check_range(diags, &ppath, &cv, spec, &mut ok);
                    out.insert(spec.name.into(), cv);
                }
                Err(want) => {
                    diags.push(ppath, format!("expected {want}, got {}", v.type_str()));
                    ok = false;
                }
            },
        }
    }
    if !ok {
        return None;
    }
    let mut table = toml::Table::new();
    for (k, v) in out {
        table.insert(k, v);
    }
    Some(table)
}

fn default_to_value(d: ParamValue) -> toml::Value {
    match d {
        ParamValue::Bool(b) => toml::Value::Boolean(b),
        ParamValue::Int(i) => toml::Value::Integer(i),
        ParamValue::Float(x) => toml::Value::Float(x),
        ParamValue::Str(s) => toml::Value::String(s.to_string()),
    }
}

fn coerce(v: &toml::Value, spec: &ParamSpec) -> Result<toml::Value, &'static str> {
    match (spec.ty, v) {
        (ParamType::Bool, toml::Value::Boolean(_))
        | (ParamType::Int, toml::Value::Integer(_))
        | (ParamType::Float, toml::Value::Float(_))
        | (ParamType::Str, toml::Value::String(_)) => Ok(v.clone()),
        // Float params accept integer literals, normalized (canonical form).
        (ParamType::Float, toml::Value::Integer(i)) => {
            let f = *i as f64;
            Ok(toml::Value::Float(f))
        }
        (ParamType::Bool, _) => Err("boolean"),
        (ParamType::Int, _) => Err("integer"),
        (ParamType::Float, _) => Err("float"),
        (ParamType::Str, _) => Err("string"),
    }
}

fn check_range(
    diags: &mut Diagnostics,
    ppath: &str,
    v: &toml::Value,
    spec: &ParamSpec,
    ok: &mut bool,
) {
    let Some((lo, hi)) = spec.range else { return };
    let x = match v {
        toml::Value::Integer(i) => *i as f64,
        toml::Value::Float(f) => *f,
        _ => return,
    };
    if !(lo..=hi).contains(&x) {
        diags.push(
            ppath.to_string(),
            format!("value {x} outside validity range [{lo}, {hi}]"),
        );
        *ok = false;
    }
}

fn is_theta_ladder_aligned(v: i64) -> bool {
    v >= 4 && v % 4 == 0 && ((v / 4) as u64).is_power_of_two()
}

fn theta_ladder_message(v: i64) -> String {
    // Nearest admissible 4·2^n below and above v.
    let mut lower: Option<i64> = None;
    let mut upper: i64 = 4;
    while upper < v {
        lower = Some(upper);
        upper *= 2;
    }
    if upper == v {
        upper *= 2; // unreachable for misaligned v; defensive
    }
    let near = match lower {
        Some(l) => format!("{l} and {upper}"),
        None => format!("{upper}"),
    };
    format!(
        "N_θ^max = {v} is not 4·2^n (n ≥ 0); the factor-2 coarsening ladder must land on \
         N_θ^guard = 4 (FND-2 §3.4). Nearest admissible: {near}"
    )
}

fn check_bindings(
    diags: &mut Diagnostics,
    bindings: &BTreeMap<String, String>,
    selected: &BTreeMap<&str, &'static Manifest>,
) {
    // Every Require port explicitly bound; candidates listed on a miss (O20).
    for (inst, m) in selected {
        for port in m.ports.iter().filter(|p| p.role == PortRole::Require) {
            let key = format!("{inst}.{}", port.name);
            let path = format!("engine.bindings.{key}");
            match bindings.get(&key) {
                None => {
                    let candidates: Vec<String> = selected
                        .iter()
                        .flat_map(|(pi, pm)| {
                            pm.ports
                                .iter()
                                .filter(|pp| pp.role == PortRole::Provide && pp.kind == port.kind)
                                .map(move |pp| format!("{pi}.{}", pp.name))
                        })
                        .collect();
                    let hint = if candidates.is_empty() {
                        "no candidate providers among selected mechanisms".to_string()
                    } else {
                        format!("candidate providers: {}", candidates.join(", "))
                    };
                    diags.push(
                        path,
                        format!(
                            "Require port {:?} ({:?}) has no binding — nothing auto-matches (O20); {hint}",
                            port.name, port.kind
                        ),
                    );
                }
                Some(target) => {
                    let Some((pi, pp_name)) = target.split_once('.') else {
                        diags.push(path, format!("malformed target {target:?}; expected \"<instance>.<provide_port>\""));
                        continue;
                    };
                    let Some(pm) = selected.get(pi) else {
                        diags.push(
                            path,
                            format!("target instance {pi:?} is not a selected mechanism"),
                        );
                        continue;
                    };
                    let Some(pp) = pm
                        .ports
                        .iter()
                        .find(|p| p.name == pp_name && p.role == PortRole::Provide)
                    else {
                        diags.push(
                            path,
                            format!("{pi:?} has no Provide port named {pp_name:?}"),
                        );
                        continue;
                    };
                    if pp.kind != port.kind {
                        diags.push(
                            path,
                            format!(
                                "kind mismatch: Require is {:?}, provider {}.{} is {:?}",
                                port.kind, pi, pp_name, pp.kind
                            ),
                        );
                    }
                }
            }
        }
    }
    // Every binding key must name an actual Require port of a selected instance.
    for key in bindings.keys() {
        let valid = key.split_once('.').is_some_and(|(inst, port)| {
            selected.get(inst).is_some_and(|m| {
                m.ports
                    .iter()
                    .any(|p| p.role == PortRole::Require && p.name == port)
            })
        });
        if !valid {
            diags.push(
                format!("engine.bindings.{key}"),
                "does not name a Require port of any selected mechanism",
            );
        }
    }
}

fn check_coupler_balance(diags: &mut Diagnostics, selected: &BTreeMap<&str, &'static Manifest>) {
    // COUP-8 §3.3-3: over the closed (kind × quantity) space, every Source
    // needs ≥1 Sink and vice-versa. Keys are Debug-formatted enum names —
    // deterministic, and the enums are closed.
    let mut sources: BTreeMap<String, Vec<String>> = BTreeMap::new();
    let mut sinks: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for (inst, m) in selected {
        for c in m.couplers {
            let key = format!("{:?}/{:?}", c.kind, c.quantity);
            if matches!(c.direction, Direction::Source | Direction::Bidirectional) {
                sources
                    .entry(key.clone())
                    .or_default()
                    .push((*inst).to_string());
            }
            if matches!(c.direction, Direction::Sink | Direction::Bidirectional) {
                sinks.entry(key).or_default().push((*inst).to_string());
            }
        }
    }
    for (key, insts) in &sources {
        if !sinks.contains_key(key) {
            diags.push(
                "couplers",
                format!(
                    "one-sided coupler edge {key}: Source(s) [{}] with no Sink (COUP-8 §3.3-3)",
                    insts.join(", ")
                ),
            );
        }
    }
    for (key, insts) in &sinks {
        if !sources.contains_key(key) {
            diags.push(
                "couplers",
                format!(
                    "one-sided coupler edge {key}: Sink(s) [{}] with no Source (COUP-8 §3.3-3)",
                    insts.join(", ")
                ),
            );
        }
    }
}

// ---------------------------------------------------------------------------
// Canonical hash & seed draw
// ---------------------------------------------------------------------------

/// §3.8: hash the canonicalized serialization (recursively sorted key order)
/// so a semantically identical config always hashes identically.
fn content_hash(doc: &toml::Value) -> String {
    let canonical = toml::to_string(&sort_value(doc.clone()))
        .expect("parsed TOML document re-serializes by construction");
    let mut h = Sha256::new();
    h.update(canonical.as_bytes());
    let bytes = h.finalize();
    let mut hex = String::with_capacity(64);
    for b in bytes {
        use std::fmt::Write;
        write!(hex, "{b:02x}").expect("writing to String cannot fail");
    }
    hex
}

fn sort_value(v: toml::Value) -> toml::Value {
    match v {
        toml::Value::Table(t) => {
            let mut sorted: Vec<(String, toml::Value)> = t.into_iter().collect();
            sorted.sort_by(|a, b| a.0.cmp(&b.0));
            let mut out = toml::Table::new();
            for (k, val) in sorted {
                out.insert(k, sort_value(val));
            }
            toml::Value::Table(out)
        }
        toml::Value::Array(a) => toml::Value::Array(a.into_iter().map(sort_value).collect()),
        other => other,
    }
}

/// Auto-draw a master seed from OS entropy (§3.1: the author may leave the
/// seed to be auto-drawn; it is recorded in the manifest, never wall-clock).
/// Masked to 63 bits so the resolved config stays TOML-representable.
fn draw_seed() -> Result<i64, getrandom::Error> {
    let mut buf = [0u8; 8];
    getrandom::fill(&mut buf)?;
    Ok((u64::from_le_bytes(buf) & (i64::MAX as u64)) as i64)
}
