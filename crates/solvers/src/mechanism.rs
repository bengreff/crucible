//! The project's first real COUP-8 registry entry (`conduction`) and the
//! config→grid→operator wiring. `from_loaded` goes through the config
//! crate's typed seam (`mechanisms_of_type` / `param_f64`) — the pattern
//! every future mechanism reuses: no hand-scanning of resolved TOML, no
//! silent pick among multiple instances, no unchecked narrowing casts.
//! Boundary conditions and sources remain caller-supplied this session —
//! their config grammar belongs to the COUP-7 boundary-object wave.

use crucible_config::{Loaded, ResolvedConfig, ResolvedGeometry};
use crucible_grid::{Grid, GridError, GridSpec};
use crucible_registry::{
    ChaoticClass, InterfaceVersion, Manifest, ParamSpec, ParamType, Regime, Registry,
};

/// SOLV-1 §3.5 diffusion-class conduction: non-chaotic in every declared
/// regime (heat diffusion has no sensitive dependence), so it composes with
/// `relaxed` determinism mode once the GPU path exists (O21).
///
/// Validity ranges (refusal bounds, sourced): kappa spans known materials
/// from aerogels (~1e-2) to graphene-class conductors (~5e3 W/(m·K)) with
/// margin ×20 ⇒ (1e-12, 1e5]; rho_cp spans gases at vacuum-adjacent density
/// through dense solids (~4e6 J/(m³·K)) with margin ⇒ (1e-12, 1e12]. The
/// 1e-12 floors exclude zero/denormal-degenerate media that make the
/// diffusion operator singular. META-3 `materials-handbook` entry pending
/// FND-7; bounds tightened per-material when the spine lands.
pub static CONDUCTION_MANIFEST: Manifest = Manifest {
    id: "conduction",
    interface_version: InterfaceVersion { major: 1, minor: 0 },
    tables: &[],
    couplers: &[],
    ports: &[],
    chaotic_class: &[
        (Regime::Steady, ChaoticClass::NonChaotic),
        (Regime::Transient, ChaoticClass::NonChaotic),
    ],
    params: &[
        ParamSpec {
            name: "kappa_w_per_m_k",
            ty: ParamType::Float,
            default: None,
            range: Some((1e-12, 1e5)),
        },
        ParamSpec {
            name: "rho_cp_j_per_m3_k",
            ty: ParamType::Float,
            default: None,
            range: Some((1e-12, 1e12)),
        },
    ],
};

/// SOLV-1's conserved-flux operator (the Euler subset, station 1). γ is the
/// gamma-law EOS parameter — pure data, the degenerate FND-7 spine occupant
/// until the general convex EOS lands (SOLV-1 §3.2).
///
/// Validity range (refusal bounds, sourced): γ ∈ [1.001, 1.667] — the
/// isothermal limit γ → 1 is a singular EOS (c² → isothermal, the gamma-law
/// energy relation degenerates) and no classical ideal gas exceeds the
/// monatomic 5/3; Sod's diatomic 1.4 sits mid-range. Chaotic class: the
/// resolved Euler subset at station-1/2 scales is non-chaotic (no
/// turbulence-resolving content); the LES-resolved tier re-declares this
/// per-regime when it lands (SOLV-1 §3.4, O21).
pub static FLOW_MANIFEST: Manifest = Manifest {
    id: "flow",
    interface_version: InterfaceVersion { major: 1, minor: 0 },
    tables: &[],
    couplers: &[],
    ports: &[],
    chaotic_class: &[
        (Regime::Steady, ChaoticClass::NonChaotic),
        (Regime::Transient, ChaoticClass::NonChaotic),
    ],
    params: &[ParamSpec {
        name: "gamma",
        ty: ParamType::Float,
        default: None,
        range: Some((1.001, 1.667)),
    }],
};

/// SOLV-1 §3.5's one wall-function heat law (station 4). Transport
/// constants are the degenerate FND-7 spine occupant (the gamma-law
/// pattern); the spine supplies per-cell transport when OFFL-5 lands.
///
/// Validity ranges (refusal bounds, sourced): c_p spans heavy combustion
/// products (~500) through hydrogen (~14 300 J/(kg·K)) with margin ⇒
/// [50, 1e5]; μ spans cold rarefied gas (~1e-6) through hot dense gas
/// (~2e-4 Pa·s) with margin ⇒ [1e-7, 1e-2] (a liquid-viscosity μ here
/// means a mis-set deck — the law is a gas-side closure); Pr for gases
/// clusters in [0.2, 1] — bounds [0.05, 5] with margin (a Pr outside that
/// is not a gas). Chaotic class: the closure is algebraic-local,
/// non-chaotic in both regimes.
pub static WALL_HEAT_MANIFEST: Manifest = Manifest {
    id: "wall_heat",
    interface_version: InterfaceVersion { major: 1, minor: 0 },
    tables: &[],
    couplers: &[],
    ports: &[],
    chaotic_class: &[
        (Regime::Steady, ChaoticClass::NonChaotic),
        (Regime::Transient, ChaoticClass::NonChaotic),
    ],
    params: &[
        ParamSpec {
            name: "cp_j_per_kg_k",
            ty: ParamType::Float,
            default: None,
            range: Some((50.0, 1e5)),
        },
        ParamSpec {
            name: "mu_pa_s",
            ty: ParamType::Float,
            default: None,
            range: Some((1e-7, 1e-2)),
        },
        ParamSpec {
            name: "pr",
            ty: ParamType::Float,
            default: None,
            range: Some((0.05, 5.0)),
        },
    ],
};

/// The production registry: one row per landed operator, sorted by id
/// (COUP-8 §3.2 — the only core touch a new mechanism makes).
static MECHANISMS: [&Manifest; 3] = [&CONDUCTION_MANIFEST, &FLOW_MANIFEST, &WALL_HEAT_MANIFEST];

pub fn registry() -> Registry {
    Registry::new(&MECHANISMS, &[])
}

/// Everything a conduction run needs from config: the built grid and the
/// operator parameters, plus the instance name it came from. Parameters are
/// unit-typed at this boundary (META-2 §4 ★); the operator extracts SI
/// `f64` at construction ([`crucible_units::si`]) — kernels stay plain.
#[derive(Debug)]
pub struct ConductionSetup {
    pub grid: Grid,
    pub instance: String,
    pub kappa: crucible_units::ThermalConductivity,
    pub rho_cp: crucible_units::VolumetricHeatCapacity,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SetupError {
    Missing(String),
    /// Multiple instances of one mechanism type where the runner supports
    /// one: refused with all names listed — never a silent first-pick
    /// (META-1 P6; review finding).
    Ambiguous {
        type_id: String,
        instances: Vec<String>,
    },
    /// A loader-blessed value does not fit the grid's index representation
    /// (defense-in-depth behind the loader's MAX_* bounds).
    OutOfRange(String),
    Grid(String),
}

impl std::fmt::Display for SetupError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Missing(m) => write!(f, "config incomplete for this mechanism: {m}"),
            Self::Ambiguous { type_id, instances } => write!(
                f,
                "{} instances of type {type_id:?} ({}) — this runner drives exactly one; \
                 refusing rather than silently picking",
                instances.len(),
                instances.join(", ")
            ),
            Self::OutOfRange(m) => write!(f, "value does not fit the grid representation: {m}"),
            Self::Grid(m) => write!(f, "grid construction failed: {m}"),
        }
    }
}

impl std::error::Error for SetupError {}

impl From<GridError> for SetupError {
    fn from(e: GridError) -> Self {
        Self::Grid(e.to_string())
    }
}

/// Build a `GridSpec` from resolved geometry — shared by every mechanism's
/// setup path (checked narrowing throughout).
pub fn grid_spec_from(geo: &ResolvedGeometry) -> Result<GridSpec, SetupError> {
    let ext = geo
        .extents
        .as_ref()
        .ok_or_else(|| SetupError::Missing("[geometry] grid extents".into()))?;
    let n_r = usize::try_from(ext.n_r)
        .map_err(|_| SetupError::OutOfRange(format!("n_r = {}", ext.n_r)))?;
    let n_z = usize::try_from(ext.n_z)
        .map_err(|_| SetupError::OutOfRange(format!("n_z = {}", ext.n_z)))?;
    let n_theta_max = u32::try_from(geo.n_theta_max)
        .map_err(|_| SetupError::OutOfRange(format!("n_theta_max = {}", geo.n_theta_max)))?;
    // Through the typed extents views (META-2 §4 ★): the seam names the
    // dimension, the kernel-facing GridSpec receives documented-SI f64.
    Ok(GridSpec {
        r_min: crucible_units::si(ext.r_min_length()),
        dr: crucible_units::si(ext.dr_length()),
        n_r,
        z_min: crucible_units::si(ext.z_min_length()),
        dz: crucible_units::si(ext.dz_length()),
        n_z,
        n_theta_max,
        axisymmetry_assertion: geo.axisymmetric,
    })
}

/// The exactly-one-instance seam shared by every single-instance runner:
/// zero instances or several are both refusals (never a silent pick).
fn sole_instance<'a>(
    loaded: &'a Loaded,
    type_id: &str,
) -> Result<(&'a str, &'a crucible_config::toml::Table), SetupError> {
    let candidates = loaded.resolved.mechanisms_of_type(type_id);
    match candidates.as_slice() {
        [] => Err(SetupError::Missing(format!(
            "a mechanism with type = {type_id:?}"
        ))),
        [one] => Ok(*one),
        many => Err(SetupError::Ambiguous {
            type_id: type_id.to_string(),
            instances: many.iter().map(|(n, _)| (*n).to_string()).collect(),
        }),
    }
}

fn resolved_geometry(loaded: &Loaded) -> Result<&ResolvedGeometry, SetupError> {
    loaded
        .resolved
        .geometry
        .as_ref()
        .ok_or_else(|| SetupError::Missing("[geometry] block".into()))
}

/// Build the grid and extract conduction parameters from a loaded config.
/// `fields` is the field set to register (config-time, FND-2 §3.6).
pub fn from_loaded(loaded: &Loaded, fields: &[&str]) -> Result<ConductionSetup, SetupError> {
    let geo = resolved_geometry(loaded)?;
    let (instance, block) = sole_instance(loaded, "conduction")?;
    let missing = |name: &str| SetupError::Missing(format!("{instance}.{name}"));
    Ok(ConductionSetup {
        grid: Grid::build(grid_spec_from(geo)?, fields)?,
        instance: instance.to_string(),
        kappa: ResolvedConfig::param_thermal_conductivity_w_per_m_k(block, "kappa_w_per_m_k")
            .ok_or_else(|| missing("kappa_w_per_m_k"))?,
        rho_cp: ResolvedConfig::param_volumetric_heat_capacity_j_per_m3_k(
            block,
            "rho_cp_j_per_m3_k",
        )
        .ok_or_else(|| missing("rho_cp_j_per_m3_k"))?,
    })
}

/// Everything a flow run needs from config: the built grid (with the `U`
/// component fields registered) and the EOS parameter.
#[derive(Debug)]
pub struct FlowSetup {
    pub grid: Grid,
    pub instance: String,
    pub gamma: f64,
}

/// Build the grid and extract flow parameters from a loaded config — the
/// same typed seam as `from_loaded` (no hand-scanning, no silent picks).
pub fn flow_from_loaded(loaded: &Loaded, fields: &[&str]) -> Result<FlowSetup, SetupError> {
    let geo = resolved_geometry(loaded)?;
    let (instance, block) = sole_instance(loaded, "flow")?;
    let gamma = ResolvedConfig::param_f64(block, "gamma")
        .ok_or_else(|| SetupError::Missing(format!("{instance}.gamma")))?;
    Ok(FlowSetup {
        grid: Grid::build(grid_spec_from(geo)?, fields)?,
        instance: instance.to_string(),
        gamma,
    })
}

/// The `wall_heat` law from config (station 4): the sole-instance seam,
/// typed accessors (META-2 §4 ★ — cp and μ dimensioned, Pr dimensionless),
/// and the law's own operand validation on top of the registry ranges. No
/// grid here — the law is algebraic; the coupler owns the geometry.
pub fn wall_law_from_loaded(
    loaded: &Loaded,
) -> Result<(String, crate::wall_heat::WallLaw), SetupError> {
    let (instance, block) = sole_instance(loaded, "wall_heat")?;
    let missing = |name: &str| SetupError::Missing(format!("{instance}.{name}"));
    let cp = ResolvedConfig::param_specific_heat_capacity_j_per_kg_k(block, "cp_j_per_kg_k")
        .ok_or_else(|| missing("cp_j_per_kg_k"))?;
    let mu = ResolvedConfig::param_dynamic_viscosity_pa_s(block, "mu_pa_s")
        .ok_or_else(|| missing("mu_pa_s"))?;
    let pr = ResolvedConfig::param_f64(block, "pr").ok_or_else(|| missing("pr"))?;
    let law = crate::wall_heat::WallLaw::new(cp, mu, pr)
        .map_err(|e| SetupError::Missing(format!("{instance}: {e}")))?;
    Ok((instance.to_string(), law))
}
