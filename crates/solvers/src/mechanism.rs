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

/// The production registry, v0: one mechanism. Grows one row per landed
/// operator (COUP-8 §3.2 — the only core touch a new mechanism makes).
static MECHANISMS: [&Manifest; 1] = [&CONDUCTION_MANIFEST];

pub fn registry() -> Registry {
    Registry::new(&MECHANISMS, &[])
}

/// Everything a conduction run needs from config: the built grid and the
/// operator parameters, plus the instance name it came from.
#[derive(Debug)]
pub struct ConductionSetup {
    pub grid: Grid,
    pub instance: String,
    pub kappa: f64,
    pub rho_cp: f64,
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
            Self::Missing(m) => write!(f, "config incomplete for a conduction run: {m}"),
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
    Ok(GridSpec {
        r_min: ext.r_min,
        dr: ext.dr,
        n_r,
        z_min: ext.z_min,
        dz: ext.dz,
        n_z,
        n_theta_max,
        axisymmetry_assertion: geo.axisymmetric,
    })
}

/// Build the grid and extract conduction parameters from a loaded config.
/// `fields` is the field set to register (config-time, FND-2 §3.6).
pub fn from_loaded(loaded: &Loaded, fields: &[&str]) -> Result<ConductionSetup, SetupError> {
    let geo = loaded
        .resolved
        .geometry
        .as_ref()
        .ok_or_else(|| SetupError::Missing("[geometry] block".into()))?;
    let candidates = loaded.resolved.mechanisms_of_type("conduction");
    let (instance, block) = match candidates.as_slice() {
        [] => {
            return Err(SetupError::Missing(
                "a mechanism with type = \"conduction\"".into(),
            ));
        }
        [one] => *one,
        many => {
            return Err(SetupError::Ambiguous {
                type_id: "conduction".into(),
                instances: many.iter().map(|(n, _)| (*n).to_string()).collect(),
            });
        }
    };

    let param = |name: &str| -> Result<f64, SetupError> {
        ResolvedConfig::param_f64(block, name)
            .ok_or_else(|| SetupError::Missing(format!("{instance}.{name}")))
    };
    Ok(ConductionSetup {
        grid: Grid::build(grid_spec_from(geo)?, fields)?,
        instance: instance.to_string(),
        kappa: param("kappa_w_per_m_k")?,
        rho_cp: param("rho_cp_j_per_m3_k")?,
    })
}
