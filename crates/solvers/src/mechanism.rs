//! The project's first real COUP-8 registry entry (`conduction`) and the
//! config→grid→operator wiring: `from_loaded` turns an FND-4 `Loaded`
//! (resolved config + manifest) into a built grid + operator parameters.
//! Boundary conditions and sources remain caller-supplied this session —
//! their config grammar belongs to the COUP-7 boundary-object wave.

use crucible_config::Loaded;
use crucible_grid::{Grid, GridError, GridSpec};
use crucible_registry::{
    ChaoticClass, InterfaceVersion, Manifest, ParamSpec, ParamType, Regime, Registry,
};

/// SOLV-1 §3.5 diffusion-class conduction: non-chaotic in every declared
/// regime (heat diffusion has no sensitive dependence), so it composes with
/// `relaxed` determinism mode once the GPU path exists (O21).
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
pub struct ConductionSetup {
    pub grid: Grid,
    pub instance: String,
    pub kappa: f64,
    pub rho_cp: f64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SetupError {
    Missing(String),
    Grid(String),
}

impl std::fmt::Display for SetupError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Missing(m) => write!(f, "config incomplete for a conduction run: {m}"),
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

/// Build the grid and extract conduction parameters from a loaded config.
/// `fields` is the field set to register (config-time, FND-2 §3.6).
pub fn from_loaded(loaded: &Loaded, fields: &[&str]) -> Result<ConductionSetup, SetupError> {
    let geo = loaded
        .resolved
        .geometry
        .as_ref()
        .ok_or_else(|| SetupError::Missing("[geometry] block".into()))?;
    let ext = geo
        .extents
        .as_ref()
        .ok_or_else(|| SetupError::Missing("[geometry] grid extents".into()))?;
    let (instance, block) = loaded
        .resolved
        .mechanisms
        .iter()
        .find(|(_, b)| b.get("type").and_then(|v| v.as_str()) == Some("conduction"))
        .ok_or_else(|| SetupError::Missing("a mechanism with type = \"conduction\"".into()))?;

    let param = |name: &str| -> Result<f64, SetupError> {
        block
            .get(name)
            .and_then(toml_float)
            .ok_or_else(|| SetupError::Missing(format!("{instance}.{name}")))
    };
    let spec = GridSpec {
        r_min: ext.r_min,
        dr: ext.dr,
        n_r: ext.n_r as usize,
        z_min: ext.z_min,
        dz: ext.dz,
        n_z: ext.n_z as usize,
        n_theta_max: geo.n_theta_max as u32,
        axisymmetry_assertion: geo.axisymmetric,
    };
    Ok(ConductionSetup {
        grid: Grid::build(spec, fields)?,
        instance: instance.clone(),
        kappa: param("kappa_w_per_m_k")?,
        rho_cp: param("rho_cp_j_per_m3_k")?,
    })
}

fn toml_float(v: &crucible_config::toml::Value) -> Option<f64> {
    v.as_float().or_else(|| v.as_integer().map(|i| i as f64))
}
