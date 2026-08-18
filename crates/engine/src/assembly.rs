//! Config → assembled engine: extract the validated mechanism/boundary
//! parameters, verify the content-addressed geometry, build the ternary
//! grid, and hand a ready [`EngineSpec`] to the runner. Everything here is
//! *reading resolved data* — the loader already validated types, ranges
//! (validity envelopes), and O20 bindings; a `None` from a param accessor
//! here means an undeclared name, i.e. a programming error surfaced loudly.

use crate::geometry::Contour;
use crucible_config::{Loaded, ResolvedConfig, parse_contour_csv};
use crucible_grid::{CellGeom, Grid, GridSpec, Region};
use crucible_solvers::euler::{EULER_FIELDS, EulerFields};
use crucible_solvers::wall_heat::WallLaw;
use crucible_units::{dynamic_viscosity_pa_s, specific_heat_capacity_j_per_kg_k};

/// Grid field names: the six `U` components + the liner temperature/rate.
pub const T_SOLID: &str = "t_solid";
pub const RATE_SOLID: &str = "rate_solid";

#[derive(Debug, Clone)]
pub struct InjectorSpec {
    pub mdot_kg_per_s: f64,
    pub mixture_ratio: f64,
    /// Z = 1/(1 + MR) — the elemental mixture fraction the surface is keyed
    /// on (OFFL-3 §3.3).
    pub z_frac: f64,
    pub h_inj_j_per_kg: f64,
}

#[derive(Debug, Clone)]
pub struct JacketSpec {
    pub h_w_per_m2_k: f64,
    pub t_coolant_k: f64,
}

#[derive(Debug, Clone)]
pub struct LinerSpec {
    pub kappa_w_per_m_k: f64,
    /// Declared steady-state continuation device when set far below the
    /// physical value (station-4 finding: the steady solution is
    /// independent of ρc_p).
    pub rho_cp_j_per_m3_k: f64,
}

/// Everything a run needs besides the table (which the caller owns — the
/// `TableEos` borrows it, so the binding happens in the runner's scope).
pub struct EngineSpec {
    pub grid: Grid,
    pub fields: EulerFields,
    pub t_solid: crucible_grid::FieldId,
    pub rate_solid: crucible_grid::FieldId,
    pub contour: Contour,
    pub injector: InjectorSpec,
    /// `None` ⇔ `liner_thickness_m = 0` (cold-flow geometry, no wall
    /// exchange, no conduction sweep).
    pub jacket: Option<JacketSpec>,
    pub liner: Option<LinerSpec>,
    pub wall_law: Option<WallLaw>,
    pub flowthroughs: f64,
    pub cfl: f64,
    pub fill_p_pa: f64,
    pub pumpdown_flowthroughs: f64,
    /// The `chem_equilibrium` pin: (file, group, data_version, digest).
    pub table_pin: (String, String, String, String),
}

/// Assemble from a loaded config. `read_file` supplies the contour CSV
/// content (same reader discipline as the loader — file access stays with
/// the caller); the content is re-verified against the recorded digest, so
/// geometry drift between load and assembly refuses.
pub fn assemble(
    loaded: &Loaded,
    read_file: &dyn Fn(&str) -> Result<String, String>,
) -> Result<EngineSpec, String> {
    let r = &loaded.resolved;
    let geo = r.geometry.as_ref().ok_or("config declares no [geometry]")?;
    let rc = geo
        .contour
        .as_ref()
        .ok_or("engine assembly needs a contour geometry (a box world is a station fixture)")?;
    let ext = geo
        .extents
        .as_ref()
        .ok_or("resolved contour config must carry derived extents")?;
    if !geo.axisymmetric || geo.n_theta_max != 1 {
        return Err(
            "this wave runs under the recorded axisymmetry assertion (n_theta_max = 1); \
             N_θ > 1 assemblies arrive with the refluxing wave (COUP-2/3)"
                .to_string(),
        );
    }

    // Content-addressed geometry: re-read and re-verify.
    let content = read_file(&rc.contour)?;
    let digest = format!("sha256:{}", crucible_config::sha256_hex(content.as_bytes()));
    if digest != rc.contour_digest {
        return Err(format!(
            "contour {} hashes to {digest} but the resolved config pinned {} — \
             geometry drifted between load and assembly; refusing (FND-5 §3.2 doctrine)",
            rc.contour, rc.contour_digest
        ));
    }
    let scale = if rc.contour_units == "in" {
        crucible_config::INCH_M
    } else {
        1.0
    };
    let stations = parse_contour_csv(&content, scale)?;
    let contour = Contour::new(stations, rc.liner_thickness_m)?;

    // Mechanism/boundary parameter extraction (loader-validated).
    let injector = {
        let (_, block) = sole(r, "injector_prior")?;
        let mdot = ResolvedConfig::param_f64(block, "mdot_kg_per_s").expect("declared");
        let mr = ResolvedConfig::param_f64(block, "mixture_ratio").expect("declared");
        let h = ResolvedConfig::param_f64(block, "h_inj_j_per_kg").expect("declared");
        InjectorSpec {
            mdot_kg_per_s: mdot,
            mixture_ratio: mr,
            z_frac: 1.0 / (1.0 + mr),
            h_inj_j_per_kg: h,
        }
    };
    let _ = sole(r, "flow_shifting")?; // the operator selection (no params)
    let table_pin = {
        let pin = r
            .tables
            .get("chem_equilibrium")
            .ok_or("flow_shifting requires the `chem_equilibrium` table pin ([tables] block)")?;
        (
            pin.file.clone(),
            pin.group.clone(),
            pin.data_version.clone(),
            pin.content_digest.clone(),
        )
    };

    let cooled = rc.liner_thickness_m > 0.0;
    let (jacket, liner, wall_law) = if cooled {
        let (_, jb) = sole(r, "jacket_coolant")?;
        let jacket = JacketSpec {
            h_w_per_m2_k: ResolvedConfig::param_f64(jb, "h_w_per_m2_k").expect("declared"),
            t_coolant_k: ResolvedConfig::param_f64(jb, "t_coolant_k").expect("declared"),
        };
        let (_, cb) = sole(r, "conduction")?;
        let liner = LinerSpec {
            kappa_w_per_m_k: ResolvedConfig::param_f64(cb, "kappa_w_per_m_k").expect("declared"),
            rho_cp_j_per_m3_k: ResolvedConfig::param_f64(cb, "rho_cp_j_per_m3_k")
                .expect("declared"),
        };
        let (_, wb) = sole(r, "wall_heat")?;
        let law = WallLaw::new(
            specific_heat_capacity_j_per_kg_k(
                ResolvedConfig::param_f64(wb, "cp_j_per_kg_k").expect("declared"),
            ),
            dynamic_viscosity_pa_s(ResolvedConfig::param_f64(wb, "mu_pa_s").expect("declared")),
            ResolvedConfig::param_f64(wb, "pr").expect("declared"),
        )
        .map_err(|e| format!("wall_heat: {e}"))?;
        (Some(jacket), Some(liner), Some(law))
    } else {
        (None, None, None)
    };

    let profile = r
        .operating_profile
        .as_ref()
        .ok_or("engine assembly needs an [operating_profile] (mode = \"steady_march\")")?;

    // The ternary grid through the FND-3 §3.3/§3.6 ingest seam: analytic
    // partial fractions + face apertures from the contour clip. Gas ⇔
    // κ > 0 (the by-center stair classification is retired); the liner
    // ring is the κ = 0 band within the declared thickness of the wall,
    // by center — grid-thickened as before. Face expressions use the
    // face-coordinate arithmetic of `Grid::face_radius`/`z_center` so a
    // shared face computes bitwise-identically from both sides.
    let spec = GridSpec {
        r_min: ext.r_min,
        dr: ext.dr,
        n_r: usize::try_from(ext.n_r).map_err(|_| "n_r overflow")?,
        z_min: ext.z_min,
        dz: ext.dz,
        n_z: usize::try_from(ext.n_z).map_err(|_| "n_z overflow")?,
        n_theta_max: 1,
        axisymmetry_assertion: true,
    };
    let mut names: Vec<&str> = EULER_FIELDS.to_vec();
    names.push(T_SOLID);
    names.push(RATE_SOLID);
    let (r_min, dr, dz, z0) = (spec.r_min, spec.dr, spec.dz, spec.z_min);
    let c = contour.clone();
    let grid = Grid::build_with_geometry(spec, &names, |i_r, i_z| {
        let r0 = r_min + i_r as f64 * dr;
        let r1 = r_min + (i_r + 1) as f64 * dr;
        let za = z0 + i_z as f64 * dz;
        let zb = z0 + (i_z + 1) as f64 * dz;
        let kappa = c.gas_volume_fraction(r0, r1, za, zb);
        if kappa > 0.0 {
            CellGeom {
                region: Region::Gas,
                kappa,
                aperture: [
                    c.r_face_aperture(r0, za, zb),
                    c.r_face_aperture(r1, za, zb),
                    c.z_face_aperture(r0, r1, za),
                    c.z_face_aperture(r0, r1, zb),
                ],
            }
        } else {
            let rc = r_min + (i_r as f64 + 0.5) * dr;
            let zc = z0 + (i_z as f64 + 0.5) * dz;
            let region = if c.liner_thickness_m > 0.0 && rc < c.r_wall(zc) + c.liner_thickness_m {
                Region::Solid
            } else {
                Region::Exterior
            };
            CellGeom {
                region,
                kappa: 0.0,
                aperture: [0.0; 4],
            }
        }
    })
    .map_err(|e| format!("grid build: {e}"))?;
    let fields = EulerFields::resolve(&grid).map_err(|e| format!("fields: {e}"))?;
    let t_solid = grid.field_id(T_SOLID).map_err(|e| format!("{e}"))?;
    let rate_solid = grid.field_id(RATE_SOLID).map_err(|e| format!("{e}"))?;

    Ok(EngineSpec {
        grid,
        fields,
        t_solid,
        rate_solid,
        contour,
        injector,
        jacket,
        liner,
        wall_law,
        flowthroughs: profile.flowthroughs,
        cfl: profile.cfl,
        fill_p_pa: profile.fill_p_pa,
        pumpdown_flowthroughs: profile.pumpdown_flowthroughs,
        table_pin,
    })
}

/// Sole-instance selection through the typed seam (session-6 convention:
/// count before you choose; never "the first of several").
fn sole<'a>(
    r: &'a ResolvedConfig,
    type_id: &str,
) -> Result<(&'a str, &'a crucible_config::toml::Table), String> {
    let of_type = r.mechanisms_of_type(type_id);
    match of_type.len() {
        0 => Err(format!("config selects no `{type_id}` mechanism")),
        1 => Ok(of_type[0]),
        n => Err(format!(
            "{n} `{type_id}` instances; this wave assembles exactly one"
        )),
    }
}
