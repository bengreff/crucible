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
    /// Declared η_c\* (labeling datum; 1.0 = full equilibrium).
    pub eta_cstar: f64,
    /// The S18 source-level knockdown realizing it (≤ 0; TableEos::h_offset).
    pub h_offset_j_per_kg: f64,
}

/// COUP-7 §3.2 closed-mode turbopump map (see TURBOPUMP_EXPANDER_MANIFEST
/// for the parameter semantics and citations).
#[derive(Debug, Clone)]
pub struct TurbopumpSpec {
    pub mdot_design_kg_per_s: f64,
    pub pump_power_design_w: f64,
    pub impedance_exponent: f64,
    pub turbine_mdot_frac: f64,
    pub turbine_eta: f64,
    pub turbine_pressure_ratio: f64,
    pub turbine_gamma: f64,
    pub coolant_cp_j_per_kg_k: f64,
    pub coolant_t_in_k: f64,
    pub mdot_envelope_lo_frac: f64,
    pub mdot_envelope_hi_frac: f64,
}

#[derive(Debug, Clone)]
pub struct JacketSpec {
    pub h_w_per_m2_k: f64,
    pub t_coolant_k: f64,
}

/// SOLV-4 §3.6 burn-progress combustion selection (S7): the two extra table
/// pins the blend needs (the burnt branch is `table_pin` — the same
/// equilibrium surface shifting mode runs on) + the declared wrinkling.
#[derive(Debug, Clone)]
pub struct BlendSpec {
    /// The `chem_unburnt` pin: (file, group, data_version, digest).
    pub unburnt_pin: (String, String, String, String),
    /// The `chem_ignition` pin.
    pub ignition_pin: (String, String, String, String),
    /// Declared `S_T/S_L` (META-3 `turbulent-flame-speed`; 1.0 = laminar).
    pub wrinkling: f64,
}

/// COUP-7 §3.3 spark-igniter placement + schedule (S7; the one cited datum
/// is `energy_j` — see SPARK_IGNITER_MANIFEST).
#[derive(Debug, Clone)]
pub struct IgniterSpec {
    pub energy_j: f64,
    pub r_m: f64,
    pub z_m: f64,
    pub half_width_m: f64,
    pub window_start_s: f64,
    pub window_s: f64,
}

/// SOLV-6 v1 margin-check operands (S7; see STRUCTURAL_MARGINS_MANIFEST).
#[derive(Debug, Clone)]
pub struct MarginsSpec {
    pub t_real_m: f64,
    pub r_shell_m: f64,
    pub p_coolant_pa: f64,
    pub e_pa: f64,
    pub alpha_per_k: f64,
    pub nu: f64,
    pub yield_cold_pa: f64,
    pub yield_hot_pa: f64,
    pub uts_cold_pa: f64,
    pub uts_hot_pa: f64,
    pub t_cold_k: f64,
    pub t_hot_k: f64,
    pub t_solidus_k: f64,
}

/// COUP-4 §3.1 Stage-1 WORKS criterion (S7), armed by commanded-profile
/// quantities in `[operating_profile]`; `flowthroughs` is then the declared
/// `T_S1_HORIZON`.
#[derive(Debug, Clone)]
pub struct VerdictSpec {
    pub commanded_p_c_pa: Option<f64>,
    pub commanded_thrust_n: Option<f64>,
    pub eps_works: f64,
    pub t_dwell_flowthroughs: f64,
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
/// FND-7 §3.3 — which occupant fills the spine's transport slot, chosen by
/// the config's mechanism id (`transport_constant` vs `transport_table` —
/// the `flow` ↔ `flow_shifting` pattern; Rule 13, data not code).
#[derive(Debug, Clone)]
pub enum TransportSpec {
    /// Fully built here: the declared-constant occupant needs no table.
    Constant(crucible_solvers::transport::ConstantTransport),
    /// The OFFL-5 §3.1a surface pin (file, group, data_version, digest) and
    /// the declared Schmidt number; opened and bound at run time, like the
    /// equilibrium surface.
    Table {
        pin: (String, String, String, String),
        schmidt: f64,
    },
}

pub struct EngineSpec {
    pub grid: Grid,
    pub fields: EulerFields,
    pub t_solid: crucible_grid::FieldId,
    pub rate_solid: crucible_grid::FieldId,
    pub contour: Contour,
    pub injector: InjectorSpec,
    /// `Some` ⇔ a `turbopump_expander` mechanism is selected: the closed
    /// (cycle-coupled) mode — delivered ṁ from the COUP-3 §3.5 fixed
    /// point. `None` = open (spec-driven) mode: the injector's declared ṁ.
    pub turbopump: Option<TurbopumpSpec>,
    /// `None` ⇔ `liner_thickness_m = 0` (cold-flow geometry, no wall
    /// exchange, no conduction sweep).
    pub jacket: Option<JacketSpec>,
    pub liner: Option<LinerSpec>,
    pub wall_law: Option<WallLaw>,
    pub flowthroughs: f64,
    pub cfl: f64,
    pub fill_p_pa: f64,
    pub pumpdown_flowthroughs: f64,
    /// Declared altitude-cell ambient floor [Pa] (see the FND-4 default's
    /// rationale; validated ≥ the pinned table's p envelope at run start).
    pub p_amb_floor_pa: f64,
    /// Injector start-ramp window (flow-through times; 0 = step start).
    pub injector_ramp_flowthroughs: f64,
    /// The `chem_equilibrium` pin: (file, group, data_version, digest).
    pub table_pin: (String, String, String, String),
    /// The FND-7 transport spine occupant (S4).
    pub transport: TransportSpec,
    /// Whether SOLV-1 §3.1's `F_visc` is scheduled — `gas_diffusion`
    /// selected in config. The operator itself is parameter-free; it is
    /// built inside `run` because it borrows its boundary closures.
    pub gas_diffusion: bool,
    /// `Some` ⇔ `combustion_blend` selected: the EOS occupant is the
    /// SOLV-4 §3.6 blend and the burn-progress source + class-R reaction
    /// are scheduled (S7).
    pub blend: Option<BlendSpec>,
    /// `Some` ⇔ `spark_igniter` selected (requires `blend`).
    pub igniter: Option<IgniterSpec>,
    /// `Some` ⇔ `structural_margins` selected (requires the cooled wall).
    pub margins: Option<MarginsSpec>,
    /// `Some` ⇔ the profile declares commanded quantities (COUP-4 S7).
    pub verdict: Option<VerdictSpec>,
    /// COUP-7 §3.2.2 compression declaration (manifest-recorded; None =
    /// no compression declared).
    pub valve_cited_timeline_s: Option<f64>,
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
            eta_cstar: ResolvedConfig::param_f64(block, "eta_cstar").expect("defaulted"),
            h_offset_j_per_kg: ResolvedConfig::param_f64(block, "h_offset_j_per_kg")
                .expect("defaulted"),
        }
    };
    let _ = sole(r, "flow_shifting")?; // the operator selection (no params)
    let turbopump = match r.mechanisms_of_type("turbopump_expander").as_slice() {
        [] => None,
        [(_, tb)] => {
            let f = |name: &str| ResolvedConfig::param_f64(tb, name).expect("declared");
            Some(TurbopumpSpec {
                mdot_design_kg_per_s: f("mdot_design_kg_per_s"),
                pump_power_design_w: f("pump_power_design_w"),
                impedance_exponent: f("impedance_exponent"),
                turbine_mdot_frac: f("turbine_mdot_frac"),
                turbine_eta: f("turbine_eta"),
                turbine_pressure_ratio: f("turbine_pressure_ratio"),
                turbine_gamma: f("turbine_gamma"),
                coolant_cp_j_per_kg_k: f("coolant_cp_j_per_kg_k"),
                coolant_t_in_k: f("coolant_t_in_k"),
                mdot_envelope_lo_frac: f("mdot_envelope_lo_frac"),
                mdot_envelope_hi_frac: f("mdot_envelope_hi_frac"),
            })
        }
        many => {
            return Err(format!(
                "{} `turbopump_expander` instances; this wave assembles at most one",
                many.len()
            ));
        }
    };
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
    if !cooled {
        // Declared cooling machinery with a zero-thickness liner is a
        // config contradiction — refuse rather than silently running an
        // adiabatic engine with cooling selected (session-12 review).
        for ty in ["jacket_coolant", "conduction", "wall_heat"] {
            if !r.mechanisms_of_type(ty).is_empty() {
                return Err(format!(
                    "`{ty}` is selected but liner_thickness_m = 0 (no solid ring): the \
                     declared cooling path cannot exist — set a liner thickness or drop \
                     the cooling mechanisms"
                ));
            }
        }
    }
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
        let (_, _wb) = sole(r, "wall_heat")?; // parameter-free since S4
        (Some(jacket), Some(liner), Some(WallLaw::new()))
    } else {
        (None, None, None)
    };

    // --- FND-7 §3.3: which occupant fills the spine's transport slot ------
    // Exactly one, always: every configuration has a medium, and leaving it
    // unstated would mean a hidden default (FND-4 §3.5).
    let n_const = r.mechanisms_of_type("transport_constant").len();
    let n_table = r.mechanisms_of_type("transport_table").len();
    let transport = match (n_const, n_table) {
        (1, 0) => {
            let (name, law) = crucible_solvers::constant_transport_from_loaded(loaded)
                .map_err(|e| e.to_string())?;
            let _ = name;
            TransportSpec::Constant(law)
        }
        (0, 1) => {
            let (_, schmidt) = crucible_solvers::table_transport_schmidt_from_loaded(loaded)
                .map_err(|e| e.to_string())?;
            let pin = r.tables.get("spine_transport").ok_or(
                "transport_table requires the `spine_transport` table pin ([tables] block)",
            )?;
            TransportSpec::Table {
                pin: (
                    pin.file.clone(),
                    pin.group.clone(),
                    pin.data_version.clone(),
                    pin.content_digest.clone(),
                ),
                schmidt,
            }
        }
        (0, 0) => {
            return Err("no transport occupant selected: declare exactly one of \
                 `transport_constant` (declared constants) or `transport_table` \
                 (the OFFL-5 §3.1a spine surface). The wall law and F_visc both \
                 read it, and neither may invent one (FND-7 §3.3)"
                .to_string());
        }
        (c, t) => {
            return Err(format!(
                "{c} `transport_constant` + {t} `transport_table` instances: the spine \
                 has exactly ONE transport slot; refusing rather than picking"
            ));
        }
    };
    // SOLV-1 §3.1's F_visc, selected as data (Rule 13).
    let gas_diffusion = match r.mechanisms_of_type("gas_diffusion").as_slice() {
        [] => false,
        [_] => true,
        many => {
            return Err(format!(
                "{} `gas_diffusion` instances; there is one gas",
                many.len()
            ));
        }
    };

    // --- SOLV-4 §3.6 combustion selection (S7) -----------------------------
    let blend = match r.mechanisms_of_type("combustion_blend").as_slice() {
        [] => None,
        [(_, cb)] => {
            let get_pin = |name: &str| -> Result<(String, String, String, String), String> {
                let pin = r.tables.get(name).ok_or(format!(
                    "combustion_blend requires the `{name}` table pin ([tables] block)"
                ))?;
                Ok((
                    pin.file.clone(),
                    pin.group.clone(),
                    pin.data_version.clone(),
                    pin.content_digest.clone(),
                ))
            };
            Some(BlendSpec {
                unburnt_pin: get_pin("chem_unburnt")?,
                ignition_pin: get_pin("chem_ignition")?,
                wrinkling: ResolvedConfig::param_f64(cb, "wrinkling").expect("declared"),
            })
        }
        many => {
            return Err(format!(
                "{} `combustion_blend` instances; there is one burn-progress field",
                many.len()
            ));
        }
    };
    let igniter = match r.mechanisms_of_type("spark_igniter").as_slice() {
        [] => None,
        [(_, ib)] => {
            if blend.is_none() {
                return Err(
                    "spark_igniter without combustion_blend: a spark with no burn-progress \
                     field is an unconsumed deposit — select combustion_blend (with its \
                     chem_unburnt/chem_ignition pins) or drop the igniter"
                        .to_string(),
                );
            }
            let f = |name: &str| ResolvedConfig::param_f64(ib, name).expect("declared");
            Some(IgniterSpec {
                energy_j: f("energy_j"),
                r_m: f("r_m"),
                z_m: f("z_m"),
                half_width_m: f("half_width_m"),
                window_start_s: f("window_start_s"),
                window_s: f("window_s"),
            })
        }
        many => {
            return Err(format!(
                "{} `spark_igniter` instances; this wave assembles at most one",
                many.len()
            ));
        }
    };
    let margins = match r.mechanisms_of_type("structural_margins").as_slice() {
        [] => None,
        [(_, mb)] => {
            if !cooled {
                return Err(
                    "structural_margins without a liner (liner_thickness_m = 0): there is \
                     no annotated shell component to check — set a liner or drop the margins"
                        .to_string(),
                );
            }
            let f = |name: &str| ResolvedConfig::param_f64(mb, name).expect("declared");
            let m = MarginsSpec {
                t_real_m: f("t_real_m"),
                r_shell_m: f("r_shell_m"),
                p_coolant_pa: f("p_coolant_pa"),
                e_pa: f("e_pa"),
                alpha_per_k: f("alpha_per_k"),
                nu: f("nu"),
                yield_cold_pa: f("yield_cold_pa"),
                yield_hot_pa: f("yield_hot_pa"),
                uts_cold_pa: f("uts_cold_pa"),
                uts_hot_pa: f("uts_hot_pa"),
                t_cold_k: f("t_cold_k"),
                t_hot_k: f("t_hot_k"),
                t_solidus_k: f("t_solidus_k"),
            };
            if m.t_cold_k >= m.t_hot_k {
                return Err(
                    "structural_margins: t_cold_k must be below t_hot_k (the two cited \
                     allowable temperatures)"
                        .to_string(),
                );
            }
            if 2.0 * m.r_shell_m / m.t_real_m <= 20.0 {
                return Err(format!(
                    "structural_margins: the declared shell (R = {} m, t = {} m) has \
                     2R/t = {:.1} <= 20 — under-idealized for the v1 thin-shell subset \
                     (SOLV-6 §5): reported, never smeared; declare the cited shell \
                     radius of the pressure-carrying member or drop the margins",
                    m.r_shell_m,
                    m.t_real_m,
                    2.0 * m.r_shell_m / m.t_real_m
                ));
            }
            Some(m)
        }
        many => {
            return Err(format!(
                "{} `structural_margins` instances; the liner is one component",
                many.len()
            ));
        }
    };

    if turbopump.is_some() && wall_law.is_none() {
        return Err(
            "turbopump_expander (closed mode) needs the cooled wall path — the jacket \
             enthalpy rise IS the drive power (COUP-7 §3.2); declare liner/wall/jacket \
             mechanisms and a liner_thickness_m > 0"
                .to_string(),
        );
    }
    let profile = r
        .operating_profile
        .as_ref()
        .ok_or("engine assembly needs an [operating_profile] (mode = \"steady_march\")")?;
    // COUP-4 §3.1 (S7): commanded quantities arm the WORKS criterion; the
    // loader materialized the criterion constants iff armed.
    let verdict = if profile.commanded_p_c_pa.is_some() || profile.commanded_thrust_n.is_some() {
        Some(VerdictSpec {
            commanded_p_c_pa: profile.commanded_p_c_pa,
            commanded_thrust_n: profile.commanded_thrust_n,
            eps_works: profile.eps_works.expect("materialized when armed"),
            t_dwell_flowthroughs: profile
                .t_dwell_flowthroughs
                .expect("materialized when armed"),
        })
    } else {
        None
    };
    if let Some(v) = &verdict
        && v.t_dwell_flowthroughs >= profile.flowthroughs
    {
        return Err(format!(
            "the dwell window ({} flow-throughs) does not fit inside the Stage-1 \
             horizon (flowthroughs = {}): the WORKS criterion can never complete — \
             raise flowthroughs or shrink t_dwell_flowthroughs (COUP-4 §3.1)",
            v.t_dwell_flowthroughs, profile.flowthroughs
        ));
    }
    if turbopump.is_some() {
        // The closed loop engages after the establishment window (pump-down
        // + start ramp); a budget that ends inside it would silently run
        // open-mode — refuse the contradiction (session-12 review).
        // Sequential establishment: ramp first, then pump-down (run.rs).
        let window = profile.pumpdown_flowthroughs + profile.injector_ramp_flowthroughs;
        if profile.flowthroughs <= window {
            return Err(format!(
                "closed mode (turbopump_expander) never engages: flowthroughs = {} does \
                 not exceed the establishment window ({} flow-throughs of pump-down/ramp); \
                 raise flowthroughs or drop the turbopump",
                profile.flowthroughs, window
            ));
        }
    }

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
        turbopump,
        jacket,
        liner,
        wall_law,
        flowthroughs: profile.flowthroughs,
        cfl: profile.cfl,
        fill_p_pa: profile.fill_p_pa,
        pumpdown_flowthroughs: profile.pumpdown_flowthroughs,
        p_amb_floor_pa: profile.p_amb_floor_pa,
        injector_ramp_flowthroughs: profile.injector_ramp_flowthroughs,
        table_pin,
        transport,
        gas_diffusion,
        blend,
        igniter,
        margins,
        verdict,
        valve_cited_timeline_s: profile.valve_cited_timeline_s,
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
