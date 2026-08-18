//! The ONE coupled stepper + the SOLV-7 readout, parameterized entirely by
//! the assembled [`EngineSpec`] — no engine-specific code (Rule 13).
//!
//! Step structure (station-4's certified pattern, generalized to contour
//! geometry and the table EOS; honest scaffolding until COUP-3's SDC-IMEX
//! class-D lands): evaluate every gas↔liner wall exchange once from the
//! pre-step state (SOLV-1 §3.5's one law — conservation by construction:
//! the same per-face q drives the solid and debits the gas), advance the
//! liner conduction (coolant Robin on its exterior faces), debit the gas,
//! advance the gas (slip-ghost stair walls about the true contour normal —
//! station-2 machinery). Δt = gas CFL, guarded against both thermal
//! stability limits (fail-loud, never silently sub-stepped).
//!
//! Readout (SOLV-7): every reported number is a plane integral of the
//! conserved `U` — thrust from exit momentum + pressure flux (vacuum,
//! SOLV-7.1), `p_c` = area-averaged stagnation pressure at the injector-end
//! plane (N11 convention), c\*/C_F/Isp per SOLV-7.2–7.4. Emergent, never
//! imposed.

use crate::assembly::EngineSpec;
use crucible_constants::G0;
use crucible_grid::{BRICK, FaceDir, Grid, InterfaceFace};
use crucible_solvers::euler::{
    Cons, EosLaw, Euler, EulerFields, FlowBc, FlowBcs, I_EN, NCOMP, TableEos, srd_neighborhood,
};
use crucible_solvers::wall_heat::NearWallGas;
use crucible_solvers::{Bcs, Conduction, Domain, FaceBc, InteriorFaces};
use crucible_tables::{Pin, Table};

/// One gas-cell wall patch — the SOLV-1 §3.5 exchange surface on the cut
/// geometry: the embedded-interface area |W| and normal from the grid's
/// closure identity (the smooth-wall area, not the stair overcount), the
/// gas↔solid grid faces it spans (the solid-side flux carriers), and the
/// SRD debit set (a small cut cell cannot absorb its own wall debit — the
/// merged control volume that stabilizes its flux update absorbs the
/// exchange too; same neighborhood rule, single owner in the solver).
struct WallPatch {
    gas: (usize, usize),
    faces: Vec<InterfaceFace>,
    /// Full grid area of each spanned face (Σ = the stair area).
    face_areas: Vec<f64>,
    /// |W| — the interface area per full ring (N_θ = 1).
    area: f64,
    /// Outward (gas→wall) unit normal (n_r, n_z).
    n_hat: (f64, f64),
    /// Primary solid partner (across the largest spanned face).
    solid: (usize, usize),
    /// Energy-debit cells (the gas cell's SRD neighborhood; `[self]` for
    /// regular cells) and Σ κV over them.
    debit_cells: Vec<(usize, usize)>,
    debit_kv_sum: f64,
}

/// Enumerate wall patches (deterministic gas-cell lexicographic order).
/// The single-valued r_wall(z) contour class cannot produce a multi-sided
/// wall (slot) inside one cell — the gas region {r < r_wall(z)} puts every
/// cell's wall on one connected outboard arc — so the closure vector is
/// always a faithful single interface here; per-side reconstruction for
/// genuine slots is the FND-3 PLIC/CSG wave. Note the stair-face sum may
/// legitimately exceed |W| by a large factor at steeply-crossing walls
/// (a covered face the wall immediately dives away from carries almost no
/// true interface): that overcount is exactly what |W| corrects.
fn build_wall_patches(g: &Grid) -> Result<Vec<WallPatch>, String> {
    let mut patches: Vec<WallPatch> = Vec::new();
    for face in g.gas_solid_faces() {
        if patches.last().map(|p| p.gas) != Some(face.gas) {
            patches.push(WallPatch {
                gas: face.gas,
                faces: Vec::new(),
                face_areas: Vec::new(),
                area: 0.0,
                n_hat: (0.0, 0.0),
                solid: face.solid,
                debit_cells: Vec::new(),
                debit_kv_sum: 0.0,
            });
        }
        let p = patches.last_mut().expect("just pushed");
        p.face_areas.push(g.interface_area_per_theta(&face, 1));
        p.faces.push(face);
    }
    for p in &mut patches {
        let (i_r, i_z) = p.gas;
        let (w_r, w_z) = g.wall_closure(i_r, i_z, 1);
        let area = (w_r * w_r + w_z * w_z).sqrt();
        if !area.is_finite() || area <= 0.0 {
            return Err(format!(
                "wall cell ({i_r}, {i_z}): zero closure interface area yet gas↔solid \
                 faces exist — geometry incoherent (a slot-class wall? per-side \
                 interface reconstruction is the FND-3 PLIC wave); refusing"
            ));
        }
        p.area = area;
        p.n_hat = (-w_r / area, -w_z / area);
        // Primary partner: across the largest spanned face (first wins ties
        // — the fixed FaceDir enumeration order).
        let mut best = 0usize;
        for (k, a) in p.face_areas.iter().enumerate() {
            if *a > p.face_areas[best] {
                best = k;
            }
        }
        p.solid = p.faces[best].solid;
        let hood = srd_neighborhood(g, i_r, i_z).map_err(|e| format!("wall patch: {e}"))?;
        match hood {
            Some(members) => {
                p.debit_kv_sum = members.iter().map(|(_, kv)| kv).sum();
                p.debit_cells = members.into_iter().map(|(c, _)| c).collect();
            }
            None => {
                p.debit_cells = vec![p.gas];
                p.debit_kv_sum = g.kappa(i_r, i_z) * g.cell_volume(i_r, 1);
            }
        }
    }
    Ok(patches)
}

/// Fraction of the solid/exchange stability limits the coupled Δt may use
/// (station-4 constants, same rationale).
pub const SOLID_DT_FRAC: f64 = 0.5;
pub const EXCHANGE_DT_FRAC: f64 = 0.5;

/// Steadiness probe cadence (steps) for the residual + progress callback.
pub const PROBE_EVERY: usize = 200;

/// Runaway guard: no coarse-to-full run on the declared tiers legitimately
/// exceeds this many steps; hitting it is a mis-sized config, not progress.
pub const MARCH_STEP_CAP: usize = 2_000_000;

#[derive(Debug, Clone)]
pub struct Progress {
    pub step: usize,
    pub t: f64,
    pub t_final: f64,
    /// max |Δρ|/ρ over the probe window (station-2-style steadiness).
    pub resid: f64,
    pub mdot_exit: f64,
    pub thrust_n: f64,
}

/// A mid-march halt with its crash artifact. META-1 P6's "fail loud, halt
/// clean" includes leaving the evidence: the dump is fault-tolerant — raw
/// conserved `U` for every active cell (the conserved state is the ground
/// truth at a halt), derived state only where the EOS projection still
/// succeeds. Empty `crash_csv` ⇔ the halt fired before a field existed.
#[derive(Debug, Clone)]
pub struct Halt {
    pub message: String,
    pub step: usize,
    pub t: f64,
    pub crash_csv: String,
}

impl std::fmt::Display for Halt {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{} (step {}, t = {:.6e} s)",
            self.message, self.step, self.t
        )
    }
}

/// The SOLV-7 §3.4 performance object (scalar subset — the p-box wrapper
/// arrives with the COUP-5 bracket wave) + run bookkeeping.
#[derive(Debug, Clone)]
pub struct Report {
    pub thrust_n: f64,
    pub isp_s: f64,
    pub v_e_m_per_s: f64,
    pub c_star_m_per_s: f64,
    pub c_f: f64,
    /// Emergent injector-end stagnation chamber pressure (N11).
    pub p_c_pa: f64,
    pub mdot_exit_kg_per_s: f64,
    pub mdot_injected_kg_per_s: f64,
    /// Total gas→liner wall heat at the final state [W] (the expander drive
    /// integrand; 0 for cold-flow geometry).
    pub jacket_watts: f64,
    pub liner_t_max_k: f64,
    pub steps: usize,
    pub t_end: f64,
    pub steady_resid: f64,
    pub active_gas_cells: usize,
    pub throat_area_m2: f64,
    /// r,z,region,rho,u_r,u_z,p,T,Z,M — one row per active cell (viz feed).
    pub fields_csv: String,
}

/// Open the pinned equilibrium table recorded in the spec (FND-5 gate:
/// version + digest enforced at open).
pub fn open_pinned_table(spec: &EngineSpec) -> Result<Table, String> {
    let (file, group, version, digest) = &spec.table_pin;
    Table::open(
        file,
        group,
        &Pin {
            data_version: version.clone(),
            content_digest: Some(digest.clone()),
        },
    )
    .map_err(|e| format!("table {file}:{group}: {e}"))
}

/// March the assembled engine to its settle budget and read out the
/// performance object. `on_progress` fires every [`PROBE_EVERY`] steps.
/// A mid-march failure returns a [`Halt`] carrying the crash artifact.
pub fn run(
    spec: &mut EngineSpec,
    table: &Table,
    on_progress: &mut dyn FnMut(&Progress),
) -> Result<Report, Halt> {
    let pre = |message: String| Halt {
        message,
        step: 0,
        t: 0.0,
        crash_csv: String::new(),
    };
    let eos = TableEos::bind(table).map_err(pre)?;

    // --- Initial fill: quiescent near-vacuum equilibrium gas + cold liner --
    let u_fill: Cons = eos
        .cons_from_phz(
            spec.fill_p_pa,
            spec.injector.h_inj_j_per_kg,
            spec.injector.z_frac,
            [0.0, 0.0, 0.0],
        )
        .map_err(|e| pre(format!("fill state: {e}")))?;
    for (k, &v) in u_fill.iter().enumerate() {
        spec.grid.fill_field(spec.fields.ids()[k], move |_, _, _| v);
    }
    let t_cool = spec.jacket.as_ref().map_or(300.0, |j| j.t_coolant_k);
    fill_solid(&mut spec.grid, spec.t_solid, t_cool);

    // --- Injector mass flux through the DISCRETE inlet plane --------------
    // The boundary object states ṁ; the flux is sized by the stair plane's
    // actual open area so the delivered ṁ is exactly the declared one.
    let a_inlet = plane_area(&spec.grid, 0, FaceDir::ZMinus);
    if a_inlet <= 0.0 {
        return Err(pre("no active gas cells on the injector plane".to_string()));
    }
    let inflow = FlowBc::MassFlowInflow {
        mdot_per_area: spec.injector.mdot_kg_per_s / a_inlet,
        h_total: spec.injector.h_inj_j_per_kg,
        c_frac: spec.injector.z_frac,
    };

    // --- Settle budget: flow-throughs of the fill-state acoustic time -----
    // (a pure function of {config, table}: deterministic, never wall-clock).
    let span = spec.contour.z_max() - spec.contour.z_min();
    let a_ref = {
        let w = eos
            .prim_checked(&u_fill)
            .map_err(|e| pre(format!("fill state sound speed: {e}")))?;
        eos.sound_speed_w(&w)
    };
    let t_final = spec.flowthroughs * span / a_ref;

    let contour = spec.contour.clone();
    let normal_fn = move |r: f64, z: f64| contour.wall_normal(r, z);
    let zero_src: fn(f64, f64, f64, f64) -> Cons = |_, _, _, _| [0.0; NCOMP];
    // Vacuum plume seam with the altitude-cell pump-down: ambient falls
    // log-linearly from the fill pressure to the table-envelope floor
    // (+25% margin) over the declared window, so the nozzle establishes
    // quasi-statically (a violent free drain shocks/starves stair corners);
    // zero upstream influence once the exit runs supersonic.
    let p_floor = 1.25 * eos.envelopes()[0].0;
    let p_fill = spec.fill_p_pa;
    let t_pump = spec.pumpdown_flowthroughs * span / a_ref;
    let pump_schedule = move |t: f64| -> f64 {
        if t_pump <= 0.0 {
            return p_floor;
        }
        let s = (t / t_pump).min(1.0);
        p_fill * (p_floor / p_fill).powf(s)
    };
    let op = Euler {
        eos: eos.clone(),
        source: &zero_src,
        bcs: FlowBcs {
            r_inner: FlowBc::Reflecting, // zero-area axis face
            r_outer: FlowBc::Reflecting, // unreached: liner/exterior interpose
            z_lo: inflow,
            z_hi: FlowBc::PressureOutflow(&pump_schedule),
        },
        wall_normal: Some(&normal_fn),
        // Slip everywhere (the station-2 certified treatment): mirror steps
        // shock the chamber off-surface even under the pump-down; the slip
        // form's residual risk is exit-lip corner starvation at hypersonic
        // Mach (the FND-3 State-Redistribution deferral class), softened by
        // the quasi-static establishment below.
        slip_wall_z_faces: true,
    };

    let patches = if spec.wall_law.is_some() {
        build_wall_patches(&spec.grid).map_err(pre)?
    } else {
        Vec::new()
    };
    let n_gas = count_active(&spec.grid);

    // --- March -------------------------------------------------------------
    let mut t = 0.0f64;
    let mut steps = 0usize;
    let mut resid = f64::INFINITY;
    let mut rho_probe = snapshot_rho(&spec.grid, &spec.fields);
    let mut jacket_watts = 0.0f64;
    while t_final - t > 1e-12 * t_final {
        let dt_cap = t_final - t;
        let dt = match coupled_step(spec, &op, &eos, &patches, t, dt_cap, &mut jacket_watts) {
            Ok(dt) => dt,
            Err(message) => {
                return Err(Halt {
                    crash_csv: crash_fields_csv(
                        &spec.grid,
                        &spec.fields,
                        &eos,
                        spec.t_solid,
                        &message,
                        steps,
                        t,
                    ),
                    message,
                    step: steps,
                    t,
                });
            }
        };
        t += dt;
        steps += 1;
        if steps >= MARCH_STEP_CAP {
            return Err(Halt {
                message: format!("runaway march: {steps} steps (mis-sized config?)"),
                step: steps,
                t,
                crash_csv: crash_fields_csv(
                    &spec.grid,
                    &spec.fields,
                    &eos,
                    spec.t_solid,
                    "runaway march",
                    steps,
                    t,
                ),
            });
        }
        if steps.is_multiple_of(PROBE_EVERY) {
            let now = snapshot_rho(&spec.grid, &spec.fields);
            resid = max_rel_change(&rho_probe, &now);
            rho_probe = now;
            let exit = exit_plane(&spec.grid);
            let mdot_exit = plane_mdot(&spec.grid, &spec.fields, exit, FaceDir::ZPlus);
            let thrust = plane_thrust(&spec.grid, &spec.fields, &eos, exit, FaceDir::ZPlus)
                .map_err(|e| halt_at(&spec.grid, &spec.fields, &eos, spec.t_solid, e, steps, t))?;
            on_progress(&Progress {
                step: steps,
                t,
                t_final,
                resid,
                mdot_exit,
                thrust_n: thrust,
            });
        }
    }

    // --- SOLV-7 readout ----------------------------------------------------
    let end = |e: String| halt_at(&spec.grid, &spec.fields, &eos, spec.t_solid, e, steps, t);
    let exit = exit_plane(&spec.grid);
    let mdot_exit = plane_mdot(&spec.grid, &spec.fields, exit, FaceDir::ZPlus);
    let thrust = plane_thrust(&spec.grid, &spec.fields, &eos, exit, FaceDir::ZPlus).map_err(end)?;
    let p_c = injector_end_stagnation_p(&spec.grid, &spec.fields, &eos, 0).map_err(end)?;
    let a_t = spec.contour.throat_area();
    let c_star = p_c * a_t / mdot_exit;
    let c_f = thrust / (p_c * a_t);
    let v_e = thrust / mdot_exit;
    let liner_t_max = if spec.wall_law.is_some() {
        max_solid(&spec.grid, spec.t_solid)
    } else {
        f64::NAN
    };
    let fields_csv = fields_csv(&spec.grid, &spec.fields, &eos, spec.t_solid).map_err(end)?;
    Ok(Report {
        thrust_n: thrust,
        isp_s: v_e / G0,
        v_e_m_per_s: v_e,
        c_star_m_per_s: c_star,
        c_f,
        p_c_pa: p_c,
        mdot_exit_kg_per_s: mdot_exit,
        mdot_injected_kg_per_s: spec.injector.mdot_kg_per_s,
        jacket_watts,
        liner_t_max_k: liner_t_max,
        steps,
        t_end: t,
        steady_resid: resid,
        active_gas_cells: n_gas,
        throat_area_m2: a_t,
        fields_csv,
    })
}

/// One coupled step (see module header). Returns the dt taken; updates the
/// running jacket-watts readout with the final-state value.
#[allow(clippy::too_many_arguments)]
fn coupled_step(
    spec: &mut EngineSpec,
    op: &Euler<'_, TableEos<'_>>,
    eos: &TableEos<'_>,
    patches: &[WallPatch],
    t: f64,
    dt_cap: f64,
    jacket_watts: &mut f64,
) -> Result<f64, String> {
    let ids = spec.fields.ids();
    let (dr, dz) = (spec.grid.spec().dr, spec.grid.spec().dz);

    // Wall exchanges from the pre-step state — SOLV-1 §3.5's one law,
    // evaluated once per wall patch: operands from the patch cell's state
    // projected on the interface normal (tangential speed incl. swirl),
    // wall distance the half-cell along the normal (first-order operands,
    // inside the law's declared band).
    let mut q = Vec::with_capacity(patches.len());
    let mut u_series_max = 0.0f64;
    if let (Some(law), Some(liner)) = (&spec.wall_law, &spec.liner) {
        for patch in patches {
            let mut u = [0.0f64; NCOMP];
            for (k, id) in ids.iter().enumerate() {
                u[k] = cell_value(&spec.grid, *id, patch.gas.0, patch.gas.1);
            }
            let w = eos
                .prim_checked(&u)
                .map_err(|e| format!("wall-patch gas state at {:?}: {e}", patch.gas))?;
            let temperature = eos
                .temperature_w(&w)
                .map_err(|e| format!("wall-patch T at {:?}: {e}", patch.gas))?;
            let (n_r, n_z) = patch.n_hat;
            let v_n = w[1] * n_r + w[3] * n_z;
            let u_t = ((w[1] * w[1] + w[3] * w[3] - v_n * v_n).max(0.0) + w[2] * w[2]).sqrt();
            let y = 0.5 * (n_r.abs() * dr + n_z.abs() * dz);
            let gas = NearWallGas {
                rho: w[0],
                u_t,
                temperature,
                y,
            };
            let t_s = cell_value(&spec.grid, spec.t_solid, patch.solid.0, patch.solid.1);
            let ex = law
                .wall_exchange(&gas, t_s, y, liner.kappa_w_per_m_k)
                .map_err(|e| format!("wall exchange at {:?}: {e}", patch.gas))?;
            u_series_max = u_series_max.max(ex.u_series);
            q.push(ex.q);
        }
    }

    // Δt: gas CFL capped by both thermal limits (fail-loud, station-4 rule).
    let dt_gas = op
        .stable_dt(&spec.grid, &spec.fields, spec.cfl)
        .map_err(|e| format!("{e}"))?;
    let dt = dt_gas.min(dt_cap);
    if let Some(liner) = &spec.liner {
        let zero = |_: f64, _: f64, _: f64, _: f64| 0.0;
        let solid_probe = Conduction {
            kappa: liner.kappa_w_per_m_k,
            rho_cp: liner.rho_cp_j_per_m3_k,
            source: &zero,
            domain: Domain::Solid,
            interior: InteriorFaces {
                gas: None,
                exterior: None,
            },
            bcs: Bcs {
                r_inner: FaceBc::HeatFlux(0.0),
                r_outer: FaceBc::HeatFlux(0.0),
                z_lo: FaceBc::HeatFlux(0.0),
                z_hi: FaceBc::HeatFlux(0.0),
            },
        };
        let dt_solid = SOLID_DT_FRAC * solid_probe.stable_dt(&spec.grid, 1.0);
        let dt_exchange = if u_series_max > 0.0 {
            EXCHANGE_DT_FRAC * liner.rho_cp_j_per_m3_k * dr.min(dz) / u_series_max
        } else {
            f64::INFINITY
        };
        if dt_solid.min(dt_exchange) < dt {
            return Err(format!(
                "thermal stability limit under the gas CFL (dt_gas {dt_gas:.3e}, \
                 dt_solid {dt_solid:.3e}, dt_exchange {dt_exchange:.3e}) — raise the liner \
                 rho_cp continuation device or refine; refusing to silently sub-step"
            ));
        }
    }

    // Solid advance (exchanges as gas-face fluxes; coolant Robin outside).
    // A patch's total watts q·|W| are carried to the solid across its
    // spanned grid faces: per-area face flux q·|W|/Σ(face areas), so the
    // solid side integrates exactly the gas side's debit (conservation by
    // construction, smooth-area heat on stair-area carriers).
    if let (Some(liner), Some(jacket)) = (&spec.liner, &spec.jacket) {
        let mut q_map: std::collections::BTreeMap<(usize, usize, u8), f64> =
            std::collections::BTreeMap::new();
        for (patch, &qp) in patches.iter().zip(&q) {
            let stair: f64 = patch.face_areas.iter().sum();
            let q_face = qp * patch.area / stair;
            for face in &patch.faces {
                q_map.insert(
                    (face.solid.0, face.solid.1, dir_code(opposite(face.dir))),
                    q_face,
                );
            }
        }
        let gas_face_q = |i_r: usize, i_z: usize, _j: u32, dir: FaceDir| -> f64 {
            *q_map
                .get(&(i_r, i_z, dir_code(dir)))
                .expect("every gas-facing solid face has an exchange")
        };
        let zero = |_: f64, _: f64, _: f64, _: f64| 0.0;
        let solid_op = Conduction {
            kappa: liner.kappa_w_per_m_k,
            rho_cp: liner.rho_cp_j_per_m3_k,
            source: &zero,
            domain: Domain::Solid,
            interior: InteriorFaces {
                gas: Some(&gas_face_q),
                exterior: Some(FaceBc::Robin {
                    h: jacket.h_w_per_m2_k,
                    t_inf: jacket.t_coolant_k,
                }),
            },
            bcs: Bcs {
                r_inner: FaceBc::HeatFlux(0.0),
                r_outer: FaceBc::Robin {
                    h: jacket.h_w_per_m2_k,
                    t_inf: jacket.t_coolant_k,
                },
                z_lo: FaceBc::HeatFlux(0.0),
                z_hi: FaceBc::HeatFlux(0.0),
            },
        };
        solid_op
            .step(&mut spec.grid, spec.t_solid, spec.rate_solid, t, dt)
            .map_err(|e| format!("liner conduction: {e}"))?;

        // Gas energy debit: the SAME q·|W|, opposite sign, distributed as a
        // uniform specific debit over the patch cell's SRD control volume
        // (Σ ΔE·κV = watts·dt exactly; a sliver cell alone cannot absorb
        // its wall debit any more than its flux update).
        let mut watts = 0.0f64;
        for (patch, &qp) in patches.iter().zip(&q) {
            let w_patch = qp * patch.area;
            let de = -w_patch * dt / patch.debit_kv_sum;
            for &(i_r, i_z) in &patch.debit_cells {
                cell_add(&mut spec.grid, ids[I_EN], i_r, i_z, de);
            }
            watts += w_patch;
        }
        *jacket_watts = watts;
    }

    // Gas advance.
    op.step(&mut spec.grid, &spec.fields, t, dt)
        .map_err(|e| format!("{e}"))?;
    Ok(dt)
}

// --- Plane diagnostics (mask-aware, EOS-threaded — SOLV-7 §3.1/§3.2) -------
// Areas are aperture-weighted on cut worlds (FND-3): the OPEN area of the
// named z-face side is what flow crosses. `side` picks which z-face of the
// cells at `i_z` the plane means: z− for upstream planes (injector end),
// z+ for the exit. Full-box worlds: aperture = 1 exactly.

fn open_area_z(g: &Grid, i_r: usize, i_z: usize, side: FaceDir) -> f64 {
    g.face_area_z(i_r, 1) * g.aperture(i_r, i_z, side)
}

/// Open flow area of z-plane `i_z` on `side`.
pub fn plane_area(g: &Grid, i_z: usize, side: FaceDir) -> f64 {
    (0..g.spec().n_r)
        .filter(|&i_r| g.is_active(i_r, i_z))
        .map(|i_r| open_area_z(g, i_r, i_z, side))
        .sum()
}

/// Mass flow through z-plane `i_z`: Σ ρu_z·a·A_z over active cells.
pub fn plane_mdot(g: &Grid, f: &EulerFields, i_z: usize, side: FaceDir) -> f64 {
    let ids = f.ids();
    (0..g.spec().n_r)
        .filter(|&i_r| g.is_active(i_r, i_z))
        .map(|i_r| cell_value(g, ids[3], i_r, i_z) * open_area_z(g, i_r, i_z, side))
        .sum()
}

/// Vacuum thrust integral over z-plane `i_z` (SOLV-7.1): Σ (ρu_z² + p)·a·A_z.
pub fn plane_thrust<E: EosLaw>(
    g: &Grid,
    f: &EulerFields,
    eos: &E,
    i_z: usize,
    side: FaceDir,
) -> Result<f64, String> {
    let ids = f.ids();
    let mut acc = 0.0f64;
    for i_r in 0..g.spec().n_r {
        if !g.is_active(i_r, i_z) {
            continue;
        }
        let mut u = [0.0f64; NCOMP];
        for (k, id) in ids.iter().enumerate() {
            u[k] = cell_value(g, *id, i_r, i_z);
        }
        let w = eos
            .prim_checked(&u)
            .map_err(|e| format!("thrust plane ({i_r},{i_z}): {e}"))?;
        acc += (w[0] * w[3] * w[3] + w[4]) * open_area_z(g, i_r, i_z, side);
    }
    Ok(acc)
}

/// N11: area-averaged stagnation pressure at the injector-end plane, per
/// cell from the local static state + Mach via the field's own EOS.
/// Weighted by the plane's open (z−) areas.
pub fn injector_end_stagnation_p<E: EosLaw>(
    g: &Grid,
    f: &EulerFields,
    eos: &E,
    i_z: usize,
) -> Result<f64, String> {
    let ids = f.ids();
    let (mut acc, mut area) = (0.0f64, 0.0f64);
    for i_r in 0..g.spec().n_r {
        if !g.is_active(i_r, i_z) {
            continue;
        }
        let mut u = [0.0f64; NCOMP];
        for (k, id) in ids.iter().enumerate() {
            u[k] = cell_value(g, *id, i_r, i_z);
        }
        let w = eos
            .prim_checked(&u)
            .map_err(|e| format!("p_c plane ({i_r},{i_z}): {e}"))?;
        let a = eos.sound_speed_w(&w);
        let m2 = (w[1] * w[1] + w[2] * w[2] + w[3] * w[3]) / (a * a);
        let g1 = w[0] * a * a / w[4]; // Γ₁ from the state itself
        let p0 = w[4] * (1.0 + 0.5 * (g1 - 1.0) * m2).powf(g1 / (g1 - 1.0));
        let a_z = open_area_z(g, i_r, i_z, FaceDir::ZMinus);
        acc += p0 * a_z;
        area += a_z;
    }
    if area <= 0.0 {
        return Err("empty injector-end plane".to_string());
    }
    Ok(acc / area)
}

/// Last z-plane with any active gas (the exit plane of the contour).
pub fn exit_plane(g: &Grid) -> usize {
    (0..g.spec().n_z)
        .rev()
        .find(|&i_z| (0..g.spec().n_r).any(|i_r| g.is_active(i_r, i_z)))
        .expect("an engine has gas somewhere")
}

// --- Small field utilities (random access at coupler rate) -----------------

fn cell_value(g: &Grid, f: crucible_grid::FieldId, i_r: usize, i_z: usize) -> f64 {
    let bi = g
        .brick_index_by_coords((i_r / BRICK) as u32, (i_z / BRICK) as u32)
        .expect("cell in an allocated brick");
    let b = g.brick(bi);
    b.field(f)[b.cell_index(0, (i_r % BRICK) * BRICK + i_z % BRICK)]
}

fn cell_add(g: &mut Grid, f: crucible_grid::FieldId, i_r: usize, i_z: usize, dv: f64) {
    let bi = g
        .brick_index_by_coords((i_r / BRICK) as u32, (i_z / BRICK) as u32)
        .expect("cell in an allocated brick");
    let idx = {
        let b = g.brick(bi);
        b.cell_index(0, (i_r % BRICK) * BRICK + i_z % BRICK)
    };
    g.brick_field_mut(bi, f)[idx] += dv;
}

fn fill_solid(g: &mut Grid, f: crucible_grid::FieldId, value: f64) {
    for bi in 0..g.n_bricks() {
        let (solid, nt) = {
            let b = g.brick(bi);
            (b.solid_mask(), b.n_theta())
        };
        if solid == 0 {
            continue;
        }
        let data = g.brick_field_mut(bi, f);
        for j in 0..nt {
            for local in 0..crucible_grid::BRICK_CELLS {
                if solid & (1u64 << local) != 0 {
                    data[j as usize * crucible_grid::BRICK_CELLS + local] = value;
                }
            }
        }
    }
}

fn dir_code(d: FaceDir) -> u8 {
    match d {
        FaceDir::RMinus => 0,
        FaceDir::RPlus => 1,
        FaceDir::ZMinus => 2,
        FaceDir::ZPlus => 3,
    }
}

fn opposite(d: FaceDir) -> FaceDir {
    match d {
        FaceDir::RMinus => FaceDir::RPlus,
        FaceDir::RPlus => FaceDir::RMinus,
        FaceDir::ZMinus => FaceDir::ZPlus,
        FaceDir::ZPlus => FaceDir::ZMinus,
    }
}

fn count_active(g: &Grid) -> usize {
    let mut n = 0usize;
    g.for_each_active_cell(|_| n += 1);
    n
}

fn snapshot_rho(g: &Grid, f: &EulerFields) -> Vec<f64> {
    let id = f.ids()[0];
    let mut out = Vec::new();
    g.for_each_active_cell(|c| out.push(g.brick(c.bi).field(id)[c.idx]));
    out
}

fn max_rel_change(before: &[f64], after: &[f64]) -> f64 {
    before
        .iter()
        .zip(after)
        .map(|(&b, &a)| ((a - b) / b.abs().max(1e-300)).abs())
        .fold(0.0, f64::max)
}

fn max_solid(g: &Grid, f: crucible_grid::FieldId) -> f64 {
    let mut m = f64::NAN;
    for bi in 0..g.n_bricks() {
        let b = g.brick(bi);
        let solid = b.solid_mask();
        if solid == 0 {
            continue;
        }
        for (local, &v) in b
            .field(f)
            .iter()
            .enumerate()
            .take(crucible_grid::BRICK_CELLS)
        {
            if solid & (1u64 << local) != 0 {
                m = if m.is_nan() { v } else { m.max(v) };
            }
        }
    }
    m
}

/// One row per active/solid cell: the viz feed (Ben's post-checkpoint
/// dataviz goal rides this surface; keep it boring and complete). T comes
/// from the equilibrium surface at each cell's projected state.
fn fields_csv(
    g: &Grid,
    f: &EulerFields,
    eos: &TableEos<'_>,
    t_solid: crucible_grid::FieldId,
) -> Result<String, String> {
    use std::fmt::Write;
    let ids = f.ids();
    let mut out = String::from("r_m,z_m,region,rho,u_r,u_z,p,T,Z,mach\n");
    for i_z in 0..g.spec().n_z {
        for i_r in 0..g.spec().n_r {
            let z = g.z_center(i_z);
            let r = g.r_center(i_r);
            if g.is_active(i_r, i_z) {
                let mut u = [0.0f64; NCOMP];
                for (k, id) in ids.iter().enumerate() {
                    u[k] = cell_value(g, *id, i_r, i_z);
                }
                let w = eos
                    .prim_checked(&u)
                    .map_err(|e| format!("csv ({i_r},{i_z}): {e}"))?;
                let a = eos.sound_speed_w(&w);
                let mach = (w[1] * w[1] + w[2] * w[2] + w[3] * w[3]).sqrt() / a;
                let temp = eos
                    .temperature_w(&w)
                    .map_err(|e| format!("csv T ({i_r},{i_z}): {e}"))?;
                writeln!(
                    out,
                    "{r:.6},{z:.6},gas,{:.6e},{:.4},{:.4},{:.6e},{temp:.2},{:.5},{mach:.4}",
                    w[0], w[1], w[3], w[4], w[5]
                )
                .expect("string write");
            } else if cell_region_is_solid(g, i_r, i_z) {
                let ts = cell_value(g, t_solid, i_r, i_z);
                writeln!(out, "{r:.6},{z:.6},solid,,,,,{ts:.2},,").expect("string write");
            }
        }
    }
    Ok(out)
}

fn halt_at(
    g: &Grid,
    f: &EulerFields,
    eos: &TableEos<'_>,
    t_solid: crucible_grid::FieldId,
    message: String,
    step: usize,
    t: f64,
) -> Halt {
    Halt {
        crash_csv: crash_fields_csv(g, f, eos, t_solid, &message, step, t),
        message,
        step,
        t,
    }
}

/// The crash artifact: one row per active/solid cell, never fails. Raw
/// conserved `U` is always emitted (it is the ground truth at a halt);
/// derived (p, T, Z, Mach) only where the equilibrium projection still
/// succeeds — `ok` = 0 marks the cells whose state left the surface, which
/// is exactly the diagnostic map a starvation/runaway inspection needs.
/// Leading `#` lines carry the halt diagnosis; (i_r, i_z) indices are
/// explicit because halts name cells by index.
pub fn crash_fields_csv(
    g: &Grid,
    f: &EulerFields,
    eos: &TableEos<'_>,
    t_solid: crucible_grid::FieldId,
    message: &str,
    step: usize,
    t: f64,
) -> String {
    use std::fmt::Write;
    let ids = f.ids();
    let mut out = format!(
        "# halt: {}\n# step: {step}  t_s: {t:.9e}\n\
         i_r,i_z,r_m,z_m,region,kappa,rho,mom_r,mom_theta,mom_z,rho_e,rho_c,ok,p,T,Z,mach\n",
        message.replace('\n', " / "),
    );
    for i_z in 0..g.spec().n_z {
        for i_r in 0..g.spec().n_r {
            let z = g.z_center(i_z);
            let r = g.r_center(i_r);
            if g.is_active(i_r, i_z) {
                let mut u = [0.0f64; NCOMP];
                for (k, id) in ids.iter().enumerate() {
                    u[k] = cell_value(g, *id, i_r, i_z);
                }
                write!(
                    out,
                    "{i_r},{i_z},{r:.6},{z:.6},gas,{:.6e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},",
                    g.kappa(i_r, i_z),
                    u[0],
                    u[1],
                    u[2],
                    u[3],
                    u[4],
                    u[5]
                )
                .expect("string write");
                match eos.prim_checked(&u) {
                    Ok(w) => {
                        let a = eos.sound_speed_w(&w);
                        let mach = (w[1] * w[1] + w[2] * w[2] + w[3] * w[3]).sqrt() / a;
                        let temp_s = match eos.temperature_w(&w) {
                            Ok(temp) => format!("{temp:.2}"),
                            Err(_) => String::new(),
                        };
                        writeln!(out, "1,{:.6e},{temp_s},{:.5},{mach:.4}", w[4], w[5])
                            .expect("string write");
                    }
                    Err(_) => writeln!(out, "0,,,,").expect("string write"),
                }
            } else if cell_region_is_solid(g, i_r, i_z) {
                let ts = cell_value(g, t_solid, i_r, i_z);
                writeln!(out, "{i_r},{i_z},{r:.6},{z:.6},solid,,,,,,,,,,{ts:.2},,")
                    .expect("string write");
            }
        }
    }
    out
}

fn cell_region_is_solid(g: &Grid, i_r: usize, i_z: usize) -> bool {
    let Some(bi) = g.brick_index_by_coords((i_r / BRICK) as u32, (i_z / BRICK) as u32) else {
        return false;
    };
    let b = g.brick(bi);
    b.solid_mask() & (1u64 << ((i_r % BRICK) * BRICK + i_z % BRICK)) != 0
}
