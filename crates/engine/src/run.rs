//! The ONE coupled stepper + the SOLV-7 readout, parameterized entirely by
//! the assembled [`EngineSpec`] — no engine-specific code (Rule 13).
//!
//! Step structure (S2 — the scaffolding retired): the ONE deterministic
//! SDC-IMEX step (COUP-3 §3.1, `crucible_solvers::sdc`): explicit
//! hyperbolic class A (the gas operator, slip-ghost/cut-cell machinery
//! unchanged), implicit class D (liner conduction, fixed-cycle CG; coolant
//! Robin on its exterior faces), the Robin-Robin wall exchange inside each
//! sweep's class-D solve (COUP-2 §3.5 — the wall-function h as the Robin
//! coefficient, gas debited exactly what the solid received), and the
//! COUP-2 conservation audit armed every step. Δt = the gas CFL alone —
//! the solid/exchange stability limits are gone by construction, which is
//! what will let S4 give the liner its physical ρc_p.
//!
//! Readout (SOLV-7): every reported number is a plane integral of the
//! conserved `U` — thrust from exit momentum + pressure flux (vacuum,
//! SOLV-7.1), `p_c` = area-averaged stagnation pressure at the injector-end
//! plane (N11 convention), c\*/C_F/Isp per SOLV-7.2–7.4. Emergent, never
//! imposed.

use crate::assembly::{EngineSpec, TransportSpec};
use crucible_constants::G0;
use crucible_grid::{BRICK, FaceDir, Grid};
use crucible_solvers::euler::{
    Cons, EosLaw, Euler, EulerFields, FlowBc, FlowBcs, NCOMP, Prim, TableEos,
};
use crucible_solvers::gas_diffusion::{
    FaceGasBc, GasDiffBcs, GasDiffusion, SpeciesBc, ThermalBc, VelocityBc,
};
use crucible_solvers::sdc::{
    AuditSpec, DiffusionClass, ExchangeClass, FlowClass, GasDiffusionClass, Sdc, build_wall_patches,
};
use crucible_solvers::transport::{TabulatedTransport, TransportProps, TransportSpine};
use crucible_solvers::{Bcs, Conduction, Domain, FaceBc, InteriorFaces};
use crucible_tables::{Pin, Table};

// --- COUP-3 §3.5 closed-mode expander constants (named, deterministic) -----

/// Fixed sweep count of the per-step ṁ fixed point. Sized (doc: "to reach
/// EPS_EXPANDER_RESID with margin") for the measured loop gain ~0.3 under
/// continuous Aitken relaxation (superlinear on the near-linear scalar
/// map; the cold engagement start converges in ~5); warm starts (every
/// subsequent step) in 2–3.
pub const N_EXPANDER_SWEEPS: usize = 8;
/// Residual acceptance |ṁ⁽ᵏ⁾−ṁ⁽ᵏ⁻¹⁾|/ṁ⁽ᵏ⁾ — orders below the pump-map
/// band, so the floor never contributes to the physics error (COUP-3 §3.5).
pub const EPS_EXPANDER_RESID: f64 = 1e-8;
/// Colburn-class jacket-pickup scaling with flow inside the frozen-field
/// loop: Q(ṁ) = Q_field·(ṁ/ṁ_field)^this (COUP-3 §3.5 "h ~ ṁ^0.8").
pub const EXPANDER_Q_EXP: f64 = 0.8;
/// Deterministic Aitken relaxation clamp (COUP-3 §3.5): for a loop gain
/// g ∈ (−1, 1) any ω ∈ (0, 2/(1−g)) keeps the relaxed map a contraction;
/// the secant-optimal ω = 1/(1−g) ≈ 1.4 at the measured RL10 gain, so the
/// clamp ceiling must sit above it.
pub const EXPANDER_OMEGA_MIN: f64 = 0.1;
pub const EXPANDER_OMEGA_MAX: f64 = 2.0;

/// Closed-mode per-step readout (the last accepted solve).
#[derive(Debug, Clone)]
pub struct ExpanderReadout {
    /// Accepted delivered ṁ [kg/s] — next step's injector inflow.
    pub mdot_kg_per_s: f64,
    /// Turbine shaft power at the accepted point [W].
    pub turbine_power_w: f64,
    /// Turbine inlet (jacket outlet) temperature at the accepted point [K].
    pub t_turbine_in_k: f64,
    /// Final sweep residual (≤ EPS_EXPANDER_RESID by acceptance).
    pub resid: f64,
}

/// One COUP-3 §3.5 fixed-point solve at frozen field state: given the
/// step's computed jacket pickup `q_field_w` (measured while the field ran
/// at `mdot_field`), find the delivered ṁ where turbine power balances the
/// pump demand along the declared impedance line. Aitken Δ² every third
/// sweep, fixed count, residual acceptance; an iterate leaving the map's
/// declared ṁ envelope is the WON'T-BOOTSTRAP physical diagnosis, distinct
/// from the numerical COUPLING_RESIDUAL failure (never conflated).
fn expander_fixed_point(
    tp: &crate::assembly::TurbopumpSpec,
    fuel_frac: f64,
    q_field_w: f64,
    mdot_field: f64,
) -> Result<ExpanderReadout, String> {
    let x_isen = 1.0
        - tp.turbine_pressure_ratio
            .powf(-(tp.turbine_gamma - 1.0) / tp.turbine_gamma);
    let lo = tp.mdot_envelope_lo_frac * tp.mdot_design_kg_per_s;
    let hi = tp.mdot_envelope_hi_frac * tp.mdot_design_kg_per_s;
    let g = |mdot: f64| -> (f64, f64, f64) {
        let q = q_field_w * (mdot / mdot_field).powf(EXPANDER_Q_EXP);
        let t_in = tp.coolant_t_in_k + q / (fuel_frac * mdot * tp.coolant_cp_j_per_kg_k);
        let p_t =
            tp.turbine_mdot_frac * mdot * tp.turbine_eta * tp.coolant_cp_j_per_kg_k * t_in * x_isen;
        let mdot_next = tp.mdot_design_kg_per_s
            * (p_t / tp.pump_power_design_w).powf(1.0 / tp.impedance_exponent);
        (mdot_next, p_t, t_in)
    };
    let mut x = mdot_field;
    if !(lo..=hi).contains(&x) {
        return Err(format!(
            "DOESN'T WORK (cycle won't bootstrap): the engagement-point ṁ {x:.4} kg/s is \
             already outside the pump map's declared envelope [{lo:.4}, {hi:.4}] — no \
             admissible starting point (physical diagnosis, COUP-3 §3.5)"
        ));
    }
    let (mut p_t, mut t_in) = (0.0f64, 0.0f64);
    let mut resid = f64::INFINITY;
    let mut prev_r: Option<(f64, f64)> = None; // (x, r) of the last sweep
    let mut omega = 1.0f64;
    for _sweep in 0..N_EXPANDER_SWEEPS {
        let (gx, p, t) = g(x);
        p_t = p;
        t_in = t;
        let r = gx - x;
        // Continuous Aitken relaxation (the scalar secant form), clamped
        // deterministically: ω_k = ω_{k−1}·r_{k−1}/(r_{k−1} − r_k).
        if let Some((_, r_prev)) = prev_r {
            let dr = r_prev - r;
            if dr != 0.0 {
                omega = (omega * r_prev / dr).clamp(EXPANDER_OMEGA_MIN, EXPANDER_OMEGA_MAX);
            }
        }
        prev_r = Some((x, r));
        let x_next = x + omega * r;
        resid = (x_next - x).abs() / x_next.abs().max(1e-300);
        x = x_next;
        if !(lo..=hi).contains(&x) {
            return Err(format!(
                "DOESN'T WORK (cycle won't bootstrap): delivered-ṁ iterate {x:.4} kg/s \
                 left the pump map's declared envelope [{lo:.4}, {hi:.4}] — an engine \
                 with no operating point inside the map's validity (physical diagnosis, \
                 COUP-3 §3.5; not a solver defect)"
            ));
        }
    }
    if !resid.is_finite() || resid > EPS_EXPANDER_RESID {
        return Err(format!(
            "COUPLING_RESIDUAL: expander fixed point residual {resid:.3e} after \
             {N_EXPANDER_SWEEPS} sweeps exceeds EPS_EXPANDER_RESID {EPS_EXPANDER_RESID:.1e} \
             (numerical — a solver defect, never an engine verdict; COUP-3 §3.5)"
        ));
    }
    Ok(ExpanderReadout {
        mdot_kg_per_s: x,
        turbine_power_w: p_t,
        t_turbine_in_k: t_in,
        resid,
    })
}

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
    /// Final injector ṁ: the declared value (open mode) or the last
    /// accepted expander solve (closed mode).
    pub mdot_injected_kg_per_s: f64,
    /// MEASURED inflow-plane mass flow at readout — the sonic-cap honesty
    /// signal (session-12 review): if the injector face is still choked,
    /// this falls short of `mdot_injected_kg_per_s` and the certificate
    /// must not claim the declared flow was delivered.
    pub mdot_inflow_plane_kg_per_s: f64,
    /// Closed-mode expander readout (None ⇔ open mode).
    pub expander: Option<ExpanderReadout>,
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

/// Open the pinned FND-7 §3.3 spine transport surface, when the config
/// selected the tabulated occupant. `None` for the declared-constant
/// occupant, which needs no table.
pub fn open_transport_table(spec: &EngineSpec) -> Result<Option<Table>, String> {
    match &spec.transport {
        TransportSpec::Constant(_) => Ok(None),
        TransportSpec::Table { pin, .. } => {
            let (file, group, version, digest) = pin;
            Table::open(
                file,
                group,
                &Pin {
                    data_version: version.clone(),
                    content_digest: Some(digest.clone()),
                },
            )
            .map(Some)
            .map_err(|e| format!("spine transport table {file}:{group}: {e}"))
        }
    }
}

/// March the assembled engine to its settle budget and read out the
/// performance object. `on_progress` fires every [`PROBE_EVERY`] steps.
/// A mid-march failure returns a [`Halt`] carrying the crash artifact.
pub fn run(
    spec: &mut EngineSpec,
    table: &Table,
    transport_table: Option<&Table>,
    on_progress: &mut dyn FnMut(&Progress),
) -> Result<Report, Halt> {
    let pre = |message: String| Halt {
        message,
        step: 0,
        t: 0.0,
        crash_csv: String::new(),
    };
    let eos = {
        let mut eos = TableEos::bind(table).map_err(pre)?;
        // S18 source-level η_c\* knockdown (0.0 = full equilibrium).
        eos.h_offset = spec.injector.h_offset_j_per_kg;
        eos
    };

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
    let mut mdot_current = spec.injector.mdot_kg_per_s;
    // Startup ramp (declared schedule — the valve-sequence class device):
    // the DELIVERED target ṁ(t) = ṁ·min(1, t/t_ramp). 0 ⇒ step start.
    let ramp_frac = |t: f64, t_ramp: f64| -> f64 {
        if t_ramp <= 0.0 {
            1.0
        } else {
            (t / t_ramp).clamp(0.0, 1.0)
        }
    };
    let inflow = FlowBc::MassFlowInflow {
        mdot_per_area: 0.0, // replaced before the first step below
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
    // Altitude-cell seam: ambient falls log-linearly from the fill pressure
    // to the DECLARED cell floor over the pump-down window, so the nozzle
    // establishes quasi-statically (a violent free drain shocks/starves
    // wall corners); zero upstream influence once the exit runs supersonic.
    // The declared floor must sit inside the pinned table's p envelope —
    // an ambient the surface cannot represent is a config error, refused.
    let p_floor = spec.p_amb_floor_pa;
    if p_floor < eos.envelopes()[0].0 {
        return Err(pre(format!(
            "declared ambient floor {p_floor} Pa is below the pinned table's p envelope \
             floor {} Pa — the plume fringe would be driven off-surface; raise \
             operating_profile.p_amb_floor_pa or regenerate a wider table",
            eos.envelopes()[0].0
        )));
    }
    if p_floor > spec.fill_p_pa {
        return Err(pre(format!(
            "declared ambient floor {p_floor} Pa exceeds the fill pressure {} Pa — the \
             pump-down schedule would pump UP; lower the floor or raise the fill",
            spec.fill_p_pa
        )));
    }
    let p_fill = spec.fill_p_pa;
    // Establishment SEQUENCE (session 12, the real altitude-start order):
    // the cell pumps down FIRST with the injector off — the modest fill
    // drains gently toward the declared floor — and the injector ramp
    // begins once the cell is at altitude, so the nozzle flows FULL as
    // chamber pressure rises: the deeply-overexpanded separation/backflow
    // regime (an ε = 61 bell held at PR ~10 mid-ramp churned h past every
    // ceiling — the dial-12 halts) never exists. Ramp-then-pump and
    // concurrent schedules both manufactured off-surface transients at
    // fine dials; pump-then-ramp is how a vacuum engine actually starts.
    let t_ramp_window = spec.injector_ramp_flowthroughs * span / a_ref;
    let t_pump_window = spec.pumpdown_flowthroughs * span / a_ref;
    let t_pump = t_pump_window + t_ramp_window; // establishment complete
    let pump_schedule = move |t: f64| -> f64 {
        if t_pump_window <= 0.0 {
            return p_floor;
        }
        let s = (t / t_pump_window).clamp(0.0, 1.0);
        p_fill * (p_floor / p_fill).powf(s)
    };
    let mut op = Euler {
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
        combustion: None, // S6 combustion is wired per-config below (mini-sim tier)
    };

    let patches = if spec.wall_law.is_some() {
        let patches = build_wall_patches(&spec.grid).map_err(pre)?;
        // Adiabatic-hole check (session-12 review): a cooled build where a
        // wall-bearing gas cell (nonzero closure vector) has no solid
        // partner would silently skip its exchange — the liner ring failed
        // to cover the contour there. Refuse with the first hole named.
        let covered: std::collections::BTreeSet<(usize, usize)> =
            patches.iter().map(|p| p.gas).collect();
        let mut hole: Option<(usize, usize)> = None;
        spec.grid.for_each_active_cell(|c| {
            if hole.is_some() || covered.contains(&(c.i_r, c.i_z)) {
                return;
            }
            let (w_r, w_z) = spec.grid.wall_closure(c.i_r, c.i_z, 1);
            // Domain-edge cells (injector/exit planes) legitimately carry
            // closure from their open boundary faces; only interior wall
            // cells (a covered r+/r− or interior z face) are holes. The
            // discriminator: a solid or exterior face-neighbor exists.
            let has_wall_nbr = [
                (c.i_r.wrapping_sub(1), c.i_z),
                (c.i_r + 1, c.i_z),
                (c.i_r, c.i_z.wrapping_sub(1)),
                (c.i_r, c.i_z + 1),
            ]
            .iter()
            .any(|&(nr, nz)| {
                nr < spec.grid.spec().n_r
                    && nz < spec.grid.spec().n_z
                    && !spec.grid.is_active(nr, nz)
            });
            if (w_r * w_r + w_z * w_z).sqrt() > 0.0 && has_wall_nbr {
                hole = Some((c.i_r, c.i_z));
            }
        });
        if let Some((i_r, i_z)) = hole {
            return Err(pre(format!(
                "cooled build with an ADIABATIC HOLE: wall cell ({i_r}, {i_z}) has a \
                 nonzero interface but no solid partner — the liner ring does not cover \
                 the contour there; raise liner_thickness_m (the declared grid-thickened \
                 ring must span ≥ 1 cell everywhere)"
            )));
        }
        patches
    } else {
        Vec::new()
    };
    let n_gas = count_active(&spec.grid);

    // --- COUP-2 audit reference scales (§3.1.1 floors): the fill state's
    // stored magnitudes — a pure function of {config, table}, deterministic.
    let audit = {
        let ids = spec.fields.ids();
        let mass_ref = spec.grid.reduce_kappa_volume_weighted(ids[0]).abs();
        let energy_gas = spec.grid.reduce_kappa_volume_weighted(ids[4]).abs();
        let energy_solid = spec.liner.as_ref().map_or(0.0, |l| {
            l.rho_cp_j_per_m3_k * spec.grid.reduce_solid_volume_weighted(spec.t_solid).abs()
        });
        AuditSpec {
            k_audit: crucible_solvers::sdc::K_AUDIT,
            ref_scale: [
                mass_ref,
                mass_ref * a_ref,
                energy_gas + energy_solid,
                mass_ref,
            ],
        }
    };
    let mut sdc = Sdc::with_audit(audit);

    // Class D (liner conduction, coolant Robin on exterior faces) — built
    // once; the exchange data arrives per sweep through the SDC step.
    let zero_heat = |_: f64, _: f64, _: f64, _: f64| 0.0;
    let solid_op = spec.liner.as_ref().map(|liner| {
        let jacket = spec.jacket.as_ref().expect("cooled build has a jacket");
        Conduction {
            kappa: liner.kappa_w_per_m_k,
            rho_cp: liner.rho_cp_j_per_m3_k,
            source: &zero_heat,
            domain: Domain::Solid,
            interior: InteriorFaces {
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
        }
    });
    let gas_temperature = |w: &Prim| -> Result<f64, &'static str> {
        eos.temperature_w(w)
            .map_err(|_| "wall-patch temperature off the pinned surface")
    };

    // --- The FND-7 §3.3 spine (S4) ----------------------------------------
    // ONE provider for the wall law and F_visc alike. The tabulated
    // occupant is interrogated at the EOS's own (p, h, Z) coordinate
    // (`interrogation_php`), so the two surfaces can never be read at
    // different states.
    let spine = match &spec.transport {
        TransportSpec::Constant(c) => TransportSpine::Constant(*c),
        TransportSpec::Table { schmidt, .. } => {
            let t = transport_table.ok_or_else(|| {
                pre(
                    "config selected `transport_table` but no spine transport table \
                     was opened — call `open_transport_table` and pass it to `run`"
                        .to_string(),
                )
            })?;
            TransportSpine::Tabulated(
                TabulatedTransport::bind(t, *schmidt)
                    .map_err(|e| pre(format!("spine transport bind: {e}")))?,
            )
        }
    };
    // The two surfaces must refuse on the same states, or a march can walk
    // off one while the other still answers (OFFL-5 §3.1a). Checked here,
    // at assembly time, against the pins actually loaded — COUP-8 §3.3(2)
    // checks each table against its consumers alone and cannot see this.
    if let TransportSpine::Tabulated(t) = &spine {
        let eos_env = eos.envelopes();
        for (i, name) in ["p", "h", "Z"].iter().enumerate() {
            let (a, b) = t.envelopes()[i];
            let (c, d) = eos_env[i];
            if a > c || b < d {
                return Err(pre(format!(
                    "spine transport {name}-envelope [{a:.6e}, {b:.6e}] does not cover the \
                     equilibrium surface's [{c:.6e}, {d:.6e}] — a march could leave the \
                     transport surface while the EOS still answers"
                )));
            }
        }
    }
    let gas_transport = |w: &Prim| -> Result<TransportProps, &'static str> {
        spine
            .eval(eos.interrogation_php(w))
            .map_err(|_| "transport off the pinned spine surface")
    };

    // --- SOLV-1 §3.1 F_visc (S4 wiring) ------------------------------------
    // Domain-edge viscous conditions, declared as data. The chamber ends are
    // OPEN planes: `Continuative` (zero normal gradient, one-sided
    // tangential stress) — a `FreeSlip` end would truncate the real τ_rz and
    // drive edge vortices (the S3 finding). The r_outer edge is the grid
    // box, not the engine wall: the contour's wall faces are gas↔solid and
    // carry no resolved diffusion at all (the wall law owns them), so the
    // box edge only ever touches exterior cells.
    let gas_op = spec.gas_diffusion.then(|| {
        GasDiffusion::new(GasDiffBcs {
            r_inner: FaceGasBc {
                velocity: VelocityBc::FreeSlip, // the r = 0 axis (zero area anyway)
                thermal: ThermalBc::Adiabatic,
                species: SpeciesBc::ZeroFlux,
            },
            r_outer: FaceGasBc {
                velocity: VelocityBc::FreeSlip,
                thermal: ThermalBc::Adiabatic,
                species: SpeciesBc::ZeroFlux,
            },
            z_lo: FaceGasBc {
                velocity: VelocityBc::Continuative,
                thermal: ThermalBc::Adiabatic,
                species: SpeciesBc::ZeroFlux,
            },
            z_hi: FaceGasBc {
                velocity: VelocityBc::Continuative,
                thermal: ThermalBc::Adiabatic,
                species: SpeciesBc::ZeroFlux,
            },
        })
    });

    // --- March -------------------------------------------------------------
    let mut t = 0.0f64;
    let mut steps = 0usize;
    let mut resid = f64::INFINITY;
    let mut rho_probe = snapshot_rho(&spec.grid, &spec.fields);
    let mut jacket_watts = 0.0f64;
    let mut expander_last: Option<ExpanderReadout> = None;
    while t_final - t > 1e-12 * t_final {
        // The step's inflow: current target ṁ (design, or the expander's
        // last accepted solve) under the declared start ramp, which begins
        // when the pump-down completes (the altitude-start order above).
        op.bcs.z_lo = FlowBc::MassFlowInflow {
            mdot_per_area: mdot_current * ramp_frac(t - t_pump_window, t_ramp_window) / a_inlet,
            h_total: spec.injector.h_inj_j_per_kg,
            c_frac: spec.injector.z_frac,
        };
        let dt_cap = t_final - t;
        // The ONE deterministic step (COUP-3 §3.1): Δt from the gas CFL
        // alone; audit armed; a failure of any kind halts with the crash
        // artifact.
        let step_result = (|| {
            let flow = FlowClass {
                op: &op,
                fields: &spec.fields,
            };
            let diffusion = solid_op.as_ref().map(|sop| DiffusionClass {
                op: sop,
                t_field: spec.t_solid,
                scratch_field: spec.rate_solid,
            });
            let exchange = spec.wall_law.as_ref().map(|law| ExchangeClass {
                patches: &patches,
                law,
                temperature: &gas_temperature,
                transport: &gas_transport,
            });
            let gas = gas_op.as_ref().map(|op| GasDiffusionClass {
                op,
                temperature: &gas_temperature,
                transport: &gas_transport,
            });
            let dt = sdc
                .stable_dt(&spec.grid, &flow, spec.cfl)
                .map_err(|e| format!("{e}"))?
                .min(dt_cap);
            let report = sdc
                .step(
                    &mut spec.grid,
                    Some(&flow),
                    diffusion.as_ref(),
                    gas.as_ref(),
                    exchange.as_ref(),
                    t,
                    dt,
                )
                .map_err(|e| format!("{e}"))?;
            Ok::<(f64, _), String>((dt, report))
        })();
        let dt = match step_result {
            Ok((dt, report)) => {
                if let Some(ex) = &report.exchange {
                    jacket_watts = ex.jacket_w;
                }
                dt
            }
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
        // COUP-3 §3.5: the closed-mode expander solve, once per step after
        // the wall exchange, at frozen field state; the accepted ṁ drives
        // the NEXT step's inflow (the declared one-step lag — identically
        // zero at steady state). Engages once the establishment window
        // (fill drain + pump-down) has run on the declared open schedule:
        // the closed mode models the steady operating point, not the
        // physical start sequence (valve schedules are out of scope).
        if let Some(tp) = &spec.turbopump
            && t >= t_pump
        {
            let ex = expander_fixed_point(tp, spec.injector.z_frac, jacket_watts, mdot_current)
                .map_err(|e| halt_at(&spec.grid, &spec.fields, &eos, spec.t_solid, e, steps, t))?;
            mdot_current = ex.mdot_kg_per_s;
            expander_last = Some(ex);
        }
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
        mdot_injected_kg_per_s: mdot_current,
        mdot_inflow_plane_kg_per_s: plane_mdot(&spec.grid, &spec.fields, 0, FaceDir::ZMinus),
        expander: expander_last,
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
