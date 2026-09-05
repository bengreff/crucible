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
use crate::eos_sel::ChemEos;
use crucible_constants::G0;
use crucible_grid::{BRICK, FaceDir, Grid};
use crucible_solvers::euler::{
    Combustion, Cons, EPS_IGNITED, EosLaw, Euler, EulerFields, FlowBc, FlowBcs, I_RB, I_RHO,
    IgnitionColumns, NCOMP, Prim, THETA_CELLS, reacting_measure,
};
use crucible_solvers::gas_diffusion::{
    FaceGasBc, GasDiffBcs, GasDiffusion, SpeciesBc, ThermalBc, VelocityBc,
};
use crucible_solvers::sdc::{
    AuditSpec, DiffusionClass, ExchangeClass, FlowClass, GasDiffusionClass, ReactionClass, Sdc,
    build_wall_patches,
};
use crucible_solvers::structural_margins::{Allowables, sigma_hoop, sigma_long, von_mises_plane};
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

/// Steadiness probe cadence (steps) for the residual + progress callback —
/// and (S7) the declared check cadence for the slow-clock halt members
/// (structural margins, the reacting-measure ignition floor, the dwell
/// tracker): those mechanisms evolve on thermal/flame timescales, orders
/// above the acoustic step, while the per-step members (audit, positivity,
/// non-finite) stay per-step inside the SDC step.
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
    /// Emergent injector-end stagnation pressure (N11) at the probe — the
    /// startup-march witness (S7: fill, light-off, choke, settle are all
    /// visible in this trace).
    pub p_c_pa: f64,
    /// Reacting measure R [kg/s] at the probe (NaN ⇔ no combustion).
    pub reacting_r: f64,
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
    /// The COUP-4 §3.3 verdict object (S7): populated for every mid-march
    /// halt — mechanism + diagnosis + location + time, the structured
    /// DOESN'T-WORK outcome. `None` ⇔ the halt fired before the march
    /// existed (assembly/pre-flight refusals are config errors, not engine
    /// verdicts).
    pub verdict: Option<Verdict>,
}

// --- COUP-4 §3.3 verdict object (S7) ----------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VerdictOutcome {
    Works,
    DoesntWork,
}

/// COUP-4 §3.2/§3.3: a physical diagnosis is an engine conclusion (the
/// mechanism is real physics — never ignited, flamed out, burst, melted,
/// failed to reach); a numerical one is a solver defect (COUPLING_RESIDUAL,
/// audit violation, non-finite field) and NEVER an engine verdict.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Diagnosis {
    Physical,
    Numerical,
}

/// The WORKS-criterion constants, recorded in the verdict (COUP-4 O12).
#[derive(Debug, Clone)]
pub struct VerdictCriterion {
    pub eps_works: f64,
    pub t_dwell_s: f64,
    pub t_s1_horizon_s: f64,
    /// The commanded quantities: (name, commanded value).
    pub commanded: Vec<(String, f64)>,
}

/// The COUP-4 §3.3 member verdict: outcome + mechanism + location + time +
/// diagnosis + the criterion it was judged against. Deterministic; a
/// DOESN'T-WORK always carries its mechanism (never a bare failure).
#[derive(Debug, Clone)]
pub struct Verdict {
    pub outcome: VerdictOutcome,
    pub mechanism: String,
    /// (i_r, i_z) cell when the mechanism is localizable.
    pub location: Option<(usize, usize)>,
    pub t_s: f64,
    pub diagnosis: Diagnosis,
    /// `None` ⇔ no commanded profile was declared (the halt is still a
    /// structured outcome; there was just no WORKS target to judge).
    pub criterion: Option<VerdictCriterion>,
}

impl std::fmt::Display for Verdict {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let outcome = match self.outcome {
            VerdictOutcome::Works => "WORKS",
            VerdictOutcome::DoesntWork => "DOESN'T WORK",
        };
        let diag = match self.diagnosis {
            Diagnosis::Physical => "physical",
            Diagnosis::Numerical => "numerical",
        };
        write!(f, "{outcome} ({}; diagnosis {diag}", self.mechanism)?;
        if let Some((i_r, i_z)) = self.location {
            write!(f, "; cell ({i_r}, {i_z})")?;
        }
        write!(f, "; t = {:.6e} s)", self.t_s)
    }
}

/// The RL10-tier reacting-measure floor scaling (SOLV-4 §3.6 / COUP-4 S7):
/// the ignition floor is `max(EPS_IGNITED, this × the delivered injector
/// ṁ)` — a flame consuming less than this fraction of the injected mass
/// rate is not an established burn (the mini-sim absolute floor does not
/// transfer to engine scale). At a settled front the b(1−b)-weighted
/// measure runs ~10⁻¹ of ṁ, three orders above this floor.
pub const EPS_IGNITED_FRAC: f64 = 1.0e-4;

/// NEVER_IGNITED establishment grace [s] (S7, measured): the check runs
/// from `window_end + this`. A just-lit kernel's reacting measure grows on
/// the FRONT-GROWTH timescale (kernel dimension / S_T — ms class at
/// laminar flame speeds), so checking at the first probe after window end
/// conflates "never lit" with "still establishing" (measured: a lit front
/// at 94% of the ṁ-scaled floor 38 µs after window end). A front that
/// lights and dies without ever crossing the floor still verdicts
/// NEVER_IGNITED once the grace expires — no outcome is masked; the grace
/// is ~1% of the Stage-1 horizon.
pub const NEVER_IGNITED_GRACE_S: f64 = 1.0e-3;

/// Igniter pulse ramp fraction (the S3 impulsive-drive discipline): the
/// deposit ramps linearly over this fraction of the firing window, then
/// holds; the delivered energy integrates to exactly `energy_j` over the
/// window (the power density accounts for the ramp).
pub const IGNITER_RAMP_FRAC: f64 = 0.3;

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
    /// Peak reacting measure R over the march [kg/s] (NaN ⇔ no combustion
    /// scheduled) — the COUP-4 ignition-floor witness.
    pub peak_reacting_measure: f64,
    /// The COUP-4 verdict (S7): `Some(WORKS …)` when a commanded profile
    /// was declared and the dwell completed; `None` for plain settle-budget
    /// marches (the stations).
    pub verdict: Option<Verdict>,
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

/// Open the blend's two extra pinned tables (`chem_unburnt` +
/// `chem_ignition`) when `combustion_blend` is selected; `None` otherwise.
pub fn open_blend_tables(spec: &EngineSpec) -> Result<Option<(Table, Table)>, String> {
    let Some(bl) = &spec.blend else {
        return Ok(None);
    };
    let open = |pin: &(String, String, String, String), what: &str| {
        let (file, group, version, digest) = pin;
        Table::open(
            file,
            group,
            &Pin {
                data_version: version.clone(),
                content_digest: Some(digest.clone()),
            },
        )
        .map_err(|e| format!("{what} table {file}:{group}: {e}"))
    };
    Ok(Some((
        open(&bl.unburnt_pin, "unburnt-reactant")?,
        open(&bl.ignition_pin, "ignition-closure")?,
    )))
}

/// March the assembled engine to its settle budget (or, with a commanded
/// profile, to a completed dwell or the Stage-1 horizon) and read out the
/// performance object. `on_progress` fires every [`PROBE_EVERY`] steps.
/// A mid-march failure returns a [`Halt`] carrying the crash artifact and
/// (S7) the COUP-4 verdict object. `blend_tables` = the `(unburnt,
/// ignition)` pair from [`open_blend_tables`] when `combustion_blend` is
/// selected.
// A Halt carries the crash artifact + verdict by design (constructed once,
// at the halt); boxing it would only obscure the one construction site.
#[allow(clippy::result_large_err)]
pub fn run(
    spec: &mut EngineSpec,
    table: &Table,
    transport_table: Option<&Table>,
    blend_tables: Option<&(Table, Table)>,
    on_progress: &mut dyn FnMut(&Progress),
) -> Result<Report, Halt> {
    let pre = |message: String| Halt {
        message,
        step: 0,
        t: 0.0,
        crash_csv: String::new(),
        verdict: None,
    };
    let eos = match (&spec.blend, blend_tables) {
        (None, _) => ChemEos::bind_table(table, spec.injector.h_offset_j_per_kg).map_err(pre)?,
        (Some(_), Some((unburnt, _))) => {
            // S18 knockdown rides the burnt branch only (blend semantics).
            ChemEos::bind_blend(unburnt, table, spec.injector.h_offset_j_per_kg).map_err(pre)?
        }
        (Some(_), None) => {
            return Err(pre(
                "config selected combustion_blend but the blend tables were not opened — \
                 call open_blend_tables and pass the pair to run"
                    .to_string(),
            ));
        }
    };
    // SOLV-4 §3.6 combustion + the class-R reaction (S7): built against the
    // SAME blend instance the flow operator clones, one closure surface.
    let comb = match (&spec.blend, blend_tables, eos.blend()) {
        (Some(bl), Some((_, ignition)), Some(blend_ref)) => Some(Combustion {
            blend: blend_ref,
            ignition: IgnitionColumns::bind(ignition).map_err(pre)?,
            wrinkling: bl.wrinkling,
            theta: THETA_CELLS,
        }),
        _ => None,
    };

    // --- Initial fill: quiescent gas + cold liner. Shifting mode fills on
    // the equilibrium surface (the inert b ≡ 0 slot); blend mode fills PURE
    // UNBURNT — the physical pre-start chamber contents (S5/S7).
    let u_fill: Cons = eos
        .fill_cons(
            spec.fill_p_pa,
            spec.injector.h_inj_j_per_kg,
            spec.injector.z_frac,
        )
        .map_err(pre)?;
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
    // COUP-7 §3.3 spark igniter (S7): a bounded, scheduled, localized energy
    // deposit — the kernel's measured κV volume converts the ONE cited datum
    // (deposited energy) to a power density that integrates to exactly
    // `energy_j` over the ramped window. No igniter ⇒ the source is zero
    // (identical arithmetic to the retired zero_src fn).
    struct IgniterKernel {
        r_lo: f64,
        r_hi: f64,
        z_lo: f64,
        z_hi: f64,
        /// `None` at N_θ = 1 (a ring — no azimuth). `Some((n_θ, j_target))`
        /// at N_θ > 1: the point spark fires in the single sector `j_target`
        /// (S11 / ◆C3), the source of the asymmetric light-off.
        theta_gate: Option<(u32, u32)>,
        t_on: f64,
        t_off: f64,
        t_ramp: f64,
        q0_w_per_m3: f64,
    }
    let igniter_kernel: Option<IgniterKernel> = match &spec.igniter {
        None => None,
        Some(ig) => {
            let (r_lo, r_hi) = (ig.r_m - ig.half_width_m, ig.r_m + ig.half_width_m);
            let (z_lo, z_hi) = (ig.z_m - ig.half_width_m, ig.z_m + ig.half_width_m);
            // Uniform, geometry-floor-pinned N_θ on the revolved contour.
            let nt_world = spec.grid.spec().n_theta_max;
            let theta_gate = if nt_world > 1 {
                let dth = std::f64::consts::TAU / f64::from(nt_world);
                let jt = ((ig.theta_rad / dth).floor() as u32) % nt_world;
                Some((nt_world, jt))
            } else {
                None
            };
            // The deposit's gas volume — summed over EXACTLY the cells the
            // source fires into (the box, and the target θ-sector at N_θ > 1),
            // each weighted by its own sector volume `κ·cell_volume(i_r, n_θ)`.
            // At N_θ = 1 this is the full-ring volume, bitwise the S7 form.
            let mut vol = 0.0f64;
            spec.grid.for_each_active_cell(|c| {
                let in_box = c.r > r_lo && c.r < r_hi && c.z > z_lo && c.z < z_hi;
                let in_theta = theta_gate.is_none_or(|(_, jt)| c.i_theta == jt);
                if in_box && in_theta {
                    let nt = brick_n_theta(&spec.grid, c.i_r, c.i_z);
                    vol += spec.grid.kappa(c.i_r, c.i_z) * spec.grid.cell_volume(c.i_r, nt);
                }
            });
            if vol <= 0.0 {
                return Err(pre(format!(
                    "spark igniter kernel at (r = {} m, z = {} m, half-width {} m) covers \
                     no gas cell — the deposit would vanish; move it inside the chamber \
                     or widen it past a cell",
                    ig.r_m, ig.z_m, ig.half_width_m
                )));
            }
            let t_ramp = IGNITER_RAMP_FRAC * ig.window_s;
            // Energy delivered = q0·vol·(window − ramp/2) with the linear ramp.
            let effective_s = ig.window_s - 0.5 * t_ramp;
            Some(IgniterKernel {
                r_lo,
                r_hi,
                z_lo,
                z_hi,
                theta_gate,
                t_on: ig.window_start_s,
                t_off: ig.window_start_s + ig.window_s,
                t_ramp,
                q0_w_per_m3: ig.energy_j / (vol * effective_s),
            })
        }
    };
    let source_fn = move |r: f64, th: f64, z: f64, t: f64| -> Cons {
        let mut src = [0.0; NCOMP];
        if let Some(k) = &igniter_kernel
            && t >= k.t_on
            && t < k.t_off
            && r > k.r_lo
            && r < k.r_hi
            && z > k.z_lo
            && z < k.z_hi
        {
            // At N_θ > 1 the deposit is confined to the target sector — the
            // continuous θ maps to its sector index (θ_center(j) = (j+½)·Δθ).
            let in_theta = match k.theta_gate {
                None => true,
                Some((nt, jt)) => {
                    let dth = std::f64::consts::TAU / f64::from(nt);
                    let j = (th / dth).floor().rem_euclid(f64::from(nt)) as u32;
                    j == jt
                }
            };
            if in_theta {
                let ramp = if k.t_ramp > 0.0 {
                    ((t - k.t_on) / k.t_ramp).min(1.0)
                } else {
                    1.0
                };
                src[4] = k.q0_w_per_m3 * ramp;
            }
        }
        src
    };
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
        source: &source_fn,
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
        combustion: comb.as_ref(),
    };
    // Class-R (COUP-3 §3.3, S7): the implicit auto-ignition occupant rides
    // the same Combustion operator the class-A source uses.
    let reaction = comb.as_ref().map(|c| ReactionClass { op: c });

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
        // The blend arm interrogates the spine over the UNION of its two
        // branches' envelopes (a pure-burnt cell roams the whole burnt
        // surface; the cold-side partition admits sub-intersection states),
        // not the projection intersection — review finding.
        let eos_env = eos.interrogation_envelopes();
        for (i, name) in ["p", "h", "Z"].iter().enumerate() {
            let (a, b) = t.envelopes()[i];
            let (c, d) = eos_env[i];
            if a > c || b < d {
                return Err(pre(format!(
                    "spine transport {name}-envelope [{a:.6e}, {b:.6e}] does not cover the \
                     EOS occupant's interrogation range [{c:.6e}, {d:.6e}] — a march could \
                     leave the transport surface while the EOS still answers"
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

    // --- COUP-4 S7: verdict machinery ----------------------------------
    // The WORKS criterion (armed by commanded quantities), the ignition
    // floor (armed by the blend+igniter), and the SOLV-6 margins (armed by
    // structural_margins) — all checked at the declared PROBE_EVERY cadence.
    let criterion = spec.verdict.as_ref().map(|v| {
        let mut commanded = Vec::new();
        if let Some(c) = v.commanded_p_c_pa {
            commanded.push(("p_c_pa".to_string(), c));
        }
        if let Some(c) = v.commanded_thrust_n {
            commanded.push(("thrust_n".to_string(), c));
        }
        VerdictCriterion {
            eps_works: v.eps_works,
            t_dwell_s: v.t_dwell_flowthroughs * span / a_ref,
            t_s1_horizon_s: t_final,
            commanded,
        }
    });
    let mut dwell_start: Option<f64> = None;
    let mut works: Option<Verdict> = None;
    let mut ignited = false;
    let mut peak_r = f64::NAN;
    let ign_window_end = spec
        .igniter
        .as_ref()
        .map(|ig| ig.window_start_s + ig.window_s);
    let margins_allow = match &spec.margins {
        None => None,
        Some(m) => Some((
            m.clone(),
            Allowables::new(
                m.yield_cold_pa,
                m.yield_hot_pa,
                m.uts_cold_pa,
                m.uts_hot_pa,
                m.t_cold_k,
                m.t_hot_k,
                m.t_solidus_k,
            )
            .map_err(|e| pre(format!("structural margins: {e}")))?,
        )),
    };
    // A physical-mechanism halt during a verdict-armed march is a
    // structured DOESN'T-WORK, never a bare failure (COUP-4 §3.3).
    let verdict_halt = |mechanism: String,
                        location: Option<(usize, usize)>,
                        t: f64,
                        diagnosis: Diagnosis,
                        criterion: &Option<VerdictCriterion>| {
        Verdict {
            outcome: VerdictOutcome::DoesntWork,
            mechanism,
            location,
            t_s: t,
            diagnosis,
            criterion: criterion.clone(),
        }
    };

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
                    reaction.as_ref(),
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
                let verdict = Some(verdict_halt(
                    step_halt_mechanism(&message),
                    None,
                    t,
                    classify_diagnosis(&message),
                    &criterion,
                ));
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
                    verdict,
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
                verdict: Some(verdict_halt(
                    "runaway march (step cap)".to_string(),
                    None,
                    t,
                    Diagnosis::Numerical,
                    &criterion,
                )),
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
            // Best-effort probe: a transient state with no stagnation
            // readout yet (mid-fill) reports NaN — a diagnostic must never
            // kill a march the physics didn't kill. The dwell comparison
            // treats NaN as outside-tolerance (the dwell resets).
            let p_c_probe =
                injector_end_stagnation_p(&spec.grid, &spec.fields, &eos, 0).unwrap_or(f64::NAN);

            // COUP-4 ignition floor (S7): the reacting measure vs the
            // ṁ-scaled floor; NEVER_IGNITED once the spark window is
            // exhausted, FLAMEOUT on collapse after establishment.
            let mut r_probe = f64::NAN;
            if let (Some(cb), Some(blend_ref)) = (&comb, eos.blend()) {
                let ids = spec.fields.ids();
                let r_meas = reacting_measure(
                    &spec.grid,
                    &ids,
                    blend_ref,
                    &cb.ignition,
                    cb.wrinkling,
                    cb.theta,
                )
                .map_err(|e| {
                    halt_at(
                        &spec.grid,
                        &spec.fields,
                        &eos,
                        spec.t_solid,
                        format!("reacting measure: {e}"),
                        steps,
                        t,
                    )
                })?;
                r_probe = r_meas;
                peak_r = if peak_r.is_nan() {
                    r_meas
                } else {
                    peak_r.max(r_meas)
                };
                let delivered = mdot_current * ramp_frac(t - t_pump_window, t_ramp_window);
                let floor = EPS_IGNITED.max(EPS_IGNITED_FRAC * delivered);
                if r_meas >= floor {
                    ignited = true;
                } else if ignited {
                    // FLAMEOUT: established burn collapsed before the dwell
                    // completed. Location = the strongest front remnant.
                    let loc = max_flame_cell(&spec.grid, &spec.fields);
                    let message = format!(
                        "FLAMEOUT: the reacting measure collapsed to {r_meas:.3e} kg/s \
                         below its floor {floor:.3e} after establishment (peak \
                         {peak_r:.3e}) — quench/flammability physics (SOLV-4 §3.6)"
                    );
                    let mut h = halt_at(
                        &spec.grid,
                        &spec.fields,
                        &eos,
                        spec.t_solid,
                        message,
                        steps,
                        t,
                    );
                    h.verdict = Some(verdict_halt(
                        "FLAMEOUT".to_string(),
                        loc,
                        t,
                        Diagnosis::Physical,
                        &criterion,
                    ));
                    return Err(h);
                }
                if let Some(w_end) = ign_window_end
                    && !ignited
                    && t > w_end + NEVER_IGNITED_GRACE_S
                {
                    let loc = spec
                        .igniter
                        .as_ref()
                        .map(|ig| nearest_cell(&spec.grid, ig.r_m, ig.z_m));
                    let message = format!(
                        "NEVER_IGNITED: the igniter schedule is exhausted (window ended \
                         at {w_end:.3e} s + the establishment grace) and the reacting \
                         measure never exceeded its \
                         floor (peak {peak_r:.3e} kg/s vs floor {floor:.3e}) — no \
                         self-sustaining burn front exists (SOLV-4 §3.6 / COUP-4 §3.2)"
                    );
                    let mut h = halt_at(
                        &spec.grid,
                        &spec.fields,
                        &eos,
                        spec.t_solid,
                        message,
                        steps,
                        t,
                    );
                    h.verdict = Some(verdict_halt(
                        "NEVER_IGNITED".to_string(),
                        loc.flatten(),
                        t,
                        Diagnosis::Physical,
                        &criterion,
                    ));
                    return Err(h);
                }
            }

            on_progress(&Progress {
                step: steps,
                t,
                t_final,
                resid,
                mdot_exit,
                thrust_n: thrust,
                p_c_pa: p_c_probe,
                reacting_r: r_probe,
            });

            // SOLV-6 margins (S7): melt/burst halt inputs on the liner as
            // the one annotated shell component.
            if let Some((m, allow)) = &margins_allow
                && let Err((message, loc, mech)) =
                    check_margins(&spec.grid, &spec.fields, &eos, spec.t_solid, m, allow)
            {
                let mut h = halt_at(
                    &spec.grid,
                    &spec.fields,
                    &eos,
                    spec.t_solid,
                    message,
                    steps,
                    t,
                );
                h.verdict = Some(verdict_halt(mech, loc, t, Diagnosis::Physical, &criterion));
                return Err(h);
            }

            // COUP-4 dwell tracker (S7): WORKS = every commanded quantity
            // held inside eps_works for a contiguous T_DWELL.
            if let (Some(v), Some(cr)) = (&spec.verdict, &criterion) {
                let p_c_now = p_c_probe;
                let mut all_in = true;
                for (name, cmd) in &cr.commanded {
                    let meas = match name.as_str() {
                        "p_c_pa" => p_c_now,
                        _ => thrust,
                    };
                    let rel = ((meas - cmd) / cmd).abs();
                    // NaN-safe: a NaN readout (no stagnation state yet) is
                    // "not holding", so the dwell resets.
                    if rel.partial_cmp(&v.eps_works) != Some(std::cmp::Ordering::Less)
                        && rel != v.eps_works
                    {
                        all_in = false;
                    }
                }
                if all_in {
                    let since = *dwell_start.get_or_insert(t);
                    if t - since >= cr.t_dwell_s {
                        works = Some(Verdict {
                            outcome: VerdictOutcome::Works,
                            mechanism: format!(
                                "dwell held: every commanded quantity inside \
                                 eps_works = {} for {:.4e} s",
                                v.eps_works, cr.t_dwell_s
                            ),
                            location: None,
                            t_s: t,
                            diagnosis: Diagnosis::Physical,
                            criterion: Some(cr.clone()),
                        });
                        break;
                    }
                } else {
                    dwell_start = None;
                }
            }
        }
    }

    // COUP-4 FAILED_TO_REACH (S7): horizon expired without a completed
    // dwell — a distinct halt naming the offending quantities, never a
    // silent timeout.
    if let (Some(v), Some(cr), None) = (&spec.verdict, &criterion, &works) {
        let end = |e: String| halt_at(&spec.grid, &spec.fields, &eos, spec.t_solid, e, steps, t);
        let exit = exit_plane(&spec.grid);
        let thrust_now =
            plane_thrust(&spec.grid, &spec.fields, &eos, exit, FaceDir::ZPlus).map_err(end)?;
        let p_c_now = injector_end_stagnation_p(&spec.grid, &spec.fields, &eos, 0).map_err(end)?;
        let mut offenders = Vec::new();
        for (name, cmd) in &cr.commanded {
            let meas = match name.as_str() {
                "p_c_pa" => p_c_now,
                _ => thrust_now,
            };
            let rel = (meas - cmd) / cmd;
            if rel.abs() > v.eps_works {
                offenders.push(format!(
                    "{name} {meas:.6e} vs commanded {cmd:.6e} ({rel:+.2}%)",
                    rel = rel * 100.0
                ));
            }
        }
        let detail = if offenders.is_empty() {
            "inside tolerance at the horizon but the dwell never ran its full span".to_string()
        } else {
            offenders.join("; ")
        };
        let message = format!(
            "FAILED_TO_REACH: the Stage-1 horizon ({t_final:.4e} s) expired without a \
             completed dwell — {detail} (COUP-4 §3.2)"
        );
        let mut h = halt_at(
            &spec.grid,
            &spec.fields,
            &eos,
            spec.t_solid,
            message,
            steps,
            t,
        );
        h.verdict = Some(verdict_halt(
            "FAILED_TO_REACH".to_string(),
            None,
            t,
            Diagnosis::Physical,
            &criterion,
        ));
        return Err(h);
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
        peak_reacting_measure: peak_r,
        verdict: works,
    })
}

/// COUP-4 §3.2 diagnosis classification of a step-failure message: the
/// numerical classes (a solver defect, never an engine verdict) are the
/// COUPLING_RESIDUAL residual-acceptance misses, conservation-audit
/// violations, and non-finite fields; everything else a step raises is a
/// physical state the operator refused (off-envelope, positivity,
/// starvation-class) — an engine conclusion.
fn classify_diagnosis(message: &str) -> Diagnosis {
    // Numerical = the solver-defect classes only (COUPLING_RESIDUAL, audit,
    // NaN). A bare "non-finite" substring match would swallow the PHYSICAL
    // "non-positive or non-finite density" starvation refusal (review
    // finding); the stringly-typed subclassification of non-finite-field
    // halts is a v1 limitation retired with FND-6's typed halt enums.
    if message.contains("COUPLING_RESIDUAL")
        || message.contains("AUDIT VIOLATION")
        || message.contains("NaN")
    {
        Diagnosis::Numerical
    } else {
        Diagnosis::Physical
    }
}

/// The mechanism label of a step-failure verdict: the message's first
/// sentence-ish head (the full text stays on `Halt.message`).
fn step_halt_mechanism(message: &str) -> String {
    let head: String = message.chars().take(120).collect();
    if head.len() < message.len() {
        format!("{head}…")
    } else {
        head
    }
}

/// The cell of strongest flame content `b(1−b)` — the FLAMEOUT location
/// witness (the extinction front's remnant). Pure field scan, no closure
/// queries.
fn max_flame_cell(g: &Grid, f: &EulerFields) -> Option<(usize, usize)> {
    let ids = f.ids();
    let mut best: Option<((usize, usize), f64)> = None;
    // Scan every θ-sector (S11): an asymmetric extinction front's remnant may
    // sit off the θ = 0 plane. At N_θ = 1 this is the θ = 0 scan unchanged.
    g.for_each_active_cell(|c| {
        let rho = cell_value_theta(g, ids[I_RHO], c.i_r, c.i_theta, c.i_z);
        let b = (cell_value_theta(g, ids[I_RB], c.i_r, c.i_theta, c.i_z) / rho).clamp(0.0, 1.0);
        let content = b * (1.0 - b);
        if best.is_none_or(|(_, m)| content > m) {
            best = Some(((c.i_r, c.i_z), content));
        }
    });
    best.map(|(loc, _)| loc)
}

/// The active cell nearest a physical (r, z) — the NEVER_IGNITED location
/// witness (the igniter kernel).
fn nearest_cell(g: &Grid, r: f64, z: f64) -> Option<(usize, usize)> {
    let mut best: Option<((usize, usize), f64)> = None;
    g.for_each_active_cell(|c| {
        let dr = g.r_center(c.i_r) - r;
        let dz = g.z_center(c.i_z) - z;
        let d2 = dr * dr + dz * dz;
        if best.is_none_or(|(_, m)| d2 < m) {
            best = Some(((c.i_r, c.i_z), d2));
        }
    });
    best.map(|(loc, _)| loc)
}

/// SOLV-6 v1 margin scan (S7): the liner as the ONE annotated shell
/// component — per z-column, R(z) from the contour of record, the REAL wall
/// thickness from config, ΔT between the column's paired inner/outer liner
/// cells, p from the column's near-wall gas cell. First violation returns
/// `(message, location, mechanism)`; melt (T ≥ solidus) checked first (a
/// melted wall has no margin to report).
/// (message, halt location, mechanism label) of a margin violation.
type MarginViolation = (String, Option<(usize, usize)>, String);

#[allow(clippy::too_many_arguments)]
fn check_margins(
    g: &Grid,
    f: &EulerFields,
    eos: &ChemEos<'_>,
    t_solid: crucible_grid::FieldId,
    m: &crate::assembly::MarginsSpec,
    allow: &Allowables,
) -> Result<(), MarginViolation> {
    let ids = f.ids();
    let n_r = g.spec().n_r;
    for i_z in 0..g.spec().n_z {
        // The column's liner cells: innermost (hot) and outermost (cold).
        let mut inner: Option<(usize, f64)> = None;
        let mut outer: Option<(usize, f64)> = None;
        for i_r in 0..n_r {
            if cell_region_is_solid(g, i_r, i_z) {
                let ts = cell_value(g, t_solid, i_r, i_z);
                if inner.is_none() {
                    inner = Some((i_r, ts));
                }
                outer = Some((i_r, ts));
            }
        }
        let (Some((ir_in, t_in)), Some((_, t_out))) = (inner, outer) else {
            continue; // no liner in this column (open planes)
        };
        if allow.melted(t_in) {
            return Err((
                format!(
                    "MELT: liner inner-surface temperature {t_in:.1} K at z = {:.4} m \
                     (column {i_z}) reached the declared solidus {:.1} K — structural \
                     material phase limit crossed (SOLV-6/COUP-4)",
                    g.z_center(i_z),
                    m.t_solidus_k
                ),
                Some((ir_in, i_z)),
                "MELT".to_string(),
            ));
        }
        // Near-wall gas pressure: the column's outermost active gas cell.
        let mut gas_p: Option<f64> = None;
        for i_r in (0..ir_in.max(1)).rev() {
            if g.is_active(i_r, i_z) {
                let mut u = [0.0f64; NCOMP];
                for (k, id) in ids.iter().enumerate() {
                    u[k] = cell_value(g, *id, i_r, i_z);
                }
                if let Ok(w) = eos.prim_checked(&u) {
                    gas_p = Some(w[4]);
                }
                break;
            }
        }
        let Some(p_gas) = gas_p else { continue };
        // COUP-4's halt input is the BURST margin on the PRIMARY
        // (load-controlled, pressure-difference) stress state — a regen
        // liner's thermal stress is strain-controlled (secondary): it
        // legitimately exceeds elastic yield locally in every real cooled
        // liner (plastic accommodation), so the elastic combined-stress
        // yield margin is a reported diagnostic, never a halt (SOLV-6 0.3,
        // the ASME primary/secondary categorization). The shell load is
        // |p_gas − p_coolant| — the declared coolant backpressure typically
        // EXCEEDS chamber pressure in an expander jacket.
        let dp = (p_gas - m.p_coolant_pa).abs();
        let sh = sigma_hoop(dp, m.r_shell_m, m.t_real_m);
        let sl = sigma_long(dp, m.r_shell_m, m.t_real_m);
        let sigma_primary = von_mises_plane(sh, sl);
        let uts = match allow.uts_at(t_in) {
            Ok(v) => v,
            Err(e) => {
                return Err((
                    format!(
                        "structural allowables interrogated beyond their cited range at \
                         z = {:.4} m (column {i_z}, liner surface {t_in:.1} K): {e}",
                        g.z_center(i_z)
                    ),
                    Some((ir_in, i_z)),
                    "MARGIN_CHECK_REFUSED".to_string(),
                ));
            }
        };
        if sigma_primary > 0.0 {
            let burst_margin = uts / (crucible_solvers::structural_margins::FS_ULT * sigma_primary);
            if burst_margin < 1.0 {
                return Err((
                    format!(
                        "BURST_MARGIN < 1 on the liner at z = {:.4} m (column {i_z}): \
                         primary von Mises {sigma_primary:.4e} Pa vs UTS {uts:.4e} Pa / \
                         FS_ULT (margin {burst_margin:.3}; |Δp| = {dp:.4e} Pa on the \
                         declared R = {} m shell; SOLV-6.5)",
                        g.z_center(i_z),
                        m.r_shell_m
                    ),
                    Some((ir_in, i_z)),
                    "BURST_MARGIN".to_string(),
                ));
            }
        }
        // The through-wall ΔT (t_in − t_out) feeds the secondary-stress
        // diagnostic when the Stage-2/reporting wave lands; the melt and
        // primary-burst checks above are the COUP-4 v1 halt inputs.
        let _ = t_out;
    }
    Ok(())
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

/// Mass flow through z-plane `i_z`: Σ ρu_z·a·A_z over active cells and θ
/// sectors. Each sector carries the per-sector annular area
/// `face_area_z(i_r, n_θ)·aperture`; at N_θ = 1 the single sector's area is
/// the full ring — bit-identical to the pre-S11 form. The z-face `aperture`
/// is the θ-plane-0 view, exact for the revolved (θ-uniform) worlds the
/// engine builds; a genuinely θ-varying wall (CSG/STL, S11) would need the
/// per-sector `aperture_at` here (latent, unreachable until that wave).
pub fn plane_mdot(g: &Grid, f: &EulerFields, i_z: usize, side: FaceDir) -> f64 {
    let ids = f.ids();
    let mut acc = 0.0f64;
    for i_r in 0..g.spec().n_r {
        if !g.is_active(i_r, i_z) {
            continue;
        }
        let nt = brick_n_theta(g, i_r, i_z);
        let a = g.face_area_z(i_r, nt) * g.aperture(i_r, i_z, side);
        for j in 0..nt {
            acc += cell_value_theta(g, ids[3], i_r, j, i_z) * a;
        }
    }
    acc
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
        let nt = brick_n_theta(g, i_r, i_z);
        let a = g.face_area_z(i_r, nt) * g.aperture(i_r, i_z, side);
        for j in 0..nt {
            let mut u = [0.0f64; NCOMP];
            for (k, id) in ids.iter().enumerate() {
                u[k] = cell_value_theta(g, *id, i_r, j, i_z);
            }
            let w = eos
                .prim_checked(&u)
                .map_err(|e| format!("thrust plane ({i_r},{j},{i_z}): {e}"))?;
            acc += (w[0] * w[3] * w[3] + w[4]) * a;
        }
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
        let nt = brick_n_theta(g, i_r, i_z);
        let a_z = g.face_area_z(i_r, nt) * g.aperture(i_r, i_z, FaceDir::ZMinus);
        for j in 0..nt {
            let mut u = [0.0f64; NCOMP];
            for (k, id) in ids.iter().enumerate() {
                u[k] = cell_value_theta(g, *id, i_r, j, i_z);
            }
            let w = eos
                .prim_checked(&u)
                .map_err(|e| format!("p_c plane ({i_r},{j},{i_z}): {e}"))?;
            let a = eos.sound_speed_w(&w);
            let m2 = (w[1] * w[1] + w[2] * w[2] + w[3] * w[3]) / (a * a);
            let g1 = w[0] * a * a / w[4]; // Γ₁ from the state itself
            let p0 = w[4] * (1.0 + 0.5 * (g1 - 1.0) * m2).powf(g1 / (g1 - 1.0));
            acc += p0 * a_z;
            area += a_z;
        }
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
    cell_value_theta(g, f, i_r, 0, i_z)
}

/// Field value in θ-sector `j` of cell (i_r, i_z). At N_θ = 1 (`j = 0`) this
/// is exactly [`cell_value`] (the certified path); the θ-sector plane
/// integrals (S11 / ◆C3) sum this over `0..n_theta`.
fn cell_value_theta(g: &Grid, f: crucible_grid::FieldId, i_r: usize, j: u32, i_z: usize) -> f64 {
    let bi = g
        .brick_index_by_coords((i_r / BRICK) as u32, (i_z / BRICK) as u32)
        .expect("cell in an allocated brick");
    let b = g.brick(bi);
    b.field(f)[b.cell_index(j, (i_r % BRICK) * BRICK + i_z % BRICK)]
}

/// The azimuthal resolution of the brick owning cell (i_r, i_z) — 1 on the
/// certified axisymmetric path, `n_theta_max` on the ◆C3 revolved world
/// (uniform, geometry-floor-pinned).
fn brick_n_theta(g: &Grid, i_r: usize, i_z: usize) -> u32 {
    let bi = g
        .brick_index_by_coords((i_r / BRICK) as u32, (i_z / BRICK) as u32)
        .expect("cell in an allocated brick");
    g.brick(bi).n_theta()
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
#[doc(hidden)]
pub fn fields_csv(
    g: &Grid,
    f: &EulerFields,
    eos: &ChemEos<'_>,
    t_solid: crucible_grid::FieldId,
) -> Result<String, String> {
    use std::fmt::Write;
    let ids = f.ids();
    // A genuinely-3-D world emits one row per (r, θ, z) with a θ column and
    // the swirl velocity (S11 / ◆C3); the certified axisymmetric feed keeps
    // its exact pre-S11 columns (no θ, no u_θ) — the two never mix.
    let world_3d = g.spec().n_theta_max > 1;
    let mut out = String::from(if world_3d {
        "r_m,theta_rad,z_m,region,rho,u_r,u_theta,u_z,p,T,Z,mach\n"
    } else {
        "r_m,z_m,region,rho,u_r,u_z,p,T,Z,mach\n"
    });
    for i_z in 0..g.spec().n_z {
        for i_r in 0..g.spec().n_r {
            let z = g.z_center(i_z);
            let r = g.r_center(i_r);
            if g.is_active(i_r, i_z) {
                let nt = if world_3d {
                    brick_n_theta(g, i_r, i_z)
                } else {
                    1
                };
                for j in 0..nt {
                    let mut u = [0.0f64; NCOMP];
                    for (k, id) in ids.iter().enumerate() {
                        u[k] = cell_value_theta(g, *id, i_r, j, i_z);
                    }
                    let w = eos
                        .prim_checked(&u)
                        .map_err(|e| format!("csv ({i_r},{j},{i_z}): {e}"))?;
                    let a = eos.sound_speed_w(&w);
                    let mach = (w[1] * w[1] + w[2] * w[2] + w[3] * w[3]).sqrt() / a;
                    let temp = eos
                        .temperature_w(&w)
                        .map_err(|e| format!("csv T ({i_r},{j},{i_z}): {e}"))?;
                    if world_3d {
                        let theta = Grid::theta_center(j, nt);
                        writeln!(
                            out,
                            "{r:.6},{theta:.6},{z:.6},gas,{:.6e},{:.4},{:.4},{:.4},{:.6e},{temp:.2},{:.5},{mach:.4}",
                            w[0], w[1], w[2], w[3], w[4], w[5]
                        )
                        .expect("string write");
                    } else {
                        writeln!(
                            out,
                            "{r:.6},{z:.6},gas,{:.6e},{:.4},{:.4},{:.6e},{temp:.2},{:.5},{mach:.4}",
                            w[0], w[1], w[3], w[4], w[5]
                        )
                        .expect("string write");
                    }
                }
            } else if cell_region_is_solid(g, i_r, i_z) {
                let ts = cell_value(g, t_solid, i_r, i_z);
                if world_3d {
                    writeln!(out, "{r:.6},,{z:.6},solid,,,,,,{ts:.2},,").expect("string write");
                } else {
                    writeln!(out, "{r:.6},{z:.6},solid,,,,,{ts:.2},,").expect("string write");
                }
            }
        }
    }
    Ok(out)
}

fn halt_at(
    g: &Grid,
    f: &EulerFields,
    eos: &ChemEos<'_>,
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
        verdict: None,
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
    eos: &ChemEos<'_>,
    t_solid: crucible_grid::FieldId,
    message: &str,
    step: usize,
    t: f64,
) -> String {
    use std::fmt::Write;
    let ids = f.ids();
    // A 3-D halt localizes in θ too (S11): the crash artifact gains an
    // `i_theta`/`theta_rad` pair and one row per sector, so an asymmetric
    // runaway/starvation is visible. N_θ = 1 keeps its exact columns.
    let world_3d = g.spec().n_theta_max > 1;
    let header = if world_3d {
        "i_r,i_z,i_theta,r_m,theta_rad,z_m,region,kappa,rho,mom_r,mom_theta,mom_z,rho_e,rho_c,rho_b,ok,p,T,Z,mach\n"
    } else {
        "i_r,i_z,r_m,z_m,region,kappa,rho,mom_r,mom_theta,mom_z,rho_e,rho_c,rho_b,ok,p,T,Z,mach\n"
    };
    let mut out = format!(
        "# halt: {}\n# step: {step}  t_s: {t:.9e}\n{header}",
        message.replace('\n', " / "),
    );
    for i_z in 0..g.spec().n_z {
        for i_r in 0..g.spec().n_r {
            let z = g.z_center(i_z);
            let r = g.r_center(i_r);
            if g.is_active(i_r, i_z) {
                let nt = if world_3d {
                    brick_n_theta(g, i_r, i_z)
                } else {
                    1
                };
                for j in 0..nt {
                    let mut u = [0.0f64; NCOMP];
                    for (k, id) in ids.iter().enumerate() {
                        u[k] = cell_value_theta(g, *id, i_r, j, i_z);
                    }
                    if world_3d {
                        let theta = Grid::theta_center(j, nt);
                        write!(
                            out,
                            "{i_r},{i_z},{j},{r:.6},{theta:.6},{z:.6},gas,{:.6e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},",
                            g.kappa(i_r, i_z),
                            u[0], u[1], u[2], u[3], u[4], u[5], u[6]
                        )
                        .expect("string write");
                    } else {
                        write!(
                            out,
                            "{i_r},{i_z},{r:.6},{z:.6},gas,{:.6e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},",
                            g.kappa(i_r, i_z),
                            u[0], u[1], u[2], u[3], u[4], u[5], u[6]
                        )
                        .expect("string write");
                    }
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
                }
            } else if cell_region_is_solid(g, i_r, i_z) {
                let ts = cell_value(g, t_solid, i_r, i_z);
                if world_3d {
                    writeln!(out, "{i_r},{i_z},,{r:.6},,{z:.6},solid,,,,,,,,,,,{ts:.2},,")
                        .expect("string write");
                } else {
                    writeln!(out, "{i_r},{i_z},{r:.6},{z:.6},solid,,,,,,,,,,,{ts:.2},,")
                        .expect("string write");
                }
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
