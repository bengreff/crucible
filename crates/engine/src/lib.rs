//! The **engine assembly** — the sandbox seam (VISION_SCOPE §4.3): one TOML
//! config = contour geometry + mechanism selections + boundary objects +
//! operating profile → one coupled run → the SOLV-7 performance report. An
//! engine is **pure data** (Rule 13): the RL10 exists in this repository
//! only as `data/anchors/rl10_contour.csv` + `configs/rl10_*.toml`; nothing
//! in this crate names an engine.
//!
//! What runs (S2): the SOLV-1 unified operator in shifting-equilibrium
//! mode (combustion in the EOS) and the one wall-function law at the
//! config-time gas↔solid faces, advanced by **the ONE deterministic
//! SDC-IMEX step** (COUP-3 §3.1, `crucible_solvers::sdc`): liner
//! conduction + the coolant Robin film live inside each sweep's class-`D`
//! implicit solve with the Robin-Robin wall exchange (COUP-2 §3.5), Δt is
//! the gas CFL alone, and the COUP-2 conservation audit is armed every
//! step. The session-11/12 explicit coupled scaffolding is retired.
//! Chamber pressure, thrust, Isp, c\*, C_F are **read out** of the field
//! (SOLV-7; emergent-quantity rule, COUP-7 §3.2) — nothing about the
//! operating point is imposed.
//!
//! Boundary objects (COUP-7 v1 subset, each a registry row whose params
//! carry validity envelopes as refusal ranges): the **prior-tier
//! injector** (premixed inflow at declared ṁ/MR/h_inj + the S18 η_c\*
//! knockdown), the **cooling-jacket coolant side** (declared film +
//! coolant state), and the **closed-mode turbopump** (COUP-3 §3.5 fixed
//! point on the jacket pickup). Deferred, loud: Bartz nozzle-envelope
//! oracle scoring + the COUP-5 band brackets (plan phases); per-quantity
//! three-field registration metadata (citation keys ride these doc
//! comments until COUP-7's registration machinery lands).

pub mod assembly;
pub mod eos_sel;
pub mod geometry;
pub mod run;

use crucible_registry::{
    ChaoticClass, InterfaceVersion, Manifest, ParamSpec, ParamType, ParamValue, PortKind, PortRole,
    PortSpec, Regime, Registry, TableReq,
};

const NON_CHAOTIC: &[(Regime, ChaoticClass)] = &[
    (Regime::Steady, ChaoticClass::NonChaotic),
    (Regime::Transient, ChaoticClass::NonChaotic),
];

/// SOLV-1 §3.4 shifting-equilibrium flow: the same conserved-flux operator
/// as `flow`, EOS occupant = the OFFL-3 (p, h, Z) equilibrium surface (the
/// `chem_equilibrium` logical table pin). No physical params — the physics
/// is in the pinned table; the operator refuses states off its envelope.
pub static FLOW_SHIFTING_MANIFEST: Manifest = Manifest {
    id: "flow_shifting",
    interface_version: InterfaceVersion { major: 1, minor: 0 },
    tables: &[TableReq {
        semantic_name: "chem_equilibrium",
        required: true,
    }],
    couplers: &[],
    ports: &[PortSpec {
        name: "inflow",
        role: PortRole::Require,
        kind: PortKind::Injector,
    }],
    chaotic_class: NON_CHAOTIC,
    params: &[],
};

/// COUP-7 §3.2.1 prior-tier injector (boundary object): premixed inflow at
/// the declared ṁ and MR, inlet enthalpy = the propellant injection
/// enthalpy (the coordinate the chem table was generated around). Citation:
/// META-3 `injector-cstar-eff` (coax-family prior band), `rl10-cycle-data`
/// (calibrated mode only). Envelopes are the refusal ranges below; the MR
/// range is exactly the table's Z envelope (MR 3–8 ⇒ Z ∈ [1/9, 1/4]) so an
/// out-of-family injector cannot silently interrogate the surface.
/// The η_c\* knockdown (source-level, S18) is calibrated in the cycle wave;
/// until then runs are full-equilibrium (η = 1) and say so.
pub static INJECTOR_PRIOR_MANIFEST: Manifest = Manifest {
    id: "injector_prior",
    interface_version: InterfaceVersion { major: 1, minor: 0 },
    tables: &[],
    couplers: &[],
    ports: &[PortSpec {
        name: "combustion_inlet",
        role: PortRole::Provide,
        kind: PortKind::Injector,
    }],
    chaotic_class: NON_CHAOTIC,
    params: &[
        ParamSpec {
            name: "mdot_kg_per_s",
            ty: ParamType::Float,
            // Sanity envelope: micro-thruster to F-1-class flow.
            range: Some((1e-4, 5e3)),
            default: None,
        },
        ParamSpec {
            name: "mixture_ratio",
            ty: ParamType::Float,
            // Family-class sanity range (LOX/LH2 practice). The MECHANICAL
            // envelope gate is the pinned table's own Z envelope at
            // interrogation (narrower per data_version — e.g. the
            // station-5 0.3.x surface covers MR ≈ 4.4–5.5); a config
            // outside it refuses at the table seam, never here
            // (session-12 review: this comment previously claimed to BE
            // the table envelope, which drifts per artifact).
            range: Some((3.0, 8.0)),
            default: None,
        },
        ParamSpec {
            name: "h_inj_j_per_kg",
            ty: ParamType::Float,
            // Cryogenic liquid-injection enthalpy sanity range (CEA
            // formation-referenced). The pinned table's h envelope is the
            // mechanical gate at interrogation (see mixture_ratio note).
            range: Some((-1.3e7, -2.0e5)),
            default: None,
        },
        ParamSpec {
            // The declared η_c\* the knockdown realizes (labeling datum for
            // the certificate: blind = coax-family band value, calibrated =
            // the anchor's fitted number). 1.0 = full equilibrium.
            name: "eta_cstar",
            ty: ParamType::Float,
            range: Some((0.5, 1.0)),
            default: Some(ParamValue::Float(1.0)),
        },
        ParamSpec {
            // S18 source-level realization of eta_cstar: the equilibrium
            // projection interrogates the surface at e + h_offset (a
            // combustion-completeness enthalpy deficit ≤ 0), calibrated so
            // delivered c\* = η_c\*·c\*_ideal at the anchor state (the
            // calibration pair is documented in the certificate). 0 = none.
            name: "h_offset_j_per_kg",
            ty: ParamType::Float,
            range: Some((-2.0e6, 0.0)),
            default: Some(ParamValue::Float(0.0)),
        },
    ],
};

/// COUP-7 §3.2/§3.3 turbopump boundary object, closed (expander) mode: the
/// `drive_power` Require port is bound to the solver's computed jacket
/// enthalpy rise; the pump/turbine machinery is the abstracted map INSIDE
/// the object (no new coupler type — Ben's Fork-2 ruling); the delivered ṁ
/// is solved once per step by COUP-3 §3.5's fixed point (Aitken, fixed
/// sweeps) and feeds the next step's injector inflow. **Emergent-quantity
/// rule:** this object states ṁ (from power balance), never p_c. Citations:
/// `rl10-cycle-data` (TM-107318 Tables 2.2.1/2.3.1/2.4.1/2.6.1 — pump
/// heads/η, turbine PR/flow, feed ΔPs; **calibrated-mode data**),
/// `huzel-huang` (pump-class envelopes for blind closed runs). The map
/// validity envelope is near-design (the declared ṁ window below); an
/// iterate leaving it is the WON'T-BOOTSTRAP physical diagnosis, never
/// clamped (COUP-3 §3.5).
pub static TURBOPUMP_EXPANDER_MANIFEST: Manifest = Manifest {
    id: "turbopump_expander",
    interface_version: InterfaceVersion { major: 1, minor: 0 },
    tables: &[],
    couplers: &[],
    ports: &[
        PortSpec {
            name: "delivered_flow",
            role: PortRole::Provide,
            kind: PortKind::Pump,
        },
        PortSpec {
            name: "drive_power",
            role: PortRole::Require,
            kind: PortKind::PowerSupply,
        },
    ],
    chaotic_class: NON_CHAOTIC,
    params: &[
        ParamSpec {
            name: "mdot_design_kg_per_s",
            ty: ParamType::Float,
            range: Some((1e-4, 5e3)),
            default: None,
        },
        ParamSpec {
            // Total shaft power all pumps demand at the design ṁ
            // (Σ ṁ_i·g·H_i/η_i over the map's stages — pre-config
            // arithmetic documented in the config).
            name: "pump_power_design_w",
            ty: ParamType::Float,
            range: Some((1.0, 1e9)),
            default: None,
        },
        ParamSpec {
            // P_pump ∝ ṁ^n along the fixed feed-system impedance line
            // through the design point (head ∝ ṁ² ⇒ n = 3).
            name: "impedance_exponent",
            ty: ParamType::Float,
            range: Some((2.0, 4.0)),
            default: None,
        },
        ParamSpec {
            name: "turbine_mdot_frac",
            ty: ParamType::Float,
            range: Some((1e-3, 1.0)),
            default: None,
        },
        ParamSpec {
            name: "turbine_eta",
            ty: ParamType::Float,
            range: Some((0.05, 0.95)),
            default: None,
        },
        ParamSpec {
            name: "turbine_pressure_ratio",
            ty: ParamType::Float,
            range: Some((1.001, 100.0)),
            default: None,
        },
        ParamSpec {
            // Effective γ of the (cold-gas) turbine working fluid.
            name: "turbine_gamma",
            ty: ParamType::Float,
            range: Some((1.05, 1.8)),
            default: None,
        },
        ParamSpec {
            name: "coolant_cp_j_per_kg_k",
            ty: ParamType::Float,
            range: Some((1e2, 2e5)),
            default: None,
        },
        ParamSpec {
            name: "coolant_t_in_k",
            ty: ParamType::Float,
            range: Some((10.0, 700.0)),
            default: None,
        },
        ParamSpec {
            // The pump map's declared validity window on delivered ṁ,
            // as fractions of design: an iterate outside is the
            // won't-bootstrap refusal (physical, COUP-3 §3.5).
            name: "mdot_envelope_lo_frac",
            ty: ParamType::Float,
            range: Some((0.05, 1.0)),
            default: None,
        },
        ParamSpec {
            name: "mdot_envelope_hi_frac",
            ty: ParamType::Float,
            range: Some((1.0, 10.0)),
            default: None,
        },
    ],
};

/// COUP-7 cooling-jacket **coolant side only** (D-C): a declared film
/// coefficient + coolant temperature (Robin) standing in for the channel
/// correlation until the jacket object's correlation set lands with the
/// cycle wave. Citation: META-3 `huzel-huang` (channel class),
/// `rl10-tm107318` Table 2.4.1/App. D (jacket data of record). The gas-side
/// h is NEVER stated here — it is SOLV-1 §3.5's one wall law. Provides the
/// `heat_pickup` power port the closed-mode turbopump's `drive_power`
/// binds to (the wall-exchange coupler's Provide, COUP-7 §4).
pub static JACKET_COOLANT_MANIFEST: Manifest = Manifest {
    id: "jacket_coolant",
    interface_version: InterfaceVersion { major: 1, minor: 0 },
    tables: &[],
    couplers: &[],
    ports: &[PortSpec {
        name: "heat_pickup",
        role: PortRole::Provide,
        kind: PortKind::PowerSupply,
    }],
    chaotic_class: NON_CHAOTIC,
    params: &[
        ParamSpec {
            name: "h_w_per_m2_k",
            ty: ParamType::Float,
            // Forced-convection liquid/supercritical films (Colburn class).
            range: Some((1e2, 1e6)),
            default: None,
        },
        ParamSpec {
            name: "t_coolant_k",
            ty: ParamType::Float,
            // Cryogenic inlet to hot-side outlet, any storable coolant.
            range: Some((15.0, 700.0)),
            default: None,
        },
    ],
};

/// SOLV-4 §3.6 burn-progress combustion (S7 config face of the S6 build):
/// the blended thermochemistry (unburnt `chem_unburnt` ⊕ the shifting
/// `chem_equilibrium` surface the flow row already pins) + the SOLV-4.4
/// bistable-Nagumo rate law with the `S_L`/`τ_ign` closures (`chem_ignition`
/// pin). Selecting this row turns the inert `ρc` slot live: the march can
/// ignite, fail to ignite, and flame out (COUP-4 halts). `wrinkling` is the
/// declared turbulent flame-speed factor `S_T/S_L` — citation META-3
/// `turbulent-flame-speed` (Zimont 2000 / Peters 2000 correlation class,
/// banded); 1.0 = laminar (the mini-sim tier). The dynamic SGS-consistent
/// closure rides the resolved-injector wave (plan S16) — until then the
/// factor is declared config data with its citation, exactly like the wall
/// law's band.
pub static COMBUSTION_BLEND_MANIFEST: Manifest = Manifest {
    id: "combustion_blend",
    interface_version: InterfaceVersion { major: 1, minor: 0 },
    tables: &[
        TableReq {
            semantic_name: "chem_unburnt",
            required: true,
        },
        TableReq {
            semantic_name: "chem_ignition",
            required: true,
        },
    ],
    couplers: &[],
    ports: &[],
    chaotic_class: NON_CHAOTIC,
    params: &[ParamSpec {
        name: "wrinkling",
        ty: ParamType::Float,
        // Laminar (1.0) to the strong-wrinkling end of the Zimont/Peters
        // correlation class at engine turbulence intensities.
        range: Some((1.0, 30.0)),
        default: None,
    }],
};

/// COUP-7 §3.3 spark igniter (S7 config face; `spark-igniter-class` pinned
/// META-3 0.8.3): a **literal electrical spark** — a scheduled, localized
/// energy deposit whose ONE cited datum is the deposited energy (H₂ MIE
/// 0.017 mJ floor → aerospace exciter class ~0.1–20 J/discharge; TM-107318:
/// "the ignition source is an electric spark"). Position, extent, and firing
/// window are config *placement* (§3.2.2 schedule class), not sourced
/// claims. Not special-cased physics: the deposit raises the kernel's `h`,
/// `T_u` sees it, and the SOLV-4 §3.6 induction term fires — ignition or
/// no-light is the FIELD's outcome, never this object's claim. The pulse is
/// **bounded** (a spark ends; sustained deposits into a burnt kernel
/// superheat mid-transition cells past the metastable-reactant validity
/// edge — the S6-close discipline).
pub static SPARK_IGNITER_MANIFEST: Manifest = Manifest {
    id: "spark_igniter",
    interface_version: InterfaceVersion { major: 1, minor: 0 },
    tables: &[],
    couplers: &[],
    ports: &[],
    chaotic_class: NON_CHAOTIC,
    params: &[
        ParamSpec {
            // The one cited datum: total deposited energy per firing
            // window. The object spans TWO cited tiers of the same class
            // (S7): the exciter SPARK (H₂ MIE ~1.7e-5 J floor → ~20 J
            // exciter top, `spark-igniter-class`) and the ASI TORCH — the
            // real RL10 augmented spark igniter is a propellant-fed torch
            // that burns CONTINUOUSLY through start into mainstage
            // (TM-107318), delivering kW-class power over the start window
            // (kJ-class energy, cited from the ASI propellant flow ×
            // heating value). The S7 ◆C2 finding that mandates the torch
            // tier: a laminar flame cannot anchor in the ~170 m/s premixed
            // fill stream (blowoff — the prior tier has no resolved
            // recirculation to hold a front), so flame-holding IS the
            // torch's job, exactly as on the real engine.
            name: "energy_j",
            ty: ParamType::Float,
            range: Some((1.0e-5, 5.0e4)),
            default: None,
        },
        ParamSpec {
            name: "r_m",
            ty: ParamType::Float,
            range: Some((0.0, 100.0)),
            default: None,
        },
        ParamSpec {
            name: "z_m",
            ty: ParamType::Float,
            range: Some((-100.0, 100.0)),
            default: None,
        },
        ParamSpec {
            // Half-width of the square deposit kernel around (r_m, z_m); at
            // N_θ = 1 this is a ring (the declared 2-D limitation — a true
            // point spark is the S8/S11 3-D capability).
            name: "half_width_m",
            ty: ParamType::Float,
            range: Some((1.0e-5, 1.0)),
            default: None,
        },
        ParamSpec {
            // Firing-window schedule [s] (COUP-7 §3.2.2 class): start time
            // in march time, window length. The deposit ramps over the
            // first `IGNITER_RAMP_FRAC` of the window (the S3
            // impulsive-drive discipline) and ends at window close.
            name: "window_start_s",
            ty: ParamType::Float,
            range: Some((0.0, 100.0)),
            default: None,
        },
        ParamSpec {
            name: "window_s",
            ty: ParamType::Float,
            range: Some((1.0e-9, 10.0)),
            default: None,
        },
    ],
};

/// SOLV-6 v1 structural margins (S7): the analytic Roark thin-shell subset —
/// hoop + longitudinal + through-wall-gradient thermal stress, von Mises,
/// margins vs T-interpolated A/B-basis allowables with the NASA-STD-5012
/// factors — evaluated on the liner as the one annotated shell component
/// (R(z) from the contour of record, REAL wall thickness `t_real_m` — the
/// grid-thickened ring is a thermal homogenization, never a stress operand).
/// A constraint check, not a structural simulation (Failure-Mode Razor):
/// margin < 1 or liner T ≥ solidus is a COUP-4 halt input (melt/burst).
/// Citations: META-3 `roark` (formulas), `nasa-std-5012` (FS), `mmpds`
/// (allowables class); the concrete allowable values are cited per config.
pub static STRUCTURAL_MARGINS_MANIFEST: Manifest = Manifest {
    id: "structural_margins",
    interface_version: InterfaceVersion { major: 1, minor: 0 },
    tables: &[],
    couplers: &[],
    ports: &[],
    chaotic_class: NON_CHAOTIC,
    params: &[
        ParamSpec {
            // REAL liner wall thickness [m] (the stress operand; the
            // grid ring is thermal-only).
            name: "t_real_m",
            ty: ParamType::Float,
            range: Some((1.0e-5, 0.1)),
            default: None,
        },
        ParamSpec {
            // The declared pressure-shell radius [m] of the annotated
            // component — for a tube-bundle liner (RL10 class) the CITED
            // tube radius, for a monocoque liner the chamber radius. The
            // v1 thin-shell subset requires 2R/t > 20 (SOLV-6 §5); an
            // under-idealized component refuses at assembly (reported,
            // never smeared).
            name: "r_shell_m",
            ty: ParamType::Float,
            range: Some((1.0e-4, 10.0)),
            default: None,
        },
        ParamSpec {
            // Declared coolant-side backpressure [Pa]: the pressure load
            // on the liner shell is |p_gas − p_coolant| (a regen liner is
            // loaded by the DIFFERENCE; the expander jacket typically
            // exceeds chamber pressure). Cited per config.
            name: "p_coolant_pa",
            ty: ParamType::Float,
            range: Some((0.0, 1.0e9)),
            default: None,
        },
        ParamSpec {
            name: "e_pa",
            ty: ParamType::Float,
            range: Some((1.0e9, 1.0e12)),
            default: None,
        },
        ParamSpec {
            name: "alpha_per_k",
            ty: ParamType::Float,
            range: Some((1.0e-7, 1.0e-4)),
            default: None,
        },
        ParamSpec {
            name: "nu",
            ty: ParamType::Float,
            range: Some((0.05, 0.49)),
            default: None,
        },
        ParamSpec {
            name: "yield_cold_pa",
            ty: ParamType::Float,
            range: Some((1.0e6, 5.0e9)),
            default: None,
        },
        ParamSpec {
            name: "yield_hot_pa",
            ty: ParamType::Float,
            range: Some((1.0e6, 5.0e9)),
            default: None,
        },
        ParamSpec {
            name: "uts_cold_pa",
            ty: ParamType::Float,
            range: Some((1.0e6, 5.0e9)),
            default: None,
        },
        ParamSpec {
            name: "uts_hot_pa",
            ty: ParamType::Float,
            range: Some((1.0e6, 5.0e9)),
            default: None,
        },
        ParamSpec {
            // The two cited temperatures the allowables are stated at;
            // interrogation outside [t_cold, t_hot] refuses (no
            // extrapolated strength).
            name: "t_cold_k",
            ty: ParamType::Float,
            range: Some((4.0, 2000.0)),
            default: None,
        },
        ParamSpec {
            name: "t_hot_k",
            ty: ParamType::Float,
            range: Some((100.0, 3000.0)),
            default: None,
        },
        ParamSpec {
            name: "t_solidus_k",
            ty: ParamType::Float,
            range: Some((200.0, 4000.0)),
            default: None,
        },
    ],
};

/// The full production registry: the landed solver rows + this crate's
/// boundary-object/assembly rows, sorted by id (COUP-8 §3.2).
static MECHANISMS: [&Manifest; 13] = [
    &COMBUSTION_BLEND_MANIFEST,
    &crucible_solvers::CONDUCTION_MANIFEST,
    &crucible_solvers::FLOW_MANIFEST,
    &FLOW_SHIFTING_MANIFEST,
    &crucible_solvers::GAS_DIFFUSION_MANIFEST,
    &INJECTOR_PRIOR_MANIFEST,
    &JACKET_COOLANT_MANIFEST,
    &SPARK_IGNITER_MANIFEST,
    &STRUCTURAL_MARGINS_MANIFEST,
    &crucible_solvers::TRANSPORT_CONSTANT_MANIFEST,
    &crucible_solvers::TRANSPORT_TABLE_MANIFEST,
    &TURBOPUMP_EXPANDER_MANIFEST,
    &crucible_solvers::WALL_HEAT_MANIFEST,
];

pub fn registry() -> Registry {
    Registry::new(&MECHANISMS, &[])
}
