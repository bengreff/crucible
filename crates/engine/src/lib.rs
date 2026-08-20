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

/// The full production registry: the three landed solver rows + this
/// crate's boundary-object/assembly rows, sorted by id (COUP-8 §3.2).
static MECHANISMS: [&Manifest; 10] = [
    &crucible_solvers::CONDUCTION_MANIFEST,
    &crucible_solvers::FLOW_MANIFEST,
    &FLOW_SHIFTING_MANIFEST,
    &crucible_solvers::GAS_DIFFUSION_MANIFEST,
    &INJECTOR_PRIOR_MANIFEST,
    &JACKET_COOLANT_MANIFEST,
    &crucible_solvers::TRANSPORT_CONSTANT_MANIFEST,
    &crucible_solvers::TRANSPORT_TABLE_MANIFEST,
    &TURBOPUMP_EXPANDER_MANIFEST,
    &crucible_solvers::WALL_HEAT_MANIFEST,
];

pub fn registry() -> Registry {
    Registry::new(&MECHANISMS, &[])
}
