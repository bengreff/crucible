//! SOLV-1 §3.1 / plan S3 — the gas-diffusion (F_visc) mini-sim battery:
//! `visc_channel` (annular Poiseuille marched deep in the stiff regime +
//! Taylor-Couette swirl — both against exact cylindrical solutions),
//! `recovery_couette` (the recovery-temperature analytic — replaces the
//! plan's flat-plate variant with an EXACT compressible solution of the
//! same physics balance: dissipation vs conduction, Pr recovery; a
//! Blasius-layer march is not a laptop mini-sim), `thermal_bl`
//! (conduction-layer erfc growth + species-layer spread at Sc ≠ Pr), the
//! whole-operator MMS with every viscous/conductive/species term active
//! (META-3 `mms`), schedule refusals, and 1-vs-N-thread bit identity
//! through the full step. The COUP-2 audit is armed on every march here —
//! closure is the standing condition, not a separate test.

use crucible_grid::{Grid, GridSpec, Region};
use crucible_solvers::euler::{
    Cons, Euler, EulerFields, FlowBc, FlowBcs, GammaLaw, NCOMP, Prim, fill_from_prim, prim6,
};
use crucible_solvers::euler_mms as emms;
use crucible_solvers::gas_diffusion::{
    FaceGasBc, GasDiffBcs, GasDiffusion, SpeciesBc, ThermalBc, VelocityBc,
};
use crucible_solvers::sdc::{
    DiffusionClass, EPS_GAS_DIFF_RESID, EPS_ROBIN_RESID, ExchangeClass, FlowClass,
    GasDiffusionClass, Sdc, SdcError, build_wall_patches,
};
use crucible_solvers::station1_sod::FIELDS;
use crucible_solvers::station4_cooled_wall::fill_solid;
use crucible_solvers::transport::{ConstantTransport, TransportProps};
use crucible_solvers::wall_heat::WallLaw;
use crucible_solvers::{Bcs, Conduction, Domain, FaceBc, InteriorFaces};
use crucible_units::{dynamic_viscosity_pa_s, specific_heat_capacity_j_per_kg_k};

// --- The one transport set of the channel fixtures (config data) -------------

const GAMMA: f64 = 1.4;
const CP: f64 = 1004.0;
const CV: f64 = CP / GAMMA;
const R_GAS: f64 = CP - CV;
const PR: f64 = 0.72;
const SC: f64 = 0.5;
const T0: f64 = 300.0;
const P0: f64 = 1.0e5;

const ZERO_SRC: fn(f64, f64, f64, f64) -> Cons = |_, _, _, _| [0.0; NCOMP];

fn temperature(w: &Prim) -> Result<f64, &'static str> {
    Ok(w[4] / (w[0] * R_GAS))
}

/// The channel fixtures' FND-7 spine reading: the declared-constant
/// occupant, exactly what a `transport_constant` config block builds. These
/// are constant-coefficient analytic fixtures — a state-varying spine would
/// destroy the exact solutions they are scored against.
fn transport(mu: f64) -> TransportProps {
    ConstantTransport::new(
        specific_heat_capacity_j_per_kg_k(CP),
        dynamic_viscosity_pa_s(mu),
        PR,
        GAMMA,
        SC,
    )
    .expect("transport set")
    .into_props()
}

/// The spine query closure a fixture hands the SDC classes.
macro_rules! spine_query {
    ($props:expr) => {{
        let p = $props;
        move |_: &Prim| -> Result<TransportProps, &'static str> { Ok(p) }
    }};
}

/// CFL-paced march of the coupled flow + gas-diffusion step; returns
/// (steps, final t, smallest dt used).
fn march(
    sdc: &mut Sdc,
    g: &mut Grid,
    flow: &FlowClass<'_, '_, GammaLaw>,
    gas: &GasDiffusionClass<'_>,
    cfl: f64,
    t_final: f64,
) -> (usize, f64, f64) {
    let mut t = 0.0f64;
    let mut steps = 0usize;
    let mut dt_min = f64::INFINITY;
    while t < t_final {
        let dt = sdc.stable_dt(g, flow, cfl).expect("dt").min(t_final - t);
        sdc.step(g, Some(flow), None, Some(gas), None, None, t, dt)
            .expect("step (audit armed)");
        dt_min = dt_min.min(dt);
        t += dt;
        steps += 1;
        assert!(steps < 200_000, "runaway march");
    }
    (steps, t, dt_min)
}

/// Cell-center values of one field along r at fixed (i_z, θ=0).
fn r_profile(g: &Grid, f: &EulerFields, slot: usize, i_z: usize) -> Vec<f64> {
    let ids = f.ids();
    (0..g.spec().n_r)
        .map(|i_r| {
            let bi = g.brick_index(i_r, i_z).expect("active");
            let b = g.brick(bi);
            let idx = b.cell_index(0, Grid::local_rz(i_r, i_z));
            b.field(ids[slot])[idx] / b.field(ids[0])[idx] // momentum → velocity
        })
        .collect()
}

fn t_profile_r(g: &Grid, f: &EulerFields, i_z: usize) -> Vec<f64> {
    let ids = f.ids();
    (0..g.spec().n_r)
        .map(|i_r| {
            let bi = g.brick_index(i_r, i_z).expect("active");
            let b = g.brick(bi);
            let idx = b.cell_index(0, Grid::local_rz(i_r, i_z));
            let u: Cons = std::array::from_fn(|k| b.field(ids[k])[idx]);
            let w = GammaLaw { gamma: GAMMA }.prim_checked(&u).expect("prim");
            w[4] / (w[0] * R_GAS)
        })
        .collect()
}

// =============================================================================
// visc_channel I — annular Poiseuille, marched ~50× beyond the explicit
// viscous bound (the stiff regime class D exists for): body-force-driven
// axial flow between no-slip walls settles onto the exact cylindrical
// profile u(r) = (f/4μ)(R1²−r²) + A·ln(r/R1), A = (f/4μ)(R2²−R1²)/ln(R2/R1).
// =============================================================================

#[test]
fn visc_channel_poiseuille_exact_profile_far_beyond_explicit_bound() {
    const MU: f64 = 3.2;
    const F_Z: f64 = 1.5e8; // declared body force [N/m³]
    let spec = GridSpec {
        r_min: 1.0e-3,
        dr: 1.0e-4,
        n_r: 16,
        z_min: 0.0,
        dz: 2.0e-4,
        n_z: 8,
        n_theta_max: 1,
        axisymmetry_assertion: true,
    };
    let (r1, r2) = (1.0e-3, 2.6e-3);
    let mut g = Grid::build(spec, FIELDS).expect("grid");
    let f = EulerFields::resolve(&g).expect("fields");
    let eos = GammaLaw { gamma: GAMMA };
    let rho0 = P0 / (R_GAS * T0);
    fill_from_prim(&mut g, &f, &eos, |_, _, _| {
        prim6(rho0, 0.0, 0.0, 0.0, P0, 0.5)
    });

    // Declared drive schedule: ramp the body force over ~5 viscous times
    // (an impulsive force at ~50x stiffness rings the truncated sweeps —
    // the walls' declared-schedule discipline applies to drives too).
    let t_ramp = 5.0 * (2.6e-3f64 - 1.0e-3).powi(2) / (MU / (P0 / (R_GAS * T0)));
    let body = move |_: f64, _: f64, _: f64, t: f64| -> Cons {
        [0.0, 0.0, 0.0, F_Z * (t / t_ramp).min(1.0), 0.0, 0.0, 0.0]
    };
    let op = Euler {
        eos,
        source: &body,
        bcs: FlowBcs {
            r_inner: FlowBc::Reflecting,
            r_outer: FlowBc::Reflecting,
            z_lo: FlowBc::Transmissive,
            z_hi: FlowBc::Transmissive,
        },
        wall_normal: None,
        slip_wall_z_faces: true,
        combustion: None,
    };
    let still = |_: f64, _: f64, _: f64, _: f64| (0.0, 0.0, 0.0);
    let wall_t = |_: f64, _: f64, _: f64, _: f64| T0;
    let tr_props = transport(MU);
    let tr_query = spine_query!(tr_props);
    let gas_op = GasDiffusion::new(GasDiffBcs {
        r_inner: FaceGasBc {
            velocity: VelocityBc::NoSlip(&still),
            thermal: ThermalBc::Isothermal(&wall_t),
            species: SpeciesBc::ZeroFlux,
        },
        r_outer: FaceGasBc {
            velocity: VelocityBc::NoSlip(&still),
            thermal: ThermalBc::Isothermal(&wall_t),
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
    });
    let flow = FlowClass {
        op: &op,
        fields: &f,
    };
    let gas = GasDiffusionClass {
        op: &gas_op,
        temperature: &temperature,
        transport: &tr_query,
    };

    let nu = MU / rho0;
    let h = r2 - r1;
    let tau = h * h / nu;
    let mut sdc = Sdc::new();
    let (steps, _, dt_min) = march(&mut sdc, &mut g, &flow, &gas, 0.4, 40.0 * tau);

    // The stiffness claim: every step sat far above the explicit viscous
    // bound (Δt from the gas CFL alone — COUP-3 §3.1's design point).
    let (drr, dzz2) = (g.spec().dr, g.spec().dz);
    let dt_visc_expl = 1.0 / (2.0 * nu * (1.0 / (drr * drr) + 1.0 / (dzz2 * dzz2)));
    let stiffness = dt_min / dt_visc_expl;
    println!("poiseuille: {steps} steps, dt/dt_visc_explicit = {stiffness:.1}");
    assert!(
        stiffness > 30.0,
        "march not in the stiff regime: {stiffness:.1}×"
    );

    // Exact annular Poiseuille profile.
    let a_log = (F_Z / (4.0 * MU)) * (r2 * r2 - r1 * r1) / (r2 / r1).ln();
    let exact = |r: f64| (F_Z / (4.0 * MU)) * (r1 * r1 - r * r) + a_log * (r / r1).ln();
    let u_max = (0..160)
        .map(|i| exact(r1 + (i as f64 + 0.5) * (h / 160.0)))
        .fold(0.0f64, f64::max);
    let uz = r_profile(&g, &f, 3, 4);
    let mut worst = 0.0f64;
    for (i_r, u) in uz.iter().enumerate() {
        let r = g.r_center(i_r);
        worst = worst.max((u - exact(r)).abs() / u_max);
    }
    println!("poiseuille: u_max {u_max:.2}, worst rel err {worst:.3e}");
    assert!(worst < 2.0e-2, "profile off: {worst:.3e}");
    // No secondary flow at steady state.
    let ur = r_profile(&g, &f, 1, 4);
    assert!(ur.iter().all(|u| u.abs() < 1.0e-2 * u_max), "u_r not ~0");
}

// =============================================================================
// visc_channel II — Taylor-Couette swirl: inner cylinder spinning (no-slip
// u_θ = W at r = R1), outer at rest. Steady u_θ(r) = A·r + B/r exactly —
// this pins the τ_rθ operator INCLUDING its geometric −u_θ/r part (a wrong
// discrete form cannot produce both this profile and the rigid-rotation
// invariant of the unit battery).
// =============================================================================

#[test]
fn visc_channel_taylor_couette_swirl_matches_exact_profile() {
    const MU: f64 = 3.2;
    const W: f64 = 20.0; // inner-wall speed [m/s]
    let spec = GridSpec {
        r_min: 1.0e-3,
        dr: 1.0e-4,
        n_r: 16,
        z_min: 0.0,
        dz: 2.0e-4,
        n_z: 8,
        n_theta_max: 1,
        axisymmetry_assertion: true,
    };
    let (r1, r2) = (1.0e-3, 2.6e-3);
    let mut g = Grid::build(spec, FIELDS).expect("grid");
    let f = EulerFields::resolve(&g).expect("fields");
    let eos = GammaLaw { gamma: GAMMA };
    let rho0 = P0 / (R_GAS * T0);
    fill_from_prim(&mut g, &f, &eos, |_, _, _| {
        prim6(rho0, 0.0, 0.0, 0.0, P0, 0.5)
    });

    let op = Euler {
        eos,
        source: &ZERO_SRC,
        bcs: FlowBcs {
            r_inner: FlowBc::Reflecting,
            r_outer: FlowBc::Reflecting,
            z_lo: FlowBc::Transmissive,
            z_hi: FlowBc::Transmissive,
        },
        wall_normal: None,
        slip_wall_z_faces: true,
        combustion: None,
    };
    // Declared wall schedule: ramp the spin over ~5 viscous times — an
    // impulsive wall at ~50x stiffness rings the truncated sweeps past
    // gas positivity (the same declared-schedule discipline as COUP-7).
    let rho0_for_ramp = P0 / (R_GAS * T0);
    let t_ramp = 5.0 * (2.6e-3f64 - 1.0e-3).powi(2) / (MU / rho0_for_ramp);
    let spin = move |_: f64, _: f64, _: f64, t: f64| (0.0, W * (t / t_ramp).min(1.0), 0.0);
    let still = |_: f64, _: f64, _: f64, _: f64| (0.0, 0.0, 0.0);
    let wall_t = |_: f64, _: f64, _: f64, _: f64| T0;
    let tr_props = transport(MU);
    let tr_query = spine_query!(tr_props);
    let gas_op = GasDiffusion::new(GasDiffBcs {
        r_inner: FaceGasBc {
            velocity: VelocityBc::NoSlip(&spin),
            thermal: ThermalBc::Isothermal(&wall_t),
            species: SpeciesBc::ZeroFlux,
        },
        r_outer: FaceGasBc {
            velocity: VelocityBc::NoSlip(&still),
            thermal: ThermalBc::Isothermal(&wall_t),
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
    });
    let flow = FlowClass {
        op: &op,
        fields: &f,
    };
    let gas = GasDiffusionClass {
        op: &gas_op,
        temperature: &temperature,
        transport: &tr_query,
    };

    let nu = MU / rho0;
    let tau = (r2 - r1) * (r2 - r1) / nu;
    let mut sdc = Sdc::new();
    let (steps, _, _) = march(&mut sdc, &mut g, &flow, &gas, 0.4, 40.0 * tau);

    // Exact: u_θ = A·r + B/r with u_θ(R1) = W, u_θ(R2) = 0.
    let a = W * r1 / (r1 * r1 - r2 * r2);
    let b = -a * r2 * r2;
    let exact = |r: f64| a * r + b / r;
    let ut = r_profile(&g, &f, 2, 4);
    let mut worst = 0.0f64;
    for (i_r, u) in ut.iter().enumerate() {
        let r = g.r_center(i_r);
        worst = worst.max((u - exact(r)).abs() / W);
    }
    println!("taylor-couette: {steps} steps, worst rel err {worst:.3e}");
    assert!(worst < 2.0e-2, "swirl profile off: {worst:.3e}");
    let ur = r_profile(&g, &f, 1, 4);
    assert!(
        ur.iter().all(|u| u.abs() < 0.02 * W),
        "secondary flow too big"
    );
}

// =============================================================================
// recovery_couette — compressible Couette with an adiabatic inner wall and
// a moving isothermal outer wall. With constant μ, k the EXACT solution
// (any Mach, exact cylindrical) is u_z = U·ln(r/R1)/ln(R2/R1) and
// T = T_w + (μU²/2k)·(1 − ln²(r/R1)/ln²(R2/R1)); the adiabatic-wall
// recovery is T(R1) − T_w = Pr·U²/(2c_p) exactly — the dissipation ↔
// conduction balance and the moving wall's ledgered work, all in one
// analytic. (The plan's flat-plate recovery mini-sim is replaced by this
// exact fixture — same physics term-balance, no smeared Blasius layer.)
// =============================================================================

#[test]
fn recovery_couette_reproduces_the_exact_recovery_temperature() {
    const MU: f64 = 3.2;
    const U_WALL: f64 = 100.0;
    let spec = GridSpec {
        r_min: 5.0e-3,
        dr: 1.0e-4,
        n_r: 16,
        z_min: 0.0,
        dz: 2.0e-4,
        n_z: 8,
        n_theta_max: 1,
        axisymmetry_assertion: true,
    };
    let (r1, r2) = (5.0e-3, 6.6e-3);
    let mut g = Grid::build(spec, FIELDS).expect("grid");
    let f = EulerFields::resolve(&g).expect("fields");
    let eos = GammaLaw { gamma: GAMMA };
    let rho0 = P0 / (R_GAS * T0);
    fill_from_prim(&mut g, &f, &eos, |_, _, _| {
        prim6(rho0, 0.0, 0.0, 0.0, P0, 0.5)
    });

    let op = Euler {
        eos,
        source: &ZERO_SRC,
        bcs: FlowBcs {
            r_inner: FlowBc::Reflecting,
            r_outer: FlowBc::Reflecting,
            z_lo: FlowBc::Transmissive,
            z_hi: FlowBc::Transmissive,
        },
        wall_normal: None,
        slip_wall_z_faces: true,
        combustion: None,
    };
    let still = |_: f64, _: f64, _: f64, _: f64| (0.0, 0.0, 0.0);
    // Declared wall schedule (see the Taylor-Couette fixture).
    let rho0_for_ramp = P0 / (R_GAS * T0);
    let t_ramp = 5.0 * (6.6e-3f64 - 5.0e-3).powi(2) / (MU / rho0_for_ramp);
    let slide = move |_: f64, _: f64, _: f64, t: f64| (0.0, 0.0, U_WALL * (t / t_ramp).min(1.0));
    let wall_t = |_: f64, _: f64, _: f64, _: f64| T0;
    let tr_props = transport(MU);
    let tr_query = spine_query!(tr_props);
    let gas_op = GasDiffusion::new(GasDiffBcs {
        r_inner: FaceGasBc {
            velocity: VelocityBc::NoSlip(&still),
            thermal: ThermalBc::Adiabatic, // the recovery wall
            species: SpeciesBc::ZeroFlux,
        },
        r_outer: FaceGasBc {
            velocity: VelocityBc::NoSlip(&slide),
            thermal: ThermalBc::Isothermal(&wall_t),
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
    });
    let flow = FlowClass {
        op: &op,
        fields: &f,
    };
    let gas = GasDiffusionClass {
        op: &gas_op,
        temperature: &temperature,
        transport: &tr_query,
    };

    let nu = MU / rho0;
    let tau = (r2 - r1) * (r2 - r1) / nu;
    let mut sdc = Sdc::new();
    let (steps, _, _) = march(&mut sdc, &mut g, &flow, &gas, 0.4, 100.0 * tau);

    let k_gas = MU * CP / PR;
    let lg = |r: f64| (r / r1).ln() / (r2 / r1).ln();
    let u_exact = |r: f64| U_WALL * lg(r);
    let dt_rec = PR * U_WALL * U_WALL / (2.0 * CP); // = μU²/2k
    assert!((dt_rec - MU * U_WALL * U_WALL / (2.0 * k_gas)).abs() < 1e-12);
    let t_exact = |r: f64| T0 + dt_rec * (1.0 - lg(r) * lg(r));

    let uz = r_profile(&g, &f, 3, 4);
    let tt = t_profile_r(&g, &f, 4);
    let mut worst_u = 0.0f64;
    let mut worst_t = 0.0f64;
    for i_r in 0..g.spec().n_r {
        let r = g.r_center(i_r);
        worst_u = worst_u.max((uz[i_r] - u_exact(r)).abs() / U_WALL);
        worst_t = worst_t.max((tt[i_r] - t_exact(r)).abs() / dt_rec);
    }
    let recovery_seen = tt[0] - T0;
    let recovery_at_cell = t_exact(g.r_center(0)) - T0;
    println!(
        "recovery couette: {steps} steps; ΔT_rec = {dt_rec:.3} K, innermost cell \
         {recovery_seen:.3} K (analytic there {recovery_at_cell:.3}); worst u \
         {worst_u:.3e}, worst T {worst_t:.3e} (of ΔT_rec)"
    );
    assert!(worst_u < 1.5e-2, "velocity profile off: {worst_u:.3e}");
    assert!(worst_t < 3.0e-2, "recovery profile off: {worst_t:.3e}");
}

// =============================================================================
// thermal_bl — transient conduction-layer growth into quiescent gas from a
// hot wall (erfc profile at the constant-pressure diffusivity) and species
// layer spread at Sc ≠ Pr in the same march (both coefficients pinned
// independently).
// =============================================================================

#[test]
fn thermal_bl_layer_growth_matches_erfc_and_species_spread() {
    const MU: f64 = 6.5e-5;
    const DELTA_T: f64 = 6.0; // 2% of T0: property variation within tolerance
    let spec = GridSpec {
        r_min: 1.0,
        dr: 2.0e-5,
        n_r: 8,
        z_min: 0.0,
        dz: 1.0e-5,
        n_z: 64,
        n_theta_max: 1,
        axisymmetry_assertion: true,
    };
    let dz = spec.dz;
    let n_z = spec.n_z;
    let mut g = Grid::build(spec, FIELDS).expect("grid");
    let f = EulerFields::resolve(&g).expect("fields");
    let eos = GammaLaw { gamma: GAMMA };
    let rho0 = P0 / (R_GAS * T0);
    let z_step = 32.0 * dz; // the species step sits at a cell face
    fill_from_prim(&mut g, &f, &eos, |_, _, z| {
        let c = if z < z_step { 1.0 } else { 0.0 };
        prim6(rho0, 0.0, 0.0, 0.0, P0, c)
    });

    let op = Euler {
        eos,
        source: &ZERO_SRC,
        bcs: FlowBcs {
            r_inner: FlowBc::Reflecting,
            r_outer: FlowBc::Reflecting,
            z_lo: FlowBc::Reflecting,
            z_hi: FlowBc::Reflecting,
        },
        wall_normal: None,
        slip_wall_z_faces: true,
        combustion: None,
    };
    let still = |_: f64, _: f64, _: f64, _: f64| (0.0, 0.0, 0.0);
    let hot = |_: f64, _: f64, _: f64, _: f64| T0 + DELTA_T;
    let tr_props = transport(MU);
    let tr_query = spine_query!(tr_props);
    let gas_op = GasDiffusion::new(GasDiffBcs {
        r_inner: FaceGasBc::free(),
        r_outer: FaceGasBc::free(),
        z_lo: FaceGasBc {
            velocity: VelocityBc::NoSlip(&still),
            thermal: ThermalBc::Isothermal(&hot),
            species: SpeciesBc::ZeroFlux,
        },
        z_hi: FaceGasBc::free(),
    });
    let flow = FlowClass {
        op: &op,
        fields: &f,
    };
    let gas = GasDiffusionClass {
        op: &gas_op,
        temperature: &temperature,
        transport: &tr_query,
    };

    // March to t*: the thermal layer reaches ~6 cells (√(α·t*) = 6Δz).
    let alpha = MU / (rho0 * PR); // k/(ρc_p): the constant-pressure diffusivity
    let d_c = MU / (SC * rho0); // species diffusivity ρD/ρ
    let target = 6.0 * dz;
    let t_star = target * target / alpha;
    let mut sdc = Sdc::new();
    let (steps, t_end, _) = march(&mut sdc, &mut g, &flow, &gas, 0.4, t_star);
    assert!((t_end - t_star).abs() < 1e-12 * t_star);

    let ids = f.ids();
    let i_r = 4usize;
    let mut worst_t = 0.0f64;
    let mut worst_c = 0.0f64;
    for i_z in 0..n_z {
        let z = g.z_center(i_z);
        let bi = g.brick_index(i_r, i_z).expect("active");
        let b = g.brick(bi);
        let idx = b.cell_index(0, Grid::local_rz(i_r, i_z));
        let u: Cons = std::array::from_fn(|k| b.field(ids[k])[idx]);
        let w = eos.prim_checked(&u).expect("prim");
        let theta_num = (w[4] / (w[0] * R_GAS) - T0) / DELTA_T;
        let theta_exact = erfc(z / (2.0 * (alpha * t_star).sqrt()));
        worst_t = worst_t.max((theta_num - theta_exact).abs());
        let c_exact = 0.5 * erfc((z - z_step) / (2.0 * (d_c * t_star).sqrt()));
        worst_c = worst_c.max((w[5] - c_exact).abs());
    }
    println!(
        "thermal_bl: {steps} steps to t* = {t_star:.3e} s; worst Θ err {worst_t:.3e}, \
         worst C err {worst_c:.3e}"
    );
    assert!(worst_t < 5.0e-2, "thermal layer off: {worst_t:.3e}");
    assert!(worst_c < 5.0e-2, "species layer off: {worst_c:.3e}");
}

/// Abramowitz–Stegun 7.1.26 complementary error function (|ε| < 1.5e-7 —
/// far below the fixture tolerances).
fn erfc(x: f64) -> f64 {
    if x < 0.0 {
        return 2.0 - erfc(-x);
    }
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let poly = t
        * (0.254829592
            + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    poly * (-x * x).exp()
}

// =============================================================================
// MMS — the whole coupled operator (Euler + F_visc + conduction + species
// diffusion, all terms active) recovers formal order against the shared
// manufactured field, with the analytic viscous residual assembled from
// the continuous stress formulas (an independent formulation of the same
// operator the code discretizes). META-3 `mms`; VAL-3 "whole PDE operator".
// =============================================================================

const MMS_MU: f64 = 0.15;
const MMS_CV: f64 = 2.5; // R = (γ−1)c_v = 1 ⇒ T = p/ρ
const MMS_CP: f64 = 1.4 * MMS_CV;
const MMS_PR: f64 = 0.7;
const MMS_SC: f64 = 0.6;

/// −D_visc(W_mms): the exact viscous/conductive/species residual of the
/// manufactured field (axisymmetric, ∂θ = 0 — the eps = 0 study).
#[allow(clippy::similar_names)]
fn mms_visc_residual(r: f64, z: f64, t: f64) -> Cons {
    let mf = emms::manufactured(r, 0.0, z, t, 0.0);
    let k_gas = MMS_MU * MMS_CP / MMS_PR;
    let rho_d = MMS_MU / MMS_SC;
    let rg = MMS_CP - MMS_CV;

    let v = |k: usize| mf.w[k];
    let d_r = |k: usize| mf.dr[k];
    let d_z = |k: usize| mf.dz[k];
    let d_rr = |k: usize| mf.drr[k];
    let d_rz = |k: usize| mf.drz[k];
    let d_zz = |k: usize| mf.dzz[k];

    let (rho, ur, ut, uz, p, c) = (v(0), v(1), v(2), v(3), v(4), v(5));
    let (rho_r, ur_r, ut_r, uz_r, p_r, _c_r) = (d_r(0), d_r(1), d_r(2), d_r(3), d_r(4), d_r(5));
    let (rho_z, ur_z, ut_z, uz_z, p_z, _c_z) = (d_z(0), d_z(1), d_z(2), d_z(3), d_z(4), d_z(5));
    let (rho_rr, ur_rr, ut_rr, uz_rr, p_rr) = (d_rr(0), d_rr(1), d_rr(2), d_rr(3), d_rr(4));
    let (ur_rz, ut_rz, uz_rz) = (d_rz(1), d_rz(2), d_rz(3));
    let (rho_zz, ur_zz, ut_zz, uz_zz, p_zz) = (d_zz(0), d_zz(1), d_zz(2), d_zz(3), d_zz(4));
    let _ = ut_rz;

    let mu = MMS_MU;
    let tau_rr = mu * ((4.0 / 3.0) * ur_r - (2.0 / 3.0) * (ur / r + uz_z));
    let tau_thth = mu * ((4.0 / 3.0) * ur / r - (2.0 / 3.0) * (ur_r + uz_z));
    let tau_zz = mu * ((4.0 / 3.0) * uz_z - (2.0 / 3.0) * (ur_r + ur / r));
    let tau_rz = mu * (ur_z + uz_r);
    let tau_rth = mu * (ut_r - ut / r);
    let tau_thz = mu * ut_z;

    let tau_rr_r = mu * ((4.0 / 3.0) * ur_rr - (2.0 / 3.0) * (ur_r / r - ur / (r * r) + uz_rz));
    let tau_rz_z = mu * (ur_zz + uz_rz);
    let tau_rz_r = mu * (ur_rz + uz_rr);
    let tau_zz_z = mu * ((4.0 / 3.0) * uz_zz - (2.0 / 3.0) * (ur_rz + ur_z / r));
    let tau_rth_r = mu * (ut_rr - ut_r / r + ut / (r * r));
    let tau_thz_z = mu * ut_zz;

    let d_mr = tau_rr / r + tau_rr_r + tau_rz_z - tau_thth / r;
    let d_mt = 2.0 * tau_rth / r + tau_rth_r + tau_thz_z;
    let d_mz = tau_rz / r + tau_rz_r + tau_zz_z;

    // T = p/(R·ρ) and its derivatives.
    let t_r = (p_r / rho - p * rho_r / (rho * rho)) / rg;
    let t_rr = (p_rr / rho - 2.0 * p_r * rho_r / (rho * rho) - p * rho_rr / (rho * rho)
        + 2.0 * p * rho_r * rho_r / (rho * rho * rho))
        / rg;
    let t_zz = (p_zz / rho - 2.0 * p_z * rho_z / (rho * rho) - p * rho_zz / (rho * rho)
        + 2.0 * p * rho_z * rho_z / (rho * rho * rho))
        / rg;

    let g_r = ur * tau_rr + ut * tau_rth + uz * tau_rz + k_gas * t_r;
    let g_r_r = ur_r * tau_rr
        + ur * tau_rr_r
        + ut_r * tau_rth
        + ut * tau_rth_r
        + uz_r * tau_rz
        + uz * tau_rz_r
        + k_gas * t_rr;
    let g_z_z = ur_z * tau_rz
        + ur * tau_rz_z
        + ut_z * tau_thz
        + ut * tau_thz_z
        + uz_z * tau_zz
        + uz * tau_zz_z
        + k_gas * t_zz;
    let d_en = g_r / r + g_r_r + g_z_z;

    let d_rc = rho_d * (d_r(5) / r + d_rr(5) + d_zz(5));

    let _ = c;
    // Burn slot carries no viscous residual (F_visc does not diffuse b).
    [0.0, -d_mr, -d_mt, -d_mz, -d_en, -d_rc, 0.0]
}

#[test]
fn mms_with_all_viscous_terms_recovers_formal_order() {
    let eos = GammaLaw {
        gamma: emms::MMS_GAMMA,
    };
    let run = |n: usize| -> [f64; NCOMP] {
        let spec = GridSpec {
            r_min: emms::MMS_R_MIN,
            dr: 1.0 / n as f64,
            n_r: n,
            z_min: 0.0,
            dz: 1.0 / n as f64,
            n_z: n,
            n_theta_max: 1,
            axisymmetry_assertion: true,
        };
        let mut g = Grid::build(spec, FIELDS).expect("grid");
        let f = EulerFields::resolve(&g).expect("fields");
        fill_from_prim(&mut g, &f, &eos, |r, th, z| {
            emms::manufactured(r, th, z, 0.0, 0.0).w
        });

        let source = move |r: f64, th: f64, z: f64, t: f64| -> Cons {
            let euler = emms::mms_source(r, th, z, t, 0.0, &eos);
            let visc = mms_visc_residual(r, z, t);
            std::array::from_fn(|k| euler[k] + visc[k])
        };
        let exact_bc =
            move |r: f64, th: f64, z: f64, t: f64| emms::manufactured(r, th, z, t, 0.0).w;
        let op = Euler {
            eos,
            source: &source,
            bcs: FlowBcs {
                r_inner: FlowBc::Prescribed(&exact_bc),
                r_outer: FlowBc::Prescribed(&exact_bc),
                z_lo: FlowBc::Prescribed(&exact_bc),
                z_hi: FlowBc::Prescribed(&exact_bc),
            },
            wall_normal: None,
            slip_wall_z_faces: true,
            combustion: None,
        };
        let wall_u = |r: f64, th: f64, z: f64, t: f64| -> (f64, f64, f64) {
            let w = emms::manufactured(r, th, z, t, 0.0).w;
            (w[1], w[2], w[3])
        };
        let wall_t = |r: f64, th: f64, z: f64, t: f64| -> f64 {
            let w = emms::manufactured(r, th, z, t, 0.0).w;
            w[4] / (w[0] * (MMS_CP - MMS_CV))
        };
        let wall_c =
            |r: f64, th: f64, z: f64, t: f64| -> f64 { emms::manufactured(r, th, z, t, 0.0).w[5] };
        let face = || FaceGasBc {
            velocity: VelocityBc::NoSlip(&wall_u),
            thermal: ThermalBc::Isothermal(&wall_t),
            species: SpeciesBc::Prescribed(&wall_c),
        };
        let tr_props = ConstantTransport::new(
            specific_heat_capacity_j_per_kg_k(MMS_CP),
            dynamic_viscosity_pa_s(MMS_MU),
            MMS_PR,
            emms::MMS_GAMMA,
            MMS_SC,
        )
        .expect("transport")
        .into_props();
        let tr_query = spine_query!(tr_props);
        let gas_op = GasDiffusion::new(GasDiffBcs {
            r_inner: face(),
            r_outer: face(),
            z_lo: face(),
            z_hi: face(),
        });
        let mms_temp =
            |w: &Prim| -> Result<f64, &'static str> { Ok(w[4] / (w[0] * (MMS_CP - MMS_CV))) };
        let flow = FlowClass {
            op: &op,
            fields: &f,
        };
        let gas = GasDiffusionClass {
            op: &gas_op,
            temperature: &mms_temp,
            transport: &tr_query,
        };
        let mut sdc = Sdc::new();
        let mut t = 0.0f64;
        while t < emms::MMS_T_FINAL {
            let dt = sdc
                .stable_dt(&g, &flow, 0.4)
                .expect("dt")
                .min(emms::MMS_T_FINAL - t);
            sdc.step(&mut g, Some(&flow), None, Some(&gas), None, None, t, dt)
                .expect("step");
            t += dt;
        }

        // Volume-weighted L1 per conserved component.
        let ids = f.ids();
        let mut num = [0.0f64; NCOMP];
        let mut den = 0.0f64;
        g.for_each_active_cell(|cell| {
            let b = g.brick(cell.bi);
            let vol = g.cell_volume(cell.i_r, b.n_theta());
            let ue = eos.prim_to_cons(
                &emms::manufactured(cell.r, cell.theta, cell.z, emms::MMS_T_FINAL, 0.0).w,
            );
            for k in 0..NCOMP {
                num[k] += vol * (b.field(ids[k])[cell.idx] - ue[k]).abs();
            }
            den += vol;
        });
        std::array::from_fn(|k| num[k] / den)
    };

    let (a, b, c) = (run(16), run(32), run(64));
    for k in 0..NCOMP {
        let o1 = (a[k] / b[k]).log2();
        let o2 = (b[k] / c[k]).log2();
        println!(
            "viscous MMS component {k}: L1 {:.3e} → {:.3e} → {:.3e}, orders {o1:.2}, {o2:.2}",
            a[k], b[k], c[k]
        );
        assert!(
            (1.6..=2.8).contains(&o2),
            "component {k} fine-pair order {o2:.2} outside [1.6, 2.8]"
        );
    }
}

// =============================================================================
// Schedule refusals + determinism.
// =============================================================================

#[test]
fn gas_diffusion_without_flow_or_at_azimuthal_resolution_refuses() {
    // Without the flow class.
    let spec = GridSpec {
        r_min: 0.5,
        dr: 0.1,
        n_r: 8,
        z_min: 0.0,
        dz: 0.1,
        n_z: 8,
        n_theta_max: 1,
        axisymmetry_assertion: true,
    };
    let mut g = Grid::build(spec, FIELDS).expect("grid");
    let tr_props = transport(0.01);
    let tr_query = spine_query!(tr_props);
    let gas_op = GasDiffusion::new(GasDiffBcs {
        r_inner: FaceGasBc::free(),
        r_outer: FaceGasBc::free(),
        z_lo: FaceGasBc::free(),
        z_hi: FaceGasBc::free(),
    });
    let gas = GasDiffusionClass {
        op: &gas_op,
        temperature: &temperature,
        transport: &tr_query,
    };
    match Sdc::new().step::<GammaLaw>(&mut g, None, None, Some(&gas), None, None, 0.0, 1e-6) {
        Err(SdcError::Config(_)) => {}
        other => panic!("expected a Config refusal, got {other:?}"),
    }

    // At N_θ > 1 (plan S8 owns the θ re-keying).
    let spec_theta = GridSpec {
        r_min: 0.5,
        dr: 0.1,
        n_r: 8,
        z_min: 0.0,
        dz: 0.1,
        n_z: 8,
        n_theta_max: 4,
        axisymmetry_assertion: false,
    };
    let mut g4 = Grid::build(spec_theta, FIELDS).expect("grid");
    let f4 = EulerFields::resolve(&g4).expect("fields");
    let eos = GammaLaw { gamma: GAMMA };
    fill_from_prim(&mut g4, &f4, &eos, |_, _, _| {
        prim6(1.0, 0.0, 0.0, 0.0, 1.0, 0.5)
    });
    let op = Euler {
        eos,
        source: &ZERO_SRC,
        bcs: FlowBcs {
            r_inner: FlowBc::Reflecting,
            r_outer: FlowBc::Reflecting,
            z_lo: FlowBc::Reflecting,
            z_hi: FlowBc::Reflecting,
        },
        wall_normal: None,
        slip_wall_z_faces: true,
        combustion: None,
    };
    let flow = FlowClass {
        op: &op,
        fields: &f4,
    };
    match Sdc::new().step(
        &mut g4,
        Some(&flow),
        None,
        Some(&gas),
        None,
        None,
        0.0,
        1e-9,
    ) {
        Err(SdcError::Config(m)) => assert!(m.contains("N_θ"), "wrong refusal: {m}"),
        other => panic!("expected the N_θ refusal, got {other:?}"),
    }
}

/// META-1 §2: the FULL step with the gas class scheduled (Euler rayon
/// sweeps + per-component CG + cross-term Picard + audit reductions) is
/// bit-identical at any thread count.
#[test]
fn coupled_step_with_gas_diffusion_is_bit_identical_across_threads() {
    const MU: f64 = 3.2;
    const U_WALL: f64 = 100.0;
    let run = |threads: usize| -> Vec<u64> {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .expect("pool");
        pool.install(|| {
            let spec = GridSpec {
                r_min: 5.0e-3,
                dr: 1.0e-4,
                n_r: 16,
                z_min: 0.0,
                dz: 2.0e-4,
                n_z: 8,
                n_theta_max: 1,
                axisymmetry_assertion: true,
            };
            let mut g = Grid::build(spec, FIELDS).expect("grid");
            let f = EulerFields::resolve(&g).expect("fields");
            let eos = GammaLaw { gamma: GAMMA };
            let rho0 = P0 / (R_GAS * T0);
            fill_from_prim(&mut g, &f, &eos, |_, _, _| {
                prim6(rho0, 0.0, 0.0, 0.0, P0, 0.5)
            });
            let op = Euler {
                eos,
                source: &ZERO_SRC,
                bcs: FlowBcs {
                    r_inner: FlowBc::Reflecting,
                    r_outer: FlowBc::Reflecting,
                    z_lo: FlowBc::Transmissive,
                    z_hi: FlowBc::Transmissive,
                },
                wall_normal: None,
                slip_wall_z_faces: true,
                combustion: None,
            };
            let still = |_: f64, _: f64, _: f64, _: f64| (0.0, 0.0, 0.0);
            let t_ramp = 5.0 * (6.6e-3f64 - 5.0e-3).powi(2) / (MU / (P0 / (R_GAS * T0)));
            let slide =
                move |_: f64, _: f64, _: f64, t: f64| (0.0, 0.0, U_WALL * (t / t_ramp).min(1.0));
            let wall_t = |_: f64, _: f64, _: f64, _: f64| T0;
            let tr_props = transport(MU);
            let tr_query = spine_query!(tr_props);
            let gas_op = GasDiffusion::new(GasDiffBcs {
                r_inner: FaceGasBc {
                    velocity: VelocityBc::NoSlip(&still),
                    thermal: ThermalBc::Adiabatic,
                    species: SpeciesBc::ZeroFlux,
                },
                r_outer: FaceGasBc {
                    velocity: VelocityBc::NoSlip(&slide),
                    thermal: ThermalBc::Isothermal(&wall_t),
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
            });
            let flow = FlowClass {
                op: &op,
                fields: &f,
            };
            let gas = GasDiffusionClass {
                op: &gas_op,
                temperature: &temperature,
                transport: &tr_query,
            };
            let mut sdc = Sdc::new();
            let mut t = 0.0f64;
            for _ in 0..150 {
                let dt = sdc.stable_dt(&g, &flow, 0.4).expect("dt");
                sdc.step(&mut g, Some(&flow), None, Some(&gas), None, None, t, dt)
                    .expect("step");
                t += dt;
            }
            let mut bits = Vec::new();
            for b in g.bricks() {
                for id in f.ids() {
                    bits.extend(b.field(id).iter().map(|v| v.to_bits()));
                }
            }
            bits
        })
    };
    let one = run(1);
    let four = run(4);
    assert_eq!(one.len(), four.len());
    let diff = one.iter().zip(&four).filter(|(a, b)| a != b).count();
    assert_eq!(
        diff, 0,
        "{diff} values differ between 1 and 4 threads through the gas-diffusion step"
    );
}

// =============================================================================
// The S4 configuration in miniature: ALL FOUR classes in one march — flow +
// gas diffusion + solid conduction + Robin-Robin exchange. Hot quiescent
// gas in an annular duct with a coolant-backed solid liner: the audit's
// combined energy row (gas-diff ledger + exchange debit + solid ledger)
// must close on every step, the wall function must own the wall (gas-side
// F_visc suppressed there — the unit battery proves the rates; here the
// coupled march proves the bookkeeping), and both fixed-sweep acceptances
// must hold live.
// =============================================================================

#[test]
fn four_class_coupled_march_audits_closed_with_wall_ownership() {
    // Air-like transport: a big μ here means a metal-like GAS conductivity
    // through k = μc_p/Pr, whose wall-law h flash-swings the debited wall
    // cell hundreds of K per step — the wall-heat registry caps μ at 1e-2
    // for exactly this reason. The stiff-diffusion regime is the
    // Poiseuille fixture's job; THIS fixture is the four-class bookkeeping.
    const MU: f64 = 5.0e-5;
    const T_GAS: f64 = 600.0;
    let spec = GridSpec {
        r_min: 1.0e-3,
        dr: 1.0e-4,
        n_r: 16,
        z_min: 0.0,
        dz: 2.0e-4,
        n_z: 8,
        n_theta_max: 1,
        axisymmetry_assertion: true,
    };
    let mut fields: Vec<&str> = FIELDS.to_vec();
    fields.push("T_solid");
    fields.push("rate_solid");
    let mut g = Grid::build_with_regions(spec, &fields, |i_r, _| {
        if i_r >= 12 {
            Region::Solid
        } else {
            Region::Gas
        }
    })
    .expect("region world");
    let f = EulerFields::resolve(&g).expect("fields");
    let t_solid = g.field_id("T_solid").expect("field");
    let rate_solid = g.field_id("rate_solid").expect("field");
    let eos = GammaLaw { gamma: GAMMA };
    let rho_hot = P0 / (R_GAS * T_GAS);
    fill_from_prim(&mut g, &f, &eos, |_, _, _| {
        prim6(rho_hot, 0.0, 0.0, 0.0, P0, 0.5)
    });
    fill_solid(&mut g, t_solid, T0); // fill_field covers gas cells only

    let op = Euler {
        eos,
        source: &ZERO_SRC,
        bcs: FlowBcs {
            r_inner: FlowBc::Reflecting,
            r_outer: FlowBc::Reflecting, // unreached: the liner interposes
            z_lo: FlowBc::Reflecting,
            z_hi: FlowBc::Reflecting,
        },
        wall_normal: None,
        slip_wall_z_faces: true,
        combustion: None,
    };
    let tr_props = transport(MU);
    let tr_query = spine_query!(tr_props);
    let gas_op = GasDiffusion::new(GasDiffBcs {
        r_inner: FaceGasBc::free(),
        r_outer: FaceGasBc::free(), // unreached: wall-law faces interpose
        z_lo: FaceGasBc::free(),
        z_hi: FaceGasBc::free(),
    });
    let zero_heat = |_: f64, _: f64, _: f64, _: f64| 0.0;
    let solid_op = Conduction {
        kappa: 20.0,
        rho_cp: 4.0e6,
        source: &zero_heat,
        domain: Domain::Solid,
        interior: InteriorFaces::refuse(),
        bcs: Bcs {
            r_inner: FaceBc::HeatFlux(0.0), // unreached: gas rings interpose
            r_outer: FaceBc::Robin {
                h: 5.0e3,
                t_inf: T0,
            },
            z_lo: FaceBc::HeatFlux(0.0),
            z_hi: FaceBc::HeatFlux(0.0),
        },
    };
    let patches = build_wall_patches(&g).expect("patches");
    assert!(!patches.is_empty(), "the duct must have wall faces");
    let flow = FlowClass {
        op: &op,
        fields: &f,
    };
    let gas = GasDiffusionClass {
        op: &gas_op,
        temperature: &temperature,
        transport: &tr_query,
    };
    let diffusion = DiffusionClass {
        op: &solid_op,
        t_field: t_solid,
        scratch_field: rate_solid,
    };
    let exchange = ExchangeClass {
        patches: &patches,
        law: &WallLaw::new(),
        temperature: &temperature,
        transport: &tr_query,
    };

    let mut sdc = Sdc::new();
    let mut t = 0.0f64;
    let mut jacket_last = 0.0f64;
    for _ in 0..200 {
        let dt = sdc.stable_dt(&g, &flow, 0.4).expect("dt");
        let report = sdc
            .step(
                &mut g,
                Some(&flow),
                Some(&diffusion),
                Some(&gas),
                Some(&exchange),
                None,
                t,
                dt,
            )
            .expect("four-class step (audit armed)");
        t += dt;
        let ex = report.exchange.expect("exchange scheduled");
        jacket_last = ex.jacket_w;
        assert!(ex.robin_resid <= EPS_ROBIN_RESID);
        assert!(report.gas_picard_resid <= EPS_GAS_DIFF_RESID);
        // The combined energy row closed (the step would have halted
        // otherwise); read it out once for the record.
        let row = report
            .audit
            .iter()
            .find(|r| r.quantity == "energy")
            .expect("energy row");
        assert!((row.delta - row.applied).abs() <= row.tol);
    }
    // Hot gas drives heat INTO the liner through the one wall law.
    assert!(
        jacket_last > 0.0,
        "hot gas must heat the liner: jacket {jacket_last:.3e} W"
    );
    // The liner warmed from its 300 K start.
    let mut t_max = 0.0f64;
    for bi in 0..g.n_bricks() {
        let b = g.brick(bi);
        for local in 0..64 {
            if b.solid_mask() & (1u64 << local) != 0 {
                t_max = t_max.max(b.field(t_solid)[b.cell_index(0, local)]);
            }
        }
    }
    println!("four-class duct: jacket {jacket_last:.3e} W, liner T_max {t_max:.2} K");
    assert!(t_max > T0 + 1.0e-3, "liner never warmed: {t_max} K");
}
