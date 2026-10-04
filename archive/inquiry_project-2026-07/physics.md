# Physics Simulation Components

Comprehensive breakdown of every physics domain the simulator must handle. Each component is categorized by computational tier, evaluated for existing tool integration, and scoped for build effort.

**Coordinate system:** Galactocentric double-double (sub-nm precision at 30 kpc). Ship octree in ship-local f64. See [determinism.md](determinism.md) for the full specification of coordinate precision, floating-point rules, and the physics/display path separation that guarantees bit-exact reproducibility.

**Trajectory integration:** Relativistic momentum formulation in the galactocentric frame. State = (position, momentum) where p = γmv. The integrator (DOP853) handles all speeds from rest to 0.99c with the same equations. Near compact objects (< 100 r_s), switches to Kerr geodesic equations in BH-centered Boyer-Lindquist coordinates. See [determinism.md](determinism.md) for regime transition hysteresis.

---

## COMPUTATIONAL TIERS

The simulation workload divides into four tiers based on when computation happens and how expensive it is:

| Tier | When | Budget | Examples |
|---|---|---|---|
| **Precomputed** | Offline, once per ship | Minutes to hours | Reaction tables, aero coefficients, shielding tables |
| **Per-Tick Core** | Every integration step | Microseconds | Gravity, trajectory, resource bookkeeping |
| **Per-Tick Conditional** | When regime applies | Microseconds-milliseconds | Atmospheric drag, relativistic corrections, GR geodesics |
| **On-Demand** | When queried | Milliseconds | Star lookups, procedural generation, minor body positions |

The tick rate target is ~1000 Hz for numerical stability during atmospheric flight, dropping to ~100 Hz in vacuum and ~10 Hz during interstellar cruise (where changes are slow). Adaptive timestep controls this automatically.

---

## 1. SHIP INTERNAL PHYSICS

The ship is a directed graph of major components (engines, tanks, reactors, pods, radiators, trusses) connected by structural links and resource lines. Each component has sub-components (valves, pumps, computers, heat exchangers) that are modeled abstractly.

---

### 1.1 Chemical Reactions

**What:** Combustion in engines (LOX/RP-1, LOX/LH2, LOX/CH4), CO2 scrubbing and O2 generation in life support, fuel cell chemistry for power, and detection of unintended reactions (propellant leaks, thermal decomposition).

**Tier:** Precomputed tables, queried per-tick.

**Approach:** Run full chemical kinetics offline across the relevant parameter space (chamber pressure, mixture ratio, inlet temperature) and store results as lookup tables. At runtime, interpolate the tables given current conditions.

**Tool evaluation:**

| Tool | What it does | License | Verdict |
|---|---|---|---|
| **Cantera** | Chemical kinetics, thermodynamics, transport. Equilibrium and time-dependent chemistry for arbitrary reaction mechanisms. | BSD-3 | **Integrate.** Use in the offline characterization tool to compute combustion products, flame temperatures, specific impulse, and exhaust composition across parameter spaces. Also use for life support chemistry (CO2 + LiOH, electrolysis, Sabatier reaction). Cantera is the standard tool for this and does exactly what we need. |
| **CoolProp** | Thermophysical properties of fluids (density, viscosity, specific heat vs. temperature and pressure). | MIT | **Integrate.** Use alongside Cantera for propellant properties at non-standard conditions. Needed for modeling propellant behavior in tanks (boiloff, pressurization, slosh thermal effects). |
| Custom | Unintended reaction detection. | -- | **Build.** Cantera could in principle do this, but runtime detection of "leak + heat = fire" is a logic problem, not a kinetics problem. A rule-based system checking temperature, pressure, and chemical proximity is simpler and faster than running Cantera at runtime. |

**Runtime cost:** Table interpolation only. Negligible.

---

### 1.2 Nuclear Reactions

**What:** Fission reactor dynamics (neutron flux, criticality, decay heat, fission product buildup), fusion reactions (D-T, D-He3 -- if future ships use fusion propulsion), antimatter annihilation (pion/muon cascades, energy deposition), and radioactive decay of activated materials.

**Tier:** Precomputed for reaction physics. Per-tick for reactor state (simplified neutronics model tracking power level, temperature feedback, control rod position).

**Tool evaluation:**

| Tool | What it does | License | Verdict |
|---|---|---|---|
| **OpenMC** | Monte Carlo neutron/photon transport. Full reactor physics. | MIT | **Use offline only.** Run OpenMC to precompute reactor behavior tables: power output vs. control rod position, temperature coefficients of reactivity, decay heat curves, neutron flux distributions for radiation exposure. Far too expensive for runtime. |
| **Geant4** | General particle physics simulation. Handles hadronic cascades, electromagnetic showers, decay chains. | Custom (free) | **Use offline only.** For antimatter annihilation energy deposition profiles and activated material decay rates. Also needed for the ISM cascade research (Section 4.1). |
| Custom | Runtime reactor model. | -- | **Build.** A point-kinetics model with temperature feedback (6 delayed neutron groups, Doppler coefficient, coolant temperature feedback) is standard nuclear engineering and runs in microseconds. This is well-understood math, not a research problem. |

**Runtime cost:** Point-kinetics ODE (6 variables) per reactor per tick. Negligible.

---

### 1.3 Thermal Management

**What:** Temperature of each major component. Fine-grained thermal zones within engines and reactors (combustion chamber, nozzle throat, turbopump, coolant channels). Heat transfer between components via conduction (structural connections), radiation (all components), and convection (atmospheric flight, internal coolant loops). External heat sources: solar flux, atmospheric heating, engine exhaust impingement, ISM particle heating.

**Tier:** Per-tick core. Thermal state is always evolving.

**Approach:** Thermal network model (resistor-capacitor analogy). Each component or sub-component is a thermal node with a heat capacity. Connections between nodes have thermal resistances. External heat loads are source terms. The system is a set of coupled ODEs: dT_i/dt = (1/C_i) * sum(Q_ij).

**Tool evaluation:**

| Tool | What it does | License | Verdict |
|---|---|---|---|
| **SfePy** / **FEniCS** | Finite element frameworks. Full 3D heat transfer with arbitrary geometry. | BSD / LGPL | **Use offline only.** Run FEA to compute thermal resistances between components and to validate the thermal network model against full 3D solutions. Too expensive for runtime. |
| Custom | Thermal network solver. | -- | **Build (Rust).** The thermal network is a sparse linear system that updates every tick. For ~50-200 thermal nodes per ship, this is a small LU solve — custom Rust implementation (thermal.md). scipy.sparse serves as the OFFLINE validation reference only; no Python at runtime. |

**Runtime cost:** Sparse matrix solve, ~50-200 nodes. Microseconds.

---

### 1.4 Resource Management

**What:** Mass and volume tracking for all consumables: propellant (fuel + oxidizer, by tank), reactor fuel, coolant, breathable air, water, food, waste, electrical power (generation, storage, consumption). Flow rates through resource lines (propellant feed, coolant loops, air circulation). Storage configuration (which tanks hold what, pressure, temperature, remaining capacity).

**Tier:** Per-tick core. Resources change every tick during powered flight.

**Approach:** Directed graph of storage nodes and flow edges. Each edge has a valve state (open/closed/throttled), maximum flow rate, and current flow. Conservation of mass enforced at every node. Power is tracked as generation minus load with battery state-of-charge.

**Tool evaluation:**

| Tool | What it does | Verdict |
|---|---|---|
| Custom | Resource flow graph. | **Build.** This is bookkeeping, not physics. A directed graph with flow constraints, mass conservation, and simple valve logic. No existing library fits because the abstraction level (ship-scale plumbing, not pipe-segment CFD) is specific to this project. Clean, simple code. |

**Runtime cost:** Graph traversal + mass updates, ~20-100 nodes. Negligible.

---

### 1.5 Structural Mechanics

**What:** Load paths through the ship under thrust, rotation, tidal forces, aerodynamic forces, and internal pressure. Structural connections between major components (how the engine is attached to the fuel tank, how the fuel tank is attached to the command pod). Failure detection: when loads exceed material limits (buckling, tensile failure, fatigue accumulation).

**Tier:** Precomputed for structural mode shapes and stiffness. Per-tick for load evaluation and failure checking.

**Approach:** Reduced-order structural model. Offline FEA computes the stiffness matrix and failure envelopes for the assembled ship. At runtime, given the current force vector (thrust + gravity gradient + aero + rotation), compute stress at critical joints by multiplying the precomputed compliance matrix. Compare against failure criteria.

**Tool evaluation:**

| Tool | What it does | License | Verdict |
|---|---|---|---|
| **FEniCS** | Finite element framework. Linear/nonlinear structural analysis. | LGPL | **Use offline.** Compute stiffness matrices, mode shapes, and failure envelopes for each ship configuration. Export reduced-order model for runtime use. |
| Custom | Runtime stress evaluator. | -- | **Build.** Matrix-vector multiply with precomputed stiffness data. Check results against precomputed failure surfaces. Standard structural engineering approach. |

**Runtime cost:** Small matrix multiply (~10-50 DOF reduced model) + comparison. Negligible.

---

### 1.6 Atmosphere and Life Support

**What:** Gas composition inside pressurized volumes (O2, N2, CO2, H2O vapor, trace contaminants like CO, NH3). Pressure regulation. Temperature regulation (HVAC). CO2 removal rate, O2 generation rate, humidity control. Leak modeling (hole size to flow rate to pressure decay).

**Tier:** Per-tick core for crewed ships. Skipped entirely for uncrewed vehicles.

**Approach:** Coupled ODE system. Each gas species has a production rate (crew metabolism, outgassing, leaks) and removal rate (scrubbers, vents, recyclers). Pressure follows from ideal gas law given total moles and volume. Temperature couples to the thermal network (Section 1.3).

**Tool evaluation:**

| Tool | What it does | Verdict |
|---|---|---|
| Custom | Atmosphere ODE system. | **Build.** This is a small system of ODEs (~10 state variables per pressurized volume). No existing library targets this specific abstraction level. The chemistry of CO2 scrubbing is simple enough to hard-code (or precompute with Cantera if multiple scrubber chemistries are needed). Standard ODE integration. |

**Runtime cost:** ~10-variable ODE per pressurized volume. Negligible.

---

### 1.7 Radiation Exposure and Degradation

**What:** Radiation field at each component's location (from reactors, cosmic rays, solar particle events, trapped radiation belts, ISM at relativistic speeds). Shielding calculation: how much radiation penetrates to each component given intervening material. Dose accumulation for crew (Sievert tracking) and electronics (total ionizing dose). Material degradation: embrittlement, insulation breakdown, solar cell efficiency loss, optics darkening.

**Tier:** Precomputed for shielding tables and degradation rates. Per-tick for dose accumulation and degradation tracking.

**Approach:** Offline: run Geant4 or OpenMC to compute dose transmission factors through the ship's shielding geometry for various radiation spectra (reactor, GCR, SPE, trapped belts, ISM). Store as lookup tables indexed by radiation type, energy spectrum, and incidence angle. Runtime: given the external radiation environment (from Section 3), look up transmission factors and accumulate dose. Degradation is a monotonic function of accumulated dose + thermal cycling count + kinetic impact history.

**Tool evaluation:**

| Tool | What it does | License | Verdict |
|---|---|---|---|
| **Geant4** | Full Monte Carlo particle transport through arbitrary geometry. | Custom (free) | **Use offline.** The gold standard for radiation transport. Compute shielding effectiveness for every radiation type and angle. Export as tables. |
| **OpenMC** | Monte Carlo neutron/photon transport, optimized for reactor geometries. | MIT | **Use offline.** Specifically for reactor radiation shielding. |
| Custom | Runtime dose accumulator. | -- | **Build.** Table lookup + integration. Track cumulative dose per component, apply degradation curves. Simple. |

**Runtime cost:** Table lookup + scalar accumulation per component. Negligible.

---

### 1.8 Electrical and Control Systems

**What:** Power bus voltage and current. Generator output (solar panels, RTGs, reactors, fuel cells). Battery state of charge. Load balancing across systems. Valve positions and actuator states. Computer/autopilot operational status. Sensor readings derived from physics state (accelerometers from force/mass, thermocouples from thermal nodes, radiation monitors from dose rate).

**Tier:** Per-tick core.

**Approach:** Power is modeled as a simple bus: total generation minus total load equals net power to/from batteries. No AC/DC distinction, no impedance matching. Valve and actuator states are discrete (open/closed/position). Sensor readings are derived quantities, not independently simulated -- the accelerometer reads whatever the physics engine computed for acceleration.

**Tool evaluation:**

| Tool | What it does | Verdict |
|---|---|---|
| Custom | Power balance + control state machine. | **Build.** This is straightforward engineering bookkeeping. No existing tool matches the abstraction level. A few dozen lines of code per subsystem. |

**Runtime cost:** Arithmetic. Negligible.

---

### 1.9 Pressure

**What:** Internal pressure of all pressurized volumes (crew cabin, propellant tanks, pressurant bottles, reactor coolant loops). External pressure (atmospheric during launch/reentry, zero in vacuum). Differential pressure across structural boundaries. Tank pressure as function of propellant mass, temperature, and ullage volume.

**Tier:** Per-tick core.

**Approach:** Ideal gas law for gaseous volumes. Propellant tanks use real-gas EOS (via CoolProp) precomputed into tables indexed by temperature and fill level. Pressure differential drives leak flow rates (Section 1.6) and contributes to structural loads (Section 1.5).

**Tool evaluation:**

| Tool | What it does | Verdict |
|---|---|---|
| **CoolProp** | Real-gas equations of state. | **Integrate (offline).** Precompute pressure-temperature-density tables for propellants. |
| Custom | Runtime pressure tracking. | **Build.** Table lookup + ideal gas law. Couples into atmosphere (1.6) and structural (1.5) modules. |

**Runtime cost:** Negligible.

---

## 2. ORBITAL AND TRAJECTORY PHYSICS

---

### 2.1 Gravity on the Ship

**What:** The total gravitational acceleration on the ship, computed uniformly everywhere from two sources:

1. **The galactic mean field** -- the self-consistent gravitational potential of the entire Milky Way, including all mass concentrations (disk, bulge, bar, spiral arms, halo, clusters, nuclear structures). See Section 3 for the complete field definition.

2. **Individual point sources above the adaptive influence threshold** -- every body (star, planet, moon, black hole) whose individual gravitational acceleration on the ship exceeds ε × |a_total|, where a_total is the ship's current total acceleration and ε ~ 10^-6.

There is one force computation rule applied everywhere: F_gravity = -∇Φ_mean_field + Σ [(GM_i / r_i^2) r̂_i - ∇Φ_smooth_i] for all active point sources. The second term in the sum subtracts each active point source's smoothed contribution to the mean field, preventing double-counting. For most stars this correction is negligible (one star out of ~10¹¹ contributes ~10^-11 of the smooth field). For Sgr A* it is critical: the bulge potential includes Sgr A*'s mass smoothly, so promoting Sgr A* to a point source without subtracting its smooth contribution would double-count ~4 million solar masses, producing ~100% gravity error at the galactic center.

In practice: for individual stars, skip the subtraction (error < 10^-10). For Sgr A*, precompute the smooth bulge contribution at the BH's location and subtract it when Sgr A* is active. For massive globular clusters (Omega Centauri at 3.5 × 10^6 M_sun), apply the same correction.

No regimes, no hierarchy, no switching between gravity models. The set of active point sources changes as the ship moves (bodies cross the threshold continuously), but the formula never changes.

**Tier:** Per-tick core. Gravity is always on.

**What affects what:**

| Object | Feels gravity from |
|---|---|
| **Ship** | Mean field + all point sources above threshold (stars, planets, Sgr A*) |
| **Stars** | Mean field + Sgr A* + bound companions (binary/triple members from same seed). No gravity from unrelated stars. |
| **Planets/moons** | Motion is ANALYTIC — Kepler + first-order secular theory (generation.md Section 4.10); sibling interactions enter through the precomputed secular rates, never runtime N-body (determinism: position is a pure function of (elements, T)). They EXERT gravity on the ship as point sources like everything else. |
| **Clusters** | Their mass is part of the mean field; cluster-scale density features orbit within the field |
| **Mean field** | Hybrid analytic potential-density pairs pinned to observed parameters + precomputed LMC-response perturbation (Section 3.1); bar/spiral rotate rigidly |

**Adaptive influence threshold:** a_threshold = ε × |a_total|. With ε = 10^-6:
- LEO: threshold ~9 × 10^-6 m/s^2. Active: Earth, Moon, Sun, Jupiter, Saturn, Venus, Mars (~8 bodies).
- Interplanetary: threshold ~6 × 10^-9 m/s^2. Active: ~5-8 planets.
- Solar neighborhood: threshold ~2 × 10^-16 m/s^2. Active: ~1,000-2,500 nearby stars (influence radius ~17 pc for a 0.4 M_sun star at 0.1 stars/pc³ — note the 1024-source hard cap BINDS in this regime; the cap's ranked truncation is the effective threshold).
- Galactic center (1 pc from Sgr A*): threshold ~5.6 × 10^-13 m/s^2. Active: Sgr A* + ~100-500 nearby stars.

Hard cap at 1024 active sources. If more bodies exceed threshold (extreme stellar density), rank by acceleration magnitude, take the top 1024.

**Solar system bodies:** The Sun, 8 planets, and ~15-20 major moons are always active point sources when the ship is within the solar system. Positions come from SPICE (DE441) inside its ±~14 kyr validity span, blended C¹ to the Kepler + secular representation beyond it (generation.md Section 4.10) — positions are pure functions of T at all epochs, with no runtime integration of solar-system bodies. Perturbations FELT BY THE SHIP include J2 oblateness for Earth, Jupiter, and Saturn, plus solar radiation pressure. Offline N-body (REBOUND) is used for validation of the representation only.

**Tool evaluation:**

| Tool | What it does | License | Verdict |
|---|---|---|---|
| **SpiceyPy / SPICE** | JPL ephemerides. Precise positions of solar system bodies at any time. | Public domain | **Integrate.** Initialization and validation of solar system body positions. |
| Custom | Adaptive N-body integrator (DOP853). | -- | **Build.** ~200 lines. The force assembly combines mean field gradient + point source sum. Validate against SPICE for solar system, against known Mercury precession for PN corrections. |
| **REBOUND** | N-body integrator. IAS15, WHFast. | GPL-3 | **Validation only.** GPL-3 prevents runtime linking. Use for offline trajectory validation. |

**Runtime cost:** Mean field evaluation (~1-3 μs, Section 3.1) + point source sum (100-1024 sources × ~10 ns = 1-10 μs) = ~2-13 μs per tick. Well within budget.

---

### 2.2 Trajectory Integration

**What:** Assemble all forces on the ship (gravity from 2.1, thrust from ship tables, aerodynamic forces from 2.3, radiation pressure, tidal forces near compact objects) and integrate the equations of motion. Adaptive timestep: small steps during atmospheric flight or close approaches, large steps during interstellar cruise.

**Tier:** Per-tick core. This is the central computation.

**Approach:** State vector: position (3), velocity (3), proper time (1), mass (1) = 8 variables minimum. Add orientation (quaternion, 4) and angular velocity (3) for attitude dynamics = 15 variables. Integrate with an embedded Runge-Kutta method (DOP853) that provides error estimation for adaptive stepping.

**Tool evaluation:**

| Tool | What it does | License | Verdict |
|---|---|---|---|
| Custom | DOP853 integrator (Rust). | -- | **Build.** Determinism requires full control: libm transcendentals, fixed evaluation order, double-double position accumulation — a wrapped scipy stepper can guarantee none of these, and the runtime is Rust (SCOPE.md). ~300 lines from the published Hairer-Nørsett-Wanner coefficients. scipy's `solve_ivp` (DOP853) is the OFFLINE cross-validation reference. |
| Custom | Force assembler + regime switcher. | -- | **Build.** The integration algorithm itself is standard. The hard part is the force assembly: which forces are active, which corrections apply, when to switch regimes. This is the connective tissue that makes the simulator novel. |

**Runtime cost:** 8-15 function evaluations per step (for DOP853). Each evaluation queries gravity, thrust tables, and any active conditional physics. Total: tens of microseconds per step.

---

### 2.3 Atmospheric Models

**What:** Density, temperature, pressure, and composition of planetary atmospheres as a function of altitude, latitude, time of day, solar activity. Needed for drag computation, heating, and aerodynamic forces during atmospheric flight.

**Tier:** Per-tick conditional. Only active during atmospheric flight.

**Tool evaluation:**

| Tool | What it does | License | Verdict |
|---|---|---|---|
| **NRLMSISE-00** (C library, FFI) | Earth atmosphere model. Standard model used by NASA and NOAA. Density, temperature, composition from 0-1000 km. Depends on solar activity indices (F10.7, Ap). | Public domain | **Integrate via Rust FFI** (one of the two C dependencies, per SCOPE.md). DETERMINISM: the driver indices F10.7/Ap must be pure functions of T — live space-weather data does not exist at T = +1 Myr. Use a synthetic solar-activity model: 11.04-yr cycle, F10.7 oscillating 70-220 sfu, plus seeded short-term fluctuations from the Sun's activity stream (same hierarchical bucket seeding as flares, generation.md 3.4); Ap derived from F10.7 with seeded storms. Validation runs against historical flights may substitute the recorded indices for that epoch. |
| **Mars-GRAM** | NASA's Mars atmosphere model. | Export-controlled | **Cannot use directly.** Mars-GRAM requires a NASA license. Implement a simplified Mars atmosphere (exponential decay with dust storm modeling) calibrated against published Mars-GRAM outputs. |
| Custom | Other planetary atmospheres. | -- | **Build.** Simple exponential or barometric models for Venus, Titan, gas giant upper atmospheres. Calibrate against published atmospheric profiles. Small amount of code per body. |

**Runtime cost:** One function evaluation per tick during atmospheric flight. Microseconds.

---

### 2.4 Special Relativity

**What:** Lorentz factor computation. Relativistic momentum and energy. Time dilation (proper time vs. coordinate time). Relativistic velocity addition. Stellar aberration (apparent star positions shift at high speed). Relativistic Doppler shift (communication frequencies, visual appearance). Relativistic beaming of ISM flux.

**Tier:** Per-tick conditional. Activates when v/c > threshold (e.g., 0.01c where relativistic effects exceed 0.005%).

**Approach:** All of this is direct application of special relativity formulas. No approximation needed -- the exact Lorentz transforms are simple closed-form expressions. The challenge is bookkeeping: maintaining both proper time and coordinate time, transforming between ship frame and coordinate frame consistently.

**Tool evaluation:**

| Tool | What it does | Verdict |
|---|---|---|
| Custom | SR kinematics module. | **Build.** This is textbook physics implemented as straightforward math functions. Lorentz factor, four-velocity, aberration formula, Doppler formula. No existing library needed -- importing one would add complexity without benefit. ~100 lines of code. Validate against known results (muon lifetime, twin paradox, published relativistic trajectory calculations). |

**Runtime cost:** A handful of floating-point operations. Negligible.

---

### 2.5 Relativistic Gravity in the Weak Field (all speeds)

**What:** Corrections to Newtonian gravity below the full-GR regime. Two distinct effects with two distinct validity domains — a fast ship in a weak field, and a slow ship near a massive body. A slow-motion PN expansion (EIH) alone is INVALID for the first case: at 0.3–0.99c the (v/c)² "corrections" are order unity, and plain −γm∇Φ underestimates ultrarelativistic deflection by up to 2× (the light-bending factor).

**Part 1 — Velocity-complete weak-field force law (always active).** The trajectory uses the geodesic equation of the static weak-field metric (g_00 = −(1+2Φ/c²+2Φ²/c⁴), g_ij = (1−2Φ/c²)δ_ij), exact in velocity to first order in Φ/c²:

- State: canonical momentum p (what the integrator already carries)
- Force: dp/dt = −γm(1 + β²)∇Φ  (+ thrust + drag + ...)
- Velocity map: dx/dt = (pc²/E)(1 + 2Φ/c²),  E = √(m²c⁴ + p²c²)

where Φ is the SAME total potential as everywhere else (mean field + active point sources − double-count corrections). One formula, all speeds. Limits: β → 0 recovers Newton exactly; β → 1 gives photon deflection 4GM/bc² and the coordinate-light-speed reduction responsible for Shapiro delay. Cost: one extra multiply per force evaluation — free.

**Part 2 — Static nonlinearity + moving-source terms (per-tick conditional, slow regime).** For precision work near massive bodies: (a) the Φ²/c⁴ term of g_00 — without it, perihelion precession comes out 4/3 of the GR value, so it is required for the Mercury 43″/century validation; (b) EIH cross-terms for source velocities (planets move) and gravitomagnetic frame-dragging (Lense-Thirring) where spin data exists, per IAU/IERS conventions and Soffel et al. (2003). Active when GM/(rc²) > 10⁻¹⁰ for any nearby body. Large-β × strong-field cross-terms are not needed: they only matter within ~100 r_s of compact objects, where full Kerr geodesic integration takes over (Section 2.6).

**Derivation note:** the canonical-momentum form of Part 1 follows from the test-particle Lagrangian L = −mc²√((1+2Φ/c²) − (1−2Φ/c²)v²/c²); a short write-up with limits and error terms is a deliverable of the trajectory-module validation (candidate research note — a velocity-complete bridging formulation is exactly the "regime bridging" claim of this project).

**Tool evaluation:**

| Tool | What it does | Verdict |
|---|---|---|
| Custom | Weak-field force law + velocity map. | **Build.** ~30 lines on top of the existing momentum integrator. Unit tests: Newtonian limit, photon-deflection limit (4GM/bc² against analytic), Shapiro delay. |
| Custom | Φ² term + EIH source-motion terms. | **Build.** ~50 lines per gravitating body. Validate against Mercury's perihelion precession (43"/century) and published post-Newtonian ephemerides. |

**Runtime cost:** Part 1 is ~free (one factor per evaluation). Part 2 is ~2x the Newtonian evaluation for active bodies. Still microseconds.

---

### 2.6 General Relativity (Compact Objects)

**What:** Full GR geodesic integration near black holes and neutron stars. Schwarzschild metric for non-rotating compact objects, Kerr metric for rotating ones. Tidal forces (Riemann tensor components in the ship's frame). Gravitational time dilation (extreme regime). Frame dragging. Gravitational lensing (for telemetry display of background star positions). Innermost stable circular orbit (ISCO), photon sphere, ergosphere geometry.

**Tier:** Per-tick conditional. Active only near compact objects (within ~100 Schwarzschild radii).

**Approach:** Replace Newtonian + PN gravity with geodesic integration in the Kerr metric. The equations of motion are the geodesic equations with a forcing term (ship thrust, expressed as a four-acceleration). This requires:
1. Metric tensor components at the ship's position
2. Christoffel symbols (or their equivalent in the chosen coordinate system)
3. Integration of the geodesic deviation equation for tidal forces
4. Coordinate transformation between Boyer-Lindquist and ship-local frames

**Tool evaluation:**

| Tool | What it does | License | Verdict |
|---|---|---|---|
| **EinsteinPy** | GR in Python. Geodesic integration for Schwarzschild and Kerr metrics. Christoffel symbol computation. Coordinate transforms. | MIT | **Evaluate for validation, likely build custom.** EinsteinPy can integrate geodesics, but it's designed for visualization and education, not for coupling into a flight simulator with thrust and ship physics. It may be too slow for per-tick use (symbolic computation under the hood). Use it to validate our custom geodesic integrator. |
| **einsteinpy-geodesics** | Optimized geodesic integration for Schwarzschild/Kerr. | MIT | **Evaluate.** May be fast enough for runtime use. If so, integrate directly. |
| Custom | Geodesic integrator with thrust. | -- | **Likely build.** Hard-code the Kerr metric Christoffel symbols (they're known analytically). Integrate geodesic equations + four-acceleration using our existing DOP853 integrator. This is the cleanest path to coupling GR with ship thrust. ~300 lines of code. The math is well-documented in Misner-Thorne-Wheeler and Chandrasekhar. |

**Runtime cost:** ~4x the cost of Newtonian integration (more state variables, more expensive force evaluation). Still sub-millisecond.

---

### 2.7 Communication Physics

**What:** Light travel time between ship and any other point (Earth, relay satellites, other ships). Doppler shift of radio signals. Signal strength vs. distance. Communication blackout during atmospheric reentry (plasma sheath). Gravitational frequency shift near compact objects.

**Tier:** On-demand. Computed when the player requests communication or checks signal status.

**Tool evaluation:**

| Tool | What it does | Verdict |
|---|---|---|
| Custom | Light delay + Doppler calculator. | **Build.** Light travel time is distance/c (with relativistic corrections if the ship is moving fast). Doppler shift is the relativistic Doppler formula from Section 2.4. Plasma blackout is a step function based on atmospheric density and velocity. ~50 lines of code. |

**Runtime cost:** One evaluation per query. Negligible.

---

## 3. THE GALACTIC MEAN FIELD

The mean field is the total gravitational potential of the Milky Way. It IS the galaxy -- all mass, all structure, all density features. It evolves under its own gravity (self-consistent). Stars orbit within it. The ship navigates through it.

**Simulation timescale:** 10 million years (sufficient to cross the galaxy at 0.01c). Over this window, the bar rotates ~22°, spirals ~13° (analytic rigid patterns; Section 3.1). Analytical stellar orbits are accurate to >99.7% everywhere, including the galactic center (0.3% velocity error from unmodeled encounters at the galactic center -- acknowledged limitation). Star-star encounters are out of scope.

---

### 3.1 Mean Field Architecture

**Hybrid: analytic structure + basis-function perturbation.** (Revised 2026-07-04, audit item
1.5. The previous pure-SCF plan had three defects: a spherical-harmonic basis needs tens of
thousands of coefficients — not 500 — to represent a 300 pc thin disk; free SCF evolution
cannot be made to end at the observed MW at T=0 (bar angle, arm loci — that is made-to-measure
research); and 10^6 tracer particles give percent-level force noise, above the accuracy
targets. The hybrid matches observations BY CONSTRUCTION and the runtime force rule is
unchanged: one Φ, one gradient.)

Φ(r, t) = Φ_axi(r) + Φ_bar(r, t) + Φ_spiral(r, t) + Φ_warp(r)
        + Σ Φ_cluster_i(|r − r_i(t)|) + Φ_NSC(r) + Φ_NSD(r) + Φ_CMZ(r) + Φ_pert(r, t)

| Term | Model | Time dependence |
|---|---|---|
| Φ_axi | McMillan 2017 best-fit potential-density pairs (thin/thick disks, HI/H2 disks, bulge, NFW halo — Section 3.2.1) | Static |
| Φ_bar | Analytic rotating bar (Ferrers-type density-potential pair), Section 3.2.2 parameters; bar mass deducted from the axisymmetric bulge to avoid double-counting | Rigid rotation, Ω_bar = 38 km/s/kpc, 27° at T=0 |
| Φ_spiral | Cox & Gomez 2002 potential-density pair, Section 3.2.3 parameters — the SAME pattern that modulates star counts | Rigid rotation, Ω_spiral = 23 km/s/kpc |
| Φ_warp | Warp/flare enter the disk DENSITY model (generation, ISM); leading-order potential correction only | Static |
| Clusters, nuclear | Plummer/King profiles on precomputed orbits; NSC/NSD/CMZ analytic at the center | Cluster orbits r_i(t) |
| Φ_pert | LMC wake + halo/disk response: offline CONSTRAINED N-body (analytic MW + LMC on its measured orbit; EXP or AGAMA), residuals fitted with a low-order basis expansion (halo-appropriate, l ≤ 4, ~150 coefficients) over the 10 Myr window | Coefficient time series |

Self-consistency by construction: every structural term is a potential-DENSITY pair, and the
density side of the same pairs is what drives star generation (Section 3.4) and the gas model.
The density that sources Φ is the density that makes stars — no SCF/analytic mismatch.

Acknowledged limitation: the bar and spiral are rigid patterns; their self-consistent
evolution over 10 Myr is neglected. The neglected effect is far smaller than the parameter
uncertainties (Ω_bar ± 3 km/s/kpc alone dominates it over this window).

**Precomputation cost:**

| Item | Offline time | Storage |
|---|---|---|
| LMC/halo response: constrained N-body + basis fit (EXP/AGAMA) | ~hours-1 day | ~12 MB (~150 coefficients × 10^4 snapshots over 10 Myr) |
| Cluster orbits (~220 clusters as test particles in the hybrid Φ) | ~minutes | ~5 MB (Chebyshev trajectory coefficients) |
| Epicyclic frequency tables Ω(R), κ(R), ν(R) from Φ_axi | ~1 minute (from galpy) | ~1 MB |
| **Total** | **~1 day** | **~18 MB** |

**Runtime cost:** analytic terms ~0.5-1 μs + perturbation expansion ~0.5-1 μs + nearby
clusters/nuclear ~0.5 μs = **~1-3 μs per potential evaluation** (cheaper than the old SCF
budget, and noise-free).

---

### 3.2 Milky Way Structure (Observed Parameters)

All values from McMillan 2017, Bland-Hawthorn & Gerhard 2016, GRAVITY Collaboration 2022, Reid+ 2019, and Gaia DR3 results.

**Fundamental constants:**

| Parameter | Value |
|---|---|
| Distance to galactic center R_0 | 8.275 ± 0.034 kpc |
| Circular velocity at Sun V_0 | 233 ± 3 km/s |
| Solar motion (U, V, W) | (11.1, 12.2, 7.25) km/s |
| Sun height above midplane z_0 | 20.8 ± 0.3 pc |
| Oort A | 15.3 ± 0.4 km/s/kpc |
| Oort B | -11.9 ± 0.4 km/s/kpc |
| Local escape velocity | 528 +16/-24 km/s |

#### 3.2.1 Axisymmetric Components (McMillan 2017 best fit)

| Component | Profile | Mass (M_sun) | Scale length | Scale height |
|---|---|---|---|---|
| Thin stellar disk | Double exponential | 3.45 × 10^10 | R_d = 2.50 kpc | h_z = 0.30 kpc |
| Thick stellar disk | Double exponential | 0.53 × 10^10 | R_d = 3.02 kpc | h_z = 0.90 kpc |
| HI gas disk | Double exponential + flare | 0.87 × 10^10 | R_d = 7.0 kpc | h_z = 0.085 kpc |
| H2 gas disk | Double exponential | 0.21 × 10^10 | R_d = 1.5 kpc | h_z = 0.045 kpc |
| Bulge | Power-law + exp cutoff: ρ ~ r^-1.8 exp(-(r/r_cut)^2) | 0.91 × 10^10 | r_cut = 2.1 kpc | 3D |
| Dark matter halo | NFW: ρ = ρ_0 / [(r/r_s)(1 + r/r_s)^2] | M_200 = 1.26 × 10^12 | r_s = 19.6 kpc | Spherical (q ~ 0.9) |

Disk density profile: ρ(R, z) = (Σ_0 / 2h_z) × exp(-R/R_d) × exp(-|z|/h_z)

Mass enclosed at key radii: 8 kpc → ~9 × 10^10 M_sun (50% DM). 50 kpc → ~4.5 × 10^11 M_sun (90% DM). 200 kpc → ~1.2 × 10^12 M_sun (97% DM).

#### 3.2.2 The Bar

The bulge is not a classical bulge but a boxy/peanut-shaped (B/P) pseudo-bulge with X-shaped morphology, formed by bar buckling instability.

| Parameter | Value |
|---|---|
| Bar stellar mass | ~1.2 × 10^10 M_sun (embedded in bulge mass) |
| Bar half-length | 5.0 ± 0.2 kpc |
| Axis ratios (x:y:z) | 1 : 0.4 : 0.3 |
| Bar angle (major axis to Sun-GC line) | 27 ± 2° |
| Pattern speed Ω_bar | 38 ± 3 km/s/kpc |
| Corotation radius | ~6 kpc |
| Rotation period | ~165 Myr |

The bar is an analytic rotating density-potential pair (Ferrers-type) pinned to these observed parameters, rotating rigidly at Ω_bar with 27° phase at T=0. Its mass is deducted from the axisymmetric bulge term to avoid double-counting (Section 3.1).

#### 3.2.3 Spiral Arms

Four-arm logarithmic spirals: R(φ) = R_ref × exp[(φ - φ_ref) × tan(ψ)]

| Arm | Pitch angle ψ | R_ref (kpc) | Notes |
|---|---|---|---|
| Scutum-Centaurus | 14.0 ± 1.5° | 5.0 | Major arm, dominates inner galaxy |
| Sagittarius-Carina | 15.5 ± 2.0° | 6.6 | Prominent in gas |
| Perseus | 9.4 ± 1.4° | 9.9 | Major arm, Sun lies interior |
| Norma/Outer | 13.8 ± 3.0° | 4.5 | Connects to outer arm |
| Local Arm (Orion Spur) | 11.4 ± 1.4° | 8.3 | Not a major arm; Sun is inside |

Pattern speed: Ω_spiral ~ 23 km/s/kpc (rotation period ~275 Myr). Arm-interarm density contrast: 2-3:1 in old stars, up to 10:1 in gas/young stars. Arm width FWHM: 0.5-1.0 kpc.

Spiral structure is the analytic Cox & Gomez (2002) potential-density pair rotating rigidly at Ω_spiral. The SAME pattern sources both the gravitational potential and the star-count modulation (self-consistent by construction; see generation.md Section 3.2 for how arm membership is encoded deterministically in star generation).

#### 3.2.4 Sgr A* and Nuclear Structures

| Structure | Mass | Scale | Profile |
|---|---|---|---|
| Sgr A* | 4.297 ± 0.013 × 10^6 M_sun | Point source (R_s = 0.08 AU) | Point mass, spin a < 0.1 |
| Nuclear Star Cluster | 2.5 ± 0.4 × 10^7 M_sun | r_eff = 4.2 pc, q = 0.7 | Plummer-like |
| Nuclear Stellar Disk | 1.0 ± 0.3 × 10^9 M_sun | R = 230 pc, h_z = 45 pc | Exponential |
| Central Molecular Zone | 3-5 × 10^7 M_sun (gas) | R < 200 pc | Twisted elliptical ring |

Sgr A* is always an active point source on the ship (and on all stars within its sphere of influence, ~2 pc). The NSC, NSD, and CMZ are analytical profiles centered at the galactic center, added to the hybrid mean field (Section 3.1).

#### 3.2.5 Disk Warp and Flare

Beyond R ~ 10 kpc, the disk warps out of the plane and the scale height increases (flares).

Warp: z_warp(R, φ) = A_w × (R - R_w)^1.7 × sin(φ - φ_LON)
- Onset R_w ~ 10 kpc
- Displacement at R = 15 kpc: 0.5-1.0 kpc. At R = 20 kpc: 3-5 kpc.
- Asymmetric: northern warp rises higher.

Flare: h_z(R) = h_z,0 × exp[(R - R_0) / R_flare], R_flare ~ 8-10 kpc. Scale height doubles by R = 12 kpc, reaches ~3-4 kpc at R = 20 kpc.

Both enter the disk DENSITY model (star generation, ISM) analytically; the warp's potential contribution is a leading-order correction term (it is second-order for dynamics over this window).

#### 3.2.6 Cluster System

~170 globular clusters (total mass ~4 × 10^7 M_sun, individual range 10^3 - 3.5 × 10^6 M_sun) plus ~50 significant open clusters (10^4 - 10^5 M_sun). Each is a localized density feature in the mean field, modeled as a Plummer or King profile orbiting within the hybrid potential.

Cluster orbits are precomputed offline as test particles in the hybrid Φ (Section 3.1) and stored as Chebyshev trajectory coefficients. At runtime, cluster positions at time T are evaluated by polynomial interpolation. The cluster's gravitational profile is added to the mean field at that position.

Stars generated near a cluster are bound to its local potential well (see Section 3.4).

Notable clusters: Omega Centauri (3.5 × 10^6 M_sun, likely stripped dwarf nucleus), 47 Tucanae (7 × 10^5 M_sun), M54 (1.5 × 10^6 M_sun, nucleus of Sgr dwarf).

#### 3.2.7 LMC and Satellite Galaxies

| Satellite | Total mass | Distance | Effect |
|---|---|---|---|
| LMC | 1.4 × 10^11 M_sun | 49.5 kpc | MW disk accelerates ~30 km/s toward LMC; DM halo wake |
| SMC | ~3 × 10^9 M_sun | 62.1 kpc | Minor |
| Sgr dSph | ~4 × 10^8 M_sun (bound remnant) | 18 kpc | Generates Sgr stream; perturbs outer disk |

The LMC is massive enough (~10-15% of M_MW) to significantly perturb the outer galaxy beyond ~30 kpc. It is included via the Φ_pert term (Section 3.1): an offline constrained N-body run (analytic MW + LMC on its measured orbit) is fitted with a low-order basis-function expansion of the halo response — wake + reflex motion — and the coefficient time series over the 10 Myr window is stored (~12 MB). The LMC itself is also a direct analytic profile term at its (precomputed) orbital position.

~60 known satellite galaxies total; most are ultra-faint dwarfs (M < 10^5 M_sun) contributing negligibly to the potential.

#### 3.2.8 Velocity Dispersions by Population

| Population | σ_R (km/s) | σ_φ (km/s) | σ_z (km/s) | Mean rotation lag |
|---|---|---|---|---|
| Thin disk (young, < 1 Gyr) | 25-30 | 18-22 | 12-15 | 0-5 km/s |
| Thin disk (old, 5-8 Gyr) | 35-45 | 25-30 | 20-25 | 5-15 km/s |
| Thick disk | 60-70 | 45-50 | 35-45 | 30-50 km/s |
| Stellar halo | 140-160 | 90-100 | 80-100 | ~200 km/s (near zero net rotation) |
| Bulge | 110-130 | 100-110 | 80-100 | Complex (bar kinematics) |

These velocity dispersions determine the epicyclic amplitudes of generated stars and are used to assign velocities during procedural generation.

#### 3.2.9 Galactic Object Census (generation targets)

The star count is DERIVED from the mass model, not chosen: the same density that sources Φ
makes the stars, so N = (M_stellar − M_remnants) / ⟨m⟩_Kroupa. With McMillan 2017 stellar
components (4.89×10¹⁰ M_☉ disk+bulge + ~1.4×10⁹ halo), remnant mass ~8×10⁹ M_☉, and
present-day Kroupa mean ⟨m⟩ ≈ 0.35 M_☉:

**N_living ≈ 1.2×10¹¹, total stellar objects ≈ 1.3×10¹¹ (~130 billion).**

(The previously used "200 billion" sits inside the popular 100-400 billion literature range
but is inconsistent with this project's own mass model by ~50% — it would require
M_stellar ≈ 7.5×10¹⁰ M_☉, above all modern estimates. Licquia & Newman 2015 meta-analysis:
6.08±1.14×10¹⁰ M_☉ total stellar; Gaia-era work trends lower, not higher.)

| Population | Count | Basis |
|---|---|---|
| Living stars (≥ 0.08 M_☉) | ~1.2×10¹¹ | Derived: mass model ÷ Kroupa PDMF mean |
| White dwarfs | ~1.0×10¹⁰ | Napiwotzki 2009; standard census |
| Neutron stars | ~10⁹ | Birth rate × Galactic age |
| Stellar-mass black holes | ~10⁸ | Birth rate × Galactic age |
| Brown dwarfs (13-80 M_jup) | 2.5×10¹⁰-10¹¹ | RCW 38 extrapolation (Muzic+ 2017) |
| Rogue planets (unbound) | ~20 per star ≈ 2-3×10¹² | MOA 9-yr microlensing (Sumi+ 2023); mass function rises steeply toward sub-Earth masses |
| Bound planets | ~1-2 per star ≈ 1.5-3×10¹¹ | Kepler/TESS occurrence (Section 3.5) |
| Globular clusters | ~170 | Harris catalog |
| Open clusters | ~10³ known; ~10⁵ estimated | Completeness-limited |
| Satellite galaxies | ~60 (most ultra-faint) | Section 3.2.7 |
| Planetary nebulae (active) | ~10⁴ estimated (~3×10³ known) | PN lifetime × death rate |
| Supernova remnants (visible) | ~10³ expected (~300 known) | SN rate × SNR lifetime |

All generation-layer normalizations (Section 3.4, generation.md) target these values;
population checks validate against them.

---

### 3.3 Stellar Motion

Every star's position at coordinate time T is a deterministic, analytical function of its orbital parameters and T. No integration history, no persistent state. Stars are generated on demand and discarded.

**Disk stars (R > 2 pc from Sgr A*):** Epicyclic orbits in the axisymmetric mean field.

- R(t) = R_guide + A_R × cos(κ × t + φ_R)
- φ(t) = Ω × t + (2Ω / κ) × A_R × sin(κ × t + φ_R) / R_guide
- z(t) = A_z × cos(ν × t + φ_z)

Where Ω(R), κ(R), ν(R) are precomputed from the potential and stored as interpolation tables. Each star's orbital parameters (R_guide, A_R, A_z, phases) are derived from its seed.

Computational cost: ~30 FLOPs per star (~10 ns).

**Galactic center stars (within Sgr A*'s sphere of influence, ~2 pc):** Keplerian orbits around Sgr A* with 1PN precession.

- Solve Kepler's equation at time T (Newton-Raphson, ~5 iterations)
- Apply precession: ω(T) = ω_0 + dω/dt × T
- Precession rate: Schwarzschild precession + Newtonian precession from extended NSC mass
- 3D coordinate rotation from orbital plane to galactic frame

Computational cost: ~150 FLOPs per star (~50 ns).

**Stars bound to clusters:** Position = cluster CM position at T (from precomputed Chebyshev trajectory) + star's orbit within the cluster's Plummer/King potential at T (epicyclic-like evaluation in the cluster potential).

Computational cost: ~200 FLOPs per star (~70 ns).

**Binary/multiple system components:** Position = system galactic position at T + component offset from the exact two-body solution at T (Kepler evaluation of the mutual orbit).

Computational cost: additional ~150 FLOPs per component.

**Halo stars (non-circular orbits):** The epicyclic approximation breaks down for highly eccentric halo orbits. These use a higher-order orbital representation precomputed with galpy: action-angle variables or Chebyshev-interpolated orbit segments. More expensive but halo stars are sparse (~10^-4 of the local density in the disk).

**Maximum system size by environment:**

| Location | R | Stellar density (pc^-3) | σ (km/s) | Max stable system (AU) |
|---|---|---|---|---|
| Near Sgr A* | 0.001 pc | 10^7 | 500 | ~0.1 |
| Inner NSC | 0.01 pc | 10^6 | 300 | ~0.5 |
| NSC edge | 0.1 pc | 10^5 | 150 | ~3 |
| 1 pc from Sgr A* | 1 pc | 10^3 | 130 | ~30 |
| Inner bulge | 100 pc | 100 | 110 | ~100 |
| Inner disk | 3 kpc | 10 | 80 | ~500 |
| Solar neighborhood | 8.2 kpc | 0.1 | 30 | ~20,000 |
| Outer disk | 12 kpc | 0.02 | 25 | ~50,000 |
| Halo | 50 kpc | 10^-4 | 150 | ~100,000 |

Systems wider than the local maximum were destroyed by encounters over the galaxy's lifetime and do not exist. The generator clips binary separations, planetary system extents, and Oort cloud radii to these limits. Binary fraction also decreases toward denser environments.

**Determinism:** The spatial index is keyed by time-invariant quantities (guiding center for disk stars, semi-major axis for galactic center stars). Queries expand the search radius by the maximum orbital excursion (epicyclic amplitude or apocenter-pericenter range). Position at T is a pure function of (orbital_params, T). Same (ship_position, T) always produces the same set of nearby stars.

**Tool evaluation:**

| Tool | What it does | License | Verdict |
|---|---|---|---|
| **galpy** | Galactic potentials, epicyclic frequencies, orbit integration, action-angle transforms. | BSD-3 | **Integrate.** Use offline to compute Ω(R), κ(R), ν(R) tables and to precompute halo star orbit representations. At runtime, use custom evaluators for speed. |
| **astropy** | Coordinate transforms, time systems. | BSD-3 | **Integrate.** Galactic ↔ ICRS ↔ ecliptic transforms throughout. |

---

### 3.4 Procedural Star Generation

~130 billion stars (derived from the mass model — Section 3.2.9). Generated on demand from deterministic seeds, evaluated at coordinate time T, discarded when no longer needed. No persistent state.

#### 3.4.1 Generation Pipeline

Given a query (ship_position, coordinate_time T, search_radius):

1. **Determine which spatial cells overlap the search volume.** Cells are defined in time-invariant orbital parameter space: (R_guide, φ_guide, z_guide) for disk stars, (a, e) bins for galactic center stars.

2. **For each cell, generate stars from the cell's deterministic seed.** The seed is derived from the cell's coordinates (hash of integer cell indices). The number of stars per cell comes from the local stellar density model (Section 3.2 axisymmetric profiles + spiral arm modulation).

3. **For each generated star, evaluate position at time T.** Epicyclic formula for disk stars, Kepler + precession for galactic center stars, cluster orbit + internal orbit for cluster members.

4. **Filter to stars actually within the search volume.** Discard stars whose evaluated position at T falls outside the query region.

5. **Check Gaia overlay.** If a real Gaia star occupies this region, suppress the corresponding number of procedural stars to avoid double-counting (Section 3.7).

6. **Return star list** with positions, velocities, masses, spectral types, and any companions.

#### 3.4.2 Stellar Density Model

The number of stars per unit volume at galactic position (R, φ, z) is:

n(R, φ, z) = n_thin(R, z) + n_thick(R, z) + n_halo(r) + n_bulge(r) + n_spiral(R, φ, z)

Each component uses the density profiles from Section 3.2. The spiral arm contribution modulates the disk density:

n_spiral = n_disk × (f_arm - 1) × S(R, φ)

where f_arm ~ 2-3 is the arm-interarm contrast and S(R, φ) is a smooth function peaked along the logarithmic spiral loci.

Total star count: ~1.3×10¹¹ — an OUTPUT of the density model + IMF, not a target (Section 3.2.9). Consistent with the 100-400 billion literature range and, by construction, with the McMillan mass model that sources the potential.

#### 3.4.3 Stellar Properties

For each generated star, draw properties from the seed:

**Mass:** Kroupa/Chabrier initial mass function (IMF):
- dN/dm ∝ m^-1.3 for 0.08 < m/M_sun < 0.5
- dN/dm ∝ m^-2.3 for 0.5 < m/M_sun

Most stars are M-dwarfs (0.08-0.5 M_sun). Massive stars (> 8 M_sun) are rare and found preferentially in spiral arms and young regions.

**Spectral type and luminosity:** From mass via main-sequence relations. L ∝ M^3.5 (approximate). Effective temperature from mass-T_eff relations. Spectral type (OBAFGKM) from T_eff.

**Age:** Drawn from the star formation history of the host population (thin disk: roughly uniform 0-10 Gyr; thick disk: 8-12 Gyr; halo: 10-13 Gyr; bulge: 8-13 Gyr).

**Metallicity:** From the age-metallicity relation and galactic gradient. [Fe/H] decreases with galactocentric radius (~-0.06 dex/kpc) and increases with age.

**Evolved stars:** Stars older than their main-sequence lifetime (t_MS ≈ 10 Gyr × (M/M_sun)^-2.5) are evolved. Depending on initial mass: red giant, white dwarf, neutron star, or black hole (see Section 3.6).

**Velocity:** Circular velocity from the rotation curve + random dispersion drawn from the population's velocity ellipsoid (Section 3.2.8). Converted to epicyclic orbital parameters.

#### 3.4.4 Binary and Multiple System Generation

When a star is generated, roll for multiplicity from the seed:

1. Draw binary fraction from the environment-dependent rate:
   - Solar neighborhood: ~46% for solar-type (Raghavan+ 2010), ~70% for massive stars (Sana+ 2012), ~26% for M-dwarfs (Winters+ 2019) — canonical sourced table: generation.md Section 3.5
   - Scales down toward galactic center per the encounter survival limit (Section 3.3)

2. If binary: draw companion mass ratio from observed distribution (Raghavan+ 2010: roughly uniform for solar-type, favoring q → 1 for massive stars). Draw orbital period from the log-normal period distribution (peak at ~10^5 days for solar-type).

3. Clip separation to the local maximum system size. If the drawn separation exceeds the limit, reject and generate a single star instead.

4. For triples (~10% of multiples): generate as hierarchical pair + distant companion. Inner pair has short period, outer companion has long period. Both clipped to local limits.

The system is one atomic generation unit from one seed. Both components get consistent orbital parameters. The ship sees individual point sources -- it doesn't know they came from the same seed.

---

### 3.5 Planetary System Generation

Every star has a probability of hosting planets, drawn from Kepler/TESS occurrence rates. (The table below is an approximate summary; generation.md Section 4.2 carries the canonical sourced values used by the generator.)

#### 3.5.1 Occurrence Rates

| Planet type | Radius | Period range | Occurrence per star |
|---|---|---|---|
| Hot Jupiters | > 6 R_earth | < 10 days | ~1% |
| Warm Jupiters | > 6 R_earth | 10-200 days | ~3% |
| Cold Jupiters | > 6 R_earth | 200 days - 10 yr | ~10% |
| Sub-Neptunes | 1.7-3.5 R_earth | < 400 days | ~30% |
| Super-Earths | 1.0-1.7 R_earth | < 400 days | ~30% |
| Earth-sized | 0.7-1.0 R_earth | < 400 days | ~15% |
| Long-period giants (RV) | > 0.3 M_jup | 1-20 yr | ~15-20% |

Overall: ~50% of Sun-like stars have at least one planet. For M-dwarfs, occurrence rates are higher for small planets, lower for giants.

#### 3.5.2 System Architecture

Generate systems following observed patterns:
- **Multiplicity:** Kepler systems show 1-7 transiting planets. Generate from the observed multiplicity distribution.
- **Spacing:** Adjacent planets are spaced by ~10-30 mutual Hill radii (Weiss+ 2018). Peas-in-a-pod pattern: adjacent planets tend to have similar sizes and regular spacing.
- **Inclination:** Small mutual inclinations (1-3°) for compact multi-planet systems.
- **Eccentricity:** Low for compact systems (e < 0.05). Higher for systems with giant planets (e ~ 0.1-0.3).

#### 3.5.3 Planet Properties

From seed: mass (from radius via mass-radius relation), bulk density, orbital elements, number and properties of moons (probabilistic, based on planet mass), ring systems (giant planets only, ~50% occurrence), atmosphere composition (from equilibrium temperature and mass).

System extent clipped to local environment maximum (Section 3.3). Oort clouds and wide debris disks truncated at the encounter survival limit for the star's galactic position.

---

### 3.6 Extreme Object Generation

Stellar remnants and compact objects are generated as the endpoints of stellar evolution.

#### 3.6.1 White Dwarfs

Stars with initial mass 0.5-8 M_sun end as white dwarfs after exhausting main-sequence lifetime. WD mass from the initial-final mass relation (Cummings+ 2018): M_WD ≈ 0.08 × M_initial + 0.48 M_sun. Surface temperature from cooling age: T_eff ∝ t^(-0.35) (Mestel's t^(-1.4) is the luminosity exponent; T follows from L ∝ T⁴). Compose ~5-10% of all stellar objects in the galaxy. ~10 billion total.

#### 3.6.2 Neutron Stars

Stars with initial mass 8-25 M_sun produce neutron stars in core-collapse supernovae. NS mass: 1.1-2.2 M_sun (peaked at ~1.4 M_sun). Radius: ~12 km. Birth rate: ~2 per century in the MW. Total population: ~10^9. Most are old, cold, and invisible. ~2000 known as radio pulsars.

Properties from seed: mass, spin period (birth: ~10 ms, slows with age), magnetic field (birth: 10^12-10^13 G, decays over ~10^7 yr). Magnetars: ~10% of NS born with B > 10^14 G, spin periods 2-12 s.

For the flight simulator: NS gravity is Newtonian at distance (point source), transitions to GR (Schwarzschild or slowly rotating Kerr) within ~100 R_s ≈ 410 km for a 1.4 M_sun NS (R_s = 4.14 km; the NS's ~12 km material radius is NOT its Schwarzschild radius). Surface gravity: ~1.3 × 10^12 m/s^2 Newtonian (~1.6 × 10^12 with the GR correction factor). Tidal forces become significant within ~1000 km.

#### 3.6.3 Stellar-Mass Black Holes

Stars with initial mass > 25 M_sun (approximately; depends on metallicity) collapse to black holes. BH mass: 3-100 M_sun. Typical: 7-15 M_sun. Birth rate: ~0.2 per century in the MW. Total population: ~10^8. Nearly all are isolated, dark, and undetectable except by gravitational lensing.

Properties from seed: mass (from initial mass via fallback prescriptions), spin (a ~ 0-0.9, uncertain distribution). Kerr metric for ship interactions.

BH in binaries: X-ray binaries (~300 known in MW) have a BH accreting from a companion star. Generate these as binary systems where the primary has evolved to a BH.

#### 3.6.4 Intermediate-Mass Black Holes

Hypothesized in the centers of some globular clusters (10^2-10^5 M_sun). Evidence is debated. For the simulator: optionally place one in the most massive globular clusters (Omega Centauri is the best candidate, possibly ~4 × 10^4 M_sun). These provide additional GR destinations within the galaxy.

#### 3.6.5 Generation Logic

When a star's age exceeds its main-sequence lifetime, it is generated as a remnant:
- M_initial < 0.5 M_sun: still on main sequence (t_MS > age of universe)
- 0.5-8 M_sun: white dwarf
- 8-25 M_sun: neutron star (with kick velocity 200-500 km/s added to galactic orbit)
- > 25 M_sun: black hole (smaller kick, ~50-100 km/s)

The original star's seed determines the remnant type, mass, and properties. Neutron star kicks displace the remnant from its birth location, reflected in modified epicyclic parameters.

---

### 3.7 Gaia DR3 Integration

~1 million real stars from the Gaia DR3 catalog provide ground truth throughout the galaxy.

#### 3.7.1 Sample Selection

Pre-download from Gaia archive (via astroquery). NOTE: "all stars within 500 pc" is NOT
feasible — at the local density 0.1 pc⁻³ that sphere holds ~52 MILLION stars (the original
"~500K" estimate was off by ~100×). The 1M-star budget is allocated instead as:

- **Complete within 100 pc** (~330K stars, GCNS-quality completeness — every star the player
  can quickly reach is real)
- Bright navigational landmarks, all-sky (V < 6, naked-eye, ~9,000)
- Stars with known exoplanets (NASA Exoplanet Archive cross-match, ~5,500 hosts)
- All O/B stars and cluster members with 6D phase space (~100K — they anchor arms and clusters)
- Representative 6D-complete sampling of the 100-500 pc shell and distant structure
  (arms, bulge, halo), ~550K, sampling fraction recorded per region so procedural
  suppression (Section 3.7.3) subtracts the CORRECT number of procedural stars per cell
  (suppress N_gaia / sampling_fraction, not N_gaia)
- Known binary systems with orbital solutions (subset of the Gaia DR3 non-single-star catalog)

Beyond 100 pc the sky is therefore a mix: real anchors + procedural fill calibrated to the
same density model. Total ~1M stars, ~500 MB with orbital parameters.

#### 3.7.2 Orbital Parameter Derivation

For each Gaia star with full 6D phase space (position + proper motion + radial velocity + parallax):
1. Convert to galactocentric coordinates using astropy
2. Compute galactic orbital parameters using galpy: guiding center radius, epicyclic amplitudes, phases at epoch
3. Store: orbital parameters + stellar properties (G magnitude, BP-RP color, T_eff, [Fe/H] where available)

For stars without radial velocity: assign from the local velocity distribution (population-dependent, using Section 3.2.8 dispersions).

#### 3.7.3 Spatial Index and Suppression

Gaia stars are stored in a separate spatial index (k-d tree over guiding center positions). When the procedural generator queries a region:

1. Query the Gaia index for real stars in the region
2. For each procedural generation cell overlapping a Gaia star, reduce the procedural star count by the number of Gaia stars present
3. Return Gaia stars alongside (reduced) procedural stars

This prevents double-counting while ensuring real stars take precedence. In the local neighborhood (< 500 pc), most visible stars are real Gaia entries. Beyond 500 pc, most are procedural with scattered Gaia landmarks.

#### 3.7.4 Known Exoplanet Systems

Stars with confirmed exoplanets (NASA Exoplanet Archive, ~5500 systems as of 2025) get their real planetary systems instead of procedurally generated ones. Orbital elements, masses, and radii from the catalog. The ship can visit real systems (TRAPPIST-1, Proxima Centauri b, etc.) with accurate planetary configurations.

---

### 3.8 Minor Bodies

**Solar system:** Orbital elements for ~10,000 known asteroids, comets, and KBOs from JPL Small-Body Database (pre-downloaded via astroquery). Positions at time T from Keplerian propagation. Kuiper belt and Oort cloud extent procedurally generated, clipped to the encounter survival limit for the solar neighborhood (~20,000 AU).

**Other star systems:** Minor body populations procedurally generated when the ship enters a system. Asteroid belt probability and extent from the system's architecture. Oort cloud extent from the local environment maximum.

**Tool evaluation:**

| Tool | What it does | License | Verdict |
|---|---|---|---|
| **SpiceyPy** | Solar system body ephemerides. | Public domain | **Integrate.** Primary source for planet/moon positions. |
| **astroquery** | Download Gaia DR3, exoplanet catalogs, small body database. | BSD-3 | **Use offline for data acquisition.** |
| **astropy** | Coordinates, units, time. | BSD-3 | **Foundation library throughout.** |
| **galpy** | Galactic potential, epicyclic frequencies, orbit integration. | BSD-3 | **Use offline for orbital parameter computation.** |

---

## 4. CROSS-CUTTING PHYSICS (RESEARCH AREAS)

These components span multiple sections and involve open research questions identified in SCOPE.md.

---

### 4.1 Relativistic ISM Interaction

**What:** At interstellar velocities, the ship collides with interstellar medium (hydrogen, helium, dust). Below ~0.3c, energy deposition follows Bethe-Bloch stopping. Above ~0.6c, hadronic cascades (pion production) produce penetrating secondary radiation. The transition between these regimes is uncharacterized.

**Tier:** Per-tick conditional (interstellar regime only).

**Coupling:** Feeds into thermal (1.3), radiation (1.7), structural degradation (1.5), and resource management (1.4 -- shield ablation consumes mass).

**Tool evaluation:**

| Tool | What it does | Verdict |
|---|---|---|
| **Geant4** | Full hadronic cascade simulation. | **Use offline (research).** Run Geant4 to compute cascade energy deposition profiles through realistic shield geometries at 0.3-0.9c. This is the core of original research paper #1 from SCOPE.md. Results become lookup tables for runtime. |
| Custom | Runtime ISM interaction model. | **Build.** Table lookup indexed by velocity and shield state. Below pion threshold: Bethe-Bloch formula (analytic). Above: precomputed cascade tables. |

---

### 4.2 Ship Characterization in Curved Spacetime

**What:** Separating ship-intrinsic physics (propulsion, thermal, structural) from environment-dependent physics (tidal stress, gravitational time dilation) near compact objects. The lookup table approach assumes forces depend on ship state, but near a black hole, forces also depend on position in the gravitational field.

**Tier:** Architectural. This affects how the entire ship model interfaces with the GR module.

**Approach (proposed):** Factor the physics into:
- **Intrinsic quantities:** Thrust, propellant flow, thermal state, structural stiffness -- these depend on ship state and are lookup table material.
- **Extrinsic quantities:** Tidal acceleration gradient across the ship, gravitational time dilation, frame dragging torque -- these depend on position in the metric and must be computed in real time by the GR module.
- **Coupled quantities:** Structural stress from tidal forces (extrinsic force applied to intrinsic stiffness matrix). Computed at runtime by combining real-time GR tidal tensor with precomputed ship structural model.

This is an architectural research area (candidate paper: "ship-systems characterization in curved spacetime"). Note: it is NOT on SCOPE.md's numbered research list — SCOPE's #2 is the laser sail; the numbering here previously contradicted it.

---

### 4.3 Regime Detection and Switching

**What:** Automatically determining which physics modules should be active based on the ship's current state and environment. This is not a physics problem per se, but a computational architecture problem that determines which of the above components run each tick.

**Tier:** Per-tick core. Evaluated every tick, but the evaluation itself is cheap.

**Approach:** Define activation thresholds for each conditional module:

| Module | Activation condition |
|---|---|
| Atmospheric drag/heating | Altitude < atmosphere_top (body-dependent) |
| Special relativity | v/c > 0.01 |
| Post-Newtonian | GM/(rc^2) > 10^-10 for any nearby body |
| General relativity | Within 100 R_s of a compact object |
| ISM interaction | Interstellar regime + v > 100 km/s |
| Life support | Crew present |
| Radiation belts | Within magnetosphere of relevant body |

Transitions must be smooth (no discontinuous jumps in forces when a module activates). Use blending functions or ensure that conditional corrections are negligible at their activation threshold.

**Tool evaluation:** Build from scratch. This is custom logic specific to the simulator's architecture.

---

## 5. EXTERNAL TOOLS SUMMARY

### Runtime dependencies (Rust — per SCOPE.md, exactly two C FFI libraries + pure-Rust crates)

| Tool | Use | License | Section |
|---|---|---|---|
| **NAIF SPICE** (C, FFI) | Solar system ephemerides inside the DE441 span | Public domain | 2.1, 3.8 |
| **NRLMSISE-00** (C, FFI) | Earth atmosphere model (deterministic synthetic solar indices) | Public domain | 2.3 |
| **libm** (Rust crate) | Bit-reproducible transcendentals | MIT/Apache | determinism.md |
| **crossterm / ratatui** (Rust crates) | Terminal interface | MIT | interface.md |

Everything else at runtime is custom Rust. No Python, no numpy/scipy/astropy in the simulator
process.

### Offline (Python — characterization, precomputation, data preparation)

| Tool | Use | License | Section |
|---|---|---|---|
| **astropy** | Coordinate transforms, time systems (data prep) | BSD-3 | 3.7 |
| **scipy / numpy** | Offline numerics + runtime-solver validation references | BSD | 1.3, 2.2 |
| **SpiceyPy** | SPICE access from the offline pipeline | MIT | 3.8 |
| **CoolProp** | Fluid property tables (precomputed) | MIT | 1.1, 1.9 |

### Use offline (precomputation only)

| Tool | Use | License | Section |
|---|---|---|---|
| **galpy** | Galactic potential definition, epicyclic frequency tables, orbit integration for halo stars | BSD-3 | 3.1, 3.3 |
| **Cantera** | Chemical kinetics, combustion, life support chemistry | BSD-3 | 1.1 |
| **Geant4** | Particle transport, radiation shielding, ISM cascades | Custom (free) | 1.2, 1.7, 4.1 |
| **OpenMC** | Neutron transport, reactor shielding | MIT | 1.2, 1.7 |
| **FEniCS** | Finite element structural/thermal analysis | LGPL | 1.3, 1.5 |
| **astroquery** | Download Gaia DR3, exoplanet catalogs, JPL SBDB | BSD-3 | 3.7, 3.8 |

### Use for validation

| Tool | What it does | License | Verdict |
|---|---|---|---|
| **EinsteinPy** | Validate GR geodesic integration | MIT | 2.6 |
| **REBOUND** | Validate N-body integration, offline trajectory optimization | GPL-3 | 2.1 |

---

## 6. BUILD-FROM-SCRATCH SUMMARY

| Component | Complexity | Lines (est.) | Section |
|---|---|---|---|
| Ship component graph (structure + resources) | Moderate | 500-800 | 1.4, 1.5 |
| Thermal network solver | Moderate | 300-500 | 1.3 |
| Resource flow manager | Low | 200-400 | 1.4 |
| Atmosphere / life support ODE | Low | 150-300 | 1.6 |
| Radiation dose accumulator | Low | 100-200 | 1.7 |
| Electrical / control state machine | Low | 200-400 | 1.8 |
| Pressure tracker | Low | 100-200 | 1.9 |
| Unified gravity (mean field + point sources) | High | 500-800 | 2.1 |
| Force assembler + regime switcher | High | 500-800 | 2.2, 4.3 |
| Special relativity module | Low | 100-200 | 2.4 |
| Post-Newtonian corrections | Low | 50-100 | 2.5 |
| GR geodesic integrator (Kerr) | High | 300-500 | 2.6 |
| Communication physics | Low | 50-100 | 2.7 |
| Mars / other atmosphere models | Low | 100-200 | 2.3 |
| Hybrid mean field evaluator (analytic terms + perturbation expansion) | High | 500-800 | 3.1 |
| LMC-response precomputation pipeline (constrained run + basis fit) | High | 800-1200 | 3.1 |
| Procedural star generator | High | 800-1200 | 3.4 |
| Stellar property assignment (IMF, evolution, remnants) | Moderate | 400-600 | 3.4, 3.6 |
| Binary/multiple system generator | Moderate | 300-500 | 3.4.4 |
| Planetary system generator | Moderate | 400-600 | 3.5 |
| Epicyclic position evaluator | Low | 100-200 | 3.3 |
| Kepler orbit evaluator (GC + binaries) | Moderate | 200-300 | 3.3 |
| Spatial index (guiding center + orbital element cells) | Moderate | 300-500 | 3.4 |
| Gaia catalog integration + suppression | Moderate | 200-400 | 3.7 |
| Unintended reaction detector | Low | 100-200 | 1.1 |
| **Total custom code** | | **~6500-11000** | |

---

## 7. COMPUTATIONAL BUDGET PER TICK

Estimated wall-clock cost per integration step on a single modern CPU core:

| Component | Cost | Frequency |
|---|---|---|
| Mean field gravity evaluation | ~1-3 us | Every tick |
| Point source gravity (100-1000 sources) | ~1-10 us | Every tick |
| Trajectory integration (DOP853) | ~10 us | Every tick |
| Thermal network (100 nodes) | ~5 us | Every tick |
| Resource flow | ~1 us | Every tick |
| Atmosphere/life support | ~2 us | Every tick (crewed) |
| Pressure/electrical | ~1 us | Every tick |
| Radiation accumulation | ~1 us | Every tick |
| Ship table lookups | ~5 us | Every tick |
| **Subtotal (always-on)** | **~30-40 us** | |
| Atmospheric model (NRLMSISE) | ~10 us | Atmospheric only |
| SR corrections | ~1 us | v > 0.01c |
| PN corrections | ~10 us | Near massive bodies |
| GR geodesics | ~50 us | Near compact objects |
| ISM interaction | ~5 us | Interstellar |
| **Subtotal (conditional)** | **~0-75 us** | |
| **Total per step** | **~30-115 us** | |
| **Max tick rate** | **~9,000-30,000 Hz** | |

This comfortably exceeds the 1000 Hz target for atmospheric flight.

**On-demand queries (not per-tick):**

| Query | Cost | When |
|---|---|---|
| Nearby star refresh (spatial index + position evaluation) | ~0.1-1 ms | Every few seconds of sim-time |
| Procedural generation for new region | ~1-10 ms | On region entry |
| Planetary system generation for destination | ~0.1 ms | On approach to star system |
| Gaia catalog lookup | ~0.01 ms | During star refresh |

---

## 8. PRECOMPUTED DATA BUDGET

| Dataset | Size | Precomputation time | Section |
|---|---|---|---|
| LMC/halo response basis coefficients (10 Myr window) | ~12 MB | ~hours-1 day | 3.1 |
| Cluster orbit trajectories (~220 clusters) | ~5 MB | ~minutes | 3.2.6 |
| Epicyclic frequency tables Ω(R), κ(R), ν(R) | ~1 MB | ~1 minute | 3.3 |
| Gaia DR3 sample (~1M stars with orbital params) | ~500 MB | ~1 hour (download + compute) | 3.7 |
| Known exoplanet catalog (~5500 systems) | ~10 MB | ~5 minutes (download) | 3.7.4 |
| JPL small body database (~10K objects) | ~5 MB | ~5 minutes (download) | 3.8 |
| Ship characterization tables (per ship) | ~100 MB | Hours (offline tool) | SCOPE.md |
| JPL SPICE kernels (DE441) | ~3 GB | Download once | 2.1 |
| **Total (excluding SPICE)** | **~630 MB** | **~1.5 hours + LMC run (~1 day, background)** | |
