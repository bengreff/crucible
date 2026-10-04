# Ship Architecture

How the ship works as an integrated system. The component graph (spatial.md) is the ship's structure. The lookup tables (characterization.md) encode its physics. This document specifies how they connect: the per-tick update loop, cascading failure model, multi-object tracking, and the three initial ship definitions.

---

## 1. SHIP STATE VECTOR

The complete state of a player-controlled ship at any instant:

### Trajectory state (integrated by DOP853)

```rust
struct TrajectoryState {
    position: DoubleDouble3,     // galactocentric meters (double-double)
    momentum: Vec3,              // relativistic 3-momentum, kg⋅m/s (f64)
    orientation: Quat,           // ship body frame → galactic frame (f64, unit quaternion)
    angular_velocity: Vec3,      // rad/s in body frame (f64)
    rest_mass: f64,              // kg (decreases as fuel burns)
    proper_time: DoubleDouble,   // seconds (double-double)
    coordinate_time: DoubleDouble, // seconds from J2000.0 (double-double)
}
```

Derived quantities (computed from state, not stored):
- velocity: v = pc²/E where E = √(m²c⁴ + p²c²)
- gamma: γ = E/(mc²)
- kinetic energy: E - mc²
- angular momentum: L = I × ω
- center of mass: from component graph

### Component state (per component)

```rust
struct ComponentState {
    temperature: f64,            // K, from thermal module
    wall_thickness_current: f64, // m, may differ from initial (ablation)
    damage_dpa: f64,             // cumulative displacement damage
    fill_fraction: f64,          // for tanks (0-1)
    failed: bool,                // structural failure flag
    failure_mode: Option<FailureMode>, // how it failed
}
```

### Subsystem state

```rust
struct SubsystemState {
    // Propulsion
    throttle: f64,               // 0-1
    engine_active: Vec<bool>,    // which engines are firing
    
    // Reactor (if present)
    reactor_power_fraction: f64, // 0-1 (fraction of rated power)
    control_rod_position: f64,   // 0-1
    fuel_burnup: f64,            // fraction of initial fuel consumed
    
    // BH drive (if present)
    bh_mass: f64,                // kg
    feed_rate: f64,              // kg/s
    
    // Resources
    valve_states: Vec<ValveState>,
    battery_charge: f64,         // joules
    
    // Atmosphere (per pressurized volume)
    gas_moles: Vec<GasState>,    // [O2, N2, CO2, H2O, ...]
    cabin_pressure: f64,         // Pa
    
    // Radiation
    crew_dose_accumulated: f64,  // Sievert
    electronics_dose: Vec<f64>,  // TID per electronics component
    cached_dose_rate: DoseRateMap, // from last radiation transport run
    
    // Electrical
    power_generation: f64,       // watts
    power_load: f64,             // watts
}
```

### Moment of inertia (recomputed on mass change)

```rust
struct RotationalProperties {
    moi_tensor: Matrix3,         // kg⋅m², in body frame
    center_of_mass: Vec3,        // m, in body frame (offset from geometric center)
    moi_dirty: bool,             // flag: recompute on next tick
}
```

MOI computed from component graph: each component contributes its shape's analytical MOI + parallel axis term for offset from CoM. Sum over all components. Recompute when any component's mass changes by >0.1% or on structural change.

---

## 2. PER-TICK UPDATE LOOP

One trajectory step. All subsystems evaluate once. Cross-subsystem coupling uses previous-tick values (explicit forward coupling, one-tick lag for circular dependencies).

```
TICK N: state at time t, compute state at time t + dt

Phase 1: FORCES (read current state, compute all forces and torques)
├── 1a. Gravity
│     Mean field gradient at position: ∇Φ(x)                    [3-5 μs]
│     Point source sum (sorted by ID): Σ GM_i/r_i²              [1-10 μs]
│     Subtract smooth contribution for massive sources (Sgr A*)  [0.1 μs]
│     Post-Newtonian corrections (if active)                     [0-10 μs]
│     OR: Kerr geodesic acceleration (if in GR mode)             [0-50 μs]
│
├── 1b. Propulsion
│     Lookup table: (throttle, ambient_pressure, fuel_state)     [1 μs]
│     → thrust_magnitude, Isp, heat_generation, exhaust_velocity
│     Thrust vector in body frame from engine positions
│     Transform to galactic frame: F_thrust = orientation × F_body
│     Relativistic bookkeeping: in the dp/dt formulation the LONGITUDINAL
│     3-force is boost-invariant (F_coord,∥ = F_proper,∥); the transverse
│     component transforms as F_coord,⊥ = F_proper,⊥/γ. No γ³ factor —
│     that belongs to the acceleration relation (a = F/γ³m), which the
│     momentum integrator already handles implicitly.
│
├── 1c. Aerodynamic forces (if in atmosphere)
│     Atmospheric density from model (NRLMSISE or equivalent)     [1-10 μs]
│     Angle of attack from velocity vs. orientation
│     Lookup table: (Mach, AoA) → (Cd, Cl, Cm)                  [1 μs]
│     (altitude enters via ρ in the force equations, not the table)
│     Drag = ½ρv²CdA, Lift = ½ρv²ClA
│     Aero heating = f(v, ρ) from characterization table
│
├── 1d. Other forces
│     Solar radiation pressure (if near a star)                   [0.1 μs]
│     ISM drag (if at relativistic speed): F = γ² ρ_ISM v² A     [0.1 μs]
│     (relativistic momentum flux, fully-stopped limit; γ² ≈ 5.3 at 0.9c —
│      partial-stopping fraction comes from the shield interaction table)
│     Tidal forces (if within activation threshold)               [0-5 μs]
│
├── 1e. Total force and torque
│     F_total = F_gravity + F_thrust + F_aero + F_radiation + F_ism + F_tidal
│     τ_total = Σ (r_i × F_i) for each force with a moment arm
│     (thrust offset from CoM, aero center of pressure offset, gravity gradient)

Phase 2: INTEGRATION (advance trajectory and attitude)
├── 2a. Trajectory (DOP853 step)
│     dp/dt = F_total
│     dx/dt = pc²/E (velocity from relativistic momentum)        [10 μs]
│     position += v × dt (double-double)
│     momentum += F × dt
│
├── 2b. Attitude (separate integrator, same dt)
│     dω/dt = I⁻¹ × (τ - ω × (I × ω))  (Euler's equation)     [2 μs]
│     dq/dt = ½ × q ⊗ (0, ω)  (quaternion derivative)
│     angular_velocity += I⁻¹ × (τ - ω × Iω) × dt
│     orientation = normalize(orientation + dq × dt)
│
├── 2c. Mass update
│     fuel_consumed = thrust / (Isp × g0) × dt
│     rest_mass -= fuel_consumed
│     Update component fill fractions
│     If mass change > 0.1%: flag MOI for recomputation

Phase 3: SHIP INTERNAL (subsystem evaluations)
├── 3a. Thermal network (implicit solve)                          [10 μs]
│     Heat sources this tick:
│       engine_heat (from 1b), solar_flux (from position),
│       aero_heating (from 1c), radiation_heating (from cached dose map),
│       reactor_heat (from reactor table), ISM_heating (from 1d)
│     (C - dt×A) × T_new = C × T_old + dt × q      [backward Euler; thermal.md]
│     (NOT "(I−dtA)T_new = T_old + dt(A·T_old + q)" — that form double-counts
│      conduction, integrating dT/dt = 2A·T + q. C = heat-capacity matrix.)
│     Sparse linear solve → new temperature per component
│
├── 3b. Structural evaluation                                     [5 μs]
│     Load vector: thrust loads + acceleration loads + 
│                  pressure differential + thermal stress(T_new) +
│                  tidal gradient (if active)
│     Displacements: u = K_reduced⁻¹ × load_vector
│     Stress: σ_e = S_e × u per element (precomputed stress-recovery matrices)
│     For each joint: if σ > yield_strength(T, dpa) → failure
│     If failure: trigger cascading failure (Section 3)
│
├── 3c. Resource bookkeeping                                      [1 μs]
│     Fuel flow: update tank fill fractions
│     Coolant flow: update coolant temperatures
│     Power balance: generation - load = battery ΔE
│     Consumables: O2 consumption, CO2 production, water usage
│
├── 3d. Atmosphere (if pressurized volumes exist)                  [2 μs]
│     dN_O2/dt = scrubber_generation - crew_consumption - leak_rate
│     dN_CO2/dt = crew_production - scrubber_removal
│     dN_H2O/dt = crew_production - condenser_removal
│     Pressure from ideal gas law: P = NkT/V
│     Leak rate: if hull breach, f(hole_size, pressure_differential)
│
├── 3e. Electrical                                                 [1 μs]
│     Generation: solar_panels(distance, orientation, degradation) +
│                 reactor_electrical + RTG_decay + fuel_cells
│     Load: life_support + computers + comms + heaters + pumps
│     Battery: charge += (generation - load) × dt
│     If battery = 0 and generation < load: brownout → systems fail by priority
│
├── 3f. Radiation dose accumulation                                [0.1 μs]
│     crew_dose += cached_dose_rate_at_crew × dt
│     electronics_dose[i] += cached_dose_rate_at_component[i] × dt
│     component.damage_dpa += cached_dpa_rate × dt

Phase 4: EVENT DETECTION                                           [1 μs]
├── 4a. Structural failure → trigger cascading failure (Section 3)
├── 4b. Lookup table boundary → pause timewarp
├── 4c. Atmosphere critical (O2 < 16%, CO2 > 4%, P < 50 kPa) → alarm
├── 4d. Temperature critical (any component > 90% melting point) → alarm
├── 4e. Radiation critical (dose rate > limit) → alarm
├── 4f. Electrical brownout → alarm, system shutdown by priority
├── 4g. Fuel exhaustion → engine shutdown
├── 4h. Collision detection with tracked objects → Section 4
├── 4i. Active source threshold crossing → update active source list
├── 4j. System entry/exit (dominant non-mean-field source changes) → Layer 2 generation

Phase 5: MOI UPDATE (if flagged)                                   [1-2 μs]
│     Recompute MOI tensor from component graph (analytical shapes + parallel axis)
│     Recompute center of mass

TOTAL PER TICK: ~35-55 μs (normal), ~50-100 μs (atmospheric flight with aero)
```

### Data flow between phases

Each phase reads from the state produced by previous phases and from PREVIOUS TICK values for circular dependencies:

| Subsystem | Reads THIS tick | Reads PREVIOUS tick |
|---|---|---|
| Gravity (1a) | position | -- |
| Propulsion (1b) | throttle, altitude, fuel | -- |
| Aero (1c) | position, velocity, orientation | -- |
| Thermal (3a) | heat sources from 1b, 1c, position | structural geometry (for conduction paths) |
| Structural (3b) | temperatures from 3a, loads from 2a | damage state |
| Resources (3c) | flow rates from 1b | -- |
| Atmosphere (3d) | leak rate from 3b, scrubber rate from 3c | -- |
| Electrical (3e) | position (solar), reactor state | -- |
| Radiation (3f) | -- | cached dose map (updated in background) |

One-tick lag for circular dependencies is physically negligible at 1000 Hz (1 ms lag for thermal ↔ structural coupling where both time constants are >> 1 ms).

---

## 3. CASCADING FAILURE MODEL

Failures propagate through the component graph physically, not as logic rules. Each failure mechanism feeds through real physics to cause downstream effects.

### Failure detection

Each tick, the structural evaluation (3b) checks every joint and component:

```
For each structural connection:
    σ = current_stress (from load vector × stiffness)
    σ_limit = yield_strength(T_component, dpa_component)
    
    if σ > σ_limit:
        connection.failed = true
        → propagate failure
```

Additionally:
- Temperature > melting_point → component destroyed
- Pressure > burst_pressure → pressure vessel rupture
- Radiation dose > electronics_failure_threshold → electronics failure

### Failure propagation

When a connection or component fails:

```
1. STRUCTURAL: Remove the connection from the graph.
   Check if the graph is still connected.
   If disconnected: ship splits (Section 4).

2. THERMAL: Failed structural connections lose thermal conduction path.
   If a heat pipe or coolant line crosses the failed connection: coolant flow stops.
   Thermal network topology changes → temperatures diverge.

3. PRESSURE: If the failed component is a pressure boundary:
   Hull breach → atmosphere leaks at rate f(hole_size, ΔP).
   Tank rupture → propellant vents. Possible explosion if propellant is hypergolic or cryogenic.

4. RESOURCE: If a resource line (fuel, coolant, power) crosses the failed connection:
   Flow stops. Downstream components lose supply.
   Engine without fuel: shutdown. Cabin without air supply: pressure drops.
   Electronics without power: off.

5. ELECTRICAL: If a power bus crosses the failed connection:
   Downstream loads lose power. Cascade through electrical tree.

6. SECONDARY STRUCTURAL: The loss of a component changes the mass distribution
   and structural load paths. Remaining connections may be overloaded by
   the redistributed forces → further failures.
```

This cascade is evaluated within a SINGLE tick. Failure detection → propagation → new failures → propagation → ... until no new failures occur. In practice, cascades resolve in 1-3 iterations because most failures don't create enough load redistribution to cause further failures immediately.

### Explosion model

Certain failures produce explosions:
- Hypergolic propellant contact (N2O4 + UDMH): instant combustion
- Cryogenic propellant release (LOX + fuel vapor): deflagration
- Battery thermal runaway: fire + toxic gas
- Reactor loss-of-coolant: core melt (no explosion, but radiation release)

Explosion effects:
- Pressure wave: damages adjacent components within blast radius
- Heat pulse: thermal damage to nearby components
- Fragmentation: shrapnel damages components in line-of-sight

These are modeled as instantaneous impulses applied to affected components, with damage proportional to distance and energy. Not a fluid dynamics simulation -- a blast model (Kingery-Bulmash or equivalent).

---

## 4. MULTI-OBJECT TRACKING

### Object types

| Type | Full sim? | Trajectory | Subsystems | Lifetime |
|---|---|---|---|---|
| Player vehicle | Yes | DOP853, full forces | All (thermal, structural, resources, atmosphere, electrical, radiation) | Permanent |
| Station | No | Keplerian + perturbations | Resources only (power, consumables depletion) | Permanent |
| Detached controllable (has power + computer + RCS) | Yes | Full | All | Until out of range or player switches |
| Debris (uncontrollable) | No | Keplerian | None | Until > threshold distance from any controllable craft |

### Ship splitting

When the component graph becomes disconnected:

```
1. Identify connected subgraphs (flood fill from each remaining node)
2. For each subgraph:
   a. Compute total mass, center of mass, velocity, angular momentum
   b. Apply conservation of momentum:
      Each piece gets v_piece = v_ship + ω × (com_piece - com_ship)
      Plus any impulse from the failure event (explosion, pressure release)
   c. Determine controllability:
      Has power source? Has computer? Has RCS/engines? Has radio?
      If yes to all: controllable. If no to any: debris.
   d. Create new tracked object with the appropriate fidelity level
   e. Build new octree for the piece (local rebuild from subgraph components)
3. The player remains attached to whichever piece contains the active crew station
4. If no piece has a crew station: game over (or switch to remote control if radio link exists)
```

### Docking

Docking requires:
- Two objects within docking range (~10 m)
- Compatible docking ports (defined as sub-components with type, position, orientation)
- Relative velocity < docking speed limit (~0.5 m/s)
- Angular alignment within tolerance (~5°)

On successful dock:
```
1. Merge component graphs (connect the two docking port sub-components)
2. Merge octrees (attach second object's octree to first)
3. Combine resource graphs (if ports support resource transfer)
4. Recompute MOI tensor for the combined object
5. The combined object is now one ship with one trajectory
```

### Collision

When bounding spheres of two tracked objects overlap:
```
1. Compute relative velocity at closest approach
2. Impact energy: E = ½ × m_reduced × v_relative²
3. Damage: distribute E among contacted components based on contact geometry
4. Apply impulse to both objects (momentum conservation)
5. If E > structural_failure_threshold for any component: trigger failure cascade
6. If either object is destroyed: remove from tracking
```

---

## 5. REGIME BRIDGING

The trajectory integrator uses different force models depending on the ship's environment. The FORCE FORMULA never changes (F = gravity + thrust + drag + ...), but HOW each term is computed varies.

### Active physics modules per regime

| Regime | Gravity | Propulsion | Drag | Thermal sources | Radiation sources |
|---|---|---|---|---|---|
| Ground/launch pad | Surface gravity | Chemical | None | Ambient | Ambient |
| Atmospheric flight | N-body + J2 | Chemical | Atmospheric (NRLMSISE + aero tables) | Aero heating + solar | Trapped belts + GCR |
| Low orbit | N-body + J2 + 1PN | Chemical/Electric | Residual atm (high alt) | Solar + albedo + IR | Trapped belts + GCR + SPE |
| Cislunar | N-body + 1PN | Chemical/Electric | None | Solar | GCR + SPE |
| Interplanetary | N-body + 1PN | Chemical/NTP/Ion | None | Solar (distance-dependent) | GCR + SPE |
| Near-star | N-body + 1PN | Any | Coronal plasma | Intense solar | Stellar radiation |
| Interstellar | Mean field + point sources | BH drive/ACMF/Sail | ISM (relativistic) | ISM heating + blueshifted CMB | GCR + ISM cascade |
| Near compact object | Kerr geodesics | Any | None (usually) | Accretion radiation | Extreme: jets, accretion disk |

### Regime detection

Each tick, evaluate conditions and activate/deactivate modules:

```rust
fn detect_regime(&self) -> ActiveModules {
    let mut modules = ActiveModules::default();
    
    // Always active
    modules.gravity_mean_field = true;
    modules.gravity_point_sources = true;
    
    // Atmospheric
    let alt = self.altitude_above_nearest_body();
    if alt < nearest_body.atmosphere_top {
        modules.atmospheric_drag = true;
        modules.aero_heating = true;
    }
    
    // Post-Newtonian
    for body in &self.active_sources {
        if body.gm / (self.distance_to(body) * C * C) > 1e-10 {
            modules.post_newtonian = true;
            break;
        }
    }
    
    // General relativity
    for body in &self.compact_objects_nearby {
        let r_s = 2.0 * body.gm / (C * C);
        if self.distance_to(body) < 100.0 * r_s && !self.gr_mode {
            modules.kerr_geodesics = true;
            // Convert state to Boyer-Lindquist
        }
        if self.distance_to(body) > 110.0 * r_s && self.gr_mode {
            modules.kerr_geodesics = false;
            // Convert state back to galactocentric
        }
    }
    
    // Special relativity (always on, but threshold for ISM interaction)
    let beta = self.speed() / C;
    if beta > 0.001 {
        modules.ism_drag = true;
        modules.relativistic_aberration = true;
    }
    
    // Tidal forces
    for body in &self.active_sources {
        let roche_dist = 10.0 * body.roche_radius(self.rest_mass);
        if self.distance_to(body) < roche_dist {
            modules.tidal_forces = true;
            break;
        }
    }
    
    modules
}
```

Transitions are smooth: each module's forces are negligible at its activation threshold, so enabling/disabling causes no discontinuity. GR transition uses hysteresis (100/110 r_s) per determinism.md.

### Ground contact

Minimal contact model — enough for launch, touchdown, and tip-over, not terramechanics:

- **Pad constraint:** while landed/pre-launch, the ship is kinematically attached to the
  body's rotating surface frame (position follows body rotation; velocity = surface velocity).
  Released when net upward acceleration (thrust − weight) > 0. This sidesteps stiff
  contact-force integration on the pad and gives exactly correct launch initial conditions
  (Earth rotation velocity included).
- **Touchdown:** each landing leg is a spring-damper (k, c per leg from the ship definition)
  against the local surface plane (body ellipsoid + Layer-4 local elevation). Contact force
  feeds the structural model; excessive impact velocity → leg failure via the normal
  stress-check path.
- **Standing stability:** static friction cone per leg (μ from surface type); tip-over check:
  gravity + thrust resultant vs. the support polygon of leg contact points.
- **Anything else** (sliding, soil deformation, crater terrain interaction) is out of scope;
  contact with terrain other than via legs or pad = collision damage (Section 4).

---

## 6. THE THREE INITIAL SHIPS

### Ship 1: Falcon 9 / Dragon

**Purpose:** Atmospheric flight through LEO. Validates against known performance data.

**Components (~60):**

First stage:
- LOX tank (cylinder, Al-2219, r=1.83m, l=16m, wall=4mm)
- RP-1 tank (cylinder, Al-2219, r=1.83m, l=8.5m, wall=4mm)
- Engine section (cone, Inconel-718, 9× Merlin 1D engines)
- Interstage (cylinder, CFRP, r=1.83m, l=4m, wall=2mm)
- Grid fins (4× panels, aluminum)
- Landing legs (4× struts, CFRP + aluminum)
- Helium pressurant bottles (spheres, COPV)

Second stage:
- LOX tank (cylinder, Al-2219)
- RP-1 tank (cylinder, Al-2219)
- Merlin Vacuum engine
- Payload fairing (2× half-shells, CFRP)

Dragon capsule:
- Pressure vessel (capsule shape, Al-2219, wall=8mm)
- Heat shield (disk, PICA-X)
- Trunk (cylinder, CFRP + solar panels)
- SuperDraco abort engines (8×)
- Crew stations (4×)
- Flight computer
- Nav sensors
- Life support system
- Docking port (forward)

**Characterization tables:**
- Merlin 1D: (throttle, altitude) → (thrust, Isp, heat) [2D, ~20×20 = 400 entries]
- Merlin Vacuum: (throttle) → (thrust, Isp, heat) [1D, ~20 entries]
- Aerodynamics: (Mach, AoA) → (Cd, Cl, Cm) [2D, ~50×20 = 1000 entries]
- Thermal: (heat_flux, duration) → (surface_temp, ablation_depth) [2D for PICA-X]

**Validation targets:**
- Max-Q at ~12 km altitude, ~80 seconds after launch
- MECO at ~160 seconds, altitude ~80 km, speed ~6000 km/h
- Second stage ignition and orbital insertion at ~8 minutes
- Payload to LEO: ~16,000 kg (expendable), ~11,000 kg (reusable)
- Entry heating: PICA-X surface temp ~1900°C during reentry

### Ship 2: Nuclear Thermal Mars Ship

**Purpose:** Interplanetary regime. NTP propulsion, reactor physics, radiation shielding.

**Components (~40):**

Propulsion section:
- NTP engine (NERVA-derivative, Isp ~900s, thrust ~110 kN)
- Reactor (U-235, graphite moderated, Cermet fuel elements)
- Shadow shield (LiH + tungsten, between reactor and crew)
- LH2 propellant tank (large, Al-2219, cryo-insulated)
- LH2 feed system (turbopump, valves)

Crew module:
- Pressure vessel (cylinder, Al-2219, wall=10mm)
- Radiation storm shelter (polyethylene-lined inner chamber)
- Life support system (CO2 scrubbers, O2 generation, water recycling)
- Crew stations (4-6×)
- Flight computer + nav sensors
- Communications array
- Solar panels (backup power, deployed in transit)

Mars entry system:
- Aeroshell (heat shield for Mars aerocapture)

**Characterization tables:**
- NTP engine: (throttle, LH2_temp, fuel_burnup) → (thrust, Isp, neutron_flux, heat) [3D]
- Reactor: (control_rods, coolant_temp, burnup) → (power, reactivity, decay_heat) [3D]
- Shadow shield: precomputed dose transmission factors (from radiation transport tool)
- Mars aero: (Mach, AoA, Mars_density) → (Cd, Cl, Cm) [3D]

**Validation targets:**
- Hohmann transfer ΔV: ~3.6 km/s per burn
- Transit time Earth-Mars: ~6-9 months
- Reactor shadow dose at crew module: < 50 mSv/year
- Mars aerocapture: peak deceleration ~5g, peak heating ~100 W/cm²

### Ship 3: BH Drive Starship

**Purpose:** Interstellar and compact object regimes. Exercises every physics module.

**Operating envelope:** BH mass is a free operational dial (hawking.md "Mass maneuvering") —
thrust from 55 MN (10⁹ kg, crewed cruise) to 21,800 MN (5×10⁷ kg, uncrewed sprint) on one
ship. The shell is sized for the minimum planned operating mass of the mission profile
(crewed band: shell thickened per the paper's muon-dose scaling; the sim enforces the
stopping-margin gate, not a hard interlock).

**Components (~30):**

Propulsion section (architecture from black_hole_paper; see hawking.md):
- SQM shell (CFL sphere, R = 2 m, 6×10⁷ kg at reference; absorbs the Hawking spectrum,
  re-emits keV pairs — there is no reflector and no aft radiation shield; the shell IS the shield)
- Anti-Helmholtz coil pair (superconducting, 1.36×10⁹ A, 50 GJ stored; quench = explosion event)
- Magnetic nozzle (820 m taper aft; exhaust exclusion axis)
- Railgun feed system + SQM-capped forward port (any-matter pellets)
- Bulk propellant store (grappled mass hopper — any matter; sized per mission, 0.7×-33× dry
  mass per the paper's Tsiolkovsky table; the sprint profile's store is asteroid-scale)
- PD confinement controller (critical electronics; loss → BH escapes in ~2.8 s → catastrophic)
- Forward shield (graded, for ISM bombardment only — from the ISM cascade research paper)

Forward shield:
- Graded shield: outer tungsten sacrificial layer + middle polyethylene + inner lead
- Designed to ablate under ISM bombardment over the mission duration
- Shield thickness and composition from the ISM cascade research paper

Crew module:
- Pressure vessel (cylinder, heavy shielding)
- Radiation storm shelter
- Full life support (closed loop for multi-year missions)
- Crew stations (2-4×)
- Redundant flight computers (radiation-hardened)
- Redundant nav sensors
- Long-range communications

Structural:
- Truss connecting propulsion to crew module (≥100 m crew standoff: muon punch-through dose)
- Radiators (sized for reactor/life-support/avionics heat only — the drive itself needs none:
  CFL firewall limits drive EM leakage to < 1 kW)

**Characterization tables:**
- BH drive: (BH_mass, shell_thickness) → (thrust, exhaust_β, system_Isp, punch_through, stopping_margin) [2D]
- BH feeding: (BH_mass, feed_rate) → (net_mass_rate, equilibrium_feed) [2D]
- ISM interaction: (velocity, shield_thickness, shield_material) → (dose_rate, ablation_rate, heating) [3D, from ISM cascade research]
- Forward shield: precomputed dose/ablation from radiation transport tool

**Validation targets:**
- Thrust = β·P_captured/c = β·0.83·P_H/c = 55.3 MN at M = 10⁹ kg (β = 0.332 Usov exhaust)
- Power budget closes: exhaust + neutrinos (~17%) + muon punch-through + EM leakage (<1 kW) = P_H
- Unfed lifetime 15.8 yr at 10⁹ kg; equilibrium feed 670 g/s
- Time dilation: proper time vs. coordinate time matches SR prediction
- S2 precession (if visiting Sgr A*): 12 arcmin/orbit

---

## 7. RADIATION TRANSPORT INTEGRATION

The radiation transport tool (radiation_transport.md) runs as a background computation, NOT per-tick. Its results are cached and used by the per-tick loop.

### When it runs

Triggers are deterministic conditions on sim state, evaluated in Phase 4 and scheduled via
the deferred-apply protocol (radiation_shielding.md): each trigger class has a fixed
sim-time apply delay, fixed primary count, and deterministic application tick.

- On system entry (new radiation environment: entering/leaving a magnetosphere, atmosphere, star system) [DRIFT class]
- On significant velocity change (ISM flux changes by >10%) [DRIFT]
- On reactor state change (startup, shutdown, power level change >10%) [DRIFT]
- On BH drive mass crossing a table band [DRIFT]
- On shield damage (octree change that affects shielding) [GEOMETRY]
- On SPE onset (new radiation source) [ACUTE]
- On an absolute sim-time grid during time acceleration (to catch gradual changes) [PERIODIC]

### What it provides

A cached dose-rate map:
- Dose rate (Gy/s and Sv/s) at each crew station
- SEE rate (events/s) at each electronics component
- dpa rate (dpa/s) at each structural component
- Energy deposition rate (W) at each thermal node (feeds into thermal as a heat source)

### Runtime cost

The full Monte Carlo run takes 10-60 wall-clock seconds (background thread). The per-tick dose accumulation is a single multiplication (rate × dt) = ~0.1 μs. The active map is replaced at each job's predetermined sim-time T_apply.

During time acceleration: if the ship is in steady-state (constant velocity, constant shielding, no events), only PERIODIC jobs fire (every 86,400 s sim time). If the environment is changing, DRIFT/GEOMETRY jobs fire as their conditions trip. At high warp, jobs fall due faster in wall-clock terms than the MC runtime and the tick loop blocks at each T_apply — transport throughput therefore sets the effective WALL-CLOCK rate in radiation-active regimes. Granted dt never depends on wall-clock (radiation_shielding.md): a slower machine runs the identical tick sequence more slowly. Wall-clock speed changes how long the loop blocks, never the trajectory.

---

## 8. SAVE/LOAD

### What is saved

```rust
struct SaveState {
    // Ship
    trajectory: TrajectoryState,
    component_graph: ComponentGraph,    // includes all component states, damage, fill fractions
    subsystem_state: SubsystemState,

    // Other tracked objects
    tracked_objects: Vec<TrackedObject>,

    // Metadata
    save_time: DoubleDouble,           // coordinate time
    universe_seed: u64,                // for regeneration
    ship_id: String,                   // which ship definition
}
```

**Damage state is per-component parameterized, never raw voxels.** The octree is always a
derived view (spatial.md), so voxel-level damage from the transport engine must round-trip
through compact per-component fields that the octree can be regenerated from EXACTLY:

```rust
struct ComponentDamage {
    ablation_depth: Grid2D<f32>,   // depth map over the exposed surface (low-res, ~32×32)
    dpa_profile: [f32; 8],         // dpa vs depth through the wall (piecewise-linear)
    dpa_mean: f64,                 // for property evaluation at component granularity
    activation: Vec<(Isotope, f64)>, // transmutation inventory summary (activity, Bq)
    holes: Vec<HoleDesc>,          // discrete breaches: position, area, shape class
}
```

The transport engine's dynamic-geometry mutations WRITE THROUGH this parameterization
(radiation_transport.md Section 7): ablating a voxel updates the depth map cell; the octree
region is then re-rasterized from the map. Quantization from the parameterization is far
below Monte Carlo statistical noise, and save/load reproduces the octree bit-exactly by
re-rasterizing from the same fields. Save file size stays ~1-10 MB.

The universe (stars, planets, biospheres, ISM) is NOT saved. It regenerates deterministically from (universe_seed, coordinate_time).

### Load procedure

```
1. Load SaveState from file
2. Set coordinate time and proper time
3. Reconstruct ship from component graph + per-component damage parameterization
   (re-rasterize octree, rebuild thermal network, structural model)
4. Regenerate the local universe (nearby stars, current system) from seed + T
5. Rebuild active source list from ship position
6. Resume simulation
```

---

## 9. ATTITUDE CONTROL AND GUIDANCE

The control layer behind interface.md's `attitude`, `maneuver`, and script commands. Scope
(year 1): attitude hold + burn execution + OPEN-LOOP scripted profiles. No closed-loop launch
guidance, powered-descent guidance, aerocapture steering, or rendezvous autopilot — those are
flown manually or by script.

### Attitude controller

Quaternion-error PD control, run in the per-tick loop (Phase 1e computes the commanded torque):

```
q_err   = q_target ⊗ q_current⁻¹          (take the sign giving the short way around)
τ_cmd   = −Kp × axis(q_err) × angle(q_err) − Kd × ω     (per-axis gains from ship definition)
```

**Actuator allocation:** reaction wheels first (τ up to wheel torque limit, momentum
accumulates in wheel state); RCS thrusters when wheels saturate or for translation; automatic
wheel desaturation burns RCS against stored momentum when |h_wheel| > 90% capacity. RCS
allocation solves a fixed (precomputed) thruster-selection table: commanded torque direction
→ thruster set + duty cycle. All gains and tables are per-ship-definition constants —
deterministic, no adaptive tuning.

**Hold modes** (target quaternion source, updated per tick): inertial (fixed q), prograde /
retrograde (±v̂ direction), normal / antinormal (±(r̂×v̂)), radial in/out, target / anti-target
(tracked object), sun-pointing, custom vector. Mode changes are player commands (inputs, so
physics-affecting by design).

### Maneuver executor

`maneuver <direction> <dv> [at <time>]`:

```
1. At t_start − t_align: command attitude to burn direction (t_align from ship's slew rate)
2. At t_start: throttle up per the ship's throttle profile
3. Integrate achieved Δv from accelerometer state (thrust/mass, gravity excluded)
4. Cutoff when Δv_achieved ≥ Δv_target (tail-off correction: subtract the known
   shutdown-transient impulse from the target)
```

Timed events are scheduled in sim time (deterministic). Burn timing for finite-burn accuracy
(splitting Δv across the node) is the player's/script's job in year 1.

### Trajectory predictor (display path)

Maneuver planning and the orbit window need predicted trajectories. The predictor runs the
SAME force assembly as the physics path, at looser tolerance (rtol ~1e-6), seeded from
current state + planned maneuvers — but in the DISPLAY path: its results render plots and
node markers and can never write physics state (determinism.md). Prediction horizon bounded
by a compute budget (~50 ms per replan).

### Scripted guidance profiles (open loop)

The script system (interface.md) drives repeatable ascent/entry profiles — the script
IS the guidance program, physics does the rest:

```
# falcon9_ascent.script (excerpt)
at T+0        throttle 1.0
wait altitude > 100m ; attitude custom(pitch_program)   # pitch vs. time table in script
wait speed > 2200 ; event MECO ; throttle 0
...
```

Script `wait`/`at` conditions are evaluated against sim state in Phase 4 — deterministic.
Profiles for the three ships (F9 ascent, NTP TMI burn sequence, BH-drive cruise program) ship
with their definitions and double as validation fixtures.

---

## 10. TIMEWARP VALIDITY LIMITS

"Unlimited time acceleration" means unlimited SIM-RATE REQUEST — the engine grants what the
active physics allows. Each subsystem declares a max dt; the integrator honors the minimum:

| Constraint | Max dt | Reason |
|---|---|---|
| Atmospheric flight | ~1-10 ms | Aero force stiffness (existing behavior) |
| Attitude dynamics, ω ≠ 0 | 0.1/\|ω\| | Explicit attitude integration accuracy; ABOVE this, if torques ≈ 0, the engine switches to ANALYTIC torque-free propagation (closed-form for axisymmetric MOI) instead of clamping dt |
| Orbit (bound, period P) | ~P/1000 | Integration accuracy per orbit (DOP853 tolerance governs; this is the practical ceiling) |
| Thermal network | unlimited (implicit) | Radiation linearization: fixed 2-iteration Newton per step (deterministic); accuracy degrades gracefully at extreme dt |
| Reactor point kinetics | unlimited (implicit) | Stiff-stable |
| Interstellar cruise | ~10⁶-10⁷ s | Mean-field smoothness; active-source refresh distance |
| Swept coverage | (refresh distance)/v | Between refreshes the trajectory segment is swept-sphere tested against star system spheres of influence — a system cannot be jumped over inside one step |
| Radiation jobs | blocks at T_apply | Transport throughput sets the effective wall-clock rate in radiation-active regimes; granted dt is never wall-clock-dependent (Section 7, radiation_shielding.md) |

Consequences worth stating honestly: in LEO the per-tick cost (~40 μs) at dt ~ 5 s caps
effective warp at ~10⁵× wall-clock on one core — warp in low orbit is COMPUTE-BOUND, not
unlimited. In interstellar cruise, dt ~ 10⁷ s makes 10 Myr ≈ 3×10⁷ ticks ≈ tens of minutes of
wall clock, punctuated by radiation-job blocks. The HUD shows requested vs. granted warp and
which constraint is binding.
