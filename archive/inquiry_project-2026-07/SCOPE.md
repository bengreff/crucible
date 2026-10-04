# Space Physics Computational Laboratory
## Scope and Research Plan

---

## WHAT IT IS

A terminal-based physics simulator for spaceflight. You fly ships through a physically accurate universe from a command line. The simulator handles everything from atmospheric launches to relativistic interstellar cruise to orbits around black holes, within a galaxy of ~130 billion procedurally generated stars (the count is DERIVED from the Milky Way mass model + IMF, not chosen — physics.md Section 3.2.9), each with planetary systems, biospheres, and minor bodies. No graphics. No game mechanics. Physics in, telemetry out.

The player launches from Earth, transfers between ships, accelerates to relativistic speeds, visits other star systems, inspects planets for life, navigates near black holes, and returns home -- finding that millions of years have passed. **The player is an operator, not a passenger:** timewarp is the operator's time machine, and the operator is not subject to mortality. Crewed missions are a bounded MODE with real constraints (life support, consumables, dose limits) -- physically honest for nearby-star trips of years to decades of proper time. Galaxy-scale voyages, where even relativistic proper time runs to millennia, fly uncrewed; it is the SHIP that returns to a changed Earth. Everything the player encounters is physically computed, not scripted. Every star has a position derived from galactic dynamics. Every planet has an atmosphere derived from formation physics. Every biosphere exists because the conditions were right and the dice (from a deterministic seed) came up in favor of life.

## WHAT IT IS NOT

- Not a game (no graphics engine, no real-time rendering, no player progression)
- Not an engineering design tool (ships are predefined, not designed by the user)
- Not a rebuild of existing software (not GMAT, not KSP, not Universe Sandbox)
- Not a detection simulator (the universe is omniscient -- all objects visible regardless of detectability)

## WHY IT DOESN'T EXIST

Every existing tool covers one or two physics domains. Trajectory tools (GMAT, REBOUND) don't do thermal or structural. Engineering tools (COMSOL, ANSYS) don't do orbital mechanics. GR codes (Einstein Toolkit) don't have spacecraft. Galaxy models (Galaxia, Besançon) don't have ships or biospheres. Nobody has built the connective tissue between these domains. The novelty is regime bridging -- one simulator that smoothly transitions between Newtonian, relativistic, and general-relativistic physics as conditions demand, within a galaxy-scale procedural universe where every object is physically grounded.

---

## CORE ARCHITECTURE

### Language

**Runtime: Rust.** Two C-library FFI dependencies: NAIF SPICE (solar system ephemerides) and NRLMSISE-00 (Earth atmosphere model). Everything else is custom Rust.

**Offline: Python/C++.** Ship characterization tool uses Cantera, Geant4, OpenMC, FEniCS, galpy, astroquery, CoolProp to precompute lookup tables. LMC halo-response run (EXP/AGAMA). Gaia catalog preparation.

### Coordinate System

Galactocentric double-double (two f64s per coordinate, ~31 significant digits). Sub-nanometer precision at 30 kpc from the galactic center. One frame for all trajectory integration -- no frame switching. Ship internal octree uses ship-local f64 (femtometer precision at 100 m). Time: double-double seconds from J2000.0 (sub-nanosecond precision at 10 million years). See [determinism.md](determinism.md).

### Determinism

Bit-exact reproducibility across runs, platforms, and query histories. All math via portable `libm` (no system math library). All summations in fixed order (sorted by permanent ID). No parallel physics computation. Physics path and display path completely separated -- player queries cannot affect the ship's trajectory. Procedural generation uses independent hash-derived seeds per object (no shared mutable RNG). See [determinism.md](determinism.md).

### Two Tools

**Ship Characterization Tool (offline):** Takes a ship definition (component graph with shapes, materials, connections) and precomputes all physics as lookup tables. Uses arbitrarily expensive methods (CFD, Monte Carlo, FEA) because it runs once per ship. See [characterization.md](characterization.md).

**Flight Simulator (real-time):** Computes trajectory, ship internals, and procedural generation in real time. ALL subsystem physics at runtime is lookup-table interpolation -- no runtime chemical kinetics, no runtime transport, no runtime CFD. The only runtime ODEs are three tiny fixed-size systems: the thermal network (implicit, unconditionally stable), the cabin atmosphere (~10 variables per volume), and reactor point kinetics (6 delayed-neutron groups + feedback -- a static table cannot produce rod-insertion transients). See [ship.md](ship.md).

---

## THE UNIVERSE

### Timescale and Epoch

10 million years from J2000.0. Sufficient to cross the galaxy at 0.01c. Time acceleration is unbounded in REQUEST, physics-limited in GRANT: every subsystem declares its maximum valid dt and the engine honors the minimum, so all subsystems always compute at full fidelity (per-regime ceilings, compute-bound warp in low orbit, and radiation-job blocking are specified in ship.md Section 10; the HUD shows requested vs. granted and the binding constraint). Timewarp pauses automatically when ship state reaches the boundary of precomputed lookup tables (structural failure, atmospheric entry, novel configuration). See [ship.md](ship.md).

### Gravity

One formula everywhere: mean field gradient + individual point sources above an adaptive threshold − double-count correction for massive sources. Same formula in LEO, interplanetary space, interstellar space, and at the galactic center. The force law is the weak-field geodesic form, exact in velocity to O(Φ/c²): dp/dt = −γm(1+β²)∇Φ with velocity map v = (pc²/E)(1+2Φ/c²) — Newtonian at rest, correct photon-limit deflection at β→1. Static Φ² and EIH source-motion corrections activate near massive bodies (Mercury precession). Full Kerr geodesic integration within 100 Schwarzschild radii of compact objects (hysteresis exit at 110 r_s; the transition is an integrator-restart event). State = (position, canonical momentum). Works at all speeds from rest to 0.99c with the same integrator. See [physics.md](physics.md).

### Galactic Mean Field

Hybrid analytic + perturbation potential of the Milky Way: McMillan 2017 axisymmetric potential-density pairs + analytic rotating bar (38 km/s/kpc, 27° at T=0) + Cox-Gomez spiral (23 km/s/kpc) + warp — all pinned to observed parameters BY CONSTRUCTION — plus a low-order basis-function expansion for the one piece analytics can't do: the LMC's halo wake and reflex response (offline constrained N-body via EXP/AGAMA, coefficients stored over the 10 Myr window). The density side of the same potential-density pairs drives star generation: self-consistent by construction. Cluster-scale density features (170 globular clusters, 50 massive open clusters) orbit within the field as analytical profiles. Nuclear structures (Sgr A* as point source, NSC, NSD, CMZ) added at the galactic center. Parameters from McMillan 2017, Bland-Hawthorn & Gerhard 2016, GRAVITY Collaboration, Gaia DR3. ~18 MB precomputed data, ~1-3 μs per evaluation. See [physics.md](physics.md) Section 3.1.

### ~130 Billion Stars (derived, not chosen)

Procedurally generated from deterministic seeds using the galactic density model (thin disk, thick disk, bulge, halo, spiral arm modulation). ~1 million real stars from Gaia DR3 provide ground truth (real stars override procedural in their vicinity). Stars get mass (Kroupa IMF), age, spectral type, luminosity, metallicity, binary companions (environment-dependent rates), and analytical orbits (epicyclic for disk, Keplerian for galactic center). Stellar events (flares, supernovae) are seed-determined. Free-floating brown dwarfs and rogue planets exist under a separate query filter. Stars are generated on demand (~55 ns each; a 10,000-star query refresh costs well under a millisecond) and discarded. No persistent state. See [generation.md](generation.md).

### Planetary Systems

Generated on demand when inspected or entered. Planet occurrence from Kepler/TESS/RV surveys (metallicity-dependent giant planet rates, peas-in-a-pod spacing). Full system architecture: planets with mass-radius relations (Chen & Kipping 2017), moons with tidal heating, asteroid belts with individual bodies > 1 km (deterministically named), ring systems as density fields, stellar wind with Parker spiral, debris disks and Oort clouds clipped to the local environment's encounter survival limit. Known exoplanet systems (~5500) are AUGMENTED -- real planets kept, procedural planets added in stable gaps. See [generation.md](generation.md) Section 4.

### Biospheres

Habitability from Kopparapu+ 2013 HZ boundaries, planet mass/atmosphere/magnetic field requirements. Life stages from a sequential Poisson model: abiogenesis (λ = 1.5 Gyr^-1, from Kipping 2020), photosynthesis, oxygenation, eukaryogenesis (λ = 0.3, Carter 2008 hard step), multicellularity, complex life. No intelligent life. Ecosystem sketches: biome types, photosynthetic pigment color (stellar-spectrum-dependent), biosignature outputs (O2, CH4, vegetation spectral edge). Expected: 2-10 billion microbial worlds, 20-200 million complex biospheres (scaled to the derived star count). Tidally heated moons with subsurface oceans are a separate habitability pathway. See [generation.md](generation.md) Section 6.

### Environmental Fields

Evaluated per-tick (~0.5 μs): ISM density and phase (5 phases from cold molecular to hot ionized, from galactic gas disk model), stellar wind (within the astropause), nebulae (HII regions around O/B stars, SNR from recent supernovae, planetary nebulae, molecular clouds in spiral arms), radiation (GCR modulated by solar wind, trapped belts near magnetic planets, stellar radiation). See [generation.md](generation.md) Section 9.

---

## THE SHIP

### Spatial Material System

Ships are defined as component graphs: lists of geometric primitives (cylinders, cones, spheres) with materials, positions, and connections (structural, thermal, resource, electrical). This human-editable definition is the single source of truth. From it, the engine derives: a voxel octree (for radiation transport and "what's at point X?" queries), a thermal network, a structural model, a resource graph, and an aerodynamic surface mesh. When the ship changes (ablation, damage, fuel depletion), the component graph updates and derived views regenerate locally. See [spatial.md](spatial.md).

### Material Database

~50-100 materials with properties as functions of temperature and radiation damage: thermal conductivity k(T, dpa), yield strength σ_y(T, dpa), etc. Interaction tables from Geant4 (element-level, thin-target mode) for radiation transport. Analytical EM physics (Bethe-Bloch, Compton, pair production) implemented in Rust. Any material defined by elemental composition and density works automatically. See [spatial.md](spatial.md), [radiation_transport.md](radiation_transport.md).

### Subsystem Physics

ALL subsystems are lookup-table-based at runtime. The characterization tool (offline) runs the expensive physics (Cantera for combustion, OpenMC for reactor neutronics, Geant4 for radiation, FEniCS for FEA) and produces tables. At runtime, each subsystem is a table interpolation (~1 μs), plus the three tiny fixed-size ODEs (thermal network — implicit, unconditionally stable, ~10 μs; cabin atmosphere; reactor point kinetics). Total per-tick cost: ~30-55 μs typical, up to ~115 μs with all conditional modules active (canonical budget: physics.md Section 7). See subsystem documents below.

### Full Rotational Dynamics

Orientation (quaternion), angular velocity (3-axis), moment of inertia tensor. Gravity gradient torque, thrust offset torque, aerodynamic torque. MOI recomputed from component graph on mass change (analytical shapes + parallel axis theorem). Attitude control via reaction wheels, RCS thrusters, or gravity gradient stabilization.

### Cascading Failure

Failures propagate physically through the component graph. Structural failure → thermal path change → pressure boundary breach → atmosphere leak → electrical disruption. Each step is computed through real physics (thermal conduction, structural stress, gas dynamics), not logic rules. Explosions (hypergolic contact, battery thermal runaway) produce blast damage to adjacent components. See [ship.md](ship.md) Section 3.

### Multi-Object Tracking

Multiple vehicles tracked simultaneously. Player-controlled ships get full subsystem simulation. Space stations get trajectory + resource depletion. Detached debris gets trajectory only, deleted when beyond threshold distance from any controllable craft. Docking via analytical contact geometry (port alignment, velocity/angular tolerance). Collisions via bounding sphere detection + impulse + damage. Ship splitting on structural disconnection: each piece becomes a tracked object, momentum conserved. See [ship.md](ship.md) Section 4.

### Radiation Transport Engine

Standalone Rust crate. Monte Carlo particle transport through the ship's voxel octree. Analytical EM physics + element-level hadronic tables from Geant4. Dynamic geometry: the octree mutates under bombardment (ablation removes voxels, sputtering thins surfaces, radiation damage degrades material properties, transmutation creates activation sources). Feedback loop: transport → damage → geometry update → transport on modified geometry. Enables novel computations: shield lifetime under sustained relativistic ISM bombardment, optimal graded shield design, reactor shadow shield degradation. Runs as a background thread; per-tick dose accumulation from cached results (~0.1 μs). See [radiation_transport.md](radiation_transport.md).

---

## INTERFACE

Terminal with command mode (type commands), flight mode (real-time keypresses for manual control), and script mode (execute command sequences from files). Telemetry HUD with trajectory, propulsion, environment, and ship status. Window system for map, orbit, thermal, radiation, structural, and system views. Query system for searching stars/bodies with filters. Alerts for critical events. Save/load: ship state + coordinate time only (~1-10 MB), universe regenerates from seed. See [interface.md](interface.md).

---

## SHIPS

### Initial 3

**Falcon 9 / Dragon:** Chemical bipropellant (LOX/RP-1). ~60 components. Atmospheric flight through LEO. Validates against known performance data (max-Q timing, payload to LEO, reentry heating). See [combustion.md](combustion.md).

**Nuclear Thermal Mars Ship:** NERVA-class NTP (Isp ~900s). ~40 components. Interplanetary regime. Reactor shadow shield degradation, Mars aerocapture. See [nuclear.md](nuclear.md).

**BH Drive Starship:** Hawking radiation propulsion. ~30 components. Interstellar and compact object regimes. Exercises every exotic physics module. Forward shield ablates under ISM bombardment. See [hawking.md](hawking.md).

### Future Ships (20+ planned)

Laser sail, antimatter-catalyzed fusion, nuclear pulse (Orion), ion drive, solar sail, Starship, space station, mass driver payload. Each novel propulsion concept is a research paper that produces a new ship for the simulator.

---

## DESIGN DOCUMENTS

| Document | Lines | Covers |
|---|---|---|
| [SCOPE.md](SCOPE.md) | 273 | This file. Vision, architecture, ships, research, canonical timeline. |
| [VISION.md](VISION.md) | 209 | The complete vision in one read + synthesis findings (F1-F7). |
| [AUDIT.md](AUDIT.md) | 359 | Audit record: all Tier-1/2/3 findings and their resolutions; strategic plan. |
| [physics.md](physics.md) | 1043 | Gravity model, hybrid mean field, Milky Way structure + object census, trajectory, velocity-complete weak-field law, GR, communication. |
| [generation.md](generation.md) | 1190 | ~130B stars (derived), planetary systems + secular ephemerides, biospheres, minor bodies, environmental fields, query system, engine loop. |
| [determinism.md](determinism.md) | 384 | Coordinate system, floating-point rules, physics/display separation, active source management, testing. |
| [spatial.md](spatial.md) | 458 | Component graph, material database, voxel octree, surface mesh. Shared by radiation tool and simulator. |
| [ship.md](ship.md) | 799 | State vector, per-tick loop, cascading failure, multi-object tracking, docking, 3 ships, radiation integration, save/load, ground contact, attitude control + guidance, timewarp limits. |
| [radiation_transport.md](radiation_transport.md) | 741 | Standalone Monte Carlo transport: dynamic octree, element tables + analytical EM, thermal-relief co-evolution loop, validation, crate structure. |
| [characterization.md](characterization.md) | 181 | Offline pipeline: Cantera, OpenMC, Geant4, FEniCS → lookup tables. Table format, decomposition strategy. |
| [combustion.md](combustion.md) | 96 | Chemical propulsion: reaction tables, nozzle tables, thermal tables. Cantera characterization. |
| [nuclear.md](nuclear.md) | 157 | Fission reactors (runtime point kinetics), fusion/ACMF engines, antimatter systems, decay heat. |
| [hawking.md](hawking.md) | 235 | BH drive per black_hole_paper: SQM shell converter, Usov exhaust, mass maneuvering, confinement failure modes, radiation tables. |
| [thermal.md](thermal.md) | 116 | Thermal network: implicit solver, heat sources, radiation nonlinearity, FEniCS conductances. |
| [structural.md](structural.md) | 98 | Reduced stiffness + stress-recovery matrices, failure modes (yield, fracture, buckling, fatigue, creep), thermal stress. |
| [resources.md](resources.md) | 95 | Resource flow: propellant, power, consumables. Conservation laws, load shedding priority. |
| [atmosphere_ship.md](atmosphere_ship.md) | 91 | Cabin gas composition ODE, leak modeling (choked/subsonic), alarm thresholds, multi-compartment exchange. |
| [electrical.md](electrical.md) | 76 | Power generation (solar, RTG, reactor, fuel cells), distribution, brownout cascade. |
| [radiation_shielding.md](radiation_shielding.md) | 156 | Ship-side radiation module: deferred-apply scheduling protocol, dose tracking (crew Sv, electronics TID, material dpa). |
| [interface.md](interface.md) | 206 | Command mode, flight mode, scripting, telemetry HUD, windows, query system, alerts. |
| **Total** | **~6,963** | **20 documents covering the complete architecture.** |

---

## ORIGINAL RESEARCH

### 1. Relativistic ISM interaction above the pion threshold

Coupled particle-cascade + thermal-transport + ablation model for shielding at 0.3-0.9c. The transition from Bethe-Bloch stopping to hadronic cascade regime is uncharacterized. Geant4 + thermal coupling. Produces BH drive forward shield tables. The radiation transport tool with dynamic octree enables the novel computation: time-resolved shield lifetime under sustained bombardment, not just static cascade profiles.

### 2. Laser sail acceleration to 0.1c

Optimal sail material and geometry for a 1000 kg payload. Three coupled problems: thermal limits under extreme flux, ride stability, beam diffraction terminating acceleration. Complete acceleration profile with material and geometry optimization. Produces a laser sail ship.

### 3. Antimatter-catalyzed fusion ignition thresholds

Minimum antiprotons per D-T/D-He3 pellet for reliable ignition. Geant4 cascade in compressed fuel + ICF ignition threshold. Specific engine design with Isp, thrust, antimatter requirements. Produces the most feasible near-term interstellar engine.

### 4. Antimatter production (angle TBD)

Laser-driven production, optimal space-based factory, or Hawking radiation harvesting. Determines whether ACMF propulsion is feasible on human timescales.

### 5. Candidates for future development

Magnetic sail braking at relativistic speeds. Relativistic gravity assist near compact objects. Dark forest stealth thermodynamics. Dyson swarm minimum seed mass.

---

## VALIDATION TARGETS

| Test | Expected result | Validates |
|---|---|---|
| Mercury perihelion precession | 43"/century | Post-Newtonian corrections |
| Hohmann transfer Earth→Mars | ΔV ~ 3.6 km/s per burn | Trajectory integration |
| Falcon 9 to LEO | Max-Q at ~80s, MECO at ~160s, payload ~16t to LEO | Chemical propulsion + atmospheric flight |
| Schwarzschild ISCO | r = 6GM/c² | GR geodesic integration |
| S2 orbital precession | ~12 arcmin/orbit | GR + galactic center dynamics |
| Solar system barycenter | Wobbles ~1 solar radius from Sun center | Analytic reflex bookkeeping + SPICE blend (generation.md §4.10) |
| Local stellar density | ~0.1 stars/pc³ | Procedural generation calibration |
| Rotation curve at R₀ | 233 ± 3 km/s | Mean field potential |
| Energy conservation (static-field test) | E conserved to integrator tolerance in a frozen mean field; Jacobi integral conserved in the bar-corotating frame | Trajectory + propulsion coupling (total energy is NOT conserved in the time-dependent mean field -- the rotating bar does work; the test must use a static or corotating configuration) |
| Photon-limit deflection | 4GM/bc² at β→1 | Velocity-complete weak-field force law |
| Determinism: cross-platform | Bit-identical trajectories on x86 and ARM | Floating-point rules |
| Determinism: query-independence | Same trajectory regardless of player queries | Physics/display separation |

---

## TIMELINE (re-baselined 2026-07-04 — gate-based, ~15 months, band 13-20)

**Planning basis:** scope is fixed, time flexes. Developer at ~20 h/wk with heavy AI-assisted
implementation. Build tasks compress ~3-5×; three things do NOT compress and set the floor:
(1) the developer's own verification hours — every physics module waits on human review,
(2) debug cycles against physical ground truth (validation targets),
(3) laptop-only offline compute (Geant4 table generation and reference runs, LMC constrained
N-body) — these are CALENDAR items started early and run in the background.

**Rules:** No phase begins until the previous gate is green. A gate slip triggers re-planning
at the gate (defer the next phase's optional items), never silent schedule slip. The
determinism CI suite runs on every commit from Gate 1 onward — nothing merges that breaks it.

### Phase 0 — Spec repair (DONE 2026-07-04)
All audit findings (7 Tier-1, 10 Tier-2, 10 Tier-3) resolved in the design docs. AUDIT.md is
the record.

### Phase 1 — Determinism kernel + gravity spine (months 1-2)
Rust scaffold, double-double math, libm policy, seed/hash infrastructure, cross-platform CI
harness. Custom DOP853 with canonical-momentum state and the velocity-complete weak-field
force law. Analytic MW potential (McMillan + bar + Cox-Gomez + warp). SPICE ingest + secular
blend. **Background compute starts now:** staged Geant4 element tables (H, C, O, Al, Fe ×
p, n first) and the LMC constrained-N-body run.
**GATE 1:** bit-identical reference trajectory on x86 + ARM; Mercury 43″/century; rotation
curve 233±3 km/s; S2 precession; photon-deflection unit test 4GM/bc²; SPICE↔secular blend
continuous at the boundary.

### Phase 2 — Ship kernel + Falcon 9 (months 3-5)
Component graph → octree → thermal network → structural (with stress recovery) → resources /
electrical / atmosphere; ground contact; attitude control + burn executor; characterization
v1 (Cantera + analytic aero correlations — no CFD); terminal interface v1 (telemetry HUD,
command mode, script mode).
**GATE 2:** Falcon 9 launch-to-LEO validation green (max-Q ~80 s, MECO ~160 s, payload
15-17 t); scripted ascent profile reproducible; determinism suite green with all subsystems
active, including query-independence.

### Phase 3 — Radiation transport engine + NTP ship (months 5-9, overlaps Phase 2 tail)
Standalone MC crate: octree transport, analytic EM physics, table-driven hadronics
(staged element set from Phase 1 background compute), variance reduction. Benchmark ladder
1-5 vs. Geant4, then dynamic geometry + thermal-relief co-evolution loop + dynamic benchmark
7. Deferred-apply integration with the simulator. NTP ship (runtime point kinetics, shadow
shield, Mars aerocapture).
**GATE 3:** benchmarks 1-5 within acceptance (10-15%); dynamic ablation benchmark within
15%; NTP shadow-shield dose scenario < 50 mSv/yr; full element table set generating or done.
**Paper 1 (shield-flux co-evolution at 0.3-0.9c) drafts from month 8, submits ~month 10-11.**

### Phase 4 — Galaxy + relativistic regime + BH ship (months 9-12)
Procedural star generation (guiding-center cells, material arms, unwinding queries), Gaia
overlay (100 pc complete + anchors + sampled shells), query system, planetary systems,
biospheres, environmental fields. SR cruise + Kerr integrator with event-boundary handoff.
BH drive ship with paper-derived tables (black_hole_paper is authoritative).
**GATE 4:** 10,000-body queries < 1 s; population checks (IMF fractions, local density
0.1/pc³, binary fraction); Schwarzschild ISCO; BH-drive cruise with time dilation matching
SR; system entry generating deterministic Layer-2 architecture.

### Phase 5 — Integration + full validation (months 12-15)
Multi-object tracking, docking, ship splitting, cascading failure with explosions. The
narrative loop end-to-end at full fidelity: launch from Earth → transfer → relativistic
cruise → another star system → biosphere inspection → compact-object approach → return to a
changed Earth. Save/load round-trip bit-exactness. Interface polish. **Paper 2 (velocity-
complete weak-field formulation or deterministic procedural galactic dynamics) drafts here.**
**GATE 5 (success bar):** every row of the validation table green; determinism suite green
end-to-end; the full mission flown start to finish without manual intervention in physics.

### Risk band and year-2 items
The ±30% band (13-20 months) is dominated by Phase 3 (transport validation debug cycles +
laptop compute) and Phase 5 unknown-unknowns (integration shakedown). Explicitly YEAR 2+:
laser sail and ACMF papers + ships, antimatter production study, the 20+ future ships, CFD-
grade aero tables, IMBHs. These are scope-preserving deferrals, not cuts — nothing in year 1
blocks them.
