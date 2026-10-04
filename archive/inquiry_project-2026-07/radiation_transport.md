# Radiation Transport Engine

Standalone Monte Carlo particle transport code in Rust. Simulates radiation propagation through spacecraft geometry -- from external environment (cosmic rays, solar particle events, trapped radiation belts) and onboard sources (fission reactors, BH drives, antimatter engines) through arbitrary shielding to crew locations, electronics, and sensitive components.

Outputs: dose maps, dose-equivalent rates (Sievert), single-event effect rates (bit flips/day), material degradation rates (dpa), energy deposition per component, time-resolved shield evolution under sustained bombardment.

This is a standalone crate with its own repository, API, and paper. The flight simulator depends on it as an external library. The simulator's spatial material system (Section 2) is shared infrastructure -- the octree and material database are used by the radiation tool AND by the simulator's thermal, structural, and resource modules.

**Determinism:** All Monte Carlo transport uses counter-based PRNG (each primary particle gets RNG stream = hash(base_seed, particle_index)). Tallies accumulated in fixed voxel-index order. Results are bit-identical across runs, platforms, and thread counts. See [determinism.md](determinism.md).

---

## 1. RESEARCH SCOPE

### What the paper presents

A Monte Carlo particle transport code with mutable octree geometry and Geant4-derived interaction tables, enabling time-resolved radiation-material co-evolution. Validated against Geant4 on spacecraft-relevant benchmarks. Applied to problems no existing tool can compute: shield lifetime under sustained bombardment, optimal graded shield design, and reactor shadow shield degradation.

### What's novel

1. **Radiation-material co-evolution.** Every existing transport code treats geometry as immutable during a simulation. Coupled transport-depletion codes (SERPENT, OpenMC) update isotopic COMPOSITION but not material SHAPE, DENSITY, or MECHANICAL PROPERTIES. This tool's octree is mutable: ablation removes material, sputtering thins surfaces, radiation damage changes density and strength, structural failure opens gaps. The transport and the geometry evolve together in a feedback loop. This extends existing coupled methods from composition-only to full material state.

2. **Octree-native geometry.** No existing transport code uses octree geometry. The octree provides O(log N) spatial queries, natural adaptive resolution for thin walls, and the voxel-level mutability needed for dynamic geometry.

3. **Geant4-derived interaction tables.** Rather than reimplementing nuclear physics models (Bertini cascade, evaporation, fragmentation), pre-run Geant4 offline for every (particle type, energy, material) combination. Store secondary particle distributions as lookup tables. Runtime Rust code samples from these tables. Physics accuracy comes from Geant4 (validated); performance comes from Rust + octree.

4. **Unified external + onboard source handling.** No existing spacecraft radiation tool handles the space radiation environment AND onboard sources (reactor, BH drive, antimatter engine) in one simulation.

### What's NOT novel (acknowledged prior work)

- Monte Carlo particle transport algorithms (established since the 1940s)
- Transport-depletion coupling for isotopic evolution (SERPENT, OpenMC, SCALE)
- Transport-thermal coupling (OpenMC + MOOSE/Cardinal)
- Space radiation environment models (Badhwar-O'Neill, AP-9/AE-9)
- Cross-section databases (ENDF, TENDL, PDG)
- Damage-property correlations (decades of experimental and theoretical work)

The contribution is extending coupled transport from composition-only to full material state (ablation, property degradation, structural failure), enabled by the mutable octree, with application to problems that existing tools cannot address.

---

## 2. SPATIAL MATERIAL SYSTEM

This is shared infrastructure between the radiation transport tool and the flight simulator. The radiation tool uses Layers 2 and 3 (material database + octree). The simulator uses all layers.

### Layer 1: Component Graph (source of truth)

The ship is defined as a list of components with geometry, material, and connections. Human-editable. Everything else is derived from it.

```
Component {
    id: "forward_shield"
    shape: Disk { radius: 5.0m, thickness: 0.3m }
    position: (0, 0, 45.0)
    orientation: (0, 0, 1)
    material: "tungsten-rhenium"
    
    connections: [
        Structural("shield_mount", joint: "bolted_ring"),
        Thermal("shield_mount", contact_area: 1.2m²),
    ]
}

Component {
    id: "crew_cabin"
    shape: Cylinder { radius: 2.5m, length: 8.0m }
    position: (0, 0, 20.0)
    orientation: (0, 0, 1)
    material: "al-2219-t87"
    wall_thickness: 8.0mm
    contents: Atmosphere { composition: "cabin_air", pressure: 101.3kPa }
    
    connections: [
        Structural("forward_truss"),
        Structural("aft_truss"),
        Thermal("forward_truss", contact_area: 0.8m²),
        Thermal("aft_truss", contact_area: 0.8m²),
        Resource("air_supply", flow: "cabin_air"),
        Electrical("main_bus"),
    ]
    
    sub_components: [
        Electronics { id: "flight_computer", position: (0.5, 0.3, 22.0), device: "rad750" },
        Electronics { id: "nav_sensor", position: (0, 0, 24.5), device: "star_tracker" },
        CrewStation { id: "pilot_seat", position: (0, 0, 23.0) },
    ]
}
```

A Falcon 9: ~50-100 components. BH drive starship: ~30-50. Each component is a geometric primitive (cylinder, cone, sphere, box, disk, torus, panel) at a known position with material and wall thickness.

Sub-components mark locations within a component for specific tallies: electronics (SEE rates), crew stations (dose), sensors (degradation).

### Layer 2: Material Database

Every material the simulator knows. Properties are functions of state (temperature T, displacement damage dpa, transmutation fraction):

```
Material {
    id: "al-2219-t87"
    name: "Aluminum 2219-T87"
    
    // Composition (for radiation cross-sections)
    composition: [(Al, 0.935), (Cu, 0.063), (Mn, 0.003)]
    nominal_density: 2.84  // g/cm³
    
    // Thermal properties (functions of T, dpa)
    thermal_conductivity: |T, dpa| 120.0 * (1.0 - 0.05 * dpa) * thermal_k_curve(T)  // W/(m·K)
    specific_heat: |T| 864.0 + 0.15 * T    // J/(kg·K)
    emissivity: |T| 0.09 + 0.0001 * T
    melting_point: 916.0   // K
    heat_of_vaporization: 10.7e6  // J/kg
    
    // Structural properties (functions of T, dpa)
    youngs_modulus: |T, dpa| 73.0e9 * (1.0 - 0.0003 * (T - 293.0)) * (1.0 - 0.02 * dpa)  // Pa
    yield_strength: |T, dpa| 352.0e6 * (1.0 - 0.001 * (T - 293.0)) * (1.0 - 0.03 * dpa)  // Pa
    fracture_toughness: |T, dpa| 26.0e6 * (1.0 - 0.05 * dpa)  // Pa·√m
    poisson_ratio: 0.33
    
    // Radiation damage response
    swelling: |dpa, T| 0.001 * dpa * exp(-0.001 * (T - 400.0).powi(2))  // fractional
    nrt_displacement_energy: 25.0  // eV
    sputter_yield_table: "al2219_sputter.bin"  // from SRIM
    
    // Interaction tables (from Geant4)
    interaction_tables: "al2219_interactions.bin"
}
```

~50-100 materials covers all spacecraft: aluminum alloys, steels, titanium, carbon composites, polyethylene, tungsten, inconel, LiH, borated PE, quartz, kapton, aerogel, regolith, SQM, propellants (LOX, RP-1, LH2), etc.

The interaction tables (from Geant4) store: for each incoming (particle_type, energy_bin), the distribution of outgoing secondaries (types, energies, angles). Pre-computed offline by running Geant4 for each material. Stored as binary files, loaded at startup.

### Layer 3: Voxel Octree (derived from component graph)

Generated by rasterizing the component graph's geometric primitives into a spatial octree.

**Per leaf voxel:**

```
struct Voxel {
    material_id: u16,     // index into material database (0 = vacuum)
    density_frac: f16,    // fraction of nominal density (1.0 pristine, <1.0 damaged, 0 void)
    damage_dpa: f16,      // cumulative displacement damage
}
// 6 bytes per leaf
```

Temperature is NOT stored per-voxel. It comes from the thermal module via the component graph: voxel position → parent component → component temperature. When a physics module needs T at a voxel, it does this lookup.

**Resolution:** Adaptive. ~1 mm at material boundaries (walls, joints), ~10 cm in bulk material, meters in empty space. For a 50 m ship: ~10-50 million leaves × 6 bytes = **60-300 MB.**

**Construction from component graph:**

```
1. Create empty octree spanning the ship's bounding box
2. For each component: rasterize its shape into the octree
   - Surface shell: material_id = component.material, density_frac = 1.0
   - Interior: material_id = component.contents (propellant, atmosphere, vacuum)
   - Wall thickness determines how many voxels deep the shell is
   - At boundaries, refine the octree to ~1 mm to capture thin walls
3. For overlapping components: inner component takes precedence (it's enclosed)
4. Mark sub-component locations (electronics, crew) in a separate index
```

**Regeneration triggers:**
- Ablation thins a wall → component wall_thickness decreases → re-rasterize that component's region
- Structural failure → component marked broken → remove its material from octree
- Propellant depletion → tank fill_fraction decreases → update interior voxel densities
- Reconfiguration (jettison a stage, deploy a shield) → remove/add components, rebuild affected region

Regeneration is LOCAL, not global. Only the affected component's octree region is rebuilt. Cost: milliseconds.

### Layer 4: Thermal Network (derived from component graph)

Each component is one or more thermal nodes. Connections define conduction paths. Not used by the radiation tool directly, but the tool's energy deposition feeds in as a heat source, and temperature from the thermal network affects material properties used by the tool.

### Layer 5: Structural Model (derived from component graph)

Components connected by structural joints. Stiffness matrix precomputed. Not used by the radiation tool directly, but radiation damage (dpa) affects structural properties, and structural failure feeds back as a geometry change.

### Layer 6: Resource Graph (subset of component graph)

Tanks, pipes, valves, engines. Resource levels affect the octree (propellant density in tanks).

### Coupling between layers

```
Component Graph (source of truth)
    │
    ├──→ Octree ←──────────────────────────────────────┐
    │     │                                              │
    │     └──→ Radiation Transport                       │
    │           │                                        │
    │           ├── dose map → Thermal (heat source)     │
    │           ├── dpa map → Material properties        │
    │           ├── ablation → Component wall update ────┘
    │           └── SEE → Electronics failure
    │
    ├──→ Thermal Network
    │     └── temperature → Material properties → Octree (density from swelling)
    │
    ├──→ Structural Model
    │     └── failure → Component graph update → Octree rebuild
    │
    └──→ Resource Graph
          └── fill levels → Octree (tank density)
```

All modifications go through the component graph. The octree is always a derived view. Consistency is guaranteed because there's one source of truth.

---

## 3. PHYSICS SCOPE

### Particle types

All particles relevant to spacecraft radiation, from 1 MeV to 10 TeV:

| Particle | Role |
|---|---|
| Protons | Primary GCR, SPE, trapped belts, ISM at relativistic speeds |
| Neutrons | Secondary from cascades, reactor source, critical for SEE and crew dose |
| Gammas | Secondary from nuclear reactions, reactor source, BH drive emission |
| Electrons/positrons | Trapped belts, pair production, BH Hawking radiation |
| Pions (π±, π⁰) | Produced above ~300 MeV/nucleon. π⁰ → 2γ. Critical above 0.3c ISM |
| Muons (μ±) | From pion decay. Deeply penetrating. Dominant dose at >0.5c ISM |
| Heavy ions (He-Fe) | ~1% of GCR flux, ~30% of dose-equivalent due to high LET |
| Kaons, hyperons | Minor contribution above ~1 GeV. Included for completeness |

### Interaction physics

Two-layer approach: analytical electromagnetic physics in Rust + element-level hadronic tables from Geant4. Any material defined by elemental composition and density works automatically.

**Layer 1: Electromagnetic processes (implemented analytically in Rust)**

These are well-documented formulas parameterized by atomic number Z and material properties. No tables needed. Covers ~99% of particle steps (EM interactions are far more frequent than nuclear).

- **Bethe-Bloch stopping power** for charged particles: dE/dx from the standard formula with density corrections, shell corrections, and Barkas effect. ~50 lines. Validated against NIST PSTAR/ASTAR/ESTAR to <1%.
- **Compton scattering:** Klein-Nishina cross-section. Exact QED formula. ~30 lines.
- **Pair production:** Bethe-Heitler cross-section with screening corrections. ~30 lines.
- **Photoelectric absorption:** parameterized cross-sections from Scofield tables (tabulated per element, ~1 KB per element). ~20 lines.
- **Bremsstrahlung:** Seltzer-Berger differential cross-sections. ~40 lines.
- **Multiple Coulomb scattering:** Highland formula for angular deflection. ~10 lines.

Total EM implementation: ~200-300 lines of core physics + ~200 lines of support code. These formulas work for ANY material composition -- just weight by atom fraction.

**Layer 2: Hadronic processes (element-level tables from Geant4)**

Nuclear reactions (spallation, fragmentation, capture) are the hard physics where cascade models are needed. Tabulated per ELEMENT, not per material.

For each (element, incident_particle_type, energy_bin):
- Nuclear interaction cross-section (barns)
- Interaction outcome: secondary particle distributions (types, energies, angles)

At runtime, when a nuclear interaction occurs:
1. Pick which element in the material is hit (probability proportional to atom fraction × cross-section)
2. Look up the outcome table for that element at the particle's energy
3. Sample secondaries from the stored distribution

**Table generation (offline, Geant4 in single-interaction mode):**

```
For each element (H, C, N, O, Al, Si, Ti, Cr, Fe, Ni, Cu, W, Pb, ...):  // ~30 elements
    For each particle type (p, n, γ, π±, μ±, α, C, O, Si, Fe):           // ~10 types
        For each energy bin (1 MeV to 10 TeV, ~100 log-spaced bins):
            Run 10,000 Geant4 single-interaction events
            (thin-target mode: record ONLY the first nuclear interaction)
            Record: cross-section and secondary (type, energy, angle) distributions
            Store as compressed binary table
```

Critical: **thin-target mode** (single interaction). Not thick-slab. The transport code handles the geometry; the table provides only the physics of one nuclear collision.

Table size: ~30 elements × 10 particle types × 100 energy bins × ~1 KB per entry = **~30 MB**
baseline. RESOLUTION RISK: ~1 KB per entry is coarse for a double-differential secondary
distribution (species × energy × angle) against the 10%-of-Geant4 acceptance bar. Plan:
staged refinement — (1) ship the 1 KB baseline, (2) after the first slab benchmark, measure
which channels dominate the dose error, (3) refine only those entries (dominant channels for
common elements may need 10-100 KB; sparse storage keeps the total ≲ 300 MB even if 10% of
entries refine). Generation: ~330 million Geant4 events at baseline. At ~0.01-0.1 s per
hadronic event, **1-10 days on a small cluster** (embarrassingly parallel across elements and
energies); a laptop-scale fallback exists (start with H, C, O, Al, Fe + protons/neutrons only
— enough for the slab benchmarks and the first ships — and grow the table set incrementally).

**Versatility:** Any material is supported by defining its elemental composition and density. Novel shield materials, exotic alloys, even lunar regolith -- specify the atom fractions and the transport works. No per-material precomputation needed.

### Decay

Unstable particles decay during transport:
- π± → μ± + ν: cτ = 7.8 m. Probability of decay before interaction sampled per step.
- μ± → e± + ν + ν̄: cτ = 660 m. Most muons decay before stopping.
- n → p + e⁻ + ν̄: τ = 880 s. Only relevant for trapped thermal neutrons.

Decay kinematics are simple (two-body or three-body) and implemented directly in Rust, not table-based.

---

## 4. RADIATION SOURCES

### External environment

**Galactic Cosmic Rays (GCR):** Badhwar-O'Neill 2020 model. Continuous, isotropic, modulated by solar cycle. ~4 particles/cm²/s at solar minimum. ~87% protons, ~12% helium, ~1% heavier.

**Solar Particle Events (SPE):** ESP or SAPPHIRE model. Sporadic, duration hours-days. 10³-10⁵ protons/cm²/s above 10 MeV. In the simulator: timing and intensity seed-determined per star.

**Trapped radiation belts:** AP-9/AE-9 models. Position-dependent (L-shell, latitude). For any planet with a magnetic field: scale from dipole moment.

**Interstellar medium at relativistic speeds:** ISM proton flux in ship frame = n_ISM × v_ship. At 0.5c: ~7.5 × 10⁹ protons/cm²/s, each at ~145 MeV. At 0.9c: ~1.2 GeV/proton. At 0.99c: ~5.7 GeV/proton.

### Onboard sources

**Fission reactor:** Watt fission spectrum (peak ~0.7 MeV, tail to ~10 MeV). Source strength proportional to thermal power.

**Black hole drive:** Hawking spectrum from hawking.md. Particle mix depends on BH mass.

**Antimatter annihilation:** ~5 pions per p-p̄ annihilation. π⁰ → 2γ (70 MeV each). π± → μ± + ν.

### Source API

```rust
trait RadiationSource {
    fn sample_particle(&self, rng: &mut Rng) -> (ParticleType, Position, Direction, Energy, Weight);
    fn total_rate(&self) -> f64;  // particles/s
}
```

All sources implement this interface. The transport code is source-agnostic.

---

## 5. TRANSPORT ALGORITHM

### Core loop

```rust
fn transport_batch(octree: &mut Octree, sources: &[&dyn RadiationSource], n_primaries: u64, tallies: &mut Tallies) {
    for _ in 0..n_primaries {
        let (ptype, pos, dir, energy, weight) = source.sample_particle(&mut rng);
        let mut stack = vec![Particle::new(ptype, pos, dir, energy, weight)];
        
        while let Some(mut p) = stack.pop() {
            while p.alive {
                let voxel = octree.lookup(p.position);
                if voxel.is_vacuum() { 
                    // Step to octree boundary or escape
                    if !octree.step_through_vacuum(&mut p) { break; }
                    continue;
                }
                
                let mat = &materials[voxel.material_id];
                let tables = &mat.interaction_tables;
                
                let mfp = tables.mean_free_path(p.ptype, p.energy);
                let dedx = tables.stopping_power(p.ptype, p.energy);
                let dist_to_exit = octree.distance_to_exit(p.position, p.direction);
                
                if dist_to_exit < mfp {
                    // Cross voxel boundary
                    let de = dedx * dist_to_exit;
                    tallies.deposit(voxel, de, p.weight);
                    tallies.deposit_dpa(voxel, de, mat.nrt_displacement_energy);
                    p.energy -= de;
                    p.position += p.direction * dist_to_exit;
                } else {
                    // Interaction occurs
                    p.position += p.direction * mfp;
                    let de_continuous = dedx * mfp;
                    tallies.deposit(voxel, de_continuous, p.weight);
                    
                    let secondaries = tables.sample_interaction(p.ptype, p.energy, &mut rng);
                    for s in secondaries {
                        if s.energy > cutoff(s.ptype) {
                            stack.push(s.with_weight(p.weight));
                        } else {
                            tallies.deposit(voxel, s.energy, p.weight);
                        }
                    }
                    p.alive = false; // primary consumed by interaction
                }
                
                // Check decay for unstable particles
                if p.ptype.is_unstable() && rng.gen::<f64>() < p.decay_probability(dist_traveled) {
                    let products = p.decay(&mut rng);
                    stack.extend(products.into_iter().map(|s| s.with_weight(p.weight)));
                    p.alive = false;
                }
                
                if p.energy < cutoff(p.ptype) {
                    tallies.deposit(voxel, p.energy, p.weight);
                    p.alive = false;
                }
            }
        }
    }
}
```

### Variance reduction

Standard techniques, applied per-voxel based on importance:
- **Weight windows:** target weight range per octree region. Auto-split (overweight) or roulette (underweight).
- **Source biasing:** preferentially emit toward regions of interest.
- **Forced collision:** in thin regions, force interaction to occur, adjust weight.

### Energy cutoffs

| Particle | Cutoff | Reason |
|---|---|---|
| Protons | 1 MeV | Range < 0.02 mm in Al |
| Neutrons | 0.025 eV (thermal) | Can still capture and produce gammas |
| Gammas | 10 keV | Photoelectric absorption dominates |
| Electrons | 100 keV | Range < 0.1 mm |
| Muons | 1 MeV | Same as protons |
| Heavy ions | 1 MeV/nucleon | Range negligible |

---

## 6. TALLIES

### Per-voxel tallies (accumulated over a batch)

- **Energy deposition** (J): total ionizing dose. Convert to Gray: D = E / m_voxel.
- **Displacement damage** (dpa): from NIEL. dpa = E_NIEL / (N_atoms × E_d) where E_d is the NRT displacement energy.
- **Nuclear reactions**: count of each reaction type, for tracking transmutation.

### Per-sub-component tallies

- **Crew dose** (Sv): dose-equivalent at crew station locations, using quality factors Q(LET) from ICRP 103.
- **SEE rates**: particle flux spectrum Φ(E) at each electronics location × device cross-section σ_SEE(E). Output: upsets/bit/day, latchups/device/day.
- **Material degradation**: cumulative dpa at structural joints, sensors, etc.

### Global tallies

- **Radiation field**: particle flux spectrum at any query point.
- **Activation inventory**: radioactive isotopes created by transmutation, their activity (Bq) and decay radiation.

---

## 7. DYNAMIC GEOMETRY

The octree is mutable. Transport modifies it. This is the core capability that enables novel physics computations.

### Modification mechanisms

After each transport batch, evaluate accumulated tallies for material changes:

**Ablation (coupled to thermal relief — energy-per-batch alone is NOT sufficient).** Whether
a surface ablates depends on the power balance, not just deposited energy: radiative cooling
(εσT⁴) and conduction into the bulk relieve the surface between impacts. The batch loop
therefore inserts a thermal-relaxation step (this is the "particle-cascade + thermal-transport
+ ablation" coupling that the research paper's novelty claim rests on):

```
1. Transport batch → volumetric heating profile q(depth) per surface patch [W/m³ at the
   batch's flux-implied power, see time mapping below]
2. 1D thermal solve per patch: T(depth) from q(depth), radiative BC at surface,
   conduction into bulk (material k(T, dpa))
3. Ablate where the energy balance fails:
   if T_surface > T_vaporization: recession rate = (q_net − q_radiated − q_conducted)
                                                    / (ρ × h_vap)   [m/s]
4. Write recession into the parent component's ablation_depth map (ship.md Section 8 —
   damage is per-component parameterized; the octree re-rasterizes from the map)
```

**Explicit fluence↔time mapping:** a batch of N primaries sampled from a source of rate R
particles/s represents Δt_macro = N_effective/R of exposure (N_effective = N/weight-sum for
variance-reduced runs). Recession, dpa, and transmutation from the batch are applied over
Δt_macro; the macro-timestep is chosen so per-step geometry change is small (recession <
1 voxel, Δdpa < 0.1). This is how microsecond-scale particle physics and month-scale erosion
coexist in one simulation.

**Sputtering.** Below ablation threshold, atomic ejection from surface impacts:

```
atoms_removed = incident_flux × sputter_yield(particle, energy, material)
mass_removed = atoms_removed × atomic_mass
voxel.density_frac -= mass_removed / voxel_mass
if voxel.density_frac < 0.01:
    voxel.material_id = VACUUM
```

**Property degradation.** Radiation damage accumulates per-voxel:

```
voxel.damage_dpa += dpa_this_batch
// Material properties now evaluate differently:
// yield_strength(T, voxel.damage_dpa) < yield_strength(T, 0)
// thermal_conductivity(T, voxel.damage_dpa) < thermal_conductivity(T, 0)
```

**Transmutation.** Nuclear reactions change elemental composition:

```
for reaction in nuclear_reactions_this_batch:
    voxel.transmutation_fraction += reaction.yield
    if transmutation_fraction > threshold:
        // Create an activation source at this voxel
        activation_sources.push(ActivationSource {
            position: voxel.center,
            isotope: reaction.product,
            activity: reaction.rate × decay_constant,
        })
```

**Structural failure (when coupled to simulator).** Radiation damage weakens material. When the structural model detects failure:

```
if yield_strength(T, dpa) < current_stress:
    component.failed = true
    // Octree removes this component's material
    // Opens a radiation pathway through the gap
```

### The feedback loop

```
1. Initialize octree from component graph
2. Define radiation sources
3. Transport batch (N primaries, ~10-60 seconds) → Δt_macro = N_effective / source_rate
4. Tally: energy deposition, dpa, nuclear reactions per voxel
5. Thermal relaxation: 1D solve per surface patch (heating profile vs. radiative +
   conductive relief) → surface temperatures, recession rates
6. Update damage state (per-component parameterization, ship.md Section 8):
   a. Ablation recession × Δt_macro → ablation_depth map (where energy balance fails)
   b. Sputtering recession × Δt_macro → same depth map (fractional; sub-voxel accumulation)
   c. dpa accumulation → dpa depth profiles
   d. Transmutation products → activation inventory (new sources)
   e. (If coupled) Check structural failure, update component graph
   Then re-rasterize affected octree regions from the updated parameterization.
7. If geometry changed significantly: go to 3 (Δt_macro chosen so per-step change is small:
   recession < 1 voxel, Δdpa < 0.1)
8. If geometry stable (or exposure target reached): report results
```

The accumulated batches give the time evolution of the shield state on the macro clock while
each batch resolves the full cascade physics — the two-timescale structure is the core of the
method.

### What this enables

**1. Shield lifetime under sustained bombardment.** How long does a shield last at 0.9c through the ISM? The answer is not (total energy / heat of vaporization) because the energy deposition PROFILE changes as the shield thins -- cascade products deposit energy deeper as the outer layer ablates, and underlayers may have been pre-damaged by penetrating secondaries. The feedback loop captures this.

**2. Optimal graded shield design.** What sequence of materials (outer sacrificial tungsten, middle neutron-moderating polyethylene, inner gamma-shielding lead) minimizes time-integrated crew dose per kg of total shield mass? The cascade physics couples the layers: neutrons from the outer layer moderated by the middle layer, but moderation efficiency changes as the outer layer ablates. Only solvable with dynamic simulation.

**3. Reactor shadow shield degradation.** Over a multi-year mission, reactor neutrons displace hydrogen from LiH shield material, reducing moderation capability. The dynamic simulation tracks this: dpa accumulates in the shield, hydrogen content decreases (transmutation), neutron transmission increases, crew dose rises. Quantifies the timeline from "safe" to "exceeds dose limits."

**4. Cascading failure under acute events.** SPE → dose spike → electronics SEE → attitude control failure → ship tumbles → new surfaces exposed → more damage → structural weakening → breach → atmosphere loss. Each step physically computed through the evolving octree.

**5. BH drive radiation environment during mass feeding cycles.** Hawking spectrum shifts as BH mass changes. Shielding response differs per spectral phase. Dynamic simulation tracks the radiation environment throughout the drive cycle.

---

## 8. PERFORMANCE

### Targets

| Metric | Geant4 | This tool |
|---|---|---|
| Setup time | Hours (C++ code, geometry, physics config) | Minutes (define octree, select sources) |
| 10⁶ primaries | ~10 minutes | ~10-60 seconds |
| Dose accuracy | 5-10% vs. experiment | 5-10% (validated against Geant4) |
| Memory | ~1 GB | ~100-500 MB |
| Dynamic geometry | Not possible | Native |

Speedup sources: octree traversal, Rust memory management, focused physics (table-based, no on-the-fly model computation), no visualization/UI overhead.

### Runtime integration with simulator

The transport tool does NOT run every tick. Scheduling follows the DETERMINISTIC
deferred-apply protocol (radiation_shielding.md): triggers are sim-state conditions, jobs
snapshot their inputs, and results become physics-visible at predetermined sim times —
never on wall-clock availability.
- Result is cached as a dose-rate map with a scheduled T_apply
- Per-tick: the simulator looks up the cached dose rate and accumulates. Microseconds.
- Dynamic geometry (ablation) is extrapolated linearly between applied maps (rate × dt)

---

## 9. VALIDATION

### Benchmarks

1. **Al slab, monoenergetic protons (50-1000 MeV).** Dose vs. depth. Compare to Geant4.
2. **Al sphere, GCR spectrum.** Dose-equivalent at center. Compare to Geant4 and OLTARIS.
3. **Multi-layer slab (Al + PE + Al).** Neutron moderation validation. Compare to Geant4 and Zeitlin+ 2006 experimental data.
4. **ISS module (simplified cylinder).** Compare to ISS dosimetry measurements.
5. **Heavy ion fragmentation (C, Fe beams in Al).** Compare fragment spectra to Geant4 and experimental data.
6. **High-energy cascade (>1 GeV protons).** Pion/muon production. Validates ISM cascade regime.
7. **Dynamic benchmark: slab ablation under sustained proton bombardment.** Run Geant4 on fixed slab, manually step through ablation, compare to the tool's automated dynamic result. Novel benchmark -- validates the feedback loop itself.

### Acceptance criteria

- Energy deposition within 10% of Geant4
- Dose-equivalent within 15% of Geant4
- Fragment spectra within 20% for heavy ions
- Dynamic simulation: ablation timeline matches manual Geant4 stepping to within 15%

---

## 10. DEVELOPMENT PLAN

(The months below are RELATIVE to this subproject's start; the canonical calendar placement
is SCOPE.md's gate-based timeline — Phase 3, months 5-9, with element-table generation
started as background compute in Phase 1.)

### Month 1: Core engine + EM physics

**Week 1-2: Octree geometry engine**
- Octree data structure with adaptive refinement
- Ray marching (step through voxels along a direction)
- Neighbor finding, distance-to-exit computation
- Construction from component graph (rasterize primitives)
- Mutability: modify individual voxels, local rebuild
- Tests: ray-trace through known geometry, verify material sequences

**Week 2-3: Photon + charged particle transport**
- Monte Carlo transport loop (the core algorithm)
- Table-based interaction sampling
- Continuous energy loss (Bethe-Bloch from tables)
- Photon processes from tables (Compton, pair production, photoelectric)
- Proton transport from tables
- Tallies: energy deposition per voxel, dpa
- Generate Geant4 interaction tables for the first ELEMENT set (H, C, O, Al — element-level
  tables per Section 3; aluminum alloys and polyethylene are then compositions of these,
  with no per-material precomputation)

**Week 3-4: Validation + dynamic geometry**
- Validate proton and photon transport against Geant4 slab benchmarks
- Implement ablation mechanism (remove voxels above vaporization threshold)
- Implement sputtering (density reduction)
- Implement the feedback loop
- First dynamic benchmark: Al slab under sustained proton bombardment

### Month 2: Full particle set + sources + paper-ready

**Week 5-6: Remaining particles + environment models**
- Neutron transport from tables (including thermal capture)
- Electron/positron transport from tables
- Pion and muon transport + decay (kinematics in Rust, interactions from tables)
- Heavy ion transport from tables (fragmentation outcomes from Geant4 tables)
- GCR, SPE, trapped belt environment models
- ISM source model (relativistic proton flux from velocity + density)
- Generate Geant4 interaction tables for all materials in the database

**Week 7-8: Integration, optimization, validation**
- Onboard source API (reactor, BH drive, antimatter -- spectrum inputs from external modules)
- SEE rate calculator
- Variance reduction (weight windows, source biasing)
- Full validation suite (all 7 benchmarks)
- Performance profiling and optimization
- Activation tracking (transmutation products as secondary sources)
- Property degradation coupling (dpa → material property updates)

**Deliverables at end of month 2:**
- Working standalone crate
- Validated against Geant4 on spacecraft-relevant benchmarks
- Dynamic geometry demonstrated (slab ablation under bombardment)
- Ready for paper writing and integration with simulator

---

## 11. CRATE STRUCTURE

```
radiation-transport/
├── Cargo.toml
├── src/
│   ├── lib.rs                    // Public API
│   │
│   ├── spatial/
│   │   ├── mod.rs
│   │   ├── octree.rs             // Octree: build, query, mutate
│   │   ├── ray_march.rs          // Step through voxels along a ray
│   │   ├── voxel.rs              // Voxel struct (material_id, density_frac, dpa)
│   │   └── component_graph.rs    // Ship definition → octree construction
│   │
│   ├── materials/
│   │   ├── mod.rs
│   │   ├── database.rs           // Material properties, constitutive relations
│   │   ├── interaction_tables.rs // Load and sample Geant4-derived tables
│   │   └── damage.rs             // Damage-property correlations
│   │
│   ├── transport/
│   │   ├── mod.rs
│   │   ├── engine.rs             // Core Monte Carlo loop
│   │   ├── particle.rs           // Particle types, kinematics, decay
│   │   └── variance.rs           // Weight windows, splitting, roulette
│   │
│   ├── sources/
│   │   ├── mod.rs
│   │   ├── gcr.rs                // Galactic cosmic rays (Badhwar-O'Neill)
│   │   ├── spe.rs                // Solar particle events
│   │   ├── trapped.rs            // Radiation belts (AP-9/AE-9)
│   │   ├── ism.rs                // Interstellar medium at relativistic speeds
│   │   └── onboard.rs            // User-defined (reactor, BH, antimatter)
│   │
│   ├── tallies/
│   │   ├── mod.rs
│   │   ├── dose.rs               // Energy deposition, Gray, Sievert
│   │   ├── dpa.rs                // Displacement damage
│   │   ├── see.rs                // Single-event effect rates
│   │   └── activation.rs         // Transmutation products, activity
│   │
│   └── dynamic/
│       ├── mod.rs
│       ├── ablation.rs           // Voxel removal from thermal overload
│       ├── sputtering.rs         // Surface mass loss
│       ├── property_update.rs    // dpa → material property changes
│       └── feedback.rs           // The outer loop: transport → update → transport
│
├── tables/                       // Geant4-generated interaction data (per ELEMENT — §3)
│   ├── generate_tables.py        // Geant4 script to produce tables
│   ├── H/                        // Per-element directories; materials are
│   ├── C/                        // runtime compositions of elements
│   ├── O/
│   ├── Al/
│   ├── Fe/
│   ├── W/
│   └── ...
│
├── data/
│   ├── materials.toml            // Material database definitions
│   └── devices/                  // SEE cross-section data for electronics
│
├── examples/
│   ├── slab_benchmark.rs         // Simple slab dose profile
│   ├── spacecraft.rs             // Full spacecraft with crew + electronics
│   ├── ism_cascade.rs            // Relativistic ISM bombardment
│   ├── dynamic_ablation.rs       // Shield ablation under sustained flux
│   └── graded_shield_opt.rs      // Find optimal graded shield design
│
├── benches/
│   └── transport_throughput.rs   // Particles/second benchmarks
│
└── tests/
    ├── geant4_comparison/        // Automated comparison against Geant4 results
    ├── octree_tests.rs
    ├── transport_tests.rs
    └── dynamic_tests.rs
```
