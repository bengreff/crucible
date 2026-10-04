# Ship Characterization Tool

Offline Python/C++ pipeline that takes a ship definition (component graph) and produces all lookup tables the runtime simulator needs. Runs once per ship design, not every frame. Can use arbitrarily expensive methods (CFD, Monte Carlo, FEA).

---

## OVERVIEW

Input: Ship definition file (component list with shapes, materials, connections -- the same format as spatial.md's component graph).

Output: A binary table package containing all precomputed physics for that ship. The runtime simulator loads this package and interpolates the tables.

The characterization tool is a SEPARATE EXECUTABLE from the simulator. It's a Python orchestrator that calls specialized tools (Cantera, Geant4, OpenMC, FEniCS) and assembles their outputs into the table package.

---

## TABLE PACKAGE CONTENTS

### Per engine (combustion.md)
- Reaction equilibrium: (mixture_ratio, chamber_pressure) → (T_flame, products, gamma)
- Nozzle performance: (expansion_ratio, ambient_pressure) → (thrust_coeff, Isp, T_exit)
- Thermal loads: (chamber_pressure, mixture_ratio) → (heat_flux at wall, throat, nozzle)
- Tool: **Cantera** + isentropic expansion solver

### Per reactor (nuclear.md)
- Reactor response: (control_rods, burnup, coolant_temp) → (power, neutron_flux, gamma_flux, temperatures)
- Decay heat curve: (time_since_shutdown) → (decay_heat_fraction)
- NTP exhaust: (reactor_outlet_temp, propellant_properties) → (thrust, Isp)
- Tool: **OpenMC** for neutronics, custom point-kinetics solver for transient response

### Per BH drive configuration (hawking.md Tables 1-4; black_hole_paper authoritative)
- BH thermodynamics: (BH_mass) → (kT_H, species-summed f, P_H, evaporation rate, lifetime, channel power fractions)
- Propulsive performance: (BH_mass, shell_thickness) → (thrust, exhaust β from the Usov equilibrium solve, exhaust/system Isp, P_captured, stopping margin x_c, validity flag). NOTE: capture fraction is NOT an input — it is 0.83 by architecture (neutrino-limited); there is no reflector.
- Radiation environment: (BH_mass, shell_thickness) → (muon punch-through power/flux, neutrino power, EM leakage, exhaust plume spec)
- Mass feeding: (BH_mass, feed_rate) → (net_mass_rate, equilibrium FEED RATE, throttle time constant, pellet capture radius)
- Tool: Custom, extending black_hole_paper/calculations/ (hawking_spectrum.py, verify_design.py)

### Aerodynamic tables (ship.md regime bridging)
- Aero coefficients: (Mach, angle_of_attack) → (Cd, Cl, Cm, Cp_center)
- Heating: (velocity, altitude, AoA) → (stagnation_heat_flux, distributed_heat_flux)
- Tool: **CFD** (OpenFOAM or SU2) for subsonic/supersonic. Analytical correlations for hypersonic (Sutton-Graves, Fay-Riddell).
- Note: generated for each distinct aerodynamic configuration (Falcon 9 with fairing, without fairing, Dragon capsule alone, etc.)

### Structural model (structural.md)
- Reduced stiffness matrix: K_reduced (10-50 DOF)
- Stress-recovery matrices S_e per monitored joint/element (map retained-DOF displacements → local stress; REQUIRED — without them runtime failure checks have no input)
- Mode shapes and natural frequencies
- Failure envelopes: (load_direction) → (max_stress, failure_mode)
- Buckling limits per component
- Tool: **FEniCS** FEA

### Radiation shielding baseline (radiation_shielding.md)
- Shield thickness map: ray-traced from component graph octree
- Baseline dose transmission factors: (source_type, source_energy, shield_config) → (transmission)
- Tool: **Geant4** (single-interaction mode, element-level tables)

### Thermal conductances (thermal.md)
- Inter-component conductances for complex geometries
- View factors for radiation exchange
- Tool: **FEniCS** for complex geometries, analytical for simple ones

---

## TABLE FORMAT

All tables stored as binary files with a common format:

```
Header:
  magic: u32           // file format identifier
  version: u32         // format version
  n_dimensions: u32    // number of input dimensions
  dimensions: [(name: String, n_points: u32, min: f64, max: f64, spacing: Spacing)]
  n_outputs: u32       // number of output values per entry
  outputs: [(name: String, unit: String)]

Data:
  values: [f64; n_points_total × n_outputs]  // row-major, dimensions vary fastest-last
```

Interpolation: multilinear (linear in each dimension). For a 3D table, this is trilinear interpolation: 8 corner lookups, 7 lerps. Cost: ~0.3 μs per interpolation.

---

## DECOMPOSITION STRATEGY

To avoid exponential blowup in high-dimensional tables, decompose the physics into independent sub-problems:

**Example: Merlin 1D engine characterization**

NOT: one 5D table (throttle × altitude × mixture_ratio × inlet_temp × fuel_remaining)

Instead:
1. Reaction table: (mixture_ratio × chamber_pressure) → products, T_flame [2D]
   - Chamber pressure is a direct function of throttle (no separate dimension needed)
2. Nozzle table: (expansion_ratio × ambient_pressure) → Isp, thrust_coeff [2D]
   - Expansion ratio is fixed per engine (not a variable)
   - So this is effectively 1D: ambient_pressure → Isp, thrust_coeff
3. Thermal table: (chamber_pressure × mixture_ratio) → heat_flux [2D]

Three 2D tables instead of one 5D table. The interdependencies (e.g., chamber pressure depends on throttle and propellant supply pressure) are resolved by simple analytical relationships at runtime, not by expanding the table dimensionality.

**Principle:** if a relationship can be expressed as an analytical equation (ideal gas law, isentropic expansion, mass conservation), use the equation at runtime. Only tabulate the parts that require expensive computation (chemical equilibrium, nuclear reactions, CFD).

---

## PIPELINE

```
ship_definition.toml
    │
    ├──→ Parse component graph
    │
    ├──→ Build octree (for radiation ray-tracing)
    │
    ├──→ For each engine:
    │      Run Cantera equilibrium across parameter grid
    │      Compute nozzle performance (isentropic expansion)
    │      Compute thermal loads (Bartz correlation + Cantera properties)
    │      Write: engine_tables.bin
    │
    ├──→ For each reactor:
    │      Run OpenMC for neutronics across parameter grid
    │      Build point-kinetics response model
    │      Compute decay heat curve
    │      Write: reactor_tables.bin
    │
    ├──→ For BH drive (if present):
    │      Compute species-summed Hawking spectrum vs. BH mass (paper Eq. 5-6)
    │      Solve the Usov shell equilibrium per (mass, shell_thickness) → β, thrust
    │      Compute punch-through/stopping-margin and feeding tables
    │      Write: bh_drive_tables.bin
    │
    ├──→ Aerodynamics:
    │      Run CFD (or analytical correlations) for each aero configuration
    │      Write: aero_tables.bin
    │
    ├──→ Structural:
    │      Build FEniCS FEA model from component graph
    │      Compute reduced stiffness matrix
    │      Run eigenvalue analysis (mode shapes)
    │      Run parametric failure study
    │      Write: structural_tables.bin
    │
    ├──→ Thermal:
    │      Compute conductances for complex geometries (FEniCS)
    │      Compute view factors (Monte Carlo ray tracing)
    │      Write: thermal_tables.bin
    │
    ├──→ Radiation baseline:
    │      Ray-trace octree for shielding thickness map
    │      Generate Geant4 element interaction tables (if not already available)
    │      Write: radiation_tables.bin
    │
    └──→ Package all tables into ship_name.tables
         Total size: ~1-100 MB per ship (depends on number of configurations)
```

## TOOLS SUMMARY

| Tool | Used for | License |
|---|---|---|
| **Cantera** | Combustion chemistry, thermodynamic properties | BSD-3 |
| **CoolProp** | Propellant and fluid thermophysical properties | MIT |
| **OpenMC** | Neutron transport, reactor neutronics, depletion | MIT |
| **Geant4** | Particle transport, radiation interaction tables | Custom (free) |
| **FEniCS** | Structural FEA, thermal conductance computation | LGPL |
| **OpenFOAM / SU2** | CFD for aerodynamic tables | GPL-3 / LGPL |
| **galpy** | Epicyclic frequency tables (not per-ship, but per-universe) | BSD-3 |
| **Python** | Orchestration, data processing, table format writing | PSF |

## RUNTIME LOADING

The simulator loads the table package at startup:
```rust
let ship = Ship::load("falcon9.tables", "falcon9.toml");
// falcon9.toml = component graph definition
// falcon9.tables = precomputed lookup tables
```

Tables are memory-mapped for fast access. Total resident memory: ~1-100 MB per ship. Interpolation cost: ~0.3 μs per table lookup.
