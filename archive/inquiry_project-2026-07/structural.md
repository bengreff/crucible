# Structural Simulation

Reduced-order structural model. The full stiffness matrix is precomputed offline (FEniCS FEA). At runtime: multiply stiffness matrix × current load vector → stress at joints. Compare against material limits (which degrade with temperature and radiation damage). Detect and propagate failures.

---

## INPUTS (per tick)

| Load | Source | Typical magnitude |
|---|---|---|
| Thrust | Propulsion module | 0 - 10 MN (Falcon 9 first stage) |
| Aerodynamic | Aero tables (drag, lift distributed across surface) | 0 - 1 MN (max-Q) |
| Inertial (acceleration) | Trajectory integrator (mass × acceleration) | All mass × current a |
| Pressure differential | Atmosphere module (cabin) vs. external | ~100 kPa (1 atm cabin in vacuum) |
| Thermal stress | Thermal module (temperature gradients → differential expansion) | Depends on ΔT and material |
| Tidal gradient | GR module (if near compact object) | 0 normally, extreme near BH |
| Rotational | Attitude dynamics (centrifugal, Coriolis) | Depends on spin rate |
| Docking/collision impulse | Multi-object module | Event-driven |

## OUTPUTS (per tick)

| Output | Destination |
|---|---|
| Stress at each structural joint | Failure detection |
| Failure flags | Cascading failure model (ship.md Section 3) |
| Structural integrity fraction per component | Telemetry display |
| Vibration estimate (from dynamic loads) | Crew comfort, instrument degradation |

## ARCHITECTURE

### Offline: Stiffness matrix from FEA

For each ship configuration, the characterization tool builds a structural FEA model:
1. Components as structural elements (beams, shells, plates)
2. Connections as joints (bolted, welded, hinged) with defined stiffness
3. Compute the REDUCED stiffness matrix: for N components with 6 DOF each (3 translation + 3 rotation), the full system is 6N × 6N. Reduce to ~10-50 DOF by retaining only the critical load paths.
4. Store: reduced stiffness matrix K, and the mapping from loads to DOF.

### Runtime: Stress evaluation

Two-step pipeline (K⁻¹ × load gives DISPLACEMENTS, not stress):

```
u       = K_reduced⁻¹ × load_vector      // displacements at retained DOF
σ_e     = S_e × u                        // stress per element/joint, via
                                          // precomputed stress-recovery matrices
```

The stress-recovery matrices S_e (one small matrix per monitored joint/element, mapping
retained-DOF displacements to local stress components) are exported by the offline FEA
alongside K_reduced — they are part of the characterization deliverable, not optional.
Both steps are small matrix multiplies: ~10-50 DOF. Cost: **~5 μs.**

For each structural connection:
```
if stress > yield_strength(T_joint, dpa_joint):
    joint.failure_mode = Yielding  → permanent deformation
if stress > ultimate_strength(T_joint, dpa_joint):
    joint.failure_mode = Fracture  → connection breaks
```

yield_strength and ultimate_strength are functions of temperature and radiation damage from the material database (spatial.md).

### Failure modes

| Mode | Detection | Consequence |
|---|---|---|
| Yielding | stress > σ_yield | Permanent deformation. Performance degradation but no break. |
| Fracture | stress > σ_ultimate | Connection breaks. Structural separation. Cascading failure. |
| Buckling | compressive stress > Euler critical load | Column/shell collapse. Component destroyed. |
| Fatigue | cumulative cycles > fatigue life | Progressive weakening over many load cycles. |
| Creep | sustained high temperature + load | Slow deformation at high T. Relevant for engine components. |

Fatigue and creep are tracked as cumulative variables per joint, incremented each tick based on current stress and temperature.

### Thermal stress computation

When adjacent components are at different temperatures, differential thermal expansion creates stress:
```
σ_thermal = E × α × ΔT / (1 - ν)
```
where E is Young's modulus, α is thermal expansion coefficient, ΔT is the temperature difference, and ν is Poisson's ratio. These come from the material database, evaluated at current T and dpa.

## CHARACTERIZATION (offline)

| Computation | Tool | Output |
|---|---|---|
| Full structural FEA | FEniCS | Reduced stiffness matrix + stress-recovery matrices per ship configuration |
| Mode shapes (vibration) | FEniCS eigenvalue solve | Natural frequencies, mode shapes |
| Failure envelopes | FEniCS parametric study | Stress limits vs. load direction |
| Buckling limits | FEniCS linear buckling analysis | Critical loads per component |

## TOOLS

| Tool | When | Purpose |
|---|---|---|
| **FEniCS** | Offline | Full 3D FEA for stiffness matrices, failure envelopes, mode shapes |
| Custom Rust matrix multiply | Runtime | u = K⁻¹F, then σ_e = S_e·u; ~10-50 DOF, ~5 μs |
