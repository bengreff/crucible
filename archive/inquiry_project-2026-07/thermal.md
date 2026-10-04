# Thermal Simulation

Thermal network model. Each component is a thermal node with heat capacity. Connections define conduction, radiation, and convection paths. One of the three runtime ODEs in the ship simulation (with cabin atmosphere and reactor point kinetics -- SCOPE.md) -- solved implicitly (backward Euler) for unconditional stability at any timestep.

---

## INPUTS (per tick)

| Source | Heat input | From |
|---|---|---|
| Engine operation | Chamber wall heat flux, nozzle heating | combustion.md tables |
| Reactor operation | Core heat, decay heat | nuclear.md tables |
| BH drive | ≈ 0 nominal (CFL shell firewall — no radiators, hawking.md); any residual deposition arrives via the cached dose map row below, not a separate source | hawking.md |
| Solar flux | Direct solar illumination on external surfaces | Position, orientation, star luminosity |
| Albedo + IR | Reflected sunlight + thermal IR from nearby bodies | Position, body properties |
| Aerodynamic heating | Convective + radiative heating during atmospheric flight | ship.md Phase 1c (aero tables) |
| ISM heating | Particle energy deposition at relativistic speeds | radiation_transport.md (cached dose map) |
| Radiation heating | Energy deposition from all radiation sources | radiation_transport.md (cached dose map) |
| Internal equipment | Electronics waste heat, life support, pumps | electrical.md load data |
| Crew metabolic | ~100 W per crew member | Constant |

## OUTPUTS (per tick)

| Output | Destination |
|---|---|
| Temperature per component (K) | structural.md (thermal stress), material property evaluation everywhere |
| Radiator performance | Determines heat rejection capability |
| Component over-temperature warning | Event detection |
| Cabin temperature | Atmosphere module, crew comfort |

## THE THERMAL NETWORK

### Nodes

Each component in the component graph is one thermal node (or 2-3 for large components with internal temperature gradients, like an engine with chamber/throat/nozzle zones).

```rust
struct ThermalNode {
    component_id: ComponentId,
    temperature: f64,           // K
    heat_capacity: f64,         // J/K = mass × specific_heat(T)
    heat_generation: f64,       // W (from all internal sources this tick)
}
```

### Connections

```rust
struct ThermalConnection {
    node_a: NodeIndex,
    node_b: NodeIndex,
    conductance: f64,           // W/K = k × A / L (material conductivity × contact area / path length)
}
```

Additionally, every external-facing node radiates to space:
```
Q_rad = ε × σ × A_surface × (T⁴ - T_env⁴)
```
Where T_env is the background temperature (2.7 K in deep space, higher near a star or planet).

### Matrix formulation

The thermal network is:
```
C × dT/dt = A × T + q
```
Where:
- C = diagonal matrix of heat capacities (N×N)
- A = conductance matrix (sparse, N×N)
- T = temperature vector (N×1)
- q = heat source vector (N×1)

### Implicit solve (backward Euler)

```
C × (T_new - T_old) / dt = A × T_new + q
(C - dt×A) × T_new = C × T_old + dt × q
```

This is a sparse linear system: (C - dt×A) is a sparse N×N matrix. Solve with LU decomposition (for N ~ 50-200, direct solve is fine, ~5-10 μs).

**Unconditionally stable** at any timestep. No sub-cycling needed. Works at 1 ms steps (atmospheric flight) and 1000 s steps (interstellar cruise at timewarp).

### Radiation nonlinearity

The T⁴ radiation term makes the system nonlinear. Linearize around the current temperature:
```
T⁴ ≈ T_old⁴ + 4×T_old³ × (T - T_old)
```
This gives a linear system that's accurate for ΔT << T. For large temperature changes (engine startup, reentry heating), iterate: solve, update linearization point, solve again. Typically converges in 1-2 iterations.

## CHARACTERIZATION (offline)

### Conductance computation (FEniCS)

For complex component geometries (engine assemblies, multi-layer insulation), the effective thermal conductance between nodes is computed offline:
- Build a 3D FEA model of the component
- Apply unit temperature difference between nodes
- Compute steady-state heat flow → effective conductance

For simple geometries (cylindrical tanks, flat panels): analytical formulas suffice.
- Cylinder wall: G = 2πkL / ln(r_outer/r_inner)
- Flat wall: G = kA/L

### View factor computation

For radiation exchange between components (not just to space): compute view factors F_ij (fraction of radiation from surface i that reaches surface j). For simple shapes, analytical or tabulated. For complex assemblies, Monte Carlo ray tracing offline.

## TOOLS

| Tool | When | Purpose |
|---|---|---|
| **FEniCS** | Offline | Complex conductance computation, model validation |
| **scipy.sparse** | Offline (validation) | Sparse matrix solve reference |
| Custom Rust sparse solver | Runtime | LU decomposition for N ~ 50-200 matrix, ~5-10 μs |
