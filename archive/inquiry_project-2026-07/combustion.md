# Chemical Propulsion

Lookup-table-based simulation of chemical rocket engines. The characterization tool (Cantera offline) computes combustion products, performance, and thermal loads across the parameter space. The runtime simulator interpolates these tables.

---

## INPUTS (per tick, from ship state)

- Throttle command (0-1)
- Ambient pressure (Pa, from atmospheric model or 0 in vacuum)
- Propellant inlet conditions (temperature from thermal module, pressure from tank state)
- Fuel remaining (affects mixture ratio if tanks drain unevenly)

## OUTPUTS (per tick, to ship systems)

| Output | Destination | Use |
|---|---|---|
| Thrust (N) | Trajectory integrator | Force on ship |
| Specific impulse (s) | Resources | Compute fuel consumption rate |
| Chamber temperature (K) | Thermal module | Heat source at engine location |
| Nozzle exit temperature (K) | Thermal module | Heat at nozzle components |
| Exhaust composition | Environment (for plume impingement if modeled) | -- |
| Engine status (nominal / off-nominal) | Event detection | Anomaly warning |

## CHARACTERIZATION (offline, Cantera)

### Table 1: Reaction equilibrium
**Tool:** Cantera equilibrium solver
**Inputs:** fuel type, oxidizer type, mixture ratio, chamber pressure
**Outputs:** flame temperature, product composition, gamma (ratio of specific heats)
**Dimensions:** 2D (mixture_ratio × chamber_pressure), ~20×20 = 400 entries per propellant combination

### Table 2: Nozzle performance
**Tool:** Isentropic expansion from Cantera equilibrium state
**Inputs:** chamber conditions (from Table 1), nozzle expansion ratio, ambient pressure
**Outputs:** thrust coefficient, specific impulse, exit temperature, exit pressure
**Dimensions:** 2D (expansion_ratio × ambient_pressure), ~10×50 = 500 entries

### Table 3: Thermal loads
**Tool:** Cantera + heat transfer correlations (Bartz equation for convective heat transfer)
**Inputs:** chamber pressure, mixture ratio, cooling channel geometry
**Outputs:** heat flux at chamber wall, nozzle throat, nozzle exit
**Dimensions:** 2D (chamber_pressure × mixture_ratio), ~20×20

### Propellant combinations

| Combo | Fuel | Oxidizer | Application |
|---|---|---|---|
| LOX/RP-1 | Kerosene | Liquid oxygen | Falcon 9 first stage (Merlin 1D) |
| LOX/RP-1 (vacuum) | Kerosene | Liquid oxygen | Falcon 9 second stage (Merlin Vac) |
| LOX/LH2 | Liquid hydrogen | Liquid oxygen | Upper stages, future ships |
| LOX/CH4 | Methane | Liquid oxygen | Raptor-class engines (Starship) |
| N2O4/UDMH | Unsymmetrical dimethylhydrazine | Nitrogen tetroxide | Hypergolic (SuperDraco, orbital maneuvering) |

Each combination gets its own set of tables. Total table size: ~5 combinations × ~1500 entries × ~32 bytes = **~240 KB.** Trivial.

## RUNTIME EVALUATION

```
fn evaluate_engine(tables: &EngineTable, throttle: f64, ambient_pressure: f64, 
                   fuel_state: &FuelState) -> EngineOutput {
    let chamber_pressure = throttle * tables.max_chamber_pressure;
    let mixture_ratio = fuel_state.current_mixture_ratio();
    
    let equilibrium = tables.reaction.interpolate(mixture_ratio, chamber_pressure);
    let nozzle = tables.nozzle.interpolate(tables.expansion_ratio, ambient_pressure);
    
    EngineOutput {
        thrust: nozzle.thrust_coefficient * chamber_pressure * tables.throat_area,
        isp: nozzle.specific_impulse,
        chamber_temp: equilibrium.flame_temperature,
        heat_flux: tables.thermal.interpolate(chamber_pressure, mixture_ratio),
    }
}
```

Cost: ~3 table interpolations × ~0.3 μs each = **~1 μs per engine per tick.**

## ENGINE-SPECIFIC PARAMETERS (not in tables)

These are fixed per engine model, not interpolated:
- Throat area (determines thrust at a given chamber pressure)
- Expansion ratio (nozzle geometry, fixed)
- Number of engines (Falcon 9: 9 first stage, 1 second stage)
- Gimbal range (±5° typically, affects thrust direction)
- Throttle range (Merlin: 40-100%. Some engines cannot throttle.)
- Ignition sequence timing
- Shutdown transient (thrust decay over ~1 second)

## TOOLS INTEGRATION

| Tool | When | Purpose |
|---|---|---|
| **Cantera** (Python) | Offline (characterization) | Equilibrium chemistry, product composition, flame temperature |
| **CoolProp** (Python) | Offline | Propellant thermophysical properties (density, viscosity, vapor pressure vs. T and P) |
| None at runtime | -- | Pure table interpolation in Rust |
