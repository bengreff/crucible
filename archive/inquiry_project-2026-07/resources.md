# Resource Flow Simulation

Bookkeeping for all consumables and energy. Not an ODE -- just conservation-law accounting evaluated once per trajectory step. Tracks mass, volume, and energy through the component graph's resource connections.

---

## TRACKED RESOURCES

| Resource | Unit | Where stored | Consumed by | Produced by |
|---|---|---|---|---|
| LOX | kg | Oxidizer tanks | Engines | None (loaded at launch) |
| RP-1 / LH2 / CH4 | kg | Fuel tanks | Engines | None |
| Helium (pressurant) | kg | COPV bottles | Tank pressurization | None |
| Reactor fuel (U-235) | kg | Reactor core | Fission | None |
| BH mass | kg | BH containment | Hawking radiation | Mass feeding |
| Coolant | kg | Coolant loops | Evaporative losses | None |
| Breathable O2 | kg | Air tanks / cabin | Crew metabolism | Electrolysis |
| N2 (cabin buffer gas) | kg | N2 tanks / cabin | Leaks | None |
| CO2 | kg | Cabin (undesired) | None | Crew metabolism |
| Water | kg | Water tanks | Crew consumption, electrolysis | Fuel cells, CO2 scrubbing |
| Food | kg | Storage | Crew consumption | None |
| Electrical energy | J | Batteries | All systems | Solar, reactor, RTG, fuel cells |

## INPUTS (per tick)

- Engine thrust level → fuel + oxidizer consumption rate (from Isp)
- Reactor power → fuel burnup rate (from nuclear.md tables)
- BH drive state → BH mass change rate (from hawking.md tables)
- Crew count → O2 consumption, CO2 production, water consumption, food consumption, metabolic heat
- System loads → electrical power consumption
- Leak rate → gas loss (from structural damage / atmosphere module)

## OUTPUTS (per tick)

| Output | Destination |
|---|---|
| Fuel remaining (per tank) | Engine tables (affects available thrust duration) |
| Tank fill fraction | Octree (propellant density in tank voxels), MOI computation |
| Power available / deficit | Electrical module |
| O2 / CO2 / H2O rates | Atmosphere module |
| Mass change | Trajectory (rest mass update), MOI recomputation flag |
| Resource depletion warnings | Event detection (fuel low, O2 low, battery low) |

## COMPUTATION

### Mass flow

For each resource connection (pipe/line) in the component graph:
```
flow_rate = min(demanded_rate, valve_max_rate × valve_position, available_in_source)
source.amount -= flow_rate × dt
sink.amount += flow_rate × dt
```

Conservation: Σ (source depletion) = Σ (sink accumulation) for each resource type. Verified every tick. Violation = bug.

### Propellant consumption

```
mass_flow_total = thrust / (Isp × g0)    // from combustion.md
fuel_flow = mass_flow_total / (1 + mixture_ratio)
oxidizer_flow = mass_flow_total × mixture_ratio / (1 + mixture_ratio)
```

### Power balance

```
generation = solar_power(distance, orientation, panel_area, degradation)
           + reactor_electrical(reactor_power, conversion_efficiency)
           + rtg_power(initial_power, age)   // Pu-238 decay: P(t) = P0 × 2^(-t/87.7yr)
           + fuel_cell_power(H2_flow, O2_flow)

load = life_support + computers + communications + heaters + pumps + actuators + instruments

battery_delta = (generation - load) × dt
battery_charge = clamp(battery_charge + battery_delta, 0, battery_capacity)

if generation < load and battery_charge == 0:
    trigger brownout → shed loads by priority
```

### Load shedding priority (lowest priority shed first)

1. Science instruments
2. Non-essential heating
3. Communications (non-emergency)
4. Lighting
5. Environmental control (partial)
6. Navigation sensors
7. Flight computers (NEVER shed -- ship is uncontrollable without)
8. Life support critical (NEVER shed)

## TOOLS

No external tools needed. Pure arithmetic bookkeeping. ~1 μs per tick.
