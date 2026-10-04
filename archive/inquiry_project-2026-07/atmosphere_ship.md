# Crew Cabin Atmosphere

Gas composition tracking inside pressurized volumes. Small ODE system (~10 state variables per volume). Evaluated once per trajectory step.

---

## STATE (per pressurized volume)

```rust
struct AtmosphereState {
    n_o2: f64,       // moles of O2
    n_n2: f64,       // moles of N2
    n_co2: f64,      // moles of CO2
    n_h2o: f64,      // moles of H2O vapor
    n_ar: f64,       // moles of Ar (trace, mostly inert)
    n_co: f64,       // moles of CO (toxic contaminant)
    n_nh3: f64,      // moles of NH3 (toxic contaminant)
    volume: f64,     // m³ (from octree enclosed volume detection)
    leak_area: f64,  // m² (from structural damage, 0 = sealed)
}
```

Derived:
- Total moles: N = Σ n_i
- Pressure: P = NRT/V (ideal gas, R = 8.314 J/(mol·K), T from thermal module)
- Partial pressures: P_i = n_i RT/V
- Mole fractions: x_i = n_i / N

## INPUTS (per tick)

| Source | Effect |
|---|---|
| Crew metabolism | O2 consumed, CO2 + H2O produced. Rate: ~1.09 mol O2/hr/person (NASA standard 0.84 kg/day — the old "0.84 mol/hr" confused kg with mol, ~23% low), ~0.95 mol CO2/hr/person (respiratory quotient ~0.87), ~0.6 mol H2O/hr/person |
| CO2 scrubber | CO2 removal rate. LiOH: irreversible, limited supply. CDRA (molecular sieve): regenerable, requires power. Rate: up to ~1 mol CO2/hr per unit. |
| O2 generation | Electrolysis of water: 2H2O → 2H2 + O2. Rate depends on power available. |
| Cabin temperature | From thermal module. Affects pressure (ideal gas law). |
| Leak rate | From structural module. If hull breached to vacuum, flow is choked: ṁ = C_d × A_hole × P₀ × √(γ_gas M_mol/(R_u T₀)) × (2/(γ_gas+1))^((γ_gas+1)/(2(γ_gas−1))). Reference value: ~0.024 kg/s per cm² for cabin air at 101.3 kPa, 293 K, C_d = 1. For subsonic leaks (downstream/upstream pressure ratio > 0.53): ṁ = C_d × A × √(2ρΔP). |
| Outgassing | Slow release of volatiles from materials. Small, constant rate per material type. |
| Equipment | Some equipment produces trace contaminants (CO from partial combustion, NH3 from waste). |

## OUTPUTS (per tick)

| Output | Destination |
|---|---|
| Pressure (Pa) | Structural (pressure loads on walls), telemetry |
| O2 partial pressure | Crew safety (alarm if < 16 kPa or > 23 kPa) |
| CO2 partial pressure | Crew safety (alarm if > 0.5 kPa, impairment > 1 kPa, lethal > 7 kPa) |
| Humidity (H2O partial pressure) | Crew comfort, condensation risk |
| Toxic contaminants (CO, NH3) | Crew safety (alarm at OSHA limits) |
| Leak mass flow rate | Resources (gas loss rate) |

## COMPUTATION

Simple explicit Euler (stable because time constants are minutes-hours, dt is << 1 second):

```
dn_o2/dt = electrolysis_rate - crew_consumption_rate - leak_rate × x_o2
dn_co2/dt = crew_production_rate - scrubber_rate - leak_rate × x_co2
dn_h2o/dt = crew_production_rate - condenser_rate - electrolysis_input - leak_rate × x_h2o
dn_n2/dt = -leak_rate × x_n2
```

Cost: ~10 multiplications + additions per gas species. **~2 μs per tick.**

## ALARM THRESHOLDS

| Condition | Threshold | Consequence |
|---|---|---|
| O2 low | < 16 kPa (~16%) | Hypoxia warning |
| O2 critical | < 12 kPa (~12%) | Incapacitation in minutes |
| CO2 elevated | > 0.5 kPa (~0.5%) | Headache, reduced performance |
| CO2 dangerous | > 4 kPa (~4%) | Loss of consciousness in minutes |
| Pressure low | < 50 kPa | Decompression warning |
| Pressure critical | < 20 kPa | Ebullism danger |
| CO elevated | > 50 ppm | Toxic exposure |
| Cabin temperature | > 35°C or < 10°C | Crew thermal stress |

## MULTIPLE PRESSURIZED VOLUMES

If the ship has multiple sealed compartments (crew cabin, airlock, laboratory), each has its own atmosphere state. Connections between compartments (hatches) allow gas exchange when open:

```
flow_rate = C_d × A_hatch × √(2 × avg_density × |P_a - P_b|)   // kg/s, subsonic orifice
// (choked-flow formula above applies if the pressure ratio exceeds 1.89 for air, γ=1.4)
```

A breached compartment can be sealed off (close hatch) to protect other compartments. This is a cascading failure mitigation action.

## TOOLS

No external tools needed at runtime. Cantera can be used offline to validate the scrubber chemistry (LiOH + CO2 → Li2CO3 + H2O) if more detailed chemistry is desired.
