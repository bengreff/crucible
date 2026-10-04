# Electrical Systems

Power generation, distribution, and storage. Bookkeeping, not physics -- evaluated once per trajectory step.

---

## STATE

```rust
struct ElectricalState {
    battery_charge: f64,         // joules
    battery_capacity: f64,       // joules (degrades with radiation dose)
    bus_voltage: f64,            // volts (simplified: one bus, one voltage)
    generator_states: Vec<GeneratorState>,
    load_states: Vec<LoadState>,
}
```

## GENERATORS

| Type | Power model | Degradation |
|---|---|---|
| Solar panels | P = η × A × S × cos(θ) × (1 - degradation). S = solar flux at distance. θ = angle to sun. η ~ 0.3. | Radiation dose: parameterized in 1 MeV ELECTRON-equivalent fluence (the industry standard — ~1% per 10^14 e/cm²; protons are ~100-1000× more damaging per particle and convert via NIEL ratios from the transport tool) |
| RTG | P = P_0 × 2^(-t/87.7yr). Pu-238 half-life decay. P_0 ~ 500 W per unit. | Age only (radioactive decay). |
| Reactor electrical | P = reactor_thermal × conversion_efficiency. η ~ 0.1-0.3 (thermoelectric or Brayton). | Fuel burnup reduces thermal power over mission. |
| Fuel cells | P = ΔG × n_dot. Gibbs energy × molar flow rate. Consumes H2 + O2, produces H2O + electricity. | Membrane degradation over time. |

## LOADS

Each sub-component (electronics, sensor, valve, pump, heater) has a power draw. Loads are categorized by priority (see resources.md load shedding list).

```rust
struct LoadState {
    component_id: SubComponentId,
    nominal_power: f64,     // watts
    is_active: bool,        // on/off
    priority: u8,           // 8 = critical (flight computer, never shed), 1 = optional
                            // (science, shed first) — matches resources.md's shed order
                            // (list position 1 sheds first)
}
```

## COMPUTATION

```
total_generation = Σ generator.current_power()
total_load = Σ active_load.nominal_power
net = total_generation - total_load

if net > 0:
    battery_charge = min(battery_charge + net × dt, battery_capacity)
else:
    battery_charge = battery_charge + net × dt  // net is negative, drains battery
    if battery_charge <= 0:
        battery_charge = 0
        trigger load_shedding(deficit = -net)
```

Cost: **~1 μs per tick.** Sum over ~10-50 generators and loads.

## BROWNOUT AND LOAD SHEDDING

When power demand exceeds supply + battery:
1. Sort loads by priority (lowest priority first)
2. Shed loads one by one until supply ≥ remaining demand
3. Shed loads are marked inactive
4. Each shed load has consequences:
   - Heaters off → temperatures drop (thermal module)
   - Pumps off → coolant flow stops → thermal consequences
   - Communications off → no contact with Earth/stations
   - Nav sensors off → reduced navigation capability
   - Flight computer off → ship is uncontrollable (catastrophic)

## TOOLS

No external tools. Pure arithmetic.
