# Nuclear Systems

Lookup-table-based simulation of fission reactors, fusion reactors, and antimatter systems. All fast nuclear physics (neutron chain reactions, fusion burns, annihilation cascades) is precomputed by the characterization tool. Runtime is table interpolation.

---

## FISSION REACTOR (NTP and power generation)

### What the table encodes

The characterization tool runs a point-kinetics reactor model offline across the full parameter space and stores the steady-state and transient responses.

**Table inputs:**
- Control rod position (0-1): determines reactivity insertion
- Fuel burnup fraction (0-1): fresh fuel to end-of-life
- Coolant inlet temperature (K): affects temperature feedback

**Table outputs:**
- Thermal power (W)
- Neutron flux at shield boundary (n/cm²/s + energy spectrum)
- Gamma flux at shield boundary (γ/cm²/s + energy spectrum)
- Fuel temperature (K)
- Coolant outlet temperature (K)
- Decay heat fraction (fraction of rated power, for shutdown cooling)

**Table dimensions:** 3D (control_rod × burnup × coolant_temp), ~20³ = 8000 entries. ~64 KB per reactor type.

### Reactor types in the simulator

| Type | Application | Fuel | Power | Isp |
|---|---|---|---|---|
| NERVA-class NTP | Mars ship propulsion | U-235 / graphite composite | 1-5 GW thermal | ~900 s |
| Kilopower-class | Surface power, small ships | U-235 / steel | 1-10 kW electric | N/A |
| SP-100 class | Large spacecraft power | UN / NbC cermet | 100 kW electric | N/A |

### NTP-specific outputs

For nuclear thermal propulsion, the reactor heats hydrogen propellant directly:
- Thrust = m_dot × v_exhaust
- v_exhaust = √(2 × c_p × (T_outlet - T_inlet))
- Isp = v_exhaust / g0

These are computed from the reactor table outputs (coolant outlet temperature) plus the hydrogen thermophysical properties (from CoolProp tables, precomputed).

### Shutdown and decay heat

After reactor shutdown, fission products continue producing decay heat:
- Immediately after shutdown: ~6% of operating power
- After 1 hour: ~1.5%
- After 1 day: ~0.5%
- After 1 year: ~0.03% (Way-Wigner gives 0.025-0.05%)

The decay heat curve is a known function (Way-Wigner or ANS 5.1 standard). Stored as a 1D table: time_since_shutdown → decay_heat_fraction. The thermal module continues to apply this heat to the reactor component even after shutdown.

### Radiation output

The reactor produces neutrons and gammas that irradiate the ship. This is handled by the radiation transport tool:
- The reactor's neutron and gamma spectra (from the table) define a radiation source
- The transport tool computes dose through the shadow shield
- Results are cached in the dose-rate map

### Offline characterization tools

| Tool | Purpose |
|---|---|
| **OpenMC** | Monte Carlo neutronics: compute neutron flux, reaction rates, depletion |
| **Cantera** (thermodynamics) | Hydrogen thermophysical properties for NTP exhaust |
| Custom point-kinetics solver | Generate the reactor response table across the parameter space |

---

## FUSION SYSTEMS (future ships)

### Fusion types considered

| Type | Fuel | Ignition method | Status in simulator |
|---|---|---|---|
| D-T ICF pellets (Daedalus-class) | Deuterium + Tritium | Electron beams or laser | Future ship |
| D-He3 ICF pellets | Deuterium + Helium-3 | Electron beams | Future ship |
| ACMF (antimatter-catalyzed) | D-T + antiprotons | Antiproton annihilation | Research paper #3 ship |

### Table structure for fusion engines

**Table inputs:**
- Pellet mass (kg)
- Pellet composition (D:T or D:He3 ratio)
- Ignition energy (J) or antiproton count (for ACMF)
- Repetition rate (Hz)

**Table outputs:**
- Thrust per pellet (N·s impulse)
- Exhaust velocity (m/s)
- Energy release per pellet (J)
- Neutron flux (for D-T: 80% of energy in 14.1 MeV neutrons)
- Radiation spectrum at containment boundary
- Waste heat to structure

**Table dimensions:** 3D for ACMF (pellet_mass × antiproton_count × compression), ~20³ = 8000 entries.

### Offline characterization tools

| Tool | Purpose |
|---|---|
| **Geant4** | Antiproton annihilation cascade in compressed fuel (ACMF ignition research) |
| Custom hydro code or published ICF tables | Pellet burn physics, energy partition |

---

## ANTIMATTER SYSTEMS

### Antimatter storage

Antihydrogen stored in Penning traps (electromagnetic confinement). Storage is a sub-component with:
- Antimatter mass (kg)
- Trap power consumption (W) -- antimatter is lost if power fails (annihilation with trap walls)
- Loss rate (kg/s) from trap imperfections

Power failure → trap collapse → antimatter annihilation → explosion (energy = 2mc²). This is a catastrophic failure mode. The energy release from even 1 microgram of antimatter: 2 × 10^-9 × (3×10^8)² = 180 MJ (equivalent to ~43 kg of TNT). Component destruction + radiation burst.

### Antimatter as propulsion

For ACMF: antiprotons are injected into fusion pellets at a controlled rate. The antiproton feed rate is a control input. The engine table maps (feed_rate, pellet_properties) → (thrust, radiation).

For pure antimatter propulsion (future, very speculative): proton-antiproton annihilation produces ~5 pions. Charged pions are magnetically directed for thrust. Neutral pions decay to gammas (lost). Efficiency ~50-60% of annihilation energy to directed thrust. Table: (annihilation_rate) → (thrust, gamma_flux, pion_flux).

---

## RUNTIME EVALUATION

Steady-state maps, flux spectra, and decay heat are table lookups. Reactor TRANSIENTS run a
point-kinetics ODE at runtime — a static steady-state table cannot produce rod-insertion
transients, startup dynamics, or feedback-driven power excursions, and these are exactly the
events the player interacts with (physics.md Section 1.2 specifies this; the "tables-only"
formulation was inconsistent with it):

```rust
fn evaluate_reactor(tables: &ReactorTable, state: &mut ReactorState, dt: f64) -> ReactorOutput {
    if state.is_shutdown {
        let decay_heat = tables.decay_heat_curve.interpolate(state.time_since_shutdown);
        return ReactorOutput { thermal_power: state.last_power * decay_heat, ..Default::default() };
    }
    // Point kinetics: 7 stiff ODEs (power + 6 delayed-neutron groups), implicit step.
    // Reactivity = rod_worth(rods, burnup) [table] + doppler(T_fuel) + coolant(T_cool) [tables]
    let rho = tables.rod_worth.interpolate(state.control_rods, state.burnup)
            + tables.feedback.interpolate(state.fuel_temp, state.coolant_temp);
    state.kinetics.step_implicit(rho, dt);          // deterministic fixed algorithm
    // Steady-state map supplies flux spectra and temperature distribution at current power:
    tables.response_at_power(state.kinetics.power, state.burnup, state.coolant_temp)
}
```

The kinetics step is implicit (stiff system; unconditionally stable at any dt, matching the
thermal network's timewarp behavior). Reactivity coefficients and rod worths come from the
OpenMC characterization tables — the expensive neutronics stays offline; only the 7-variable
transient integrator runs live.

Cost: **~2 μs per tick.** One 7×7 implicit solve + table interpolations.
