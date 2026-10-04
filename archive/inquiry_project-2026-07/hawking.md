# Black Hole Drive

Lookup-table-based simulation of the CFL-shell black hole drive from
[black_hole_paper/](black_hole_paper/) (Greff 2026, draft 6). The paper is the authoritative
source for all drive physics; this document specifies how that physics becomes runtime tables.
The characterization tool extends the paper's calculation scripts
(`black_hole_paper/calculations/`). Runtime is table interpolation.

---

## PHYSICS SUMMARY (from the paper)

A micro black hole of mass M radiates the full Standard Model spectrum at

    kT_H = ℏc³/(8πGMk_B)          (10.57 GeV at M = 10⁹ kg)
    P_H  = ℏc⁶ f(M) / (G²M²)      (60,200 TW at M = 10⁹ kg)

where f(M) is the species-summed emission factor (Page greybody coefficients × DOF ×
finite-mass suppression S(m/kT); f = 3.503×10⁻³ at 10.6 GeV, rising to ~4.1×10⁻³ above the
top/W/Z/Higgs thresholds). Unfed lifetime τ = G²M³/(3ℏc⁴f) = 15.8 yr at 10⁹ kg. Equilibrium
fuel rate dM/dt = P_H/c² = 670 g/s at 10⁹ kg. Fuel is **any matter** (mc² per kg regardless
of composition).

**The engine is not a photon rocket and there is no reflector.** A color-flavor-locked SQM
shell (R = 2 m) absorbs ~83% of P_H (irreducible losses: 7% primary neutrinos, ~7.7% secondary
neutrinos from stopped-hadron decay chains — budgeted 10% — and 0.1% gravitons) and re-emits
the absorbed power as keV e± pairs via the Usov electrosphere mechanism at a self-regulating
equilibrium T_eq = 0.356 GK. An anti-Helmholtz magnetic bottle (227 T/m axial gradient,
227–390 T at the shell, 1000 T at the throat) confines the magnetically charged BH
(g = 1.182×10⁶ A·m) and funnels the pair exhaust through a 0.3 mm aft bore into an 820 m
tapered magnetic nozzle (divergence 0.010°, loss ~10⁻⁸).

    Exhaust velocity:  β = 0.332 (Usov equilibrium; ±3% over 100× uncertainty in J(ζ))
    Thrust:            F = β · P_captured / c = β · 0.83 · P_H / c   (55.3 MN at 10⁹ kg)
    Exhaust Isp:       1.01×10⁷ s
    System Isp:        8.4×10⁶ s   (v_eff = 0.83·β·c = 0.276c, after 17% neutrino loss)

**Thermal management is architectural, not a ship subsystem.** The CFL gap blocks conduction
into the shell bulk (suppression ~10⁻²⁸²⁸); the three-layer outer firewall limits off-axis EM
leakage to < 1 kW. The drive needs **no radiators and no cooling loop**. Nominal drive heat
load on the rest of the ship ≈ 0. Non-exhaust emissions: neutrinos (~17% of P_H, free-streaming,
no deposition), muon punch-through (84 MW at 3000 fm shell; free-streaming through structure
but a severe biological dose — see radiation tables), EM leakage < 1 kW.

**There is no off switch.** Hawking emission cannot be stopped; the shell always converts;
the exhaust always flows. Thrust is set by M alone (F ∝ f(M)/M²) and is adjusted only by
letting M drift via the feed balance dM/dt = feed − P_H/c². Throttle time constants are
months at 10⁹ kg (Δ M of 1% ≈ 173 days at max differential feed) and hours–days at sprint
masses. Drive termination options: (a) stop feeding and let the BH evaporate through
end-of-life (Section: failure modes), or (b) structurally eject the drive assembly.

### Design space

| M (kg) | kT (GeV) | P_H (W) | F (MN) | Shell t (fm) | M_shell (kg) | Fuel (kg/s) | τ_unfed |
|---|---|---|---|---|---|---|---|
| 10⁹ (reference) | 10.6 | 6.02×10¹⁶ | 55.3 | 3000 | 6.0×10⁷ | 0.67 | 15.8 yr |
| 2×10⁸ | 53 | 1.5×10¹⁸ | 1,375 | 1.5×10⁴ | 3×10⁸ | 17 | ~46 d |
| 5×10⁷ (sprint) | 212 | 2.4×10¹⁹ | 21,800 | 6×10⁴ | 1.2×10⁹ | 264 | ~17 h |

Shell thickness scales ∝ kT ∝ 1/M to hold the stopping budget x_c = E_c/kT fixed
(E_c = 80 MeV/fm × t). Crewed configurations thicken the shell (5000 fm at reference) to cut
muon punch-through from 84 MW → 58 W (crew dose at 100 m: 4×10⁴ Sv/yr → 0.4 mSv/yr).

---

## TABLE STRUCTURE

### Table 1: BH thermodynamics
**Inputs:** BH mass (kg)
**Outputs:** kT_H (GeV), emission factor f (species-summed, with all SM thresholds),
P_H (W), evaporation rate P_H/c² (kg/s), unfed lifetime (s), species power fractions
(charged / neutral-hadron / e± / μ+τ channel / neutrino / graviton)
**Dimensions:** 1D (mass), ~200 log-spaced points from 10⁵ to 10¹² kg. The wide range covers
end-of-life runaway and bootstrap studies; the drive's *operational* validity is flagged by
Table 2.
**Source:** `hawking_spectrum.py` / `verify_design.py` species sum (Eq. 5–6 of the paper).

### Table 2: Propulsive performance
**Inputs:** BH mass (kg), shell thickness (fm)
**Outputs:** thrust (N), exhaust β (from the per-point Usov equilibrium solve), exhaust Isp (s),
system Isp (s), P_captured (W), stopping-budget margin x_c = E_c/kT (dimensionless),
validity flag (x_c > x_min)
**Dimensions:** 2D (~200 mass × 10 shell-thickness points).
**Note:** capture fraction is NOT a free parameter — it is 0.83 (neutrino-limited) by
architecture. The magnetic-bottle first-pass fraction only sets the shell's internal Usov
equilibrium (energy recycles); it does not appear in net thrust. Verified value: ~8.1%, not
the paper's 12.1% — the paper's 227 T equatorial field is a linear-gradient extrapolation
invalid at r = R_c; exact Biot-Savart gives 155.5 T (paper erratum, AUDIT.md). The
characterization solve uses the exact field map; β and thrust are logarithmically robust to
this (±3% over 100× J), so headline performance stands.

### Table 3: Radiation environment
**Inputs:** BH mass (kg), shell thickness (fm)
**Outputs:** muon punch-through power P_punch = 0.093·P_H·ε(x_c) (W) and muon flux at 1 m
(scale as 1/r²), neutrino power (W, non-interacting, listed for bookkeeping), off-axis EM
leakage (< 1 kW, constant), exhaust plume power and β (for the nozzle-axis exclusion zone)
**Dimensions:** 2D (mass × shell thickness).
**Note:** feeds the radiation transport tool as an onboard source. Punch-through muons
free-stream (verified exit energy ~11.5 GeV, γ ≈ 110, decay length ~72 km — still far beyond
ship scale): negligible energy deposition in structure, dominant crew dose term. The dose is
NUMBER-flux driven (MIPs), so the paper's flux/dose pair stands; the "84 MW" figure is the
above-budget incident CHANNEL power (P_punch = 0.093·P_H·ε(x_c)) while actual escaping power
is ~4 MW — the tables carry both. ε(x_c) = e^(−x_c)(x_c² + 4x_c + 6)/(7π⁴/120).

### Table 4: Mass feeding
**Inputs:** BH mass (kg), feed rate (kg/s)
**Outputs:** net dM/dt (feed − P_H/c²), equilibrium feed rate, throttle time constant
(seconds per 1% mass change at max feed differential), pellet capture radius vs injection
velocity (magnetic capture by the monopole field; ~11 cm at β = 0.33, larger at lower speed)
**Dimensions:** 2D (~200 mass × 50 feed-rate points).
**Feed ceiling:** set by bore aperture and railgun cadence (ship definition parameter).

---

## RUNTIME EVALUATION

```rust
fn evaluate_bh_drive(tables: &BHDriveTables, state: &BHDriveState) -> BHDriveOutput {
    let thermo = tables.thermodynamics.interpolate(state.bh_mass);
    let prop   = tables.performance.interpolate(state.bh_mass, state.shell_thickness);
    let rad    = tables.radiation.interpolate(state.bh_mass, state.shell_thickness);
    let feed   = tables.feeding.interpolate(state.bh_mass, state.feed_rate);

    BHDriveOutput {
        thrust: prop.thrust,                    // along coil axis; no gimbal — steer the ship
        system_isp: prop.system_isp,
        bh_mass_rate: feed.net_mass_rate,       // feed − evaporation; thrust throttles via M
        radiation_source: RadiationSource::bh_drive(rad), // muons + EM leak + plume
        thermal_load: 0.0,                      // architectural: CFL firewall, no radiators
        stopping_margin: prop.stopping_margin,  // x_c; event detection watches this
        confinement_ok: state.coil_power_ok && state.pd_control_ok,
    }
}
```

Cost: **~1 μs per tick.** Four table interpolations.

There is no `drive_on` flag and no `capture_fraction` control. Control inputs are:
- **feed_rate** (kg/s): player/autopilot; any-matter pellets via railgun through the forward port
- **shell configuration**: fixed per ship definition (thickness chosen crewed/uncrewed)
- **coil + PD controller power**: electrical loads; losing them is a failure mode, not a throttle

## MASS MANEUVERING (thrust is a free dial with slow dynamics)

The BH mass — and therefore thrust, F ∝ f(M)/M² — is fully commandable in operation. One
ship covers the whole design space (no fixed variants needed). But "at will" has physics:

**Thrust-UP (shrink the BH):** stop/reduce feeding; the BH evaporates toward the target mass.
Time = τ(M_start) − τ(M_target) with τ ∝ M³/f — dominated by the high-mass phase:

| Maneuver | Duration |
|---|---|
| 10⁹ → 5×10⁸ kg (F: 55 → 220 MN) | ~14 yr |
| 10⁹ → 2×10⁸ kg (F: 55 → 1,375 MN) | ~15.7 yr |
| 2×10⁸ → 10⁸ kg (F: 1,375 → 5,500 MN) | ~40 days |
| 10⁸ → 5×10⁷ kg (F → 21,800 MN) | ~5 days |

At the reference mass the drive is effectively fixed-thrust on mission timescales; agility
lives at the low-mass end. Mission planning therefore chooses the CRUISE mass band up front.

**Thrust-DOWN (grow the BH):** feed above equilibrium; time = ΔM / (feed_max − P_H/c²),
limited by the railgun ceiling (ship definition parameter). Growing 10⁸ → 10⁹ kg requires
injecting ≥ 9×10⁸ kg of matter — the bulk store must carry it.

**The shell is the envelope boundary (fixed hardware):** shell thickness must be sized for
the MINIMUM planned operating mass, since kT ∝ 1/M and the stopping budget is
E_c = 80 MeV/fm × t. Operating below M_min collapses the stopping margin x_c = E_c/kT and
punch-through grows exponentially (failure table below). Ship 3's shell is therefore chosen
for its intended mass band: sized at M_min with crewed thickening (paper Table 4/5 scaling,
M_shell ∝ 1/M_min), and the drive refuses no command — it just reports the margin.

## FAILURE MODES (event detection)

| Mode | Trigger | Physics | Severity |
|---|---|---|---|
| Stopping-budget breach | x_c = 22.7·(M/10⁹ kg) at 3000 fm: punch-through hits 1% of P_H at M ≈ 2×10⁸ kg and grows exponentially below (kT equals the full 240 GeV budget already at M = 4.4×10⁷ kg — the paper's §6.2 "M < 10⁷ kg" understates the threshold ~20×; see errata) | Punch-through power grows as ε(x_c); thrust efficiency degrades; shell not structurally harmed | Graceful degradation; controlled shutdown = stop feeding |
| Confinement loss | Coil quench, PD controller power loss, or actuator failure | Lateral instability e-folds in τ ≈ 2.8 s; BH drifts off-center, exits shell; ship is then beside an unshielded isotropic multi-10¹⁶ W GeV source | Catastrophic within seconds; the defining emergency of this ship |
| Evaporation end-of-life | Unfed for ~τ(M) | M shrinks, kT rises, punch-through failure occurs long before final burst; final evaporation (~10⁶ kg → seconds) only reached if deliberately run down | Catastrophic only if ignored through the degradation phase |
| Feed starvation | Fuel exhausted | Same as unfed evaporation; τ gives the clock (15.8 yr at reference — ample warning) | Slow |

The per-tick event detector (ship.md Phase 4) watches `stopping_margin`, `confinement_ok`,
and fuel reserves against τ(M).

## SHIP INTEGRATION (component graph)

The drive assembly is real geometry, not a point:

- **SQM shell**: sphere component, R = 2 m, material `sqm-cfl`, mass 6×10⁷–1.2×10⁹ kg
  (mass scales with design point — it dominates the ship's MOI). Radiation-opaque to
  everything except muons above E_c and neutrinos.
- **Anti-Helmholtz coils**: two torus components at z = ±2 m, superconducting,
  I = 1.36×10⁹ A, hoop stress 6.7 GPa. Stored energy: the paper's "~50 GJ" is low —
  the field map alone (≥227 T/m over the r ≤ 2 m cavity) integrates to ~800 GJ, with the
  near-conductor region pushing toward TJ (paper erratum; recompute at characterization).
  Electrical load + quench failure mode. Stored-energy release on quench is an explosion
  event (ship.md §3) — at ~TJ scale, quench is hull-loss, not a component failure.
- **Magnetic nozzle**: 820 m truss-mounted taper aft. Defines the exhaust exclusion axis.
- **Railgun feed + forward port**: SQM-capped during cruise; feed system draws electrical
  power per pellet.
- **PD confinement controller**: electronics sub-component, critical=true. Magnetometer
  position readout: the BH's monopole field at 2 m is 0.030 T (μ₀g/4πr², g = 1.182×10⁶ A·m),
  NOT the paper's "30 T" (factor-10³ erratum); at 1 nT noise the readout precision is
  ~30 nm — still 6 orders below the 100 mm equilibrium offset, so control closes fine, but
  the spec number and the ~11 cm pellet-capture radius derive from the field and need
  re-derivation at characterization. Coil-current modulation actuator (~1%).
- **The BH itself**: sub-voxel point source at the bottle center; its gravity on ship
  components (6.7×10⁻⁶ m/s² at 100 m for 10⁹ kg) is a static structural load, negligible
  for trajectory (the drive is not self-gravitating in any meaningful sense).
- **Crew standoff**: crewed configs put the cabin ≥ 100 m forward with the 5000 fm shell
  (0.4 mSv/yr from punch-through — below deep-space background).

## CHARACTERIZATION (offline)

| Computation | Tool | Notes |
|---|---|---|
| f(M) species sum with thresholds | `hawking_spectrum.py` (extend) | Paper Eq. 5–6; optionally cross-check with BlackHawk+PYTHIA (paper §7.5) |
| Usov equilibrium per (M, t_shell) | `verify_design.py` solver (extract) | T_eq, β, J(ζ) self-consistency |
| Punch-through / dose tables | `sputter_bound.py`, `secondary_neutrinos.py` + Geant4 muon transport | Dose at crew stations via transport tool |
| Feed capture radii | `mission_profiles.py` + bore geometry | vs pellet velocity |
| Validation | `verify_design.py` end-to-end | Table below |

## VALIDATION TARGETS (from the paper)

| Quantity at M = 10⁹ kg, t = 3000 fm | Value |
|---|---|
| kT_H | 10.57 GeV |
| P_H | 60,200 TW |
| Thrust | 55.3 MN |
| Equilibrium feed | 670 g/s (21,140 t/yr) |
| Unfed lifetime | 15.8 yr |
| Exhaust β / v_eff | 0.332 / 0.276c |
| Muon punch-through | 84 MW (3000 fm), 58 W (5000 fm) |
| Off-axis EM leakage | < 1 kW |
| Sprint point (5×10⁷ kg) | 21,800 MN, 0.45c peak, 13 yr Alpha Cen |

## OPEN RESEARCH (tracked in the paper, §7.5)

- Usov pair-emission validity below 10⁹ K (T_eq = 3.56×10⁸ K is an extrapolation)
- BlackHawk+PYTHIA spectrum at kT ~ 10 GeV (replaces the analytic species sum)
- ISM interaction at >0.3c (forward shield co-evolution — research paper #1, feeds the
  forward shield tables; the SQM shell handles only the *drive's own* radiation)
- Bootstrap/seed phase is out of simulator scope (the ship spawns with an operating drive)
