# C2 step 2: a freely propagating laminar flame against Cantera's free flame

Status: done (5 October 2026). All four criteria pass at 10 um (see Verdict). The 10 um run was made on backhouse after the open-face fix, the 40 and 20 um runs on the Mac before it. Criteria were stated 4 October 2026 in the header of `tests/flame_study.cpp`, before any engine flame result, and have not changed. The case and the outlet were changed twice after runs, and the header records both changes. Raw outputs are in `c2/`.

## The test

- **Case.** H2/air at equivalence ratio 1, 300 K, 1 atm, mechanism h2o2.yaml (verification use only).
  - Duct: radius 1 mm, length 6 mm, one cell high.
  - The injector end is a closed adiabatic slip plate; the far end is the subsonic outlet at 101325 Pa.
  - At t = 0 the gas is at rest: premix for z < 4 mm, its constant-(h, p) equilibrium beyond.
  - The flame runs toward the closed end into gas at rest. There is no igniter.
- **Reference.** Cantera Flow1D free flame, same mechanism, mixture-averaged transport, refined three times. S_L = 2.330470 m/s (716 points). Its two finest levels differ by 0.089% (precondition < 0.2%). Its own S_c against its S_L is -0.029% (precondition < 0.2%). Both preconditions pass.
- **Measurements.**
  - S_c: the H2 balance of the region ahead of x_f + 1 mm, including the diffusive flux there.
  - S_d = -dx_f/dt over the last 0.2 ms.
- **Criteria.**
  1. At dz = 10 um, S_c within 1% of S_L.
  2. |S_c/S_L - 1| falls at each refinement, 40, 20, 10 um.
  3. At 10 um, S_d within 1% of S_L.
  4. Settled: the S_c means over the two halves of the last 0.2 ms differ by less than 0.2%.

## Run 1: fixed-pressure outlet (40 um) — the duct resonates

Measured, 1 ms, 581 s wall: S_c 2.806549 m/s (+20.428%), S_d 2.860950 m/s (+22.763%), drift over the last 0.2 ms 3.115%. Mass budget 7.4e-15, energy budget 4.7e-15.

- The flame did not settle. S_c swung between 2.0 and 3.0 m/s through the run, and the pressure swung between 92 and 114 kPa (`c2/reflecting_history.png`).
- Mechanism (derived, consistent with the plots): with a fixed-pressure outlet the duct is a closed-open resonator. The plate is closed and the outlet is a pressure node, which gives a quarter wave near 20 kHz in the hot gas. The flame's heat release drives it.
- The reference is a flame in open space, so this run cannot meet the criteria by construction. It is kept as evidence; it is not a verdict on the flame physics.
- A 20 um run with the same outlet was started and stopped before its first output.

## The outlet change

The outlet was made partially non-reflecting (Poinsot and Lele 1992), as an option on every subsonic outlet (`Definition::outletRelaxation`, sigma). The flame study uses sigma = 0.25. Zero keeps the old fixed-pressure node, bit for bit. The ambient became the burned gas, so backflow at the outlet cannot draw cold gas.

Theory: the incoming characteristic relaxes toward the back pressure at the rate K = sigma (1 - M^2) a / L. The reflection is -1 / (1 + 2 i omega / K), and the impulse response is -(K/2) exp(-K t / 2). The derived reflection for sigma = 0.25 at 20 kHz with burned gas at the outlet is about 0.13.

Core check (`outletChecks` in `tests/core_tests.cpp`, criteria stated first). A Gaussian pulse in a 1 m duct reflects off the outlet; the reference signal comes from a 2 m duct.

| Check | Criterion | Measured |
|---|---|---|
| sigma 0, reflected/reference peak | -1 within 0.02 | -1.00266 |
| sigma 100, reflected peak against the reference convolved with the theory impulse response | within 10% | -0.823 against -0.846 |
| sigma 0.25, largest reflected excursion | < 0.02 | 0.00688 |
| 100 K hot spot leaving at 15 m/s (sigma 0.25), largest pressure change | < 1e-3 of ambient | 5.8e-4 |

Raw output: `core_verification_outlet_2026-10-04.txt`.

### How the outlet got there (in order, measured)

1. **Backflow switching.** In the first version sigma 100 reflected a peak of -1.018 against the theory's -0.846 (waveform difference 0.456). At zero mean flow the outlet velocity changed sign every step, and the backflow branch took a different pressure.
   - Fix: backflow keeps the same p and u and takes only the ambient entropy and composition.
2. **Extrapolated incoming invariant.** Relaxing the face's extrapolated incoming invariant gave unbounded linear growth.
   - Fix: the incoming invariant comes from the last cell's own state.
3. **Rate four times too fast.** Fitting an exponential to the reflected signal gave K_eff/K = 3.57 to 3.96 for sigma 0.05 to 100 and nz 800 and 1600.
   - Mechanism (derived): the cell lags the face by half a cell, and Poinsot and Lele write the relaxation rate as K/2. That makes the per-step weight beta = 0.25 sigma (1 - M^2) dz / L, not sigma (1 - M^2) dz / L.
   - With the factor, K_eff/K is 0.981 (sigma 0.25), 0.964 (sigma 10) and 0.906 (sigma 100) at nz 800. Committed in 3197bab.
4. **Entropy waves made sound.** The first non-reflecting 40 um flame run blew up. Diagnostic stops measured the pressure range in the duct: 101.5 to 102.4 kPa at 45 us, 102.1 to 103.5 kPa at 52 us, 103.3 to 121.9 kPa at 55 us. By 0.1 ms T_max had fallen to 190 K, so the gas had been expanded far below the 300 K fill (inferred). After 0.2 ms CVODES failed with flag -4 (`c2/f40_isentropic_outlet_log.txt`).
   - Mechanism (measured in a perfect-gas duct): the outlet used isentropic Riemann invariants, u +- 2a/(gamma-1), mixing the cell's and the face's sound speeds. An entropy gradient at the outlet then counts as an acoustic wave, and with gamma near 1.2 the factor 2/(gamma-1), about 10, amplifies it. A 1 K hot spot convected out at 15 m/s grew exponentially to about 870 kPa.
   - Fix: the characteristics are W = p +- rho a u_z with the face's impedance. The entropy is the face's own. 1 K and 100 K hot spots now stay at the start-up level (about 1.4 kPa), and the reflection checks are unchanged. Committed in 32adf4e, with the hot spot as a core check.

## Run 2: non-reflecting outlet (sigma 0.25)

| dz | S_c (m/s) | S_c / S_L - 1 | S_d (m/s) | S_d / S_L - 1 | drift, last 0.2 ms | wall |
|---|---|---|---|---|---|---|
| 40 um | 2.361664 | +1.339% | 2.352410 | +0.941% | 0.031% | 1057 s |
| 20 um | 2.341939 | +0.492% | 2.335204 | +0.203% | 0.023% | 5369 s |
| 10 um | 2.336968 | +0.279% | 2.331505 | +0.044% | 0.018% | 6084 s (backhouse) |

40 um, measured: mass budget -5.2e-15, energy budget 5.5e-15. 20 um: -3.6e-16 and -1.4e-14 (`c2/f20_log.txt`). 10 um: 8.1e-14 and 1.6e-13 (`c2/f10_log.txt`; history and profiles in `c2/f10_history.png` and `c2/f10_profiles.png`). The 10 um run: backhouse (WSL2 Ubuntu, gcc 13, Cantera 3.2.0), commit 4e2ee5c, after the open-face fix (1e0fd9e), 16 threads beside other users' jobs, 5 October 03:02 to 04:46.
- The start-up transient dies out by about 0.3 ms. After that S_c sits within 0.02 m/s of its final value (`c2/f40_history.png`).
- The pressure spread in the duct falls from 1.2 kPa at start-up to 78 to 91 Pa after 0.5 ms. The resonance of run 1 is gone.
- The temperature, H2, H and OH profiles lie on the free flame's (`c2/f40_profiles.png`). T_max (2273 K at 1 ms) stays below T_ad (2388 K) because the products are still recombining behind the flame, as in the free flame.

**The outlet pressure drifts above the back pressure.** Measured: the duct sat at 101.72 kPa at 0.5 ms and 101.58 kPa at 1 ms, against 101.325 kPa. At 1 ms the gas ahead of the flame was at 300.24 K and 101.607 kPa.
- Mechanism (inferred): the products are still recombining as they reach the outlet, so the velocity still rises there (du/dz about 104 /s over the last two cells). A partially non-reflecting outlet only reaches back pressure when the outgoing wave vanishes. In steady flow the outlet settles at about p_inf + rho a L (du/dz) / sigma. That is 318 Pa with rho a estimated at about 127 kg/(m^2 s) (guessed molar mass 24.5 g/mol, gamma 1.2), against 254 Pa measured. The offset scales as 1/sigma and does not depend on the grid.
- Effect on the comparison (derived with Cantera, same free-flame setup, at the engine's unburned state at 1 ms): S_L = 2.333328 m/s, +0.123% above the stated reference.
  - S_d is measured against gas at rest, so an exact engine would show S_d +0.12%.
  - S_c divides the mass burning rate by the reference's unburned density (300 K, 1 atm), as stated. The engine's unburned gas is 0.198% denser, so an exact engine would show S_c +0.32%.
  - The criteria and the reference are unchanged. These offsets are part of the gap and are reported, not subtracted.

## Verdict (5 October 2026)

| Criterion | Measured | Limit | Result |
|---|---|---|---|
| 1. S_c at 10 um | +0.279% | within 1% | pass |
| 2. \|S_c/S_L - 1\| falls at each refinement | 1.339%, 0.492%, 0.279% | falls | pass |
| 3. S_d at 10 um | +0.044% | within 1% | pass |
| 4. Settled at 10 um: S_c means over the two halves of the last 0.2 ms | 0.018% | < 0.2% | pass |

- At 10 um the temperature, H2, H and OH profiles lie on the free flame's (`c2/f10_profiles.png`). After the start-up transient (about 0.3 ms), S_c sits within 0.02 m/s of its final value (`c2/f10_history.png`).
- Criterion 2 compares across a code change: the 10 um run has the open-face fix, the 40 and 20 um runs do not. The fix changes only the end columns of the duct. Its effect on these flames was not measured.
- The outlet offset derived above (an exact engine would show S_c +0.32% and S_d +0.12%) is as large as the 10 um gaps. So at 10 um this check cannot separate the engine's own error from the outlet's offset; it bounds the sum. The gap's fall from 20 to 10 um (a ratio of 1.76, against 2.72 from 40 to 20; derived) is consistent with the offset dominating (inferred).
- What this verifies: flow, mixture-averaged transport and CVODES chemistry, coupled through the Strang split, reproduce a laminar flame speed within 0.3% at 10 um. It does not test turbulence, the closure or walls.
