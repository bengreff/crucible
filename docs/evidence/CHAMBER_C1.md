# C1: time-evolving chamber and nozzle ("valves open, ignite")

Status: criteria stated 4 October 2026, before any settled run. Results are appended below the line, they do not change the criteria.

## What C1 adds to the engine

- **Wall contour as data** (`Definition::contour`): (z, r) points sampled at the mesh stations by linear interpolation. It replaces the built-in cosine nozzle when given.
- **Chamber case** (`Case::Chamber`): the chamber starts as ambient gas at rest (composition, `backPressure`, `ambientTemperature`). Its exit draws ambient gas in if the flow reverses during start-up.
- **Supplies** (`Supply`): propellant streams through rings of the injector face. Each has an imposed mass flow, a frozen total enthalpy and a composition. The face pressure follows from the interior along the outgoing characteristic, and a face that would need supersonic inflow carries the choked (sonic) state. The valve ramps open linearly. Plate faces between supplies are slip walls, and their pressure force goes into the thrust ledger.
- **Igniter** (`Igniter`): a bounded energy deposit at constant power into a region of cells, shared by volume.
- **Time steps land on schedule events** (valve opening and full opening, igniter on and off). The reacting split shortens its step to the same events.
- **Chemistry modes** (`thermo::Chemistry`):
  - `FiniteRate`: kinetics.
  - `LocalEquilibrium`: constant-(u, v) equilibrium after every flow step, the verification limit.
  - `Frozen`.
- **Equilibrium trace-element cutoff**: elements below a mass fraction of 1e-14 in a cell are left out of that cell's equilibrium problem, and the species carrying them keep their amounts. Cells whose remaining equilibrium is trivial (for example pure N2) are skipped.
  - Why it is needed: advection leaves traces of 1e-50 and below in every cell, and Cantera's solvers took about 0.5 s per such cell (measured: one coarse step, 5.9 s of 6 s).
  - With the cutoff the same state costs 1 ms in total (measured).
  - The energy those traces could release is below 1e-13 of the cell's internal energy.

## Fast checks (core_verification, every build)

- The contour table reproduces the cosine nozzle at its own stations, also when the table's z origin is not the injector face.
- A shut chamber at ambient stays at rest (max |u| about 1.6e-12 m/s after 300 steps, measured).
- The igniter delivers exactly its energy, steps land on its switching times, and the mass and energy budgets close to 1e-12.
- The supply delivers exactly its scheduled flow, both mid-ramp (step average) and fully open, through the choked start and into subsonic filling, with budgets closed to 1e-12.

## The verification case (tests/chamber_study.cpp)

**Propellant and geometry**
- Premixed gaseous H2/O2 at O/F 5, 0.4 kg/s, total temperature 300 K, over the whole injector face. The valve opens at t = 0 with a 0.5 ms ramp.
- Mechanism h2o2.yaml (verification only).
- The chamber starts full of N2 at 1 kPa and 300 K.
- Throat radius 10 mm, chamber radius 25 mm, contraction ratio 6.25.
- Circular throat arcs of 20 mm on both sides (R = r_c/r_t = 2).
- 30 degree convergent cone, 15 degree conical exit, area ratio 10.
- The throat sits at mid-length, so an even nz puts a mesh station exactly on it.

**Igniter:** 0.5 J over 0.2 ms from t = 0.2 ms.
- In LocalEquilibrium mode the igniter has no role, because local equilibrium burns any premixed gas the moment it enters. Declared.
- Without molecular transport (C2), premixed flame propagation in FiniteRate mode is numerical diffusion. Declared.

**Reference:** the one-dimensional ideal rocket (the CEA rocket problem, `Mixture::idealRocket`, same h2o2.yaml thermo; the thermo itself matches CEA to 0.131%, see THERMO_VERIFICATION.md), evaluated at the simulated chamber pressure:
- **Chamber pressure p0:** the mass-flux-averaged equilibrium isentropic stagnation pressure of the last cylinder column.
- **c\*:** the two-dimensional throat passes Cd times the one-dimensional flow, so c\*_2D = c\*_1D(p0)/Cd.
  - Cd comes from the Kliegel-Levine series (AIAA J. 7(7), 1969), with coefficients as tabulated in Johnson and Wright, J. Fluids Eng. 130, 071202 (2008), Table 4.
  - Cd = 1 − a2/L² + a3/L³ − a4/L⁴, L = 1 + R.
  - At R = 2 and gamma about 1.21, Cd is about 0.9965.
  - The last retained term is 0.04%. This is the series uncertainty.
- **Vacuum Isp:** the conical divergence factor lambda = (1 + cos 15°)/2 applies to the momentum part only.
- **Vacuum thrust:** device thrust plus ambient pressure times exit area. This is valid only with a supersonic exit, which is checked.

## Criteria (LocalEquilibrium mode, stated before the settled runs)

1. **Settled:** relative drift of outlet mass flow, injector pressure and vacuum thrust below 1e-3 over the last 1 ms of the run.
2. **Mass:** settled outlet mass flow within 1e-3 of the supply. Mass and energy budgets below 1e-11 over the whole run.
3. **c\*:** within 0.5% of c\*_1D(p0)/Cd on the finest grid, with the error falling under refinement. The bracket is built from:
   - Cd series 0.04%;
   - gamma choice in Cd about 0.01%;
   - p0 averaging definition about 0.1%;
   - the rest is discretization.
4. **Vacuum Isp:** within 1.0% of the lambda-corrected shifting-equilibrium value on the finest grid, with the error falling under refinement. The bracket is built from:
   - the conical lambda approximation after a circular-arc throat, about 0.3%;
   - exit non-uniformity;
   - discretization.
5. **Grids:** 32x6, 64x12 and 128x24, if the finest fits the time budget. Otherwise two grids, with the observed trend stated and no order assumed.

FiniteRate mode is reported, not judged against CEA. Its vacuum Isp is expected to fall between the frozen and the shifting values. If it does not, that is a finding to explain.

---

## Results

### FiniteRate with the declared 0.5 J igniter: no light (32x6, 4 October)

- Measured: the peak temperature during and after the igniter window (0.2 to 0.4 ms) stayed below 600 K. It then fell to 300 K. The chamber settled into cold flow: injector pressure 0.98 MPa, 0.40 kg/s through the exit at 7 ms.
- Mechanism (derived): the igniter region (r < 10 mm) carries 16% of the face flow. During the window that is about 5e-6 kg of premix at 0.16 to 0.32 kg/s supply. Heating it from 300 K to 1500 K takes about 4.1 MJ/kg (mean cp about 3.4 kJ/(kg K) for the O/F 5 premix). That is about 20 to 30 J, against 0.5 J deposited. The 0.5 J was my case input and was never sized. This is a case-design error, not a solver fault.
- The engine predicts no light for an undersized igniter. That outcome is kept as evidence that ignition is not automatic in FiniteRate mode.

### Igniter sized for the FiniteRate rerun (stated before the rerun)

- Rule: deposit enough power to heat the full-flow throughput of the igniter region to 1500 K: 0.16 x 0.4 kg/s x 4.1 MJ/kg = about 260 kW (derived).
- Declared: 300 J over 1.0 ms from 0.2 ms (300 kW, margin 1.15 at full flow and more during the ramp). The region is unchanged.
- This is an input sized by a stated physical rule, set before the rerun. It is not adjusted to any observable.
- Expectation: the kernel lights. Without molecular transport (C2) and recirculation, a premixed flame in uniform inflow can only propagate by numerical diffusion. Whether it holds at the face, stands off, or blows out is a grid-dependent outcome, reported and not judged.

### LocalEquilibrium grid study (criteria 1 to 5)

Runs: `crucible_chamber_study eq 32 6` and `eq 64 12`, each to 8 ms, 5 threads, this Mac; `eq 128 24` to 5 ms (the time budget), built from c1aaa6c plus the shared-counter chemistry scheduling of daf408c. Between the two codes every change to the flow is behind `hasTransport()`, which is off in this case, so the inviscid chamber physics is the same (read from the diff, not rerun). Logs, history plots and contact sheets are in `c1/`; the raw CSVs were not kept in git.

| Grid | Wall time | p0 (MPa) | c\* sim (m/s) | c\* predicted (m/s) | c\* error | Isp_vac sim (s) | Isp_vac predicted (s) | Isp error | F_vac (N) |
|---|---|---|---|---|---|---|---|---|---|
| 32x6 | 341 s | 3.1558 | 2478.55 | 2438.06 | +1.66% | 427.99 | 429.09 | −0.26% | 1678.85 |
| 64x12 | 2375 s | 3.1164 | 2447.60 | 2437.83 | +0.40% | 428.59 | 429.07 | −0.11% | 1681.21 |
| 128x24 | 11457 s | 3.1067 | 2439.98 | 2437.77 | +0.09% | 430.01 | 429.06 | +0.22% | 1686.77 |

All values are measured, except the predicted columns, which are derived from the 1-D ideal rocket at the simulated p0 with Cd = 0.99645 and lambda = 0.98296. The frozen-flow Isp_vac bound is 412.7 s.

**Rerun after the axis-row clamp (5 October 2026, 32x6, Mac, built from the tree committed as 03d56e3; [log](c1/eq32_after_axis_clamp_2026-10-05.txt)).** The clamp (f09b077) bounds the axis-row r^2 curvature of the mass fractions by the cell value. Every printed result equals the recorded 32x6 run to the printed digit: p0, c\*, vacuum Isp, thrust, the drift and the outlet mass flow. Only the budgets differ, at round-off: mass −4.05e-13 (was −6.74e-13) and energy −8.35e-12 (was 1.25e-11). The run is no longer bit-identical. The cause may be the clamp or another change since c1aaa6c; they are not separated. 64x12 and 128x24 were not rerun.

Verdicts (finest grid 128x24; the 64x12 verdicts written before it are kept in git history):
1. **Settled: pass.** Drift over the last 1 ms on 128x24: outlet mass flow 1.02e-4, injector pressure 8.3e-5, vacuum thrust 8.5e-5 (limit 1e-3). It is larger than on 64x12 (5e-8 at 8 ms) because the run stopped at 5 ms; the thrust still rose 0.1 N per 0.5 ms at the end, about 6e-5 of it, far below the grid changes below.
2. **Mass: pass. Energy budget: pass on 128x24, fail as stated on the two coarser grids.**
   - Outlet mass flow 0.399996 kg/s against the supply 0.400000 (1e-5). On 64x12 it matched to 1e-6 at 8 ms; the 1e-5 here is the unfinished settling of verdict 1. Mass budgets −6.7e-13 (32x6), −1.4e-12 (64x12) and −1.10e-12 (128x24), limit 1e-11.
   - Energy budgets 1.25e-11 (32x6), −1.16e-11 (64x12) and 5.34e-12 (128x24), limit 1e-11.
   - Mechanism of the two failures (derived): the budget is divided by the energy of the initial fill. That fill is N2 at 1 kPa and 300 K in 2.5e-4 m^3, about −0.25 J. The same absolute error, about 3e-12 J, is about 1e-14 of the energy of the gas in the running engine (estimated at several hundred joules). The normalization, not a leak, is what fails the 1e-11 limit.
   - The criterion stays as written. Future chamber runs will also report the budget against the gas content.
3. **c\*: pass.** +0.09% on 128x24 against the 0.5% bracket, falling from +1.66% and +0.40%.
   - The simulated c\* differences between grids are 30.95 and 7.62 m/s (ratio 4.06, observed order 2.02). The Richardson limit, 2437.44 m/s, is 0.014% below the prediction at the 128x24 p0 (derived).
4. **Vacuum Isp: within the bracket, but the error did not fall, so the criterion fails as stated.**
   - +0.22% on 128x24 against the 1.0% bracket. The error was −0.26% on 32x6 and −0.11% on 64x12, so its size grew from 0.11% to 0.22% and its sign changed.
   - The simulated Isp rises with refinement: 427.99, 428.59 and 430.01 s. The increments are 0.60 s, then 1.42 s, so they are growing. The grids are not in the asymptotic range for thrust, and no limit or order is claimed.
   - The prediction is not exact. It carries the conical lambda after a circular-arc throat, about 0.3% (criterion 4's own bracket). So the difference from it is not a pure discretization error, and once the simulation passes the prediction, further convergence can raise it. That is a reading of the bracket, not a demonstration.
   - What would settle it (derived costs):
     - a 256x48 grid, about 25 hours to 5 ms (8 times the 128x24 cost); or
     - an independent inviscid 2-D reference for this contour (method of characteristics for the divergent section) in place of the conical lambda.
   - Neither has run.
5. **Grids: three.** c\* converges at an observed order of 2.02. The vacuum Isp does not yet show convergence (verdict 4).

In the frames, local equilibrium burns the premix as it enters ("valves open" without an igniter role, as declared). A starting shock crosses the nozzle by 0.1 ms. The chamber reaches about 3.3 MPa at the face by 2 ms.

### FiniteRate with the sized 300 J igniter (32x6, reported, not judged)

- **The kernel lit.** At about 0.22 ms the face pressure spiked to 1.2 MPa and the vacuum thrust swung between about −950 N and +1600 N. That is the deflagration of the premix that had filled the chamber since the valve opened (measured, history plot).
- The burning front then filled the chamber by 0.3 ms, and the flame held at the injector face. Without molecular transport the face-anchoring is grid-dependent numerical diffusion, as declared.
- **Settled at 8 ms:** drift below 2.4e-7. c\* 2477.20 m/s (+1.61% against the equilibrium prediction at the same grid). Vacuum Isp 424.31 s.
- **The expectation holds.** 424.31 s sits between the frozen 412.77 s and the shifting 429.09 s predictions. Against the equilibrium run on the same grid (427.99 s), finite-rate recombination in the nozzle costs 3.7 s, or 0.86%.
- Budgets: mass −1.6e-12. Energy 2.2e-11, which fails 1e-11 by the same normalization as above.
