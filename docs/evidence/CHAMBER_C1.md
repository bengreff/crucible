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

(appended after the runs)
