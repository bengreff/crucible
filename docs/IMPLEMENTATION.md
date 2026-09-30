# First gas-flow implementation

17 September 2026. This records implemented behavior and its limits, not the full product scope.

## What runs

One C++ core powers both the native Qt/VTK application and `crucible_run`. The default experiment is a smooth converging-diverging nozzle, 0.60 m long, inlet/exit radii 0.035 m and throat radius 0.020 m at 36% of length. Gas is calorically perfect (gamma 1.4, R 287.05 J/kg/K). Inlet reservoir pressure/temperature are 300 kPa / 300 K; ambient pressure is 15 kPa. The UI permits reservoir pressure changes from 240–360 kPa. This is a numerical experiment, not an engine dataset.

The mesh is a fixed body-fitted grid in axial/radial coordinates. Every control volume represents a complete annular frustum. The displayed lower half mirrors the calculated upper half; it is not a second calculation. The prototype has no swirl, viscosity, heat conduction, turbulent transport, reactions, magnetic fields, particles, wall heating, or liquid injection. The nozzle preset establishes gas numerics without narrowing the project's combustion-chamber scope.

## Discretization

Evolved variables are mass, axial momentum, radial momentum, and total gas energy per volume. Face fluxes use an HLLC Riemann solver with HLLE fallback when the candidate star state is inadmissible. Piecewise linear primitive reconstruction uses minmod slopes in logical mesh directions. Boundary rows/columns (wall, inlet, outlet) use the limited one-sided slope minmod(c-n1, n1-n2), with density/pressure face values kept at least half the cell value; the axis row keeps a zero slope (measured second order there). Inlet/outlet characteristic models receive the reconstructed interior face state. (Before 30 September 2026 boundary rows used zero slopes, which made wall-row errors first order.) Time advancement is two-stage SSPRK2. A first-order spatial option is available in the headless runner.

The finite-volume update uses exact axisymmetric volumes and integrated face-area vectors. The radial momentum source is pressure times the integral of 1/r over each volume. Its geometric cancellation with pressure fluxes preserves stationary gas in the curved mesh. Internal face exchanges are computed once and applied with opposite signs. The timestep is CFL times volume divided by the sum of face acoustic/advection rates. A nonpositive density/internal-energy stage causes timestep halving, up to 14 attempts; no cell energy/density floors repair an invalid accepted solution. Unsupported boundary states stop the run with an explanation.

Measured second-order convergence on the default nozzle (below) is evidence for this smooth body-fitted grid, not for arbitrary curved grids. Upwind (HLLC) dissipation scales with sound speed, so low-Mach observables are sensitive: at exit Mach ~0.15 mass flow amplifies total-pressure error ~1/(gamma M^2) ~ 33x (venturi evidence below). Combustion chambers run at Mach 0.1-0.3; decide on a low-Mach correction with chamber evidence, not before.

## Boundaries and initial state

- Axis: zero-area radial face and the axisymmetric momentum source.
- Outer wall: stationary reflected slip state; exactly zero mass/energy transfer.
- Inlet: subsonic reservoir with prescribed total pressure/temperature and zero radial velocity; outgoing acoustic invariant determines the axial speed. Reverse/supersonic inlet states are unsupported.
- Outlet: extrapolation for supersonic axial outflow; otherwise specified static pressure with an acoustic/isentropic reconstruction. Backflow is unsupported (apart from roundoff tolerance around rest).

The initial nozzle state is a **prepared quasi-one-dimensional isentropic flow**, with radial velocity estimated from wall slope. It is not a simulated startup from stationary gas, and it is not an exact 2D steady solution. Its total enthalpy has a small radial-kinetic contribution beyond the 1D prescription. The actual solver evolves this field and its transients. A nearby steady mass-flow comparison alone is therefore insufficient verification; the separate shock-tube and pressure-change tests matter.

Uniform-duct and shock-tube definitions are verification presets using this same Euler update, with extrapolated axial boundaries. They are not separate simulation engines. Before broader experiment editing, replace preset-specific initialization/boundary selection with explicit experiment data.

## Measurements and ownership

Mass/energy balance errors are `(current inventory - initial inventory - integrated net boundary transfer) / initial inventory`. Integrated transfers use the exact numerical boundary fluxes and RK weights of accepted updates. Inlet/outlet rates and all axial forces are averages over the two RK stages of the most recent accepted timestep; before the first step these are zero/unmeasured. Exit Mach is area weighted from the last cell column.

Axial momentum / device thrust (item 3, 30 September 2026). Forces on the gas are positive in +z (the exhaust direction); thrusts are forces on the device, positive against the exhaust. The gas obeys `dP_z/dt = inletMomentumFlux - outletMomentumFlux + wallAxialForce + bodyAxialForce`, each term being the exact numerical flux (area times HLLC axial momentum flux) summed over the boundary faces with the RK weights. `momentumBalanceError` is `(P_z - P_z0 - integrated net source) / (|P_z0| + integrated gross exchange)`. Device thrust is the reaction of every force the device applies to the gas, minus ambient pressure on the closed exterior: `deviceThrust = inletMomentumFlux + wallAxialForce + bodyAxialForce - backPressure * exitArea`. This equals `exitPlaneThrust + dP_z/dt`, so the two agree only in steady state. The supply plane is treated as part of the device (a reservoir face). `bodyAxialForce` is the slot for a future volumetric force such as the Lorentz force `J x B`; its volume integral enters the gas budget and its reaction acts on the coils, so a magnetic nozzle's thrust enters through the same ledger without new accounting. It is zero today.

Measured (80x12, 8 ms, second order): C_F device 1.37282 vs quasi-1D ideal 1.37351 (ratio - 1 = -5.0e-4); (device - exit plane)/device 3.6e-6; momentum residual 1.4e-16. At 160x24 the components are supply plane 1183.66 N, wall force on gas -608.16 N, ambient 57.73 N, device thrust 517.77 N. Grid study at 8 ms (`docs/evidence/nozzle_thrust_convergence.csv`): C_F/ideal - 1 = -1.76e-3, -5.0e-4, -5.7e-5, +1.4e-4 for 40x6 to 320x48. Successive device-thrust changes shrink by 2.8 and 2.2 (observed order about 1.5 then 1.1). On the finest grid device and exit-plane thrust differ by 0.107 N (2e-4), which is dP_z/dt from the transverse ringing still present at 8 ms; the exit-plane value converges faster (changes shrink by 3.3 then 6.2). Inferred: the 2e-4 level is not yet steady, so thrust claims at that level need a longer run or time averaging.

The worker alone owns `Flow`. Controls are accepted between timesteps and recorded with sequence, applied value and physical time. Pending pressure requests may coalesce; applied events remain recorded. Immutable snapshots cross a short mutex-protected pointer exchange. The GUI owns all Qt and VTK objects. Pause acknowledges only after publishing the final accepted state. Restart discards this prototype's current state and restores the default initial experiment; there is no history comparison or checkpoint recovery yet.

Snapshots currently copy every cell, nominally at 20 Hz. This is appropriate for the default 3,840 cells, not a completed million-cell display architecture. CSV exports a single accepted field with time and SI units; it is neither an experiment definition nor a restart file.

## Evidence from this session

Tested locally with AppleClang 17, C++20, Release and Debug ASan/UBSan; Qt 6.11.2 and VTK 9.7.0 from Homebrew. These are the tested versions, not a lockfile or a redistributable dependency bundle. Windows/Linux builds and packaging remain untested.

Automated core checks cover cylinder volume, state conversion, identical-state flux, straight-duct uniform flow, curved-wall rest preservation, independent analytic Sod shock density convergence, radial acoustic mode convergence, nozzle mass flow and conservation, physical response to changed reservoir pressure, snapshot immutability, control acceptance, pause stability and replay from the accepted control times.

- A small-amplitude cylindrical acoustic eigenmode uses J0(k r) pressure with J1(k R)=0 at the wall. At a quarter period, area-weighted radial-velocity error normalized by acoustic velocity amplitude falls from 0.008800 (16 radial cells) to 0.0009632 (48). This independently exercises radial transport and the cylindrical source.
- Sod density mean absolute error at t=0.15: 0.008686 at 80 axial cells, 0.003451 at 240.
- 64×8 nozzle at 4 ms: outlet 0.878703 kg/s against ideal choked 0.879654 kg/s. Acceptance is within 2%; this checks the gradual-nozzle limit, not measured propulsion performance.
- After increasing reservoir pressure by 10%, the same nozzle settles within 2% of the correspondingly scaled ideal mass flow; mass/energy residuals remain below 1e-11.
- Separate 160×24, 20 ms run: outlet 0.879545 kg/s, exit Mach 2.65898; mass/energy relative budget residuals below 1e-14, no rejected steps.
- Native smoke test runs the actual worker/rendering path, changes pressure, waits for acceptance, pauses, and saves a screenshot. Screenshot was visually inspected. This is not exhaustive manual interaction testing.

Performance observations on this development Mac (single numerical worker, Release; not controlled hardware benchmarks): 3,840 cells / 54,960 steps / 20 ms physical duration took 26.6 s. One million cells / 10 steps took 1.24 s compute plus 0.043 s initialization (after moving repeated column initialization out of the radial loop; previously 3.58 s), advancing only 0.199 microseconds. Both yield about 8 million cell updates/s. The million-cell/minutes aspiration is **not achieved for useful flow durations**. Priorities are profiling, scalable storage/display, parallel execution, and justified timestep strategies; increasing throughput alone does not remove acoustic timestep restrictions.

## Verification review, 30 September 2026 (branch `claude/verify-core`)

Re-measured on this Mac: every number above reproduced (Sod, acoustic mode, nozzle mass flow, residuals, 26.6 → 27.1 s, 8.2 M cell-updates/s, smoke screenshot). The review found no wrong formula: face area vectors close exactly per cell, the p/r source equals the exact meridional-area integral, HLLC star states and inlet/outlet invariants are standard. It found first-order boundary reconstruction (fixed above) and untested subsonic-outlet, contact and strong-rarefaction behavior (now tested).

Nozzle grid study (`crucible_convergence`, 8 ms, nr = 0.15 nz; CSVs in `docs/evidence/`). The exact steady solution is isentropic with uniform total enthalpy:

| Grid | Entropy L1 | order | Wall-row entropy | order | Mass flow vs quasi-1D |
|---|---|---|---|---|---|
| 40x6 | 8.96e-4 | | 8.90e-4 | | -1.25e-3 |
| 80x12 | 2.11e-4 | 2.09 | 2.19e-4 | 2.02 | -3.68e-4 |
| 160x24 | 5.14e-5 | 2.04 | 5.48e-5 | 2.00 | -7.46e-5 |
| 320x48 | 1.28e-5 | 2.01 | 1.37e-5 | 2.00 | -3.30e-5 |
| 640x96 | 3.21e-6 | 1.99 | 3.47e-6 | 1.98 | -2.15e-5 |

First-order scheme: entropy order 1.02, 1.01 (confirms the measure). Before the boundary fix the wall row converged at order 1.23, 1.12. Total enthalpy and mass flow stop improving at ~2e-5 on fine grids: the prepared initial state excites transverse acoustic modes (measured 5-7.6 kHz; first radial mode estimate ~6 kHz) that decay slowly, more slowly on finer grids. Time-averaged (4-24 ms) mass-flow error: -3.67e-4, -7.49e-5, -3.18e-5 (80, 160, 320); Richardson extrapolation lands ~2.4e-5 below quasi-1D, consistent with Hall's leading 2D discharge correction (~2.5e-5 for the upstream throat curvature, Rc/rt = 31.5). Consistent, not proven, given the residual ringing.

New tests: exact Riemann solver cross-checked against the Sod constants; planar Sod in a 5-row duct (rows identical, radial velocity zero; L1 0.00737 → 0.00227); Toro 1-2-3 rarefaction (0.0135 → 0.00451, no failures); 1e5 pressure-ratio shock (0.170 → 0.060); HLLC antisymmetry, rotation invariance, supersonic upwinding and exact stationary contact; nozzle entropy order > 1.8 including the wall row; choked nozzle with internal normal shock (back pressure 0.7 p0, from rest): exit Mach 0.26773 vs quasi-1D 0.26803, mass flow -0.11%. Slow-labelled (`ctest -L slow`, 37 s): subsonic venturi (0.985 p0 back pressure, from rest, 160 ms) mass flow -9.49% (40x6) → -2.63% (80x12), order 1.85; the error is total-pressure loss after the throat and a spurious gain upstream, both converging. Acoustic mode error with the new boundary slopes: 0.00422 (16) → 0.000401 (48), order 2.14. ASan/UBSan: pass (fast suite).

## Next work

1. Broaden genuinely multidimensional verification, boundary applicability/error tests, and full axial momentum/force accounting. Avoid adding chemistry before the gas substrate is understood.
2. Add experiment save/load, accepted-control history and time traces; make restart semantics and control history clearer. Preserve one engine and the native workflow.
3. Select a published propulsion geometry/performance dataset with adequate inputs and uncertainties. Establish equilibrium/frozen gas thermodynamics and the path to an in-domain reacting chamber. No real engine is validated yet.
4. Resolve AMReX versus body-fitted infrastructure with the planned RZ/curved-wall/field feasibility case before committing to a general mesh architecture. This prototype does not settle that decision.
