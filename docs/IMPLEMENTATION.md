# First gas-flow implementation

17 September 2026. This records implemented behavior and its limits, not the full product scope.

## What runs

One C++ core powers both the native Qt/VTK application and `crucible_run`. The default experiment is a smooth converging-diverging nozzle, 0.60 m long, inlet/exit radii 0.035 m and throat radius 0.020 m at 36% of length. Gas is calorically perfect (gamma 1.4, R 287.05 J/kg/K). Inlet reservoir pressure/temperature are 300 kPa / 300 K; ambient pressure is 15 kPa. The UI permits reservoir pressure changes from 240–360 kPa. This is a numerical experiment, not an engine dataset.

The mesh is a fixed body-fitted grid in axial/radial coordinates. Every control volume represents a complete annular frustum. The displayed lower half mirrors the calculated upper half; it is not a second calculation. The prototype has no swirl, viscosity, heat conduction, turbulent transport, reactions, magnetic fields, particles, wall heating, or liquid injection. The nozzle preset establishes gas numerics without narrowing the project's combustion-chamber scope.

## Discretization

Evolved variables are mass, axial momentum, radial momentum, and total gas energy per volume. Face fluxes use an HLLC Riemann solver with HLLE fallback when the candidate star state is inadmissible. Piecewise linear primitive reconstruction uses minmod slopes in logical mesh directions; boundaries use unreconstructed adjacent states. Time advancement is two-stage SSPRK2. A first-order spatial option is available in the headless runner.

The finite-volume update uses exact axisymmetric volumes and integrated face-area vectors. The radial momentum source is pressure times the integral of 1/r over each volume. Its geometric cancellation with pressure fluxes preserves stationary gas in the curved mesh. Internal face exchanges are computed once and applied with opposite signs. The timestep is CFL times volume divided by the sum of face acoustic/advection rates. A nonpositive density/internal-energy stage causes timestep halving, up to 14 attempts; no cell energy/density floors repair an invalid accepted solution. Unsupported boundary states stop the run with an explanation.

These choices do not establish global second-order accuracy on arbitrary curved grids. Boundary treatment, axis treatment, reconstruction geometry, and multidimensional convergence still need systematic study. Shock convergence and stationary balance are narrower evidence.

## Boundaries and initial state

- Axis: zero-area radial face and the axisymmetric momentum source.
- Outer wall: stationary reflected slip state; exactly zero mass/energy transfer.
- Inlet: subsonic reservoir with prescribed total pressure/temperature and zero radial velocity; outgoing acoustic invariant determines the axial speed. Reverse/supersonic inlet states are unsupported.
- Outlet: extrapolation for supersonic axial outflow; otherwise specified static pressure with an acoustic/isentropic reconstruction. Backflow is unsupported (apart from roundoff tolerance around rest).

The initial nozzle state is a **prepared quasi-one-dimensional isentropic flow**, with radial velocity estimated from wall slope. It is not a simulated startup from stationary gas, and it is not an exact 2D steady solution. Its total enthalpy has a small radial-kinetic contribution beyond the 1D prescription. The actual solver evolves this field and its transients. A nearby steady mass-flow comparison alone is therefore insufficient verification; the separate shock-tube and pressure-change tests matter.

Uniform-duct and shock-tube definitions are verification presets using this same Euler update, with extrapolated axial boundaries. They are not separate simulation engines. Before broader experiment editing, replace preset-specific initialization/boundary selection with explicit experiment data.

## Measurements and ownership

Mass/energy balance errors are `(current inventory - initial inventory - integrated net boundary transfer) / initial inventory`. Integrated transfers use the exact numerical boundary fluxes and RK weights of accepted updates. Inlet/outlet rates and outlet force are averages over the most recent accepted timestep; before the first step these are zero/unmeasured. Exit Mach is area weighted from the last cell column. Outlet-plane force is axial numerical momentum flux minus ambient pressure times exit area. **It is not total device thrust**: inlet and wall forces/storage must be included in a full device momentum budget.

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

## Next work

1. Broaden genuinely multidimensional verification, boundary applicability/error tests, and full axial momentum/force accounting. Avoid adding chemistry before the gas substrate is understood.
2. Add experiment save/load, accepted-control history and time traces; make restart semantics and control history clearer. Preserve one engine and the native workflow.
3. Select a published propulsion geometry/performance dataset with adequate inputs and uncertainties. Establish equilibrium/frozen gas thermodynamics and the path to an in-domain reacting chamber. No real engine is validated yet.
4. Resolve AMReX versus body-fitted infrastructure with the planned RZ/curved-wall/field feasibility case before committing to a general mesh architecture. This prototype does not settle that decision.
