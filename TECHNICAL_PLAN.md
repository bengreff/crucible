# CRUCIBLE — Technical plan

17 September 2026; current-state section, chemical plan and dated plan to January updated 4 October 2026. Target architecture; most capabilities remain planned. Implemented behavior and its verification are indexed in [docs/evidence/README.md](docs/evidence/README.md); the live state is in [docs/SESSION_HANDOFF.md](docs/SESSION_HANDOFF.md). Read [VISION_SCOPE.md](VISION_SCOPE.md) for purpose and scope and [RESEARCH.md](RESEARCH.md) for evidence and unresolved science. Current user instructions take precedence. Architecture changes must preserve those physical and human requirements.

**How to read decisions:** product requirements are fixed by the vision and Ben's instructions; selected architecture is the current implementation direction; algorithm candidates require the stated tests. Do not promote a candidate to supported physics by implementing its interface. Keep one authoritative home for each decision: purpose in the vision, model evidence in research, implementation here, historical lessons in RESEARCH.md section 7.

## Current state and the plan to January (4 October 2026)

**Implemented and verified** (evidence index: `docs/evidence/README.md`): a fixed body-fitted axisymmetric finite-volume mesh with uniform axial spacing; HLLC with limited primitive reconstruction and SSPRK2; a thermally perfect multi-species mixture (NASA-7 data from the mechanism) with species carried by the mass flux; Cantera thermochemistry behind `adapters/`; per-cell stiff kinetics (CVODES BDF over Cantera rates) joined to the flow by symmetric Strang splitting; a device-thrust momentum ledger; a Qt/VTK app operating a gas nozzle live (Mac only so far). Release is the default build. C1 (4 October): wall contours from data, gas supply faces with valve ramps, an igniter, local-equilibrium and frozen chemistry modes, and a "valves open, ignite" chamber case checked against the ideal rocket (`docs/evidence/CHAMBER_C1.md`). C2 under way: mixture-averaged transport properties match Cantera; the transport operator (viscous stress, conduction, diffusion; no-slip, slip, isothermal and adiabatic walls) passes its solution and budget checks, and two of its four truncation-order criteria still fail (`docs/evidence/TRANSPORT_C2.md`). Wall-clustered rings and SST-2003 are in the engine; the SST wall distance, homogeneous decay and turbulent operator checks pass, and the fully developed pipe check is queued (`docs/evidence/turbulence_verification_2026-10-04.txt`). AMReX is not adopted; the fixed custom mesh has been sufficient so far, and the decision stays conditional (section 6).

**Not yet implemented:** the FreeFlame check of the transport operator (running), turbulent inflow conditions on supplies and the nozzle inlet (so turbulence runs on ducts only), the turbulence-chemistry closure (PaSR), interleaved fuel and oxidizer rings, liquid injection, breakup and evaporation, persistence, experiment editing, other platforms.

**Direction** (Ben, 4 October 2026; recorded in VISION_SCOPE *Changes since the pivot*): one time-evolving reacting engine; shifting equilibrium as the validated baseline produced by that engine; liquid LOX injection with breakup and evaporation in the shared engine before RL10 ("Add evaporation first"; no pre-vaporized RL10); every chemical run starts "valves open, ignite"; RL10A-3-3A next; then CEA sweeps for trends; then a database of hundreds of engines. By January all four regimes work in the same engine: chemical, magnetic nozzle, pulsed fusion (ICF post-burn pulse) and antimatter. Chemical stays lean and correct; the parked MHD spike returns right after it; the first draft of the laboratory sandbox app follows the chemical milestone. The dated order is in *Dated plan to January* below. The steps here are the selected architecture for the chemical work. Each lands with its verification.

1. **Always a time-dependent calculation.** Every chemical experiment marches the unsteady reacting equations from a declared initial state. No steady-state solver, no pseudo-time acceleration that alters the transient. Steady performance is the average over a late, settled window of the history; the settling test and window are recorded with the result.
   - **Start: "valves open, ignite"** (Ben). The chamber and nozzle start filled with ambient gas at rest (declared composition, pressure, temperature). The propellant valves open on a declared ramp schedule (valve travel is an input, not modeled equipment), the igniter fires, the flame establishes, the chamber pressurizes, and the flow settles.
   - A run that cannot afford the full start says so with the measured cost; a declared partly developed initial state that still evolves in time is then a development fallback, never the reported RL10 case.
2. **Chemistry.** Finite-rate kinetics in the existing Strang slot is the production model. Local equilibrium per cell (constant-u,v equilibrium of the cell's own elements, in the same slot) is its fast-chemistry limit and is used to verify against CEA shifting equilibrium; frozen chemistry is the other limit. Use a hydrogen/oxygen mechanism validated to rocket chamber pressures (the Burke et al. 2012 or Li et al. 2004 class) and record its validity range; Cantera's `h2o2.yaml` is a GRI-3.0 subset suitable for verification only.
3. **Geometry as data.** The wall is r(z) from a published contour; for an area-against-station table r = r_t sqrt(A/A*), so no shape is assumed. Uniform axial spacing first; non-uniform axial spacing if measured cost requires it.
4. **Gas supplies.** The injector face is a boundary object. Each radial face carries a stream with imposed mass flux, total enthalpy and composition; static pressure comes from the interior (one outgoing characteristic), so the inlet is well posed; backflow is flagged, not silently reversed. A coaxial-element injector becomes annular averages: either premixed at the face (a declared limit with no mixing loss, used for verification) or interleaved oxidizer and fuel rings, whose mixing is left to the transport closure. A valve ramp scales each stream's mass flux in time.
5. **Liquid propellant: injection, breakup and evaporation** (Ben, "Add evaporation first"; required before RL10). Liquid enters as stochastic Lagrangian parcels in the RZ domain: each parcel carries position, three velocity components, drop diameter, temperature and a statistical weight, and exchanges mass, momentum, energy and species with the gas through sources inside the shared step, normalised to physical ring volumes. Primary atomization of a coaxial element is not resolved in 2-D axisymmetry: the initial drop size distribution comes from a published coaxial-injector correlation and is carried as a declared bracket, the leading uncertainty of the model. Secondary breakup uses a KH-RT class model; evaporation uses a film-theory model (Abramzon-Sirignano class) with Cantera gas properties and cited liquid O2 properties. The liquid's enthalpy is the real feed enthalpy, so the latent heat is honoured without a separate correction. Applicability is checked as the state evolves: subcritical evaporation requires chamber pressure below the O2 critical pressure. The same particle container later carries energetic antimatter products with their own laws. Verification: single-drop evaporation against a reference solution, mass and energy budgets closing with parcels present, and a published subcritical LOX/GH2 spray flame as a component case.
6. **Igniter.** A bounded energy deposit (declared energy, volume and duration) or a defined hot-gas supply, counted in the energy budget once.
7. **Unresolved transport, applied uniformly.** Molecular transport first (mixture-averaged, Cantera's fits, verified against Cantera and a laminar flame speed from Cantera's FreeFlame). Then SST-2003 URANS with a declared wall treatment and turbulent Prandtl/Schmidt numbers, then PaSR for turbulence-chemistry interaction. The same closure acts on every stream, the chamber and the nozzle; there is no combustion-only mixing shortcut. Until it exists, results that depend on mixing are not claimed, and premixed runs isolate everything else.
   - **SST-2003 as built** (equations from the NASA Turbulence Modeling Resource `sst.html`, extracted 4 October 2026). The model is the standard SST with the 2003 changes and nothing else:
     - production limited to min(P, 10 beta* rho omega k) in both equations;
     - mu_t = rho a1 k / max(a1 omega, S F2) with the strain invariant S = sqrt(2 S_ij S_ij), which includes the hoop strain u_r / r;
     - gamma1 = 5/9, gamma2 = 0.44, CD_komega floor 1e-10;
     - other constants: sigma_k 0.85 / 1.0, sigma_omega 0.5 / 0.856, beta 0.075 / 0.0828, beta* 0.09, a1 0.31, blended by F1.
     - P = tau_ij du_i/dx_j with the full Boussinesq stress, including -2/3 rho k delta_ij.
   - **Energy accounting** (RESEARCH *Gas turbulence*: account consistently): k is part of the total energy, and the thermodynamic state uses e = E/rho - |u|^2/2 - k.
     - The Reynolds stress includes -2/3 rho k delta_ij, carried in the transport flux; its work u . (-2/3 rho k I) is therefore in the energy flux.
     - The energy flux carries (mu + sigma_k mu_t) grad k.
     - Turbulent heat and species fluxes use Pr_t and Sc_t (declared inputs; 0.9 and 0.7 first, sensitivity reported). The species flux adds mu_t / (rho Sc_t) to every D_km in the existing mixture-averaged form, which with the correction velocity is exactly -(mu_t / Sc_t) grad Y_k.
   - **State and numerics.**
     - rho k and rho omega are two more conserved fields, carried by the mass flux like the species and reconstructed like them.
     - Their sources are integrated in a Strang split inside each step (half step, transport step, half step), with the mean-flow coefficients frozen per half step.
     - The sources are integrated per cell at fixed rho and E (decayed k becomes heat) by MPRK22, the second-order positive modified Patankar Runge-Kutta scheme (Kopecz and Meister, BIT 58, 2018): production explicit, destruction weighted by the new over the old value. The strain, divergence, grad k . grad omega, nu and wall distance are frozen per half step; F1, F2 and mu_t are evaluated at each stage's k and omega. Near walls beta omega dt is O(1) or larger, so an explicit source is not an option.
     - Measured (check 2): positive and second order (error 2.2e-5 at beta2 omega dt = 0.025). At beta2 omega dt = 300 it stays positive but is not accurate (k 2.9 against an exact 8.3 after 10 steps), so a run whose source is that stiff in a region that matters reports it.
   - **Walls and boundaries.**
     - No-slip walls: k = 0, mu_t = 0 and omega = 10 * 6 nu / (beta1 d1^2) as the wall-face value, with d1 the wall distance of the adjacent cell centroid (declared; sensitivity checked). The wall treatment integrates to the wall and needs y+ about 1, reported per run. No wall functions.
     - The k flux through a no-slip wall face stays in the cell as heat, so the energy through the wall is the conduction alone. In the continuum it is zero (k ~ y^2 at the wall); on a grid it is the discretisation's, and this keeps it out of the wall heat ledger.
     - Slip walls pass no k or omega flux and are not walls for the wall distance; their normal stress includes the cell's -2/3 rho k.
     - Wall distance: exact distance from each centroid to the no-slip wall segments.
     - Supplies and the nozzle inlet carry a declared turbulence intensity and viscosity ratio, inputs with provenance.
     - The ambient fill and backflow carry declared ambient values, with the viscosity ratio in TMR's freestream range 1e-5 to 1e-2.
   - **Mesh.** A near-wall resolution of y+ about 1 needs rings clustered at the wall. `Definition::radialStretching` (a tanh distribution of the ring fractions; zero keeps equal rings bit for bit) comes first, verified on the existing transport and core checks.
   - **Verification, with criteria fixed in the test headers before any run.**
     1. Wall distance exact on straight and conical walls.
     2. Decaying homogeneous turbulence (at rest, slip walls, so F1 = 0) against the exact k(t) and omega(t).
     3. Truncation order of the turbulent transport operator (stress, heat, species, k and omega fluxes) on a straight slip-wall duct, where F1 = F2 = 0 and mu_t = rho k / omega, as in `transport_verification`; radial momentum is judged against the laminar order on the same field.
     4. Fully developed turbulent pipe flow (a short duct with transmissive ends, which keep an axially uniform state uniform) driven by a body force, against an independent 1-D solution of the same SST-2003 equations on the same rings (`tools/sst_pipe_1d.py`, `docs/evidence/SST_REFERENCE.md`). Agreement measures the implementation; the grid convergence of c_f and bulk velocity is the model's and is reported. Budgets with turbulence on.
     - TMR's flat-plate data are for SST-V, not SST-2003, so they are a cross-check only; validation against pipe DNS and the TMR axisymmetric subsonic jet (SST-Vm results, PIV data) comes after verification.
8. **Walls.** Inviscid slip first (declared: no boundary-layer loss). Then no-slip under the turbulence closure with a declared wall temperature or heat-flux condition; a regeneratively cooled wall's temperature is an input or a bracket.
9. **Measurements.** Vacuum thrust from the device-thrust ledger with zero ambient; Isp from thrust over the total supplied mass flow (every propellant stream, liquid included, and any igniter flow); c* = Pc A_t / mdot; chamber pressure at the injector face and at the measured tap location; each loss as the delta between two runs of the same engine.
10. **Verification toward RL10, each against CEA:** premixed reactants with local equilibrium give CEA shifting c* and Isp after accounting for finite chamber area and divergence; finite-rate kinetics lands between shifting and frozen; the "valves open, ignite" start reaches the same end state as a partly developed start; grid and time-step refinement of the settled observables.
11. **RL10A-3-3A.** Geometry from NASA TM-107318 (area against axial station); operating inputs (propellant flows or O/F, feed temperatures) from the same or primary sources. LOX enters as liquid through step 5. The predicted vacuum Isp and thrust, with brackets for declared missing inputs, are committed to git before any measured performance is read. Compare, then explain the gap with the loss ledger (divergence, kinetics, boundary layer, mixing, atomization and evaporation). Missing inputs get declared brackets, never tuned values.
12. **Engines as data.** An engine definition holds its contour (a table, or a parametric shape from published dimensions), propellants and feed states, flows or O/F and chamber pressure, ambient pressure, and a provenance label on each field (measured, derived or assumed). A headless batch runner takes a database of such definitions and reports Isp, thrust and c* with their trends against CEA and measurements. There are no per-engine efficiency factors. A fast ideal-rocket tier (the engine's own thermochemistry, as in the CEA check) and the full time-evolving tier read the same definition.
13. **Output.** Ben wants the data representation before the UI: each run writes its time history (CSV) and field frames, with contact sheets and a single exported video per run.
14. **Laboratory sandbox app, first draft** (Ben; after the chemical milestone). "Draw and run": edit the cross-section (wall contour points), place injectors (supply faces and liquid injection sites), run, and watch the fields evolve. "Probe and plot": click anywhere for time histories of the local fields, with thrust and Isp traces for the whole device. It drives the same engine and run records as the headless runner.
15. **Cost.** Reaction integration dominates the step. Measure the cost per cell-step for each chemistry mode, use the worker threads, and run grid convergence on the Windows GPU PC or long Mac runs.

### Dated plan to January

Ben, 4 October 2026: all four regimes work by January. Weeks start Monday. Each milestone exits only with its verification evidence committed and the earlier regression cases still passing; a slipped milestone is reported with re-estimated dates at once, so Ben can choose between the date and the scope of the later regimes. Acceptance standards do not move.

| Dates | Milestone | Done means |
|---|---|---|
| 5–11 Oct | **C1** Time-evolving chamber and nozzle | Contour from data; gas supply faces with valve ramps; igniter; local-equilibrium option. Premixed gaseous H2/O2 "valves open, ignite" run settles; its end state matches CEA after the declared corrections |
| 12–18 Oct | **C2** Transport, walls, mixing closure | Mixture-averaged molecular transport verified against Cantera and FreeFlame; SST-2003 URANS and PaSR; no-slip and thermal walls; interleaved rings mix through the closure |
| 19–25 Oct | **C3** Liquid LOX | Parcels, breakup, evaporation and two-way coupling verified (step 5); a subcritical LOX/GH2 component case |
| 26 Oct–1 Nov | **C4** RL10A-3-3A | Inputs extracted, prediction pre-registered in git, run, comparison table and loss ledger. The chemical milestone |
| 2–8 Nov | **App draft** | "Draw and run" and "Probe and plot" on the chemical engine |
| 9–29 Nov | **M** Magnetic nozzle | The parked spike rejoins the shared engine: resistive MHD with constrained transport on the body-fitted mesh, applied coil fields and a declared vacuum-region treatment, field stresses in the thrust ledger. The zero-field gas cases are unchanged; MHD Riemann and magnetic-nozzle limit checks pass; a published magnetic-nozzle measurement is compared under the blind protocol |
| 30 Nov–13 Dec | **P** ICF post-burn pulse | A declared supplied pulse (mass, species, energy partition; from a published source, or labelled synthetic) expands through a pulsed magnetic nozzle; integrated impulse per pulse and coupling efficiency with the coil energy in the budget; compared with a published pulsed-nozzle result |
| 14–27 Dec | **A** Antimatter | Annihilation products from published yields as particle populations in the shared container: relativistic pusher in the coil field, decay, escape and deposition into the medium, thrust from escaping momentum plus coil reaction; orbit checks and comparison with the published Geant4 antimatter-nozzle study |
| 28 Dec–3 Jan | **All four** | One regression suite runs all four regimes in the same engine and app; buffer |

Lean choices that keep these dates possible: the energetic products in A use tabulated published yields and a tested pusher before Geant4 is adopted (Geant4 enters when material interactions of the products matter); the pulse in P is a supplied state, not a simulated implosion; plasma closures in M and P are the simplest that pass their applicability checks in the chosen cases. None of these is a claim about devices outside the checked cases.

## 1. The application we are building

One native desktop application owns an open experiment: its geometry, physical assumptions, controls, evolving calculation, measurements, and saved history. The person operates and observes that experiment in place. Starting a calculation does not submit a report or leave the editor for a separate analysis product.

Confirmed requirements:

- One shared simulation engine; added capabilities apply across compatible experiments, with no application-specific extensions or separately maintained physics engines.
- Mac, Windows, and Linux each run the app and simulation locally.
- Feed, heating, and current can change during operation when the selected equipment model supports them.
- Geometry changes require pausing and starting a new calculation with an explicit initial state.
- Viewing and controls remain responsive; physical evolution may take minutes.
- Primary development/use is on Ben's Mac, with a Windows PC for larger local runs; the published product must work across hardware classes. Plan for million-cell meshes and measure the target of minutes for short runs.

Use a single executable with an in-process simulation library on a worker thread. Parallel numerical work may use additional CPU threads. Keep a headless entry point for automated verification and later sweeps, using the same engine. Networking, remote execution, MPI distribution, and a browser runtime are not initial dependencies.

## 2. Dependency decisions

| Dependency | Responsibility | Adoption |
|---|---|---|
| C++20, CMake, CTest | Application/core language, build, test orchestration | Baseline; verify compiler/library compatibility in the first build |
| Qt 6 Widgets | Desktop windows, controls, object/property views, input, undo, file dialogs | Selected |
| VTK C++ | Field rendering, picking, contours, vectors, scientific plots, optional revolved view | Selected; Qt integration prototype first |
| Cantera C++ | Supported gas thermodynamics, reaction rates, molecular transport | Selected for chemical physics |
| SUNDIALS | CVODE/BDF for initial stiff chemistry; ARKODE if tighter IMEX coupling is justified | Selected library; candidate methods require coupling tests |
| AMReX | Grid storage, parallel loops, geometry facilities, applicable field operators | Preferred, conditional on the combined RZ/curved-wall/field experiment |
| HDF5 | Numerical history and complete restart data | Selected; dedicated persistence owner |
| Geant4 | Applicable annihilation/product/decay/transport models | Introduced into the shared engine with energetic-product physics |

Python with NumPy/SciPy/Matplotlib supports independent calculations and development analysis. CEA, PlasmaPy, and selected atomic datasets support research and comparisons. They are not required user-facing application processes. Do not add a second rendering engine, general plugin loader, or custom numerical-library replacement without a demonstrated need.

Use Qt Widgets rather than Qt Quick for the initial workbench: desktop selection, property editing, keyboard interaction, and docked scientific views are the immediate needs. Widgets supports a direct VTK integration without introducing QML/render-thread ownership as another application layer. This is an engineering judgment, not a claim that Qt Quick is incapable. [Qt main-window facilities](https://doc.qt.io/qt-6/qmainwindow.html), [VTK widget](https://vtk.org/doc/nightly/html/classQVTKOpenGLNativeWidget.html), [VTK Quick ownership constraints](https://vtk.org/doc/nightly/html/classQQuickVTKItem.html).

Choose and pin released dependency versions as a tested set during the first build; do not use floating development branches. Record source revisions, build options, compiler, platform, and required data files. Begin with CMake package discovery and reproducible dependency build instructions. Decide packaging from that working set; users eventually install an app without assembling a Python or scientific-library environment themselves.

## 3. Ownership and module boundaries

Proposed source layout; create modules when their first behavior is implemented:

```text
app/             Qt shell, tools, property panels, session controller
core/            experiment types, units, commands, session state
geometry/        axisymmetric shapes, physical regions, mesh construction
physics/         equations, material/reaction/transport models, equipment
numerics/        spatial operators, integration, particles, conservation
adapters/        Cantera, SUNDIALS, AMReX, Geant4 integrations
visualization/   VTK presentation, picking, probes, display data conversion
persistence/     experiment serialization, HDF5 history and checkpoints
tests/           focused numerical, interaction, restart, portability tests
examples/        versioned, explained experiments and reference inputs
```

Dependencies point inward: the app calls core/session APIs; numerical execution uses physical models and adapters; visualization consumes published data. Physics and numerics do not include Qt widgets or VTK objects. Keep physical and numerical types explicit; avoid inventing a universal graph of arbitrary physics plugins.

Modules are internal organization, not user-installed extensions. Build the evolution from physical state, materials, reactions and equipment, not a top-level `chemical`/`fusion`/`antimatter` solver switch. Avoid computing irrelevant processes when a physical population is absent, but make implemented capabilities available wherever their assumptions hold. Different closures and integration methods can operate within this one engine; do not force one formula to cover incompatible regimes merely to make the source code look uniform.

A capability is integrated only when its exchanges enter the shared time advancement and budgets, its state is saved/restarted, its results are inspectable, and existing reference cases still pass within declared tolerances. Test a limiting case and a coupled case: adding magnetic physics must preserve the zero-field gas case; adding particle deposition must heat the existing medium rather than a private replacement for it. Application presets contain ordinary experiment data and cannot introduce hidden physical laws.

| State | Owner and rule |
|---|---|
| `ExperimentDefinition` | Editable scene, equipment, initial conditions, model/data selection, numerical settings; versioned when a calculation starts |
| `SimulationState` | Worker alone mutates fields, closure variables such as k/omega, particles, equipment state, time, integrator history, random state, and exchange totals |
| `ControlEvent` | Requested physical input change; worker records its accepted value and effective simulation time |
| `DisplaySnapshot` | Immutable, bounded representation of one accepted state, with time, experiment revision, applied-event sequence, fields and measurements |
| `Checkpoint` | Complete restartable state plus its definition, model/data identities and control history; never confused with a display snapshot |

Use SI internally and explicit dimensions at public interfaces. Stable identifiers connect visible objects, mesh regions, equipment, probes, events, and errors. The scene definition is physical geometry; a VTK actor or grid cell number is not an authoritative physical object.

Provide versioned model descriptions containing evolved variables, closure/equation variant, energy convention, data identities, coefficient defaults, applicability checks and reference cases. A supported configuration selects a compatible set. The interface explains the assumptions and supplies suitable defaults; ordinary users should not need to assemble numerical methods. Expert changes create an explicit experiment revision.

## 4. Interaction and execution

The main thread owns Qt widgets and their VTK render objects. A worker owns the engine and checks a thread-safe command queue. Long numerical routines never run inside a UI callback. Qt's [threading rules](https://doc.qt.io/qt-6/threads-qobject.html) govern object affinity and cross-thread delivery.

Execution states are `Ready`, `Running`, `PauseRequested`, `Paused`, `Completed`, and `Failed`. Editing an incompatible definition marks it as needing initialization. An edited geometry retains the previous calculation for inspection/comparison but cannot resume its state on the changed mesh.

Live-control sequence:

1. The UI validates dimensions/range and immediately displays the requested setting as pending.
2. The worker accepts the command at a safe step boundary, or rejects it with a physical explanation. It records the value and actual simulation time applied.
3. The next advance uses that input, including any specified equipment response or ramp. Prescribed current is not a substitute for induced-field energy accounting.
4. A published snapshot identifies the applied event. The UI then shows the setting as active and places an event marker on the history.

Commands requesting the next safe time and commands scheduled at a specified physical time are distinct. Schedule discontinuities so a step does not cross them unnoticed. Log actual accepted events for replay; wall-clock mouse timing is not reproducible physics. Coalesce unsent slider motion, but never erase already-applied control history. Undoing a live change issues another physical event; it does not erase the medium's response.

Pause takes effect at a consistent accepted state. Poll cancellation within bounded work units; reject incomplete trial advances before publishing. Long library calls may delay safe pause, so report pending status immediately and measure the delay. Avoid relying on a worker event loop that cannot process commands during its own long-running solve.

A rejected advance must also roll back closure/equipment state, random streams, event application and exchange tallies. Never redraw stochastic trials merely to select a more favorable accepted outcome. Record control discontinuities once, and reset multistep integration history when required. Model/data changes and incompatible numerical changes restart or branch the calculation; they are not live controls.

The renderer receives the newest available display snapshot through a bounded mailbox. Dropping obsolete display frames is allowed; dropping physical integration steps or accepted control events is not. No view holds a pointer to mutable solver memory. An additional persistence queue owns saved output; its backpressure may slow simulation but must not freeze interaction or silently discard requested checkpoints.

## 5. Physics representation and accounting

Use finite-volume evolution on the axisymmetric meridional domain. Geometry supplies physical ring volumes, face areas, wall normals, and consistent source metrics. Treat the axis with symmetry/parity and regularity conditions. Retain azimuthal vector components when required, without implying resolved azimuthal variation.

Start with compressible reacting gas: species mass, momentum, and a thermodynamically consistent energy variable. Later plasma models add their required thermal and electromagnetic variables; particles carry position, momentum, species and statistical weight. Do not force all materials into the gas equation of state or equate an energetic particle population with bulk temperature.

The initial plasma candidate is a conducting-fluid model in a justified collisional regime, allowing separate electron/ion energies. Its field formulation, closures and low-density treatment must pass the research decision before becoming a supported mode. A future hybrid description can share geometry, controls, accounting and observation without pretending to share every evolution equation.

The electromagnetic calculation may extend beyond the material domain. Distinguish external coil fields from plasma-generated response using a consistent total-field formulation and force/work budget. A changing coil current requires the associated induction and external power; changing only B in a particle pusher is not adequate. Initial prescribed-field test problems must be labeled as neglecting backreaction. Empty regions do not become conducting fluid through an unreported density floor. Field/vacuum boundaries and any fluid/kinetic handoff must be specified before magnetic-device predictions.

Axisymmetry applies to the ensemble and fields, not to individual particle motion confined to a plane. Retain three momentum components and the required trajectory geometry; normalize particle weights and deposition to physical ring volumes. Test trajectories crossing the axis and radial mesh boundaries.

Each model declares required state and its valid inputs. Assemble a supported model combination at initialization; reject incompatible combinations with an explanation. Presets supply arrangements and parameters, not private engine-specific laws.

Recheck model applicability as the state evolves, using available quantities such as collisionality, magnetization, temperature/composition ranges and optical depth. Distinguish a computable out-of-domain extrapolation from a numerical failure. Preserve the last inspectable state and identify affected observables. Passing these checks does not prove that all missing physics is insignificant.

Boundary objects include injector hardware, supplies, walls/cooling, igniters, coils/drivers, and environment. They may be distributed within the domain and may have response state. Their interfaces return physical exchanges or constraints and their derivatives where needed; they do not overwrite arbitrary solution fields. Combustion and internal transport remain part of the medium calculation.

Translate equipment into well-posed mathematical conditions: do not independently prescribe incompatible pressure, velocity and mass flow at an inlet. Supply reservoir/orifice response or a declared imposed-flow assumption, composition, enthalpy and inlet turbulence as needed. Treat backflow, ambient pressure and wave exit consistently with the equations. An axisymmetric inlet represents an annulus or justified averaged injector distribution; it does not resolve individual circumferential injector holes. Walls declare thermal response, slip/no-slip and relevant particle/sheath exchanges.

Each transfer has one numerical owner. Species formation energy, ionization/binding energy, particle rest energy and field energy require explicit reference conventions before coupling models. Cantera reaction rates must not be combined with an additional heat source that counts the same chemical release again. Product deposition removes the corresponding particle energy/momentum. External power, stored energy and escaping fluxes close the experiment's budget.

The particle/medium/field coupling must also conserve charge consistently with the selected field equations. Specify whether energetic-product charge/current affects those fields or is neglected under a tested dilute-population approximation. Energy/momentum deposition alone is not a complete electromagnetic coupling. Include cylindrical tensor/curvature terms in averaged stresses and field operators even though azimuthal derivatives vanish.

Compute thrust from a declared control volume and device force convention, accounting for pressure, material momentum, field stresses and energetic products as applicable. Cross-check transfers to walls/coils; do not sum duplicate force estimates. Record the mass-flow denominator used for specific impulse, including which supplied streams it counts. Pulses need integrated impulse and consumed mass, not a misleading instantaneous ratio near zero flow.

Radiation emission/absorption and unresolved kinetic energy have explicit owners in that budget. Start with optically thin thermal losses only in an assessed regime, and transport energetic photons through the selected particle treatment when applicable. A general optically thick radiation model is unresolved; do not hide that gap behind an unrestricted cooling coefficient. Restart files include all state required by these models.

## 6. Spatial and temporal methods

Begin with a fixed 2D mesh. AMReX adoption depends on proving the required RZ geometry, curved-wall discretization, small-cell treatment, near-axis behavior and field operators together. Adopt no AMR or embedded-boundary capability solely because it appears in a library feature list. A failed experiment requires revising this decision before building applications around it.

Use the following starting algorithms for focused prototypes. They are candidates within the stated regimes, not an obligation to implement every row before a useful app exists:

| Process | Starting algorithm and acceptance condition |
|---|---|
| Gas advection/shocks | Second-order limited piecewise-linear reconstruction, HLLC flux, SSP-RK2 explicit transport; test an HLLE fallback with suitable wave speeds where admissibility/shock behavior requires it |
| Chemical stiffness | CVODE variable-order BDF with a consistent species/energy system and scale-aware tolerances; verify elemental conservation and thermodynamic recovery |
| Gas turbulence | Favre URANS with documented SST-2003 equations; evolve rho-k and rho-omega, declare turbulent pressure/energy conventions and inlet/wall conditions |
| Turbulent reaction | PaSR plus finite-rate chemistry as the first mixing-controlled candidate; validate mixing/chemical time definitions and held-out reacting profiles |
| Molecular/turbulent diffusion | Conservative face fluxes; explicit on early small cases, implicit when measured diffusion timestep limits require it; conserve net species mass and associated enthalpy transport |
| Collisional plasma transport | Applicable Braginskii-type anisotropic coefficients with distinct electron/ion energies; test diffusion along oblique fields and energy exchange |
| Magnetic evolution | Face-centered magnetic flux with constrained transport is preferred; prove the cylindrical/curved-wall combination before selecting its implementation |
| Energetic trajectories | Use Geant4's supported field integration for segments it owns; if a separate plasma tracker is required, test a relativistic Boris/Vay-class pusher against reference orbits and electric-field work |

HLLC is a gas solver; an ideal-MHD prototype may assess HLLD/HLLE only with an appropriate EOS and field formulation. These choices are not interchangeable with arbitrary two-temperature/nonideal plasma equations. [Athena solver distinctions](https://princetonuniversity.github.io/Athena-Cversion/AthenaDocsUGRiemann.html), [general-EOS considerations](https://arxiv.org/abs/1909.05274), [CVODE methods](https://sundials.readthedocs.io/en/latest/cvode/Mathematics_link.html), [relativistic particle-pusher comparisons](https://warpx.readthedocs.io/en/latest/theory/kinetic_particles.html). Closure rationale and primary references are in [RESEARCH.md](RESEARCH.md#how-unresolved-transport-will-be-approximated).

An admissible update needs reconstruction, flux, source and timestep treatment together; selecting HLLE does not guarantee positivity for every EOS or cut cell. Prefer conservative limiting or step rejection to independently clipping species or pressure. Log any remaining repair and quantify its contribution. Use the same face exchanges in evolution and diagnostic balances. Numerical diffusion is not an accepted turbulent-transport model.

For each closure implementation, record one compact specification alongside its tests: exact equations/variant, units, state, boundary data, reference coefficients, validity checks, discretization and conserved transfers. Reproduce a reference calculation before changing coefficients. Keep fitted parameters and the data used to fit them separate from the held-out comparisons. These specifications become the implementation authority for each supported model; library names and this algorithm table alone are insufficient.

Use Cantera for local constitutive/rate evaluations with independent contexts per executing thread. The CRUCIBLE integrator owns global time advancement; do not have both an independent Cantera reactor network and SUNDIALS advance the same state. Begin with documented symmetric reaction/transport splitting and stiff reaction integration; measure splitting error and replace it with tighter coupling if the target regimes require it. Nominal scheme order is not evidence of observed accuracy.

SUNDIALS supplies integration machinery. Our code supplies state mapping, residuals, tolerances, Jacobian/preconditioner strategy where needed, admissibility handling and acceptance criteria. Magnetic evolution additionally requires a documented divergence-control strategy compatible with the geometry. Do not claim magnetic support from a field visualization alone.

Coupling acceptance checks all relevant timescales: wave propagation, reaction, diffusion, electron/ion exchange, particle interaction and equipment response. An error-controlled chemical substep does not control the error of the complete split calculation. Strong feedback may require iterations or an IMEX treatment; establish convergence rather than doing one ordered pass and assuming the loop is closed. Guiding-center particle reductions, reduced chemistry and precomputed response tables are later optimizations with separate validity/convergence tests.

For energetic products, batch interaction/transport work over a tested coupling interval, return conservative exchanges, and update the receiving medium. Refine both interval and particle count. Geant4 may own applicable transport segments, but its material models and external fields must match the experiment; it is not automatically a hot-plasma closure. No second tracker may simultaneously own the same particle segment.

## 7. Rendering, saving and performance

The main view is the physical cross-section with selectable geometry/equipment, field overlays, vectors, units, probes, and a timeline. Provide contextual property editing, clear calculated-versus-prescribed quantities, and simple saved comparisons with matching scales. Add a revolved 3D view only as another presentation of the axisymmetric state. Invalid values must remain visibly invalid rather than being smoothed into attractive contours.

Label modeled mean fields, sampled particle trajectories and prescribed inputs correctly. Offer transport/reaction/escape budget views and comparison across uncertain parameters. Keep model maturity (verified, experimentally compared) separate from applicability to the current state. Do not synthesize turbulent eddies for visual realism.

Use `QVTKOpenGLNativeWidget` with its supported render-window type and required surface setup. VTK owns presentation data; expensive derived-field preparation happens outside UI callbacks. Initial plots can use VTK's chart facilities, avoiding a second plotting dependency.

Live viewing uses memory snapshots, not files. Keep a versioned experiment manifest plus HDF5 numerical data in the saved project. One persistence worker owns HDF5 access. Select dataset layout/chunks for expected time/field reads and benchmark them; [HDF5 chunking guidance](https://portal.hdfgroup.org/documentation/hdf5/latest/hdf5_chunking.html) explains why access patterns matter. Commit checkpoints as complete records and recover the last complete record after interruption. History scrubbing is playback; resuming an earlier state requires a checkpoint and a new continuation branch.

Portable restart means agreement to declared tolerances, not serializing opaque library pointers or promising identical internal integrator histories. Save physical/closure state, random state, control schedules and required numerical state; reinitialize library solvers through supported APIs and test split-run versus uninterrupted results. Version the saved schema and preserve provenance during migration. Curated examples include redistributable inputs or explicit acquisition instructions; files on a developer's machine are not an implicit dependency.

Initial interaction targets, to measure on named hardware: ordinary control acknowledgement within 100 ms, useful manipulation of the displayed view around 30 frames/s, bounded memory over a long session, and progress/safe-pause status even when a step is slow. These are design targets, not measured performance promises. Physical time and wall time are displayed separately; never animate fabricated intermediate physics to conceal a slow solver.

CPU execution is the baseline on all platforms. Leave compute headroom for the UI and avoid nested thread oversubscription. Profile reaction evaluation, field solves, transport, particles, copying and rendering separately before optimizing. Optional CUDA acceleration follows evidence; do not promise Apple GPU acceleration through AMReX. Host-only Cantera calls do not become GPU kernels automatically. Read back only requested display data at a bounded rate.

Performance is an early architecture check. Benchmark roughly 100,000, 1 million and several million cells as memory permits, reporting actual hardware, species/state count, timestep, simulated interval and accuracy settings. Use contiguous field storage and bounded scratch/snapshot buffers; avoid per-cell heap objects, dynamic dispatch in hot loops, and one heavyweight library/integrator instance per cell. Reuse execution contexts and evaluate batching. One million cells holding 20 float64 fields already require about 160 MB for one state copy, before mesh data, ghosts, integration stages, chemistry workspaces and display buffers. Track peak memory as well as throughput.

A transport-only throughput result cannot establish reacting-flow speed. Benchmark the coupled workload early, including stiff chemistry and any implicit solves. Build optimized configurations for timing. Adopt validated reduced chemistry, batching, adaptive refinement or accelerator kernels when profiling justifies them; preserve the physical model and verify approximation errors. If the CPU path misses the Mac target, reassess the compute strategy explicitly rather than assuming a Windows GPU will satisfy local-Mac performance.

## 8. Implementation sequence and completion evidence

1. **Native and numerical foundation:** build the same small Qt/VTK app on Mac, Windows and Linux; run a worker calculation with real field snapshots; exercise controls, pause, save/reopen and shutdown. Pin the working dependencies and verify RZ gas flow. Audit experimental inputs before constructing a large chemical example. Synthetic fields are explicitly development fixtures.
2. **Chemical validation milestone:** operate a reacting chamber/nozzle with gas supplies and liquid LOX, change feed live, inspect the response, restart a geometry revision, and compare. Verify chemistry, conservation and numerical sensitivity. Assess SST/PaSR against appropriate flow/flame data and compare declared chemical-engine observables against real measurements. A laminar prototype, CEA match or wall-heat comparison alone does not establish all turbulent-engine performance quantities.
3. **Magnetic and pulse physics in the shared engine:** resolve the plasma regime and field/vacuum treatment, then demonstrate magnetic expansion and justified post-burn pulse input using the established scene, controls, budgets and observation tools. Preserve the chemical regression cases.
4. **Energetic reactions and products in the shared engine:** add selected antimatter interactions, applicable fusion reactions, particle transport and deposition as their model/data evidence permits. Test existing and newly coupled cases. Remaining burning inside a supplied pulse is included when significant; the upstream implosion remains excluded.
5. **Research use and delivery:** test with unfamiliar users, package desktop builds, and support substantive investigations selected by Ben. Usability, packaging checks and validation develop throughout the earlier stages, not only here. MCF configuration searches remain a stretch dependent on suitable transport and stability evidence.

Appropriate checks include conservation and convergence, command ordering and rejected steps, stable snapshot ownership, complete restart after interruption, repeated start/pause/close, cross-platform numerical agreement to declared tolerances, and measured interaction under load. A passing UI demonstration does not validate physics; a passing solver test does not establish usability.

Before broad implementation, produce a working desktop/worker loop, a measured RZ flow prototype, and an audited chemical-validation case with defined observables and a cost estimate. Also perform a small geometry/field compatibility check to prevent an unsuitable spatial-library commitment; it must not become a competing plasma-development project. Chemical validation is the first delivery priority. The desired horizon is a few months; weekly time availability and measured development throughput are not yet known.

The provisional 12–16-week envelope of 17 September is replaced by the *Dated plan to January* at the top of this document (Ben, 4 October 2026). Failure to reproduce the chemical evidence changes the schedule before it changes acceptance standards. MCF optimization remains a stretch beyond January. Avoid polishing a separate chemical-only product or adding independent application solvers to meet dates.

Agents should implement one demonstrable increment at a time, update the relevant decision with evidence, and preserve a working app. Do not scaffold every proposed module, invent missing datasets, tune a case to desired performance, revive archived gates, or claim unimplemented support. [RESEARCH.md section 7](RESEARCH.md#7-lessons-from-the-pre-pivot-attempt) identifies reusable evidence and code; the archive stays historical.

## 9. Coding and review workflow

Use one active implementation owner per change, with bounded tasks and an independent review of consequential work. Given Ben's current subscriptions, the recommended operating pattern is Claude Code for most implementation and Codex for selected architecture, numerical-physics and regression reviews. This is a budget allocation, not an assumption that either model can certify the other's physics. Roles can be swapped without changing project authority.

For each increment: state the behavior and independent acceptance evidence, implement it on `main` (Ben works solo; use a branch only for a risky experiment and merge or abandon it promptly), run the relevant checks, inspect the actual result, review material risks, then commit the coherent change. Avoid concurrent edits to the same checkout; parallel work needs separate worktrees and clear file ownership. Review a stable diff or commit, not files changing underneath the reviewer. Ben retains direction and consequential model choices; agents handle routine implementation within the agreed task.

Use fresh implementation/review conversations at natural task boundaries. Carry context through the active documents, code, tests and a short handoff containing current commit, completed behavior, checks run, unresolved issue and next task. Do not repeatedly ingest the archive or require both agents to reread the full conversation. Evidence should include executed test commands, numerical comparisons and actual UI inspection where appropriate. Agreement between models does not replace an independent reference or measurement.
