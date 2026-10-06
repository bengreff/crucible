# CRUCIBLE — Technical plan

17 September 2026; current-state section, chemical plan and dated plan to January updated 5 October 2026. Target architecture; most capabilities remain planned. Implemented behavior and its verification are indexed in [docs/evidence/README.md](docs/evidence/README.md); the live state is in [docs/SESSION_HANDOFF.md](docs/SESSION_HANDOFF.md). Read [VISION_SCOPE.md](VISION_SCOPE.md) for purpose and scope and [RESEARCH.md](RESEARCH.md) for evidence and unresolved science. Current user instructions take precedence. Architecture changes must preserve those physical and human requirements.

**How to read decisions:** product requirements are fixed by the vision and Ben's instructions; selected architecture is the current implementation direction; algorithm candidates require the stated tests. Do not promote a candidate to supported physics by implementing its interface. Keep one authoritative home for each decision: purpose in the vision, model evidence in research, implementation here, historical lessons in RESEARCH.md section 7.

## Current state and the plan to January (4 October 2026)

**Implemented and verified** (evidence index: `docs/evidence/README.md`): a fixed body-fitted axisymmetric finite-volume mesh with uniform axial spacing; HLLC with limited primitive reconstruction and SSPRK2; a thermally perfect multi-species mixture (NASA-7 data from the mechanism) with species carried by the mass flux; Cantera thermochemistry behind `adapters/`; per-cell stiff kinetics (CVODES BDF over Cantera rates) joined to the flow by symmetric Strang splitting; a device-thrust momentum ledger; a Qt/VTK app operating a gas nozzle live (Mac only so far). Release is the default build. C1 (4 October): wall contours from data, gas supply faces with valve ramps, an igniter, local-equilibrium and frozen chemistry modes, and a "valves open, ignite" chamber case checked against the ideal rocket (`docs/evidence/CHAMBER_C1.md`; on 128x24 c* is within 0.09% and converging at order 2.0, while the vacuum Isp is within 0.22% but not yet converging under refinement). C2 under way: mixture-averaged transport properties match Cantera; the transport operator (viscous stress, conduction, diffusion; no-slip, slip, isothermal and adiabatic walls) passes its solution and budget checks, and two of its four truncation-order criteria still fail (`docs/evidence/TRANSPORT_C2.md`). Wall-clustered rings and SST-2003 are in the engine, with turbulent Chamber supplies. The wall distance, homogeneous decay, turbulent operator and supply checks pass. The fully developed pipe matches the independent 1-D reference (with the SSTs production of 5 October: u_b and c_f within 1.0e-4 on nr 16 and 2.1e-5 on nr 32; nr 64 was 3.2e-5 before the change, on backhouse, and was not rerun), and fails only its 1e-12 mass-budget criterion, at 6e-12 to 6e-11; the cause is inferred to be rounding of the stored cell densities over millions of steps (`docs/evidence/TURBULENCE_C2.md`). A freely propagating H2/air flame through the full engine matches Cantera's free-flame speed within 0.28% (S_c) and 0.04% (S_d) at 10 um (`docs/evidence/FLAME_C2.md`). The PaSR closure passes its criteria 1 to 3, and criterion 4 (a turbulent reacting chamber from a cold start, with and without the closure) on 32x6; its smoke runs found and fixed two k-omega faults (`docs/evidence/PASR_C2.md`). AMReX is not adopted; the fixed custom mesh has been sufficient so far, and the decision stays conditional (section 6).

**Not yet implemented:** turbulent inflow on the nozzle inlet (Chamber supplies have it), the PaSR closure on a non-premixed flame (the criterion 4 chamber is premixed), interleaved fuel and oxidizer rings, liquid injection, breakup and evaporation, persistence, experiment editing, other platforms.

**Direction** (Ben, 4 October 2026; recorded in VISION_SCOPE *Changes since the pivot*): one time-evolving reacting engine; shifting equilibrium as the validated baseline produced by that engine; liquid LOX injection with breakup and evaporation in the shared engine before RL10 ("Add evaporation first"; no pre-vaporized RL10); every chemical run starts "valves open, ignite"; RL10A-3-3A next; then CEA sweeps for trends; then a database of hundreds of engines. By January all four regimes work in the same engine: chemical, magnetic nozzle, pulsed fusion (ICF post-burn pulse) and antimatter. Chemical stays lean and correct; the parked MHD spike returns right after it; the first draft of the laboratory sandbox app follows the chemical milestone. The dated order is in *Dated plan to January* below. The steps here are the selected architecture for the chemical work. Each lands with its verification.

1. **Always a time-dependent calculation.** Every chemical experiment marches the unsteady reacting equations from a declared initial state. No steady-state solver, no pseudo-time acceleration that alters the transient. Steady performance is the average over a late, settled window of the history; the settling test and window are recorded with the result.
   - **Start: "valves open, ignite"** (Ben). The chamber and nozzle start filled with ambient gas at rest (declared composition, pressure, temperature). The propellant valves open on a declared ramp schedule (valve travel is an input, not modeled equipment), the igniter fires, the flame establishes, the chamber pressurizes, and the flow settles.
   - A run that cannot afford the full start says so with the measured cost; a declared partly developed initial state that still evolves in time is then a development fallback, never the reported RL10 case.
2. **Chemistry.** Finite-rate kinetics in the existing Strang slot is the production model. Local equilibrium per cell (constant-u,v equilibrium of the cell's own elements, in the same slot) is its fast-chemistry limit and is used to verify against CEA shifting equilibrium; frozen chemistry is the other limit. Use a hydrogen/oxygen mechanism validated to rocket chamber pressures (the Burke et al. 2012 or Li et al. 2004 class) and record its validity range; Cantera's `h2o2.yaml` is a GRI-3.0 subset suitable for verification only. Since 5 October 2026 (*Lightweight engine*), the runtime chemistry in both modes comes from tables built offline with Cantera: Table A for local equilibrium, then Table B for finite rate. Direct Cantera equilibrium and CVODES build the tables and remain the reference checks.
3. **Geometry as data.** The wall is r(z) from a published contour; for an area-against-station table r = r_t sqrt(A/A*), so no shape is assumed. Uniform axial spacing first. Since 5 October 2026, resolution is adaptive: block-structured refinement in the mapping's index space (*Lightweight engine*).
4. **Gas supplies.** The injector face is a boundary object. Each radial face carries a stream with imposed mass flux, total enthalpy and composition; static pressure comes from the interior (one outgoing characteristic), so the inlet is well posed; backflow is flagged, not silently reversed. A coaxial-element injector becomes annular averages: either premixed at the face (a declared limit with no mixing loss, used for verification) or interleaved oxidizer and fuel rings, whose mixing is left to the transport closure. A valve ramp scales each stream's mass flux in time.
5. **Liquid propellant: injection, breakup and evaporation** (Ben, "Add evaporation first"; required before RL10). Liquid enters as stochastic Lagrangian parcels in the RZ domain: each parcel carries position, three velocity components, drop diameter, temperature and a statistical weight, and exchanges mass, momentum, energy and species with the gas through sources inside the shared step, normalised to physical ring volumes. Primary atomization of a coaxial element is not resolved in 2-D axisymmetry: the initial drop size distribution comes from a published coaxial-injector correlation and is carried as a declared bracket, the leading uncertainty of the model. Secondary breakup uses a KH-RT class model; evaporation uses a film-theory model (Abramzon-Sirignano class) with Cantera gas properties and cited liquid O2 properties. The liquid's enthalpy is the real feed enthalpy, so the latent heat is honoured without a separate correction. Applicability is checked as the state evolves: subcritical evaporation requires chamber pressure below the O2 critical pressure. The same particle container later carries energetic antimatter products with their own laws. Verification: single-drop evaporation against a reference solution, mass and energy budgets closing with parcels present, and a published subcritical LOX/GH2 spray flame as a component case.
6. **Igniter.** A bounded energy deposit (declared energy, volume and duration) or a defined hot-gas supply, counted in the energy budget once.
7. **Unresolved transport, applied uniformly.** Molecular transport first (mixture-averaged, Cantera's fits, verified against Cantera and a laminar flame speed from Cantera's FreeFlame). Then SST-2003 URANS with a declared wall treatment and turbulent Prandtl/Schmidt numbers, then PaSR for turbulence-chemistry interaction. The same closure acts on every stream, the chamber and the nozzle; there is no combustion-only mixing shortcut. Until it exists, results that depend on mixing are not claimed, and premixed runs isolate everything else.
   - **Turbulence is modelled, never resolved.** The engine evolves Reynolds- and Favre-averaged mean fields (URANS). SST-2003 supplies the turbulent stresses and heat and species fluxes, and PaSR supplies the turbulence-chemistry interaction. No eddy is resolved on any grid. Grid refinement converges the model's mean solution; it does not resolve turbulence. The RZ domain could not hold resolved turbulence anyway, because eddies are three-dimensional. The time dependence is that of the mean flow: startup, valve transients and acoustics. LES and DNS are not planned.
   - **SST-2003 as built** (equations from the NASA Turbulence Modeling Resource `sst.html`, extracted 4 October 2026). The model is the standard SST with the 2003 changes and nothing else:
     - production limited to min(P, 10 beta* rho omega k) in both equations;
     - mu_t = rho a1 k / max(a1 omega, S F2) with the strain invariant S = sqrt(2 S_ij S_ij), which includes the hoop strain u_r / r;
     - gamma1 = 5/9, gamma2 = 0.44, CD_komega floor 1e-10;
     - other constants: sigma_k 0.85 / 1.0, sigma_omega 0.5 / 0.856, beta 0.075 / 0.0828, beta* 0.09, a1 0.31, blended by F1.
     - Production in TMR's SSTs form, P = mu_t S^2 in both equations ([NASA TMR](https://tmbwg.github.io/turbmodels/sst.html); adopted 5 October 2026, Ben's default of published over bespoke). The -2/3 rho k delta_ij term stays in the stress for the momentum and energy fluxes. It replaces the exact P = tau_ij du_i/dx_j, whose dilatation part, -(2/3) rho k div in omega's production where the limiter is active, collapsed omega to 1e-138 in the turbulent chamber's divergent nozzle (PASR_C2, criterion 4), and a CRUCIBLE-only fix of that term used for a few hours. Both productions are now non-negative. The real dilatation term no longer enters k (the turbulent pressure's expansion work stays in the energy flux).
   - **Energy accounting** (RESEARCH *Gas turbulence*: account consistently): k is part of the total energy, and the thermodynamic state uses e = E/rho - |u|^2/2 - k.
     - The Reynolds stress includes -2/3 rho k delta_ij, carried in the transport flux; its work u . (-2/3 rho k I) is therefore in the energy flux.
     - The energy flux carries (mu + sigma_k mu_t) grad k.
     - Turbulent heat and species fluxes use Pr_t and Sc_t (declared inputs; 0.9 and 0.7 first, sensitivity reported). The species flux adds mu_t / (rho Sc_t) to every D_km in the existing mixture-averaged form, which with the correction velocity is exactly -(mu_t / Sc_t) grad Y_k.
   - **State and numerics.**
     - rho k and rho omega are two more conserved fields, carried by the mass flux like the species and reconstructed like them.
     - Their sources are integrated in a Strang split inside each step (half step, transport step, half step), with the mean-flow coefficients frozen per half step.
     - The sources are integrated per cell at fixed rho and E (decayed k becomes heat) by MPRK22, the second-order positive modified Patankar Runge-Kutta scheme (Kopecz and Meister, BIT 58, 2018): production explicit, destruction weighted by the new over the old value. The strain, divergence, grad k . grad omega, nu and wall distance are frozen per half step; F1, F2 and mu_t are evaluated at each stage's k and omega. Near walls beta omega dt is O(1) or larger, so an explicit source is not an option.
     - Measured (check 2): positive and second order (error 2.2e-5 at beta2 omega dt = 0.025). At beta2 omega dt = 300 it stays positive but is not accurate (k 2.9 against an exact 8.3 after 10 steps), so a run whose source is that stiff in a region that matters reports it.
     - The k and omega diffusion is not sign-preserving on contoured rings. The tangential part of the face gradient can draw k out of a cell at k = 0, and next to quiescent gas it left rho k at -1.3e-198, which no step reduction cures (check 5, run 1). After each RK stage rho k is therefore set to max(rho k, 0). E is untouched (k is part of it), so the energy budget stays exact. The rho k V added is summed in `Measurements::clippedTurbulentEnergy`; every run reports it against the domain integral of rho k, and a run where it is not negligible says so (check 5: 0 in every judged case).
   - **Walls and boundaries.**
     - No-slip walls: k = 0, mu_t = 0 and omega = 10 * 6 nu / (beta1 d1^2) as the wall-face value, with d1 the wall distance of the adjacent cell centroid (declared; sensitivity checked). The wall treatment integrates to the wall and needs y+ about 1, reported per run. Since 5 October 2026 this wall-resolved treatment is the verification reference, and wall functions are the default (step 17).
     - The k flux through a no-slip wall face stays in the cell as heat, so the energy through the wall is the conduction alone. In the continuum it is zero (k ~ y^2 at the wall); on a grid it is the discretisation's, and this keeps it out of the wall heat ledger.
     - Slip walls pass no k or omega flux and are not walls for the wall distance; their normal stress includes the cell's -2/3 rho k.
     - Wall distance: exact distance from each centroid to the no-slip wall segments.
     - Chamber supplies carry a declared turbulence intensity I and viscosity ratio, inputs with provenance (built, check 5): k = 3/2 (I u)^2 and omega = rho k / (ratio mu) at the face; the stream's total enthalpy satisfies h(T) + u^2/2 + 5/3 k = h0 (k plus the work of the normal stress 2/3 rho k); the momentum flux is g u (1 + I^2) + p. The nozzle inlet will carry the same pair (not built).
     - The ambient fill and backflow carry declared ambient values, with the viscosity ratio in TMR's freestream range 1e-5 to 1e-2. The ambient omega must not be degenerate. With omega near zero, the negative cross-diffusion term -2 (1 - F1) rho sigma_w2 |grad k . grad omega| / omega drives omega through zero once a k front arrives (check 5, run 2: omega 1 1/s failed at 134 us). Default: the Spalart and Rumsey free-stream values (AIAA J 45, 2544, 2007), k = 9e-9 a^2 and omega = 1e-6 rho a^2 / mu (mu_t / mu about 0.009).
   - **Mesh.** A near-wall resolution of y+ about 1 needs rings clustered at the wall. `Definition::radialStretching` (a tanh distribution of the ring fractions; zero keeps equal rings bit for bit) comes first, verified on the existing transport and core checks.
   - **Verification, with criteria fixed in the test headers before any run.**
     1. Wall distance exact on straight and conical walls.
     2. Decaying homogeneous turbulence (at rest, slip walls, so F1 = 0) against the exact k(t) and omega(t).
     3. Truncation order of the turbulent transport operator (stress, heat, species, k and omega fluxes) on a straight slip-wall duct, where F1 = F2 = 0 and mu_t = rho k / omega, as in `transport_verification`; radial momentum is judged against the laminar order on the same field.
     4. Fully developed turbulent pipe flow (a short duct with transmissive ends, which keep an axially uniform state uniform) driven by a body force, against an independent 1-D solution of the same SST-2003 equations on the same rings (`tools/sst_pipe_1d.py`, `docs/evidence/SST_REFERENCE.md`). Agreement measures the implementation; the grid convergence of c_f and bulk velocity is the model's and is reported. Budgets with turbulence on.
     5. Turbulent supplies: the face routine against the energy, characteristic and turbulence relations to round-off; inflow turbulence decaying in plug flow against the exact decay along the streamline; budgets and positivity in a no-slip chamber with two supply rings.
     - TMR's flat-plate data are for SST-V, not SST-2003, so they are a cross-check only; validation against pipe DNS and the TMR axisymmetric subsonic jet (SST-Vm results, PIV data) comes after verification.
   - **PaSR (selected design, 4 October 2026; implemented 5 October, `docs/evidence/PASR_C2.md`: criteria 1, 2 and 3 pass, criterion 2's test-power check after a dated restatement; criterion 4, a turbulent reacting chamber against the same run without the closure, is stated in `tests/pasr_tests.cpp` and not yet run).** One turbulence-chemistry closure, applied in every cell alike, inside the existing Strang chemistry slot. Algorithm reference: the OpenFOAM `PaSR` model (RESEARCH *Turbulent combustion*), with the departures below stated and justified.
     - **Rate.** The mean source is kappa_eff times the finite-rate source at the cell's mean state. CVODES integrates dz/dt = kappa_eff(z) f(z) at fixed rho and e, with z = [T, Y] and f the mechanism's rates, so heat release and species are scaled together.
       - The mixing time and the segregation (below) depend on flow quantities and are frozen per half step, as the SST sources freeze theirs. The chemical time is re-evaluated from z at every right-hand-side call.
       - OpenFOAM instead multiplies the laminar change over the step, (Y_lam(dt) - Y0) / dt, by kappa. When the laminar chemistry reaches equilibrium inside the step, that rate depends on dt. Integrating the scaled equation is the model itself and does not.
       - With kappa_eff = 1 the code path is the laminar integration, bit for bit.
     - **Reacting fraction.** kappa = tau_c / (tau_c + tau_mix) (Chomiak; the OpenFOAM form).
     - **Mixing time.** tau_mix = C_mix sqrt(nu_eff / epsilon), as RESEARCH names for the first test.
       - nu_eff = (mu + mu_t) / rho, and epsilon = beta* k omega, the SST model's own dissipation.
       - C_mix is a declared input. The first runs use 1 (the Kolmogorov time) and report 0.1 as sensitivity; it is never fitted to an engine result.
     - **Chemical time.** Formation-rate time scales (Ferrarotti et al., Energy Fuels 32, 10228, 2018, eq. 4, after Li et al. 2017): tau_c,i = Y_i / |dY_i/dt| over a declared species set S.
       - tau_c is the slowest active one; species with Y_i <= 0 or no net rate are dormant and excluded. If every species of S is dormant, kappa = 1.
       - For hydrogen and oxygen, S = {H2, O2, H2O} (the reactants and the major product, as Ferrarotti's major-species set).
       - Known concern, to examine in validation rather than assume away: in a non-premixed flame the abundant reactant's time is long, so the slowest-species rule can leave kappa near 1 off stoichiometric.
     - **Zero-turbulence limit (a CRUCIBLE choice, not in the reference).**
       - The problem: as k goes to 0, tau_mix grows like k^(-1/2), so the reference form sends kappa to 0 and switches reaction off in laminar gas. That contradicts the verified laminar flame and would stop ignition in quiet gas.
       - Why the reference form must be modified: a Reynolds average of a flow with no fluctuations is the flow itself, so the closure must return the mean-state rate.
       - The modification: the reacting fraction is weighted by the segregation the turbulence produces. kappa_eff = 1 - s (1 - kappa). A fraction s of the cell is segregated and reacts through the PaSR fine structures; the rest reacts at the mean state.
       - How s is computed: from the algebraic scalar variance at production-dissipation equilibrium, with production 2 D_t |grad X|^2 and dissipation 2 (epsilon / k) X''^2. That gives X_i''^2 = (nu_t / Sc_t) |grad X_i|^2 / (beta* omega), using the engine's mole-fraction gradients. Then s_i = min(1, X_i''^2 / (X_i (1 - X_i))) for each i in S with 0 < X_i < 1, and s is the largest.
       - Limits: as k goes to 0, s falls like k and kappa_eff goes to 1. With strong segregation, s = 1 and the form is the reference PaSR. A homogeneous field has no mean gradient, so s = 0 there, as the model says.
       - The quiet ambient is not exactly quiet. With the Spalart-Rumsey values, s = 9e3 nu^2 |grad X|^2 / (Sc_t beta* a^2 X (1 - X)). For a 100 um layer that is below 1e-3 at 3 MPa but can reach about 0.5 in hot gas at 1 atm (inferred, order of magnitude). The equilibrium variance also assumes the gas has spent about 1/(beta* omega) in the layer, which in a quiet ambient is milliseconds, so there it overstates s. The first build keeps the algebraic form, reports the s field per run, and names a transported variance equation as the upgrade if s outside the turbulent region is not negligible.
     - **Where it applies.** PaSR requires FiniteRate chemistry and turbulence; validation rejects it otherwise. LocalEquilibrium and Frozen are limits and take no closure.
     - **Verification, criteria in the test header before the code.**
       1. kappa_eff from the engine against the formula evaluated in the test from Cantera's rates, and s on a field whose mole fractions are linear in z (least-squares gradients exact there) against the analytic value: both within 1e-12.
       2. One chemistry substep of an auto-igniting H2/O2 cell with frozen tau_mix and s against an independent explicit integration in the test (Cantera rates, step refined until converged), within the integrator tolerance.
       3. Limits: k = 0 gives the laminar substep exactly; C_mix -> 0 approaches it; an inert mixture is unchanged exactly.
       4. A turbulent reacting chamber with PaSR (stated in full in `tests/pasr_tests.cpp`, 5 October): the C1 chamber made viscous and turbulent, with and without the closure, on 32x6 and 64x12 to 8 ms. Judged: budgets, positivity, settling, and admissible kappa_eff and s fields. Reported: the fields and the differences from the control.
8. **Walls.** Inviscid slip first (declared: no boundary-layer loss). Then no-slip under the turbulence closure with a declared wall temperature or heat-flux condition; a regeneratively cooled wall's temperature is an input or a bracket.
9. **Measurements.** Vacuum thrust from the device-thrust ledger with zero ambient; Isp from thrust over the total supplied mass flow (every propellant stream, liquid included, and any igniter flow); c* = Pc A_t / mdot; chamber pressure at the injector face and at the measured tap location; each loss as the delta between two runs of the same engine.
10. **Verification toward RL10, each against CEA:** premixed reactants with local equilibrium give CEA shifting c* and Isp after accounting for finite chamber area and divergence; finite-rate kinetics lands between shifting and frozen; the "valves open, ignite" start reaches the same end state as a partly developed start; grid and time-step refinement of the settled observables.
11. **RL10A-3-3A.** Geometry from NASA TM-107318 (area against axial station); operating inputs (propellant flows or O/F, feed temperatures) from the same or primary sources. LOX enters as liquid through step 5. The predicted vacuum Isp and thrust, with brackets for declared missing inputs, are committed to git before any measured performance is read. Compare, then explain the gap with the loss ledger (divergence, kinetics, boundary layer, mixing, atomization and evaporation). Missing inputs get declared brackets, never tuned values.
12. **Engines as data.** An engine definition holds its contour (a table, or a parametric shape from published dimensions), propellants and feed states, flows or O/F and chamber pressure, ambient pressure, and a provenance label on each field (measured, derived or assumed). A headless batch runner takes a database of such definitions and reports Isp, thrust and c* with their trends against CEA and measurements. There are no per-engine efficiency factors. A fast ideal-rocket tier (the engine's own thermochemistry, as in the CEA check) and the full time-evolving tier read the same definition.
13. **Output.** Ben wants the data representation before the UI: each run writes its time history (CSV) and field frames, with contact sheets and a single exported video per run.
14. **Laboratory sandbox app, first draft** (Ben; after the chemical milestone). "Draw and run": edit the cross-section (wall contour points), place injectors (supply faces and liquid injection sites), run, and watch the fields evolve. "Probe and plot": click anywhere for time histories of the local fields, with thrust and Isp traces for the whole device. It drives the same engine and run records as the headless runner.
15. **Cost.** Reaction integration dominates the step (96% of it on one thread, measured). Measure the cost per cell-step for each chemistry mode, use the worker threads, and run grid convergence on the Windows GPU PC or long Mac runs. The measured costs and the projection for a 1 s RL10 run are in *Compute budget for a 1 s RL10 run* below. The target since 5 October 2026 is a cold start to full thrust in about 5 minutes on the Mac; its budget is in *Lightweight engine*.
16. **Dual time stepping** (parked at 21:45 on 5 October 2026 by *Lightweight engine*; it is built only if that section's decision rule calls for it. It was selected earlier that day by the ruling in *Compute budget*. The design is kept here; no code exists). It removes the acoustic limit of the wall cells (about 7e-11 s on the RL10 grid) from the physical step and keeps the march physical (step 1). The inner pseudo-time iterations converge every physical step and are discarded, so they do not alter the history.
   - **Physical time: a one-step implicit method.** The candidate is the implicit part of Kennedy and Carpenter's ARK3(2)4L[2]SA: an ESDIRK with 3 implicit stages, third order, L-stable and stiffly accurate, with a second-order embedded estimate.
     - Because the method is one-step, the Strang split stays as it is: react(dt/2), the implicit flow step, react(dt/2). The CVODES chemistry is unchanged.
     - BDF2 is the alternative: one implicit solve per step instead of three. But its second history level is not defined for the flow sub-problem inside a split step.
     - Which of the two is cheaper at our tolerances is measured on the chamber before the choice is fixed.
   - **Step size.** The physical step comes from the embedded error estimate, with a declared relative tolerance on p, T and u. It is capped so that the convective CFL of the core flow stays below a declared value; the acoustic CFL of the wall cells may be large. Each run records its tolerance. Refining the time step is part of every settled result's verification (step 10).
   - **Inner iterations.** Each stage solves (U − U_known)/(a_ii dt) + R(U) = 0. R is the existing spatial operator, unchanged: HLLC with MUSCL, transport, SST, sources and boundaries.
     - Pseudo-time steps are local to each cell. That is allowed because only the converged stage value is used.
     - The implicit operator is approximate: the first-order (Rusanov) Jacobian of the inviscid flux, the thin-layer viscous Jacobian normal to the wall, and the SST destruction terms. The residual is exact, so the iterations converge to the second-order discretization (defect correction).
     - The stiff direction is wall-normal, along each column of rings at fixed z. It is solved as a block-tridiagonal system per column, with symmetric Gauss-Seidel sweeps in z.
     - Blocks are 4 + ns + 2 unknowns, 16 for h2o2.yaml.
     - A candidate to speed up convergence at chamber Mach 0.1 to 0.2 is low-Mach preconditioning applied in pseudo time only. It does not change the converged solution.
   - **Conservation.** The step's update is formed from the stage residuals, U(n+1) = U(n) + dt sum b_i R(U_i), and the boundary fluxes are accumulated with the same weights. So the budgets close to round-off whatever the inner tolerance. The difference from the last stage (equal in exact arithmetic for a stiffly accurate method) is logged as the inner-convergence error. The inner tolerance is declared per run and recorded.
   - **Admissibility.** The engine checks that rho, p, T, Y_k, k and omega are positive at every inner iteration (using a damped update) and at the end of the step. A failed step is rejected and retried at dt/2, as the explicit steps already do. No clipping.
   - **Threads.** Within each sweep colour (red-black in z), the residual and the column solves are independent across columns, so the flow step becomes threaded. Today the serial flow step is the second-largest cost (measured).
   - **Verification**, stated as criteria before the first run:
     1. Temporal order on a smooth unsteady case (an acoustic wave in the duct, or pipe decay) against an explicit run at a tiny step.
     2. The settled pipe of TURBULENCE_C2 check 4 matches the explicit result to the inner tolerance.
     3. The C1 "valves open, ignite" history (light-off time, injector pressure, outlet mass flow) converges to the explicit run as the physical step is refined, at three steps.
     4. Budgets close to round-off at any inner tolerance.
     5. The measured cost per simulated microsecond on a wall-clustered chamber with RL10-like wall-cell aspect ratios, against explicit.
   - **Open, settled by measurement:**
     - ESDIRK against BDF2, and the third-order ESDIRK against Kennedy and Carpenter's fourth-order one (ARK4(3)6L[2]SA, five implicit stages): Bijl et al. found the fourth order the most efficient below 10% error (read below).
     - The step update: the residual-weighted form above (budgets to round-off) against the last stage value (Bijl et al., p. 4, who call the weighted form "potentially damaging"). Inferred reason: R at a stage converged only to the inner tolerance carries that error times the stiffness, here up to the wall cells' acoustic CFL. Measured on the wall-clustered chamber. If the last stage is chosen, verification criterion 4 is restated before any run.
     - The inner tolerance at which the transient stops depending on it. Starting points from the papers: 1/10 of the temporal tolerance (Bijl et al., measured) to 0.005 of it (Kennedy and Carpenter).
     - Point LU-SGS (scalar diagonal, Yoon and Jameson) against wall-normal block-tridiagonal lines, by convergence per unit cost on a wall-clustered chamber.
     - Whether the chemistry must move inside the stage residual (point-implicit) at physical steps of 1e-7 to 1e-6 s. The Strang splitting error at those steps is measured on the C1 light-off (criterion 3) first.
   - **Sources.** The titles, venues and volumes were checked on 5 October 2026 against index pages. Four have been read (5 October, below); Jameson 1991 has not. The design above is checked against them; the tableau is entered in code from the TM's appendix, not retyped from here.
   - **Read on 5 October 2026** (a Sonnet agent read the full text; page numbers are the papers').
     - Jameson, "Application of dual time stepping to fully implicit Runge Kutta schemes for unsteady flow calculations" ([PDF](http://aero-comlab.stanford.edu/Papers/jameson_dts_irk.pdf); no year or venue printed on this copy).
       - Dual time with BDF2 or with fully coupled implicit RK (Gauss, Radau IIA), smoothed by one LU-SGS sweep per direction on a first-order Roe Jacobian (p. 8).
       - A plain pseudo-time march on fully coupled RK stages can be unstable in pseudo time at small physical steps (p. 5). The fix is preconditioning by the inverse RK matrix (pp. 6-7). The ESDIRK here solves one stage at a time with a scalar a_ii, so that mechanism does not carry over (inferred). The 1/(a_ii dt) term still goes on the implicit operator's diagonal, so the physical-time term is never explicit in pseudo time (inferred; the standard remedy).
       - The inner tolerance has sharp thresholds: one 2-stage Gauss airfoil case fails at 14 inner iterations and runs at 15 (p. 8). So the inner tolerance is measured, not assumed (already an open item above).
       - DIRK schemes are noted as needing many stages: Kennedy and Carpenter's fourth-order scheme has one explicit and five implicit stages (p. 4). The third-order ARK3(2)4L[2]SA has three implicit stages. Its cost against BDF2 stays a measured choice.
     - Yoon and Jameson 1988 (AIAA J. 26(9), 1025-1026; a synoptic of the full paper, presented as AIAA Paper 87-0600).
       - The Jacobian is split by the spectral radius, A± = (A ± r_A I)/2. This gives a scalar diagonal, so no block inversions are needed (p. 1025). This supports the Rusanov-type Jacobian above.
       - It argues against line solves: avoiding the block inversions of line Gauss-Seidel is its point (p. 1025). The design's wall-normal block-tridiagonal lines go the other way, for the RL10 wall cells' aspect ratios (inferred to need them; not shown). So point LU-SGS against wall-normal lines is added to the open items, settled by measured convergence on a wall-clustered chamber.
     - Kennedy and Carpenter, NASA/TM-2001-211038 (July 2001; the report version of the 2003 paper; [PDF](https://www.cs.odu.edu/~mln/ltrs-pdfs/NASA-2001-tm211038.pdf)).
       - ARK3(2)4L[2]SA's implicit part (Appendix D, p. 47; properties Appendix B, p. 45): gamma = a_ii = 1767732205903/4055673282236, about 0.43587; c = (0, 2 gamma, 3/5, 1); stiffly accurate (b is the last row), so A-stable and L-stable; third order; stage order 2. The embedded second-order method is not L-stable (Appendix C, p. 46). This matches the design.
       - Step size: a PID controller is recommended (pp. 16-17, eqs. 39-43), with k_I 0.25, k_P 0.14, k_D 0.10 and a safety factor of about 0.9. The design's "embedded error estimate with a declared tolerance" takes this controller.
       - Stage predictor: dense-output extrapolation, with the previous stage as the fallback when the step ratio is large (p. 16).
       - Inner tolerance: residual and displacement tolerances about 0.005 times the temporal tolerance (p. 15, eq. 36).
       - Order reduction (section 8.3, pp. 30-34; Table 15, p. 33): with stage order 2, stiff (algebraic-like) components converge as dt^3 + eps dt^2. Criterion 1's temporal order is therefore measured on a smooth case where that does not apply, and the chamber's observed order is reported, not assumed.
       - Operator splitting is not discussed: the paper presents additive RK as an alternative to splitting, with the stiff terms in the implicit part. That is the route if the chemistry has to move inside the stage residual (open item above).
     - Bijl, Carpenter and Vatsa, AIAA Paper 2001-2612, "Time integration schemes for the unsteady Navier-Stokes equations" (the conference precursor of the 2002 J. Comput. Phys. paper, which is not open; [PDF](https://ntrs.nasa.gov/api/citations/20010066914/downloads/20010066914.pdf)).
       - Laminar cylinder at Re 1200, Mach 0.3: BDF2 is 2.5 times less efficient than ESDIRK4 at 10% error in lift; at 1% ESDIRK4 needs 1.5% of BDF1's work; at 0.1% it needs 10% of BDF2's (p. 9). Fourth order is recommended over fifth for robustness and storage (pp. 9, 11).
       - The inner (multigrid) iterations must be converged to at least 1/10 of the desired accuracy; 1/2 degrades it, and 1/20 and 1/200 are equivalent (p. 10).
       - Stiff accuracy is used to take U(n+1) as the last stage, removing "the potentially damaging explicit update" U(n) + dt sum b_j R_j (p. 4). Added to the open items above.
   - **Citations** (titles, venues and volumes checked against index pages on 5 October 2026).
     - Jameson 1991, "Time dependent calculations using multigrid, with applications to unsteady flows past airfoils and wings", AIAA Paper 91-1596: dual time stepping ([cited in Jameson's later paper](http://aero-comlab.stanford.edu/Papers/jameson_dts_irk.pdf)).
     - Kennedy and Carpenter 2003, "Additive Runge-Kutta schemes for convection-diffusion-reaction equations", Appl. Numer. Math. 44(1), 139-181: ARK3(2)4L[2]SA and its ESDIRK part ([Semantic Scholar](https://www.semanticscholar.org/paper/Additive-Runge-Kutta-Schemes-for-Equations-Kennedy-Carpenter/ad463b85089ac66ae41dad57e06523403acb11e6)).
     - Bijl, Carpenter, Vatsa and Kennedy 2002, "Implicit time integration schemes for the unsteady compressible Navier-Stokes equations: laminar flow", J. Comput. Phys. 179, 1-17: ESDIRK against BDF2 for unsteady flow.
     - Yoon and Jameson 1988, "Lower-upper symmetric-Gauss-Seidel method for the Euler and Navier-Stokes equations", AIAA J. 26, 1025-1026: LU-SGS ([PDF](http://aero-comlab.stanford.edu/Papers/AIAA-10007-471.pdf)).
17. **Wall functions: the default wall treatment** (Ben's direction of 5 October 2026, 21:45, which reverses that morning's ruling "a declared option, never the default"). An SST wall treatment with no resolved sublayer, designed in *Lightweight engine*. The wall-resolved treatment of step 7 is kept as its verification reference. The method is Nichols and Nelson's; its verification criteria were stated before code on 5 October (`docs/evidence/WALL_FUNCTIONS.md`, summarised in *Lightweight engine*), and every result labels its wall treatment.

### Compute budget for a 1 s RL10 run (5 October 2026)

Ben asked what a full 1 s RL10 run costs. Every number is labelled: **measured** (on the named machine), **derived** (stated arithmetic from measured or published values) or **guessed** (an input not yet known).

**Measured cost of the reacting step.** Measured with `tests/step_cost.cpp`; the output of four runs (and one on backhouse) is in `docs/evidence/step_cost_2026-10-05.txt`. Another project's job shared the Mac during run 4; its rates match runs 1 to 3.
- Machine: Mac M2 Pro (6 performance and 4 efficiency cores).
- Case: the C1 chamber on 64x12 (768 cells) with h2o2.yaml (10 species). It is marched in local equilibrium to 0.6 ms, when the chamber burns at 1.9 MPa and the nozzle flows.
- Timing: 200 steps per configuration at that grid's flow step, about 5.7e-8 s.
- Settings: FiniteRate is CVODES at rtol 1e-6 and atol 1e-12, as in C1. "Viscous" means mixture-averaged transport, SST-2003 and no-slip walls at 600 K.

| Step (cell updates per second, range over four runs) | 1 thread | 10 threads |
|---|---|---|
| Inviscid, frozen chemistry (the flow step alone) | 7.2e5 to 7.4e5 | not run (serial) |
| Inviscid, finite rate | 8.7e3 to 9.6e3 | 4.4e4 to 5.2e4 |
| Viscous (SST), frozen chemistry | 2.3e5 to 2.4e5 | not run (serial) |
| Viscous (SST), finite rate | 8.3e3 to 9.3e3 | 3.7e4 to 4.3e4 |
| Inviscid, local equilibrium (one run) | 1.9e4 | 6.3e4 |
| Viscous (SST), local equilibrium (one run) | 1.8e4 | 5.4e4 |

- Chemistry is 96% of the single-thread viscous step (derived). It costs about 1.0e-4 s per cell per step.
- **Local (shifting) equilibrium needs no CVODES, but it is not much cheaper** (measured with the `equilibrium` option of `step_cost`, 5 October 2026).
  - Viscous, per step: 4.2e-2 s on one thread, against 8.3e-2 s for finite rate. That is 2.0 times cheaper.
  - On 10 threads: 1.4e-2 s against 1.8e-2 s, only 1.25 times cheaper. Equilibrium gains 3.0x from 10 threads where finite rate gains 4.6x; the cause is not measured.
  - One cell's constant-(u, v) equilibrium costs about 5.1e-5 s (derived from the step). That agrees with the 62 us of an isolated near-equilibrium cell (measured, `crucible_chemistry_cost`).
  - It is also a different physical model (no finite-rate kinetics, no ignition delay), so it cannot stand in for finite rate on the reference path.
- **A chemistry cache would buy about 2x on this chamber, and at most 5x** (the `cache` option of `step_cost`, 5 October 2026; output in `docs/evidence/chemistry_cache_2026-10-05.txt`).
  - The stream (measured): every cell of every reaction call in 200 viscous finite-rate steps on 64x12, 464,640 calls. The steps start from the inviscid equilibrium setup state, so the cells are adjusting, not settled.
  - The table: calls are binned by T, Y_k, ln rho and ln dt. A call whose bin already holds an earlier call is a hit. It is answered from that record either as the record's increment (constant) or plus a first-order correction from the record's Jacobian (linear, as ISAT, [Pope 1997](https://iopscience.iop.org/article/10.1088/1364-7830/1/1/006)).
  - Accuracy is measured against each call's own CVODES result, on CVODES's weighted scale (1 = its local tolerance). On that scale the march's own integration error, against rtol 1e-10, is 0.74 at the median and 4.4 at the 99th percentile (measured).
  - Measured, bins of 1 K, 1e-3 in Y, 0.1% in rho and 1% in dt:

    | | Hits | Hits within weighted error 10, linear | Same, constant |
    |---|---|---|---|
    | Records from any cell | 85% of calls | 84% | 24% |
    | Records from other cells only | 28% | 28% | 6% |

    - Most reuse is a cell reusing its own previous call. Across cells only 28% of calls could be answered within 10.
    - Wider bins hit more often but answer worse: at 10 K bins only 31% of calls are within 10. Narrower bins hit less: at 0.1 K, 56%.
    - The worst linear answer at 1 K bins has a weighted error of 2400 (|dY| 1.5e-6, against 8e-8 for the integration). So the cache needs error control on every retrieval, as ISAT has.
  - Cost:
    - One retrieval costs 3.1e-7 s (measured), about 1% of a direct call (3.4e-5 s, derived above).
    - Each record needs a Jacobian of the reaction map. Here it took 14 CVODES calls (forward differences). With CVODES sensitivities it is guessed at 2 to 5 calls.
  - What it buys (derived, with 96% of the step in chemistry):
    - With free records: a 6x saving on the chemistry and 5x on the step.
    - With records at 2 direct calls: 2x on the step. At 5 calls: nothing.
    - From reuse across cells alone: at most 1.4x.
  - Pope reports a thousand-fold saving, on a statistically stationary reactor ([abstract](https://iopscience.iop.org/article/10.1088/1364-7830/1/1/006)). Here each record serves only about 5 later calls. A light-off, or the larger state changes per step of dual time stepping, should give less reuse (inferred).
  - Inferred: tabulation is not the large lever. The larger one is the cost per call. A CVODES call costs as much as about 15 to 40 rate evaluations, whatever the substep (derived from the measured costs above). A lean stiff integrator for short substeps, with an analytic Jacobian, is the candidate, and it would also suit the GPU port. Its gain is not measured.
- 10 threads give about 4.6x. The flow step is serial, and 4 of the 10 cores are efficiency cores.
- Each step makes 3.0 reaction substeps rather than 2 (measured). Once per step the first half-substep lowers the CFL step and is redone ("replans", 1.01 per step).
- One CVODES substep costs about 3.4e-5 s per cell at the measured half step (derived from the above).
- **Shortening the substep barely lowers its cost.**
  - The state: after 20 finite-rate steps with the valve open. 20 substeps were repeated on that one state.
  - Measured cost per cell of one substep:
    - 1.7e-5 s at 2.8e-8 s (half the flow step);
    - 1.35e-5 s at 5e-10 s;
    - 1.35e-5 s at 5e-11 s.
  - The repeated substeps relax the state toward equilibrium, so these are about half the in-march cost.
  - So CVODES costs about the same per call whatever the substep length: the cost is per call, not per simulated time. At the RL10 step, chemistry by CVODES costs about as much per step as it does here.
- One Cantera rate evaluation costs 0.78 to 0.84 us (measured, `crucible_chemistry_cost`, two runs). A single isolated cell costs 100 to 144 us per 50 ns CVODES call at rtol 1e-6 (measured).

**RL10 grid and time step** (derived from a guessed resolution; the real contour is extracted at C4).
- The domain is 1.37 m long (TM-107318 Table E1, stations -12 in to +41.84 in), with area ratio 61. The throat radius is about 6 cm (guessed until the contour is extracted).
- Axial cells: uniform 1 mm, so nz is about 1400 (guessed resolution, about 60 cells along the throat radius).
- Radial rings: about 120, clustered to a first wall cell of 0.4 um.
  - That is y+ about 1 at the throat, which SST-2003 needs to integrate to the wall.
  - It is derived from a skin-friction coefficient of about 2e-3 and the gas viscosity at a 600 K wall, at a chamber pressure of a few MPa (guessed). It is uncertain by about 3x.
- About 1.7e5 cells.
- The explicit step is set by sound crossing the wall cell: dt = 0.2 dr / a at the engine's CFL of 0.4. That is about 7e-11 s for a = 1200 m/s (derived).
  - The diffusion bound, 0.1 dr^2 / D with D up to about 3e-5 m^2/s for H2 (guessed), is about 5e-10 s or more, so it does not bind.
- 1 s is therefore about 1.5e10 steps and 2.5e15 cell updates (derived).

**Projected wall time for 1 s of simulated time** (explicit, wall-resolved, as the engine is built today)

| Machine | Cell updates per second | Wall time for 1 s |
|---|---|---|
| Mac, 10 threads | 4e4 (measured on 768 cells) | 6e10 s, about 2000 years (derived) |
| backhouse CPU (i7-14700K: 8 performance and 12 efficiency cores, 28 threads) | 7.3e4 (measured, one run, 28 threads, with family use on the box; 1.2e4 on one thread) | 3.4e10 s, about 1100 years (derived) |
| GPU port (RTX 4070 Ti SUPER, f64) | 1e6 to 1e7 (guessed; a consumer card runs f64 at 1/64 of its FP32 rate) | 8 to 80 years (derived from the guess) |

A 1 s explicit run at wall-resolved resolution is out of reach on every machine here. Taking a week as acceptable, it is too slow by about 400 times on the GPU (guessed) and about 1e5 times on the Mac.

**What would cut it** (factors multiply only where the rows are independent)

| Cut | Factor | Cost to the physics |
|---|---|---|
| Simulate to the settled state instead of 1 s | About 50 (guessed: settled by about 20 ms; the settling test measures it) | None. Step 1 already takes performance from a late settled window. |
| Plan each step slightly below the CFL limit, so the first half-substep is not redone | 1.5 on the chemistry (derived from the measured 3.0 substeps per step) | None. |
| Integrate the chemistry inside the explicit step | About 10 to 20 on the step (guessed). At steps of 1e-10 s the kinetics are not stiff. Explicit integration then needs 2 to 4 rate evaluations of 0.84 us each, against about 1.0e-4 s per cell of CVODES chemistry per step (measured). The flow step, 4.2 us per cell viscous on one thread (measured), then dominates. | None while the step is shorter than the fastest chemical time. Each cell must check this and fall back to CVODES when it fails. |
| Thread the flow step | Up to about 8 on 10 cores, once chemistry is no longer dominant (guessed). Measured on the inviscid C1 with Table A: 3.0 on 4 threads and 4.2 on 6, bit-identical to serial (5 October, `docs/evidence/THREAD_POOL.md`). Transport is still serial | None. |
| Wall functions: first cell at y+ 30 to 100 instead of 1 | 30 to 100 on the step (derived from the cell size), and fewer rings | A declared wall model replaces the resolved sublayer. Wall heat flux and friction, which are entries in the RL10 loss ledger, become model outputs. |
| Implicit or dual time stepping for the flow | 10 to 700 (guessed) | None if the physical step resolves the transient. The inner pseudo-time iterations converge each physical step and do not alter the history, so this is consistent with step 1. Assumes a physical step of 1e-7 to 1e-6 s with 10 to 30 inner iterations, each costing 2 to 5 explicit steps. Needs a new solver and its verification, about a week (guessed). |
| Chemistry tabulation (ISAT, Pope 1997) | About 2 on the step, at most 5 (derived from the measured reuse, 5 October; see *A chemistry cache* above). It replaces the guess of 5 to 20. | An approximation with an error tolerance that must be verified, and it needs error control on every retrieval. |
| Local time stepping once settled | 1e4 or more on that part (guessed) | The history after the switch is pseudo-time. That conflicts with step 1 as written ("no pseudo-time acceleration"), so it needs Ben's ruling. |
| GPU port | 20 to 200 over the Mac (guessed) | None in f64. A large port of the flow and the chemistry. |

**What is feasible** (derived from the rows above; the factors are guessed).
- **1 s:** not feasible on any of these machines with any combination short of local time stepping.
- **A settled run of about 20 ms needs two structural changes.** One is a wall treatment or time stepping that removes the 7e-11 s wall step: wall functions, or implicit or dual time stepping. The other is cheaper chemistry or a GPU.
- **Example: wall functions plus in-step chemistry plus a threaded flow step.**
  - About 1e7 steps of 2e-9 s on about 1e5 cells, so about 1e12 cell updates.
  - That is about 10 days on the Mac at a guessed 1e6 cell updates per second, a few days on backhouse, and a few hours to a day on a GPU port.
- **Example: dual time stepping on the wall-resolved grid.**
  - About 2e5 physical steps of 1e-7 s on 1.7e5 cells, each costing about 3e-4 s per cell on one thread.
  - That is about 4 weeks on the Mac with 10 threads (at the measured 4.6x) and about 2 weeks on backhouse.
  - Both are 10 times less with a physical step of 1e-6 s.

**Ruling** (Ben's decision procedure, through the Director, 5 October 2026). Superseded at 21:45 the same day by Ben's direction in *Lightweight engine* below; kept as the record.
- **Main line: dual time stepping on the wall-resolved grid**, plus the cuts that cost no physics (plan below the CFL limit, chemistry inside the explicit step, a threaded flow step). The reason: the chemical baseline is a time-evolving chamber at full fidelity. For an expander-cycle engine such as the RL10, wall heat flux drives the cycle, so on the reference path it must be resolved, not a model output. The design is step 16.
- **Parallel experiment: wall functions**, a declared, switchable option and never the default. Before any use, it is checked against the wall-resolved dual-time result on one short case. Database sweeps (hundreds of engines) may use it only after that comparison is recorded, and their results are labelled with it.
- **Local time stepping is rejected**, per Ben's earlier ruling against pseudo-time acceleration (step 1).
- **Chemistry first.** Chemistry is 96% of the step, and cheaper chemistry helps both lines and the GPU question. Before solver code is written, two things are measured and labelled: the step cost in the local-equilibrium mode (no CVODES), and what a chemistry cache or tabulation (ISAT class) would buy on the chamber.
- **C4 slips by up to one week** for the dual-time solver and its verification. The dated plan below carries the re-estimate.

### Lightweight engine (5 October 2026, from 21:45)

**Direction** (Ben, through the Director, 5 October 2026, 21:45, verbatim): "I would like this engine simulator to not be crazy compute heavy. It shouldn't be resolving needless turbulence, adaptive resolution is key. After all, it is 2d. Ideally I should be able to do a run from cold to full thrust in minutes, heavy chemistry for lookup tables done OFFLINE." It supersedes the *Ruling* above (dual time stepping on the wall-resolved grid as the main line).

The new main line:
- Runtime chemistry comes from tables built offline with Cantera (on backhouse, or the Mac when small). The tables are validated against Cantera and CEA, with labelled errors. CVODES stays only for building tables and for reference checks.
- Wall functions are the wall treatment, with no resolved sublayer (step 17, restated).
- 2-D adaptive resolution where the physics needs it: flame front, shocks, throat, shear layers.
- Implicit or dual time stepping only if the measured budget still needs it (step 16 is parked; decision rule below).
- The validation bar is unchanged: c*, Isp and thrust against CEA and RL10 data, plus trends. The march stays physical from a cold start (step 1).

**Target.** A cold start to full thrust in about 5 minutes on the Mac (M2 Pro).

**Budget.** Labels: measured (on the named machine), derived (stated arithmetic), guessed (an input not yet known).

| Quantity | C1 chamber, 64x12 | RL10 (contour guessed until C4) |
|---|---|---|
| Cells | 768 on the fixed mesh. No adaptivity needed (measured: c* within 0.40% of the corrected ideal rocket on this grid) | Base level 140 x 16 = 2,240 (10 mm axial, rings to the wall). Level 1 (ratio 2) on 30% of the base: 2,700. Level 2 (ratio 2) on 10%: 3,600. Total about 8,500 (guessed fractions, derived counts). The same finest resolution everywhere would be 560 x 64 = 35,800 (derived) |
| First wall cell | Inviscid (slip) in C1. The turbulent C1 runs have first-cell y+ 56 to 105 (laminar estimate, measured) | 0.94 mm at level 2 at the throat, centre at y+ about 1,000 to 2,000. The viscous length is 0.24 to 0.4 um (derived from a guessed c_f 2e-3, a 600 K wall and a chamber pressure of a few MPa). That is inside the log layer if the throat boundary layer is 2 to 4 mm thick (guessed) |
| Step | 5.9e-8 s (measured) | Level 0 6e-7 s, level 1 3e-7 s, level 2 1.5e-7 s, each level subcycled (derived from the cell sizes, a throat sound speed of about 1,500 m/s and the engine's CFL of 0.4; the inputs are guessed) |
| Cell updates per simulated ms | 1.3e7 (derived) | 3.7e7 (derived: the sum over levels of cells over step) |
| Physical time, cold to full thrust | About 4 ms. Vacuum thrust is within 0.1% of its settled value at 3.5 to 4 ms (measured, equilibrium 64x12, `docs/evidence/c1/eq64_log.txt`). The judged C1 protocol runs to 8 ms | 10 to 20 ms after the supplies open (guessed). Derived from a fill time constant L* / (Gamma^2 c*) of about 1 ms for a guessed L* of 0.8 to 1 m, seven time constants to settle to 0.1%, and a guessed 2 to 3 ms for the gas to cross the chamber. The cycle's own start (turbopumps, the expander) is not modelled: the supply faces impose the propellant flow schedule (step 4) |
| Cost per cell update, one thread | Inviscid flow step 1.35 to 1.39 us (measured, `step_cost`), plus a table lookup of 0.07 to 0.08 us for a coherent query stream or 0.33 to 0.35 us for a random one (measured, prototype `tests/table_lookup_cost.cpp`, `docs/evidence/table_lookup_cost_2026-10-05.txt`). About 1.5 us (derived) | Viscous SST flow step 4.2 to 4.3 us (measured) plus the lookup, wall functions (wall cells only, guessed negligible) and AMR overhead (ghost cells, refluxing, regridding; guessed 25%). About 5.5 us (derived) |
| Today's chemistry cost, for comparison | Cantera equilibrium about 51 us per cell, CVODES finite rate about 100 us per cell (measured) | the same |
| Wall time, one thread | 78 s to 4 ms, 156 s to 8 ms (derived). With Table A: 84 s to 4 ms and 168 s to 8 ms (measured, 5 October). Cantera equilibrium: 2,375 s to 8 ms on 5 threads (measured, 4 October; 2,048 s on 5 October) | 2,000 to 4,100 s (derived) |
| Wall time with a threaded flow step (4 to 5 times on 10 cores, guessed; the chemistry measured 4.6 times) | 16 to 40 s (derived). Measured 5 October with Table A: 28 s to 4 ms on 4 threads and 20 s on 6, 3.0 and 4.2 times the 1-thread run (`docs/evidence/THREAD_POOL.md`) | 400 to 800 s, about 7 to 14 min (derived). Transport, most of the viscous step, is not yet threaded |

What the budget says (derived from the rows above):
- **C1 meets the target with tables alone**, on one thread. That is tonight's build and its measurement (*Table A* below). The turbulent C1 (viscous SST, 4.3 us per cell update) would take about 230 s to 4 ms on one thread once its finite-rate chemistry is tabulated (*Table B*).
- **RL10 needs every lever and is still about 1.5 to 3 times over 5 minutes** on these guesses. The candidates for the rest, in order of expected yield per unit work:
  1. Tabulated transport properties. The viscous step costs 3.1 times the inviscid one (measured); the share spent evaluating mixture-averaged properties is not measured. A profile comes first.
  2. The GPU port (guessed 20 to 200 times over the Mac, *Compute budget*).
  3. A wall-normal implicit step for the wall strip only, if the measured wall-strip step binds.
- **Wall functions pay only if they hold far out.** At the y+ 30 to 100 of the old step 17, the throat wall cell is 7 to 40 um, and the explicit step about 1e-9 to 5e-9 s (derived). The design therefore uses the wall function up to the top of the log layer, which is y+ about 1,000 to 2,000 at the throat (guessed boundary-layer thickness). That is beyond the range tested in the papers read: grid independence of skin friction up to y+ 60 to 80 with SST's automatic treatment (Alrutz and Knopp, DLR). Verification decides (*Wall functions* below). If the wall function holds only to y+ about 300, the throat strip needs one or two more levels. Two rows of cells along the throat region at levels 3 and 4 would add about 25% to 100% (guessed).
- **Implicit time stepping: the decision rule.** Measure the RL10-like chamber's cost per simulated ms after Table A, threading, wall functions and AMR exist. Design an implicit step only if the measured cold start exceeds 15 minutes on the Mac (three times the target) and the binding step is the wall strip's acoustic limit. Step 16's dual-time design is kept for that case.

**Table method** (literature: a Sonnet agent's search on 5 October 2026; full text read where stated).
- **Candidates.**
  - Flamelet/progress-variable (FPV; Pierce and Moin, J. Fluid Mech. 504, 2004, abstract read) and FGM (van Oijen and de Goey, Combust. Sci. Technol. 161, 2000, abstract read): two scalars, mixture fraction and a progress variable, map onto precomputed flamelets.
  - Compressible FPV (Saghafian, Terrapon and Pitsch, Combust. Flame 162, 2015, equations read as reproduced in Baumgart, Yao and Blanquart 2025, LLNL-CONF-2002007, full text read): temperature from energy through a linearised gamma(T). Baumgart et al. find that linearisation errs by several hundred K, and the progress-variable source by up to about 10%, far from the flamelet's reference state (their test: detonation states).
  - FPV in rocket chambers: a ten-injector H2/O2 FPV study (arXiv 2211.06594, full text read) tabulates on (Z, Z variance, C), applies pressure by a power law rather than a table axis, starts from a hot chamber because "no explicit ignition model is included", and reports no comparison of c* or Isp with CEA or finite rate. The agent found no rocket use of RCCE.
  - ISAT-class caching was measured here to buy about 2 times on the chamber (*Compute budget*).
- **Choice, in two tables built in order.** Both store all species mass fractions (the composition is part of the medium's state, VISION_SCOPE, *The physical scope*). Temperature always comes from the core's exact e-to-T inversion of the NASA-7 mixture, so the energy is conserved exactly and the gamma-linearisation error above does not arise. Multilinear interpolation on axes in which the element mass fractions are linear (or bilinear) reproduces the query's elements exactly, so element conservation is exact by construction.
  - **Table A: shifting equilibrium.** The composition is a function of the element fractions, e and rho. It is exactly today's LocalEquilibrium model, so it can be verified to interpolation error against Cantera with no change of model, and it is the CEA shifting-equilibrium limit that C1's c* and Isp are judged against. It cannot represent ignition delay or kinetic lag (it burns cold premixed gas at once, as LocalEquilibrium does).
    - Axes as built (5 October, about 22:20; the first design, with axes f_H, Z_N, normalised e and ln rho, missed criterion 2(a) by about 100 times, see `docs/evidence/TABLE_A.md`): f_H = Z_H / (Z_H + Z_O), the temperature T and ln rho_r, with rho_r = rho (Z_H + Z_O) the partial density of the reacting H/O species.
    - N2 and Ar drop out exactly. Each is the only carrier of its element in `h2o2.yaml`, so its amount is fixed. In an ideal gas the H/O species' equilibrium at a given T depends only on their own elements and partial density. The builder checks the single-carrier condition, so a mechanism with NOx chemistry would be refused.
    - The lookup finds the T at which the interpolated energy, the reacting mass's plus the diluents', equals the cell's e (it increases with T), then interpolates Y there.
    - The f axis is clustered geometrically at the stoichiometric f and at both ends, where the composition has a near-kink.
    - Built and verified on 5 October (`docs/evidence/TABLE_A.md`). The table has 131 x 129 x 97 nodes (144 MB, 11 s to build on the Mac). Its 99th-percentile |dT| against direct Cantera is 1.06 K on random states, and its largest is 0.83 K on C1's own states.
  - **Table B: a finite-rate progress manifold** (next; design and research pass before code). The candidate: trajectories of 0-D constant-volume reactors integrated offline by CVODES on the full mechanism, parameterized by (f_H, Z_N, e, ln rho, c), with c a progress variable. The table stores Y_k(c) and the rate dc/dt. The reaction step advances each cell's c and sets its composition from the manifold.
    - Why a 0-D reactor manifold rather than flamelets: PaSR's fine structure is a homogeneous reactor at the cell's state, which is what these trajectories are. A flamelet table brings its own turbulence-chemistry closure (a presumed PDF over Z) that would duplicate PaSR.
    - Open, settled before code:
      - a progress variable that resolves the induction period (radical growth with little H2O), so that the ignition delay survives;
      - the table size (five axes: about 2 GB at Table A's resolution with 24 nodes in c, so the axes must be coarser or adaptive);
      - mixing states that lie off the manifold.
    - Verification against CVODES: 0-D ignition delays across the C1 and RL10 range, the FLAME_C2 free flame, and the C1 light-off time of the finite-rate runs.

**Wall functions** (step 17, restated).
- Read in full: Menter, Carregal Ferreira, Esch and Konno, IGTC2003-TS-059 (SST's automatic wall treatment). The velocity and omega blend between the sublayer and log-layer solutions (eqs. 15 to 18). The wall heat flux follows Kader's thermal law, Theta+ = (T_w - T) rho c_p u_tau / q_w as a function of Pr and y+ (eqs. 11 to 14). It switches "gradually ... from a classical low-Re formulation on fine grids to a log-wall function formulation on coarser meshes" (p. 2).
- The chamber's wall layer is strongly compressible and heated: a wall at about 600 K under gas near 3,000 K is a density ratio of about 5. Kader's law is incompressible. A related NASA compressible derivation (De Chant and Tattar, CR-191185, 1994, full text read) reports 12.5% error for adiabatic flow and 18.5% with heat transfer against experiments.
- **Choice: Nichols and Nelson** (AIAA J. 42(6), 2004). Read in full from the author's own chapter of the method (Nichols, *Turbulence Models and Their Application to Complex Flows*, Rev. 4.01, ch. 10). The journal article is paywalled and unread. Spalding's single formula from the wall to the log layer, with White and Christoph's compressible, heated outer law. Crocco-Busemann gives the temperature, and the first cell's k and omega are prescribed. One formula holds at every y+, so there is no switch to the resolved wall. Two printed equations (10.9, 10.13) disagree with their derivation; the corrected forms are used (`docs/evidence/WALL_FUNCTIONS.md`). The author caps the first point at y+ 100 (a rule of thumb of 50), against this plan's design point of 1,000 to 2,000 at the throat.
- Verification: criteria stated before code on 5 October in `docs/evidence/WALL_FUNCTIONS.md`.
  0. The algebra (limits, inversion, the derivative), bit identity on 1 and 4 threads, and closed budgets.
  1. A fully developed pipe at Re_tau about 10,000 (about the RL10 throat's, derived), not the Re_tau 182 pipe of TURBULENCE_C2. Wall functions at first-cell y+ 1, 30, 100, 300 and 1,000 against the wall-resolved SST reference: c_f within 5%.
  2. A heated compressible pipe at the RL10's density ratio of about 5: c_f within 5% and the Stanton number within 10%.
  3. C1 turbulent (Table A, SST) on 12 rings against the wall-resolved run: c*, Isp and thrust within 0.2%, wall heat flow within 10%, wall axial force within 5%.
  4. At most 10% extra cost per cell update.
  - The largest y+ passing 1 and 2 is y+_max, which sets the wall strip's AMR level. Below 300, the budget in *Lightweight engine* is re-derived before AMR is built.
- Every result labels its wall treatment. The wall heat flux and friction become model outputs, with error bands in the RL10 loss ledger from criteria 2 and 3.

**Adaptive resolution** (Berger and Colella, J. Comput. Phys. 82, 1989, full text read; Berger and Oliger 1984 by summary).
- **Choice: block-structured AMR in the logical (i, j) index space of the body-fitted mapping**, with refinement ratio 2 in space and time.
  - Refined patches are rectangles of index space, so each patch is a small structured mesh and the existing kernels (HLLC, MUSCL, transport, SST, sources) run on it unchanged.
  - Time-accurate subcycling: each level takes two steps per step of its parent, so the same explicit scheme is stable on every level (Berger and Colella, p. 67). This is not local time stepping: every level marches in physical time.
  - Conservation: coarse cells under a fine patch are replaced by the conservative average of the fine cells. At coarse-fine faces, the coarse flux is replaced by the sum of the fine fluxes over the fine substeps (refluxing, their eqs. 1 to 3, pp. 69-70).
  - Geometry from the finest level. The mapping (contour and ring distribution) is evaluated once at the finest allowed level. A coarser cell is the exact union of its fine cells, and its faces are the polylines of the fine faces, with area vectors and volumes summed. Coarse and fine geometry then nest exactly; refluxing and free-stream preservation stay exact (the meridional area equals the closed sum of r n dl for any polygon). The wall and the throat area are as accurate on every level as on the finest. If the coarse cells were straight-edged instead, a 10 mm chord on the throat arc (radius about 6 cm) would err by 0.2 mm, and the throat area by about 0.7% (derived).
  - Tags, re-evaluated every few base steps:
    - static: the throat region, and a wall strip held at a target y+ from the wall function's u_tau;
    - dynamic: the flame front (heat release or the gradient of the progress variable), shocks (pressure jump above a threshold; Berger and Colella treat shocks apart from their error estimate because shock capturing is zeroth order there, pp. 75-76), shear layers (vorticity) and mixing layers (the gradient of f_H).
- **Not chosen.**
  - AMReX: its complex-geometry path is embedded boundaries on Cartesian grids (from its documentation, not re-read tonight). Adopting it would replace the body-fitted engine.
  - Chombo's mapped AMR (known by summary only): a large dependency with its own fourth-order operators.
  - Moving-mesh r-adaptivity (Tang and Tang 2003; applied to shock-induced combustion by Yuan and Tang 2007, abstracts only): the cell count is fixed, the smallest cell sets everyone's step (no subcycling), the wall and the flame compete for the same nodes, and every remap adds diffusion.
  - A cell-based quadtree: it loses the structured kernels.
- **Verification, criteria stated before code:**
  1. Free-stream preservation on a refined curved mesh, to round-off.
  2. Budgets closed to round-off across coarse-fine faces under subcycling.
  3. A moving shock and a contact through a patch boundary, converging at the scheme's order, with no spurious reflection above a stated level.
  4. C1 on a base mesh plus levels against C1 on the uniform finest mesh (c*, Isp, light-off).

**Order of work.** Each item exits with its evidence.
1. Table A, and its measured speedup on the C1 cold start. **Done 5 October, 22:40** (`docs/evidence/TABLE_A.md`, all four criteria pass). On 64x12, the C1 cold start to 8 ms takes 163 s on 5 threads against Cantera's 2,048 s, a 12.6 times speedup. It reaches full thrust (4 ms) in 82 s, or 84 s on 1 thread. c* is within 0.0053% of Cantera's, Isp within the printed digit and vacuum thrust within 0.0006%. The flow step is now 85% to 90% of the step.
2. Thread the flow step. **Done 5 October, 23:25** (`docs/evidence/THREAD_POOL.md`). One persistent pool runs the flow step and the reaction call. The result is the same bit for bit on any thread count and equal to the serial code's. The C1 64x12 cold start reaches full thrust (4 ms) in 28 s on 4 threads (77 s before) and 20 s on 6; to 8 ms, 56 s and 41 s. Transport (`core/transport.cpp`) is still serial; threading it is the lever for the turbulent cases.
3. Wall functions, with the verification above.
4. Table B, the finite-rate manifold, through the PaSR closure.
5. AMR.
6. Measure an RL10-like chamber's cold start. Apply the decision rule for implicit time, and the transport-table and GPU levers.

### Dated plan to January

Ben, 4 October 2026: all four regimes work by January. Weeks start Monday. Each milestone exits only with its verification evidence committed and the earlier regression cases still passing; a slipped milestone is reported with re-estimated dates at once, so Ben can choose between the date and the scope of the later regimes. Acceptance standards do not move.

**Revised 5 October 2026, 21:45** (reason: Ben's direction in *Lightweight engine*). T and L are added. Dual time stepping is parked, so its week in C4 goes to L. Wall functions move into C2. End dates from C4 on are unchanged. Risks (guessed): AMR may take two weeks rather than one, and Table B's ignition behaviour is untested; either would slip C3 and C4 by a week, reported at once.

| Dates | Milestone | Done means |
|---|---|---|
| 5–11 Oct | **C1** Time-evolving chamber and nozzle | Contour from data; gas supply faces with valve ramps; igniter; local-equilibrium option. Premixed gaseous H2/O2 "valves open, ignite" run settles; its end state matches CEA after the declared corrections |
| 6–11 Oct | **T** Tables and threads (added 5 October, 21:45) | Table A (equilibrium) verified against direct Cantera, with a measured speedup on the C1 cold start; the flow step threaded; Table B (finite rate) designed, then verified against CVODES |
| 12–18 Oct | **C2** Transport, walls, mixing closure | Mixture-averaged molecular transport verified against Cantera and FreeFlame; SST-2003 URANS and PaSR (PaSR acting on the tabulated source); no-slip and thermal walls through wall functions, verified against the wall-resolved SST; interleaved rings mix through the closure |
| 19–25 Oct | **L** Lightweight engine (added 5 October, 21:45) | Block-structured AMR verified to its four criteria; a cold start to full thrust in about 5 minutes on the Mac, measured on C1 and on an RL10-like chamber; the decision rule for implicit time applied |
| 26 Oct–1 Nov | **C3** Liquid LOX | Parcels, breakup, evaporation and two-way coupling verified (step 5); a subcritical LOX/GH2 component case |
| 2–8 Nov | **C4** RL10A-3-3A | Inputs extracted, prediction pre-registered in git, run, comparison table and loss ledger. The chemical milestone. Revised on 5 October at 21:45: the morning's one-week slip for dual time stepping became milestone L, so C3 moved a week later and C4 keeps the slipped end date |
| 9–15 Nov | **App draft** | "Draw and run" and "Probe and plot" on the chemical engine |
| 16 Nov–6 Dec | **M** Magnetic nozzle | The parked spike rejoins the shared engine: resistive MHD with constrained transport on the body-fitted mesh, applied coil fields and a declared vacuum-region treatment, field stresses in the thrust ledger. The zero-field gas cases are unchanged; MHD Riemann and magnetic-nozzle limit checks pass; a published magnetic-nozzle measurement is compared under the blind protocol |
| 7–20 Dec | **P** ICF post-burn pulse | A declared supplied pulse (mass, species, energy partition; from a published source, or labelled synthetic) expands through a pulsed magnetic nozzle; integrated impulse per pulse and coupling efficiency with the coil energy in the budget; compared with a published pulsed-nozzle result |
| 21 Dec–3 Jan | **A** Antimatter | Annihilation products from published yields as particle populations in the shared container: relativistic pusher in the coil field, decay, escape and deposition into the medium, thrust from escaping momentum plus coil reaction; orbit checks and comparison with the published Geant4 antimatter-nozzle study |
| by 3 Jan | **All four** | One regression suite runs all four regimes in the same engine and app. The week added on 5 October (first dual time, now L) used this milestone's buffer week, so the suite now runs inside A's last week, with no buffer left |

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
