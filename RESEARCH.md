# CRUCIBLE — Research and evidence plan

17 September 2026; updated 4 October 2026. This document supports [VISION_SCOPE.md](VISION_SCOPE.md); [TECHNICAL_PLAN.md](TECHNICAL_PLAN.md) turns the decisions into an implementation plan. Sources below were consulted for this plan. A library capability is not evidence that CRUCIBLE implements or validates it. Proposed experiments and acceptance criteria below are project decisions, not claims from the cited sources.

**Decision status.** The product requirements are commitments. The algorithm candidates below are a starting research program, not established capability. A model becomes supported only after its equations, data, validity limits, numerical checks and relevant comparisons are recorded. Distinguish verification, experimental validation and conditional exploration in the app; validation applies to particular observables and regimes, not an entire application name.

Ben has selected chemical propulsion with real-data validation as the first substantive milestone. There is one engine: every added capability enters its shared calculation and is usable in any compatible experiment. Differences in evidence or applicability are properties of physical models and results, not separate supported/exploratory products or installable extensions. The proposal is complete; the target development horizon is a few months.

## 1. What we need to establish

Can one approachable instrument let a person investigate how geometry, operating inputs, reacting matter, and fields together determine useful propulsion? The proposed contribution is a shared, inspectable experimental workflow with credible coupling. Neither assembling libraries nor offering several engine presets proves originality.

Compare representative tasks against existing tools before claiming a gap:

| Existing capability | What to learn or reuse | What CRUCIBLE must demonstrate beyond it |
|---|---|---|
| [NASA CEA](https://github.com/nasa/CEA): equilibrium composition and idealized rocket performance | Thermochemical limiting cases and consistent reference assumptions | Spatial, time-dependent consequences of chamber/nozzle geometry and operating changes |
| [PeleC](https://github.com/Pele-Suite/PeleC): compressible reacting flow on adaptive meshes | Numerical methods and practical reacting-flow infrastructure | A usable experimental instrument that also supports justified plasma/product coupling |
| [Geant4 antimatter nozzle study](https://arxiv.org/abs/1205.2281): particle tracking and magnetic-nozzle optimization | A precedent and particle-level comparison problem | Feedback between the changing medium, reaction/product transport, fields, and operating inputs where it matters |

Reproduce a small comparable task, record what already works, and identify the additional question our shared workflow actually answers. Do not describe magnetic-nozzle simulation itself as new. A minimum-size antimatter engine study is one possible use of the instrument, not its definition.

Two claims in the proposal need qualification in subsequent explanations. Spatial chemical-engine simulation is not limited to turbulence-resolving cluster calculations: published axisymmetric averaged-flow work already exists. Magnetic nozzles for pulsed fusion also have computational and experimental precedents, including [NASA's pulsed-fusion nozzle program](https://www.nasa.gov/directorates/stmd/space-tech-research-grants/experimental-and-computational-validation-and-scaling-of-power-generating-magnetic-nozzles-for-pulsed-fusion-propulsion/). The specific added coupling, accessible experimental workflow and resulting investigations must establish the contribution. Do not claim that an entire application has never been simulated without a narrower literature finding.

## 2. Product decisions and their consequences

Ben has selected local operation on Mac, Windows, and Linux; live changes to feed, heating, and current; pause and restart for geometry edits; responsive controls and viewing even when physical results take minutes.

Use a native C++ app with Qt Widgets and VTK. Qt supplies established desktop controls and layout facilities; VTK has a direct Qt widget integration. This is a practical choice for editing a cross-section, selecting physical objects, probing fields, and arranging comparison views. It is not a claim that native software automatically calculates faster. [Qt Widgets](https://doc.qt.io/qt-6/qtwidgets-index.html), [Qt main windows](https://doc.qt.io/qt-6/qmainwindow.html), [VTK integration](https://vtk.org/doc/nightly/html/classQVTKOpenGLNativeWidget.html).

Start with a portable CPU calculation. Optional acceleration follows measurement. AMReX documents CPU execution and CUDA/HIP/SYCL GPU paths; that is not a supported Apple Metal solver backend. Mac-local operation must therefore succeed without assuming GPU acceleration. [AMReX execution support](https://amrex-codes.github.io/amrex/docs_html/GPU.html).

The first usability experiment is concrete: a newcomer opens an example, predicts what increasing feed will do, changes it, locates the actual change on the time history, inspects the response, and compares a geometry revision. Record misunderstandings as well as task completion. Smooth animation is insufficient if users confuse calculated pressure with prescribed pressure.

## 3. Preserve feedback without resolving everything

Every model must state its inputs, evolving state, outputs, validity range, and uncertain parameters. The common rule is physical exchange and accounting, not one universal equation set.

For a reacting chamber, temperature and composition determine reaction rates; reactions change composition and energy; pressure and motion then alter supply, mixing, and subsequent reaction. A prescribed igniter contributes external energy to that loop. For an ICF handoff, the supplied event defines an initial state, after which the included interactions evolve. Both fit the same experiment without pretending that both energy sources are independent.

Approximations must respond to relevant local conditions. A model that gets thrust right through compensating errors is inadequate if it gets wall loads or energy partition wrong. Select observations that expose those differences.

### Reacting gas

Begin with a small gaseous hydrogen/oxygen example, with a documented mechanism and validity range. Use a prepared reacting state if needed to separate the first flow/coupling tests from ignition and flame stabilization. This is a development case, not an RL10 reproduction or the project's permanent fuel choice.

The primary validation targets are multiple real engines with published chamber/nozzle geometry, operating conditions and measured thrust/Isp or equivalent performance data. Ben chose RL10A-3-3A as the first engine on 4 October 2026 (NASA TM-107318 gives area against axial station; inputs and the pre-registered prediction live in `docs/validation/`). The TUM single-element GOX/GCH4 chamber, chosen on 30 September, stays a component case; its inputs are incomplete (`docs/validation/VALIDATION_TUM_ROUND.md`). After RL10 come CEA sweeps for trends and then a database of hundreds of engines from public sources. Audit phase/injection assumptions against the initial gaseous-medium support; a liquid-fed engine is not validated by silently treating liquid injection as an equivalent gas supply. Record missing information and whether its uncertainty could control the intended comparison before selecting each case.

The Penn State preburner combustor (RCM-1) is only a possible component benchmark for wall heat flux, not the primary engine-performance target. A [DLR study](https://elib.dlr.de/99741/1/jpp%20zhukov%202015.pdf) demonstrates axisymmetric SST-based modeling on a workstation and discusses uncertainty/inconsistency in supplied inlet conditions versus measured pressure. It is evidence that this kind of reduced simulation already exists and that boundary-data auditing matters. If used, recover the original geometry, preburner-stream compositions/enthalpies, wall conditions, measurements and uncertainties. Do not replace the streams with pure cold hydrogen/oxygen or tune the unknown inputs solely to obtain the expected result.

Use complementary evidence for distinct observables: flow/flame profiles for mixing and chemistry, measured chamber wall heat for heat transfer, and a separately audited nozzle/engine dataset for thrust or discharge performance. RL10A-3-3A is the first thrust/Isp dataset. A successful RCM-1 heat-flux comparison does not validate thrust/Isp. Compare chemical equilibrium and idealized nozzle limits against CEA for verification/reference purposes, keeping experimental validation distinct.

Test the desired low-detail setup explicitly. Begin with chamber/nozzle contours, propellant composition and inlet state, delivery rates/distribution, ambient pressure and a declared wall model. Use supported abstract injector descriptions rather than demanding resolved injector hardware. Vary uncertain mixing and heat-transfer inputs; if geometry/performance conclusions depend strongly on them, report that dependence and seek the missing data. Do not calibrate a private efficiency factor for every engine to make a common simulator appear predictive.

[Cantera](https://www.cantera.org/stable/cxx/index.html) supplies thermodynamic, reaction-rate, and molecular transport calculations. It does not by itself supply resolved injector mixing or a turbulent combustion closure. First verify laminar chemistry/transport; then choose an averaged mixing and combustion model against an appropriate dataset before making turbulent-chamber claims. Test ignition delay, flame propagation, species, wall heat, and thrust as applicable; a single equilibrium temperature cannot validate a combustion chamber.

### Shifting-equilibrium baseline and the loss ledger (Ben, 4 October 2026)

Shifting equilibrium is the validated chemistry baseline, and it must come out of the general reacting machinery: finite-rate kinetics integrated in time, with local equilibrium as its fast-chemistry limit. Local equilibrium is an option for verification, not the production model. The reasoning is that later regimes (plasma, fusion products) have no equilibrium calculator to fall back on, so the chemical case must prove the time-evolving machinery itself.

The evidence chain, each step a measured delta in the same engine:

1. The engine's thermochemistry reproduces the CEA ideal rocket (done: 18 points, worst 0.131%, `docs/evidence/THERMO_VERIFICATION.md`).
2. A time-evolving chamber and nozzle in the local-equilibrium limit, run to its steady end state, reproduces CEA shifting-equilibrium c* and Isp once the differences in problem definition are accounted for: finite chamber area (CEA's finite-area-combustor option), two-dimensional nozzle divergence, and the injected enthalpy.
3. Finite-rate kinetics instead of local equilibrium: the difference is the kinetic (recombination-freezing) loss. It must lie between the shifting and frozen CEA limits.
4. Mixing through the unresolved-transport closure, then wall friction and heat loss: each adds its loss.
5. The real engine: the remaining gap to measurement is explained, not fitted.

This mirrors the classic loss breakdown used for liquid rocket performance prediction (equilibrium ideal, then divergence, kinetics, boundary layer and energy-release efficiency, as in the JANNAF performance prediction methodology and two-dimensional kinetics codes). The difference is that every term here comes from one time-dependent field calculation rather than separate correction codes. A match on Isp alone does not validate the terms individually; report each one.

Trends matter as much as values. Ben's medium-term target is "pasting in hundreds of chemical engine params and matching CEA/real values and trends": the same engine, with no per-engine efficiency factor, should reproduce how Isp and thrust move with O/F, chamber pressure, area ratio and propellant.

### How unresolved transport will be approximated

The unifying idea is to calculate the **net transport caused by motion we do not resolve**, using the evolving medium as input. The same accounting rules apply across applications; the physical formulas need not be identical. There are three distinct tasks:

| Task | First candidate | What it actually estimates |
|---|---|---|
| Turbulent gas motion | Time-dependent, density-weighted averaged flow (Favre URANS) with a specified SST k–omega model | Unresolved momentum transport, with modeled heat/species transport |
| Turbulence affecting chemical reaction | Finite-rate chemistry with a partially stirred reactor (PaSR) closure | How limited small-scale mixing changes the average reaction rate |
| Transport in magnetized plasma | Classical anisotropic transport within its domain, plus a separately identified anomalous-transport hypothesis if needed | Heat, particles and momentum moving along/across the field; classical collisions and turbulence are different contributions |

**Gas turbulence.** Use SST as the first candidate for the chamber/nozzle because wall layers and separation matter. Its two extra fields describe turbulent kinetic energy, k, and a characteristic dissipation frequency, omega; they evolve with the flow and determine an effective viscosity. Specify one published variant, initially SST-2003, rather than mixing constants and corrections from different implementations. The [Turbulence Modeling Resource](https://tmbwg.github.io/turbmodels/sst.html) documents those distinctions. This is our candidate selection, not evidence that SST is best for every injector jet.

Model turbulent heat and species fluxes using explicit turbulent Prandtl and Schmidt numbers and test sensitivity to them. Supply inlet turbulence intensity and length scale with provenance. Initially integrate the mean wall layer down to the wall with an appropriate near-wall mesh; this does not resolve turbulent eddies. Verify wall distance, near-wall resolution and heat transfer. Wall functions can be a later, separately tested approximation, not an automatic cure for a coarse mesh. Account consistently for energy stored in unresolved motion, its production and its conversion to heat; [compressible RANS energy conventions](https://tmbwg.github.io/turbmodels/implementrans.html) matter here.

**Turbulent combustion.** Averaging a nonlinear reaction rate is not equivalent to evaluating it at average temperature/composition. More turbulent diffusion alone does not fix that. The initial PaSR candidate compares a mixing time with a chemical time to estimate the reacting fraction. One documented formulation weights the finite-rate source by `tau_chem / (tau_chem + tau_mix)`; both time definitions and their coefficients are part of the model. [OpenFOAM PaSR reference implementation](https://api.openfoam.com/2306/PaSR_8C_source.html). This is an algorithm reference, not a decision to import OpenFOAM code or add that dependency.

PaSR is a testable starting point for a gaseous mixing-controlled example, not a promise of correct flame speed, ignition, extinction or reignition. Compare scalar and velocity profiles using a documented dataset such as the [Sandia/ETH hydrogen jet flames](https://tnfworkshop.org/data-archives/simplejet/h2he/), with the actual oxidizer, pressure and boundary conditions. That dataset does not validate a high-pressure hydrogen/oxygen rocket. If PaSR cannot reproduce the needed reaction/mixing behavior, assess EDC or an appropriate flamelet/PDF treatment for the target regime before adding features. Do not introduce several combustion closures at once or tune only the final thrust.

For the first PaSR test, use the reference mixing-time form `C_mix * sqrt((nu + nu_t) / epsilon)`, with dissipation derived consistently from the selected SST model, and pin a documented chemical-time definition. The mixing coefficient, species used in the chemical-time estimate, near-zero-rate handling and pressure range are explicit research choices. Test the fast-mixing limit, reaction-free mixing, and zero-turbulence behavior; do not obtain a laminar limit by an undocumented numerical floor.

**Plasma transport.** Braginskii-type coefficients are a first classical baseline for sufficiently collisional, near-Maxwellian, weakly coupled ionized plasma. They distinguish transport directions relative to the magnetic field and electron/ion behavior; they are not a turbulence model. [Classical transport documentation](https://docs.plasmapy.org/en/stable/formulary/braginskii.html). Check the selected formula's composition and magnetization assumptions. Heat-flux limiting can bound a fluid approximation, but does not restore missing nonlocal electron physics; compare with a justified nonlocal or kinetic treatment if that controls the answer. [Nonlocal transport research](https://www.osti.gov/servlets/purl/1305830).

For an exploratory torch, expose a small, named family of additional cross-field particle and thermal transport laws, with explicit coefficients and domains. These may respond to local density, temperatures, field and gradients, but that responsiveness alone does not validate them. Keep particle diffusivity, electron/ion thermal diffusivity, momentum transport and resistivity distinct. Do not silently add a gas eddy viscosity or universal Bohm multiplier to every plasma. First ask whether the conclusion survives a physically justified range of transport assumptions. A surviving trend is a conditional finding; if the trend reverses, determining the transport becomes necessary research.

MCF transport models require their own geometry and equilibrium assumptions. For example, [TGLF](https://gafusion.github.io/doc/tglf.html) combines linear mode calculations with a saturation model and has tokamak-specific development/validation. It is a reference for what an informed reduced model entails, not a ready-made closure for an open magnetic nozzle. No MCF turbulence package is selected yet.

The displayed URANS fields are evolving averages, not a movie of real turbulent eddies. Neither finer axisymmetric cells nor an added diffusion term recovers missing azimuthal instabilities, discrete injector asymmetries, or three-dimensional magnetic disruption. Such omissions must remain visible when they could determine a proposed limit.

### Magnetized plasma: the earliest scientific decision

Use collisional, conducting-fluid expansion as a bounded first plasma case. Assess separate electron/ion energies, anisotropic transport, ionization, and nonideal field response against the chosen conditions. Do not freeze a universal ideal-MHD model into the architecture.

An expanding magnetic nozzle may become weakly collisional. Electron distributions and cooling then affect acceleration and detachment; a conventional fluid closure is not automatically adequate. This is a documented issue in [magnetic-nozzle electron-cooling research](https://arxiv.org/abs/2212.07161).

Before extending to a torch, map representative density, temperature, field, and length scales. Compare mean free paths, gyroradii, collision times, and flow times. Decide whether the target experiment is covered by a fluid description with justified closures or needs a hybrid kinetic description. Record the equations and the evidence behind that choice. Use [PlasmaPy's formulary](https://docs.plasmapy.org/en/stable/formulary/index.html) for independent scale/formula checks and selected [OPEN-ADAS data](https://open.adas.ac.uk/) where applicable; neither supplies a universal plasma model.

*Parked with the MHD spike (Ben, 30 September 2026) until the chemical milestone is done; the two decisions below stand for when it resumes.*

**Director default, 30 September 2026 (reversible by Ben): resistive MHD first, Hall next.** The first field physics will be resistive MHD, with a conductivity and magnetic Reynolds number taken from the plasma state. Ideal MHD freezes the field into the flow, so plasma can never detach from a magnetic nozzle's field lines; detachment is the central question for a magnetic nozzle, so ideal MHD cannot be the target model (the spike also measured the coil field being swept out under ideal MHD, `docs/parked/MHD_SPIKE.md`). The Hall term is designed in as the next increment, because most electric-propulsion regimes need it; the induction equation and its boundary conditions should leave room for it. Neither model resolves the weakly collisional expansion described above. An outside check on the model class has been requested by the Director. Ben may reverse this default.

**Director decision, 30 September 2026 (reversible by Ben): Maeno et al. 2013 uses a prescribed plume.** The first validation case (`docs/parked/VALIDATION_MAENO2013.md`) starts from a supplied plume state. The laser-ablation event lies outside the simulated window, and VISION_SCOPE treats such an event as a supplied state. The plume's mass, velocity distribution, ionisation and temperature are constrained from independent published ablation data, not from Maeno's own fitted model. Each input is labelled measured, derived or assumed. Plausible ranges are propagated through the run. No value is picked because it matches the impulse. Consequence: the 1/2/3-omega wavelength trend becomes an input to the case, not a prediction.

**Expansion into vacuum is a model decision.** A conducting fluid cannot be extended into empty space merely by maintaining an arbitrary minimum density. Test whether the measurement can be taken before fluid validity fails, or whether a kinetic/hybrid continuation and conservative coupling are needed. Magnetic fields also occupy regions outside the material; external field geometry, plasma response and open boundaries need a coherent treatment. Wall sheaths and ambipolar escape may require reduced interface models even when their microscopic scales are not resolved. A model of turbulence cannot substitute for these effects.

The early plasma feasibility case must examine the intended measurement region and exit conditions, not only an easy dense interior. Otherwise it cannot answer whether the shared design can predict useful exhaust or pulse impulse. Record the parameter envelope for that case before choosing the fluid or hybrid implementation.

### Antimatter and energetic products

First select the reactants and target state explicitly. An electron/positron example cannot establish antiproton-engine support. The initial research recommendation is a hydrogen/antiproton case because it exercises charged and neutral products and coupling; this remains a model-selection decision until its interaction data and compute cost have been assessed.

Keep direct beamed-core exhaust distinct from a thermalized plasma-core torch when defining an experiment. The former can obtain thrust directly from guided relativistic products; the latter relies substantially on energy transferred to a working medium. The shared engine should represent the relevant populations and transfers, including intermediate arrangements, rather than silently treating these as equivalent configurations.

[Geant4 physics lists](https://geant4-internal.web.cern.ch/support/physics_lists) are constituent models requiring validation for the application. Evaluate annihilation channels, spectra, decay, magnetic deflection, material interactions, and escape independently. Establish which interactions are valid for the hot, ionized medium; add justified plasma stopping/collision models where necessary. Never substitute a cold material solely because it has the desired density.

Separate event generation from subsequent transport conceptually, even if Geant4 handles both in a supported case. Deposit each transfer once, update the receiving medium, and remove the corresponding particle energy/momentum. Annihilating reactants must also be depleted consistently. Account for escaping photons, neutral products, and neutrinos where relevant.

Unresolved turbulent transport may affect residence time, mixing, and confinement in a torch. It cannot be dismissed because particle trajectories are calculated. Vary plausible transport models/coefficients and compare the resulting loss and deposition budgets. If the proposed conclusion changes, transport evidence is necessary before that conclusion is supported.

**Radiation is an explicit part of the energy budget.** First estimate the importance of thermal/atomic emission, energetic photons, absorption and re-emission for the selected regime. An optically thin loss model is a useful starting point only where reabsorption is negligible; it still needs valid emissivity and ionization data. If absorption or radiation pressure matters, select a transport model and opacity data before claiming that regime. Diffusion and flux-limited diffusion are not universal solutions for directed, weakly interacting radiation; [Castro's radiation documentation](https://amrex-astro.github.io/Castro/docs/radiation.html) illustrates a different, explicitly specified radiation-hydrodynamic approximation. We have not selected a general radiation solver. Photons represented by particle transport must not also be counted as a separate escaped thermal loss.

### Post-ICF and MCF

There is no universal number of nanoseconds at which to inject a post-ICF state. A usable handoff is after the excluded implosion/burn and before significant omitted interaction with the engine, with the supplied matter inside the receiving model's validity. Required data include mass/species, spatial density, velocity, thermal partition, energetic-product distributions, timing, and already-escaped energy. Yield alone is insufficient. Locate a published or collaborator-provided dataset satisfying those requirements; no such dataset has yet been selected here. Until then, synthetic states are numerical demonstrations, labeled accordingly.

MCF remains a stretch. First consider an open axisymmetric configuration; fueling and deposition, thermal transport, fusion rates, product heating, depletion, and exhaust must form a closed feedback calculation. A prescribed temperature profile tests a reaction calculation, not predictive confinement. Select a geometry-specific transport model and evidence before adding an MCF preset.

## 4. Numerical and library feasibility experiments

| Decision | Small experiment needed | Consequence |
|---|---|---|
| AMReX as spatial infrastructure | Axisymmetric volumes/areas, curved wall, near-axis behavior, small cut cells, conservation, and required magnetic-field operators in one small problem | Adopt only the demonstrated combination; a list of individually supported features is insufficient |
| Stiff integration | Couple compression and chemical reaction; compare reduced timesteps and reference chemistry; exercise step rejection | Select the SUNDIALS method and coupling order using measured error and cost |
| Particle coupling | Known orbit/decay plus conservative deposition into a changing medium | Determine sampling noise, coupling interval, and whether batching preserves feedback |
| Native portability | Same small numerical case and field view on all three desktop platforms | Establish the dependency/compiler combination before broad implementation |
| Averaged gas closure | Wall-flow and jet mixing checks, then measured reacting-jet profiles with held-out conditions | Determine whether SST plus PaSR covers the intended chemical case |
| Plasma closure and losses | Collisionality/scale map, conduction and escape checks, radiation significance estimate | Bound the first plasma model and identify any required kinetic/radiative capability |

AMReX's [embedded-boundary facilities](https://amrex-codes.github.io/amrex/docs_html/EB_Chapter.html) provide geometry infrastructure, not every application discretization. Its [time-integration interface](https://amrex-codes.github.io/amrex/docs_html/TimeIntegration_Chapter.html) connects to SUNDIALS, whose [ARKODE methods](https://sundials.readthedocs.io/en/latest/arkode/Introduction_link.html) include explicit, implicit, and mixed approaches. These reduce infrastructure work; they do not settle our conservation, field, or coupling scheme.

Start on a fixed grid. Adaptive refinement and multiple accelerator backends are later decisions. If the AMReX geometry/field combination fails, investigate a fixed body-fitted finite-volume mesh for the supported geometry family and revise the technical decision before expanding either implementation.

## 5. Evidence that earns a scientific claim

Follow a hierarchy: mathematical/analytical checks, solution convergence, independent implementation comparisons, and experiments with documented measurement uncertainty. These answer different questions. NASA's [verification and validation tutorial](https://www.grc.nasa.gov/www/wind/valid/tutorial/tutorial.html) and [validation assessment](https://www.grc.nasa.gov/www/wind/valid/tutorial/valassess) support this building-block approach.

For each claimed capability retain: the question and regime, input/model/data versions, reference source, measured observables, numerical sensitivity, uncertain assumptions, and conclusion with its limits. These are developer evidence and inspectable experiment metadata; the normal user workflow does not submit reports.

Initial evidence ladder:

1. Axisymmetric uniform/rest states, shocks, advection, and known nozzle limits; quantify conservation and near-axis error.
2. Constant-volume/pressure chemistry against a reference; then reacting flow with compression and diffusion. Distinguish chemistry agreement from experimental validation.
3. Magnetic-field identities, waves/diffusion appropriate to the selected equations, particle orbits, and force/energy transfer.
4. Integrated chamber and magnetic-expansion cases, measuring local profiles and loss partition as well as propulsion.
5. An independently documented experimental comparison for each supported physical regime before presenting it as validated.

Use spatial and temporal refinement separately, along with particle-count/coupling-interval refinement where applicable. Estimate convergence only when the asymptotic behavior supports it; three grids alone do not prove accuracy. [NASA spatial convergence guidance](https://www.grc.nasa.gov/www/wind/valid/tutorial/spatconv.html).

Choose tolerances from the observable and intended conclusion before examining pass/fail outcomes. Explore uncertain transport, wall, supply, and source parameters separately from mesh error. A nominal optimum smaller than those uncertainties is not an established optimum.

Keep calibration data separate from evaluation data. Start with a small designed parameter sweep; extend to global sensitivity analysis if interactions prevent interpretation. Display uncertainty in observations as well as assumptions: for example, compare nozzle shapes across the same plausible transport range and identify whether their performance ordering persists. Numerical refinement, Monte Carlo sampling error and physical-model uncertainty are separate questions.

The first chemical closure comparison will use TMR wall/jet cases and the TNF hydrogen-flame data above where their actual conditions fit. Specific magnetic-expansion, antimatter-interaction and post-ICF datasets remain to be selected. Whole-device experimental validation may be unavailable for speculative engines; constituent validation and conditional predictions are still valuable, but must not be described as validation of a complete torch.

## 6. What to resolve first

The next implementation should establish the native app, local CPU execution, safe live controls, and a chemical case with audited experimental data. Chemical validation is the first substantive result. Small geometry/field and state-representation checks should prevent architectural dead ends while that work proceeds; substantial plasma development follows the chemical milestone. All later physics enters the same engine and is checked against the established chemical cases.

Unresolved items are specific: a compatible dependency build on each platform; the curved-wall/field discretization; the first plasma closure and vacuum interface; turbulent combustion evidence; radiative losses; antimatter species/interaction validity; and a usable ICF handoff dataset. None is resolved by adding another generic plugin interface. Update this document with measured decisions as they become available.

The foundation is ready for bounded feasibility work, not a commitment that all primary applications will be predictive in the first release. Before promising a delivery date, select a representative spatial resolution, physical duration, chemistry size and particle count; measure their cost on named local hardware. The first demonstration is chemical validation, and the requested horizon is a few months. Re-estimate the staged schedule in [TECHNICAL_PLAN.md](TECHNICAL_PLAN.md) from measured progress rather than relaxing the evidence requirements.

Computational feasibility must include million-cell workloads on the main Mac and later the Windows PC, not only small demonstrations. Measure memory, time per accepted step and time per simulated physical interval for transport alone, chemistry, then their coupling. Establish the physical interval from the observable: gas residence/acoustic behavior, ignition, equipment heating and repeated-pulse recovery may have different timescales. Steady or prepared initial states can shorten an investigation only when they preserve its question. The minutes-per-short-run target remains to be demonstrated on stated hardware and model settings.

Translate the proposal's later investigations into measurable questions. Pulse repetition requires repeated events with residual plasma, heating, coil/circuit response and declared recovery constraints; a single-pulse impulse does not determine a maximum rate. Compare engine types under stated input-power, supplied-mass, geometry and equipment constraints. A minimum size is conditional on those constraints and model validity. An MCF shape search needs a defined parameter family, objective, transport assumptions and relevant stability limits; finding the best sampled candidate is not a comprehensive search of every physical design. These remain investigations enabled by the instrument, not separate implementations or promised discoveries.

## Implemented evidence boundary (4 October 2026)

Verification evidence for what is implemented is indexed in [docs/evidence/README.md](docs/evidence/README.md): the axisymmetric gas core and its device-thrust ledger, thermochemistry against CEA, stiff reaction integration against Cantera, the thermally perfect mixture core, Strang-coupled reacting flow (a detonation against the equilibrium Hugoniot), and the low-Mach study. None of it validates a turbulence or mixing closure, a plasma model, or any real engine. The million-cell throughput number in the gas-core record is ten timesteps of single-gas flow, not evidence that a reacting run completes in minutes.

## 7. Lessons from the pre-pivot attempt

Folded in from the former `REBUILD_NOTES.md` on 4 October 2026, text unchanged except where marked.

17 September 2026. This is a selective extraction from the old project, not a renewed implementation plan. `VISION_SCOPE.md` is authoritative. Archive-relative references below point into `archive/pre-pivot-2026-09-17/`; their contents were preserved unchanged. Line numbers refer to that frozen snapshot.

The extraction inspected repository structure, selected numerical routines, certificate machinery, and relevant session records. It did not independently rerun or validate the legacy solver. Statements about historical results below are attributed to those records.

### Lessons worth carrying forward

1. **Build the human experiment loop early.** The previous plan explicitly deferred visualization until a large certification run. That conflicts with the new product. The new instrument must let a person see the arrangement, observe the calculated state, ask where energy goes, and compare a change. An inaccessible solver is an incomplete sandbox. This is a product decision from the pivot, not a measured failure of an interface that did not yet exist.

2. **Keep the physical scope distinct from resolution.** The reacting chamber belongs inside the simulator. Turbulence and chemistry can be represented by responsive models. This does not justify assuming away their influence, nor does inclusion require resolving every eddy or collision. Avoid inheriting the old escalation from a physical question to mandatory three-dimensional startup simulation.

3. **A good headline number can hide the wrong internal behavior.** The [session record](archive/pre-pivot-2026-09-17/SESSION_LOG.md#L703) reports that adding/fixing transport changed scored performance by at most about 0.32% in one coarse study while wall heat changed substantially. An apparent 3.27% heat change became 19.8% after correcting a driving potential; the earlier small difference involved cancellation. Inspect fields, local budgets, heat loads, and intermediate quantities, not only thrust/Isp. This is a reason visualization belongs in the scientific workflow.

4. **Identify exactly what a certificate proves.** The [RL10 certificate generator](archive/pre-pivot-2026-09-17/crates/engine/src/bin/station5_rl10_certificate.rs#L1) explicitly recomputes scoring from recorded readouts rather than rerunning the underlying flow. Its constants include different model assumptions and calibrated versus prior-band cases. A green certificate diff establishes that score generation is unchanged; it does not establish that the current executable reproduces the run or that combustion startup is validated.

5. **A numerically completed run is not a physical feasibility verdict.** The [later GPU startup record](archive/pre-pivot-2026-09-17/SESSION_LOG.md#L2240) reports a conservation-audited trajectory that ignited a small sustained kernel but reached only about 22.5 kN against a 75.6 kN target. The log attributes the shortfall to the flame-spreading/modeling gap at that setup. Preserve this as a limitation of that run and model, not evidence that the real engine cannot work. The new interface must distinguish numerical failure, an unmet target, and missing physical capability.

6. **Different processes settle on different clocks.** The [wall-response discussion](archive/pre-pivot-2026-09-17/SESSION_LOG.md#L719) compares an approximately 37 ms wall thermal timescale with an 11 ms run. Settled gas did not mean a settled wall. Show what has equilibrated, what is still changing, and why a run ended.

7. **Geometry and conservation need shared definitions.** [Grid metrics](archive/pre-pivot-2026-09-17/crates/grid/src/lib.rs#L1026) assign each shared face one radius calculation. The comments record that algebraically equivalent radius formulas rounded differently and broke exact cancellation. Preserve this ownership principle in axisymmetric volumes, face areas, source terms, and diagnostics; do not inherit the entire adaptive-azimuthal grid just to keep it.

8. **Table validity and thermodynamic consistency are part of the physics.** Keep explicit units, data provenance, interpolation coordinates, and applicability checks. The [old blend model](archive/pre-pivot-2026-09-17/crates/solvers/src/euler/blend_eos.rs#L1) contains branch thresholds, pressure-root selection, and warm-start behavior that deserve a fresh physical review before reuse. Continuous interpolation or selecting the first root is not, by itself, proof that a closure describes the new medium correctly.

9. **Measure the expensive operation.** The [GPU profiling record](archive/pre-pivot-2026-09-17/SESSION_LOG.md#L2263) reports that reaction-node projections consumed 79% of kernel time in that case; flux sweeps were about 5%. Reusing pressure hints reduced the measured composed step from roughly 0.137 to 0.028 seconds. These are historical case-specific timings. The reusable lesson is to profile a representative coupled experiment rather than optimize an assumed bottleneck.

10. **Save enough state to explain and reproduce a run.** Configuration and table identities, model version, controls, field state, numerical caches needed for continuation, and diagnostic histories all matter. Preserve the restart concept. Replace opaque failure labels with useful explanations and retain the last inspectable state. Reproducibility does not require universal cross-device bit identity.

**Specific code to consult**

These are references for extraction or redesign, not approved drop-in components. Read their tests and assumptions before carrying code into the new implementation.

| Need | Specific archived chunk | Evidence and caution |
|---|---|---|
| Axisymmetric geometric metrics | [grid/src/lib.rs](archive/pre-pivot-2026-09-17/crates/grid/src/lib.rs#L1036): `face_radius`, `ring_radii`, `cell_volume`, `face_area_r`, `face_area_z` | Shared-face ownership and cylindrical measures. Extract from the larger 3D/brick architecture. |
| Compressible-gas numerical flux | [euler/hllc.rs](archive/pre-pivot-2026-09-17/crates/solvers/src/euler/hllc.rs#L18): `physical_flux`, `hllc_flux`; [recon.rs](archive/pre-pivot-2026-09-17/crates/solvers/src/euler/recon.rs#L34): `ppm_faces` | Review EOS and fixed state-component assumptions. This is gas-dynamics code, not an electromagnetic/plasma solver. |
| Primitive-state recovery from thermochemistry | [euler/table_eos.rs](archive/pre-pivot-2026-09-17/crates/solvers/src/euler/table_eos.rs#L119): `TableEos::bind`, `cons_from_phz` | Useful lookup/projection pattern; medium coordinates and admissible states must fit the new model. |
| Fast checked table access | [tables/src/bound.rs](archive/pre-pivot-2026-09-17/crates/tables/src/bound.rs#L177): `BoundColumn::interpolate`; [interp.rs](archive/pre-pivot-2026-09-17/crates/tables/src/interp.rs#L60): reference API | Checks table domain and narrower validity envelope, then interpolates in declared coordinates. Do not assume table interpolation error equals physical-model uncertainty. |
| Content-addressed data | [tables/src/digest.rs](archive/pre-pivot-2026-09-17/crates/tables/src/digest.rs#L1): v3 encoding; [Python tables.py](archive/pre-pivot-2026-09-17/offline/crucible_offl/tables.py#L139): `content_digest` | Length-prefixed metadata, units, data, and provenance; check [cross-language fixture tests](archive/pre-pivot-2026-09-17/crates/tables/tests/fnd5_python_seam.rs#L26). Old logs mention earlier digest versions; consult current archived code. |
| Chemistry-library interface | [chemistry.py](archive/pre-pivot-2026-09-17/offline/crucible_offl/chemistry.py#L135): `EquilibriumEngine`, `state_php`, `chamber`, `performance`; `FrozenReactantEngine` at line 404 | Explicit SI conversions and reactant/equilibrium distinction. Equilibrium performance is a reference, not validation of mixing or finite-rate combustion. |
| Audit of coupled transfers | [sdc.rs](archive/pre-pivot-2026-09-17/crates/solvers/src/sdc.rs#L1928): `audit_stored`, `audit_check`, `audit_check_n` | Stored changes compared with transfers using the integration weights actually applied. Extend accounting for fields, reaction products, and external work; do not inherit all SDC interfaces. |
| Save/resume | [gpu_engine_run.rs](archive/pre-pivot-2026-09-17/crates/gpu/src/gpu_engine_run.rs#L94): `write_checkpoint`, `read_checkpoint` | Includes physical state, hints, clock, trackers, and input digests. Temporary-write/rename is useful; inspect crash-durability and format compatibility requirements afresh. |
| Measurements and spatial export | [engine/run.rs](archive/pre-pivot-2026-09-17/crates/engine/src/run.rs#L1515): `plane_area`, `plane_mdot`, `plane_thrust`; `fields_csv` at line 1716 | Useful reference for measurements and first field views. Gas-plane thrust alone does not supply complete magnetic, radiation, or transient force accounting. Avoid coupling the new UI to the old run controller. |
| Small numerical reference problems | [Sod tests](archive/pre-pivot-2026-09-17/crates/solvers/tests/solv1_station1_sod.rs#L19); [manufactured-solution tests](archive/pre-pivot-2026-09-17/crates/solvers/tests/solv1_euler_mms.rs#L11); [table holdouts](archive/pre-pivot-2026-09-17/offline/tests/test_station3_surfaces.py#L59) | Preserve the physical questions and independently expected answers. Reconsider exact-output and architecture-specific assertions when the model changes. |
| GPU implementation references | [gpu/src/engine_host.rs](archive/pre-pivot-2026-09-17/crates/gpu/src/engine_host.rs); [residency_diffusion.cu](archive/pre-pivot-2026-09-17/crates/gpu/cuda/residency_diffusion.cu); [residency_engine_xcheck.rs](archive/pre-pivot-2026-09-17/crates/gpu/src/residency_engine_xcheck.rs) | Device-resident state and CPU/device comparisons are worth studying. These depend on the old state layout and coupling; they are not the foundation of the new architecture by default. |

### Three small patterns, extracted verbatim

Shared radial-face arithmetic, from `crates/grid/src/lib.rs:1036`:

```rust
pub fn face_radius(&self, f: usize) -> f64 {
    self.spec.r_min + f as f64 * self.spec.dr
}
```

The value is that neighboring cells refer to this same definition, not that the formula is elaborate. Preserve that invariant when redesigning the grid.

Unambiguous string hashing, from `crates/tables/src/digest.rs`:

```rust
fn put_str(h: &mut Sha256, s: &str) {
    h.update((s.len() as u64).to_le_bytes());
    h.update(s.as_bytes());
}
```

This avoids using a separator that can also appear inside metadata. The full schema and ordering must remain explicit; copying this helper alone does not define a reproducible file format.

The energy flux in the old gas solver, from `crates/solvers/src/euler/hllc.rs:28`:

```rust
f[n] += p;
f[I_EN] = u_n * (e_tot + p);
```

This compactly exposes pressure's contribution to momentum and energy transport. Preserve the governing accounting when redesigning interfaces; do not reuse this gas-only expression as the total energy flux for electromagnetic fields and energetic particles.

### Design decisions to retain, reconsider, and retire

Retain explicit units and material/data provenance; one owner for each exchange; conservative spatial measures; model-domain checks; independently specified reference problems; reproducible experiments; useful saved state; and performance work guided by measurements. Expose these through understandable displays and diagnostics rather than requiring users to read internal logs.

Reconsider the state representation, EOS/reaction closure, grid and time integrator, input format, data storage, and CPU/GPU division against the supported physics and the human workflow. The archived Rust/Python/CUDA work is available, but the language split and HDF5-only interface are not governing constraints. Prefer a shared experiment and result representation for the interface and batch work so they cannot silently run different physics.

Retire mandatory full 3D, automatic axisymmetric-to-3D ambition, universal single-law blending, mandatory physical cold startup of an entire engine (equipment such as pumps; the chamber contents still evolve in time and a cold start of the medium is the aim, Ben 2026-10-04, VISION_SCOPE), numerical-verdict-as-real-feasibility language, unconditional certificate byte equality, per-document crate architecture, and the requirement to finish physics before visualization. Detailed structures, service lifetime, pumps, whole-vehicle simulation, and arbitrary liquid injection are not automatic obligations of the new medium/reaction scope.

### What remains genuinely unresolved

The plasma/fluid/particle treatment, appropriate transport closures, first antimatter species and configuration, defensible post-ICF data/handoff, and specific validation anchors need evidence. The sections above identify starting closure candidates and tests; [TECHNICAL_PLAN.md](TECHNICAL_PLAN.md) owns current implementation choices. This section remains the historical extraction. Chemical reacting-chamber models must be assessed in their own right; the old burn-progress closure is not accepted just because code exists. None of these is settled by writing the new vision.

The new instrument's usefulness also needs a literature comparison and human testing. The first implementation should let a person define, run, inspect, and compare a small supported experiment, while magnetic and particle cases test the physical design early. A chosen investigation can then test the tool's research value without becoming the only question the instrument can ask.

### Archive and reset verification

The archive manifest records the pre-pivot revision, branch, original paths, hashes, and modes for 236 tracked files and four untracked planning documents. The archive operation preserved them; no legacy source was rewritten. Ignored local artifacts were moved but are not included in Git. The prior active CI and agent-specific settings are archived, so they cannot impose the previous roadmap on the root project.

This reset was checked for archive integrity and active-document consistency. Legacy solver tests, physical benchmarks, and GPU runs were not repeated, because this change neither alters nor adopts that implementation. The old environment may need recreation after relocation. No claim of newly validated physics follows from the reset.
