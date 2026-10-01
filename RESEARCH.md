# CRUCIBLE — Research and evidence plan

17 September 2026. This document supports [VISION_SCOPE.md](VISION_SCOPE.md); [TECHNICAL_PLAN.md](TECHNICAL_PLAN.md) turns the decisions into an implementation plan. Sources below were consulted for this plan. A library capability is not evidence that CRUCIBLE implements or validates it. Proposed experiments and acceptance criteria below are project decisions, not claims from the cited sources.

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

The primary validation targets are multiple real engines with published chamber/nozzle geometry, operating conditions and measured thrust/Isp or equivalent performance data. RL10 remains a candidate based on prior interest, not an inherited acceptance target. Select well-documented engines first. Audit phase/injection assumptions against the initial gaseous-medium support; a liquid-fed engine is not validated by silently treating liquid injection as an equivalent gas supply. Record missing information and whether its uncertainty could control the intended comparison before selecting each case.

The Penn State preburner combustor (RCM-1) is only a possible component benchmark for wall heat flux, not the primary engine-performance target. A [DLR study](https://elib.dlr.de/99741/1/jpp%20zhukov%202015.pdf) demonstrates axisymmetric SST-based modeling on a workstation and discusses uncertainty/inconsistency in supplied inlet conditions versus measured pressure. It is evidence that this kind of reduced simulation already exists and that boundary-data auditing matters. If used, recover the original geometry, preburner-stream compositions/enthalpies, wall conditions, measurements and uncertainties. Do not replace the streams with pure cold hydrogen/oxygen or tune the unknown inputs solely to obtain the expected result.

Use complementary evidence for distinct observables: flow/flame profiles for mixing and chemistry, measured chamber wall heat for heat transfer, and a separately audited nozzle/engine dataset for thrust or discharge performance. No thrust dataset has yet been selected. A successful RCM-1 heat-flux comparison does not validate thrust/Isp. Compare chemical equilibrium and idealized nozzle limits against CEA for verification/reference purposes, keeping experimental validation distinct.

Test the desired low-detail setup explicitly. Begin with chamber/nozzle contours, propellant composition and inlet state, delivery rates/distribution, ambient pressure and a declared wall model. Use supported abstract injector descriptions rather than demanding resolved injector hardware. Vary uncertain mixing and heat-transfer inputs; if geometry/performance conclusions depend strongly on them, report that dependence and seek the missing data. Do not calibrate a private efficiency factor for every engine to make a common simulator appear predictive.

[Cantera](https://www.cantera.org/stable/cxx/index.html) supplies thermodynamic, reaction-rate, and molecular transport calculations. It does not by itself supply resolved injector mixing or a turbulent combustion closure. First verify laminar chemistry/transport; then choose an averaged mixing and combustion model against an appropriate dataset before making turbulent-chamber claims. Test ignition delay, flame propagation, species, wall heat, and thrust as applicable; a single equilibrium temperature cannot validate a combustion chamber.

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

**Director default, 30 September 2026 (reversible by Ben): resistive MHD first, Hall next.** The first field physics will be resistive MHD, with a conductivity and magnetic Reynolds number taken from the plasma state. Ideal MHD freezes the field into the flow, so plasma can never detach from a magnetic nozzle's field lines; detachment is the central question for a magnetic nozzle, so ideal MHD cannot be the target model (the spike also measured the coil field being swept out under ideal MHD, `docs/MHD_SPIKE.md`). The Hall term is designed in as the next increment, because most electric-propulsion regimes need it; the induction equation and its boundary conditions should leave room for it. Neither model resolves the weakly collisional expansion described above. An outside check on the model class has been requested by the Director. Ben may reverse this default.

**Director decision, 30 September 2026 (reversible by Ben): Maeno et al. 2013 uses a prescribed plume.** The first validation case (`docs/VALIDATION_MAENO2013.md`) starts from a supplied plume state. The laser-ablation event lies outside the simulated window, and VISION_SCOPE treats such an event as a supplied state. The plume's mass, velocity distribution, ionisation and temperature are constrained from independent published ablation data, not from Maeno's own fitted model. Each input is labelled measured, derived or assumed. Plausible ranges are propagated through the run. No value is picked because it matches the impulse. Consequence: the 1/2/3-omega wavelength trend becomes an input to the case, not a prediction.

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

## Implemented evidence boundary — first gas increment

The ideal-gas prototype now has analytic shock, geometric balance, nozzle limiting-case and live-pressure-response checks; see [the implementation record](docs/IMPLEMENTATION.md) for methods, errors and limitations. It does not validate any candidate turbulence/chemistry/plasma model or any published engine. Its million-cell throughput benchmark is only ten timesteps, not evidence of a useful physical simulation completing in minutes.
