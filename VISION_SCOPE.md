# CRUCIBLE — Vision and scope

Authoritative pivot, 17 September 2026; amended 4 October 2026 (see *Changes since the pivot* at the end). This document replaces the previous vision, proposals, implementation roadmap, and architecture rules. Current user instructions take precedence. Everything under `archive/` is historical reference, not an active requirement. This is the intended product and scope, not a claim of implemented or validated capabilities.

**The problem and purpose**

The problem is determining which engine arrangements can turn a reaction's available energy into useful propulsion, and what physically limits their performance. Geometry, reaction rates, confinement, and losses affect one another. The project investigates the operating conditions that emerge from those interactions and how changing a design changes what is achievable.

CRUCIBLE will be a technical sandbox that a human can use to construct, operate, observe, and compare such experiments. The user describes the device and its supplies or controls. The simulator calculates the medium's coupled behavior and makes the causes of the resulting performance inspectable.

There is one simulation engine. Every added physical capability becomes part of its shared calculation and is available to any compatible experiment. Chemical, fusion and antimatter arrangements are applications of that engine, not separate solvers or installable extensions. Internal numerical methods and physical closures can differ where the physics requires them; their exchanges, control history and measurements remain coupled consistently.

The project is the reusable instrument. The minimum size of an antimatter engine, comparisons of magnetic nozzles, and the effect of pulse timing are possible investigations, not permanent definitions of the project. Neither broad applicability nor an attractive interface establishes scientific novelty by itself; an early comparison with existing tools must identify the useful experiments and human workflow this project adds.

**Who it is for**

The intended user has an intuition for physics and is willing to learn the instrument for a few hours, but is not required to be a numerical-methods specialist. The target is approachable use with visible technical meaning. A guided example should explain the difference between a control, an initial condition, a calculated quantity, and a model assumption.

The normal workflow is visual: open or create an experiment; shape its axisymmetric geometry; place supplies and equipment; choose supported materials and physical models; specify operation; run; inspect; change something; compare. Everyday use must not require code edits, handwritten solver configuration, or knowledge of mesh algorithms. Expert configuration and batch runs may expose the same underlying experiment.

The product is one native desktop app with local simulation on Mac, Windows, and Linux. Editing, operation, and observation belong to the same open experiment. The user operates an instrument rather than submitting reports. Feed, heating, and current can change live through supported equipment controls; changing geometry requires pausing and restarting with an explicit initial state. Controls and viewing must stay responsive even when physical results take minutes.

Ben intends to publish the tool. His Mac is the main development/use machine, with a Windows PC available for larger calculations. Design for multiple hardware classes and meshes with millions of cells. The performance ambition is short, physically credible runs in minutes; actual cell counts, physical durations and achieved accuracy must be established by measurement and convergence. High fidelity means adequacy for the claimed observables within the supported models, not simply a large mesh.

Ordinary setup should center on chamber/nozzle shape, supplies, injection conditions and equipment controls. Supported injector and wall descriptions provide explicit, inspectable assumptions where full hardware detail is unnecessary. Composition, supply thermodynamic state, spatial delivery and thermal boundary behavior still matter when they affect the answer; an easy interface must not silently invent them.

**The physical scope**

The simulated subject is the working medium and its reactions. The combustion chamber's reacting contents are included. Bulk flow, species/composition, appropriate temperatures, reactions, energy transport, and relevant electromagnetic response evolve together. Energetic particles and radiation are represented separately when treating them as local heat would lose important behavior.

The spatial scope is two-dimensional axisymmetry. Vector components around the axis may be retained, but azimuthally varying structures are not resolved. Turbulent transport and other unresolved processes use explicit, responsive models rather than resolved eddies. Axisymmetry and model applicability remain limits on conclusions about real devices.

The unifying principle is consistent physical exchanges and feedback, not one universal turbulence formula. The display distinguishes modeled averages from resolved motion. Results should reveal whether a design comparison survives plausible uncertainty in mixing, heat transport and losses. A minimum size inferred from this model is conditional on its scope; it does not by itself establish stability, manufacturability or a complete working engine.

The initial bulk-medium target is nonrelativistic gas/plasma; energetic products can require relativistic particle transport. This does not automatically support every relativistic plasma or arbitrary material phase. Liquid propellant enters as liquid: breakup and evaporation are part of the shared engine and are added before the first liquid-fed engine is run (Ben, 4 October 2026). No liquid-fed engine is run as pre-vaporized gas.

Generalization means rearranging supported physics through geometry and operating inputs. A new arrangement must not need private engine-specific physics code. A new mechanism or an unsupported regime may require new models. Different physical descriptions are allowed where justified; numerical smoothness alone does not prove a physically valid transition.

Fields and radiation may extend beyond the material, and expanding matter may leave a fluid model's domain of validity. Account for those exchanges and transitions explicitly. Boundary-object classification does not justify omitting a physical effect that controls the requested result.

**The system boundary**

Boundary objects are anything outside the working medium and reaction: injector hardware, igniters, coils and their drivers, walls and cooling equipment, upstream supplies, and the surrounding environment. They specify physical exchanges and, when necessary, responses to the medium. A boundary object is not necessarily a fixed source or a geometrical exterior surface.

For example, an injector supplies composition, energy, momentum, and a spatial flow distribution; its delivered flow can respond to chamber pressure. An igniter supplies a bounded energy deposit or a defined hot-gas supply. Coils impose operating conditions while the plasma's electromagnetic response remains part of the calculation. Walls exchange heat and momentum; thermal or circuit response models may retain internal state where that feedback matters. Maintaining a prescribed temperature, current, or flow is an explicit assumption about the external equipment.

Approximating turbulence, chemistry, or particle deposition is an internal physical model, not a way to move those processes outside the system. Included does not mean microscopically resolved. If an omitted response could control the requested measurement, represent it, establish its insignificance, or report the result as conditional or unresolved.

**Primary applications**

| Application | Inputs and modeled subject | Intended questions |
|---|---|---|
| Chemical propulsion | Supplies and igniter interfaces; reacting chamber, mixing, modeled combustion, and nozzle expansion | How chamber/nozzle arrangement and operation affect the evolving state and propulsion performance |
| Antimatter propulsion | Specified matter and antimatter supply; supported annihilation reactions, product transport, deposition, and coupled medium/field response | How geometry and operation affect usable energy, force, losses, and physical limits |
| Post-ICF propulsion | A justified prepared post-burn state or time-dependent event output; subsequent plasma/field interaction | How a prescribed event becomes impulse, exhaust, and loads |

For ICF and other supplied pellet events such as the AMCF events in Ben's proposal, implosion and production of the supplied pulse are outside the simulated time window. The supplied state must describe mass, composition, spatial state, motion, and energy partition—not only total yield. Choose a handoff after the excluded burn and before significant omitted nozzle interaction, within the receiving model's validity. If no defensible handoff exists, the experiment needs compatible upstream output or is unsupported. Remaining reactions cannot silently be discarded. Energy already escaped must be accounted for without depositing it again.

The first antimatter species and physical configuration must be selected during model design; support for one does not establish support for all. Interaction libraries supply constituent models, not automatic validation of the coupled engine.

Magnetic confinement fusion is a stretch application where it fits the physical descriptions. Fuel delivery/deposition, electron and ion energy, reaction rates, product heating, depletion, and transport must respond to one another if fusion power is predicted. A geometry-specific transport model is not a universal confinement model. General tokamak performance and three-dimensional stability prediction are not initial commitments.

**What the user sees and learns**

Visualization is a central product requirement and a scientific observation tool. It develops alongside the solver, beginning with the first useful experiment. It is not deferred until every physics feature is complete.

The interface must make these actions possible:

- Construct or modify a cross-section and see what physical volume it represents. Any revolved 3D view is labeled as an axisymmetric representation.
- See supplies, equipment, initial conditions, and assumptions separately from calculated fields; inspect units and the meaning of controls.
- Observe spatial fields, flow and field directions, reactions, energy deposition, and losses when supported. Use probes and histories to connect local changes with thrust, impulse, and loads.
- Pause and inspect a running calculation; scrub saved time states; save, reopen, and compare experiments using consistent scales. Scrubbing recorded output is distinct from rerunning the physics.
- Understand a refused input, numerical failure, unmet operating target, or unsupported physical regime through a localized explanation and retained diagnostic state.
- Export the experiment, measurements, model/data versions, and figures so another person can reproduce and assess the result.

Changes during operation must be recorded as controls or new experiment versions; the interface must not silently rewrite the inputs behind an existing result. Example cases and presets are transparent starting points, not hidden sources of desired performance. Validation status and uncertainty are accessible where they affect interpretation, with detail available on demand.

The usability target is that, after a few hours of guided learning, a person with basic physical intuition can modify an example, run it, inspect the response, make a controlled comparison, and explain an important assumption without developer intervention. Test that with people unfamiliar with the code. Attractive plots alone do not meet the target.

**Scientific and computational commitments**

Each physical transfer has one accounting path. Sources, exchanges, stored energy/momentum, and escaping fluxes must balance within justified numerical tolerances. Reaction energy and product deposition must not be counted twice. Force includes the relevant transfer to walls and coils and is checked against a consistent momentum budget.

Use analytical/reference problems, convergence studies, and appropriate experimental comparisons. Distinguish verification of numerics, comparison between models, and validation against measurements. Retain input and data provenance. Report numerical sensitivity and model uncertainty honestly; there is no universal promised error bar, pedigree score that substitutes for evidence, or unrestricted hardware-feasibility verdict.

Results must be reproducible to declared numerical tolerances. Bitwise identity can be useful for debugging but is not a universal product requirement. Choose the implementation language, data format, numerical schemes, and CPU/GPU split for demonstrated physical and usability needs, not inherited doctrine. Fast interaction with the editor and results is required; real-time physical simulation is not. Long runs need clear progress, cancellation, and recoverable saved state.

**How the rebuild proceeds**

1. Specify representative experiments and the shared physical state; compare existing capabilities and choose defensible fluid, particle, reaction, and unresolved-transport models. Use the archived lessons as evidence, not requirements.
2. Make chemical propulsion the first substantive milestone: an operable reacting chamber/nozzle, the visual define/run/inspect/change/compare workflow, and validation using published geometry, operating conditions and measured performance for multiple real engines. The first engine is RL10A-3-3A from its published geometry; then many engines (see *Changes since the pivot*). Component experiments support particular model checks. The medium always evolves in time, and a chemical run starts "valves open, ignite": the chamber holds ambient gas, the propellant valves open on a declared schedule, the igniter fires, and the flow goes through its transient to steady operation. A declared partly developed state that still evolves in time is only a development fallback, with the reason stated. Valve travel is a declared input; pumps, valves and other equipment are not simulated.
3. Add magnetic-plasma evolution and energetic-product coupling to the same engine after establishing the chemical milestone. Ben's target is all four regimes working by January 2027 (see *Changes since the pivot*). Use small early architecture checks where they prevent costly redesign, without displacing the chemical-validation priority. Distinguish synthetic or prescribed demonstrations from validated device predictions.
4. Verify and validate the supported cases, measure computational cost, and test usability with unfamiliar users. Reuse archived code only after its assumptions and behavior fit the new model.
5. Conduct substantive investigations chosen by Ben using the same instrument. Use the findings to improve the supported models and observation tools without turning each engine into its own implementation.

Success requires both scientific credibility and human usability. The primary demonstrations must use a shared experiment workflow, changes in geometry/operation must produce inspectable physical consequences, and at least one investigation must explain a useful mechanism or limit with its assumptions. The schedule and detailed acceptance tolerances follow the initial model and compute feasibility work. A solver-only deliverable or a visual shell without supported physics is incomplete.

The primary applications describe the intended reach of the instrument, not simultaneous first-release coverage. Chemical propulsion with real-data validation is the first demonstration. Ben's proposal is complete, and the desired development horizon is a few months. Treat that as a planning target to test against measured progress, not evidence that all proposed investigations can be completed in that time. Later physical capabilities belong to the same engine and must earn their credibility while preserving earlier results within declared tolerances.

**Authority and changes**

The pre-pivot repository is preserved under `archive/pre-pivot-2026-09-17/`. Lessons, specific references, and unresolved issues from it are in [RESEARCH.md section 7](RESEARCH.md#7-lessons-from-the-pre-pivot-attempt). No old milestone, prohibition on visualization, all-regime promise, mandatory three-dimensional simulation, source-code language split, or universal determinism requirement carries over unless restated here.

**Changes since the pivot**

4 October 2026, Ben:

- **One project.** "I want to make sure there is ONE project, coherent vision." The post-pivot work is the main line; the 17 September pivot is the scope; the old 3-D plan is retired.
- **Chemistry baseline.** Shifting equilibrium is the validated chemistry baseline, produced by the general reacting machinery in the engine (finite-rate kinetics, with local equilibrium as its fast-chemistry limit), not by a separate equilibrium calculator, "because this needs to be the baseline for all future regimes".
- **Never a static solution.** "It cannot be pure static equilibrium, the gas inside the chamber must evolve over time, ideally it actually starts up from cold but the combustion mixing just gets the same approximation as everything else ... it cannot be a static image because engines are not static and especially when we get to fusion everything will be evolving over various timescales." Mixing uses the same unresolved-transport closure as the rest of the flow, not a special combustion shortcut.
- **Validation order.** RL10A-3-3A from published geometry next, matching vacuum Isp and thrust. Then CEA sweeps for trends. Then a database of hundreds of chemical engines from public sources, matching CEA and measured Isp/thrust and their trends.
- **Evaporation first.** "Add evaporation first": liquid LOX injection with breakup and evaporation is built into the shared engine before RL10. There is no pre-vaporized RL10.
- **Valves open, ignite.** The start of a chemical run: an empty chamber at ambient, the propellant valves open, the igniter fires, and the flow goes through its transient to steady state.
- **All four regimes by January.** Chemical, magnetic nozzle, pulsed fusion (an ICF post-burn pulse) and antimatter all work in the same engine by January 2027. Chemical stays lean and correct; the parked MHD spike comes back right after it. TECHNICAL_PLAN.md holds the order and dates. This restates the all-regime goal for this window; it does not revive the pre-pivot promises.
- **App draft after chemical.** The first draft of the laboratory sandbox app: "Draw and run" (edit the cross-section, place injectors, run, watch the fields evolve) and "Probe and plot" (click anywhere for time histories; thrust and Isp traces).

Keep this vision short enough for a human to understand and challenge. Record consequential scope changes here with their reasons. [RESEARCH.md](RESEARCH.md) records evidence and open model decisions; [TECHNICAL_PLAN.md](TECHNICAL_PLAN.md) defines the implementation architecture. They explain how the current scope is achieved and do not silently enlarge it.
