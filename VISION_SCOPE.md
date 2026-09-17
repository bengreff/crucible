# CRUCIBLE — Vision and scope

Authoritative pivot, 17 September 2026. This document replaces the previous vision, proposals, implementation roadmap, and architecture rules. Current user instructions take precedence. Everything under `archive/` is historical reference, not an active requirement. This is the intended product and scope, not a claim of implemented or validated capabilities.

**The problem and purpose**

The problem is determining which engine arrangements can turn a reaction's available energy into useful propulsion, and what physically limits their performance. Geometry, reaction rates, confinement, and losses affect one another. The project investigates the operating conditions that emerge from those interactions and how changing a design changes what is achievable.

CRUCIBLE will be a technical sandbox that a human can use to construct, operate, observe, and compare such experiments. The user describes the device and its supplies or controls. The simulator calculates the medium's coupled behavior and makes the causes of the resulting performance inspectable.

The project is the reusable instrument. The minimum size of an antimatter engine, comparisons of magnetic nozzles, and the effect of pulse timing are possible investigations, not permanent definitions of the project. Neither broad applicability nor an attractive interface establishes scientific novelty by itself; an early comparison with existing tools must identify the useful experiments and human workflow this project adds.

**Who it is for**

The intended user has an intuition for physics and is willing to learn the instrument for a few hours, but is not required to be a numerical-methods specialist. The target is approachable use with visible technical meaning. A guided example should explain the difference between a control, an initial condition, a calculated quantity, and a model assumption.

The normal workflow is visual: open or create an experiment; shape its axisymmetric geometry; place supplies and equipment; choose supported materials and physical models; specify operation; run; inspect; change something; compare. Everyday use must not require code edits, handwritten solver configuration, or knowledge of mesh algorithms. Expert configuration and batch runs may expose the same underlying experiment.

**The physical scope**

The simulated subject is the working medium and its reactions. The combustion chamber's reacting contents are included. Bulk flow, species/composition, appropriate temperatures, reactions, energy transport, and relevant electromagnetic response evolve together. Energetic particles and radiation are represented separately when treating them as local heat would lose important behavior.

The spatial scope is two-dimensional axisymmetry. Vector components around the axis may be retained, but azimuthally varying structures are not resolved. Turbulent transport and other unresolved processes use explicit, responsive models rather than resolved eddies. Axisymmetry and model applicability remain limits on conclusions about real devices.

The initial bulk-medium target is nonrelativistic gas/plasma; energetic products can require relativistic particle transport. This does not automatically support every relativistic plasma or arbitrary material phase. The first chemical implementation will use gaseous or already-vaporized reactants, with modeled mixing, combustion, and expansion. Liquid breakup and evaporation require a supported extension before raw liquid injection is claimed.

Generalization means rearranging supported physics through geometry and operating inputs. A new arrangement must not need private engine-specific physics code. A new mechanism or an unsupported regime may require new models. Different physical descriptions are allowed where justified; numerical smoothness alone does not prove a physically valid transition.

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

For ICF, implosion and production of the supplied pulse are outside the simulated time window. The supplied state must describe mass, composition, spatial state, motion, and energy partition—not only total yield. Choose a handoff after the excluded burn and before significant omitted nozzle interaction, within the receiving model's validity. If no defensible handoff exists, the experiment needs compatible upstream output or is unsupported. Remaining reactions cannot silently be discarded. Energy already escaped must be accounted for without depositing it again.

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
2. Establish a small human workflow through a supported reacting-gas experiment: define, run, inspect, change, compare. A prepared operating state is acceptable; complete cold startup of all equipment is not required.
3. Exercise magnetic-plasma expansion and then energetic-product coupling early enough to test the shared design before polishing the chemical application. Distinguish synthetic or prescribed demonstrations from validated device predictions.
4. Verify and validate the supported cases, measure computational cost, and test usability with unfamiliar users. Reuse archived code only after its assumptions and behavior fit the new model.
5. Conduct substantive investigations chosen by Ben using the same instrument. Use the findings to improve the supported models and observation tools without turning each engine into its own implementation.

Success requires both scientific credibility and human usability. The primary demonstrations must use a shared experiment workflow, changes in geometry/operation must produce inspectable physical consequences, and at least one investigation must explain a useful mechanism or limit with its assumptions. The schedule and detailed acceptance tolerances follow the initial model and compute feasibility work. A solver-only deliverable or a visual shell without supported physics is incomplete.

**Authority and changes**

The pre-pivot repository is preserved under `archive/pre-pivot-2026-09-17/`. Consult `REBUILD_NOTES.md` for lessons, specific references, and unresolved issues. No old milestone, prohibition on visualization, all-regime promise, mandatory three-dimensional simulation, source-code language split, or universal determinism requirement carries over unless restated here.

Keep this vision short enough for a human to understand and challenge. Record consequential scope changes here with their reasons. Implementation documents explain how the current scope is achieved; they do not silently enlarge it.
