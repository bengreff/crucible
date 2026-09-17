# CRUCIBLE — proposal foundation

Working draft, 16 September 2026. Recommended scope for discussion; this does not replace the existing plans or claim that the capabilities below are implemented. Based on the vision and high-level planning material and the current conversation, without an implementation review.

**The project in ordinary language**

I want to experiment with engines whose operation depends on hot gas or plasma interacting with walls, magnetic fields, and reactions. I will build a virtual laboratory where I describe an arrangement and how it is supplied or driven, then calculate what happens inside it and what performance results. The purpose is to investigate unfamiliar arrangements without having to assume the internal conditions I am trying to discover.

The simulator follows both directions of cause and effect: reactions change the material, and the changing material affects reactions; fields influence the plasma, and plasma currents influence the fields. Some experiments begin with a supplied energy pulse. Others calculate energy release as part of this evolution. The distinction is explicit in the experiment.

The project is the reusable laboratory. Questions such as the minimum size of an antimatter engine are investigations conducted with it, not its permanent organizing question.

**What counts as an arrangement**

The first version represents a cross-section revolved around an axis. Chambers, nozzles, coils, and supply openings share that symmetry. Flow and fields may have components around the axis, but the calculated state has no variation around it. Discrete injectors or nonuniform coils require a justified symmetric representation; drawing a symmetric outline does not establish that real operation remains symmetric.

An arrangement consists of stationary walls, coils, material supplies, and open boundaries surrounding gas or plasma. It can contain multiple connected regions. It can operate continuously or through time-dependent inputs and pulses. Coil currents are prescribed initially; the plasma's electromagnetic response is calculated where the supported model requires it. Predicting the response of an external power circuit is a separate capability.

The starting physical scope is nonrelativistic bulk gas/plasma with separately represented energetic particles, including relativistic reaction products, and radiation where these transport significant energy or momentum. Direct particle exhaust can be supported through that particle description; a relativistic bulk pair plasma is not automatically covered by the same fluid model.

An accepted arrangement must use supported materials, reactions, and transport models within their stated domains. Changing dimensions, supply rates, or coil placement should be an experiment configuration. Introducing a new physical mechanism requires a model extension. This is the boundary of generality.

**What makes an experiment complete**

Every experiment specifies four things:

| Part | What must be supplied | Example |
|---|---|---|
| Arrangement | Geometry, materials, coils, and connections to the surroundings | Chamber contour, coil positions, inlet and exhaust openings |
| Preparation and controls | Initial state and externally imposed inputs, including their timing | Initial density and temperatures, fuel supply, current history |
| Physical description | Included reactions, energy transport, unresolved-transport models, and prescribed processes | Annihilation products transported; turbulent transport modeled |
| Measurements | Quantities to compare and any definition of successful operation | Thrust and heat load at a fixed feed rate; impulse from one pulse |

An energy pulse needs a location, spatial extent, duration or initial post-pulse state, affected mass, composition, and energy partition. Total energy alone does not specify the resulting physical experiment.

An experiment need not begin from cold startup. It can ask about an operating state or a prepared transient. Its conclusion is conditional on that preparation; it does not also establish startability.

**The rule that keeps the physics connected**

For each influential process, the experiment uses one of three explicit treatments:

| Treatment | Meaning | Example |
|---|---|---|
| Calculated evolution | The process changes with the evolving state | Fuel depletion and reaction rate respond to local conditions |
| Responsive approximation | An unresolved process is replaced by a model that responds to relevant conditions | Heat transport depends on temperatures, gradients, and magnetic conditions |
| Prescribed condition | The process is an external input or a deliberate experimental assumption | A supplied post-fusion state; heating imposed independently of plasma response |

These treatments apply to processes, not labels such as “chemical” or “fusion.” A single experiment can contain all three. They are not necessarily interchangeable accuracy settings.

The decision rule is: if a response could materially change the requested measurement, represent that response or state that the answer is conditional on fixing it. A prescribed pulse can answer what happens after a burn. It cannot answer how the surrounding field changes that burn. A fusion-core experiment must connect reaction rates, product deposition, temperature, density, and losses if it claims to predict fusion power.

Modeling turbulence preserves only the feedback represented by the chosen model. It does not confer an ability to predict every confinement transition or instability. Where that missing behavior could decide the result, the answer is unresolved by the present model.

**What the calculation actually does**

The calculation evolves the bulk material, electromagnetic fields, and populations that carry energy independently of the bulk. These descriptions exchange energy, momentum, and charge consistently. Reactions consume and produce species; their released energy is accounted for consistently with the reaction model. Energy carried by products is not also deposited immediately as heat.

For example, annihilation products may heat the plasma downstream from where they were born. The plasma then expands, changing both its interaction with the field and the environment encountered by later products. The model must preserve this chain when it affects the experiment.

One physical transfer has one accounting path, regardless of the application. This does not require one equation for every material state. A fluid description and a particle description can coexist if their responsibilities, exchanges, and domains are explicit. Hybrid fluid/kinetic magnetic-nozzle research provides an example of why such distinctions are physically useful. [Fluid–kinetic nozzle model](https://e-archivo.uc3m.es/entities/publication/c6518003-4bea-4f7b-af8f-bd051dcf771b)

The field model must represent the electric and magnetic effects needed for the selected experiment; following magnetic field lines alone is not sufficient. Choosing an adequate fluid/kinetic treatment for nozzle expansion is an early research and feasibility task. Electron cooling in nearly collisionless nozzles is a documented modeling difficulty. [Experiments and modeling](https://arxiv.org/abs/2212.07161)

Existing libraries can supply interaction models and reference calculations. They do not automatically supply the coupled plasma model. Geant4, for example, provides particle transport through matter and electromagnetic fields; deciding how that transport exchanges energy and momentum with an evolving plasma remains part of this project. [Geant4 transportation](https://geant4.web.cern.ch/documentation/pipelines/master/prm_html/PhysicsReferenceManual/generalities/particletransport/transportation.html)

**Where the physical system ends**

An inlet represents a connection to a specified supply. A prescribed coil current represents an external driver. A thermal boundary represents a stated surrounding temperature, cooling response, or heat-removal capacity. These are physical experimental conditions, not unmentioned infinite resources.

Wall or coil temperature can be calculated with an appropriate thermal-response model when it affects the answer. Calculated heat and force loads can also be compared with specified equipment limits. A result obtained with fixed wall temperature is conditional on maintaining that temperature. If cooling capacity controls the proposed engine size, it must enter the size analysis.

This first scope does not include arbitrary moving machinery, detailed injector atomization, pellet implosion, resolved turbulence, or unrestricted three-dimensional stability. These require additional capabilities or explicit prepared/input conditions. Mechanical construction, antimatter production/storage, and whole-vehicle operation are not predicted by a successful core-flow calculation.

**What an experiment returns**

The result is the evolving state together with measurements: forces on the device, exhaust mass and momentum, energy deposition, radiation escape, and thermal loads. Steady cases report thrust and consistently defined specific impulse; pulses report impulse and the time history. Average performance requires a specified repetition schedule and treatment of any interaction between pulses.

Forces must account for transfer to walls and coils as applicable, and be consistent with momentum carried by matter and radiation and stored within the domain. A visually narrow plume alone is not a thrust measurement.

Each comparison states its prescribed assumptions, applicable model limits, and numerical and model sensitivity. Unsupported extrapolation is identified. No universal accuracy percentage or hardware “works” verdict is promised. A geometry ranking is a result only when the modeled difference can be distinguished from relevant uncertainty.

**Three demonstrations test the shared scope**

| Demonstration | Shared capability it exercises | What would show useful generality |
|---|---|---|
| Chemical nozzle | Compressible expansion, thermochemistry, wall interaction, force measurement | Changing a contour changes spatial flow and performance without new engine-specific physics code |
| Prescribed post-fusion expansion | Transient plasma/field interaction and impulse transfer | Changing source extent or coil arrangement changes capture and losses through the same state evolution |
| Antimatter interaction region | State-dependent reaction and product transport coupled to plasma/fields | Changing supply or geometry changes deposition and exhaust without prescribing an energy-utilization efficiency |

These are development demonstrations, not three required full engine designs. Chemical validation checks relevant flow and thermochemistry; it does not validate plasma confinement or annihilation coupling. Magnetic-nozzle benchmarks and particle-transport references are needed for those capabilities.

The proposed minimal experiments are concrete:

1. **Two solid nozzles, one supply.** Supply the same specified hot combustion-product mixture at the same inlet stagnation pressure and temperature to two contours with the same throat and exit areas. Hold the wall treatment and ambient conditions fixed. Calculate flow, mass flow, force, and exit energy partition. Compare the ideal limiting case with a reference expansion calculation, then assess contour-dependent effects against appropriate nozzle data or independently established solutions. This tests whether geometry produces a meaningful difference while upstream combustion is a declared input.
2. **One prepared plasma pulse, two coil arrangements.** Place the same axisymmetric volume of initially ionized hydrogen, with specified density, electron/ion temperatures, and initial velocity, inside two coil arrangements. Evolve the plasma and its electromagnetic response, measuring impulse on the device, escaping energy, and peak loads. This is the post-event expansion test; the first version need not pretend to represent a particular fusion pellet. A later post-fusion source must supply the appropriate composition and energy partition. Compare elementary expansion/field limits and then a suitable published or experimental magnetic-nozzle case before interpreting novel arrangements.
3. **One antimatter supply, different receiving plasma conditions.** For the provisional proton–antiproton case, specify the antiproton beam distribution and supply history, hydrogen supply/preparation, walls, and coil controls. Change the initial hydrogen density while holding the other selected controls fixed. Calculate annihilation location and rate, product transport, deposition, plasma response, and force. Compare a dilute limit against a particle-transport reference, then establish the additional evidence needed for the coupled case. The intended demonstration is a deposition fraction and flow response that emerge together, rather than an efficiency selected in advance.

The first antimatter species and engine configuration must be selected explicitly during model specification. Proton–antiproton heating of a secondary plasma is a reasonable initial candidate, not a commitment that every antimatter configuration is already covered. A particle-emission test alone does not complete this demonstration; the requested plasma feedback must also be exercised.

Magnetic confinement fusion remains a stretch capability. An axisymmetric confinement experiment belongs when its geometry, reaction feedback, and transport can be represented adequately. General tokamak performance prediction is not implied by adding a fusion reaction rate.

**How to build without returning to the old scope**

First, specify one minimal experiment in each primary demonstration, including measurements and evidence needed to trust them. Identify the common physical descriptions and the assumptions that need separate support. Select the plasma and energetic-particle models before committing to implementation details.

Next, establish the connected path from supplied matter or energy through state evolution to measured force. Begin with tractable reference cases and introduce magnetic response and nonlocal energy transport in tested increments. Exercise different configurations early enough to expose hidden application-specific assumptions.

Then demonstrate feedback explicitly: vary the plasma conditions and check the predicted change in deposition; vary deposition and check the resulting plasma response. Check conservation, numerical convergence, and relevant analytic, experimental, or independently computed references. A different numerical model is a cross-check, not experimental validation.

Finally, use the instrument for an exploratory campaign chosen by Ben. Its result should identify a mechanism, limit, or design relationship and explain the dependence on assumptions. The campaign can change without changing the laboratory's identity.

Success is a reusable, evidenced capability to conduct the stated class of experiments, demonstrated across the three primary cases and through at least one substantive investigation. The exact schedule follows a feasibility check of the hardest coupled model and available compute; existing repository claims are not treated as independently verified here.

**What could distinguish this project**

The intended contribution is making coupled virtual experiments accessible across this bounded class of devices: geometry and controls are changed while important internal relationships remain calculated. Reuse matters because it permits new questions with comparable measurements and visible assumptions.

Neither axisymmetry nor combining established models is sufficient evidence of originality. A targeted literature comparison must establish which experiments existing tools already support, what must presently be supplied by hand, and what this instrument would add. That comparison, model selection, and selection of specific reference cases are the next planning deliverables, before a detailed implementation roadmap.
