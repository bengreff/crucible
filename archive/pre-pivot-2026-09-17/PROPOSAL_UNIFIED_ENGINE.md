# CRUCIBLE — Unified Engine Proposal

**Status:** DRAFT v0.3 for Ben's review, 2026-09-09. Written to stand alone: no prior document is
assumed. If adopted it becomes VISION_SCOPE v2.0 and replaces PLAN_CHEMICAL_SANDBOX.md; it amends
neither yet. Appendix C records what changed between versions and why, including the v0.2 review.

**One sentence.** CRUCIBLE is a computational test stand: you describe an engine's geometry and the
things you would physically control (coils, injectors, heating, electrodes, a spark), the instrument
marches the working medium forward in time under one set of physical laws, and it returns thrust,
specific impulse, efficiency, power flows, and heat loads, each as a range with a stated pedigree,
together with a picture of every field in every cell, for any engine whose working medium heats
itself by known laws.

**The honest headline.** For chemical engines the instrument predicts absolute performance blind
(demonstrated). For magnetic fusion engines it predicts **which control geometry is better, over a
stated fraction of the closure band**, and it computes the core temperature as an emergent quantity
with that band attached. It does not predict absolute fusion gain to better than the band, and it
says so on every result.

**What makes it unique.** Its primary outputs are the rocket engineer's numbers: **thrust, specific
impulse, and propulsive efficiency, each with a band**, for a geometry you drew; the fusion physics
inside is in service of those three numbers. Thrust from a plasma code is not itself new (Princeton's
UEDGE run of the direct-fusion-drive edge reported thrust and specific impulse as its primary result,
and electric-propulsion models now ship thrust with calibrated ~10 % bands for Hall thrusters). What
no existing tool does is run the whole chain in one domain, from the coils and injectors through a
self-consistent plasma-made field and a continuous neutral-gas-to-plasma medium to the exhaust plane,
with a per-cell validity and pedigree on every number and a rank-robustness statement between
geometries (§13). That is why the same solver that reproduced a chemical rocket's thrust blind is the
right foundation rather than a plasma code adapted after the fact.

**The first program's deliverable, stated so it cannot drift.** One year of solo work ends at the
magnetic-nozzle and gas-dynamic-mirror tiers (milestones 1–7): thrust and specific impulse for drawn
magnetic-nozzle thrusters validated against measured devices, and the first self-heated engine whose
core temperature emerges from control geometry. Fusion-core outputs in that program are **rank-only**:
the instrument orders control geometries over the declared band and refuses where it cannot. Absolute
fusion thrust to better than the band is not a deliverable of the first program, for a physics
reason (§16), and no milestone claims it.

---

## 1. What the instrument is for

A rocket engine turns stored energy into directed momentum through a working medium. In a chemical
engine the energy is released by reactions inside the gas. In a fusion engine it is released by
nuclear reactions inside the plasma. In a magnetic-nozzle plasma thruster it is deposited by fields
into the plasma. In every case the physics that decides performance is the same problem: a
compressible, conducting, reacting medium moving through a shaped domain under forces, exchanging
heat with walls, and leaving through a nozzle. The differences between engines are differences in
the medium's state, not in the laws.

CRUCIBLE takes that literally. There is one medium-state vector, one law per phenomenon evaluated
against that vector, and no branch in the code that asks which kind of engine this is. A cold
un-ionized gas with no magnetic field is one corner of the state space; a hot magnetized fusion
plasma is another; the laws are continuous between them, and a law that does not apply in a corner
evaluates to zero there by the state, never by a switch. Geometry, boundary inputs, and operating
schedules are pure data. That is what makes it a sandbox: you can change the shape of the machine,
not only its numbers, and the answer changes for physical reasons.

The workflow contract: **define geometry → choose boundary objects from a library → set an
operating profile → run → receive distributions and fields.** One configuration file, one command,
no mesh sessions, no per-concept code.

What exists today (seven weeks of work; about 47k lines of Rust, 4.4k of CUDA, 4.7k of Python; 314
tests) is the chemical corner: a 2-D axisymmetric (swirl-carrying, optionally 3-D) cylindrical
finite-volume solver with tabulated thermochemistry, viscous and thermal transport, conjugate wall
conduction, an implicit-explicit integrator, cut-cell geometry from contours or STL, a conservation
audit every step, bit-reproducible on CPU and GPU. It reproduces the RL10 rocket's thrust and
specific impulse blind against NASA test data. Everything below is a widening of that solver, not a
replacement.

**Three terms used throughout.** A *tier* is a selectable level of physics for one quantity: the
chemical injector, for example, has a *prior tier* (mixing efficiency supplied as a validated band
for a known injector family) and a *resolved tier* (the unmixed streams enter and mixing is
computed). A *closure* is a law for a quantity the grid does not resolve, carrying a declared band.
A *side solver* is a lower-dimensional problem posed on the computed fields (along a field line, or
an eigenproblem on a snapshot) whose result feeds back as a source or a verdict.

---

## 2. Scope

**Criterion.** An engine is in scope when the energy that produces its thrust is released inside the
working medium by a law whose rate depends on the medium's own local state (density, temperatures,
composition, field), and the release is continuous or repetitively pulsed. Chemical reaction rates,
fusion reactivities, Joule heating of a computed current, and resonant field deposition at a
computed field strength all qualify. Energy deposited at a rate independent of the medium (a laser
pulse, a reactor at fixed power) is admitted as a source term. Energy released by a chain reaction
whose rate depends on the state everywhere at once (fission criticality) or by a burn the firewall
forbids (pellet implosion) is not simulated; its products enter as boundary objects.

| Regime | Status | Computed | Boundary objects |
|---|---|---|---|
| Chemical (liquid bipropellant; expander and gas-generator cycles) | built, validated blind on RL10 | chamber and nozzle flow, combustion, wall heat, emergent chamber pressure | pumps, valves, tank heads, igniter pulse |
| Magnetic fusion (field-reversed, mirror, gas-dynamic trap, hybrids) | new | the confined plasma including its core, its edge layer, propellant mixing, magnetic nozzle, thrust; the core temperature emergent | coils (geometry + current schedule), heating deposition, fuel and propellant injectors, biasing electrodes |
| Magnetic-nozzle plasma thrusters (helicon, ECR, VASIMR-class) | new; same physics as the fusion exhaust; the validation anchors for it | expansion, detachment, thrust | RF or ECR deposition object |

**Deferred to Appendix A:** pulsed engines (magneto-inertial, laser inertial, antimatter-catalyzed).

**Out of scope, permanently:** solid propellants; nuclear-thermal engines and the neutronics pipeline
(shelved in the design record, not built); Hall and ion thruster discharge physics; fusion-electric
power conversion; combustion-instability certification; 3-D turbulence resolution; trajectories; the
rest of the spacecraft.

---

## 3. Principles

1. **One medium, one law per phenomenon.** Every law is written once and evaluated against the local
   medium-state vector; no `if(material)`, no `if(regime)`. A configuration that would need bespoke
   physics is a defect in the core. *Reason:* the code that reproduced the RL10 is the code that runs
   the fusion exhaust; that is what makes an unanchored result credible.

2. **Physical march only; reach-and-hold is the verdict.** Every result is a transient march from a
   declared initial state under commanded inputs. Stage 1 asks whether the engine reaches and holds
   its operating point. No steady-state shortcuts. *Change from the previous doctrine:* the
   requirement that the march begin from a cold chamber with a spark is withdrawn as a certifying
   requirement; the measured ◆C4 run showed a spark-to-steady verdict at desktop resolution cannot
   distinguish "the engine does not light" from "the instrument cannot resolve flame-holding."
   Spark startability stays in the code as a flagged, lower-pedigree capability.

3. **Every result is a distribution with a pedigree.** No bare numbers. Closures carry declared bands
   and are swept; the comparison of two configurations is a rank-robustness statement (§9.4), and the
   instrument refuses to rank where the band swamps the geometry signal. **The band on an output is
   the band on the closure multiplied by the loop gain of the operating point (§9.3)**; near a
   self-heated operating point that multiplier can be large, and it is reported, never absorbed.

4. **Closures are calibrated once, globally, and frozen; no device appears in both the calibration
   set and the validation set.** Fit once to the declared calibration set (§10.3), recorded with a
   band, never retuned per geometry or per case. Calibrated-once is an instrument; retuned-per-case is
   a fudge. This is the mechanical form of "same physics, no anchor fitting."

5. **The reduced law has three terms; the instrument keeps one, models one, and drops one, and says
   which.** By the Mori-Zwanzig theorem the exact evolution of a cell's coarse variables is a term in
   the current coarse state, a memory of its own history, and a noise term carrying the unresolved
   detail. The instrument computes the first, models the second with a stated first-order kernel
   (§7.4), and carries only the *statistics* of the third: its mean enters the transport coefficients,
   its variance is what the declared band and the dwell-averaged readout stand for. The noise's path
   is not simulated. Every run reports per cell how well the preconditions of this reduction hold
   (§8).

6. **Determinism: bit-exact per build, per device.** Gather-formulated kernels, fixed-topology
   reductions, no atomics in physics. Marker populations deposit through per-cell sorted gathers, not
   scatters, so they obey the same rule; a fixed seed makes marker noise reproducible, not small, and
   the noise level is a reported validity item. A run is one reproducible realization; the physical
   result is the statistic over many decorrelation times.

7. **Double precision throughout the certified spine.** Offline calibration solvers and external
   oracle codes run in whatever precision they use; they produce tables and cross-checks, not
   certified results.

8. **2-D axisymmetric with swirl is the default.** The magnetic state (poloidal flux plus azimuthal
   field, §4) is complete only in axisymmetry; the retained 3-D capability applies to the gas corner
   and is not extended to the magnetized state in this proposal. Non-axisymmetric stability is
   evaluated as a labelled diagnostic on the axisymmetric state (§6, row 27).

9. **Firewall, no amendment path.** No export-controlled or restricted codes or data. Pulsed-burn
   physics stops at published gains, yields, and partitions.

10. **Verification by parts, three rungs:** analytic or manufactured solution → published benchmark
    or open reference code → measured hardware. Every operator climbs all three before its output
    counts in a pedigree.

11. **Fail loud, never clamp.** Out of envelope means refuse or flag. A cell outside a closure's
    validity keeps its computed value, receives a lowered pedigree, and is flagged. An implicit solve
    that does not converge to its declared tolerance is a refusal, never a silently capped answer.

12. **Visualization is a deliverable, now.**

13. **Compute: one laptop and one desktop GPU** (RTX 4070 Ti SUPER, 16 GB). Minutes on the laptop at
    coarse resolution; hours to a few days on the GPU at fine resolution.

---

## 4. The medium-state vector

Every cell carries the following, averaged over the cell. Each component is present everywhere;
components that do not apply in a given corner of state space are inert there by continuity. The
azimuthal coordinate is θ throughout; φ is reserved for the electrostatic potential.

| Component | Symbol | Notes |
|---|---|---|
| Mass density | ρ | |
| Momentum (r, θ, z) | ρu | swirl carried at N_θ = 1 |
| Total energy | E | |
| Electron internal energy | e_e | identically zero where there are no free electrons; collisional equilibration collapses T_e → T_i as collisionality rises |
| Species vector | ρX_k | the single owner of composition: chemical reactants and products; neutral atoms and molecules; ions of each fuel and propellant; fusion products. The ionization fraction is derived from it, not stored separately |
| Reaction progress | c | the chemical burn-progress *statistic* (the fraction of a cell's mass that has crossed the flame), not a composition; it selects between the unburnt and burnt thermochemical branches and its rate law self-zeroes where there is nothing to burn |
| Poloidal magnetic flux | ψ | B_pol = ∇ψ × ∇θ/(2πr); divergence-free by construction; zero in the chemical corner |
| Azimuthal field | rB_θ | |
| Electrostatic potential | φ | a derived, constrained field: the curl-free part of the poloidal electric field given by the one Ohm's law (§6, row 14); identically zero where there are no carriers |
| Turbulence energy | k_t | fluctuation energy of unresolved eddies (§7) |
| Turbulence decorrelation rate | ω_t | |
| Fast-particle moments | (n_f, p_f, E_f) | minority population with orbits too large for the fluid; carried by counter-seeded δf markers (§6, row 10) |
| Validity diagnostics | ρ_i/L, s*, k_⊥ρ_i, ν*, Kn, d_i/L, … | derived, not evolved; drive flags and pedigree (§8). Defined in §8 |

Temperature is derived in every cell from (ρ, E, e_e, X_k) through the constitutive spine. There is
no assumed region. Where the field and the carriers are zero the magnetic and electrostatic
components are exactly zero and the equations reduce to the validated gas equations.

**Decision recorded: how the electric field and potential are handled.** *v0.2 proposed* an
independent potential solved from current continuity with the density kept inside the polarization
term. *The v0.2 review showed* that operator's coefficient diverges as the field vanishes and the
equation degenerates in the ideal limit, so "dormant by continuity" was asserted, not designed.
*v0.3 formulation:* the electric field is given everywhere by the one generalized Ohm's law (row 14);
current continuity is then automatic. The potential is the curl-free part of that field, obtained by
a plain Poisson solve (an operator that is well-conditioned everywhere and has nothing to diverge),
and it exists as a state component because three things need it: the E×B shearing rate that feeds
the turbulence closure, electrode boundary conditions (a biased plate imposes φ and the sheath's
current-voltage law, row 23, supplies the current), and consistency with the ambipolar potential the
kinetic tier computes along open lines. *For:* one Ohm's law in both the induction equation and the
potential; general; no blocker. *Against:* one extra Poisson solve per step.

---

## 5. Geometry, fields, and the three global couplings

- **Grid.** Cylindrical-structured (r, θ, z) about a declared axis; exact metrics; adaptive azimuthal
  resolution with N_θ = 1 as the default; static (r, z) refinement zones.
- **Geometry in.** Constructive solid geometry for sweepable designs; watertight STL for CAD
  components. Both voxelized to cut cells with per-cell volume fractions and face apertures;
  conservation exact through State Redistribution.
- **Control geometry versus magnetic geometry.** The user controls coils, walls, injectors, and
  electrodes: the control geometry, and the design variables of any search. The plasma's own
  currents make the magnetic geometry: in a high-pressure bottle (a field-reversed configuration) the
  plasma current is what reverses the field and the coils only set the boundary; in a low-pressure
  mirror the coils make the geometry and the plasma barely bends it; plasma pressure over magnetic
  pressure is the continuous dial between. The magnetic geometry (separatrix radius and length, null
  position, effective mirror ratio, connection lengths) is an **output**, reported and drawn.
- **Coils.** Geometry objects with current schedules I_j(t). Biot-Savart is linear in current, so
  each coil's unit field is computed once at configuration time (tested against the analytic
  on-axis loop field and the magpylib oracle) and the vacuum flux at any instant is Σ I_j(t) ψ_j.
  A changing coil current is a loop voltage: ∂ψ_vac/∂t enters the induction equation as a source
  (row 12). Field lines are the contours of ψ.
- **Walls.** Solid regions carry conjugate heat conduction. The gas–wall interface is one law
  blended by ionization (row 23).

**The three global couplings.** A cell's energy source depends on the cell's own state; the cell's
state is delivered by three couplings to everything else. Each is handled differently, and this is
where the difficulty and the innovation live.

1. **Transport of heat and particles.** Parabolic, neighbor to neighbor, slow (milliseconds). Closed
   by the turbulence-energy field with its calibrated band (§7). This is where the approximation lives.
2. **The electromagnetic field.** Elliptic, everywhere at once on the Alfvén time (microseconds):
   the poloidal flux and azimuthal field by induction, the potential by the Poisson solve. Its only
   closure is the resistivity, whose anomalous part inherits the turbulence band. The plasma-made
   geometry is a computed output.
3. **Nonlocal energy deposition.** Fast fusion products, neutrons, and radiation depositing far from
   where they were born. Carried as bounded tiers with declared prompt-loss, escape, and reflection
   fractions.

The timescales are ordered: the field responds in microseconds, turbulence decorrelates in tens of
microseconds, profiles evolve in milliseconds. Fast couplings settle inside slow ones, so the slow
fixed point is what the dwell finds. The same ordering is why the field and the stiff transport are
treated implicitly (§12).

---

## 6. Every law, as one law

R = resolved on the grid; C = closure with a declared band; S = side solver on the computed fields;
B = boundary object with citation, envelope, and band. "Chemical corner" describes what the law does
in a neutral reacting gas; "plasma corner" what it does in a magnetized plasma; the same code runs
both.

| # | Law | Chemical corner | Plasma corner | How | Rung-one test |
|---|---|---|---|---|---|
| 1 | Mass conservation | as built | same | R | Sod, MMS |
| 2 | Momentum with pressure, viscous stress, Lorentz force J×B, and interspecies friction (ion–neutral drag from charge exchange and elastic collisions) | J = 0, one fluid: Euler + viscous as built | J×B from the computed field; friction from tabulated rates | R | Brio-Wu; Bennett pinch held static |
| 3 | Total energy, with Poynting flux and Joule heating | as built | adds field energy terms | R | MHD energy audit |
| 4 | Electron internal energy with ion–electron equilibration | e_e ≡ 0 | evolves; equilibration ∝ n²/T^{3/2} | R | two-temperature relaxation |
| 5 | Species conservation with diffusion and reaction | mixture fraction + burn progress, as built | species vector | R | flame_1d grid independence |
| 6 | Ionization balance: ionization, recombination, charge exchange, as species reactions | rates vanish below ~3000 K | tabulated collisional-radiative rates | R + tables | Saha equilibrium recovery |
| 7 | Equation of state | tabulated (p, h, Z) surfaces, as built | two-temperature plasma with tabulated ionization energy; blended by ionization fraction | R + tables | table round-trips |
| 8 | Chemical reaction rate; flame speed and induction closures | as built, validated; the flame-speed wrinkling factor now fed by row 20 | fuel ions are not chemical reactants; rate ≡ 0 | R + C | flame_1d |
| 9 | Fusion reaction rate n_a n_b ⟨σv⟩(T_i) per fuel pair | no fusion fuel; 0 | Bosch-Hale tables (D-T, D-D, D-³He, p-¹¹B) | R + tables | reactivity vs published |
| 10 | Charged fusion-product birth, orbits, slowing-down to ions and electrons | 0 | counter-seeded δf markers deposited by per-cell sorted gather; Trubnikov-Stix slowing; prompt-loss fraction keyed to gyroradius over minor radius | S + C | slowing-down distribution |
| 11 | Neutron escape and wall deposition (D-T, D-D branches) | 0 | birth from row 9; straight-line escape with view factors and tabulated attenuation; wall heating and a fluence tally | S + C (escape and deposition bands) | point-source view-factor identity |
| 12 | Radiation loss: bremsstrahlung, line, synchrotron | negligible; declared magnitude band | analytic bremsstrahlung; tabulated line loss; synchrotron with wall-reflection band | R + C | loss-function tables |
| 13 | Induction: ∂ψ/∂t = −u·∇ψ + (η/μ₀)Δ*ψ + Hall and electron-pressure terms from the one Ohm's law (row 14) + ∂ψ_vac/∂t from the coil schedules; the azimuthal induction equation for rB_θ; ∇·B = 0 by construction | ψ ≡ 0 stays 0 exactly | evolved | R (resistive and Hall terms implicit) | Brio-Wu; Solov'ev equilibrium; resistive decay; whistler dispersion |
| 14 | The one generalized Ohm's law: E = −u×B + ηJ + (J×B − ∇p_e)/(e n_e), electron inertia neglected; the potential φ as the curl-free part of the poloidal E by a Poisson solve | no carriers: every term is zero; E ≡ 0, φ ≡ 0 | radial electric field, E×B shear, biasing, ambipolar potential | R (Poisson) | radial force balance with the swirl velocity; ambipolar potential analytic limit |
| 15 | Resistivity | dormant (no current) | Spitzer plus an anomalous term from k_t | C (spine) | Spitzer table |
| 16 | Viscous stress | Chapman-Enskog μ, as built | Braginskii ion viscosity, anisotropic in the computed field; blended by ionization | R + spine | Couette; parallel viscosity |
| 17 | Heat conduction, anisotropic; parallel flux blended continuously between the collisional (Braginskii) value and the kinetic value from row 22 by the harmonic form q = q_B /(1 + q_B/q_kin), applied as a face flux | isotropic Fourier, as built | k_∥ along field lines, k_⊥ across | R + S | Günter ring; the harmonic blend's two limits |
| 18 | Species diffusion, with a Knudsen validity flag for neutrals | Fick with ρD from the spine | ambipolar cross-field plus classical for ions; fluid neutrals flagged where their mean free path exceeds the cell (a deterministic neutral-transport tier is the recorded upgrade) | R + spine | binary diffusion |
| 19 | Turbulence-energy transport: ∂k_t/∂t + ∇·(uk_t − D_k∇k_t) = Γk_t − C_ε k_t ω_t − S_shear, with a companion equation for ω_t | Γ is the mean-flow strain-rate production of a k-ω model (§7.1); the chemical regime acquires a RANS turbulence model with a declared wall treatment | Γ adds the drift/interchange growth rate from row 21; S_shear from the E×B and diamagnetic shearing rates (row 14) | C (calibrated coefficients, §7) | gyro-Bohm scaling; decaying homogeneous turbulence; wall-heat-flux and mixing metrics (§10) |
| 20 | Turbulent transport coefficients: eddy viscosity ν_t = c k_t/ω_t; heat and species diffusivity ν_t/Pr_t, ν_t/Sc_t with Pr_t, Sc_t named constants with bands; floored at laminar or classical | a selectable tier alongside the injector prior and the wall law, so regressions are attributable | supersedes the "anomalous multiplier" | C | jet mixing station; L-mode critical gradient |
| 21 | Local linear rate Γ_lin in each cell: the growth rate of the fastest drift, interchange, lower-hybrid, or temperature-gradient mode, with finite-Larmor-radius terms, from nine geometric and state inputs (§7.2) | zero at B = 0 (no drift branches); the neutral production is the strain-rate term of row 19, not this row | computed | R (algebraic, per cell) | slab ITG; LHDI; FLR stabilization |
| 22 | Kinetic transport along open field lines, both species: a conservative low-rank Fokker-Planck solver in (s, v_∥, v_⊥) along each computed open line | no field lines; dormant | ion end loss; electron cooling and trapped populations; parallel heat flux (into row 17); ambipolar potential (consistent with row 14) | S | Pastukhov; gas-dynamic outflow; Ahedo-Merino kinetic nozzle solutions |
| 23 | Wall interface: heat transfer and current | turbulent wall function with declared y⁺ treatment | sheath: Bohm flux, heat transmission coefficient, and the sheath current-voltage law that closes electrode circuits | R + C | Colburn; Stangeby sheath |
| 24 | Wall material response: conjugate conduction; surface recombination of ions to neutrals | conduction as built | recombination recycles ions as neutral gas | R | station-4 cooled wall |
| 25 | Detachment diagnostic: the surface where the flow kinetic energy density exceeds the magnetic energy density (ρu² > B²/μ₀), and the electron magnetization on the expansion scale, both drawn and reported; **labelled a diagnostic, not a criterion**, because the 2025 review of propulsive magnetic nozzles finds the detachment mechanism (resistive, electron-inertia, induced-field) unresolved and the three candidates divergent | not applicable | drawn as surfaces; not a halt | R (derived) | Ahedo-Merino detachment analysis |
| 26 | Heating deposition: neutral beams (attenuation computed from the species and rate tables), RF and ECR (resonance located from the computed field; absorbed fraction a declared band) | igniter energy deposit as built | as stated | R + B | beam attenuation analytic |
| 27 | Non-axisymmetric stability diagnostic: ideal-MHD eigenproblem for the first two non-axisymmetric modes on the checkpointed state, plus the published kinetic criteria (S*/E for tilt, average-minimum-B for flutes) evaluated on the same state | trivially stable | reported with growth time; **labelled** "ideal, no flow, no FLR"; not a halt | S | a known ideal-MHD kink benchmark; the published FRC tilt boundary as a *comparison*, not a pass/fail |
| 28 | Structural margins; coil limits (field, current density, heat) | as built | plus coil objects | B/S | Roark |
| 29 | Boundary objects: injectors, pumps, tank heads, igniter pulse, heating objects, electrodes, radiator | as built | plus coils, heating, electrodes | B | envelope refusal |
| 30 | Conservation audit every step: mass, momentum, energy, species, magnetic energy and flux | as built | extended with the field | R | identity to round-off |
| 31 | Predictability and validity diagnostics | Knudsen and e-folding only | all of §8 | R (derived) | flagged, never clamped |

Three rows are where the unification is real rather than cosmetic, and each carries its own honesty
note.

**Rows 13 and 14 are one Ohm's law.** The Hall term is not small in the engines in scope: the ion
skin depth is about 3 cm at 10²⁰ m⁻³ and 30 cm at 10¹⁸ m⁻³, comparable to a cell and to the device.
Dropping it in the induction equation while keeping it in the potential would be the regime branch
principle 1 forbids, so it is kept in both and treated implicitly (it introduces a whistler wave whose
explicit time step scales with the cell size squared). Where the ion skin depth exceeds the gradient
length the fluid description itself is failing, and that ratio (d_i/L) is a validity flag.

**Row 19 changes the chemical regime, and the change is priced.** Today the chamber has no
turbulence model: mixing is a validated prior and wall heat is a wall law. One turbulence-energy law
evaluated in a neutral gas is a standard k-ω model, and published k-ω predictions of rocket-chamber
wall heat flux scatter from 1.4 % to 28.5 % against measurements depending on wall treatment, mesh,
and compressibility corrections. So the new tier's wall-heat band may be wider than the wall law's
±20–30 %, and it depends on the mesh. The old tiers stay selectable as controls. The RL10 thrust and
specific impulse are weakly sensitive to wall heat and mixing, so a re-earned overlap is weak
evidence; the gate for the unified law is the **wall-heat-flux profile and a mixing metric** against
anchors that are not in the calibration set (§10.4, milestone 6).

**Row 22 is the largest single fidelity gain, and it still assumes five things.** Four exhaust
closures (the nozzle electron-cooling exponent; the parallel heat-flux limiter; the mirror end-loss
coefficient; the ambipolar drop) are one piece of physics: what the velocity distribution does along
an open field line. The cooling exponent is well measured on the devices that exist (1.15 ± 0.03 and
1.17 ± 0.02 on two thrusters; 1.2 ± 0.1 is the usual fit; the wider spread quoted in earlier drafts
came from its variation along a single field line, not between devices), so the tier's value is not
fixing a wide band there; it is computing the cooling and the end loss for geometries and regimes
(fusion-heated layers, strong mirrors, biased expanders) where no measured exponent exists. One deterministic solver along each computed line, for both species, replaces all four with a
computation whose precedent is the Merino-Ahedo fluid-kinetic nozzle model (free, reflected, and
doubly trapped electron populations). It still assumes: the line geometry is quasi-static over the
solver's cadence (the lag is a reported validity item); lines are decoupled except through the
cross-field source terms taken from rows 17–20; the magnetic moment is conserved (flagged near a
null and for fast ions, which row 10 carries instead); a model collision operator; and a low-rank
truncation that must be mass- and energy-conserving (the conservative dynamical-low-rank schemes of
Einkemmer et al.). "Low-rank" means the distribution along a line is represented as a short sum of
products of a function of position and a function of velocity, typically five to twenty terms, which
is why it costs megabytes and seconds rather than a kinetic grid.

---

## 7. The turbulence closure

### 7.1 The law

The cross-field leak in every real magnetic bottle is dominated by turbulence at the ion gyroradius
scale, 10 to 100 times the classical value; in a rocket chamber the mixing and the wall heat are
dominated by hydrodynamic turbulence. Neither is resolved on this hardware. Both are carried by one
transport law (row 19) whose production rate Γ is a sum of two self-zeroing terms:

Γ = c_s |S|² / ω_t  +  C_p Γ_lin,

where |S| is the mean-flow strain rate (this is exactly the k-ω production P = ν_t|S|² written as a
rate on k_t, so the neutral corner is a real k-ω model, wall-bounded shear included) and Γ_lin is the
drift-wave growth rate of row 21, zero where there is no field. The dissipation is C_ε k_t ω_t and the
shear suppression uses the E×B and diamagnetic shearing rates. Below the critical gradient Γ_lin → 0
and, away from mean shear, the field decays to laminar or classical transport: critical gradients
emerge, they are not imposed. The transport coefficients of row 20 follow. The coefficients
(c_s, C_p, C_ε, D_k, c, Pr_t, Sc_t, and the wall treatment) are calibrated once (§10.3) and frozen.
Declared band after calibration: about ±50 % in tokamak-like and neutral-gas regimes, a factor of
three either way when extrapolated to field-reversed and mirror cores.

Why this is defensible in a magnetized plasma: the eddies are nearly two-dimensional (constrained to
the plane across the field), two-dimensional turbulence conserves both energy and enstrophy and
self-organizes into large-scale flows, and its decorrelation time (tens of microseconds) is a
hundred times shorter than the profile evolution time, which is the separation the reduction needs.

**The named precedent that failed, and what is different here.** The one published attempt at a
k-ε-type turbulence closure in an edge-plasma transport code (SOLEDGE2D-EIRENE) found its
coefficients "vary substantially depending on the machine, the type of experiment, and even the
location inside the device," which destroyed its predictive power for a new configuration. That is
exactly the portability §10.3 needs. The difference claimed here, and it must be tested rather than
assumed: that model's production term was a constant; here the production rate is computed per cell
from the local linear physics (§7.2), so the machine-to-machine and place-to-place variation that
SOLEDGE2D had to absorb into its coefficients is carried by the computed rate, and only the
saturation level is a constant. Milestone 6's calibration campaign tests this directly: if the fitted
coefficients drift across the box states by more than the declared band, the closure has failed in
the same way, the rank-robustness readout refuses, and the first program's fusion outputs remain
rank-only by construction. This is the single biggest scientific risk in the proposal.

### 7.2 How far geometry reaches into the closure

Geometry enters the closure through the local linear rate and a mixing-length assumption (the
transport coefficient is proportional to the growth rate over the wavenumber squared). That is the
quasilinear approximation every reduced transport model makes; what is computed is the linear rate,
and it takes nine inputs from the geometry and the state, each a term in a known linear theory with
its own test and a physical case where it is the whole story:

| Input | Physics | Where it is the whole story |
|---|---|---|
| Gradient scale lengths of n, T_i, T_e | drive for every branch | profile stiffness: extra heating steepens the gradient to critical and no further |
| Magnetic curvature projected on ∇p | bad curvature drives interchange; good curvature stabilizes | a simple mirror's central cell is flute-unstable; minimum-B or biased ends fix it |
| Magnetic shear (rotation of the field direction across flux surfaces) | localizes modes, reduces growth | reversed-shear transport barriers |
| Connection length along the line | sets the parallel wavenumber; ballooning (modes localized where curvature is bad) and sheath-connected thresholds | scrape-off-layer blobs (filaments of plasma that detach from the edge and travel outward) scale with the length to the target |
| Flux-tube expansion and mirror ratio | trapped fraction; trapped-electron modes; loss cone | the gas-dynamic trap's confinement and expander cooling |
| E×B shearing rate (row 14) | suppresses turbulence when it exceeds the growth rate | vortex confinement in the Budker mirror: transverse losses fall to near classical as bias rises |
| Finite-Larmor-radius ratio ρ_i/L | large orbits average out short modes | TAE C-2 core: ion-scale turbulence measured suppressed |
| Collisionality and beta | branch switching; electromagnetic stabilization | high-beta FRC core |
| Open versus closed field lines | finite connection length; sheath in the dispersion | intermittent blob transport at the edge |

The finite-Larmor-radius term is mandatory: without it the closure would predict tokamak-like
turbulence in a field-reversed core and be wrong by an order of magnitude. Magnetic shear, flux-tube
expansion, and connection length are ill-defined near a field null, which is the headline case for
the field-reversed core; those cells are flagged on those terms (§8), and the closure there rests on
the curvature, gradient, and FLR terms alone. Only the saturation amplitude, the spreading rate, and
the dissipation are calibrated.

### 7.3 Zonal flows

Turbulence can generate its own sheared flow, which then suppresses the turbulence, a predator-prey
pair whose oscillations trigger confinement transitions in tokamaks. *Decision (revised in v0.3):*
this is not built in the first program. The scalar form carries no radial scale, its damping
coefficients would have no mirror or field-reversed data to calibrate against, and its known payoff
is a tokamak phenomenon. The mean-flow shear that matters for the engines in scope (biased-electrode
vortex confinement in mirrors) is already carried by row 14. A state slot is reserved; the offline
calibration boxes record zonal-flow statistics so the extension can be calibrated later without a
new campaign.

### 7.4 Memory and neighbors: how far the cell-and-history model goes

The turbulence-energy equation is a first-order memory: a cell's turbulence responds to a change in
its drive over a decorrelation time and leaks into neighbors over a correlation length (a few
gyroradii to a few centimeters). That is an exponential kernel with diffusive spreading; no memory
integral is evaluated, and the noise term is carried only through its statistics (principle 5).

- **Enough** where memory decays faster than profiles change and transport is diffusive on average:
  tokamak cores, mirror central cells, most of a nozzle. Reduced models of this class reach 5 to 25 %
  on measured tokamak core profiles when the edge boundary is taken from experiment; that is the
  ceiling for any local-plus-spreading model, and it is a profile accuracy, not a thrust accuracy.
- **Needs a measured kernel** where transport arrives in bursts across many cells (avalanches,
  blobs). The time-averaged flux is still predictable, the instantaneous is not, and the effective
  kernel is not exponential. The calibration boxes measure the actual kernel; if it is not
  exponential, the ω_t equation is fit to the measured kernel rather than assumed.
- **Stops** near a field null, where orbits are the size of the region and no local dispersion
  exists, and for three-dimensional mode coupling, which row 27 reports as a labelled growth time.
  Those cells are flagged.

---

## 8. Closure honesty: validity, pedigree, and predictability

Every cell computes, from its own state, the quantities that decide whether each closure applies
there. Defined once here:

- **ρ_i/L**: ion gyroradius over the local gradient length. Above about 0.3 the fluid and local
  closures fail (large orbits).
- **s\***: the number of ion gyroradii across the plasma (s* = r_s/ρ_i); small s* is the kinetic
  field-reversed regime.
- **k_⊥ρ_i**: the dominant mode's wavelength relative to the gyroradius; the dispersion solve is
  valid below about 1.
- **ν\***: collisionality, the ratio of the collision rate to the transit rate along the line; it
  selects collisional versus collisionless branches and drives the row-17 blend.
- **d_i/L**: ion skin depth over gradient length; the Hall regime and the limit of the fluid
  description.
- **Kn**: Knudsen number for neutrals, mean free path over cell size.
- **Kn_flight**: the flight length (flow speed plus thermal speed, over the collision frequency)
  over the local gradient length: collisions per transit. This, not the thermal Knudsen number, is
  the condition for local equilibrium to be reachable in a fast flow; a fast, thin flow can have a
  small thermal Knudsen number and still transit without a single collision. (Added 2026-09-10
  after the solver's counterexample to the closure theorem, which exploited exactly that gap.)
- **Inlet non-equilibrium**: the relative distance of each injected velocity distribution from a
  Maxwellian with the same moments. A coarse inlet condition carries only moments; a beam-like
  inlet must enter as its own object (row 26), never as a moment boundary condition.
- **Loss-to-collision ratio**: radiative and escape losses per collision time. Where it is not
  small, the collisional transport coefficients describe a relaxation that is not occurring.
- **Cadence lag**: the time since the row-22 solution on a line was refreshed, over the local
  transport time.

Each of the nine inputs of §7.2, each side-solver tier, and each blend carries its own validity hull
in these quantities. A cell can be outside the hull of one term and inside the others; its pedigree is
the weakest of its terms, and the flags name which term failed. Cells outside a hull keep their
computed values. A run whose thrust-carrying cells are mostly flagged says so in the verdict.

Every run reports its own predictability horizon: the e-folding time of the fastest local
instability and the decorrelation time of the turbulence field, against the dwell over which readouts
are averaged. Readouts are statistics over many decorrelation times.

---

## 9. Execution model and outputs

### 9.1 Stage 1: reach and hold

From a declared initial state (a filled chamber, or a formed plasma at declared parameters) under
commanded boundary inputs and declared schedules, march until the readouts hold within tolerance for
a dwell, or until a halt.

### 9.2 Halts

Checked every step: melt or vaporization of a structural material; burst margin below one; coil
quench; radiative collapse; recombination of the exhaust; failure to reach the operating point within
the horizon; conservation-audit failure; an implicit solve that misses its tolerance. Evaluated per
checkpoint and reported, not halted on: the detachment surfaces (row 25) and the stability diagnostic
(row 27). First essential loss ends the run with mechanism, location, and time.

### 9.3 Outputs

Thrust; specific impulse; propulsive efficiency (jet power over power delivered to the exhaust);
fusion power; gain Q; wall, coil, and neutron loads; the emergent magnetic geometry (§5); **loop
gain** g = ∂(P_heat + P_fusion,retained − P_loss)/∂P_loss at the operating point, dimensionless:
g < 0 is self-stabilizing, g → 0⁻ is marginal, g > 0 means a perturbation grows (runaway or quench).
The output band is the closure band times |1/(1 − g)| and is reported as such. Also: the stability
diagnostic; the predictability horizon; the full field set per cell with validity and pedigree; every
scalar as a band from the closure sweep.

### 9.4 Rank-robustness

For any two configurations the instrument reports the fraction of the closure band over which the
ordering of each figure of merit holds, and declines to rank where that fraction is below a declared
threshold. This is the primary output of a geometry search. Geometry variables split into those
inside the feedback loop (field shape, mirror ratio, where fuel, heating, and bias go in), which need
the full loop and the band sweep, and those downstream of it (nozzle flare, expansion ratio), which
can be swept cheaply once a loop solution exists.

### 9.5 Visualization

Field lines (ψ contours), potential contours, density, ion and electron temperature, ionization,
turbulence energy, the detachment surface (row 25), the emergent separatrix, the validity map, and
the band plots, from the existing per-cell field output. Built in milestone 1.

---

## 10. Verification and calibration

### 10.1 Rung one: analytic and manufactured, per operator

Sod and the whole-operator manufactured solution (existing); Brio-Wu; the whistler dispersion
relation; the on-axis current-loop field; a Solov'ev equilibrium and a Bennett pinch held static for
100 Alfvén times; the Günter anisotropic-conduction ring, on cut cells, with the perpendicular
numerical leakage measured; the harmonic blend's two limits; the Ahedo-Merino isothermal and
polytropic nozzle solutions and their kinetic electron extension; the Pastukhov loss rate and the
gas-dynamic outflow limit; the ambipolar potential in the collisionless limit; slab
ion-temperature-gradient and lower-hybrid-drift growth rates; finite-Larmor-radius stabilization;
the gyro-Bohm scaling of χ; decaying homogeneous turbulence; the slowing-down distribution; the
point-source view-factor identity.

### 10.2 Rung two: open reference codes (oracles, never linked)

DIMAGNO (two-fluid nozzle); the Merino-Ahedo fluid-kinetic nozzle model; TokaMaker (axisymmetric
equilibria); PLUTO and Athena++ (MHD and Hall-MHD); GX (flux-tube gyrokinetics, tokamak-like limit
only, on the same GPU); QuaLiKiz (reduced transport, overlap region); Hermes-3 (drift-reduced
multi-species edge transport, the closest cousin); the 2016 UEDGE direct-fusion-drive edge case at
its own resolution; magpylib (coil fields); an OpenFOAM k-ω case of a rocket chamber.

### 10.3 Calibration set (used to fit the frozen coefficients; nothing here is also an anchor)

Embedded two-dimensional drift-fluid turbulence boxes run offline in batch over the instrument's own
computed states (a separate small solver whose only job is to measure the closure's saturation and
kernel; a thousand 64² boxes need about 160 MB of state and minutes to hours to saturate and
flux-average, not one minute as v0.2 said); GX in the regime where both apply; the published
tokamak L-mode critical-gradient database and QuaLiKiz in the overlap region; a canonical
turbulent-jet and backward-facing-step dataset for the neutral corner; the published sheath
heat-transmission range; the published fusion-product prompt-loss studies. The coefficients are
recorded with their band and frozen.

### 10.4 Rung three: experimental anchors (validation only; none used in §10.3)

Acceptance at every anchor is declared in advance **relative to the measurement's own uncertainty**:
magnetic-nozzle thrust-balance measurements carry about ±25 %, and published 2-D models of the same
devices miss by up to +59 % at high power. A pass is an overlap of the instrument's band with the
measured band; a claim of agreement tighter than the measurement is not made.

| Anchor | Validates | Status |
|---|---|---|
| RL10A-3-3A (NASA TM-107318): thrust, Isp, and the jacket heat pickup as the wall-heat-flux check | chemical corner end to end; the unified turbulence law in the neutral corner (milestone 6, gated on wall heat flux and a mixing metric with acceptance declared in advance) | earned; to re-earn under row 19 |
| VASIMR VX-200 (72 ± 9 % at 4900 ± 300 s) | magnetic nozzle, kinetic electron cooling | planned |
| Helicon and ECR magnetic-nozzle thrusters | nozzle, detachment, facility effects | planned |
| PFRC-2 scrape-off-layer measurements | edge layer, propellant heating, blow-off | planned |
| Budker gas-dynamic trap, incl. vortex confinement vs bias | mirror confinement, expander, the E×B shear term, the sheath current law | planned |
| TAE C-2 / C-2W: FLR suppression and global confinement | the FLR term; the turbulence band as extrapolated | planned |

---

## 11. Regime by regime: computed, assumed, coupled

| Engine | Assumed (band) | Computed | Backward coupling |
|---|---|---|---|
| Chemical | turbulence coefficients (global); wall treatment; or the injector prior at the prior tier | everything else | none |
| Field-reversed drive | turbulence coefficients; heating deposition; the fluid description itself near the null (flagged) | core temperature and fusion rate, plasma-made field and separatrix, edge layer, propellant heating, nozzle, thrust, gain, loop gain, stability diagnostic | neutral penetration into the core (row 18, flagged by Knudsen), edge temperature, layer blow-off: computed and checked |
| Mirror / gas-dynamic trap | turbulence coefficients (near-classical here); beam deposition | the whole plasma including throat loss (row 22) and expander; temperature emergent from control geometry | two-way; both sides computed |
| Magnetic-nozzle plasma thruster | source deposition | expansion, kinetic cooling, detachment, thrust | none |
| Sheared-flow Z-pinch | — | — | **refused:** the exhaust flow stabilizes a 3-D kink; 2-D cannot compute it |

The core temperature for magnetic engines is computed through rows 4, 17, 19–22. Its uncertainty is
the turbulence band times the loop gain, and because the coefficients are global, geometry
sensitivity is real: two control geometries are compared under the same coefficients across the
same band.

---

## 12. Tractability

Explicit marching is limited by the Alfvén speed. At 1 cm cells, CFL 0.4, deuterium:

| Density, field | Alfvén speed | Steps per ms of physics | Wall time per ms at ~5 ms/step | Same, at the 2–4× priced below (before the preconditioner) |
|---|---|---|---|---|
| 10²⁰ m⁻³, 0.3 T | 4×10⁵ m/s | 10⁵ | ~8 min | ~20–35 min |
| 10²⁰ m⁻³, 1 T | 1.4×10⁶ m/s | 3.5×10⁵ | ~30 min | ~1–2 h |
| 10¹⁸ m⁻³, 4 T (the dilute plume) | 5.5×10⁷ m/s | 1.4×10⁷ | ~20 h | ~40–80 h |

A 2-D device mesh is 6k cells at 1 cm and 600k at 1 mm. The dense core is minutes to an hour per
millisecond of physics; a quasi-steady approach is an overnight. The dilute plume is the explicit
bottleneck and the decision between a semi-implicit fast-wave treatment and a declared density floor
with its own convergence test is made in milestone 3, not deferred.

The additions are priced honestly. The per-cell dispersion solve is comparable to a hyperbolic
sweep, about 5 % of a step. The two elliptic tiers dominate: the anisotropic conduction (parallel
over perpendicular conductivity up to 10⁸) and the implicit Hall and resistive induction will take
thousands of conjugate-gradient iterations on the present Jacobi preconditioner, an estimated 2 to
4× the current step, until a field-aligned or geometric-multigrid preconditioner exists. That
preconditioner, with a measured iteration-count-versus-anisotropy test, is in milestone 3. An implicit
solve converges to its declared tolerance or the step is refused; there is no fixed iteration cap
standing in for convergence.

For comparison, measured in 2026: a fully kinetic 3-D field-reversed formation (79M cells, 190M
particles, 5 µs of plasma, electron mass inflated 15×, light speed cut 10×) took 3.7 days on four
V100 GPUs; a flux-tube gyrokinetic box (GX) takes under four minutes on one GPU but only in
tokamak-like geometry; a gyrokinetic eigenvalue solve now takes 0.01 to 0.1 s. These are oracles, not
the engine.

---

## 13. Similar projects, and what none of them does

Hermes-3 (UKAEA/LLNL) is the closest cousin: one code, configurable 1-D to 3-D, transport or
turbulence, arbitrary ion and neutral species, drift-reduced fluid equations on curvilinear grids. It
has no neutral-gas combustion corner, no coils-to-nozzle-to-thrust chain, and no declared-band
pedigree. SOLPS-ITER and UEDGE are multi-fluid edge transport codes with prescribed cross-field
transport. GBS and TOKAM3X are drift-reduced turbulence codes. QuaLiKiz and TGLF are reduced
transport models for tokamak cores. DIMAGNO and the Merino-Ahedo kinetic model cover magnetic nozzles
in isolation. The one-dimensional Miki-Diamond model carries turbulence intensity, flow shear, and
profiles as coupled fields, the precedent for §7. Engineering RANS codes carry the neutral-gas half.
None spans neutral gas to magnetized plasma continuously in one operator, none puts the plasma-made
field, the exhaust, and the thrust in the same domain, and none reports a per-cell validity and a
rank-robustness statement. That intersection is the instrument.

---

## 14. Relationship to the existing code and documents

### 14.1 Code: widen in place

Grid, tables, config, registry, units, and constants carry over unchanged. The conduction module
hosts the ψ diffusion and the Poisson solve. PPM reconstruction is width-agnostic. The
conserved-state width sits behind one type alias (283 references, 36 array positions recompile
untouched); the real edits are the full-slot writers, the flux assembly, the Riemann solver, and the
conservation audit's hardcoded seven-name table with two special-cased slots, about 25 to 35
functions plus a few CUDA defines. Two closed enums (the transport spine and the diffusion component
selector) become seams first. The GPU conjugate gradient's per-iteration host round-trip and the
duplicated deterministic reduction are fixed before any physics lands. Staging: widen with inert
slots (certificates byte-identical) → decide the Ohm's-law scope (done: row 14) and then swap the
Riemann solver to HLLD for the ideal part in one reviewed commit (chemical certificates move at the
10⁻¹²–10⁻¹⁰ level; gate 5 becomes a declared-tolerance diff) → induction with the field audit →
potential → two temperatures and ionization → turbulence field → kinetic tier and markers →
stability diagnostic. A live coarse RL10 march joins the gate before the state is touched.

### 14.2 Documents

VISION_SCOPE v1.6 and PLAN_CHEMICAL_SANDBOX are replaced by this proposal's adopted form. Of the 35
Layer-2 documents: the foundations (FND-1…7), the coupling set (COUP-2…8), SOLV-1 (widened), SOLV-4
(widened to fusion), SOLV-6/7, OFFL-3, OFFL-5, and VAL-1…3 stay and are amended. SOLV-3 (energetic
particles) is revised to rows 10 and 11. SOLV-2 (radiation operator) and SOLV-5 (pulsed events) go
dormant with Appendix A. OFFL-1/2 (neutronics), SOLV-8 (degradation clocks), OFFL-4 (annihilation
beyond published envelopes), and COUP-1 go dormant. New documents: the electromagnetic tier (rows
13–14); the turbulence closure and its calibration protocol (§7, §10.3); the kinetic-along-field-
lines tier (row 22); the stability diagnostic (row 27); validity, pedigree, loop gain, and
rank-robustness (§8, §9).

---

## 15. Build sequence

Each milestone ends with its tests green, the chemical certificate regenerating, and a picture.

1. **Consolidation.** Adopt this proposal; retire the startup-as-verdict and full-3-D-product
   rulings; lift the visualization gate; shelve the nuclear leg. Build the field viewer. Fix the two
   enum seams, the CG round-trip, the duplicated reduction; add the live RL10 march to the gate.
2. **Month-one de-risk, before solver code.** Compute the collisionality and Hall maps for the
   published direct-fusion-drive edge; reproduce the 2016 UEDGE point at its own resolution; run a
   sensitivity screen on a cheap surrogate over closure knobs versus control-geometry knobs,
   including loop gain. Where closure variance times loop gain swamps geometry variance for a
   candidate question, change the question now.
3. **3a. Field and flux.** Widen the state; the field audit; coil unit fields and the loop-voltage
   source; HLLD for the ideal part; Brio-Wu; Solov'ev and Bennett held static; the plume fast-wave
   decision (semi-implicit versus floor) made and tested. **3b. The electromagnetic tier.** Implicit
   Hall and resistive induction; the whistler test; the Poisson potential with its two analytic
   limits; the elliptic preconditioner with the measured iteration-count-versus-anisotropy test.
4. **Two temperatures, ionization, anisotropic transport, and the nozzle oracle.** Rate tables;
   Braginskii spine occupant; the Günter ring on cut cells with the leakage measured; the harmonic
   blend; DIMAGNO and Ahedo-Merino fluid solutions reproduced.
5. **The nozzle, fluid tier, anchored.** Thrust and divergence versus coil geometry against the
   helicon and ECR data: the first anchored result where control geometry changes the answer.
6. **The turbulence-energy field, both corners.** Dispersion solve with its growth-rate tests; the
   transport equation; the jet and backward-step calibration first; **then the RL10 re-earned under
   the unified law, gated on the jacket heat pickup and a mixing metric with acceptance declared in
   advance and the wall-law tier retained as a control**; calibration-box campaign; GX cross-check;
   the validity field and the rank-robustness readout.
7. **The kinetic tier along open field lines.** Pastukhov and gas-dynamic limits; the kinetic
   nozzle solutions; VX-200 and helicon anchors with computed electron cooling; the gas-dynamic trap
   as the first end-to-end fusion-heated engine whose temperature is emergent from control geometry,
   validated against the suppression-versus-bias curve it was not calibrated on.
8. **The edge layer and propellant mixing.** Neutral gas with its Knudsen flag, charge exchange,
   blow-off; PFRC-2 anchor; the first campaign: propulsive efficiency and specific impulse of a
   direct-fusion-drive-class exhaust as a surface over layer thickness, mirror ratio, injection
   location, and expansion ratio, under the declared band, with loop gain reported.
9. **The stability diagnostic.** Ideal-MHD kink benchmark; the published FRC tilt boundary as a
   comparison; the diagnostic wired into the readout.
10. **Ensembles and certificates.** Latin-hypercube sweeps over the declared bands; rank-robustness
    across a control-geometry family; the fusion-exhaust certificate.

---

## 16. What the instrument cannot do

- Predict the specific state of a turbulent region beyond its predictability horizon. It reports
  the horizon.
- Compute anomalous transport from first principles. It carries a calibrated turbulence field with a
  factor-of-three band outside the tokamak and neutral-gas regimes, amplified by loop gain.
- Watch a three-dimensional instability evolve. It reports an ideal-MHD growth time labelled as such.
- Rank confinement schemes by absolute fusion power when their difference is inside the band. It
  says so.
- Describe a plasma as a fluid where the ion orbit or the ion skin depth exceeds the gradient length.
  It flags those cells.
- Simulate pellet implosion, ignition, or burn; fission chain reactions; Hall or ion thruster
  discharges; solid motors.

---

## Appendix A. Pulsed engines: a feeder project

Magneto-inertial, laser-inertial, and antimatter-catalyzed engines release their energy in a
millimeter-scale burn over nanoseconds. The burn is radiation-hydrodynamics plus nuclear burn, the
class of code the firewall forbids, and for antimatter-catalyzed microfission it has not been
computed since the 1990s Penn State studies (one 2025 Geant4 fast-ignition study aside). What this
instrument can do for them is everything after the burn: a burnt-plasma initial-state object (yield,
partition into neutrons, photons, and debris, from published envelopes) expanding against the coils'
field in the same domain as any other engine, returning impulse per shot, field recovery, wall and
coil loads, and the rep-rate thermal state. That needs no new law beyond those in §6. It is scheduled
after the continuous engines because it does not exercise the self-heating coupling that is the
instrument's innovation, and it is kept open as the interface to a separate burn-physics project.

## Appendix B. Terms

- **Δ\*ψ**: the Grad-Shafranov operator, the axisymmetric form of the Laplacian acting on the flux
  function; Δ*ψ/(μ₀ r) is the toroidal current density.
- **Boussinesq**: two unrelated usages. In plasma vorticity equations it means replacing the density
  in the polarization term by a constant (rejected here, and now moot: §4). In turbulence modelling
  it is the eddy-viscosity hypothesis of row 20 (adopted here, as in every k-ω model).
- **Gyro-Bohm scaling**: the expected size of turbulent diffusivity, ρ_i² v_th/L, the reference
  against which the closure is checked.
- **Critical gradient**: the gradient below which drift-wave growth vanishes; above it transport
  rises steeply, which is why heated profiles are "stiff."
- **S\*/E**: the field-reversed tilt-stability parameter, the number of ion gyroradii across the
  plasma divided by its elongation; experiments find stability improving above about 3 to 3.5, a
  kinetic result the ideal-MHD probe of row 27 does not reproduce and is not asked to.

## Appendix C. Change log

| Version | Change | Why |
|---|---|---|
| v0.1 (2026-09-08) | first no-context draft | — |
| v0.2 (2026-09-09) | complete law table; potential as a state component; kinetic tier along open lines; nine geometric inputs; control vs magnetic geometry; three couplings; loop gain; per-term validity; similar projects | see the v0.2 log in the previous revision |
| v0.3 (2026-09-09) | **Electromagnetic tier rebuilt:** one generalized Ohm's law with Hall and electron-pressure terms used by both the induction equation and the potential; the potential is the curl-free part of that field by a Poisson solve (no divergent operator, no degenerate limit); coil schedules drive a loop voltage; μ₀ restored; Hall treated implicitly with the whistler CFL priced; d_i/L added as a validity flag | the v0.2 review found two Ohm's laws (a regime branch) and a potential operator that diverges as B → 0 |
| v0.3 | **Calibration and validation sets separated:** no device in both; RL10 and the Budker bias curve moved to validation only; the neutral corner calibrates on jet and backward-step data | anchor fitting by construction, contradicting principle 4 |
| v0.3 | **Turbulence production redefined** as strain-rate production plus drift-wave growth (a real k-ω model in the neutral corner, wall-bounded shear included); turbulent Prandtl and Schmidt numbers and the wall treatment named as constants; old wall-law and injector-prior tiers kept selectable as controls; milestone 6 re-gated on the wall-heat-flux profile and a mixing metric, not thrust | the v0.2 form gave zero production in wall-bounded shear and an unattributable gate |
| v0.3 | Row 17/22 parallel flux blended continuously (harmonic form) and applied as a face flux | a collisionality switch and an audit break |
| v0.3 | Row 22's five remaining assumptions listed; conservative low-rank scheme required; cadence lag made a validity item | overclaimed as written |
| v0.3 | Laws added: interspecies friction (row 2), neutron escape and deposition (row 11), fluid neutrals with a Knudsen flag (row 18), sheath current-voltage law (row 23), detachment criterion (row 25), heating deposition (row 26) | promised outputs or anchors with no law behind them |
| v0.3 | Zonal flows removed from the first program (slot reserved); stability probe demoted from halt to labelled diagnostic | no mirror/FRC calibration data; the tilt boundary is kinetic, not ideal-MHD |
| v0.3 | Species vector made the single owner of composition; burn progress clarified as a statistic; θ/φ symbol collision fixed; halt cadence reconciled; timing arithmetic corrected; "under 30 %" replaced by the priced 2–4× until the preconditioner exists; box campaign cost corrected | consistency and honesty items from the review |
| v0.3 | Honest headline added; output band = closure band × loop gain stated as a principle; Mori-Zwanzig framing restated as keep/model/drop; "temperature in every cell" restated as "core temperature emergent"; 3-D retained for the gas corner only; marker determinism via sorted gathers | overstatements named by the review |
| v0.3 | Appendix B (terms) added | readability for a non-specialist decision-maker |
| v0.3b (2026-09-09) | After the independent yes/no assessment: uniqueness claim corrected (UEDGE already reports DFD thrust/Isp; EP models ship calibrated bands), differentiators restated; the first program's deliverable pinned to milestones 1–7 with fusion outputs rank-only; electron-cooling spread corrected to the measured 1.15–1.2 ± 0.1 and row 22's payoff restated; detachment demoted to a diagnostic (mechanism unresolved per the 2025 review); §12 wall times shown at the priced 2–4×; §10.4 acceptance declared relative to measurement uncertainty (±25 % thrust balance, +59 % model misses); the SOLEDGE2D k-ε portability failure named as the biggest scientific risk with the test that decides it | five factual corrections and the verdict's condition |
