# VISION & SCOPE — CRUCIBLE
### A Generalized Physics Sandbox for In-Vacuum Propulsion, from Chemical to Antimatter

*(**CRUCIBLE** = **C**ross-**R**egime **U**nified **C**omparator **I**ntegrating **B**ounded-**L**ifetime **E**ngines. Working name; the backronym is descriptive, not binding.)*

**Document status:** Source of truth for all design decisions, 12-month inquiry project, July 2026 – June 2027.
**Change control:** This document outranks any design doc, code comment, or conversation. If a design decision contradicts this document, either the decision is wrong or this document gets amended *first*, with a dated entry in the Amendment Log (§15). Amendments to §4 (Scope Doctrine), §9 (Validation), and §10 (Firewall) require explicit justification of why the original boundary failed, not just preference.
**Owner:** Ben. **Code authorship:** Claude writes all code (~3 h/day budget). **Hardware:** one high-end desktop (assume 16–32 cores, 64–128 GB RAM, **and a high-end GPU used for the runtime 3-D field solver** — v1.3; see §8). *Offline* nuclear/table generation remains CPU-bound.

---

## 1. Vision

One simulator, one consistent physics framework, every propulsion regime that operates in vacuum — chemical, nuclear thermal, nuclear pulsed, fission fragment, fusion, antimatter — evaluated with the *same* thermal, radiation, structural, and uncertainty machinery, so that results across regimes are actually comparable and every output carries a defensible confidence range.

No such tool has ever existed. The 2026 landscape survey (three independent research sweeps, results embedded throughout this doc) found three disconnected strata: per-regime physics codes that don't compose (CEA, Hall2De/WarpX, Griffin+BISON, Geant4), cycle-level system tools that top out at ~2.5 regimes at 0-D fidelity (NPSS), and mission tools that reduce every engine to an Isp number (GMAT/Copernicus/STK). Fusion and antimatter propulsion have never been in *any* reusable engine-design tool. Meanwhile DRACO was cancelled (2025) and FY2026 zeroed the NTP/NEP budget lines: the field's center of gravity has shifted from hardware to modeling and concept down-selection — exactly the niche this tool occupies.

**The inquiry question:** *When every propulsion regime is evaluated under one consistent, uncertainty-quantified physics framework, which fast-Mars propulsion concepts survive — and what operating points exist in the gaps between regimes that no one has proposed?*

**The backbone of the project is the inquiry question itself: full-regime integration and search inside one consistent system.** That capability — not any single pre-chosen paper topic — is the quest. The research campaigns (§3) are a portfolio of candidate expressions of it, selected opportunistically as the instrument matures; the strongest paper may well be a topic not yet on the list. The workflow contract (§4.3) is what keeps 90% of effort on reaction physics and out of engineering rabbit holes. The fidelity doctrine in one line: **utter completeness for the reaction, zero simulation for what we don't care about.**

---

## 2. Success Criteria (measurable, end-of-project)

| # | Criterion | Threshold |
|---|-----------|-----------|
| S1 | Chemical regime validated | The predicted **p-box** for vacuum Isp and c* **overlaps the RL10A-3-3A reference band** (NASA TM-107318: 73.4 kN, Isp 440–446 s, Pc 475 psia, MR 5.0, c\* 7824 in/s), reported with the **Ferson area metric**, both **blind** (open-mode feed) and **calibrated** (closed-mode expander) — VAL-2 §3.1. **≤ 2% is the reported target, not a binary gate:** a wider honest, overlapping band is a success if it is the tightest the stated assumptions allow (cf. S8). c\* and C_F are validated separately (v_e = c\*·C_F) |
| S2 | Nuclear thermal regime validated | Chamber temperature, flow, and Isp within **≤ 10%** of NERVA test state points (Pewee, NRX-A6, XE-Prime; NASA-CR-184270); k-eff within **≤ 300 pcm** of the KRUSTY ICSBEP benchmark (HEU-MET-FAST-101) |
| S3 | Every headline result is a distribution | No bare-number outputs; all campaign results ship with propagated 68/95% intervals and a pedigree score (§7.8) |
| S4 | **The backbone works:** all §5 regimes integrated and cross-comparable in one system | The cross-regime sweep (§3.0) runs end-to-end over the mechanism library by M12 |
| S4b | ≥ 2 research campaigns completed to paper-draft quality | Chosen opportunistically from the §3 portfolio, or from topics that emerge during the backbone search |
| S5 | The two-stage execution model works end-to-end | Any configured engine returns WORKS / DOESN'T WORK (+ reason), and if WORKS, an operational-lifetime estimate with limiting mechanism (§7.6) |
| S6 | Reproducibility | Any result regenerable from config + pinned tables + seed; tables regenerable from their pipeline scripts. Determinism relaxed to a **negligible tolerance** (deviations orders of magnitude below the physics error bound, provably non-spiraling); regimes with genuine chaotic/turbulent sensitivity use **fixed-order deterministic reductions** (accepting the cost) so instantaneous state cannot butterfly-diverge and results never conflict |
| S7 | Class deliverables | A demoable milestone every ~2 months; thin web viewer by month 12 if schedule allows (§12) |
| S8 | **Assumption-bounded prediction in untested regimes** | For every regime — including those with no experimental anchor — a run states its **explicit, minimal, research-grounded assumptions** and reports a **single honest error bound** under them. Target ≤10% vs reality where validated constituent physics + feasible resolution allow; a wider bound (up to order-of-magnitude) is a **success if it is the tightest the narrowest reasonable assumptions yield** — provided it is stated, not hidden. The instrument exists to **resolve** untested-regime uncertainty as far as the physics permits (hierarchical validation-by-parts of each constituent law), not to accept it |

Explicit non-goals: publishing the code as a polished open-source product, real-time performance, and photorealistic visualization. Nice if they happen; never traded against S1–S8. **In scope as of v1.3: full-resolution 3-D field simulation at runtime** (adaptive-dimension — §5/§7); reduced-dimensional solving is an *optimization applied only where geometry and state are symmetric*, never a fidelity ceiling.

---

## 3. The Backbone and the Campaign Portfolio

### 3.0 The Backbone — Full-Regime Integration & Search *(the quest itself; not optional, not a campaign)*
The non-negotiable deliverable is the **capability**: every regime in §5 running inside one consistent, uncertainty-quantified system, cross-comparable on the same axes — culminating in the systematic sweep of (energy source × energy-to-propellant coupling × propellant) across the 2,000–20,000 s Isp / 10–100 kN fast-Mars gap. No such sweep has ever been done; every published hybrid concept (wave-rotor cycles, LANTR, bimodal NTP/NEP) was found by human intuition, not search. Pedigree scoring (§7.8) partitions the resulting map into "validated-physics finds" and "extrapolated leads." This is where the anticipated not-yet-known core insight is most likely to surface; the architecture, roadmap, and descope ladder all serve it first.

Below the backbone sits a **portfolio of candidate campaigns — possibilities, not commitments.** Each is publishable regardless of which way the answer falls (a deliberate selection criterion), each is listed with the module capabilities it needs, and none is the backbone. Campaign slots in the roadmap (§12) are filled opportunistically as the instrument matures; topics that emerge from the backbone search itself take precedence over anything pre-listed here.

### C1 — Nuclear Salt-Water Rocket: settling the order-of-magnitude dispute *(candidate — default first slot)*
**Question:** Does the NSWR deliver Zubrin's claimed ~10,000 s-class Isp at high thrust (Zubrin 1991), or does the water→steam transition throttle the neutron flux and crater performance, as a 2024 1-D neutronics/thermal-hydraulics study found? And under what feed geometry does the propellant *not* pre-detonate upstream?
**Why it's the default first pick:** The least-simulated concept in the entire field — its central feasibility question has never had a published coupled 3-D-informed treatment. The physics needed (transient neutronics + two-phase flow) is the smallest superset of the validated NTP module. An order-of-magnitude dispute means any rigorous answer is a paper.
**Method:** Transient point/few-group kinetics with reactivity coefficients from offline OpenMC sweeps over (salt species, enrichment, void fraction, geometry); quasi-1D two-phase reacting flow through chamber and nozzle; feed-system criticality margin maps from dedicated OpenMC geometry sweeps; full nuclear-data + correlation UQ.
**Output:** Isp/thrust distributions vs design parameters; a feasibility map (stable-burn region vs predetonation region); comparison against both prior claims.

### C2 — Centrifugal Nuclear Thermal Rocket: does vapor entrainment kill the Isp? *(candidate)*
**Question:** The CNTR (rotating liquid-uranium fuel, ~1800 s Isp claimed) lives or dies on uranium vapor/droplet loss into the hydrogen exhaust — the acknowledged biggest open question in advanced NTP (active UAH/OSU literature through 2025: Acta Astronautica 2025, AIAA 2025-0377). What is the actual entrainment rate across the operating envelope, and what Isp survives it?
**Method:** Reduced-order rotating annular-film model (film stability, bubble transport, evaporation at the liquid-gas interface) coupled to the validated NTP flow/kinetics machinery; entrainment closures from published correlations with explicit error bands; bound the answer with limiting-case analyses where the reduced model's validity is thin.
**Why it matters:** Results land in a live, funded research conversation; the propellant-mass-loss economics of CNTR flip on this number.

### C3 — ICAN-II Antimatter-Catalyzed Microfission, Re-Optimized *(candidate — weakest anchors)*
**Question:** The 1990s Penn State design needs ~140 ng of antiprotons for a Mars-class mission — orders of magnitude beyond world production. Modern optimization over target geometry, antiproton delivery, and pellet composition, within *published* energy-partition envelopes: how far can the antiproton requirement actually be pushed down?
**Method:** 1-D Lagrangian radiation-hydrodynamics of pellet response; antiproton annihilation/fission source terms from digitized PS177 data + one-time Geant4 runs; magnetic nozzle from the plasma module; parameter-space optimization under UQ. **Hard constraint:** stays inside the sensitive-physics firewall (§10) — published pellet gains and energy partitions only, no implosion-design physics beyond the open literature.
**Why:** Untouched since the 1990s; unique territory. Weakest validation anchors of the portfolio — hence candidate-only, and hence the widest (honestly reported) uncertainty bands.

*(The former "C4 hybrid map" is not a campaign — it is the backbone, §3.0. Additional candidates — fission-fragment thrust envelopes, magnetic-nozzle efficiency limits, hybrid duty-cycle concepts surfaced by the sweep — may be added to this portfolio via the Amendment Log without ceremony; the portfolio is a menu, not a contract.)*

---

## 4. Scope Doctrine

### 4.1 The Two Razors (test every feature request against these)
1. **The Reaction Razor:** *Simulate a phenomenon only if energy release or thrust generation flows through it. Everything else is a boundary object with literature specs and an error band.*
2. **The Failure-Mode Razor:** *If a subsystem's failure mode is mechanical or electrical engineering, it is a constraint check, not a simulation.*

Worked consequences (precedents, not exhaustive):
- Turbopumps: boundary objects (ΔP, ṁ, η from literature). CNTR's rotating fuel: **in scope** — the rotation *is* the reaction chamber physics. Wave rotors: **excluded from v1** — the moving hardware *is* the concept; better eight concepts honest than nine with one faked.
- Circuits/power electronics/power conversion: spec-sheet boundary objects (a pulsed-power driver is a current waveform with a rise time and jitter). **But magnetostatics is core physics**: coil geometry → Biot-Savart → static B-field is computed, because magnetic nozzles, confinement, and fragment collimation live on field geometry.
- Bearing life, seal leakage, insulation breakdown: pass/fail limit checks against literature values.

### 4.2 Hard Boundaries
| Boundary | Ruling |
|---|---|
| Environment | Vacuum only. No atmosphere, no launch, no reentry. |
| Thermal | Full thermal model inside the engine; heat crossing the radiator port exits the simulation. The radiator interface is a boundary object: a heat sink with stated capacity and coolant return temperature. |
| Ship | Not simulated. Fuel, propellant, and electrical power enter through **defined ports** with stated conditions. Thrust, torque, heat, and dose exit as reported outputs. |
| Mechanics | General Newtonian bookkeeping: net thrust vector, torque, momentum flux. No trajectory propagation — export thrust/Isp/mass curves for GMAT et al. |
| Radiation | Fully modeled, deterministic — **not** runtime Monte Carlo. *(v1.3)* It is **one** emit-once/transport-once transport operator over an energy-partitioned photon field (M1 two-moment in the thermal/optically-thick band; deterministic multigroup Sₙ in the streaming/nuclear band — a numerical split by **energy**, not a seam by **origin**; SOLV-2, radiation-partition invariant in COUP-2). Neutrons share the same Sₙ operator with different interaction data. Nuclear vs thermal are *source spectra* into one operator, not separate solvers. |
| Failure | Binary WORKS / DOESN'T WORK plus lifetime estimate. No cascading-failure simulation: the moment something essential to engine function is lost, the sim halts and reports (§7.6). |
| Degradation | Ablation/erosion/corrosion/burnup run as a **separate second-stage simulation** after the function sim passes (§7.6). |

### 4.3 The Workflow Contract
Everything in the architecture serves this loop, and any feature that doesn't serve it is out of scope:

> **Define geometry → pick mechanisms from the library → set an operating profile → run → receive distributions.**

A user (Ben, or the batch sweeper in C4) writes one config file, runs one command, and gets: verdict, performance distributions, field snapshots, lifetime estimate, pedigree score. No mesh generation sessions, no solver babysitting, no per-concept custom code.

**The sandbox doctrine (companion to the fidelity doctrine).** This is a *sandbox*: dozens of geometries and novel reactor/engine configurations are **pure data** over **one seamless physics core**. Two commitments make that real and are non-negotiable:
1. **No physics seams.** Each physical law is written once and evaluated against a **local medium-state vector** (density, temperature, ionization, composition, B-field, flow). "Cold solid," "warm dense matter," and "hot magnetized plasma" are corners of that vector, blended by continuous functions — never selected by an `if(material)`/`if(regime)` branch. A fast particle slowing in plasma and in a wall, a neutron from fusion and from fission, obey the *same* operator with different local state / source spectrum. (This is the operative form of §5's classification and is enforced at design time.)
2. **No per-concept code.** Defining a new geometry or reactor configuration requires no new physics code and introduces no new seam; a configuration that would need bespoke physics is a defect to be fixed in the core, not the config. Reaction modules (fusion, fission, annihilation) are *sources* that emit particles into the shared transport; they never carry their own private transport.

---

## 5. Regime & Physics Classification

The three-way classification: **RUNTIME** (computed every step, in Rust), **OFFLINE** (precomputed once into HDF5 tables by the Python pipeline, interpolated at runtime), **IGNORED** (replaced by literature closures, boundary objects, or nothing).

> **v1.3 doctrine override (adaptive-dimension unified grid).** The runtime default is a **single 3-D multiphysics field simulation on the world-state grid** — the grid itself *evolves* flow / plasma / fields / reactions over the medium-state vector `M` every tick (Rules 12/13 taken to their end: one representation, one set of conservation laws, reactions as sources). **Dimensional reduction is an adaptive optimization, not the primary model:** the solver collapses a region to 2-D-axi or 1-D **only where geometry *and* state are symmetric** (a **conservative, controlled-error** reduction via **adaptive azimuthal resolution** — conservative coarsening/refinement of the ring grid, v1.4, superseding the v1.3 spectral-truncation formulation — conserved integrals preserved exactly, discarded azimuthal content bounded by the collapse threshold; *not* lossless), spending 3-D cost only where symmetry breaks. The per-row "quasi-1D / 2-D-axi" fidelity labels below are **superseded** by this and are rewritten in the Layer-2 cascade (FND-2 + the SOLV/OFFL re-scope). Kinetic (PIC) fidelity is **not** a runtime default — it is offline calibration or an occasional one-time high-compute run (§8).

### 5.1 RUNTIME — Rust, every timestep
| Physics | Fidelity | Build/Source | Reference |
|---|---|---|---|
| Compressible reacting flow (nozzles, channels, chambers) | Quasi-1D + 2D-axisymmetric, tabulated equilibrium chemistry, HLLC Godunov | **BUILD** | Toro, *Riemann Solvers*; Anderson's MacCormack nozzle chapter as didactic baseline |
| Two-phase flow extension (NSWR, boiling/void) | Quasi-1D drift-flux / homogeneous-relaxation | **BUILD** | Standard two-phase closures w/ error bands |
| Heat conduction in solids | 3-D on the voxel world-state | **BUILD** | Explicit/implicit FV, textbook |
| Thermal radiation exchange | The **one radiation operator's thermal band** (M1 two-moment, SOLV-2); optically-thin cavity/surface exchange enters the *same* operator as a surface-emission/radiosity term with config-time Monte Carlo ray-cast view factors | **BUILD** (~200-line view-factor MC + a page of linear algebra) | pyviewfactor as cross-check oracle |
| Reactor kinetics | Point / few-group kinetics ODEs; reactivity feedback ρ(T_fuel, T_mod, density, drum angle) interpolated from offline tables | **BUILD** | Standard PKE; coefficients from OpenMC (§8) |
| Nuclear heating & dose | Deposition/dose kernels from tables applied on voxel grid — for neutrons & photons of **any origin** (fission, fusion, activation, annihilation): one transport operator, only the source spectrum differs | **BUILD** (application layer) | Kernels precomputed offline |
| Plasma energy balance + magnetic nozzle | **One uniform** 0-D/quasi-1D two-fluid flux-tube energy equation; source/sink menu self-zeroing outside its regime (RF/ohmic in; fusion self-heat from the transport module; ionization/excitation, wall/sheath, charge-exchange, bremsstrahlung, synchrotron, line, conduction; PdV→directed KE = thrust); electron closure = piecewise-polytropic γ_e tied to the detachment plane, **calibrated to published kinetic results** | **BUILD** | Ahedo/Merino & Little-Choueiri magnetic-nozzle literature; γ_e and detachment fraction are declared uncertainties |
| Pulsed events (pellets, plasma slugs) | 1-D Lagrangian radiation-hydrodynamics: von Neumann–Richtmyer staggered grid, tabulated EOS/opacity, radiation via the **same M1 closure as the one radiation operator** (grey → multigroup; a private flux-limited-diffusion model is a Rule-12 closure seam — SOLV-5/SOLV-2); runs as a sub-simulation "event" returning impulse + energy partition | **BUILD-WITH-REFERENCE** | Bowers & Wilson; Castor, *Radiation Hydrodynamics*; MULTI-IFE (CPC/Mendeley) and SNEC sources as validation oracles — reference codes, never dependencies |
| Fusion reactions & product sourcing | Bosch-Hale ⟨σv⟩ (D-T, D-D, D-³He); Sikora-Weller 2016 / Putvinski 2019 for p-¹¹B. Reactions **emit birth-spectrum source terms** (species, energies, branching) into the shared transport (charged) and neutron operator (neutrons) — they carry no private transport | **BUILD** (pages of coefficients, implemented directly) | Bosch & Hale, *Nucl. Fusion* 32 (1992) |
| Magnetostatics | Biot-Savart from declared coil geometry, computed at config time, static during run | **BUILD** (an afternoon) | magpylib used offline as cross-validation oracle |
| Structural margins | Analytic quasi-static checks: thin-shell hoop stress, thermal stress, pressure limits → pass/fail + margin | **BUILD** | Roark's formulas; voxel-hex linear-elastic FEM is a *stretch* upgrade, not v1 |
| Ablation / degradation clocks (stage 2 only) | Energy-balance surface recession with B′ tables; H₂ corrosion correlations; burnup; fluence/DPA limits; decay heat via ANS standard correlation | **BUILD** | B′ tables generated offline from equilibrium chemistry |
| Energetic-particle transport & interaction (ALL fast particles: fission fragments, fusion α/p, annihilation pions/muons, sputtered ions, fast e⁻) | **One** friction-&-diffusion law over a local medium-state vector (cold↔plasma continuum, no seam): slowing-down = stopping = self-heat source; coupled guiding-center↔Boris orbits in sampled B; self-zeroing radiative losses (brems ∝1/m², synchrotron ∝1/m⁴); Z_eff(v) charge state; optional in-flight secondary reactions | **BUILD** | Unified per §4.3 "no seams" doctrine; stopping from offline cold→plasma data; charge-state band declared |
| Electric-propulsion performance | 0-D validated performance tables (power → thrust/Isp/η) from published throttle data; discharge physics **not** simulated — EP is mature, measured technology and needs no rediscovery for the backbone map | **SOURCE (data)** | NSTAR throttle table (Goebel & Katz), Hall/VASIMR published envelopes |
| Steady magnetic-confinement fusion source (DFD-class) | Parameterized energy/particle source from published confinement scalings; confinement physics **not** simulated (it is the fusion field's own open question) | **BUILD (parameterization)** | Published FRC/mirror scalings, with declared bands |
| Newtonian outputs | Thrust/torque/momentum-flux integration over exit planes | **BUILD** | trivial |
| UQ engine | Latin-hypercube ensembles over table covariances + correlation bands + BC tolerances; Sobol indices; rayon-parallel across ensemble members | **BUILD** | Saltelli; SALib as offline cross-check |

### 5.2 OFFLINE — Python pipeline, run once per table version
| Physics | Tool | Product |
|---|---|---|
| Neutron/photon transport: k-eff, deposition kernels, dose kernels, kinetics parameters, reactivity coefficient sweeps | **OpenMC** (MIT license — embed/redistribute freely; **MCNP is export-controlled: never used, never requested**) | HDF5 kernel & coefficient tables per geometry class |
| Nuclear data | **ENDF/B-VIII.1, TENDL-2023, FENDL-3.2c** via OpenMC data API / NJOY / ENDFtk | Pointwise + multigroup libraries |
| Nuclear-data uncertainty | **SANDY** perturbed-library sampling (~300 samples × 3–6 dominant nuclides: U-235, H, C, Be/BeO, Mo, Zr) + OpenMC re-runs | Covariance-propagated coefficient distributions |
| Equilibrium chemistry & rocket thermo | **Cantera + NASA CEA/RocketCEA** (CEA now open-source on GitHub, full ~2000-species Glenn database) | Property tables T(p, MR, species), c*, expansion tables; B′ ablation tables |
| Antiproton annihilation | **Geant4** one-time runs across **multiple physics lists** (FTFP/INCL++/CHIPS — the inter-list spread *is* the uncertainty) + **digitized PS177 data** (antiproton-induced fission probabilities on U-238/Bi/Pb/Au) + AEgIS-era prong-multiplicity data (arXiv:2407.06721) | Product spectra, energy-partition, fission-probability tables. **Conversion efficiency and waste-heat load are computed downstream** by transporting these products (neutrinos escape; γ / charged pions·muons deposit or escape per geometry — SOLV-2/3), never assumed. **Model-form uncertainty = the measured per-quantity physics-list spread** (tight where PS177/AEgIS pins it, e.g. prong multiplicity; a factor of several where lists diverge, e.g. U-238 fission probability), carried per-quantity — not a single blanket % |
| Plasma EOS & opacity | Analytic Saha/QEOS **default**; SESAME (request to LANL) and TOPS opacities (web request) as upgrades if granted | EOS/opacity tables |
| Magnetic-nozzle closures | Optional: a few expensive WarpX kinetic runs to calibrate the electron-cooling exponent | Closure coefficient + band |
| Particle–medium interaction data (cold→plasma) | Cold: NIST ASTAR/PSTAR + SRIM-class. Plasma: Li-Petrasso / BPS / RPA-LDA dielectric. Combined into **one continuous** S(E, species; medium-state) with no join discontinuity | Unified stopping/friction + Z_eff + scattering for all fast species (fragments, fusion products, pions) |
| **3-D verification oracles** | **OpenFOAM** (3-D compressible flow) and **Athena++** (3-D MHD) — established open codes, run offline on selected cases only | Spot-checks of reduced models where 3-D/non-axisymmetric effects are suspected (CNTR film stability, swirl, plume asymmetry, pulsed-plasma/nozzle interaction). Discrepancy → widened model-form band. Never coupled at runtime. |

### 5.3 IGNORED — with the replacement stated
| Excluded | Replaced by |
|---|---|
| Kinetic plasma (PIC) as the *runtime default* | Continuum (two-fluid / MHD) with validated transport closures; kinetic PIC reserved for **offline calibration / occasional one-time high-compute runs** (§8), not every-tick runtime |
| ~~3-D flow/plasma at runtime~~ — **removed v1.3: now the runtime default** (§5.1 banner, §7.1) | — (offline OpenFOAM/Athena++ runs remain a cross-check + calibration source, no longer the only 3-D) |
| 3-D turbulent combustion CFD | **Tiered mixing integration (R2, 2026-08-13):** **prior tier** — equilibrium chemistry + a cited injector c\*-efficiency **prior**, legitimate *only inside* the validity envelope of the injector family that measured it (mechanically enforced by the standard envelope refusal); **resolved tier** — the actual unmixed streams enter at the injection plane and mixing is *computed* on the grid (same operator, resolution/mode-ceiling lift), η_c\* an **output** with per-quantity model-form on universal closures. Architected now (COUP-7/SOLV-1/OFFL-3); resolved tier built when a blind prediction needs it |
| EP discharge physics (Hall/ion/RF) | 0-D validated performance tables from published throttle data (§5.1) |
| Fusion confinement physics (FRC stability, RF heating) | Published confinement-scaling parameterizations with declared bands (§5.1) |
| Runtime Monte Carlo transport | Precomputed deterministic kernels |
| Rotordynamics, circuits, power-conversion internals | Boundary objects (§7.5) |
| Fatigue, fracture, vibration, cascading failures | Static margins + halt-on-essential-loss (§7.6) |
| Trajectories, GNC, crew, shielding optimization, cost | Exported thrust/Isp/mass curves; downstream tools' problem |
| Wave-rotor concepts | Out of v1 entirely (Reaction Razor casualty) |
| Detailed fission-burn/implosion design physics | Published pellet gain & energy-partition envelopes only (§10) |

---

## 6. Two-Language Rule

**The runtime simulator is 100% Rust.** The offline table-generation pipeline is Python **only because it must be**: OpenMC's and Cantera's first-class APIs are Python-only and non-negotiable for nuclear/chemistry table generation. The seam between the two is exactly one artifact type: **versioned HDF5 table files** (Rust reads via the `hdf5` crate). Python never runs at simulation time; Rust never generates tables. No third language, no exceptions.

Core crates (defaults, swappable without doc amendment): `ndarray`, `rayon`, `serde` + TOML for configs, `hdf5`, `rand`/`rand_distr`. Viewer (late, if schedule allows): static HTML/JS reports or a small `wgpu`/`egui` tool — viewer choices never constrain the core.

Reference codes consulted but **never linked**: MULTI-IFE, SNEC (rad-hydro oracles); magpylib, pyviewfactor, SALib (cross-check oracles run in Python during validation only).

---

## 7. Architecture

### 7.1 Overview
*(v1.3 — adaptive-dimension unified grid; see the §5.1 banner.)* The runtime is a **single 3-D multiphysics field simulation on the world-state grid**: the grid **evolves** matter and fields — flow, plasma, conduction, radiation, transport, and reactions-as-sources — over `M` every tick, and **audits its own conservation**. Layers: **(1)** the voxel world-state (now also the field solver); **(2)** an *adaptive dimensional-reduction* pass that collapses regions to 2-D-axi / 1-D where geometry and state are symmetric (a compute optimization, not a separate physics); **(3)** the boundary-object library (reaction-decoupled engineering only); **(4)** the orchestrator (two-stage execution + UQ ensemble loop). The earlier **"grid is the accountant; reduced solvers are the physicists" split is superseded** — the accountant *is* the physicist; carrying two discretizations reconciled by gather/scatter was itself a seam that Rules 12/13 forbid. The reduced-dimension native-mesh + conservative-transfer machinery (former COUP-1) is absorbed into the adaptive-reduction pass and re-scoped in the Layer-2 cascade.

### 7.2 World-State Grid
- **Structure** *(v1.4)*: hand-rolled sparse block grid (shallow-wide, OpenVDB-style indexing; ~500–1500 lines) on a **cylindrical-structured index space** (i_r, i_θ, i_z about a **config-declared symmetry axis**; cell volumes ∝ r; reflecting axis treatment). **Azimuthal resolution N_θ is adaptive per region** — conservative ring coarsening/refinement; **N_θ=1 *is* the axisymmetric-with-swirl calculation**; regions eligible for adaptive collapse keep a guard resolution (N_θ ≥ guard) so growing asymmetry remains detectable, and full N_θ=1 is a recorded config assertion, never an adaptive decision. Refinement in (r,z) near material interfaces and steep gradients; data layout GPU-portable (SoA) from day one. Static (r,z) topology per run (geometry doesn't move; recession updates fractions, not topology).
- **Per-cell state:** material ID + volume fractions, temperature, density, nuclear heating rate, dose rates (n/γ), stress-margin flags, ablation recession depth, degradation-clock states; **reserved sharp-interface fields** (geometric liquid-surface capture — dormant capability, designed into the state layout now, built with the liquid-interface wave; v1.4).
- **Geometry in:** two paths, both voxelized at config time. **(a)** Config-file CSG (primitives + booleans + revolved profiles) for parametric/sweepable geometry — the backbone sweep (§3.0) requires this path. **(b) Watertight surface-mesh import (STL, plus OBJ/3MF)** per component, voxelized by ray-parity testing (~a few hundred lines, BUILD) — design components in any CAD tool and export STL; STEP→STL conversion in FreeCAD or similar is the supported route for parametric CAD files. An in-simulator CAD/B-rep kernel is permanently out of scope (§13): voxelization fidelity is controlled by grid refinement, not by carrying exact surfaces.

### 7.3 Regime Solvers & Native Meshes
Solvers declare **flow paths** (1-D station sequences), **axisymmetric zones** (2-D r-z patches, carrying azimuthal momentum — "2.5-D swirl" — where rotation is reaction physics, e.g. CNTR), **flux tubes** (magnetic nozzle), or **pellet lines** (1-D Lagrangian), each bound to world-state regions in config. Nearly every engine of interest is axisymmetric — this is exploited ruthlessly **as an adaptive, conservative, controlled-error reduction** (adaptive azimuthal resolution, v1.4 — conservative ring coarsening; collapse toward N_θ=1 only where geometry *and* state are symmetric — conserved integrals exact, discarded angular content bounded, *not* lossless). *(v1.3: the former "full 3-D solvers permanently out of the runtime" ruling is **superseded** — 3-D is the default, §5.1/§7.1. These declared flow-path / axi-zone / flux-tube / pellet-line forms survive as the **reduction targets** of the adaptive pass, not as the primary model.)* Offline 3-D oracle runs (§5.2) remain a cross-check and calibration source, no longer the *only* 3-D.

### 7.4 Couplers (the exhaustive list — a new coupler type requires a doc amendment)
1. **Deposition:** table kernels × current reactor power → volumetric heat sources on voxels (and dose accumulation).
2. **Wall exchange:** flow/plasma solver computes h(T_wall); voxel conduction returns T_wall — Robin BC, operator-split Picard iteration each step.
3. **Recession (stage 2 only):** ablation rates → voxel material fractions → flow-path area updates.
4. **Field sampling:** precomputed B-field interrogated by plasma/fragment solvers.
5. **Kinetics feedback:** voxel/flow temperatures + densities → table-interpolated reactivity → PKE → power → coupler 1.
6. **Event injection:** pulsed rad-hydro sub-simulation consumes pellet spec + driver boundary object, returns impulse, energy partition, and wall loading to the world-state.
7. **Port accounting:** boundary objects inject/extract mass, power, heat at declared ports; global conservation audited every step (violation > tolerance = bug, halt).

### 7.5 Boundary-Object Library (v1 set)
Pumps/pressurization (ΔP, ṁ, η), power supplies & pulsed drivers (waveform + jitter), lasers/particle beams (energy, pulse shape, spot), coil sets (geometry in → field computed), reactivity control schedules (drum angle vs t), antiproton delivery (rate, stopping position spread), injectors (tiered, §5.3 — envelope-bounded η_c\* prior or resolved unmixed-stream conditions), radiator interface (capacity, return temperature). Every boundary object carries: a literature citation, a validity envelope, and an error band. **An object without all three cannot be registered.**

### 7.6 Two-Stage Execution & the Failure Model
- **Stage 1 — FUNCTION:** Does the engine reach and hold its commanded operating profile? Transient march (or pseudo-transient to steady state; pulsed concepts run event sequences). **Halt conditions** (checked every step): melt/vaporization of a structural material, burst margin < 1, loss of criticality control (PKE excursion beyond commanded envelope), quench/extinction, choking/starvation at a port, conservation-audit failure. First essential-function loss → **halt immediately** → verdict `DOESN'T WORK (mechanism, location, time)`. Otherwise `WORKS` + performance report (thrust, Isp, power flows, field snapshots). **No cascading-failure modeling, ever** — one loss ends the run.
- **Stage 2 — LIFETIME:** Runs only after Stage 1 passes; separate simulation as specified. Freeze the operating point; march the slow degradation clocks (ablation recession, hydrogen corrosion of carbides, fuel burnup, neutron-fluence/DPA material limits, decay-heat accumulation where relevant) with large timesteps; re-check Stage-1 margins as geometry/materials degrade. First essential margin crossed → stop → **operational lifetime estimate + limiting mechanism + confidence interval**. If nothing crosses within a configured horizon, report `lifetime > horizon`.

### 7.7 Config & Results
One TOML config = geometry + materials + mechanism selections + operating profile + UQ settings + table-version pins. One results bundle = verdict, distributions (HDF5 + summary JSON), field snapshots, generated plots, pedigree report, and the exact config + seeds for reproduction (S6).

### 7.8 UQ Engine & Pedigree (first-class, not a post-processor)
- Every table carries uncertainty (nuclear covariances via SANDY-propagated coefficient distributions; Geant4/PS177 model-form bands; correlation error bands; boundary-object tolerances).
- Every campaign result is an **ensemble** (Latin hypercube over declared uncertainties, rayon-parallel members): outputs are distributions, never numbers. Sobol indices identify what drives the spread.
- **Validity enforcement:** interrogating any table/correlation outside its declared envelope flags the run (or refuses, per config). This is the structural defense against the tool becoming a speculation generator.
- **Pedigree** on every result (v1.3): a **PCMM-style predictive-maturity vector** (physics-model fidelity, code verification, solution verification, validation, UQ, input pedigree, results robustness — the seven axes of META-1 §4.1), each graded by rigor + independence, **reported by its weakest-link minimum — never averaged into a single "fraction of validated path" scalar** (COUP-6; META-1 §4.1). Reported alongside every headline result; the backbone map (§3.0) is partitioned by pedigree. *(Supersedes the earlier per-module-status × usage-weight scalar.)*

---

## 8. Offline Pipeline & Compute Budget

**Compute budget (v1.3): ~5000 GPU-desktop-hours across the project, to be *spent*, not conserved.** 3-D continuum multiphysics as the default is feasible in this envelope with adaptive resolution + GPU acceleration; the budget is sized to run the backbone sweep on the 3-D solver **at adaptive (mostly-reduced) dimensionality**, with full-3-D reserved for anchors/validation and UQ-at-full-3-D done **multi-fidelity** (META-1 §9; COUP-5). Occasional one-time high-fidelity / kinetic runs are reserved for an **HPC-grant allocation (e.g. NSF ACCESS) if secured** — the architecture must keep that door open without depending on it. Measured/published throughput on desktop-class hardware (from the 2026 compute survey):
- **OpenMC k-eigenvalue**, compact HEU fast-spectrum core (KRUSTY/NERVA-class): ~5–15 min/run to <10 pcm statistics on ~32 cores; ~1 min coarse. Fixed-source heating tallies to ~1%: ~10 min–2 h/run.
- **Table sweeps:** 100 points ≈ a day; 1,000 points ≈ ~10 days; 1,000 points with heavy heating tallies ≈ ~6 weeks. **Tier the fidelity:** coarse statistics across the full grid, high statistics at anchor points. Tables are regenerable dozens of times across the year — this budget is comfortable.
- **SANDY UQ:** ~300 samples/nuclide ≈ 1–2 days per dominant nuclide per model. Dominant nuclides only; all-nuclide simultaneous TMC is out of scope.
- **Geant4 antiproton tables:** ~1e6 events per configuration in minutes-to-hours. Effectively free; accuracy (not cost) is the declared limitation.
- **GPU: the runtime field solver's primary compute (v1.3).** The unified 3-D solver is GPU-accelerated (the ~5000-hour budget is GPU-desktop-hours) — this supersedes the pre-pivot "GPU unavailable for physics" ruling. *Offline* OpenMC stays CPU (its GPU port is branch-only/data-center-targeted), so nuclear table-gen is CPU-bound. GPU determinism is handled per META-1 §2 and S6: declared relaxed-reduction GPU paths reproduce within a negligible non-spiraling tolerance in **non-chaotic** regimes; **chaotic/turbulent** regimes mandate fixed-order deterministic reductions (bit-exact per build, accepting the cost); cross-platform comparison is statistical (ECT).
- Coupled steady-state OpenMC↔thermal Picard table points cost ~10× a bare run — feasible for ~100-point tables where needed.

Storage/versioning: every table file carries provenance metadata (generator script hash, library versions, date) and a semantic version; configs pin table versions (S6).

---

## 9. Validation Ladder

Every runtime module climbs three rungs before its outputs count as validated in pedigree scoring: **(i)** analytic/manufactured solutions, **(ii)** published benchmarks or reference-code cross-checks, **(iii)** real hardware data. No module skips rungs; modules stuck at rung (ii) (plasma, antimatter) are *permanently flagged* in pedigree and the papers say so.

**Blind-mode rule (v1.4, Ben 2026-08-14):** a blind validation run consumes **no quantity measured on the engine under test** — its inputs are exactly what a designer of a never-tested engine would have: the design specification (geometry, commanded profile), universal physics and closures, and technology-class data measured on *other* hardware, each with a cited band ("we assume we know exactly what we know for a fusion engine"). Anchor-derived values (a fitted η_c\*, a measured pump map) bind only in **calibrated** runs, labeled as such. Comparison against the anchor's measured values happens only *after* the blind run.

| Anchor | Data | Source | Target |
|---|---|---|---|
| RL10A-3-3A | 73.4 kN, Isp ≈ 444 s, Pc 475 psia, MR 5.0 + component maps | NASA TM-107318 (NTRS 19970010379); 1966 design report | ≤ 2% Isp, c* |
| NERVA Pewee | >500 MW, fuel-exit 2556 K, chamber 1833 K @ 18.6 kg/s, 4275 kPa; peak ideal vac Isp 901 s | NASA-CR-184270 (NTRS 19920005899); LA-4217-MS | ≤ 10% state points |
| NRX-A6 | 1120 MW × 60 min, Pc 4089 kPa, 32.7 kg/s | ibid.; WANL-TNR-223/224 | ≤ 10% |
| XE-Prime | 1140 MW, 244.75 kN thrust (thrust-stand), Isp ~710 s, 2272 K, 35.8 kg/s + operating map | ibid.; Aerojet RN-S-510 | ≤ 10%, map shape |
| KRUSTY | Criticality benchmark + warm criticals + full-power coupled run | ICSBEP HEU-MET-FAST-101 / KRUSTY-SPACE-EXP-001; MOOSE VTB KRUSTY model (open, incl. geometry/materials) | ≤ 300 pcm k-eff; coupled-transient shape |
| NSTAR | 16-level throttle table: 0.58–2.57 kW, 19–92 mN, Isp 1950–3120 s | Goebel & Katz Ch. 9 (JPL Descanso, free) | Within published η bands |
| VASIMR VX-200 | 5.8 ± 0.4 N, Isp 4900 ± 300 s @ 200 kW (plume-derived — treat error bars seriously) | Longmier et al. | Within stated error bars |
| Rad-hydro | Sod/Noh/Sedov analytics; MULTI-IFE & SNEC cross-runs | open sources | Standard convergence |
| Antimatter | Prong multiplicities & annihilation spectra | PS177; arXiv:2407.06721 | Computed source physics brackets PS177/AEgIS within the Geant4 physics-list spread (per-quantity, not a blanket %) |
| Coupling machinery | BEAVRS (not space-like; shakes down coupling only) | MIT open benchmark | qualitative |

The paper-credibility argument this ladder buys: *the same physics stack that reproduces RL10 and NERVA within stated bounds, applied to NSWR/CNTR/ICAN-II with propagated uncertainties, yields X.*

---

## 10. Sensitive-Physics Firewall (non-negotiable, no amendment path)

1. **No export-controlled or restricted codes or data**: no MCNP, no RSICC codes, no restricted evaluations. OpenMC + public libraries cover every need (verified, §5.2/§8).
2. **Pulsed-fission modeling stops at published envelopes**: pellet gains, energy partitions, and yields come from the open literature (ICAN-II/AIMStar papers, PS177, open ICF scaling laws). The simulator never performs implosion-to-criticality or fission-burn design optimization beyond parameter ranges already published. C3 optimizes *antiproton dose, target geometry at parameter level, and delivery* — not device physics.
3. **FLUKA licensing** (non-profit/institutional terms) is checked before any use; if unlicensable for an independent student, Geant4-only, with the wider uncertainty band that implies.
4. Anything that feels like it's drifting toward weapons-adjacent design detail gets cut, even at cost to a campaign. C3 is a candidate, never a commitment, precisely so this cut is always affordable.

---

## 11. Toolchain Summary

| Role | Tool | License / status |
|---|---|---|
| Runtime core | Rust (`ndarray`, `rayon`, `serde`+TOML, `hdf5`, `rand`) | all permissive |
| Offline transport | OpenMC | MIT |
| Offline nuclear data | ENDF/B-VIII.1, TENDL-2023, FENDL-3.2c + NJOY/ENDFtk; SANDY for covariance sampling | public / open |
| Offline chemistry | Cantera, NASA CEA (github.com/nasa/cea) + RocketCEA | BSD / open / GPLv3 (pipeline-only) |
| Offline annihilation | Geant4 (+ FLUKA only if licensable) | Geant4 license (permissive) |
| EOS/opacity | Saha/QEOS analytic default; SESAME + TOPS on request | request-based upgrades |
| Fusion reactivities | Bosch-Hale; Sikora-Weller/Putvinski coefficients | published |
| Oracles (never dependencies) | MULTI-IFE, SNEC, magpylib, pyviewfactor, SALib; OpenFOAM + Athena++ as offline 3-D verification oracles | reference/cross-check only |
| Validation data | NTRS/OSTI reports, ICSBEP/IRPhEP papers, JPL Descanso, MOOSE VTB | public (§9) |
| Mission handoff | GMAT-compatible thrust/Isp/mass curve export | out-of-scope consumer |

---

## 12. Twelve-Month Roadmap

Cadence rule for the 3 h/day Claude budget: **every session ends with tests green and work committed**; solvers are built test-first against analytic solutions; no session leaves a solver half-refactored. A demoable milestone every ~2 months (class deliverable, S7).

| Months | Phase | Exit criterion (demo) |
|---|---|---|
| 1–2 | Skeleton: config schema, world-state grid, conduction + view-factor radiation, quasi-1D HLLC flow, CEA/Cantera table pipeline | **M2 demo:** chemical engine validated vs RL10 (S1 met) |
| 3–4 | Nuclear leg: OpenMC pipeline, deposition kernels, PKE + reactivity tables, NTP module; two-stage execution + halt conditions + degradation clocks | **M4 demo:** NERVA/KRUSTY validation (S2 met); an engine that WORKS with a lifetime estimate |
| 5–7 | UQ engine + SANDY; two-phase flow; **campaign slot 1** (default candidate: C1 NSWR, incl. feed predetonation maps) | **M7 demo:** first campaign results with confidence bands → paper 1 draft data |
| 7–9 | Swirl/rotating-film capability; charged-particle transport; **campaign slot 2** (default candidate: C2 CNTR) | **M9 demo:** second campaign envelope → paper 2 draft data |
| 9–11 | Pulsed/plasma leg: 1-D Lagrangian rad-hydro, Geant4/PS177 tables, magnetic nozzle, fusion reactivities + confinement parameterizations, EP tables — **all §5 regimes integrated**; optional **campaign slot 3** (candidate: C3 ICAN-II) | **M11 demo:** every regime runs in one system |
| 11–12 | **The backbone sweep (§3.0)** — batch ensembles over the mechanism library; thin web viewer if schedule allows; writing | **M12:** the cross-regime map (S4) + ≥ 2 paper drafts (S4b) |

Slack policy: campaign slots absorb all schedule slips — fill them with shallower studies, different candidates, or leave slot 3 empty. The backbone integration thread (the module build order above) and validation/UQ never absorb slips: cutting those deletes the product. Campaign choice within each slot is free and decided at slot time, not now.

---

## 13. Risk Register & Descope Ladder

| Risk | Likelihood | Mitigation / descope |
|---|---|---|
| CNTR reduced-model validity too thin to defend | Medium | Bound with limiting-case analyses; widen declared bands; worst case C2 becomes "envelope + sensitivity" paper (still novel — nobody has published even that) |
| 1-D rad-hydro numerics fight back | Medium | 60 years of literature; MULTI-IFE/SNEC oracles; grey diffusion before multigroup; slack absorbed by C3 |
| Geant4 model-form uncertainty dominates C3 | High (per-quantity physics-list spread: ~few % where PS177/AEgIS pins a quantity, ×several where lists diverge — R3) | Per-quantity spread by design; C3's claim is *relative* optimization (dose reduction factor), which is robust to correlated model error |
| SESAME/TOPS requests denied/slow | Medium | Analytic Saha/QEOS default is the plan of record; tables are upgrades |
| Two-phase NSWR closures disputed | Medium | Multiple closure sets run as model-form ensemble members — the dispute becomes UQ, and the paper reports it |
| 3-D multiphysics cost & complexity *(v1.3 makes 3-D the intended default — the former "avoid 3-D" trap is now a budgeted capability)* | High | Adaptive dimensional reduction (§5.1/§7.1) spends 3-D cost only where symmetry breaks; ~5000 h budget + optional ACCESS grant (§8); the descope ladder coarsens 3-D resolution, it never deletes 3-D |
| Claude-session fragmentation (context loss across days) | High | This doc + per-module DESIGN.md files as session-start context; tests-green-per-session rule |
| Table/config version drift breaks reproducibility | Medium | Provenance metadata + pinned versions (S6) enforced by the loader, not by discipline |

**Descope ladder** (apply top-down under schedule pressure): drop viewer → leave campaign slot 3 empty → reduce campaigns to envelope/sensitivity studies → coarsen the backbone sweep (fewer axes, coarser grids — it shrinks, it never disappears) → reduce 2D-axi zones to quasi-1D everywhere. **Never descope:** the backbone integration S4, validation anchors S1/S2, UQ-on-every-result S3, two-stage execution S5, the firewall §10.

**Permanently out of scope (no amendment without a new project):** runtime Monte Carlo *transport* as a primary method *(clarified v1.4: the ban targets statistical transport — noise, irreproducibility, cost — as a primary path; SOLV-3's bounded, counter-seeded, fully-reproducible δf-marker fallback is permitted)*, cascading-failure modeling, trajectories/GNC, atmospheres, the rest of the spacecraft, multi-user/product engineering. *(v1.3 removed from this list — now in scope: **3-D flow/plasma solving in the runtime** as the adaptive-dimension default (§5.1/§7.1); an **in-simulator CAD/B-rep geometry kernel** as a permitted future capability if a workflow needs it. Kinetic **PIC** / high-fidelity kinetic transport is not a runtime default but is permitted as offline or one-time-grant calibration — no longer banned outright.)*

---

## 14. Prior-Work Integration Note
Ben's earlier design for a general whole-ship spaceflight simulator bridging regimes is superseded by this narrowed scope; the follow-on design agent decides what survives. Standing rule for that agent: anything imported must pass the two razors (§4.1) and the classification table (§5) as-is — the earlier design's whole-ship ambitions do not reopen §4.2 boundaries.

## 15. Amendment Log
| Date | Section | Change | Why the original boundary failed |
|---|---|---|---|
| 2026-07-14 | — | v1.0 baseline | — |
| 2026-07-14 | §1–3, 5, 7, 11–13 | v1.1: backbone reframed — full-regime integration & search is the quest, campaigns demoted to candidate portfolio (former C4 promoted to §3.0); STL/mesh geometry import added alongside CSG; offline 3-D verification oracles (OpenFOAM/Athena++) sanctioned; 2.5-D swirl in axisymmetric zones; regime audit added charged-particle/fragment transport, EP performance tables, fusion-confinement parameterizations | Owner directive: papers are possibilities not requirements; geometry authoring needed a CAD path; regime coverage audited for the backbone map |
| 2026-07-14 | §4.3, 5.1, 5.2 | v1.2: **sandbox / no-seams doctrine** made explicit (dozens of geometries & reactor configs as pure data over one seamless physics core; each law evaluated against a local medium-state vector, no `if(regime)` branch). Energetic-particle transport **unified** into one medium-agnostic friction law spanning cold matter↔plasma (fusion products, fragments, annihilation products, fast electrons all share it — slowing-down = stopping = self-heat). Reaction modules reframed as particle **sources** feeding shared transport (one neutron operator, one charged-particle operator). Stopping data extended to the cold→plasma continuum | Owner directive: the tool must be a sandbox with no physics seams (Rule 12). Research 2026-07-14 confirmed one uniform transport law spans cold↔plasma (RPA-LDA/Li-Petrasso/BPS) and that reactions cleanly factor into phase-space source terms feeding shared transport |
| 2026-07-19 | §2, §5, §7, §8, §13 | **v1.3: adaptive-dimension unified-grid pivot.** Runtime default becomes a **single 3-D multiphysics field simulation on the world-state grid** — the grid *is* the physicist; the accountant/physicist split and the reduced-dimension native-mesh + gather/scatter machinery are **superseded** (they were a Rule-12/13 seam). Dimensional reduction demoted to a **conservative, controlled-error adaptive optimization** (spectral azimuthal-mode truncation, *not* lossless) used only where geometry+state are symmetric. **3-D runtime solving and an in-sim CAD/B-rep kernel removed from §13 permanent-exclusion**; kinetic/PIC reserved for offline / one-time-grant calibration (not a runtime default, not banned). New goal **S8** (assumption-bounded prediction in untested regimes: minimal explicit research-grounded assumptions → one honest error bound; **resolve, don't accept**). Determinism (**S6**) relaxed to a negligible non-spiraling tolerance, with fixed-order deterministic reductions retained in chaotic/turbulent regimes. Compute budget stated: **~5000 GPU-desktop-hours + optional HPC (ACCESS) grant**. Materials → one **universal constitutive spine** over `M` (folds FND-7 with plasma EOS/opacity/stopping; Layer-2). Layer-2 docs (FND-2, FND-7, COUP-1/2, catalog) re-scoped in the cascade; the in-doc "needs-Ben" registers are purged (open questions asked before writing / during review, never parked in-doc) | Owner directive (2026-07-19): the point of the project is to **resolve** untested-regime uncertainty with a first-principles high-fidelity 3-D simulator, not to accept it; Rules 12/13 are supreme and the reduced-dim / per-material-table design violated their spirit |
| 2026-07-20 | §0, §4.2, §7.8, §8 | **v1.3.1: consistency-review reconciliations.** (a) §4.2 radiation reworded — one emit-once/transport-once operator over an energy-partitioned photon field (M1 + Sₙ split by energy, not by origin); (b) §7.8 pedigree corrected from a single "fraction of validated path × usage weight" scalar to the **PCMM weakest-link maturity vector** (aligns with the v1.3 predictive-V&V doctrine, META-1 §4.1); (c) §0/§8 **GPU reconciled** — the runtime 3-D field solver is GPU-accelerated (budget = GPU-desktop-hours), superseding the pre-pivot "GPU unavailable for physics" ruling (offline OpenMC stays CPU); (d) §8 "run the sweep in 3-D" reworded to "at adaptive/mostly-reduced dimensionality, full-3-D for anchors, UQ multi-fidelity." | Consistency-review pass (2026-07-20) found these Layer-1 statements were left un-reconciled by the v1.3 amendment |
| 2026-07-21 | §2 (S1) | **S1 overlap-band reframing.** S1's chemical-validation criterion changed from a hard **≤2% cutoff** to **p-box overlap with the RL10 reference band + reported Ferson area metric** (blind open-mode + calibrated closed-mode); ≤2% retained as the reported *target*, not a binary gate. Aligns S1 with the v1.3 predictive-V&V/p-box doctrine (S8; META-1 §4.1; VAL-1/VAL-2). The RL10 source (TM-107318) is itself a *validated model with a published error distribution* (Isp 440.3 s chamber sub-model vs 445.6 s cycle), so a single-point cutoff presumed a reference value that does not exist and collapsed an honest interval into a cliff. Turbopump confirmed a boundary object (delivers mass-flow; chamber pressure emerges); radiation left out of the chemical stage-1 by declared magnitude, not a seam. |
| 2026-08-13 | §5.3 | **Injector mixing tiered (R2).** The "injector mixing-efficiency boundary parameter" replacing 3-D turbulent combustion CFD is re-scoped as the **prior tier** of a tiered integration. A supplied η_c\* prior is measurement, not assumption, *only inside* the envelope of the injector family it was measured on — so its use is bounded by the existing COUP-7/FND-5 validity-envelope refusal (novel injectors mechanically cannot use it). Outside that envelope, mixing must be **computed**: the injection plane states the actual unmixed streams (composition, per-stream ṁ, momentum), the grid mixes them (SOLV-1 — same operator, resolution/azimuthal-mode lift; species advection `ρX_k` already carries the state), and η_c\* becomes an **output** whose residual uncertainty is per-quantity model-form on universal atomization/mixing closures. Resolved tier **architected now** (COUP-7 §3.2.1, SOLV-1 §3.4, OFFL-3 local-composition envelope), **built later** (rides the W4 two-phase extension for liquid injection). | Owner ruling (R2): a supplied injector η_c\* is the same category error as the supplied annihilation partition retired by R3 — it asserts an answer the instrument exists to compute, and a novel engine cannot supply one (blind-prediction contract §4.3). Tiering keeps the validated-family prior honest while making the extrapolation path computed, not assumed (Rule 12 / S8) |
| 2026-08-13 | §5.2, §9 | **Antimatter energy-partition reframed (R3).** Conversion efficiency and waste-heat load are **computed** by transporting Geant4 annihilation products (neutrinos escape; γ / charged pions·muons deposit or escape per geometry — SOLV-2/3), never assigned. The blanket **"15–30% model-form band" is retired**; model-form uncertainty is the **measured per-quantity Geant4 physics-list spread** (FTFP/INCL++/CHIPS ensemble members), validated against PS177/AEgIS. OFFL-4's headline "factor-of-several band" likewise demoted to the per-quantity spread it already computes. | Owner ruling (R3): a supplied efficiency/partition band is the same category error as a supplied injector η_c\* — it asserts the answer instead of computing it. A single blanket both **over-states** where data pins the quantity (prong multiplicity, ~few %) and **under-states** where lists diverge (U-238 antiproton fission probability, ~×several); only a per-quantity spread propagated through product transport is honest (Rule 12 / S8) |
| 2026-08-13 | §5.1, §7.5, §7.8, §8, §13, §15 | **v1.3.2: consistency-sweep reconciliations.** (a) §5.1 thermal-radiation row reworded to the one-operator form (M1 thermal band, SOLV-2; view-factor/radiosity = the *same* operator's optically-thin surface term); (b) §5.1 pulsed row: the private "flux-limited diffusion" replaced by the reused M1 closure (SOLV-5's Rule-12 ruling — a second radiation closure is a seam); (c) §7.5 injector entry updated to the R2 tiers; (d) §7.8 pedigree parenthetical expanded to the authoritative seven axes (META-1 §4.1); (e) §8 GPU-determinism parenthetical corrected to the S6 contract — it had **inverted** the mapping (S6: relaxed tolerance is the default, fixed-order is the *chaotic-regime* mandate; the same inversion is fixed in META-1 §2.1/§2.4 and FND-6 §3.6); (f) §13 Geant4 risk row's "~15%+" blanket updated to the R3 per-quantity spread; (g) §15 rows restored to date order | Consistency sweep (2026-08-13, two-agent corpus audit) found these Layer-1 statements stale against the reconciled Layer-2 corpus (SOLV-2/SOLV-5 one-operator radiation, R2/R3, META-1 §4.1) or self-contradictory (S6 vs §8) |
| 2026-08-14 | §5.1, §5.3, §7.2, §7.3, §9, §13 | **v1.4: review rulings (D-A…D-I, REVIEW_FINDINGS.md §1).** (a) World grid declared **natively cylindrical-structured** (config-declared symmetry axis) and adaptive dimensionality re-realized as **adaptive azimuthal resolution N_θ(r,z)** — conservative finite-volume ring coarsening with a mandatory guard resolution in adaptively-collapsible regions — superseding the v1.3 spectral-mode-truncation formulation (§5.1 banner, §7.2, §7.3); GPU-portable SoA layout mandated; CPU fixed-order reference solver precedes the GPU port (D-H). (b) **Sharp liquid-interface capture** added to the per-cell state design as a dormant capability (§7.2, D-B). (c) **One local wall-heat law** for all engines (Bartz demoted to a nozzle-envelope validation oracle) and **one universal offline-calibrated subgrid mixing closure** for the resolved injector tier (D-C, D-D) — capabilities are always general, never engine-specific (the tool/sandbox doctrine). (d) **Blind-mode rule** added to §9 (D-G, sharpened by Ben: no measured values of the engine under test, period). (e) §13 runtime-MC ban clarified to primary-method statistical transport; deterministic δf-marker fallback permitted (D-I). (f) P(WORKS)/reliability is a first-class reported result (D-F; COUP-4/COUP-5). | Review pass (2026-08-13, REVIEW_FINDINGS.md) found the v1.3 spectral mechanism unimplementable for shocked flow on the specified grid (S3–S5: Cartesian/cylindrical mismatch, unfireable re-expansion trigger, no shock-capturing story — precedents are smooth-field PIC codes); the ring-FV realization keeps the same physics and cost lever with proven shock-capable machinery. Liquid-interface, wall-law, mixing-closure, and blind-rule gaps were latent contract violations (§1.3 promises, Rule 12 seams, circular validation) surfaced by the same pass |
