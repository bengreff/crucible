# META-0 — Master Catalog

| Field | Value |
|---|---|
| **ID** | META-0 |
| **Family** | META |
| **Status** | Draft |
| **Depends on** | VISION_SCOPE.md (Layer 1, now v1.3) |
| **Version** | 0.5 (v1.3 unified-grid pivot) |

---

## 0. Purpose

This is the **table of contents for the entire simulator** and the authority on *which* design docs
exist and *what each one owns*. It is Layer 2's index and its change-control anchor. Layer 1
(`VISION_SCOPE.md`) says *what we build and why*; Layer 2 (these docs) says *how*, at spec level;
Layer 3 is code.

Rules:
- **No orphan design work.** If something needs deciding and doesn't fit an existing doc's scope, a
  doc is added *here first* (with an ID) before it is written.
- **This catalog outranks the individual docs on scope boundaries between docs.** If two docs both
  claim a decision, this file adjudicates who owns it.
- **VISION_SCOPE.md still outranks this file** on scope/doctrine (Layer 1 > Layer 2). A Layer-2 doc
  may not quietly reopen a Layer-1 boundary; it must trigger a Layer-1 Amendment (§15 there).

## 1. ID scheme

`<FAMILY>-<n>` where `n` is a **stable identifier, not an ordering**. New docs append the next free
number in their family; they never force a renumber. Filenames: `<ID>_<slug>.md`. Ordering is not
tracked as a formal dependency graph — docs are written **as needed, roughly following the writing
waves (§5)**; each doc's header notes the handful of docs worth reading first.

Families:

| Code | Family | Role |
|---|---|---|
| **META** | Meta | Catalog, philosophy, conventions, sources. Ambient to everything. |
| **FND** | Foundations / spine | Data model, grid+solver, geometry, config, tables, results, constitutive spine. Everything inherits these. |
| **COUP** | Coupling & orchestration | Time integration / operator split, conservation audit, two-stage execution, UQ, pedigree, boundary objects, the operator registry. |
| **SOLV** | Runtime physics operators | One doc per **operator** in the unified 3-D field solver (Rust, runtime). |
| **OFFL** | Offline pipeline | One doc per Python table-generation pipeline (OpenMC, Cantera, Geant4, DFT average-atom, …). |
| **VAL** | Validation & test | Predictive-V&V hierarchy, per-anchor specs, test/CI strategy. |

Status lifecycle (see META-2): **Draft → In Review → Reviewed → Frozen → Superseded**. A doc must be
**Reviewed** before any Layer-3 code implements it; **Frozen** means changes require the same ceremony
as a VISION_SCOPE amendment.

## 2. Reading order for a fresh session

`VISION_SCOPE.md` (v1.3) → this file (META-0) → META-1 (philosophy; Rules 12/13 are supreme) →
META-2 (conventions) → the specific module doc(s) you're working on (each header notes what to read
first). META-3 (sources) is a reference, consulted per datum.

---

## 3. The Catalog

> Scope = the one thing each doc owns. META-1/2/3 are ambient (every doc follows them). The **v1.3
> unified-grid pivot** collapsed the runtime into one 3-D field solver over the medium-state vector
> `M`; the old reduced-dimension "solver zoo" and the accountant/physicist split are retired
> (§6 log). SOLV is now a small set of **operators** on the one grid; couplers became terms of the
> grid's own coupled update.

### META — Meta (4)

| ID | Title | Scope (owns) |
|---|---|---|
| META-0 | Master Catalog | The doc set, IDs, writing waves, optional add-ons |
| META-1 | Design Philosophy & Principles | Governing principles; **Rules 12/13 supremacy**; determinism/reproducibility mandate (relaxed, v1.3); predictive-validation doctrine; precision & error philosophy; fail-loud doctrine |
| META-2 | Conventions: Docs, Notation & Code | Doc template, symbol/units notation, cross-ref syntax, forward-looking code conventions |
| META-3 | Sources, Data Provenance & Research Ledger | Paper-ready ledger of every constant/equation/table/tool with citation + derivation; data-hygiene schema |

### FND — Foundations / spine (7)

| ID | Title | Scope (owns) |
|---|---|---|
| FND-1 | Data & Uncertainty Model | The value-with-provenance-and-uncertainty type; units; `M`; how "a result is a distribution/p-box" is represented end to end |
| FND-2 | World-State Grid **& Unified Field Solver** | *(v1.3 rescoped)* The single 3-D grid that **evolves** all matter/fields over `M` and audits itself; cell state; adaptive **azimuthal resolution N_θ** (v1.4 ring-FV); static-topology sparse structure; determinism |
| FND-3 | Geometry & Voxelization | CSG (first-class) + STL import; material/volume-fraction assignment; interface refinement |
| FND-4 | Config Schema & Run Manifest | TOML schema; registry-driven sub-schemas; config↔manifest split; the regeneration key |
| FND-5 | Table Format, Loader & Interpolation | HDF5 table schema; provenance; validity-envelope enforcement; N-D interpolation + interpolation-error budget |
| FND-6 | Results Bundle & Reproducibility | Output bundle (HDF5+JSON); p-box serialization; provenance chain; build fingerprint; the S6 regeneration + verification contract |
| FND-7 | **Constitutive Spine** *(was Materials Property Library)* | *(v1.3 rescoped)* One continuous provider of EOS + transport + opacity + stopping over `M` for **arbitrary materials, cold solid → hot plasma**, no per-material seam; hybrid data calibration. Spec; generated offline by OFFL-5 |

### COUP — Coupling & orchestration (8; COUP-1 retired)

| ID | Title | Scope (owns) |
|---|---|---|
| COUP-1 | ~~Solver–Grid Binding & Conservative Source-Term Mapping~~ **RETIRED (v1.3)** | Absorbed — the unified grid *is* the solver, so there are no reduced-dimension native meshes to bind. Adaptive azimuthal resolution (FND-2 §3.4) replaces it. ID retained, retired, not reused |
| COUP-2 | Conservation Audit & Operator Coupling *(rescoped)* | The every-step global conservation/consistency audit over the one grid; the operator-coupling contract (which operator reads/writes which field of `U`); the **radiation-partition invariant** (emit-once/transport-once) |
| COUP-3 | Time Integration & Orchestrator | Global clock; **SDC-coupled IMEX** four-class operator split (explicit hydro/MHD + spatially-coupled implicit diffusion + cell-local stiff sources + declared global-ODE systems); pulsed-event sequencing; stiffness/convergence *(sub-cycling cut from v1, 2026-08-14; pseudo-transient mode DELETED 2026-08-19, VISION_SCOPE v1.5 — physical march only)* |
| COUP-4 | Two-Stage Execution & Halt/Failure Model | Stage-1 FUNCTION / Stage-2 LIFETIME control flow; halt conditions; verdict object |
| COUP-5 | UQ Engine | Latin-hypercube ensembles; correlation; Sobol; **multi-fidelity UQ** (cheap reduced-dim members anchored by few full-3-D runs, v1.3); parallel ensemble execution & determinism |
| COUP-6 | Pedigree Scoring | Per-result **PCMM-style maturity vector, reported by weakest-link min** (v1.3); PIRT assumption register; partition of the backbone map |
| COUP-7 | Boundary-Object Library & Registration Contract | The v1 boundary-object set (incl. EP performance tables & fusion-confinement source as source parameterizations); citation + validity-envelope + error-band registration |
| COUP-8 | Solver Interface & Operator Registry | The contract every **operator** implements to be config-composable with no per-concept code: declared reads/writes on `U`, required tables/boundary-objects, halt conditions, config-selectable registration |

### SOLV — Runtime unified-grid operators (8)

*(v1.3 — collapsed from the prior 16-module reduced-dimension "solver zoo"; merge mapping in §6. All
are operators/source-terms in the one field solver of FND-2, coupled by COUP-3, audited by COUP-2.)*

| ID | Title | Scope (owns) |
|---|---|---|
| SOLV-1 | Unified Field Operator | The conserved-variable system on the grid: compressible **reacting flow + two-phase + two-fluid/MHD + conduction**, one hyperbolic+parabolic update; HLLD/HLLC fluxes; **constrained transport** for ∇·B; config-time magnetostatics setup |
| SOLV-2 | Radiation-Transport Operator | Thermal radiation (**M1** two-moment) **+** deterministic **multigroup Sₙ** nuclear (neutron/photon) deposition — **one** transport operator; only the source spectrum differs (fission/fusion/activation/annihilation); no runtime Monte Carlo |
| SOLV-3 | Energetic-Particle-Transport Operator | The **one** friction/diffusion law for **all** fast charged particles (fission fragments, fusion α/p, annihilation π/μ, sputtered ions, fast e⁻) over `M`; slowing-down = stopping = self-heat; guiding-center↔Boris orbits; self-zeroing radiative losses |
| SOLV-4 | Reaction Sources | Fission point/few-group kinetics (incl. source-driven subcritical), fusion ⟨σv⟩ product sourcing, annihilation — **emit phase-space particle/energy sources** into SOLV-2/3; carry no private transport |
| SOLV-5 | Pulsed-Event Mode | The unified solver in a transient **event** configuration (Lagrangian rad-hydro for pellets/plasma slugs); returns impulse + energy partition + wall loading |
| SOLV-6 | Structural Margins | Analytic quasi-static hoop/thermal/burst checks → pass/fail + margin |
| SOLV-7 | Newtonian Outputs | Thrust/torque/momentum-flux integration over exit planes |
| SOLV-8 | Degradation Clocks (Stage 2) | Ablation recession, H₂ corrosion, burnup, fluence/DPA, decay heat — slow second-stage clocks |

### OFFL — Offline pipeline (6)

*(v1.3 — collapsed from 8; the stopping pipeline folds into the constitutive-spine pipeline, and the
two verification harnesses merge.)*

| ID | Title | Scope (owns) |
|---|---|---|
| OFFL-1 | OpenMC Transport Pipeline | Geometry classes; k-eff, deposition/dose kernels, kinetics params, reactivity sweeps; multigroup cross-sections for SOLV-2's Sₙ; fidelity tiering |
| OFFL-2 | Nuclear Data & SANDY UQ | ENDF/TENDL/FENDL handling; SANDY perturbed-library covariance sampling |
| OFFL-3 | Equilibrium Chemistry Pipeline | Cantera + CEA/RocketCEA property tables; B′ ablation tables |
| OFFL-4 | Annihilation Source Pipeline | Geant4 annihilation **source terms only** (products moved by SOLV-2/3); per-quantity physics-list spread (model-form, R3) |
| OFFL-5 | Constitutive-Spine Pipeline | **DFT average-atom backbone** (atoMEC-class) → consistent EOS + transport + opacity + stopping tables over `M` for arbitrary materials; QEOS/FEOS + Lee-More-Desjarlais + Stanton-Murillo closures; **GP data-calibration**. Generates the FND-7 spine |
| OFFL-6 | Verification-Oracle Harness | OpenFOAM (flow) & Athena++ (MHD) offline 3-D oracles + optional WarpX **kinetic calibration** (closure coefficients + bands); discrepancy → widened model-form band |

### VAL — Validation & test (3)

| ID | Title | Scope (owns) |
|---|---|---|
| VAL-1 | Predictive-V&V & Validation-Ladder Framework | *(v1.3)* The **4-tier hierarchy** (unit → benchmark → subsystem → system); **solution verification** (GCI); **predictive model-form UQ** (area metric → extrapolated band → p-box); how maturity maps to the PCMM pedigree (COUP-6) |
| VAL-2 | Anchor Specifications | Full spec of every anchor: **unit-physics** anchors (EOS Hugoniot, transport, radiation analytics, stopping) **+** engine anchors (RL10, NERVA/KRUSTY, NRX-A6, XE-Prime, NSTAR, VASIMR, rad-hydro analytics, antimatter, BEAVRS): data, targets, retrieval |
| VAL-3 | Test Strategy, CI & Conservation Harness | Test-first workflow; CI gates; conservation-audit harness; **determinism + ensemble-consistency** tests (v1.3); manufactured-solution suite |

**Total: META 4 + FND 7 + COUP 8 (7 active, COUP-1 retired) + SOLV 8 + OFFL 6 + VAL 3 = 36 entries (35 active).**

*Every SOLV operator implements the COUP-8 interface contract and is composed from config as pure data (Rule 13).*

---

## 4. Optional add-ons (side quests)

Ideas that are **not part of the main plan, the success criteria (S1–S8), or the writing waves.**
Pursued only if the main instrument is healthy and Ben wants to — parked here so the idea isn't lost
and, if taken up, slots in cleanly under the same razors and firewall. Adding one does not reopen any
VISION_SCOPE boundary.

| ID | Title | Scope | Status |
|---|---|---|---|
| SIDE-1 | Black-Hole-Drive Regime Add-on | Integrate Ben's `black_hole_paper/` BH-drive concept as an optional exotic "energy source" regime, reusing the existing field solver (SOLV-1), charged-particle operator (SOLV-3), and Newtonian outputs (SOLV-7); the Hawking/annihilation source enters as a boundary object / reaction source specified by the paper | Optional — not scheduled |

**SIDE-1 caveat.** The BH drive has **no experimental anchor** — it can never rise above validation-
ladder tier (ii), so it carries a permanently flagged pedigree and the widest declared bands, and sits
outside the "applied / validatable" spirit of the main campaigns. That's fine for a side quest; it
just means SIDE-1 never counts toward S1–S8 and never competes with backbone/validation work for
schedule. If pursued, it gets its own full DESIGN doc (and its ID graduates into the catalog proper).

---

## 5. Writing waves (rough execution order)

Sequencing per Ben's directive: **meta → spine → chemical vertical slice → fan out.** A loose guide to
what gets written when, not a rigid schedule.

| Wave | Docs | Rationale / exit |
|---|---|---|
| **W0 — Meta** *(done)* | META-0, META-1, META-2, META-3 | Conventions & principles fixed (now v1.3). |
| **W1 — Spine** *(done; rescoped for v1.3)* | FND-1, FND-2, FND-3, FND-5, FND-7, FND-6, FND-4, COUP-8 | Data model, the grid+solver, materials→constitutive spine, the table seam, the operator contract. FND-2 (grid=solver + adaptive dim) and FND-7 (constitutive spine) rewritten for the pivot; FND-4/6/COUP-8 (config/results/registry) are dimension-agnostic and survive. |
| **W2 — Chemical vertical slice** *(done 2026-07-21)* | COUP-2, COUP-3, COUP-4, COUP-7, SOLV-1 (reacting-flow config), SOLV-2 (thermal), SOLV-6, SOLV-7, OFFL-3, OFFL-5 (chemical/cold corner of the spine), VAL-1, VAL-2 (RL10 + unit anchors), VAL-3 | Everything to run the chemical engine end-to-end to the RL10 anchor (S1) on the unified solver. Proves the architecture on the most-validated regime. |
| **W3 — Nuclear leg** *(docs done 2026-07-21)* | OFFL-1, OFFL-2, SOLV-4 (fission kinetics), SOLV-2 (Sₙ nuclear deposition), SOLV-8, COUP-5, COUP-6, VAL-2 (NERVA/KRUSTY) | The NTP path + UQ + pedigree; validated vs NERVA/KRUSTY (S2). |
| **W4 — Advanced regimes fan-out** *(docs done 2026-07-21)* | SOLV-1 (two-phase/MHD extensions), SOLV-3, SOLV-4 (fusion/annihilation), SOLV-5; OFFL-4, OFFL-6; the rest of the constitutive spine (WDM/plasma corner); VAL-2 (remaining) | Two-phase, plasma/MHD, energetic-particle transport, fusion, pulsed, charged-particle — the rest of §5. *(Every catalog doc now has a v0.1 Draft; the SOLV-1 two-phase/MHD extensions and the FND-7/OFFL-5 WDM-plasma spine corner are explicitly deferred content **within** those written docs, to be filled as those regimes are built.)* |

These are *design-doc* waves, distinct from (though aligned with) the VISION_SCOPE §12 *build* roadmap.

---

## 6. Change control for the catalog

Adding, splitting, merging, or retiring a doc is an edit to this file with a dated line in the log
below. **Renumbering or reusing an ID is forbidden for any *instantiated* doc (one that has a file or
code).** IDs that were only ever catalog entries may be restructured under a Layer-1 amendment.

| Date | Change | Why |
|---|---|---|
| 2026-07-14 | v0.1 baseline: 44-doc catalog, dependency DAG, 5 writing waves | Layer-2 kickoff |
| 2026-07-14 | v0.2: added FND-7 and COUP-8 → 46 docs; tightened COUP-1/COUP-2 boundary; removed status dashboard and per-doc Owner field | Coverage/cleanup pass |
| 2026-07-14 | v0.3: removed the dependency graph and per-doc "Depends on" columns; added §4 optional add-ons with SIDE-1 | Ben directive |
| 2026-07-14 | v0.4: Rule-12 (no-seams) pass on calc-file scopes; SOLV-12 rescoped to medium-agnostic particle transport, SOLV-9 to product-sourcing, OFFL-6 to cold→plasma stopping, OFFL-4 to annihilation sources, SOLV-6 to any-origin deposition, SOLV-5 gains source-driven subcritical, SOLV-7 to one uniform energy-balance + polytropic nozzle; COUP-2 gains the radiation-partition invariant | Ben directive: sandbox, no seams; research 2026-07-14 |
| 2026-08-13 | v0.5.1: OFFL-4 scope cell updated to R3's per-quantity physics-list spread (was "model-form band") | Consistency sweep post-R3 |
| 2026-08-14 | **v0.6: review cycle closed — 29 docs promoted Draft → Reviewed (2026-08-14).** Full trail: 3-agent critical-path review (68 findings, REVIEW_FINDINGS.md) → Ben's rulings D-A…D-I (VISION_SCOPE v1.4: cylindrical grid + adaptive N_θ, liquid-interface capability, one wall law, universal SGS closure, staged spine, P(WORKS), blind rule, CPU-first, δf blessed) → 4-agent fix wave (all findings discharged) → verification pass (V1–V7 closed). Catalog scope rows aligned to v1.4 (FND-2 N_θ; COUP-3 four-class split, sub-cycling cut; COUP-1 wording). **Not promoted:** SOLV-5, OFFL-4 (deferred set — review before their build wave), COUP-1 (retired), META-0..3 (ambient). The META-0 §1 gate is now open: Layer-3 code may implement any Reviewed doc | Review cycle 2026-08-13/14 |
| 2026-07-19/20 | **v0.5: unified-grid pivot (VISION_SCOPE v1.3).** Runtime collapsed to one 3-D field solver over `M` (FND-2 = grid+solver); the reduced-dimension solver zoo and the accountant/physicist split retired. **SOLV 16 → 8 operators**, **OFFL 8 → 6**, **COUP-1 retired** (binding absorbed by adaptive azimuthal-mode reduction), COUP-2 rescoped (audit + operator coupling), COUP-3 central (SDC-IMEX), COUP-5 gains multi-fidelity UQ, COUP-6 → PCMM weakest-link pedigree. **FND-7 rescoped** materials→constitutive spine; **OFFL-5** rescoped EOS/opacity→DFT-average-atom spine generator (absorbs old OFFL-6 stopping); **OFFL-6** = merged verification/kinetic-oracle harness. **VAL-1** → predictive-V&V 4-tier framework; **VAL-2** gains unit-physics anchors. Prior SOLV/OFFL entries were never instantiated → restructured (ID-reuse rule waived for them, §1). **Merge mapping below.** | Ben directive (2026-07-19): resolve untested-regime uncertainty with a first-principles high-fidelity 3-D simulator; Rules 12/13 supreme; research 2026-07-19/20 |

### 6.1 v0.5 merge mapping (old reduced-dimension modules → new operators)

| Old (never instantiated) | New home |
|---|---|
| SOLV-1 Compressible Reacting Flow; SOLV-2 Two-Phase; SOLV-3 Conduction; SOLV-7 Plasma Energy Balance/Nozzle; SOLV-10 Magnetostatics | **SOLV-1 Unified Field Operator** (one conserved system; magnetostatics = config-time setup) |
| SOLV-4 Thermal Radiation; SOLV-6 Nuclear Heating/Dose | **SOLV-2 Radiation-Transport Operator** (M1 + multigroup Sₙ) |
| SOLV-12 Energetic-Particle Transport | **SOLV-3** (unchanged in spirit) |
| SOLV-5 Reactor Kinetics; SOLV-9 Fusion Reactions & Product Sourcing; (annihilation) | **SOLV-4 Reaction Sources** |
| SOLV-8 Pulsed Rad-Hydro | **SOLV-5 Pulsed-Event Mode** |
| SOLV-11 Structural; SOLV-15 Newtonian; SOLV-16 Degradation | **SOLV-6 / SOLV-7 / SOLV-8** |
| SOLV-13 EP Performance Tables; SOLV-14 Fusion-Confinement Source | **COUP-7 boundary objects** (source parameterizations) |
| OFFL-1/2/3/4 | **OFFL-1/2/3/4** (retained) |
| OFFL-5 EOS/Opacity; OFFL-6 Particle–Medium Stopping | **OFFL-5 Constitutive-Spine Pipeline** (one DFT-AA backbone) |
| OFFL-7 3-D Oracle Harness; OFFL-8 Mag-Nozzle Closure Calibration | **OFFL-6 Verification-Oracle Harness** |
