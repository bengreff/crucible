# REVIEW_FINDINGS — Critical-Path Review Results (2026-08-13)

| Field | Value |
|---|---|
| **Status** | **CLOSED 2026-08-14.** 68 findings → rulings D-A…D-I (§1, all adopted 2026-08-14) → 4-agent fix wave → verification pass: **68/68 discharged**, all mirror pairs verified in full text (wall-law, η_c\* knockdown, (p,h,Z) chemistry, P(WORKS), kernel-mode partition, blind rule, determinism mapping), key closure clean (157 keys resolve), V1–V7 punch list closed (incl. the class-`G` global-ODE slot in COUP-3 §3.1 and the O1 mirrors). **29 docs promoted Draft → Reviewed (2026-08-14)** (META-0 §6 v0.6). The coding gate is open. |
| **Scope** | Per REVIEW_PREP §3 (R1 ruling): spine + field solver (S#), orchestration + UQ (O#), nuclear leg + validation (N#). Deferred set (SOLV-5, OFFL-4) not reviewed. |
| **Baseline reviewed** | Post-sweep corpus of 2026-08-13 (R2/R3 applied; v1.3.2; determinism reconciliation) |

## 0. Verdict table

| Doc | Verdict | Findings |
|---|---|---|
| FND-1 | PASS-WITH-FIXES | S1, S2 |
| FND-2 | **HOLD** | **S3, S4, S5** + S6, S7, S8, S9 |
| FND-3 | PASS-WITH-FIXES | S10 |
| FND-4 | PASS-WITH-FIXES | O20, O21 |
| FND-5 | PASS-WITH-FIXES | S11 |
| FND-6 | PASS-WITH-FIXES | O22, O23 |
| FND-7 | PASS-WITH-FIXES | S11–S14 (+ S15 pending Ben) |
| COUP-2 | PASS-WITH-FIXES | O7–O11 |
| COUP-3 | **HOLD** | **O1** + O2–O6 |
| COUP-4 | PASS-WITH-FIXES | O4, O12, O13, O14 |
| COUP-5 | **HOLD** | **O14** + O15–O17 |
| COUP-6 | PASS-WITH-FIXES | O18 |
| COUP-7 | PASS-WITH-FIXES | O10, O19 |
| COUP-8 | PASS-WITH-FIXES | O20 |
| SOLV-1 | PASS-WITH-FIXES | S16, S18, S19, S21 (+ S5, S17, S20 pending Ben) |
| SOLV-2 | **HOLD** | **N2** + N1, N3, N4 |
| SOLV-4 | PASS-WITH-FIXES | N5, N6, N7 |
| SOLV-6 | **HOLD** | **N8** + N9, N10 |
| SOLV-7 | PASS-WITH-FIXES | N11 |
| OFFL-1 | **HOLD** | N2, N12, N13, N14 |
| OFFL-3 | **HOLD** | **S22**, S23 |
| OFFL-5 | PASS-WITH-FIXES | S12, S14, S24 (+ S15 pending Ben) |
| VAL-1 | PASS-WITH-FIXES | N15 |
| VAL-2 | PASS-WITH-FIXES | N16, N17, N19 (+ N18 pending Ben) |
| VAL-3 | PASS-WITH-FIXES | N20, N21 |
| SOLV-3 / SOLV-8 / OFFL-2 / OFFL-6 | LIGHT-OK | — |

**Reviewer consensus:** the W2 chemical-slice architecture (emergent-p_c, UQ/p-box plumbing, m=0 solver, registry/config/results skeletons, COUP-2 audit concept, COUP-7 emergent-quantity design) is sound. The HOLDs are concentrated where the v1.3 pivot's mechanism was asserted rather than designed (FND-2 §3.4), where quasi-1-D thinking survived (OFFL-3), and where one-sentence stubs sit on the critical path (COUP-3's owned mechanisms, SOLV-2 Sₙ numerics, SOLV-4 kinetics closure). None requires abandoning the architecture; all require writing what is implicit — plus the §1 decisions.

---

## 1. DECISION REGISTER — **RULED by Ben, 2026-08-14. All recommendations adopted, with two sharpenings:**

> **D-A** ✅ cylindrical grid + adaptive ring-FV N_θ. **D-B** ✅ interface tracking added to `U` (dormant). **D-C** ✅ one local wall-heat law for all engines — *and generalized: duplicated definitions/ownership anywhere are defects; redundancy makes implementation harder, not safer.* **D-D** ✅ one universal closure, calibrated offline, "make it smart" (dynamic-coefficient class) — *and: never implement engine-specific features; every capability is general. CRUCIBLE is a TOOL — the codebase builds the instrument, the research uses it. Sandbox is key.* **D-E** ✅ Saha/QEOS analytic backbone in W2/W3. **D-F** ✅ P(WORKS)/reliability is a first-class result. **D-G** ✅ **sharpened beyond the rec:** a blind run consumes **no measured values of the engine under test at all** — inputs are what we would have for a fusion engine (design spec + universal physics + technology-class data from *other* hardware, banded); comparison against measurement happens only *after*. **D-H** ✅ CPU reference first. **D-I** ✅ δf markers blessed via §13 clarification. **E-1/E-2/E-3** unvetoed → adopted.

*(Original decision table retained below for the record.)*

| # | Decision | Findings | Recommendation |
|---|---|---|---|
| **D-A** | **Grid representation + m≥1 mechanism** (the FND-2 blockers). §3.2 specifies a Cartesian sparse brick tree; §3.4's spectral azimuthal modes require cylindrical — irreconcilable, coordinate system never declared. And no shock-capturing algorithm exists for nonlinear fluxes on spectral θ-modes (FBPIC/QPAD precedents are smooth-field PIC). | S3, S5 (+S4, S6, S7 follow) | Declare the world grid **natively cylindrical-structured** (config-declared symmetry axis per run, cell volumes ∝ r; FND-3 restated in that metric), and realize adaptive dimensionality as **adaptive per-region θ-resolution N_θ(r,z)** (conservative ring-FV — shock-capable, one operator, mode-count = cell-count; spectral truncation survives as the low-N_θ limit with a mandatory guard band ≥ 1 azimuthal DOF so re-expansion can trigger, S4). Alternative: keep true spectral θ + specified pseudo-spectral flux + validity restriction to smooth azimuthal content (research risk). Layer-1 §7.2 amendment either way. |
| **D-B** | **Liquid interfaces.** §1.3 promises NSWR fuel-in-water and the CNTR rotating film as fully-simulated multi-material cells, but α_k is static-geometry and the deferred two-phase extension is drift-flux (interpenetrating continua — cannot represent a sharp rotating film). | S9 | Add **dormant VOF/PLIC-style α_k interface transport** to the `U` design now (the honest cost of the §1.3 promise). Alternative: amend §1.3/VISION_SCOPE to mixture-fidelity + declared model-form band for free liquid interfaces. |
| **D-C** | **Wall-heat closure.** Bartz ownership is smeared across SOLV-1/COUP-2/COUP-7 (double-count risk on the expander-drive quantity), and Bartz is global-operand (throat D, p_c) — NTP channels will need a second correlation ⇒ prospective Rule-12 seam. | S16, S17, O10 | **One local wall-function closure** over `M` + local geometry (declared band; Bartz demoted to VAL-2 oracle for the nozzle envelope); wall faces selected **geometrically at config time** (data, not an `if(material)`), resolved viscous flux suppressed at those faces; delineation written into COUP-2 §3.5 (who evaluates h, which faces, which sub-step). |
| **D-D** | **Resolved-mixing tier closure.** Computed η_c\* at engine Re (~10⁷) requires a turbulence/SGS closure; none exists in the corpus — laminar-resolved mixing under-mixes by orders of magnitude, so the R2 "computed, not assumed" promise cannot be discharged as architected. | S20 | Adopt **one universal LES-class SGS closure** (e.g. dynamic Smagorinsky + turbulent species diffusivity) over `M`+grid-scale as the resolved tier's declared closure, per-quantity model-form band. Alternative: restrict the tier's validity envelope to genuinely resolvable mixing (stated in COUP-7 §3.2.1). |
| **D-E** | **Milestone spine coverage.** The W2 data-first spine's envelope = the data envelope = exactly the validity cliff FND-7 forbids; the NSWR stretch (ionizing steam) sits outside it with no scheduled provider (the §5.2 Saha/QEOS default appears in neither FND-7 nor OFFL-5's plan). | S15 | Pull the cheap **Saha/QEOS analytic tier into OFFL-5's W2/W3 deliverable** as the wide-envelope backbone (data-corner GP-calibrated on top); full DFT-AA later for the WDM valley; per-regime bands stated honestly. |
| **D-F** | **Ensemble verdict semantics.** MFMC is undefined over halted members (no QoI value; correlation premise fails across the halt discontinuity); what a reported distribution means when 30% of members halt is unspecified. | O14 | **P(WORKS) is its own functional** with outer-loop epistemic bounds; performance functionals **conditional-on-WORKS** (stated on every plot); P(WORKS) estimated by plain/stratified MC on the survey level, verdict-discordant members promoted to HF and excluded from control variates, induced-bias bound recorded. |
| **D-G** | **What "blind" means for RL10.** Blind open-mode currently consumes RL10's own measured η_c\* = 0.9892 (η_c\* is *defined* from the anchor's own c\* data) — partial circularity in the headline blind claim. | N18 | Blind mode uses the **coaxial-injector-family prior band** (±1–3%); RL10's own fitted value only in the **calibrated** closed-mode run; both scores labeled honestly in VAL-2 §3.1. |
| **D-H** | **GPU sequencing (ratify the standing rec).** No GPU execution model exists in FND-2/SOLV-1/COUP-3; a coding agent will build CPU-only or GPU-first unvalidated. | S8, O21 follows | **CPU fixed-order reference solver first** (the Tier-1 oracle), GPU port later as a **declared relaxed-reduction path** validated against it; FND-2 constrains the data layout to GPU-portable (SoA) now; config gains the determinism-mode block (O21). |
| **D-I** | **SOLV-3 runtime δf markers** vs §13's "no runtime Monte Carlo transport" (escalated by the sweep). | sweep | **Bless** bounded, deterministic counter-seeded markers via a §13 clarifying amendment (the ban targets statistical transport as primary method); alternative: demote rung 3 to offline-only. Not milestone-critical. |

**Resolved as engineering unless vetoed (per META-1 §6 escalation rule):**
- **E-1 (N2, kernel-vs-Sₙ partition):** precomputed kernels = the *precomputed solution mode of the same transport operator*, valid only within a tabulated geometry class (milestone NTP path); runtime Sₙ = the general mode; **exactly one mode per particle-class+band per run**, config-declared, COUP-2-asserted. FND-2 §3.4.1 / SOLV-2 §3.2 / OFFL-1 reworded to match.
- **E-2 (N7, NERVA kernel fidelity):** deposition kernels gain the **drum-angle axis** (cheap at tiered statistics); steady-state sufficiency argument written into SOLV-4 §5 leaning on OFFL-1's sparse coupled-Picard anchors.
- **E-3 (O2, radiation CFL):** **RSLA** (reduced-speed-of-light) with declared per-regime ĉ + PIRT-recorded band, matching the Quokka oracle; per-band Δt rule + sub-stepping in COUP-3 §3.4.

---

## 2. Findings — Spine + field solver (reviewer 1)

[S1] MAJOR MECH FND-1 §3.3 — `M` fields "present-or-null" force per-law null-branches (an `if(regime)` in disguise). Fix: all fields always populated with physically-degenerate values (cold neutral: n_e=0, ⟨Z⟩=0, T_e=T_i=T, B=0); laws are total functions of `M`.
[S2] MINOR MECH FND-1 §3.5 — sample-set selection "`k mod N` (or a seeded permutation)" unresolved. Fix: pin the seeded permutation keyed on the §3.5 RNG tuple.
[S3] BLOCKING JUDG FND-2 §3.2/§3.4 — Cartesian brick tree vs cylindrical spectral modes: irreconcilable discretizations; coordinate system, axis metric, mode storage never declared. → **D-A**.
[S4] BLOCKING MECH FND-2 §3.4 — re-expansion indicator A = Σ|û_m≥1|²/|û_0|² ≡ 0 in a collapsed region (modes not carried): can never fire; collapsed regions suppress symmetry-breaking forever. Fix: mandatory guard band (collapse = truncate to m_g ≥ 1, never to m=0 alone); trigger on guard-mode growth; error bound restated on guard amplitude. (Realized inside D-A.)
[S5] BLOCKING JUDG FND-2 §3.4 / SOLV-1 §3.3 — no nonlinear-flux/shock algorithm for m≥1 (PPM/HLLC/EMF not per-mode evaluable; precedents are smooth-field PIC; θ-FV fallback pre-rejected). → **D-A**.
[S6] MAJOR MECH FND-2 §3.4 — collapse semantics unspecified (conserved-projection thermalizes discarded KE; primitive-projection trips the audit). Fix: conserved-variable projection, thermalized ΔKE logged per event, included in declared truncation bound.
[S7] MAJOR MECH FND-2 §3.4 — symmetry indicator unimplementable: fields entering A, ε value/units (blows up as |û_0|→0), τ values, region granularity, cadence/dwell all unstated. Fix: normalized energy norm over all conserved fields, per-brick granularity, absolute+relative floor, named τ defaults with rationale, config-documented cadence.
[S8] MAJOR JUDG FND-2 §3.7/§3.8 — no GPU execution model anywhere (all thread-count language); GPU port becomes a rewrite. → **D-H**.
[S9] MAJOR JUDG FND-2 §3.3/§1.3 — α_k static vs promised flowing-liquid multi-material cells (NSWR/CNTR); no doc owns interface transport; drift-flux can't do a sharp rotating film. → **D-B**.
[S10] MINOR MECH FND-3 §3.1 — voxelization sampling has no sample count/convergence bound procedure. Fix: fixed jittered pattern (named N), fraction error derived into the manifest.
[S11] MAJOR MECH FND-5 §3.3 / FND-7 §3.2 — independently interpolating p,e,s,a breaks thermodynamic consistency/convexity (imaginary sound speeds near vapor dome). Fix: spine EOS tables store a Helmholtz potential F(ρ,T) per species; p/e/s/a by consistent differentiation of the interpolant; offline convexity/positivity audit in metadata.
[S12] MAJOR MECH FND-7 §3.7 / OFFL-5 §3.2 — actual table coordinates never stated (full `M` untabulatable); additive-volume mixture rule named without its algorithm. Fix: per-species F(ρ_s,T) + transport vs (ρ,T,⟨Z⟩); exact isobaric additive-volume iteration (fixed count); mix/derive order stated.
[S13] MAJOR MECH FND-7 §2 — spine contract "returns value + uncertainty band" violates pure-outer-loop. Fix: runtime spine(M) returns member-sampled scalars; bands are declarations consumed by COUP-5.
[S14] MAJOR MECH FND-7 §4 / OFFL-5 §3.4 — both still RSS interpolation error into the band (sweep fixed FND-5, missed these mirrors). Fix: FND-5 §3.4 wording (epistemic interval, never RSS'd).
[S15] MAJOR JUDG FND-7 §5 / OFFL-5 §3.1 — milestone spine = data-envelope cliff; NSWR stretch unprovided. → **D-E**.
[S16] MAJOR MECH SOLV-1 §3.5 — wall-heat ownership smeared (Bartz vs resolved F_visc vs COUP-2 CHT vs COUP-7 jacket). Fix per **D-C** delineation.
[S17] MAJOR JUDG SOLV-1 §3.5 — Bartz global-operand ⇒ per-geometry correlation branch looming. → **D-C**.
[S18] MAJOR MECH SOLV-1 §3.4 — η_c\* "applied as a c\* deficit" has no mechanism; output-side multiplication leaves field state inconsistent. Fix: in-solver combustion-completeness knockdown at source-term level (delivered c\* = η_c\*·c\*_ideal at anchor state); output-side multiplication forbidden. (= O19; both reviewers converge.)
[S19] MAJOR MECH SOLV-1 §3.4/§3.1 — advected composition under shifting equilibrium unspecified (advect-then-overwrite violates element conservation). Fix: shifting mode advects **elemental** fractions (equilibrium projection each step); frozen mode advects full species; frozen↔shifting = discrete epistemic dimension in COUP-5.
[S20] MAJOR JUDG SOLV-1 §3.4 — resolved-mixing tier has no turbulence closure; laminar-resolved mixing at Re~10⁷ is orders under-mixed. → **D-D**.
[S21] MAJOR MECH SOLV-1 §5/§6 — no resolution/cost statement for the RL10 anchor; no steady-state accelerator owned. Fix: stated anchor-run budget (grid, steps, cell-updates/s, ensemble size) + pseudo-transient continuation as an owned COUP-3 mode (with O6).
[S22] BLOCKING MECH OFFL-3 §2/§3.3 — expansion tables keyed by area/pressure ratio: a quasi-1-D coordinate that doesn't exist for a field-solver cell; naive wiring double-counts the expansion physics. Fix: re-parameterize runtime tables to local state (p, h, elemental composition/mixture fraction) consistent with S19; area-ratio tables retained only as SOLV-7/VAL-2 oracles.
[S23] MAJOR MECH OFFL-3 §2 — OFFL-3 transport table (μ,k,Pr) vs FND-7 spine: two providers of the same quantities (seam test fails at spec level). Fix: Cantera transport becomes an input feed to OFFL-5's spine assembly; OFFL-3→SOLV-1 transport row deleted.
[S24] MINOR MECH OFFL-5 §3 — no compute-cost/table-size estimate. Fix: add the cost paragraph (AA points × per-solve × elements; hours-to-days on 32 cores).

## 3. Findings — Orchestration + UQ (reviewer 2)

[O1] BLOCKING MECH COUP-3 §3.1/§3.3 — conduction/viscous fluxes classed "cell-local implicit," but ∇·(k∇T) is neighbor-coupled: a cell-local implicit diffusion solve does not exist; the Jacobi fallback loses stability exactly where conduction is stiff (fine wall cells — the RL10 closed-mode path). Fix: third operator class — spatially-coupled implicit diffusion via deterministic fixed-cycle solver (multigrid/CG) — or rule diffusion explicit with a stated diffusion-CFL; SOLV-1 §3.2/FND-2 §3.4 updated to match.
[O2] MAJOR MECH COUP-3 §3.4 — (E_rad,F_rad) explicit in U ⇒ light-speed CFL (~10³–10⁴× below fluid) whenever M1 is active; no RSLA/implicit ruling anywhere. Fix → **E-3** (RSLA, declared ĉ + band, sub-stepping).
[O3] MAJOR MECH COUP-3 §3.4 / COUP-7 §3.2 — closed-mode expander solve: two docs point at each other; iteration variable, placement in SDC, sweep count, relaxation unspecified. Fix: COUP-3 owns it — fixed-point on ṁ with Aitken acceleration, fixed sweep count, once per step after wall-exchange; ordering stated.
[O4] MAJOR MECH COUP-3 §2 / COUP-4 §3.2 — "converged (fixed-sweep)" self-contradictory; no acceptance test; no COUPLING_RESIDUAL halt; numerical non-convergence indistinguishable from physical won't-bootstrap. Fix: named residual-acceptance constant post-sweeps + new halt variant (diagnosis=numerical); won't-bootstrap defined as divergence/inconsistency test; contraction expectation stated for the RL10 loop.
[O5] MAJOR MECH COUP-3 §3.4 — sub-cycling is one sentence: no region rule, ratio, refluxing/audit interaction, SDC interaction. Fix: cut region sub-cycling from v1 (single global Δt) or fully specify the multi-rate scheme. Recommend cut.
[O6] MAJOR MECH COUP-3 §3.1 / COUP-4 §3.1 — pseudo-transient mode named but undefined; acoustic-CFL march to thermal steady state ≈ 10⁷–10⁸ steps. Fix: COUP-3 defines deterministic local-Δt pseudo-transient continuation, its audit semantics (audit applies to converged state), and cost vs §8 budget.
[O7] MAJOR MECH COUP-2 §2/§5 — audit tolerance: no number, no scaling law; §6.2 "~1e-10 relative" contradicts "absolute." Fix: named constant with derivation (ε_machine·√N_cells·‖U‖-scaled, absolute floor), manifest-recorded.
[O8] MAJOR MECH COUP-2 §3.1/§3.4 — Σ(sources) undefined: independent recomputation can't close under SDC (O(Δt²) splitting structure ≫ 1e-10); kernel deposition needs explicit escape-fraction ledger + offline normalization. Fix: audit = bookkeeping identity over **integrator-applied increments**; kernel escape as boundary-ledger term; normalization verified offline.
[O9] MAJOR MECH COUP-2 §3.1 — momentum audit can't close: no wall-force/mount-reaction ledger term; anchored-solid momentum rule absent. Fix: per-step impulse-on-anchored-material ledger term + SOLV-7 thrust cross-check.
[O10] MAJOR JUDG COUP-2 §3.5 (+SOLV-1/COUP-7) — wall-face closure ownership + face identification without if(material). → **D-C**.
[O11] MINOR MECH COUP-2 §3.3 — "runtime enforcement... caught in test" conflates mechanisms. Fix: test-time property checks + debug-build write asserts; release relies on Manifest-granted typed field access.
[O12] MAJOR MECH COUP-4 §3.1 — WORKS has no criterion (tolerance, dwell, horizon, steady-state test). Fix: named constants — per-quantity tolerance on commanded profile, dwell window, Stage-1 horizon with FAILED_TO_REACH halt; recorded in verdict.
[O13] MINOR MECH COUP-4 §3.2 — nuclear halts deferred but milestone includes NTP; excursion check needs the commanded envelope FND-4 also defers. Fix: tracked obligation — PKE-excursion halt + envelope source land before the nuclear coding wave.
[O14] BLOCKING JUDG COUP-5 §3.2 / COUP-4 §5 — MFMC undefined over halted members; reporting semantics of a 30%-halt ensemble unspecified. → **D-F**.
[O15] MAJOR MECH COUP-5 §3.3 — epistemic box omits interval-family `UncertainInput`s (boundary objects/config); the dominant Bartz ±20–30% band read as aleatory ⇒ silently narrowed p-box (false-confidence failure). Fix: epistemic box = union of ALL interval-typed inputs; COUP-7 §3.4 disambiguated (epistemic-vs-aleatory by FND-1 family).
[O16] MAJOR MECH COUP-5 §3.3 — outer-loop cost unbounded (2^m inner ensembles; optimizer unspecified; no fallback when enclosure uncertified). Fix: epistemic-dimension budget + screening/grouping (pedigree-recorded), deterministic bound-constrained optimizer (fixed structure, manifest-frozen), conservative box-corner widening when not certified.
[O17] MINOR MECH COUP-5 §3.2 — MFMC levels must be *pinned* mode ceilings (adaptivity capped below), else pilot ρ doesn't govern production. Fix: level = config-pinned global mode ceiling.
[O18] MAJOR MECH COUP-6 §3.1/§3.2 — grading function not computable: no maturity-level table, no rigor/independence rubric, θ unstated. Fix: per-axis PCMM 0–3 descriptors bound to concrete V&V product thresholds; θ named with rationale; independence as recorded enum (self/oracle/external).
[O19] MAJOR MECH COUP-7 §3.2.1 / SOLV-1 §3.4 — η_c\* application operator unspecified (= S18). Fix: upstream combustion-enthalpy efficiency so p_c, jacket ΔH, cycle closure all see the deficit; stated in both docs.
[O20] MAJOR MECH COUP-8 §3.3 / FND-4 [engine] — Require↔Provide "matching" undefined (name? kind? explicit binding?); two same-kind Provides auto-match ambiguously. Fix: explicit named bindings in [engine] (Require → provider instance); ambiguity/missing = load error listing candidates; [engine] grammar added to FND-4's deferred list with this semantics fixed now.
[O21] MAJOR MECH FND-4 §3.2 (+FND-6, META-1 §2) — no declaration surface for relaxed-vs-fixed-order paths or chaotic classification; nothing refuses relaxed+chaotic. Fix: config determinism block (fixed-order | relaxed) in manifest; per-mechanism chaotic-classification in COUP-8 Manifest; load-time refusal of chaotic+relaxed.
[O22] MINOR MECH FND-6 §3.1 — HDF5 summary sketch is scalar-valued (collapses the p-box the JSON tier keeps); `/rng_keys [N]` can't hold within-member key components. Fix: interval-valued [7]×[2] summary; store {master_seed, algorithm, keying-rule id} (keys derivable).
[O23] MINOR MECH FND-6 §3.2 — "aleatory ensemble at one epistemic setting" — which setting, what's retained per outer evaluation, where Sobol' draws sit: unstated. Fix: full members at nominal (id recorded); per-outer-evaluation summaries retained; Sobol' at nominal.

## 4. Findings — Nuclear leg + validation (reviewer 3)

[N1] MAJOR MECH SOLV-2 §3.2 — Sₙ numerics absent: no quadrature, no acceleration (source iteration stalls at scattering ratio→1: graphite/H₂ cores), no negative-flux fixup, no solve cadence/cost bound. Fix: level-symmetric S₈-class quadrature, DSA/CMFD acceleration, fixup named, quasi-static change-triggered cadence, ray-effect item in §6.
[N2] BLOCKING JUDG SOLV-2 §3.2 / OFFL-1 / FND-2 §3.4.1 — two transport owners for the nuclear band (runtime Sₙ vs precomputed kernels), no partition rule: Rule-12 seam-test failure + double-count/gap risk. → **E-1** (one operator, two solution modes, one per run, config-declared, COUP-2-asserted).
[N3] MAJOR MECH SOLV-2 §3.3 — M1↔Sₙ band edge: no value/selection rule; below-edge down-scatter disposal unspecified (energy-leak risk). Fix: edge (or per-run rule) stated; sub-edge out-scatter deposits to matter (KERMA discipline), re-enters only via thermal re-emission; §6.4 conservation test extended.
[N4] MAJOR MECH SOLV-2 §3.1/§5 — radiosity-vs-M1 ownership criterion unspecified; M1 crossing-beam failure untied to the NTP transparent-cavity case (2500–2800 K walls across H₂ — the pathological regime the milestone exercises). Fix: continuous optical-depth-based ownership (emit once, τ-blended per Rule 12); NTP-cavity PIRT entry + two-beam benchmark in §6.
[N5] MAJOR MECH SOLV-4 §3.2 / OFFL-1 §2 — neutron side of kinetics unwritten: no amplitude equation, no φ(r,t) representation, no adjoint/importance weighting for advected precursors (THE flowing-fuel β_eff physics of NSWR); OFFL-1 ships no shape/adjoint fields. Fix: explicit quasi-static closure (amplitude ODE + tabulated ψ(r) + ψ†(r) weighting for delayed source and drifted-precursor loss); OFFL-1 gains per-geometry-class shape+adjoint products.
[N6] MAJOR MECH SOLV-4 §3.2/§3.5 — kinetics is a GLOBAL stiff system (ρ = integral over core fields) but declared cell-local for COUP-3; no Δt control for fast transients. Fix: declared global-stiff-ODE slot in the SDC step (like the pulsed-event carve-out) + kinetics-driven Δt limiter.
[N7] MAJOR JUDG SOLV-4 §3.2/§5 — kernels-shape-static vs drum/feedback-moved power shape: ≤10% NERVA + XE-Prime map-shape risk undeclared. → **E-2** (drum-angle axis + sufficiency argument).
[N8] BLOCKING MECH SOLV-6 §3 — σ_therm = E·α·ΔT is uniaxial; constrained liner is biaxial EαΔT/(1−ν) (40–100% systematic error); scalar-summing orthogonal components is not a valid stress state. Margins drive halts ⇒ spurious DOESN'T-WORK verdicts. Fix: per-direction superposition with (1−ν) biaxial forms (Roark), margin vs a stated equivalent-stress criterion (von Mises).
[N9] MAJOR MECH SOLV-6 §2/§3 — D, t, "load-bearing cell" not derivable from per-cell state; no structural-component abstraction exists. Fix: config-time structural annotation (FND-3/FND-4): shell primitives R(z), t(z) from CSG, mapped to member cells; ΔT between paired surfaces of the same component.
[N10] MINOR MECH SOLV-6 §3 — safety factors as unsourced ranges. Fix: single sourced values (NASA-STD-5012-class) as named constants.
[N11] MAJOR MECH SOLV-7 §3.2/§3.3 — p_c convention unpinned (static vs stagnation, which plane) — moves the c\*/C_F split by ~1–2% (= the whole S1 target); "p_c = ṁ·c\*/A_t" circular as definition. Fix: area-averaged stagnation pressure at a declared reference plane matching TM-107318's convention (state which); equivalence demoted to consistency check.
[N12] MAJOR MECH OFFL-1 §3 — per-cell MT-301/MT-901 mixing double-counts/gaps photon energy at every regime boundary (non-conservative kernels). Fix: one globally consistent decomposition (neutron KERMA excluding photon production + transported-photon heating everywhere); range test reserved for the runtime product split (FND-2 §3.4.1).
[N13] MAJOR MECH OFFL-1 §2/§5 — "one re-run per realization" unbounded (300 × 1000-point sweep ≈ years); the separability assumption nobody declared. Fix: factorization specified — per-realization re-runs at nominal + extreme state points; state-indexed offset with declared separability check, PIRT-recorded; budget arithmetic shown.
[N14] MAJOR MECH OFFL-1 §2 — no photon multigroup library / (n,γ) production matrices, yet SOLV-2 transports γ of any origin. Fix: coupled n-γ multigroup product (photoatomic XS, production matrices, γ KERMA) + γ-heating verification check.
[N15] MAJOR MECH VAL-1 §3.3 — area-metric→regression→95% PI undefined at milestone anchor counts (n=1–4). Fix: small-n fallback ladder — n=1–2: ASME V&V20 u_val + declared sourced extrapolation multiplier; n≥3: the regression; crossover stated; fallback recorded in pedigree.
[N16] MAJOR MECH VAL-2 §3 — nuclear anchor specs absent (KRUSTY config/retrieval, NERVA state-point definitions, pass criteria) though R1 puts them in milestone-1; S2-thresholds-vs-overlap-doctrine interaction unstated. Fix: write them now; S2 thresholds as reported targets over a p-box (mirroring S1) — flag to Ben only if he wants binary gates instead.
[N17] MAJOR MECH VAL-2 §3.1 — the S1 pass criterion is not computable: no recipe for the reference p-box, no overlap/area-metric definition between two p-boxes. Fix: reference p-box = reference interval convolved with published error distribution; d = ∫max(0, F̲_pred−F̄_ref, F̲_ref−F̄_pred)dq; pass = nonzero overlap per quantile, d reported.
[N18] MAJOR JUDG VAL-2 §3.1/§3.2 — blind mode consumes RL10's own η_c\*. → **D-G**.
[N19] MINOR MECH VAL-2 §3.2 — geometry-of-record missing (contours, contraction ratio) — blind config can't be authored. Fix: add geometry retrieval (design-report contour tables), cached.
[N20] MAJOR MECH VAL-3 §3.4/§2 — determinism harness still encodes the pre-reconciliation inverted mapping (ECT for intra-build chaotic; GPU relaxed-path gate missing). Fix: rewrite to the reconciled contract; add the relaxed-path tolerance + non-spiraling gate.
[N21] MINOR MECH VAL-3 §3.2 — GCI + area-metric appear in no gate tier (doctrine that never runs). Fix: add to the milestone gate tier.

---

## 5. Fix-wave plan

**Stage 1 (independent of §1 decisions — can start immediately):** S1, S2, S10, S11, S12, S13, S14, S18/O19, S19, S21, S22, S23, S24, O1, O3–O9, O11–O13, O15–O18, O20, O22, O23, N1, N3, N5, N6, N8–N17, N19–N21, E-1, E-2, E-3.
**Stage 2 (gated on D-A…D-I):** FND-2 §3.2–§3.4 rewrite (D-A: S3–S7), α_k transport (D-B: S9), wall-heat closure rewrite (D-C: S16, S17, O10), resolved-tier closure (D-D: S20), OFFL-5/FND-7 spine plan (D-E: S15), COUP-5 §3.2 verdict semantics (D-F: O14), VAL-2 blind definition (D-G: N18), GPU/determinism declaration surface (D-H: S8 + O21), §13 amendment (D-I).
**Then:** re-verify the HOLD docs (FND-2, OFFL-3, COUP-3, COUP-5, SOLV-2, SOLV-6, OFFL-1), promote the critical-path set to **Reviewed**, and open the coding wave (M1–M2 skeleton per VISION_SCOPE §12).
