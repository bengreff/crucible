# META-1 — Design Philosophy & Principles

| Field | Value |
|---|---|
| **ID** | META-1 |
| **Family** | META |
| **Status** | Draft |
| **Depends on** | VISION_SCOPE.md |
| **Version** | 0.4 (2026-08-19: §2.5 GPU determinism policy — per-device bit-exactness confirmed by Ben, gather kernels/fixed-topology reductions/no physics atomics, cross-device = tolerance/ECT, >30%-cost escape hatch. Prior: 0.3 v1.3 pivot; 2026-08-13 S6 reconciliation — Principle 3, §2.1, §2.4) |

---

## 0. Purpose

The standing engineering principles every design doc and every line of code must satisfy. Where
VISION_SCOPE.md gives *scope doctrine* (what to simulate), this gives *build doctrine* (how to
simulate it correctly, reproducibly, and defensibly). When a design decision is genuinely
underdetermined by the physics, these principles break the tie. They are ranked: when two conflict,
the lower number wins.

---

## 1. The principles, ranked

> **Rules 12 and 13 are supreme and inviolable (Ben, 2026-07-19).** The no-seams mandate (Rule 12)
> and the sandbox / no-per-concept-code mandate (Rule 13) are **not** tiebreakers — they are hard
> architectural constraints on *every* design. A design that requires a physics seam or per-concept
> code is **wrong by definition** and is fixed in the core, never accommodated in the config. The
> numbered ranking below breaks ties among the *remaining* engineering principles; it never licenses
> a Rule-12/13 violation. Every design doc is audited against them (§the seam test, Rule 12).

1. **Correctness before completeness before speed.** A result that is wrong is worse than absent.
   A narrow result that is right and honestly bounded beats a broad one that is fudged. Performance
   is optimized only after correctness is demonstrated and only where a profiler says it matters.
   (Non-goal reaffirmed from VISION_SCOPE: no real-time requirement.)

2. **Every result is a distribution with a pedigree, never a bare number** (S3). Uncertainty is a
   first-class data type (FND-1), propagated, never bolted on afterward. A headline number without a
   band is a bug, not a simplification.

3. **Determinism & reproducibility by construction** (S6). See §2 — this is large enough to own a
   section. On a fixed build/target, any result is reproducible from `{config + pinned table versions
   + seed}`: **bit-exact on fixed-order paths** (mandatory in chaotic regimes; the default wherever
   the cost is acceptable), **within a negligible, provably non-spiraling tolerance** on declared
   relaxed-reduction paths (non-chaotic only); any table is regenerable from its pinned pipeline
   script + inputs.

4. **Fidelity doctrine: utter completeness for the reaction, zero simulation for what we don't care
   about.** The two razors (VISION_SCOPE §4.1) are applied at *design-doc* time, not discovered at
   code time. Every doc states, up front, what it simulates and what it treats as a boundary object,
   and cites the razor that drew the line.

5. **The grid is the solver** (VISION_SCOPE §7.1, v1.3). One unified 3-D field solver evolves all
   matter and fields over the medium-state vector `M` on the world-state grid **and** audits its own
   conservation. The earlier "grid is the accountant; reduced-dimension solvers are the physicists"
   split is **retired** — it was a Rule-12 seam (two discretizations reconciled by transfer
   machinery). Reduced dimensionality survives only as an *adaptive, conservative* projection
   (adaptive azimuthal resolution — conservative ring coarsening, v1.4) applied where geometry and state are symmetric: a
   controlled-error compute optimization within the one solver, never a separate physics.

6. **Fail loud, halt clean, never guess.** Out-of-envelope table access, conservation-audit
   violation, or loss of an essential function **halts with a diagnosis** (mechanism, location, time)
   — it never silently clamps, extrapolates, or continues on a fabricated value. No cascading-failure
   modeling: first essential loss ends the run (VISION_SCOPE §7.6).

7. **Provenance is not optional** (data hygiene, META-3). Every constant, correlation, table, and
   closure carries {value, units, uncertainty, source citation, retrieval date, derivation note}.
   An input a paper cannot cite is an input we do not use. **Open sources are archived to a local
   read-only cache** so retrieval stays reproducible years on; paywalled sources are cited by DOI
   *(Ben, 2026-07-14)*.

8. **Validity envelopes are enforced by machinery, not discipline** (VISION_SCOPE §7.8). Every table
   and correlation declares its domain of validity; interrogating outside it flags or refuses per
   config. This is the structural defense against the tool degenerating into a speculation generator.

9. **3-D by default; dimensionality follows the physics** (v1.3). The runtime solves in full 3-D and
   *adaptively, conservatively* collapses to 2-D-axisymmetric / 1-D (adaptive azimuthal resolution N_θ, v1.4) only
   where geometry and state are symmetric — a controlled-error optimization, never a fidelity
   ceiling. Compute reality (VISION_SCOPE §8): the backbone **sweep** runs reduced; full-3-D is
   reserved for anchors, validation, and symmetry-breaking cases, and UQ at full-3-D is **multi-
   fidelity** (cheap reduced-dim members anchored by a few full-3-D runs). Offline 3-D oracle runs
   (OFFL-6) remain a cross-check and calibration source, not the only 3-D.

10. **Test-first, tests-green-per-session.** Solvers are built against analytic / manufactured
    solutions before benchmarks before hardware (the validation ladder, VAL-1). No session ends with
    a solver half-refactored or tests red (VISION_SCOPE §12). This is also the antidote to
    Claude-session context fragmentation.

11. **Two-language purity.** Runtime = 100% Rust; offline = Python only where OpenMC/Cantera/Geant4
    APIs force it; the *only* seam is versioned HDF5 tables (VISION_SCOPE §6). No third language;
    Python never runs at simulation time; Rust never generates tables.

12. **One uniform rule per physical law, everywhere it applies — no seams.** *(Standing user
    directive; a defining property of the sandbox — reaffirmed and sharpened 2026-07-14.)* A physical
    law is written **once** and evaluated against a **local medium-state vector** `M` = (species &
    number densities, nuclear charge, free- & bound-electron density, T_e, T_i, degeneracy, B-field,
    flow, …). "Cold solid," "warm dense matter," and "hot magnetized plasma" are **not different
    models** — they are different corners of `M`. This forces, at design time:
    - **No `if(material)` / `if(regime)` branch ever selects the physical law.** Where closures must
      change with regime (e.g. bound- vs free-electron stopping), they blend by **continuous switching
      functions of `M`**, never a discrete label. A discontinuity in a computed quantity across a
      material or regime boundary is a bug, not an approximation.
    - **One quantity, one owner, one code path.** A fast alpha slowing in plasma and in a wall obey the
      *same* friction law; stopping (dE/dx in space) and slowing-down (dE/dt in time) are the *same*
      friction (`dE/dt = v·dE/dx`); a neutron born from fusion and from fission is the *same* transport
      operator with a different source spectrum. Terms that should vanish in a regime vanish by
      **physics** (mass/field scaling, e.g. ion synchrotron ∝ 1/m⁴), not by a branch.
    - **Seam test (apply to every solver scope):** if two code paths could compute the same physical
      quantity, for the same particle or cell, under different labels, they must be merged into one law
      over `M`. This test is run against every SOLV/OFFL scope in META-0.

13. **The simulator is a sandbox.** Geometries, reactor/engine configurations, and mechanism
    selections are **pure data** composed through one uniform path — CSG/STL geometry (FND-3) + TOML
    config (FND-4) + the mechanism registry (COUP-8). Defining a new geometry or a novel reactor
    configuration must require **no new physics code and introduce no new seam** (VISION_SCOPE §4.3,
    "no per-concept custom code"). Dozens of configurations share one physics core; a configuration
    that needs bespoke physics code is a design failure to be fixed in the **core**, not papered over
    in the config. This principle and Rule 12 are two faces of the same commitment.

---

## 2. The determinism & reproducibility mandate

**Statement.** *Everything is fully deterministic and bit-for-bit reproducible wherever it is
physically and computationally reasonable. Where bit-exactness is not reasonable, the residual
nondeterminism is (a) named, (b) bounded, and (c) prevented from reaching a headline result
un-quantified.* Determinism is a load-bearing feature of this project, not a nicety.

### 2.1 Two tiers of reproducibility

- **Tier 1 — Reproducibility within a negligible, non-spiraling tolerance** *(v1.3, required for the
  runtime).* Same binary/target/`{config, pinned tables, seed}` ⇒ results reproducible to a tolerance
  **orders of magnitude below the physics error bound**, with **no chaotic (butterfly) divergence** and
  **never a conflicting verdict**. The regime→guarantee mapping (S6): **fixed-order reductions are the
  default wherever their cost is acceptable** (always on CPU), giving **byte-identical** output at any
  thread count; a **declared relaxed-reduction path** (GPU throughput) may replace byte-identity with
  the negligible-tolerance contract **only in non-chaotic regimes**, where a round-off perturbation
  provably damps rather than amplifies. In **chaotic/turbulent regimes** — where round-off genuinely
  butterflies (real sensitive dependence, not a bug) — relaxed reductions are **forbidden**: fixed-order
  deterministic reductions are **mandatory, accepting the cost**, so that *within a build* even chaotic
  results are byte-identical and a run can neither spiral nor contradict itself. Only
  **cross-build/cross-platform** comparison of chaotic results falls back to **ensemble/statistical
  consistency** (ECT-style, §2.4). (Ben, 2026-07-19.)
- **Tier 2 — Scientific reproducibility** *(required where Tier 1 is unreasonable — primarily the
  offline Monte-Carlo pipelines).* Same `{script, seed, code+data versions, rank count}` ⇒ results
  identical within declared Monte-Carlo statistics. The *runtime never depends on Tier-2
  reproducibility*, because it consumes only the **frozen, versioned table artifacts** produced by
  the offline step — not live MC.

### 2.2 Enumerated sources of nondeterminism and the ruling on each

| Source | Ruling |
|---|---|
| **Floating-point non-associativity under parallel reduction** (the big one: `rayon` sums, ensemble reductions) | Reductions in numeric paths use a **fixed, order-independent reduction structure** (deterministic tree reduction / fixed chunking, or serial accumulation where cost allows). Parallelism must not change the arithmetic order. Tested by "N-thread vs 1-thread ⇒ identical bytes." |
| **RNG in parallel** (LHC sampling, MC view factors, any stochastic step) | **Counter-based / splittable PRNG** (e.g. seed → per-stream sub-seed by a fixed rule keyed on ensemble-member index, cell index, etc.). A member's draws depend only on its index and the master seed, never on execution order or thread assignment. RNG algorithm + seed recorded in the results bundle. |
| **Hash-map / set iteration order** | Forbidden in any numeric or output-ordering path. Iterate ordered containers or explicitly sorted keys. (A Rust lint/convention in META-2.) |
| **Wall-clock, PID, addresses, environment** | Never enter physics or output ordering. Timestamps are recorded as *metadata*, never used in computation. |
| **Compiler/target FP variance** (fma contraction, `-ffast-math`-equivalents, SIMD reassociation, x87) | Fast-math-style reassociation is **off** for physics. FMA-contraction policy is fixed and recorded. Bit-exactness is guaranteed *per build/target*; cross-target bit-exactness is **not** promised (declared Tier-1 boundary). The build fingerprint (compiler version, target triple, relevant flags) is recorded in every results bundle. |
| **Transcendental function libraries** (`sin`, `exp`, … differ across libm) | The libm/build is part of the recorded build fingerprint. Where a result is sensitive, prefer a pinned implementation. |
| **Iterative-solver convergence** (Picard/implicit steps) | Convergence criteria are **absolute + deterministic** (fixed tolerance, fixed max-iters, fixed iteration order). "Converged in a variable number of steps by wall-clock budget" is forbidden. |
| **Offline Monte Carlo (OpenMC/Geant4)** | Tier-2: seeded, fixed rank count, versions pinned. The *table* it emits carries its seed and statistics; downstream runtime reproducibility rests on the table hash, not on re-running MC. |
| **Table regeneration** | Deterministic given pinned inputs; provenance metadata (generator hash, library versions, seed) travels with the table (VISION_SCOPE §8, FND-5). |

### 2.3 What "reasonable" excludes (the honest boundary)

- Cross-CPU-architecture / cross-compiler **bit**-exactness is not required (only per-build). A
  result bundle records enough build fingerprint to reproduce Tier-1 on a matching build.
- Offline MC bit-exactness across differing MPI rank counts is not required (Tier-2 covers it).
- Wall-clock timing/throughput is never part of "a result."

### 2.4 Enforcement

Determinism is **tested, not trusted**: a standing CI check (VAL-3) reruns representative configs
at 1 vs N threads and across repeat invocations. Payloads on **fixed-order paths — which includes
every chaotic-regime run** — (physics + provenance, excluding recorded timestamps) must be
**byte-identical**; payloads on **declared relaxed-reduction paths** (non-chaotic regimes only) must
agree **within the declared negligible tolerance** and pass the **non-spiraling check**; **no
comparison may ever show verdict divergence**. Cross-platform jobs compare statistically — an
**ensemble-consistency test** (ECT-style — statistically indistinguishable from the accepted
ensemble) for chaotic regimes. A determinism regression — bitwise where required, a failed
tolerance / non-spiraling / consistency check elsewhere — is a build-breaking failure, ranked with
correctness.

### 2.5 GPU policy *(2026-08-19, Ben — plan-of-record confirmation)*

Bit-exact-per-build applies **per device**: the CPU fixed-order build and the GPU build are each internally
bit-exact (same binary + config + tables + seed ⇒ identical bytes at any thread/block count), including in
chaotic regimes — the §2.1 chaotic-regime mandate binds the GPU too. Implementation consequences, stated
here so kernels are designed to them rather than retrofitted: **gather-formulated kernels** (one writer per
cell — scatter-style accumulation via atomics is forbidden in physics paths), **fixed-topology tree
reductions** for all grid/ensemble statistics, no vendor library calls with nondeterministic internals.
**Cross-device (CPU↔GPU) bit-identity is not promised** (different FMA/libm) — cross-device verification is
tolerance-based on non-chaotic fixtures and ensemble-statistical in chaotic regimes (§2.4, ECT).
Golden-byte artifact diffs remain pinned to the dev host; GPU regression = same-build rerun identity.
**Escape hatch:** a specific kernel where determinism is measured to cost > ~30% goes to Ben for an explicit
ruling (non-chaotic paths only, per §2.1); it is never relaxed silently.

---

## 3. Precision & numerical-error philosophy

- **f64 for all physics state and accumulation.** f32 only ever for bulk storage/visualization where
  a doc explicitly justifies it and bounds the error. No implicit narrowing in numeric paths.
- **Every solver states its error budget.** Discretization order, expected convergence rate, and the
  grid/step at which truncation error drops below the physics uncertainty are stated in the doc and
  demonstrated by a convergence test (manufactured solution, VAL-1 rung i). A solver whose numerical
  error is comparable to its physics uncertainty must say so.
- **Conservation is checked to a stated tolerance every step** (mass, momentum, energy, and charge
  where relevant). The tolerance is an absolute number in the doc; exceeding it is a bug → halt
  (Principle 6). "Roughly conserved" is not acceptable phrasing anywhere.
- **Interpolation error is budgeted into UQ, not ignored** (FND-5). Table interpolation carries a
  quantified error contribution that flows into the result's band; it is never assumed negligible
  without a bound.
- **Tolerances are named constants with rationale, never magic numbers.** Each appears once, sourced,
  in config or a documented default (META-2, META-3).
- **Units are explicit and consistent.** SI internally; any non-SI at a boundary (e.g. psia in a
  citation) is converted at the edge with the conversion recorded (META-2 units convention).
- **Prefer formulations that are well-conditioned in the regime of interest.** Where a standard form
  loses significance (e.g. differences of large near-equal quantities, stiff source terms), the doc
  chooses and justifies a conditioned reformulation rather than discovering the blow-up in code.

---

## 4. Uncertainty & pedigree philosophy

- **Sources of uncertainty are typed** (FND-1): nuclear-data covariance (SANDY), model-form band
  (e.g. Geant4 lists, reduced-model closures), correlation/BC tolerance, interpolation error. Each is
  tracked with its type so Sobol attribution (COUP-5) can say *what* drives a spread, not just how big.
- **UQ is an outer loop over deterministic inner runs** *(architectural commitment).* Each ensemble
  member is a fully independent, deterministic simulation; all uncertainty lives in the ensemble
  sampling (Latin hypercube) and post-processing. The inner solver carries **no distribution
  arithmetic** — a major simplification the whole runtime relies on. (Confirmed reading of
  VISION_SCOPE §7.8; formalized in FND-1/COUP-5.)
- **Pedigree is reported with every headline result** (COUP-6): a **PCMM-style predictive-maturity
  vector reported by its weakest-link minimum — not a single averaged "fraction of validated path"
  scalar** (see §4.1). The backbone map is explicitly partitioned into "validated-physics finds" vs
  "extrapolated leads." Modules stuck at ladder rung (ii) (plasma, antimatter) are *permanently
  flagged*; the papers say so.
- **Model-form disputes become UQ, not hidden choices.** Where the literature disagrees (e.g. NSWR
  two-phase closures), competing closures run as model-form ensemble members and the spread is
  reported — the dispute is a result, not a decision we bury.

### 4.1 Predictive validation & honest bounds in untested regimes *(v1.3 — the basis of S8)*

The project's purpose is to **resolve** untested-regime uncertainty as far as the physics permits,
not accept it. That is earned, not asserted:

- **Hierarchical validation-by-parts.** Each *constituent* law (EOS, transport, opacity, stopping,
  and each operator) is validated against fundamental benchmarks on a **4-tier ladder** — unit →
  benchmark → subsystem → system (VAL-1). A composite prediction at an untested operating point is
  trusted only because every constituent is validated **and** the composition introduces no
  un-validated physics (which must be argued, and checked at the system tier wherever any data exist).
  This is *why* Rule 12 is load-bearing: **first-principles laws extrapolate defensibly; calibrated
  curve-fits do not** (Oberkampf-Trucano-Hirsch). A regime patchwork could never make this claim.
- **Numerical error is bounded and separated** (solution verification). Every headline result carries
  a grid-convergence (GCI) numerical band with the **observed order checked against the formal order**;
  a result outside the asymptotic range is flagged *unverified*, never reported with a false bound.
- **Model-form error is measured and honestly extrapolated.** An area-metric discrepancy is measured
  where data exist and regressed against a physical coordinate; its **95% prediction interval at the
  application condition** is the model-form band — which *widens with distance from data* by
  construction. A *calibrated* discrepancy term is never carried into extrapolation (it confounds with
  parameters and feigns confidence far from data).
- **The single honest bound is a p-box, never a lone distribution.** Aleatory inputs → distributions;
  epistemic (model-form, numerical) → intervals; propagated nested (epistemic outer, aleatory inner)
  into an interval-valued CDF. Collapsing to one probability invites the false-confidence theorem. A
  ≤10% bound and an order-of-magnitude bound are **equally successful** *if* each is the tightest the
  stated minimal assumptions allow, and is stated not hidden (S8).
- **Pedigree is a maturity vector, reported by weakest-link, never averaged** (PCMM / NASA-STD-7009
  style). Independent axes — physics-model fidelity, code verification, solution verification,
  validation, UQ, input pedigree, results robustness — each graded by *rigor + independence of
  assessment*. The headline is the **minimum** axis, not a mean (averaging incommensurable axes is the
  chain-strength fallacy). Every result carries a PIRT-style **assumption register** (phenomena ranked
  importance × state-of-knowledge), which is where the minimal explicit assumptions of S8 are recorded.

---

## 5. Interface & data philosophy

- **Boundary objects carry a contract or they don't exist** (VISION_SCOPE §7.5): citation + validity
  envelope + error band, all three, enforced at registration (COUP-7).
- **The HDF5 table is the only cross-language contract** (Principle 11). Table schema, provenance,
  and validity envelopes are specified once (FND-5) and honored by both sides.
- **Config is the single source of a run's intent** (FND-4): one TOML file fully determines a run
  (plus pinned table versions + seed). No hidden defaults that aren't documented; no behavior that
  depends on anything outside `{config, tables, seed, build fingerprint}`.
- **Conservation audit is a global invariant, not a solver's private business** (COUP-2): the
  orchestrator audits port-level mass/energy/momentum every step; a violation beyond tolerance is a
  bug and halts.

---

## 6. Process philosophy

- **Docs are the durable memory.** VISION_SCOPE + these Layer-2 docs are the session-start context
  that survives context fragmentation. A decision that isn't written in a doc did not happen.
- **Resolve at design time; escalate only genuine design/creativity forks.** *(Standing directive,
  2026-07-14 — supersedes the earlier "ask on any underdetermined decision.")* Every question a coder
  can answer better — units, data structures, numerics, algorithms, libraries, uncertainty
  representation, table formats, any **engineering / implementation / representation** choice — is
  **resolved by Claude** with sourced reasoning and written into the body as a decision. Only questions
  that **genuinely need Ben** are escalated: **scientific direction & scope** (which campaigns / regimes
  / concepts to prioritize), **risk appetite** on speculative physics, **product intent**,
  **aesthetic / naming**, and tradeoffs between competing *goals* (not competing implementations). When
  escalation is needed it happens **before the doc is written, batched** — never draft-then-ask. A rare
  item that needs an **experiment/analysis** to close is recorded with the analysis that will close it.
  Escalation happens **before the doc is written or during a review pass**, batched — never draft-then-ask,
  and **never parked in-doc as a tagged open-questions register** (Ben, 2026-07-19). A doc ships with its
  design forks already resolved or explicitly raised to Ben, not with a standing question table.
- **Resolved questions are deleted, not archived.** When an item is resolved, its decision + rationale
  goes into the relevant body section and the question is **removed**. Docs never accumulate
  resolved-Q&A tables, and **design forks are never parked in-doc** — they go to Ben before writing or
  during a review pass (META-2 §5, v1.3). The only in-doc open item is a `needs-analysis` TODO; a doc
  with none omits any register entirely.
- **Descope down the ladder, never through the firewall.** Under schedule pressure, apply the
  VISION_SCOPE §13 descope ladder; never descope the backbone integration (S4), validation anchors
  (S1/S2), UQ-on-every-result (S3), two-stage execution (S5), or the sensitive-physics firewall (§10).
- **Amendments are cheap but logged.** Changing a Reviewed/Frozen doc is fine when warranted, via its
  own change log; changing a VISION_SCOPE boundary requires the §15 ceremony.
