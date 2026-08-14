# VAL-3 — Test Strategy, CI & Conservation Harness

| Field | Value |
|---|---|
| **ID** | VAL-3 |
| **Family** | VAL (Validation & test) |
| **Status** | Reviewed (2026-08-14) |
| **Depends on** | VAL-1, VAL-2, COUP-2, FND-6; META-1 §2 |
| **Version** | 0.2 (2026-08-14 review fixes: N20, N21) |

---

## 0. Purpose

VAL-3 owns the **test-first workflow and the CI gates** that make the whole architecture's guarantees
*tested, not trusted*: the **method-of-manufactured-solutions** suite, the **conservation-audit harness**, the
**determinism + ensemble-consistency** tests, and the **differential-oracle harness** (the correctness engine —
port algorithms from open reference codes and keep them as oracles). It operationalizes the cadence rule:
**every session ends with tests green and work committed** (VISION_SCOPE §12).

Read after VAL-1 (the framework it runs) and META-1 §2 (the determinism mandate it enforces).

## 1. Scope & razor ruling
**Owns:** the test-first workflow; the CI gate set; the MMS suite; the conservation/determinism harnesses; the
**oracle-per-module** harness + tolerances. **Defers:** the framework recipes (GCI/area-metric) → **VAL-1**;
the anchor data → **VAL-2**; the regeneration gates → **FND-6 §3.6** (VAL-3 runs them); the conservation audit
itself → **COUP-2** (VAL-3 tests it).

**Razor ruling:** infrastructure. Its obligation: **no guarantee is believed until a standing test enforces
it** — determinism, conservation, order-of-accuracy, and oracle agreement are all build-breaking gates.

## 2. Interfaces & contracts
| Interface | Consumed by | Contract |
|---|---|---|
| **CI gate set** (§3.2) | the build | a fixed, deterministic-ordered battery; any red gate breaks the build |
| **Oracle harness** (§3.3) | every runtime operator | a module is "not done" until it reproduces its oracle to tolerance |
| **Determinism harness** (§3.4) | all numeric paths | intra-build fixed-order → byte-identity (incl. all chaotic-regime runs); declared relaxed paths (non-chaotic only) → declared tolerance + non-spiraling; ECT cross-build/platform only |

**Invariant:** a determinism/conservation/order regression is **ranked with correctness** (build-breaking,
META-1 §2.4).

## 3. Method

### 3.1 Test-first workflow
Solvers are built against **analytic / manufactured solutions before benchmarks before hardware** (the VAL-1
ladder). No session leaves a solver half-refactored or tests red (the antidote to context fragmentation,
VISION_SCOPE §12). Tests are co-located and named to the doc ID + version they cover (META-2 §4).

### 3.2 CI gate set
Fast gates every commit; heavy gates at milestones (§3.3). Per-commit gates: **MMS order-of-accuracy**;
**conservation audit** (closed-box + port closure, COUP-2); **determinism** (§3.4); **envelope/fail-loud**
(out-of-envelope refuses; injected NaN halts with diagnosis); **analytic anchors** (Sod, Su-Olson, Hugoniot,
PSTAR/ASTAR — VAL-2); the **FND-6 §3.6 regeneration gates** (determinism, thread-invariance, full-regeneration,
provenance-closure). **Milestone gate tier (N21):** the heavy code oracles (§3.3) **plus the 3-grid GCI
study** on the milestone anchor's headline QoIs (≥3 refined grids, r ≥ 1.3, **observed order checked against
formal** — out-of-asymptotic flags *unverified*; VAL-1 §3.2) **and the anchor area-metric/overlap computation**
(the §3.1 `d` score per VAL-2, blind + calibrated) — so the VAL-1 doctrine is a gate that actually runs, not
prose. All diagnoses in deterministic order. [META-3: `gci-roache`, `area-metric`]

### 3.3 Differential-oracle harness — tiered *(default; Ben's Q3 unanswered → tiered)*
For every runtime operator, port the algorithm from an open reference code and keep it as a **differential
oracle**. **Tiered execution:** cheap **analytic oracles run on every commit** (MMS, exact Sod, Su-Olson,
NIST tables — they catch most regressions fast); heavyweight **code oracles run offline/at milestones** to
cover regimes with no closed-form answer and to cross-check the analytic gates. (This is the tiered option;
switch to full-CI or one-time on Ben's call — it only changes *when* the heavy oracles run.)

| Module | Oracle | Checks | Tolerance |
|---|---|---|---|
| SOLV-1 hydro | exact Riemann + **Castro**/Athena++ | shock/contact/rarefaction, geometry | L1→0 at formal order; within oracle truncation |
| SOLV-1 + COUP-3 reacting coupling | **Castro** (SDC path) + **PeleC** (real chem) | stiff advection-reaction coupling | match reference to few % |
| whole PDE operator | **MMS** (self-oracle) | every discretized term | observed order within ~5–10% of formal |
| SOLV-2 radiation | **Su-Olson/Marshak** + **Quokka** | non-equilibrium rad-matter coupling | T_rad,T_mat within ~1–3% |
| FND-7/OFFL-5 EOS | principal Hugoniot + SESAME | shock-compression EOS | within experimental bars |
| SOLV-3 stopping | **NIST PSTAR/ASTAR** | dE/dx, range | ~1–5% (Bethe regime) |
| OFFL-3 chemistry | **CEA/RP-1311** cases | equilibrium performance | within manual tolerances |
| system | **RL10** (VAL-2) | end-to-end QoI | p-box overlap + area metric |

[META-3: `mms`, `sod-shock`, `su-olson`, `castro-source`, `radiation-m1`, `hugoniot-anchor`, `stopping-astar`]

### 3.4 Determinism harness — the reconciled S6 mapping *(N20)*
Three gates, mapping one-to-one onto META-1 §2.1/§2.4 (the pre-reconciliation inverted mapping is retired):
- **Intra-build, fixed-order paths — which includes ALL chaotic-regime runs** (fixed-order reductions are
  *mandatory* in chaotic regimes): rerun representative configs at **1 vs N threads** and across
  invocations; payloads must be **byte-identical** after canonicalization (excluding the FND-6 §3.7
  metadata zone). There is **no intra-build ECT** — within a build even chaotic results are byte-identical
  by construction.
- **Declared relaxed-reduction paths (non-chaotic regimes only) — the GPU gate (D-H):** the relaxed run
  (GPU throughput path) is compared against the **CPU fixed-order reference** on the same config: agreement
  within the **declared negligible tolerance**; the **non-spiraling gate** — the reference↔relaxed
  trajectory difference is monitored over the run and must damp or stay bounded, never trend growing; and
  **no verdict divergence, ever**. The harness also verifies the **load-time refusal** of a
  relaxed+chaotic config fires (the O21 declaration surface, FND-4/COUP-8).
- **Cross-build / cross-platform jobs only:** the **ensemble-consistency test (ECT)** — statistically
  indistinguishable from the accepted ensemble — never used intra-build.
A regression in any gate is build-breaking (META-1 §2.4). [META-3: `gamer2-determinism`, `ect-consistency`]

### 3.5 Determinism of the harness itself
All gates iterate ordered/sorted keys, emit diagnoses in deterministic order, and use tolerance-based
(ULP/relative) comparison against code oracles (which are not bit-identical to CRUCIBLE), reserving
byte-identity for the intra-build determinism gate.

## 4. Coupling relationships
- **VAL-1** provides GCI/area-metric/order recipes the gates apply; **VAL-2** provides anchor data; **COUP-2**
  provides the conservation audit under test; **FND-6** provides the regeneration gates VAL-3 runs; **COUP-6**
  consumes gate results into pedigree (a module below its oracle tolerance cannot claim its ladder tier).

## 5. Uncertainty & validity
VAL-3 originates no physics uncertainty — it *enforces* the guarantees others rely on. A failed gate is a
build break, not a caveat. Validation-ladder status: N/A (it is the test layer).

## 6. Validation plan
Self-consistency: (1) a seeded regression on a known-good solver stays green; (2) an injected order-reduction,
conservation leak, non-deterministic reduction, or oracle-tolerance breach each **turns the correct gate red**
(the harness is tested by deliberate faults); (3) diagnostic ordering is deterministic.

## 7. References
META-3 keys: `mms`, `sod-shock`, `su-olson`, `castro-source`, `radiation-m1`, `hugoniot-anchor`,
`stopping-astar`, `gci-roache`, `gamer2-determinism`, `ect-consistency`. Depends on VAL-1, VAL-2, COUP-2,
FND-6, META-1 §2.

*(No open questions — oracle-harness cadence is the tiered default pending Ben's Q3; switching tiers changes
only when heavy oracles run.)*

## 8. Change log
| Date | Version | Change |
|---|---|---|
| 2026-08-14 | 0.2 | Review fixes (N20, N21). **N20:** determinism harness rewritten to the reconciled S6 mapping — intra-build fixed-order paths (incl. ALL chaotic-regime runs) → byte-identical; declared relaxed-reduction paths (non-chaotic only, the GPU gate, D-H) → declared tolerance + non-spiraling gate vs the CPU fixed-order reference + verified relaxed+chaotic load refusal (O21); ECT reserved for cross-build/cross-platform jobs only (§3.4, §2). **N21:** milestone gate tier gains the 3-grid GCI study (observed-vs-formal order; out-of-asymptotic → unverified) and the anchor area-metric/overlap (`d`) computation (§3.2). |
| 2026-07-21 | 0.1 | Initial draft. Test-first workflow + CI gate set; MMS order-of-accuracy suite; conservation-audit and determinism/ensemble-consistency harnesses; the **tiered differential-oracle harness** (analytic oracles every commit, heavy code oracles Castro/PeleC/Quokka/Athena++ at milestones — Ben's default) with an oracle-per-module tolerance table; runs the FND-6 regeneration gates; self-tested by deliberate-fault injection. |
