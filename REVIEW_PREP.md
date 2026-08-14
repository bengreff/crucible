# REVIEW_PREP — Pre-Review Checklist & Rubric

| Field | Value |
|---|---|
| **Status** | **CYCLE CLOSED 2026-08-14.** R1–R3 applied → sweep (~34 fixes) → 3-agent review (68 findings) → Ben's rulings D-A…D-I (all recs adopted + sharpenings; VISION_SCOPE v1.4) → 4-agent fix wave (68/68 discharged) → verification pass (V1–V7 closed) → **29 docs promoted Draft → Reviewed (2026-08-14)**. The META-0 gate to Layer-3 code is **OPEN**. Record: `REVIEW_FINDINGS.md`. Not promoted: SOLV-5/OFFL-4 (deferred — review before their build wave), META set (ambient). Residual low-stakes items: S1 §2 wording sign-off; oracle cadence (tiered default) confirm. |
| **Purpose** | Everything that must be true before an independent review pass on the Layer-2 corpus is worth launching |
| **Milestone driver** | Working simulator in **5–6 months** with **≥1 novel regime** |
| **Gate significance** | META-0: a doc must be **Reviewed** before any Layer-3 code implements it. This review is the gate to coding. |
| **Written** | 2026-08-13 (Ben away; prepared on best judgment, no authoritative doc altered) |

---

## 0. Why this exists

An independent review pass is expensive attention. It pays off only if reviewers find **new** problems
instead of (a) re-discovering issues already logged, or (b) reviewing a design that is about to change
under them. This file does three things so the pass is productive:

1. Names what must be **locked before** the review (baseline-changing decisions — §2).
2. Gives reviewers a **rubric** to check against your invariants, not their taste (§4).
3. Points review effort at the docs actually on the **5–6 month critical path** (§3), and hands over a
   **known-open-questions register** so reviewers hunt forward (§5).

Nothing in this file is a design decision I made unilaterally. The three rulings in §2 are yours; my
recommendation is recorded next to each.

---

## 1. Timeline reframe → what the milestone actually is

**Milestone-1 (5–6 months):** a running unified-grid solver that
- reproduces the **RL10 chemical** anchor as a distribution overlapping the reference p-box (VAL-2), and
- produces a pedigreed prediction for **one novel regime** (recommended: NTP/NERVA-class — validated,
  low-risk, exercises the whole nuclear leg), with
- **every result a distribution + pedigree** (COUP-5/COUP-6), never a bare number.

This is a subset of the full 35-doc program. Review depth must follow the split in §3 — deep review on
the ~20 critical-path docs, light or deferred on the advanced fan-out (pulsed events, annihilation,
plasma/magnetic-nozzle corners, full nuclear-data UQ).

**Recommended novel-regime strategy:** commit to NTP/NERVA as the guaranteed regime (it was tested → it
can be *validated*, de-risking the deadline), then take one cheap shot at a genuinely **untested** regime
(NSWR or fission-fragment) as a blind-prediction showcase riding the same machinery. That satisfies "≥1
novel regime" with a validated result and demonstrates the extrapolation thesis as upside.

---

## 2. Three rulings required before the review (baseline-changing — yours to make)

These change the authoritative baseline. Reviewing the affected docs before they are settled wastes the pass.

| # | Decision | Docs it changes | My recommendation |
|---|---|---|---|
| **R1** | **Novel-regime target** for milestone-1 | Sets review priority across the nuclear-leg docs | **NTP/NERVA floor + untested (NSWR/fission-fragment) stretch** |
| **R2** ✅ | **Injector η_c\* fork** — resolve mixing in-sim vs. boundary prior | VISION_SCOPE §5.3 (razor), COUP-7, SOLV-1, OFFL-3 | **RULED & APPLIED (2026-08-13): tiered integration.** Prior tier envelope-bounded (measurement, not assumption); resolved tier architected (COUP-7 §3.2.1 / SOLV-1 §3.4 / OFFL-3 §3.3), built later. |
| **R3** ✅ | **Antimatter partition** — blanket band | VISION_SCOPE §5.2/§9 + OFFL-4/META-3/SOLV-4 | **RULED & APPLIED: both blankets retired.** Efficiency/waste-heat *computed* by product transport; model-form = per-quantity physics-list spread. §15 amended 2026-08-13. |

R1–R3 are ruled (2026-08-13). **R3 is applied**; **R2's doc cascade is the one remaining baseline-changing write**
before review. R1 decides review depth: NTP/NERVA floor + NSWR stretch → the nuclear-leg docs get the deep pass.

---

## 3. Critical-path scoping (review-depth map)

**DEEP review (on the milestone-1 critical path):**

- Foundations: **FND-1** (data/uncertainty model), **FND-2** (symmetry/homogenization/geometry),
  **FND-7** (materials/constitutive spine) — plus the remaining FND foundations they inherit from.
- Solvers: **SOLV-1** (unified field operator), **SOLV-4** (reaction sources / fission), **SOLV-7**
  (Newtonian outputs), **SOLV-6** (structural margins), **SOLV-2** (radiation — matters for NTP).
- Coupling: **COUP-2** (conservation audit), **COUP-3** (SDC-IMEX time integration), **COUP-4**
  (two-stage/halt), **COUP-7** (boundary objects), **COUP-5** (UQ engine), **COUP-6** (pedigree).
- Offline: **OFFL-3** (equilibrium chemistry), **OFFL-5** (constitutive spine), **OFFL-1** (OpenMC — NTP).
- Validation: **VAL-1** (predictive V&V), **VAL-2** (anchors: RL10 + NERVA), **VAL-3** (test/CI/oracle).

**LIGHT review (needed but reducible for milestone-1):**

- **SOLV-3** (energetic-particle transport) — NTP power shape can come from OFFL-1 tables + SOLV-4; full
  charged-particle transport is not on the milestone-1 critical path.
- **SOLV-8** (degradation clocks) — needed for NTP lifetime (hot-H₂ corrosion, burnup, dpa); stage-2, so
  reviewable after the stage-1 core.
- **OFFL-2** (nuclear-data / SANDY UQ) — can be simplified (point-value bands) for the first milestone.
- **OFFL-6** (verification-oracle harness) — the OpenFOAM/Athena++ verification path is useful; the
  WarpX/plasma calibration path is not milestone-1.

**DEFERRED (not milestone-1):**

- **SOLV-5** (pulsed-event mode) — antimatter/fusion.
- **OFFL-4** (annihilation source) — antimatter.
- Deferred *content within* written docs: SOLV-1 two-phase/MHD extensions (except the injector resolved
  rung if R2 = integration), FND-7/OFFL-5 WDM-plasma spine corner.

---

## 4. Review rubric (what reviewers check against)

Reviewers grade every doc against **your invariants**, in priority order:

1. **Rule 12 (no seams).** One physical law over the local medium-state vector `M`; no `if(material)` /
   `if(regime)` branch. Physics self-zeroes by magnitude, never by a code path. **This is paramount** —
   any seam is a blocking finding.
2. **Rule 13 (sandbox).** Configs are pure data; no per-concept physics code. Reactions are sources.
3. **Reaction Razor.** Simulated only if energy/thrust flows through it; else a boundary object with
   {citation, validity envelope, error band}. Flag anything simulated that shouldn't be, or asserted
   that should be predicted.
4. **Failure-Mode Razor.** Mechanical/electrical failure = constraint check, not simulated physics.
5. **Distribution + pedigree.** No result is a bare number. Every closure coefficient is an
   `UncertainInput`; the p-box is never silently collapsed; pedigree headline = weakest-link min.
6. **Blind-prediction contract & computed-partition rule (§4.3 workflow).** Geometry + operating profile
   in → performance out. Flag any *per-engine measured input* a novel engine could not supply. **Sharper
   form (R2/R3):** never supply an efficiency or energy-partition *band* that could be computed by
   transporting the actual products/flow through the actual geometry — a supplied injector η_c\* and a
   supplied annihilation efficiency are the *same* category error. Efficiency is an **output**; its
   residual uncertainty is model-form on universal closures, carried **per-quantity, not as a blanket**.
7. **Interface contracts.** Every "consumes from X" must have a matching "provides" in X. Dangling
   interfaces are cheap to find now, a rewrite to find at implementation. (Mechanized in §6.)
8. **No misdirection of the next agent.** Ambiguity that could be built wrong is a finding, even if the
   doc is internally consistent.

Reviewers should report **new** issues; §5 is what is already known — don't re-report it.

---

## 5. Known-open-questions register (already logged — hunt forward, don't re-report)

| Item | State | Owner |
|---|---|---|
| Injector η_c\* — in-sim vs boundary prior | **Ruled & applied R2** (tiered integration; cascade written 2026-08-13) | — |
| Antimatter partition — blanket band | **Ruled & applied R3** (computed by transport; per-quantity spread); §15 amended 2026-08-13 | — |
| S1 pass criterion "≤2%" → reported target (not a binary gate) | Needs VISION_SCOPE §2 wording sign-off | Ben |
| Oracle cadence | Defaulted **tiered** (cheap analytic every commit, heavy codes at milestones); switchable | Confirm |
| Deferred content: SOLV-1 two-phase/MHD; FND-7/OFFL-5 WDM-plasma spine | Intentionally deferred, not milestone-1 | Noted |
| SOLV-3 §3.1 rung-3 **runtime δf markers** vs VISION_SCOPE §13's permanent exclusion of "runtime Monte Carlo transport" — are deterministic-seeded, bounded marker fallbacks banned MC or permitted? | **JUDGMENT — escalated to Ben** (sweep, 2026-08-13). Reviewers: don't re-report | Ben |
| Nothing is Reviewed/Frozen yet | All 35 docs are Draft; this pass moves the critical-path set toward Reviewed | — |

---

## 6. Pre-review execution checklist

- [x] **R1/R2/R3 ruled** (§2, 2026-08-13). **R3 applied** via §15 (VISION_SCOPE §5.2/§9, OFFL-4, META-3, SOLV-4).
- [x] **R2 doc cascade written** (2026-08-13) — §5.3 razor amendment (tiered mixing) + §15 log; COUP-7 v0.2 (tiered injector §3.2.1, sharpened emergent-quantity rule); SOLV-1 v0.2 (resolved-mixing rung); OFFL-3 v0.2 (local-composition envelope); META-3 `injector-cstar-eff`.
- [x] **Consistency + interface-contract sweep** (2026-08-13) — two agents over all 37 files; 38 raw
      findings (~32 unique); all mechanical fixes applied (dangling OFFL-1↔SOLV-8 products, FND-6↔COUP-5
      rng-key mismatch, FND-5 RSS'd-return vs pure-outer-loop, the **S6 determinism-mapping inversion**
      in META-1 §2.1/§2.4 + FND-6 §3.6 + VISION_SCOPE §8, stale IDs/versions/pointers, pre-R2/R3
      phrasing, VISION_SCOPE v1.3.2 batch). 1 judgment item escalated to Ben (§5: SOLV-3 δf markers).
- [x] **Rubric (§4) + register (§5) handed to reviewers** up front (2026-08-13).
- [x] **Critical-path scoping (§3) applied** — deep vs light vs deferred set per the R1 ruling (3 agents:
      spine+field-solver / orchestration+UQ / nuclear-leg+validation).
- [ ] **S1 wording + oracle cadence** confirmed (can ride inside the review — §5).

Once R1–R3 are ruled I can execute the sweep and apply the R2/R3 doc edits in one batched pass, then the
review launches against a stable, scoped baseline.
