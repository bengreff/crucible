# COUP-5 — UQ Engine

| Field | Value |
|---|---|
| **ID** | COUP-5 |
| **Family** | COUP (Coupling & orchestration) |
| **Status** | Reviewed (2026-08-14) |
| **Depends on** | FND-1 (sampling contract), VAL-1 (bands), FND-2/COUP-3 (fidelity levels), COUP-4 (verdicts), FND-6 |
| **Version** | 0.2 |

---

## 0. Purpose

COUP-5 is the **ensemble outer loop** that turns FND-1 `UncertainInput`s into a `ResultDistribution`/p-box:
correlated joint sampling, the **multi-fidelity estimator** (cheap reduced-dim members anchored by full-3-D),
Sobol' sensitivity, deterministic parallel execution, and the envelope-hit coverage report. It **implements
FND-1's contract** and **propagates VAL-1's bands** — it originates no uncertainty of its own (the pure-outer-loop
commitment, META-1 §4).

Read after FND-1 (the sampling contract + p-box types it executes) and VAL-1 (the model-form/numerical bands
it propagates).

## 1. Scope & razor ruling
**Owns:** the ensemble loop; joint correlated sampling; the MFMC multi-fidelity estimator; Sobol' sensitivity;
deterministic parallel execution; the coverage/envelope-hit report. **Defers:** the sampling *contract* + RNG
key + interval/p-box *types* → **FND-1**; the band *widths* (model-form, GCI, interpolation) → **VAL-1**/FND-5;
the reduced-dim/full-3-D *models* → **FND-2/COUP-3**; serialization → **FND-6**; pedigree → **COUP-6**.

**Razor ruling:** pure infrastructure — its representational choices (inverse-CDF, Iman–Conover, MFMC control
variates, interval enclosure) are exact/standard estimators, not model-form.

## 2. Interfaces & contracts
| Interface | Consumed by | Contract |
|---|---|---|
| **`ResultDistribution`/p-box** | FND-6, COUP-6 | MFMC-corrected summaries {mean, std, p2.5/16/50/84/97.5} as **intervals** (p-box) + fidelity-tagged member vector + RNG record — **conditional-on-WORKS, labeled** (§3.2.1) |
| **`P(WORKS)`** (§3.2.1) | FND-6, COUP-6, reporting | first-class reliability functional: survey-level MC estimate + binomial CI + epistemic interval `[P_lo, P_hi]` + dominant limiting mechanism (COUP-4 per-member verdicts) |
| **Sobol' indices** | COUP-6, reporting | first-order + total-effect, **grouped by uncertainty type**, with bootstrap-CI ranking |
| **Coverage report** | COUP-6 | per-member envelope hits → honest partial-ensemble reporting (FND-1 §3.7) |

**Invariant:** every draw is a pure function of its counter-based RNG key (Tier-1 determinism); the p-box is
**never collapsed** to one distribution; MFMC is applied **per functional** (a pooled histogram would be biased
by the low-fidelity shape); performance functionals are conditional-on-WORKS and every output says so (§3.2.1).

## 3. Method

### 3.1 Sampling (implements FND-1 §3.5)
**LHS** (McKay-Beckman-Conover) for propagation; inverse-CDF mapping; empirical sets by **joint index**;
**Iman–Conover** rank-correlation for parametric correlation groups (empirical sets carry correlation
intrinsically). All draws from one **counter-based Philox**, keyed on
`(master_seed, input_id, member_index, dimension, outer_epistemic_index, purpose_tag)` — a strict superset of
FND-1's key tuple (the extra fields index the double loop, Sobol' columns, and bootstrap replicates).
[META-3: `latin-hypercube`, `iman-conover`, `random123-philox`]

### 3.2 Multi-fidelity aleatory inner loop — MFMC
A fidelity level is a **config-PINNED global `N_θ` ceiling** (O17, v1.4): the level pins `N_θ^max` for every
member run at that level, with the FND-2 §3.4 adaptivity free to coarsen **below** the ceiling but never rise
above it — levels are **N_θ ceilings, not spectral mode counts**. Pinning is what makes the pilot ρ govern
production (an adaptively-drifting model is not the model the allocation was computed for). Levels are ordered
by decreasing correlation ρ_{1,i} to the full-`N_θ^max` HF (index 1) — **not different physics**, so no
Rule-12 seam, and the coarsening hierarchy gives naturally high ρ. [META-3: `mfmc-estimator`,
`multifidelity-uq`]
- **Estimator:** `ŝ_MFMC = ŝ₁^{(m₁)} + Σ_{i≥2} α_i(ŝ_i^{(m_i)} − ŝ_i^{(m_{i-1})})`, nested `m₁⊂…⊂m_k`.
- **Control coeff:** `α_i* = ρ_{1,i} σ₁/σ_i`.
- **Sample ratios:** `r_i = √[w₁(ρ²_{1,i}−ρ²_{1,i+1}) / (w_i(1−ρ²_{1,2}))]`, `ρ_{1,k+1}≔0`.
- **Admissibility filter:** keep the model set only if `1=|ρ_{1,1}|>|ρ_{1,2}|>…` **and**
  `w_{i-1}/w_i > (ρ²_{1,i-1}−ρ²_{1,i})/(ρ²_{1,i}−ρ²_{1,i+1})`; drop violators (a filter, not hand-tuning).
- **Pilot** (~15–20% budget) estimates ρ, σ; the allocation is computed by fixed-order reduction, integer
  `m_i` by a deterministic largest-remainder rule, then **frozen into the manifest** (reproducible allocation).
- **Per functional:** run MFMC separately for mean, variance, and each quantile (control-variate quantile
  estimator); the serialized member vector is fidelity-tagged; summaries are the corrected functionals — never
  a raw pooled histogram. **MFMC not MLMC** (the reduced-dim model is a projection, not a grid-refinement
  level); **ACV** is a future variance-squeeze behind the same control-variate interface.

### 3.2.1 Verdict semantics over halted members — P(WORKS) *(O14/D-F — owned here)*
MFMC's control-variate premise fails across the halt discontinuity (a halted member has **no QoI value**, and
LF↔HF correlation is undefined through the survival boundary). Ruling:
- **P(WORKS) is its own functional** — a first-class reported result (VISION_SCOPE v1.4, D-F) — estimated by
  **plain or stratified MC on the survey fidelity level** (the cheapest pinned `N_θ` ceiling), never by MFMC
  control variates; its MC sampling error carries a binomial (Wilson) interval [META-3: `binomial-ci`], and
  the **epistemic outer loop (§3.3) brackets it**: the report is the interval `[P_lo, P_hi](WORKS)` over the
  epistemic box.
- **Performance functionals are CONDITIONAL-ON-WORKS**: every mean/variance/quantile/p-box is computed over
  surviving members only and **every output (plot, summary, serialized field) is labeled**
  `conditional-on-WORKS` — an unconditional performance number over a partially-halted ensemble is
  meaningless and forbidden.
- **Verdict-discordant members** — where fidelity levels disagree on survival — are **promoted to the high-
  fidelity model** (HF verdict governs) and **excluded from the control variates** (their LF evaluations are
  not correlation-valid); the **induced-bias bound** of this exclusion (discordant fraction × functional
  range over discordant members) is computed and **recorded in the bundle/pedigree** — never silently absorbed.

### 3.3 Epistemic outer loop → p-box (double loop)
**The epistemic box is the union of ALL interval-typed inputs (O15):** VAL-1 model-form intervals + the GCI
numerical interval + FND-5 interpolation intervals **+ every boundary-object/config `UncertainInput` of the
FND-1 interval family** (COUP-7 §3.4) — all **consumed**, carried as **intervals**, never RSS'd-Gaussian.
The epistemic-vs-aleatory assignment **follows the FND-1 family**: interval-typed ⇒ this outer loop;
distribution-family ⇒ the §3.2 aleatory inner loop — never re-classified here. In particular the
**wall-function ±20–30% band (SOLV-1 §3.5) is EPISTEMIC** — treating the dominant expander-drive band as
aleatory would silently narrow the p-box (the false-confidence failure O15 names).
- **Monotone dims → vertex enclosure** (2^m corners; exact when the QoI is monotone). Monotonicity is
  **checked, not assumed** (physics argument + finite-difference sign check at vertices).
- **Non-monotone dims → interval-optimization enclosure** — the optimizer is a **deterministic
  bound-constrained global method with a fixed iteration structure** (DIRECT-class: fixed division depth,
  fixed evaluation budget, fixed tie-breaking order), its structure **frozen into the manifest** — never a
  stochastic or wall-clock-bounded search. LHS-outer serves only as a screening seed; raw LHS-outer does
  **not** enclose and would give a too-narrow p-box (false-confidence violation). [META-3: `direct-optimizer`]
- **Cost control (O16):** the epistemic dimension count is capped by the named budget `M_EPI_MAX` (default 12).
  Above it, a **screening/grouping rule** applies: dims whose screened influence (elementary-effects pass
  [META-3: `morris-screening`], or a cited negligibility argument) is provably below a stated fraction of the
  leading dim are **frozen at nominal**, and physically-kin dims may be **grouped** to one interval; every
  freeze/group is **pedigree-recorded** (COUP-6) — screening is a recorded modeling act, never silent.
- **Uncertified enclosure ⇒ conservative box-corner widening (O16):** if enclosure is not certified for a
  dim (the monotonicity check fails *and* the optimizer budget ends without its optimality gap closing), the
  reported envelope over that dim is the hull of **all** evaluations (box corners + every optimizer iterate)
  **widened outward by the named factor `K_EPI_WIDEN`** (default 1.5× the observed evaluation range beyond
  the hull), flagged in pedigree — never the possibly-too-narrow point result of an uncertified search.
- **p-box = per-quantile min/max envelope** over all outer evaluations of the inner CDF; never collapsed.
  [META-3: `pbox`, `false-confidence`]

### 3.4 Sobol' sensitivity
Dedicated **Saltelli-2010 radial design** on a scrambled-Sobol' base, cost `N(k+2)`; first-order S_i +
**Jansen-1999** total-effect S_Ti; **grouped by uncertainty type** (nuclear-data / model-form / correlation-BC
/ interpolation). **Convergence = bootstrap-CI on the ranking** with explicit tied-group handling (never a
spurious order). [META-3: `saltelli-jansen`]

### 3.5 Deterministic parallel execution
`rayon` `par_iter` over members (Philox keyed per §3.1); **fixed-order tree reduction** for all ensemble
statistics — retained even in chaotic regimes so thread count can't butterfly (META-1 §2). Per-member
**envelope-hit logging** → coverage report. [META-3: `gamer2-determinism`]

## 4. Coupling relationships
- **FND-1** contract implemented; **VAL-1** bands consumed (never recomputed); **FND-2/COUP-3** supply the
  fidelity ladder; **FND-6** serializes the p-box + the RNG record (FND-6 §3.2 — keys derivable, O22);
  **COUP-4** supplies the per-member verdicts §3.2.1 consumes; **COUP-6** consumes coverage + Sobol'
  attribution; **OFFL-2/COUP-7** are the empirical-set / boundary-object `UncertainInput` producers.

## 5. Uncertainty & validity
Originates none. Correctness (unbiased MFMC, enclosing p-box, determinism) is CI-gated (VAL-3). Validation-
ladder status: infrastructure.

## 6. Validation plan
1. Sampling determinism (1-vs-N-thread byte-identical, FND-1 test 1); marginal + Iman–Conover correlation
   reproduction.
2. **MFMC unbiasedness + variance reduction** on a synthetic 2-fidelity model; **admissibility filter** drops
   a planted non-qualifying LF model.
3. **p-box enclosure:** on a synthetic non-monotone band, vertex-only under-covers and the optimization
   enclosure recovers the true interval.
4. **Sobol' self-test** (Ishigami/Sobol-g analytic indices) + ranking stability + tied-group behavior.
5. Allocation determinism: frozen `m_i` reproduce across reruns.
6. **Verdict semantics (O14):** on a synthetic mixed-halt 2-fidelity ensemble, P(WORKS) is unbiased; planted
   verdict-discordant members are promoted to HF, excluded from control variates, and the recorded
   induced-bias bound covers the exclusion; all performance outputs carry the conditional-on-WORKS label.
7. **Epistemic hygiene (O15/O16):** an interval-typed boundary-object input lands in the outer box (never the
   inner ensemble); a planted negligible dim is frozen with a pedigree record; an uncertified-enclosure
   fixture reports the widened hull, flagged.

## 7. References
META-3 keys: `latin-hypercube`, `iman-conover`, `saltelli-jansen`, `mfmc-estimator`, `multifidelity-uq`,
`random123-philox`, `pbox`, `false-confidence`, `gamer2-determinism`, `sandy-samples`, `binomial-ci`,
`direct-optimizer`, `morris-screening`. Depends on FND-1, VAL-1, FND-2/COUP-3 (pinned `N_θ` ceilings),
COUP-4 (per-member verdicts), COUP-7 (interval-typed inputs), FND-6.

*(No open questions — LHS-for-propagation / Saltelli-radial-Sobol'-for-sensitivity, MFMC per-functional,
vertex+optimization p-box enclosure, frozen allocation all resolved 2026-07-21; FND-1 §3.5 clarified to match.)*

## 8. Change log
| Date | Version | Change |
|---|---|---|
| 2026-08-14 | 0.2 | **Review fix wave (O14/D-F, O15, O16, O17).** New §3.2.1: **verdict semantics over halted members** — P(WORKS) is its own first-class functional (survey-level plain/stratified MC + binomial CI + epistemic interval; never MFMC'd across the halt discontinuity); performance functionals conditional-on-WORKS, labeled on every output; verdict-discordant members promoted to HF and excluded from control variates with the induced-bias bound recorded (O14/D-F). §3.3: epistemic box = **union of ALL interval-typed inputs** incl. boundary-object/config intervals; assignment follows the FND-1 family; wall-function ±20–30% band epistemic (O15); `M_EPI_MAX` dimension budget + Morris-screened freeze/group rule (pedigree-recorded); deterministic DIRECT-class bound-constrained optimizer, fixed structure, manifest-frozen; uncertified enclosure ⇒ `K_EPI_WIDEN` box-corner widening, flagged (O16). §3.2: a fidelity level = a **config-pinned global N_θ ceiling** (adaptivity capped below; v1.4 — ceilings, not spectral mode counts) (O17). §2 P(WORKS) interface row + conditional-on-WORKS invariant; §6 items 6–7. |
| 2026-07-21 | 0.1 | Initial draft. Implements FND-1's sampling contract (LHS propagation + Iman–Conover + joint-index empirical sets); MFMC multi-fidelity estimator over FND-2 azimuthal-mode fidelity levels (formulas, admissibility filter, frozen allocation, per-functional); double-loop p-box (vertex enclosure for checked-monotone dims, interval-optimization for non-monotone, per-quantile envelope); Saltelli-radial-Sobol' + Jansen total-effect grouped by uncertainty type with bootstrap-CI ranking; Philox + fixed-order-reduction determinism; envelope-hit coverage report. |
