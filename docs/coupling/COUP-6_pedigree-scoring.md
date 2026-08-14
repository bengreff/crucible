# COUP-6 — Pedigree Scoring

| Field | Value |
|---|---|
| **ID** | COUP-6 |
| **Family** | COUP (Coupling & orchestration) |
| **Status** | Reviewed (2026-08-14) |
| **Depends on** | META-1 §4.1 (the seven axes), VAL-1 (V&V products), COUP-5 (UQ), FND-6 |
| **Version** | 0.3 |

---

## 0. Purpose

COUP-6 attaches to every headline result a **PCMM-style maturity vector over META-1 §4.1's seven axes,
headlined by the weakest-link minimum**, and partitions the backbone/cross-regime map into
**validated-physics finds** vs **extrapolated leads**. It **scores** the products VAL-1/COUP-5 compute — it
runs **no V&V machinery itself** and never averages axes.

Read after META-1 §4.1 (the authoritative axis set + weakest-link rule) and VAL-1 (whose GCI/area-metric/PIRT
outputs it consumes).

## 1. Scope & razor ruling
**Owns:** the seven-axis maturity vector + weakest-link headline; the finds-vs-leads derived predicate; the
per-result assumption register (a *reference* into VAL-1's PIRT). **Defers:** all band computation →
VAL-1/COUP-5; PIRT *authoring* → VAL-1; serialization → FND-6; validation data → VAL-2.

**Razor ruling:** pure infrastructure — no physics and **no V&V machinery**. Its obligation is an honest,
non-collapsible maturity report (the false-confidence guard).

## 2. Interfaces & contracts
| Interface | Consumed by | Contract |
|---|---|---|
| **Pedigree object** | FND-6, reporting | seven-axis vector [{axis, grade, evidence_handle, independence}] + headline **min** + finds/leads label + assumption register |
| **Finds/leads label** | backbone-map reporting | derived predicate per result/path (§3) |

**Invariant:** headline = **minimum** axis (never a mean); each axis is fed by **one distinct** V&V product
(so nothing is averaged); a rung-(ii) module is **hard-wired to lead** regardless of a favorable numeric.

## 3. Method

### 3.1 The seven axes, each fed by one product
[META-3: `pcmm`, `nasa-std-7009`]

| Axis (META-1 §4.1) | Graded from | Key |
|---|---|---|
| physics-model fidelity | Rule-12 first-principles vs calibrated-fit status of each constituent law on the path + physics completeness + PIRT-weighted coverage | `vv-oberkampf-roy`, `pirt` |
| code verification | MMS / order-of-accuracy pass | `mms` |
| solution verification | GCI band + observed-vs-formal order within the asymptotic range | `gci-roache` |
| validation | 4-tier ladder tier achieved + area-metric within data range + `u_val` | `area-metric`, `asme-vv20` |
| UQ | COUP-5 completeness: all typed sources sampled, p-box non-collapsed, Sobol attribution present | `pbox` |
| input pedigree | `ProvenanceRef` + META-3 star-field completeness + covariance-vs-point data quality + OFFL-2's non-Gaussianity diagnostic & dominant-nuclide truncation justification | (FND-1/META-3/OFFL-2) |
| results robustness | p-box width vs extrapolation distance + model-form ensemble spread + envelope-hit fraction + false-confidence guard satisfied | `false-confidence` |

Each axis → a maturity level **graded by rigor + independence of assessment** (the NASA-STD-7009 grading
discipline atop META-1's seven axes). **Headline = min** (chain-strength fallacy forbids a mean); the full
vector is still reported (FND-6 serializes it). **Distinctness guard:** physics-model fidelity is kept
separate from validation — a first-principles model can be high-fidelity yet poorly validated (exactly the
plasma/antimatter case, permanently flagged).

**Per-axis maturity levels (O18) — PCMM 0–3, each level bound to a concrete V&V product threshold** (the
grading is *computable* from the products above, never judgment-per-result):

| Axis | 0 | 1 | 2 | 3 |
|---|---|---|---|---|
| physics-model fidelity | judgment/heuristic closure on the path | calibrated-fit closures, bands declared | first-principles laws over `M`; every high-importance PIRT phenomenon covered with declared model-form | 2 + per-quantity model-form measured against independent data (physics-list-spread style) |
| code verification | none/smoke only | regression/benchmark tests pass | MMS order-of-accuracy verified for every operator term class in use | 2 + coupled-operator MMS, standing CI |
| solution verification | single grid, no estimate | qualitative grid check | GCI band formed **in the asymptotic range** (observed order within 0.5 of formal) | 2 + numerical band ≤ a stated fraction (default 1/3) of the physics band |
| validation | ladder rung (i) only | rung (ii): reference-code/benchmark agreement within stated bounds | rung (iii): area metric vs a hardware anchor, formed **within data range**, nonzero overlap | 2 + ≥2 independent anchors + `u_val` quantified (ASME V&V20) |
| UQ | point values | partial sources sampled | all typed sources sampled; p-box non-collapsed; Sobol' attribution present | 2 + verdict/discordance handling recorded + clean coverage report |
| input pedigree | uncited value present (registration forbids — flag state) | citations complete, point values | full `ProvenanceRef` + META-3 star fields + bands on all inputs | 2 + covariance-grade data with non-Gaussianity diagnostic + truncation justification (OFFL-2) |
| results robustness | single run | ensemble spread reported | p-box stable vs extrapolation distance; envelope-hit fraction below the FND-1 §3.7 threshold; false-confidence guard satisfied | 2 + model-form ensemble spread bracketed by independent closures |

**Independence is a recorded enum per axis (O18):** `independence ∈ {self, oracle, external}` — `self` =
graded from the project's own tests; `oracle` = an independent reference code/dataset cross-check (OFFL-6,
VAL-2 oracles); `external` = independently published benchmark or external review. The pedigree object
carries `(level, independence, evidence_handle)` per axis; levels 2–3 on the validation and solution-
verification axes require `independence ≥ oracle` by construction of their products.

### 3.2 Finds-vs-leads (derived predicate)
`find ⟺ (pedigree_min ≥ θ) ∧ (no high-importance/low-knowledge PIRT phenomenon on the path) ∧ (area-metric
formed within data range, not extrapolated)`; else **lead**. **`θ = THETA_FIND = 2`** (named constant).
Rationale: level 2 is the first grade at which every axis rests on quantified, application-domain evidence
(GCI in the asymptotic range, in-data-range area metric, complete sampled UQ — the §3.1 table); requiring 3
would empty the finds set outside the anchored regimes, and 1 admits calibrated-but-unvalidated paths —
exactly what "extrapolated lead" exists to label.
**Honesty override:** any module at ladder rung (ii) (plasma, antimatter — "permanently flagged", META-1
§4.1) is hard-wired to **lead** regardless of a favorable numeric (the triple guard already forces this;
recorded clause-by-clause so the label is auditable). The predicate is *computed*, not hand-labeled (Rule 13).

### 3.3 Assumption register & determinism
Per result, COUP-6 emits the S8 minimal-assumptions record by **referencing** the VAL-1 PIRT entries on the
path (importance × knowledge) — it does **not** re-author the PIRT. The whole object is a pure deterministic
function of its inputs (no RNG) → CI-checkable.

## 4. Coupling relationships
- **META-1 §4.1** fixes the seven axes; **VAL-1** supplies ladder tier / GCI / area-metric / u_val / PIRT (all
  consumed, never recomputed); **COUP-5** supplies coverage + Sobol' + p-box completeness; **FND-1/META-3**
  supply provenance completeness; **OFFL-2** supplies the non-Gaussianity diagnostic + dominant-nuclide
  truncation justification for the input-pedigree axis; **FND-6** serializes; **VAL-2** supplies the anchor
  results.

## 5. Uncertainty & validity
Originates none; it is the frame others are graded in. Determinism by construction. Infrastructure rung.

## 6. Validation plan
1. **Weakest-link:** headline equals the min axis on a constructed vector, never the mean.
2. **Distinctness:** a high-fidelity/low-validation fixture (antimatter-like) scores high physics-fidelity, low
   validation, classifies as **lead**.
3. **Finds/leads predicate:** fixtures crossing θ / PIRT / data-range boundaries flip correctly; a rung-(ii)
   override forces lead despite good numerics.
4. **No-recompute:** given VAL-1/COUP-5 outputs, COUP-6 reproduces the pedigree without invoking any
   GCI/area-metric routine (guards the boundary).

## 7. References
META-3 keys: `pcmm`, `nasa-std-7009`, `pirt`, `vv-oberkampf-roy`, `gci-roache`, `area-metric`, `asme-vv20`,
`pbox`, `false-confidence`, `mms`. Depends on META-1 §4.1, VAL-1, COUP-5, FND-6.

*(No open questions — seven axes per META-1 §4.1 (NASA's eight rejected; its grading discipline kept);
PIRT owned by VAL-1 and referenced here; finds/leads a computed predicate with a rung-(ii) override. Resolved
2026-07-21.)*

## 8. Change log
| Date | Version | Change |
|---|---|---|
| 2026-08-14 | 0.3 | **Review fix wave (O18).** §3.1: the grading function made computable — **per-axis PCMM 0–3 maturity-level table**, each level bound to a concrete V&V product threshold (asymptotic-range GCI / in-data-range area metric + overlap / ladder tier / sampled-source completeness / provenance completeness / envelope-hit coverage); **independence = recorded enum {self, oracle, external}** per axis, carried as `(level, independence, evidence_handle)`. §3.2: **`THETA_FIND = 2`** named with rationale (first all-quantified grade; 3 empties finds outside anchors, 1 admits calibrated-unvalidated paths). |
| 2026-08-13 | 0.2 | Consistency sweep: OFFL-2's non-Gaussianity diagnostic + dominant-nuclide truncation justification wired into the input-pedigree axis (§3.1, §4) — OFFL-2 already declared COUP-6 the consumer. |
| 2026-07-21 | 0.1 | Initial draft. Seven-axis PCMM maturity vector (META-1 §4.1), each axis fed by one distinct V&V product (no averaging), headline = weakest-link min, graded by rigor + independence (NASA-STD-7009 discipline). Finds-vs-leads as a computed triple-guard predicate with a permanent rung-(ii)→lead override. Assumption register references VAL-1's PIRT. Consumes VAL-1/COUP-5 outputs; runs no V&V machinery. |
