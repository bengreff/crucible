# COUP-4 — Two-Stage Execution & Halt/Failure Model

| Field | Value |
|---|---|
| **ID** | COUP-4 |
| **Family** | COUP (Coupling & orchestration) |
| **Status** | Reviewed (2026-08-14) |
| **Depends on** | COUP-3, COUP-2, SOLV-6, SOLV-8; COUP-8 (`halts[]`), COUP-5 (§3.2.1 P(WORKS)) |
| **Version** | 0.3 (2026-08-19: physical-march-only Stage 1; NEVER_IGNITED/FLAMEOUT halts — VISION_SCOPE v1.5) |

---

## 0. Purpose

COUP-4 owns the **control flow of a run**: the two-stage execution model (**Stage 1 FUNCTION** → does the
engine reach and hold its commanded operating profile; **Stage 2 LIFETIME** → how long until a slow-degradation
margin is crossed), the **closed set of halt conditions**, and the **verdict object**. The governing rule
(VISION_SCOPE §7.6): **first loss of an essential function halts the run immediately — no cascading-failure
modeling, ever.**

Read after COUP-3 (the step it wraps), COUP-2 (whose audit failure is a halt), and SOLV-6 (whose margins are
halt inputs).

## 1. Scope & razor ruling
**Owns:** the Stage-1/Stage-2 control flow; the halt-condition **enum** and the check schedule; the verdict
object; the Stage-1→Stage-2 handoff. **Defers:** the per-step advance → **COUP-3**; conservation-audit
detection → **COUP-2**; structural margins → **SOLV-6**; degradation clocks → **SOLV-8**; which halts an
operator can raise → its **COUP-8 `Manifest.halts[]`**.

**Razor ruling (Failure-Mode Razor):** mechanical/electrical failure modes are **constraint checks, not
simulations**; COUP-4 checks margins and halts — it does not model the failure's evolution. One loss ends the
run.

## 2. Interfaces & contracts
| Interface | Consumed by | Contract |
|---|---|---|
| **Stage controller** | the run entry point | run Stage 1 to a verdict; if WORKS, run Stage 2 to a lifetime estimate |
| **Halt check** (§3.2) | COUP-3 (every step) | evaluate the closed halt set; first trip → immediate halt + diagnosis |
| **Verdict object** (§3.3) | FND-6, COUP-6 | `WORKS` + performance (SOLV-7) **or** `DOESN'T WORK` (mechanism, location, time) |

**Invariant:** the run is a **pure function** of `{config, tables, seed}` up to the first halt; the verdict is
**deterministic** and carries a diagnosis, never a bare failure.

## 3. Method

### 3.1 Two stages
- **Stage 1 — FUNCTION.** A **physical transient march** to the commanded operating profile — always
  (VISION_SCOPE §7.6 v1.5: accelerated convergence deleted; start-up is part of what Stage 1 verifies).
  External-system timelines may be compressed as declared boundary-object schedules (COUP-7 §3.2.2). For
  the chemical slice this is the march from declared fill state through ignition to the RL10 operating
  point (COUP-7 flow in → emergent `p_c`/thrust out). Reaches and holds ⇒ `WORKS` + performance report;
  otherwise a halt (§3.2).

  **The WORKS criterion is three named constants (O12), recorded in the verdict object (§3.3):**
  - `EPS_WORKS[q]` — per-quantity **relative tolerance on the commanded profile** (default 0.02 — the S1
    reported-target scale; per-quantity config override). "Reaches" = every commanded quantity within
    `EPS_WORKS[q]` of its commanded value.
  - `T_DWELL` — the **dwell window**: WORKS requires *holding* every commanded quantity inside `EPS_WORKS[q]`
    for a contiguous physical span `T_DWELL` (default 20 chamber flow-through times; config-overridable),
    as the tail of the one physical march.
  - `T_S1_HORIZON` — the **Stage-1 horizon** (config-required). Horizon expiry without a
    completed dwell ⇒ the **distinct `FAILED_TO_REACH` halt** (§3.2) — never a silent timeout, never
    conflated with a physical mechanism halt.
- **Stage 2 — LIFETIME.** Runs **only after Stage 1 passes**; freezes the operating point and marches the slow
  degradation clocks (SOLV-8) with large deterministic timesteps, re-checking Stage-1 margins as geometry/
  materials degrade. First essential margin crossed ⇒ **lifetime estimate + limiting mechanism + CI**; nothing
  crossed within the horizon ⇒ `lifetime > horizon`.

### 3.2 The closed halt set
Checked every step (COUP-3); the enum a mechanism may raise is its `Manifest.halts[]` (COUP-8). Chemical-slice
members:
- **Melt / vaporization** of a structural material (spine state crosses phase limit at a load-bearing cell).
- **Burst margin < 1** (SOLV-6 hoop+thermal+pressure vs allowable).
- **Choking / starvation at a port** (mass flow cannot be sustained; e.g. the expander cycle cannot close —
  divergent iterates or an envelope-refused fixed point, COUP-3 §3.5 — no consistent operating point →
  **won't-bootstrap → `DOESN'T WORK`**; diagnosis = **physical**).
- **`NEVER_IGNITED`** *(v0.3, VISION_SCOPE v1.5)*: the commanded ignition sequence completed (igniter
  schedule exhausted, COUP-7 §3.3) and no self-sustaining burn front exists — the burn-progress field's
  reacting measure never exceeds its named threshold (SOLV-4 §3.6). Diagnosis = **physical**
  (`DOESN'T WORK (never ignited)`).
- **`FLAMEOUT`** *(v0.3)*: a previously established burn extinguishes before the dwell completes (the
  reacting measure collapses — quench/flammability physics, SOLV-4 §3.6). Diagnosis = **physical**, with
  location + time of the extinction front.
- **`COUPLING_RESIDUAL`** (O4): a fixed-sweep coupling solve failed its named residual-acceptance test —
  defined in **COUP-3 §3.5** (`EPS_EXPANDER_RESID` et al.). Diagnosis = **numerical** (a solver defect, never
  an engine verdict); **never conflated** with the physical won't-bootstrap above.
- **`FAILED_TO_REACH`** (O12): `T_S1_HORIZON` expired without a
  completed dwell — diagnosis names the quantities still outside `EPS_WORKS[q]` and by how much.
- **Conservation-audit failure** (COUP-2 beyond `TOL_AUDIT` — a bug, not a physical result).
- **Non-finite field** (NaN/Inf — FND-1 §3.8).
- *(Nuclear/plasma/pulsed halts — loss of criticality control, quench/extinction — added with those legs.
  **Tracked obligation (O13):** the PKE-excursion halt spec **and** its commanded-envelope source (the FND-4
  operating-profile grammar entry it checks against) must land **before the nuclear coding wave** — the
  milestone includes NTP; this deferral has a due date, not an open end. The nuclear mechanisms' COUP-8
  `halts[]` extension is noted there.)*

First trip → **halt immediately**, emit the verdict, stop. No subsequent-failure modeling.

### 3.3 The verdict object
Per member: `{ outcome ∈ {WORKS, DOESNT_WORK}, mechanism, location (cell/port), time, diagnosis ∈ {physical,
numerical}, criterion { EPS_WORKS[q], T_DWELL, T_S1_HORIZON }, performance? (SOLV-7 object if WORKS),
lifetime? (SOLV-8 estimate + limiting mechanism + CI if Stage 2 ran) }` — the WORKS-criterion constants are
**recorded in the verdict** (O12). Serialized by FND-6; pedigree-scored by COUP-6. Every field is
deterministic and diagnosis-carrying (mechanism+location+time on a `DOESN'T WORK`).

**Ensemble level (D-F):** the ensemble verdict object records the **per-member verdicts**, and **P(WORKS) is
a first-class reported result** — its estimator, epistemic interval bounds, and the conditional-on-WORKS rule
for performance functionals are owned by **COUP-5 §3.2.1** (cross-reference; COUP-4 only defines what a
member verdict is).

### 3.4 Determinism
The stage controller, halt-check order, and Stage-2 clock schedule are fixed; the verdict is a pure function of
the run inputs. A halt is a **structured, reported outcome** (not a `panic!`) — the one exception to the
"solvers return typed errors" rule (META-2 §4).

## 4. Coupling relationships
- **COUP-3** calls the halt check each step and runs the Stage-2 slow-clock march; **COUP-2** audit failure is
  a halt; **SOLV-6** margins and **SOLV-8** clocks are halt/lifetime inputs; **COUP-7** port starvation (failed
  expander closure) is a halt.
- **SOLV-7** supplies the WORKS performance report; **FND-6/COUP-6** consume the verdict + pedigree.

## 5. Uncertainty & validity
COUP-4 originates no physics uncertainty; the **verdict itself is an ensemble outcome** — over the UQ members
some may WORK and some DOESN'T WORK. **P(WORKS)** is reported as its own first-class functional with
epistemic bounds (estimator owned by COUP-5 §3.2.1), with the dominant limiting mechanism; performance is
**conditional-on-WORKS**, labeled; never collapsed to a single bit. Validation-ladder status: infrastructure,
CI-gated.

## 6. Validation plan
1. **Halt-on-melt / burst:** a deliberately over-driven case halts at the correct mechanism/location/time.
2. **Starvation:** an expander config that cannot close halts as `DOESN'T WORK (cycle won't bootstrap)`
   (physical); a planted non-contractive solve halts `COUPLING_RESIDUAL` (numerical) — the diagnoses never
   conflate (mirrors COUP-3 §6.7).
3. **No cascading:** exactly one halt fires; nothing runs after it.
4. **Stage handoff:** Stage 2 runs only if Stage 1 = WORKS; lifetime reported with limiting mechanism.
5. **Ensemble verdict:** a mixed ensemble carries per-member verdicts and reports P(WORKS) +
   dominant failure mechanism, not a scalar; performance outputs are labeled conditional-on-WORKS.
6. **Determinism:** identical verdict/timing at 1 vs N threads.
7. **WORKS criterion:** a fixture holding just inside / drifting just outside `EPS_WORKS[q]` across `T_DWELL`
   verdicts WORKS / not; horizon expiry trips `FAILED_TO_REACH` with the offending quantities named.

## 7. References
META-3 keys: `modelica-connector` (port starvation). Depends on COUP-3 (step), COUP-2 (audit halt), SOLV-6
(margins), SOLV-8 (clocks), COUP-8 (`halts[]`), SOLV-7 (performance), FND-6/COUP-6 (verdict/pedigree).

*(No open questions.)*

## 8. Change log
| Date | Version | Change |
|---|---|---|
| 2026-08-19 | 0.3 | **VISION_SCOPE v1.5 (Ben).** §3.1: Stage 1 = **physical march only** (pseudo-transient references removed; dwell = the tail of the one march; compressed external schedules per COUP-7 §3.2.2). §3.2: chemical halt set gains **`NEVER_IGNITED`** and **`FLAMEOUT`** (physical diagnoses, thresholds owned by SOLV-4 §3.6). |
| 2026-08-14 | 0.2 | **Review fix wave (O12, O4, O13, D-F).** §3.1: the WORKS criterion made concrete — named constants `EPS_WORKS[q]` (per-quantity relative tolerance on the commanded profile, default 0.02), `T_DWELL` (dwell window, default 20 flow-through times), `T_S1_HORIZON`; horizon expiry = the distinct `FAILED_TO_REACH` halt; constants recorded in the verdict (O12). §3.2: **`COUPLING_RESIDUAL`** halt variant added (diagnosis = numerical; acceptance test defined in COUP-3 §3.5) — never conflated with physical won't-bootstrap (O4); nuclear deferral upgraded to a **tracked obligation** — PKE-excursion halt spec + commanded-envelope source due before the nuclear coding wave, COUP-8 `halts[]` extension noted (O13). §3.3/§5: verdict object records **per-member verdicts** + criterion + physical/numerical diagnosis; **P(WORKS) is a first-class result** with estimator/bounds owned by COUP-5 §3.2.1; performance conditional-on-WORKS, labeled (D-F). §6 items 2, 5, 7 updated/added. |
| 2026-07-21 | 0.1 | Initial draft. Stage-1 FUNCTION / Stage-2 LIFETIME control flow; the closed chemical-slice halt set (melt/vaporization, burst<1, choking/starvation incl. failed expander closure, conservation-audit failure, non-finite); halt = immediate structured verdict (mechanism/location/time), no cascading failure; ensemble verdict as a WORKS-fraction + dominant mechanism; deterministic control flow. |
