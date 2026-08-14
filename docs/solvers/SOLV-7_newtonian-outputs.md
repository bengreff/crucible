# SOLV-7 — Newtonian Outputs

| Field | Value |
|---|---|
| **ID** | SOLV-7 |
| **Family** | SOLV (Runtime unified-grid operators) |
| **Status** | Reviewed (2026-08-14) |
| **Depends on** | FND-2, SOLV-1; FND-1, VAL-2 |
| **Version** | 0.3 (2026-08-14 review fix: N11) |

---

## 0. Purpose

SOLV-7 turns the field solution into the **reported engine performance**: **thrust, specific impulse (Isp),
characteristic velocity (c\*), thrust coefficient (C_F), torque, and momentum flux**, computed by
**integrating the conserved fluxes over exit/reference planes** of the grid. It is the operator that produces
the numbers the RL10 anchor is validated against (S1). It is deliberately thin — a diagnostic integration over
`U`, no new physics.

Crucially, SOLV-7 is where the design decision of COUP-7 §3.2 lands: **chamber pressure and thrust are
emergent** — SOLV-7 *reads out* what the reacting flow produced, it never imposes an operating point.

Read after SOLV-1 (the field it integrates) and COUP-7 (the boundary flow that drives the reaction).

## 1. Scope & razor ruling
**Owns:** the exit-plane / reference-plane flux integrals and the standard performance identities (§3);
the reported **performance object** (§3.4). **Defers:** the field solution itself → **SOLV-1**; trajectory
propagation → **out of scope** (VISION_SCOPE §4.2 — thrust/Isp/mass curves are *exported*, not propagated);
the anchor targets and pass criterion → **VAL-2**; boundary inflow (ṁ, MR) → **COUP-7**.

**Razor ruling:** pure diagnostic — no energy release flows through SOLV-7, so it simulates nothing; it is the
Newtonian bookkeeping the fidelity doctrine explicitly keeps (net thrust vector, torque, momentum flux —
VISION_SCOPE §4.2 "General Newtonian bookkeeping"). It introduces **no model-form uncertainty**; it only
propagates what SOLV-1 and the boundary objects already carry.

## 2. Interfaces & contracts
| Interface | Consumed by | Contract |
|---|---|---|
| **Plane-flux integral** (§3.1) | reporting, VAL-2 | momentum + pressure flux integrated over a declared exit/reference plane of the grid |
| **Performance object** (§3.4) | FND-6 (results bundle), COUP-6 | {thrust, Isp, c\*, C_F, v_e, torque, momentum flux}, each an ensemble quantity (p-box), with the exit state it was derived from |
| **Emergent-quantity read** (§3.2) | — | reads `p_c` (from the field), never writes it |
| **Pulsed-impulse intake** *(deferred — W4, with SOLV-5)* | reporting, COUP-4 | consumes SOLV-5's per-event return {impulse, energy partition, wall loading}; folds event sequences into time-averaged thrust/Isp alongside the plane integrals |

**Invariant:** every reported number is a **plane integral of the conserved `U`**, so it inherits the solver's
conservation and determinism; it is an **ensemble** (never a bare scalar, S3); and it is **derived, never
imposed** — SOLV-7 asserts no operating point.

## 3. Method & governing definitions

### 3.1 Exit-plane integration
Thrust is the momentum-plus-pressure flux over the exit plane (vacuum ⇒ p_a = 0):
- **(SOLV-7.1)** `F = ∮_exit (ρ u_x² + (p − p_a)) dA` — recovers both the momentum term `ṁ v_e` and the
  pressure term `(p_e − p_a) A_e` directly from the field (no separate one-dimensional assumption).
- Torque and net momentum flux are the corresponding moment/vector integrals over the same plane set.

### 3.2 Emergent chamber pressure *(COUP-7 §3.2 realized; convention pinned, N11)*
`p_c` is **not** an input. The boundary turbopump injects ṁ at the injector port; the reacting flow combusts
and chokes at the throat, and `p_c` is read from the chamber field. SOLV-7 **reports** this emergent `p_c`
alongside thrust/Isp — it is a *predicted* quantity, which is what makes the RL10 comparison a genuine test
rather than a tautology.

**The p_c convention (N11):** `p_c` ≡ the **area-averaged stagnation pressure at the injector-end reference
plane** (a config-declared plane just downstream of the injector face), computed per cell from the local
static state + Mach via the field's own EOS and area-averaged over the plane. **Why this one:** it is the
convention of the RL10 data-of-record — TM-107318's cycle-table chamber-pressure station and the JANNAF c\*
convention are both **injector-end stagnation** — so S1's comparison is convention-matched by construction.
The choice is not cosmetic: at chamber Mach ~0.2–0.3 the **static-vs-stagnation difference (and the
plane-location drift toward the converging section) is ~1–2%** — the size of the **entire S1 target** — so
an unpinned convention would silently spend the whole error budget on a definition mismatch.
[META-3: `jannaf-pc-convention`, `rl10-tm107318`]

**`p_c = ṁ·c*/A_t` is demoted to a consistency check, not a definition** (it is circular against SOLV-7.2 —
each would define the other): the harness verifies the reported `p_c`, `ṁ`, `c*`, `A_t` close this identity
to within the plane-averaging + truncation tolerance, flagging a convention or integration bug if not.

### 3.3 Performance identities *(the validation-friendly split)*
[META-3: `cstar-cf-defs`]
- **(SOLV-7.2)** c\* = `p_c · A_t / ṁ` — combustion/chamber performance (upstream of the throat), with `p_c`
  per the §3.2 convention (injector-end stagnation — the JANNAF c\* convention).
- **(SOLV-7.3)** C_F = `F / (p_c · A_t)` — nozzle performance (downstream of the throat).
- **(SOLV-7.4)** Isp = `C_F · c\* / g₀` = `v_e / g₀`, with `v_e` the effective exhaust velocity and g₀ =
  9.80665 (META-2 §2.1). Reported in **seconds and as `v_e` (m/s)**.

The identity `v_e = c\*·C_F` **factors Isp into a combustion term (c\*) and a nozzle term (C_F)**, which VAL-2
validates **separately** against RL10's published c\* (7824 in/s, η_c\* 0.989) and C_F — so a discrepancy is
localized to combustion vs expansion rather than smeared into one Isp number.

### 3.4 The performance object
A single reported object carries {F, Isp, v_e, c\*, C_F, torque, momentum flux, emergent `p_c`, MR}, each a
`ResultDistribution`/p-box (FND-1 §3.4c) over the UQ ensemble, plus the exit-plane state it was integrated
from and a `ProvenanceRef`. Serialized by FND-6; pedigree-scored by COUP-6.

### 3.5 Determinism
Plane integrals use the grid's canonical Morton-ordered traversal + fixed-order reductions (FND-2 §3.7), so
every reported number is bit-reproducible at any thread count.

## 4. Coupling relationships
- **SOLV-1** provides the field `U`; SOLV-7 integrates its fluxes over declared planes — no back-coupling
  (diagnostic only).
- **COUP-7** sets the inflow (ṁ, MR); SOLV-7 reads the emergent `p_c`/thrust the reaction produced.
- **COUP-5** turns each per-member integral into the p-box; **FND-6** serializes the performance object;
  **VAL-2** compares it to the RL10 anchor under the overlap-band criterion (Ben, 2026-07-21 — no hard cutoff).
- **COUP-6** scores pedigree; the c\*/C_F split lets the pedigree note which factor (combustion vs nozzle)
  limits confidence.
- **SOLV-5** *(deferred, W4)*: per-event impulse returns fold into time-averaged thrust/Isp for pulsed
  concepts (SOLV-5 names SOLV-7 as this consumer).

## 5. Uncertainty & validity
SOLV-7 **originates no uncertainty** — it propagates SOLV-1's discretization band and the boundary objects'
bands (esp. the turbopump band, COUP-7, and the wall-function band, SOLV-1 §3.5, that drive the `p_c`
spread) into the reported p-box. Validity
is a statement about the exit-plane resolution (a plane thinner-resolved than the exit gradient is flagged,
FND-2 §5). Validation-ladder status: infrastructure (rung i); its *numbers* are validated at the system tier
by VAL-2 (RL10).

## 6. Validation plan
1. **Analytic nozzle:** for an ideal quasi-1-D nozzle, (SOLV-7.1)–(SOLV-7.4) reproduce the closed-form
   C_F(γ, p_e/p_c, ε) and c\* within solver truncation error.
2. **Identity consistency:** `v_e = c\*·C_F` holds to round-off on any run (guards the plane integrals);
   `p_c·A_t ≈ ṁ·c*` closes within the plane-averaging + truncation tolerance (the §3.2 consistency check —
   a check, never the definition).
3. **RL10 system anchor:** integrated F, Isp, c\*, C_F and **emergent `p_c`** overlap the RL10 reference
   p-box (440–446 s Isp, 73.4 kN, 475 psia, c\* 7824 in/s), reported with the Ferson area metric — **blind**
   (open-mode turbopump) and **calibrated** (closed-mode), per VAL-2. [META-3: `rl10-cycle-data`,
   `cstar-cf-defs`]
4. **Determinism:** identical reported numbers at 1 vs N threads.

## 7. References
META-3 keys: `cstar-cf-defs`, `rl10-cycle-data`, `rl10-tm107318`, `jannaf-pc-convention` *(new)*. Depends on
FND-2 (grid/`U`, traversal), SOLV-1 (field), COUP-7 (inflow), FND-1 (p-box), VAL-2 (anchor + overlap
criterion).

*(No open questions.)*

## 8. Change log
| Date | Version | Change |
|---|---|---|
| 2026-08-14 | 0.3 | Review fix (N11). `p_c` convention pinned: **area-averaged stagnation pressure at the declared injector-end reference plane**, matching the RL10 data-of-record (TM-107318 cycle-table station) and the JANNAF c\* convention — the ~1–2% static-vs-stagnation/plane difference is the size of the whole S1 target, stated. `p_c = ṁ·c*/A_t` demoted from definition to a harness consistency check (§3.2, §3.3, §6.2). |
| 2026-08-13 | 0.2 | Consistency sweep: added the deferred pulsed-impulse intake (SOLV-5's per-event {impulse, partition, wall loading} → time-averaged thrust/Isp, W4) — SOLV-5 already named SOLV-7 as this consumer. |
| 2026-07-21 | 0.1 | Initial draft. Exit-plane flux integration for thrust/torque/momentum; the c\*/C_F/Isp identities with the `v_e = c\*·C_F` split enabling separate combustion/nozzle validation; **emergent chamber pressure** (read, never imposed — realizes COUP-7 §3.2); the ensemble performance object (p-box); determinism via fixed-order plane reductions; RL10 validated under the overlap-band criterion (blind + calibrated), no hard cutoff (Ben 2026-07-21). |
