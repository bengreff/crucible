# SOLV-6 — Structural Margins

| Field | Value |
|---|---|
| **ID** | SOLV-6 |
| **Family** | SOLV (Runtime unified-grid operators) |
| **Status** | Reviewed (2026-08-14) |
| **Depends on** | FND-2, FND-3/FND-4 (structural annotation), FND-7 (allowables), COUP-4 (halt input); FND-1 |
| **Version** | 0.2 (2026-08-14 review fixes: N8, N9, N10) |

---

## 0. Purpose

SOLV-6 computes **analytic quasi-static structural margins** — thin-shell hoop, longitudinal, thermal, and
pressure/burst stresses — as **pass/fail + margin** at load-bearing cells, feeding the `burst margin < 1` halt
(COUP-4). It is a **constraint check, not a structural simulation** (Failure-Mode Razor): no FEM, no fatigue,
no fracture in v1 (voxel-hex linear-elastic FEM is a stretch upgrade, VISION_SCOPE §5.1).

Read after FND-2 (the wall p/T fields it reads) and FND-7 (the material allowables).

## 1. Scope & razor ruling
**Owns:** the analytic stress formulas (Roark) and the margin/pass-fail output. **Defers:** the pressure/
temperature fields → **SOLV-1/conduction**; material strength allowables (static handbook limits) → **FND-7**;
the halt decision → **COUP-4**; slow degradation of margins → **SOLV-8** (stage 2 re-checks these formulas as
geometry/material degrade).

**Razor ruling:** Failure-Mode Razor — structural failure is mechanical engineering, so SOLV-6 checks margins
against literature allowables and reports pass/fail; it does not simulate crack growth or vibration.

## 2. Interfaces & contracts
| Interface | Consumed by | Contract |
|---|---|---|
| **Margin field** | COUP-4 | per structural component (and its member cells): {σ_eq, allowable, margin, pass/fail} |
| **Structural annotation** (§3.1) | FND-3/FND-4 | config-time declaration of load-bearing components as **shell primitives** — R(z), t(z) derived from the CSG geometry (FND-3), cell membership mapped at voxelization; the geometry/config machinery is FND-3/FND-4's, SOLV-6 only consumes the annotation |
| **Allowables query** | FND-7 | `allowable(T, degradation-state)` — A/B-basis yield/UTS as a function of T **and** accumulated dpa/burnup (irradiation knockdown, FND-7); the Stage-1 call passes dpa=0, the Stage-2 re-check passes SOLV-8's accumulated state. One allowable law over the degradation state, no branch (Rule 12) |

**Invariant:** a margin is a **derived diagnostic** of the current p/T field + geometry; it writes no physics
state; a `margin < 1` is a halt input, not a silent flag.

## 3. Method
[META-3: `roark`, `nasa-std-5012`, `ntp-fuels`, `mmpds`]

### 3.1 Structural annotation (config-time, N9)
Per-cell state alone cannot yield D, t, or "load-bearing" — those are **component** properties. At config
time, each load-bearing component is declared a **shell primitive**: a `(R(z), t(z))` profile **derived from
the same CSG geometry FND-3 voxelizes** (the revolved-profile primitives carry it analytically; no second
geometry description), with its member cells mapped at voxelization and recorded in the manifest (FND-4).
The annotation supplies each formula's operands: local radius `R`, thickness `t`, internal `p` (from the
bounding fluid cells), and **ΔT taken between the paired inner/outer surfaces of the same component**
(paired at equal z along the shell normal — never between cells of different components). The
geometry/config machinery is FND-3/FND-4's; SOLV-6 consumes the annotation. A load-bearing component
*without* an annotation is a config error (fail loud), not a guessed shell.

### 3.2 Stress state & margin (N8)
Per-direction superposition of pressure and thermal terms — never a scalar sum of orthogonal components:
- **(SOLV-6.1)** Hoop (thin shell, 2R/t > 20): `σ_hoop = p·R/t`.
- **(SOLV-6.2)** Longitudinal: `σ_long = p·R/(2t)`.
- **(SOLV-6.3)** Thermal, **biaxial constrained-liner form** (equal in hoop and longitudinal — a flat/shell
  wall constrained in both in-plane directions): `σ_therm = E·α·ΔT/(1−ν)`; for a **through-wall linear
  gradient** the surface stress is `σ_therm = E·α·ΔT/(2(1−ν))` (Roark). The annotation declares which case
  applies per component (uniform-constrained vs linear-gradient); the uniaxial `E·α·ΔT` is **retired** — it
  under-predicts the constrained case by the factor (1−ν) and mis-states the gradient case entirely.
  For a cooled liner σ_therm often **dominates** the pressure stress.
- **(SOLV-6.4)** Per-direction totals: `σ_θ = σ_hoop + σ_therm`, `σ_z = σ_long + σ_therm` (signs per load
  sense; radial ≈ 0, thin shell) — then the margin is taken against a **stated equivalent-stress criterion,
  von Mises** (plane stress): `σ_eq = √(σ_θ² − σ_θ·σ_z + σ_z²)`.
- **(SOLV-6.5)** Margin = `allowable / (FS · σ_eq)` against **yield** (allowable = A/B-basis yield,
  `FS_YIELD`) and **burst** (allowable = UTS, `FS_ULT`). Margin < 1 on any component → COUP-4 halt.
- **Safety factors (N10):** single sourced named constants — **`FS_YIELD` = 1.1**, **`FS_ULT` = 1.4**
  (NASA-STD-5012-class factors for liquid-fueled propulsion engine structures; `nasa-std-5012`) — replacing
  the earlier unsourced ~1.5 / 2.5–4 ranges. They appear once, in config defaults, with this citation.
- Material allowables (Inconel 625, SS, Cu alloys) are **T-dependent** from FND-7's static-limit table
  (A/B-basis), interrogated at the component's surface temperature.

The dominant uncertainty is **ΔT**, which rides the same **±20–30% wall-function band** as the wall heat
flux (SOLV-1 §3.5, the one local wall-heat law; Bartz is a VAL-2 oracle only) — so the structural margin's
band and the performance band share a source, and Sobol will show it.

## 4. Coupling relationships
- **SOLV-1/conduction** provide wall p and T; **FND-3/FND-4** provide the structural annotation (§3.1);
  **FND-7** provides allowables; **COUP-4** consumes the margin as a halt condition; **SOLV-8** (stage 2)
  re-evaluates (SOLV-6.1)–(SOLV-6.5) as recession thins `t` and burnup/DPA degrade the material, producing
  the lifetime-limiting margin crossing.

## 5. Uncertainty & validity
Introduces **model-form error of the analytic idealization** (thin-shell, elastic, quasi-static, von Mises
equivalent stress) plus the allowable and ΔT bands (FND-1). Validity: `2R/t > 20`, elastic, quasi-static;
outside these it flags under-idealized rather than silently applying the formula. A component where the
thin-shell assumption fails is reported, not smeared. Ladder rung: analytic (i), with the *allowables*
validated at their handbook source.

## 6. Validation plan
1. **Analytic pressure vessel:** (SOLV-6.1)/(SOLV-6.2) reproduce the closed-form thin-shell stresses.
2. **Thermal stress:** (SOLV-6.3) matches the Roark biaxial constrained-plate and through-wall-gradient
   closed forms — including the (1−ν) factor vs the retired uniaxial form.
3. **Equivalent stress:** (SOLV-6.4) matches hand-computed von Mises for a combined hoop+thermal case.
4. **Annotation:** a revolved-profile CSG chamber yields the analytic R(z), t(z), paired inner/outer
   surfaces, and correct cell membership (§3.1).
5. **Margin/halt:** an over-pressurized/over-heated wall yields margin < 1 at the right component (COUP-4).
6. **Band propagation:** the ΔT (wall-function) band produces a margin distribution; Sobol attributes it.
7. **Determinism:** identical margins at 1 vs N threads.

## 7. References
META-3 keys: `roark` (biaxial thermal-stress + combined-stress forms), `nasa-std-5012` *(new)*, `mmpds`,
`ntp-fuels`, `wall-function-heat` (ΔT band). Depends on FND-2 (p/T fields), FND-3/FND-4 (structural
annotation), FND-7 (allowables), COUP-4 (halt), SOLV-8 (stage-2 re-check).

*(No open questions.)*

## 8. Change log
| Date | Version | Change |
|---|---|---|
| 2026-08-14 | 0.2 | Review fixes (N8, N9, N10). **N8:** uniaxial `E·α·ΔT` retired — biaxial constrained-liner `EαΔT/(1−ν)` + through-wall-gradient `EαΔT/(2(1−ν))` (Roark); per-direction superposition (σ_θ, σ_z); margin against the **stated von Mises equivalent-stress criterion** (§3.2). **N9:** config-time structural annotation — load-bearing components as shell primitives R(z), t(z) derived from the FND-3 CSG, member-cell mapping, ΔT between paired inner/outer surfaces of the same component (§3.1, §2). **N10:** safety factors as single sourced named constants `FS_YIELD` = 1.1 / `FS_ULT` = 1.4 (`nasa-std-5012`), replacing unsourced ranges. ΔT-band source updated Bartz → the one wall-function law (SOLV-1 §3.5, D-C). |
| 2026-07-21 | 0.1 | Initial draft. Analytic quasi-static hoop/longitudinal/thermal/burst margins (Roark) → pass/fail + margin as a COUP-4 halt input; T-dependent A/B-basis allowables from FND-7; ΔT (Bartz) as the dominant band, shared with the performance spread; validity gated to thin-shell/elastic/quasi-static; stage-2 re-check hook. FEM/fatigue/fracture out of v1 (Failure-Mode Razor). |
