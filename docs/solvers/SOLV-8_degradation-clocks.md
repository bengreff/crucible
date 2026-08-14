# SOLV-8 — Degradation Clocks (Stage 2)

| Field | Value |
|---|---|
| **ID** | SOLV-8 |
| **Family** | SOLV (Runtime unified-grid operators) |
| **Status** | Reviewed (2026-08-14) |
| **Depends on** | FND-2 (recession), FND-7 (degraded allowables), SOLV-6 (margin re-check), COUP-4 (Stage-2/verdict), OFFL-1/OFFL-3 |
| **Version** | 0.1 |

---

## 0. Purpose

SOLV-8 owns the **slow Stage-2 degradation clocks** — ablation/recession, hot-H₂ corrosion, fuel burnup,
neutron fluence/DPA, and decay heat — marched with large deterministic timesteps after Stage-1 passes,
re-checking the Stage-1 margins as geometry and material degrade. The **first essential margin crossed** →
operational lifetime + limiting mechanism + confidence interval (VISION_SCOPE §7.6). It computes **rates** and
drives the owners (FND-2 geometry, SOLV-1 thermal, SOLV-6 structure); it carries no distribution arithmetic
(UQ is the outer loop).

Read after COUP-4 (Stage-2 control flow), SOLV-6 (the re-check hook), and FND-7 (degraded allowables).

## 1. Scope & razor ruling
**Owns:** the degradation-rate **closures** and their competition → lifetime + limiting mechanism + CI.
**Defers:** geometry mutation → **FND-2**; the thermal solve → **SOLV-1**; the structural formulas →
**SOLV-6**; the control flow/verdict → **COUP-4/COUP-3**; degraded allowables + damage-conversion (E_d, ξ) →
**FND-7**; B′ tables → **OFFL-3**; depletion Δρ(B) + damage-energy kernel → **OFFL-1**.

**Razor ruling (Failure-Mode Razor):** corrosion/DPA/burnup are **constraint checks against literature limits
+ bands** (uncertainty in the outer p-box) — not simulated crack/microstructure physics.

## 2. Interfaces & contracts
| Interface | Consumed by | Contract |
|---|---|---|
| **Recession rate** ṡ | FND-2 | `dκ/dt = −ṡ·A/V` monotone solid-fraction decrease + conservative inward spill; area from **PLIC** |
| **Degraded state** | SOLV-6 (re-check), FND-7 | accumulated dpa/burnup/corrosion → `allowable(T, degradation-state)` re-evaluation |
| **Lifetime verdict** | COUP-4 | `{lifetime, limiting mechanism, CI}` (min over clocks) or `lifetime > horizon` |
| **Clock coefficients** | COUP-5 | each rate coefficient is an `UncertainInput` → outer-loop lifetime p-box |

**Invariant:** clocks march **slow, large-Δt, fixed reduction order**; the limiting mechanism is the **MIN**
crossing over all clocks + re-checked margins; decay heat is **not** a monotone clock (§3.5); no cascading
failure (first crossing ends the stage, COUP-4).

## 3. Method

### 3.1 Ablation / recession
Enthalpy/VOF **monotone solid-fraction κ decrease on the fixed grid**; ṡ from the **equilibrium B′ surface
energy balance** (`bprime`/OFFL-3), with **Extended-B′** (`ablation-recession`) as the *assumption ledger*
(flags where equilibrium/finite-rate-pyrolysis fails); flux area from **PLIC** (`vof-plic`), never
marching-cubes (~8% biased, `mc33-dc`); `dκ/dt=−ṡ·A/V` conservative inward spill → FND-2. Validated vs analytic
Stefan + **Lachaud Ablation Test-Case #2** (this closes FND-2-Q1). [META-3: `lachaud-ablation-tc2`]

### 3.2 Hot-H₂ corrosion
Two physical regimes, blended by a **continuous switching function of the coating-integrity state** (not a
boolean `if(breached)` — Rule 12): parabolic diffusion-limited `k_p(T)√t` (Arrhenius; protective layer intact)
→ linear reaction-limited substrate attack (layer consumed). Coefficients are literature closures with
citation + validity envelope + band (`carbide-h2-corrosion`, `ntp-fuels`). Breach-onset uncertainty is
**epistemic** (an `UncertainInput` in the outer p-box); each member runs a fixed threshold deterministically →
the spread gives the lifetime CI (bimodal-safe — reported as a distribution, never collapsed).

### 3.3 Burnup — tiered, often self-zeroing
Time-integrated %FIMA; reactivity swing **tiered**: analytic `Δρ≈−αB` default → full **OpenMC depletion**
(OFFL-1) for long-burn NEP/fission-fragment → **self-zeroing** for short-burn NTP (minutes → `burnup-limited
lifetime ≫ horizon`). Re-checks the **criticality-control** margin (drum authority exhausted → COUP-4 nuclear
halt).

### 3.4 Neutron fluence & DPA
OpenMC `damage-energy` **MT-444 is RAW damage energy, not dpa** — unit-scaled in the loader/OFFL-1
(`nrt-arc-dpa`). SOLV-8 integrates fluence and applies **both NRT and arc-dpa** (ξ per material; **W-specific,
never Fe-scaled**); the **NRT↔arc spread is the model-form band**. Checks vs the FND-7 dpa limit (~20–100 dpa)
+ **He-appm co-limit**; accumulated dpa indexes `allowable(T,dpa)` for the SOLV-6 re-check (embrittlement).

### 3.5 Decay heat — a shutdown transient, not a clock
Decay heat **decreases** with time, so a "crossing time" is a category error; the hazard is a **peak
temperature after shutdown** with reduced/lost active cooling. Modeled as an **ANS-5.1-2014 afterheat P(t)
source into SOLV-1 conduction, re-checking the COUP-4 melt margin** — a separate pass/fail transient
(triggered by a scram in the operating profile or a worst-case end-of-life scram), **not** part of the MIN
competition. **Fission legs only** (fusion/fission-fragment afterheat is activation decay via SOLV-2/OFFL-1).
LWR-derivation → **conservative over-prediction band** declared. [META-3: `ans51-decayheat`]

### 3.6 Multi-clock competition
Monotone clocks marched in parallel with large deterministic timesteps; after each step, **re-evaluate the
SOLV-6 margins + COUP-4 halts** against the degraded state. Limiting mechanism = **MIN** over {clock-limit
crossings, re-checked margin < 1, reactivity-margin loss}; the UQ ensemble turns the crossing time into a
**CI**; no crossing in the horizon → `lifetime > horizon`.

## 4. Coupling relationships
- **FND-2** applies the recession fraction update (owns geometry); **SOLV-1** runs the decay-heat transient;
  **SOLV-6** re-checks margins on the degraded state; **FND-7** supplies `allowable(T,dpa,burnup)` + E_d/ξ/He;
  **OFFL-1** supplies depletion Δρ(B) + damage-energy; **OFFL-3** supplies B′; **COUP-4** consumes the lifetime
  verdict; **COUP-5** samples the clock coefficients.

## 5. Uncertainty & validity
Every closure coefficient (B′, k_p, breach onset, ξ, α/Δρ) is an `UncertainInput` → outer LHC → lifetime
p-box/CI (bimodal-safe); NRT↔arc = declared model-form band; corrosion/DPA/burnup closures are **rung (ii)**
(permanently pedigree-flagged). Validity envelopes on all rate closures; out-of-envelope flags/refuses.

## 6. Validation plan
1. **Recession** vs analytic Stefan + **Lachaud Ablation Test-Case #2** (closes FND-2-Q1).
2. **Corrosion** vs published `k_p(T)√t` data.
3. **DPA** vs `nrt-arc-dpa` reference for Fe **and** W (guards against Fe-scaling W).
4. **Decay heat** P(t) vs ANS-5.1 tabulated curve.
5. **Determinism:** identical limiting mechanism/CI at 1 vs N threads.

## 7. References
META-3 keys: `bprime`, `ablation-recession`, `lachaud-ablation-tc2`, `vof-plic`, `mc33-dc`,
`carbide-h2-corrosion`, `ntp-fuels`, `nrt-arc-dpa`, `ans51-decayheat`. Depends on FND-2, FND-7, SOLV-1, SOLV-6,
COUP-4, OFFL-1, OFFL-3.

*(No open questions — decay heat is a shutdown-transient re-check (not a monotone clock); coating-breach is a
continuous state switch with epistemic onset in the outer p-box; burnup tiered/self-zeroing; W never
Fe-scaled. Resolved 2026-07-21. FND-2-Q1 closes at test §6.1.)*

## 8. Change log
| Date | Version | Change |
|---|---|---|
| 2026-07-21 | 0.1 | Initial draft. Stage-2 slow clocks: enthalpy/VOF monotone recession (PLIC area, Extended-B' ledger); two-regime hot-H₂ corrosion (continuous coating-state switch, epistemic breach onset); tiered/self-zeroing burnup; NRT+arc-dpa (W-specific, MT-444 unit-scaled); decay heat as a post-shutdown melt-margin transient (not a clock; fission-only; conservative LWR band). Multi-clock MIN-crossing competition → lifetime + limiting mechanism + CI; re-checks SOLV-6 margins on degraded state via FND-7 allowable(T,dpa). |
