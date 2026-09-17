# VAL-1 — Predictive-V&V & Validation-Ladder Framework

| Field | Value |
|---|---|
| **ID** | VAL-1 |
| **Family** | VAL (Validation & test) |
| **Status** | Reviewed (2026-08-14) |
| **Depends on** | META-1 §4.1, FND-1, COUP-6 |
| **Version** | 0.2 (2026-08-14 review fixes: N15, D-G) |

---

## 0. Purpose

VAL-1 owns the **predictive-validation framework** that makes S8 real: the **4-tier validation hierarchy**
(unit → benchmark → subsystem → system), **solution verification** (GCI), **predictive model-form UQ** (area
metric → extrapolated band → p-box), and how ladder maturity maps to the **PCMM weakest-link pedigree**
(COUP-6). It is the doctrine of META-1 §4.1 made operational.

Read after META-1 §4.1 (the doctrine), FND-1 (the p-box/uncertainty types), and before VAL-2/VAL-3 (which
instantiate anchors and tests against this framework).

## 1. Scope & razor ruling
**Owns:** the tier hierarchy; the GCI solution-verification recipe; the area-metric model-form recipe (incl.
the small-n fallback ladder, §3.3); the p-box construction; the PIRT assumption register; the maturity→PCMM
mapping; the **operational blind-mode rule** (§3.6). **Defers:** per-anchor data + per-anchor blind/calibrated
input lists → **VAL-2**; the CI/test harness that runs it → **VAL-3**; the pedigree **score** → **COUP-6**;
the uncertainty **representation** → **FND-1**.

**Razor ruling:** infrastructure — no physics. Its obligation is that every headline result carry an
**honest, non-collapsible** bound and a **weakest-link** maturity, per the false-confidence guard.

## 2. Interfaces & contracts
| Interface | Consumed by | Contract |
|---|---|---|
| **Tier assignment** | every SOLV/OFFL/FND module | each constituent law states its highest achieved ladder tier |
| **GCI band** | FND-6, COUP-6 | a numerical (solution-verification) interval on every headline QoI, with observed vs formal order |
| **Model-form band** | FND-6, COUP-6 | the area-metric 95% prediction interval at the application condition (widens with extrapolation) |
| **PIRT register** | COUP-6 | phenomena ranked importance × knowledge — the S8 assumption register |

**Invariant:** a headline result is a **p-box** (aleatory distribution + epistemic model-form/numerical
intervals), never a lone distribution; its pedigree is the **minimum** maturity axis, never a mean.

## 3. Method

### 3.1 The 4-tier hierarchy
**Unit** (a single constitutive law/operator vs analytic or fundamental data) → **benchmark** (a published
benchmark or reference-code cross-check) → **subsystem** (a coupled subset) → **system** (a full engine
anchor). Validate **bottom-up, predict top-down**: a system-level extrapolated prediction is trusted only
because every constituent is validated **and** the composition introduces no un-validated physics (argued,
and checked at the system tier wherever data exist). This is *why* Rule 12 is load-bearing — **first-principles
laws extrapolate defensibly; calibrated curve-fits do not**. [META-3: `vv-oberkampf-roy`]

### 3.2 Solution verification — GCI
Every headline QoI carries a **Grid Convergence Index** from ≥3 systematically refined grids (r ≥ 1.3):
compute the **observed order** `p` and **check it against the formal order**; GCI = `Fs·|ε|/(rᵖ−1)` with
**Fs = 1.25** (3 grids, in the asymptotic range) or **3** otherwise. A result outside the asymptotic range is
flagged **unverified**, never reported with a false bound. [META-3: `gci-roache`, `asme-vv20`]

### 3.3 Model-form UQ — area metric → extrapolated band
Where validation data exist, measure the **Ferson area metric** (area between predicted and empirical CDFs,
in the QoI's units) with **u-pooling** across heterogeneous conditions; **regress it against a physical
coordinate** and take the **95% prediction interval at the application condition** as the model-form band —
which **widens with distance from data** by construction. A *calibrated* discrepancy term is **never** carried
into extrapolation (it feigns confidence far from data). [META-3: `area-metric`]

**Small-n fallback ladder (N15).** A regression + 95% PI is undefined at milestone anchor counts. The recipe
by anchor count `n` (per QoI, along the declared physical coordinate):
- **n = 1–2:** no regression — at each anchor compute an **ASME V&V20-style validation uncertainty**
  `u_val = √(u_num² + u_input² + u_D²)` and take the model-form band at the anchor as `|E| + u_val`
  (E = prediction − data). Extrapolating beyond the anchor's regime multiplies this band by a **declared,
  sourced extrapolation multiplier** — default **×3**, the Roache conservative-factor precedent for
  out-of-asymptotic/underdetermined evidence (`gci-roache`; `asme-vv20`) — never presented as a prediction
  interval it isn't.
- **n ≥ 3:** the regression recipe above (area metric vs physical coordinate → 95% PI).
**Crossover is exactly n = 3**; which branch produced a band is **recorded in the pedigree validation axis**
(COUP-6) — a ×3-multiplied point estimate and a regressed PI are different maturities and must grade
differently. [META-3: `asme-vv20`, `gci-roache`]

### 3.4 The single honest bound — a p-box
Aleatory inputs → distributions; epistemic (model-form, numerical) → **intervals**; propagate nested
(epistemic outer, aleatory inner) into an **interval-valued CDF**. A ≤10% band and an order-of-magnitude band
are **equally successful** if each is the tightest the stated minimal assumptions allow and is **stated, not
hidden** (S8). Collapsing to one probability invites the false-confidence theorem. [META-3: `pbox`,
`false-confidence`]

### 3.5 PIRT & the maturity→PCMM mapping
Each regime opens with a **PIRT** (phenomena ranked importance × state-of-knowledge) — the S8 assumption
register and the driver of *where* the p-box must widen. Maturity is reported as the **PCMM vector**
(physics-model fidelity, code verification, solution verification, validation, UQ, input pedigree, robustness),
each graded by rigor + independence, headlined by the **weakest-link minimum** (averaging incommensurable axes
is the chain-strength fallacy). Modules stuck at rung (ii) (plasma, antimatter) are **permanently flagged**.
[META-3: `pirt`, `pcmm`]

### 3.6 The operational blind-mode rule *(D-G; VISION_SCOPE §9 v1.4, Ben 2026-08-14)*
A **blind** validation run consumes **no quantity measured on the engine under test — period**. Its config may
reference exactly what the designer of a never-tested engine would have:
- the **design specification**: geometry-of-record, materials, commanded operating profile;
- **universal physics and closures**: the first-principles laws and their offline-calibrated universal
  closures (wall function, subgrid mixing, EOS/opacity spine) with their declared bands;
- **technology-class data measured on *other* hardware**, each with a cited band (e.g. an injector-family
  η_c\* prior band, a pump-class efficiency envelope) — "we assume we know exactly what we know for a fusion
  engine."
**Anchor-derived values** (a fitted η_c\*, the anchor's own measured pump map) bind only in **calibrated**
runs, and every reported score is **labeled blind or calibrated**. Comparison against the anchor's measured
values happens only *after* the blind run. Enforcement is mechanical, not disciplinary: blind configs carry a
manifest flag and their input provenance (META-1 Principle 7) is checked against the anchor's source keys —
an anchor-measured quantity in a blind config is a load error. VAL-2 instantiates this per anchor.

## 4. Coupling relationships
- **VAL-2** supplies the anchor data each tier validates against; **VAL-3** runs the GCI/area-metric/MMS
  machinery as CI; **COUP-6** consumes the tier + PCMM vector into the pedigree; **FND-6** serializes the GCI +
  model-form intervals into the result p-box; **FND-1** provides the p-box type.

## 5. Uncertainty & validity
VAL-1 originates no physics uncertainty — it is the framework others are graded in. Its correctness is that
the recipes are applied faithfully (GCI observed-order check, non-collapsed p-box, weakest-link pedigree),
CI-enforced (VAL-3). Validation-ladder status: N/A (it *is* the ladder).

## 6. Validation plan
1. **GCI self-test:** on a manufactured solution of known order, the recipe recovers the formal order and a
   correct GCI. [META-3: `mms`, `gci-roache`]
2. **Area-metric self-test:** on synthetic data with known discrepancy, the metric + 95% PI reproduce it and
   widen correctly under extrapolation.
3. **P-box non-collapse:** a result with epistemic inputs yields an interval-valued CDF, never a scalar.
4. **Weakest-link:** the PCMM headline equals the minimum axis, not the mean.

## 7. References
META-3 keys: `vv-oberkampf-roy`, `gci-roache`, `asme-vv20`, `area-metric`, `pbox`, `false-confidence`,
`pcmm`, `pirt`, `mms`. Depends on META-1 §4.1, FND-1 (p-box), COUP-6 (pedigree).

*(No open questions.)*

## 8. Change log
| Date | Version | Change |
|---|---|---|
| 2026-08-14 | 0.2 | Review fixes (N15, D-G). **N15:** small-n fallback ladder for the model-form band — n=1–2: ASME V&V20-style u_val at the anchor + declared sourced ×3 extrapolation multiplier; n≥3: the area-metric regression; crossover at n=3; branch recorded in the pedigree validation axis (§3.3). **D-G:** operational blind-mode rule added (§3.6): blind = design spec + universal physics/closures + technology-class data from other hardware with cited bands; nothing measured on the engine under test; anchor-derived values calibrated-only, labeled; mechanical provenance enforcement (VISION_SCOPE §9 v1.4). |
| 2026-07-21 | 0.1 | Initial draft. 4-tier validate-bottom-up/predict-top-down hierarchy; GCI solution verification (observed vs formal order, Fs 1.25/3, out-of-asymptotic → unverified); Ferson area-metric + u-pooling → regressed 95% prediction-interval model-form band (widens with extrapolation; no calibrated term extrapolated); p-box as the single honest non-collapsible bound; PIRT assumption register; maturity → PCMM weakest-link mapping. |
