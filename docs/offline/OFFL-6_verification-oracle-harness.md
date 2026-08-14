# OFFL-6 — Verification-Oracle Harness

| Field | Value |
|---|---|
| **ID** | OFFL-6 |
| **Family** | OFFL (Offline pipeline) |
| **Status** | Reviewed (2026-08-14) |
| **Depends on** | VAL-1 (area-metric/GCI recipes), FND-2 (symmetry indicator), FND-5; feeds COUP-5/COUP-6, SOLV-1 |
| **Version** | 0.1 |

---

## 0. Purpose

OFFL-6 is the **offline oracle harness** that (a) **verifies** the reduced/adaptive-dimension runtime solver
against full-3-D oracles → a **model-form band**, and (b) **calibrates** the magnetic-nozzle electron closure
via kinetic runs → a **coefficient + band**. It runs established open codes **offline, never linked at
runtime** (VISION_SCOPE §5.2/§10), and feeds frozen bands/coefficients into the UQ (COUP-5) and pedigree
(COUP-6). It **owns running the oracle and measuring discrepancy**; **VAL-1 owns turning discrepancy into a
band** (no duplicated machinery).

Read after VAL-1 (the area-metric/GCI recipes it calls) and FND-2 §3.4 (the symmetry indicator it reuses as a
trigger).

## 1. Scope & razor ruling
**Owns:** the oracle runs (OpenFOAM, Athena++, WarpX), their GCI gating, the area-metric discrepancy
*measurement*, and the trigger logic. **Defers:** the area-metric→95%-PI→p-box *recipe* → **VAL-1 §3.3**; the
symmetry *indicator* → **FND-2 §3.4**; **CI code-verification** (matched-setup differential tests) → **VAL-3**
(a different role of possibly the same tool — verification, not validation); the constitutive spine → FND-7.

**Razor ruling:** infrastructure/verification — no new physics; open codes only (firewall-clean). Products are
frozen versioned tables/bands (Tier-2).

## 2. Interfaces & contracts
| Interface | Consumed by | Contract |
|---|---|---|
| **Model-form discrepancy** | VAL-1 → COUP-5/COUP-6 | a GCI-gated area-metric datum at the checked condition; VAL-1 regresses it into the widening band |
| **γ_e calibration** | SOLV-1 (plasma closure) | γ_e(M) as a function with a band that widens off-condition |
| **Trigger flag** | the run (deterministic) | fires on symmetry-break; auto-flags + applies a conservative band; the oracle run itself is human-gated |

**Invariant:** an oracle run is **GCI-gated first** (so discrepancy = model-form, not oracle mesh error); the
trigger *flags and widens conservatively* on every run but **never auto-launches** an oracle (budget); a
calibrated *parameter* (γ_e) is legitimate but a calibrated *discrepancy term* is never extrapolated (META-1
§4.1).

## 3. Method
- **Verification oracles:** **OpenFOAM (`rhoCentralFoam`)** for 3-D compressible flow, **Athena++** for 3-D
  MHD. GCI-gate both (`gci-roache`, VAL-1 §3.2), then measure the **area-metric** discrepancy (VAL-1 §3.3
  recipe) between the reduced-model CDF and the oracle CDF; hand it to VAL-1, which regresses it vs the
  extrapolation coordinate (swirl number, magnetization, Reynolds) into the 95%-PI band → COUP-5's p-box.
  [META-3: `openfoam`, `athena++`, `area-metric`, `gci-roache`]
- **Calibration oracle:** **WarpX** PIC (offline-only, §5.3) calibrates the magnetic-nozzle electron polytropic
  **γ_e** (piecewise 1→5/3 tied to the detachment plane, `mn-polytropic`) as a **function γ_e(M)** with a band
  that widens away from the calibration condition (VAL-1 regression) — a legitimate calibrated *parameter*, not
  a carried discrepancy fudge. Resampling-threshold / particles-per-cell convergence caveat in the band.
  [META-3: `warpx`, `warpx-gammae`, `mn-polytropic`]
- **Trigger (four items, all queue a human-gated batch candidate — none auto-runs):** (1) azimuthal
  resolution **saturated at N_θ^max with the asymmetry indicator still high** (closure/kinetic inadequacy the
  adaptive refinement cannot fix — distinct from FND-2's τ_collapse/τ_expand); (2) a published **azimuthal
  (m=1-class) instability critical parameter** crossed; (3) a **new pedigree corner** (no prior oracle
  coverage); (4) **audit sampling**. Meanwhile the default is a **conservative band derived from FND-2's
  τ_collapse discarded-variance bound** (a provable upper bound); an oracle run **replaces** it with the
  measured (tighter) band. [META-3: `symmetry-indicator`]
- **Determinism/budget:** the trigger flag is a deterministic function of the run's own field (reproducible);
  the oracle campaign is **human-gated + batched**, so a sweep never auto-burns the ~5000 GPU-h budget and
  never silently produces an over-confident result (fail-loud via the conservative default).

## 4. Coupling relationships
- **VAL-1** owns the banding recipe OFFL-6 calls; **FND-2** owns the symmetry indicator OFFL-6 triggers on;
  **VAL-3** runs CI code-verification with the same codes (verification ≠ this doc's validation);
  **COUP-5/COUP-6** consume the band into UQ/pedigree; **SOLV-1** consumes the γ_e calibration.

## 5. Uncertainty & validity
VERIFICATION output = a widened, extrapolation-aware model-form band (via VAL-1). CALIBRATION output = γ_e(M)
with a band widening away from data. Firewall-clean (all open codes; WarpX kinetic offline per §5.3). Ladder:
the oracles are the 3-D/kinetic truth against which reduced-model model-form is bounded.

## 6. Validation plan
1. **Oracle self-checks:** OpenFOAM/Athena++ reproduce their canonical benchmarks before use; GCI observed-order
   check on each run.
2. **γ_e** recovers `mn-polytropic` published values (Little-Choueiri) at the measured condition.
3. **Trigger/band path:** a deliberately symmetry-breaking config fires the flag, applies the τ-derived
   conservative band, and does **not** auto-run an oracle.

## 7. References
META-3 keys: `openfoam`, `athena++`, `warpx`, `warpx-gammae`, `mn-polytropic`, `area-metric`, `gci-roache`,
`pbox`, `symmetry-indicator`. Depends on VAL-1 (recipes), FND-2 (indicator), FND-5; feeds COUP-5, COUP-6,
SOLV-1.

*(No open questions — OFFL-6 measures discrepancy, VAL-1 bands it, VAL-3 does CI verification (distinct roles);
trigger auto-flags + conservative τ-band, oracle runs human-gated/batched; γ_e is a calibrated parameter with a
widening band, not a fudge. Resolved 2026-07-21.)*

## 8. Change log
| Date | Version | Change |
|---|---|---|
| 2026-07-21 | 0.1 | Initial draft. OpenFOAM/Athena++ 3-D verification oracles (GCI-gated → area-metric discrepancy → VAL-1 bands it → COUP-5 p-box); WarpX offline calibration of magnetic-nozzle γ_e(M) as a widening-band parameter; trigger = modes-saturated-and-still-high OR regime-boundary (auto-flag + conservative τ-derived band; oracle runs human-gated/batched to protect the compute budget). Ownership lines drawn vs VAL-1 (banding) and VAL-3 (CI verification). |
