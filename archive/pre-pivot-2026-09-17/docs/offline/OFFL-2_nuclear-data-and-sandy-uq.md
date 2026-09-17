# OFFL-2 — Nuclear Data & SANDY UQ

| Field | Value |
|---|---|
| **ID** | OFFL-2 |
| **Family** | OFFL (Offline pipeline) |
| **Status** | Reviewed (2026-08-14) |
| **Depends on** | FND-1 (empirical sample-set), FND-5 (schema); feeds OFFL-1, COUP-5 |
| **Version** | 0.1 |

---

## 0. Purpose

OFFL-2 handles the **evaluated nuclear-data libraries** and the **SANDY covariance sampling** that turns
nuclear-data uncertainty into the **empirical sample-set** the outer-loop UQ consumes (FND-1 §3.9). It is the
origin of the nuclear-data covariance uncertainty type. Read after FND-1 (the empirical sample-set contract)
and FND-5 (how the sample-set is stored).

## 1. Scope & razor ruling
**Owns:** ENDF/TENDL/FENDL handling (via OpenMC data API / NJOY / ENDFtk); **SANDY perturbed-library
covariance sampling**; the frozen joint sample-set. **Defers:** the OpenMC re-runs that map realizations →
coefficients → tables → **OFFL-1**; the *runtime* sampling (joint index) → **COUP-5** (per FND-1); the
uncertainty *type* → **FND-1**.

**Razor ruling:** infrastructure feeding the Reaction Razor's data path. **Firewall (§10.1):** open libraries
only; the IAEA-patched NJOY build is pinned in provenance (some FENDL-3.2c data needs it).

## 2. Interfaces & contracts
| Interface | Consumed by | Contract |
|---|---|---|
| **Perturbed libraries** | OFFL-1 | N joint realizations, one coherent perturbed set across all dominant nuclides per index |
| **Empirical sample-set** | COUP-5 (joint index), FND-5 | one row per realization, joint across correlated nuclides — the FND-1 §3.9 carrier |
| **Non-Gaussianity diagnostic** | COUP-6 | an *offline* fitted covariance + area-metric departure from raw — pedigree input, **never a runtime draw source** |

**Invariant:** the **raw joint sample-set is the sole runtime UQ carrier** — consumed by joint index (member
`k` → realization `π(k mod N)`, the seeded permutation of FND-1 §3.5, S2), preserving full
cross-nuclide/cross-coefficient correlation **and** non-Gaussianity. Per-nuclide marginals are forbidden as a runtime carrier.

## 3. Method
- **Libraries:** default transport = **ENDF/B-VIII.1**; **TENDL-2023** for nuclides lacking ENDF/B
  covariances (its TMC-native covariances); **FENDL-3.2c** for fusion neutronics (pin the IAEA-patched
  NJOY2016 in provenance). [META-3: `endf-b-viii.1`, `tendl-2023`, `fendl-3.2c`]
- **SANDY sampling:** sample **MF31/33/34/35** (ν̄, cross-section, angular, χ — not MF33 alone; ν̄ and
  fission-spectrum covariance dominate k/kinetics uncertainty) → **JOINT realizations** across the dominant
  nuclides with **one master seed per realization index**, so realization `i` is reproducible and internally
  consistent (Tier-2). ~300 samples (`sandy-samples`, VISION_SCOPE §8). [META-3: `sandy-mf-coverage`]
- **Dominant-nuclide selection:** frozen by a **one-time OpenMC sensitivity screen** per geometry class (rank
  by uncertainty-weighted sensitivity |dk/dσ|·σ_unc; keep ~95% of response variance) — the pedigree
  justification for truncating to `{U-235, H, C, Be/BeO, Mo, Zr}` (`sandy-samples`); all-nuclide simultaneous
  TMC stays out of scope (§8).
- **If N=300 is insufficient:** raise N in SANDY offline (the §8 budget is comfortable) — **never** fit a
  Gaussian and draw from it at runtime (that is a second sampling seam + a false-confidence Gaussian
  assumption that discards the non-Gaussianity which is the whole reason to keep raw samples).

## 4. Coupling relationships
- **OFFL-1** re-runs OpenMC once per realization → the joint coefficient sample-set; **COUP-5** consumes it by
  joint index; **FND-5** stores it (HDF5 dataset); **FND-1** provides the empirical-set type; **COUP-6** reads
  the truncation justification + non-Gaussianity diagnostic into pedigree.

## 5. Uncertainty & validity
Originates the **nuclear-data covariance** uncertainty type; must preserve cross-nuclide correlation (joint
rows) and non-Gaussianity (raw samples). The dominant-nuclide truncation is a declared PIRT/pedigree entry.
Ladder: infrastructure feeding benchmark-tier transport.

## 6. Validation plan
1. **Sample-set fidelity:** the frozen set reproduces the source ENDF covariance (FND-1 §6 test 3).
2. **Screen coverage:** the dominant-nuclide screen captures ≥ a stated fraction of total response variance
   (recorded).
3. **Non-Gaussianity:** the offline Gaussian-fit diagnostic reports its area-metric departure from raw (so a
   reviewer sees when Gaussianity would have lied).
4. **Provenance:** libraries + NJOY build + seeds pinned; regenerable.

## 7. References
META-3 keys: `endf-b-viii.1`, `tendl-2023`, `fendl-3.2c`, `sandy-samples`, `sandy-mf-coverage`,
`false-confidence`. Depends on FND-1 (empirical sample-set + joint-index), FND-5 (schema); feeds OFFL-1,
COUP-5, COUP-6.

*(No open questions — raw joint sample-set as sole runtime carrier resolved against a fit-and-draw path per
FND-1, 2026-07-21.)*

## 8. Change log
| Date | Version | Change |
|---|---|---|
| 2026-07-21 | 0.1 | Initial draft. ENDF/B-VIII.1 default + TENDL-2023/FENDL-3.2c (IAEA-patched NJOY pinned); SANDY MF31/33/34/35 joint sampling with one seed per realization index → the FND-1 empirical sample-set (raw joint set is the sole runtime carrier; fitted covariance demoted to an offline non-Gaussianity diagnostic); dominant nuclides fixed by a one-time OpenMC sensitivity screen; raise-N-not-fit-and-draw. |
