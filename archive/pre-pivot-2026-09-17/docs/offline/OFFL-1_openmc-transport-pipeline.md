# OFFL-1 — OpenMC Transport Pipeline

| Field | Value |
|---|---|
| **ID** | OFFL-1 |
| **Family** | OFFL (Offline pipeline) |
| **Status** | Reviewed (2026-08-14) |
| **Depends on** | FND-5 (table schema), OFFL-2 (perturbed libraries); consumed by SOLV-2, SOLV-4 |
| **Version** | 0.3 (2026-08-14 review fixes: N12, N13, N14, E-1, E-2) |

---

## 0. Purpose

OFFL-1 is the **offline OpenMC pipeline** that produces every nuclear datum the runtime consumes: k-eff,
deposition/dose kernels, few-group kinetics parameters + flux/adjoint shapes, reactivity-coefficient tables,
and the **coupled neutron-photon multigroup libraries that drive the runtime deterministic Sₙ operator
(SOLV-2)**. Its kernel products are **the precomputed solution mode of the one SOLV-2 transport operator**
(E-1) — the same operator solved offline within a tabulated geometry class, never a second transport owner.
It runs once per table version and emits versioned HDF5 (FND-5). Python (OpenMC) never runs at simulation
time; the runtime consumes frozen kernels — there is **no runtime Monte Carlo** (VISION_SCOPE §4.2).

Read after FND-5 (the table contract) and alongside SOLV-2 (the Sₙ consumer) and SOLV-4 (the kinetics
consumer).

## 1. Scope & razor ruling
**Owns:** OpenMC geometry classes; k-eff eigenvalue calculations; KERMA deposition/dose kernels (the
**precomputed solution mode** of SOLV-2's operator, E-1); **flux-shape ψ and adjoint-shape ψ† products**
(SOLV-4 §3.2); kinetics parameters (β_i, λ_i, Λ); reactivity-coefficient sweeps ρ(state); **coupled
neutron-photon multigroup + multi-delayed-group XS generation** for SOLV-2's Sₙ (incl. the photon library,
N14); birth spectra; **Stage-2 nuclear products** (depletion reactivity decrement Δρ(B); MT-444
damage-energy kernels, unit-scaled to dpa in the loader; activation/decay-heat sources); the **factorized**
OpenMC re-runs that turn each OFFL-2 nuclear-data realization into a coefficient row (§3, N13); fidelity
tiering. **Defers:** covariance *sampling* → **OFFL-2**; the Sₙ *solve* → **SOLV-2**; kinetics *evolution*
→ **SOLV-4**; table schema/interpolation → **FND-5**.

**Razor ruling (Reaction Razor):** nuclear energy release flows through the neutron/photon field, so transport
is simulated to completeness — but **offline, at its correct scale** (fine lattice calc → homogenized
constants → coarse runtime cell, FND-2 §3.4.1). **Firewall (§10.1):** OpenMC + open libraries only; no MCNP,
no RSICC, no restricted evaluations, ever.

## 2. Interfaces & contracts
| Interface | Consumed by | Contract |
|---|---|---|
| **Reactivity table** ρ(state) | SOLV-4 | ρ vs (T_fuel, T_mod, density, drum-angle **and, for NSWR, void fraction + salt species/enrichment**); `interp_rule` monotone-cubic (it is differenced) |
| **Kinetics parameters** | SOLV-4 | β_i, λ_i, Λ (adjoint-weighted) |
| **Flux/adjoint shapes** ψ(𝐫,E), ψ†(𝐫,E) | SOLV-4 | per geometry class, on the **same state + control axes as ρ** (E-2) — the quasi-static factorization + importance weighting of SOLV-4 §3.2 |
| **Coupled n-γ multigroup library** | SOLV-2 | neutron group XS + scattering matrices + χ_p + per-delayed-group χ_d + ν; **photon group structure, photoatomic group XS + scattering, neutron→photon production matrices, γ KERMA factors** (N14) — the object SOLV-2 §3.2 sweeps, both particle classes |
| **KERMA deposition/dose kernels** | SOLV-2 | multilinear positivity-preserving; **globally consistent n/γ decomposition** (§3, N12); per geometry class on the state + **control-state** axes (E-2); valid only as SOLV-2's declared precomputed mode (E-1) |
| **Birth spectra** | SOLV-4 | prompt Watt + prompt-γ + delayed-group spectra |
| **Depletion reactivity** Δρ(B) | SOLV-8 | burnup-dependent reactivity decrement vs burnup, per geometry class (long-burn legs; self-zeroing for short-burn NTP) |
| **Damage-energy kernel** (MT-444) | SOLV-8 | raw damage-energy tally, unit-scaled to dpa in the loader (MT-444 is damage *energy*, not dpa) |
| **Activation/decay source** | SOLV-2, SOLV-8 | activation inventories → decay-heat + delayed-photon source spectra (fusion/fission-fragment afterheat path) |

**Invariant:** every table carries FND-5 mandatory metadata (provenance, versions, envelope, per-column
units/interp-rule/uncertainty); offline MC is **Tier-2 reproducible** (seed + statistics recorded in table
provenance). k-eff validation ≤300 pcm is *experimental* (S2); MGXS↔CE ~50 pcm is *code* verification — kept
distinct.

## 3. Method
- **Geometry:** OpenMC-native **CSG** (surfaces/cells/universes/lattices) — offline nuclear models are clean
  axisymmetric reactors CSG expresses best; **DAGMC/STL is skipped** (avoids a MOAB/embree dependency and mesh
  healing; CRUCIBLE's own voxelizer serves the runtime, not OpenMC).
- **Mode partition (E-1):** every kernel product is the **precomputed solution mode of the one SOLV-2
  transport operator** — the same operator solved offline, valid only within its tabulated geometry class +
  envelope (FND-5 refusal outside). **Exactly one mode per particle-class + band per run**, config-declared,
  COUP-2-asserted; **SOLV-2 §3.2 owns the partition rule** (cross-reference, not restated here).
- **k-eff:** eigenvalue mode (inactive/active batch split), validated ≤300 pcm vs **KRUSTY / HEU-MET-FAST-101 /
  MOOSE-VTB** (`krusty`, S2) and ICSBEP fast-HEU-metal benchmarks.
- **Kinetics:** β_eff, Λ via **Iterated Fission Probability** (adjoint-weighted; run IFP over a generation
  window with high active-batch counts — adjoint quantities converge slower than k). [META-3: `ifp-kinetics`]
- **Reactivity sweeps:** **windowed-multipole on-the-fly Doppler** across a structured state grid, with
  **pair-correlated sampling (shared seed across the T+ΔT pair)** so the small Δk survives MC noise in the
  reactivity difference. [META-3: `table-precedent`]
- **Multigroup XS — coupled n-γ (N14):** `openmc.mgxs`/`mdgxs` on a problem-adapted few-group structure,
  verified to reproduce CE k within **~50 pcm**; delayed-group MGXS so SOLV-2's Sₙ and SOLV-4's PKE share
  consistent delayed data. The library is **coupled neutron-photon**: a declared **photon group structure**,
  **photoatomic group XS + scattering matrices** (ENDF/B-VIII.1 photoatomic sublibrary), **neutron→photon
  production matrices** (prompt fission γ, capture γ, inelastic γ), and **γ KERMA factors** — SOLV-2
  transports γ of any origin, so the γ side of the library is not optional. [META-3: `openmc-mgxs`,
  `photon-mgxs`, `endf-b-viii.1`]
- **Flux/adjoint shapes (E-2, N5-intake):** per geometry class, on the **same state + control axes as ρ**:
  forward shape ψ(𝐫,E) from mesh flux tallies; adjoint shape ψ†(𝐫,E) from the **multigroup adjoint solve**
  (transposed scattering/fission operator on the same MGXS — OpenMOC or an offline SOLV-2 adjoint sweep).
  These are the ψ/ψ† products SOLV-4 §3.2's quasi-static closure consumes.
- **Deposition — one globally consistent decomposition (N12):** in a **coupled n-γ run**, heating everywhere
  = **neutron KERMA excluding photon production** + **transported-photon heating** (γ KERMA where photons
  deposit). The former per-cell MT-301/MT-901 choice by range test is **retired**: mixing `heating` and
  `heating-local` across cells double-counts or drops the photon-production energy at every boundary between
  the two regimes — non-conservative kernels by construction. The range test (R vs cell size ℓ, META-3 §6.1)
  survives **only** as the *runtime* local-deposit-vs-transport split of **reaction products** (FND-2
  §3.4.1); it never selects tally scores. Dose via ICRP-116 fluence-to-effective-dose. Kernels carry the
  general **control-state axis** (any declared control coordinate — drum angle at minimum; E-2) so
  commanded power-shape movement is a table dimension, not a frozen configuration. [META-3:
  `kerma-decomposition`]
- **Nuclear-data-perturbation factorization (N13):** naïvely, propagating OFFL-2's covariance sample-set
  means **~300 realizations × the full ~1000-point state sweep ≈ ~8 years** (a 1000-point sweep ≈ ~10 days,
  VISION_SCOPE §8) — unbuildable. The factorization that makes it feasible: the **nominal library runs the
  full state sweep once**; each realization re-runs **only the nominal state point + a small extreme-state
  set** (~4–8 envelope corners of the dominant axes). The realization's effect enters the tables as a
  **state-indexed offset** `Δc(realization; state)` interpolated from {nominal + extremes}, under a
  **declared separability assumption** (the data perturbation varies slowly with state) that is **checked at
  the extremes**: offset spread across the extreme set beyond a declared tolerance = separability failure →
  PIRT-recorded, band widened (never silently averaged). Budget: 300 × (1 + ~6) ≈ **~2100 runs ≈ 1–3 weeks
  at tiered statistics** — the difference between feasible covariance propagation and none.
- **Fidelity tiering:** coarse statistics across the sweep grid, high statistics at anchors, **FW-CADIS**
  weight windows only for deep dose/heating tallies. [META-3: `fw-cadis`]
- **Coupled Picard** (OpenMC↔thermal, ~10× a bare run) is used **only** at a sparse anchor set — to build the
  KRUSTY coupled steady-state anchor and to *verify* the runtime PKE feedback loop reproduces the coupled
  equilibrium; the bulk ρ-grid is uncoupled and the **runtime PKE closes the feedback** (SOLV-4).

## 4. Coupling relationships
- **OFFL-2** supplies perturbed libraries; OFFL-1 re-runs OpenMC per realization **under the §3
  factorization** (nominal + extreme-state set, N13) → the coefficient sample-set. **SOLV-2** sweeps the
  coupled n-γ MGXS + applies the KERMA kernels (its precomputed mode, E-1); **SOLV-4** evolves the kinetics
  from β_i/λ_i/Λ/ρ/birth-spectra + the **ψ/ψ† shapes** (§3, E-2); **SOLV-8** consumes the Stage-2 products (Δρ(B), damage-energy/dpa, activation
  decay heat — activation photons transported by SOLV-2); **FND-5** stores/interpolates; **VAL-2** holds the
  `krusty` anchor.

## 5. Uncertainty & validity
Originates the nuclear-transport statistical uncertainty (Tier-2 MC) and carries the nuclear-data covariance
via the OFFL-2 joint sample-set under the **§3 factorization** (correlated sampling or high history counts so
MC noise ≪ data spread). The factorization's **separability assumption is a declared, checked, PIRT-recorded
epistemic item** (N13) — its extreme-state residual widens the coefficient band, never disappears into it.
Interpolation error budgeted (FND-5 §3.4). Validity envelopes per axis; refuse out-of-envelope. Ladder:
benchmark (ii) vs KRUSTY/ICSBEP.

## 6. Validation plan
1. **KRUSTY** ≤300 pcm k-eff + warm criticals + coupled full-power shape (S2). [META-3: `krusty`]
2. **MGXS↔CE ~50 pcm**, and **OpenMOC** (MIT, open) as the offline deterministic MGXS-consistency oracle.
   [META-3: `openmc-mgxs`, `openmoc`]
3. **Deposition** vs analytic attenuation + an OpenMC multigroup cross-check; **conservation of the n/γ
   decomposition (N12):** summed (neutron-KERMA-excl-γ-production + transported-γ heating) over a closed
   problem equals total energy release — no double-count/gap at any regime boundary.
4. **γ-heating check (N14):** a γ-emitting benchmark (capture-γ-dominated slab) reproduces the photon
   heating profile vs a CE OpenMC reference within the multigroup discretization band.
5. **Separability check (N13):** per-realization offsets at the extreme-state set agree with the nominal
   offset within the declared tolerance, or the run is PIRT-flagged and the band widened.
6. **Provenance/regeneration:** every table regenerable from generator hash + library/tool versions + seed.

## 7. References
META-3 keys: `krusty`, `ifp-kinetics`, `openmc-mgxs`, `photon-mgxs` *(new)*, `kerma-decomposition` *(new)*,
`fw-cadis`, `openmoc`, `table-precedent`, `watt-spectrum`, `sn-deposition`, `radiation-partition`,
`endf-b-viii.1`, `sandy-samples`. Depends on FND-5, OFFL-2; consumed by SOLV-2 (§3.2 mode partition),
SOLV-4 (§3.2 ψ/ψ† closure), VAL-2.

*(No open questions.)*

## 8. Change log
| Date | Version | Change |
|---|---|---|
| 2026-08-14 | 0.3 | Review fixes (N12, N13, N14, E-1, E-2). **N12:** per-cell MT-301/MT-901 range-test tally mixing retired (non-conservative at regime boundaries) → one globally consistent decomposition: neutron KERMA excluding photon production + transported-photon heating everywhere, in a coupled n-γ run; range test survives only as the runtime reaction-product split (FND-2 §3.4.1); conservation test added (§6.3). **N13:** nuclear-data-perturbation factorization specified — nominal full sweep once; per-realization re-runs at nominal + extreme-state set; state-indexed offset with declared, checked, PIRT-recorded separability; budget arithmetic shown (~8 years naïve → ~1–3 weeks). **N14:** coupled n-γ multigroup product added (photon groups, photoatomic XS + scattering, n→γ production matrices, γ KERMA) + γ-heating verification (§6.4). **E-1:** kernels labeled the precomputed solution mode of the one SOLV-2 operator; one mode per particle-class+band per run; SOLV-2 §3.2 owns the partition. **E-2:** kernels + ψ/ψ† gain the general control-state axis (drum angle at minimum); ψ/ψ† flux+adjoint shape products added on ρ's state+control axes (SOLV-4 §3.2 intake). |
| 2026-08-13 | 0.2 | Consistency sweep: Stage-2 nuclear products (depletion Δρ(B), MT-444 damage-energy/dpa kernels, activation/decay sources) added to Owns/§2/§4 — SOLV-8 already named OFFL-1 their supplier; §3 range-test reference qualified to META-3 §6.1. |
| 2026-07-21 | 0.1 | Initial draft. OpenMC-native CSG (DAGMC skipped); k-eff ≤300 pcm vs KRUSTY; IFP β_eff/Λ; windowed-multipole pair-correlated reactivity sweeps (NSWR gets a void-fraction axis); openmc.mgxs/mdgxs MGXS for SOLV-2's Sₙ (~50 pcm CE verification); KERMA MT-301/901 by the range test; FW-CADIS tiering; coupled Picard only at sparse anchors (runtime PKE closes feedback); OpenMOC offline oracle. |
