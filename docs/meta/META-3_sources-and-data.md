# META-3 — Sources, Data Provenance & Research Ledger

| Field | Value |
|---|---|
| **ID** | META-3 |
| **Family** | META |
| **Status** | Draft — living document (grows for the whole project) |
| **Depends on** | META-1, META-2 |
| **Version** | 0.3 (v1.3 keys added, §6.5) |

---

## 0. Purpose

The **single ledger** of every hard number, equation, correlation, table, tool, and validation
datum the project uses — each sourced and derivation-traceable so it drops straight into a research
paper's methods/references with no re-hunting. This is the project's bibliography *and* its data-of-record.

Two hard rules (Principles 7 & 8, META-1):
- **Every datum a doc uses has an entry here.** Design docs cite `META-3:<key>`, never a bare paper.
- **An input a paper cannot cite is an input we do not use.** No entry ⇒ not usable.

This doc is **living**: it is seeded below from VISION_SCOPE.md and grows every time a design doc
pulls in a constant, correlation, or dataset. Items still to be gathered are marked
**`[RESEARCH-PENDING]`** with the doc that will need them.

---

## 1. Data-hygiene schema

Every ledger entry carries these fields (a datum missing any starred field cannot be used in a
result):

| Field | Meaning |
|---|---|
| **key** | stable citation handle, e.g. `bosch-hale-1992`, `rl10-tm107318` |
| ★ **value / expression** | the number, table description, or equation |
| ★ **units** | SI, or native + SI conversion (META-2 §2.1) |
| ★ **uncertainty** | value ± (meaning), or band, or "model-form X%", or "exact-by-definition" |
| ★ **source** | full citation (authors, title, venue/report no., year, identifier: DOI/NTRS/arXiv/ICSBEP) |
| **retrieval** | date obtained + where (URL / archive / repo), for reproducibility |
| **derivation note** | how a *derived* value was obtained (which equation, which upstream keys) |
| **consumers** | doc IDs that use it (back-reference; keeps the ledger and catalog in sync) |
| **license / export status** | for tools & data: license + export-control clearance (firewall §10) |

The same provenance travels into every generated HDF5 table (FND-5 metadata: generator hash,
library versions, source keys) and into every results bundle (FND-6), so a result → table → source
chain is fully walkable. **This is the data hygiene that makes the whole thing paper-defensible.**

**Archival rule** *(Ben, 2026-07-14):* every **open** source (PDF/dataset/report) is copied into a
local **read-only source cache** keyed by its `key`, so "retrieval" is reproducible years later even
if the upstream moves or is retracted; **paywalled** sources are cited by DOI and not cached. The
cache location and integrity (hash per file) are recorded; the cache is the retrieval-of-record.

### 1.1 Flow into papers

A campaign paper's methods section is assembled by walking, for each headline result, its pedigree
path (COUP-6) → the tables it touched (FND-5) → their META-3 source keys. The references list is the
union of those keys. Nothing is cited that isn't here; nothing here that a result used is omitted.

---

## 2. Tool inventory

| key | Tool | Role | Version (pin) | License | Export status | Consumers |
|---|---|---|---|---|---|---|
| `rust` | Rust toolchain | Runtime core | `[RESEARCH-PENDING: pin]` | MIT/Apache-2.0 | clear | all runtime |
| `openmc` | OpenMC | Offline neutron/photon transport, k-eff, kernels, reactivity sweeps | `[pin]` | MIT | clear (embeds/redistributes freely) | OFFL-1/2 |
| `cantera` | Cantera | Offline equilibrium chemistry / thermo | `[pin]` | BSD-3 | clear | OFFL-3 |
| `nasa-cea` | NASA CEA + RocketCEA | Rocket equilibrium performance, ~2000-species Glenn DB | `[pin]` | open (github.com/nasa/cea) / GPLv3 wrapper (pipeline-only) | clear | OFFL-3 |
| `geant4` | Geant4 | Offline antiproton annihilation/fission (FTFP/INCL/CHIPS) | `[pin]` | Geant4 license (permissive) | clear | OFFL-4 |
| `njoy`/`endftk` | NJOY / ENDFtk | Nuclear data processing | `[pin]` | open | clear | OFFL-2 |
| `sandy` | SANDY | Perturbed-library covariance sampling | `[pin]` | open | clear | OFFL-2 |
| `hdf5` | HDF5 (+ Rust `hdf5` crate) | The cross-language table seam | `[pin]` | BSD-style | clear | FND-5, all OFFL |
| `openfoam` | OpenFOAM | Offline 3-D compressible-flow oracle | `[pin]` | GPL (oracle only, never linked) | clear | OFFL-6 |
| `athena++` | Athena++ | Offline 3-D MHD oracle | `[pin]` | open (oracle only) | clear | OFFL-6 |
| `warpx` | WarpX | Optional kinetic runs to calibrate mag-nozzle closure | `[pin]` | open | clear | OFFL-6 |
| **Oracles (never linked)** | MULTI-IFE, SNEC, magpylib, pyviewfactor, SALib | Cross-check only, run in Python during validation | — | various | clear | VAL, SOLV cross-checks |
| **FORBIDDEN** | MCNP, RSICC codes, restricted evaluations | — | — | **export-controlled: never used, never requested** (firewall §10.1) | — |
| **Conditional** | FLUKA | Antiproton, *only if licensable* for an independent student; else Geant4-only with wider band | — | non-profit/institutional terms | check before any use (§10.3) | OFFL-4 |

## 3. Nuclear-data libraries

| key | Library | Role | Notes | Consumers |
|---|---|---|---|---|
| `endf-b-viii.1` | ENDF/B-VIII.1 | Primary evaluated nuclear data | pointwise + multigroup via OpenMC data API | OFFL-1/2 |
| `tendl-2023` | TENDL-2023 | Supplementary evaluations | for nuclides sparse in ENDF | OFFL-1/2 |
| `fendl-3.2c` | FENDL-3.2c | Fusion-relevant evaluations | | OFFL-1/2 |
| `sandy-samples` | SANDY perturbed sets | Covariance UQ | ~300 samples × dominant nuclides (U-235, H, C, Be/BeO, Mo, Zr) | OFFL-2, COUP-5 |

## 4. Physical constants

Single pinned source: **the latest official CODATA recommended values as of project start**, with the
exact set version recorded here and read by the FND-1 constants module. (Each value carries its CODATA
uncertainty; exact-by-definition constants like g₀, c are flagged as such.) `[RESEARCH-PENDING: stamp
the exact CODATA release + retrieval date when the constants module is built.]`

| key | Symbol | Value | Units | Note |
|---|---|---|---|---|
| `g0` | g₀ | 9.806 65 | m/s² | standard gravity, **exact by definition** (Isp convention) |
| `c` | c | 2.997 924 58 ×10⁸ | m/s | speed of light, exact by definition |
| `kB` | k_B | `[CODATA]` | J/K | Boltzmann |
| `NA` | N_A | `[CODATA]` | 1/mol | Avogadro |
| `e` | e | `[CODATA]` | C | elementary charge |
| `me`,`mp`,`mn`,`mu` | masses | `[CODATA]` | kg | electron/proton/neutron/atomic-mass-unit |
| `eV` | eV | `[CODATA]` | J | eV→J conversion |
| `barn` | b | 1×10⁻²⁸ | m² | exact |

(Full CODATA set filled in when FND-1 fixes the constants module; each value carries its CODATA
uncertainty.)

## 5. Validation-anchor data (seed from VISION_SCOPE §9)

Expanded per-anchor specs live in **VAL-2**; this is the *sourced data-of-record* each anchor rests on.

> *(v1.3 note, 2026-07-20: the `consumers` back-reference columns in §5 and §6/§6.1–§6.3 still use
> **pre-v0.5 SOLV/OFFL IDs**; map old→new via META-0 §6.1's merge table until each consuming doc is
> written and updates them. §6.4/§6.5 already use current IDs. Physics content is unaffected.)*

| key | Anchor | Data-of-record | Source | Target | Consumers |
|---|---|---|---|---|---|
| `rl10-tm107318` | RL10A-3-3A | 73.4 kN vac thrust; Isp ≈ 444 s; Pc = 475 psia (→ 3.275 MPa); MR = 5.0; + component maps | NASA TM-107318 (NTRS 19970010379); 1966 design report | ≤ 2% Isp, c* (S1) | SOLV-1, OFFL-3, VAL-2 |
| `nerva-pewee` | NERVA Pewee | >500 MW; fuel-exit 2556 K; chamber 1833 K @ 18.6 kg/s; 4275 kPa; peak ideal vac Isp 901 s | NASA-CR-184270 (NTRS 19920005899); LA-4217-MS | ≤ 10% state points (S2) | SOLV-5/6, OFFL-1, VAL-2 |
| `nrx-a6` | NRX-A6 | 1120 MW × 60 min; Pc 4089 kPa; 32.7 kg/s | NASA-CR-184270; WANL-TNR-223/224 | ≤ 10% | SOLV-1/5, VAL-2 |
| `xe-prime` | XE-Prime | 1140 MW; 244.75 kN (thrust-stand); Isp ~710 s; 2272 K; 35.8 kg/s; + operating map | NASA-CR-184270; Aerojet RN-S-510 | ≤ 10%, map shape | SOLV-1/5, VAL-2 |
| `krusty` | KRUSTY | criticality benchmark + warm criticals + full-power coupled run | ICSBEP HEU-MET-FAST-101 / KRUSTY-SPACE-EXP-001; MOOSE VTB KRUSTY model (open) | ≤ 300 pcm k-eff; coupled-transient shape (S2) | OFFL-1, SOLV-5, VAL-2 |
| `nstar` | NSTAR | 16-level throttle table: 0.58–2.57 kW, 19–92 mN, Isp 1950–3120 s | Goebel & Katz, *Fundamentals of Electric Propulsion* Ch. 9 (JPL Descanso, free) | within published η bands | SOLV-13, VAL-2 |
| `vasimr-vx200` | VASIMR VX-200 | 5.8 ± 0.4 N; Isp 4900 ± 300 s @ 200 kW (plume-derived) | Longmier et al. | within stated error bars | SOLV-7, VAL-2 |
| `radhydro-analytics` | Rad-hydro analytics | Sod / Noh / Sedov exact solutions; + MULTI-IFE & SNEC cross-runs | open (Sod 1978; Noh 1987; Sedov 1959) | standard convergence | SOLV-8, VAL-2 |
| `antimatter-ps177` | Antimatter | antiproton-induced fission probabilities (U-238/Bi/Pb/Au); prong multiplicities; annihilation spectra | PS177; AEgIS-era arXiv:2407.06721 | within declared 15–30% model-form band | OFFL-4, SOLV-8, VAL-2 |
| `beavrs` | Coupling shakedown | LWR benchmark (not space-like; shakes down coupling machinery only) | MIT BEAVRS open benchmark | qualitative | COUP-2, VAL-2 |

## 6. Equation & correlation ledger

Governing equations and closures, each cited so a paper can reproduce the method. Seeded from
VISION_SCOPE §5.1; grows as each SOLV/OFFL doc is written. `[RP]` = `[RESEARCH-PENDING]` — the
exact reference/coefficients get pinned when the consuming doc is written.

| key | Equation / correlation | Source | Uncertainty | Consumers |
|---|---|---|---|---|
| `hllc` | HLLC approximate Riemann solver (Godunov flux) | Toro, *Riemann Solvers and Numerical Methods for Fluid Dynamics* (3rd ed.) | numerical (scheme-order) | SOLV-1/2 |
| `maccormack-nozzle` | MacCormack quasi-1D nozzle (didactic baseline) | Anderson, *Modern Compressible Flow* / *CFD* | numerical | SOLV-1 |
| `drift-flux` | Quasi-1D drift-flux / homogeneous-relaxation two-phase closures | `[RP: pick closure set(s); NSWR dispute → run as model-form ensemble]` | model-form (disputed → UQ) | SOLV-2 |
| `pke` | Point/few-group reactor kinetics equations + delayed-neutron params | `[RP: standard PKE; β, Λ per fuel]`; coefficients from OpenMC | data + model-form | SOLV-5 |
| `bosch-hale-1992` | Fusion reactivity ⟨σv⟩ parameterizations (D-T, D-D, D-³He) | Bosch & Hale, *Nucl. Fusion* **32** (1992) 611 | fit-level (stated in ref) | SOLV-9 |
| `pb11-sikora-2016` | p-¹¹B reactivity | Sikora & Weller (2016) | stated | SOLV-9 |
| `pb11-putvinski-2019` | p-¹¹B reactivity (alt.) | Putvinski et al. (2019) | stated | SOLV-9 |
| `biot-savart` | Biot–Savart law for coil B-field | standard EM (Jackson) | numerical (quadrature) | SOLV-10 |
| `view-factor-mc` | Monte-Carlo ray-cast view factors + radiosity | standard (Modest, *Radiative Heat Transfer*); pyviewfactor as oracle | MC statistics | SOLV-4 |
| `vn-richtmyer` | von Neumann–Richtmyer artificial-viscosity staggered Lagrangian hydro | von Neumann & Richtmyer (1950); Bowers & Wilson | numerical | SOLV-8 |
| `radhydro-fld` | Grey → few-group flux-limited diffusion | Castor, *Radiation Hydrodynamics*; Levermore-Pomraning limiter | model-form | SOLV-8 |
| `roark` | Thin-shell hoop/thermal stress, pressure/burst formulas | Roark's *Formulas for Stress and Strain* | analytic | SOLV-11 |
| `bprime` | B′ transfer-coefficient ablation tables | generated offline from equilibrium chemistry | model-form | SOLV-16, OFFL-3 |
| `saha-qeos` | Saha ionization + QEOS-style analytic EOS/opacity (default) | Saha; More et al. QEOS (1988) | model-form (upgradeable to SESAME/TOPS) | OFFL-5, SOLV-8 |
| `stopping-astar` | ASTAR/PSTAR/SRIM-class stopping powers & ranges | NIST ASTAR/PSTAR; SRIM (Ziegler) | tabulated | OFFL-6, SOLV-12 |
| `guiding-center` | Guiding-center / gyro-tracking in sampled B-field | standard plasma (Northrop) | model-form (fragment charge-state: wide band) | SOLV-12 |
| `saltelli-sobol` | Saltelli sampling + Sobol sensitivity indices | Saltelli et al.; SALib as oracle | statistical | COUP-5 |
| `latin-hypercube` | Latin-hypercube ensemble sampling | McKay, Beckman & Conover (1979) | statistical | COUP-5 |

### 6.1 Energetic-particle transport, plasma & reaction sourcing (researched 2026-07-14)

Gathered to inform the Rule-12 unification of SOLV-12 / SOLV-7 / SOLV-9 / OFFL-4 / OFFL-6 (see
META-0 v0.4). `[RP]` = coefficients/tables to pin when the consuming doc is written. Confidence flags
carried from the research pass.

**Unified stopping / friction (cold → plasma continuum) — one dielectric-response law, `S_elec + S_nuclear + S_ion`:**
| key | Content | Source | Consumers |
|---|---|---|---|
| `stopping-rpa-lda` | Enhanced RPA-LDA average-atom dielectric electronic stopping (+ nuclear + ionic), continuous cold↔plasma; agrees with PSTAR/IAEA across the table + reproduces plasma data; open-source (GitHub) — target reference architecture | *An Enhanced RPA-LDA Model…*, arXiv:2606.30978 (1 Jul 2026) *(verified 2026-07-14)* | OFFL-6, SOLV-12 |
| `stopping-li-petrasso` | Plasma charged-particle stopping (fast ions on Maxwellian) | Li & Petrasso, PRL **70**, 3059 (1993); Erratum PRL **114**, 199901 (2015); Zylstra et al., Phys. Plasmas **26**, 122703 (2019) | OFFL-6, SOLV-12 |
| `stopping-bps` | Brown-Preston-Singleton exact-order plasma stopping | Brown, Preston & Singleton, Phys. Rep. **410**, 237 (2005), arXiv:physics/0501084 | OFFL-6 |
| `stopping-cold` | Cold-matter end: Bethe-Bloch + ZBL nuclear + tables | NIST ASTAR/PSTAR/ESTAR; SRIM (Ziegler et al., NIM B **268**, 1818, 2010) | OFFL-6 |
| `zeff-betz` | Effective-charge closure Z_eff(v;M) (heals dressed-ion↔bare-ion seam) | Betz, Rev. Mod. Phys. **44**, 465 (1972); Sigmund, NIM B **174**, 535 (2001) | OFFL-6, SOLV-12 |

**Fast-ion slowing-down / self-heating:**
| key | Content | Source | Consumers |
|---|---|---|---|
| `spitzer-slowing` | Spitzer slowing-down time τ_s ∝ T_e^{3/2}/(n_e Z²); friction shared with stopping (`dE/dt=v·dE/dx`) | Spitzer (1962); Trubnikov (1965); NRL Plasma Formulary (2019) | SOLV-12, SOLV-7 |
| `stix-critical-energy` | Critical energy E_c (electron vs ion heating split); E_c≈33·T_e for α in D-T | Stix, Plasma Phys. **14**, 367 (1972) | SOLV-12, SOLV-7 |

**Orbits & radiative cooling (self-zeroing terms):**
| key | Content | Source | Consumers |
|---|---|---|---|
| `gc-boris` | Coupled guiding-center↔Boris pusher, switched on adiabaticity ε_gc=max(ρ_L/L, 1/ω_cτ) — the no-seam orbit integrator | Bacchini et al., ApJS **251**, 10 (2020); Northrop (1963) | SOLV-12, SOLV-7 |
| `brems-emc` | Bremsstrahlung as radiation-reaction on f(v); Bethe-Heitler + Elwert + relativistic e-e | Embréus, Stahl & Fülöp, New J. Phys. **18**, 093023 (2016) | SOLV-12, SOLV-7 |
| `synchrotron-larmor` | P_sync ∝ q⁴B²γ²β_⊥²/m⁴ (ion term self-zeros by 1/m⁴) | Larmor/synchrotron standard; Rossi (1952) | SOLV-12, SOLV-7 |

**Fusion product sourcing & neutrons:**
| key | Content | Source | Consumers |
|---|---|---|---|
| `bosch-hale-1992` | ⟨σv⟩ + σ(E) Padé fits, D-T/D-D/D-³He; reactivity fit ~0.25%, σ fit ~2% (do not conflate); B_G, m_r c², C₁…C₇ tables VII–VIII | Bosch & Hale, Nucl. Fusion **32**, 611 (1992) + Erratum NF **33**, 1919 (1993) | SOLV-9 |
| `pb11-reactivity` | p-¹¹B ⟨σv⟩; σ resonance ~600 keV vs reactivity peak ~250 keV (Ti); modern +20–30% | Nevins & Swain, NF **40**, 865 (2000); Sikora & Weller, JFE **35**, 538 (2016); Putvinski et al., NF **59**, 076018 (2019); Tentori & Belloni, NF **63**, 086001 (2023) | SOLV-9 |
| `dt-neutron-broadening` | Thermal broadening of the 14.1 MeV D-T neutron peak: **FWHM(keV) ≈ 177.2·√T_i** (T_i in keV); the D-D 2.45 MeV peak is FWHM ≈ 82.6·√T_i *(verified 2026-07-14 — do not swap these two coefficients)* | Ballabio, Källne & Gorini, NF **38**, 1723 (1998); Brysk (1973) | SOLV-9, OFFL-1 |
| `watt-spectrum` | Fission source χ(E)=C·exp(−E/a)·sinh(√(bE)); ²³⁵U-th a=0.988, b=2.249 MeV⁻¹. **≥3 incompatible (a,b) conventions — verify sampler** | standard; ENDF/B | OFFL-1 |
| `nrt-arc-dpa` | Damage: OpenMC `damage-energy` MT-444 (raw T_dam, NOT dpa) → NRT N_d=0.8T_dam/2E_d, E_d≈40 eV (Fe); arc-dpa ξ≈0.286 (Fe, fails for W) | Norgett-Robinson-Torrens; Nordlund et al., Nat. Commun. **9**, 1084 (2018) | SOLV-16, OFFL-1/2 |

**Antiproton annihilation source (factor-of-several model-form band):**
| key | Content | Source | Consumers |
|---|---|---|---|
| `pbar-annihilation` | ⟨n_π⟩≈5 (~3 charged), mean π KE ~230 MeV; local deposit ~150–300 MeV/annihilation; ~16 neutrons/annihilation on U (~52% via fission); fragments ~160 MeV | von Egidy; PRC **45**, 2332 (1992); PRC **63**, 034616 (2001) | OFFL-4 |
| `pbar-fission-prob` | Prompt fission probabilities U-238/Bi/Pb/Au (per-nucleus values paywalled) | Bocquet et al., Z. Phys. A **342**, 183 (1992) | OFFL-4 |
| `pbar-modern-data` | Charged-track/fragment multiplicities, multi-target — model constraints | ASACUSA arXiv:2407.06721 (2024); AEgIS/ASACUSA arXiv:2503.04868 (2025) | OFFL-4 |
| `pbar-list-spread` | Geant4-list disagreement: FLUKA closest; FTFP underestimates heavy fragments ~12×; neutron-dose spread up to ~50% → carry factor-of-several systematic | ASACUSA benchmark; Galoyan & Uzhinsky arXiv:1610.08341 | OFFL-4 |
| `ican-ii` | Antiproton-catalyzed microfission driven-subcritical pellet; ~10¹⁷ (~140 ng) p̄ Mars-class; annihilations to start DT burn ≈10²¹/k² | Smith/Lewis/Gaidos, AIP CP **324**, 555 (1995); Gsponer & Hurni, arXiv:physics/0507125 | SOLV-5, SOLV-8, OFFL-4 |

**Plasma energy balance, magnetic nozzle & radiation partition:**
| key | Content | Source | Consumers |
|---|---|---|---|
| `mn-polytropic` | Electron closure T_e/T_e0=(n_e/n_e0)^{γ_e−1}; γ_e≈1.15±0.03 (measured); piecewise 1.13(magnetized)→1.30(detached), break = detachment plane | Little & Choueiri, PRL **117**, 225003 (2016); Andrews et al., arXiv:2111.11698 (2022) | SOLV-7, OFFL-8 |
| `mn-thrust` | Dual electrothermal/electromagnetic thrust; detachment (ion demagnetization) = **largest declared unknown** | Ahedo & Merino, Phys. Plasmas **17**, 073501 (2010); Merino & Ahedo (2017) | SOLV-7 |
| `plasma-balance-terms` | Minimal-complete source/sink menu (self-zeroing by regime) | Emoto, Takahashi & Takao, Front. Phys. **9**, 779204 (2021); Shubov, arXiv:2104.06251 (2021) | SOLV-7 |
| `radiation-partition` | Photon "emit-once/transport-once" architecture; single thermal/suprathermal boundary; aligned MC-floor/opacity-ceiling energy edges (HYDRA 3-package pattern) | Marinak et al. (HYDRA, 2013); Wollaeger et al., ApJS **209**, 36 (2013); Embréus et al. (2016) | COUP-2, SOLV-4/6/7/8 |
| `pb11-brems-limit` | p-¹¹B: brems ≥ fusion at T_e=T_i (no ignition); narrow window only for T_e<T_i; relativistic + e-e brems mandatory, no tail double-count | Rider, Phys. Plasmas **4**, 1039 (1997); Putvinski et al. (2019) | SOLV-7, SOLV-9 |
| `plasma-eos-opacity` | QEOS analytic default (TF electrons + Cowan ions); Saha ⟨Z⟩ w/ McWhirter validity; TOPS/OPLIB Planck/Rosseland means + ⟨Z⟩ | More et al., Phys. Fluids **31**, 3059 (1988); TOPS/OPLIB (LANL) | OFFL-5, SOLV-7/8 |

**Reaction↔transport interface contract:** every reaction/annihilation module emits a phase-space
source `S(r,E,Ω,t)` + species + strength (particles/s); transport owns everything after birth. Handoff
invariants: Σ birth KE = Q; energy transported must be subtracted from local deposit (OpenMC KERMA
MT-301/901 discipline); per-event multiplicities correct; **range test** R vs cell size ℓ decides
local-deposit vs transport. Mirrors OpenMC `IndependentSource` / R2S operator-split (NF **64**, 076023,
2024). Consumers: COUP-1/2, SOLV-9, OFFL-4, SOLV-5.

### 6.2 World-state grid: structure, coupling & geometry (researched 2026-07-14)

Method references informing FND-2 (grid), FND-3 (geometry/voxelization), COUP-1 (mapping), COUP-2
(conservation audit). To be pinned per-doc as those are written; captured here so the research isn't lost.

**Sparse grid data structure & layout:**
| key | Content | Source |
|---|---|---|
| `vdb` | Shallow-wide fixed-depth tree, 8³ leaves, bitmask topology, ValueAccessor node-chain cache (~3× on coherent access) | Museth, *ACM TOG* **32**(3):27 (2013) |
| `nanovdb` | Static-topology linearized VDB (pointer-less, computed offsets) — matches "topology frozen at config time" | NanoVDB (Academy Software Foundation / NVIDIA) |
| `amrex` | Block-structured AMR, `FArrayBox` component-major SoA, ghost-cell halos | Zhang et al., arXiv:2009.12009 (2021) |
| `p4est` | Linear octree + Morton/Z-order (borrow for iteration/partition order, not random access) | Burstedde, Wilcox & Ghattas, *SIAM J. Sci. Comput.* **33** (2011) |
| `aosoa-cabana` | AoSoA (tile width = SIMD lanes) as the SoA/cache compromise | Cabana (ECP); Wikipedia AoS/SoA |

**Determinism & conservation reductions:**
| key | Content | Source |
|---|---|---|
| `repro-sum` | Reproducible/binned summation — order-independent totals | Demmel & Nguyen, *ACM TOMS* **46**(3) (2020) |
| `gamer2-determinism` | Bitwise reproducibility across thread/rank counts via forced deterministic FP order | Schive et al., arXiv:1712.07070 |
| `fp-nonassoc` | FP non-associativity, Kahan/compensated summation | Goldberg, *ACM Comput. Surv.* (1991) |

**Refinement (config-time, frozen at run):**
| key | Content | Source |
|---|---|---|
| `lohner` | Normalized-second-difference error indicator (no time history — ideal one-shot) | Löhner, *CMAME* **61** (1987) |
| `berger-amr` | Structured AMR, Richardson/truncation tagging, 2:1 balance, proper nesting, Berger-Rigoutsos clustering | Berger & Oliger, *JCP* **53** (1984); Berger & Colella, *JCP* **82** (1989) |

**Accountant↔solver coupling & conservative transfer:**
| key | Content | Source |
|---|---|---|
| `moose-multiapp` | Main-app orchestrates outer loop + Picard; sub-apps are black boxes; directional Transfers; `MultiAppConservativeTransfer` rescales to preserve an integral | Permann et al., *SoftwareX* **11**:100430 (2020); Gaston et al., *Nucl. Eng. Des.* **239** (2009); Giudicelli et al., *Front. Nucl. Eng.* (2025) |
| `precice` | Peer coupling, black-box adapters, consistent vs conservative maps, acceleration | Bungartz et al., *Comput. Fluids* **141** (2016); Chourdakis et al., *ORE* **2**:51 (2022) |
| `consistent-vs-conservative` | Conservative map = transpose of consistent map; can't have both on non-matching meshes | de Boer, van Zuijlen & Bijl, *CMAME* **197** (2008) |
| `supermesh` | Exact conservative Galerkin/L2 projection via supermesh (conservative but not bounded) | Farrell et al., *CMAME* **198** (2009); Farrell & Maddison, *CMAME* **200** (2011) |
| `donor-cell` | Exact-intersection-volume transfer — conservative AND bounded, first-order (the practical minimum) | (standard FV; de Boer 2008) |
| `peaceman-well` | Point/1-D ↔ 3-D block coupling via a precomputed geometric weight (well index) | Peaceman, *SPE J.* (1978, 1983) |
| `koch-line-source` | Kernel-distributed (tube) line source restores regularity vs singular Dirac line (D'Angelo) | Koch, Schneider, Helmig & Jenny, *JCP* (2020); DuMuX 3 |
| `axi-2pir` | 2-D↔3-D maps must carry the 2πr revolved-volume factor or mass is silently lost | Chen, Heger et al., arXiv:1112.3033 |

**Operator-split / Picard wall coupling & audit:**
| key | Content | Source |
|---|---|---|
| `robin-robin-cht` | Dirichlet-Neumann CHT unstable for Biot>1; Robin-Robin with optimal coefficient → unconditional stability | Giles, *IJNMF* **25** (1997); Errera et al., *IJHMT*/*JCP* (2018/2019) |
| `aitken-iqnils` | Fixed-point acceleration: Aitken dynamic relaxation; IQN-ILS quasi-Newton | Küttler & Wall, *Comput. Mech.* **43** (2008); Degroote (preCICE) |
| `fv-telescoping` | Flux-form FV telescopes → global audit = boundary/port flux accounting; Lax-Wendroff justifies conservative form | LeVeque, *Finite Volume Methods for Hyperbolic Problems* (2002) |
| `modelica-connector` | Port = flow connector; sum of through-flows at a node = 0 (the port-accounting invariant) | Modelica connector spec |

**Geometry, sub-cell materials & recession:**
| key | Content | Source |
|---|---|---|
| `eb-cutcell` | Embedded-boundary: per-cell volume fraction κ + face apertures α; EB centroid/normal/area | AMReX-EB docs; EBChombo (Colella et al.) |
| `state-redistribution` | Small-cut-cell stability without merging — conservative + 2nd-order + stable | Berger & Giuliani, arXiv:2005.05734 (2020) |
| `gen-winding-number` | Robust inside/outside on imperfect STL (threshold 0.5); fast Barnes-Hut eval → partial volume fractions by sampling | Jacobson, Kavan & Sorkine-Hornung, *SIGGRAPH* (2013); Barill et al., *SIGGRAPH* (2018) |
| `vof-plic` | Volume-of-fluid fractions (mass-conservative) + PLIC plane reconstruction (volume-exact, analytic offset) for surface areas | Hirt & Nichols (1981); Youngs (1982); Scardovelli & Zaleski |
| `mc33-dc` | Deterministic manifold isosurface: Marching Cubes 33 (asymptotic decider) / Manifold Dual Contouring; MC area is ~8% biased → derive area from PLIC | Chernyaev (1995); Nielson & Hamann (1991); Ju et al. (2002); Schaefer et al. (2007) |
| `multimat-closure` | Mixed-cell closure: keep segregated per-material ρ/e/T/EOS, reconcile to one (p,u) by rate-based pressure relaxation (Tipton); MOF sub-cell reconstruction | Kamm & Shashkov (LANL); Tipton pressure-relaxation; Dyadechko & Shashkov (MOF); Baer-Nunziato / Kapila |
| `ablation-recession` | Recession as monotone solid-fraction decrease on a fixed grid (enthalpy/VOF style); B′ blowing-rate tables; conservative fraction spill | Extended B′, arXiv:2310.07080; NASA CHAR/PATO/ITRAC (NTRS 20140008557) |

### 6.3 Tables/interpolation & materials (researched 2026-07-14)

Informing FND-5 (table format/interpolation), FND-3 (geometry, prior keys in §6.2), FND-7 (materials).

**Table schema & interpolation:**
| key | Content | Source |
|---|---|---|
| `hdf5-schema` | Group-per-table; explicit coordinate-axis datasets (CF-style); chunk-aligned to stencil; mandatory attrs = provenance (tool+version, input hash, library id, seed), schema+data semantic versions, per-axis validity envelope, per-column units/uncertainty/interp-rule | CF conventions; `h5rdmtoolbox`; FAIR provenance arXiv:2604.25944 |
| `interp-monotone` | Multilinear (positivity/bound-preserving, cheapest, bit-stable) as default for conservation/positivity-critical quantities; tensor **monotone cubic (Steffen 1990 / Fritsch-Carlson PCHIP)** for C¹ needs; natural cubic/RBF overshoot — avoid, pre-resample scattered→regular offline | Steffen, *A&A* **239** (1990); Fritsch & Carlson, *SIAM JNA* **17** (1980) |
| `interp-error-budget` | Interp error is knowable (multilinear ≈ h²/8·max|∂²f|; cubic O(h⁴)); measure offline by refinement/holdout, store `interp_error_bound`, RSS into the physical uncertainty → propagates in UQ; interpolate in linearizing space (log-log XS, log-T opacity); **refuse** out-of-envelope | standard numerical analysis |
| `table-precedent` | Tabulate-offline/eval-online precedent: OpenMC union-grid lin-lin/log-log + windowed multipole; CoolProp/Cantera bicubic on (log p, log T) with cached coeffs. **Avoid** OpenMC's *stochastic* temperature interpolation (RNG) in a deterministic runtime | OpenMC methods docs; Bell et al., CoolProp, *IECR* **53** (2014) |

**Materials property data (FND-7):**
| key | Content | Source | Access |
|---|---|---|---|
| `tprc-cindas` | Gold standard T-dependent k(T), cp(T), CTE(T), emissivity + accuracy classes | Touloukian et al., *Thermophysical Properties of Matter* (TPRC, 1970s); CINDAS TPMD/ASMD | printed citable / DB paywalled |
| `nist-janaf` | cp, enthalpy, S(T), phase/vapor data → melting/vaporization | NIST-JANAF (Chase, 1998); NIST WebBook; NIST cryo fits | open |
| `mmpds` | Statistically-based A/B-basis design allowables (Inconel, stainless, Cu alloys) | MMPDS (successor to MIL-HDBK-5) | paywalled |
| `iter-mph` | Refractory/nuclear W, CuCrZr, SS316, Be incl. irradiation effects + uncertainty | ITER Material Properties Handbook; many in *J. Nucl. Mater.*/*Fusion Eng. Des.* | restricted/partly open |
| `ntp-fuels` | (U,Zr)C & carbide fuels, W-UO₂/W-UN CERMET, hot-H₂ carbide corrosion, DPA/burnup limits | Rover/NERVA LASL reports; Benensky (2018) NTP review; Zinkle & Was, *Acta Mater.* **61** (2013); Was, *Fund. of Radiation Materials Science* | NTRS/OSTI open; JNM paywalled |
| `mat-representation` | Store tabulated (T, value, ±band) + citation + validity range + state; interpolate k/cp/CTE with PCHIP; reproduce published fits exactly (NIST cryo, JANAF); phase transitions as explicit breakpoints; A/B-basis for allowables. **(v1.3: superseded as the *runtime* mechanism by the FND-7 constitutive spine — these handbook forms now serve as unit-physics *validation anchors* and cold-solid calibration data for the spine, not a per-material runtime lookup.)** | (synthesis of above) | FND-7, OFFL-5, VAL-2 |

### 6.4 Config, results-bundle, registry & reproducibility architecture (researched 2026-07-19)

Method/architecture references informing FND-4 (config schema & run manifest), FND-6 (results bundle
& reproducibility), COUP-8 (solver interface & mechanism registry). Pinned per-doc as those are
written; captured here so the research isn't lost. Reused keys (not repeated below): `moose-multiapp`
(orchestrator/Transfer pattern), `precice` (declarative composition + `precice-config-validate`),
`modelica-connector` (balanced-model / connector-balance wiring check), `gamer2-determinism` (bitwise
repro via fixed reduction order + position-keyed RNG).

**Config schema, run manifest & load-time validation (FND-4):**
| key | Content | Source | Consumers |
|---|---|---|---|
| `moose-inputparams` | Static pre-construction `validParams()`/`InputParameters` declaring params + data deps (coupled vars, material props by *semantic name*); `registerMooseObject` + Factory build objects from a `type` string; Action system expands input blocks; `deprecateParam`/`renameParam` for versioned migration | MOOSE docs (mooseframework.inl.gov); Permann et al., *SoftwareX* **11**:100430 (2020) | FND-4, COUP-8 |
| `openfoam-runtimeselection` | `declareRunTimeSelectionTable` + `addToRunTimeSelectionTable`: distributed self-registration into a string-keyed constructor table; `New(dict)` selects by name; new models in separate libs, core unedited (**caveat: HashTable/static-init order is unordered — do not replicate the ordering**) | OpenFOAM Foundation/ESI API guides (runTimeSelectionTables.H) | FND-4, COUP-8 |
| `cantera-yaml` | Fixed top-level mapping + per-entry `type` selecting a published per-parameterization subschema; import-time validation of duplicates/discontinuities | Cantera YAML Input Reference (cantera.org, v3.x) | FND-4 |
| `cargo-lock` | Intent (`Cargo.toml`, version ranges) vs machine-resolved record (`Cargo.lock`: exact versions + source + checksums, own format version for forward-compat) | doc.rust-lang.org/cargo ("Cargo.toml vs Cargo.lock") | FND-4 |
| `dvc-lock` | `dvc.lock` content-hashes every dependency, param, command **and output**; hash set = the reproducibility proof and change-detector | dvc.org (project structure / internal files) | FND-4, FND-6 |
| `conda-lock` | Fully-solved, **per-platform** lockfiles pinning every transitive dep to version + hash; a pin is valid only for the target it was solved on | github.com/conda/conda-lock; Snakemake deployment docs | FND-4 |
| `sacred-config` | Auto-emits the **final resolved** `config.json` (incl. the auto seed) + `run.json` (source, git, deps, host) as run provenance | Sacred, IDSIA (github.com/IDSIA/sacred) | FND-4, FND-6 |
| `mlflow-run` | Two-tier split: light queryable backend (params/metrics/times) vs heavy artifact store; run links to an immutable code identity (git commit / build) + logged env, seeds, data hashes | Zaharia et al., Databricks (2018); mlflow.org | FND-4, FND-6 |
| `garde-validation` | serde/toml stops at the **first** error and can't express cross-field/registry constraints → add an **error-accumulating** semantic pass (garde/validator) with source-span diagnostics (miette/ariadne) | docs.rs/garde; entropicdrift.com "Premortem vs Figment" | FND-4 |
| `k8s-crd-versioning` | Config outlives any binary → mandatory `apiVersion`; stable schemas immutable, breaking changes = new version; **incremental** conversion (v1→v2→v3), not N-to-latest | kubernetes.io (CRD versioning) | FND-4 |

**Results bundle, provenance & reproducibility (FND-6):**
| key | Content | Source | Consumers |
|---|---|---|---|
| `nexus-hdf5` | Self-describing HDF5: groups carry class attrs, one mandatory root entry, units/axes via reserved attributes → walkable without bespoke reader code | Könnecke et al., *J. Appl. Cryst.* **48**:301 (2015); nexusformat.org | FND-6 |
| `ro-crate` | Single root JSON-LD manifest describing every file + contextual entities (people/software) by typed links + content hash; a pragmatic **PROV subset** over Schema.org | Soiland-Reyes et al., *Data Science* **5**:97 (2022); researchobject.org/ro-crate | FND-6 |
| `w3c-prov` | PROV-O: Entity/Activity/Agent + `used`/`wasGeneratedBy`/`wasDerivedFrom`; the derivation subset makes "reference list = transitive closure to leaf source keys" a reachability query | Lebo, Sahoo & McGuinness, W3C Recommendation (2013) | FND-6 |
| `fair-digital-object` | Content + PID + typed metadata + operations; give each table version and result a resolvable typed handle (semver + content hash *is* the PID) so provenance edges point at stable ids | De Smedt, Lannom & Schwardmann, *Publications* **8**(2):21 (2020) | FND-6 |
| `codemeta` | Crosswalk vocabulary for software identity (`softwareVersion`, `runtimePlatform`, `buildInstructions`) → exportable Agent/software nodes for Zenodo/Software-Heritage | Boettiger, Jones et al. (2016); codemeta.github.io | FND-6 |
| `hdf5-chunking` | Compression requires chunked layout (per-chunk filters); store each field `[member, …]`, chunk on the member axis, shuffle + fast codec (zstd/gzip), ~256 KB–2 MB chunks; thin huge sweeps but keep all summaries | HDF Group guidance; parallel-HDF5 concatenation arXiv:2205.01168 | FND-6 |
| `random123-philox` | Stateless counter-based RNG: Nth value = mix(counter, key); distinct keys → independent streams, O(1) seek → maps onto `{master_seed, input_id, member_index, dimension}`; each member's draws a pure function of coordinates (no shared stream) | Salmon, Moraes, Dror & Shaw, *SC'11* (2011) | FND-6, FND-1, COUP-5 |
| `rust-fp-rfc3514` | Rust follows IEEE-754; by default **no FMA contraction** (only explicit `mul_add`) and **no reassociation** → per-build bit-exactness achievable, but contingent on target-cpu/`target-feature` (FMA) and on never inspecting NaN sign/payload bits (nondeterministic across arch) | Rust RFC 3514 "float-semantics" (rust-lang/rfcs, 2023); `-Cextra-fp-precision` | FND-6, FND-4 |
| `flit-fp-determinism` | FP results vary across compilers/flags/ISA (FMA, reassociation) — motivates recording the full FP/FMA/libm policy in the build fingerprint and testing bit-equality | Sawaya et al. (FLiT/PRUNERS), *IISWC* (2017) | FND-6 |
| `reproducible-builds` | "Same source + build environment + instructions → bit-for-bit identical artifacts"; embedded timestamps are the canonical nondeterminism (fix: `SOURCE_DATE_EPOCH`); diffoscope = recursive structural diff of divergent artifacts | reproducible-builds.org (Debian/community) | FND-6 |
| `guix-repro` | A build = pure function of its full closure (toolchain+flags+deps), content-hashed; **env pinning ≠ reproducible output** — only a deterministic toolchain + a tested byte-equality gate guarantees it | Courtès, *Nature Sci. Data* **9**:597 (2022); arXiv:2501.15919 | FND-6 |
| `esm-bitwise-repro` | Bitwise reproducibility as a continuously-tested contract; cross-platform bit-identity is generally unattainable → bit-exact within a build, statistical agreement across platforms | Baker et al. / Milroy et al., *Geosci. Model Dev.* (2015) | FND-6 |

**Solver interface, registry & wiring validation (COUP-8):**
| key | Content | Source | Consumers |
|---|---|---|---|
| `petsc-ts` | Classes as a struct of function pointers (`ops`) + string-keyed registry (`SNESRegister`/`SNESSetType`); everything selectable by an options-database string; uniform `ops` contract every impl fills | Abhyankar et al., "PETSc/TS," arXiv:1806.01437 (2018) | COUP-8 |
| `sundials-generic` | Generic base classes split **required vs optional** methods + a `GetType()` self-description so the caller checks formulation compatibility up front | Gardner et al., *ACM TOMS* **48**(3) (2022); arXiv:2011.10073 | COUP-8 |
| `ufl-fenics` | Physics (weak form) expressed as **symbolic data** over a closed, versioned vocabulary; a form compiler generates the kernel — "what" fully separated from "how" | Alnæs et al., "Unified Form Language," *ACM TOMS* **40**(2) (2014) | COUP-8 |
| `rust-registry-dispatch` | `inventory`/`linkme`/`typetag` give **no guaranteed iteration order** (life-before-main / link-order / unspecified map) — forbidden by the no-HashMap-ordering doctrine; `enum_dispatch` = static dispatch (~10× over `Box<dyn>`), fixed source order | dtolnay/{inventory,linkme,typetag}; enum_dispatch (docs.rs) | COUP-8 |
| `di-validate-on-build` | Validate the whole wiring graph at build/load time and report **all** failures together; blind spot: anything resolved dynamically (service-locator) escapes validation → forbid undeclared deps | Andrew Lock (ASP.NET Core service-provider validation); MS Learn DI | COUP-8, FND-4 |

### 6.5 Unified-grid solver, constitutive spine & predictive V&V (researched 2026-07-19/20)

Architecture references informing the **v1.3 pivot** (VISION_SCOPE v1.3; META-0 v0.5): the unified
single-grid field solver (FND-2, SOLV-1/2/3/4, COUP-3), adaptive dimensionality (FND-2 §3.4), the
constitutive spine (FND-7, OFFL-5), and predictive validation/pedigree (VAL-1, COUP-6). Reused keys
(not repeated): `gamer2-determinism` (GPU bitwise repro via fixed reduction order), `stopping-rpa-lda`
/`stopping-li-petrasso`/`stopping-bps`/`zeff-betz` (§6.1, spine stopping), `tprc-cindas`/`nist-janaf`
/`mmpds`/`iter-mph`/`ntp-fuels` (§6.3, spine unit-physics data anchors).

**Unified single-grid solver & radiation (FND-2, SOLV-1/2, COUP-3):**
| key | Content | Source | Consumers |
|---|---|---|---|
| `castro-source` | Whole system written once as `U_t = A(U) + R(U)` (advection + source); conserved state `(ρ, ρX_k, ρ𝐮, ρE, …)` — the Rule-12 state-vector algebra | Almgren et al., *ApJ* **715**:1221 (2010); Zingale et al., *JOSS* **5**:2513 (2020) | FND-2, SOLV-1 |
| `athena-ct` | Constrained transport: 𝐁 stored as face-area averages holds ∇·𝐁=0 to machine precision (structural invariant, not a runtime clean-up) | Stone et al., *ApJS* **249**:4 (2020) | FND-2, SOLV-1 |
| `radiation-m1` | M1 two-moment radiation: (E_rad, 𝐅_rad) carried *in* the conserved vector; correct thick+thin limits; single local closure; GPU-proven | Rosdahl et al., *MNRAS* **449**:4380 (2013); Wibking & Krumholz (Quokka), *MNRAS* **512**:1430 (2022) | FND-2, SOLV-2 |
| `sdc-imex` | SDC-coupled IMEX: explicit hyperbolic hydro/MHD + cell-local implicit stiff sources (reactions, radiation, conduction), iterated to 2nd order; energy-conserving, far below Lie-Trotter error at CFL | Zingale et al., arXiv:2411.12491 (2024); PeleLMeX (Esclapez et al., *JOSS* **8**:5450, 2023); IMEX rad-hydro, *JCP* (2024) | COUP-3, FND-2 |
| `stiff-reactions` | Strang/Lie-Trotter splitting **loses order and conservation** when reactions are stiff/energetic (the nuclear/antimatter end) → use SDC/IMEX there | Zingale et al., "Challenges of Modeling Astrophysical Reacting Flows," arXiv:2411.12491 (2024) | COUP-3, SOLV-4 |
| `sn-deposition` | Deterministic multigroup discrete-ordinates (Sₙ) transport swept over the grid = **one** MC-free neutron/photon deposition operator; cost = groups × ordinates | Lewis & Miller, *Computational Methods of Neutron Transport*; sweep kernels, EPJ Web Conf. (2024) | SOLV-2, OFFL-1 |

**Adaptive dimensionality (FND-2 §3.4):**
| key | Content | Source | Consumers |
|---|---|---|---|
| `spectral-azimuthal` | Azimuthal Fourier-mode truncation: keep m=0 + adaptive m≥1 per region → conservative collapse/re-expand of near-axisymmetric fields at ~2-D cost (quasi-3-D) | Lehe et al. (FBPIC), *CPC* **203**:66 (2016); Li et al. (QPAD), arXiv:2002.08494 (2020) | FND-2 |
| `symmetry-indicator` | Azimuthal-energy indicator `Σ_{m≥1}|û_m|²/(|û_0|²+ε)` (Löhner-analog, ε-normalized) + **dual-threshold hysteresis + dwell** to gate collapse/expand without thrashing | Löhner, *CMAME* **61**:323 (1987); AMReX-Astro AMR hysteresis (`n_error_buf`) | FND-2 |
| `dim-hetero-coupling` | Stitched heterogeneous-dimension (3-D↔1-D) coupling → defective BCs + spurious wave reflections, or iterative coupling, or an interface multiplier — **a seam; rejected** in favor of spectral collapse | Formaggia, Gerbeau, Nobile & Quarteroni, *CMAME* **191**:561 (2001); Boon, Nordbotten & Vatne, arXiv:1705.06876 (2018) | FND-2 |
| `anisotropic-amr` | Per-direction anisotropic refinement ratios + conservative flux-register (Berger-Colella) interfaces; static "2-D-axi = 3-D one-cell-thick azimuthal wedge + rotation BC" | AMReX (Zhang et al., *IJHPCA* **35**(6), 2021); OpenFOAM wedge (CFD Direct, 2022) | FND-2 |

**Constitutive spine (FND-7, OFFL-5):**
| key | Content | Source | Consumers |
|---|---|---|---|
| `dft-avg-atom` | Kohn-Sham DFT average-atom: EOS + transport + opacity + stopping from **one** self-consistent solve, continuous cold→plasma (no phase cliff) | Callow, Hansen, Kraisler & Cangi, arXiv:2103.09928 / *Phil. Trans. R. Soc. A* (2022); Sterne et al. (Purgatorio), *HEDP* **3**:278 (2007); Starrett et al. (Tartarus), arXiv:1804.01613 (2018) | FND-7, OFFL-5 |
| `atomec` | Open-source **Python** KS-DFT average-atom code (matches the offline stack; reference/fork target) | Callow et al., SciPy Proc. / arXiv:2206.01074 (2022); github atomec-project/atoMEC | OFFL-5 |
| `purgatorio-kg` | Electrical/thermal conductivity via Kubo-Greenwood over the **same** AA orbitals used for EOS | Starrett, *HEDP* / arXiv:1603.06925 (2016) | FND-7, OFFL-5 |
| `qeos` | Wide-range EOS: Thomas-Fermi electrons + Cowan ions + **bonding/cold-curve correction** anchoring ρ₀ & bulk modulus (smooth, no cliff) | More, Warren, Young & Zimmerman, *Phys. Fluids* **31**:3059 (1988) | FND-7, OFFL-5 |
| `feos` | FEOS/MPQEOS adds **liquid-vapor coexistence + critical point** (vapor dome) to QEOS — needed for expansion/vaporization | Faik, Tauschwitz & Iosilevskiy, *CPC* **227** (2018) | FND-7, OFFL-5 |
| `lee-more-desjarlais` | One continuous electrical/thermal conductivity cold metal → plasma; **Saha↔TF ionization blend** across the metal-insulator transition (the no-cliff template) | Desjarlais, *Contrib. Plasma Phys.* **41**:267 (2001); Lee & More, *Phys. Fluids* **27**:1273 (1984) | FND-7, OFFL-5 |
| `stanton-murillo` | Algebraic ion transport (self-diffusion, interdiffusion, viscosity, ion conductivity) weak→strong coupling from ⟨Z⟩ + screening length | Stanton & Murillo, *Phys. Rev. E* **93**:043203 (2016) | FND-7, OFFL-5 |
| `aa-opacity` | Planck/Rosseland opacity (free-free + bound-free + bound-bound) from the **same** AA orbitals | Piron et al., *Phil. Trans. R. Soc. A* **381**:20220212 (2023); Piron & Blenski, *HEDP* (2017) | FND-7, OFFL-5 |
| `tops-oplib` | Detailed-line opacity tables — **calibration/validation anchor** for AA mean opacities (AA less accurate on Planck/line-dominated) | Colgan et al. (OPLIB), *ApJ* **968**:56 (2024); TOPS (LANL) | FND-7, OFFL-5 |
| `ko-discrepancy` | Kennedy-O'Hagan GP model-discrepancy; a stationary-kernel GP **reverts to prior (zero) far from data** → smooth calibration with **no cliff**; must be physics-constrained (limits/positivity) so extrapolation stays physical | Kennedy & O'Hagan, *JRSS-B* **63**:425 (2001); Brynjarsdóttir & O'Hagan, *Inverse Problems* **30**:114007 (2014) | FND-7, OFFL-5 |

**Predictive V&V, uncertainty representation & pedigree (VAL-1, COUP-5/6, FND-1/6, META-1):**
| key | Content | Source | Consumers |
|---|---|---|---|
| `vv-oberkampf-roy` | Hierarchical validation (unit→benchmark→subsystem→system, validate bottom-up, predict top); **physics-based laws extrapolate defensibly, calibrated curve-fits do not** (the external justification for Rule 12) | Oberkampf & Roy, *V&V in Scientific Computing* (Cambridge, 2010); Oberkampf, Trucano & Hirsch, *Appl. Mech. Rev.* **57**:345 (2004) | VAL-1, META-1 |
| `asme-vv20` | Validation uncertainty `u_val = √(u_num² + u_input² + u_D²)`; model-form error bracketed as `E ± u_val`; a hard floor propagated into the composite budget | ASME V&V 20-2009 (R2016) | VAL-1 |
| `gci-roache` | Grid Convergence Index from Richardson extrapolation; **observed order checked vs formal**; Fs=1.25 (≥3 grids) / 3.0 (2 grids); out-of-asymptotic-range ⇒ flag unverified | Roache, *J. Fluids Eng.* **116**:405 (1994); Roy, *JCP* **205**:131 (2005) | VAL-1, VAL-3 |
| `area-metric` | Ferson area-metric model-form discrepancy + u-pooling; **regress vs a physical coordinate, take the 95% prediction-interval at the application condition** = a model-form band that widens with extrapolation distance | Ferson, Oberkampf & Ginzburg, *CMAME* **197**:2408 (2008); Roy & Oberkampf, *CMAME* **200**:2131 (2011) | VAL-1, COUP-6 |
| `pbox` | Probability box (interval-valued CDF): epistemic = width, aleatory = slope; nested epistemic-outer/aleatory-inner; the **single honest bound** — never collapse to one distribution | Ferson et al.; Roy & Oberkampf, *CMAME* **200** (2011); Helton et al., *RESS* **85**:39 (2004) | FND-1, FND-6, COUP-6 |
| `false-confidence` | Additive/probabilistic epistemic reps can assign high belief to false propositions + **probability dilution** (more ignorance → smaller apparent risk) → use widening intervals, not shrinking scalars | Balch, Martin & Ferson, *Proc. R. Soc. A* **475**:20180565 (2019) | COUP-6 |
| `pcmm` | Predictive Capability Maturity Model: a **vector** of maturity axes graded by rigor + independence; **report the weakest-link minimum, never a mean** (averaging incommensurable axes = the chain-strength fallacy) | Oberkampf, Pilch & Trucano, SAND2007-5948 (2007); NASA-STD-7009 (min-of-eight rule) | COUP-6 |
| `pirt` | Phenomena Identification & Ranking Table: rank each phenomenon **importance × state-of-knowledge** → focuses validation and *is* the S8 assumption register | Wilson & Boyack, *Nucl. Eng. Design* **186**:23 (1998); Boyack et al. (CSAU), *NED* **119**:1 (1990) | COUP-6, VAL-1 |
| `ect-consistency` | Ensemble-consistency test: a rerun is *statistically indistinguishable* from the accepted ensemble — the rigorous bit-repro replacement under chaos | Baker et al. (CESM-ECT), *GMD* (2015 / 2025) | META-1, VAL-3 |
| `multifidelity-uq` | Multi-fidelity UQ: cheap low-fidelity (reduced-dim) ensemble members anchored by a few high-fidelity (full-3-D) runs — how UQ stays affordable at 3-D | Peherstorfer, Willcox & Gunzburger, *SIAM Review* **60**:550 (2018) | COUP-5 |

## 7. Compute-budget data (seed from VISION_SCOPE §8)

Recorded so schedule/feasibility claims in the roadmap are sourced (the "2026 compute survey"):
OpenMC k-eig compact HEU core ~5–15 min/run to <10 pcm on ~32 cores; heating tallies to ~1%
~10 min–2 h; table sweeps ~1 day/100 pts, ~10 days/1000 pts; SANDY ~1–2 days/nuclide; Geant4
~1e6 events in minutes–hours; coupled Picard table points ~10× bare. `[RESEARCH-PENDING: re-verify
throughput on Ben's actual desktop at project start; record measured numbers here.]`
