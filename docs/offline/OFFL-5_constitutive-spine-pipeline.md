# OFFL-5 — Constitutive-Spine Pipeline

| Field | Value |
|---|---|
| **ID** | OFFL-5 |
| **Family** | OFFL (Offline pipeline) |
| **Status** | Reviewed (2026-08-14) |
| **Depends on** | FND-7 (the spec it generates), FND-5 (table schema), FND-1 |
| **Version** | 0.3 (2026-08-20: §3.1a the **chemical-regime transport stage** — the first stage of the slot to ship; landed with plan S4) |

---

## 0. Purpose

OFFL-5 is the **offline pipeline that generates the constitutive spine** (FND-7): continuous EOS + transport +
opacity + stopping over the medium-state vector `M`, for arbitrary materials cold solid → hot plasma, into
FND-5 HDF5 tables keyed on **state**, not `material_id`. The standing scope decision (Ben, 2026-07-21):
**this is a table-generation pipeline, never a runtime DFT solver.** The staging is **superseded (Ben,
2026-08-14, D-E/S15):** the W2/W3 backbone is the **Saha/QEOS analytic model over the full (ρ,T) range from
day one** — smooth everywhere, honest ~10–20% bands — with reference data **GP-blended on top where data
exist** (never a data-first corner whose envelope cliffs at the data edge); the DFT average-atom backbone is
the later-wave upgrade of the same slot (the WDM valley).

Read after FND-7 (the spec/invariants) and FND-5 (the table contract).

## 1. Scope & razor ruling
**Owns:** the offline generation of the spine tables — the W2/W3 Saha/QEOS analytic backbone with its
data GP-blend (D-E), the DFT-average-atom backbone upgrade (later waves), the closure blends, and the
Gaussian-process no-cliff calibration.
**Defers:** the spine **spec/invariants** (continuity, no-cliff mandate, per-regime bands) → **FND-7**; the
table schema/interpolation → **FND-5**; the runtime use → **SOLV-1/2/3** via `M`; the uncertainty type →
**FND-1**.

**Razor ruling:** infrastructure, load-bearing for S8 (validation-by-parts of each constitutive law). The
fidelity obligation is a **continuous** law over all of `M` (no phase/regime cliff) with a **declared** band.

## 2. Interfaces & contracts
| Interface | Consumed by | Contract |
|---|---|---|
| **Spine tables** | FND-5 → SOLV-1/2/3 (via `M`) | EOS, transport, opacity, stopping vs the §3.2 coordinates (per-species Helmholtz F_s(ρ_s,T); transport/opacity vs (ρ,T,⟨Z⟩)); declared bands (COUP-5 declarations); continuous across phase/regime |
| **Material definition** | config (FND-4) | fundamental data (composition, cohesive energy, ρ₀, bulk modulus, ionization potentials) + optional GP correction — **not** a handbook table |
| **Transport feed (input, S23)** | OFFL-5 ← OFFL-3 | Cantera mixture `μ, k, c_p, c_v, ∂h/∂Z\|_{p,T}` vs (T, p, Z) — the chemical regime's backbone evaluation (§3.1a) and, in the plasma/WDM regimes, GP calibration data for the Saha/QEOS assembly. Either way the spine is the **sole runtime provider** of transport |
| **Chemical-regime transport surface** *(0.3, plan S4)* | FND-5 → SOLV-1 (via `M`) | μ, k, c_p, c_v, ∂h/∂Z vs the **local state (p, h, Z)**; axes/envelope matched to the equilibrium surface's; declared **10–20%** band; measured interp-error bounds (§3.1a) |

**Invariant:** every generated property is a **continuous** function of `M` (a discontinuity across melt/
ionization/metal-insulator is a bug, FND-7); it carries a **declared band**; the GP correction **reverts to
the physics model away from data** (no calibration cliff).

## 3. Method

### 3.1 The W2/W3 backbone — Saha/QEOS analytic spine + data GP-blend *(D-E, S15)*
The milestone spine is the **Saha/QEOS analytic backbone over the full (ρ,T) range from day one**: QEOS
(Thomas-Fermi electrons + Cowan ions) with the FND-7 §3.2 cold-curve/vapor-dome anchoring, Saha ionization
for ⟨Z⟩, Lee-More-Desjarlais conductivity, Stanton-Murillo ion transport, Kramers-class analytic opacity —
smooth everywhere, honest **~10–20% bands** (FND-7 §5). **Tabulated experimental data are GP-blended on top
where data exist** (§3.3): **Cantera** hot-combustion-gas properties (H₂O/H₂/OH/O₂) — including OFFL-3's
Cantera mixture-transport feed (S23: μ, k enter *here* as calibration data; the spine is the **sole runtime
provider** of transport); **CoolProp/REFPROP** cryogens (parahydrogen, oxygen); **TPRC/CINDAS + NIST** cold
structural metals (Cu alloys, stainless). The data are **calibration on an everywhere-defined backbone,
never the backbone itself** — coverage-weighted, reverting to the analytic model outside data (the no-cliff
guarantee): the spine's envelope is the model's, not the data's, so an off-data query (the NSWR
ionizing-steam stretch) gets the analytic backbone with its declared band, not a refusal cliff. These are
the data anchors FND-7 §3.6 names. [META-3: `saha-qeos`, `qeos`, `cantera`, `coolprop`, `tprc-cindas`,
`nist-janaf`, `lee-more-desjarlais`, `stanton-murillo`]

### 3.1a The chemical-regime transport stage *(0.3, plan S4 — the stage that ships first)*
The same transport slot, in the regime the chemical sandbox marches in. Backbone: **mixture-averaged
Chapman-Enskog kinetic theory** (FND-7 §3.3) evaluated by **Cantera** on the pinned species set and its LJ
transport data — the chemical sibling of Saha/QEOS, and, exactly like it, **a physics model defined over the
whole state space, not a data corner**. In this regime the §3.3 GP discrepancy is therefore **identity**: the
backbone *is* the evaluation, there is no independent measured set to nudge it toward, and pretending
otherwise would fabricate a correction. That is a stage property, recorded, not a weakening of §3.3 — when
high-temperature H₂O/H₂ transport measurements are brought in, they enter as a coverage-weighted GP on top of
this backbone with no change to the runtime contract.

**Emitted product:** one FND-5 `regular` table, keyed on the **local state (p, h, Z)** (FND-7 §3.7 — the S22
rule; the Cantera evaluation happens in its natural (T, p, Z) coordinate at each grid node and is written out
in the runtime coordinate), carrying **μ, k, c_p, c_v, ∂h/∂Z|_{p,T}** with per-column measured
interpolation-error bounds (absolute + rule-space, FND-5 §3.4) and the **declared 10–20% band**. Its axes and
envelope are generated to **match the equilibrium surface's**, so the two runtime surfaces a cell interrogates
refuse and accept on the same set of states — an envelope gap between them would be a mid-march surprise the
COUP-8 §3.3(2) coverage check cannot see across two independently declared tables. Provenance stamps the
Cantera version, the mechanism/transport-data deck, **and the transport model** (mixture-averaged vs
multicomponent) — the model form is part of the input-deck hash, per the §3.4 lesson that two different
physics choices must never share a deck hash. [META-3: `chapman-enskog-mixavg`, `h2o2-transport-data`,
`cantera`]

### 3.2 The DFT average-atom backbone upgrade *(later waves)* — and the emitted coordinates *(S12)*
The WDM-valley upgrade of the same slot: a **Kohn-Sham DFT average-atom** solve (prototyped from open-source
**atoMEC**) from which EOS (free-energy derivatives), transport (Kubo-Greenwood on the same orbitals), and
opacity (same orbitals) are all derived consistently; **QEOS/FEOS** cold-curve + vapor-dome anchoring;
**Lee-More-Desjarlais** conductivity; **Stanton-Murillo** ion transport; **RPA-LDA + Li-Petrasso + BPS**
stopping. All run **offline** into tables — no AA solve at runtime. [META-3: `dft-avg-atom`, `atomec`,
`qeos`, `feos`, `lee-more-desjarlais`, `stanton-murillo`, `aa-opacity`, `stopping-rpa-lda`]

**Emitted table coordinates (S12, every stage — Saha/QEOS and DFT-AA alike):** per-species Helmholtz
**F_s(ρ_s, T)** with tabulated derivatives, plus the FND-5 **`thermo_audit`** convexity/sound-speed-positivity
audit (S11 — a generation obligation of this pipeline); transport and opacity vs **(ρ, T, ⟨Z⟩)**;
stopping/collision bundles vs (E, projectile species; ρ, T, ⟨Z⟩). The runtime mixture combination (isobaric
additive-volume pressure equilibration, fixed iteration count) is **FND-7 §3.7's** — this pipeline emits
per-species sub-EOS in exactly those coordinates and never pre-mixes engine-specific compositions
(capabilities are general — the tool doctrine). [META-3: `helmholtz-table`]

### 3.3 The no-cliff GP calibration
Where data exist, the physics model is **smoothly nudged toward it** by a **coverage-weighted Kennedy-O'Hagan
GP discrepancy** in `M`-coordinates: **multiplicative** for positive-definite quantities (κ, opacity, p),
**additive** for signed/energy quantities; physics-constrained (log-space, monotonicity/thermo limits). A
stationary-kernel GP **reverts to zero correction away from data** — the mathematical guarantee against a
calibration cliff — and its **posterior variance is the data-coverage-weighted band**. (Caveat internalized:
K&O is identifiable in *prediction* mode, which is how it is used here.) [META-3: `ko-discrepancy`]

### 3.4 Uncertainty & regime bands
The pipeline emits FND-7's declared per-regime bands: cold solid (EOS good if cold-curve-anchored; opacity
weakest), WDM (the hard valley — EOS few–15%, transport 20–50%, opacity ~2×; ~10–20% overall at the
Saha/QEOS stage), hot plasma (converges to Spitzer/Saha/Kramers — tightest). Each is typed (FND-1); the
FND-5 interpolation bound is **declared table metadata consumed by COUP-5 as an epistemic interval** (FND-5
§3.4) — never RSS'd with the physical band (S14).

### 3.5 Compute cost & table size *(S24)*
Cost = state points × per-solve time × elements. A per-element (ρ,T) grid of ~50×50 ≈ 2.5×10³ points at
~10–60 s per AA solve (atoMEC KS solve + Kubo-Greenwood response) is ~7–40 CPU-h per element; ~10–20
elements per wave ⇒ ~10²–10³ CPU-h ⇒ **hours-to-days wall on 32 cores** — comfortably inside the
VISION_SCOPE §8 CPU budget, and regenerable dozens of times. The Saha/QEOS analytic stage is orders cheaper
(ms per point). Table size is trivial: 2.5×10³ points × O(10) quantities × f64 (+ derivative datasets for
the Helmholtz tables) ≈ **MB-class per material** — negligible against the nuclear tables.

## 4. Coupling relationships
- **FND-7** is the spec OFFL-5 realizes; **FND-5** stores/interpolates the output and enforces the (wide)
  envelope; **SOLV-1/2/3** consume EOS/transport, opacity, stopping via `M`; **VAL-2** validates at the
  unit-physics tier (Hugoniot, transport, opacity, stopping); **COUP-6** scores data-anchored regions higher.

## 5. Uncertainty & validity
Originates the spine's **model-form band** (regime-dependent, §3.4) — often the dominant spread for advanced
concepts (Sobol will say so). Validity is **as wide as the model supports** (the spine's raison d'être);
out-of-envelope still refuses/flags (FND-5). Ladder: unit-physics tier (VAL-1) — the constituent-law
validation that underwrites S8's extrapolation claim.

## 6. Validation plan
1. **Data-anchored regions:** hot-gas transport reproduces the Cantera reference (incl. the OFFL-3 feed);
   cryo props match REFPROP; metal `k(T)` matches CINDAS/NIST within class accuracy — and off-data the
   backbone reverts to Saha/QEOS within its declared ~10–20% band (no data-edge cliff, D-E).
1b. **Thermo audit (S11):** every Helmholtz table passes the convexity/sound-speed-positivity audit over its
   declared envelope before shipping (FND-5 `thermo_audit`).
2. **EOS:** principal-Hugoniot reproduced; cold-curve matches ρ₀/bulk modulus; no kink across pressure
   ionization. [META-3: `hugoniot-anchor`]
3. **No-cliff:** every property + first derivative continuous across melt/vaporization/ionization/MIT (the
   Rule-12 seam test applied to the spine).
4. **GP calibration:** correction reverts to the physics model outside data coverage; posterior variance
   tracks data density.
5. **Provenance:** every table regenerable from generator hash + tool/data versions.

## 7. References
META-3 keys: `saha-qeos`, `cantera`, `chapman-enskog-mixavg`, `h2o2-transport-data`, `coolprop`,
`tprc-cindas`, `nist-janaf`, `dft-avg-atom`, `atomec`,
`qeos`, `feos`, `helmholtz-table`, `lee-more-desjarlais`, `stanton-murillo`, `aa-opacity`,
`stopping-rpa-lda`, `ko-discrepancy`, `hugoniot-anchor`. Depends on FND-7 (spec), FND-5 (schema), FND-1
(uncertainty); consumes OFFL-3's transport feed (S23).

*(No open questions — W2/W3 backbone = Saha/QEOS analytic over the full (ρ,T) range + coverage-weighted GP
data blend, Ben 2026-08-14 (D-E), superseding the 2026-07-21 data-first corner; DFT-AA backbone remains the
later-wave WDM upgrade.)*

## 8. Change log
| Date | Version | Change |
|---|---|---|
| 2026-08-20 | 0.3 | **The chemical-regime transport stage (plan S4 — landed with the code).** New §3.1a: the transport slot's chemical stage, mixture-averaged Chapman-Enskog via Cantera on the pinned species/LJ deck, **the first stage of the slot to ship**. Records that in this regime the §3.3 GP discrepancy is **identity** — the backbone *is* the evaluation and there is no independent measured set to nudge toward; a fabricated correction would be worse than none, and measured high-T transport data enter later with no runtime-contract change. Emitted product: one FND-5 `regular` table on the **local state (p, h, Z)** (FND-7 §3.7/S22; Cantera evaluates in (T, p, Z) at each node and the table is written in the runtime coordinate) carrying **μ, k, c_p, c_v, ∂h/∂Z\|_{p,T}**, axes/envelope **matched to the equilibrium surface's** (a gap between two independently declared surfaces is invisible to the COUP-8 §3.3(2) coverage check), measured abs + rule-space interp bounds, declared **10–20%** band, and provenance stamping the transport model into the deck hash. §2 gains the emitted-surface row and widens the OFFL-3 feed's quantity set to the caloric companions. |
| 2026-08-14 | 0.2 | **Post-review fix wave (S12, S14, S15/D-E, S23, S24).** §0/§3.1: W2/W3 deliverable re-staged per D-E — **Saha/QEOS analytic backbone over the full (ρ,T) range from day one** (~10–20% bands; QEOS+Saha ⟨Z⟩ feeding LMD/Stanton-Murillo/Kramers-class closures) with reference data **GP-blended on top** (data = calibration on an everywhere-defined backbone, never the backbone — no data-edge cliff; NSWR ionizing-steam covered); DFT-AA = later-wave upgrade of the same slot. §3.2: **emitted table coordinates stated** — per-species Helmholtz F_s(ρ_s,T) + `thermo_audit` obligation (S11), transport/opacity vs (ρ,T,⟨Z⟩); mixture rule cross-ref'd to FND-7 §3.7, never pre-mixed (S12). §3.4: interpolation bound = declared metadata → COUP-5 epistemic interval, never RSS'd (S14). New §3.5: **compute-cost/table-size budget** — points × per-solve × elements ⇒ hours-to-days on 32 cores, MB-class tables (S24). §2: OFFL-3 Cantera **transport feed named as input**; spine = sole runtime transport provider (S23). |
| 2026-07-21 | 0.1 | Initial draft. Table-generation pipeline (never a runtime AA solver). W2 chemical/cold corner built data-first (Cantera + CoolProp/REFPROP + TPRC/CINDAS + NIST); DFT-average-atom backbone (atoMEC/QEOS/FEOS/LMD/Stanton-Murillo/RPA-LDA) deferred to WDM/plasma waves; coverage-weighted Kennedy-O'Hagan GP no-cliff calibration; per-regime declared bands; unit-physics validation tier. |
