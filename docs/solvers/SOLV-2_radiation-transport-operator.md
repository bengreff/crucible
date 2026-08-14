# SOLV-2 — Radiation-Transport Operator

| Field | Value |
|---|---|
| **ID** | SOLV-2 |
| **Family** | SOLV (Runtime unified-grid operators) |
| **Status** | Reviewed (2026-08-14) |
| **Depends on** | FND-2, FND-7 (opacity), COUP-2 (radiation-partition invariant), COUP-3 (stiff-source slot); FND-1 |
| **Version** | 0.2 (2026-08-14 review fixes: N1, N2/E-1, N3, N4) |

---

## 0. Purpose

SOLV-2 is the **one** radiation/particle-field transport operator on the unified grid: **thermal radiation
(M1 two-moment)** and **deterministic multigroup Sₙ nuclear (neutron/photon) transport** are the **same
operator over the medium-state vector `M`**, differing only by the **source spectrum** and the energy band
they act in — not by origin, not by a regime branch (VISION_SCOPE §4.2, radiation-partition invariant;
Rule 12). Reactions (SOLV-4) emit birth spectra *into* this operator; it owns everything after birth
(deposition, dose, radiative energy exchange). There is no runtime Monte Carlo.

This doc fixes the operator's governing form, its coupling into the conserved state `U`, and — the decision
that shapes the **chemical vertical slice** — **how radiation is handled when it is physically negligible
without ever becoming a seam** (§3.4, Ben 2026-07-21).

Read after FND-2 (the grid/`U` it writes into), FND-7 (the opacity it reads), and COUP-2 (the emit-once/
transport-once invariant it must satisfy).

## 1. Scope & razor ruling

**Owns:** the radiation/particle-field transport update on the grid — the **M1 two-moment** thermal band and
the **multigroup Sₙ** streaming/nuclear band, as one operator selected by **energy**, not origin; the
operator's **two solution modes** (runtime Sₙ = general mode; OFFL-1 precomputed kernels = precomputed mode —
§3.2, E-1); the mapping from reaction-emitted source spectra to deposited energy/dose; the range test
(local-deposit vs transport).

**Defers (owner named):** the **opacity / absorption-emission coefficients** over `M` → **FND-7**
(constitutive spine); **reaction birth spectra** (fission/fusion/annihilation source terms) → **SOLV-4**;
the **energetic-charged-particle** friction/slowing-down law → **SOLV-3**; **surface-to-surface view
factors** for optically-thin cavity exchange → computed at config time (`view-factor-mc`) and consumed here;
the **time-integration** of the (stiff) radiation source → **COUP-3**; the **conservation audit** →
**COUP-2**.

**Razor ruling (Reaction Razor + fidelity doctrine):** radiation is simulated **wherever it carries energy or
dose**; where it does not, it **contributes zero by physics** (small opacity → negligible emission/transport),
never by an `if(regime)` switch. The energy split (M1 vs Sₙ) is **numerical, by photon energy band** (thermal/
optically-thick vs streaming/nuclear), with **aligned band edges** — a discretization choice, **not** a seam
by origin (`radiation-partition`).

## 2. Interfaces & contracts

| Interface | Consumed by | Contract |
|---|---|---|
| **Radiation moments in `U`** | FND-2, COUP-2/3 | (E_rad, 𝐅_rad) for the M1 band are carried **in the conserved vector `U`** and transported by the same finite-volume update as hydro; Sₙ group fluxes are swept over the grid |
| **Deposition/dose write** | FND-2 (field-physics state) | volumetric heating rate + n/γ dose rates written per cell, for photons/neutrons of **any origin** |
| **Source-spectrum intake** | SOLV-4 | accepts a phase-space birth source `S(𝐫,E,Ω,t)` + species; owns all transport after birth |
| **Opacity query** | FND-7 | absorption/emission + scattering coefficients as continuous functions of `M` (never a material lookup) |

**Invariant promised to everyone:** radiation energy + matter energy is conserved to the COUP-2 tolerance
(emit-once/transport-once — a photon emitted is subtracted from the local deposit, OpenMC KERMA discipline);
the M1↔Sₙ band split conserves total radiative energy at the aligned edge; and **no result changes by a
discrete regime label** — turning a regime "on" or "off" happens through `M`-driven coefficients or a
declared config magnitude justification (§3.4), never a branch in the law.

## 3. Method & governing equations

### 3.1 Thermal band — M1 two-moment
Carry the first two angular moments (E_rad, 𝐅_rad) **inside `U`**, closed by the **M1 (maximum-entropy)
Eddington factor**, which gives the correct **optically-thick (diffusion) and free-streaming limits** with a
single local closure and is GPU-proven (Quokka). The radiation four-force couples to matter energy/momentum
as a **stiff source** handled cell-locally implicitly by COUP-3's IMEX step. [META-3: `radiation-m1`,
`sdc-imex`]

**Surface-exchange vs volumetric ownership — the τ-continuous rule (N4).** The optically-thin cavity/
surface-exchange (view-factor/radiosity) term and volumetric M1 are **one emission budget**: every emitting
cell or surface emits **once**, apportioned by a smooth function of local optical depth (Rule 12 — no
cavity-vs-participating branch). With `τ_c = κ_P(M)·L_c` (Planck-mean opacity over the config-derived
characteristic cavity path L_c, from the same geometry the view factors are cast on):
- a **wall surface** emits into the **radiosity system with weight e^(−τ_c)** (transparent limit) and into the
  **M1 field with weight 1 − e^(−τ_c)** (participating limit); absorbed radiosity flux deposits at the
  receiving surface;
- a **gas cell** always emits into M1 — its emission is ∝ κ and **self-zeroes by physics** in the transparent
  limit, so no weight is needed.
Both limits are exact: τ_c → 0 gives pure view-factor exchange (precisely the regime where M1's two-moment
closure fails for crossing beams); τ_c → ∞ gives pure M1 (the radiosity weight vanishes). The blend is a
continuous function of `M`; the blend-region residual is a declared model-form band (§5).

### 3.2 Streaming/nuclear band — multigroup Sₙ
Deterministic **discrete-ordinates (Sₙ)** multigroup transport, swept over the grid; cost = groups ×
ordinates. This is the **one MC-free neutron/photon deposition operator**; fission, fusion, activation and
annihilation differ only in the **group source spectrum** SOLV-4 injects. Multigroup cross-sections come from
OFFL-1 tables (FND-5). [META-3: `sn-deposition`]

**Two solution modes of the one operator (E-1, Ben 2026-08-14).** The nuclear band has one operator and two
ways of solving it:
- **General mode — runtime Sₙ:** the sweep specified below; valid for any geometry.
- **Precomputed mode — deposition kernels (OFFL-1):** the **same transport operator solved offline** at
  tabulated state/control points **within a declared geometry class**; runtime application is interpolation ×
  source strength. Valid **only inside the tabulated class + envelope** (FND-5 refusal outside). This is the
  milestone NTP path.
**Exactly one mode per particle-class + band per run**, config-declared and **asserted by COUP-2** (§3.4
partition invariant): the modes never co-own a particle class in a run — no double-count, no gap. The former
kernel-vs-Sₙ ambiguity was a seam-test failure (N2); this partition is its resolution. The runtime
range test (FND-2 §3.4.1) splits reaction *products* between local deposit and transport — it never selects
between modes.

**Sₙ numerics (N1):**
- **Quadrature:** **level-symmetric S₈** default (LQ_N family; 80 ordinates in 3-D), order config-selectable
  within the family. [META-3: `sn-numerics`]
- **Spatial scheme:** diamond-difference on the grid's FV cells, with the classic **set-to-zero negative-flux
  fixup** (negative angular fluxes zeroed, cell balance renormalized — conservation preserved).
- **Acceleration:** **CMFD** (coarse-mesh finite-difference) low-order acceleration of the scattering/fission
  source iteration — chosen over DSA because (a) the structured cylindrical grid supplies the coarse mesh for
  free, (b) CMFD accelerates the outer (fission/upscatter) iterations as well as within-group scattering, and
  (c) it avoids DSA's transport-consistent-discretization requirement, whose violation degrades exactly as the
  scattering ratio → 1 — graphite/H₂ cores, the case that matters. [META-3: `cmfd-acceleration`]
- **Cadence (quasi-static, change-triggered):** the flux shape is re-solved only when the **trigger metric**
  fires: max relative change since the last solve in any Sₙ input field (macroscopic group XS over `M` —
  density, temperature; control state; source *shape*) exceeding **`EPS_SN_RESOLVE` = 1%** (named constant;
  below it, group-XS movement is second-order in deposition). Between solves, deposition scales linearly with
  source amplitude.
- **Cost bound (documented):** one solve ≈ N_sweeps × N_cells × N_groups × N_ordinates updates; at S₈
  (80 ordinates), ~20 groups, ~10⁶ active cells, ≲30 CMFD-accelerated sweeps ⇒ ≲5×10¹⁰ cell-angle-group
  updates ≈ seconds-to-minutes per solve, amortized by the change trigger to ≪ the hydro cost of the
  intervening steps. A config whose trigger fires every step is mis-matched to the general mode and is
  flagged.

### 3.3 One operator, split by energy not origin
The thermal and nuclear bands are **two energy windows of one transport operator** with **aligned band
edges** so total radiative energy is conserved across the split (the HYDRA 3-package pattern). "Nuclear vs
thermal" is a **source spectrum**, not a solver; a photon from a fission γ and from a hot wall obey the same
transport, entering whichever band their energy falls in. This is the operative form of the §4.2 radiation-
partition invariant, audited by COUP-2. [META-3: `radiation-partition`]

**Band edge (N3).** Default **E_edge = 1 keV**, subject to the per-run selection rule: **E_edge ≥ 30·k_B·T_max**
(T_max the run's declared temperature-envelope maximum — the Planck tail above 30 kT carries <10⁻⁹ of thermal
emission, so no thermal photon meaningfully crosses up) **and** E_edge ≤ the lowest OFFL-1 photon group
boundary (so no nuclear group straddles the edge). A config violating either bound refuses at load. **Sub-edge
disposal:** a nuclear-band photon that down-scatters below E_edge **deposits its energy to matter at the
interaction site** (KERMA discipline — it joins matter internal energy) and re-enters radiation **only via
thermal re-emission in the M1 band**; there is **no direct cross-edge photon hand-off** (no Sₙ→M1 transfer
term). Energy crosses the edge through the matter-energy ledger only, which COUP-2 audits (§6.4).

### 3.4 Negligible-but-not-seamed: the chemical slice *(Ben, 2026-07-21)*
For an RL10-class LH₂/LOX engine, gas radiation transports **~0.1–0.3% of enthalpy** (≈11% of *wall* flux,
which is ~1–3% of enthalpy and **near-fully recuperated** in an expander cycle) — ≈10× below the 2% Isp bar —
and **gray/WSGG models run >100% wrong** vs line-by-line for H₂O/CO₂ (`rocket-gas-radiation`). So radiation
must be *absent from the chemical stage-1 performance path* **without** creating a seam. The rule:

- **The operator is always architecturally present.** Its strength is set by the **spine opacity over `M`**
  (§3.1). For combustion-gas `M` the opacity is small, so the operator **self-zeroes by physics** — the
  correct Rule-12 behavior (terms vanish by physics, not by a branch, META-1 Rule 12).
- **Omission is a declared magnitude justification, not a code branch.** A stage-1 chemical config may leave
  the operator **un-instantiated** as an **optimization**, but only with a **recorded PIRT entry** (COUP-6):
  *radiative enthalpy transport <0.3%, below the model-form band; excluded by magnitude.* This is a config-
  level, envelope-checked, revisitable decision (FND-1 §3.7) — **there is no `if(chemical)` in the law.**
- **Nothing the reaction touches is magic-boundaried.** Wall heating, conduction, combustion and jacket
  pickup remain full grid physics; only the *radiative transport channel* — genuinely negligible here — is
  the omitted term, and it is omitted honestly (declared, magnitude-justified), which is exactly the
  "not excluding something arbitrarily" test.
- **Retained where it matters:** the ~11%-of-wall-flux radiative term **is kept for stage-2 wall/lifetime**
  (liner/injector-face temperature) via this **same operator** (a surface-emission/radiosity term on the
  wall-thermal model), so radiation is excluded only from the stage where it is negligible, not from the
  physics. This respects the stage-1 (function/Isp) vs stage-2 (lifetime) split (VISION_SCOPE §7.6).

The general doctrine SOLV-2 encodes: **include a radiation operator wherever it carries power (nuclear/plasma/
pulsed regimes); elsewhere let it self-zero, and if omitted for cost, declare the magnitude.** A bad model of
a 0.2% effect is worse than a declared omission.

## 4. Coupling relationships
- **COUP-3** advances the M1 radiation source as a cell-local implicit (stiff) term in the IMEX step; the
  transport of (E_rad,𝐅_rad) rides the same flux update as hydro (FND-2 §3.4).
- **COUP-2** audits radiation↔matter energy exchange every step and enforces the emit-once/transport-once
  radiation-partition invariant; a violation halts.
- **FND-7** supplies opacity over `M` (the self-zeroing knob, §3.4); **SOLV-4** injects birth spectra;
  **SOLV-3** owns charged-particle slowing (a distinct operator, coupled by shared deposition).
- **FND-2** stores the deposited heating/dose per cell; **SOLV-8** (stage-2) consumes wall radiative loading
  for lifetime.

## 5. Uncertainty & validity
SOLV-2 introduces **transport/closure model-form uncertainty**: the M1 closure error (anisotropic streaming),
Sₙ angular/group discretization error (incl. ray effects, §6.7), and — where included — opacity band error
from FND-7. For the chemical slice the **declared magnitude omission** (§3.4) contributes a bounded,
PIRT-recorded ≤0.3%-of-enthalpy epistemic interval to the p-box (COUP-6), not a hidden zero.

**Transparent-cavity PIRT entry (N4):** M1's known model-form failure is **crossing beams in optically-thin
media** — two hot walls exchanging across transparent gas collapse into a single flux direction. That is
exactly the **NTP flow-channel case** (2500–2800 K walls across optically-thin H₂), which the milestone
exercises. Recorded as a standing PIRT entry (importance high, knowledge medium); the τ-blend (§3.1) routes
that regime to the radiosity term, and the residual blend-region model-form band is carried per-quantity in
the p-box.

Validation-ladder rungs: analytic (Su-Olson/Marshak, Sod-radiation), reference-code (Quokka), then hardware
where dose/heating anchors exist. Modules that stay at rung (ii) are pedigree-flagged.

## 6. Validation plan
1. **M1 analytic limits:** recovers the diffusion and free-streaming limits; **Su-Olson non-equilibrium
   Marshak wave** reproduced within ~1–3% (T_rad, T_mat). [META-3: `su-olson`]
2. **Radiating shock / advecting-diffusion:** matches Quokka reference profiles to a few % (differential
   oracle). [META-3: `radiation-m1`]
3. **Sₙ deposition:** multigroup deposition reproduces an analytic/attenuation benchmark and cross-checks an
   OpenMC multigroup case. [META-3: `sn-deposition`]
4. **Band-split conservation (extended, N3):** total radiative + matter energy is conserved across the M1↔Sₙ
   edge **including sub-edge down-scatter**: a hard-γ slab source thermalizing through matter deposit and M1
   re-emission closes the ledger to the COUP-2 audit tolerance, with zero direct Sₙ→M1 photon transfer.
5. **No-seam / self-zero test:** sweeping `M` from combustion-gas to optically-thick, the operator's
   contribution varies **continuously** (no jump), and its combustion-gas contribution is ≤0.3% of enthalpy —
   confirming §3.4's omission is magnitude-justified, not arbitrary.
6. **Emit-once/transport-once:** energy transported out of a cell equals energy subtracted from its deposit.
7. **Ray effects / streaming (N1):** a localized source in a transparent region exhibits the expected S₈
   ray-effect artifacts; quantified against a refined-quadrature (S₁₆) reference and recorded as the angular-
   discretization band (§5).
8. **Two-beam / cavity exchange (N4):** a concentric-cylinder cavity (hot wall / cold wall across transparent
   gas) reproduces the analytic view-factor solution through the τ-blend; sweeping κ from transparent to
   optically thick, the computed exchange transitions **continuously** between the radiosity and M1 answers —
   no jump at any τ.
9. **Mode-partition assertion (E-1):** a config declaring both modes for one particle-class+band refuses at
   load (COUP-2 §3.4); kernel-mode deposition matches general-mode Sₙ on a tabulated-class case to the
   kernel's declared interpolation band.

## 7. References
META-3 keys: `radiation-m1`, `sn-deposition`, `sn-numerics` *(new)*, `cmfd-acceleration` *(new)*,
`radiation-partition`, `su-olson`, `rocket-gas-radiation`, `view-factor-mc`, `sdc-imex`. Depends on FND-2
(grid/`U`), FND-7 (opacity over `M`), COUP-2 (partition invariant + audit + mode assertion), COUP-3
(stiff-source IMEX), SOLV-4 (birth spectra), OFFL-1 (MGXS + precomputed-mode kernels).

*(No open questions — chemical-slice radiation resolved by Ben, 2026-07-21: one unified operator, self-zeroing
by physics, stage-1 omission only by declared magnitude justification, ~11%-of-wall-flux term retained for
stage-2; no regime branch, no magic boundary.)*

## 8. Change log
| Date | Version | Change |
|---|---|---|
| 2026-08-14 | 0.2 | Review fixes. **N2/E-1:** precomputed kernels = the **precomputed solution mode of the same operator** (valid only within a tabulated geometry class); runtime Sₙ = general mode; exactly one mode per particle-class+band per run, config-declared, COUP-2-asserted (§3.2, §6.9). **N1:** Sₙ numerics specified — level-symmetric S₈ quadrature, diamond-difference + set-to-zero fixup, **CMFD** acceleration (rationale stated), quasi-static change-triggered cadence (`EPS_SN_RESOLVE` = 1%), documented cost bound; ray-effect item added (§6.7). **N3:** band edge pinned (default 1 keV; rule E_edge ≥ 30·k_B·T_max ∧ ≤ lowest photon group boundary); sub-edge down-scatter deposits to matter (KERMA discipline), re-enters only via M1 thermal re-emission — no cross-edge hand-off; conservation test extended (§6.4). **N4:** τ-continuous emit-once ownership rule between radiosity and M1 (§3.1); transparent-cavity (NTP-channel) PIRT entry (§5); two-beam/cavity benchmark (§6.8). |
| 2026-07-21 | 0.1 | Initial draft. One radiation operator = M1 two-moment (thermal, moments in `U`) + multigroup Sₙ (nuclear), split by **energy band not origin** with aligned edges (radiation-partition invariant); source-spectrum intake from SOLV-4; opacity from FND-7 over `M`. **§3.4 chemical-slice ruling (Ben 2026-07-21):** self-zero by physics, stage-1 omission only as a declared PIRT magnitude justification (radiative enthalpy <0.3%, below model-form band), ~11%-of-wall-flux term retained for stage-2 wall/lifetime — no `if(regime)` branch, no magic boundary. Validation via Su-Olson/Marshak + Quokka oracle. |
