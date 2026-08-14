# SOLV-4 — Reaction Sources

| Field | Value |
|---|---|
| **ID** | SOLV-4 |
| **Family** | SOLV (Runtime unified-grid operators) |
| **Status** | Reviewed (2026-08-14) |
| **Depends on** | FND-2, OFFL-1 (fission data), OFFL-4 (annihilation), COUP-3 (integration), SOLV-2/SOLV-3 (transport); FND-1 |
| **Version** | 0.3 (2026-08-14 review fixes: N5, N6, E-2) |

---

## 0. Purpose

SOLV-4 is the **reaction-source operator**: fission, fusion, and annihilation, unified as **one
`ReactionSource`** that emits **phase-space birth sources** `S(𝐫,E,Ω,t)` + species + strength into the shared
transport operators (charged → SOLV-3, neutrons/photons → SOLV-2). It **carries no private transport** and
**owns no time integrator** — the three reaction families differ *only* in (rate law, birth spectrum), so
there is no fusion-vs-fission-vs-annihilation seam (Rule 12).

Read after SOLV-2/SOLV-3 (the transport it feeds), COUP-3 (the SDC-IMEX step that advances its stiff terms),
and OFFL-1/OFFL-4 (the data it interpolates).

## 1. Scope & razor ruling
**Owns:** the rate laws and birth-spectrum emission for all reactions — **fission** point/few-group kinetics
(with advected precursors + source-driven subcritical), **fusion** ⟨σv⟩ product sourcing, **annihilation**
source emission. **Defers:** *all transport after birth* → SOLV-2 (neutrons/γ), SOLV-3 (charged); *time
integration* → **COUP-3** (SOLV-4 hands stiff cell-local `R(U)` terms; it owns no integrator); *data* →
OFFL-1 (fission), OFFL-4 (annihilation), META-3 (fusion coefficients); *depletion/burnup* → SOLV-8.

**Razor ruling:** the reaction is the energy origin — fully in scope — but only as a **source**; the transport
it feeds is elsewhere (VISION_SCOPE §4.1). Turbulent combustion/implosion CFD is out (equilibrium tables /
published envelopes).

## 2. Interfaces & contracts
| Interface | Consumed by | Contract |
|---|---|---|
| **Birth source** `S(𝐫,E,Ω,t)` | SOLV-2 (neutrons/γ), SOLV-3 (charged) | phase-space spectrum + species + strength; Σ birth KE = Q; carries no transport |
| **Stiff source** `R(U)` | COUP-3 | cell-local reaction/precursor terms for the SDC-IMEX implicit stage |
| **Reaction-data intake** | OFFL-1 (β_i,λ_i,Λ,ρ,spectra, **ψ/ψ† flux+adjoint shapes** — §3.2), OFFL-4 (annihilation), META-3 (⟨σv⟩) | `M`-indexed table lookups; shapes per geometry class on the same state+control axes as ρ |
| **Distribution intake** `f(E)` | SOLV-3 | fast-projectile distribution for in-flight secondary rates (§3.3) |

**Invariant:** SOLV-4 **emits and stops** — it never transports a product; energy handed to transported
particles is subtracted from the local deposit (emit-once/transport-once, COUP-2 §3.4); every reaction is one
`ReactionSource` differing only by (rate, spectrum).

## 3. Method & governing relations

### 3.1 The one rate functional
Every reaction rate is `R = n_target · ∫ σ(E) v f_projectile(E) dE`. "Thermal fusion," "beam-target," and
"in-flight secondary burn" are the **same functional** on different `f`: Maxwellian(T_i) recovers a tabulated
⟨σv⟩; a slowing-down `f` gives secondaries (§3.3); an injected beam gives driven systems. No thermal-vs-beam
branch.

### 3.2 Fission kinetics — advected precursors, source-driven subcritical
Delayed-neutron precursors are **advected scalar fields** on the grid:
- **(SOLV-4.1)** `∂C_i/∂t + ∇·(𝐮 C_i) = β_i νΣ_f φ − λ_i C_i`.
Classical **point kinetics is the static-fuel reduction** (recovered where fuel is fixed, e.g. NTP solid
core); for **flowing fuel (NSWR)** precursors drift downstream and reduce effective β — a first-order effect a
0-D PKE would miss (and a 0-D ODE amid grid fields would itself be a Rule-12 seam). **Source-driven
subcritical** (inhomogeneous) kinetics carry an external source `S(t)` for NSWR/ICAN. Reactivity feedback
ρ(T_fuel,T_mod,ρ,void,drum) is interpolated from OFFL-1 and captured **inside the SDC sweeps** (coupler 5) —
no separate feedback iteration. The stiff system is handed to **COUP-3**; SOLV-4 owns **no integrator** (a
private Rosenbrock/CRAM would be a temporal seam — and CRAM is a *depletion* method, SOLV-8, not PKE).
[META-3: `pke`, `precursor-advection`, `source-driven-pke`]

**Quasi-static closure (N5).** The neutron side factorizes **φ(𝐫,E,t) ≈ A(t)·ψ(𝐫,E)**, with the **flux shape
ψ and adjoint (importance) shape ψ† tabulated per geometry class by OFFL-1** on the same state + control axes
as ρ. The amplitude obeys the adjoint-weighted point-kinetics ODE:
- **(SOLV-4.2)** `dA/dt = [(ρ − β_eff)/Λ]·A + Σ_i λ_i c̃_i + s̃`,
with every term an **adjoint-weighted functional of the grid fields**:
- **(SOLV-4.3)** `c̃_i(t) = ⟨ψ†, C_i⟩ / ⟨ψ†, (1/v)ψ⟩` — the delayed source from the **advected precursor
  fields** (SOLV-4.1). A precursor advected into a region of low importance (ψ† → 0) contributes little; one
  advected **out of the domain** enters an explicit **adjoint-weighted loss ledger** at the outflow ports
  (COUP-2 port accounting). The **flowing-fuel β_eff reduction therefore emerges from advection + importance
  weighting** — for *any* advected-precursor configuration (flowing aqueous fuel, liquid/gas cores,
  circulating salt), not as a concept-specific feature.
- β_eff and Λ are OFFL-1's adjoint-weighted (IFP) values; s̃ = ⟨ψ†, S_ext⟩-weighted external source
  (source-driven subcritical). Static fuel recovers classical PKE identically (the weighting reduces to
  constants). **Validity:** the factorization holds while the shape is quasi-static within the tabulated
  (state × control) envelope; shape movement beyond it is flagged model-form (§5), with runtime Sₙ (SOLV-2
  general mode) as the escape hatch. [META-3: `quasi-static-kinetics`]

**Global-stiff-ODE slot (N6).** ρ, c̃_i, and s̃ are **integral functionals over core fields** — the kinetics is
one **global stiff ODE system**, not a cell-local `R(U)`. It occupies the **declared global-stiff-ODE slot in
COUP-3's SDC step** (the pulsed-event carve-out precedent; slot added to COUP-3 in this fix wave): evaluated
once per SDC sweep in fixed order (functionals = fixed-order grid reductions, FND-2 §3.7), its power feeding
the deposition sources within the same sweep.

**Kinetics-driven Δt limiter (N6):** per step, **|Δρ| ≤ `EPS_RHO_STEP` = 0.1·β_eff** and **|ΔP|/P ≤
`EPS_POWER_STEP` = 0.05** (named constants; rationale: keeps the amplitude ODE inside the SDC sweep's
contraction regime and resolves the prompt-jump timescale). Exceeding either bound halves Δt through COUP-3's
limiter interface; the limiter is part of the declared Δt rule, not an ad-hoc clamp.

### 3.3 Fusion product sourcing
⟨σv⟩ × n_i n_j → volumetric rate, emitting **birth spectra, not single energies**: Bosch-Hale D-T/D-D/D-³He
(`bosch-hale-1992`; use the 0.25% reactivity fit and 2% σ fit distinctly); p-¹¹B → a **3α continuum** peaked
~4 MeV (`pb11-reactivity`, `pb11-branch-spectrum`); D-D dual branch with **secondary** D-T/D-³He burn as the
same functional (§3.1) on SOLV-3's slowing-down `f` (interface I-SEC — the D-T ~64 keV resonance is folded
against the analytic Gaffey `f`, resolving it). Neutron peaks thermally broadened: **FWHM ≈ 177.2·√T_i**
(D-T 14.1 MeV), **82.6·√T_i** (D-D 2.45 MeV) — do not swap. Charged → SOLV-3, neutrons → SOLV-2 Sₙ.
[META-3: `bosch-hale-1992`, `pb11-reactivity`, `dt-neutron-broadening`]

### 3.4 Annihilation sourcing
Consumes OFFL-4's Geant4 source-term tables and emits the pion/charged-fragment (→ SOLV-3) and
annihilation-neutron (→ SOLV-2) birth spectra + energy partition. Carries OFFL-4's **per-quantity
physics-list spread** (not a blanket band); conversion efficiency and waste-heat are **computed downstream**
by SOLV-2/3 transport (neutrinos escape), never emitted here.

### 3.5 Determinism & halts
All rates are fixed-order `M`-indexed lookups; the global kinetics functionals (§3.2) are fixed-order grid
reductions (FND-2 §3.7); the SOLV-3↔SOLV-4 secondary cycle converges inside the SDC sweep (COUP-3). A **PKE
excursion beyond the commanded envelope** is a COUP-4 halt (loss of criticality control); a Δt driven below
the configured floor by the kinetics limiter (§3.2) halts with diagnosis rather than under-resolving a
transient.

## 4. Coupling relationships
- **COUP-3** advances the stiff `R(U)`; **SOLV-2/SOLV-3** receive birth spectra; **COUP-2** enforces
  Σ birth KE = Q + emit-once/transport-once; **OFFL-1/OFFL-4** supply data; **SOLV-3** provides `f(E)` for
  secondaries; **COUP-7** antiproton-delivery/beam boundary objects drive external sources (later wave).

## 5. Uncertainty & validity
Reactivity/⟨σv⟩ fit-level bands (stated in refs) + the annihilation **per-quantity physics-list spread** + the
point-kinetics/quasi-static model-form (flagged where the flux shape moves — fast transients, strong void
gradients). p-¹¹B net-energy caveat (`pb11-brems-limit`: brems ≥ fusion at T_e=T_i) cross-referenced to the
plasma-balance corner. Ladder: unit (analytic PKE / Bosch-Hale reproduction) → system (KRUSTY coupled shape).

**Deposition-kernel steady-state sufficiency (E-2, Ben 2026-08-14).** The precomputed deposition kernels
(SOLV-2 precomputed mode) and the ψ/ψ† shapes carry a **general control-state axis** — any declared control
coordinate, drum angle at minimum (OFFL-1 §2) — so the *commanded* power-shape dimension is captured by the
table, not frozen at one configuration. The residual is **feedback-induced shape movement at fixed control
state**. Its magnitude is measured against **OFFL-1's sparse coupled-Picard anchors**: kernel-predicted
deposition vs the coupled OpenMC↔thermal equilibrium at each anchor point, the discrepancy carried as the
kernels' model-form band. For quasi-static operation (the NTP milestone: slow throttle, shape lag second-order
behind amplitude) the anchors bound the residual — that is the sufficiency argument for the ≤10% NERVA
state-point and XE-Prime map-shape targets. Transients that move the shape faster than quasi-static are
flagged out-of-envelope; the remedy is the general mode (runtime Sₙ), never a silent extrapolation.

## 6. Validation plan
1. **Fission:** prompt-jump analytic + inverse-kinetics round-trip; reactivity-insertion vs analytic PKE;
   precursor-drift vs a flowing-fuel β_eff-reduction analytic. [META-3: `source-driven-pke`]
2. **Fusion:** reproduce Bosch-Hale ⟨σv⟩/σ(E); recover the Maxwellian limit of the rate functional;
   secondary-burn yield vs published beam-target data; neutron-peak width vs analytic.
3. **Emit-not-transport:** Σ birth KE = Q; no product advected inside SOLV-4 (COUP-2 audit).
4. **Determinism:** identical sources at 1 vs N threads.

## 7. References
META-3 keys: `pke`, `quasi-static-kinetics` *(new)*, `precursor-advection`, `source-driven-pke`, `watt-spectrum`, `bosch-hale-1992`,
`pb11-reactivity`, `pb11-branch-spectrum`, `dt-neutron-broadening`, `pb11-brems-limit`, `pbar-list-band`,
`radiation-partition`, `sdc-imex`, `stiff-reactions`. Depends on FND-2, OFFL-1, OFFL-4, COUP-3, SOLV-2, SOLV-3.

*(No open questions — SOLV-4 owns no integrator (COUP-3 does); precursors advect (point kinetics = static
reduction); reactions emit spectra not single energies; secondaries are the one rate functional on a
slowing-down `f`. Resolved 2026-07-21.)*

## 8. Change log
| Date | Version | Change |
|---|---|---|
| 2026-07-21 | 0.1 | Initial draft. One `ReactionSource` (fission/fusion/annihilation differ only in rate law + birth spectrum). Advected-precursor fission kinetics with point kinetics as the static-fuel reduction + source-driven subcritical (NSWR/ICAN); no private integrator (COUP-3 owns it; CRAM≠PKE). Fusion ⟨σv⟩ sourcing emitting birth spectra (p-¹¹B 3α continuum, D-D secondaries as the one rate functional on SOLV-3's slowing-down f, thermal-broadened neutron peaks). Annihilation via OFFL-4 with the factor-of-several band. Emit-once/transport-once; no private transport. |
| 2026-08-13 | 0.2 | R3 terminology: annihilation model-form reframed from a blanket "factor-of-several band" to OFFL-4's **per-quantity physics-list spread**; efficiency/waste-heat noted as computed downstream by SOLV-2/3 transport, not emitted here. |
| 2026-08-14 | 0.3 | **Review fix wave (N5, N6, E-2).** §3.2: quasi-static closure written explicitly — amplitude ODE (SOLV-4.2) + tabulated flux/adjoint shapes ψ/ψ† from OFFL-1 + adjoint-weighted delayed source (SOLV-4.3) and drifted-precursor loss ledger (flowing-fuel β_eff emerges from advection + importance weighting, configuration-general); kinetics declared a **class-`G` global stiff ODE** in COUP-3's SDC step + `EPS_RHO_STEP`/`EPS_POWER_STEP` Δt limiters; kernels' control-state axis noted (E-2). |
