# SOLV-4 — Reaction Sources

| Field | Value |
|---|---|
| **ID** | SOLV-4 |
| **Family** | SOLV (Runtime unified-grid operators) |
| **Status** | Reviewed (2026-08-14) |
| **Depends on** | FND-2, OFFL-1 (fission data), OFFL-4 (annihilation), COUP-3 (integration), SOLV-2/SOLV-3 (transport); FND-1 |
| **Version** | 0.4 (2026-08-19: §3.6 chemical burn-progress source — VISION_SCOPE v1.5; prior: 0.3 review fixes N5, N6, E-2) |

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
source emission, and the **chemical burn-progress source** (§3.6, v0.4 — the rate law that ignites,
propagates, and extinguishes tabulated-chemistry combustion; the same pattern: a rate law + a "spectrum"
that is heat/species release into `U`, no private transport). **Defers:** *all transport after birth* → SOLV-2 (neutrons/γ), SOLV-3 (charged); *time
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
| **Burn-progress source** `dc/dt` (§3.6) | SOLV-1 (advects c, evaluates the blended EOS), COUP-4 (`NEVER_IGNITED`/`FLAMEOUT`) | the one continuous chemical rate law: flame propagation + auto-ignition over local `M`; closures from OFFL-3 with {citation, envelope, band} |

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

### 3.6 Chemical burn-progress source *(v0.4, VISION_SCOPE v1.5 — ignition as computed physics)*

Equilibrium chemistry cannot express ignition: mixed propellants are "already burnt" by construction, so
neither start-up, failure-to-light, nor flameout exists in it — while the frozen mode never burns at all.
The chemical regime therefore carries a **burn-progress field c ∈ [0,1]** on `U` (SOLV-1 §3.4): the cell's
**burnt mass fraction** — sub-cell combustion state represented statistically, per the fidelity contract
(`PLAN_CHEMICAL_SANDBOX.md` §2.2). Thermochemistry is the c-blend of two OFFL-3 branches: **unburnt**
(frozen reactant mixture, valid to cryogenic temperatures) and **burnt** (the shifting-equilibrium surface).
**Sub-cell partition (the blend rule):** the two sub-states share the cell pressure `p`; the cell's
specific enthalpy splits mass-weighted, `h = (1−c)·h_u + c·h_b`. The **unburnt enthalpy closure** is the
standard premixed-SGS assumption, declared: unburnt gas rides its injection-state isentrope to the local
pressure (`h_u = h_u(p; h_inj, Z)` from the unburnt-reactant surface), giving `T_u(p, Z)` — the coordinate
`S_L` and `τ_ign` are keyed on — without carrying a second energy field; `h_b` is then `(h − (1−c)h_u)/c`
and the burnt branch is interrogated/projected at `(p, h_b, Z)` exactly as in pure shifting mode. `c = 1`
recovers shifting mode identically; `c = 0` recovers the frozen reactant branch. The closure is data with a
band (it neglects unburnt preheat by the flame — bounded at build against the offline 1-D flame solution),
never a branch.

**The one continuous rate law (Rule 12 — no ignition `if`), TFC form:**
- **(SOLV-4.4)** `∂(ρc)/∂t|_source = ρ_u·S_T·|∇c| + ρ·(1−c)/τ_ign(p, T_u, Z)` — the standard
  turbulent-flame-closure **propagation term** (unburnt density × turbulent flame speed × progress-gradient
  magnitude; `S_T` = laminar `S_L(p, T_u, Z)` from the OFFL-3 surface × the wrinkling factor,
  `turbulent-flame-speed`) plus an **auto-ignition term** (induction-time surface `τ_ign`, Arrhenius-class
  fits). **The grid-independence mechanism is the TFC pairing, stated explicitly:** the source is
  accompanied by a **matched front-thickening diffusion of c** (`D_c` sized so the front spans a fixed
  Θ ≈ 2–4 cells at every resolution; the KPP front-speed identity of the paired diffusion–source system
  holds the propagation speed at `S_T`, not at whatever numerical diffusion provides) — a bare source with
  no matched diffusion self-sharpens to the grid scale and its front speed becomes grid-set, which is
  exactly what the `flame_1d` gate exists to catch. `D_c` is the front-carrier device (declared, like SRD),
  not a physics claim about flame thickness. **On the two terms adding:** in the deflagration regime
  `τ_ign` is orders longer than the flame's passage time, so the sum is propagation-dominated (no
  double-count in practice); where compressed end-gas approaches autoignitive states the added induction
  consumption **is** the physics (pre-ignition/knock class), and the two closures are anchored separately
  (`flame_1d` vs the ignition-delay reproduction) so each is right where it dominates.
  Both closures carry **{citation, envelope, band}** (META-3 `h2-flame-speed`, `h2-ignition-delay`,
  `turbulent-flame-speed`, `h2-kinetics-mech`) and blend continuously over `M`: outside flammability /
  below quench conditions, `S_T → 0` and `τ_ign → ∞` **by the closure surfaces themselves**, never by a
  regime branch. Heat release is *implicit in the blend* — as c rises, the cell's EOS reads the burnt branch
  (the S18 η_c\* knockdown composes unchanged: it offsets the burnt branch's interrogation coordinate).
- **The igniter is a boundary object** (COUP-7 §3.3): a declared, scheduled, localized energy deposit
  (position, radius, duration, energy). It is not special-cased physics — it merely creates a state whose
  induction term fires. Spark → kernel → propagating burn front → pressure rise → downstream compression →
  (possibly) auto-ignition ahead of the front: **events cause events**, every link a field evolution.
- **Grid-independence invariant (THE acceptance test):** the front's propagation *speed* is closure-set,
  never grid-set — the `flame_1d` fixture must reproduce the same front speed on coarse and fine grids
  (refinement sharpens *where*, never changes *what*). A model whose burn timeline moves with resolution
  fails this gate.
- **Halt thresholds (single owner):** the run's **reacting measure** `R(t) = ∫ c(1−c)·rate dV` (the live
  flame content). `NEVER_IGNITED` = igniter schedule exhausted and `R` never exceeded `EPS_IGNITED` (named
  constant, fixed at build); `FLAMEOUT` = `R` falls below `EPS_IGNITED` after having exceeded it, before the
  dwell completes (COUP-4 §3.2 consumes both).
- **Integration & determinism:** `dc/dt` is a genuinely **cell-local class-`R` stiff source** (COUP-3 §3.1)
  — fixed-order `M`-indexed lookups, no private integrator, bit-reproducible. **Runtime finite-rate
  networks remain out of scope**: `S_L` and `τ_ign` are *offline-computed* surfaces (Cantera 1-D flames +
  0-D reactors on a cited mechanism, OFFL-3 §3.3), validated at the unit tier against measured flame-speed
  and shock-tube ignition-delay data (VAL-2 §3.3).
- **What c is not:** not a species (elements still advect exactly; Z is untouched), not an efficiency knob
  (η_c\* stays S18's separate, declared knockdown), and not resolved flame structure (a real H₂/O₂ front is
  ~10 μm; the model's front is the closure speed carried on grid scales — declared model form, banded).

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
5. **Burn progress (§3.6):** `flame_1d` — front speed matches the closure and is **grid-independent**
   (coarse vs fine, THE gate); `spark_box` lights; `lean_no_light` refuses (NEVER_IGNITED); `quench_box`
   extinguishes (FLAMEOUT); τ_ign surface reproduces shock-tube ignition delays and S_L the measured flame
   speeds within declared bands (VAL-2 §3.3 anchors).

## 7. References
META-3 keys: `pke`, `quasi-static-kinetics` *(new)*, `precursor-advection`, `source-driven-pke`, `watt-spectrum`, `bosch-hale-1992`,
`pb11-reactivity`, `pb11-branch-spectrum`, `dt-neutron-broadening`, `pb11-brems-limit`, `pbar-list-band`,
`radiation-partition`, `sdc-imex`, `stiff-reactions`; §3.6: `progress-variable`, `h2-flame-speed`,
`h2-ignition-delay`, `turbulent-flame-speed`, `h2-kinetics-mech`. Depends on FND-2, OFFL-1, OFFL-3 (§3.6
closures), OFFL-4, COUP-3, SOLV-2, SOLV-3.

*(No open questions — SOLV-4 owns no integrator (COUP-3 does); precursors advect (point kinetics = static
reduction); reactions emit spectra not single energies; secondaries are the one rate functional on a
slowing-down `f`. Resolved 2026-07-21.)*

## 8. Change log
| Date | Version | Change |
|---|---|---|
| 2026-07-21 | 0.1 | Initial draft. One `ReactionSource` (fission/fusion/annihilation differ only in rate law + birth spectrum). Advected-precursor fission kinetics with point kinetics as the static-fuel reduction + source-driven subcritical (NSWR/ICAN); no private integrator (COUP-3 owns it; CRAM≠PKE). Fusion ⟨σv⟩ sourcing emitting birth spectra (p-¹¹B 3α continuum, D-D secondaries as the one rate functional on SOLV-3's slowing-down f, thermal-broadened neutron peaks). Annihilation via OFFL-4 with the factor-of-several band. Emit-once/transport-once; no private transport. |
| 2026-08-13 | 0.2 | R3 terminology: annihilation model-form reframed from a blanket "factor-of-several band" to OFFL-4's **per-quantity physics-list spread**; efficiency/waste-heat noted as computed downstream by SOLV-2/3 transport, not emitted here. |
| 2026-08-19 | 0.4 | **VISION_SCOPE v1.5 (Ben): the chemical burn-progress source (§3.6).** c ∈ [0,1] burnt-fraction field on `U`; blended unburnt↔equilibrium thermochemistry; one continuous rate law (SOLV-4.4: flame propagation via `S_L`×wrinkling + auto-ignition via `τ_ign`), closures offline-generated (OFFL-3) and unit-anchored (VAL-2); igniter = COUP-7 energy-deposit object; grid-independent front speed = the acceptance gate; reacting-measure thresholds feed the new COUP-4 `NEVER_IGNITED`/`FLAMEOUT` halts; class-`R` cell-local stiff source; runtime finite-rate networks remain out. §1/§2/§6 updated. *Same-day adversarial-review fixes:* SOLV-4.4 restated in **TFC form** (`ρ_u·S_T·|∇c|` + the **matched front-thickening diffusion `D_c`** as the explicit grid-independence mechanism — a bare source self-sharpens and its speed goes grid-set); **sub-cell partition rule stated** (common p, mass-weighted h, unburnt-isentrope closure for `h_u`/`T_u` — no second energy field); the flame/auto-ignition overlap addressed (deflagration regime is propagation-dominated; autoignitive end-gas addition is physics, each closure anchored where it dominates). |
| 2026-08-14 | 0.3 | **Review fix wave (N5, N6, E-2).** §3.2: quasi-static closure written explicitly — amplitude ODE (SOLV-4.2) + tabulated flux/adjoint shapes ψ/ψ† from OFFL-1 + adjoint-weighted delayed source (SOLV-4.3) and drifted-precursor loss ledger (flowing-fuel β_eff emerges from advection + importance weighting, configuration-general); kinetics declared a **class-`G` global stiff ODE** in COUP-3's SDC step + `EPS_RHO_STEP`/`EPS_POWER_STEP` Δt limiters; kernels' control-state axis noted (E-2). |
