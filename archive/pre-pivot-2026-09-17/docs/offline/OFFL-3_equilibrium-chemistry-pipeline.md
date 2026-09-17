# OFFL-3 — Equilibrium Chemistry Pipeline

| Field | Value |
|---|---|
| **ID** | OFFL-3 |
| **Family** | OFFL (Offline pipeline) |
| **Status** | Reviewed (2026-08-14) |
| **Depends on** | FND-5 (table schema), FND-1 (uncertainty); consumed by SOLV-1/7/8, OFFL-5 (transport feed) |
| **Version** | 0.6.4 (2026-08-24 plan S7: **unburnt v0.3.0 — a declared RE-GRIDDED VARIANT** (§3.3 R2 envelope doctrine, NOT a strict extension): the Z band narrows to the burnt surface's own grid band (0.145–0.195, envelope 0.155–0.185 — the premixed design-window class the blend intersects to anyway) so the rectangular h-envelope's floor, which binds at the O₂-rich Z edge, stops pinning the mid-Z design line's effective floor at ~158 K; the cold face re-derives at a **75 K binding edge** (measured: T at the envelope floor on the mid-Z line ≈ 96 K, covering the 120 K wall-cooled fill states — declared metastable supersaturated-vapor overhang toward the ~50 kPa+ corner only); sub-floor overhang = the base axis's own ~2-cell margin, no deeper cold prepends (deepest binding-corner node ~40 K, γ ≈ 1.15); h ceiling rises to 8.63e6; the ignition surface's T_u envelope (75–2934 K) covers the new envelope at BOTH ends — the envelope-consistency contract now stated two-faced; the wide-Z 0.2.0 grid stays on disk as the S5b species-work base). 0.6.3 (2026-08-24 plan S7: the **cold-side closure floor** — ignition v0.3.0 extends the T_u axis *down* (graded prepended nodes, old nodes bit-exact — the strict-extension discipline; envelope floor ~229.5 → **75 K**, covering the marched flame front's cold-reactant queries at T_u ≈ 150–230 K, and now covering the **unburnt surface's own validity floor** — the cold face of the 0.6.2 envelope-consistency contract); below the widened envelope the runtime carries a **declared non-reactive floor** (the cold analogue of SOLV-4 §3.6's BURN_COMPLETE domain guard): where p or T_u sits under the ignition surface's declared envelope floor both rate terms read **zero by declaration** — near-vacuum fill (~10²–10³ Pa) and sub-cryo corners are declared no-burn, never a refusal mid-march and never an interpolated guess. Generation gains an **isolated-zero refusal** (an S_L = 0 node with flammable neighbors on both T_u sides is a failed solve, not a flammability limit — refuse, don't record)). 0.6.2 (2026-08-24 plan S6 close: the **ignition-headroom envelope set** — burnt v0.4.0 ceiling +4.0e6 → +1.225e7 J/kg (rule: burnt ≫ unburnt, since the SOLV-4 §3.6 blend interrogates both branches at one cell enthalpy); unburnt v0.2.0 ceiling t_ceil 2200 → 2900 K (confined-ignition blast compression measured past the old ceiling on mid-transition cells); ignition v0.2.0 T_u ceiling ~1900 → ~3000 K (the **envelope-consistency contract**: the ignition surface's T_u envelope must cover T_u at the unburnt surface's h-ceiling). All three are **strict extensions** — existing nodes bit-exact, new nodes appended — so every in-old-envelope certified number is unchanged). 0.6.1 (2026-08-21: the **ignition/flame closure surface shipped** — the `S_L(p,T_u,Z)` + `τ_ign(p,T_u,Z)` two-column surface via Cantera `FreeFlame`/const-`P` reactor on `h2o2.yaml`, extinction as the columns' own values (`S_L→0`, `τ_ign` capped `τ_max`) — plan S6). 0.6 (2026-08-20: the **unburnt-reactant surface shipped** — the c = 0 cold branch, a frozen gas-phase ideal-gas reactant mixture via CEA, cryo-valid; the frozen↔shifting **bracket** promoted from a validation check to a shipped declared model-form band; the frozen-*advection* (p, h, {X_k}) surface + per-species diffusion re-deferred to plan **S5b** with the species-vector state — plan S5) |

---

## 0. Purpose

OFFL-3 is the **offline Python pipeline** that generates the chemical-slice tables the runtime reacting-flow
operator (SOLV-1) interpolates: **local-state equilibrium/frozen composition-thermo surfaces** (S22), the
chamber performance reference, quasi-1-D expansion **oracles** (SOLV-7/VAL-2 only), and B′ ablation tables.
Its Cantera **transport** computations are an **input feed to OFFL-5's spine assembly** — the spine (FND-7)
is the sole runtime provider of μ, k (S23). It runs **once per table version** and emits versioned HDF5 to
the FND-5 schema — the only cross-language seam (VISION_SCOPE §6). Python never runs at simulation time.

Read after FND-5 (the table contract it must satisfy) and SOLV-1 §3.4 (the consumer).

## 1. Scope & razor ruling
**Owns:** the CEA/Cantera equilibrium computations and the products they yield (local-state surfaces,
performance reference, expansion oracles, B′, the transport feed). **Defers:** the on-disk
schema/interpolation/error-budget → **FND-5**; the runtime use → **SOLV-1**; **runtime transport provision →
the spine (FND-7), assembled by OFFL-5** — OFFL-3 never ships a runtime transport table (S23); the
cold-solid/plasma constitutive corners → **OFFL-5** (the spine); the injector mixing tiers (prior η_c\* /
resolved stream conditions) → **COUP-7** (a boundary object, not a table).

**Razor ruling:** the reaction thermochemistry is done **completely, offline, at its correct scale**
(equilibrium free-energy minimization over the full Glenn species set), then homogenized into tables — the
standard tabulate-offline / interpolate-online method (FND-2 §3.4.1).

## 2. Interfaces & contracts
| Interface | Consumed by | Contract |
|---|---|---|
| **Equilibrium surface** *(S22)* | SOLV-1 | equilibrium composition + thermo `(T, ρ, {X_k}^eq, γ_eff, a, heat release)` vs **(p, h, Z)** — pressure, specific enthalpy, elemental composition/mixture fraction; the shifting-mode per-step projection lookup (SOLV-1 §3.4, S19) |
| **Frozen-path surface** *(S22; deferred → plan S5b)* | SOLV-1 | frozen-composition thermo `(T, ρ, γ, a)` vs `(p, h, {X_k})` for the frozen **advection** mode. Its axis set is the frozen-advection consumer's, so it lands **with the species-vector state** (plan S5b — SOLV-1 §3.1's `{ρX_k}` block); until then the delivered-performance frozen end is carried by the §3.2 **bracket**, not a field surface |
| **Unburnt-reactant surface** *(v0.6, SHIPPED — plan S5)* | SOLV-1/SOLV-4 §3.6 | frozen **reactant-mixture** thermo `(T, ρ, γ_fr, a)` vs (p, h, Z), valid to cryogenic T — the burn-progress blend's **c = 0** branch. Generated as the **gas-phase ideal-gas frozen mixture** of the propellant's reactant species at the elemental proportions Z (§3.3), by the **same CEA engine** as the equilibrium surface so the two share one enthalpy reference (the blend `h = (1−c)h_u + c·h_b` is only meaningful on one reference — SOLV-4 §3.6). Same FND-5 schema as the equilibrium surface, so `TableEos` binds either (SOLV-1 §3.4); the two-phase liquid-injection/vaporization coupling is SOLV-1's W4 drift-flux extension (plan S15), not this surface |
| **Ignition/flame closures** *(v0.6.1, SHIPPED — plan S6)* | SOLV-4 §3.6 | laminar flame speed `S_L(p, T_u, Z)` + induction time `τ_ign(p, T_u, Z)` on **one (p, T_u, Z) surface, two columns**, computed **offline** on the cited `h2-kinetics-mech` (`h2o2.yaml`): `S_L` from Cantera 1-D freely-propagating flames (`FreeFlame`, mixture-averaged transport), `τ_ign` from 0-D constant-pressure reactors (`IdealGasConstPressureReactor`). Extinction is the **columns' own values**, so a runtime cell never branches: `S_L` is **linear-valued** and **→ 0** where the flame fails to propagate (non-flammable ⇒ the flame solve does not light ⇒ `S_L = 0`); `τ_ign` is **log-valued** and **capped at a declared `τ_max`** (no `∞`) so the auto-ignition source vanishes at the quench/low-`T` corner. Holdout error bounds per FND-5; unit-anchored (VAL-2 §3.3, `h2-flame-speed`/`h2-ignition-delay`) |
| **Performance reference** | SOLV-7, SOLV-1 (§3.4 anchor knockdown) | `c*_ideal, T_c, γ, M̄` vs `(p_c, MR)` — a chamber-stagnation performance functional, not a field lookup |
| **Expansion oracles** *(quasi-1-D, demoted S22)* | SOLV-7, VAL-2 | shifting **and** frozen state `(T,p,ρ,h,s,γ_eff,a)` + species vs area/pressure ratio — **oracle cross-checks only, never interpolated by the field solver**; labeled `oracle` in FND-5 metadata |
| **Transport feed** *(S23; quantity set widened v0.5)* | OFFL-5 | Cantera mixture-averaged `μ, k` **plus the caloric companions `c_p, c_v, ∂h/∂Z\|_{p,T}`** vs `(T, p, Z)` — all five from the **one** Cantera evaluation (FND-7 §3.3: a second model for c_v would be a seam inside the spine). **Input to the spine assembly**; OFFL-3 still ships no runtime transport table — the spine (FND-7) is the **sole runtime provider**, and OFFL-5 §3.1a owns the emitted surface and its runtime (p, h, Z) coordinate |
| **B′ ablation table** | SOLV-8 | blowing-rate/recession coefficients vs surface state |

**Invariant:** every table carries FND-5 mandatory metadata (provenance, versions, envelope, per-column
units/interp-rule/error-bound); every value is reproducible from `{generator script hash, tool versions,
inputs}` (Tier-2 determinism, META-1 §2.1).

## 3. Method

### 3.1 Tools
**NASA CEA** (`nasa/cea`, Apache-2.0, Python bindings, >2000-species Glenn database) as the equilibrium/rocket-
performance engine, with **Cantera** (BSD-3) as an independent equilibrium cross-check and the source of
**transport** properties. **RocketCEA** (GPLv3) is a validation oracle only (kept out of the redistributed
pipeline). [META-3: `nasa-cea`, `cantera`, `nist-janaf`]

### 3.2 Frozen, equilibrium, and the kinetic middle *(the fidelity decision)*
Real delivered performance sits between the **frozen** and **shifting-equilibrium** limits; the
**frozen↔shifting gap is a discrete epistemic model-form dimension** sampled by COUP-5 (S19) and the JANNAF
**kinetic-efficiency** factor is tabulated so the bracket is data-anchored — the disagreement between models
**becomes the error bar** (S3), not a hidden choice. For LOX/LH₂ at MR≈5 the gap is ~0.8–1%. [META-3:
`jannaf-eff`]

**The bracket is a shipped declared band, not just a validation check** *(v0.6, plan S5).* The delivered-
performance frozen end is the **frozen-from-chamber rocket solution** (composition frozen at the chamber,
expanded frozen) against the shifting rocket solution — the pair a `RocketSolver` yields at one area ratio
(`n_frz = 1` vs shifting). Their vacuum-Isp gap is the model-form band, declared with its `jannaf-eff`
anchor and recorded in the run's pedigree, so a run reports the bracket rather than a point. **Two numbers,
not one:** this *raw bracket width* is a few percent and **widens with the area ratio** (the frozen limb
leaves more recombination energy unclaimed the longer it expands — **~3.7–4.8% at ε = 61 over the MR 5.0–5.5
window**, ≈ 4.3% at the RP-1311 example-8 anchor); the **~0.8–1% above is the `jannaf-eff`-anchored
*delivered* position inside that bracket** (H/O kinetics are fast, so the real engine hugs the shifting end).
The bracket is the outer model-form interval; the JANNAF figure says where in it the delivered value lands.
This is
**distinct from the frozen-*advection* field surface** (§2, deferred to plan S5b): the bracket is a
delivered-performance interval carried without advecting a species vector, so it is available now; the field
surface is the *mechanism* for marching the frozen limb as a live field, and it rides the species-vector
state. **The c = 0 unburnt branch is a third, orthogonal thing** (§2, §3.3): it selects *whether burnt*
(the burn-progress blend, SOLV-4 §3.6), while frozen↔shifting selects *how the burnt products behave*.

### 3.3 What is tabulated (and why offline)
**Runtime surfaces are keyed on local state (S22)** — a field-solver cell has no area ratio. The
**equilibrium surface** (composition, T, ρ, γ_eff, a, heat release) and the **frozen-path surface** are
tabulated vs **(p, h, Z)** — pressure, specific enthalpy, elemental composition/mixture fraction —
consistent with SOLV-1's advected-element shifting mode (S19). The **performance reference**
`(c*_ideal, T_c, γ, M̄)` vs `(p_c, MR)` is a chamber-stagnation functional read by SOLV-7 and by SOLV-1's
§3.4 anchor knockdown. **B′ ablation coefficients** vs surface state. The **area-ratio/pressure-ratio
expansion tables survive only as quasi-1-D oracles** — SOLV-7's C_F cross-check and VAL-2's nozzle-envelope
checks — labeled `oracle` in their FND-5 metadata and **never interpolated by the field solver** (naive
wiring would double-count the expansion physics the grid already computes). **Transport** (Cantera
mixture-averaged `μ, k, c_p, c_v, ∂h/∂Z|_{p,T}` vs (T, p, Z) — Cantera's natural coordinate) is computed here
but shipped **only as OFFL-5's input feed** (S23); OFFL-5 §3.1a converts it to the runtime local-state
(p, h, Z) coordinate and owns the shipped surface. Grid spacing is chosen offline to meet the FND-5 interpolation-error budget (interpolate in
linearizing space where curvature warrants). Runtime reads via multilinear/monotone-cubic interpolation only
— never a runtime equilibrium solve.

**The unburnt-reactant (c = 0) surface** *(v0.6, plan S5)* is tabulated on the **same (p, h, Z) coordinate
and the same FND-5 schema as the equilibrium surface** (`temperature, density, gamma_eff, sound_speed,
mbar`), so `TableEos` binds either with no new occupant (SOLV-1 §3.4). Its physics is the **gas-phase
ideal-gas frozen mixture** of the propellant's reactant species (for LOX/LH₂: gaseous H₂ + O₂ at the mass
proportions Z sets), evaluated by the **same CEA engine** as the equilibrium surface — `calc_property` on
the reactant `Mixture` (the `calc_property(cea.ENTHALPY, …)` mechanism `injection_enthalpy` uses, on the
gasified reactant mixture — same absolute formation reference, gas rather than liquid species) — so the
unburnt and burnt branches
carry **one enthalpy reference** (the S6 blend `h = (1−c)h_u + c·h_b` is a category error on two references).
`T(h, Z)` is the (pressure-independent) inverse of the frozen mixture's monotone enthalpy; `ρ = pM̄/(R̄T)`
(ideal gas — the same assumption CEA makes for the products, so the two branches are *consistently* ideal,
differing only in composition); `c_p,fr = ∂h/∂T`, `c_v = c_p,fr − R̄/M̄`, `γ_fr = c_p/c_v`,
`a = √(γ_fr R̄T/M̄)` (CEA's own sound-speed algebra). Its **envelope floor (v0.3.0, plan S7) is ~96 K on
the mid-Z design line / 75 K at the binding O₂-rich edge** (the rectangular h-envelope's floor binds at the
O₂-rich Z corner), over a **Z envelope of 0.155–0.185**: v0.3.0 is the **narrow-Z re-gridded variant**
matched to the burnt surface's own grid band (0.145–0.195) — the premixed design-window class the blend
intersects to anyway — while the **wide-Z 0.2.0 grid stays on disk as the S5b species-work base**. The
reactant elements (zero formation enthalpy) give a still-physical γ ≈ 1.47 at the floor, while the
equilibrium surface (whose products want to condense) already refuses. CEA *converges* below the floor (down to ~35 K), which is what
lets the grid overhang the envelope for interpolation, but that is convergence, not validity: below ~60 K
the extrapolated NASA-polynomial thermo degrades (γ → 1.14 at 40 K) and `state_php` refuses where `c_p ≤ 0`.
The sub-floor nodes are the declared metastable-model overhang, gated off at runtime (the same pattern as
the equilibrium surface's gas-only overhang). Because the p-invariant columns depend
on p only through the exact ideal-gas `ρ ∝ p`, a coarse log-p axis is interpolation-exact and the measured
holdout bound lives on the h-axis (the c_p(T) curvature). The **cryogenic-liquid injection state** and its
vaporization isentrope are the two-phase drift-flux extension (SOLV-1 W4, plan S15); this surface is the
gas-phase branch.

**The ignition/flame closure surface** *(v0.6.1, plan S6)* is tabulated vs **(p, T_u, Z)** — the coordinate
the SOLV-4 §3.6 rate law is keyed on (`T_u` from the unburnt-isentrope closure, not the burnt temperature).
Its physics is **finite-rate `h2-kinetics-mech`**, done at its correct scale offline and homogenised into
the surface exactly as the equilibrium thermochemistry is: the reactant mixture at (p, T_u, Z) — the same
`{H₂: Z, O₂: 1−Z}` mass mapping the frozen-reactant surface uses (§3.3, `gasify`d streams) — is handed to
**Cantera** on `h2o2.yaml`. `S_L` is the `FreeFlame` burning velocity `velocity[0]` (mixture-averaged
transport, declared fixed refine criteria; a non-lighting solve ⇒ `S_L = 0`, the flammability-limit value).
`τ_ign` is the `IdealGasConstPressureReactor` induction time (declared temperature-rise trigger; a
non-igniting integration ⇒ `τ_ign = τ_max`, the declared cap). Both are **interpolated at runtime, never
re-solved** — no runtime kinetics network (SOLV-4 §3.6). The declared band absorbs the mechanism-vs-
measurement gap; the surface is validated to reproduce a fresh Cantera solve within its FND-5 holdout bound
(interpolation fidelity) **and** the mechanism to reproduce the VAL-2 measured anchors within a declared
model-form band (`h2-flame-speed`, `h2-ignition-delay`). Because the Cantera pin is part of the record, the
provenance deck stamps `cantera {version}` and the mechanism file (the OFFL-5 transport-feed pattern).

**The cold-side closure floor** *(v0.6.3, plan S7).* The startup march legitimately queries this surface
below any envelope a Cantera flame solve can honestly cover: the pre-ignition fill sits at ~10² Pa and the
front, once lit, propagates into injector-cold reactants at T_u ≈ 150–230 K. The two regions get the two
honest treatments. (1) **Where flames really burn, the surface has real rows:** the T_u axis extends *down*
with graded prepended nodes (v0.3.0; old nodes bit-exact — strict extension) so the declared envelope floor
(~75 K) covers every unburnt-branch temperature the blend can produce — the **cold face** of the 0.6.2
envelope-consistency contract (the ignition envelope must cover the unburnt surface's own validity range at
both ends, or a legal blend state becomes a closure refusal). (2) **Below the declared envelope floor the
medium is declared non-reactive:** the runtime rate law carries the cold analogue of the SOLV-4 §3.6
BURN_COMPLETE domain guard — where the cell's p or T_u sits under the surface's declared envelope floor,
both rate terms read **zero by declaration** (no query, no refusal, no clamp of a mid-range value). This is
a *declared model floor*, honest because cm-scale flames genuinely do not propagate at the ~10²–10³ Pa fill
pressures it covers (quenching distance ∝ 1/p reaches the hardware scale) and the certified light-off must
be *measured* to occur well inside the live envelope (recorded per run). Generation arms an
**isolated-zero refusal**: along the T_u axis, an `S_L = 0` node with flammable (nonzero) neighbors on both
sides is a failed cold solve masquerading as a flammability limit — the generator refuses rather than
records it (extinction is contiguous from an edge, never an interior hole).

**Local-composition envelope (R2, architect-now):** at the resolved-mixing tier (SOLV-1 §3.4, COUP-7
§3.2.1) local composition leaves the global-MR line, so the (p, h, Z) surfaces must be generatable over the
**full local mixture-fraction range** (deep fuel-rich → deep ox-rich), not just the design-MR window. This
is a **table-envelope setting of the same pipeline, not a new pipeline**: W2 tables grid the design window
densely; a resolved-tier table version widens the Z axis (dense core, coarse wings, within the FND-5
interpolation budget). The declared envelope always states which range a given table version actually covers
— the FND-5 refusal then keeps a resolved-tier run from interpolating a design-window-only table.

### 3.4 Uncertainty typing
Thermochemical model-form (equilibrium assumption; the frozen↔shifting band), thermodynamic-data uncertainty
(Glenn/JANAF), and the interpolation-error bound (FND-5) are each typed (FND-1) so Sobol can attribute a
performance spread to its source.

## 4. Coupling relationships
- **FND-5** defines the schema/loader/interpolation and enforces the envelope; **SOLV-1** interpolates the
  (p, h, Z) surfaces at runtime; **SOLV-7** reads c\*_ideal from the performance reference and uses the
  expansion oracles for its C_F cross-check; **SOLV-8** reads B′; **OFFL-5** consumes the transport feed —
  the spine (FND-7) is the sole runtime μ,k provider (S23); **COUP-7** supplies the injector tiers
  separately (prior η_c\* or resolved stream conditions — a boundary object, not a table).
- **VAL-2** validates the tables against RP-1311 example cases and the RL10 anchor (expansion oracles serve
  its nozzle-envelope checks).

## 5. Uncertainty & validity
Originates the **thermochemical model-form band** (incl. the frozen↔shifting epistemic dimension) and
carries the thermodynamic-data band; the envelope is the `(p, h, Z)` ranges of the runtime surfaces plus
`(p_c, MR)` for the performance reference (area-ratio ranges apply to the oracles only), enforced by FND-5
(refuse out-of-envelope). Ladder rung:
benchmark (ii) against RP-1311/CEA reference cases; the *delivered-Isp* claim is validated at the system tier
via RL10 (VAL-2).

## 6. Validation plan
1. **RP-1311 reproduction:** CEA example cases reproduced within the manuals' tolerances (self-verification of
   the pipeline). [META-3: `nasa-cea`]
2. **CEA↔Cantera cross-check:** the two equilibrium solvers agree within tolerance on chamber state.
3. **Frozen/shifting bracket:** the two tables bracket the JANNAF-corrected delivered value.
4. **Table round-trip + error bound:** interpolated values meet the stored FND-5 `interp_error_bound`.
5. **Provenance:** every table carries generator hash + tool/data versions; regenerable.
6. **Coordinate consistency (S22):** the (p, h, Z) equilibrium surface reproduces the (p_c, MR) chamber
   solution along the design line, and the quasi-1-D expansion oracles agree with a field-solver quasi-1-D
   nozzle run (the two parameterizations agree where both are defined).

## 7. References
META-3 keys: `nasa-cea`, `cantera`, `nist-janaf`, `jannaf-eff`, `bprime`, `cstar-cf-defs`, `h2-kinetics-mech`, `h2-flame-speed`, `h2-ignition-delay` (§2 ignition/flame closures). Depends on FND-5
(schema/interp), FND-1 (uncertainty); consumed by SOLV-1/7/8, OFFL-5 (transport feed), COUP-7 (injector
separate), VAL-2 (oracles).

*(No open questions — chemistry fidelity (both frozen+shifting, gap as UQ band) and engine choice (`nasa/cea`
Apache-2.0 + Cantera cross-check) resolved 2026-07-21.)*

## 8. Change log
| Date | Version | Change |
|---|---|---|
| 2026-08-24 | 0.6.4 | **Plan S7 — unburnt v0.3.0: a declared RE-GRIDDED VARIANT (§3.3's R2 envelope doctrine — a table-envelope setting of the same pipeline), NOT a strict extension.** The **Z band narrows to the burnt surface's own grid band** (0.145–0.195, envelope 0.155–0.185 — the premixed design-window class the blend intersects to anyway), because the **rectangular h-envelope's floor binds at the O₂-rich Z edge**: the wide S5 band left the mid-Z design line's *effective* floor at ~158 K — while the S7 startup march holds REAL gas-phase states below it (fill gas cooled by the 120 K coolant wall at 1–20 kPa sits 45+ K above its own O₂ saturation, T_sat ≈ 61–82 K at the band's O₂ partial pressures — honestly gas, wrongly refused). The cold face re-derives the box floor at a **75 K binding edge** (`t_floor_ext` = 75 K; declared metastable supersaturated-vapor overhang toward the ~50 kPa+ corner only; measured: T at the envelope floor on the mid-Z line ≈ **96 K**, covering the 120 K wall-cooled states). **Sub-floor overhang = the base axis's own ~2-cell (3 %) margin, no deeper cold prepends** — the deeper NASA-polynomial region degrades toward `c_p ≤ 0` (~35 K) and stays out of the artifact; the deepest binding-corner node measures ~40 K at γ ≈ 1.15 (grid min 1.151); the offline physicality gate relaxed 1.1 → 1.05 as declared headroom, its units-bug purpose intact. The **wide-Z 0.2.0 grid remains on disk as the S5b species-work base** (narrow-Z is the blend's runtime binding, not a replacement of the species program's coordinate). The **h ceiling rises to 8.63e6** (narrow Z lifts min-Z h(2900 K)); the ignition surface's T_u envelope (**75–2934 K**) covers the new unburnt envelope **at BOTH ends** — the 0.6.2/0.6.3 envelope-consistency contract now stated two-faced (hot face: T_u at the unburnt h-ceiling; cold face: the unburnt validity floor). §3.3's unburnt paragraph updated (envelope claim + variant/base split). |
| 2026-08-24 | 0.6.3 | **Plan S7 — the cold-side closure floor (ignition v0.3.0 + the declared non-reactive floor).** The S6 review wave flagged that `Combustion` queries `S_L`/`τ_ign` on every non-burnt cell while the ignition surface floored at T_u ≈ 229.5 K / p ≈ 6.75 kPa — a cryo-fill or near-vacuum cell refuse-halts at first contact, and a pure guard-to-zero *at those floors* would falsely stall the front into injector-cold gas (the march's flame consumes reactants at T_u ≈ 150–230 K at chamber pressure — real S_L territory, not extinction). Cure is the composite in §3.3: **(1) ignition v0.3.0** prepends graded cold T_u nodes (60/90/120 K at Δ30 under the 150 K base — the FND-5 loader requires monotone axes, not uniform spacing) with real Cantera solves; every pre-existing node bit-exact (strict extension, verified column-by-column); envelope floor 229.5 → **75 K**, now covering the unburnt branch's own validity floor (the contract's cold face). The p-axis **floor** is unchanged (FreeFlame below ~3 kPa does not honestly converge; that region belongs to the declared floor), but the **top p-cell is refined**: the coarse 6-node log axis put the declared p-envelope *ceiling* at 4.44 MPa purely through the half-log-cell inset — *below* the RL10 chamber class (p_c ≈ 3.2 MPa + ignition transients), so a mid-transition cell compressed past it turned a legal blend state into a closure refusal (the S7 stiff mini-sim caught it — the p-ceiling face of the 0.6.2 envelope-consistency contract). The inserted geometric-midpoint node (~4.44 MPa, real Cantera solves — the 100-bar row already converges) moves the envelope ceiling to **6.67 MPa**; refinement, not extension — old nodes bit-exact, in-cell interpolants legitimately improve. **(2) The declared non-reactive floor**: below the surface's declared envelope floor in p or T_u, both rate terms are zero by declaration — the cold analogue of the BURN_COMPLETE guard, owned by SOLV-4 §3.6 at runtime; the light-off pressure is measured and recorded per run to sit well inside the live envelope. **(3) Isolated-zero generation refusal** armed (an interior S_L = 0 hole along T_u = a failed solve, refused). **KNOWN LIMIT (S7 review):** three sub-envelope cold-corner nodes (p = 3 kPa row, Z = 0.252/0.30, T_u = 60/90 K) carry failed-solve-scale S_L (137/155/5.2 m/s vs the real ≤ 4.3 m/s cold maximum); they sit below the p-envelope floor but bracket legal fuel-rich queries just above it, which would interpolate a spurious tens-of-m/s S_L. Unreachable at the premixed startup band's uniform Z = 1/6; the cure is a max-S_L cold-row generation gate (the isolated-zero refusal's spurious-NONZERO sibling), riding the next regen. |
| 2026-08-24 | 0.6.2 | **Plan S6 close — the ignition-headroom envelope set (burnt v0.4.0, unburnt v0.2.0, ignition v0.2.0).** The S6 igniter mini-sims (`spark_box`/`quench_box`) measured the S5/S6 surfaces' ceilings as **ignition-incompatible**, in two layers. (1) **Envelope inversion:** the v0.3.2 products ceiling (+4.0e6 grid / +3.8e6 envelope, sized for *station* transients) sat **below** the unburnt surface's (+5.0e6) — but the SOLV-4 §3.6 blend interrogates **both** branches at the same cell enthalpy, so an igniting kernel (auto-ignition needs `T_u` ~ 1000–1300 K ⇒ h ~ 2.5–4e6 J/kg) died on the **burnt** branch while still valid cold gas — a table-envelope defect, not an igniter-sizing constraint (Ben ruling: raise the ceiling, never throttle the spark to fit the table). **Burnt v0.4.0:** h-grid +4.0e6 → **+1.225e7** J/kg (envelope +1.2e7; T ~ 5000–6000 K class, trivially CEA-convergent, still chemistry), with the standing rule **burnt ceiling ≫ unburnt ceiling** so a burning cell can never refuse where the same cold gas was fine. (2) **Blast headroom + the envelope-consistency contract:** with (1) fixed, confined-ignition blast compression drove **mid-transition (0 < c < 1) cells** to h ≈ 5.3e6 — past the unburnt ceiling while still reading that branch. **Unburnt v0.2.0:** `t_ceil` 2200 → **2900 K** (h-ceiling ≈ +6.9e6 envelope / +7.06e6 grid; superheated reactants are the metastable branch's self-consistent transient — τ_ign there is sub-µs). **Ignition v0.2.0:** T_u ceiling ~1900 → **~3000 K**, honoring the new **contract: the ignition surface's T_u envelope must cover T_u at the unburnt surface's h-ceiling** (the rate law queries S_L/τ_ign at T_u(p, h, Z), so a mismatch turns a legal blend state into a closure refusal); the new rows are the auto-ignitive regime (S_L reads 0 via the ceiling guard, τ_ign sub-µs — the columns' own extinction/instant-burn values, Rule 12). Discipline throughout: **strict extensions** — every pre-existing node reproduced bit-exact, new nodes appended at the same spacing, verified column-by-column at regeneration. For the burnt table even the stored interp bounds are bit-identical (only pins/digests move); the unburnt/ignition bounds legitimately re-measure over their new rows (unburnt temperature abs bound 1.66 → 2.26 K). Certificates: three stations + convergence byte-identical; station 1 gains only the new ρc MMS column (orders 1.92–2.17, all pre-existing values byte-unchanged — the NCOMP widening's certificate face) and station 5 only its provenance text (the v0.4.0 note). Envelope widening remains a table-version setting of the same pipeline (§3.3 R2 doctrine), never a new pipeline. |
| 2026-08-21 | 0.6.1 | **The ignition/flame closure surface shipped (plan S6 — landed with the code).** §2/§3.3: the `S_L(p,T_u,Z)` + `τ_ign(p,T_u,Z)` closures move from a contract row to a **shipped two-column (p, T_u, Z) surface**, generated on the cited `h2-kinetics-mech` (`h2o2.yaml`): `S_L` = Cantera `FreeFlame` burning velocity (a non-lighting solve ⇒ `S_L = 0`, **linear-valued** so extinction is representable), `τ_ign` = `IdealGasConstPressureReactor` induction time (**log-valued**, capped at a declared `τ_max` so there is no `∞` and the auto-ignition source vanishes smoothly at quench). Extinction is thus the surfaces' own values, never a runtime branch (SOLV-4 §3.6, Rule 12). Cantera-based provenance (deck stamps `cantera {version}` + mechanism, per the OFFL-5 transport-feed pattern). Validated: fresh-Cantera interpolation-fidelity holdout (FND-5 bound) + mechanism-vs-measured VAL-2 anchors within a declared model-form band. No FND-5 schema change; no existing table moved (a new artifact, no station consumer until S7). |
| 2026-08-20 | 0.6 | **The cold/unburnt branch + the shipped bracket (plan S5 — landed with the code).** §2: the **unburnt-reactant surface** moves from a contract row to **shipped** — a gas-phase ideal-gas frozen reactant mixture vs (p, h, Z), valid to cryogenic T, generated by the **same CEA engine and enthalpy reference** as the equilibrium surface (§3.3) so the SOLV-4 §3.6 blend is on one reference; it uses the equilibrium surface's FND-5 schema, so `TableEos` binds it with no new occupant. §3.2: the **frozen↔shifting bracket** is promoted from a validation check to a **shipped declared model-form band** (the frozen-from-chamber vs shifting rocket-Isp gap, `jannaf-eff`-anchored, recorded in the pedigree). §2: the frozen-**advection** (p, h, {X_k}) field surface and per-species diffusion are **re-deferred to plan S5b** with the species-vector state (SOLV-1 §3.1's `{ρX_k}` block) — the delivered-performance frozen end is the bracket until then, and per-species D needs the per-species gradients the widened state carries (FND-7 §3.3). No FND-5 schema change (the unburnt surface reuses the equilibrium schema). |
| 2026-08-20 | 0.5 | **Transport-feed quantity set (plan S4 — landed with the code).** §2/§3.3: the S23 Cantera transport feed widens from `μ, k` to **`μ, k, c_p, c_v, ∂h/∂Z\|_{p,T}`** — all five from the one Cantera evaluation, because the resolved diffusion operator needs the caloric companions (the implicit T-solve's linearization slope and the species-enthalpy diffusion-flux coefficient) and a second model for them would be a seam inside the spine (FND-7 §3.3). The S23 ruling is otherwise unchanged and restated: OFFL-3 ships **no** runtime transport table; OFFL-5 §3.1a converts the (T, p, Z) feed to the runtime local-state (p, h, Z) coordinate and owns the shipped surface; the spine (FND-7) is the sole runtime provider. |
| 2026-08-19 | 0.4 | **VISION_SCOPE v1.5 / SOLV-4 §3.6 products.** §2 gains the **unburnt-reactant surface** (the burn-progress c = 0 branch, cryo-valid) and the **ignition/flame closure surfaces** `S_L(p,T_u,Z)` + `τ_ign(p,T_u,Z)` — offline finite-rate chemistry (Cantera flames/reactors on the cited `h2-kinetics-mech`) is in scope *for table generation only*; runtime stays interpolation. Closure envelopes carry the flammability/quench limits so extinction is the closure's own smooth behavior. Unit-anchored per VAL-2 §3.3 (`h2-flame-speed`, `h2-ignition-delay`). |
| 2026-08-14 | 0.3 | **Post-review fix wave (S22, S23).** Runtime tables **re-parameterized to local state**: equilibrium + frozen-path composition/thermo surfaces vs **(p, h, Z)** (consistent with SOLV-1's S19 advected-element mode); the (p_c, MR) chamber product survives as the **performance reference** (a stagnation functional for SOLV-7 and the §3.4 anchor knockdown); **area-ratio expansion tables demoted to SOLV-7/VAL-2 quasi-1-D oracles** (`oracle`-labeled, never interpolated by the field solver). **OFFL-3→SOLV-1 transport interface deleted** — two providers of μ,k was a seam: Cantera transport becomes an **input feed to OFFL-5's spine assembly**; the spine (FND-7) is the sole runtime provider; the transport axis uses the same local (T, p, Z) coordinate. §2/§3.3/§4/§5 reworded to match; coordinate-consistency validation item added. |
| 2026-08-13 | 0.2 | **R2 applied — local-composition envelope.** §3.3 gains the resolved-mixing-tier requirement: property surfaces generatable over the full local mixture-fraction range as a table-envelope setting (dense design core, coarse wings; FND-5 refusal guards a resolved-tier run against design-window-only tables). Defer/coupling wording updated to the tiered injector (COUP-7 §3.2.1). |
| 2026-07-21 | 0.1 | Initial draft. `nasa/cea` (Apache-2.0) + Cantera cross-check pipeline; chamber/expansion/transport/B′ tables to FND-5; **both frozen and shifting-equilibrium** emitted with the frozen↔shifting gap + JANNAF kinetic efficiency as an explicit UQ band; offline-tabulate/online-interpolate (no runtime equilibrium solve); provenance + envelope per FND-5; validated vs RP-1311 and RL10. |
