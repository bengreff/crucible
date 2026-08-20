# OFFL-3 — Equilibrium Chemistry Pipeline

| Field | Value |
|---|---|
| **ID** | OFFL-3 |
| **Family** | OFFL (Offline pipeline) |
| **Status** | Reviewed (2026-08-14) |
| **Depends on** | FND-5 (table schema), FND-1 (uncertainty); consumed by SOLV-1/7/8, OFFL-5 (transport feed) |
| **Version** | 0.5 (2026-08-20: transport-feed quantity set widened to the caloric companions c_p/c_v/∂h∂Z — plan S4) |

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
| **Frozen-path surface** *(S22)* | SOLV-1 | frozen-composition thermo `(T, ρ, γ, a)` vs `(p, h, {X_k})` for the frozen advection mode |
| **Unburnt-reactant surface** *(v0.4)* | SOLV-1/SOLV-4 §3.6 | frozen **reactant-mixture** thermo vs (p, h, Z), valid to cryogenic T — the burn-progress blend's c = 0 branch |
| **Ignition/flame closures** *(v0.4)* | SOLV-4 §3.6 | laminar flame speed `S_L(p, T_u, Z)` + induction time `τ_ign(p, T_u, Z)` surfaces, computed **offline** (Cantera 1-D freely-propagating flames + 0-D constant-pressure reactors on the cited `h2-kinetics-mech`); measured-data holdout bounds; envelopes = flammability/quench limits (outside: S_L → 0, τ_ign → ∞ smoothly — the closure carries its own extinction) |
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
Emit **both** the frozen-path and **shifting-equilibrium** surfaces. Real delivered performance sits between
them; the **frozen↔shifting gap is a discrete epistemic model-form dimension** sampled by COUP-5 (S19) and
the JANNAF **kinetic-efficiency** factor is tabulated so the bracket is data-anchored — the disagreement
between models **becomes the error bar** (S3), not a hidden choice. For LOX/LH₂ at MR≈5 the gap is ~0.8–1%.
[META-3: `jannaf-eff`]

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
| 2026-08-20 | 0.5 | **Transport-feed quantity set (plan S4 — landed with the code).** §2/§3.3: the S23 Cantera transport feed widens from `μ, k` to **`μ, k, c_p, c_v, ∂h/∂Z\|_{p,T}`** — all five from the one Cantera evaluation, because the resolved diffusion operator needs the caloric companions (the implicit T-solve's linearization slope and the species-enthalpy diffusion-flux coefficient) and a second model for them would be a seam inside the spine (FND-7 §3.3). The S23 ruling is otherwise unchanged and restated: OFFL-3 ships **no** runtime transport table; OFFL-5 §3.1a converts the (T, p, Z) feed to the runtime local-state (p, h, Z) coordinate and owns the shipped surface; the spine (FND-7) is the sole runtime provider. |
| 2026-08-19 | 0.4 | **VISION_SCOPE v1.5 / SOLV-4 §3.6 products.** §2 gains the **unburnt-reactant surface** (the burn-progress c = 0 branch, cryo-valid) and the **ignition/flame closure surfaces** `S_L(p,T_u,Z)` + `τ_ign(p,T_u,Z)` — offline finite-rate chemistry (Cantera flames/reactors on the cited `h2-kinetics-mech`) is in scope *for table generation only*; runtime stays interpolation. Closure envelopes carry the flammability/quench limits so extinction is the closure's own smooth behavior. Unit-anchored per VAL-2 §3.3 (`h2-flame-speed`, `h2-ignition-delay`). |
| 2026-08-14 | 0.3 | **Post-review fix wave (S22, S23).** Runtime tables **re-parameterized to local state**: equilibrium + frozen-path composition/thermo surfaces vs **(p, h, Z)** (consistent with SOLV-1's S19 advected-element mode); the (p_c, MR) chamber product survives as the **performance reference** (a stagnation functional for SOLV-7 and the §3.4 anchor knockdown); **area-ratio expansion tables demoted to SOLV-7/VAL-2 quasi-1-D oracles** (`oracle`-labeled, never interpolated by the field solver). **OFFL-3→SOLV-1 transport interface deleted** — two providers of μ,k was a seam: Cantera transport becomes an **input feed to OFFL-5's spine assembly**; the spine (FND-7) is the sole runtime provider; the transport axis uses the same local (T, p, Z) coordinate. §2/§3.3/§4/§5 reworded to match; coordinate-consistency validation item added. |
| 2026-08-13 | 0.2 | **R2 applied — local-composition envelope.** §3.3 gains the resolved-mixing-tier requirement: property surfaces generatable over the full local mixture-fraction range as a table-envelope setting (dense design core, coarse wings; FND-5 refusal guards a resolved-tier run against design-window-only tables). Defer/coupling wording updated to the tiered injector (COUP-7 §3.2.1). |
| 2026-07-21 | 0.1 | Initial draft. `nasa/cea` (Apache-2.0) + Cantera cross-check pipeline; chamber/expansion/transport/B′ tables to FND-5; **both frozen and shifting-equilibrium** emitted with the frozen↔shifting gap + JANNAF kinetic efficiency as an explicit UQ band; offline-tabulate/online-interpolate (no runtime equilibrium solve); provenance + envelope per FND-5; validated vs RP-1311 and RL10. |
