# SOLV-4 — Reaction Sources

| Field | Value |
|---|---|
| **ID** | SOLV-4 |
| **Family** | SOLV (Runtime unified-grid operators) |
| **Status** | Reviewed (2026-08-14) |
| **Depends on** | FND-2, OFFL-1 (fission data), OFFL-4 (annihilation), COUP-3 (integration), SOLV-2/SOLV-3 (transport); FND-1 |
| **Version** | 0.4.5 (2026-08-24 plan S7 ◆C2 shake-out: the sub-cell partition extension made **two-sided** — reactant sub-state pins at the NEAR edge of its envelope under heat loss (floor) or spark/blast superheat (ceiling), products absorb the balance; rate-law coordinates read the pinned state so closures stay live and the class-R reaction completes superheated mid-transition cells). 0.4.4 (2026-08-24 plan S7: §3.6 integration split of record — propagation reaction + front-thickening diffusion stay **explicit class `A`** (non-stiff by construction), the auto-ignition term becomes the **cell-local implicit class-`R` occupant** (COUP-3 §3.3: fixed-structure backward-Euler node solve, linear in `ρc` at frozen `τ_ign`, `N_TAU_REFREEZE = 2`; the burnt fixed point honored **exactly** — a stiff update parks AT `ρ(1−BURN_COMPLETE)`, the exact integral of the declared law, not a clamp; measured: `dt/τ = 3.04` marches where the explicit tier halted); the **cold-side non-reactive floor** declared (two faces: the ignition surface's own envelope floor + the unburnt branch's h-floor) and the **cold-side partition extension** (reactant sub-state pinned at the unburnt h-floor under real heat loss, `B_PARTITION_MIN = 0.01`); `quench_box` **REALIZED** (§6.5: cold-wall FLAMEOUT with the full blend↔class-D coupling, R 2.06e-4 → 0 at burned = 0.150, adiabatic control alive); `turbulent-flame-speed` **PINNED**, first consumer the S7 startup config's declared wrinkling factor). 0.4.3 (2026-08-24 plan S6 close: the blend's pure-limit thresholds made **asymmetric** — pure burnt from `1 − 2·BURN_COMPLETE` (one owner; strictly contains the pinning attractor — `c` asymptotes ~2e-7 *below* the fixed point, so an equal threshold is a knife edge), pure unburnt at the original `1e-9` (no attractor on the cold side, and the cold-side crossing step scales by `v_b/v_u` ≈ 7.8×, so a 2e-3 threshold there would step density ~1.5 % — outside the unburnt density bound); the same domain guard added to the reacting-measure/consumption-rate diagnostics; igniter mini-sims obey the bounded-pulse spark discipline, COUP-7 0.4.1). 0.4.2 (2026-08-21: §3.6 **built — plan S6**; the blend is the energy-conserving flamelet form (both branches at the common `(p, h, Z)`, mass-weighted specific volume — the projection generalized; supersedes the injection-isentrope closure so the igniter's enthalpy deposit fires `τ_ign`), the SOLV-4.4 operator-class decomposition + S6-explicit/S7-class-`R` scoping stated, extinction pinned to the surfaces' own behaviour (`S_L→0`, `τ_ign` capped)) |

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
**Sub-cell partition (the blend rule) — the energy-conserving flamelet form (v0.4.2, plan S6, made
explicit).** Both sub-states share the cell pressure `p` **and the cell's actual static specific enthalpy**
`h = e + p/ρ`. That the two share `h` is not an assumption but the **adiabatic constant-pressure combustion
identity on S5's shared CEA formation reference**: burning at constant `p` with no heat loss conserves
enthalpy, so `h_b = h_u = h` (the burnt products at the *same* absolute `h` as the reactants — this is the
definition of the adiabatic flame temperature, and it is exactly why the shared reference S5 shipped is a
correctness requirement). The mass-weighted split `h = (1−c)h_u + c·h_b` is then satisfied identically for
every `c`, with **no `h_b = (…)/c` division and so no `c → 0` singularity**. The two sub-states differ only
in *composition* — unburnt reactants vs equilibrium products — read from the two OFFL-3 branches at the
**same** `(p, h, Z)`:
- **Unburnt** `ρ_u, T_u = ` unburnt-reactant surface at `(p, h, Z)`; **burnt** `ρ_b, T_b = ` equilibrium
  surface at `(p, h, Z)` (the latter *is* the pure shifting-mode interrogation).
- **Density closure = mass-weighted specific volume:** `1/ρ = (1−c)/ρ_u(p, h, Z) + c/ρ_b(p, h, Z)`, and the
  pressure `p` is the **root of this identity** over the admissible bracket — the SOLV-1 §3.4 equilibrium
  projection **generalized to interrogate both branches** (the same deterministic root find). It recovers
  the limits exactly: `c = 1` ⇒ `1/ρ = 1/ρ_b(p, h, Z)` (the pure burnt projection, **bit-for-bit** shifting
  mode) and `c = 0` ⇒ `1/ρ = 1/ρ_u(p, h, Z)` (the pure unburnt projection).
- **`T_u = T_unburnt(p, h, Z)` is the reaction coordinate for both `S_L` and `τ_ign`.** Because a flame is
  energy-conserving, `h` is ≈ uniform across the front, so `T_unburnt(p, h, Z)` reads the **cold** reactant
  temperature there (the flame consumes cold reactants — right for `S_L`); where an **igniter energy deposit
  or adiabatic compression raises `h`**, `T_u` rises with it, so `τ_ign` drops and the induction term
  fires — the igniter works **because** `T_u` sees the deposited enthalpy (the earlier injection-isentrope
  phrasing `T_u(p, Z)` is the *unheated* special case `h = h_isentrope(p)` of this general
  `h`-dependent form, and it could not fire on a non-isentropic deposit; superseded). No second energy field
  is carried — `h` is the cell's own conserved enthalpy.
- The **blended acoustic speed** is the mass-weighted `a² = (1−c)a_u² + c·a_b²`, and the diagnostic cell
  temperature the mass-weighted `T = (1−c)T_u + c·T_b` (declared model-forms recovering each pure limit; the
  front speed is set by `S_T`/`D_c`, not the acoustics, so the blend rides the closure band). The whole
  partition is one continuous formula over `c ∈ [0,1]` — no regime branch (Rule 12).

**Cold-side partition extension (S7, v0.4.4).** The shared-`h` flamelet identity presumes *adiabatic*
burning; under real heat loss (cold-wall quench) a mid-transition cell's mixture enthalpy legitimately
falls **below any representable reactant state**. Below the unburnt branch's h-floor the declared extension
**pins the reactant sub-state AT the floor and the products absorb the balance**,
`h_b = (h − (1−c)·h_floor)/c` — continuous at the floor, mass-consistent, and self-limiting (`h_b` walks
down the burnt branch until **its** envelope refuses — the blend's true cold edge). Below a declared
trace-weight floor `B_PARTITION_MIN = 0.01` the products sub-state reads the cell `h` directly (the balance
form amplifies sub-floor deficits by `1/c`; at < 1 % weight the products' state is immaterial within the
density columns' own bounds). Both are **declared model forms, never clamps**.

**The one continuous rate law (Rule 12 — no ignition `if`), the reaction-diffusion (thickened-flame) form:**
- **(SOLV-4.4)** `∂(ρc)/∂t|_source = ∇·(ρ D_c ∇c) + ρ_u·K·c(1−c)(c−a) + ρ·(1−c)/τ_ign(p, T_u, Z)` — a
  **matched front-thickening diffusion** `D_c`, a **bistable (Nagumo/Allen-Cahn) propagation reaction**
  `ρ_u·K·c(1−c)(c−a)` (unburnt density × a rate coefficient `K` × the cubic progress shape with a declared
  sub-cell **ignition threshold** `a ∈ (0, ½)`), and an **auto-ignition term** (induction-time surface
  `τ_ign`, Arrhenius-class). **The grid-independence mechanism is the bistable pushed-front identity, sized
  explicitly:** with `S_T = S_L(p, T_u, Z)·(`wrinkling`)`, a local cell length `Δ`, and a declared front
  width `w = Θ·Δ` (`Θ` the width parameter in cells, ≈ 1.5), set `D_c = w·S_T/(1−2a)` and
  `K = 2·S_T/((1−2a)·w)`. The Nagumo travelling wave `c(ξ) = [1 + exp(ξ/w)]⁻¹` of the paired `(D_c, K)`
  system then travels at exactly `√(D_c·K/2)·(1−2a) = S_T` with width scale `√(2D_c/K) = w = Θ·Δ` (a fixed
  `Θ` cells at **every** resolution — both coefficients scale with `Δ`, so the speed is `S_T` and the width
  `Θ` cells independent of `Δ`). `D_c` (and `K`) are the front-carrier device (declared, like SRD), not a
  physics claim about flame thickness — the *speed* is the closure, the *width* is numerical; the exact
  `tanh` profile lets `flame_1d` initialize the settled front directly, so there is no pulled-front
  relaxation transient. **Why bistable and not the monostable KPP `c(1−c)` or the FSD `ρ_u·S_T·|∇c|` form:**
  (i) the FSD `|∇c|` source paired with *linear* diffusion is analytically degenerate — the linearized
  travelling-wave `D c'' + (V−S_T)c' = 0` admits **no bounded front**, so the speed is not closure-set.
  (ii) The monostable Fisher-KPP `c(1−c)` is a **pulled** front — its speed is set by the *unstable* `c → 0`
  leading edge, which is both **numerically pathological** (grid-sensitive, slow-converging) and
  **unphysical**: `c = 0` is linearly *unstable*, so any noise `c = ε` grows and the whole unburnt domain
  self-ignites. The **bistable** cubic makes `c = 0` **metastable** (`c < a` decays back to unburnt —
  physical: reactants do not spontaneously burn without a trigger) and gives a **pushed** front whose speed
  is set by the front *core*, so it is robust and grid-independent. It is the artificially-thickened-flame
  (ATF) realization of the 0.4 draft's "propagation speed at `S_T`" intent. A uniform unburnt region
  (`c ≡ 0`) is a **stable** fixed point and **cannot self-ignite** by propagation (the auto-ignition term
  seeds it where `τ_ign` is finite: the igniter kernel, compressed end-gas — pushing `c` past `a`); a uniform
  burnt region (`c ≡ 1`) is the other stable fixed point. **On the two terms adding:** in the deflagration
  regime
  `τ_ign` is orders longer than the flame's passage time, so the sum is propagation-dominated (no
  double-count in practice); where compressed end-gas approaches autoignitive states the added induction
  consumption **is** the physics (pre-ignition/knock class), and the two closures are anchored separately
  (`flame_1d` vs the ignition-delay reproduction) so each is right where it dominates.
  Both closures carry **{citation, envelope, band}** (META-3 `h2-flame-speed`, `h2-ignition-delay`,
  `turbulent-flame-speed`, `h2-kinetics-mech`) and blend continuously over `M`: outside flammability /
  below quench conditions, `S_T → 0` and `τ_ign → ∞` **by the closure surfaces themselves**, never by a
  regime branch. Heat release is *implicit in the blend* — as c rises, the cell's EOS reads the burnt branch
  (the S18 η_c\* knockdown composes unchanged: it offsets the burnt branch's interrogation coordinate).
  **The `turbulent-flame-speed` citation is now PINNED (v0.4.4; META-3: the Zimont 2000 / Peters 2000
  correlation class), and its first consumer is the S7 startup config's declared wrinkling factor** —
  config data with citation + band; the dynamic SGS-consistent closure rides the resolved-injector wave
  (plan S16).
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
- **The cold-side non-reactive floor (S7):** the cold analogue of the `BURN_COMPLETE` domain guard, with
  **two faces**. **(i) Closure-envelope face:** where the cell's `p` or `T_u` sits below the ignition
  surface's **own declared envelope floor** (cached at bind from the artifact's envelope attrs — one owner,
  the table), both rate terms read **zero by declaration**. Near-vacuum fill (~10²–10³ Pa) is thereby
  **declared no-burn** — cm-scale flames genuinely cannot propagate there — never a refusal mid-march. The
  surface's cold rows (ignition v0.3.0, envelope floor `T_u` 75 K) carry the *real* `S_L` down to the
  reactant branch's own validity floor, so the front's consumption of injector-cold (150–230 K) reactants
  is **data, not this guard**. **(ii) Reactant-representability face:** where the cell's enthalpy sits
  below the **unburnt branch's h-envelope floor** (colder than any representable reactant — its
  condensation edge), the cell is likewise non-reactive; reached in practice by **wall cooling of
  nearly-burnt cells** whose `c` parks just under the `1−BURN_COMPLETE` fixed point and therefore keeps
  closure queries live (the quench trajectory). the three terms of SOLV-4.4 decompose across COUP-3's operator classes and
  are all **fixed-order `M`-indexed lookups, no private integrator, bit-reproducible**. The **auto-ignition**
  term `ρ(1−c)/τ_ign` and the **bistable propagation reaction** `ρ_u·K·c(1−c)(c−a)` are **cell-local**
  (each reads only the cell's own state — the 0.4.2 recast retired the neighbour-reading FSD `|∇c|` form);
  the **front-thickening diffusion** `∇·(ρ D_c∇c)` is the neighbour-coupled piece — all ride the explicit
  advective composition (class `A` of the same step), which is why `∂(ρc)/∂t` as a whole is *not* purely
  cell-local. **S7 class split (of record, v0.4.4):** the **propagation reaction + front-thickening
  diffusion stay explicit class `A`** — non-stiff *by construction*: the matched `(D_c, K)` clock is
  `Δ/S_T ≫` the acoustic Δt (`S_T` far below the sound speed), so explicit is exact-tier there, behind the
  same loud positivity guard (`c` leaving `[0,1]` beyond round-off is a halt diagnosis, never a silent
  clamp). The **auto-ignition term is now the cell-local implicit class-`R` occupant** (COUP-3 §3.3),
  advanced per SDC sweep by a **fixed-structure backward-Euler node solve**: the equation is *linear in
  `ρc` at frozen `τ_ign`* (closed form, unconditionally stable), and `τ`'s weak dependence on the unknown
  (through the blend's projected `p`) is converged by a fixed `N_TAU_REFREEZE = 2` re-evaluations. The
  guarded law's burnt fixed point is honored **EXACTLY**: the source vanishes at `c ≥ 1−BURN_COMPLETE`, so
  a stiff update **parks AT `ρ(1−BURN_COMPLETE)`** — the exact integral of the declared discontinuous law,
  not a clamp. Realized (applied-increment) rates ride the SDC trapezoid quadrature and the COUP-2 audit's
  `burn_progress` row; the node-0 quadrature rate is **capped at the realizable `(cap−ρc)/Δt`** (the raw
  `ρ/τ` at extreme stiffness would overshoot the composition — measured hazard, fixed); the mild limit
  reduces to the explicit value. **Measured acceptance:** a superheated 15-bar open tube at
  `T_u = 1888 K`, `τ_ign = 1.758e-7 s`, `dt/τ = 3.04` (past the explicit positivity bound, where the S6
  tier halted) marches to `burned = 0.999` with `b` parking **bit-exactly at `1−BURN_COMPLETE`**.
  Only `ρc` is sourced — mass, momentum, and energy are untouched; the heat release is *implicit in the
  blend* (as `c` rises the same conserved `h` reads hotter/higher-`p` off the burnt branch, the two branches
  sharing one CEA formation reference, S5), so the conservation audit's energy row is unperturbed by the
  reaction. **Runtime finite-rate networks remain out of scope**: `S_L` and `τ_ign` are *offline-computed*
  surfaces (Cantera 1-D freely-propagating flames + 0-D constant-pressure reactors on the cited
  `h2-kinetics-mech`, OFFL-3 §3.3), tabulated so extinction is the surfaces' **own** smooth behaviour
  (`S_L → 0` at the flammability limits; `τ_ign` capped at a declared `τ_max` so the source vanishes at the
  quench/low-`T` corner — no `τ_ign = ∞`), validated at the unit tier against measured flame-speed and
  shock-tube ignition-delay data (VAL-2 §3.3).
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
   (coarse vs fine, THE gate); `spark_box` lights; `lean_no_light` refuses (NEVER_IGNITED); τ_ign surface
   reproduces shock-tube ignition delays and S_L the measured flame speeds within declared bands (VAL-2
   §3.3 anchors). **`quench_box` (FLAMEOUT) — REALIZED at S7 (v0.4.4; re-scoped from S6, 0.4.3 finding:
   quenching is a heat-loss phenomenon — a lit closed adiabatic box cannot flame out,
   `adiabatic_box_cannot_flame_out`):** a lit front in the 0.8 mm-gap tube between cold NoSlip+Isothermal
   walls, marched with the **full blend↔class-D gas-diffusion coupling** (flow + gas class-D + class-`R`
   scheduled together), at 0.1 atm and 420 K walls — the coldest wall the GAS-PHASE model honestly
   supports: the products surface's own envelope floor is ~407 K (envelope; grid floor ~332 K) at these pressures (H₂O condenses below
   it, measured from the v0.4.0 artifact); colder walls are the S15 two-phase wave — **DIES**: `R`
   collapses 2.06e-4 → **exactly 0 kg/s** with `burned = 0.150` (the front froze mid-tube) — the FLAMEOUT
   R-trajectory COUP-4 consumes — while the **adiabatic control on the same fixture holds `R` at 2.20e-4**
   (front alive; at 0.1 atm the laminar front crosses ~half a cell in the march window, so the
   discriminator is the live-vs-dead reacting measure, not consumption).

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
| 2026-08-24 | 0.4.5 | **Plan S7, ◆C2 shake-out — the partition extension is TWO-SIDED.** The certifying startup march exposed the hot mirror of the 0.4.4 cold-side case: a spark-kernel cell at c = 0.91 with its enthalpy past the reactant branch's h-CEILING (the deposit + burn drove h ≈ 8.7e6 vs the 8.63e6 ceiling) is physically sub-µs from fully burnt — τ_ign at the ceiling temperature says so — but the shared-h projection could not represent the state to let the reaction finish it (the S6-close "blast headroom" issue one rung up; chasing it with ever-hotter metastable ceilings is unbounded). The declared extension now pins the reactant sub-state at the NEAR EDGE of its envelope (floor under heat loss, ceiling under spark/blast superheat) with the products absorbing the balance `h_b = (h − (1−c)·h_pin)/c` — continuous at both edges, mass-consistent, self-limiting at the products' own envelope, `B_PARTITION_MIN` regularizing trace weights on both sides. The rate-law coordinates `T_u`/`ρ_u` read the PINNED reactant state, so a superheated mid-transition cell's closures stay live at the ceiling values (τ_ign sub-µs) and the class-R solve completes it within a step — the model self-heals through the reaction, which is the physically true statement. **Addendum (◆C2 cold-purge finding):** the balance form applies only where its 1/c-amplified result is representable on the products surface; otherwise the deficit is declared unattributed (the trace form) — measured: pre-light bell gas purged at ~105 K with c ≈ 0.01 amplified its sub-floor deficit ×100 past the products ceiling; the fallback's bookkeeping error is bounded by reactant-weight × off-envelope depth (sub-percent, transient purge states only). |
| 2026-08-24 | 0.4.4 | **Plan S7 — the class-`R` slot fills, the cold side gets its floors, `quench_box` flames out for real.** Five amendments, all as-built. (1) **§3.6 integration split of record:** the propagation reaction + front-thickening diffusion stay **explicit class `A`** — non-stiff *by construction*, the matched `(D_c, K)` clock is `Δ/S_T ≫` the acoustic Δt — while the **auto-ignition term becomes the cell-local implicit class-`R` occupant** (COUP-3 §3.3), advanced per SDC sweep by a fixed-structure backward-Euler node solve: linear in `ρc` at frozen `τ_ign` (closed form, unconditionally stable), `τ`'s weak dependence on the unknown (through the blend's projected `p`) converged by a fixed `N_TAU_REFREEZE = 2` re-evaluations. The guarded law's burnt fixed point is honored **exactly** — the source vanishes at `c ≥ 1−BURN_COMPLETE`, so a stiff update parks AT `ρ(1−BURN_COMPLETE)`, the exact integral of the declared discontinuous law, not a clamp. Realized rates ride the SDC trapezoid quadrature + the COUP-2 `burn_progress` audit row; the node-0 quadrature rate is capped at the realizable `(cap−ρc)/Δt` (raw `ρ/τ` at extreme stiffness overshoots the composition — measured hazard, fixed); mild limit = the explicit value. Measured: a superheated 15-bar open tube at `T_u = 1888 K`, `τ_ign = 1.758e-7 s`, `dt/τ = 3.04` — past the explicit positivity bound where the S6 tier halted — marches to `burned = 0.999`, `b` parking bit-exactly at `1−BURN_COMPLETE`. (2) **The cold-side non-reactive floor**, two faces: below the ignition surface's own declared envelope floor (cached at bind from the artifact's envelope attrs — one owner, the table) both rate terms read zero by declaration — near-vacuum fill (~10²–10³ Pa) is declared no-burn, never a mid-march refusal, while the surface's cold rows (ignition v0.3.0, floor `T_u` 75 K) carry the *real* `S_L` down to the reactant branch's validity floor so injector-cold (150–230 K) consumption is data, not the guard; and below the unburnt branch's h-envelope floor (colder than any representable reactant) the cell is likewise non-reactive — reached by wall cooling of nearly-burnt cells whose `c` parks under the fixed point and keeps closure queries live (the quench trajectory). (3) **Cold-side partition extension:** the shared-`h` flamelet identity presumes adiabatic burning; under real heat loss the declared extension pins the reactant sub-state AT the unburnt h-floor with the products absorbing the balance `h_b = (h − (1−c)·h_floor)/c` — continuous, mass-consistent, self-limiting (the burnt envelope is the true cold edge) — with a trace-weight floor `B_PARTITION_MIN = 0.01` below which the products read the cell `h` directly (the balance form amplifies sub-floor deficits by `1/c`). Declared model forms, never clamps. (4) **§6.5 `quench_box` REALIZED:** 0.8 mm gap, cold NoSlip+Isothermal 420 K walls (the products surface's ~407 K (envelope; grid floor ~332 K) envelope floor at these pressures is the honest gas-phase limit — H₂O condensation below it is S15), 0.1 atm, full blend↔class-D coupling: `R` 2.06e-4 → exactly 0 kg/s, `burned = 0.150` (front froze mid-tube) — THE FLAMEOUT trajectory — vs the adiabatic control alive at 2.20e-4. (5) **`turbulent-flame-speed` PINNED** (META-3: Zimont 2000/Peters 2000 correlation class); first consumer = the S7 startup config's declared wrinkling factor (config data, cited, banded; dynamic SGS-consistent closure rides plan S16). |
| 2026-08-24 | 0.4.3 | **Plan S6 close — the burnt fixed point owns ONE threshold; the spark is a bounded pulse.** Three S6-close consistency fixes, found by the `spark_box`/`quench_box` mini-sims and landed with the code. (1) **The blend's pure-BURNT threshold is derived from the reaction's burnt fixed point** — `EPS_B_PURE_BURNT = 2·BURN_COMPLETE` (one owner, `BURN_COMPLETE` = 10⁻³). The reaction zeroes its source (and its `S_L`/`τ_ign` queries — the §3.6 domain guard: post-flame gas is not a reactant) at `c ≥ 1 − 10⁻³`, so `c` **asymptotes toward that point from below and never crosses it** (measured: pinned ~2×10⁻⁷ under it); the blend's pure-burnt skip sat at `1 − 10⁻⁹`, so every burnt cell kept interrogating the unburnt branch at ~10⁻³ weight forever — and hot burnt gas, whose enthalpy legitimately exceeds the unburnt surface's ceiling (the OFFL-3 0.6.2 headroom ruling's other face), refused on a branch describing 0.1 % of its mass, gas the reaction itself had declared complete. An **equal** threshold is a knife edge (the attractor lands epsilon on its wrong side — measured before the factor 2 was added), so the pure-burnt region **strictly contains** the attractor. **The pure-UNBURNT threshold stays at the original `1e-9` — deliberately asymmetric (review-wave finding):** the reaction pins an attractor only at the burnt end (the metastable fringe decays smoothly through every value to 0), and the crossing step of each skip scales with the *dropped branch's* specific volume — burnt-side `EPS·(v_u/v_b)` ≈ 2.6×10⁻⁴ (the dropped unburnt branch is dense: fine), but cold-side `EPS·(v_b/v_u)` ≈ 7.8·EPS at flame states, so a symmetric 2×10⁻³ threshold would step density ~1.5 % — *outside* the unburnt density column's own ~0.8 % bound. Both declared steps now sit far inside their columns' bounds. (2) The **reacting-measure and consumption-rate diagnostics carry the same domain guard** (they queried `T_u`/closures on fully-burnt cells; their `c(1−c)`-weighted contribution there is ≤ 10⁻³ and is now exactly 0). (3) **The igniter mini-sims obey the spark discipline** (COUP-7 0.4.1, Ben ruling): the deposit is a **bounded pulse ending at about the ignition time** (~0.7–8 J, 20 µs, `spark-igniter-class`) — a sustained deposit into an already-burnt kernel superheats the mid-transition cells past the unburnt/ignition envelopes (the metastable-reactant validity edge; the hot ASI-torch regime is declared S7 hardening territory, with the stiff class-`R` treatment it already owns). **(4) `quench_box` re-scoped to S7 (finding):** a lit **closed adiabatic** box cannot flame out — with no loss channel the pressure rise compression-heats the reactants and auto-ignition completes the burn, so the mini-sim as drafted asked for unphysical behavior and the model refused to fake it. §6.5 updated: S6 asserts the converse (`adiabatic_box_cannot_flame_out`); the flameout demonstration = conductive loss to cold walls via the **blend↔class-D coupling S7's startup march builds regardless** (the fixture's 0.8 mm gap is at the quench-distance scale — though H₂/O₂'s own quenching distance is a few× smaller than H₂/air's ~0.6 mm, so the S7 box may need lower pressure or a narrower gap to actually flame out). Until then no test exercises the FLAMEOUT (R-collapse-after-exceed) trajectory — its consumer, the COUP-4 verdict object, is also S7, so coverage lands with its consumer. |
| 2026-08-21 | 0.4.2 | **Plan S6 — §3.6 built (blend + rate law + igniter + halts), landed with the code.** Three amendments, code-driven: (1) **the sub-cell partition is now the energy-conserving flamelet form** — both sub-states share the cell pressure `p` **and enthalpy** `h = e + p/ρ` (the adiabatic constant-`p` identity `h_b = h_u = h` on S5's shared CEA reference: burning conserves enthalpy, so the two branches are read at the *same* `(p, h, Z)`, differing only in composition), density partitioned by **mass-weighted specific volume** `1/ρ = (1−c)/ρ_u + c/ρ_b` — the SOLV-1 §3.4 projection generalized to both branches, recovering `c=1` (bit-for-bit shifting) and `c=0` (pure unburnt) exactly. This **supersedes the 0.4 injection-isentrope `h_u(p; h_inj, Z)` closure**, which is the *unheated* special case `h = h_isentrope(p)`: it could not fire `τ_ign` on a non-isentropic igniter deposit (the deposit raises `h` but not the isentrope). There is now **no `h_b = (…)/c` division and so no `c→0` singularity** (the earlier draft's `EPS_C_SPLIT` regularization is unnecessary and was dropped). `T_u = T_unburnt(p, h, Z)` sees the deposited enthalpy, so the igniter works; declared mass-weighted blended sound speed + diagnostic temperature. (2) **SOLV-4.4 operator-class decomposition + integration scoping** — the auto-ignition term is the cell-local class-`R` piece, propagation/diffusion are neighbour-coupled; **S6 advances all three explicitly** (matched non-stiff `D_c`, long deflagration `τ_ign`) behind a **loud positivity guard** (never a clamp), with the stiff class-`R` implicit auto-ignition a **declared S7 hardening**; only `ρc` is sourced (heat release EOS-implicit via S5's shared reference, energy audit unperturbed). (3) **extinction pinned to the surfaces** — `S_L→0` at flammability limits, `τ_ign` capped at a declared `τ_max` (no `∞`), so no-light/flameout are the closures' own smooth behaviour feeding the COUP-4 halts. The igniter is the existing scheduled `S(r,θ,z,t)` energy deposit (COUP-7 §3.3), ramped. **SOLV-4.4 also re-cast from the FSD `ρ_u·S_T·|∇c|` propagation form to the well-posed bistable (Nagumo) reaction `ρ_u·K·c(1−c)(c−a)` with the matched `(D_c, K)`** — the `|∇c|` form paired with linear diffusion is analytically degenerate (no bounded travelling-wave front ⇒ the speed is *not* closure-set, the exact failure `flame_1d` catches), and the **monostable Fisher-KPP `c(1−c)` is a pulled front** (speed set by the *unstable* `c→0` edge — grid-pathological AND unphysical, since noise self-ignites the whole unburnt domain). The **bistable cubic** makes unburnt `c=0` **metastable** (physical ignition threshold `a`) and gives a **pushed** front at `√(D_c·K/2)·(1−2a)=S_T`, width `Θ·Δ`, robust + **provably grid-independent** — the ATF realization of the 0.4 draft's "propagation speed at `S_T`" intent (initialized directly at its exact `tanh` profile, so no relaxation transient). |
| 2026-08-20 | 0.4.1 | **Plan S5 — the c = 0 branch is shipped (dependency note, no design change).** §3.6's blend `h = (1−c)h_u + c·h_b` reads `h_u` from the **unburnt-reactant surface**, which OFFL-3 §3.3 now ships (a gas-phase ideal-gas frozen reactant mixture vs (p, h, Z), cryo-valid, on the equilibrium surface's FND-5 schema and enthalpy reference). The blend, the TFC rate law, the sub-cell partition, and the igniter object are unchanged and remain S6 build content; only the surface `h_u`/`T_u` are keyed on is now a real artifact rather than a contract row. |
| 2026-07-21 | 0.1 | Initial draft. One `ReactionSource` (fission/fusion/annihilation differ only in rate law + birth spectrum). Advected-precursor fission kinetics with point kinetics as the static-fuel reduction + source-driven subcritical (NSWR/ICAN); no private integrator (COUP-3 owns it; CRAM≠PKE). Fusion ⟨σv⟩ sourcing emitting birth spectra (p-¹¹B 3α continuum, D-D secondaries as the one rate functional on SOLV-3's slowing-down f, thermal-broadened neutron peaks). Annihilation via OFFL-4 with the factor-of-several band. Emit-once/transport-once; no private transport. |
| 2026-08-13 | 0.2 | R3 terminology: annihilation model-form reframed from a blanket "factor-of-several band" to OFFL-4's **per-quantity physics-list spread**; efficiency/waste-heat noted as computed downstream by SOLV-2/3 transport, not emitted here. |
| 2026-08-19 | 0.4 | **VISION_SCOPE v1.5 (Ben): the chemical burn-progress source (§3.6).** c ∈ [0,1] burnt-fraction field on `U`; blended unburnt↔equilibrium thermochemistry; one continuous rate law (SOLV-4.4: flame propagation via `S_L`×wrinkling + auto-ignition via `τ_ign`), closures offline-generated (OFFL-3) and unit-anchored (VAL-2); igniter = COUP-7 energy-deposit object; grid-independent front speed = the acceptance gate; reacting-measure thresholds feed the new COUP-4 `NEVER_IGNITED`/`FLAMEOUT` halts; class-`R` cell-local stiff source; runtime finite-rate networks remain out. §1/§2/§6 updated. *Same-day adversarial-review fixes:* SOLV-4.4 restated in **TFC form** (`ρ_u·S_T·|∇c|` + the **matched front-thickening diffusion `D_c`** as the explicit grid-independence mechanism — a bare source self-sharpens and its speed goes grid-set); **sub-cell partition rule stated** (common p, mass-weighted h, unburnt-isentrope closure for `h_u`/`T_u` — no second energy field); the flame/auto-ignition overlap addressed (deflagration regime is propagation-dominated; autoignitive end-gas addition is physics, each closure anchored where it dominates). |
| 2026-08-14 | 0.3 | **Review fix wave (N5, N6, E-2).** §3.2: quasi-static closure written explicitly — amplitude ODE (SOLV-4.2) + tabulated flux/adjoint shapes ψ/ψ† from OFFL-1 + adjoint-weighted delayed source (SOLV-4.3) and drifted-precursor loss ledger (flowing-fuel β_eff emerges from advection + importance weighting, configuration-general); kinetics declared a **class-`G` global stiff ODE** in COUP-3's SDC step + `EPS_RHO_STEP`/`EPS_POWER_STEP` Δt limiters; kernels' control-state axis noted (E-2). |
