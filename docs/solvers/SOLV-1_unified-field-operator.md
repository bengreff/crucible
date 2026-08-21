# SOLV-1 — Unified Field Operator

| Field | Value |
|---|---|
| **ID** | SOLV-1 |
| **Family** | SOLV (Runtime unified-grid operators) |
| **Status** | Reviewed (2026-08-14) |
| **Depends on** | FND-2, FND-7 (EOS/transport), COUP-3 (time integration), COUP-2 (audit); FND-1, COUP-7 |
| **Version** | 0.4.3 (2026-08-20 plan S5: frozen-mode `{ρX_k}` state widening explicitly rides **plan S5b**; the c = 0 unburnt (p, h, Z) branch is shipped and bound by the existing `TableEos` — no widening — §3.4). 0.4.2 (S4 review: the film-mean-vs-local c_p ruling; the reference-state band claim withdrawn; catalycity stated). 0.4.1: `F_visc`'s species-enthalpy flux `Σ h_k j_k` named; per-cell spine transport + face-averaging rule; the wall law's transport operands are spine queries |
| **Skeleton/complete split** | **Fixed now:** the conserved system, the reacting-flow update, the N_θ=1 axisymmetric-with-swirl corner, EOS/transport from the spine, wall-heat coupling (the one local wall-function law, §3.5), halts, oracles. **Deferred (W4):** the two-phase (drift-flux) and two-fluid/MHD (HLLD + constrained-transport) extensions of the *same* operator; the **resolved-mixing rung** (R2, §3.4) — unmixed multi-stream injection with grid-computed mixing (a resolution/N_θ lift, riding the two-phase extension for liquid injection) under the **one universal LES-class subgrid mixing closure** (D-D, §3.4). |

---

## 0. Purpose

SOLV-1 is the **one conserved-variable field operator** on the grid: it advances mass, species, momentum,
energy (and, in later waves, magnetic field and radiation moments) as a single hyperbolic-plus-parabolic
update over the medium-state vector `M` (Rule 12). The W2 build is its **reacting-flow corner at N_θ=1** —
2-D axisymmetric with swirl — which is exactly the N_θ=1 ring of the adaptive-azimuthal-resolution solver
(FND-2 §3.4, D-A), and the vehicle for the RL10 anchor. Two-phase and MHD are **extensions of this same
operator**, not new solvers (§1.2).

Read after FND-2 (the grid/`U` it evolves), FND-7 (EOS/transport over `M`), COUP-3 (the SDC-IMEX schedule it
runs inside), and COUP-7 (the boundary inflow that drives it).

## 1. Scope & razor ruling

### 1.1 Owns
The conserved-flux update on the grid: reconstruction, the Riemann flux, the geometric (axisymmetric) source
terms, the EOS/transport closure calls into the spine, and the reaction-source intake. It writes `U` in place.

### 1.2 Rule-12/13 compliance
- **One operator, extended by `M`, not branched.** Reacting flow, two-phase, and two-fluid/MHD are the **same
  conserved system** with more active components of `U`/`M`; they are added by widening the state and the
  spine, **never** by an `if(regime)` solver switch (VISION_SCOPE §5.1 banner). W2 activates the
  reacting-gas subset; the drift-flux and HLLD/CT terms are dormant (their `M`-components null), not absent.
- **No per-species code.** Species are advected scalars `ρX_k`; the EOS is a general convex law over `M`
  from the spine (FND-7) — a new propellant is data (composition + spine), not code (Rule 13).

### 1.3 Defers (owner named)
Time integration / operator split → **COUP-3**; the conservation audit + operator-coupling contract →
**COUP-2**; EOS/transport/opacity/stopping over `M` → **FND-7**; reaction birth/heat sources → **SOLV-4**
(chemical: tabulated equilibrium composition from **OFFL-3**, §3.4); radiation transport → **SOLV-2**;
geometry/apertures → **FND-3/FND-2**; boundary inflow (ṁ, MR) → **COUP-7**.

**Razor ruling:** SOLV-1 *is* the reaction — anything the reaction touches (combustion, wall heating, phase
change, ablation) is full grid physics here (VISION_SCOPE §4.1). Turbulent combustion CFD is **not** resolved
(replaced by tabulated equilibrium chemistry + the tiered injector integration — envelope-bounded η_c\*
prior, or grid-resolved mixing under the one universal subgrid closure, §3.4).

## 2. Interfaces & contracts
| Interface | Consumed by | Contract |
|---|---|---|
| **Conserved update** `U ← U + Δt·(flux + sources)` | COUP-3 | given a state and Δt, return the flux-divergence + geometric source contribution; stiff sources are handed to COUP-3's implicit stage |
| **Flux at a face** | COUP-2 | the HLLC numerical flux per face — the quantity the conservation audit telescopes |
| **EOS/transport query** | FND-7 | p, T, sound speed, μ, k from `M` — no material label read; also the wall-function operands (§3.5) |
| **Reaction-source intake** | SOLV-4 / OFFL-3 tables | composition/heat-release as **local-state (p, h, elemental-composition) lookups** (OFFL-3 §3.3, S22) |
| **Wall-function flux at wall faces** | COUP-2 (exchange sweep), COUP-7 (coolant side) | h from local near-wall state per §3.5; `F_visc` suppressed at those faces (no double-count) |

**Invariant:** the update is **flux-form** (conservative by construction, telescoping to port fluxes for
COUP-2); it reads `M` and writes `U` **in place** on the one grid; and it is **deterministic** (fixed
reconstruction/flux/reduction order, FND-2 §3.7).

## 3. Method & governing equations

### 3.1 Conserved system
Per cell, `U = (ρ, {ρC_j}, ρ𝐮, ρE)` for the W2 reacting-gas subset (the full `U` also carries `𝐁_face,
E_rad, 𝐅_rad`, dormant here — FND-2 §3.4). The composition block `{ρC_j}` is **mode-dependent (S19)**: in
**shifting-equilibrium** mode it is the **elemental mass fractions** (advected, source-free — element
conservation is exact by construction); in **frozen** mode it is the full species set `{ρX_k}` with zero
reaction source (§3.4). Governing law (compressible reacting Navier-Stokes, source form):
- **(SOLV-1.1)** `∂U/∂t + ∇·F(U) = ∇·F_visc(U,∇U) + S_geom + S_react + S_wall`, with `F` the hyperbolic flux,
  `F_visc` the viscous/conductive flux (transport from the spine), `S_geom` the axisymmetric geometric source
  (§3.3), `S_react` the reaction source (§3.4), `S_wall` the wall-heat coupling (§3.5).

**`F_visc`'s energy flux is `τ·u + k∇T + Σ_k h_k j_k`** *(0.4.1, plan S4)* — the third term is the enthalpy
the diffusing species carry, and omitting it is not a small error when Lewis ≠ 1: it is the difference between
conducting heat and *transporting* it with the mixture. Under the one composition coordinate of §3.4's
shifting mode it reduces exactly to `(∂h/∂Z)|_{p,T}·j_Z` with the coefficient a spine output (FND-7 §3.3), so
it costs one more face term and no new model. It vanishes identically on a single-composition gas — which is
why the S3 constant-transport occupant could defer it honestly. **Transport is a per-cell spine query on the
cell's own local state** (FND-7 §3.7's (p, h, Z)); face coefficients are the arithmetic mean of the two
cells', which is exact in the constant-coefficient limit and second-order for the smooth transport fields a
gas has — there is no material discontinuity to cross, because a wall face carries no resolved diffusion at
all (§3.5).

### 3.2 Discretization — Godunov FV, PPM + HLLC
Finite-volume Godunov: **PPM (or PLM) reconstruction with characteristic tracing**, **HLLC** numerical flux
with **Batten wavespeed estimates** (restores the contact + species-advection waves HLLE smears — mandatory
with multiple species). General convex EOS in the Riemann solve (no ideal-gas assumption). Explicit
hyperbolic part; viscous/conductive and reaction parts are the stiff pieces COUP-3 treats implicitly.
[META-3: `hllc`, `castro-source`, `maccormack-nozzle` (didactic baseline)]

### 3.3 Axisymmetric with swirl — the N_θ=1 corner
Solve on the cylindrical grid at azimuthal resolution N_θ=1 with **well-balanced geometric source terms** so a
uniform/steady state is preserved exactly at the axis; carry the **azimuthal (swirl) momentum as an advected
out-of-plane component** ("2.5-D") with its centrifugal term in `S_geom`. This is precisely the **N_θ=1
axisymmetric-with-swirl corner** of the adaptive azimuthal resolution (FND-2 §3.4, D-A) — no separate 2-D
solver, just the one operator with its ring resolution pinned at N_θ=1. The **r=0 axis** uses a reflecting BC
(Castro-style 2.5-D), checked by an axis-symmetry conservation test (§6). Higher azimuthal resolution
(N_θ>1) is a capability lift of the *same* finite-volume code — more rings of the same cells and fluxes,
never spectral modes. [META-3: `castro-source`]

### 3.4 Reactions for the chemical slice — tabulated equilibrium, not finite-rate
Combustion is **tabulated shifting-equilibrium chemistry** (OFFL-3), interrogated at **local state** — the
(p, h, elemental-composition/mixture-fraction) surfaces of OFFL-3 §3.3 (S22) — **not** a runtime finite-rate
network (turbulent combustion CFD is out, VISION_SCOPE §5.3).

**Composition modes (S19).** Config-selected, mirrored in `U` (§3.1): **shifting-equilibrium** mode advects
**elemental mass fractions** and applies an **equilibrium projection each step** — species, temperature, and
heat release are the OFFL-3 equilibrium surface evaluated at the cell's local (p, h, Z); element conservation
is exact because only elements are advected. **Frozen** mode advects the full species set with zero reaction
source. The frozen↔shifting gap is a **discrete epistemic model-form dimension** sampled by COUP-5 (the
JANNAF kinetic-efficiency deviation becomes an error bar, not a hidden choice).

**Frozen-mode state widening rides plan S5b** *(0.4.3, plan S5).* The frozen mode's `{ρX_k}` composition
block is a **wider `U`** than shifting mode's single advected element — a change to the conserved-state width
`NCOMP`, which the fixed-array runtime carries as a compile-time constant. Ben split that widening into its
own session (S5b): it is bit-identity-critical for the shifting stations (the `c ≡ 1` / single-element
corner must stay byte-for-byte), has **no consumer before the COUP-5 ensemble wave** (S18 — the delivered-
performance frozen end is carried by the OFFL-3 §3.2 *bracket* until then), and is orthogonal to S5's actual
mission, the cold/unburnt branch. **The c = 0 unburnt branch needs no widening**: it is a single-element
(p, h, Z) surface (a two-stream reactant mixture is set by Z), shipped at S5 (OFFL-3 §3.3) and bound by the
**existing** `TableEos` occupant with no new state — so a cold non-reacting flow is representable now, and
S6's burn-progress blend (§3.4 below) composes the two (p, h, Z) branches without a species vector.

**Burn progress (v0.4, VISION_SCOPE v1.5).** `U` additionally carries the **burn-progress field c ∈ [0,1]**
(cell burnt mass fraction): the cell's thermochemistry is the SOLV-4 §3.6 blend of the **unburnt**
(frozen-reactant) and **burnt** (equilibrium) branches, and c's source — flame propagation + auto-ignition,
via universal banded closures, with the igniter as a COUP-7 energy-deposit object — is **owned by SOLV-4
§3.6** (a reaction source, exactly the §1.3 deferral); SOLV-1 only advects c and evaluates the blended EOS.
**With c active, the per-step equilibrium projection above applies to the burnt fraction only** — the cell
state is the c-blend of the unburnt-branch state and the projected equilibrium state (partition rule owned
by SOLV-4 §3.6); `c = 1` recovers this section's shifting mode identically, so pure-shifting configs (the
stations, the session-12 certificate) are the c ≡ 1 corner, bit-unchanged. Ignition, failure-to-ignite,
and flameout thereby become computed field behavior (COUP-4 halts).

**Injector tiers (R2, COUP-7 §3.2.1).** Injector-scale mixing/vaporization enters at the config-selected tier:
- **Prior tier** (W2/milestone-1): premixed inflow at the declared MR with the cited η_c\* prior applied as an
  **in-solver combustion-completeness knockdown at the source-term level** (S18): the effective heat release
  of the equilibrium projection is scaled so that delivered c\* = η_c\*·c\*_ideal at the anchor state — so T_c,
  wall heat flux, jacket ΔH, and the emergent p_c are all self-consistent with the deficit. **Output-side
  multiplication of c\*/Isp by η_c\* is forbidden** (it leaves the field state inconsistent with the reported
  performance). COUP-7 §3.2.1 carries the mirror statement on the boundary-object side. Legitimate only
  inside the cited injector-family envelope.
- **Resolved tier** (deferred lift of the *same* operator): the injection plane states the actual unmixed
  streams — the composition advection already carries them, nothing new enters `U` — mixing is computed on
  the grid (azimuthal resolution N_θ>1 and/or interface refinement where the mixing indicator demands it,
  FND-2 §3.4; liquid streams ride the two-phase extension for atomization/evaporation via universal
  closures), and η_c\* is an **output** (SOLV-7 reads c\*_delivered/c\*_ideal).

**Declared closure of the resolved tier (D-D, S20).** Mixing at engine Reynolds numbers (~10⁷) is turbulent;
a laminar-resolved calculation **under-mixes by orders of magnitude** and would discharge the R2
"computed, not assumed" promise falsely. The resolved tier therefore carries **one universal LES-class
subgrid mixing model** — dynamic-coefficient eddy-viscosity + gradient-diffusion species mixing with a
turbulent Schmidt number, coefficients from the dynamic procedure (**no hand-tuned constants**) — as
deferred-build content of this same operator: an added `F_visc`-class subgrid flux, never a second solver.
It carries a **per-quantity model-form band** and is calibrated **offline against canonical turbulence data**
— universal physics, never per-engine (the tool doctrine). [META-3: `les-dynamic-sgs`,
`canonical-turbulence-data`]

The tier is config data, never a code branch (Rule 13). At the resolved tier local composition leaves the
global-MR line — covered natively by the (p, h, Z) table coordinates (OFFL-3 §3.3).
[META-3: `jannaf-eff`, `injector-cstar-eff`]

### 3.5 Wall heat and emergent chamber pressure
**One local wall-function heat-flux law — every wall, every engine (D-C).** The gas-side convective wall flux
is a **local Reynolds-analogy/Colburn-class wall function**: h computed from the **local near-wall state
only** — ρ, tangential velocity, T, T_wall, μ, k, Pr (spine transport, FND-7), and wall distance — with a
declared **±20–30% band**. **The transport operands are spine queries at the near-wall gas cell's own state**
*(0.4.1, plan S4)*, so the law and the resolved `F_visc` next to it read the **same one provider** — the wall
law states no transport constant of its own.

**Which c_p, and the reference-state question** *(0.4.2, S4 review — this replaces a claim that was measured
false).* The Colburn analogy transports **enthalpy**: `q_w = St·ρu·(h_aw − h_w)`. Writing that as `h·ΔT`
requires the **film-mean** slope `(h_aw − h_w)/(T_aw − T_w)`, not the local one. For a dissociating gas the
local equilibrium c_p is the *peak* of a strongly-peaked curve — ~7970 J/(kg·K) against a film mean of ~4130
at the RL10 chamber — so driving on it overpredicts `q_w` by **1.7–2.4×**, one-signed, right through this
band. The law therefore drives its convective limb on the spine's **c_p,frozen** (a measured 0.94–1.06 proxy
for the film mean at chamber/throat states, 0.76–0.86 below ~0.2 MPa) and its **recovery** term on the local
c_p, which is an edge-state property. *Recorded deferral:* the exact form drives on `h_aw − h_w` directly —
the spine is already keyed on h and the EOS can invert `T(p,h,Z)` at the wall temperature, so it is a
root-solve per wall patch, not a new model; it rides the mount-reaction wave with the skin-friction debit.

The earlier 0.4.1 wording said that evaluating the operands at the gas state rather than at a film or Eckert
reference temperature was "a declared choice **inside** the ±20–30% band". **That is withdrawn**: with
constant properties it was vacuously true, but with state-dependent properties the measured lever is
h(film)/h(gas) ≈ 0.24–0.40 — a factor 2.5–4, not a ±30%. The reference-state choice is now a **named
model-form limit** of this closure, not a band member: at a resolved near-wall tier the wall-adjacent cell's
state approaches T_wall and the operand set shifts systematically. It is declared here rather than folded
into the band, and it is the reason the ±20–30% must not be read as covering resolution changes.

**Wall catalycity** *(0.4.2)*: this law's operands are the equilibrium (fully-catalytic) ones, which is the
upper bound on `q_w`; the resolved `F_visc` next to it declares a **non-catalytic** species wall
(zero-flux). They never meet on the same face — `F_visc` is suppressed where the wall function runs — so
there is no contradiction in the discretization, but the two assumptions are different and both are now
stated rather than inferred. Because its operands are purely local `M` + local geometry, the same law is valid
at every wall face in every engine — nozzle, reactor channel, duct — with no per-geometry correlation branch
(Rule 12; capabilities are general, never engine-specific). **Bartz is demoted to a VAL-2 nozzle-envelope
oracle** — a cross-check on the integrated nozzle heat load, no longer a runtime closure. [META-3:
`wall-function-heat`, `bartz` (oracle)]

**Ownership delineation (S16, one owner per quantity):** **SOLV-1** evaluates the wall function at
**config-time-identified wall faces** — geometric identification from apertures/volume fractions (data, not
an `if(material)`); the resolved viscous/conductive flux `F_visc` is **suppressed at those faces** (the wall
function replaces, never adds to, the resolved flux — no double-count). **COUP-2** owns the exchange-sweep
placement (which sub-step, the Robin-BC iteration); **COUP-7**'s cooling-jacket object owns only the
**coolant side**.

This wall flux (a) sets the liner thermal state and (b) **is the expander-cycle drive**: its integral is the
jacket enthalpy rise that COUP-7's closed-mode turbopump converts to delivered ṁ. **Chamber pressure is
emergent** — the injected ṁ combusts and chokes at the throat, and `p_c` is whatever the field produces (read
out by SOLV-7), never imposed.

### 3.6 Error budget & determinism
2nd-order in smooth regions (PPM/PLM), 1st-order at captured shocks; the grid/step at which truncation drops
below the physics band is demonstrated by a convergence test (VAL-1 rung i). f64 throughout; fixed
reconstruction/flux/reduction order → bit-reproducible at any thread count (FND-2 §3.7). Small-cut-cell
stability is handled by **State Redistribution** (FND-3/FND-2). [META-3: `state-redistribution`]

## 4. Coupling relationships
- **COUP-3** wraps the update in the SDC-IMEX step (explicit hyperbolic here; conduction/viscous fluxes in
  the **spatially-coupled implicit diffusion class**, COUP-3 §3.1 — never a cell-local conduction solve;
  cell-local implicit reserved for reactions/M1 source coupling); anchor runs are **physical marches**
  (VISION_SCOPE v1.5 — the former pseudo-transient accelerator is deleted, COUP-3 §3.6 tombstone);
  **COUP-2** audits the flux-form conservation and enforces
  operator-coupling (incl. the wall-exchange sweep placement, §3.5).
- **FND-7** supplies EOS/transport over `M`; **OFFL-3** supplies the equilibrium/frozen composition tables;
  **SOLV-2** (radiation) is dormant in the chemical slice (SOLV-2 §3.4).
- **COUP-7** injects ṁ/MR and consumes the jacket ΔH; **SOLV-7** integrates the exit plane for
  thrust/Isp/c\*; **SOLV-6** reads wall p/T for structural margins; **SOLV-8** (stage 2) consumes recession.

## 5. Uncertainty & validity
Introduces **discretization (truncation) error** (budgeted, feeds UQ) and consumes the spine's EOS/transport
band, the wall-function band (§3.5), and the reaction frozen↔shifting dimension (§3.4). Validity:
resolution-limited (a feature thinner than the finest cell is flagged, FND-2 §5); the reacting-gas subset is
valid for single-phase combustion flow — the two-phase/MHD corners are dormant, not silently approximated.
Ladder rungs: analytic (Sod, nozzle), reference code (Castro/PeleC), hardware (RL10 via SOLV-7/VAL-2).

**Anchor-run budget (S21, rewritten 2026-08-19 — VISION_SCOPE v1.5).** The certifying RL10 anchor is a
**full-3-D physical march** from declared fill state through ignition to settled steady state: grid class
**~5×10⁶ cells** (static-refined (r,z) × adaptive N_θ ≤ 64), Δt ~3×10⁻⁷ s (acoustic CFL), a **declared
compressed start window of 100–200 ms** (COUP-7 §3.2.2 — external schedules only) ⇒ **~10⁶ steps ≈ 10¹³
cell-updates**. Throughput: the CPU fixed-order build is the bit-exact reference and the mini-sim tier
(~10⁷–5×10⁷ cell-updates/s); certification runs live on the **GPU port at ~1–5×10⁸ f64 cell-updates/s
(RTX 4080 class, bandwidth-bound; measured, not assumed — plan S12)** ⇒ **8–30 h, inside the declared 24-h
cap** with checkpoint/restart (FND-6). Ensembles are **multi-fidelity** (COUP-5 §3.2: pinned-N_θ-ceiling
ladder — the bulk of members at coarse-3-D/axisymmetric rungs, few at full resolution), never
N × full-resolution. Budget details + laptop/desktop tiers: `PLAN_CHEMICAL_SANDBOX.md` §3.

## 6. Validation plan
1. **Sod / Sedov:** exact-Riemann shock/contact/rarefaction reproduced; L1 error → 0 at formal order.
   [META-3: `sod-shock`]
2. **MMS:** manufactured-solution source recovers the formal order over every term (VAL-3). [META-3: `mms`]
3. **Axis symmetry:** a symmetric transient stays symmetric across r=0; conservation holds at the axis.
4. **Quasi-1-D nozzle:** reproduces the analytic isentropic/area-Mach relation; feeds SOLV-7's C_F check.
5. **Differential oracle:** matches Castro (hydro + reaction coupling) and Athena++ (flux/geometry) to their
   truncation error on shared cases; PeleC for a real H₂/O₂ case. [META-3: `castro-source`]
6. **Determinism:** bit-identical at 1 vs N threads and across runs.
7. **Anchor budget held (S21):** the RL10 anchor meets the §5 budget (physical-march step count under the
   declared compressed schedule; measured GPU throughput) — recorded in the results bundle.

## 7. References
META-3 keys: `hllc`, `castro-source`, `maccormack-nozzle`, `state-redistribution`, `jannaf-eff`,
`injector-cstar-eff`, `wall-function-heat`, `bartz` (VAL-2 oracle), `les-dynamic-sgs`,
`canonical-turbulence-data`, `sod-shock`, `mms`. Depends on FND-2, FND-7, COUP-3, COUP-2, COUP-7,
SOLV-4, OFFL-3.

*(No open questions — chemical-slice reaction handling (tabulated equilibrium + frozen↔shifting UQ band) and
emergent-`p_c` design resolved with COUP-7/SOLV-7, Ben 2026-07-21; injector mixing tiered per R2, Ben
2026-08-13.)*

## 8. Change log
| Date | Version | Change |
|---|---|---|
| 2026-08-20 | 0.4.3 | **Plan S5 — the cold/unburnt branch (landed with the code).** §3.4: the frozen mode's `{ρX_k}` composition-block **state widening** (a `NCOMP` change to the fixed-array conserved state) is stated to ride **plan S5b**, the session Ben split it into — bit-identity-critical for the shifting stations, no consumer before the COUP-5 wave (S18; the frozen delivered-performance end is the OFFL-3 §3.2 bracket until then), orthogonal to S5's mission. The **c = 0 unburnt branch needs no widening**: a single-element (p, h, Z) frozen-reactant surface (OFFL-3 §3.3) bound by the **existing** `TableEos` occupant, so a cold non-reacting flow is representable at S5 and S6's blend composes the two (p, h, Z) branches. No operator/flux/EOS change — the shifting stations are byte-for-byte unchanged. |
| 2026-08-20 | 0.4.2 | **S4 review corrections (same session).** §3.5: the wall law's convective limb must drive on the **film-mean** enthalpy slope, not the local equilibrium c_p — measured, the latter overpredicts `q_w` by 1.7–2.4× one-signed through the ±20–30% band (the pre-S4 constant c_p = 5000 was accidentally *inside* it at 1.21×, so S4 made the property more accurate and the flux less so until this landed). The law now takes the spine's c_p,frozen for that limb and the local c_p for recovery, with the exact `h_aw − h_w` form a named deferral. **The 0.4.1 sentence claiming the gas-state-vs-reference-temperature choice sits inside the ±20–30% band is WITHDRAWN** — measured h(film)/h(gas) ≈ 0.24–0.40; it is now a declared model-form limit, and the band explicitly does not cover resolution changes. Wall catalycity stated (equilibrium operands in the law, non-catalytic species wall in `F_visc`; they never share a face). |
| 2026-08-20 | 0.4.1 | **Clarification wave, landed with plan S4's code (the COUP-3 0.4.1/0.4.2 pattern).** §3.1: `F_visc`'s energy flux stated in full as `τ·u + k∇T + Σ_k h_k j_k` — the species-enthalpy term is the Lewis ≠ 1 content, reduces exactly to `(∂h/∂Z)\|_{p,T}·j_Z` under §3.4's one composition coordinate with the coefficient a spine output, and vanishes on a single-composition gas (which is what made the S3 deferral honest); **transport is a per-cell spine query on the local state**, face coefficients the arithmetic mean (exact in the constant limit, 2nd order for smooth gas transport — no material discontinuity crosses a resolved face because wall faces carry no resolved diffusion). §3.5: the wall function's μ/k/Pr are **spine queries at the near-wall gas state** — one provider shared with `F_visc`, no transport constant stated by the law; the gas-state-vs-reference-temperature choice is declared inside the existing ±20–30% band, refinement deferred. |
| 2026-08-19 | 0.4 | **VISION_SCOPE v1.5 (Ben).** §3.4: `U` gains the **burn-progress field c** (blended unburnt↔equilibrium thermochemistry; rate laws owned by SOLV-4 §3.6; igniter = COUP-7 object). §4: pseudo-transient reference removed (COUP-3 §3.6 tombstone). §5 (S21): anchor budget rewritten — full-3-D **physical march** under a declared compressed start window, GPU-ported, ≤ 24 h; multi-fidelity ensembles; plan pointer. §6.7 updated to match. |
| 2026-08-14 | 0.3 | **Post-review fix wave (S16, S17, S18, S19, S20, S21; rulings D-A, D-C, D-D).** §3.5 rewritten to the **one local wall-function heat-flux law** for all engines (local near-wall operands, ±20–30% band; Bartz demoted to VAL-2 nozzle-envelope oracle) with the S16 ownership delineation (SOLV-1 evaluates at config-time-identified wall faces, `F_visc` suppressed there; COUP-2 owns sweep placement; COUP-7 owns coolant side only). §3.4: η_c\* prior applied as **in-solver source-term combustion-completeness knockdown** (output-side multiplication forbidden, S18); **shifting mode advects elemental fractions with per-step equilibrium projection / frozen mode advects full species** (frozen↔shifting = discrete epistemic dimension, S19; §3.1 updated); resolved tier gains the **one universal LES-class dynamic-coefficient subgrid mixing closure**, offline-calibrated on canonical turbulence data, per-quantity band (D-D/S20 — laminar-resolved mixing at Re~10⁷ under-mixes by orders of magnitude). §5/§6: **anchor-run budget** stated (grid class, pseudo-step count, cell-updates/s target, ensemble size) with steady-state acceleration = COUP-3's pseudo-transient mode (S21). D-A wording: "m=0 corner" → "N_θ=1 axisymmetric-with-swirl corner"; N_θ>1 is a finite-volume capability lift, no spectral modes (§0, §3.3, §3.4). Reaction lookups re-keyed to local (p, h, Z) state (S22 cross-ref). |
| 2026-08-13 | 0.2 | **R2 applied — resolved-mixing rung architected.** §3.4 mixing entry re-specified as the config-selected COUP-7 tier: prior tier (premixed + envelope-bounded η_c\* prior, milestone-1) vs resolved tier (unmixed streams on the existing `ρX_k` state; mixing via mode-ceiling/refinement lift; η_c\* emergent; per-quantity closure model-form; liquid injection rides the W4 two-phase extension). Added to the deferred split. OFFL-3 local-mixture-fraction table-envelope dependency noted. |
| 2026-07-21 | 0.1 | Initial draft (reacting-flow m=0 corner). Conserved system + PPM/HLLC-Batten Godunov FV; 2.5-D axisymmetric-with-swirl as the m=0 truncation of the adaptive solver; tabulated shifting/frozen equilibrium chemistry with the frozen↔shifting gap as a UQ band (finite-rate CFD out; injector c\*-efficiency boundary); Bartz wall-heat coupling as the expander drive; emergent chamber pressure; State-Redistribution cut-cell stability; oracles Castro/PeleC/Athena++. Two-phase/MHD extensions of the same operator deferred to W4. |
