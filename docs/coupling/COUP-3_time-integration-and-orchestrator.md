# COUP-3 — Time Integration & Orchestrator

| Field | Value |
|---|---|
| **ID** | COUP-3 |
| **Family** | COUP (Coupling & orchestration) |
| **Status** | Reviewed (2026-08-14) |
| **Depends on** | FND-2, SOLV-1, COUP-2; COUP-8, COUP-4 |
| **Version** | 0.4 (2026-08-19: §3.6 pseudo-transient DELETED — VISION_SCOPE v1.5) |

---

## 0. Purpose

COUP-3 owns the **global clock and the per-step operator schedule**: how the explicit hyperbolic update
(SOLV-1), the **spatially-coupled implicit diffusion** (conduction/viscous), and the **genuinely cell-local
implicit stiff sources** (reactions, M1 source coupling) are advanced **together** without the order-loss and
conservation error that operator splitting incurs for energetic stiff reactions. The method is **SDC-coupled
IMEX** (spectral deferred corrections): explicit hydro/MHD + implicit stiff solves, iterated a fixed number of
sweeps to 2nd order. It also owns the Δt rule (single global Δt, v1), radiation (RSLA) sub-stepping, the
**closed-mode expander consistency solve** (§3.5), and pulsed-event sequencing. *(The former §3.6
pseudo-transient steady-state mode is **DELETED** — VISION_SCOPE v1.5: physical march only.)*

Read after SOLV-1 (the explicit update it drives), COUP-2 (the audit it triggers each step), and META-1 §2
(the determinism mandate it must satisfy).

## 1. Scope & razor ruling
**Owns:** the time-integration scheme (SDC-IMEX), the three operator classes and their solvers (§3.1), the
operator-split *schedule* (call order of the SOLV operators inside a step), Δt/CFL control incl. the RSLA
radiation sub-step rule (§3.4), the closed-mode expander consistency solve (§3.5), pulsed-event sequencing
(SOLV-5 hook), and the per-step→audit handoff. **Defers:** which fields each operator reads/writes and the conservation audit →
**COUP-2**; the operators themselves → **SOLV-1…8**; the two-stage (function/lifetime) control flow and halt
handling → **COUP-4**; the mechanism registry → **COUP-8**.

**Razor ruling:** pure infrastructure — no physics. Its fidelity obligation is that the *coupling* introduce
no spurious order-loss or conservation error (a splitting seam would be a Rule-12 violation in time).

## 2. Interfaces & contracts
| Interface | Consumed by | Contract |
|---|---|---|
| **Step advance** | the orchestrator/COUP-4 | advance `U` by Δt with a stated order; every fixed-sweep coupling solve passes its named residual-acceptance check (§3.5) or raises `COUPLING_RESIDUAL` |
| **Operator schedule** | SOLV operators | a deterministic call order per step; explicit fluxes then implicit stiff solves, iterated |
| **CFL / Δt** | SOLV-1 | a deterministic timestep from a fixed CFL rule over the wave speeds |

**Invariant:** the advance is **deterministic** (fixed sweep count, fixed sub-cycle structure, fixed reduction
order) and **2nd-order** where the SDC sweeps converge; it never uses a wall-clock or variable-iteration stop
(META-1 §2.2).

## 3. Method

### 3.1 SDC-coupled IMEX *(the core)* — four operator classes
Write the system as `U_t = A(U) + D(U) + R(U)`; each class is advanced by the solver its coupling structure
demands (O1):
- **`A` — explicit hyperbolic** (advection/MHD fluxes, SOLV-1): explicit at the fluid CFL.
- **`D` — spatially-coupled implicit diffusion** (conduction `∇·(k∇T)`, viscous fluxes): neighbor-coupled, so
  a *cell-local* implicit solve **does not exist** for it. Solved per SDC node by a **deterministic
  fixed-cycle solver**: geometric multigrid (fixed V-cycle count, fixed smoother type and sweep order) or
  preconditioned CG with a **fixed iteration structure** (fixed count, fixed-order reductions). Implicit
  treatment removes the diffusion CFL `Δt ∝ Δx²`, which binds exactly where conduction is stiff (fine wall
  cells — the RL10 closed-mode path); a Jacobi/cell-local fallback loses unconditional stability there and is
  forbidden. [META-3: `geometric-multigrid`]
- **`R` — genuinely cell-local stiff sources** (reactions, M1 radiation–matter *source* coupling): the
  cell-local implicit ODE solve of §3.3. "Cell-local implicit" is **reserved** for this class.
- **`G` — declared global stiff ODE systems** (N6; the pulsed-event carve-out generalized): a system whose
  right-hand side is an **integral functional over grid fields** — the reaction-kinetics amplitude system
  (SOLV-4 §3.2: ρ, the adjoint-weighted delayed source c̃_i, s̃) is the first occupant. A class-`G` system is
  evaluated **once per SDC sweep in fixed order** (its functionals are fixed-order grid reductions, FND-2
  §3.7) against the sweep-current fields, its outputs (e.g. power) feeding the deposition sources within the
  same sweep. It is neither cell-local (`R`) nor neighbor-diffusive (`D`); a private integrator inside an
  operator remains forbidden — class-`G` membership is declared in the COUP-8 Manifest.

Integrate with **spectral deferred corrections**: the stiff solves carry the explicit advection source
evaluated at the time node, and the coupled system is iterated **2–3 SDC sweeps** to lift to 2nd (optionally
4th) order. This is **energy-conserving and far below Lie-Trotter error at CFL** for stiff/energetic reactions.
**Stability/order consequence:** the step is unconditionally stable in `D` and `R` (Δt is set by the
hyperbolic CFL alone), and formal 2nd order is retained provided each implicit sub-solve's terminal residual
sits below the SDC truncation — guaranteed by sizing the fixed cycle/iteration counts in the
temporal-convergence test (§6) and guarded at runtime by the residual-acceptance checks (§3.5).
[META-3: `sdc-imex`, `castro-source`]

### 3.2 Why not Strang/Lie-Trotter for stiff reactions
Operator splitting decouples advection and reaction each substep, so a **stiff, energetic** reaction never
sees the advective change over the step → it **loses formal order and mis-partitions energy** unless Δt is cut
hard, and Strang has **no clean extension above 2nd order**. SDC-IMEX avoids this and evaluates the RHS
**fewer** times than splitting for vigorous burning. The chemical slice's combustion is mild, but the scheme
is chosen once for the whole chemical→antimatter roadmap (the nuclear/antimatter end is where splitting fails
outright). [META-3: `stiff-reactions`]

### 3.3 Stiff inner integrator *(class `R` only)*
The cell-local implicit solve uses a **variable-order BDF (VODE-class)** integrator with an analytic Jacobian
(numeric retry), a fixed max-step cap, and a tabulated-equilibrium hand-off at the high-`T` end — a clean,
**deterministic** Rust port (fixed iteration structure for Tier-1 determinism). For the chemical slice the
`R` set is tabulated heat-release; conduction/viscous terms belong to class `D` (§3.1) and are solved by the
fixed-cycle diffusion solver, never by this integrator. Full reaction networks arrive with the nuclear leg.

### 3.4 Timestep, radiation sub-stepping (RSLA), pulsed events
**One global Δt (v1) — region sub-cycling is cut (O5).** Δt from a fixed CFL over the hyperbolic wave speeds,
**global over the grid**; stiff classes `D`/`R` need no CFL (implicit). **Class-`G` Δt-limiter interface
(N6):** a class-`G` system may declare per-step accuracy limiters (named constants in its own doc — e.g.
SOLV-4's `EPS_RHO_STEP`, `EPS_POWER_STEP`); an exceeded limiter halves Δt for the next step through this
interface — part of the declared Δt rule, deterministic, never an ad-hoc clamp. Per-region sub-cycling is **out of
v1**: the refluxing bookkeeping at region interfaces, the audit interaction (COUP-2's telescoping assumes one
flux ledger per step), and the multi-rate-SDC coupling it would require are not worth their complexity for
milestone-1 physics. The door stays open — a later doc amendment must specify refluxing + audit interaction +
multi-rate SDC *together*, never one alone. [META-3: `berger-amr`]

**Radiation band — RSLA + sub-stepping (E-3).** When M1 radiation is active, `(E_rad, 𝐅_rad)` propagate at a
**reduced speed of light ĉ**, **declared per regime in config**, with the induced model-form error a
**PIRT-recorded band per regime** (VAL-1) — matching the Quokka oracle lineage. **Per-band Δt rule:** the
fluid Δt stands; within each fluid step the radiation band advances by `n_rad = ⌈CFL_rad·ĉ/max|λ_hyp|⌉`
**fixed sub-steps** of `Δt/n_rad` (a deterministic integer from the same fixed-order wave-speed reduction),
so the radiation CFL at ĉ is honored without collapsing the global Δt by the 10³–10⁴× light-speed ratio.
ĉ-sufficiency is validated by the standard RSLA insensitivity check (§6). [META-3: `rsla`, `radiation-m1`]

**Pulsed events** (SOLV-5, later waves) are sequenced here as bounded sub-simulations returning
impulse/energy-partition to `U`. A pulsed event is spatially extended and violently scale-separated (ns–µs
internal vs ms global), so it is **not** a cell-local scalar solve: it is a **non-cell-local, internally
sub-cycled stiff sub-integration** whose *net* deposition over the SDC node interval is the stiff-source
contribution the sweeps iterate (a multi-rate SDC term, admitted via the event-injection coupler,
VISION_SCOPE §7.4 #6) — kept inside SDC, never a Lie-Trotter split seam (§3.2). An event's *internal*
sub-cycling is a bounded sub-simulation, not the region sub-cycling cut above.

### 3.5 Closed-mode expander consistency solve *(owned here — COUP-7 §3.2 defers to this section)*
The closed expander loop (`jacket ΔH → turbine power → pump map → ṁ → chamber conditions → jacket ΔH`) is
solved **once per step, after the wall-exchange (Robin-Robin) solve of COUP-2 §3.5**, as a **fixed-point
iteration on the delivered ṁ** with **Aitken Δ² acceleration** and a **fixed sweep count `N_EXPANDER_SWEEPS`**
(named constant, default 4). [META-3: `aitken-iqnils`]
- **Ordering & SDC interaction:** the loop runs at **frozen field state** — the step's post-sweep `U` and its
  converged wall-heat integral; each sweep re-evaluates only the boundary-object algebra (wall-function heat
  pickup as a function of ṁ, turbine power, pump map → ṁ′). The accepted ṁ sets the injector-port inflow for
  the **next** step's explicit stage; the loop is never nested inside the SDC sweeps. The one-step lag is
  first-order in Δt, below the SDC truncation near the operating point, and identically zero at a converged
  steady state (§3.6).
- **Residual acceptance (O4):** after the fixed sweeps, require `|ṁ⁽ᵏ⁾−ṁ⁽ᵏ⁻¹⁾|/ṁ⁽ᵏ⁾ ≤ EPS_EXPANDER_RESID`
  (named constant, default 1e-8 relative — orders below the pump-map band, so the residual floor never
  contributes to the physics error). Failure ⇒ the **`COUPLING_RESIDUAL` halt** (COUP-4 §3.2; diagnosis =
  **numerical** — a solver defect, never an engine verdict).
- **Won't-bootstrap is a distinct, physical test:** monotone divergence of the iterates toward zero/maximum
  flow, or a fixed point reachable only by interrogating the pump map **outside its validity envelope** (the
  standard envelope refusal fires *inside* the loop): diagnosis = **physical**, verdict `DOESN'T WORK (cycle
  won't bootstrap)` (COUP-4 §3.2). A small-but-unaccepted residual is a numerical defect; a divergent or
  envelope-refused iteration is an engine with no operating point — never conflated.
- **Expected contraction:** the loop gain `|∂ṁ′/∂ṁ| < 1` in the expander regime (jacket pickup sub-linear in
  ṁ — h ∼ ṁ^0.8-class scaling — and delivered ṁ sub-linear in drive power through the pump map), and Aitken's
  adaptive relaxation is clamped to a deterministic `[ω_min, ω_max]` so the relaxed gain stays a contraction;
  expected factor ~0.5–0.8/sweep at the RL10 point, which sizes `N_EXPANDER_SWEEPS` to reach
  `EPS_EXPANDER_RESID` with margin. Clamp bounds and constants are manifest-recorded.

### 3.6 ~~Pseudo-transient mode~~ — DELETED (VISION_SCOPE v1.5, Ben 2026-08-19)
**Tombstone.** The pseudo-transient/local-Δt steady-state continuation mode (and dual-time stepping, its
per-step form) is **deleted, not deferred**: local Δt is not physical evolution, and reaching the operating
point — start-up included — is itself part of what Stage 1 verifies (S5). **Every Stage-1 run is a physical
transient march at the global Δt of §3.4.** The declared cost levers are: **compressed external schedules**
(boundary-object timelines only — COUP-7 §3.2.2), adaptive resolution (FND-2 §3.4/§3.5), and GPU throughput
(`PLAN_CHEMICAL_SANDBOX.md` §3). Grid-sequenced restart (a finer grid initialized from a coarser grid's
settled *physical* state, FND-6) remains legal — every step of every march is still physical. The section
number is retained to keep §3.7 stable; the `pseudo-transient` META-3 key is retired from §7. *(Forward
note: stepping far over the acoustic CFL during genuinely slow phases — e.g. an NTR heat-up over minutes —
would require an implicit/all-speed treatment of the hyperbolic class, a possible amendment for the nuclear
wave; it is not licensed by this section.)*

### 3.7 Determinism
Fixed sweep count, fixed diffusion cycle count (§3.1), fixed sub-step structure (§3.4), fixed CFL rule, fixed
Aitken relaxation clamp (§3.5), fixed-order reductions (FND-2 §3.7). Convergence
criteria are **absolute + deterministic** (fixed tolerance, fixed max-iters) — never "converged by wall-clock
budget" (META-1 §2.2). The schedule is a `match` over the COUP-8 operator enum in fixed source order (no
`HashMap`, no distributed-registration order — COUP-8 §3.2).

## 4. Coupling relationships
- **SOLV-1** provides `A(U)` (explicit fluxes) and the stiff terms; **COUP-2** audits after each advance and
  owns which operator touches which field; **COUP-4** wraps the step in the function/lifetime march and
  handles halts; **COUP-8** fixes the deterministic operator order.
- **SOLV-2** radiation-matter exchange and **SOLV-4** reaction sources enter as implicit `R` terms; **SOLV-8**
  (stage 2) runs with large deterministic timesteps (slow clocks), a different Δt regime COUP-3 also owns.

## 5. Uncertainty & validity
Introduces **temporal discretization error** (budgeted, feeds UQ) and a **coupling-residual** floor (bounded
by the fixed sweep/cycle counts and guarded by the §3.5 acceptance constants). The RSLA ĉ is the one
model-form declaration it carries — per-regime, PIRT-recorded (§3.4). Validity: the SDC order is demonstrated
by a temporal convergence test; a stiff case outside the integrator's convergence flags/halts, never silently
under-resolves. CI-gated (VAL-3).

## 6. Validation plan
1. **Temporal order:** SDC-IMEX recovers 2nd order on a manufactured advection-reaction problem. [META-3: `mms`]
2. **Strang failure fixture:** reproduce the order-loss of Strang splitting on a stiff reacting problem as a
   regression guard (and confirm SDC does not exhibit it). [META-3: `stiff-reactions`]
3. **Conservation:** a closed advance conserves mass/energy to the COUP-2 tolerance across the coupling.
4. **Differential oracle:** matches Castro's SDC path on a shared reacting case. [META-3: `castro-source`]
5. **Determinism:** identical trajectory at 1 vs N threads and across runs; fixed sweep count verified.
6. **Diffusion class:** manufactured conduction/viscous solution — the fixed-cycle multigrid/CG solve reaches
   design order and remains stable at wall-cell stiffness where a Jacobi/cell-local treatment fails (O1).
7. **Expander fixed point:** contraction on an RL10-class fixture; a planted non-contractive case trips
   `COUPLING_RESIDUAL` (numerical) while a no-fixed-point/envelope-refused case trips won't-bootstrap
   (physical) — the two diagnoses never conflate (§3.5).
8. **RSLA:** result insensitive to ĉ above the declared per-regime value; `n_rad` deterministic (§3.4).
9. ~~Pseudo-transient~~ — deleted with §3.6 (VISION_SCOPE v1.5).

## 7. References
META-3 keys: `sdc-imex`, `castro-source`, `stiff-reactions`, `mms`, `gamer2-determinism`,
`geometric-multigrid`, `rsla`, `radiation-m1`, `aitken-iqnils`, `berger-amr` (`pseudo-transient` retired
2026-08-19). Depends on
FND-2 (grid/traversal), SOLV-1 (explicit update), COUP-2 (audit/coupling/wall-exchange), COUP-7 (the expander
boundary object whose consistency §3.5 solves), COUP-8 (operator order), COUP-4 (halts).

*(No open questions — SDC-IMEX over Strang is settled by the stiff-reaction literature; determinism policy
(fixed sweep count) per META-1 §2.)*

## 8. Change log
| Date | Version | Change |
|---|---|---|
| 2026-08-19 | 0.4 | **VISION_SCOPE v1.5 (Ben): §3.6 pseudo-transient mode DELETED** (tombstoned in place; §3.7 numbering kept). Stage 1 = physical march only; cost levers = compressed external schedules (COUP-7 §3.2.2) + adaptive resolution + GPU; grid-sequenced restart from settled physical states stays legal (FND-6). `pseudo-transient` key retired; §3.7/§6.9 references cleaned; forward note on implicit/all-speed acoustics for slow-phase (nuclear-wave) profiles. |
| 2026-08-14 | 0.3 | **Verification pass (V2/N6).** Class-`G` declared global-stiff-ODE slot added to §3.1 (integral-functional systems — reaction kinetics first occupant; evaluated once per SDC sweep in fixed order; Manifest-declared) + the class-`G` Δt-limiter interface in §3.4 — the slot SOLV-4 §3.2 binds to. |
| 2026-08-14 | 0.2 | **Review fix wave (O1, O3, O4, O5, O6, E-3).** §3.1: third operator class — **spatially-coupled implicit diffusion** via deterministic fixed-cycle multigrid/CG; "cell-local implicit" reserved for genuinely local stiff sources (reactions, M1 source coupling); stability/order consequence stated. §3.4: region **sub-cycling cut from v1** (single global Δt; rationale + amendment door); **RSLA** radiation with per-regime declared ĉ, PIRT-recorded band, fixed integer sub-stepping (per-band Δt rule). New §3.5: COUP-3 **owns the closed-mode expander solve** — fixed-point on ṁ, Aitken, `N_EXPANDER_SWEEPS`, once per step after wall-exchange, frozen-field ordering + SDC interaction; `EPS_EXPANDER_RESID` acceptance → `COUPLING_RESIDUAL` (numerical) vs won't-bootstrap divergence/envelope-refusal (physical); contraction expectation + clamped Aitken. New §3.6: **pseudo-transient mode** defined (local-Δt/SER, deterministic; audit applies to the converged state; cost vs §8 budget). Determinism renumbered §3.7. |
| 2026-07-21 | 0.1 | Initial draft. SDC-coupled IMEX (explicit hyperbolic + cell-local implicit stiff sources, 2–3 sweeps to 2nd order); rationale vs Strang/Lie-Trotter for stiff energetic reactions; deterministic BDF/VODE-class inner integrator; CFL/sub-cycling/pulsed-event sequencing; fixed-sweep determinism and fixed-order operator schedule; Castro SDC differential oracle. |
