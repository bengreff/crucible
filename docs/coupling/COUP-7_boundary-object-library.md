# COUP-7 — Boundary-Object Library & Registration Contract

| Field | Value |
|---|---|
| **ID** | COUP-7 |
| **Family** | COUP (Coupling & orchestration) |
| **Status** | Reviewed (2026-08-14) |
| **Depends on** | FND-1, COUP-8; FND-4 (config blocks), COUP-2 (ports), SOLV-1 (§3.4–3.5) |
| **Version** | 0.4 (2026-08-19: §3.2.2 compressed external schedules; igniter + valve/start-sequence objects — VISION_SCOPE v1.5) |

---

## 0. Purpose

A **boundary object** is a piece of the engine that is **fully decoupled from the reaction** and therefore
*not simulated as matter on the grid* — it enters the physics only as a **port** carrying stated conditions
(VISION_SCOPE §7.5; FND-2 §1.3). COUP-7 owns the **v1 boundary-object set** and the **registration
contract** every boundary object must satisfy: **{citation, validity envelope, error band}** — all three, or
it cannot be registered. This is the structural guarantee that the sandbox never smuggles an un-sourced
number into a result.

The governing design decision for the chemical slice (Ben, 2026-07-21): **the turbopump is a boundary object
that delivers a propellant mass-flow; chamber pressure is never imposed — it *emerges* from the reacting-flow
solution.** §3.2 makes this concrete and shows why it needs **no new coupler type**.

Read after FND-1 (`UncertainInput` — a boundary object is the archetypal producer of one) and COUP-8 (the
`ports[]`/`Manifest` contract a boundary object registers through).

## 1. Scope & razor ruling

**Owns:** the boundary-object abstraction and its **three-field registration contract** (§3.1); the **v1
chemical-slice set** (§3.3); how a boundary object exposes uncertainty (§3.4); the **emergent-quantity rule**
(a boundary object states inputs, never outputs the reaction should determine — §3.2).

**Defers (owner named):** the `ports[]` declaration + wiring validation → **COUP-8**; the config block a
boundary object is selected/parameterized by → **FND-4**; the every-step port-flux conservation audit →
**COUP-2**; the uncertainty *type* a band is expressed in → **FND-1**; the radiator/thermal port's
downstream heat handling → out of scope (heat crossing the radiator port exits the simulation, VISION_SCOPE
§4.2).

**Razor ruling (Reaction Razor, VISION_SCOPE §4.1):** a subsystem is a boundary object **iff no energy
release or thrust generation flows *through its internals***. The turbopump qualifies — its blade
aerodynamics, bearings and gear train are solved mechanical engineering — **even though the energy that
*drives* it (regen-jacket heat pickup) is fully simulated**. The line is drawn at the object's casing: what
crosses its ports is physics; what happens inside is a literature performance map. A subsystem whose
*internal* process carries the reaction (CNTR rotating fuel, ablating throat) is **not** a boundary object —
it is grid matter (FND-2 §1.3).

## 2. Interfaces & contracts

| Interface | Consumed by | Contract |
|---|---|---|
| **Boundary-object descriptor** (§3.1) | COUP-8 registry, FND-4 loader | id + `ports[]` (role/kind) + citation + validity envelope + error band + config params |
| **Port emission/consumption** (§3.2) | COUP-2 (audit), the operators on the grid | injects/extracts mass, momentum, power, or heat at a declared port with a stated condition; participates in the every-step port balance |
| **`UncertainInput` production** (§3.4) | COUP-5 (sampler) | every stated condition + map coefficient is an `UncertainInput` (nominal + family/interval + provenance + envelope) |

**Invariant promised to everyone:** a boundary object **states conditions at its ports and nothing else** —
it never writes a field the grid physics should determine (§3.2), and it exposes **no quantity without a
{citation, envelope, band}** (§3.1). An object missing any of the three is a registration error and the run
**halts at load** (COUP-8 §3.3).

## 3. Method

### 3.1 The registration contract — three fields, non-negotiable

Every boundary object carries, per stated quantity:
- **Citation** → a META-3 `key` (the literature source for the value/map). No key ⇒ cannot register.
- **Validity envelope** → the `ValidityEnvelope` (FND-1 §3.7) over which the value/map is declared valid
  (e.g. a pump map valid near its design speed/flow-coefficient). Interrogation outside it flags/refuses
  per config (FND-5 §3.5).
- **Error band** → an `UncertainInput` (FND-1 §3.4a): a distribution family or an epistemic **interval**,
  propagated by COUP-5 into the result's p-box.

Enforced **mechanically** at registration (COUP-8 §3.3, structural-resolve check), not by discipline — this
is the boundary-object analog of "provenance is not optional" (META-1 Principle 7).

### 3.2 The emergent-quantity rule *(the load-bearing decision)*

A boundary object **states its inputs to the reaction and lets the reaction determine everything downstream.**
It must **never inject a quantity the grid physics is responsible for producing.** Applied to the turbopump:

- The turbopump **emits a propellant mass-flow ṁ (and inlet enthalpy/composition) at the injector port** —
  the condition at which propellant *enters* the reaction.
- **Chamber pressure `p_c` is NOT a boundary condition.** It **emerges** from the reacting-flow solution:
  the injected ṁ combusts and chokes at the throat, and `p_c ≈ ṁ·c*(p_c,MR,T_c)/A_t` is whatever the field
  solver produces (SOLV-1/SOLV-7). Imposing `p_c` would feed in the very answer the RL10 anchor is meant to
  *predict* (circular validation) — forbidden.
- The turbopump's **delivered ṁ is set by its performance map + the drive power available to it.** For an
  **expander cycle** that drive power is the **regen-jacket enthalpy rise the solver computes** (`expander-cycle`,
  `rl10-cycle-data`): the jacket heat is real reaction energy through the wall→coolant path — **computed**,
  on the gas side, through the one local wall-function closure (SOLV-1 §3.5) with its declared **±20–30%
  band** (the honest dominant band on the expander drive, not a solved-exactly path), plus the coolant-side
  channel correlation this object carries (§3.3); only the pump/turbine *machinery* is the abstracted map.

**Drive-power input as a port — one object, two modes, no seam.** The turbopump descriptor exposes a
`drive_power` **Require** port. It may be satisfied two ways, **selected in config, not by a code branch**:
1. **Open (spec-driven):** `drive_power` = a declared literature constant → the pump delivers its rated ṁ.
   Used to validate the **reacting flow in isolation** (given the right ṁ and MR, does `p_c`, thrust, Isp
   emerge at the RL10 values?).
2. **Closed (cycle-coupled):** `drive_power` bound to the solver's jacket enthalpy rise (supplied by the
   **wall-exchange coupler**, COUP-2/VISION_SCOPE §7.4 coupler 2) → the full expander loop
   `wall heat → turbine power → pump map → ṁ → chamber conditions → wall heat` closes, and the pump
   operating point is solved to consistency with the downstream impedance within the step.

Both modes are the **same boundary object** reading the same port; the power balance lives **inside the
object** (its turbine+pump map), fed by existing couplers (wall-exchange in, port-accounting out). **No new
coupler type is introduced** (the §7.4 list is unchanged) — the earlier "cycle-power coupler" amendment
considered under Fork 2 is **not taken**: Ben's ruling keeps the turbopump squarely outside the reaction,
so its balance is internal to the boundary object, not a first-class coupler. The closed-mode within-step
consistency solve uses the **fixed-sweep-count** relaxation of COUP-3 (determinism, META-1 §2).

*Consequence:* RL10 validation proceeds **open-mode first** (isolates reacting-flow error) then **closed-mode**
(adds the cycle-balance error) — the blind/calibrated pair of VAL-2's overlap-band criterion.

**The rule's sharp form (R2+R3, 2026-08-13):** a boundary object must never state an **efficiency or
energy-partition** the grid could compute by transporting the actual flow/products through the actual
geometry. Efficiency is an *output*; residual uncertainty is model-form on *universal* closures, carried
per-quantity, never as a blanket band. A supplied prior is admissible only where it is **measurement, not
assumption** — i.e. strictly inside the validity envelope of the hardware family it was measured on — and
the envelope refusal (§3.1, FND-5 §3.5) enforces that mechanically: a novel design *cannot* silently ride a
measured family's number.

### 3.2.1 The tiered injector (R2) — prior tier vs resolved tier

Like the turbopump's open/closed `drive_power` modes, the injector is **one boundary object with two
config-selected tiers** — data, never a code branch (Rule 13):

1. **Prior tier** *(W2 / milestone-1)*: the object states the injection-plane inflow **premixed at the
   declared MR**, plus a cited **η_c\* prior** (`injector-cstar-eff`). **Application operator (O19):** the
   prior is applied as SOLV-1 §3.4's **in-solver combustion-completeness knockdown at the source-term level**
   — the equilibrium heat release is scaled so delivered c\* = η_c\*·c\*_ideal at the anchor state, keeping
   T_c, wall heat, jacket ΔH, and emergent `p_c` self-consistent with the deficit; **output-side
   multiplication of c\*/Isp is forbidden** (SOLV-1 §3.4 carries the mechanism; this is its boundary-object
   mirror). Legitimate *only inside* the cited injector-family/MR/Pc envelope. **Blind vs calibrated (D-G,
   VISION_SCOPE §9 v1.4):** a **blind** run consumes **no value measured on the engine under test** — its
   prior is the **injector-class band measured on other hardware** (coaxial-element family, ±1–3%); the
   RL10-fitted **η_c\* = 0.9892 (`rl10-cycle-data`) binds only in CALIBRATED (closed-mode) runs**, labeled as
   such. A config whose injector lies outside every cited family envelope **fails envelope interrogation**
   and must use the resolved tier: the blind-prediction contract is enforced by the registration machinery,
   not by discipline.
2. **Resolved tier** *(architected now; built as a later lift — rides SOLV-1's W4 two-phase extension for
   liquid injection, needs only the reacting-gas subset at sufficient resolution for gas–gas injection)*:
   the object states the **actual unmixed injection-plane conditions** per element class — stream
   composition, per-stream ṁ and inlet enthalpy, injection momentum/angle, and (liquid streams) droplet
   size/velocity distributions from **universal atomization closures** with their own {citation, envelope,
   band}. Mixing then happens **on the grid** (SOLV-1 §3.4); η_c\* becomes an **emergent output** (read by
   SOLV-7 as c\*_delivered/c\*_ideal) and the residual uncertainty is the per-quantity model-form of the
   atomization/mixing closures. The object states *conditions*, the reaction determines *performance* —
   the emergent-quantity rule with no exception.

**Tier-consistency check** (validation obligation, deferred with the build): the resolved tier, run on an
anchored injector (RL10 coax), must reproduce that family's measured η_c\* within the declared closure
bands — the same blind→calibrated laddering as open/closed turbopump modes (§6 item 7).

### 3.2.2 Compressed external schedules *(VISION_SCOPE §7.6 v1.5, Ben 2026-08-19)*

A boundary object whose stated condition is a **timeline** (valve travel, pump-rotor spin-up, tank-head
sequence, igniter firing window) may declare that timeline **time-compressed** relative to the hardware it
cites — the cost lever that keeps physical-march start-ups affordable now that accelerated convergence is
deleted. Rules: (a) the compression is **config data, recorded in the manifest** (schedule + the cited
physical timeline it compresses); (b) only *externally imposed* timelines qualify — **nothing that feeds
back through the reaction may be compressed** (the expander's jacket pickup, the liner's thermal mass, any
flame/acoustic/thermal timescale the grid evolves — those are physics); (c) a certificate states the
compression next to its labels. The declared-schedule class (valve/tank-head/igniter objects, §3.3) all
carry it.

### 3.3 The v1 boundary-object set (chemical slice)

| Object | Ports (role) | Stated quantity + map | Citation | Envelope | Band |
|---|---|---|---|---|---|
| **Turbopump / feed** | injector-flow (Provide), `drive_power` (Require) | delivered ṁ, inlet enthalpy, from pump/turbine η-maps + drive power (§3.2) | `rl10-cycle-data`, `huzel-huang` | near design speed/flow-coeff | η ±3–5 pts → propagates to `p_c` |
| **Injector** *(tiered, §3.2.1)* | combustion-inlet (Provide) | **Prior tier:** premixed inflow @ MR + cited η_c\* prior (a c\* deficit). **Resolved tier:** actual unmixed stream conditions (composition, per-stream ṁ, momentum, droplet distributions); η_c\* an *output* | `injector-cstar-eff`, `rl10-cycle-data` | **prior: cited injector family/MR/Pc only** (outside ⇒ resolved tier, mechanically); resolved: closure envelopes | prior: ±1–3% mature families; resolved: per-quantity closure model-form |
| **Cooling-jacket / regen interface** | wall-heat (Bidirectional), coolant-return (Provide) | **coolant side ONLY (D-C):** channel correlation + coolant state → jacket ΔH uptake. Gas-side h is **SOLV-1 §3.5's one wall-function law** (`wall-function-heat`; Bartz demoted to a VAL-2 nozzle oracle) — never stated here | `huzel-huang` (coolant channels) | coolant channel design range | coolant-side correlation band (declared per source); the gas-side ±20–30% band is the wall function's (SOLV-1 §3.5, epistemic) |
| **Radiator interface** | heat-sink (Require) | capacity + coolant-return T; heat crossing it **exits the sim** (VISION_SCOPE §4.2) | *(per-config)* | stated capacity | stated |
| **Pressurization / tank feed** | inlet (Provide) | inlet pressure/temperature/composition | *(per-config)* | tank state range | stated |
| **Igniter** *(v0.4)* | ignition-energy (Provide) | a **scheduled, localized energy deposit**: position, radius, duration, energy, firing window (§3.2.2 schedule class). Not special-cased physics — it creates a state whose SOLV-4 §3.6 induction term fires; ignition/no-light is the field's outcome, never this object's claim | `spark-igniter-class` | deposit within grid + table envelopes; energy in the cited igniter class | stated per source |
| **Valve / start-sequence schedules** *(v0.4)* | commanded-profile (Provide) | declared timelines for valve open fractions, tank-head/pump spin-up ṁ availability — the compressed-schedule class (§3.2.2): compression declared + manifest-recorded | `rl10-cycle-data` (start sequence), *(per-config)* | within the cited hardware's sequence class | stated |

*(Nuclear/plasma/pulsed boundary objects — pulsed-power drivers, beam/laser sources, coil sets, reactivity-
control schedules, antiproton delivery, EP performance tables, fusion-confinement source — are added to this
table in later waves under the same three-field contract.)*

### 3.4 How a boundary object exposes uncertainty

Every stated condition and every map coefficient is an `UncertainInput` (FND-1 §3.4a): a nominal plus either
a distribution family or an epistemic **interval**, with its `ProvenanceRef` and `ValidityEnvelope`.
**Which COUP-5 loop it lands in follows the FND-1 family (O15):** interval-typed inputs enter the
**epistemic outer loop** (the box of COUP-5 §3.3); distribution-family inputs enter the **aleatory inner
ensemble** — the assignment is the type, never a per-object choice. Sobol (COUP-5) attributes how much of a
result's spread each boundary object drives — for the expander cycle, the **wall-function ±20–30% band
(SOLV-1 §3.5, an epistemic interval — not this object's)** is expected to dominate the `p_c` spread, and the
report will say so.

### 3.5 Determinism

A boundary object holds **no state that depends on schedule** — its maps are pure functions of port inputs +
config coefficients, evaluated in fixed order. The closed-mode consistency solve (§3.2) uses fixed-sweep-count
relaxation (COUP-3), so the operating point is a pure function of `{config, computed jacket ΔH, seed}` — bit-
reproducible at any thread count (META-1 §2).

## 4. Coupling relationships
- **COUP-8** registers the descriptor and validates its `ports[]` wiring (every `Require` has a `Provide`);
  the `drive_power` Require of a closed-mode turbopump is satisfied by the wall-exchange coupler's Provide.
- **COUP-2** audits port fluxes every step (mass/energy/momentum in = stored + out); a boundary object is a
  port term in that balance, never an interior source that escapes it.
- **FND-1** supplies the `UncertainInput`/`ValidityEnvelope` types; **COUP-5** samples them; **COUP-6** scores
  a result's pedigree partly from boundary-object band width and citation independence.
- **SOLV-1/SOLV-7** consume the injected ṁ/inlet condition and produce the emergent `p_c`, thrust, Isp, c\*.
- **FND-4** carries each object's `[couplers]`/boundary block; **VISION_SCOPE §7.4** coupler list is
  unchanged (§3.2).

## 5. Uncertainty & validity
COUP-7 originates **no physics model-form uncertainty of its own** — it *carries* each boundary object's
declared band into the ensemble. Its correctness obligation is that **no un-banded, un-cited, or
un-enveloped quantity can enter a result** (registration contract, §3.1), CI-gated (VAL-3). The physical
risk it governs is the boundary-object bands themselves, which are first-class UQ inputs (the dominant
gas-side wall-heat band is SOLV-1 §3.5's, not this doc's).
Validation-ladder status: infrastructure (rung i for its tests); the *values* it carries are validated where
their source anchor is (e.g. the RL10 pump/jacket data, VAL-2).

## 6. Validation plan
1. **Registration completeness:** an object missing citation, envelope, or band fails registration (COUP-8) —
   no run starts.
2. **Emergent-quantity guard:** a config that tries to have a boundary object *impose* `p_c` (a reaction-
   determined quantity) is rejected — boundary objects may only write declared input ports.
3. **Open vs closed turbopump:** open-mode delivers rated ṁ; closed-mode's operating point converges (fixed
   sweeps) and, fed the RL10 jacket heat, reproduces the RL10 pump discharge/ṁ within the declared band.
4. **Port conservation:** injected mass/energy at all ports closes the COUP-2 audit to tolerance.
5. **Uncertainty propagation:** the wall-function ±20–30% band (SOLV-1 §3.5) produces a `p_c` spread via the
   jacket drive, in the epistemic outer loop; Sobol attributes it.
6. **Determinism:** identical operating point at 1 vs N threads and across runs.
7. **Tier consistency** *(deferred with the resolved-tier build)*: the resolved injector tier on the RL10
   coax element reproduces the measured η_c\* = 0.9892 within the declared atomization/mixing closure bands;
   a config outside every cited family envelope is refused at the prior tier (blind-prediction guard).

## 7. References
META-3 keys: `expander-cycle`, `rl10-cycle-data`, `huzel-huang`, `injector-cstar-eff` (also the
other-hardware family band for blind mode), `wall-function-heat` (SOLV-1's law, referenced), `bartz` (VAL-2
oracle only),
`modelica-connector` (port balance), `spark-igniter-class` (igniter object, §3.3). Depends on FND-1 (`UncertainInput`, `ValidityEnvelope`), COUP-8
(`ports[]`, registration), COUP-2 (port audit), SOLV-1 §3.4–3.5 (knockdown mechanism, wall function),
FND-4 (config block).

*(No open questions — the turbopump/emergent-`p_c` fork resolved by Ben 2026-07-21; the injector η_c\* fork
resolved by Ben 2026-08-13 (R2): tiered integration, prior tier envelope-bounded, resolved tier architected
now / built later.)*

## 8. Change log
| Date | Version | Change |
|---|---|---|
| 2026-08-19 | 0.4 | **VISION_SCOPE v1.5 (Ben).** New §3.2.2 **compressed external schedules** (declared, manifest-recorded, external-timelines-only — the physical-march cost lever; certificates state compression). §3.3 gains the **igniter** (scheduled localized energy deposit — outcome belongs to the field, SOLV-4 §3.6) and **valve/start-sequence schedule** objects. |
| 2026-08-14 | 0.3 | **Review fix wave (O19, D-C, D-G).** §3.2.1 prior tier: the η_c\* application operator stated — SOLV-1 §3.4's **in-solver combustion-completeness knockdown** (delivered c\* = η_c\*·c\*_ideal at the anchor state; output-side multiplication forbidden; this doc is the boundary-object mirror) (O19); **blind vs calibrated split** — blind runs use the injector-class band measured on *other* hardware (±1–3%); the RL10-fitted η_c\* = 0.9892 binds in calibrated mode only (D-G, VISION_SCOPE §9 v1.4). §3.3 cooling-jacket row + §3.2 wording: the jacket object owns **only the coolant side**; gas-side h is SOLV-1 §3.5's one wall-function law (Bartz → VAL-2 nozzle oracle); "simulated in full" softened to name the wall-function ±20–30% closure band honestly (D-C). §3.4: epistemic-vs-aleatory assignment follows the FND-1 family (interval ⇒ outer loop) — COUP-5 §3.3's O15 disambiguation mirrored; §5/§6.5 Bartz references re-pointed to the wall function. |
| 2026-08-13 | 0.2 | **R2 applied — tiered injector.** Injector re-specified as one object with two config-selected tiers (§3.2.1): envelope-bounded η_c\* **prior** tier (milestone-1) vs **resolved** tier (actual unmixed stream conditions; mixing computed on the grid; η_c\* an output; per-quantity closure model-form). Emergent-quantity rule sharpened to its R2+R3 form (never supply a computable efficiency/partition; priors only inside their measurement envelope, enforced mechanically). v1 set table updated; tier-consistency validation item added (deferred with build). |
| 2026-07-21 | 0.1 | Initial draft. Three-field registration contract (citation+envelope+band, enforced at registration); the **emergent-quantity rule** and the turbopump-as-flow-source model with a config-selected open/closed `drive_power` port (chamber pressure emergent, no imposed `p_c`, no new coupler type — Ben 2026-07-21); the v1 chemical-slice boundary set (turbopump, injector c\*-efficiency, cooling-jacket/Bartz, radiator, pressurization) with per-object bands; uncertainty exposure via `UncertainInput`; determinism via fixed-sweep closed-mode solve. |
