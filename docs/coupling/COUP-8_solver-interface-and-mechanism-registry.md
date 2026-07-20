# COUP-8 — Solver Interface & Operator Registry

| Field | Value |
|---|---|
| **ID** | COUP-8 |
| **Family** | COUP (Coupling & orchestration) |
| **Status** | Draft — **skeleton** (per META-0 §5, W1) |
| **Depends on** | FND-1, FND-2, FND-5 |
| **Version** | 0.2 (v1.3 unified-grid pivot) |
| **Skeleton/complete split** | **Fixed now:** the `Manifest` capability-declaration contract; the closed-enum vocabularies (grid-field reads/writes on `U`, coupler kinds, exchanged quantities, halt conditions); the deterministic registry + dispatch decision; the four config-time wiring-validation checks; the interface-version policy; the *declared-or-forbidden* invariant. **Deferred** (filled as dependents land, W2+): the concrete Rust trait method signatures; each of the 8 SOLV operators' actual manifests (each SOLV doc supplies its own); the per-step operator-split call sequence (→ COUP-3). |

---

## 0. Purpose

COUP-8 is the **plugin contract every solver/mechanism implements** so that an engine is composed **from config as pure data, with no per-concept custom code** (VISION_SCOPE §4.3; META-1 Rule 13). It is the concrete enforcement machinery for the sandbox mandate: defining a new engine or reactor selects mechanisms from a **registry**, and the orchestrator wires them — a configuration that would need bespoke physics is a defect fixed in the *core*, never papered over in the config.

This doc fixes *what a mechanism must declare about itself* (§3.1), *how it is registered and dispatched deterministically* (§3.2), and *how a selected set of mechanisms is validated as fully wired before any solve* (§3.3). It is the contract all 8 SOLV operators bind to; getting its extension surface right is why it is written (as a skeleton) in the spine wave.

Read after META-1 (Rules 12/13), FND-2 (the grid = unified solver). It is consumed by **FND-4** (the loader that invokes §3.3) and by the orchestrator (COUP-3/COUP-4). *(v1.3: "mechanisms" are the unified-grid **operators** of SOLV-1…8; the former reduced-dimension mesh-binding is retired with COUP-1.)*

## 1. Scope & razor ruling

**Owns:** the mechanism **capability declaration** (`Manifest`); the registry data structure and its **deterministic** construction/iteration; the **static-vs-dynamic dispatch** ruling; the **config-time wiring-validation** checks (their *definitions*); the **interface-version** contract; the invariant that a mechanism may depend only on what it declares.

**Defers (owner named):**
- The **config document structure**, the loader pipeline, and the run manifest → **FND-4** (FND-4 *invokes* the §3.3 checks; COUP-8 defines them).
- The **adaptive azimuthal-mode reduction** (how 3-D collapses to 2-D-axi where symmetric) → **FND-2 §3.4** (COUP-8 operators are dimension-agnostic; they declare only which fields of `U` they read/write on the one grid).
- The **coupler catalog + conservation audit** (the 7 coupler types, the every-step global audit) → **COUP-2** (COUP-8 only records *which* coupler ports a mechanism participates in).
- The **boundary-object registration contract** (citation + envelope + band) → **COUP-7** (COUP-8 only records *which* ports a mechanism needs/provides).
- **Table schema/loader/interpolation** → **FND-5**; the **halt/verdict control flow** → **COUP-4**; **time integration / operator split** → **COUP-3**.

**Razor ruling:** pure infrastructure — no physics is simulated here. Its fidelity obligation is that the declaration be *complete enough* that no mechanism is forced to reach for an undeclared dependency (which would be a hidden seam, violating §3.3's soundness). It is the load-bearing wall for Rule 13: if a mechanism can declare a need the registry cannot resolve and validate, that is a core defect, not a config problem.

## 2. Interfaces & contracts

| Interface | Consumed by | Contract |
|---|---|---|
| **`Manifest`** (§3.1) | FND-4 loader, orchestrator (COUP-3/4) | a static, pre-construction descriptor: a mechanism's id, interface version, mesh kind, required tables (by semantic name + envelope), coupler ports, boundary ports, halt conditions it can raise, and its config-parameter schema |
| **Registry** (§3.2) | FND-4 loader | a fixed, source-ordered map `mechanism-id → {Manifest, constructor}`; deterministic to iterate; the sole way config strings become mechanisms |
| **Wiring-validation checks** (§3.3) | FND-4 loader | four total, all-errors-collected checks that a selected mechanism set is fully wired; each returns diagnoses in a deterministic order or the run halts at load |
| **Interface version** (§3.4) | every SOLV module, the core | a SemVer contract split into required + optional capabilities; modules migrate independently within a supported range |

**Invariant promised to everyone:** a mechanism may **read or emit only what its `Manifest` declares** (the *declared-or-forbidden* rule, §3.3). This is what makes config-time validation *sound* — an undeclared, dynamically-resolved dependency is invisible to the checks and is therefore prohibited.

## 3. Method

### 3.1 The `Manifest` — capability declaration *(the load-bearing contract)*

Every mechanism exposes a **static, const-evaluable `Manifest`**, returned **before construction**, so the loader can validate a config block without instantiating anything (the MOOSE `validParams()` pattern, but *data*, not a code side effect). [META-3: `moose-inputparams`] The fields (spec level; Layer 3 owns the Rust types):

| Field | Meaning |
|---|---|
| `id` | stable registry key; config selects the mechanism by this string |
| `interface_version` | the COUP-8 contract SemVer this module was written against (§3.4) |
| `grid_fields` | which components of the conserved state `U` (FND-2 §3.4) the operator **reads** and **writes** — the unified-grid analog of a mesh declaration. Every operator runs on the one 3-D grid; adaptive dimensionality (azimuthal-mode reduction) is FND-2's, not per-operator |
| `tables[]` | required offline tables, each `{semantic_name, envelope, required}` — **named**, never a pointer/index; the loader resolves the supplier (FND-5) |
| `couplers[]` | participation in the 7 coupler types, each `{coupler_kind, quantity, direction}` — `coupler_kind` and `quantity` are **closed enums**, `direction ∈ {Source, Sink, Bidirectional}` |
| `ports[]` | boundary-object ports, each `{name, role ∈ {Require, Provide}, kind}` (COUP-7) |
| `halts[]` | the closed set of halt conditions this mechanism can raise (COUP-4) |
| `params` | typed config schema: per parameter `{name, type, default, required, validity-range}` — the sub-schema FND-4 dispatches a config block against |

Design rules, each from precedent:
- **Semantic names, resolved by the orchestrator**, never pointers — a mechanism asks for `"eos.pressure"`, not a table handle; the registry binds supplier↔consumer (MOOSE `getMaterialProperty`, preCICE `data`). [META-3: `moose-inputparams`, `precice`]
- **Closed enums for couplers and quantities** (the 7 couplers of VISION_SCOPE §7.4; the conserved-quantity set) — so wiring validation (§3.3) is an *exhaustive match*, not a string compare, and a typo is a compile error, not a silent mis-wire (Modelica's closed connector vocabulary). [META-3: `modelica-connector`]
- **Required-vs-optional split** on tables/ports/capabilities — lets the loader separate a hard failure from a degraded-mode allowance (SUNDIALS generic base classes carry required + optional methods + a self-describing type query). [META-3: `sundials-generic`]
- **The table requirement carries the envelope it needs covered** — this is the project's extension beyond the frameworks surveyed: it lets §3.3 check a loaded table's *coverage* against every consumer at config time, rather than discovering a gap as a mid-solve out-of-envelope halt (FND-5 §3.5). Nothing off-the-shelf does this; it is the accountant doctrine applied to composition.
- **No hidden dependencies:** anything not in the `Manifest` may not be used. The dependency-injection lesson — a container's build-time validation is *only* as complete as the declared graph; anything resolved dynamically (service-locator style) escapes it. [META-3: `di-validate-on-build`]

### 3.2 Registry & dispatch — deterministic by construction

**Ruling: compile-time composition over a closed `enum` of the ~8 operators, plus a source-ordered `const` registry table. Do *not* use distributed-registration crates (`inventory`/`linkme`/`typetag`).**

Rationale — determinism first (META-1 §2):
- The distributed-registration crates give **no guaranteed iteration order**: `inventory` visits plugins "in any order," `linkme`'s order is link-order/platform dependent, `typetag` iterates an unspecified-order registry. That is exactly the non-deterministic ordering META-1 §2.2 forbids in any numeric or output-ordering path. A closed `enum` + `match` has a fixed, source-defined order that is bit-stable across thread counts, builds, and platforms. [META-3: `rust-registry-dispatch`, `gamer2-determinism`]
- **~8 is a closed, owned set** — the textbook case for enums over `dyn`. There is no runtime third-party plugin loading (VISION_SCOPE §6: 100% Rust runtime, no third language, no `.so` operator modules).
- **No `dyn` in hot kernels.** The orchestrator matches the mechanism enum **once per solver step** at the outer loop and then runs a concrete, monomorphized kernel; dispatch never happens per voxel/per cell (`dyn` in a tight loop is un-inlinable and ~5× slower). `enum_dispatch` is the sanctioned convenience if hand-written matches become unwieldy. [META-3: `rust-registry-dispatch`]

The registry is a single explicit table `id → {Manifest, constructor}`, either hand-written in one core file or **emitted by a build script that scans the mechanisms directory and sorts by `id` before emitting** (never trusting filesystem read order). Adding a mechanism = one new enum variant + one new row — the *only* core touch, and it is mechanical registration, not physics. This is the OpenFOAM `runTimeSelectionTable` / PETSc `SNESRegister` pattern (string-keyed constructor table, input-driven `New()`), with the ordering non-determinism engineered out. [META-3: `openfoam-runtimeselection`, `petsc-ts`]

### 3.3 Config-time wiring validation *(the fail-loud composition check)*

When FND-4's loader has parsed a config and dispatched each block to the registry, it runs these **four checks before any solve**, collecting **all** failures (never bailing on the first — the DI `ValidateOnBuild` lesson) and reporting them in a **deterministic order** (sorted by `(mechanism-id, field)`; never HashMap order). Any failure **halts at load with a diagnosis** (META-1 Principle 6). The checks are *defined* here and *invoked* by FND-4 §3.4.

1. **Structural resolve.** Every `Manifest` semantic name resolves: every required table name exists in the pinned table set; every `Require` port has a matching `Provide`; every config parameter satisfies its schema. Unknown mechanism → error listing the valid ids. (preCICE `precice-config-validate`: undefined-data / missing-mesh detection, run without launching the solve.) [META-3: `precice`]
2. **Envelope coverage.** For each required table, intersect the loaded table's covered domain (FND-5 metadata) with the **union of consumer envelopes**; a gap fails loud naming the mechanism, table, and uncovered region — *not* a silent mid-solve extrapolation. (Project extension; no framework does this for you.)
3. **Coupler source/sink balance** *(the one most easily skipped)*. Build the directed graph of coupler edges over the closed `(coupler_kind × quantity)` space; for every `Source` require ≥1 `Sink` and vice-versa; flag dangling or one-sided edges. This is the Modelica "locally balanced model" analog — the compile-time twin of the grid's every-step conservation audit (COUP-2). preCICE explicitly notes its *basic* validator does **not** check "is all necessary data actually exchanged," so this pass is built deliberately, not assumed free. [META-3: `modelica-connector`, `precice`]
4. **Interface-version.** Each `Manifest.interface_version` must lie in the core's supported range (§3.4); otherwise fail loud (`mechanism X built against interface 2.x; core requires 3.x`).

The **soundness precondition** for all four is the *declared-or-forbidden* invariant (§2): because a mechanism may use only what it declares, "fully wired" as certified here is the whole truth — there is no dynamic resolution to slip past the check.

### 3.4 Interface versioning

- `interface_version: SemVer` on every `Manifest`. The contract is split into a **stable required core** + **optional capabilities** (SUNDIALS pattern): additive capabilities are queried, so a module that doesn't provide an optional one still composes. [META-3: `sundials-generic`]
- **Additive changes bump minor**; the core supports a *range* of minors, so the 8 operators migrate **independently**, not in lockstep. **Breaking changes bump major** and are gated by check §3.3(4).
- **Never silently reinterpret a semantic name.** A renamed table/param/quantity goes through an explicit **deprecation alias with a removal date** (MOOSE `deprecateParam`/`renameParam`), emitting a warning, not a break. [META-3: `moose-inputparams`]

### 3.5 Determinism

Registry order is source-/id-defined (§3.2). Every validation pass iterates enum-variant order or an explicitly sorted key list — never a `HashMap` (config-time hash maps are permitted for construction, then baked into sorted arrays, exactly as FND-2 §3.7 does). Diagnostic ordering is deterministic. No part of composition depends on link order, wall-clock, or hash seed.

## 4. Coupling relationships

- **FND-4** dispatches each config block against a `Manifest` sub-schema and **invokes the §3.3 checks**; COUP-8 supplies the registry and the check definitions. This is the single shared seam — FND-4 owns *when/where* validation runs, COUP-8 owns *what* it means.
- **FND-2/COUP-3** run each operator's declared `grid_fields` reads/writes on the one grid inside the SDC-IMEX schedule; **COUP-2** consumes `couplers[]` to reconcile operator coupling with the conservation audit (§3.3(3)); **COUP-7** consumes `ports[]`; **FND-5** supplies the tables named in `tables[]` and their envelopes; **COUP-4** owns the halt conditions enumerated in `halts[]`.
- **Every SOLV module** supplies its concrete `Manifest` when its doc is written (deferred per the skeleton split) and implements the required interface capabilities.

## 5. Uncertainty & validity

COUP-8 introduces **no physical uncertainty** — it is composition infrastructure. Its correctness obligations (declaration completeness, sound wiring validation, deterministic registry) are **CI-gated** (VAL-3), not physically validated. Validation-ladder status: N/A (infrastructure, rung i for its tests).

## 6. Validation plan

Unit/property tests (test-first, VAL-3):
1. **Unknown mechanism / missing required param** fails loud, and the diagnosis lists the valid ids / names.
2. **Envelope gap:** a required table whose coverage misses a consumer's envelope fails at load, naming the uncovered region — no run starts.
3. **Under-wired coupler:** a `Source` with no `Sink` (and vice-versa) fails the balance pass with a precise edge diagnosis.
4. **Interface-version mismatch** fails loud with the required range.
5. **All-errors-collected:** a config with several independent faults reports them all at once, in deterministic order (not first-only).
6. **Registry determinism:** the registry's iteration order is identical across repeated builds and across platforms (guards against any distributed-registration regression).
7. **Sandbox smoke test:** a new (toy) mechanism composes from config with **no core edit beyond its one registration row** — the Rule-13 acceptance test.

## 7. References

META-3 keys: `moose-inputparams`, `openfoam-runtimeselection`, `petsc-ts`, `sundials-generic`, `precice`, `modelica-connector`, `ufl-fenics`, `rust-registry-dispatch`, `di-validate-on-build`, `gamer2-determinism`. Depends on FND-1 (`M`, `ValidityEnvelope`), FND-2 (the grid + conserved state `U`), FND-5 (tables + envelopes).

*(No open questions — all engineering decisions resolved into the body per META-1 §6.)*

## 8. Change log

| Date | Version | Change |
|---|---|---|
| 2026-07-19 | 0.1 | Initial **skeleton**. Fixed the `Manifest` capability-declaration contract (id/interface-version/mesh-kind/tables+envelopes/couplers/ports/halts/params; semantic names resolved by orchestrator; closed enums; required-vs-optional; declared-or-forbidden invariant); ruled compile-time enum dispatch + deterministic source-ordered `const` registry over distributed-registration crates (determinism); defined the four config-time wiring-validation checks (structural resolve, envelope coverage, coupler source/sink balance, interface-version) with all-errors-collected fail-loud; interface-version policy (SemVer, required+optional split, deprecation aliases). Concrete trait signatures + per-mechanism manifests deferred to the SOLV docs (W2+). |
| 2026-07-19 | 0.2 | **v1.3 unified-grid pivot.** "Mechanisms" are now the unified-grid **operators** (SOLV-1…8, 8 not 16); `mesh_kind` retired → `grid_fields` (which components of the conserved state `U` an operator reads/writes on the one grid); COUP-1 binding dependency removed (adaptive azimuthal-mode reduction is FND-2 §3.4, not a per-operator binding); coupling relationships rewired to FND-2/COUP-3 (SDC-IMEX schedule) + COUP-2 (audit). Registry/dispatch/validation contracts unchanged (dimension-agnostic). |
