# FND-4 — Config Schema & Run Manifest

| Field | Value |
|---|---|
| **ID** | FND-4 |
| **Family** | FND (Foundations / spine) |
| **Status** | Reviewed (2026-08-14) — **skeleton** (per META-0 §5, W1) |
| **Depends on** | FND-1, FND-3, FND-5, COUP-8; FND-2 (§3.4 θ-ladder) |
| **Version** | 0.3 (2026-08-24 plan S7: `[operating_profile]` gains the COUP-4 §3.1 commanded-profile grammar + the COUP-7 §3.2.2 `valve_cited_timeline_s` compression declaration; three new registry mechanisms (`combustion_blend`, `spark_igniter`, `structural_margins`) ride the type-keyed grammar unchanged). 0.2 (2026-08-14 review fix wave) |
| **Skeleton/complete split** | **Fixed now:** the top-level TOML document structure; the **config-vs-manifest split**; the loader pipeline and where validation runs; the config-derived part of the run-manifest field list; the no-hidden-defaults policy; schema versioning + migration; the registry-driven sub-block dispatch that makes new mechanisms pure data. **Deferred** (filled as dependents land, W2+): the per-mechanism / per-material / per-boundary-object parameter blocks (each supplied by its own `Manifest`, COUP-8 §3.1, and its SOLV/OFFL/COUP-7 doc); the complete operating-profile grammar; the complete UQ-settings block (co-designed with COUP-5); the complete **`[engine]` composition grammar** — with its **explicit Require→provider binding semantics fixed now** (§3.2, O20; deferred is the full grammar, never the matching rule). |

---

## 0. Purpose

**Config is the single source of a run's intent** (META-1 §5): one TOML file, plus pinned table versions and a seed, **fully determines a run** — nothing depends on anything outside `{config, tables, seed, build fingerprint}`, and there are no hidden defaults. This doc fixes the **config schema** the author writes and the **run manifest** the loader emits: the machine-captured, fully-resolved, content-addressed record that is the actual regeneration key for S6.

FND-4 is the sandbox's data-entry surface (Rule 13): geometries, materials, mechanism selections, operating profile, and UQ settings are **pure data** composed through this one path. Its load-bearing decisions are (a) the **config→manifest split** and (b) the **registry-driven sub-schema** mechanism by which a new mechanism adds a config block with *zero* edits to the top-level schema.

Read after FND-1 (`UncertainInput`, `ValidityEnvelope`), COUP-8 (the registry it dispatches to), FND-5 (table pins). Its output (the manifest) is consumed by **FND-6**, which adds the build fingerprint and bundles it with results.

## 1. Scope & razor ruling

**Owns:** the top-level TOML **document structure**; the **loader pipeline** (parse → migrate → structural deserialize → registry dispatch → wiring validation → default resolution → manifest emission); the **config-derived portion** of the run manifest; the **no-hidden-defaults** policy; **schema versioning & migration**; the single-source-of-intent rule.

**Defers (owner named):**
- The **per-mechanism/material/boundary-object parameter schemas** → each is a **`Manifest`** (COUP-8 §3.1) supplied by the owning SOLV/OFFL/FND-7/COUP-7 doc. FND-4 owns only the *document* the blocks live in and the *dispatch* to those schemas.
- The **wiring-validation checks** (unknown mechanism, missing table, unsatisfied coupler) → **COUP-8 §3.3** *defines* them; the FND-4 loader *invokes* them.
- The **geometry grammar** (CSG tree, STL reference) → **FND-3** (config carries it; FND-3 owns its syntax and voxelization).
- The **uncertainty representation** (`UncertainInput`, distribution families) → **FND-1**; the **UQ ensemble/sampler** it feeds → **COUP-5**.
- The **table schema, version-pin resolution and content-hashing** → **FND-5** (config states the pin; the loader records the resolved hash via FND-5).
- The **build fingerprint** and the **results bundle** the manifest is embedded in → **FND-6**.

**Razor ruling:** pure infrastructure — no physics. Its fidelity obligation is that config be *complete enough* that a run's entire intent is captured in one file (no ambient state, no hidden default), which is what makes S6 achievable.

## 2. Interfaces & contracts

| Interface | Consumed by | Contract |
|---|---|---|
| **Author config** (TOML) | the loader | one human-authored file = a run's full intent; fixed top-level grammar + registry-dispatched blocks; carries `schema_version` |
| **Loader pipeline** (§3.4) | the orchestrator (COUP-3/4) | parse → migrate → validate (COUP-8) → resolve defaults → produce a fully-specified in-memory config with **zero implicit values** |
| **Run manifest (config part)** (§3.6) | FND-6 | the resolved config + table version-pins & content-hashes + master seed + RNG algorithm + schema versions — the config-derived slice of the S6 regeneration key |

**Invariant promised to everyone:** a run's behavior is a pure function of `{resolved config, pinned+hashed tables, seed, build fingerprint}`. No environment variable, CLI overlay, wall-clock, or ambient default may change a result — the file is the single source of intent, and every effective value it resolves to is recorded in the manifest.

## 3. Method

### 3.1 Two artifacts: author config vs run manifest *(the load-bearing split)*

There are **two** artifacts, and conflating them breaks reproducibility:
- **Author config** = *intent*: human-authored, human-diffable, intentionally underspecified (it may omit defaulted parameters, say `table = "eos@latest"`, or leave the seed to be auto-drawn), versioned by `schema_version`. This is what Ben edits and what the batch sweeper templates.
- **Run manifest** = *the regeneration key*: the **fully-resolved** config (every effective default materialized — §3.5), with the table pins resolved to explicit versions **and content hashes**, the master seed and RNG algorithm fixed, and (added by FND-6) the build fingerprint. You **regenerate from the manifest**, never from the author config alone.

This mirrors every serious reproducibility tool: Cargo's `Cargo.toml` (ranges) vs `Cargo.lock` (exact, hashed, format-versioned); DVC's `dvc.lock` (content-hashed deps *and* outputs); conda-lock/Snakemake's platform-stamped fully-solved pins; Sacred's auto-emitted final `config.json` (with the seed as part of the record). [META-3: `cargo-lock`, `dvc-lock`, `conda-lock`, `sacred-config`, `mlflow-run`]

### 3.2 Top-level TOML structure *(fixed grammar; bodies are registry-dispatched)*

The parser knows only these **section kinds** and the `type = "<registry id>"` key convention; it **never enumerates mechanisms**:

```
schema_version = 1

[meta]                  # run name, description (non-load-bearing)
[geometry]              # CSG tree OR { stl = "ref@version" } — grammar owned by FND-3
[materials.<name>]      # type = "<FND-7 id>", + that material's params
[mechanisms.<name>]     # type = "<COUP-8 registry id>", params, table pins used
[couplers.<name>]       # inter-mechanism coupling selections (7 kinds, COUP-2)
[engine]                # composition: which geometry/materials/mechanisms/couplers wire together,
                        #   incl. [engine.bindings]: explicit Require -> provider-instance port
                        #   bindings (semantics fixed now, §below; full grammar deferred)
[operating_profile]     # time / throttle / boundary schedule (grammar deferred)
[determinism]           # mode = "fixed-order" | "relaxed" (O21; checked at load, recorded in manifest)
[tables]                # logical_name -> version pin (content hash resolved into manifest)
[uq]                    # sampler, N, declared distributions over knobs (block deferred, COUP-5)
[rng]                   # master_seed, algorithm (seed may be auto-filled into the manifest)
```

**`[engine]` binding semantics (O20 — fixed now, ahead of the full deferred grammar):** every boundary/coupler
`Require` port of every selected mechanism instance is satisfied by an **explicit named binding** in
`[engine.bindings]` — `"<instance>.<require_port>" = "<provider instance>.<provide_port>"`. **Nothing
auto-matches by kind or name:** a missing or ambiguous binding (two same-kind `Provide`s, no binding) is a
**load error listing the candidate providers** (COUP-8 §3.3 check 1 defines the resolution rule; this block
is where the bindings live).

**`[determinism]` block (O21):** `mode ∈ {"fixed-order", "relaxed"}` — the S6/META-1 §2.1 regime→guarantee
declaration surface. `fixed-order` (the default) = byte-identical fixed-order reductions everywhere.
`relaxed` = the declared relaxed-reduction (GPU) path, admissible **only for non-chaotic runs**: at load the
checker consults every selected mechanism's `Manifest.chaotic_class` (COUP-8 §3.1) for the configured
regime, and **REFUSES the run** — load error naming the chaotic-classified mechanisms — if any is `chaotic`
while `mode = "relaxed"`. The resolved mode is recorded in the manifest (§3.6); the refusal is what makes
S6's "fixed-order mandatory in chaotic regimes" mechanical rather than disciplinary.

The top level uses typed structs with `deny_unknown_fields`, so an unknown top-level key is a **hard error**, not a silent ignore. This is the MOOSE/OpenFOAM/Cantera pattern — a fixed generic grammar + a `type` key dispatched to a per-module subschema, with unknown options rejected loudly rather than ignored. [META-3: `moose-inputparams`, `openfoam-runtimeselection`, `cantera-yaml`]

### 3.3 Registry-driven sub-schemas *(the Rule-13 extension surface)*

Each `[mechanisms.*]`, `[materials.*]`, and boundary-object block declares a `type`; the loader looks that id up in the **COUP-8 registry** and validates the block against **that mechanism's `Manifest.params`** (COUP-8 §3.1) — and only that mechanism's params. There is **no monolithic hand-edited schema**: the effective schema is the union of what registered mechanisms declare. Adding a mechanism therefore adds its `Manifest` (one registration in COUP-8) and *nothing* in FND-4 — no top-level grammar edit, no config-path change. This is precisely what makes "define a new engine = pure data, no new code" true.

### 3.4 The loader pipeline

Ordered, and it either produces a fully-specified config or **halts at load with all diagnoses** (never a partial run):

1. **Parse** the TOML (span-annotated errors via `toml`/`miette`).
2. **Migrate** if `schema_version` < current, via the ordered migration chain (§3.7).
3. **Structural deserialize** into typed top-level structs (`deny_unknown_fields`) — the compile-time-checked layer.
4. **Dispatch** each `type`-keyed block to its registry `Manifest` sub-schema (§3.3).
5. **Wiring validation** — invoke COUP-8 §3.3's four checks (structural resolve — incl. the explicit `[engine.bindings]` resolution, O20 — envelope coverage, coupler source/sink balance, interface-version). **Collect all errors** and report together in deterministic order. The same pass runs the **FND-4-owned load-time value checks**: (a) **determinism refusal (O21)** — `[determinism].mode = "relaxed"` with any selected mechanism chaotic-classified for the configured regime (`Manifest.chaotic_class`, COUP-8 §3.1) ⇒ **refuse to load**, naming the mechanisms (§3.2); (b) **θ-ladder alignment** — the declared **`N_θ^max` must equal `4·2^n`, `n ≥ 0`**, so the factor-2 coarsening ladder lands **exactly on the guard resolution `N_θ^guard = 4`** (FND-2 §3.4); any other value ⇒ load error stating the nearest admissible values.
6. **Resolve defaults** (§3.5) → an in-memory config with zero `None`.
7. **Emit the manifest** (§3.6).

Two-layer validation is deliberate: **serde/toml stops at the first error** and cannot express cross-field or registry constraints, so structural checks (steps 3–4) use serde while **semantic checks (step 5) use an error-accumulating pass** (`garde`-style) so the author sees *every* problem at once, each annotated to its TOML span. Single-file `serde + toml` is the whole stack; **layered/overlay config loaders (`figment`/`config-rs`) are rejected** — their merge/override precedence reintroduces exactly the hidden-precedence ambient state that "one file = single source of intent" forbids. [META-3: `garde-validation`]

### 3.5 No hidden defaults

Defaults live **in the registry descriptors** (each `Manifest.params` carries its default), never as scattered `unwrap_or` in code. Step 6 resolves them into a config where **nothing is implicit**, and that resolved config is serialized into the manifest (§3.6). The doctrine is enforced *mechanically*: a test asserts the resolved config **round-trips and contains zero implicit values**; a `--print-resolved-config` command emits it without running. (MOOSE `--show-input`/`--json` dumps every param's effective default; Sacred writes the final resolved config.) [META-3: `moose-inputparams`, `sacred-config`]

### 3.6 Run manifest — the config-derived regeneration key

The manifest records (FND-4's portion; FND-6 adds the build/environment portion and does the bundling):

- `manifest_schema_version`;
- `config_content_hash` (hash of the canonicalized author TOML);
- `resolved_config` (the full config, all defaults materialized — §3.5);
- per pinned table: `{logical_name, resolved_version, content_hash, producing_generator_version}` (hashes resolved via FND-5's loader);
- `master_seed` + `rng_algorithm` (+ crate/version + the seed-derivation/stream policy — the counter-based key of FND-1 §3.5);
- `determinism_mode` (`fixed-order` | `relaxed`, O21) + the per-mechanism `chaotic_class` values in force for the configured regime (the record behind the load-time refusal);
- `schema_version` (author) and the migration chain applied, if any.

*(The build fingerprint — rustc/target/FP-FMA policy/libm/lockfile hash — is **FND-6 §3.4**'s field list; the complete manifest is FND-4-part + FND-6-part, emitted into every results bundle.)* Table pins carry a **content hash, not just a version label**: a re-uploaded "v3" with the same label but different bytes must fail (DVC/conda-lock lesson). [META-3: `dvc-lock`, `conda-lock`]

### 3.7 Schema versioning & migration

A mandatory `schema_version` heads both the author config and the manifest — config outlives any one binary. The loader runs an **ordered chain of pure migration passes** (`v1→v2→v3`, incremental, each small and tested, never skipped) to bring an old config to current *before* validation. It **refuses loudly** a version newer than the build understands, or older than the oldest supported (Kubernetes CRD versioning; Cargo.lock's internal format version). [META-3: `k8s-crd-versioning`, `cargo-lock`]

### 3.8 Determinism

Config hashing is over a **canonicalized** serialization (stable key order), so a semantically identical config always hashes identically. The seed is carried in the manifest, never drawn from wall-clock at run time. No env var or CLI overlay may alter a result (§3.4 rejection of overlay loaders) — the file is the sole intent, which is what lets FND-6 assert bit-identical regeneration from the manifest.

## 4. Coupling relationships

- **COUP-8** supplies the registry and the wiring-validation check definitions; FND-4's loader dispatches blocks to `Manifest.params` and invokes the checks (§3.4 step 5). Single shared seam: FND-4 owns *when/where*, COUP-8 owns *what*.
- **FND-5** resolves each `[tables]` pin to a version + content hash and enforces envelopes; FND-4 records the result in the manifest.
- **FND-3** owns the `[geometry]` grammar; **FND-7** the `[materials.*]` types; **COUP-5** the `[uq]` block; **COUP-7** the boundary-object blocks — each a deferred sub-schema per the skeleton split.
- **FND-6** consumes the manifest, adds the build fingerprint, and embeds it in the results bundle (the S6 chain).
- **The orchestrator (COUP-3/4)** consumes the resolved config to construct and wire the selected mechanisms.

## 5. Uncertainty & validity

FND-4 introduces **no physics uncertainty**; it carries the *declarations* of uncertainty (the `[uq]` block feeds FND-1/COUP-5) but originates none. Its correctness (no hidden defaults, sound migration, canonical hashing, single-source-of-intent) is **CI-gated** (VAL-3). Validation-ladder status: N/A (infrastructure).

## 6. Validation plan

Unit/property tests (test-first, VAL-3):
1. **Config↔manifest round-trip:** the resolved config round-trips and contains **zero implicit values**.
2. **Unknown key / unknown mechanism** fails loud (`deny_unknown_fields`; registry miss lists valid ids).
3. **All-errors surfaced:** a multi-fault config reports every fault at once, span-annotated, in deterministic order (guards against serde first-error truncation).
4. **Table-pin content-hash mismatch:** a same-label/different-bytes table fails the manifest check.
5. **Migration:** an old-`schema_version` config migrates to current and validates; a future version is refused loudly.
6. **Single source of intent:** setting an environment variable or CLI overlay does **not** change a result — only the file does.
7. **Determinism refusal (O21):** a `relaxed`-mode config selecting a chaotic-classified mechanism is refused at load, naming the mechanisms; the same config in `fixed-order` mode loads.
8. **θ-ladder alignment:** `N_θ^max = 24` is refused (names 16 and 32); `N_θ^max = 32` loads.
9. **Binding errors (O20):** an unbound `Require` with two same-kind `Provide`s is refused listing both candidates (exercises COUP-8 §3.3(1) through the loader).

## 7. References

META-3 keys: `moose-inputparams`, `openfoam-runtimeselection`, `cantera-yaml`, `cargo-lock`, `dvc-lock`, `conda-lock`, `sacred-config`, `mlflow-run`, `garde-validation`, `k8s-crd-versioning`. Depends on FND-1 (`UncertainInput`, `ValidityEnvelope`), FND-5 (table pins/hashes), COUP-8 (registry + validation), FND-3 (geometry grammar).

*(No open questions — all engineering decisions resolved into the body per META-1 §6.)*

## 8. Change log

| Date | Version | Change |
|---|---|---|
| 2026-08-24 | 0.3 | **Plan S7 (landed with the code) — the startup config faces.** The `[operating_profile]` subset gains the **COUP-4 §3.1 commanded-profile grammar**: `commanded_p_c_pa` / `commanded_thrust_n` arm the WORKS criterion; `eps_works` (default 0.02) and `t_dwell_flowthroughs` (default 20) **materialize only when armed** — absent-for-absent replay preserved; declaring criterion constants without a commanded quantity is a **load error**. Plus the **COUP-7 §3.2.2 compression declaration `valve_cited_timeline_s`** (declaration-only, manifest-recorded — the cited physical timeline the compressed schedule stands in for). Three new registry mechanisms — **`combustion_blend`** (with `chem_unburnt` + `chem_ignition` table pins and the wrinkling param), **`spark_igniter`**, **`structural_margins`** — ride the existing type-keyed grammar **unchanged** (new mechanism = pure data, the Rule-13 extension surface working as designed). All new blocks are manifest-recorded for free through the embedded resolved config. |
| 2026-08-14 | 0.2 | **Review fix wave (O20, O21) + θ-ladder check.** §3.2/header split: **`[engine]` composition grammar** added to the deferred list with its **explicit Require→provider binding semantics fixed now** — `[engine.bindings]` named bindings, no auto-matching, missing/ambiguous = load error listing candidates (O20). §3.2/§3.4(5): **`[determinism]` block** (`mode = fixed-order \| relaxed`), manifest-recorded; load-time check **refuses** relaxed-mode runs selecting any chaotic-classified mechanism (`Manifest.chaotic_class`, COUP-8 §3.1) — the S6 guarantee made mechanical (O21). §3.4(5): **θ-ladder alignment check** — `N_θ^max = 4·2^n` enforced at load so the coarsening ladder lands exactly on `N_θ^guard = 4` (FND-2 §3.4). §3.6: `determinism_mode` + in-force `chaotic_class` values recorded in the manifest. §6 items 7–9. |
| 2026-07-19 | 0.1 | Initial **skeleton**. Fixed the config-vs-manifest split (author intent vs fully-resolved, content-addressed regeneration key); the fixed top-level TOML grammar with `type`-keyed registry dispatch (`deny_unknown_fields`); registry-driven sub-schemas as the Rule-13 extension surface (new mechanism = pure data, zero top-level edit); the loader pipeline with two-layer validation (serde structural + error-accumulating semantic pass invoking COUP-8 §3.3), rejecting overlay loaders to preserve single-source-of-intent; no-hidden-defaults (defaults in registry descriptors, resolved into the manifest, mechanically tested); the config-derived manifest field list (content-hashed table pins, seed, RNG algo); schema versioning + ordered incremental migration. Per-module parameter blocks, operating-profile & UQ grammars deferred (W2+). |
