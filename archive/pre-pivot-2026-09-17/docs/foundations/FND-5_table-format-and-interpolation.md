# FND-5 — Table Format, Loader & Interpolation

| Field | Value |
|---|---|
| **ID** | FND-5 |
| **Family** | FND (Foundations / spine) |
| **Status** | Reviewed (2026-08-14) |
| **Depends on** | FND-1 |
| **Version** | 0.3 |

---

## 0. Purpose

The **only cross-language seam** in the project (VISION_SCOPE §6): versioned HDF5 tables produced offline
(Python: OpenMC/Cantera/Geant4/…) and read at runtime (Rust, `hdf5` crate). This doc fixes the **table
schema**, the **loader** (provenance, versioning, validity-envelope enforcement), and the **N-D
interpolation** with a quantified error budget. Everything the runtime knows about nuclear, chemical,
EOS/opacity, stopping, and kernel data comes through here.

## 1. Scope & razor ruling
**Owns:** the HDF5 group/attribute schema; the load-time contract (version pins, provenance, envelope); the
runtime N-D interpolation and its error budget. **Defers:** how each table is *generated* → the OFFL docs;
the uncertainty *type* it carries → FND-1; the material tables' content → FND-7. One uniform table contract
for all data (Rule 12): a cross-section table and an EOS table obey the same schema and loader.

## 2. Interfaces & contracts
| Interface | Consumed by | Contract |
|---|---|---|
| **Table schema** (§3.1) | all OFFL producers, the loader | group-per-table + explicit axes + mandatory metadata |
| **Loader** (§3.2) | orchestrator (config-time) | verify version pins + provenance + envelope; refuse on mismatch |
| **Interpolate(query)** (§3.3) | every solver reading a table | a deterministic **scalar** value, or a hard `OutOfEnvelope` — never a distribution (META-1 §4 pure outer loop) |
| **Provenance handle** | FND-6, COUP-6 | table semantic-version + content hash → walkable to META-3 keys |

**Invariant:** interpolation returns **scalars only** — the inner run carries no distribution arithmetic
(META-1 §4). The table's physical uncertainty enters the ensemble by **outer-loop sampling** (FND-1: a
sampled realization / perturbed sample-set member per ensemble member), and the declared
`interp_error_bound` is carried by COUP-5 as an **epistemic interval** into the p-box (COUP-5 §3.3) — never
RSS'd into a runtime return value. No value is *ever* returned outside the declared validity envelope.

## 3. Method

### 3.1 HDF5 table schema
- **One group per physical table** (e.g. `/xs/U235/fission`), holding **explicit coordinate-axis datasets**
  (1-D arrays per parameter — never an implied linspace; robust to non-uniform grids) + one value dataset
  per quantity + a sibling `sigma_<name>` uncertainty dataset (or an attribute for a scalar band).
- **Chunk + compress** value datasets, chunk shape aligned to the interpolation stencil (a query touches ≤1–2
  chunks).
- **Mandatory attributes** (a table missing any cannot load):
  - **Provenance:** producer tool + version, input-deck hash, source-library id (e.g. `ENDF/B-VIII.1`),
    generator git commit, RNG seed (if stochastic) — the same set that flows into META-3 and the results
    bundle (data-hygiene, META-1 Principle 7).
  - **Versions:** `schema_version` and `data_version` (semantic). The loader **refuses a major mismatch**.
  - **Validity envelope:** per-axis `[min,max]` (+ optional convex-hull/mask for scattered data).
  - **Per-column semantics:** `units` (SI or tagged), `interp_rule` (e.g. `log-log`, `lin-lin`, `log-T`) so
    the runtime interpolates in the space the producer intended, and `interp_error_bound` (§3.4).
- **Empirical sample-sets** (SANDY-perturbed nuclear-data realizations, FND-1 §3.9) are carried as an HDF5
  dataset: one row per realization, joint across correlated nuclides. [META-3: `hdf5-schema`]

### 3.2 Loader & version pinning
Config pins each table's `data_version` (+ content hash) (S6). At load the loader: verifies the pin and hash;
checks schema major-version; records the full provenance into the run manifest (FND-4/FND-6); and validates
the envelope metadata is present. Any failure **halts** with a diagnosis — the loader enforces
reproducibility by construction, not by discipline (META-1 §5).

### 3.3 Interpolation
- **Default = multilinear on regular grids** for any conservation- or positivity-critical quantity
  (cross-sections, opacities, densities, reactivity): it **cannot overshoot** (stays within cell bounds), is
  cheapest, and is trivially bit-reproducible.
- **Monotone cubic (Steffen / Fritsch–Carlson PCHIP)**, tensor-product, where C¹ smoothness is needed
  (reactivity coefficients and EOS quantities that get differenced) — monotonicity/positivity preserved per
  interval, unlike natural cubic splines (which overshoot) or RBF (ill-conditioned, non-deterministic).
- **Scattered data is pre-resampled to a regular grid offline**; RBF is not used at runtime.
- **Thermodynamic-potential tables (spine EOS, S11):** a table declared `kind = thermo_potential` stores a
  **Helmholtz free energy F(ρ,T) per species**, with its tabulated partial derivatives, interpolated by
  **tensor-product Hermite (biquintic-class) interpolation** (the Timmes–Swesty precedent). p, e, s, and
  sound speed are obtained by **fixed-order analytic differentiation of that one interpolant** — never
  independently interpolated — so Maxwell-relation consistency holds by construction. The generator writes a
  mandatory **`thermo_audit` metadata block**: the offline convexity / sound-speed-positivity audit over the
  declared envelope; a potential table without a passing audit **cannot load**. (The physical requirement and
  the quantity set are FND-7 §3.2's; this bullet owns only the schema + interpolation rule.) [META-3:
  `helmholtz-table`]
- **Interpolate in the linearizing space** named by `interp_rule` (log-log for XS, log-T for opacity) to
  shrink curvature and hence error.
- **Coefficients are precomputed offline** and shipped in the table; the runtime does stencil-fetch + fixed-
  order Horner evaluation only. [META-3: `interp-monotone`, `table-precedent`]

### 3.4 Interpolation-error budget (fed to UQ, never assumed negligible)
Interpolation error is bounded and knowable (multilinear ≈ (h²/8)·max|∂²f|; cubic ≈ O(h⁴)). **Offline**, the
generator measures the actual error per table (grid refinement or holdout) and writes `interp_error_bound`;
the grid spacing is *chosen* offline to hit a target budget. The bound is **declared table metadata, not a
runtime return value**: COUP-5 consumes it as an **epistemic interval** in the outer loop (§3.3 there), so it
**propagates into the p-box** without ever being RSS'd with the physical band or attached to an inner-loop
scalar. This is the structural defense against silent interpolation error corrupting a headline result. [META-3:
`interp-error-budget`]

### 3.5 Validity-envelope enforcement
Interrogating outside the envelope returns a hard **`OutOfEnvelope`** — never a clamped or extrapolated value.
Policy is `{flag | refuse}` per config, default **refuse** for headline runs (META-1 Principle 8); COUP-5 records
per-member envelope hits so a partially-out-of-envelope ensemble is reported honestly.

### 3.6 Determinism
Fixed FP evaluation order; FMA policy pinned (META-1 §2); no `-ffast-math`. Multilinear/Horner in fixed order
is bit-reproducible on a target. **The runtime never uses stochastic table interpolation** (e.g. OpenMC's RNG-
sampled temperature interpolation is for MC transport, not a deterministic runtime — use true weighted
interpolation). [META-3: `table-precedent`]

## 4. Coupling relationships
- **OFFL-1…6** produce tables to this schema (incl. OFFL-5, the constitutive-spine generator); a producer
  that can't supply provenance+envelope+uncertainty cannot ship a table.
- **Every table-reading solver** calls §3.3; it receives a scalar value or a halt, never a guess —
  uncertainty rides the outer loop (COUP-5).
- **FND-6** records the pinned versions/hashes for regeneration (S6); **COUP-6** scores pedigree from the
  provenance chain.
- **FND-1** supplies the uncertainty type; **FND-7** material tables are a client of this schema.

## 5. Uncertainty & validity
FND-5 introduces **interpolation error** (budgeted, §3.4) and enforces validity (§3.5); it carries but does
not originate physical uncertainty (that is the producing OFFL doc's). Infrastructure rung; CI-gated.

## 6. Validation plan
1. **Round-trip:** write→load→read reproduces values and metadata exactly.
2. **Interpolation accuracy & error bound:** interpolated values meet the stored `interp_error_bound` on
   analytic test functions; the bound is not exceeded.
3. **Monotonicity/positivity:** multilinear and PCHIP never overshoot on steep-gradient test tables.
4. **Envelope refusal:** out-of-envelope queries return `OutOfEnvelope`; flag/refuse honored per config.
5. **Version/provenance:** major-version mismatch and hash mismatch both halt; provenance round-trips to the
   results bundle.
6. **Determinism:** identical interpolation results across threads/runs.
7. **Potential-table consistency (S11):** p/e/s/a derived from a `thermo_potential` table satisfy the Maxwell
   relations to round-off; a² > 0 everywhere the `thermo_audit` asserts; a table with a failing audit is
   refused at load.

## 7. References
META-3 keys: `hdf5-schema`, `interp-monotone`, `interp-error-budget`, `table-precedent`, `helmholtz-table`.
Depends on FND-1 (uncertainty type, empirical sample-sets).

## 8. Change log
| Date | Version | Change |
|---|---|---|
| 2026-08-18 | 0.3.1 | **Rule-space error bound (session-12 review).** §3.4's measured bound gains a companion for LOG-valued columns: the producer additionally measures the holdout error in the value's own rule space (`interp_error_bound_log`, \|Δ ln value\| ≈ relative) — an absolute bound attained at the dense end of a log-valued range is scale-blind at the sparse end (found: the equilibrium density acceptance admitted states 10⁴× the local ρ off-surface at the plume fringe). Holdout sweeps now include envelope-EDGE points (a declared envelope cutting through an axis cell previously left that cell's in-envelope slab structurally unsampled). Consumers use the rule-space bound when present, the absolute bound otherwise (pre-existing artifacts unchanged). **Recorded deferral:** the new attr rides outside digest v3; digest v4 folds it in. |
| 2026-08-14 | 0.3 | **Post-review fix wave (S11).** §3.3 gains the **thermodynamic-potential table rule**: `kind = thermo_potential` stores Helmholtz F(ρ,T) per species + derivatives, tensor-product Hermite (biquintic-class, Timmes–Swesty) interpolation, p/e/s/sound-speed by fixed-order analytic differentiation of the one interpolant (Maxwell consistency by construction); mandatory `thermo_audit` (offline convexity/sound-speed-positivity) metadata — no passing audit, no load. §6 consistency-validation item added. Physical requirement owned by FND-7 §3.2 (clean split, no restating). |
| 2026-08-13 | 0.2 | **Consistency sweep: pure-outer-loop reconciliation.** Interpolate() returns scalars only; the table's physical band enters by outer-loop sampling and `interp_error_bound` is consumed by COUP-5 as an epistemic interval, never RSS'd into a runtime return value (§2 row + invariant, §3.4, §4 reworded to match META-1 §4 / FND-1 §3.4b / COUP-5 §3.3). |
| 2026-07-14 | 0.1 | Initial draft. HDF5 group-per-table schema with mandatory provenance/version/envelope/units/uncertainty; loader with version-pin + hash + envelope enforcement (halt on mismatch); multilinear-default / monotone-cubic interpolation with offline-measured error bound RSS'd into UQ; refuse-out-of-envelope; deterministic fixed-order evaluation; scattered→regular offline (no runtime RBF/stochastic interpolation). |
