# META-2 — Conventions: Docs, Notation & Code

| Field | Value |
|---|---|
| **ID** | META-2 |
| **Family** | META |
| **Status** | Draft |
| **Depends on** | META-1 |
| **Version** | 0.2 (v1.3) |

---

## 0. Purpose

The mechanical standards that make the ~36-doc catalog (and later the codebase) read as one artifact:
the doc template, symbol/units notation, cross-reference syntax, and forward-looking code
conventions. Consistency here is what lets a fresh session pick up any doc without re-learning
local dialect.

---

## 1. Doc template

Every Layer-2 doc begins with the **header table**, then the **standard sections**. Sections that
don't apply are kept with an explicit "N/A — reason" rather than deleted (so the reader knows it was
considered).

### 1.1 Header table (required fields)

| Field | Meaning |
|---|---|
| **ID** | Stable `<FAMILY>-<n>` (META-0 §1). Never renumbered. |
| **Family** | META / FND / COUP / SOLV / OFFL / VAL |
| **Status** | Draft → In Review → Reviewed → Frozen → Superseded |
| **Depends on** | IDs whose contracts this doc consumes (META-1/2/3 omitted — ambient) |
| **Version** | `major.minor`; bump minor on edits, major on contract change |
| **Skeleton/complete split** | *(only if the doc is written in two passes per META-0 §5)* what's fixed now vs later |

Owner is always Ben (project-wide), so it isn't repeated in headers. The **Depends on** field is a
light "read these first" hint, **not** a formal dependency graph — CRUCIBLE tracks no such graph
(docs are written as needed, roughly by the META-0 §5 waves).

### 1.2 Standard sections (in this order)

0. **Purpose** — the one thing this doc owns; who reads it and when.
1. **Scope & razor ruling** — what is simulated here vs treated as a boundary object, and *which
   razor* drew each line (Principle 4). Explicit in/out lists.
2. **Interfaces & contracts** — inputs consumed (from which docs/tables/ports), outputs produced,
   and the invariants promised to downstream docs. This is the part other docs bind to; it is the
   most change-controlled part.
3. **Method & governing equations** — the physics/numerics at spec level: governing equations,
   discretization scheme *named* (not coded), boundary conditions, closure choices with error bands.
   Equations are numbered and each carries a META-3 source reference.
4. **Coupling relationships** — precisely how this doc's module exchanges source terms / BCs with
   others (which of the 7 couplers, what quantity, what direction, what conservation statement).
5. **Uncertainty & validity** — the uncertainty sources this module introduces (typed per FND-1),
   its validity envelope, and its rung on the validation ladder.
6. **Validation plan** — the specific analytic/benchmark/hardware cases this module must pass
   (pointer into VAL-1/VAL-2), with target tolerances.
7. **Open questions** *(present only if there are open items)* — the format in §5 below. Resolved
   items are **deleted and integrated** into the body, never left here (META-1 §6); a doc with none
   omits this section.
8. **References** — pointers into META-3 (not free-floating citations; META-3 is the ledger).
9. **Change log** — dated lines; `major` bumps note the contract change.

### 1.3 Status lifecycle

`Draft` (being written) → `In Review` (handed to Ben) → `Reviewed` (approved; **Layer-3 may
implement**) → `Frozen` (contract stable; change needs ceremony) → `Superseded` (replaced; kept for
history, header points to successor). Downstream code must cite the doc **version** it implements.

---

## 2. Notation & symbols

- **Symbols follow the dominant field convention** for each module, declared in a per-doc symbol
  table on first heavy use. Where fields collide (e.g. `ρ` = density vs reactivity), the doc
  disambiguates locally (`ρ_m` density, `ρ_k` reactivity) and notes it.
- **Vectors** bold or arrow (`𝐁`, `\vec{v}`); **tensors** bold uppercase; **scalars** italic.
  Consistency within a doc matters more than a global choice.
- **Equations are numbered** `(SOLV-1.3)` = doc ID + local number, so other docs can cite them.
- **Uncertain quantities** are written `value ± 1σ (dist)` or `[lo, hi] @95%`, and the *distribution
  family* is named (normal, lognormal, uniform, empirical). A number with no band in a results
  context is a defect (Principle 2). "±" without a stated meaning is forbidden — always say 1σ / 95% / etc.
- **Significant figures** reflect the uncertainty: don't print digits the band doesn't support. Raw
  constants keep source precision (META-3).

### 2.1 Units convention

- **SI internally, everywhere.** kg, m, s, K, mol, A, and coherent derived units (Pa, W, J, N, T…).
- **Non-SI only at boundaries** (a citation in psia, an Isp in seconds). Converted at the edge; the
  conversion factor and its source recorded in META-3. Isp is reported in **seconds** (with g₀ =
  9.80665 m/s² stated) *and* effective exhaust velocity **`v_e`** in m/s, since the field uses both.
  Symbol convention: **`v_e`** = effective exhaust velocity; **`c`** is reserved for the speed of light;
  **`c*`** = characteristic velocity (standard rocketry).
- **eV / MeV, barns, and nuclear-data units** are permitted in the *offline nuclear* docs where they
  are the native convention, but tables handed to the runtime are converted to SI (or carry an
  explicit unit tag in FND-5 metadata). Every table column has a unit in its metadata — no bare arrays.
- **Angles in radians** internally; degrees only at config/report edges (e.g. drum angle), converted
  and labeled.
- **Temperatures in K**; never °C in computation.

### 2.2 Constants

Physical constants come from **one pinned source** (CODATA, version recorded in META-3), referenced
by name, never re-typed per doc. g₀, k_B, N_A, c, e, etc. live in META-3's constants table with
their CODATA values and are cited by ID.

---

## 3. Cross-reference syntax

- Reference another doc by **ID**: "see FND-5 §2" or the wiki-link form `[[FND-5]]` in prose.
- Reference an equation: `(SOLV-1.3)`. A datum/constant/correlation: `META-3:<entry-key>`.
- A validation anchor: `VAL-2:RL10`. A coupler: "coupler 2 (wall exchange)" matching VISION_SCOPE §7.4.
- Never cite a bare external paper inline in a design doc — cite the **META-3 entry** that holds the
  full reference. META-3 is the single bibliography, so a paper's citation lives in exactly one place.

---

## 4. Code conventions (forward-looking, for Layer 3)

Recorded now so Layer-2 interface specs are written against a known code style. Refined when coding
starts; **not** binding until then, but the determinism/units items (★) are hard requirements traceable
to META-1.

- **Language split** (Principle 11): Rust runtime, Python offline, HDF5 seam. No third language.
- **Module layout mirrors the doc families.** A Rust crate/module per FND/COUP/SOLV area; a Python
  package per OFFL pipeline. A code unit names the doc ID + version it implements.
- ★ **Determinism rules in code** (META-1 §2): deterministic reductions in numeric paths; no
  `HashMap`/`HashSet` iteration in numeric or output-ordering paths (use ordered/sorted); counter-based
  seeded RNG keyed on stable indices; no wall-clock/env/address in physics; fast-math reassociation off;
  FMA policy fixed and recorded. These are CI-enforced (VAL-3), not style suggestions.
- ★ **f64 for physics** (META-1 §3); f32 only with a doc-justified bound. **Units are type-checked at
  interface boundaries** (config parse, boundary-object declarations, table columns, results, public
  signatures) via a units type — `uom` (zero-cost for f64, base-SI normalized) as the starting choice —
  and **plain `f64` / `ndarray<f64>` inside hot kernels** working in documented SI base units (FND-1 §3.2).
- **Error handling:** `Result` at boundaries; solver hot paths return typed errors, they don't
  `panic!` — except a **halt condition**, which is a structured, reported outcome (COUP-4), not a crash.
- **No magic numbers:** every tolerance/constant is a named item sourced to config or META-3.
- **Tests co-located and test-first** (Principle 10): analytic/manufactured-solution tests before
  benchmark tests; determinism test in CI; conservation-audit test in CI (VAL-3).
- **`unsafe` requires a written justification** in-line and is avoided in physics paths.
- **Provenance in outputs** (Principle 7): every results bundle embeds config, seed, table versions,
  and build fingerprint (FND-6).
- **Precision-critical idioms:** prefer compensated/ordered summation for large reductions where the
  error budget demands it; prefer well-conditioned formulas (META-1 §3); document any place where
  catastrophic cancellation was avoided.

---

## 5. Open items

*(v1.3, 2026-07-19 — the in-doc "needs-Ben" register is purged.)* **Design / creativity forks that
need Ben are never parked in a doc.** They are raised to him **before the doc is written or during a
review pass**, batched, and integrated into the body once resolved (META-1 §6). There is no in-doc
question register for them and no `needs-Ben` tag anywhere. Engineering / implementation /
representation choices are resolved by Claude and written straight into the body.

The **only** open item a doc may carry is a **`needs-analysis`** TODO — something whose closure
requires an experiment or offline study — recorded inline where it applies (or in a short **Open
items** section) as `<DOCID>-Q<n>` · the question · the analysis that will close it · what it blocks.
On resolution it is deleted and the finding written into the body.

---

## 6. File & repo conventions

- Docs live under `docs/<family>/<ID>_<slug>.md` (lowercase-kebab slug).
- One doc = one file. Diagrams inline as Mermaid (with a textual fallback, since terminal markdown
  may not render Mermaid — see META-0 §4).
- Dates are absolute ISO `YYYY-MM-DD` (never "today"/"last week").
- Tables/data destined for the runtime are **not** committed as loose CSVs in docs; they are described
  here and generated by the OFFL pipelines into versioned HDF5 (FND-5). Small illustrative numbers in
  a doc are fine but are cited to META-3.

