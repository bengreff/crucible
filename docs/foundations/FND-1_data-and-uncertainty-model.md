# FND-1 — Data & Uncertainty Model

| Field | Value |
|---|---|
| **ID** | FND-1 |
| **Family** | FND (Foundations / spine) |
| **Status** | Draft (no open questions) |
| **Depends on** | — (root of the spine; META-1/2/3 ambient) |
| **Version** | 0.3 (v1.3 p-box) |

---

## 0. Purpose

FND-1 defines the **atoms** every other doc uses: how a physical quantity carries its **units**,
its **uncertainty**, and its **provenance**, and how the project's headline promise — *"every result
is a distribution, never a bare number"* (S3) — is represented from input declaration through to the
results bundle. It is the contract that makes the blessed architectural commitment concrete:
**UQ is a pure outer loop; the inner simulation carries no distribution arithmetic** (META-1 §4).

Read this before FND-2 (which stores these quantities per cell), FND-5 (tables that expose them),
COUP-5 (which samples them), and FND-6 (which reports them).

---

## 1. Scope & razor ruling

**FND-1 owns (the types and their semantics):**
- The **numeric substrate** for physics (§3.1).
- The **units strategy** and the canonical internal unit system (§3.2).
- The **`Quantity`** notion (a physical value + its unit + optional provenance).
- The **`MediumState` vector `M`** — the canonical local-state argument every constitutive law takes
  (§3.3). This is the data-model expression of Rule 12 (one law over `M`, no `if(material)` branch).
- The three uncertainty representations: **`UncertainInput`** (declaration), **runtime scalar**,
  **`ResultDistribution`** (output) (§3.4).
- The **sampling contract** (how a declaration becomes one member's concrete value) and the
  **deterministic RNG contract** (§3.5) — *defined here, executed by COUP-5*.
- The **`ProvenanceRef`** and **`ValidityEnvelope`** types (§3.6, §3.7).
- The **finiteness / NaN policy** for physics values (§3.8).

**FND-1 does *not* own (deferred, with the owner named):**
- Storage of quantities on the grid → **FND-2** (FND-1 says what a per-cell value *is*; FND-2 says how the grid holds it).
- On-disk table schema & interpolation → **FND-5** (FND-5 *exposes* `UncertainInput`s and enforces `ValidityEnvelope`s).
- The ensemble loop, correlated joint sampling, Sobol → **COUP-5** (FND-1 defines the contract it obeys).
- The results bundle file format → **FND-6** (FND-1 defines the `ResultDistribution` *content*).
- Pedigree computation → **COUP-6** (FND-1 provides the provenance handle it scores).

Razor ruling: FND-1 is pure infrastructure — no physics is simulated here. Its "completeness" obligation
(fidelity doctrine) is that the *representation* be complete enough that no downstream module is forced
to invent an ad-hoc uncertainty or provenance path (which would be a seam).

## 2. Interfaces & contracts (what downstream binds to)

The types below are described at **spec level** (fields + semantics + invariants), not as Rust code
(Layer 3 owns the concrete `struct`s). Names are indicative.

| Type | Consumed by | Produced by | One-line contract |
|---|---|---|---|
| `Quantity` | all solvers | config, tables, boundary objects | a finite f64 in a known unit; unit handling per §3.2 |
| `MediumState M` | every operator (SOLV-1…4) via the constitutive spine (FND-7) | FND-2 (per cell, reconstructed from `U`) | the *only* state argument a material/medium law may read — no law reads a material label |
| `UncertainInput` | COUP-5 (sampler) | tables (FND-5), boundary objects (COUP-7), config (FND-4) | declares how one scalar varies + its provenance + validity |
| runtime scalar (f64) | inner solvers | COUP-5 per member | one deterministic realization; no ± attached |
| `ResultDistribution` | FND-6, COUP-6, plots | COUP-5 (ensemble reduce) | empirical member values + summaries + provenance + pedigree handle |
| `ProvenanceRef` | COUP-6, FND-6 | META-3 keys, table loader | walkable link result → table → source (§3.6) |
| `ValidityEnvelope` | FND-5, COUP-5 | tables, boundary objects | domain of validity + flag/refuse policy |

**Invariant promised to everyone:** any value that reaches a `ResultDistribution` is traceable, via
its `ProvenanceRef`s, to META-3 source keys and pinned table versions (data-hygiene, META-1 Principle 7); and
any headline output is an ensemble, never a scalar (S3).

## 3. Method & governing definitions

### 3.1 Numeric substrate
- **f64 for all physical state and accumulation** (META-1 §3). f32 only for bulk storage/visualization
  with a doc-justified error bound. No implicit narrowing in numeric paths.
- **Integers** (`u32`/`u64`/`usize`) for counts, indices, cell/table addresses — never f64 for these.
- **No rational/fixed-point** in physics (rejected: breaks with irrational constants and transcendental
  functions; f64 + the determinism rules of META-1 §2 give reproducibility without it).
- **Finiteness is mandatory** (§3.8).

### 3.2 Units strategy
- **Internal unit system: strict SI base + coherent derived**, normalized to base units (m, kg, s, K,
  mol, A, and Pa/W/J/N/T/…). Angles in **radians**; temperature in **K**; energy in **J** internally.
- **Decision — units-typed at boundaries, `f64` in kernels:** physical quantities are **unit-typed at
  interfaces** (config parse, boundary-object declarations, table columns, results reporting, public
  function signatures) — where a wrong unit is catastrophic and the path is cold — and **plain `f64` /
  `ndarray<f64>` inside hot solver kernels**, which operate in documented SI base units. Rationale
  (researched 2026-07-14): a units crate (`uom`) is *zero-cost for f64* and normalizes to base SI, but
  typed `Quantity` arrays fight `ndarray`/`hdf5`/`rayon` in kernels; uom's own guidance is "units only
  at interface boundaries." `uom` is the starting choice; a hand-rolled newtype layer is the fallback
  if its compile-time/interop cost bites (a Layer-3 call).
- **Boundary conversion policy:** any non-SI quantity (a citation in psia, drum angle in degrees, Isp in
  s) is converted **at the edge**, and the conversion factor is recorded via META-3 (META-2 §2.1).
- **Reporting convention:** Isp in **seconds** (with g₀ = 9.80665 stated) *and* effective exhaust
  velocity `v_e` in m/s; `c` reserved for speed of light, `c*` for characteristic velocity.
- **Constants** come from the pinned CODATA set via META-3; a typed constants module exposes them
  (the latest official CODATA set as of project start, version-stamped in META-3). No re-typed constants per solver.

### 3.3 The medium-state vector `M`  *(the data-model form of Rule 12)*
Every constitutive/transport law (stopping, reactivity, EOS, opacity, conduction closure, …) takes its
local physical state **only** through a single `MediumState` bundle `M`. No law may branch on a material
ID or a regime label. `M` carries (fields present-or-null as the cell warrants):

- per-species number densities `n_s` and nuclear charge `Z_s`; mass density `ρ`;
- **free**-electron density `n_e^free`, **bound**-electron density `n_e^bound`, mean ionization `⟨Z⟩`;
- electron and ion temperatures `T_e`, `T_i` (and a single `T` when in equilibrium);
- degeneracy parameter `θ = T/E_F`;
- magnetic field vector `B`; bulk flow velocity `u`;
- optional derived cache (plasma frequency, Coulomb log, Debye length) — computed from the above, never
  an independent input.

"Cold solid," "warm dense matter," "hot magnetized plasma" are **corners of `M`**, not types. A law's
regime-dependence is expressed as continuous functions of `M` (e.g. bound↔free stopping via a continuous
weight), enforced by the **seam test** (META-1 Rule 12). FND-2 stores `M` per cell; FND-1 fixes its shape so
every solver reads the *same* `M`.

### 3.4 The three uncertainty representations
Because the inner run is deterministic and scalar (no distribution arithmetic), uncertainty lives at the
*ends*, not the middle:

**(a) `UncertainInput` — declaration.** How a single scalar input varies. Fields:
- `nominal` (Quantity), `units`;
- **either** a *parametric* law — family ∈ {normal, lognormal, uniform, triangular, delta} + parameters
  — **or** an *empirical sample-set* reference (a pinned set of pre-generated realizations, e.g. SANDY
  nuclear-data samples; carried as a pinned HDF5 dataset per FND-5, one row per realization, joint
  across correlated nuclides) (§3.5, §3.9);
- `correlation_group` (id) and/or membership in a parametric correlation matrix;
- `ProvenanceRef` (§3.6);
- `ValidityEnvelope` (§3.7).
A boundary object or table column that cannot supply a family-or-samples + provenance + envelope **cannot
be registered** (COUP-7 contract).

**(b) Runtime scalar — one realization.** For ensemble member `k`, each `UncertainInput` is resolved to a
single `f64` (in SI base units) by the sampling contract (§3.5). The solver sees only these numbers.
`delta` inputs resolve to their nominal (deterministic inputs are the trivial case — one code path, no
seam between "certain" and "uncertain").

**(c) `ResultDistribution` / `PBox` — output.** For each reported quantity: the **full vector of member
values** (retained by default — Sobol and arbitrary quantiles need it; config may thin for huge sweeps),
plus a **summary** {mean, std, p2.5/p16/p50/p84/p97.5}, plus a `ProvenanceRef` set and a pedigree handle
(COUP-6). **(v1.3) The member vector is the *aleatory* ensemble at one epistemic setting; epistemic
uncertainty — model-form band and numerical/interpolation error (FND-5 §3.4) — is carried as *intervals*
and propagated in a nested outer loop (COUP-5), so the reported result is a `p-box` (an interval-valued
CDF: the summary quantiles become intervals `[lo, hi]`), never a single collapsed distribution (META-1
§4.1; the false-confidence guard).** This p-box is the object FND-6 serializes and plots. [META-3: `pbox`,
`false-confidence`]

### 3.5 Sampling & determinism contract  *(defined here, executed by COUP-5)*
- **Method: Latin-hypercube over declared uncertainties** (COUP-5); FND-1 fixes the per-input mapping.
- **Parametric inputs** are sampled by **inverse-CDF** of the declared family applied to the member's
  stratified uniform draw (inverse-CDF is exact, monotone, and reproducible).
- **Empirical sample-sets** are sampled by **joint index**: member `k` selects realization `k mod N`
  (or a seeded permutation), so a SANDY set's *internal correlations are preserved automatically*.
- **Correlation** (sample-sets + parametric matrix): sample-sets carry correlation
  intrinsically; parametric inputs sharing a correlation matrix are correlated by **Iman–Conover rank
  correlation** over the LHC sample (rank-based, distribution-free, preserves the marginals) — pinned as
  the standard method here.
- **Deterministic RNG contract (META-1 §2):** every draw is produced by a **counter-based / splittable
  PRNG** keyed on `(master_seed, input_id, member_index, dimension)`. A draw depends *only* on those
  keys — never on thread, iteration order, or wall-clock. The RNG algorithm + master seed are recorded
  in the results bundle. This guarantees the same ensemble on 1 or N threads (Tier-1 determinism).

### 3.6 Provenance
- **`ProvenanceRef`** = {META-3 source `key`(s), table semantic-version + content hash (FND-5), and for a
  sampled value the RNG key of §3.5}. It is a *handle*, not heavy data.
- **Where provenance lives:** on `UncertainInput` declarations and on `ResultDistribution` outputs —
  **never carried on inner-loop `f64` scalars** (that would wreck the plain-f64 kernels and cost more than
  it's worth). The chain result→table→source is reconstructed from the input registry + the results
  bundle, not from per-scalar tags.
- Guarantee: a paper's reference list is the union of the `ProvenanceRef` source keys along a result's
  pedigree path (META-3 §1.1).

### 3.7 Validity envelopes
- **`ValidityEnvelope`** = the domain over which a table/correlation/input is declared valid (an
  axis-aligned box in its parameter space at minimum; richer regions allowed).
- **Policy hook:** interrogation outside the envelope triggers `{flag | refuse}` per config (default
  **refuse** for headline runs — the structural defense against speculation, META-1 Principle 8). FND-5 enforces
  for tables; COUP-5 records envelope hits per member so a partially-out-of-envelope ensemble is reported
  honestly, not silently clamped.

### 3.8 Finiteness & fail-loud
- Every physics `f64` must be **finite**. A NaN or ±Inf in a physics field is a **bug**, not a signal:
  the conservation/consistency audit (COUP-2) detects it and the run **halts with a diagnosis**
  (mechanism, location, member index) — never propagated (META-1 Principle 6). "Missing" is represented by an
  explicit optional/absent field, never by a sentinel NaN.

### 3.9 Distribution families (v1 set)
normal, lognormal (for positive-definite quantities like rates/cross-sections), uniform, triangular,
delta (deterministic), **empirical sample-set** (the general/nuclear-data carrier), and **interval**
(epistemic — a bound with *no* assumed distribution, propagated in the outer loop and yielding the
p-box of §3.4c). Extensible via a declared family enum; adding one is a doc amendment here, not
per-solver code.

## 4. Coupling relationships
- **FND-2** stores `Quantity`/`MediumState` per cell; FND-1 fixes their shape. FND-2 never redefines a field.
- **FND-5** table columns are typed `Quantity`s carrying `UncertainInput` (covariance-derived) + `ValidityEnvelope`; the loader enforces §3.7.
- **COUP-5** executes §3.4(a)→(b) sampling under the §3.5 contract and reduces (b)→(c); FND-1 is the contract it implements. Sobol consumes the retained member vectors of §3.4(c).
- **COUP-6** reads `ProvenanceRef` + per-module validation status → pedigree handle on each `ResultDistribution`.
- **COUP-7** boundary objects are the archetypal `UncertainInput` producers (citation+envelope+band = the three required fields).
- **FND-6** serializes `ResultDistribution` (members + summary + provenance) and enforces the regeneration contract (S6).

## 5. Uncertainty & validity (of this doc)
FND-1 introduces **no physics uncertainty of its own** — it is the framework others declare uncertainty
*in*. Its only "closure" choices are representational (inverse-CDF sampling, Iman–Conover correlation),
which are exact/standard, not model-form. Validation is therefore about *faithful representation*, not
physical accuracy (§6). Validation-ladder status: N/A (infrastructure), but its determinism and sampling
correctness are CI-gated (VAL-3).

## 6. Validation plan
Unit/property tests (test-first, VAL-3):
1. **Sampling determinism:** identical `{seed, declarations}` → byte-identical member draws at 1 vs N threads and across repeat runs (Tier-1).
2. **Marginal reproduction:** large-N ensembles reproduce each declared family's mean/variance/quantiles within sampling error.
3. **Correlation reproduction:** Iman–Conover induces the target rank-correlation (within tolerance) while preserving marginals; empirical sample-sets reproduce the source covariance.
4. **Inverse-CDF exactness:** known analytic quantiles recovered.
5. **Provenance round-trip:** a `ResultDistribution` reconstructs the correct META-3 keys + table versions.
6. **Envelope policy:** out-of-envelope interrogation flags/refuses per config; partial-ensemble envelope hits are reported.
7. **Finiteness guard:** an injected NaN halts with the correct diagnosis.

## 7. References
META-3 keys: `latin-hypercube` (LHC), `saltelli-sobol` (Sobol, downstream), `sandy-samples` (empirical
sample-sets), CODATA constants (§4 META-3). Iman & Conover (1982) rank-correlation — *to add to META-3
when COUP-5 is written* (`iman-conover`). Units crate: `uom` (docs.rs/uom), zero-cost dimensional
analysis, verified 2026-07-14.

*(No open questions — all engineering decisions resolved into the body per META-1 §6.)*

## 8. Change log
| Date | Version | Change |
|---|---|---|
| 2026-07-14 | 0.1 | Initial draft. Three-representation uncertainty model; medium-state vector `M` as the Rule-12 data form; deterministic sampling & provenance contracts. |
| 2026-07-14 | 0.2 | Resolved all six questions as engineering calls (units→types-at-boundaries+f64-kernels; ensemble-only; sample-sets+parametric correlation via Iman–Conover; CODATA-latest; `uom` start; HDF5 sample-sets) and integrated into the body; register removed per the new escalation criteria. |
| 2026-07-20 | 0.3 | **v1.3 p-box.** `ResultDistribution` extended to a **p-box** output (aleatory member vector at one epistemic setting + epistemic model-form/numerical *intervals* propagated in a nested outer loop → interval-valued CDF; never collapsed to one distribution, META-1 §4.1). Added the **interval** (epistemic) distribution family (§3.9). Updated the `M`-consumer interface row to the unified operators + constitutive spine (v1.3). (Header was stale at 0.1; corrected.) |
