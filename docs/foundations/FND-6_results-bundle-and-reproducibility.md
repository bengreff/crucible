# FND-6 — Results Bundle & Reproducibility

| Field | Value |
|---|---|
| **ID** | FND-6 |
| **Family** | FND (Foundations / spine) |
| **Status** | Reviewed (2026-08-14) — **skeleton** (per META-0 §5, W1) |
| **Depends on** | FND-1, FND-4, FND-5, COUP-6 |
| **Version** | 0.4 |
| **Skeleton/complete split** | **Fixed now:** the two-tier bundle layout (JSON summary + HDF5 payload); the `ResultDistribution` serialization; the provenance model (a PROV 3-verb subset, bespoke JSON); the **build-fingerprint** field list; the **S6 regeneration + verification contract**. **Deferred** (filled as dependents land, W2+): the exact per-result field schemas; plot/report specifics; the pedigree-report layout (COUP-6 owns the *score*; FND-6 serializes it); the thin web viewer (M12, if schedule allows). |

---

## 0. Purpose

FND-6 owns the **output bundle** and the **reproducibility contract S6**: *any result is regenerable from `{config file + pinned table versions + seed}` on a fixed build/target; any table is regenerable from its offline pipeline scripts.* It makes concrete the two headline promises — **every result is a distribution with a pedigree, never a bare number** (S3), and **the provenance chain result→table→source is fully walkable** so a campaign paper's reference list is assembled mechanically (META-3 §1.1).

This doc fixes *how a `ResultDistribution` (FND-1 §3.4c) is serialized*, *what the results bundle contains and how it self-describes*, *the build fingerprint that completes the run manifest*, and *how regeneration is verified*. Read after FND-1 (the result/provenance types), FND-4 (the manifest it embeds and completes), and META-1 §2 (the determinism mandate it operationalizes).

## 1. Scope & razor ruling

**Owns:** the results-bundle layout (light JSON tier + heavy HDF5 tier); `ResultDistribution` on-disk serialization; the **provenance model** (nodes + edges, and the source-key-union rule); the **build fingerprint** (the build/environment slice of the run manifest); the **regeneration & verification contract** (the tests that prove S6).

**Defers (owner named):**
- The **config-derived manifest** (resolved config, table pins, seed) → **FND-4** (FND-6 *embeds* it and *adds* the fingerprint; the complete manifest = FND-4-part + FND-6-part).
- The **`ResultDistribution` content/semantics** (member vector, summary set, `ProvenanceRef`) → **FND-1** (FND-6 owns only its *serialization*).
- The **pedigree score** computation and the backbone-map partition → **COUP-6** (FND-6 *serializes* the score + handle it produces).
- The **ensemble loop / Sobol / member reduction** → **COUP-5**; **table provenance metadata & hashing** → **FND-5**; the **CI harness** that runs the verification gates → **VAL-3**.

**Razor ruling:** pure infrastructure — no physics. Its fidelity obligation is that the bundle be *self-describing and complete enough* that a result survives loss of the reader code and can be regenerated years on from the manifest alone. It **carries** uncertainty and provenance; it originates neither.

## 2. Interfaces & contracts

| Interface | Consumed by | Contract |
|---|---|---|
| **Results bundle** (§3.1) | Ben, campaign papers, the viewer, COUP-6 | a directory/HDF5 with a JSON summary at the root indexing a self-describing HDF5 payload; every payload file referenced by relative path + content hash |
| **`ResultDistribution` serialization** (§3.2) | plots, COUP-6, re-analysis | member vectors (HDF5, chunked/thinnable) + summary {mean, std, p2.5/16/50/84/97.5} (JSON) + `ProvenanceRef` set + pedigree handle |
| **Provenance graph** (§3.3) | COUP-6, paper assembly | a typed node/edge graph whose transitive closure from a result to leaf source-keys **is** that result's reference list |
| **Complete run manifest** (§3.4–3.5) | regeneration, VAL-3 | FND-4 config part + build fingerprint = the exact S6 regeneration key, embedded in the bundle |

**Invariant promised to everyone:** any value in a bundle is traceable, via its `ProvenanceRef`s and the provenance graph, to META-3 source keys and pinned+hashed table versions; and the bundle carries enough to **regenerate the run bit-for-bit on a matching build** (or, cross-platform, within declared statistics — §3.6).

## 3. Method

### 3.1 Two-tier bundle *(light JSON index + heavy self-describing HDF5)*

The bundle splits by access pattern (MLflow's queryable-backend vs heavy-artifact split): **anything you'd grep, filter, diff, or cite across a campaign goes in JSON; anything you'd plot or resample goes in HDF5.** [META-3: `mlflow-run`]

**Light tier — `bundle.json`** (the walkable index and the campaign-paper source): `schema_version`; `result_id` (content hash); `run_hash` (hash of config + table hashes + build fingerprint); the **complete run manifest** (§3.4); `tables_used` (`[{semantic_version, content_hash, source_keys}]`); the **provenance graph** (§3.3) with a precomputed `source_keys_union` per result; the **pedigree score** (COUP-6); per output field a `summary` block; `rng` (`{algorithm, master_seed, keying_rule_id}` — O22); and `payload_ref` (relative path + content hash of the HDF5).

**Heavy tier — `payload.h5`**, laid out **NeXus-style self-describing**: every group/dataset carries reserved attributes (`@schema_version`, `@units`, `@semantic`), there is one mandatory root result entry, and axes/signals are tagged so the file is walkable **without FND-6's reader**. A mirror of the provenance graph is embedded so the payload is self-contained. [META-3: `nexus-hdf5`, `ro-crate`]

Sketch:
```
bundle.json                         # light index + manifest + provenance + summaries + pedigree
payload.h5
  /                                 @schema_version, @created (excluded from compare, §3.7)
  /result                          @result_id
    /config                        resolved-config bytes (or hash-linked to bundle.json)
    /ensemble                      @epistemic_setting_id = "nominal" (O23)
      /<field>                     dataset [member, …field_shape]; chunked on member axis;
                                   shuffle + zstd; @units @semantic @source_keys
                                   @conditional_on = "WORKS" where applicable (COUP-5 §3.2.1)
      /member_index                [N] int
      /rng                         {master_seed, algorithm, keying_rule_id} — every draw's key
                                   derivable (O22); no per-member key table
    /summary/<field>               [7]×[2] lo/hi intervals over
                                   {mean,std,p2.5,p16,p50,p84,p97.5} (the p-box, O22)
    /outer_summary/<field>         [n_outer]×[7]×[2] per-outer-evaluation summaries + setting ids
                                   (the envelope source, O23)
    /provenance                    serialized graph (mirror of bundle.json)
```

### 3.2 `ResultDistribution` serialization

FND-1 §3.4c fixes the *content*: the aleatory member vector, the summary set, the **epistemic intervals / p-box** (model-form + numerical), a `ProvenanceRef` set, and a pedigree handle. On disk:
- **Member vectors → HDF5**, one 2-D dataset per field `[member, …field_shape]`, **chunked along the member axis** (a chunk = a slab of members) so members append cheaply and per-member slices read cheaply; **shuffle + a fast codec (zstd/gzip)** per chunk; chunk size ~256 KB–2 MB. For huge sweeps, **thinning** keeps full member vectors for a decimated subset while **always** retaining the summary for *all* members. Any thinning is recorded, never silent (META-1 fail-loud spirit). [META-3: `hdf5-chunking`]
- **Summary + p-box → JSON** (the queryable/citable tier, §3.1). *(v1.3)* The member vector is the aleatory ensemble **stored at the NOMINAL epistemic setting (O23)** — `epistemic_setting_id` (= `nominal`) recorded on `/ensemble`; **per-outer-evaluation summary sets are retained** (`/outer_summary`, §3.1 — the source of the per-quantile envelope); **Sobol' columns are drawn at the nominal setting**. The epistemic **intervals** (model-form band; numerical/GCI error) and the resulting **interval-valued CDF** (lower/upper bound curves per field) are serialized alongside, so the reported bound is the **p-box** META-1 §4.1 requires — summary quantiles are stored as `[lo, hi]` intervals in **both tiers** (the JSON summary *and* the `[7]×[2]` HDF5 `/summary` datasets, O22), never collapsed to scalars, so the p-box survives the self-describing tier alone. [META-3: `pbox`, `false-confidence`]
- **RNG record → HDF5 (O22):** `{master_seed, algorithm, keying_rule_id}` — **not** a per-member key table. The versioned keying rule derives every draw's counter-based key from its coordinates `{input_id, member_index, dimension, outer_epistemic_index, purpose_tag}` (COUP-5 §3.1, a strict superset of FND-1 §3.5's base tuple), so **every key is derivable** and any draw replays from record + coordinates; storing N key tuples was redundant with the rule and could not hold within-member draws anyway. [META-3: `random123-philox`]

### 3.3 Provenance model *(a PROV 3-verb subset, bespoke JSON)*

Adopt the **semantics** of exactly three W3C PROV verbs — `used`, `wasGeneratedBy`, `wasDerivedFrom` — over typed nodes (a **result** Entity; an **inner run** Activity; **config/table/source** Entities; the **build** Agent), serialized as a small typed JSON graph, **not** RDF/OWL/SPARQL. This is the sweet spot the whole research-object community (RO-Crate) converged on. [META-3: `w3c-prov`, `ro-crate`]

The three verbs are the minimal set that makes the headline claim *mechanically* true: a result `wasGeneratedBy` an inner run that `used` config + table Entities; each table `wasDerivedFrom` its source-key Entities. Then **"a campaign paper's reference list = the union of source keys along a result's pedigree path" is exactly the transitive closure of `used`/`wasDerivedFrom` from the result down to leaf source Entities** (META-3 §1.1) — a graph-reachability query, precomputed and cached per result as `source_keys_union`. Every node is a **resolvable typed handle** (its semver + content hash *is* its identifier, FAIR-Digital-Object style), so the graph is portable to Zenodo/Software-Heritage later with no rework, and Agent/software nodes use CodeMeta vocabulary. [META-3: `fair-digital-object`, `codemeta`]

### 3.4 Build fingerprint *(the build/environment slice of the manifest)*

The exact fields that must match for bit-exact regeneration (completing FND-4's config-derived manifest, §3.5). **Software identity:** `rustc_version` (channel + version + commit); `llvm_version`; `target_triple`; **`target_cpu` + resolved `target_feature` set** (FMA-relevant — *mandatory*); `opt_level` + `codegen_units` (=1 for determinism); **FP policy** (contraction / `extra-fp-precision` off, fast-math-equivalents **off**, FMA-contraction policy, panic strategy); `libm_impl` + version (governs transcendentals); `Cargo.lock` content hash (full dependency closure) + workspace source commit; build profile. **Numerics/RNG:** `rng_algorithm` (e.g. `philox4x64`) + `master_seed`; reduction-ordering mode (fixed/tree — the thread-count-invariance guarantor). **Data seam:** for every table touched, `{semantic_version, content_hash}` (shared with FND-4/FND-5).

Rationale, load-bearing for the bit-exact claim: Rust follows IEEE-754 and by default **does not contract `a*b+c` into an FMA** and **does not reassociate** — which is *why* per-build bit-exactness is achievable — but this is contingent on `target-feature`/`target-cpu` (enabling FMA features changes emitted fused ops) and on **never inspecting NaN sign/payload bits** (non-deterministic across arch/flags). Hence target-cpu/features and libm are mandatory fingerprint fields, and cross-*target* bit-exactness is not promised (META-1 §2.3). [META-3: `rust-fp-rfc3514`, `flit-fp-determinism`, `reproducible-builds`, `guix-repro`]

### 3.5 The complete run manifest

`complete manifest = FND-4 config part (§3.6 there) + the §3.4 build fingerprint`, embedded verbatim in `bundle.json`. **You regenerate from this manifest against a matching build — never from the author config alone.** A pinned environment does *not* by itself guarantee bit-identical output (the Nix/Guix at-scale lesson: only a *deterministic* toolchain **plus** a *tested* byte-equality gate does), which is why §3.6 exists. [META-3: `guix-repro`]

### 3.6 Regeneration & verification contract *(the tests that prove S6)*

Payloads are compared **after stripping the excluded-from-compare zone** (§3.7). The gates (run by VAL-3):
1. **Determinism gate** (per build): run a result twice on the same build + same thread count → payloads on **fixed-order paths — which includes every chaotic-regime run** — **byte-identical** after canonicalization; payloads on **declared relaxed-reduction paths** (GPU throughput, non-chaotic regimes only) agree **within the declared negligible tolerance**, pass the **non-spiraling check**, and show **no verdict divergence** (the S6 regime→guarantee mapping, META-1 §2.1/§2.4). [META-3: `ect-consistency`]
2. **Thread-count-invariance gate:** run at thread counts {1, 4, N} → **byte-identical on fixed-order paths** (mandatory in chaotic regimes precisely so thread count can't seed a butterfly divergence — the reductions are deterministic exactly where the physics is sensitive); **tolerance-bounded + non-spiraling** on declared relaxed-reduction paths. This catches **non-deterministic parallel reductions**, the thing most likely to break reproducibility. Cross-*build*/cross-*platform* comparison relaxes to ensemble-consistency (gate 4). [META-3: `gamer2-determinism`, `ect-consistency`]
3. **Full-regeneration gate:** from `{config + pinned table versions + master_seed}` alone, on a matching build, regenerate and compare against the archived payload — **byte-compare on fixed-order paths; tolerance-compare (against the recorded bound) on relaxed paths** — exercises table pinning + RNG-key replay end-to-end.
4. **Cross-platform gate** (statistical, not bitwise): on a different target, assert summary stats agree within ensemble sampling error (`|Δmean| ≤ k·std/√N`) — matching the reality (earth-system-model and GAMER experience) that cross-platform bit-identity is unattainable; hence the Tier-1/Tier-2 split (META-1 §2.1). [META-3: `esm-bitwise-repro`, `gamer2-determinism`]
5. **Provenance-closure gate:** the cached `source_keys_union` equals a freshly computed transitive closure over the graph (guards a stale reference list).
6. **Failure diagnostics:** on any byte-mismatch, emit an **HDF5-aware structural diff** (group/dataset granularity, diffoscope-style) to localize the divergent field — raw byte-equality stays the pass/fail gate, but the *report* is structural. [META-3: `reproducible-builds`]

### 3.7 Determinism & the excluded-from-compare zone

Recorded but **ignored in byte-equality**: wall-clock timestamps, hostname, username, absolute paths, thread count, elapsed times (the `SOURCE_DATE_EPOCH` discipline). The HDF5 payload is **canonicalized** before hashing: HDF5 internal object timestamps disabled, dataset creation order fixed, volatile attributes zeroed. Everything that affects a *number* is in the fingerprint (§3.4) and thus *inside* the compare; everything that is mere metadata is outside it. [META-3: `reproducible-builds`]

## 4. Coupling relationships

- **FND-4** produces the config-derived manifest; FND-6 embeds it and adds the build fingerprint (§3.5) — the single shared seam (FND-4 owns config/seed/table-pins; FND-6 owns fingerprint + bundling).
- **FND-1** defines `ResultDistribution` and `ProvenanceRef`; FND-6 serializes them (§3.2–3.3).
- **COUP-5** produces the ensemble (member reduction → `ResultDistribution`) FND-6 writes; **COUP-6** produces the pedigree score + backbone-map partition FND-6 serializes.
- **FND-5** supplies table provenance + content hashes that populate `tables_used` and the provenance leaves.
- **VAL-3** runs the §3.6 gates as standing CI (a determinism regression is build-breaking, ranked with correctness — META-1 §2.4).

## 5. Uncertainty & validity

FND-6 introduces **no physics uncertainty** — it *carries* the ensemble and its bands, and *records* the interpolation/model-form/nuclear-data uncertainties originated upstream (FND-1/FND-5/COUP-5). Its own correctness (self-describing layout, sound provenance closure, exact regeneration) is **CI-gated** (VAL-3). Validation-ladder status: N/A (infrastructure).

## 6. Validation plan

Handled by the §3.6 gates, run as VAL-3 CI: (1) determinism, (2) thread-count invariance, (3) full regeneration from manifest, (4) cross-platform statistical agreement, (5) provenance-closure equality, (6) structural-diff failure reporting. Plus: (7) **bundle self-description** — the HDF5 payload is readable and walkable by a generic HDF5/NeXus reader with no FND-6 code; (8) **thinning honesty** — a thinned sweep records the thinning and still carries all-member summaries.

## 7. References

META-3 keys: `mlflow-run`, `nexus-hdf5`, `ro-crate`, `w3c-prov`, `fair-digital-object`, `codemeta`, `hdf5-chunking`, `random123-philox`, `rust-fp-rfc3514`, `flit-fp-determinism`, `reproducible-builds`, `guix-repro`, `esm-bitwise-repro`, `gamer2-determinism`. Depends on FND-1 (`ResultDistribution`, `ProvenanceRef`, RNG key), FND-4 (config manifest), FND-5 (table provenance/hashes), COUP-6 (pedigree).

*(No open questions — all engineering decisions resolved into the body per META-1 §6.)*

## 8. Change log

| Date | Version | Change |
|---|---|---|
| 2026-07-19 | 0.1 | Initial **skeleton**. Fixed the two-tier bundle (light JSON index/campaign-source + heavy NeXus-style self-describing HDF5); `ResultDistribution` serialization (member vectors chunked on the member axis with shuffle+zstd, thinnable-with-record, `rng_keys` for replay; summaries in JSON); the provenance model (a W3C-PROV 3-verb subset — `used`/`wasGeneratedBy`/`wasDerivedFrom` — as bespoke typed JSON, with `source_keys_union` = transitive closure = the paper reference list; resolvable typed handles); the build-fingerprint field list (rustc/llvm/target-triple/target-cpu+features/FP-FMA policy/libm/lockfile hash/rng algo+seed/reduction-order) completing the FND-4 manifest; the S6 regeneration + verification contract (determinism, thread-count-invariance, full-regeneration, cross-platform-statistical, provenance-closure gates; excluded-from-compare zone + canonicalized HDF5; HDF5-aware structural diff on failure). Per-result schemas, plot/report layout, and viewer deferred. |
| 2026-08-14 | 0.4 | **Review fix wave (O22, O23).** §3.1/§3.2: the HDF5 `/summary` datasets become **interval-valued `[7]×[2]` lo/hi** so the p-box survives the self-describing tier alone; the `/rng_keys [N]` dataset **replaced by `{master_seed, algorithm, keying_rule_id}`** (every draw's key derivable from the versioned keying rule + coordinates; a key table was redundant and couldn't hold within-member draws) (O22). §3.1/§3.2: member vectors stored at the **NOMINAL epistemic setting** with `epistemic_setting_id` recorded; **`/outer_summary` `[n_outer]×[7]×[2]`** per-outer-evaluation summaries retained as the envelope source; **Sobol' drawn at nominal** — stated in §3.2 (O23). Light-tier `rng` record gains `keying_rule_id`; `/ensemble` fields carry `@conditional_on = "WORKS"` where applicable (COUP-5 §3.2.1). |
| 2026-08-13 | 0.3 | **Consistency sweep.** §3.6 gates 1–3 aligned to the S6 regime→guarantee mapping (fixed-order paths — including *all* chaotic-regime runs — byte-identical; declared relaxed-reduction paths, non-chaotic only, tolerance + non-spiraling; ECT = cross-platform only): the previous gate text had inverted the chaotic/fixed-order assignment (same fix as META-1 §2.1/§2.4, VISION_SCOPE §8). `rng_keys` serialization (§3.1 sketch, §3.2) extended to COUP-5 §3.1's full key set (`outer_epistemic_index`, `purpose_tag`) so double-loop/Sobol'/bootstrap draws are replayable. |
| 2026-07-20 | 0.2 | **v1.3 consistency-review fixes.** §3.2 now serializes the **p-box** (aleatory members + epistemic model-form/numerical intervals + interval-valued CDF; JSON quantiles as `[lo,hi]`), matching FND-1 v0.3 / META-1 §4.1 (was a bare distribution — the review's top coherence gap). §3.6 determinism gates gained the **chaotic-regime carve-out**: byte-identical/ECT per the v1.3 relaxed contract (fixed-order reductions retained in chaotic regimes so thread count can't butterfly). |
