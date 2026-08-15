# CLAUDE.md — Session-Start Context for CRUCIBLE

CRUCIBLE is an in-vacuum propulsion physics sandbox (chemical → antimatter): one unified 3-D
adaptive-dimension field simulator, every result a distribution + pedigree. **This codebase builds the
instrument; the research uses the instrument.** No engine-specific features, ever.

## State (2026-08-14)

- **Design phase complete. The coding gate is OPEN.** All 29 critical-path Layer-2 docs are
  **Reviewed (2026-08-14)** after a full review cycle (68 findings found, fixed, verified —
  `REVIEW_FINDINGS.md` is the record; `REVIEW_PREP.md` the process).
- **Coding started 2026-08-14.** Session 1: Rust 1.93.1 pinned, workspace scaffolded (crate per doc
  area under `crates/`), `crucible-constants` (CODATA 2022 stamped in META-3 §4), `scripts/check.sh`
  = the VAL-3 fast gate battery. Session 2: **FND-4 loader done** (`crucible-config`: §3.4 pipeline,
  §3.5 no-hidden-defaults fixed point, §3.6 manifest, O20/O21/θ-ladder checks; §6 tests 1–3, 5–9
  green) + `crucible-registry` (COUP-8 §3.1–§3.4 loader-facing subset). Known deferrals are listed
  in `crates/config/src/lib.rs` (FND-5 pin resolution + §6-4; deferred grammars refuse non-empty;
  span-annotated diagnostics pending — paths only).
- Session 3 (2026-08-15): **FND-5 table loader done** (`crucible-tables`: §3.1 schema on statically
  pinned libhdf5 2.2.0, §3.2 pin/digest/provenance gates, §3.3 multilinear in `interp_rule` space,
  §3.5 Refuse/Flag envelope policy, §6 tests 1–6 green + golden digest vector). The canonical
  content-digest algorithm in `crates/tables/src/digest.rs` is the **cross-language contract** —
  the Python OFFL side must mirror it exactly. Deferred with loud load refusals: `thermo_potential`
  (S11) + `sample_set` kinds, `pchip` method, config→tables pin wiring (FND-4 §6-4).
- Session 4 (2026-08-15): **FND-2 grid core done** (`crucible-grid`: cylindrical index space with
  exact ring metrics + zero-area axis faces + θ↔θ+π parity pairing, Morton-ordered 8×8 brick arena
  with SoA fields and per-brick N_θ, conservative θ-coarsen/refine with thermalized-ΔKE ledger,
  symmetry indicator + guard/hysteresis/dwell controller, fixed-shape tree reductions; §6 items
  3/7/8/9 green). Deferred with owners (in `crates/grid/src/lib.rs` header): §3.3 full cell model
  (FND-1/FND-7 wave), §3.5 tiles, §3.6 FND-3 ingest, §3.3(7) PLIC fields.
- Session 5 (2026-08-15): **GOAL A COMPLETE — the CONVERGENCE CERTIFICATE is earned and committed**
  (`certificates/convergence_certificate.md`, regenerable via
  `cargo run --bin convergence_certificate`; criteria CI-enforced in
  `crates/solvers/tests/goal_a_certificate.rs`). `crucible-solvers`: flux-form conduction on the
  exact cylindrical metric (one law, no branches; axis + N_θ=1 handled by geometry, not code),
  first real registry mechanism (`conduction`), config→grid `from_loaded` wiring, `[geometry]`
  extents + `axisymmetric` assertion grammar in FND-4. **Data: MMS orders 2.001/2.000 (2-D) and
  2.006/2.002 (3-D m=2 θ-mode); annulus 4.3e-4 rel; Bessel cylinder 6.3e-3 K on 100 K (crosses
  r=0); conservation drift 1.1e-16; byte-identical rerun.** Honest scaffolding note: explicit
  fixed-order reference integrator drives the certificate; superseded (not extended) by COUP-3's
  SDC-IMEX class-D implicit path when it lands. Uniform-N_θ sweeps only (refluxing = COUP-2/3 wave).
- **Next goal (Ben to confirm scope): the SOD CERTIFICATE** — SOLV-1 quasi-1D HLLC/PPM at N_θ=1 vs
  the exact Riemann solution, L1 convergence at formal order + conservation to round-off. Then
  OFFL-3 chemistry tables (the seam certificate) → blind RL10 (M2).
- Not yet reviewed (do **not** implement without a review pass first): SOLV-5, OFFL-4 (deferred
  pulsed/antimatter set). COUP-1 is a retired tombstone — never build it.
- `VISION_SCOPE.md` is at **v1.4** and outranks everything (§15 = amendment log).

## Reading order (fresh session)

`VISION_SCOPE.md` → `docs/meta/META-0` (catalog + gate rules) → `docs/meta/META-1` (build doctrine;
**Rules 12/13 are supreme**: one law per phenomenon over the medium-state vector `M`, no
`if(material)`/`if(regime)` branch; configs are pure data; reactions are sources) → `docs/meta/META-2`
(conventions) → the specific doc(s) for what you're building. `docs/meta/META-3` = source ledger,
consulted per datum. One-owner rule: docs cross-reference, never restate — trust the owner doc.

## The coding wave (next work — VISION_SCOPE §12, months 1–2 skeleton)

Target: the **blind RL10 chemical anchor** (S1) on the unified solver's N_θ=1 corner. Build order:

1. **Rust workspace + config schema** (FND-4; `[engine]` explicit bindings, determinism block) and the
   **table loader** (FND-5; scalar returns only — uncertainty rides the outer loop).
2. **Cylindrical world grid** (FND-2 v0.6): sparse bricks on (i_r, i_θ, i_z), adaptive N_θ with guard = 4,
   SoA layout. **CPU fixed-order reference first** — the GPU port comes later as a declared
   relaxed-reduction path validated against it (D-H).
3. **Conduction + the one wall-function heat law** (SOLV-1 §3.5 / COUP-2 §3.5 / COUP-3 class-`D`).
4. **N_θ=1 reacting flow** (SOLV-1: PPM/HLLC-Batten, elemental-composition advection, emergent p_c).
5. **Chemistry table pipeline** (OFFL-3: CEA + Cantera, (p, h, Z) local-state surfaces).
6. RL10 blind run per the §9 blind-mode rule: **no quantity measured on the engine under test** —
   design spec + universal closures + technology-class data from other hardware only.

First-session chores (META-3 flags these `RESEARCH-PENDING`): pin the Rust toolchain version, stamp
CODATA constants at first use, re-benchmark OpenMC/solver throughput on the actual desktop.

## Non-negotiables while coding

- **Two-language rule:** runtime = 100% Rust; offline pipelines = Python; versioned HDF5 tables are the
  only seam (VISION_SCOPE §6). No third language.
- **Test-first, tests-green-per-session:** solvers are built against analytic/manufactured solutions
  before benchmarks before hardware (VAL-1/VAL-3); no session ends with red tests or a half-refactored
  solver. Every session ends committed.
- **Determinism** per META-1 §2 (S6): fixed-order paths bit-exact (mandatory in chaotic regimes);
  declared relaxed paths tolerance-bounded + non-spiraling; chaotic + relaxed = load refusal.
- **Fail loud, halt clean, never guess** (META-1 P6): out-of-envelope ⇒ refuse/flag, never clamp.
- **Every result is a distribution + pedigree** (S3); the p-box is never collapsed.
- **Firewall** (VISION_SCOPE §10): no MCNP, no restricted codes/data, pulsed-fission stops at
  published envelopes. No amendment path.
- A design change that contradicts a Reviewed doc ⇒ amend the doc first (its change log; Layer-1
  boundaries need a VISION_SCOPE §15 entry). Never code around a doc.
