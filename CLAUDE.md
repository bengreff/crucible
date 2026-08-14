# CLAUDE.md — Session-Start Context for CRUCIBLE

CRUCIBLE is an in-vacuum propulsion physics sandbox (chemical → antimatter): one unified 3-D
adaptive-dimension field simulator, every result a distribution + pedigree. **This codebase builds the
instrument; the research uses the instrument.** No engine-specific features, ever.

## State (2026-08-14)

- **Design phase complete. The coding gate is OPEN.** All 29 critical-path Layer-2 docs are
  **Reviewed (2026-08-14)** after a full review cycle (68 findings found, fixed, verified —
  `REVIEW_FINDINGS.md` is the record; `REVIEW_PREP.md` the process). No code exists yet.
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
