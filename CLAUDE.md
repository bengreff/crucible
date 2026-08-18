# CLAUDE.md — Session-Start Context for CRUCIBLE

CRUCIBLE is an in-vacuum propulsion physics sandbox (chemical → antimatter): one unified 3-D
adaptive-dimension field simulator, every result a distribution + pedigree. **This codebase builds
the instrument; the research uses the instrument.** No engine-specific features, ever.
`VISION_SCOPE.md` (v1.4) outranks everything; its §15 is the amendment log.

## Reading order (fresh session)

`VISION_SCOPE.md` → `docs/meta/META-0` (catalog + gate rules) → `docs/meta/META-1` (build doctrine;
**Rules 12/13 are supreme**: one law per phenomenon over the medium-state vector `M`, no
`if(material)`/`if(regime)` branch; configs are pure data; reactions are sources) →
`docs/meta/META-2` (conventions) → **the owning doc's §3 before coding any operator**.
`docs/meta/META-3` = source ledger, consulted per datum. One-owner rule: docs cross-reference,
never restate. `SESSION_LOG.md` holds the detailed per-session history (measured data, findings,
review waves) — consult it for the story behind a surface; this file carries only current state.

## State (2026-08-17)

Design complete: all 29 critical-path Layer-2 docs **Reviewed 2026-08-14** (`REVIEW_FINDINGS.md`).
Ten coding sessions, every one gates-green and committed; two multi-agent code reviews (sessions
6 and 10) with every confirmed finding fixed. `scripts/check.sh` = the 5-gate battery (fmt,
clippy, cargo test, offline pytest, certificate regen + diff). **104 Rust + 25 Python tests.**

**Goal A ✓** — conduction convergence certificate (`certificates/convergence_certificate.md`).

**Goal B — the BLIND RL10 (M2). Five certificate stations, strict ladder order, four done:**

| # | Physical system | Certificate | Status |
|---|---|---|---|
| 1 | Bursting diaphragm (Sod) vs exact Riemann | `certificates/station1_sod_certificate.md` | ✓ |
| 2 | De Laval nozzle vs isentropic theory; emergent p_c | `certificates/station2_nozzle_certificate.md` | ✓ |
| 3 | Flame seam: CEA → (p,h,Z) HDF5 vs RP-1311 | `certificates/station3_flame_certificate.md` | ✓ |
| 4 | Cooled wall: one wall law + conjugate liner | `certificates/station4_cooled_wall_certificate.md` | ✓ |
| 5 | **Blind RL10 assembly** (predict 73.4 kN / Isp ≈ 444 s, §9 blind rule) | — | **← NEXT** |

## What exists (by area — deferral owners live in each module's header)

- `crates/constants` — CODATA 2022, provenance-typed.
- `crates/units` — sole owner of pinned uom 0.38 (META-2 §4 ★): typed quantities at interface
  boundaries, `si()` extraction, documented-SI `f64` inside kernels.
- `crates/config` + `crates/registry` — FND-4 loader (no-hidden-defaults fixed point,
  resolved-config replay, dimensioned param accessors) + COUP-8 subset; 3 mechanism rows
  (`conduction`, `flow`, `wall_heat`) with validity-range refusals.
- `crates/tables` — FND-5 loader/interp on static libhdf5: pin/digest/envelope gates, multilinear
  in `interp_rule` space, `expect_units` bind gate. **`digest.rs` = digest v3, THE cross-language
  contract** (schema_version and exactly-one-sigma-form are pinned; golden vector asserted in both
  languages). Deferred kinds refuse loudly.
- `crates/grid` — FND-2 core: exact cylindrical metrics (`face_radius` single owner), Morton 8×8
  brick arena, SoA fields, per-brick N_θ with θ-coarsen/refine + symmetry controller, **ternary
  regions** (gas/solid/exterior) through the §3.6 ingest seam, `gas_solid_faces()` wall-face
  enumeration; surface sealed.
- `crates/solvers` — conduction (domain-selected, interface-aware, Robin faces; Goal-A certified);
  Euler = SOLV-1 §3.1–3.3 subset (PPM/HLLC-Batten on the exact metric, mask-aware sweeps,
  slip-ghost walls, whole-operator MMS at formal order); `wall_heat` = the one Colburn-class wall
  law (**±20–30% declared band**); station fixtures + SOLV-7 plane diagnostics. All explicit
  integrators are **honest scaffolding**, superseded by COUP-3's SDC-IMEX when it lands.
- `offline/` — `crucible_offl` (Python 3.13, exact pins incl. `cea==3.3.2`): digest-v3 mirror +
  h5py writer, NASA-CEA engine behind SI boundaries, Cantera cross-check, surface generators with
  measured interp-error bounds (×1.5 declared sampling margin; fresh-holdout CI gate). Production
  tables: `tables/chem/lox_lh2_v0.1.0.h5`; **`…pins.toml` sidecar is the single pin owner**
  (Rust tests parse it; the FND-4 §6-4 grammar will too).
- `data/anchors/` — **TM-107318 cached** (sha256 2d25422c…, META-3 `rl10-tm107318`): geometry
  Table 2.5.1 p.6 + App. E p.135; jacket Table 2.4.1 + App. D; 16-station cycle Table 6.1.1 p.37;
  pump/turbine maps App. B/C are **calibrated-mode only** (blind = coax-family η_c\* ±1–3% +
  pump-class envelopes, VAL-2 N18/D-G).

## Station 5 — the blind RL10 (2–3 sessions; wave (a) committable alone)

- **(a) EOS seam:** `BoundColumn` handle on crucible-tables FIRST (interpolate() re-parses rule
  strings + allocates per query; SOLV-1 makes ~10⁸ calls/run); then SOLV-1 §3.4 shifting mode —
  Z advects in the existing ρC slot, per-cell equilibrium projection at (p, h, Z) via fixed-count
  p-iteration (h = e + p/ρ), gamma-law retired to table data; FND-4 §6-4 `[tables]` pin grammar
  consuming the pins.toml sidecar. Add the cross-process table-regeneration determinism test.
- **(b) Assembly:** digitize the RL10 contour stations (Table 2.5.1/App. E → cached CSV, cited);
  `build_with_regions` from the contour (gas + liner + exterior); rebuild the coupler as ONE
  stepper parameterized by BC/ledger config, with per-face brick/cell indices + areas precomputed
  (retires station-4's duplicate stepper, BTreeMap-per-step, and inline grid indexing); COUP-7
  injector prior tier (premixed at MR; η_c\* band applied as SOLV-1 §3.4's **source-level**
  heat-release knockdown — output-side multiplication forbidden) + cooling-jacket coolant side;
  Bartz nozzle-envelope oracle rides here. Shared mask-aware, eos-threaded plane-diagnostics
  module (station-4's copies retire).
- **(c) Cycle + blind run:** COUP-3 §3.5 expander fixed point on ṁ (Aitken, fixed sweep count,
  pump-class envelopes); p_c emerges; SOLV-7 reads thrust/Isp/c\*/C_F vs the TM-107318 p-box
  (overlap metric d reported; every score labeled blind/calibrated). Certificate = the blind
  headline. Thrust from exit momentum flux, never assumed; uncertainty from declared closure
  bands (full COUP-5 UQ is a later wave).

## Cross-cutting deferrals (module headers carry the per-module lists)

- **COUP-3** SDC-IMEX class-D supersedes every explicit scaffolding integrator and brings
  Robin-Robin Picard/Aitken wall coupling.
- **COUP-7** boundary objects supersede `StagnationInflow`, caller-supplied BC closures, and the
  coolant Robin film.
- **FND-3** voxelization: partial apertures/cut cells retire stair walls + transpiration; brings
  mask-disjointness validation and aperture-aware interface classification.
- **OFFL-3** remaining products: frozen-path surface (SOLV-1 frozen mode), expansion oracles
  (station 5), B′ (SOLV-8), transport feed (OFFL-5, S23); per-point sigma columns (COUP-5).
- **FND-7 spine**: per-cell transport/material data; retires the sole-instance two-material
  ceiling and degenerate constant-transport occupants (gamma-law pattern).
- libm/powf note: certificate byte-diffs are valid on the pinned dev host only.

## Working rules & non-negotiables

- **Read the owning doc §3 before coding its operator.** Explain results to Ben as physical
  systems. Every session ends gates-green and committed.
- **Two-language rule:** runtime = 100% Rust; offline = Python; versioned HDF5 the only seam.
- **Test-first:** analytic/manufactured before benchmark before hardware (VAL-1); no red tests,
  no half-refactored solver at session end.
- **Determinism** (META-1 §2): fixed-order paths bit-exact; declared relaxed paths
  tolerance-bounded + non-spiraling; chaotic + relaxed = load refusal.
- **Fail loud, halt clean, never guess** (META-1 P6): out-of-envelope ⇒ refuse/flag, never clamp.
- **Every result is a distribution + pedigree** (S3); the p-box is never collapsed.
- **Firewall** (VISION_SCOPE §10): no MCNP, no restricted codes/data; pulsed-fission stops at
  published envelopes. No amendment path.
- A change contradicting a Reviewed doc ⇒ **amend the doc first** (change log; Layer-1 needs a
  VISION_SCOPE §15 entry). Never code around a doc.
- **Do not build:** SOLV-5, OFFL-4 (unreviewed); COUP-1 is a retired tombstone.
- **Ben's post-checkpoint goal:** data visualizations of CRUCIBLE results once the ladder yields
  real data — its own wave after the stations; keep FND-6 viz-friendly meanwhile.
