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

## State (2026-08-18)

Design complete: all 29 critical-path Layer-2 docs **Reviewed 2026-08-14** (`REVIEW_FINDINGS.md`).
Eleven coding sessions, every one gates-green and committed; two multi-agent code reviews (sessions
6 and 10) with every confirmed finding fixed. `scripts/check.sh` = the 5-gate battery (fmt,
clippy, cargo test, offline pytest, certificate regen + diff). **122 Rust + 26 Python tests.**
**Blind rule v1.4.1** (Ben, session 11): blind = mechanical input-blindness; every certificate
declares `development-observed: yes/no`; the RL10 campaign is declared **open development**.

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
  resolved-config replay, dimensioned param accessors) + COUP-8 subset; **§6-4 `[tables]` pin
  grammar live** (explicit pair or pins-sidecar via `load_str_with_sidecars`; resolved replays
  purely; manifest `table_pins`); **FND-3 contour grammar** (`contour` CSV content-addressed like
  a pin + the fidelity dial `cells_across_throat` → derived extents, manifest-recorded);
  `[operating_profile]` steady-march subset (flowthroughs/cfl/fill_p_pa/pumpdown).
- `crates/tables` — FND-5 loader/interp on static libhdf5: pin/digest/envelope gates, multilinear
  in `interp_rule` space, `expect_units` bind gate; **`BoundColumn`** (bind-once units gate +
  rule parse + ln-hoist; allocation-free queries bit-identical to `interpolate()`). **`digest.rs`
  = digest v3, THE cross-language contract** (schema_version and exactly-one-sigma-form are
  pinned; golden vector asserted in both languages). Deferred kinds refuse loudly.
- `crates/grid` — FND-2 core: exact cylindrical metrics (`face_radius` single owner), Morton 8×8
  brick arena, SoA fields, per-brick N_θ with θ-coarsen/refine + symmetry controller, **ternary
  regions** (gas/solid/exterior) through the §3.6 ingest seam, `gas_solid_faces()` wall-face
  enumeration; surface sealed.
- `crates/solvers` — conduction (domain-selected, interface-aware, Robin faces; Goal-A certified);
  Euler = SOLV-1 §3.1–3.4 (PPM/HLLC-Batten on the exact metric, mask-aware sweeps, slip-ghost
  walls w/ per-op `slip_wall_z_faces` policy, MMS at formal order; **`EosLaw` seam, NPRIM=8 aux
  (e, Γ₁) slots**); **`TableEos` = shifting-equilibrium mode** (combustion in the EOS: per-cell
  (p,h,Z) projection, Illinois + declared-bound slow path; η_c\* `h_offset` knockdown hook,
  S18; mass-flow injector inflow w/ sonic startup cap; scheduled `PressureOutflow`); `wall_heat`
  = the one Colburn-class law (**±20–30% band**). Explicit integrators = **honest scaffolding**
  until COUP-3's SDC-IMEX.
- **`crates/engine` — the sandbox seam (new, session 11)**: config → assembly (content-verified
  contour → ternary grid) → the ONE coupled stepper (wall law + liner conduction + coolant
  Robin, flux-matched) → SOLV-7 readout (N11 p_c, exit thrust integral) + fields-CSV viz feed.
  COUP-7 subset rows: `flow_shifting`, `injector_prior`, `jacket_coolant` (coolant side only,
  D-C). **An engine is pure data**: `crucible run configs/rl10_coarse.toml` (34 s laptop, dial 5)
  / `rl10_full.toml` (dial 16, desktop). O20 bindings exercised for real.
- `offline/` — `crucible_offl` (Python 3.13, exact pins incl. `cea==3.3.2`): digest-v3 mirror +
  h5py writer, NASA-CEA engine behind SI boundaries, Cantera cross-check, surface generators with
  measured interp-error bounds (×1.5 declared sampling margin; fresh-holdout CI gate). Production
  tables: `lox_lh2_v0.1.0.h5` (station-3 certificate) + **`v0.2.0` (station-5 envelope: p floor
  10 Pa for the vacuum-plume fringe)**; **`…pins.toml` sidecars are the single pin owners** (Rust
  tests + the §6-4 grammar parse them). Regen ≈ 5 s (`make_station5_tables.py`).
- `data/anchors/` — **TM-107318 cached** (sha256 2d25422c…, META-3 `rl10-tm107318`) + the
  **digitized geometry-of-record `rl10_contour.csv`** (Table E1 + Table 2.5.1 + Fig. E1 planes;
  closures declared in-header). **ERRATUM (session 11): Table 2.5.1's "Diameter" values are
  RADII** — r_throat = 2.47 in (c\* identity, ε=61 exit ≈ 40 in bell, Fig. E1 axis; VAL-2 0.2.3).
  Pump/turbine maps App. B/C stay **calibrated-mode only** (blind = coax-family η_c\* ±1–3% +
  pump-class envelopes, VAL-2 N18/D-G).

## Station 5 — the blind RL10 (waves (a)+(b) DONE, session 11)

- **(a) EOS seam ✓** and **(b) assembly ✓** — see SESSION_LOG session 11. The coarse preset
  (`configs/rl10_coarse.toml`, fidelity dial = 5) runs the establishment phase in **34 s on the
  laptop**; the dial is continuous (full preset = 16). Deferred out of (b), recorded: station-4
  fixture rewire onto the engine stepper (its certificate stays on its own stepper until then);
  Bartz oracle scoring; the per-face precompute perf item.
- **(b′) FIRST ITEM NEXT SESSION — the stair-corner blocker:** past ~2.9 flow-throughs an
  exit-lip corner cell starves (ρ runaway → CFL collapse; slip z-faces transpire, mirror z-faces
  shock — both are the small-cut-cell class). Designated cure: **FND-3 partial apertures + State
  Redistribution** (SOLV-1 §3.6 names it). Instrument first: dump fields at halt (crash-artifact
  writer), inspect the corner, then implement. Until fixed, no steady coarse readout and no
  settle for the full tier.
- **(c) Cycle + blind-config run:** COUP-3 §3.5 expander fixed point on ṁ (Aitken, fixed sweep
  count, pump-class envelopes); η_c\* knockdown calibration (the `h_offset` hook exists); p_c
  emerges; SOLV-7 vs the TM-107318 p-box (overlap metric d; every score labeled blind/calibrated
  **+ `development-observed: yes`**). Certificate = the headline. Thrust from exit momentum flux,
  never assumed; uncertainty from declared closure bands (full COUP-5 UQ later).

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
