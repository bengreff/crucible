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
Twelve coding sessions, every one gates-green and committed; three multi-agent code reviews
(sessions 6, 10, 12) with every confirmed finding fixed. `scripts/check.sh` = the 5-gate battery
(fmt, clippy, cargo test, offline pytest, certificate regen + diff). **136 Rust + 28 Python
tests.** **Blind rule v1.4.1**: blind = mechanical input-blindness; every certificate declares
`development-observed: yes/no`; the RL10 campaign is declared **open development**.

**Goal A ✓** — conduction convergence certificate (`certificates/convergence_certificate.md`).

**Goal B — the BLIND RL10 (M2). Five certificate stations, ALL FIVE EARNED (session 12):**

| # | Physical system | Certificate | Status |
|---|---|---|---|
| 1 | Bursting diaphragm (Sod) vs exact Riemann | `certificates/station1_sod_certificate.md` | ✓ |
| 2 | De Laval nozzle vs isentropic theory; emergent p_c | `certificates/station2_nozzle_certificate.md` | ✓ |
| 3 | Flame seam: CEA → (p,h,Z) HDF5 vs RP-1311 | `certificates/station3_flame_certificate.md` | ✓ |
| 4 | Cooled wall: one wall law + conjugate liner | `certificates/station4_cooled_wall_certificate.md` | ✓ |
| 5 | **RL10 assembly vs the TM-107318 p-box** | `certificates/station5_rl10_certificate.md` | ✓ (coarse tier) |

**Station-5 headline (all scores `development-observed: yes`):** BLIND (coax-family η_c\* band ×
wall-law band corners): **F and Isp OVERLAP the record** (Ferson d = 0); p_c/c\*/C_F miss
coherently ~5.5–5.9% — one coarse-tier discretization signature (→ ~1.4% under the indicative
refinement). CALIBRATED (closed expander, TM component data): **F, Isp, AND the emergent p_c all
OVERLAP** (nominal 463.8 psia vs 475–482; the wall-law band sweeps delivered ṁ 16.74→18.56 kg/s
while Isp self-regulates flat at ~440 s — real expander behavior, reproduced not imposed).
**KNOWN LIMIT (owner named):** dials ≥ 8 cannot ESTABLISH by physical march (startup transients
genuinely leave the equilibrium surface; five schedules probed, all refuse loudly — see
`configs/rl10_full.toml`); the cure is **COUP-3 §3.6 pseudo-transient continuation** (the doc's
own default route to steady points) — the prime next-session candidate, unlocking the certified
refinement sweep. **Ben's dataviz checkpoint is REACHED**: the ladder now yields real data
(settled RL10 fields in `runs/*/fields.csv`, r/z/ρ/u/p/T/Z/M per cell + solid liner T).

## What exists (by area — deferral owners live in each module's header)

- `crates/constants` — CODATA 2022, provenance-typed.
- `crates/units` — sole owner of pinned uom 0.38 (META-2 §4 ★): typed quantities at interface
  boundaries, `si()` extraction, documented-SI `f64` inside kernels.
- `crates/config` + `crates/registry` — FND-4 loader (no-hidden-defaults fixed point,
  resolved-config replay, dimensioned param accessors) + COUP-8 subset; **§6-4 `[tables]` pin
  grammar live** (explicit pair or pins-sidecar via `load_str_with_sidecars`; resolved replays
  purely; manifest `table_pins`); **FND-3 contour grammar** (`contour` CSV content-addressed like
  a pin + the fidelity dial `cells_across_throat` → derived extents, manifest-recorded);
  `[operating_profile]` steady-march subset (flowthroughs/cfl/fill_p_pa/pumpdown/
  **p_amb_floor_pa** (declared altitude-cell floor)/**injector_ramp_flowthroughs** — session 12).
- `crates/tables` — FND-5 loader/interp on static libhdf5: pin/digest/envelope gates, multilinear
  in `interp_rule` space, `expect_units` bind gate; **`BoundColumn`** (bind-once units gate +
  rule parse + ln-hoist; allocation-free queries bit-identical to `interpolate()`). **`digest.rs`
  = digest v3, THE cross-language contract** (schema_version and exactly-one-sigma-form are
  pinned; golden vector asserted in both languages). Deferred kinds refuse loudly. Session 12:
  optional **`interp_error_bound_log`** (rule-space/relative bound for log-valued columns) —
  RECORDED DEFERRAL: rides outside digest v3; digest v4 folds it in.
- `crates/grid` — FND-2 core: exact cylindrical metrics (`face_radius` single owner), Morton 8×8
  brick arena, SoA fields, per-brick N_θ with θ-coarsen/refine + symmetry controller, **ternary
  regions** (gas/solid/exterior) through the §3.6 ingest seam, `gas_solid_faces()` wall-face
  enumeration; **FND-3 cut geometry (session 12): `build_with_geometry`** (per-cell κ + 4 face
  apertures, validated: Gas ⇔ κ>0, bitwise shared faces, covered-vs-κ=0) + **`wall_closure`** =
  THE discrete interface identity (well-balance-defined wall vector; |W| = smooth wall area).
- `crates/solvers` — conduction (domain-selected, interface-aware, Robin faces; Goal-A certified);
  Euler = SOLV-1 §3.1–3.4 + §3.6 (PPM/HLLC-Batten on the exact metric; **aperture-weighted
  sweeps + Berger–Giuliani State Redistribution** (κ < 0.5, conservation exact, slivers at the
  UNCUT CFL — session 12; full-box worlds bit-identical by arithmetic-identity defaults);
  **rayon-parallel by brick-row/column ownership partition — bit-exact at any thread count,
  asserted**; `EosLaw` seam, NPRIM=8 aux slots, datum-free Roe-averaged c² wavespeeds);
  **`TableEos`** (per-cell (p,h,Z) projection: warm-started + uniqueness-guarded fast path,
  rule-space slow-path acceptance vs the density column's own log bound; **S18 `h_offset`
  knockdown FIXED session 12** — store true energy, interrogate at h+δ; measured slope −0.847%
  c\* per −3e5 J/kg); `wall_heat` = the one Colburn-class law (**±20–30% band**). Explicit
  integrators = **honest scaffolding** until COUP-3's SDC-IMEX; **COUP-3 §3.6 pseudo-transient
  continuation = the named cure for fine-dial establishment (see Station 5)**.
- **`crates/engine` — the sandbox seam**: config → assembly (content-verified contour →
  FND-3 cut-geometry grid; refusals: cooling-with-zero-liner, closed-mode-never-engages,
  adiabatic liner holes) → the ONE coupled stepper (wall law on closure-vector patches with the
  SRD-neighborhood debit + liner conduction + coolant Robin + **the COUP-3 §3.5 closed-mode
  expander fixed point** — session 12: `turbopump_expander` boundary object, drive_power ←
  jacket heat_pickup, Aitken ≤ 1e-8, engages post-establishment) → SOLV-7 readout (+ measured
  inflow-plane ṁ honesty signal) + fields-CSV viz feed + **fault-tolerant crash artifact on
  halt** (`runs/<name>/crash_fields.csv`). Presets: `rl10_coarse.toml` (dial 5, ~60 s laptop,
  the certified tier), `rl10_calibrated.toml` (closed mode), `rl10_full.toml` (dial 16 —
  KNOWN LIMIT: awaits pseudo-transient). Certificate regen: `station5_rl10_certificate` bin
  (recorded readouts + Ferson rescoring, gate 5).
- `offline/` — `crucible_offl` (Python 3.13, exact pins incl. `cea==3.3.2`): digest-v3 mirror +
  h5py writer, NASA-CEA engine behind SI boundaries (**`gas_only` metastable mode**, deck-stamped,
  condensed-suffix filter), Cantera cross-check, surface generators with measured interp-error
  bounds (×1.5 margin; abs + **rule-space log bounds**; **envelope-EDGE holdout** — session-12
  review fix; fresh-holdout CI gates on BOTH pinned artifacts). Production tables:
  `lox_lh2_v0.1.0.h5` (station-3 certificate — untouched) + **`v0.3.2`** (station-5: gas-only
  metastable, Z narrowed to the premixed class MR ≈ 4.4–5.5, p ∈ [10 Pa, 7 MPa],
  h ∈ [−1.23e7, +3.8e6] with transient-sized ceiling); sidecars = single pin owners. Regen ≈
  2 min (`make_station5_tables.py`).
- `data/anchors/` — **TM-107318 cached** (sha256 2d25422c…, META-3 `rl10-tm107318`) + the
  **digitized geometry-of-record `rl10_contour.csv`** (Table E1 + Table 2.5.1 + Fig. E1 planes;
  closures declared in-header). **ERRATUM (session 11): Table 2.5.1's "Diameter" values are
  RADII** — r_throat = 2.47 in (c\* identity, ε=61 exit ≈ 40 in bell, Fig. E1 axis; VAL-2 0.2.3).
  Pump/turbine maps App. B/C stay **calibrated-mode only** (blind = coax-family η_c\* ±1–3% +
  pump-class envelopes, VAL-2 N18/D-G).

## Station 5 — DONE at the coarse tier (session 12); next steps

- All five ladder stations earned; the station-5 certificate scores blind + calibrated boxes vs
  the TM-107318 reference p-box (details in the certificate + SESSION_LOG session 12).
- **NEXT-SESSION PRIME CANDIDATE — COUP-3 §3.6 pseudo-transient continuation** (local-Δt SER):
  the named cure for the fine-dial establishment KNOWN LIMIT (see `configs/rl10_full.toml`),
  unlocking the certified refinement sweep (dial 8/12/16) that collapses the coarse-tier
  p_c/c\*/C_F miss (~5.9% → ~1.4% indicated). Grid-sequenced restart (FND-6) is the alternative.
- **Ben's dataviz checkpoint is REACHED** — real settled RL10 fields exist (`runs/*/fields.csv`);
  the dataviz wave is now unblocked as its own wave (Ben's call on priority vs pseudo-transient).
- Deferred, recorded: station-4 fixture rewire onto the engine stepper; Bartz nozzle-envelope
  oracle scoring; digest v4 (folds `interp_error_bound_log`); full COUP-5 UQ ensembles (the
  certificate's boxes are declared-band corner brackets); the multi-sided-wall (slot) interface
  class (FND-3 PLIC wave).

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
