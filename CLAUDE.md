# CLAUDE.md — Session-Start Context for CRUCIBLE

CRUCIBLE is an in-vacuum propulsion physics sandbox (chemical → antimatter): one unified 3-D
adaptive-dimension field simulator, every result a distribution + pedigree. **This codebase builds
the instrument; the research uses the instrument.** No engine-specific features, ever.
`VISION_SCOPE.md` (v1.5) outranks everything; its §15 is the amendment log.

## Reading order (fresh session)

`VISION_SCOPE.md` → **`PLAN_CHEMICAL_SANDBOX.md` (THE PLAN OF RECORD — session-by-session build
plan S1–S20; read §1 rulings + the current session's §5 entry before anything else)** →
`docs/meta/META-0` (catalog + gate rules) → `docs/meta/META-1` (build doctrine;
**Rules 12/13 are supreme**: one law per phenomenon over the medium-state vector `M`, no
`if(material)`/`if(regime)` branch; configs are pure data; reactions are sources) →
`docs/meta/META-2` (conventions) → **the owning doc's §3 before coding any operator**.
`docs/meta/META-3` = source ledger, consulted per datum. One-owner rule: docs cross-reference,
never restate. `SESSION_LOG.md` holds the detailed per-session history (measured data, findings,
review waves) — consult it for the story behind a surface; this file carries only current state.

## State (2026-08-19 — session 15 = plan S3 complete)

Design complete: all 29 critical-path Layer-2 docs **Reviewed 2026-08-14** (the review register was
closed 68/68-discharged and deleted 2026-08-19; findings live in doc change logs + git history).
Fifteen sessions, every one gates-green and committed; three multi-agent code reviews
(sessions 6, 10, 12) with every confirmed finding fixed. `scripts/check.sh` = the 5-gate battery
(fmt, clippy, cargo test, offline pytest, certificate regen + diff). **154 Rust + 28 Python
tests.** **Blind rule v1.4.1**: blind = mechanical input-blindness; every certificate declares
`development-observed: yes/no`; the RL10 campaign is declared **open development**.

**Session 14 = plan S2: THE REAL INTEGRATOR.** `crucible-solvers::sdc` = the ONE deterministic
SDC-IMEX step (COUP-3 §3.1: IMEX-Euler predictor + 2 fixed trapezoid correction sweeps; explicit
hyperbolic class A = `Euler::eval_rhs`; class D = spatially-coupled implicit diffusion via a
deterministic fixed-structure Jacobi-CG solve — matrix-free on the conduction assembly, δ-form
warm start, `EPS_CG_RESID`/`N_CG_ITERS_MAX` acceptance → `COUPLING_RESIDUAL`); **Robin-Robin
gas–wall exchange inside each sweep's class-D solve** (COUP-2 §3.5: wall-function h as the Robin
coefficient — linear in T_solid, unconditionally stable at Biot > 1; fixed Picard sweeps +
clamped Aitken; the gas debits exactly the per-face heats the accepted solid solve received);
**the COUP-2 §3.1 conservation audit armed EVERY step** (flux-telescoping port ledger in the
sweeps + applied-increment source ledger, `TOL_AUDIT[q] = K_AUDIT·ε·√N·S[q]` per §3.1.1,
violation = halt with diagnosis — closed at ≤ round-off through every certificate march incl.
15.8k-step RL10 members). Every scaffolding integrator retired (session-5 explicit conduction,
session-7 SSP-RK2, session-10/11 flux-matched coupled splitting); Δt = the gas CFL alone — the
solid/exchange stability guards are gone by construction (the S4 physical-ρc_p prerequisite).
Goal-A anchors now march 4–32× past the explicit bound (a new test holds at 512×); temporal
order 2 verified by dt-Richardson (flow 2.0; linear diffusion superconverges ~3). Certified
numbers moved as expected and were regenerated (headline: Sod star-plateau u* 7.3e-5 → 8.8e-4 —
a dissipation-profile shift at the captured shock, global L1/orders unchanged; station-5 boxes
shifted < 0.1%, conclusions identical). COUP-3 0.4.1 clarification amendment (CG = the §3.7
fixed-tolerance/fixed-cap rule). Deferrals: TOL_AUDIT constants manifest-recording rides FND-6;
mount-reaction vs SOLV-7 thrust cross-check (COUP-2 §6-7) with the verdict wave; class-D
assembly is serial (perf; GPU wave brings multigrid/parallel).

**VISION_SCOPE v1.5 (Ben, 2026-08-19) — the session-13 rulings, all doc-amended:** accelerated
convergence (pseudo-transient/local-Δt) is **DELETED** — every certified result is a **physical
march incl. start-up** (COUP-3 §3.6 tombstone); external-system timelines may be **compressed as
declared schedules** (COUP-7 §3.2.2), internal physics never; ignition = the **burn-progress
field c** (SOLV-4 §3.6: unburnt↔equilibrium blend, flame-speed + induction closures, igniter =
energy-deposit object, `NEVER_IGNITED`/`FLAMEOUT` halts); **full 3-D is the product tier**;
GPU (RTX 4080, 24-h cap) on the critical path with **bit-exact-per-device determinism**
(META-1 §2.5); no solids; no viz until THE RUN (plan S19).

**Session 15 = plan S3: THE MISSING FORCES.** `crucible-solvers::gas_diffusion` = SOLV-1 §3.1's
`F_visc` — compressible viscous stress + Fourier conduction + species diffusion on the exact
cylindrical metric (swirl included), the **gas occupant of COUP-3 class D** (0.4.2 amendment
landed with the code): per-component symmetric fixed-structure CG cores (u_r; **u_θ solved as
ω = u_θ/r in the angular-momentum form** — rigid rotation discretely stress-free, angular
momentum telescopes; **T in total-energy flux form** — dissipation from the KE bookkeeping;
C constant-ρD Fickian), cross-stress couplings converged by the same fixed Picard sweeps as the
wall exchange (`EPS_GAS_DIFF_RESID` = a declared **contraction guard**, gain ≲ 1/12 structural;
accuracy owned by the order gates). Transport (μ, Pr→k, Sc→ρD, c_p/c_v) = pure config data from
the ONE owner (`WallLaw`'s constants). **Wall ownership:** resolved diffusion flows only through
gas↔gas faces (aperture-aware); wall-law faces contribute nothing (unit-proven bitwise) — no
double count. Declared viscous BCs incl. **Continuative** (zero-normal-gradient open plane) +
wall-velocity schedules (a moving wall does ledgered work). The S2 seam refusal retired. Battery:
Poiseuille at 33.5× the explicit viscous bound (0.73%), Taylor-Couette swirl (0.80%), **exact
recovery-Couette** 3.584 vs 3.581 K analytic (replaces the flat-plate mini-sim — same physics
balance, exact solution), thermal_bl erfc + species at Sc≠Pr, whole-operator MMS **orders
1.92–2.21 on all six components**, a four-class march (the S4 configuration in miniature), N_θ>1
refusals, thread bit-identity. **Certificates byte-identical** (gas diffusion is opt-in config;
stations don't schedule it — S4's ◆C1 turns it on with real transport tables). Deferrals:
skin-friction gas debit (mount-reaction wave); species-enthalpy flux + TableEos-T refresh (S4
spine); COUP-8 row + config grammar (S4); serial assembly (GPU wave); θ-diffusion (S8). Finding:
impulsive walls/drives at stiffness ring the truncated sweeps past gas positivity — declared
ramp schedules (the COUP-7 discipline) are the cure; RL10 already carries its injector ramp.

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
`configs/rl10_full.toml`); the cure (v1.5: physical march only) is the **plan's Phase 1–2 physics**
(gas diffusion + implicit integrator + the cold/unburnt branch), which makes those states
representable instead of refused. **Viz gate (Ben):** no dataviz until THE RUN (plan S19 — the
full-3-D spark-to-steady certification overnight); `runs/*/fields.csv` keeps accumulating as the
future feed (r/z/ρ/u/p/T/Z/M per cell + solid liner T).

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
  c\* per −3e5 J/kg); `wall_heat` = the one Colburn-class law (**±20–30% band**), the ONE
  transport-constant owner; **`gas_diffusion` = F_visc (S3)** — the class-D gas occupant
  (per-component symmetric CG + fixed-Picard cross terms; ω-form swirl; total-energy T-solve;
  suppressed at wall-law faces; declared viscous BCs incl. Continuative; N_θ = 1, S8 re-keys).
  Fine-dial establishment is cured by the plan's Phase 1–2 physics (pseudo-transient DELETED,
  v1.5); S4 wires gas transport into the engine config + real tables.
- **`crates/engine` — the sandbox seam**: config → assembly (content-verified contour →
  FND-3 cut-geometry grid; refusals: cooling-with-zero-liner, closed-mode-never-engages,
  adiabatic liner holes) → **the ONE SDC-IMEX step (S2)** — wall law on closure-vector patches
  with the SRD-neighborhood debit, liner conduction + coolant Robin inside the class-D implicit
  solve, audit armed, Δt = gas CFL alone — plus **the COUP-3 §3.5 closed-mode expander fixed
  point** (session 12: `turbopump_expander` boundary object, drive_power ← jacket heat_pickup,
  Aitken ≤ 1e-8, engages post-establishment; consumes the post-sweep converged wall-heat
  integral) → SOLV-7 readout (+ measured inflow-plane ṁ honesty signal) + fields-CSV viz feed +
  **fault-tolerant crash artifact on halt** (`runs/<name>/crash_fields.csv`). Presets:
  `rl10_coarse.toml` (dial 5, ~155 s laptop on the S2 spine, the certified tier),
  `rl10_calibrated.toml` (closed mode), `rl10_full.toml` (dial 16 — KNOWN LIMIT: awaits the
  plan's Phase 1–2 physics). Certificate regen: `station5_rl10_certificate` bin (recorded
  readouts + Ferson rescoring, gate 5).
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

## Where we are in the plan (PLAN_CHEMICAL_SANDBOX.md — the authority on "what next")

- All five ladder stations earned (session 12); the station-5 certificate scores blind +
  calibrated boxes vs the TM-107318 reference p-box (details: certificate + SESSION_LOG).
- **Session 13 = plan S1 DONE** (the v1.5 amendment wave); **session 14 = plan S2 DONE** (the
  real integrator); **session 15 = plan S3 DONE** (the missing forces: gas F_visc as the
  class-D gas occupant, wall-law suppression, exact-analytic mini-sim battery). **NEXT = plan
  S4**: real properties + ◆C1 — transport tables μ/k/c_p(T, p, Z) through the FND-7 seam
  (OFFL-3 Cantera feed, ~10–20% declared bands); liner gets physical ρc_p (continuation device
  retired); **◆C1: 2-D RL10 with full diffusion physics re-settles with NO schedule tuning**;
  fine-dial 2-D establishment attempted — read the plan's §5 S4 entry + FND-7 §3 before coding.
- The old per-item deferral list (station-4 fixture rewire; Bartz oracle scoring; digest v4;
  COUP-5 ensembles; FND-3 PLIC/slot class) is absorbed into the plan's phases: §4 maps each to
  its session.

## Cross-cutting deferrals (module headers carry the per-module lists)

- **COUP-3** SDC-IMEX class-D supersedes every explicit scaffolding integrator and brings
  Robin-Robin Picard/Aitken wall coupling.
- **COUP-7** boundary objects supersede `StagnationInflow`, caller-supplied BC closures, and the
  coolant Robin film.
- **FND-3** voxelization: partial apertures/cut cells retire stair walls + transpiration; brings
  mask-disjointness validation and aperture-aware interface classification.
- **OFFL-3** remaining products: frozen-path + unburnt-reactant surfaces (plan S5), S_L/τ_ign
  ignition closures (plan S6), expansion oracles, B′ (SOLV-8), transport feed (OFFL-5, S23 —
  plan S4); per-point sigma columns (COUP-5).
- **FND-7 spine**: per-cell transport/material data; retires the sole-instance two-material
  ceiling and degenerate constant-transport occupants (gamma-law pattern).
- libm/powf note: certificate byte-diffs are valid on the pinned dev host only.

## Working rules & non-negotiables

- **Read the owning doc §3 before coding its operator.** Explain results to Ben as physical
  systems. Every session ends gates-green and committed.
- **Two-language rule:** runtime = 100% Rust; offline = Python; versioned HDF5 the only seam.
- **Test-first:** analytic/manufactured before benchmark before hardware (VAL-1); no red tests,
  no half-refactored solver at session end.
- **Determinism** (META-1 §2, GPU policy §2.5 — Ben-confirmed 2026-08-19): bit-exact per build
  **per device**, CPU and GPU (gather kernels, fixed-topology reductions, no physics atomics);
  chaotic + relaxed = load refusal; cross-device = tolerance/ECT.
- **Fail loud, halt clean, never guess** (META-1 P6): out-of-envelope ⇒ refuse/flag, never clamp.
- **Every result is a distribution + pedigree** (S3); the p-box is never collapsed.
- **Firewall** (VISION_SCOPE §10): no MCNP, no restricted codes/data; pulsed-fission stops at
  published envelopes. No amendment path.
- A change contradicting a Reviewed doc ⇒ **amend the doc first** (change log; Layer-1 needs a
  VISION_SCOPE §15 entry). Never code around a doc.
- **Do not build:** SOLV-5, OFFL-4 (unreviewed); COUP-1 is a retired tombstone.
- **Ben's viz gate (2026-08-19):** NO dataviz work until THE RUN exists (plan S19 — the full-3-D
  spark-to-steady certification overnight); keep the fields-CSV feed boring and complete meanwhile.
- **Compute strategy:** laptop = per-session mini-sims (the test battery); desktop (RTX 4080,
  24-h cap) = the ◆ checkpoint overnights (plan §3/§6).
