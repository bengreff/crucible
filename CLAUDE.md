# CLAUDE.md — Session-Start Context for CRUCIBLE

CRUCIBLE is an in-vacuum propulsion physics sandbox (chemical → antimatter): one unified 3-D
adaptive-dimension field simulator, every result a distribution + pedigree. **This codebase builds
the instrument; the research uses the instrument.** No engine-specific features, ever.
`VISION_SCOPE.md` (v1.6) outranks everything; its §15 is the amendment log.

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

## State (2026-08-28 — session 24 = plan S12 done; chemical sandbox at S12/20, Phase 4)

**This section is current state only.** Per-session build history (measured deltas, findings,
review waves) lives in `SESSION_LOG.md`; the by-area capability detail is under "What exists" below;
the standing rulings are in the plan (`PLAN_CHEMICAL_SANDBOX.md` §1) and VISION_SCOPE §15.

The instrument today runs the **chemical regime end-to-end on the one unified 3-D field solver**: a
deterministic **SDC-IMEX step** (explicit hyperbolic + implicit diffusion + cell-local stiff sources
+ Robin-Robin gas–wall coupling, **conservation-audited every step**); gas **viscous/thermal/species
transport** as the class-D occupant; real **Chapman-Enskog transport properties** over the FND-7
spine; **burn-progress ignition chemistry** (bistable-Nagumo front + stiff auto-ignition + blended
thermochemistry); full **3-D azimuthal** capability (r=0 axis, mixed-N_θ reflux, θ-CFL + controller);
**CSG/STL geometry** with cut cells + PLIC; static (r,z) refinement; and typed **COUP-4 verdicts**
(WORKS / DOESN'T-WORK + mechanism/location/time). Five certificate stations earned; RL10 F and Isp
**overlap the record blind**.

**Session 24 added the GPU tier** (compute-strategy rule below + `docs/gpu-box.md`): the two hot
kernels measured on the RTX 4070 Ti SUPER — hyperbolic sweep **~1×10⁸ cell-updates/s**, EOS
projection **~4.9×10⁸/s**, STREAM **401 GB/s**; determinism validated (bit-identical reruns; CPU↔GPU
ECT ~5×10⁻¹⁰); the `crates/gpu` Rust→CUDA binding built. Finding: the sweep is **occupancy-bound
(~25% occ), not f64-bound** — the throughput target is reachable via register reduction at S13/S14;
THE-RUN envelope re-sized to ~1.7–14 h tuned (§3).

**Doctrine (VISION_SCOPE v1.6 / plan ruling #14, 2026-08-27):** the **torch/ASI flame-holder object is
DELETED** — flame-holding is **emergent** (resolved recirculation), the only start boundary-inputs are
the **bounded spark + injected fuel**, so the startability verdict is not circular; **refinement is
re-scoped to turbulence (S16)**, not the flame front. The S7/S8 torch-tier code + `rl10_startup*.toml`
are **DEPRECATED** (removal rides S16 so ◆C2/◆C3 keep a starter).

Design complete: all 29 critical-path Layer-2 docs **Reviewed 2026-08-14** (register closed 68/68,
deleted; findings live in doc change logs + git). Every session gates-green + committed;
`scripts/check.sh` = the 5-gate battery (fmt, clippy, cargo test, offline pytest, certificate regen +
diff); **253 Rust + 50 Python tests**. **Blind rule v1.4.1**: blind = mechanical input-blindness;
every certificate declares `development-observed: yes/no`; the RL10 campaign is **open development**.
**NEXT = plan S13 (full residency).**

**VISION_SCOPE v1.5 (Ben, 2026-08-19) — the session-13 rulings, all doc-amended:** accelerated
convergence (pseudo-transient/local-Δt) is **DELETED** — every certified result is a **physical
march incl. start-up** (COUP-3 §3.6 tombstone); external-system timelines may be **compressed as
declared schedules** (COUP-7 §3.2.2), internal physics never; ignition = the **burn-progress
field c** (SOLV-4 §3.6: unburnt↔equilibrium blend, flame-speed + induction closures, igniter =
energy-deposit object, `NEVER_IGNITED`/`FLAMEOUT` halts); **full 3-D is the product tier**;
GPU (RTX 4080, 24-h cap) on the critical path with **bit-exact-per-device determinism**
(META-1 §2.5); no solids; no viz until THE RUN (plan S19).

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
  c\* per −3e5 J/kg; **also binds the S5 c = 0 unburnt-reactant surface with no new code** — same
  (p, h, Z) schema, so the cold branch is a table, not an occupant); **`BurnBlendEos` + `Combustion`
  (S6)** — the burn-progress blend (energy-conserving flamelet: both branches at the cell's (p, h, Z);
  pure-limit threshold = 2·`BURN_COMPLETE`, one owner) + the SOLV-4.4 bistable-Nagumo source riding
  `Euler::eval_rhs` as class A, with `IgnitionColumns` on the (p, T_u, Z) closure surface and
  `reacting_measure`/`consumption_rate` = the COUP-4 R-floor diagnostics (N_θ = 1; explicit tier —
  stiff class-R and the blend↔class-D coupling are S7); **`transport` = THE FND-7 §3.3 spine seam (S4)** — the ONE provider of
  (μ, k, c_p, c_v, ρD, ∂h/∂Z, Pr) over the medium state, two occupants selected as config data
  (`transport_constant` = the declared set the stations use, bit-identical to sessions 7–15;
  `transport_table` = the OFFL-5 §3.1a surface on the local (p, h, Z) state); `wall_heat` = the
  one Colburn-class law (**±20–30% band**, now a declared `band_factor` on h), **stateless** —
  transport is a spine operand; **`gas_diffusion` = F_visc (S3, per-cell at S4)** — the class-D
  gas occupant (per-component symmetric CG + fixed-Picard cross terms; ω-form swirl;
  total-energy T-solve with the spine's own c_v as the slope; **the species-enthalpy flux
  `Σ h_k j_k` = ρD·∂h/∂Z·∇C**; per-cell transport, faces = the two-cell arithmetic average;
  suppressed at wall-law faces; declared viscous BCs incl. Continuative; N_θ = 1, S8 re-keys).
  Fine-dial establishment is cured by the plan's Phase 1–2 physics (pseudo-transient DELETED,
  v1.5); the cold/unburnt branch (S5–S6) is the remaining leg.
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
  review fix; fresh-holdout CI gates on BOTH pinned artifacts); **`FrozenReactantEngine` (S5)** =
  the gas-phase ideal-gas frozen reactant mixture (the c = 0 unburnt branch, same CEA enthalpy
  reference as the equilibrium surface); **`KineticEfficiencyBand` (S5)** = the shipped
  frozen↔shifting model-form band; **`IgnitionEngine` (S6)** = the Cantera `FreeFlame`/const-`P`
  reactor closure generator (`h2o2.yaml`). Production tables (S6-close versions — the OFFL-3 0.6.2
  ignition-headroom envelope set, all strict extensions of their predecessors, old nodes bit-exact):
  `lox_lh2_v0.1.0.h5` (station-3 certificate — untouched) + **`v0.4.0`** (station-5 + the blend's
  burnt branch: gas-only metastable, Z ≈ MR 4.4–5.5, p ∈ [10 Pa, 7 MPa], h ∈ [−1.23e7, **+1.2e7**];
  rule: burnt ceiling ≫ unburnt) + **`lox_lh2_unburnt_v0.2.0.h5`** = the c = 0 branch, cryo-~100 K
  to **2900 K** + **`lox_lh2_ignition_v0.2.0.h5`** = `S_L`/`τ_ign` (p, T_u, Z), T_u 150–**3013 K**
  (contract: covers T_u at the unburnt ceiling) + `tables/spine/lox_lh2_transport_v0.2.0.h5`
  (re-derived from the v0.4.0 envelope; hot-side c_p two-fit tolerance 4% = measured ×1.6);
  sidecars = single pin owners. Regen ≈ 3 min (`make_station5_tables.py`) / ≈ 8 s
  (`make_unburnt_tables.py`) / **≈ 30 min** (`make_ignition_tables.py`, streams per-row progress;
  NOT in the fast gate — only the tiny regen-determinism probe is) / ≈ 20 s
  (`make_spine_transport_tables.py`).
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
  class-D gas occupant, wall-law suppression, exact-analytic mini-sim battery); **session 16 =
  plan S4 DONE** (real properties: the FND-7 transport spine filled for the chemical regime,
  physical liner thermal mass, **◆C1 met**); **session 17 = plan S5 DONE (split)** — the cold/unburnt
  branch: the c = 0 unburnt-reactant (p, h, Z) surface shipped (frozen gas H₂/O₂ via CEA, cryo-valid,
  bound by the existing `TableEos` with no new code), the frozen↔shifting bracket shipped as a declared
  band; the frozen-mode `{ρX_k}` **field-advection** widening + the (p, h, {X_k}) advection surface +
  per-species diffusion (superseding the single Sc) split to **S5b** (may ride S8), no consumer before
  S18; **session 18 = plan S6 DONE** (ignition: the c-field + blended thermochemistry + the
  bistable-Nagumo rate law + the pinned spark; the S6-close envelope set — burnt v0.4.0 / unburnt
  v0.2.0 / ignition v0.2.0 / transport v0.2.0, all strict extensions; battery 5/5 with `flame_1d`
  grid-independent at 2.3%; the adiabatic-flameout finding re-scoped `quench_box` to S7).
  **Session 20 = plan S8 DONE** (azimuthal capability — axis crossing, mixed-N_θ reflux, θ-CFL +
  controller, N_θ > 1 combustion, the S_T-CFL + first turbulent consumer; the F_visc θ-extension
  split → S9 by the improvisation rule, recorded in §8 v1.8).
  **Session 21 = plan S9 DONE** (full geometry: the geom3d kernel — CSG SDF trees, exact-winding
  STL import, jittered voxelization with the measured C_jitter bound, PLIC; 3-D apertures incl.
  θ-faces at uniform N_θ with per-sector SRD and the `stl_toy_chamber` end-to-end gate; the
  F_visc θ-stress tensor completed at uniform N_θ — ◆C3 unblocked; the three S8 carries incl.
  the N_TAU_REFREEZE finding and its restated contract; plan §8 v1.9).
  **Session 22 = plan S10 DONE** (refinement + the AMR gate: the conservative static (r,z)
  level-interface primitive built + gated — the meridional sibling of the ring reflux, annular-metric
  single-difference well-balance, FND-2 §3.6.1; the AMR go/no-go **measured NO-GO** on dynamic
  front-tracking — timeline grid-independent, refinement buys 1.93× sharpness at 4× cost/level, so
  S10b/S10c not inserted and the static-plus-closure-set-front fallback stands measured; two S9
  carries — the coarsest-reproducing-N_θ floor + the face-order unification; plan §8 v1.10).
  **Session 23 = plan S11 DONE** (◆C3 — the first genuinely 3-D engine run + review wave A: the
  engine N_θ > 1 assembly (`build_with_geometry_theta` on the revolved contour, N_θ = 1 bit-identical),
  θ-aware `run.rs` readouts + the **point-in-θ spark**, and the **combustion cut-θ D_c stencil made
  per-sector** so the revolved cut RL10 lights at N_θ = 8 (θ-varying still refused via
  `geometry_is_theta_uniform`); ◆C3 = ADIABATIC flow+combustion — the cooled-wall/F_visc/conduction
  couplings stay the typed S11/S9 solver refusals; the laptop-mini verdict = **DOESN'T WORK
  (FAILED_TO_REACH): p_c −92.2% / F −94.8% of commanded** — the 3-D START MACHINERY proven, "doesn't
  fully light" the honest outcome; review wave A over phases 1–3 found no confirmed bug, the transport
  Pr-rail LOW finding fixed; certificates byte-identical; plan §8 v1.11).
  **NEXT = plan S13 (full residency)** *(S12 done — see the session-24 State block at the top)*: entire SDC step on-device (flow, diffusion cycles, projection,
  progress, conduction, closures); CPU orchestrates; same-build rerun bit-identity; the CPU↔GPU
  tolerance cross-check on real fixtures; FND-6 checkpoint/restart (the 24-h-cap instrument). The S12
  spike proved the pieces (measured throughput + determinism + the `crates/gpu` Rust→CUDA binding); S13
  is residency + the register reduction toward the throughput target (the sweep is occupancy-bound at
  146 regs/thread — see the S12 finding). Read the S13 §5 entry + META-1 §2.5 + the owning docs' §3.
  **Owed to a later CPU wave (S11's typed refusals still standing, owners named):** the cooled-3-D
  wall-patch/conduction wave (per-θ wall patches + coupled flow+conduction — `sdc::build_wall_patches`
  and the step guards refuse N_θ > 1); F_visc on cut θ-faces; the CSG/STL config grammar + geom3d
  production consumer; the cross-pencil (r,z) flux-register integration (the S10 Option A/B decision,
  still Ben's — unbuilt); per-sector activity masks; combustion on θ-VARYING cut geometry. Read the
  S12 §5 entry + the owning docs' §3 before coding.
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
- **Compute strategy:** laptop = per-session mini-sims (the test battery); desktop (`backhouse` =
  RTX 4070 Ti SUPER, 24-h cap; reached over Tailscale SSH — WSL2 Ubuntu + CUDA 13.3) = the new GPU
  measurement tier + the ◆ checkpoint overnights (plan §3/§6). **How to reach + drive the box (SSH/WSL
  mechanics, long-job persistence, sharing): `docs/gpu-box.md` — read it before any GPU session.**
