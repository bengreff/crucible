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
- Session 6 (2026-08-16): **review fix wave** — 8-angle multi-agent review of sessions 1–5 (41
  candidates, independently verified; 10 confirmed findings) then all fixes landed: resolved-config
  replay fixed for extents configs (§6-1 fixed point now flat-serialized + regression-tested);
  loader refuses NaN/inf extents, huge worlds (MAX_AXIS_CELLS/MAX_N_THETA, named), and the θ-ladder
  diagnostic no longer overflows/hangs; digest **v2** (length-prefixed strings — Python mirror must
  implement v2), reader/writer refuse `sigma_` orphans + unaccounted datasets; symmetry indicator
  fails loud on NaN/bad floors (now `Result`); `from_loaded` refuses multiple same-type instances
  and uses checked narrowing via the new typed seam (`mechanisms_of_type`/`param_f64` +
  `grid_spec_from`); conduction sweep restructured (≤4 Morton lookups per BRICK, none per cell —
  certified bit-identical, the SOLV-1 template); **grid surface sealed** (private bricks/masks/
  n_theta; `for_each_active_cell` visitor, `fill_field`, `r_center`/`global_rz` single owners);
  certificate criteria are named shared constants and **check.sh gate 4 regenerates + diffs the
  committed artifact**; .gitignore restored (Python rules back), README/check.sh headers current;
  META-3 `bessel-j0-zeros` entry added. True test count: **58** (session-5 commit message said 54 —
  wrong, was 48).
- Session 7 (2026-08-16): **GOAL-B STATION 1 COMPLETE — the bursting diaphragm (Sod)**
  (`certificates/station1_sod_certificate.md`, regenerable via
  `cargo run --release --bin station1_sod_certificate`; criteria CI-enforced in
  `crates/solvers/tests/solv1_station1_sod.rs`; check.sh gate 4 now diffs both certificates).
  `crucible-solvers::euler` = SOLV-1 §3.1–3.3 reacting-gas subset: U = (ρ, ρu_r, ρu_θ, ρu_z, ρE,
  ρC) on the exact cylindrical metric, PPM (CW84, MOL form — no characteristic tracing; Castro's
  SDC path uses the same, so COUP-3 consumes this operator unchanged) + HLLC with Batten
  wavespeeds + well-balanced geometric sources; gamma-law EOS struct = first degenerate FND-7
  spine occupant (γ is config data); one passive composition ρC advected from day one (the
  Batten contact wave). Exact Riemann oracle (Toro; star state vs Table 4.2 to 5 figs). Registry
  row 2: `flow` mechanism (γ validity [1.001, 1.667], load-refused outside; `flow_from_loaded`
  through the typed seam). MOL SSP-RK2 = honest scaffolding (session-5 pattern), superseded by
  COUP-3 SDC-IMEX. **Data (n_z=800 finest): shock position 0.16 cells; contact 1.11 cells; star
  plateaus ≤ 7.3e-5; star-window L1 order 2.2–2.6 (formal order in smooth regions); fan interior
  ~1.0 (known centered-fan startup behavior, reported honestly); smooth advection: C at textbook
  2.0, ρ mean 2.09; closed-tube mass/energy drift ≤ 4e-16 through wall reflections; uniform gas
  at rest = BITWISE fixed point (axis N_θ=1 + annulus N_θ=8); radially-uniform tube stays
  radially uniform bitwise through the whole shock.** Two upstream fixes the certificate forced:
  grid **`face_radius` single owner** (ring i's outer face ≡ ring i+1's inner face bitwise — the
  old `r_i + dr` rounding broke flux telescoping by 1 ulp at particular radii, latent under
  conduction's tolerances; Goal-A artifact unchanged at printed precision) and z/θ sweeps use
  exact metric ratios (A_z/V = 1/dz) instead of per-ring A·F products. `[profile.test]
  opt-level = 2` (real marches in tests; f64 bit-identical across opt levels — no fast-math);
  the test profile's overflow checks caught a θ-wrap underflow that release wrapping had masked.
  Euler deferrals (loud refusals, owners in `euler/mod.rs` header): r_min=0 with N_θ>1
  (cross-axis parity gather), mixed per-brick N_θ (COUP-2/3 refluxing), apertures/cut cells
  (FND-3 not yet consumed — worlds are full boxes).
  **Same session, second wave: SOLV-1 §6-2 whole-operator MMS done** (`euler_mms.rs` +
  `tests/solv1_euler_mms.rs`; new artifact section) — manufactured field with radial flow, swirl,
  and axial flow all active (one shared mode, analytic residual assembled exactly by product rule
  — the terms Sod structurally cannot exercise since its u_r ≡ 0): **all six components of U at
  observed order 1.88–2.22 on levels [16,32,64], both 2-D-axisym-with-swirl (N_θ=1) and 3-D m=2
  θ-mode (N_θ=n).** Norm finding, recorded in the artifact: limiter clipping at the θ-mode's
  smooth extrema is locally 1st-order over O(h) measure — L2 shows a ~O(h^1.6) tail while **L1
  (the shock-capturing verification norm, SOLV-1 §6-1's own choice) retains formal order**;
  MMS gates therefore assert L1 per component. True test count: **79**.
- Session 8 (2026-08-17): **GOAL-B STATION 2 COMPLETE — the De Laval nozzle**
  (`certificates/station2_nozzle_certificate.md`; criteria CI-enforced in
  `crates/solvers/tests/solv1_station2_nozzle.rs`; gate 4 diffs all three certificates).
  Machinery this station forced, all doc-seamed: grid gains **`build_with_activity`** — the
  FND-2 §3.6 ingest surface in binary degenerate form (FND-3 fractions/apertures arrive through
  the same seam) + `is_active`; Euler sweeps are **mask-aware** (pencil lines decompose into
  maximal active runs; interior run boundary = stair wall face, domain edge = configured BC);
  **`FlowBc::StagnationInflow`** (reservoir isentrope at interior-extrapolated u; vacuum-limit
  refusal); **`Euler::wall_normal` slip-ghost walls** — stair faces mirror about the true contour
  normal (ghost-cell immersed boundary), cutting spurious waves from O(slope) to O(h·curvature);
  SOLV-7 §3.1/§3.2-subset plane diagnostics (ṁ, momentum+pressure thrust integral, emergent p_c
  as inlet-plane stagnation avg per N11). Fixture: parabolic-radius contour, 9.5° max slope,
  area ratio 1.5625, reservoir-fed, exit supersonic. **Findings the physics forced:** (1)
  Anderson's `1+2.2(z−1.5)²` profile revolved literally is a 44° wall — the 2-D flow genuinely
  departs ~2× from quasi-1-D (steep-wall shocks); the station contour must sit inside the
  oracle's validity envelope. (2) Naive stair-mirror walls in supersonic flow shed a wave train
  (a shock crossed the axis at z≈5); slip ghosts fixed it (steady resid 4.5e-2 → 7.6e-4 at
  finest, shock gone). **Data: Cd = ṁ/ṁ_ideal = 1.0060/1.0022/1.0020 across the ladder (choked
  mass flow to 0.2% at finest, |Cd−1| non-increasing); emergent p_c/p0 = 0.99999; centerline
  area–Mach within 3.6–4.2% past the entrance band (compound of stair + the oracle's own 2-D
  centerline floor); C_F within 4.3% of ideal vacuum C_F; masked-cavity uniform gas at rest =
  bitwise fixed point under BOTH wall treatments; declared slip-wall transpiration (plane-ṁ
  spread) 7.2→5.2→3.8% shrinking with h.** True test count: **83**. Deferrals now owned by
  FND-3/COUP-7 waves: partial apertures + cut cells + State Redistribution retire the stair
  wall + transpiration; the COUP-7 injector object retires StagnationInflow; check.sh notes the
  station-2 ladder as the battery's heavy item (split to milestone tier if it grows).
- Session 9 (2026-08-17): **GOAL-B STATION 3 COMPLETE — the flame seam**
  (`certificates/station3_flame_certificate.md`, regenerable via
  `offline/scripts/station3_flame_certificate.py`; criteria CI-enforced in `offline/tests/`
  pytest + `crates/tables/tests/{fnd5_python_seam,station3_tables}.rs`; **check.sh is now 5
  gates** — gate 4 = the offline pytest battery, gate 5 diffs all four certificates). **The
  Python side begins**: `offline/` = `crucible_offl` on Python 3.13.7 (`offline/.venv`,
  bootstrap line in check.sh), exact pins in `pyproject.toml` (h5py 3.16.0, numpy 2.5.2,
  cantera 3.2.0, cea 3.3.2 — all stamped in META-3). **Digest v2 mirrored byte-for-byte**
  (`tables.py`: stdlib-pure digest half + h5py writer to the reader.rs schema; golden vector
  asserted in both languages); committed cross-language fixture (UTF-8 byte-length probe, all
  three sigma markers, log rules, rng_seed) written by Python, opened by the Rust loader under
  the Python-stamped pin. **NASA CEA is pip-installable as `cea` 3.3.2** — the modern
  Apache-2.0 github.com/nasa/cea re-implementation, native Python bindings, bundled Glenn DB,
  RP-1311 samples in-tree (RocketCEA never needed). `chemistry.py` wraps it behind SI
  boundaries ((p,h,Z) HP solves; (p_c,MR) IAC rocket; Z = 1/(1+MR); h_SI = h/R·R_CEA pinned by
  round-trip). **Production tables committed** (`tables/chem/lox_lh2_v0.1.0.h5` + pins
  sidecar): (p,h,Z) equilibrium surface 41×31×13 log-p (T, ρ, γ_eff, a, M̄, condensed_fraction,
  6 X_k) + (p_c,MR) performance reference 11×11 (c*_ideal ±2.7 m/s, T_c ±10.4 K stored
  bounds). Error-bound honesty loop the §6-4 gate forced twice: midpoint-only holdout
  understated a minor-species onset (→ midpoints + ¼-offsets), then 1.28× sampling variance on
  a coarse grid (→ declared ×1.5 margin); enforcement = fresh disjoint ⅜-offset CI sweep on the
  committed artifact — violation means refine the grid, never relax the gate. **Data: RP-1311
  example 8 reproduced through the pipeline to ≤2.7e-4 rel across 14 quantities (T_c 3383.845 K,
  c* 2332.34 m/s, compositions, Isp ladder); CEA↔Cantera two-solver agreement ≤0.13% T, ≤0.05%
  M̄ (no shared code, different thermo fits); frozen↔shifting bracket 3.91% at ε=25 — contains
  the JANNAF 0.8–1% kinetic band, sampled by COUP-5 as the S19 epistemic dimension;
  (p,h,Z)↔(p_c,MR) coordinate consistency ≤2.5 K along the design line; the RL10 chamber
  through the whole seam: Rust pin-verifies and interpolates T within 0.42 K and c* within
  0.05 m/s of CEA direct, Rust≡Python interpolation to 1e-6.** Honest note: the surface is the
  equilibrium *including* condensed H₂O in deep-cold rectangle corners (`condensed_fraction`
  column reports it; no engine trajectory goes there; full-envelope bounds are kink-dominated,
  gas-region ~50× tighter). Deferrals owned (certificate + `chemistry.py` header): frozen-path
  surface (SOLV-1 frozen-advection wave), expansion oracles (station 5), B′ (SOLV-8), transport
  feed (OFFL-5, S23), config→tables pin wiring (FND-4 §6-4, sidecar TOML is the interim
  record), per-point sigma columns (COUP-5 wave). True test count: **89 Rust + 23 Python**.
- **Ben's post-checkpoint BIG goal (2026-08-17): data visualizations** of CRUCIBLE results once
  the station ladder yields meaningful data — its own wave after the stations; keep FND-6
  results-bundle design viz-friendly meanwhile (auto-memory `project-dataviz-goal`).
- Session 10 (2026-08-17): **units wave + GOAL-B STATION 4 COMPLETE — the cooled wall.**
  **Wave 1 — units (META-2 §4 ★, review finding 9 CLOSED; own commit):** `crucible-units` =
  sole owner of pinned uom 0.38.0 (SI-base constructors, `si()` kernel-boundary extraction,
  the TemperatureKind ΔT guard — `delta_t` → `TemperatureInterval`); config gains
  dimensioned param accessors (named for their unit, paired with the FND-4 name-suffix
  convention) + typed extents views; `ConductionSetup` carries typed κ/ρc_p; kernels stay
  documented-SI f64; the serialized resolved-config grammar is unchanged (§6-1 replay).
  **Wave 2 — station 4** (`certificates/station4_cooled_wall_certificate.md`; criteria
  CI-enforced in `crates/solvers/tests/solv1_station4_cooled_wall.rs`; gate 5 diffs all five
  certificates). Machinery, all doc-seamed: **grid regions** (FND-2 §3.6 ingest widened
  binary→ternary gas/solid/exterior as config-time data; `solid_mask`,
  `build_with_regions`, `gas_solid_faces()` deterministic wall-face enumeration,
  `interface_area_per_theta`); **conduction generalized** (Domain selector,
  `InteriorFaces` {gas-exchange closure, exterior BC}, `FaceBc::Robin` as film + half-cell
  series, unallocated-neighbor handling, `validate()` fail-loud — Goal-A behavior
  preserved, artifact unchanged); **`wall_heat.rs` = SOLV-1 §3.5's one law** (Colburn-class
  wall function: C_f/2 = 0.0225·Re_y^−1/4, St = (C_f/2)·Pr^−2/3, molecular floor k_g/y —
  smooth to u_t = 0, no regime branch; recovery T_aw = T + Pr^{1/3}u²/2c_p; **±20–30%
  declared band**; registry row 3 + typed `wall_law_from_loaded`); **coupled exchange**:
  one wall-function evaluation per face per step applied to both sides (conservation by
  construction; the same-face-area guarantee is session 7's `face_radius` single owner),
  explicit flux-matched splitting at the gas CFL dt, guarded against both thermal
  stability limits (fail-loud) — honest scaffolding, superseded by COUP-3's class-D
  Robin-Robin Picard/Aitken. Fixture: M = 2 / 800 K gas (recovery ≈ 1374 K) in a straight
  duct, steel-class liner, coolant Robin film. **Data: the coupled system reproduces the
  cylindrical series-resistance conjugate solution to 9.4e-4 worst-case (oracle built from
  gas state + declared coolant data only — the simulated solid never enters); three-way
  energy ledger (gas enthalpy deficit / wall exchange / coolant extraction) closes to
  3.2e-3 and 1.4e-7; q ≈ 1.06 MW/m² rocket-scale, liner ΔT 762→400 K; stepped-cavity
  stair interface (r- and z-faces) conserves gas→liner energy to 4.2e-12; Robin annulus
  analytic anchor ≤ 2e-3; uniform rest at coolant T = bitwise fixed point of the full
  coupled step; reruns bit-identical.** Findings: (1) near-wall sampling — the
  wall-adjacent cell is itself cooled, reading recovery T_aw 6% under free-stream at
  2.5 mm cells; declared, inside the band, shrinks with FND-3 + finer wall cells.
  (2) liner ρc_p = 200 is a declared steady-state continuation device (steady solution
  independent of ρc_p). Deferrals owned: Bartz nozzle-envelope cross-check → station 5
  (the nozzle+liner assembly exists there; its stair machinery is certified here);
  Robin-Robin-in-class-D → COUP-3; coolant closure → COUP-7 jacket object; per-cell
  transport → OFFL-5 spine. True test count: **102 Rust + 23 Python**.
- **GOAL B (confirmed by Ben 2026-08-16): the BLIND RL10 — M2 itself.** Five certificate stations,
  each a physical system, ladder order strict (VAL-1: analytic → benchmark → hardware):
  **(1) Bursting diaphragm (Sod)** — SOLV-1 HLLC/PPM Euler at N_θ=1 vs exact Riemann; L1 at formal
  order in smooth regions, conservation to round-off. **✓ COMPLETE (session 7, 2026-08-16).**
  **(2) De Laval nozzle** — cold-gas choked flow vs isentropic area–Mach relations; emergent p_c.
  **✓ COMPLETE (session 8, 2026-08-17).**
  **(3) Flame (seam)** — OFFL-3 Python CEA/Cantera → HDF5 (p,h,Z) surfaces vs CEA manual cases;
  h5py mirrors `crates/tables/src/digest.rs` byte-for-byte. **✓ COMPLETE (session 9, 2026-08-17).**
  **(4) Cooled wall** — the one wall-function heat law (SOLV-1 §3.5/COUP-3 class-D) + certified
  conduction; declared ±20–30% band. **✓ COMPLETE (session 10, 2026-08-17).**
  **← NEXT: (5) the BLIND RL10 assembly.** Big pieces before the blind run, each its own
  wave: (p,h,Z) surface wired into the Euler operator as its EOS (SOLV-1 §3.4 shifting
  mode — the FND-4 §6-4 config→tables pin grammar lands here); nozzle+liner assembly
  (Bartz cross-check rides along); COUP-7 injector + cooling-jacket boundary objects;
  the expander-cycle closure so p_c emerges; SOLV-7 thrust/Isp readout vs TM-107318
  under the VAL-2 §9 blind rule.
  **(5) Blind RL10** — 1966 design geometry + universal closures + technology-class data only
  (§9 blind rule: nothing measured on the engine under test); predict 73.4 kN / Isp ≈ 444 s ≤ 2%;
  thrust measured from exit momentum flux (SOLV-7), never assumed; uncertainty from declared
  closure bands (full COUP-5 UQ is a later wave).
- Working rules for sessions: read the owning doc §3 before coding its operator (SOLV-1 for
  stations 1–2); explain results to Ben as physical systems (auto-memory has the preference);
  every session ends gates-green and committed.
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
