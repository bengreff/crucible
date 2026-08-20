# SESSION_LOG.md — CRUCIBLE coding-session history

The detailed per-session record: what was built, the measured data behind each
certificate, findings the physics forced, honest-scaffolding declarations, and
review waves. `CLAUDE.md` carries only current state and points here; consult
this file when you need the story behind a surface (why a deferral exists, what
a certificate number was when first earned, which review found what).

Newest entries at the bottom. Every session listed here ended gates-green and
committed; certificate data quoted here is the value at the time it was earned —
the committed artifacts in `certificates/` are the living record.

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
- Session 10 review wave (2026-08-17): **10-angle multi-agent review of sessions 7–10 + gap
  sweep** (~40 verified candidates; every confirmed correctness bug fixed same-session, gates
  green; certificates stayed byte-identical through the conduction fix — no certified number
  stood on a bug). Headline fixes: (1) **conduction rate-buffer corruption** — pass 1 wrote
  only domain cells but the reused buffer was copied whole per brick and the unmasked pass 2
  integrated stale rates into out-of-domain cells' T (latent under full-box grids, activated
  by Domain::Solid; found independently by three angles; regression test pins T_solid ≡ 0 on
  gas cells). (2) **Digest v3** — v2 left `schema_version` (reader-gated!) and the
  both-sigma-forms corner unpinned: exactly one sigma form is now refused by writer AND reader
  in both languages, schema_version joined the hashed trailer, shape/σ-length validated before
  any I/O (Python `assert`s → real exceptions, `-O`-proof); all pins re-stamped (golden
  a3b0bc59…, fixture 5f91b7ab…, production via sidecar — the fresh-process regeneration also
  reproduced bit-identical table data, exercising cold-start CEA determinism). (3)
  **pins.toml = single owner**: Rust station-3 tests parse the sidecar (`toml` dev-dep), the
  certificate script uses `tomllib` — hand-copied digest constants are gone. (4) **`cea==3.3.2`
  was missing from pyproject pins** (venv-only install; the engine's version was unrecorded).
  (5) all-Exterior worlds refuse at the build seam (was: Ok(0 bricks) then `brick(0)` panic).
  (6) **release profile gains overflow-checks** — certificates had regenerated unguarded since
  session 7 moved gate 5 to `--release` (the θ-wrap-underflow bug class). (7) smaller fixes:
  coolant ledger reads pre-step T (transient closure was biased); `face_geometry` sized from
  `grid.spec()` not fixture constants (station-5 landmine); station-2 plane diagnostics
  eos-threaded (hardcoded γ=1.4 decode); station-4 criteria → named shared constants in src
  (test + certificate bin, session-6 convention restored); `InteriorFaces::gas` under
  `Domain::FlowActive` refused loudly (was silently dead); coupled-march panics → typed
  `CoupledError`s; symmetry indicator skips solid-only bricks (0/0-NaN read as "corrupted
  data"); `EquilibriumEngine` validates Propellant coherence (overlapping streams gave NaN
  weights); station-1/2 doc-vs-measurement mismatches corrected in place (star-span includes
  the contact deliberately — u,p continuous there; no THROAT_EXCLUSION exists); check.sh
  refuses untracked certificates; VAL-2 §3.2 erratum c\* `in/s`→`ft/s` (changelog 0.2.1,
  mirrored in META-3). **Perf items deliberately deferred INTO the station-5 waves that
  rebuild those surfaces** (recorded in the NEXT block): tables `BoundColumn` handle
  (`interpolate` re-parses rule strings + allocates per query — ~10⁸ calls/run at station 5),
  Euler per-step scratch/allocation churn, coupler face-cache + one-stepper merge, shared
  mask-aware plane-diagnostics module. True test count: **104 Rust + 25 Python**.
- Session 11 (2026-08-18): **STATION-5 WAVES (a)+(b) — the EOS seam and the FIRST
  CONFIG-DRIVEN ENGINE.** Rulings first (own commit): VISION_SCOPE **v1.4.1** — blind =
  mechanical input-blindness, `development-observed: yes/no` declared per certificate; RL10
  campaign = **open development** (Ben); VAL-2 0.2.2 mirrors.
  **Wave A (commit 2):** `BoundColumn` (FND-5 §3.3) — bind-time units gate/rule parse/ln-hoist,
  allocation-free queries **bit-identical** to `interpolate()` on the production surface
  (asserted); **SOLV-1 §3.4 shifting mode** — `EosLaw` seam (monomorphized), NPRIM=8 with aux
  (e, Γ₁) slots (Castro/PeleC general-convex-EOS treatment: Batten wavespeeds from local/
  Roe-averaged Γ₁); `TableEos` = per-cell equilibrium projection at (p, h=e+p/ρ, Z), Illinois
  regula-falsi (fixed tol/iters); η_c\* knockdown hook = `h_offset` (S18, source-level;
  calibration = cycle wave). **Combustion lives in the EOS**: RL10 chamber T 3225.4 K emerges
  from (ρ,e,Z) alone; uniform equilibrium rest = bitwise fixed point of the full step; closed
  hot/cold tube conserves to <1e-12; **all four gamma-law certificates byte-identical through
  the seam** (gate 5). FND-4 §6-4 `[tables]` grammar (explicit pair or pins-sidecar via
  `load_str_with_sidecars`; resolved form replays purely; manifest `table_pins`); cross-process
  regen-determinism pytest.
  **Wave B (commit 3):** geometry-of-record digitized → `data/anchors/rl10_contour.csv` (Table
  E1 verbatim + Table 2.5.1 scalars + Fig. E1 planes ±0.5 in; closures declared in-header).
  **FINDING (anchor erratum):** TM Table 2.5.1's "Diameter 2.47/5.13 in" are **RADII** — c\*
  identity (Ø-reading → 600 m/s, 4× off; radius → 2381 vs 2385 record), ε=61 exit Ø 38.6 in ≈
  the known ~40 in bell, and Fig. E1's radius axis all agree (VAL-2 0.2.3 erratum; META-3
  updated). FND-3 contour grammar in `[geometry]` (content-addressed CSV like a table pin;
  **fidelity dial `cells_across_throat`** → derived isotropic extents, manifest-recorded;
  replay pure); `[operating_profile]` steady-march subset (flowthroughs/cfl/fill/pump-down).
  **`crates/engine`**: the ONE config-parameterized coupled stepper (wall law at
  `gas_solid_faces`, liner conduction + coolant Robin, flux-matched debit — honest scaffolding
  until COUP-3 class-D), COUP-7 subset rows (`flow_shifting`, `injector_prior` w/ mass-flow
  inflow BC + sonic startup cap, `jacket_coolant` coolant-side-only per D-C), O20 bindings
  exercised for real, SOLV-7 plane readout (N11 stagnation p_c, exit momentum+pressure thrust),
  fields-CSV viz feed; **`crucible run <config>`** CLI. The RL10 exists ONLY as data (CSV + 2
  preset TOMLs + pins). Table **v0.2.0** (station5 envelope: p widened to 10 Pa floor for the
  vacuum-plume fringe — an OFFL-3 R2 envelope setting, 5 s regen; v0.1.0 + station-3 certificate
  untouched). **Findings the physics forced:** (1) near-vacuum (ρ,e)→p inversion is
  ill-conditioned (constraint line grazes the ρ-contour) → projection slow path: fixed log-scan
  + golden-section, acceptance = the density column's own declared interp-error bound (the
  projection IS the per-step re-equilibration; a miss within the surface's declared error is the
  surface). (2) Low-p fill startup = injector piston shock heats past the table h-ceiling →
  quiescent fill near the operating class + **altitude-cell pump-down schedule**
  (`FlowBc::PressureOutflow` with declared p_amb(t): log-linear fill→floor; supersonic exit =
  pure extrapolation, zero upstream influence). (3) Slip-ghost stair z-faces transpire; at
  M≈4.4 the exit-lip corner cell STARVES (ρ runaway → CFL collapse → envelope exit) — mirror
  z-faces instead shock the chamber off-surface: **both are the small-cut-cell class whose
  designated cure is FND-3 partial apertures + State Redistribution (SOLV-1 §3.6)** —
  `slip_wall_z_faces` operator policy flag added (stations keep certified behavior
  bit-for-bit). **Measured: the coarse preset (dial 5, 41×119, 2584 gas cells) marches 2935
  steps / 2.5 flow-throughs in 34 s wall on the laptop**, trending F 106.6 kN falling toward
  the 73 kN class, p_c 4.2 MPa falling, jacket 11 MW (record 8.4), honestly labeled non-steady
  (resid 0.9); engine smoke test (dial 2.5, 2 FT, full pipeline incl. MPa-class emergent
  chamber + liner heating) rides the fast battery. **KNOWN LIMIT recorded in the preset:**
  settle past ~2.9 FT blocked by the exit-lip starvation — first item of the next wave, before
  the expander cycle. True test count: **122 Rust + 26 Python**.
- Session 12 (2026-08-18): **STATION-5 WAVES (b′)+(c) — cut cells, the closed cycle, and the
  FIRST RL10 CERTIFICATE.** Instrument first: `run()` returns a structured `Halt` with a
  fault-tolerant crash CSV (raw conserved `U` always; derived state where the projection holds;
  κ column) — the artifact that diagnosed every pathology this session.
  **Wave A (commit 1): the stair-corner blocker retired.** The crash artifact pinned the
  session-11 death: the shadow column behind each stair step starving to ρ ~ 10⁻¹² (dial-8
  check: same class, earlier — not resolution-curable). Cure as designated: **FND-3 §3.3
  analytic partial fractions + apertures** for the revolved contour (CSG revolved-profile
  leaf, zero sampling error, exact cylindrical measure) through a validated
  `build_with_geometry` seam (Gas ⇔ κ>0, bitwise shared faces, covered-vs-κ=0 rule);
  `wall_closure` = THE discrete interface identity (well-balance defines the wall vector; |W| =
  the smooth wall area, not the stair overcount); aperture-weighted sweeps with
  arithmetic-identity defaults (full-box worlds bit-identical — gate 5 held); **Berger–Giuliani
  SRD** after each RK stage (κ < 0.5, fixed lexicographic neighborhoods, conservation exact;
  `srd_neighborhood` public — the wall debit deposits into the same merged volume). Engine wall
  exchange rewired to per-cell closure patches. Tests: frustum volume exact 1e-13; uniform-rest
  well-balance < 1e-11 in a cut cone; closed-domain conservation < 1e-12 through shocks;
  κ ~ 4e-3 slivers at the UNCUT CFL. The 6-FT march that died at 2.9 FT completed same-day.
  **Wave B (commit 2):** table **v0.3.x** — the settled fringe PINNED the v0.2 h floor
  (T = 698.55 K × 35 cells: binding physics, quasi-clamp); v0.3.2 = **gas-only METASTABLE
  products** (declared plume model; the only CEA-convergent branch below the condensed cliff;
  deck-stamped), **Z narrowed to the premixed class** (the prior-tier field holds Z = Z_inj
  everywhere — what buys the cold floor on a rectangular grid), h ∈ [−1.23e7, **+3.8e6**]
  (ceiling sized by measured transient overshoots: piston → −1e5, backflow recompression →
  +1.33e6). `[operating_profile] p_amb_floor_pa` (declared ~1 mbar altitude cell; refused below
  the table floor) + `injector_ramp_flowthroughs` (declared valve-sequence class). **COUP-3
  §3.5 expander closed mode**: `turbopump_expander` boundary object (drive_power ← jacket
  heat_pickup, O20; TM-107318 component data = calibrated-only; algebra internal per Fork-2);
  per-step Aitken fixed point (residual ≤ 1e-8, machine-0 at steady), engages post-
  establishment. **S18 FIXED**: the session-11 h_offset wiring was a pure gauge relabeling
  (calibration trial returned the baseline BIT-IDENTICAL) — correct asymmetry: store true
  energy, interrogate at h+δ; measured slope −0.847% c\* per −3e5 J/kg. **Optimization**
  (profiler-led: interpolate = 70%): warm-started projection (uniqueness-guarded) + rayon by
  brick-row/column ownership partition — **bit-exact at any thread count (asserted)**, coarse
  12 FT 204 s → 60 s (3.4× on 6 P-cores).
  **Review wave (Ben's session-6/10 pattern; finder fan-out, verifier fan-out cut for usage —
  Sonnet for mechanical checks per Ben):** 22 findings / 8 dimensions; confirmed + fixed:
  **`TableEos::roe_sound_speed` silently ZERO at every face** (γ-law enthalpy identity under
  the CEA datum + `.max(0.0)` — now datum-free Roe-averaged c²); **warm-start root hysteresis**
  (non-monotone corner — now gated on the cold fast path's own straddle precondition);
  **vacuous slow-path acceptance** (absolute bound at fringe ρ = de-facto clamp — now
  rule-space `interp_error_bound_log`, producer-measured incl. **envelope-EDGE holdout** whose
  structural hole the review proved against live CEA; digest-v4 deferral recorded; FND-5 0.3.1
  change log); station-5 fresh-holdout CI gate (was never armed on the pinned artifact — now
  28 offline tests); perf-deck gas_only stamp; condensed-suffix gas filter; assembly refusals
  (cooling-declared-with-zero-liner, closed-mode-never-engages, adiabatic liner holes);
  inflow-plane ṁ honesty signal; pin pair-vs-sidecar precedence documented.
  **Wave C (commit 3): `certificates/station5_rl10_certificate.md`** — recorded-readout
  rescoring bin in gate 5; VAL-2 §3.1 reference box (published intervals ⊗ TM's own 0.598%
  model scatter); Ferson d by fixed quadrature; every score labeled + `development-observed:
  yes`. **BLIND (coarse, η-family × wall-band corners): F and Isp OVERLAP the record**
  (d = 0); p_c/c\*/C_F miss coherently ~5.5–5.9% (one discretization signature; →1.4–1.5%
  under the indicative dial-8 deltas). **CALIBRATED (closed): F, Isp, AND the emergent p_c
  ALL OVERLAP** — nominal p_c 463.8 psia vs 475–482; the wall-law ±25% band sweeps delivered
  ṁ 16.74→18.56 kg/s while **Isp self-regulates flat (439.0–441.0 s)** — real expander
  behavior, reproduced not imposed. **KNOWN LIMIT (recorded, owner named):** dials ≥ 8 cannot
  ESTABLISH by physical march under the honest acceptance (five schedules probed — hot/low
  fill, ramp/pump orders, altitude start — each halting loudly at a different envelope edge;
  thin gas against the 120 K liner physically equilibrates below any CEA floor); the
  pre-review dial-8 "success" rode the vacuous bound. Cure = **COUP-3 §3.6 pseudo-transient
  continuation** (the doc's own default route to steady points); grid-sequenced restart
  (FND-6) alternative. True counts: **136 Rust + 28 Python**; all gates green.
- Session 13 (2026-08-19): **PLAN S1 — the v1.5 amendment wave (docs only; no solver code).**
  Context: Ben's post-session-12 design review (full-conversation, this session's chat) produced
  the rulings now in `PLAN_CHEMICAL_SANDBOX.md` §1 — the plan of record for finishing the
  chemical-regime sandbox (S1–S20, six phases, THE RUN at S19 = full-3-D spark-to-steady RL10 on
  the RTX 4080 ≤ 24 h). Landed this session: **VISION_SCOPE v1.5** (§7.6 physical-march-only —
  accelerated convergence DELETED; compressed external schedules; ignition-as-physics; full-3-D
  product tier; §15 entry with rationale); **COUP-3 0.4** (§3.6 tombstoned, `pseudo-transient`
  key retired, forward note on all-speed acoustics for the nuclear wave); **COUP-4 0.3** (physical
  march wording; new `NEVER_IGNITED`/`FLAMEOUT` halts); **SOLV-1 0.4** (c field in `U`; §5 anchor
  budget rewritten to the full-3-D physical march + GPU envelope); **SOLV-4 0.4 §3.6** — THE new
  design: the chemical burn-progress source (c ∈ [0,1] burnt fraction; SOLV-4.4 rate law = S_L ×
  wrinkling flame propagation + τ_ign auto-ignition; unburnt↔equilibrium blended thermochemistry;
  igniter = energy-deposit object; grid-independent front speed = the acceptance gate; reacting
  measure → the halt thresholds; runtime finite-rate networks stay out — closures are offline
  surfaces); **OFFL-3 0.4** (unburnt-reactant surface + S_L/τ_ign closure products — offline
  Cantera finite-rate is in scope for table generation only); **COUP-7 0.4** (§3.2.2 compressed
  external schedules; igniter + valve/start-sequence objects); **VAL-2 0.2.4** (H₂/O₂ flame-speed
  + shock-tube ignition-delay unit anchors); **META-1 0.4 §2.5** (GPU determinism policy,
  Ben-confirmed: bit-exact per device, gather kernels, no physics atomics, fixed-topology
  reductions, >30% escape hatch); **META-3 0.7 §6.9** (burn-progress keys, `[RP]` pins due at
  S6); **META-0** COUP-3 row; **FND-2 0.5.1** (refinement staging note; §3.9 GPU-relaxed wording
  superseded by META-1 §2.5). Cleanup: `rl10_full.toml` stale comment corrected (the certified
  tier is dial 5 ONLY — the comment predated the review's retirement of the vacuous bound);
  station-5 certificate cure text regenerated (2-line diff, scores untouched);
  **REVIEW_FINDINGS.md + REVIEW_PREP.md deleted** (closed 68/68-discharged register + its
  prep checklist; findings live in doc change logs + git history; README updated); CLAUDE.md
  restructured around the plan. Adversarial review wave: two
  independent reviewers — consistency (14 findings: 8 unbumped version headers, 3 incomplete §7
  key lists, §15 date order, 2 stale register references → REVIEW_PREP.md also deleted) and
  plan-coherence (12 findings; the load-bearing four: SOLV-4.4 needed the TFC |∇c| form + matched
  front-thickening diffusion or the front speed is grid-set; the equilibrium projection had to be
  gated on c or every flammable cell burns regardless of ignition; the θ-CFL at inner rings is
  handled by the N_θ(r) profile with a centered-spark-is-near-axisymmetric argument; the 100–300 ms
  window is physically sufficient because the slow bootstrap clocks are external/compressible while
  the physical 0.33 mm liner's thermal time is ~tens of ms) — all 26 fixed in-session. Test counts unchanged (**136 Rust + 28 Python**); all gates green.
- Session 14 (2026-08-19): **PLAN S2 — THE REAL INTEGRATOR: the one deterministic SDC-IMEX
  step; every scaffolding integrator retired; the COUP-2 audit armed every step; all
  certificates re-earned on the new spine.**
  **The step (`crucible-solvers::sdc`, COUP-3 §3.1):** 2 Lobatto nodes, IMEX-Euler predictor +
  N_SDC_CORRECTIONS = 2 fixed trapezoid correction sweeps (fixed point = trapezoid, 2nd order
  both classes; explicit-part stability polynomial 1 + z + z²/2 + z³/4 — imaginary-axis stable
  to |z| ≤ 2, upwind-stable at CFL 1; stiff-part sweeps damp from the L-stable BE predictor,
  |R(−∞)| = 7/8). Class A = `Euler::eval_rhs` (the session-7 spatial operator unchanged; SRD
  applied stagewise after every node-state composition). Class D = the conduction operator
  refactored into ONE affine heat assembly (`assemble_heat`: Affine/Linear modes, diagonal
  extraction, COUP-2 ledger lines, per-face exchange-heat records) solved matrix-free by
  Jacobi-preconditioned CG in δ-form (warm start = the field itself: a settled state costs zero
  iterations; EPS_CG_RESID = 1e-12 rel, N_CG_ITERS_MAX = 512, NaN-safe acceptance →
  COUPLING_RESIDUAL). COUP-3 0.4.1 amendment landed with the code: §3.1's CG "fixed iteration
  structure" bound to the §3.7 fixed-tolerance/fixed-cap rule.
  **Robin-Robin exchange (COUP-2 §3.5), placed inside each sweep's class-D solve:** the
  wall-function h is the Robin coefficient (linear in T_solid — stable at Biot > 1);
  N_ROBIN_SWEEPS = 3 Picard sweeps with clamped Aitken (ω ∈ [0.1, 2]); the gas debits exactly
  the per-face heats of the accepted assembly (one owner — conservation by construction);
  `WallPatch`/`build_wall_patches` moved into the solvers crate, generalized to both world
  classes (cut: |W|-closure patches + SRD debit sets; box: one patch per face — station 4 and
  the engine share ONE exchange law). EPS_ROBIN_RESID = 1e-6, sized on the §3.5 principle
  (orders below the ±20–30% wall band; the residual measures only operand staleness — measured
  1.1e-8/sweep spikes during the RL10 establishment transient at the old 1e-8 gate).
  **The audit (COUP-2 §3.1), armed on EVERY step of every march:** port ledger accumulated at
  the sweeps' run boundaries (domain BCs, wall faces incl. declared stair transpiration — Σ κV
  telescoping), applied-increment source ledger (geometric/closure/external, exact κV-weighted
  increments), solid heats via the assembly's ledger lines, exchange pair cancelling;
  Δ(stored) = Σ(weighted ports + sources) checked against TOL_AUDIT[q] = K_AUDIT·ε·√N·S[q]
  (§3.1.1; S = |stored| + gross throughput; declared floors optional) — violation = halt with
  diagnosis. Closed at round-off everywhere, including 15.8k-step RL10 members. Deferrals
  (recorded): TOL constants' manifest recording rides FND-6 (build fingerprint pins them);
  the §3.1.2 mount-reaction vs SOLV-7 thrust cross-check (§6-7) arrives with the verdict wave;
  the class-D assembly is serial (perf — GPU wave brings multigrid/parallelism).
  **Retired:** `Euler::step/step_ws/advance` (SSP-RK2), `Conduction::step/advance` (explicit),
  station-4's fixture stepper and the engine's `coupled_step` + BOTH thermal-stability dt
  guards (Δt = the gas CFL alone — the point of class D; `stable_dt` survives as the explicit
  BOUND the stiffness tests measure against). InteriorFaces lost its frozen-flux gas closure
  (exchange data now flows per-sweep through the assembly).
  **S2 mini-sims (`coup3_sdc.rs`, 6 tests):** flow temporal order 2.0 by dt-Richardson (smooth
  acoustic tube); class-D temporal order ~3.0 measured (linear fixture superconverges — the
  2-sweep composition matches trapezoid through z³; gate [1.7, 3.5], the flow test is the
  strict one); class-D stable + steady-accurate at 512× the explicit bound; audit rows close
  on a live transient (θ-momentum tol legitimately 0 on a no-swirl fixture: 0 = 0); planted
  K_AUDIT = 1e-6 trips AuditViolation (the halt path); 1-vs-4-thread bit identity through the
  FULL coupled step (Euler rayon sweeps + CG + exchange + audit reductions).
  **Certified numbers, re-earned (the expected drift class — spatial operators untouched):**
  Goal-A MMS orders 1.99–2.00 (errors shift in the 3rd digit); Bessel 6.3e-3 → 9.5e-3 K — now
  marched at 4× the explicit bound (annulus at 32×, drift fixture at 8×; the anchors CERTIFY
  stiff stability now); closed-sweep drift 2.3e-16 → 7.9e-16 (CG conserves to round-off).
  Station 1: global L1 ladder and smooth orders unchanged (adv ρ mean order 2.09 → 2.57);
  shock position 0.155 → 0.845 cells (criterion ≤ 1); star-plateau u* 7.3e-5 → 8.8e-4 — the
  less-dissipative stage structure rings more behind the captured shock (STAR_PLATEAU_TOL
  re-pinned 5e-4 → 2.5e-3 with the shift recorded); closed-tube conservation at round-off.
  Station 2: Cd 1.0060/1.0022/1.0020 → 1.0050/1.0023/1.0020; steadiness residuals shift in
  kind. Station 4 (now marching the production coupled step): oracle agreement 9.4e-4 class
  held, ledger closure 3.4e-3/1.1e-6, stepped-cavity mismatch 4.2e-12 → 3.4e-12; the
  ThermalLimitUnderCfl refusal class deleted (premise gone). Station 5: all 9 members re-run
  (~155 s each on the new spine vs ~60 s — 3 rhs evals + implicit solid + audit); readouts
  within ~0.1% of the session-12 records (blind: F/Isp still OVERLAP, p_c/c*/C_F same coherent
  coarse-tier signature; calibrated: F, Isp, AND p_c still OVERLAP; expander self-regulation
  reproduced); certificate prose updated (the honest-scaffolding integrator paragraph RETIRED;
  audit line added). **Review wave (same session, Ben's pattern):** self-review (δ-form CG
  RHS + audit composition re-derived by hand) + two independent agents (SDC/audit algebra;
  assembly/ledger/refits) — **no confirmed correctness bugs on any reachable path**; sweep
  weights, roll timing, ledger telescoping, exchange single-owner pairing, determinism, and
  the bitwise-unchanged rate arithmetic all verified against COUP-2/COUP-3. Latent hazards
  hardened same-session: `ensure_ws` staleness rebuild; station-4 `coolant_joules` read the
  exterior ledger line and recorded −0.0 (the duct's coolant Robin fires on the r_outer
  DOMAIN-EDGE line — `ExchangeStepReport.applied_bc_j` added, fixture reads both);
  `rel_resid` NaN propagation (f64::max drops NaN — the NaN acceptance was dead code) +
  NaN-‖b‖ → COUPLING_RESIDUAL; `march_flow` sub-ulp-dt stall refusal; `build_wall_patches`
  N_θ > 1 refusal (S8 re-keys per θ); scratch==t_field refusal; Linear-mode ledger
  debug_assert; N_SDC_CORRECTIONS ≥ 1 compile guard. Certificates byte-identical through
  the fix wave (no physics touched). True test count: **142 Rust + 28 Python**; all gates
  green.
