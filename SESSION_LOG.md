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
- Session 15 (2026-08-19): **PLAN S3 — THE MISSING FORCES: gas-phase F_visc (compressible
  viscous stress + Fourier conduction + species diffusion) on the exact cylindrical metric,
  swirl included, as the gas occupant of COUP-3's class D; suppressed at wall-law faces;
  certificates byte-identical.**
  **The operator (`crucible-solvers::gas_diffusion`, SOLV-1 §3.1):** one flux-form assembly
  over (u_r, ω = u_θ/r, u_z, T, C) — transport (μ, Pr→k, Sc→ρD, c_p/c_v) is pure config data
  derived from THE one owner (`WallLaw`'s constant set; new accessors — nothing restated).
  Per-component **symmetric two-point implicit cores** solved by the same fixed-structure
  Jacobi-CG as the solid class-D (shared `EPS_CG_RESID`/`N_CG_ITERS_MAX`, δ-form warm start,
  module-owned buffers): u_r gets the 4/3-μ radial core + the negative-definite −(4/3)μu_r/r̄²
  geometric diagonal; **u_θ is solved as angular velocity ω in the angular-momentum form** —
  the whole τ_rθ/τ_θz operator collapses to a pure symmetric diffusion with r³-class face
  weights, rigid rotation is DISCRETELY stress-free, angular momentum telescopes exactly, and
  the linear θ-momentum it induces is ledgered as an applied source (the flow operator's swirl
  pattern); **T is solved in total-energy flux form** (k∇T + τ·u work fluxes; dissipation
  emerges from the KE bookkeeping in the RHS — exact for the gamma-law class); C is constant-ρD
  Fickian. The cross-stress couplings (τ_rz cross-derivatives + the −⅔μ∇·u compressible
  corrections; ω has NO lagged remainder) converge by the SAME fixed Picard sweeps as the
  Robin-Robin exchange, riding one loop; **COUP-3 0.4.2 landed with the code**: the gas
  occupant named in §3.1, and `EPS_GAS_DIFF_RESID` (0.25) declared a **contraction guard**
  (the lagged remainder's structural gain is ≲ 1/12 at any Δt by AM-GM over the implicit
  diagonals; the truncated Picard's remainder is a temporal-truncation term of the same order
  class as the truncated SDC sweeps — accuracy is owned by the order gates), measured as the
  rate-staleness STATE effect per step against per-component conserved scales (momentum-joint —
  a near-zero component's own rate scale would read fp noise as divergence; found live).
  **Ownership at walls (SOLV-1 §3.5/COUP-2 §3.5):** resolved diffusion flows ONLY through
  gas↔gas faces (aperture-weighted); gas↔solid and gas↔exterior faces contribute nothing —
  unit-proven (gas rates bitwise independent of solid-cell operand garbage). Declared viscous
  BCs (data, COUP-7-closure style): NoSlip with wall-velocity schedules (a moving wall does
  ledgered work — the Couette drive), FreeSlip (τ·n̂ = 0), **Continuative** (zero-normal-
  gradient open plane — one-sided tangential stress; a FreeSlip channel end would truncate the
  real τ_rz and drive edge vortices), Isothermal/Adiabatic, species ZeroFlux (non-catalytic;
  Prescribed exists for MMS only). **Retired:** the S2 seam refusal ("a gas-domain diffusion
  class alongside the flow class is the S3 viscous wave") — `GasDiffusionClass` is the real
  thing; scalar gas-domain conduction beside flow stays refused (superseded).
  **S3 mini-sims (`solv3_gas_diffusion.rs`, 8 tests + 4 module unit tests):** annular
  Poiseuille marched at **33.5× the explicit viscous bound** (Δt = the gas CFL alone), exact
  cylindrical profile to 0.73%; Taylor-Couette swirl to 0.80% (+ rigid-rotation/uniform-state
  exactness and every face coefficient checked by hand at unit level); **recovery Couette
  EXACT analytic replaces the plan's flat-plate mini-sim** (same physics balance —
  dissipation vs conduction vs moving-wall work; a Blasius march is not a laptop mini-sim;
  plan §1 improvisation rule): adiabatic-wall recovery **3.584 K vs 3.581 K analytic at the
  cell** (ΔT_rec = Pr·U²/2c_p = 3.586 K, worst profile err 0.07% of ΔT_rec); thermal_bl erfc
  layer 4.9% (the fixture's declared isobaric-limit class) with the species layer at Sc ≠ Pr
  to 0.29%; **whole-operator MMS with every viscous/conductive/species term active: all six
  components at order 1.92–2.21** (the analytic residual assembled from the continuous stress
  formulas — an independent formulation; `euler_mms` gained public second-derivative bundles);
  a **four-class march** (flow + gas diffusion + solid conduction + exchange — the S4
  configuration in miniature) with the combined energy row closing every step and the liner
  warming through the one wall law; refusals (gas without flow; N_θ > 1 → plan S8);
  1-vs-4-thread bit identity through the full gas-scheduled step.
  **Findings (recorded):** impulsive wall/drive starts at ~50× stiffness ring the truncated
  trapezoid sweeps (bounded 7/8-damped oscillation — harmless on the unconstrained solid,
  positivity-fatal on gas) — fixtures use declared ramp schedules, the COUP-7 discipline;
  S4's RL10 already carries the injector ramp. A metal-like fixture μ gives the gas a
  metal-like k through μc_p/Pr and a flash-swinging wall-law h — the wall-heat registry's
  μ ≤ 1e-2 cap is confirmed physics, not caution.
  **Deferrals (owners named, module header):** wall-function skin-friction momentum debit
  (COUP-2 §3.1.2 mount-reaction ledger, verdict wave); species-enthalpy diffusion flux
  Σh_k·j_k + TableEos-consistent T refresh in the Picard (S4 spine — constant-c_v is exact
  for the gamma-law class); COUP-8 registry row + config grammar (S4 engine wiring); the gas
  assembly is serial like the solid one (perf — GPU wave); θ-diffusion fluxes + per-θ operand
  keying (S8).
  **Certificates: byte-identical through gate 5** — gas diffusion is opt-in config (Rule 13);
  no station schedules it (S4's ◆C1 turns it on with real transport tables); the
  no-gas-diffusion arithmetic is untouched. True test count: **154 Rust + 28 Python**; all
  gates green.
  **Review wave (same session, Ben's pattern):** two independent agents — SDC/audit-algebra
  lens and continuum-to-discrete lens — on the committed S3 diff. **No confirmed correctness
  bugs on any reachable path.** The physics reviewer built an INDEPENDENT continuum oracle
  (own analytic field, own algebra from τ = 2μe − ⅔μΔ) and measured the discrete rates
  against it: D_mr 1.99/2.00, D_mt 1.98/1.99, D_mz 1.97/1.99, D_en 1.99/2.00, D_rc 2.00/2.00;
  then finite-differenced the true Jacobian of `assemble_rates` and confirmed `apply_linear`
  reproduces it to round-off (≤ 3.4e-17 on velocities, 8e-14 relative on T) — the CG solves
  exactly the operator the composition applies, boundary arms and geometric diagonal included.
  The angular-momentum reduction, the z-face r̄ weight, the ρr̄²κV mass, the energy-flux slot
  mapping, `mms_visc_residual` term by term, and the drag/work signs at all four edges (6/6
  configurations) were each re-derived and confirmed. The SDC reviewer hand-traced all three
  sweeps' buffer generations (d0/dprev/dlag/dcur), the audit's weight correspondence, the
  T-solve KE bookkeeping, and the ΣκV·d[k] = port_net + src_net per-assembly identity
  (θ-momentum deliberately non-telescoping, ledgered wholly as a source).
  **The one real gap — a TEST hole, fixed:** the S3 battery could not see the compressible
  (dilatation) terms. `MMS_AMP[1] == MMS_AMP[3]` and one shared mode make manufactured u_r and
  u_z have IDENTICAL gradient fields, and `face_coefficients_by_hand_at_one_cell` drives one
  velocity at a time (the other's gradients identically zero) — so swapping `duz_dz` for
  `dur_dz` in the −⅔μ∇·u corrections was invisible. **Verified by planting exactly that
  mutation: all 12 tests stayed green.** New unit test
  `dilatation_and_cross_shear_discriminate_independent_gradients` — a field with BILINEAR
  cross terms making all four lag gradients pairwise distinct AND varying across the stencil.
  The bilinearity is load-bearing, and finding out why was itself a result: a spatially
  UNIFORM dilatation error cancels exactly out of the r-momentum (the τ_rr face term carries
  it with weight (A_out−A_in)/V and the −τ_θθ/r source with `geo` — the same number), which is
  correct physics (a uniform isotropic stress exerts no net force) and is why a simpler
  independent-shape field was ALSO blind. Mutation-tested: all four gradient-confusion classes
  now caught, pristine code passes.
  **Latent hazards hardened same-session:** T ≤ 0 now refuses (was `is_finite`-only — a
  non-positive temperature would have flowed into k∇T unremarked, META-1 P6); the gas ledger's
  `port_abs` no longer double-counts boundary ports (it added `tot` — which already contains
  them — on top of the per-port lines, silently LOOSENING `TOL_AUDIT` on every edge-touching
  cell); `rate_resid`'s scale loop NaN-checks explicitly (`f64::max` drops NaN — the S2
  `rel_resid` finding class); `compose_gas`'s masked-slot write now debug-asserts the
  zero-fill contract it depends on (κ = 0 would hide any drift from both the audit and the
  reductions); `ensure_gb`'s brick-count-only staleness test documents the full-rewrite
  invariant that makes it sound; `EPS_GAS_DIFF_RESID` records its dependence on
  N_ROBIN_SWEEPS ≥ 2; `validate` no longer indexes `brick(0)` unconditionally; the
  zero-gradient `(None, None)` case is named honestly (correct for the symmetric/quasi-1-D
  cases it reaches today; a genuinely under-resolved one-cell gas island belongs to the FND-3
  refinement wave). Accepted and recorded, not cured: boundary-cell local truncation is O(1)
  and mesh-independent (the expected half-cell Dirichlet closure — the solution order survives
  by supraconvergence, measured; a sign error would instead grow as 1/h, and does not).
  Certificates byte-identical through the fix wave. True test count: **155 Rust + 28 Python**.
- Session 16 (2026-08-20): **PLAN S4 — REAL PROPERTIES: the FND-7 constitutive spine's
  chemical-regime transport slot, filled; the wall law and F_visc read ONE owner; the liner's
  thermal mass made physical; ◆C1.**

  **The doc wave first (the working rule).** FND-7 0.5 gains the transport slot's
  **chemical-regime stage**: where ⟨Z⟩ ≡ 0 and Lee-More-Desjarlais/Stanton-Murillo degenerate,
  the everywhere-defined analytic backbone is **mixture-averaged Chapman-Enskog kinetic
  theory** — the chemical sibling of Saha/QEOS, *a physics model, not a data table*, so §3.6's
  no-cliff guarantee and §5's envelope-is-the-model's rule carry unchanged; declared band
  **10–20%**. It states, once, the quantity set the slot returns (μ, k, c_p, **c_v**,
  **∂h/∂Z|_{p,T}**), that **Pr = μc_p/k is derived, never a second datum**, and that the
  runtime coordinate is the **local state (p, h, Z)** — the S22 rule — because keying transport
  on T would make every per-cell query a chained `T(p,h,Z) → μ(T,p,Z)` interpolation,
  compounding two `interp_error_bound`s into a quantity neither describes. OFFL-5 0.3 §3.1a
  owns the emitted surface (and records that in this regime the §3.3 GP discrepancy is
  **identity** — the backbone *is* the evaluation; a fabricated correction would be worse than
  none). OFFL-3 0.5 widens the S23 transport feed to the caloric companions. SOLV-1 0.4.1 names
  `F_visc`'s third energy limb and makes the wall law's operands spine queries. COUP-3 0.4.3
  gives class `D` its **variable-coefficient rule**. META-3 0.8 §6.10 adds the four keys.

  **The one flux, decomposed (the design decision the session turned on).** A dissociated gas
  moves far more energy down a temperature gradient by recombination than by collisions. The
  tempting move is to tabulate an equilibrium conductivity — but that would be a second model
  of a flux the operator already carries. On the equilibrium manifold `Y_k(T,p,Z)` the chain
  rule is exact:
  `Σ_k h_k j_k = −ρD[(c_p,eq − c_p,fr)∇T + (∂h/∂Z)|_{p,T}∇Z]`, so its ∇T limb folds into the
  Fourier flux as `k_eff = k_fr + ρD(c_p,eq − c_p,fr)` — precisely the classical equilibrium
  conductivity — and its ∇Z limb is the resolved species-enthalpy flux. **One flux, two limbs,
  never counted twice**; the table therefore ships the *molecular* pieces and the *caloric*
  pieces and the runtime composes them with the one declared Schmidt number. The payoff is
  measurable: at the dissociated low-pressure corner `k_eff/k_fr` exceeds 3, and **Pr stays a
  gas Prandtl number (0.2–1.5) at both ends** — had k been left frozen while c_p went
  equilibrium, Pr would have blown up with the dissociation and the Colburn analogy would have
  silently misfired.

  **The offline product (`crucible_offl::transport`, OFFL-5 §3.1a).** Two engines, one state:
  CEA (the pinned equilibrium engine — the *same* solver, `gas_only` mode, and (p,h,Z)
  coordinate as the shipped EOS surface) supplies the caloric columns; **Cantera evaluates
  mixture-averaged Chapman-Enskog transport on that composition at that (T,p)** and is never
  asked to equilibrate anything. Using a second thermochemistry would have put a seam between
  the EOS the runtime projects onto and the heat capacity it linearizes with. Six columns
  (`viscosity`, `conductivity_frozen`, `cp_frozen`, `cp_equilibrium`, `cv_equilibrium`,
  `dh_dz`), five of them log-valued with rule-space bounds. Two refusals armed at generation:
  a CEA species absent from the Cantera set above 1e-6 mole fraction (O₃ appears at ≤ 1e-8),
  and a **two-fit c_p cross-check** (Cantera NASA-7 vs CEA NASA-9 on the identical state, 2%)
  — which is what *bounds* Cantera's above-3500 K extrapolation rather than hoping about it.
  `∂h/∂Z|_{p,T}` is a central difference over two CEA **TP** solves, step-independent to
  ~8e-5 relative over δZ ∈ [5e-5, 2e-4].

  **Grid density is a measured choice**: 57 × 37 × 7 = 14,763 nodes spanning **exactly the
  equilibrium surface's declared envelope** (so every node is a state a run may legally reach,
  and the two surfaces refuse and accept on the same set — the engine cross-checks that at
  assembly, because COUP-8 §3.3(2) checks each table against its consumers alone and cannot see
  a gap *between* two tables). Measured rule-space bounds: μ 3.0%, k 3.3%, c_p,eq 1.6%,
  c_v,eq 1.7%, c_p,fr 0.11% — a factor ≥ 3 inside the declared 10–20% physics band, which is
  what FND-5 §3.4's grid-sizing rule asks for. Artifact `tables/spine/
  lox_lh2_transport_v0.1.0.h5` (727 KB) + sidecar; regen ≈ 10 s.

  **The runtime seam (`crucible-solvers::transport`, FND-7 §3.3).** ONE provider, two occupants
  selected by config id (the `flow` ↔ `flow_shifting` pattern): `transport_constant` (the
  declared c_p/μ/Pr/γ/Sc set the stations have used since session 7 — still the right occupant
  for every analytic fixture, where a constant-coefficient exact solution is the point) and
  `transport_table`. **`WallLaw` is now stateless** — its private (c_p, μ, Pr) *was* the
  degenerate spine occupant, and it is now that occupant, stated once; transport arrives as a
  `TransportProps` operand at the near-wall cell's own state, the same provider the resolved
  `F_visc` next door reads. `gas_diffusion` likewise holds no transport: the caller refreshes a
  per-cell `GasTransportField` from the spine **once per Picard iterate, from that iterate's lag
  state** (COUP-3 0.4.3 — never inside the CG, which must stay the solve of one fixed linear
  operator). Interior face coefficients are the **arithmetic two-cell average** — exact in the
  constant-coefficient limit (which is what keeps every fixture and certificate bit-identical)
  and commutative, so the coefficient seen from either side is bit-identical and the CG stays
  SPD while the flux telescopes.

  **What the constant occupant preserves, bitwise:** `k = μc_p/Pr`, `c_v = c_p/γ`, `ρD = μ/Sc`,
  `∂h/∂Z ≡ 0` — the retired expressions verbatim, asserted by test, which is why stations 1–5
  regenerate byte-identically.

  **The wall-law band is now named.** The p-box corners used to be reached by scaling the wall
  law's private `cp_j_per_kg_k` at fixed Pr and μ — a proxy that happened to scale h
  proportionally but also moved the recovery temperature and (after S3) the resolved viscous
  fluxes. `wall_heat` gains a declared `band_factor` (default 1.0, range [0.5, 2.0]) that
  multiplies **h and nothing else** — the direct realization of SOLV-1 §3.5's ±20–30% band, and
  a labeled band coordinate rather than a tuning dial.

  **Liner thermal mass made physical (plan ruling #4).** The ~10³-fast ρc_p continuation device
  is deleted. Below the dial at which the liner is resolved at physical thickness the modeled
  ring is thicker than the metal it stands for, so it carries **two** declared homogenizations,
  not one: κ resistance-preserving (since session 11) and now ρc_p **capacitance-preserving** —
  `ρc_p_metal · t_real/t_model`, because what sets the expander bootstrap clock is thermal mass
  *per unit wall area*. 0.33 mm of SS-347 (ρc_p ≈ 3.95e6, META-3 `ss347-liner-thermal`) on the
  20 mm coarse ring ⇒ 6.52e4 J/(m³·K) (65× the retired device); on the 8 mm dial-16 ring ⇒
  1.63e5. Both retire together when the liner is resolved (THE RUN). The silver braze girdle is
  excluded — declared, a faster wall clock.

  **New tests, chosen for what nothing else could see.** Every S3 fixture runs a UNIFORM
  transport field, so a face coefficient that read only the visiting cell would have been
  invisible to all of them *and* would have broken the CG's symmetry silently.
  `face_coefficients_are_the_two_cell_average_of_a_varying_spine` gives the spine a
  two-directionally varying reading and checks the assembly, `apply_linear`, and the
  either-side flux antisymmetry by hand — **mutation-proven**: swapping `FaceTr::between` for
  `FaceTr::at` leaves the entire S3 battery (8 tests, 91 s) green and fails only this one.
  `species_enthalpy_flux_carries_energy_down_a_composition_gradient` does the same for the new
  term (zero on a single-composition gas — which is what made the S3 deferral honest);
  `temperature_solve_mass_uses_the_per_cell_spine_slope` for the c_v mass.
  `crates/solvers/tests/fnd7_transport_spine.rs` (7) scores the production artifact: units-gate
  refusals per column, the derived-group identities, envelope refusal instead of extrapolation,
  the interpolation bounds against the physics band, the equilibrium-conductivity behaviour
  above, and — the session's motivating claim as a test — that the tabulated occupant genuinely
  **departs from the constants it replaces** where it should (within ~2× at the chamber, where
  those constants were sized; far below at the recombined fringe; c_v spanning > 3× between the
  fringe and the dissociated core, which is exactly the deferral the S3 header recorded).
  `engine_smoke` gains the S4 configuration end to end and a proof that a surface without the
  transport columns **refuses at bind**, not somewhere downstream.

  **Deferrals (owners named).** Per-species diffusion coefficients (the single effective Sc is
  a declared one-composition-coordinate closure) and Soret/Dufour/pressure diffusion ride plan
  S5's species-vector state; the wall law's film/Eckert reference-temperature refinement is
  declared inside its existing ±20–30% band; measured high-temperature H₂O/H₂ transport data
  enter later as a coverage-weighted GP on this backbone with no runtime-contract change
  (OFFL-5 §3.1a); `interp_error_bound_log` still rides outside digest v3 (digest v4);
  θ-diffusion at N_θ > 1 refuses (plan S8); the gas assembly and the per-cell spine query are
  serial (perf — the GPU wave), and the spine query is now a measurable share of step cost.

  **◆C1 — MET (open mode, coarse tier).** `configs/rl10_coarse.toml` with the tabulated spine
  and F_visc scheduled marched **15 819 steps to a settled readout with no halt, no refusal and
  no schedule tuning** — the COUP-2 audit armed on every one of them, zero violations. Numbers
  below are the **post-review** march (the review's E1/E2 both changed the wall term, so the
  first measurement was re-run rather than published):

  | | S4 (full diffusion + real transport) | S3 spine (constant transport) | Δ |
  |---|---|---|---|
  | thrust | 74 382.2 N | 74 162.9 N | +0.30% |
  | Isp | 447.37 s | 446.25 s | +0.25% |
  | c\* | 2255.5 m/s | 2257.0 m/s | −0.07% |
  | C_F | 1.9451 | 1.9389 | +0.32% |
  | p_c | 3.0926 MPa | 3.0933 MPa | −0.02% |
  | **jacket heat** | **7.812 MW** | **9.737 MW** | **−19.8%** |
  | liner T_max | 410.1 K | — | — |
  | ṁ_inj | 16.947 kg/s | 16.947 kg/s | 0 |

  **The finding: the missing forces barely move plane-integrated performance at this tier, and
  move the wall term by two orders of magnitude more.** Every scored quantity lands within
  ±0.32% — the wall law already owned the wall and the boundary layers are not resolved at
  dial 5, so resolved viscous stress and conduction have little to add to a plane integral. What
  changes is the wall: real transport, and the corrected driving potential, drop the jacket
  pickup by a fifth and the liner's peak temperature from 535.7 K to 410.1 K.

  **The pre-review measurement was a trap, and the reviewer called it in advance.** Before the
  E1 fix the same march reported jacket heat −3.27%, which read as "the new physics is harmless
  here". The physics reviewer flagged that reading as *"a coarse-tier coincidence — most
  near-wall cells sit near ~2200–2700 K where the c_p rise and the μ/k changes happen to
  cancel"*, and predicted the sign would move at a resolved tier. Fixing the potential turned
  −3.27% into **−19.8%**: the cancellation was real and it was hiding a 1.9× error. The lesson is
  recorded because it is about method, not about this number — a small measured delta after a
  large model change is evidence to *investigate*, not evidence of harmlessness.

  **What did NOT settle: the liner.** Steadiness residual 7.06e-3 (better than the pre-fix run's
  1.38e-2 — the extra Robin sweeps converge the coupling harder) but liner T_max is still
  climbing at cutoff. With physical areal capacitance the wall's thermal time constant is
  ~C/(h_gas + h_cool) ≈ 37 ms against an 11 ms march — the gas field is settled and the wall is
  not. This is plan ruling #4 doing exactly what it was for (the retired device made the solid
  clock ~10³ fast, which is what hid it), and the cure is run length, not tuning: the certified
  budget must grow past the wall clock. **Cost measured:** 1625 s solver wall clock at dial 5 vs
  ~155 s on the S2 spine — **~10.5×**, of which roughly 7× is F_visc + the per-cell spine query
  and the rest is the review's Robin-sweep fix (3 → 5). The spine is queried for every gas cell
  on every Picard iterate of every SDC sweep; that and the serial gas assembly are the recorded
  perf deferrals the GPU wave owns.

  **New tests, chosen for what nothing else could see.** Every S3 fixture runs a UNIFORM
  transport field, so a face coefficient that read only the visiting cell would have been
  invisible to all of them *and* would have broken the CG's symmetry silently.
  `face_coefficients_are_the_two_cell_average_of_a_varying_spine` gives the spine a
  two-directionally varying reading and checks the assembly, `apply_linear`, and the
  either-side flux antisymmetry by hand — **mutation-proven**: swapping `FaceTr::between` for
  `FaceTr::at` leaves the entire S3 battery (8 tests, 91 s) green and fails only this one.
  `species_enthalpy_flux_carries_energy_down_a_composition_gradient` does the same for the new
  term (zero on a single-composition gas — which is what made the S3 deferral honest);
  `temperature_solve_mass_uses_the_per_cell_spine_slope` for the c_v mass.
  `crates/solvers/tests/fnd7_transport_spine.rs` (7) scores the production artifact: units-gate
  refusals per column, the derived-group identities, envelope refusal instead of extrapolation,
  the interpolation bounds against the physics band, the equilibrium-conductivity behaviour
  above, and — the session's motivating claim as a test — that the tabulated occupant genuinely
  **departs from the constants it replaces** where it should (within ~2× at the chamber, where
  those constants were sized; far below at the recombined fringe; c_v spanning > 3× between the
  fringe and the dissociated core, which is exactly the deferral the S3 header recorded).
  `engine_smoke` gains the S4 configuration end to end and a proof that a surface without the
  transport columns **refuses at bind**, not somewhere downstream.

  **Deferrals (owners named).** Per-species diffusion coefficients (the single effective Sc is
  a declared one-composition-coordinate closure) and Soret/Dufour/pressure diffusion ride plan
  S5's species-vector state; the wall law's film/Eckert reference-temperature refinement is
  declared inside its existing ±20–30% band; measured high-temperature H₂O/H₂ transport data
  enter later as a coverage-weighted GP on this backbone with no runtime-contract change
  (OFFL-5 §3.1a); `interp_error_bound_log` still rides outside digest v3 (digest v4);
  θ-diffusion at N_θ > 1 refuses (plan S8); the gas assembly and the per-cell spine query are
  serial (perf — the GPU wave), and the spine query is now a measurable share of step cost.

  **◆C1 — MET (open mode, coarse tier).** `configs/rl10_coarse.toml` with the tabulated spine
  and F_visc scheduled marched **15,778 steps to a settled readout with no halt, no refusal and
  no schedule tuning** — the COUP-2 audit armed on every one of them, zero violations. Measured
  against the S3-spine recorded baseline (`R_COARSE_ETA1`):

  | | S4 (full diffusion + real transport) | S3 spine (constant transport) | Δ |
  |---|---|---|---|
  | thrust | 74 080.8 N | 74 162.9 N | −0.11% |
  | Isp | 445.72 s | 446.25 s | −0.12% |
  | c\* | 2248.7 m/s | 2257.0 m/s | −0.37% |
  | C_F | 1.9437 | 1.9389 | +0.25% |
  | p_c | 3.0822 MPa | 3.0933 MPa | −0.36% |
  | jacket heat | 9.419 MW | 9.737 MW | **−3.27%** |
  | ṁ_inj | 16.947 kg/s | 16.947 kg/s | 0 |

  **The finding: the missing forces barely move integrated performance at this tier, and move
  the wall heat by ten times as much.** That is the physically right shape of the answer — the
  wall law already owned the wall, the boundary layers are not resolved at dial 5, so resolved
  viscous stress and conduction have little to add to a plane-integrated thrust; what changes is
  the *wall* term, where real (equilibrium) transport replaced a chamber-fitted constant. It is
  also the honest reading of the previous tier: the S3 numbers were not wrong because the
  physics was missing, they were right *for the quantities they reported* and blind to the one
  the expander cycle depends on.

  **What did NOT settle: the liner.** Steadiness residual 1.38e-2 vs the S3 spine's 6.16e-3, and
  liner T_max 535.7 K still climbing at cutoff. With physical areal capacitance the wall's
  thermal time constant is ~C/(h_gas + h_cool) ≈ 1.3e3/3.5e4 ≈ 37 ms against an 11 ms march —
  so the gas field is settled and the wall is not. This is plan ruling #4 doing exactly what it
  was for (the retired device made the solid clock ~10³ fast, which is what hid this), and the
  cure is run length, not tuning: the certified budget must grow past the wall clock. **Cost
  measured:** 1075.6 s solver wall clock at dial 5 vs ~155 s on the S2 spine — ~7× (this run
  shared the machine with the test battery and the review agents; the clean figure is lower).
  The per-cell spine query (six interpolations per gas cell per Picard iterate) is a real share
  of that and is a recorded perf deferral alongside the serial gas assembly.

  **Gates green.** True test count: **172 Rust + 34 Python** (155 + 28 at S3). Certificates:
  stations 1, 2, 4 and the convergence certificate regenerate **byte-identically** — the
  constant occupant reproduces the retired expressions bitwise, which is the whole point of
  keeping it. Stations 3 and 5 carry text-only amendments: station 3's transport-feed deferral
  is discharged, and station 5 now declares **S2/S3 as its spine of record** with the measured
  S4 delta and the reason the eight members are not re-scored here (~5 h of laptop march, and
  the COUP-5 ensemble wave replaces these hand-run corner brackets anyway — plan S18/S19).

  **Review wave (same session, two independent agents — Ben's pattern).** One on the discrete
  algebra (conservation, the SDC/Picard structure, the CG operator identity, determinism), one
  as an independent continuum/thermodynamics oracle. Both worked from first principles: the
  algebra reviewer built the full dense Jacobian column-by-column from `apply_linear` and
  finite-differenced `assemble_rates` against it (worst deviation **1.0e-10** across all five
  components on a 100×-contrast spine, boundary arms and geometric diagonal included — the S3
  identity survives variable coefficients), asserted `L[i][j] == L[j][i]` **bitwise**, proved
  masked transport slots unreadable (two fillings differing by 1.0 vs 9.87e11 give bit-identical
  rates), and confirmed the constant occupant's bit identity down to the sign of a zero. The
  physics reviewer re-derived the chain-rule decomposition against CEA + Cantera partial-molar
  enthalpies (`Σ h_k ∂Y_k/∂T` vs `c_p,eq − c_p,fr`: 0.19–0.32%; `∂h/∂Z` two ways: 0.05–0.07%)
  and measured the >3500 K extrapolation directly by rebuilding the mechanism with stretched fit
  ranges (**μ ≤ 0.29%, k ≤ 2.87%**).

  **Four confirmed findings, all fixed:**

  1. **The wall law's driving potential.** The Colburn analogy transports *enthalpy*, so
     `h·ΔT` needs the **film-mean** slope `(h_aw − h_w)/(T_aw − T_w)`. S4 fed it the local
     equilibrium c_p — the *peak* of a strongly-peaked curve (~7970 vs a film mean of ~4130 at
     the chamber) — **overpredicting q_w by 1.7–2.4×**, one-signed, through a ±20–30% band. The
     sting: the pre-S4 constant c_p = 5000 was accidentally *inside* that band at 1.21×, so S4
     made the property more accurate and the flux less so. Fixed: the spine returns **both**
     c_p's, the convective limb drives on `cp_film` (frozen c_p, a measured 0.94–1.06 proxy at
     chamber/throat) and recovery on the local one; exact `h_aw − h_w` is a named deferral.
     **SOLV-1 0.4.2 withdraws** the 0.4.1 sentence claiming the gas-state-vs-reference-temperature
     choice sat inside the band — measured h(film)/h(gas) ≈ 0.24–0.40, a factor 2.5–4, now a
     declared model-form limit instead.
  2. **The Robin-Robin Picard lost its margin.** With `h` a config constant the exchange map was
     **affine** in T_gas and three sweeps converged it to round-off; per-cell transport makes it
     nonlinear. Constructed failing case: a four-class duct at a near-wall gradient of ~1.7e5 K/m
     (milder than an RL10 chamber wall) **halts** on `EPS_ROBIN_RESID` at 3.8e-6. The contraction
     survives — 4 sweeps give 5.5e-10, 5 give 1.6e-11 — so it is a margin loss, not a divergence.
     `N_ROBIN_SWEEPS` raised **3 → 5** (four would clear it; five because `dh/dT_gas` grows with
     dissociation and the shipped surface reaches `k_eff/k_fr` = 22.7), with the measurement in
     the constant's doc and the `EPS_ROBIN_RESID` rationale corrected — the acceptance itself is
     deliberately unchanged, since loosening it would hide the staleness it exists to catch.
  3. **Sc = 0.5 was 25% off, and META-3's justification was factually wrong.** Measured directly
     (impose a pure ∇Z at fixed (T,p), contract the mixture-averaged fluxes onto the elemental
     fuel fraction): **Sc_eff = 0.402–0.412 across four decades of pressure**. The per-species
     numbers are H 0.21, H₂ 0.27, O 0.68, OH 0.70, O₂ 1.02, **H₂O 1.79** — the 90%-mass species
     is the *slowest* diffuser, not the 0.8 the entry claimed. Nominal → **0.40**, which moves
     `k_eff` from 21–32% below the mixture-averaged truth to ~11%, inside the declared band.
  4. **A claimed guard that did not exist.** Both reviewers found it: the `cp_eq ≥ cp_fr` comment
     asserted the positivity loop would catch a violation, but that loop tests each column alone
     — 2.3% of shipped nodes carry a tiny negative Δ (CEA round-off, ≤ 7e-10, in undissociated
     corners). The data are left alone (doctoring a table to satisfy an inequality is worse than
     the 1e-13 W/(m·K) it would fix); an explicit `k ≥ k_fr(1 − K_ORDERING_SLACK)` floor now
     guards the real failure, and the **tabulated occupant gained the same declared rails the
     constant one's manifest carries** — a table is not more trustworthy than a config value.

  **The test hole they found, closed and mutation-proven.** The species-enthalpy term had no
  coverage the battery could feel: deleting the **boundary-arm** enthalpy port left every test
  green (the omission is conservation-neutral — it transports the wrong physics without breaking
  the ledger), and a one-sided `dh_dz` read at interior faces left all 27 lib tests green while
  breaking conservation. Two new tests, each verified to fail on exactly its own mutation and
  nothing else: `the_boundary_species_port_carries_its_enthalpy_too` (hand-computed half-cell
  flux + its ledger line) and `a_varying_spine_still_telescopes_exactly` (the
  `Σ rate·κV == ports + sources` identity on independently-varying μ/k/ρD/∂h∂Z fields — a shared
  scale factor would let one coefficient's one-sided read hide behind another's correct average).

  **Also recorded:** `h_offset` now shifts the transport interrogation as well as the EOS one —
  the consistent choice (one coordinate) but a widened reach for a calibration dial, declared in
  `table_eos.rs`. The ∇p limb's neglect was re-justified: it is **not** small (22–45% of the
  retained ∇T limb) but vanishes where diffusion matters (∂p/∂n ≈ 0 in a boundary layer) and is
  irrelevant where it doesn't (Pe ≫ 1) — falsifiable as stated, unlike "third-order small".
  Wall catalycity is now stated (equilibrium operands in the law, non-catalytic species wall in
  `F_visc`; they never share a face). Dead `BadCoefficient` refusal removed; station-4's
  transport bundle hoisted out of its per-call path. Perf deferral widened: the spine is queried
  for every gas cell on every Picard iterate of every SDC sweep — 9 full-field refreshes per step
  at 6 interpolations each — which is a real share of the ~7× cost.

  **Certificate consequence of the Robin-sweep fix — checked, not assumed.** I expected station 4
  to stay byte-identical (with the constant occupant `h` is independent of T_gas, so the map the
  sweeps relax is affine) and **it moved**. The reason is that the affine argument covers only
  the gas half: the Robin-Robin fixed point also relaxes the *solid* surface temperature the
  exchange debits against, so more sweeps genuinely converge it further. The shift is
  noise-level and in the right direction — ledger closure **improves** 3.36e-3 → 3.18e-3, the
  analytic series-resistance oracle moves 9.44e-4 → 9.49e-4 against its 5e-3 gate, mid-duct h
  2298.9 → 2298.8 W/(m²·K), gas enthalpy deficit 133295.0 → 133319.1 W. Regenerated and
  committed. Recorded because the prediction was wrong and gate 5 is what caught it: a
  correctness fix in a shared constant reaches every fixture that runs the coupled step, and the
  right response is to re-measure rather than to reason about which ones "shouldn't" care.


- Session 17 (2026-08-20): **PLAN S5 — THE COLD/UNBURNT BRANCH: the burn-progress
  c = 0 branch shipped as a (p, h, Z) surface `TableEos` binds with no new occupant; the
  frozen↔shifting bracket shipped as a declared model-form band; the species-vector FIELD
  widening split out to S5b.**

  **The scope decision, made first (Ben).** Plan S5's headline goal — "species-vector state
  widening" — is a `NCOMP` change to the flat, compile-time-sized conserved state
  (`Cons = [f64; NCOMP]`, ~125 references across `sdc`/`gas_diffusion`/`euler` + the engine
  consumers), and stable Rust blocks the clean const-generic because `NPRIM = NCOMP + 2` needs
  `generic_const_exprs` (nightly). It is bit-identity-critical for the shifting stations and has
  **no consumer before the COUP-5 ensemble wave (S18)**. The physics that de-risked the split:
  the session's *named* mission — the cold/unburnt branch and the fine-dial establishment cure —
  needs none of it. The unburnt-reactant branch is a **(p, h, Z)** surface (a two-stream reactant
  mixture is set by Z), and S6's ignition blend is unburnt(p,h,Z)↔burnt(p,h,Z); the full
  frozen-species *field advection* is the orthogonal frozen end of the frozen↔shifting bracket,
  and the bracket carries that delivered-performance content as a declared band without advecting
  a species vector. So the widening is both the riskiest and the least urgent of the five goals.
  Ben chose to split it into **S5b** (may ride S8's state-layout/refluxing wave); per-species
  diffusion coefficients + Soret/Dufour + Stefan-Maxwell baro-diffusion ride S5b with it
  (per-species D is meaningless without per-species gradients).

  **The doc wave first (the working rule).** OFFL-3 0.6 promotes the **unburnt-reactant surface**
  from a contract row to shipped and the **frozen↔shifting bracket** from a validation check to a
  shipped declared band; FND-7 0.5.2 retargets the frozen-chemistry-ceiling + per-species-D owners
  to S5b (with the S5 bracket recorded as the interim delivered-performance mitigation); SOLV-1
  0.4.3 states the `{ρX_k}` widening rides S5b and the c=0 branch needs none of it; SOLV-4 0.4.1
  records that its §3.6 blend's `h_u`/`T_u` surface is now a real artifact; META-3 0.8.1 fixes the
  `schmidt-combustion-gas` per-species-D handoff (S5 → S5b); plan §5/§8 record the split and the
  two flagged rulings.

  **The offline product (`FrozenReactantEngine`, OFFL-3 §3.3).** The unburnt branch is the
  **gas-phase ideal-gas frozen reactant mixture** — gaseous H₂ + O₂ at the mass proportions Z
  sets — evaluated by the **same CEA reactant-`Mixture` machinery** the equilibrium surface's
  reactants use (the `calc_property(cea.ENTHALPY, …)` mechanism, on a *gasified* reactant mixture
  rather than the liquid one `injection_enthalpy` reads — same absolute formation reference, different
  species state), so the two branches carry **one enthalpy reference**. That is a correctness
  requirement, not tidiness: the SOLV-4
  §3.6 blend `h = (1−c)h_u + c·h_b` is a category error on two references. `T(h, Z)` is the
  (pressure-independent) monotone inverse of the frozen mixture enthalpy by fixed-count bisection;
  `ρ = pM̄/(R̄T)`, `c_p,fr = ∂h/∂T`, `c_v = c_p,fr − R̄/M̄`, `γ_fr = c_p/c_v`,
  `a = √(γ_fr R̄T/M̄)` — CEA's own sound-speed algebra. The mixture is ideal-gas exactly as CEA's
  products are, so the two branches are *consistently* ideal, differing only in composition. Its
  **envelope floor is 100 K** — where the reactant elements (zero formation enthalpy) still give a
  physical γ ≈ 1.47 while the equilibrium surface's condensing products already refuse, the *point*
  of the branch. CEA converges below the floor (to ~35 K), which is what lets the grid overhang for
  interpolation, but convergence is not validity: the review corrected an early "~50 K valid" phrasing —
  below ~60 K the extrapolated NASA thermo degrades (γ → 1.14 at 40 K) and `state_php` refuses where
  `c_p ≤ 0`; those sub-floor nodes are the declared gated-off overhang. **One M̄ bug caught in the
  probe:** `of_ratio_to_weights` returns UN-normalized
  mass weights (sum ≈ 34), so M̄ = Σ(mass)/Σ(moles), never `1/Σ(moles)` — the first version gave
  c_v < 0 and refused loudly, exactly as it should.

  **The surface artifact.** `tables/chem/lox_lh2_unburnt_v0.1.0.h5` (670 KB) + sidecar, the same
  FND-5 schema as the equilibrium surface (density/sound_speed/temperature/gamma_eff/mbar keyed
  (p, h, Z)), regen ≈ 4 s. Grid 9 × 121 × 15 = 16 335 nodes is a measured choice: the p-axis is
  coarse on purpose (ideal-gas ρ is exactly log-linear in p and the other columns are p-invariant,
  so p resolution costs and buys nothing beyond bracketing the envelope), the h-axis carries the
  T(h) curvature and the Z-axis the M̄/mixture-c_p curvature (worst at the H₂-rich hot corner where
  H₂'s vibrational modes activate — found by locating the worst holdout point, not guessed).
  Measured holdout bounds: temperature 1.66 K abs (~0.08 % hot, ~1.7 % at the 100 K floor), density
  0.81 % rule-space, sound speed 0.90 m/s, γ 1.2e-3 — a factor ≥ 10 inside any reasonable model
  band. The rectangular h-envelope is Z-coupled (the frozen enthalpy is Z-dependent), so it is
  derived as `h_env_lo = max_Z h(t_floor, Z)`, `h_env_hi = min_Z h(t_ceil, Z)` — a hot floor at the
  H₂-poor edge, a cold floor at the H₂-rich edge, every node convergent. `t_floor = 100 K` is the
  declared cold floor of the gas-phase branch; below the liquefaction line the ideal-gas frozen
  mixture is a declared metastable model and the real two-phase state is SOLV-1's W4 drift-flux
  extension (plan S15).

  **The frozen↔shifting bracket shipped (`KineticEfficiencyBand`, OFFL-3 §3.2).** Promoted from a
  §6-3 validation check to a shipped datum a run records in its pedigree, so performance is reported
  as the `[frozen, shifting]` interval, not a point. `relative_gap` is the raw bracket width
  `(shift−frozen)/shift` — a few percent, widening with the area ratio (~3.7–4.8% at ε = 61 over
  MR 5.0–5.5, ≈ 4.3% at the RP-1311 anchor) as the frozen limb leaves more recombination energy
  unclaimed; the JANNAF kinetic-efficiency knockdown
  (~0.8–1 % of shifting Isp, `jannaf-eff`) is the data-anchored *delivered estimate inside* that
  bracket (H/O kinetics are fast, so the real engine hugs the equilibrium end). The ordering is
  guarded (frozen can never exceed shifting — an inverted bracket refuses).

  **The runtime seam — zero new production code (SOLV-1 §3.4).** The unburnt surface uses the
  equilibrium surface's FND-5 schema, so the **existing `TableEos` occupant binds it unchanged** —
  the whole runtime cost of the cold branch is "bind the same occupant to a different table," no
  `if(unburnt)`, no new state. `crates/solvers/tests/solv1_unburnt_branch.rs` (4) proves it against
  the production artifact through the full FND-5 pin gate: `TableEos` binds and projects a cold state
  (~168 K, frozen Γ₁ ≈ 1.4 — genuinely cold reactant gas, not the ~3000 K a burnt surface would
  read); uniform cold rest is a **bitwise fixed point** (well-balance is geometric, not
  EOS-specific); a cold/warm transient **marches through the projection every stage without
  refusing and telescopes to round-off** with elemental Z untouched (the establishment-cure core —
  a startup cell has an honest home); and a sub-floor state **refuses rather than extrapolates**
  (the gas-phase branch is a declared model down to its floor; liquid/vapor is S15).

  **What "retires the envelope refusals" actually means here — the honest scoping.** The cure is
  *enabled* by the shipped branch and *realized* at S6. S5 ships the branch that makes cold unburnt
  states representable and proves it is marchable; the *routing* of a transient cell to the unburnt
  branch is S6's burn-progress c-blend. The pinned v0.3.2 burnt surface and the ◆C1 config are
  untouched, so the shifting stations stay byte-identical — the log does **not** claim the RL10
  dial-16 establishment refusals are gone (they are cured when S6 wires the blend). A correction of
  the reviewer's own reflex: I first wrote a test asserting the equilibrium surface *refuses* at the
  cold h; it does not — at the same h the burnt branch gives a *hot* state (the two branches map
  h→T differently by the heat of reaction). The honest, robust claim the test now makes is that the
  branches are *physically distinct* (Σ several MJ/kg apart at the same p,T,Z), which is *why* the
  c=0 branch is needed.

  **Two flagged rulings (mine).** (1) **Settle budget** — grow the certified budget to cover the
  ~37 ms liner clock as *run length* (ruling #4) at the coarse tier; ◆C2 is S7, so the exact
  flow-through count is set when S7 needs it, not paid now. (2) **Station-5 re-score** — left to the
  COUP-5 ensemble wave (S18); S2/S3 stay spine-of-record with the S4 delta declared (the hand-run
  brackets are replaced there anyway). Neither costs a march this session.

  **Deferrals (owners named), all to S5b unless noted.** The frozen-mode `{ρX_k}` field advection
  (the `NCOMP` widening); the frozen-composition (p, h, {X_k}) advection surface; per-species
  diffusion coefficients on the transport spine (superseding the single Sc); Soret/Dufour + Stefan-
  Maxwell baro-diffusion. The liquid-injection / vaporization coupling of the unburnt branch is
  SOLV-1's W4 drift-flux extension (plan S15).

  **Gates green.** True test count: **179 Rust + 45 Python** (**175** + 34 at S4 — the CLAUDE.md "172"
  was a stale harness undercount, corrected here): +4 Rust (`solv1_unburnt_branch`) and +11 Python
  (9 unburnt-surface incl. the Cantera frozen-c_p cross-check and the in-envelope diatomic gate, 1
  shipped-bracket, 1 unburnt-path regen-determinism probe — the frozen-reactant build path earns its
  own cross-process digest match, since its bisection T-inversion is a new pipeline). Certificates
  regenerate **byte-identically** — no production Rust changed and no station-pinned table moved (the
  unburnt surface is a new artifact with no station consumer yet).

  **Review wave (two independent agents, before commit — Ben's pattern; one pointed at the doc
  amendments, per the S4 finding).** The code/physics reviewer found **no correctness bug** and
  verified every load-bearing claim empirically against CEA — the shared enthalpy reference
  (`h_u − h_b` = the physical heat of reaction, +12.4→13.5 MJ/kg across MR 4–6), M̄, units, cross-process
  determinism, and the envelope bounds (in-envelope error 0.67–0.84× the stamped bounds). Its hardening,
  all applied: an explicit fail-loud `c_p ≤ 0` guard in `state_php` (the bisection's monotonicity
  assumption, checked at the root — META-1 P6); a **Cantera NASA-7 cross-check** bounding the cold-corner
  CEA extrapolation (worst 2.8% at the 100 K floor, < 0.3% above 200 K — the OFFL-3 §6-2 two-fit
  discipline applied to the frozen branch, which had no built-in thermo cross-check); tightened Python
  physicality gates (γ > 1.1 all-nodes + an in-envelope γ ∈ [1.3, 1.5] / c_p gate — a units-collapse
  bug gives a self-consistent γ ≈ 1.0003 that a `> 1.0` gate misses); a tightened reference backstop.
  The doc reviewer caught the exact S4-class defect: the plan §5 line's **trailing "Retires the envelope
  refusals" sentence stood uncorrected** by the split parenthetical — read standalone it claimed a cure
  S5 only *enables* (§8 and this log had it right); fixed. It also corrected the test count (179, not
  the 172-baseline undercount) and flagged the ε = 61 band figure as condition-specific (~3.7–4.8% over
  MR 5.0–5.5, ≈ 4.3% at the RP-1311 anchor — qualified in all three places). No finding survived to the
  commit unfixed.

- Session 18 (2026-08-21 build / 2026-08-24 close — one commit): **plan S6, IGNITION.** The chemical
  regime can now light: a cell carries the burn-progress field c (`NCOMP` 6→7, the fixed `+1` `ρc` slot —
  inert by default, so every shifting station is untouched), its thermochemistry is the **energy-conserving
  flamelet blend** of the S5 unburnt and the equilibrium surfaces (`BurnBlendEos`: both branches read at
  the cell's own (p, h, Z) on S5's shared CEA reference — burning is the same conserved energy re-read as
  hotter gas, no explicit heat-release term; density = mass-weighted specific volume, pressure = one
  deterministic Illinois root; c = 1 recovers shifting mode bit-for-bit), and c evolves by the SOLV-4.4
  **bistable (Nagumo) reaction-diffusion law** (`Combustion`): matched `(D_c, K)` give a **pushed front** at
  exactly `S_T` with width Θ·Δ at every resolution — the FSD `|∇c|` form is analytically degenerate and
  monostable KPP is a pulled front (grid-pathological + noise self-ignites), both rejected during the build.
  Closures `S_L(p, T_u, Z)`/`τ_ign(p, T_u, Z)` are one offline Cantera surface (OFFL-3 §3.3, `h2o2.yaml`;
  extinction = the columns' own values). All three rate terms explicit at this tier behind the loud
  positivity guard; the stiff class-`R` implicit auto-ignition is the declared S7 hardening.

  **The close found the v1.5 build's tables ignition-incompatible — the S6-close envelope set (OFFL-3
  0.6.2; Ben ruling: expand the tables, never throttle the spark to fit them).** Two layers, both caught by
  `spark_box`/`quench_box` halting mid-march. (1) **Envelope inversion:** the burnt ceiling (+4.0e6 J/kg,
  sized for station transients) sat *below* the unburnt's (+5.0e6) — but the blend interrogates both
  branches at one enthalpy, so an igniting kernel (auto-ignition wants T_u ~ 1000–1300 K ⇒ h ~ 2.5–4e6)
  died on the *products* surface while still valid cold gas. (2) With that fixed, confined-ignition blast
  compression drove mid-transition cells to h ≈ 5.3e6 — past the unburnt ceiling while still burning.
  Cure: **burnt v0.4.0** to +1.225e7 J/kg (standing rule: burnt ≫ unburnt — a burning cell can never
  refuse where the same cold gas was fine), **unburnt v0.2.0** to t_ceil 2900 K, **ignition v0.2.0** T_u
  to ~3000 K under the new **envelope-consistency contract** (the ignition surface's T_u range must cover
  T_u at the unburnt h-ceiling, or a legal blend state becomes a closure refusal), and **transport
  v0.2.0** re-derived from the new equilibrium envelope — the engine's own armed transport↔EOS consistency
  refusal caught that dependency, and the Cantera↔CEA two-fit c_p cross-check got a **measured** hot-side
  tolerance (4% above 3500 K; measured max drift 2.44% at the 4276 K corner — equilibrium dissociation
  buffers the hot edge to ~4300 K; the 2% mapping-bug gate stands below). Every extension is **strict** —
  old nodes bit-exact, verified column-by-column at regeneration — so all five certificates are
  byte-identical except station 1 (gains only the new ρc MMS column, order 1.92–2.17 — the widening's
  passive-scalar proof made visible) and station 5 (provenance text). Transport bounds *tightened*
  (N_H 37→59, spacing ~6 % finer than 0.1.0 over the widened span).

  **Two code fixes the mini-sims forced.** (1) The blend's pure-limit threshold vs the reaction's burnt
  fixed point: the source zeroes (domain guard) at c ≥ 1−BURN_COMPLETE, so c **asymptotes from below and
  never crosses** (measured: pinned ~2e-7 under it) — the blend's old 1−1e-9 pure-burnt threshold left
  every burnt cell interrogating the unburnt branch at ~1e-3 weight forever, refusing on hot states that
  branch cannot describe. The cure is **asymmetric**
  (SOLV-4 0.4.3, sharpened by the review wave): `EPS_B_PURE_BURNT = 2·BURN_COMPLETE` (one owner; the
  pure-burnt region strictly contains the attractor) while `EPS_B_PURE_UNBURNT` stays at the original
  1e-9 — the reaction pins an attractor only at the burnt end, and each skip's crossing step scales
  with the *dropped branch's* specific volume: burnt-side `EPS·(v_u/v_b)` ≈ 2.6e-4 (dropping the dense
  branch — fine), cold-side `EPS·(v_b/v_u)` ≈ 7.8·EPS (dropping the light one), so a symmetric 2e-3
  cold threshold would have stepped density ~1.5 % — *outside* the unburnt density column's own 0.8 %
  bound (the review wave measured this; the first fix had it symmetric and mis-declared). The same
  domain guard added to the reacting-measure/consumption-rate diagnostics. (2) **The spark is a literal electrical energy deposit** (Ben ruling):
  its one cited datum is the deposited energy — `spark-igniter-class` **pinned** (META-3 0.8.3: H₂ MIE
  ≈ 0.017 mJ floor per Lewis & von Elbe; aerospace exciter class ~0.1–20 J/discharge; TM-107318: "the
  ignition source is an electric spark", ASI at the injector-face center) — delivered as a **bounded pulse
  ending ~at the ignition time** (battery: ~0.7–8 J in 20 µs; a sustained deposit into an already-burnt
  kernel superheats mid-transition cells past the metastable-reactant validity edge — the hot ASI-torch
  regime is S7 hardening territory). COUP-7 0.4.1: position/extent/window are config *placement*, not
  sourced claims.

  **Battery green (5/5, ~1 min release).** `flame_1d` — THE gate — S_c 5.738 (coarse) vs 5.625 (fine)
  m/s = 2.0% with the front a fixed Θ-cells wide on both; `spark_box` lights (peak R = 7.4e-2 ≫ floor,
  burned 0.999); `lean_no_light` refuses (peak R = 4.7e-12 — NEVER_IGNITED class); `ignition_delay`
  fires on the surface timescale (t_ign/τ = 0.12 — thermal runaway *shortens* the naive (1−c)/τ clock,
  as it should); and **`adiabatic_box_cannot_flame_out` pins the close's physics finding**: a lit closed
  adiabatic box *must* burn to completion (one-cell kernel, ~0.7 J ⇒ burned 0.999) because quenching is a
  heat-loss phenomenon — FLAMEOUT is unreachable without a loss channel, so the model refused to fake the
  originally-drafted `quench_box`. That test rides **S7** as conductive loss to cold isothermal walls via
  the blend↔class-D coupling S7's startup march builds regardless (the tube's 0.8 mm gap is at the
  quench-distance scale, though H₂/O₂'s own quenching distance is a few× smaller than H₂/air's ~0.6 mm,
  so the S7 box may need lower pressure or a narrower gap); its consumer — the COUP-4 FLAMEOUT verdict
  object — is S7 too.

  **Deferrals (owners named).** Config-grammar + `run.rs` igniter wiring → S7/◆C2 (no unconsumed
  manifests); stiff class-`R` implicit auto-ignition → S7; `quench_box`/FLAMEOUT demonstration → S7 (above);
  `turbulent-flame-speed` pin → S7 (first turbulent consumer; the battery is laminar, wrinkling = 1);
  near-vacuum tangency acceptance for the blend projection → S7 (the RL10-plume feature, as `TableEos`);
  **the cold-side/low-p closure-envelope guard → S7** (review-wave flag: `Combustion::accumulate` queries
  `S_L`/`τ_ign` on every non-burnt gas cell, and the ignition surface floors at T_u ≈ 230 K / p ≈ 6.8 kPa —
  a cryo-fill or near-vacuum RL10 cell will refuse-halt at S7 wiring unless the cold analogue of the
  BURN_COMPLETE domain guard, or wider closure floors, lands with it);
  N_θ > 1 combustion → S8. Session interrupted once by a machine restart (the first ignition-table regen
  died silently to a laptop sleep; the generator now streams per-row progress — a 30-minute silent
  pipeline is undiagnosable by design).

  **Gates green: 184 Rust + 50 Python** (+5 each side: the `solv4_combustion` battery; the
  `test_ignition_surface` VAL-2 anchor suite). Gate 5: certificates regenerate deterministically with
  exactly two intended diffs committed this session — station 1 (+ the ρc MMS column, every old value
  bit-identical) and station 5 (provenance text: the v0.4.0 strict-extension note). A two-agent review
  wave (code/physics + doc-claims) ran before commit. **Code/physics:** no behavioral correctness bug in
  the marched physics (audit-armed battery re-run independently; strict extensions re-verified
  column-by-column; igniter arithmetic reproduced in exact f64; ledger-closure argument checked) — but it
  caught the one real declaration bug: the first threshold fix was symmetric at 2e-3 and its declared
  ≤1e-3 step neglected the v_b/v_u ≈ 7.8 volume-ratio scaling of the cold-side crossing (measured 1.55 %
  at the flame state — outside the 0.81 % unburnt density bound). Fixed as the asymmetric-threshold
  design above. Its hardenings, applied: the one-cell adiabatic kernel moved off an exact cell-face ulp
  edge (half 0.5 → 0.6 cells); the regen-determinism ignition probe repointed to a measured dual-active
  (flammable AND auto-ignitive) corner box — 50–100 bar × 1000–1100 K, all 8 corners measured — because
  the dual-active region is a thin diagonal band in (p, T_u) that a naive probe box misses, and given
  n_tu_ext=0 so it no longer inherits the production grid's appended rows; generator progress moved to
  stderr (the probe's stdout is the digest contract); the stale transport-comment sentence handed the
  >3500 K claim to the hot constant; the S7 cold-side envelope trap recorded above. **Doc-claims:**
  verified the artifact numerology, strict extensions, envelope contract (2900.0 ≤ 2934.1 K measured),
  battery numbers, pins, TM-107318 quotes verbatim, and META-3 mechanism facts against the bundled file;
  confirmed four consistency defects (the owning-doc SOLV-4 0.4.3 row carrying the pre-fix threshold
  wording; the §3.6 integration bullet still describing the retired FSD `|∇c|` term as the propagation
  form; a stale ~1.6 J spark-energy floor; OFFL-3's certificate parenthetical contradicting the correct
  PLAN/CLAUDE statements) — all fixed, plus the minor wording risks (quench-distance phrasing, spacing
  claim, order range, wrinkling-in-fixtures). The final full gate-4 rerun then caught one more:
  the VAL-2 fresh-holdout test filtered activity at the *point* (truth + estimate both active) while the
  stored bound's declared domain is the builder's *fully-active-cell* criterion — on the extended surface
  a ⅜-offset point landed in a τ-cap-straddling cell (est 3.1e-3 vs truth 8.1e-6 at 6.75 kPa / 2894 K,
  the declared sharp feature, not covered interpolation) and correctly failed; the test now applies the
  builder's own 8-corner mask. No finding survived to the commit unfixed.

- Session 19 (2026-08-24): **plan S7, FIRST STARTUP VERIFICATION + ◆C2.** "DOESN'T START" became a
  computed outcome: the march now ends in a typed **COUP-4 verdict object** (WORKS / DOESN'T-WORK with
  mechanism, location, time, physical-vs-numerical diagnosis, and the criterion it was judged against),
  attached to every mid-march halt and to the WORKS report, written to `runs/<name>/verdict.txt`. The
  WORKS criterion is config grammar (`[operating_profile]`: commanded p_c/thrust arm it; `eps_works`
  default 0.02, `t_dwell_flowthroughs` default 20; `flowthroughs` IS the declared T_S1_HORIZON; horizon
  expiry without a completed dwell = the distinct FAILED_TO_REACH halt naming the offending quantities).
  NEVER_IGNITED/FLAMEOUT consume the reacting measure R against the **ṁ-scaled floor**
  `max(EPS_IGNITED, 1e-4 × delivered ṁ)` (the mini-sim absolute floor does not transfer to engine scale);
  ignition/margin/dwell checks run at the declared PROBE_EVERY = 200-step cadence (slow-clock members;
  the per-step members — audit, positivity, non-finite — stay per-step).

  **The class split of record (SOLV-4 0.4.4 / COUP-3 0.4.4):** the auto-ignition term moved from the S6
  explicit tier (retired) to the **cell-local implicit class-R occupant** — a fixed-structure
  backward-Euler node solve per SDC sweep (linear in ρc at frozen τ_ign; τ's weak state dependence
  converged by N_TAU_REFREEZE = 2; realized applied-increment rates ride the trapezoid quadrature and the
  audit's burn_progress row; the node-0 rate is capped at the realizable (cap−ρc)/Δt — the raw ρ/τ at
  extreme stiffness overshoots the composition, a measured hazard). Applied after each sweep's accepted
  composition + SRD, so the next sweep's class-A propagation sees the auto-ignited b (seed→propagate).
  ONE treatment across the whole regime, no stiffness branch. Acceptance: a superheated 15-bar open tube
  at T_u = 1888 K, τ = 1.758e-7 s, **dt/τ = 3.04** — past the explicit positivity bound where S6 halted —
  marches to burned 0.999 with b parking **bit-exactly** at 1−BURN_COMPLETE; a unit test drives the node
  solve at w/τ = 10⁶ (parks exactly) and the mild limit (reduces to the explicit rate within the measured
  ~1% τ-drift of BE semantics).

  **The real quench_box — FLAMEOUT demonstrated (SOLV-4 §6.5 realized):** a lit front in the 0.8 mm-gap
  tube between cold NoSlip+isothermal walls, marched with the full blend↔class-D coupling (flow + gas
  diffusion + class-R in one audited step — the S7 coupling built for the startup). At 0.1 atm / 400 K
  walls the front DIES: R collapses 2.06e-4 → exactly 0 kg/s, burned 0.150 — the FLAMEOUT R-trajectory
  the verdict consumes — while the adiabatic control on the same fixture holds R at 2.20e-4 (alive).
  420 K walls are the coldest the GAS-PHASE model honestly supports: the products surface's own envelope
  floor is ~407 K at the fixture state (declared envelope; grid floor ~332 K, H₂O condensation ~310 K — measured from the artifact); colder walls are
  S15 two-phase territory. Getting here forced the **cold-side partition extension** (SOLV-4 0.4.4):
  wall cooling breaks the adiabatic shared-h flamelet identity, so below the reactant branch's h-floor
  the reactant sub-state pins AT the floor and the products absorb the balance (continuous, mass-
  consistent, self-limiting at the burnt branch's own envelope; trace-weight floor B_PARTITION_MIN = 0.01
  regularizes the 1/b amplification), plus the **below-unburnt-floor face of the non-reactive guard**
  (cells colder than any representable reactant are declared no-burn — reached because b parks just under
  the 1−BURN_COMPLETE fixed point and keeps closure queries live forever, the second knife-edge of the
  S6-close class).

  **The cold-side closure floor (OFFL-3 0.6.3/0.6.4) — the S6 review-wave flag cured as the composite:**
  ignition v0.3.0 (strict extension, verified bit-exact column-by-column): graded cold T_u rows 60/90/120 K
  (real Cantera solves — 120 K carries real S_L 0.5–4.3 m/s across the row; 60/90 K fail-to-zero across
  the live envelope at the O₂ condensation edge, honest contiguous extinction) putting the envelope floor at 75 K, plus a p-axis top-cell midpoint
  node at 4.443 MPa lifting the p-envelope ceiling 4.44 → 6.67 MPa — the S7 stiff mini-sim caught
  mid-transition cells compressing past the old ceiling (a coarse-p-axis inset artifact, not physics; the
  p-ceiling face of the envelope-consistency contract), and an armed **isolated-zero generation refusal**
  (an interior S_L = 0 hole along T_u = a failed solve, refused not recorded). Runtime: the declared
  non-reactive floor below the surface's own envelope floors (the cold analogue of the BURN_COMPLETE
  domain guard; one owner — the bound artifact's envelope attrs). Unburnt v0.3.0 (a declared RE-GRIDDED
  variant per the R2 envelope doctrine, NOT a strict extension — no certificate consumes it): Z narrowed
  to the burnt surface's own band so the rectangular h-floor stops binding ~60 K hot of the design line;
  cold face at a 75 K binding edge = mid-Z validity to ~96 K, covering the 120 K wall-cooled fill states
  the startup march actually holds (45+ K above their real O₂ saturation); sub-floor overhang limited to
  the base axis's own ~2-cell margin (deepest binding-corner node ~40 K, γ ≈ 1.15 — no deeper cold
  prepends into the degraded-polynomial region); the wide-Z 0.2.0 grid stays as the S5b base.

  **SOLV-6 v1 (0.3) + margins wiring:** `structural_margins.rs` (SOLV-6.1–6.5, FS_YIELD 1.1 / FS_ULT 1.4,
  7 closed-form unit tests). The COUP-4 halt inputs are MELT (surface ≥ solidus) and **BURST_MARGIN on the
  primary (pressure-difference) stress state** — the ASME primary/secondary categorization: a regen
  liner's thermal stress is strain-controlled and legitimately exceeds elastic yield locally, so the
  combined-stress elastic margins are reported diagnostics, never halts (an elastic yield-halt would kill
  every real cooled liner). The shell is the DECLARED pressure-carrying member (the RL10's brazed tube:
  r 3.5 mm / wall 0.33 mm, 2R/t = 21.2; loaded by |p_gas − p_coolant| with the declared ~6.9 MPa jacket
  backpressure); 2R/t ≤ 20 refuses at assembly (reported, never smeared); allowables = two-point A-basis
  pair, cold anchor at 77 K with RT strengths (declared conservative flattening — the chilled liner
  interrogates at ~120 K at start).

  **Config faces (FND-4 0.3 / COUP-7 0.4.2):** three new type-keyed mechanisms — `combustion_blend`
  (chem_unburnt + chem_ignition pins + the declared wrinkling, now consuming the **pinned**
  `turbulent-flame-speed` key: Zimont 2000 / Peters 2000 class; the startup config declares 1.0 — the
  laminar tier — after the measured finding that wrinkling > 1 multiplied into the S_L crossover band
  drives the front-carrier S_T toward the sound speed and trips the explicit class-A guard (the S_T-CFL
  hardening rides S8/S16 with the dynamic closure), `spark_igniter` (energy_j = the ONE cited datum,
  exciter-class range gate 1e-5–20 J; placement + §3.2.2 firing window; the deposit ramps over 0.3 of the
  window and integrates to exactly energy_j over the kernel's measured κV volume; at N_θ = 1 the kernel
  is a declared one-cell RING — the declared 10 J discharge heats it past auto-ignition at the measured
  fill densities),
  `structural_margins` (13 params). `valve_cited_timeline_s` records the compressed timeline's citation
  (the ~2 s RL10 sequence) in the manifest. The engine grew the ONE config-selected EOS seam
  (`ChemEos: Table | Blend` — the dispatch lives at the seam once), blend warm-started through the
  pure-limit delegation (which also brought the near-vacuum tangency acceptance to b = 0/1 for free, one
  owner; mid-transition gains its own mass-weighted-bound acceptance).

  **Startup shake-out findings (each a recorded cure, found by the dial-3 probe marches):** (1) the
  Robin-Robin exchange acceptance failed on the near-vacuum cold fill — the wall-adjacent gas cell's
  thermal mass is ~10⁴× smaller than the stations' and the fixed sweeps land at ~3e-5 relative on a ~10 W exchange
  (0.3 mW of staleness): EPS_ROBIN_RESID relaxed 1e-6 → 1e-4 (a halt-gate, not a solution modifier — no
  accepted number moves, certificates byte-identical) + an absolute EPS_ROBIN_Q_FLOOR_W = 1e-2 W
  (COUP-2 0.3). (2) The injector-face fixed point transiently dips below the reactant h-floor when
  choking into vacuum (the first iterate evaluates sound speed at the total state) — cured by the
  declared gas-phase injection state at −2.5e5 J/kg (~215 K premixed; the real cryo-liquid state is
  S15's, recorded limitation). (3) The pre-ignition cold jet over-expands below the reactant model's
  condensation edge at cell backpressures under ~10 mbar — the declared altitude cell is set at the
  10–20 mbar ejector-cell class and the spark fires early in the fill (as a real sequence does; chamber
  crosses the ignition surface's 6.8 kPa floor at a measured ~0.9 ms). (4) Margin allowables interrogated
  at the 120 K chilled liner → the 77 K cold anchor. (5) An un-ignited fill march ends in the CORRECT
  refusal: 14 ms of accumulating unlit propellant compresses off the reactant model — the hazard state a
  real engine cannot sit in either.

  **◆C2 — the certifying startup march (`configs/rl10_startup.toml`) and what it computed.** The
  shake-out was the season's physics in miniature — thirteen attempts, every halt a distinct finding,
  each cured and battery-verified before relaunch: (1–2) NEVER_IGNITED verdicts taught the spark's
  form (the fill chamber is a 200–400 m/s stream; a deposit longer than the ~60 µs kernel residence
  pays the ignition enthalpy once per gas replacement, and firing late into a densified chamber buys
  only ms-class τ_ign parcels that advect out — the cure is the exciter's own short intense burst);
  (3–4) the positivity blowup exposed the class-R quadrature-transient (the invariant-set projection
  fix) and the S_T-CFL limit of the explicit front-carrier in the S_L crossover band (wrinkling > 1
  there drives S_T toward the sound speed — the laminar tier certifies; the turbulent consumer rides
  the S8 guard); (5–7) the blend's projection learned its full validity structure the hard way — the
  two-sided partition (spark superheat past the reactant ceiling), the pulse-cut-at-light discipline
  (the tail of a 50 µs window superheats the blast-rarefied kernel wisp ×100 faster than b escapes to
  pure-burnt), and the establishment grace (a just-lit front at 94% of the ṁ-scaled floor 38 µs after
  window end is establishing, not absent); (8–13) the TORCH pivot and its consequences — a lit kernel
  CANNOT anchor as a flame in the fill stream (blowoff: laminar S_T ~5 m/s vs 170 m/s flow; the prior
  tier's uniform plane inflow has no recirculation), so flame-holding is the ASI torch's job exactly
  as on the real engine (TM-107318: the ASI burns continuously through start into mainstage) — the
  igniter object gained the cited kJ-class torch tier — and the supersonic cold-purge exposed the
  mid-transition projection's bracket (now the products branch's whole h-window, with a fixed
  log-scan-first root isolation and a deterministic bisection backstop in the ONE root finder —
  a strictly-fewer-halts change, previously-converging projections bit-unchanged).

  **The certified march itself:** valve opens into the 10–20 mbar declared altitude cell at t = 0;
  chamber crosses the ignition surface's 6.8 kPa floor at 0.9 ms; the torch (21 kJ declared window,
  ~300 kW post-ramp) lights the kernel at ~3.04 ms; **chamber-scale light-off at ~24 ms** (R jumps
  6.9e-7 → 1.3e-2 kg/s in ~4 ms; p_c through 0.86 MPa, full 16.9 kg/s delivered); the pre-light cold
  bell gas purges supersonically (u_z ~ 1060 m/s — the projection-hardening states); then a stable
  torch-anchored flame (R ≈ 1.6e-2, 10× the FLAMEOUT floor) with p_c PLATEAUED at ~0.90 MPa and
  F ≈ 22 kN — implied c\* ≈ 660 m/s: **the chamber runs mostly cold.** The laminar front cannot
  spread across a face swept at 100+ m/s; the spreading agents — turbulent flame speed, the 216
  distributed coax elements, 3-D recirculation — are precisely the plan's S8/S16 physics.
  **The verdict (55,658 steps, the full 84.06 ms horizon marched, audit clean throughout):**
  `FAILED_TO_REACH: the Stage-1 horizon (8.4065e-2 s) expired without a completed dwell — p_c_pa
  9.085571e5 vs commanded 3.150000e6 (−71.16%); thrust_n 2.206685e4 vs commanded 7.560000e4
  (−70.81%)` — verdict object `DOESN'T WORK (FAILED_TO_REACH; diagnosis physical)`, written to
  `runs/rl10-startup/verdict.txt` with the crash artifact beside it. Final state: R = 1.7e-2 kg/s
  (the torch-held flame never wavered — 10× the FLAMEOUT floor for 60 ms), ṁ delivered 16.97 kg/s,
  no margin/melt/audit trip over the whole march.
  "DOESN'T START" is now a computed outcome with a true mechanism: at the laminar-coarse tier with a
  single torch, the RL10 does not reach its commanded operating point — and the instrument says so,
  names the shortfall, and the diagnosis maps one-to-one onto the physics the plan already schedules.

  **Dial-16 establishment:** the cure is demonstrated in kind at dial 5 (the blend startup march IS
  the establishment path the KNOWN LIMIT awaited); a dial-16 blend startup costs ~33× the dial-5
  march (≈ 10⁵ s laptop) — recorded as the desktop-tier command (same config, dial 16), not run here.

  **Gates: all five green.** fmt/clippy clean; gate 3 workspace tests pass (the battery grew to 9 —
  stiff_auto_ignition, implicit_node_solve, cold_floor, quench_box joined the five — plus
  structural_margins' 7 in-module tests); gate 4 offline suite green (unburnt/ignition suites
  repointed to the v0.3.0 artifacts, regen probes parameterized); **gate 5: all six certificates
  regenerate byte-identical** — the entire S7 surface (blend partition, class-R, the shared root
  finder's bisection backstop) never moved a station number, exactly as designed (the backstop runs
  only where the old code refused). Review wave: two agents (code/physics: NO confirmed correctness
  bugs, three plausible findings all fixed — the diagnosis mislabel, the blend transport-gate
  union, the volume-weighted tangency bound; doc-claims: ~15 number drifts corrected against the
  artifacts, four flags all resolved incl. the 420 K quench-wall correction and the OFFL-3
  rich-Z cold-corner KNOWN LIMIT). ◆C2's own shake-out then added the run.rs establishment grace,
  the two-sided-plus-fallback partition, the products-window bracket + scan-first + bisection
  backstop — each battery-verified before its relaunch.**

- Session 20 (2026-08-25): **plan S8, AZIMUTHAL CAPABILITY — Phase 3 opens.** The flow operator went
  genuinely 3-D: the r = 0 axis at N_θ > 1 and mixed per-brick N_θ marched, gated, and audited.

  **Docs first** (the working rule): FND-2 0.5.2 (§3.4's "conservatively aggregated/subdivided" clause
  made concrete — 2:1 ladder adjacency, THE FINE SIDE OWNS THE FLUX, aggregate applied coarse-side,
  proper nesting ≥ NGHOST, cut geometry stays uniform; §3.2 axis parity pairing recorded as built);
  COUP-2 0.3.1 (interface faces interior to the one ledger; the §3.1.1 S[q] stored term = the GROSS
  Σκv|q| at both endpoints — the finding below); COUP-3 0.4.5 (the class-A wave-speed members stated:
  meridional + per-brick azimuthal + SOLV-4's σ_front; the gas class-D θ-design fixed with the build
  split to S9/S11); SOLV-1 0.4.6; SOLV-4 0.4.6 then 0.4.7 (the review wave); META-3 0.8.5.

  **The axis crossing (FND-2 §3.2 as built):** the innermost ring's cross-axis ghosts gather from the
  θ+π partner with u_r AND u_θ negated (the basis flip); at N_θ = 1 the partner is the cell itself, so
  the certified axisymmetric corner is arithmetically identical (gate 5's byte-identity is the proof).
  Gates: `axis_pulse_3d` — an off-axis 3-D pulse crosses the axis with the armed audit every step,
  reaches the far side (θ ≈ π) through r = 0, preserves θ → −θ mirror symmetry to 1e-10, and conserves
  mass/energy to 1e-12; the exact gate — an AXISYMMETRIC pulse at N_θ = 8 reproduces the N_θ = 1 march
  **bitwise per θ-plane** (the metric factors differ by exact powers of two, so the arithmetic cancels
  exactly); and the review-mandated sign discriminator — a uniform TRANSVERSE Cartesian flow
  (u_r = U·cosθ, u_θ = −U·sinθ) holds steady through the axis (worst u_θ drift ≪ 0.05·U; axisymmetric
  and mirror-symmetric gates are structurally blind to a wrong u_θ flip — this one is not).

  **Mixed-N_θ refluxing (FND-2 §3.4 as built):** the r/z sweeps decompose each pencil into maximal
  uniform-N_θ segments and process each run in TWO PASSES — reconstruct+flux everything (fine ghosts
  prolonged piecewise-constant from the coarse side, coarse ghosts pair-mean-restricted from the fine
  primitives), then replace each coarse side's jump-face entry with the exact aggregate of the fine
  children (plain sum where af carries area; the exact ×½ area ratio in the z-sweep's metric-ratio
  form), then accumulate as the SINGLE well-balanced difference (af[q]−af[q+1])/(κV). **The first cut
  used push-style split accumulation and the gates caught it at step 0:** splitting the difference into
  separate roundings breaks the exact cancellation against the geometric pressure source — the uniform
  state stopped being a bitwise fixed point at the interface. The two-pass form restores it exactly.
  Gates: uniform state bitwise fixed point on a mixed-N_θ axis world; an axisymmetric transient on a
  mixed world reproduces the N_θ = 1 march bitwise (refluxing active every step, BOTH fix-up branches —
  the right-finer and the review-flagged left-finer orientations — covered in r and z); genuinely-3-D
  data across r- and z-interfaces conserves to 1e-12 with the audit armed; 1-vs-4-thread bit identity.
  Typed refusals: ladder membership, 2:1 face adjacency, ≥ NGHOST proper nesting; combustion at mixed
  N_θ → S11; cut geometry/SRD stay N_θ = 1 (FND-3 3-D wave).

  **The audit-tolerance finding (COUP-2 0.3.1):** the first genuinely 3-D mirror-symmetric world halted
  the audit on its own arithmetic — momentum_theta's stored total cancels to ~0 by symmetry, so the old
  |net-total| term in S[q] collapsed the tolerance below the reduction's cancellation rounding
  (measured: delta 4e-22 vs tol 7e-24 on a healthy step), and the θ-sweep's large equal-and-opposite
  increments never enter the ledger's gross throughput (interior faces telescope, unledgered). Cure per
  §3.1.1's own rationale ("scaled by what was actually summed"): S[q]'s stored term is the gross Σκv|q|
  at BOTH step endpoints (`reduce_kappa_volume_weighted_abs`). Halt-gate only; certificates
  byte-identical.

  **θ-CFL + controller in anger:** the Δt rule's azimuthal member runs at each brick's own N_θ; the
  plan-§3 gate holds (a coarse-inner N_θ(r) profile lifts Δt > 1.5× vs uniform-fine on an axis world).
  The symmetry controller ran mid-march: an axisymmetric flow collapses every brick to N_θ^guard = 4
  (ΔKE ledgered), a seeded m = 1 field re-expands to N_θ^max, the Euler workspace re-keying per brick
  (staleness rebuild) with the audit armed throughout. Deferral: a production 2:1-enforcement pass over
  controller decisions (bricks adapt independently; ladders spanning ≥ 3 levels can transiently violate
  adjacency) rides S11 with the profile controller's first engine use.

  **N_θ > 1 combustion (◆C3's prerequisite):** the whole SOLV-4 §3.6 operator re-keyed per θ-plane —
  rate accumulation, the class-R node solve + quadrature buffers (per-brick plane sizes), the
  diagnostics — and ∇·(ρD_c∇c) gained its θ-direction faces (periodic within-brick ring stencil,
  face ρD_c = two-cell mean, exactly zero cost at N_θ = 1). Gate: a spark kernel that is a POINT in θ
  (sector 2 of 8) lights its sector and the front reaches the ADJACENT sectors while the opposite
  sector is still dark (ordering asserted at first adjacent light-off), through the real SDC step with
  class R and the audit armed.

  **The S_T-CFL (SOLV-4 0.4.6/0.4.7) — wrinkling > 1 becomes marchable:** the front-carrier's
  stability is no longer a premise ("S_T ≪ a") but a Δt-rule member: σ_front = 2·D_c·Σ_d 1/Δ_d² +
  C_NAGUMO_SLOPE·(ρ_u/ρ)·K joins the class-A wave-speed reduction wherever the rate law is live (both
  matched coefficients scale with Δ, so σ_front ~ S_T/Δ is CFL-class, not parabolic). The
  scale-separation guard is stated on MODEL-FORM quantities — S_T > (2/3)·c refuses, typed and
  cell-named. That form is itself a review catch: the first cut compared mesh rates (σ_front vs the
  acoustic signal), and the reviewer measured it resolution-dependent — the θ-arc's 1/arc² carrier
  rate outgrows the 1/arc acoustic rate, so a rate-ratio guard tightens linearly in N_θ at the
  innermost rings until it refuses mild flames the certified tier marches (S_T/c thresholds 0.41 →
  0.043 from N_θ 8 → 64). Velocities are resolution-independent; mesh rates belong to the Δt member.
  Gates: dt(wrinkling 8) < dt(1.0) on a mid-flame state; a sonic-class wrinkling refuses typed; a
  mild flame near the axis at N_θ = 16 does NOT refuse.

  **Review wave (three agents) + the fix wave — every CONFIRMED finding fixed in-session:**
  *Axis/reflux math (the plan-§7 risk):* the mesh-rate guard above (CONFIRMED, fixed); left-finer
  fix-up branch untested (fixed: two new gates); u_θ sign undiscriminated (fixed: the transverse-flow
  gate); brick(0)-sampled N_θ guards panic on disconnected mixed worlds (fixed: all-brick scans);
  stable_dt +inf on an empty active set (fixed: typed refusal); conservation bookkeeping at jumps,
  ghost-index walks, θ-ring bitwise telescoping, determinism partitions all verified sound.
  *S7-carry fresh-eyes (the mandated re-exam):* the blend partition's balance↔trace switch was
  DISCONTINUOUS inside the projection's own scan bracket (a pseudo-root site) and the B_PARTITION_MIN
  crossing stepped the mixture volume by (1−b)·δ·∂v_b/∂h — the 1/b amplification cancels the b weight,
  so the step was NOT trace-bounded (the ◆C2 purge cells sat at the measured knife edge b = 0.010);
  CURED by the continuous window form h_b = clamp(balance, products window) (SOLV-4 0.4.7) — the
  trace form and B_PARTITION_MIN retired. The ONE root finder accepted on bracket width alone — a
  pseudo-root could return with a finite unchecked density mismatch; CURED: the volume-weighted
  rule-space acceptance is armed on EVERY accepted root, and the blend's true cold edge is enforced
  there. The class-R base projection was asymmetric (negative quadrature artifacts integrated through
  the BE) and shaved advected ρc ≥ cap excursions silently; CURED: symmetric [0, cap] base projection,
  and advected at/past-cap content takes the zero-source exact trajectory (owned by the loud
  EPS_BURN_BOUND guard). Also corrected: the S7 log's "previously-converging projections bit-unchanged"
  claim was too broad — it holds for the bisection backstop itself, not for the mid-transition bracket
  re-derivation that landed with it (recorded here; stations were never affected — they schedule no
  blend). Recorded deferrals (owner: the S9 review wave): a class-R dt-Richardson order gate for the
  w/τ ~ 1 mid-stiffness regime; a mid-transition root uniqueness/branch-residency guard (the scan
  takes the first crossing); a measured N_TAU_REFREEZE convergence check in that band.
  *Doc-claims:* five mismatches found and fixed (the "conservative restriction" overclaim → primitive
  pair-mean stated; the S[q] both-endpoints wording; the META-3 consumer-provenance sentence — the
  declared 2.0 is the pin's band factor, not the raw Zimont estimate ~6; a stale plan-S8 pointer; the
  refusal-payload naming); everything else verified accurate against the code.

  **SPLIT (improvisation rule, plan §8 v1.8):** the gas-diffusion F_visc θ-extension. The DESIGN is
  fixed (COUP-3 0.4.5: θ-θ implicit cores — μ for u_r/u_z, (4/3)μ·r̄² for the ω angular-momentum form,
  k for T, ρD for C — with the curvature/cross couplings on the same fixed Picard lag as the S3
  cross-stress terms); the BUILD rides S9. Reason: a PARTIAL θ-stress tensor is wrong physics — the
  θ-θ core alone spuriously damps m = 1 translation modes whose curvature-coupling partners are absent
  (u_r = U cosθ, u_θ = −U sinθ has zero true stress; the partial operator sees −μU cosθ/r²) — and the
  full tensor + total-energy θ-work bookkeeping + a θ-MMS order gate is a session-scale block that
  would have crowded the §7-risk axis/reflux work this session exists to get right. The refusal stands
  loud (N_θ > 1 gas diffusion names S9; mixed-N_θ class-D and per-θ wall patches name S11). S9/S10
  must land the uniform θ-tensor before ◆C3 marches 3-D with diffusion.

  **Gates: all five green** on the final reviewed code — fmt/clippy clean; gate 3 grew the Rust
  battery 184 → 209 (the 11 azimuthal gates in `solv1_azimuthal.rs` + the point-spark, S_T-CFL,
  and near-axis-guard gates in `solv4_combustion.rs`, plus review-wave additions); gate 4 = the
  50 offline tests; **gate 5: every certificate regenerates byte-identical** — the axis parity
  gather, the two-pass sweep restructure, the audit-tolerance change, and the blend/class-R review
  cures never moved a certified bit (the N_θ = 1 arithmetic-identity claims, proven).

  **The ◆C2 rerun (wrinkling 2.0) marches as this entry is written** — launched on the final
  reviewed build (`runs/rl10-startup-s8`; the first launch was killed and restarted when the
  review wave's blend/class-R cures landed — a rerun on a stale binary would not be reproducible
  from this commit). Early trajectory vs S7's laminar march: kernel lights on the same ~3 ms
  schedule; p_c at 9.2 ms already 0.50 MPa (the laminar run crossed ~0.5 MPa only after its
  ~24 ms chamber-scale light-off) — the 2× front is spreading measurably faster. The verdict +
  measured delta land in a follow-up commit when the ~3 h march completes.

- Session 21 (2026-08-25): **plan S9, FULL GEOMETRY + the S8 carries.** The FND-3 geometry kernel
  exists (CSG, STL, PLIC, voxelization), the F_visc θ-stress tensor is complete at uniform N_θ
  (the ◆C3 blocker), cut geometry runs at N_θ > 1, and the three S8-recorded carries landed —
  one of them a genuine finding.

  **Docs first** (the working rule): FND-3 0.4 (the as-built block: exact winding evaluation with
  Barnes-Hut a recorded perf deferral; PLIC normal sources per authoring path; the S-Z offset by
  deterministic fixed-count bracketed bisection on the monotone V(d); per-face-class sampling
  measures; the SplitMix64 counter key scheme with the θ-sector deliberately excluded; the S9
  conservative N_θ^geom floor; measured C_jitter + rate recorded); FND-2 0.5.3 (§3.4(iv)
  superseded — cut geometry legal at UNIFORM N_θ, six-aperture per-θ-sector storage, the
  wall-closure θ-limb W_θ = (a_θ+−a_θ−)·A_θ, mixed-N_θ cut → S11); SOLV-4 0.4.8 (the carry
  designs, then the finding disposition); COUP-3 0.4.6 + SOLV-1 0.4.7 (as-built); the FND-4
  schema comment re-pointed (the CSG/STL config GRAMMAR rides the engine 3-D assembly wave,
  S10/S11 — its first consumer).

  **The geometry kernel (`crucible_grid::geom3d` — a subagent build, math-reviewed in-session):**
  the analytic SDF tree (sphere/box/cylinder-z/cone-z/torus exact; the revolved-(s,z)-polygon
  leaf the engine contour corresponds to; min/max booleans with the §3.1 sign-exact/magnitude-
  bound caveat; every node 1-Lipschitz ⇒ the |field| > R_circum pure-cell test is sound); STL
  import (binary + ASCII, the 84+50n size identity disambiguates "solid"-headed binaries,
  degenerate triangles counted never silently dropped) with the EXACT van Oosterom–Strackee
  winding sum — a tetrahedron with a deleted face still classifies where ray-parity provably
  flips (gated); the jittered-stratified voxelizer in the cylindrical measure (uniform in
  (r², θ, z); apertures in each face's own measure — r: (θ,z), z: (r²,θ), θ: (r,z) planar);
  PLIC via the Scardovelli–Zaleski corner-sum V(d) (verified against hand values) with
  area = dV/dd exact and the volume match by 80 fixed bisections. **Measured:** C_jitter 0.907
  → declared 1.4 (×1.5 margin, recorded in FND-3); sampling-rate exponent −0.723 vs the derived
  −2/3; sphere PLIC area sum 2.77% (first-order-in-curvature class, as declared). **Two build
  findings, doc-amended:** a watertight mesh's winding number is piecewise CONSTANT, so its
  gradient cannot supply the PLIC normal — mesh cells take the area-weighted outward facet
  normal (CSG keeps the SDF gradient); and the jitter key excludes the θ-sector index so an
  axisymmetric solid voxelizes bit-identically per sector (up to the classifier's rotational
  round-off — the honest form of the §6.6 zero-variance gate) with canonical face arrays giving
  shared-face bitwise coherence by construction. Ten gates green (convergence rate, exact ring
  forms, imperfect-STL robustness, CSG↔STL cross-check, planar-cut exactness, sphere area,
  determinism/seed-liveness, θ-congruence + floors, STL round-trip).

  **The F_visc θ-extension (the S8 split repaid — COUP-3 0.4.5's design, built):** every solved
  component gains its θ-θ implicit core on the periodic within-brick ring stencil (μ for u_r/u_z,
  (4/3)μ·r̄² for ω — both ring cells share r̄ exactly, so the coefficient is symmetric and the CG
  stays SPD; k for T; ρD for C); the curvature/cross couplings ride the same fixed Picard lag as
  the S3 cross terms (the τ_rθ θ-limb (1/r)∂u_r/∂θ at meridional faces + r∂ω/∂r at θ-faces, the
  τ_θz mirror, and the ∂ω/∂θ dilatation limb extending e_θθ in τ_rr, τ_zz, AND the −τ_θθ/r
  volume source — the m = 1 curvature partners); θ work/conduction/species fluxes enter the same
  total-energy bookkeeping (dissipation still emerges from the KE ledger). Buffers went
  per-brick-θ-plane through sdc.rs (the ensure_rb S8 pattern); the CG dot's partials NEST per
  (brick, θ-plane) into the fixed tree — at N_θ = 1 exactly the prior per-brick list, and on an
  axisymmetric N_θ = 2^k world the equal per-plane partials combine EXACTLY (pairwise doubling;
  the power-of-two θ-ladder is what makes this work), so the CG's accept/iterate decisions are
  bit-identical per plane to the N_θ = 1 march. **Gates:** the m = 2 θ-MMS with the full
  continuous θ-stress residual (an independent formulation, every ∂θ limb) — fine-pair L1
  orders 2.03–2.61 on all seven components; **the trap gate**: uniform transverse flow
  (u_r = U cosθ, u_θ = −U sinθ — zero true stress, the state a partial tensor damps at μU/r²)
  measures interior residual 0.9% of that scale at N_θ = 8 and 0.09% at N_θ = 32 (Δθ²
  convergence; species stays an exact zero); an axisymmetric no-swirl coupled march at N_θ = 8
  reproduces N_θ = 1 BITWISE per θ-plane through the full flow+diffusion step; mixed-N_θ and
  cut-θ class-D refuse typed (→ S11). N_θ = 1 arithmetic identity held everywhere (the whole
  existing battery + gate 5).

  **The S8-recorded carries (a subagent build) — landed, one FINDING:**
  - *Root uniqueness (SOLV-4 0.4.8):* the blend's mid-transition scan now completes its full
    fixed N_P_SCAN sweep and COUNTS crossings — more than one admissible root is a typed refusal
    (branch residency undecidable), the TableEos warm-path straddle analogue the mid-b path
    lacked. Single-root brackets are bit-identical to the old first-hit path (verified in-suite);
    a genuine two-root fixture is not constructible on the production surfaces (monotone v(p)),
    so the counting logic is pinned on synthetic closures.
  - *Class-R dt-Richardson (the mid-stiffness blind spot):* measured temporal order **1.80** at
    dt/τ ≈ 0.4–3 on a marched superheated tube (2.09/1.84 on horizon/step variants; the parked
    3.2τ variant collapses order exactly as the doc predicts and is guarded out). The composed
    step holds its order where the class-R truncation dominates.
  - *N_TAU_REFREEZE (the FINDING):* the 0.4.4 sizing premise — "τ depends on the unknown only
    weakly (Δp/p per node is CFL-bounded)" — is FALSE in the reaction-driven-compression band:
    within ONE node solve the burn's constant-volume compression heating drives T_u 1057 → 1512 K
    and τ down ×1/93 at the fixture state; the 2-refreeze solve lags the converged frozen-τ fixed
    point by rel Δ 3.0e-1 at w/τ = 0.3 (2.2e-2 at 1, 2.1e-3 at 3; the reference converges to
    machine precision by pass ~13, so the witness is well-posed). **Disposition (this session,
    recorded in SOLV-4 0.4.8 + the constant's doc comment): the contract is restated, the
    constant is NOT resized** — the fixed count is structure (COUP-3 §3.7), the per-node lag is a
    temporal-truncation-class term the SDC sweeps' own re-evaluations absorb (which is exactly
    what the order gate measures through the same band), the same split as the truncated
    gas-diffusion Picard; the witness ships ARMED as a measured-envelope regression pin
    ({4e-1, 3e-2, 3e-3}) so lag GROWTH fails loudly. The two results are coherent: the isolated
    sizing claim was false while the composed step was always fine — the doc now says what is
    actually true.

  **3-D apertures incl. θ-faces (a subagent build, diff-reviewed in-session; the FND-2 0.5.3
  wave):** `build_with_geometry`'s N_θ = 1 restriction is retired for genuinely 3-D worlds —
  `BrickGeom` went θ-plane-major with the FULL six-face aperture set (FaceDir gained the θ pair;
  every pre-S9 `[_;4]` consumer verified 4-long), `build_with_geometry_theta` validates per
  sector (bitwise shared-face coherence incl. the θ-pair rule aperture[θ+][j] ≡ aperture[θ−][j+1];
  the covered-face rule; the S9 SCOPE RULE — a θ-sector fully covered inside a gas ring is a
  typed refusal, per-sector activity masks ride S10/S11; thin walls legal), the wall closure
  gained its θ-limb W_θ = (a_θ+−a_θ−)·A_θ (per-sector `wall_closure_cell`), the sweeps
  aperture-weight θ-faces with ONE canonical side per face (face fi reads θ+ of ring cell
  (fi−1) mod n — both accumulation directions identical bits, ring telescoping exact), SRD went
  per-sector with θ-neighbors as flow-connected candidates (tie order r−, r+, θ−, θ+, z−, z+; at
  N_θ = 1 the θ-candidates are structurally absent), and cut-geometry bricks pin
  n_theta_geom_floor at the built N_θ (coarsen/refine/assert refuse on geometry bricks — adaptive
  N_θ on cut worlds rides S11). The certified N_θ = 1 path held bit-identity ON THE FIRST RUN of
  the station batteries; both new arms are mutation-proven load-bearing (neutering the θ-limb or
  forcing θ-apertures to 1.0 fails the well-balance gate). Extra refusals beyond the spec, all
  typed: combustion + cut θ > 1 (its D_c stencil reads the (r,z) geometry view → S11); the
  (r,z)-keyed SRD debit view refuses θ-recruited members (per-θ wall patches → S11). Gates (7):
  the revolved cut world at N_θ = 8 ≡ N_θ = 1 **bitwise per plane** through both builders (the
  gate-5-in-miniature); well-balance on a θ-varying world (ρ/ρE/ρC/ρb bitwise, momenta ≤ 1e-14 —
  the separately-rounded (A·ap)·p associations, the S12 round-off class, documented); audited
  shock-transient conservation; per-sector SRD (a sliver in ONE sector merges through its θ-faces
  while the other sectors stay bitwise untouched); the plane contract; the analytic-supplier toy
  chamber; and **`stl_toy_chamber` END TO END on the real sampled path** — a CSG chamber (big
  cylinder − revolved cavity + a shallow axial rib protruding into the gas over one θ-side)
  voxelized by geom3d, ingested by `build_with_geometry_theta` with its per-sector validation
  accepting the sampled output AS-IS (canonical face arrays ⇒ bitwise coherence by construction),
  θ-varying gas rings confirmed, 55 audited cold-flow steps, mass/energy drift < 1e-12.
  **Recorded seam wart (review-wave flag):** the voxelizer emits FND-2 §3.3(1) face order
  {r−,r+,θ−,θ+,z−,z+}; the grid's `CellGeomTheta` uses FaceDir-index order {r−,r+,z−,z+,θ−,θ+} —
  mapped explicitly at the one seam where they meet.

  **Gates: all five green** on the final merged tree — fmt/clippy clean; gate 3 grew the Rust
  battery 209 → 243 (the 10-gate geom3d battery, the 7-gate θ-geometry battery incl. both toy
  chambers, the θ-MMS + transverse-flow + bitwise-per-plane gas gates, the three carry gates, the
  blend uniqueness units, + geometry-seam additions); gate 4 = the 50 offline tests, untouched;
  **gate 5: every certificate regenerates byte-identical** — the θ-plane-major geometry storage,
  the six-aperture restructure, the θ-stress tensor, the per-plane CG partials, and the blend
  scan-completion never moved a certified bit (the N_θ = 1 arithmetic-identity claims, proven
  again). Certified numbers stay S2/S3-spine of record.

  **Review wave (two agents: voxelization/PLIC math + θ-stress derivation; doc-claims) — no
  CONFIRMED correctness defect.** The math reviewer re-derived the continuous cylindrical stress
  tensor and its divergence independently, replicated the code's exact discretizations in Python
  (the m = 1 residual decays ~×4/×8 per N_θ doubling with no O(1) term — the tensor inventory is
  complete), Monte-Carlo-verified the S-Z corner-sum V(d) + the area = dV/dd identity + the
  signed-normal mirroring algebra, verified the van Oosterom–Strackee winding form numerically,
  and confirmed every sampling measure, the SPD/Jacobian-identity contracts, the ring-telescoping
  with apertures ≠ 1, the wall-closure θ-limb fixed-point sign, and the per-(brick,plane)
  reduction-scaling argument. Four low/nano findings, all cured as doc/comment precision: the
  mesh pure-cell hint's threshold-stability claim scoped to WATERTIGHT meshes (an imperfect
  mesh's w is merely harmonic away from triangles and can cross 0.5 in a triangle-free cell near
  a hole — recorded limit; a watertightness gate rides the production-STL wave); the θ-congruence
  "exact, not statistical" claim softened to the honest ulp-bounded form (kernel + FND-3); the
  z-face torque-arm "exactly" comment corrected to the consistent second-order lumped form (the
  exact moment is (r̄²+Δr²/12)ΔrΔθ; the r̄² form matches the ω inertia and telescopes exactly);
  a ±0.0 corner in the bit-identity argument noted unreachable. The doc-claims reviewer verified
  every S9 doc addition against the code and found six mismatches, all fixed: the battery count
  (209 → 243, +34 — the entry above now says it right; CLAUDE.md's stale 184 fixed with it);
  COUP-3's "six components" → seven; the SOLV-4 change-log "zero new evaluations" overclaim →
  the bounded-cost truth; the FND-3 0.4 change-log row moved to its oldest-first position;
  FND-2 0.5.3's face-order phrasing disambiguated (grid FaceDir order vs the voxelizer's §3.3(1)
  emission order) and its bit-identity wording tightened; the geometry floor stated in its hard
  shipped form (EVERY geometry brick pins at built N_θ, stricter than the kernel's θ-varying
  minimum). Everything else verified accurate against the code, including all measured numbers.

  **Deferrals recorded (owners named):** the CSG/STL CONFIG GRAMMAR + engine 3-D assembly (the
  kernel's first config consumer) → S10/S11 with refinement tiles + ◆C3 (FND-4 schema comment
  re-pointed); per-sector activity masks (full sector coverage inside gas rings) → S10/S11;
  mixed-N_θ cut worlds, mixed-N_θ class-D, cut-θ class-D, per-θ wall patches, adaptive N_θ on
  cut worlds → S11 (all typed refusals); Barnes-Hut winding evaluation → the production-STL
  wave; the FND-2 §3.3(7) sharp-interface grid fields stay dormant (the voxelizer now emits
  their initial condition — PLIC planes live in the VoxelWorld product); the voxelizer↔grid
  face-order unification → the S9 review wave's cleanup list. **The ◆C2-rerun verdict
  (wrinkling 2.0, S8's follow-up) has NOT landed as of this entry** — the parked S8 session owns
  that follow-up commit; this entry deliberately does not restate its early trajectory.
  **◆C2-rerun appendix (the follow-up commit; the march outlived the session commit):** the full
  84.06 ms horizon marched clean — 90,656 steps (vs S7's 55,658: the σ_front Δt member + the hotter
  held flame buy ~1.6× more steps for the same horizon), audit green throughout, verdict
  `DOESN'T WORK (FAILED_TO_REACH; diagnosis physical)` at `runs/rl10-startup-s8/verdict.txt`:
  p_c 0.9218 MPa vs 3.15 commanded (**−70.74%**, vs S7's −71.16%); F 22.06 kN vs 75.6 (**−70.82%**,
  vs −70.81%). THE MEASURED DELTA of the first turbulent consumer: a 2.0× front holds **2.2× the
  flame content** (R steady at 3.8e-2 vs the laminar 1.7e-2, both ~an order above the FLAMEOUT
  floor) and deepens the torch-anchored burning zone — p_c +1.5% relative — but the settled point
  is otherwise THE SAME plateau: light-off timing is schedule-set (~25 vs ~24 ms), and a ~10 m/s
  wrinkled front against a 100+ m/s swept face still cannot spread laterally, so the chamber still
  runs mostly cold. Reported as measured, not tuned (the plan's own instruction): the mechanism
  named by the verdict is unchanged, and the spreading agents remain exactly the scheduled physics —
  S11's 3-D recirculation and S16's distributed elements + dynamic S_T. The laminar-vs-wrinkled
  pair now brackets the tier honestly: front speed alone, even doubled, is not the RL10's missing
  starter.

- Session 22 (2026-08-26): **plan S10, REFINEMENT + THE AMR GATE (+ two S9 carries).** The headline is
  a **measurement, not an assumption**: the AMR go/no-go (ruling #7 / plan §5 / §7's "dynamic AMR
  complexity spiral") is **NO-GO on dynamic front-tracking**, measured. The conservative static (r,z)
  level-interface primitive is built and gated (the meridional sibling of the S8 ring reflux), two S9
  carries are repaid, and the heavy production integration is recorded to S11 — **justified by the
  measurement**, not by budget.

  **Docs first** (the working rule): FND-2 0.5.4 (new §3.6.1 — static (r,z) refinement level interfaces:
  the fine-owns-flux area-weighted aggregate as the meridional sibling of §3.4's ring rule, 2:1 (r,z)
  balance + proper nesting, tiles = value representation; the annular-metric single-difference
  well-balance; the AMR-gate MEASURED NO-GO recorded); COUP-2 0.3.2 (§3.1 (r,z) level-interface faces are
  interior to the one flux ledger — the S8 mixed-N_θ clause's sibling); FND-3 0.5 (§3.4 the
  coarsest-reproducing-N_θ floor built; the face-order seam wart retired — FaceDir is the single owner).

  **THE AMR GATE — a measured go/no-go (`solv10_amr_front_study`).** The closure-set pushed front
  (SOLV-4 §3.6) travels at S_T with a width of a fixed Θ cells at every resolution — coarse ⇒ physically
  WIDE (smeared), fine ⇒ SHARP — so refinement sharpens *where*, not *when*. The startup mission's COUP-4
  verdict is **timeline-driven** (when the chamber fills / pressure rises / the flame reaches a plane), so
  the operative question is whether sharpening the front moves the timeline. The study marched the same
  flame tube at coarse h = 2e-4 m and fine h/2, tracking the front position (b = ½ crossing) at three
  instants. **Measured:** the front timeline agreed to **< 0.5 coarse cells**, a bounded sub-cell
  registration offset that *shrinks* over the march (0.50 → 0.25 cells) — the closure-set front carries
  the same speed on both grids (the separate `flame_1d` gate pins the consumption speed grid-independent
  to ~2%), while refinement bought only a **1.93× sharper** front (thickness ratio, ≈ the
  fixed-Θ-cells prediction) at **4.0× cost per level** (2× cells × 2× steps from the halved CFL Δt).
  **VERDICT — NO-GO:** dynamic front-tracking cannot move the timeline-driven verdict; the uniform-coarse
  vs uniform-fine pair *brackets* static front refinement from above and bounds dynamic front-tracking
  (which can buy at most what uniform-fine buys, paying fine cost only in the band), so both are measured
  not-worth-it for the front. The plan's declared fallback — **static refinement + closure-set speed =
  blurry front, correct timeline (§7)** — stands *measured*. Static refinement's purpose is genuine
  geometry/wall/throat gradient resolution, never front-chasing. **Ben-visible call:** dynamic
  front-tracking (an S10b/S10c that would need a frozen-topology amendment, FND-2 §3.2/§0.5.1) is **NOT
  inserted** — the measurement says the sharpness it buys is invisible to the verdict.

  **The conservative (r,z) level-interface primitive (`solv10_rz_level_interface`, the meridional sibling
  of the S8 θ reflux).** Ruling #7 wants static declared refinement zones first, and the conservative
  level interface is the reusable machinery they ride. Built + gated on the **real HLLC flux + real
  cylindrical `Grid` metric** at N_θ = 1, class-A, uncut: (a) **fine-owns-flux telescopes BITWISE** — the
  conserved flux leaving a coarse cell equals the sum entering its two fine children, to the bit; (b) the
  naive **coarse-owns-flux scheme leaks** (measured 0.8% relative — mutation-proof that fine-owns-flux is
  load-bearing); (c) the **new content vs the θ sibling is the metric** — a z-interface's fine children
  carry *unequal ANNULAR* z-face areas ½(r²_{k+1}−r²_k)·Δθ, so the well-balanced uniform fixed point turns
  on the **single-difference form** (the coarse interface area is the children-sum, not an independently
  metricked ½(r²_hi−r²_lo)·Δθ — the S8 trap restated for the cylindrical metric); the annular
  reconstruction residual is the S12 round-off class (≤ 1e-14, measured 0 on the gate fixture), and the
  **r-interface** (equal-area children) holds the fixed point BITWISE. **Scope, recorded to S11 (the
  improvisation-split pattern, exactly as S8 split the F_visc θ-tensor):** threading this primitive through
  the parallel-pencil `sweep_r`/`sweep_z` + the brick-arena refinement topology is a **cross-pencil flux
  register** (Berger–Colella — genuinely more than the within-pencil ring reflux), and it lands with its
  consumer (a refined 3-D RL10). This is not a budget deferral: the measured NO-GO means front-refinement
  is not the RUN's need, so the primitive is built and proven while the heavy integration waits for a real
  consumer. Mixed level × N_θ, level × cut geometry, level × class-D refuse there, typed.

  **Two S9 carries repaid (both safe — no certified path touched).**
  - *The coarsest-reproducing-N_θ geometry floor (FND-3 §3.4, `voxel.rs`):* the S9 binary form pinned
    `N_θ^geom = N_θ^max` on ANY θ-variation; the built S10 form returns the **coarsest power-of-two ladder
    rung whose θ-coarsening reproduces the cell's κ and all six apertures within ε_α** (per-group mean
    deviation ≤ ε_α) — 1 for θ-uniform, N_θ^max for a sharply localized feature, an intermediate rung for
    gentle variation. Unit-gated at every rung (`coarsest_reproducing_finds_the_ladder_rung`: uniform→1,
    half/half→2, quadrants→4, single sector→8, sub-ε ripple→1); the existing off-axis-sphere floor test
    updated to assert ladder-validity + the strongly-cut row still needs θ (it floors at 8 for that sharp
    fixture — the form is a refinement, not a weakening). `theta_geom_floor` has no production consumer, so
    certificates are untouched.
  - *The voxelizer↔grid face-order unification (the recorded S9 seam wart):* the voxelizer now emits
    `CellCut.aperture` in the grid's **`FaceDir::index` order** {r−,r+,z−,z+,θ−,θ+} (was §3.3(1)'s
    {r,θ,z}), so FaceDir is the single face-order owner and the `CellCut → CellGeomTheta` ingest is a plain
    identity copy — no permutation to get wrong. The two order-dependent consumers updated (the fnd3_geom3d
    slot test's θ↔z swap; the `stl_toy_chamber` bridge → identity); the toy chamber still marches 55 audited
    steps with drift < 1e-12.

  **Review wave (two agents: level-interface conservation math + AMR-study methodology; doc-claims).**
  Both confirmed the load-bearing math sound — the level-interface *design* telescopes by construction, the
  naive-coarse-owns-flux mutation genuinely leaks (0.8%, load-bearing), the coarsest-reproducing floor is a
  genuine safe floor (a level that truly needs 8 cannot falsely pass at 2), the face-order reorder is
  consistent across every consumer, and the AMR bracketing argument (uniform-fine upper-bounds any
  front-only refinement) is valid. Every CONFIRMED finding fixed in-session: (1) **two vacuous `x == x`
  assertions** in the level-interface gate whose comments overclaimed (a `af_coarse == af0+af1` "telescoping"
  check and a `(a0+a1)·f == (a0+a1)·f` "fixed point" check) — rewritten to the honest discriminator
  (conservation-by-aggregate is stated as the by-construction theorem; the *load-bearing* assertion is that
  the naive scheme leaks under the SAME measure while the aggregate conserves; the fixed-point test now
  measures the independent-area residual the single-difference form eliminates). (2) **The "annular is
  ulp-not-bitwise" narrative was asserted but never observed** — the residual is bitwise 0 on the tested 2:1
  off-axis interface; corrected everywhere to the honest form: **the r²-band differences are Sterbenz-exact
  on a 2:1 off-axis interface, so the annular reconstruction is BITWISE (measured 0), and ≤ 1e-14 is the
  declared *safety bound* for the general/near-axis case** (test header + §3.6.1 + COUP-2). (3) The AMR
  "fronts travel at the same speed to ~2%" overclaim (the sampled tracks imply ~7% lab-frame window speed;
  the ~2% is the separate `flame_1d` consumption-speed gate) — restated as the measured **shrinking sub-cell
  offset (0.50 → 0.25 cells)**; the "1.96×" sharpness corrected to the measured **1.93×**. (4) Two **stale
  face-order references** surviving the seam-wart retirement (a `stl_toy_chamber` comment and FND-2 §3.4(iv),
  which had begun to contradict the sibling FND-3 0.5) — both updated to the FaceDir-single-owner truth.

  **Recorded to S11 (owners named):** the CSG/STL **config grammar + engine 3-D assembly** (the geom3d
  kernel's first production consumer — the schema slot exists, FND-4 §6-4; rides ◆C3's refined 3-D RL10);
  the **cross-pencil (r,z) flux-register production integration** + dynamic-AMR-if-ever (needs a
  frozen-topology amendment first, and the measurement says it isn't needed for the front); per-sector
  activity masks; mixed level × N_θ / cut / class-D (all typed refusals when they arrive).

- Session 23 (2026-08-26): **plan S11, ◆C3 — THE FIRST GENUINELY 3-D ENGINE RUN + review wave A.** The
  azimuthal machinery built and gated in the SOLVER at S8/S9 now reaches the ENGINE: a coarse full-3-D
  RL10 lights from a **true point spark**, shows **asymmetric light-off**, and marches spark→settle to a
  typed COUP-4 verdict. ◆C3 proves the 3-D START MACHINERY, not a working engine — and the honest
  measured outcome is **DOESN'T WORK**, exactly as the 2-D ◆C2 (S7) and the plan predicted.

  **The engine N_θ > 1 assembly (`crates/engine/src/assembly.rs`).** The `!axisymmetric || n_theta_max
  != 1` guard (the S9-map's ~line-202 refusal) is lifted. `n_theta_max` from `[geometry]` now threads
  into `GridSpec`; the N_θ = 1 path keeps `build_with_geometry` **bitwise unchanged** (the certified
  stations + ◆C2), and N_θ > 1 builds via `build_with_geometry_theta` on the revolved contour (a
  `_j`-ignoring closure — every sector the identical analytic clip, θ-face apertures = κ, the S9
  `stl_toy_chamber` convention). **Scope, refused loudly at assembly:** the cooled-wall Robin exchange
  (per-θ wall patches), coupled flow+solid conduction, and F_visc on cut θ-faces are the typed S11/S9
  refusals still standing IN THE SOLVER (`sdc::build_wall_patches` refuses N_θ > 1; the `Sdc::step`
  guards refuse coupled flow+conduction and gas-diffusion-on-cut-θ) — so a cooled or viscous 3-D RL10
  lands with that wall-patch/conduction wave. ◆C3 therefore runs the **ADIABATIC flow+combustion start
  core** (`liner_thickness_m = 0`, slip-wall nozzle); the two new assembly refusals name the deferral
  in kind (gated by `s11_azimuthal_3d`).

  **`run.rs` made θ-aware — N_θ = 1 bit-identical (certificates byte-identical, gate 5).** The
  field-weighted plane integrals (`plane_mdot`/`plane_thrust`/`injector_end_stagnation_p`) now sum over
  θ sectors, each carrying its per-sector annular area `face_area_z(i_r, n_θ)·aperture`; at N_θ = 1 the
  single sector's area is the full ring and `cell_value_theta(…,0,…)` = `cell_value`, so the arithmetic
  is bitwise the pre-S11 form (the geometry `aperture` reads the θ-plane-0 view — exact for revolved,
  latent for a θ-varying wall). `max_flame_cell` scans every sector; `fields_csv`/`crash_fields_csv`
  emit one row per (r, θ, z) with a θ column + swirl velocity for N_θ > 1 (the certified N_θ = 1 columns
  unchanged). The quiescent fill + audit reductions were already θ-complete (`fill_field`,
  `reduce_kappa_volume_weighted`).

  **The point-in-θ spark (`IgniterSpec.theta_rad`, default 0 so the certified configs are untouched).**
  The igniter deposit is confined to the SINGLE θ-sector nearest `theta_rad` — a true point, not the
  N_θ = 1 ring: the `IgniterKernel.theta_gate` computes the target sector, the volume-normalization
  loop selects exactly that sector's cells (per-sector `cell_volume(i_r, n_θ)`), and `source_fn` maps
  the continuous θ back to the sector index (`(θ/Δθ).floor()`) — the same cell set, so the deposited
  energy still integrates to exactly `energy_j`. At N_θ = 1 `theta_gate = None` ⇒ the ring, bitwise the
  S7 form.

  **The combustion cut-θ D_c stencil built (the S9 "combustion on cut θ > 1 worlds" refusal RETIRED for
  revolved walls — landed with its gate, exactly the brief's sanction).** The exploratory ◆C3 march
  surfaced the last blocker: on the CUT RL10 contour the combustion front-diffusion refused at N_θ > 1
  (`euler/mod.rs` `validate`), because its D_c stencil read the θ-uniform (r,z) geometry view. Fixed
  per-sector (`combustion.rs`): `kv` uses `kappa_cell(j, local)`, the meridional faces use
  `aperture_at(i_r, j, i_z, dir)`, and the θ-direction faces now **weight by the θ-face aperture**
  `aperture_at(…, ThetaMinus/ThetaPlus)` — a revolved cut cell's constant-θ plane is partly blocked, and
  full dr·dz would over-diffuse the front across a wall-clipped sector; the shared θ-face aperture is
  bit-coherent from both sides (build validation) so the flux still telescopes. **Every change is
  bit-identical at N_θ = 1 and on box worlds** (`kappa_cell(0,·) == kappa_rz`, `aperture_at(·,0,·,·) ==
  aperture`, θ-block skipped at N_θ = 1, ×1.0 exact on box worlds) — the whole N_θ = 1 combustion
  battery (`flame_1d`, `spark_box`, …) stays green. The refusal is **narrowed** to genuinely θ-VARYING
  cut geometry via a build-time `Grid::geometry_is_theta_uniform()` flag (computed once in
  `build_with_geometry_theta` by bit-comparing every sector's κ + 6 apertures to sector 0) — a CSG/STL
  θ-varying wall still refuses, typed, until that wave. `reacting_measure`/`consumption_rate` aligned to
  per-sector κ too (review-wave consistency; bit-identical on θ-uniform).

  **Gates (all new, all green).** `solv4_combustion::s11_revolved_cut_combustion_n_theta_8_is_theta_
  symmetric_and_reduces`: an axisymmetric burn on a revolved CUT world at N_θ = 8 stays **exactly
  θ-symmetric** (every sector bit-identical to sector 0 through the SDC step — the per-sector D_c
  stencil + zero θ-flux introduce NO spurious azimuthal asymmetry) and reduces to the certified N_θ = 1
  cut march (tight tolerance early; finding: unlike the flow sweep, the per-θ face-area accumulation is
  NOT power-of-two bit-reducible, so the sub-ULP round-off amplifies chaotically — different N_θ
  discretizations need only converge, not bit-match); `s11_theta_varying_cut_combustion_still_refuses`
  (a θ-varying κ world refuses, typed). `engine/s11_azimuthal_3d`: a coarse 3-D RL10 **assembles,
  marches audited steps, audit closes**, reads out θ-summed (the fields CSV carries the θ column); the
  cooled-wall and gas_diffusion 3-D configs **refuse at assembly**.

  **Uniform N_θ is structural, not a shortcut.** `build_with_geometry_theta` pins every geometry-bearing
  brick's `n_theta_geom_floor = n_theta_max`, and the controller refuses to regrid a geometry-bearing
  brick — so a revolved-contour world is **uniform N_θ by construction** (the θ-CFL near the axis sets
  the step count; the N_θ(r) controller that would coarsen near-axis rings is a later-wave item, moot on
  the contour's geometry floor). No controller wiring was needed or possible.

  **◆C3 — the run (verdict AS MEASURED).** `configs/rl10_startup_3d.toml` = the desktop-confirm spec
  (N_θ = 8, dial 5, grid 40×119, 84 ms horizon; assembles + marches — Ben's to run to completion). The
  **laptop MINI** (N_θ = 8, dial 3, grid 24×71, commanded p_c 3.15 MPa / F 75.6 kN) marched end to end:
  spark fires at 2 ms into the target sector → an asymmetric reacting kernel grows (R climbs 0 →
  ~1.2e-7 kg/s) but **cannot anchor** (three orders below the ṁ-scaled ignition floor — the laminar
  front can't spread across the swept face, exactly the ◆C2/S16 physics) → the Stage-1 horizon expires
  and the march ends in the typed verdict **DOESN'T WORK (FAILED_TO_REACH; physical): p_c 0.246 MPa
  (−92.2%), thrust 3.93 kN (−94.8%) of commanded** (2610 steps, audit green, verdict + crash artifacts
  written). "DOESN'T FULLY LIGHT" remains the honest expected outcome until the spreading physics
  (distributed injection S16, resolved 3-D recirculation) exists — ◆C3 proves the 3-D START MACHINERY
  (point spark, asymmetric light-off, θ-summed readout, typed verdict), not a working engine.

  **Review wave A (multi-agent review of phases 1–3 + the S11 diff; 4 finder agents by subsystem).**
  Verdict: **no confirmed correctness bug** across the S2 integrator, S3 F_visc, S4 transport, S5 cold
  branch, S6/S7 ignition, S8 azimuthal, S9 geometry, S10 refinement, or the S11 3-D engine path. The
  agents re-derived the load-bearing invariants by hand: the SDC audit genuinely closes and catches a
  violation; F_visc's metric/CG-symmetry/wall-suppression are correct; heat release is EOS-implicit and
  conserved by construction; the Nagumo front speed is genuinely closure-set; the winding number,
  cylindrical voxel measures, PLIC bisection, and level-interface telescoping are sound; and the S11
  N_θ = 1 bit-identity + D_c telescoping + point-θ mapping + θ-uniformity detection all hold. **One
  confirmed finding (LOW), fixed in-session:** the tabulated transport occupant did not enforce the
  Prandtl rail its own doc claimed the constant occupant carries (`[0.05, 5]`, `TRANSPORT_CONSTANT_
  MANIFEST`) — a metal-like-`k` node would pass the column rails yet deliver a sub-0.05 Pr; now railed
  (`PR_RANGE`, defense-in-depth — shipped Pr is 0.2–1.5, no data rejected, S4 spine march unchanged).
  **Recorded to their owning waves (latent, not bugs):** the θ-summed readouts + `reacting_measure` read
  the θ-plane-0 aperture/κ view (exact on revolved, needs per-sector on a CSG/STL θ-varying wall); the
  igniter's sector index derives from `n_theta_max` (agrees with the brick N_θ only via the geometry
  floor); the blend `partition_h` clamp trades the enthalpy identity for the density identity at the
  cold/hot edges (a bounded fidelity edge — worth a mid-`b` deep-quench unit test); the ignition↔unburnt
  envelope-consistency contract is convention-enforced (an offline cross-check would harden it, like the
  transport↔EOS refusal that caught the S6 inversion).

  **Docs:** SESSION_LOG (this entry); CLAUDE.md State; PLAN §8 v1.11. No Reviewed Layer-2 doc contract
  changed (the combustion cut-θ extension realizes COUP-3 0.4.6 / SOLV-4 §3.6 at uniform N_θ, already
  the design of record; the engine 3-D assembly is FND-4/engine-header scope). Certificates
  byte-identical (gate 5). **NEXT = plan S12 (GPU spike).**

  **OPEN DECISION FOR BEN (carried from S10 close, still unanswered): static (r,z) refinement wiring —
  Option A (wire the S10 flux-register primitive into the engine now, for wall/throat/injector gradient
  resolution) vs Option B (defer until the S19 RUN needs the cell budget — my recommendation, justified
  by the measured AMR NO-GO). NOT built this session; the primitive stays proven-but-unwired. New angle
  from ◆C3: refining the injector/shear-layer region may help the RESOLVED turbulent mixing that spreads
  the flame — but only after distributed injection (S16), so even under A it is a measure-first item.**

- Session 24 (2026-08-27): **plan S12 — THE GPU SPIKE (measured on real hardware) + the S11 doctrine
  amendment + remote-GPU bring-up.** Three things landed.
  **(1) Remote GPU box brought up** (the new logistics). Ben had nothing set up; from a standing start we
  established: Tailscale on both machines (PC = `backhouse`, RTX **4070 Ti SUPER** — NOT the assumed 4080;
  16 GB, 256-bit, ~672 GB/s, driver 591.86); Windows OpenSSH server (installed from Microsoft's GitHub
  after the Windows-Update FoD route returned `NotPresent`); key-based login; and the dev environment
  **inside WSL2 Ubuntu 24.04** (already present, GPU visible in WSL) — base toolchain, **CUDA 13.3**
  toolkit, **Rust 1.93.1**. The full CRUCIBLE workspace **builds on the box in ~27 s incl. hdf5-from-source**.
  Working pattern of record: edit on the laptop → `rsync -e ssh --rsync-path="wsl rsync"` into WSL → build/run
  over ssh (Windows sshd lands in cmd.exe, so bash runs via `wsl bash -l` with the script piped over stdin;
  long jobs held in the FOREGROUND of a laptop-side background ssh session — WSL kills nohup'd detached jobs
  when the launching session exits). Sudo needed no persistence hack: `wsl -u root` runs apt as root
  per-command (the NOPASSWD-sudoers route was correctly refused by the auto-mode guard as unrequested
  persistence).
  **(2) The S11 doctrine amendment (doc-first, ruling #14).** The **torch / ASI flame-holder object is
  DELETED** (reverses the S7/S8 "torch tier"): the only boundary inputs to a start are the **bounded initial
  spark** (a deposit that ends) and the **injected fuel**; flame-holding / recirculation / turbulent mixing
  are **emergent**, resolved on the grid — so the startability verdict cannot be circular. Refinement is
  re-scoped to **turbulence, not the flame front** (the S10 NO-GO measured front sharpness invisible to the
  timeline verdict): the S16 solution-adaptive subsystem (Pope-80% / vorticity / Q, 3-D incl. θ, GPU-resident,
  LES-commutation-aware), with the archived-dormant S10 (r,z) flux-register revived there as its
  interface-conservation kernel. Landed as PLAN §1 ruling #14 + §2.2/§2.3 + §5 S16 + §8 v1.12; VISION_SCOPE
  §4.1 reaction-razor worked consequence + §15 v1.6; META-3 `spark-igniter-class` 0.8.6 DEPRECATED note. The
  torch-tier **code + configs stay in place, marked DEPRECATED, until S16** replaces them with emergent
  flame-holding (ripping them out now would orphan ◆C2/◆C3 with no starter) — `rl10_startup*.toml` flagged,
  values untouched.
  **(3) The spike (STEP 2).** Ported the two hot kernels to CUDA, matching the CPU math (indices, GammaLaw
  closures, PPM edge algebra, HLLC-Batten, the Illinois root-find + its constants). **Measured f64 throughput
  (RTX 4070 Ti SUPER, CUDA 13.3, arch sm_89):**
  · STREAM triad = **401 GB/s** sustained (60% of the 672 peak) — the bandwidth ceiling.
  · PPM+HLLC hyperbolic z-sweep (NCOMP=7) = **~0.96×10⁸ cell-updates/s**, run-twice **BIT-IDENTICAL**.
  · Illinois (p,h,Z) EOS projection = **~4.9×10⁸ cell-projections/s**, root recovery 2.8e-15, run-twice
    **BIT-IDENTICAL**.
  **Determinism strategy VALIDATED** (META-1 §2.5, ruling #13): gather-only, fixed-order, no atomics ⇒
  same-build reruns bit-identical. The **Rust→CUDA binding skeleton** = `crates/gpu` (`crucible-gpu`,
  raw-FFI + nvcc `build.rs`, its own empty `[workspace]` so the laptop `check.sh` never builds it — CPU
  reference untouched, gate 5 byte-identical). Its `gpu_spike` bin cross-checks the GPU HLLC against the
  bit-exact `crucible_solvers::euler::hllc_flux` over 1,048,576 face pairs × 3 directions: **worst rel diff
  5.0×10⁻¹⁰** (the declared cross-device ECT tolerance — FMA contraction differs CPU↔GPU; bit-identity
  cross-device is impossible and not promised), GPU rerun bit-identical, PASS.
  **THE KEY FINDING (governs the envelope):** the sweep's ~10⁸ cups is **occupancy/latency-bound, NOT
  f64-ALU-bound and NOT bandwidth-saturated.** Evidence: the shared-memory formulation barely beat the naive
  fused one (9.6 vs 9.3×10⁷); the f32 twin is only **~2×** faster (would be 20–60× if f64 arithmetic were the
  wall); effective DRAM traffic is ~10 GB/s vs the 401 GB/s ceiling; `ptxas` reports **146 registers/thread ⇒
  ~25% occupancy**. So ruling #9's "bandwidth-bound" premise is **reachable** — the lever is register reduction
  (staged predictor/correct kernels, fewer live stencil arrays) at the S13 residency + S14 profiling, not
  blocked by the consumer f64 1/64 rate. **§3 THE-RUN envelope RE-SIZED with data:** ~3×10¹²–1.5×10¹³
  cell-updates ÷ the measured sweep floor ~1×10⁸ cups = **~8–42 h** (low end in-cap; high end rides the legal
  multi-night checkpointed continuation), dropping to **~1.7–14 h** at a realistic tuned ~3–5×10⁸ cups.
  **Spike simplifications (S13 generalizes, recorded in `crates/gpu/cuda/bench/README.md`):** z-sweep on a
  uniform (r,z) field N_θ=1 (z-face annular areas cancel — metric-exact for z; r-sweep face_radius + cut
  apertures + mixed N_θ are S13); the EOS surface is a synthetic **multilinear** (p,h,Z) grid exercising the
  exact trilinear gather + Illinois flow (the real HDF5 surface + the SOLV-4 0.4.8 multi-root scan guard are
  S13); fused single-kernel sweep, not the staged SDC pipeline; class-D diffusion CG unspiked. Bench kernels
  archived under `crates/gpu/cuda/bench/`. **Docs:** this entry, CLAUDE.md State, PLAN §8 v1.12+v1.13 + ruling
  #9 + §3, VISION_SCOPE §15 v1.6, META-3 0.8.6. Local `check.sh` green (the detached `crucible-gpu` is
  invisible to it); certificates byte-identical. **NEXT = S13 (full residency): entire step on-device,
  FND-6 checkpoint/restart, the CPU↔GPU tolerance cross-check on real fixtures, register-reduction toward
  the throughput target.**

- Session 25 (2026-08-28): **plan S13 — GPU RESIDENCY (the class-A step on-device) + determinism
  cross-check + the FND-6 checkpoint primitive + tuning; SPLIT — class-D/combustion/real-EOS/geometry →
  S13b.** S13 as briefed ("the ENTIRE SDC step on-device") is much larger than one ssh-remote session; per
  the improvisation rule I drew the split at the **class-A (explicit hyperbolic) subset** — a complete,
  self-consistent resident step — and recorded the rest to S13b. Branch `s13-gpu-residency` off
  `s12-gpu-spike` (the S12→main merge left to Ben — it was never pushed; S13 contains all of S12).
  **What landed, all measured on the RTX 4070 Ti SUPER (CUDA 13.3, sm_89):**
  **(1) The residency core (STEP 1).** The class-A RHS ported to a device-resident kernel set
  (`crates/gpu/cuda/residency.cu`), generalizing two of the four S12 spike simplifications: **(a)** the
  sweep is now the REAL 2-direction operator — the exact `face_radius` cylindrical metric (r-sweep
  area-weighted `af = A·F`, z-sweep metric-ratio `(F_l−F_r)/dz`) + the SOLV-1 §3.3 geometric sources
  (radial pressure + centrifugal `ρu_θ²/r̄` + swirl `−ρu_ru_θ/r̄`) — not the S12 z-only-uniform pass where
  the annular z-face areas cancel; **(c)** the kernels are **staged** (`fill_prims` → per-direction `rate`
  → `compose`), not one fused kernel. PPM (CW84 + van-Leer MC limiter) + HLLC-Batten are bit-for-formula
  copies of `recon.rs`/`hllc.rs`. GammaLaw EOS on-device (the real HDF5 `TableEos` (p,h,Z) projection is
  S13b). **Cross-check vs the bit-exact CPU `Euler::eval_rhs`** on a smooth subsonic N_θ=1 box fixture,
  INTERIOR cells only (≥ NGHOST from every edge, so the compact PPM stencil is all real interior data —
  the BC/reflux/axis machinery is out of the compared set, S13b): **worst rel 1.2×10⁻¹⁰** over 26,460
  scalar comparisons (FMA-order, inside the declared ECT ~5×10⁻¹⁰), **same-build rerun bit-identical**.
  **(2) The marched resident SDC step + determinism (STEP 1+2).** 2-node Lobatto IMEX-SDC, explicit-only
  (1 predictor + `N_SDC_CORRECTIONS`=2 corrections; `compose_gas` flow-only composition `U = u0 +
  we0·A(u0) + we1·A(U)`) — the state **lives on-device across the whole internal step loop**, CPU
  orchestrates (uploads once, downloads once). Fixed Δt fed to CPU and GPU alike (stable_dt-on-device is
  the reduction, S13b). **Cross-check vs the CPU `Sdc::step_flow`** over a 5-step march on the ≥3·M
  interior (neither the CPU's BCs nor the GPU's frozen boundary can reach it in M compact steps): **worst
  rel 3.1×10⁻¹¹ — ECT does NOT grow over the march** (phase-1 was 1.2×10⁻¹⁰); **GPU resident-march rerun
  bit-identical**. (META-1 §2.5: gather-only, one writer per cell, fixed control flow, no physics atomics.)
  **(3) FND-6 checkpoint/restart (STEP 3).** The device-state host round-trip proven **bit-faithful**:
  `march(2) ▸ resume ▸ march(3)` == continuous `march(5)` **byte-for-byte** on the same device/build — the
  §7 "one physical trajectory across wall-clock segments" made real. **FND-6 amended first** (0.4→0.5,
  §3.8): checkpoints at a **step boundary** (only the field state is live — the SDC node scratch is
  step-local, rebuilt deterministically, so not serialized); byte-identical-continuation contract
  extending §3.6 gate 1 across a wall-clock seam; fail-loud on manifest/build/table mismatch. The
  overnight auto-checkpoint/resume/halt-artifact **harness + profiling-to-target is S14 (◆C4)** — this
  session fixes the contract + the primitive.
  **(4) Tuning (STEP 4) — the register-reduction finding.** Staged `fill_prims`/`compose` out of the
  sweep (28 / 12–14 regs, cheap), then split the rate kernel **per-direction** (each holds ONE 7-cell
  pencil). `ptxas -v`: the **fused monolithic** rate kernel = **228 regs** + 288 B stack (WORSE than the
  S12 fused sweep's 146 — it fuses both directions + sources + the whole PPM pencil into local memory);
  the **per-direction split = 206 regs each** (occupancy ~21%). Throughput on 256×512 = 131,072 cells:
  **8.4×10⁷ full-RHS-evals/s ≈ 1.7×10⁸ cell-direction-updates/s** — the **same ~1×10⁸ cups regime** S12
  measured. **FINDING: naive fused-RHS residency is register-heavy, and the per-direction split buys only
  a marginal win (228→206).** The real register-reduction lever — compute each face **once** into a device
  flux buffer (retiring the current 2×-per-cell face recompute) + localize the PPM temporaries + the
  staged predictor/correct pipeline — is **confirmed as the S14 target**, not reached by the naive
  staging. **§3 THE-RUN envelope UNCHANGED:** the ~8–42 h un-tuned floor stands (now reproduced on the
  resident path), the ~1.7–14 h tuned target still gated on the S14 register work.
  **CPU reference untouched** except one read-only accessor (`EulerWorkspace::rates()`, for the
  cross-check — no physics path reads it); local `check.sh` green (the detached `crucible-gpu` is invisible
  to it); **certificates byte-identical** (gate 5). Build: `crates/gpu` gained `cuda/residency.cu` (the
  resident class-A kernels + the march + throughput-bench FFIs) + `src/residency_xcheck.rs` (the 3-phase
  harness); `build.rs` generalized to archive both `.cu`; `crucible-grid` added as a dep.
  **SPLIT → S13b (improvisation rule):** class-D implicit diffusion CG residency (gas `F_visc` +
  solid conduction + Robin-Robin Picard); combustion (Nagumo + front-thickening + class-R auto-ignition)
  residency; the real HDF5 `TableEos` (p,h,Z) projection + multi-root scan guard on-device; cut apertures
  + mixed-N_θ + SRD + BC + `stable_dt` residency (the whole-step generality). ◆C4's overnight harness +
  profiling-to-target = **S14**. **Docs:** this entry, CLAUDE.md State, PLAN §8 v1.14 + §5, FND-6 0.5.
  **NEXT = S13b (finish residency), then S14 (hardening + ◆C4).**

- Session 26 (2026-08-31): **plan S13b — GPU RESIDENCY, the class-D IMPLICIT DIFFUSION per-component
  symmetric CG on device + the fixed-topology-reduction determinism primitive; SPLIT —
  b-assembly / Robin / solid-conduction / combustion / real-EOS / whole-step geometry → S13c.** S13b as
  briefed (class-D diffusion + combustion + real `TableEos` + geometry generality) is far larger than one
  ssh-remote session; per the improvisation rule I drew the split at the **class-D per-component symmetric
  CG** — the brief's own "single biggest piece" — a complete, self-consistent resident solver, and
  recorded the rest to S13c. Branch `s13-gpu-residency` (continues from S13; the S12→main merge still
  Ben's, deferred). **Why the CG first:** the S13 class-A residency had **zero reductions** (pure gather
  stencils); the class-D CG introduces the first on-device REDUCTIONS (the dot products), the genuinely
  new and only risky determinism surface — a **fixed-topology tree reduction** (META-1 §2.5:
  "fixed-topology tree reductions for all grid/ensemble statistics"; COUP-3 §3.1: "per-component
  symmetric fixed-structure CG ... fixed-order reductions ... reduction topology coefficient-independent").
  Nail the primitive here against a bit-exact CPU reference and every downstream residency piece (Robin
  coupling, `stable_dt`, the combustion node solve) reuses it.
  **What landed, measured on the RTX 4070 Ti SUPER (CUDA 13.3, sm_89, first box trip clean):**
  **(1) The resident class-D CG (`crates/gpu/cuda/residency_diffusion.cu`).** `apply_linear` ported as a
  gather kernel — the exact cylindrical-metric two-point face stencil `coef·A·(x_nbr−x_c)/d` with the
  transport read as the two-cell **face average** (`FaceTr::between`, so a one-sided coefficient would
  show), the Uᵣ negative-definite geometric diagonal `−(4/3)μ·u_r/r̄·geo·κV` (SOLV-1 §3.3), and the
  per-component `face_coef` (Uᵣ/U_z the μ vs (4/3)μ split, ω the μ·r_face² angular-momentum form, T the
  k, C the ρD) — bit-for-formula from `gas_diffusion.rs`. `fill_mass` (ρκV; ρr̄²κV for ω; ρc_vκV for T)
  as a kernel. The Jacobi-preconditioned CG driver runs the **fields resident on-device across the whole
  solve loop** (only the O(1) scalars — α,β + the termination predicates — round-trip to host per
  iteration, exactly as a real resident CG does). **The fixed-topology reduction** = block-level pairwise
  shared-mem tree → fixed-stride block-partial pre-sum → single-block pairwise tree; shape a pure function
  of (ncell, TPB), never of scheduling. **All 5 components (Uᵣ,U_z,ω,T,C) are the ONE component-generic
  `cg_solve`** — the Picard-lagged cross terms live in the RHS `b` (`assemble_rates`), never in the CG
  matrix (module doc), so validating the five validates the whole class-D CG.
  **(2) Cross-check vs the bit-exact CPU `GasDiffusion::cg_solve`** (`crates/gpu/src/residency_diffusion_xcheck.rs`,
  via the additive doc-hidden `xcheck_cg_dense` accessor — runs the real solve, changes no production
  number, the S13 `EulerWorkspace::rates()` pattern) on a real N_θ=1 (r,z) box fixture (48×96), varying
  nowhere-symmetric transport, free (zero-flux Neumann) BCs, a nonzero smooth RHS. Driven for the SAME
  iteration count the CPU measured (identical work; only arithmetic order differs), per component
  (iters, resid → worst rel on x / on δ=x−x0): **Uᵣ 270, 8.9e-13 → 3.4e-12 / 7.6e-12; U_z 293,
  9.3e-13 → 9.6e-13 / 8.1e-13; ω 271, 9.5e-13 → 2.2e-13 / 3.4e-12; T 20, 4.9e-13 → 7.0e-12 / 2.3e-13;
  C 290, 9.3e-13 → 3.8e-12 / 2.1e-12** — **overall worst 7.6×10⁻¹²**, well inside the declared
  converged-solve ECT **1×10⁻⁸** (looser than the class-A per-op 1e-9 because the reduction-shape/FMA
  difference accumulates over the CG iterations to the residual floor). **The GPU's OWN data-dependent
  termination lands on the SAME iteration count and the same residual as the CPU for all five components
  (270/293/271/20/290)** — the reduction is tight enough that the shape difference never flips a
  termination decision (the light velocity/species masses give the stiff ~270–293-iter solves; T's
  c_v-weighted mass gives the easy 20-iter solve — the reduction is exercised hundreds of times). **All
  same-build GPU reruns bit-identical** (gather-only, one writer per cell, no physics atomics — META-1
  §2.5). The S13 class-A `residency_xcheck` re-ran **unchanged** on the same build (single-RHS 1.2e-10,
  marched 3.07e-11, FND-6 checkpoint byte-identical, throughput 1.71×10⁸ cell-updates/s).
  **(3) SPLIT → S13c (improvisation rule):** the RHS `b`-assembly residency (the Picard-lagged
  cross-stress + species-enthalpy `Σ h_k j_k` flux — `assemble_rates`); the **Robin-Robin fixed-Picard
  wall coupling** + the **solid-conduction CG** (the other two legs of brief item 1); combustion residency
  (bistable-Nagumo propagation + front-thickening diffusion + the class-R BE-with-τ-refreeze auto-ignition
  node solve); the real HDF5 `TableEos` (p,h,Z) Illinois projection + the SOLV-4 0.4.8 multi-root scan
  guard on-device (replacing the GammaLaw closure the class-A path carries); cut apertures + mixed-N_θ
  refluxing + SRD + domain BCs + `stable_dt`-on-device (the CFL reduction — the reduction primitive's next
  consumer). Acceptance for S13c stays the full-physics resident step CPU↔GPU cross-check on a coarse 3-D
  RL10 fixture. **S14** (profiling to the throughput target via the flux-buffer register reduction + the
  overnight auto-checkpoint/resume harness + ◆C4) is unchanged.
  **CPU reference untouched** except the additive `xcheck_cg_dense` accessor + one guard test
  (`xcheck_cg_dense_converges_for_all_components`), both doc-hidden/test-only, unused by production; local
  `check.sh` green (the detached `crucible-gpu` is invisible to it); **certificates byte-identical**
  (gate 5). Build: `crates/gpu` gained `cuda/residency_diffusion.cu` + the `residency_diffusion_xcheck`
  bin; `build.rs` archives the third `.cu`. **Docs:** this entry, CLAUDE.md State, PLAN §8 v1.15 + §5
  S13c, COUP-3 0.4.7 (device-resident CG breadcrumb). **NEXT = S13c (finish residency: b-assembly + Robin
  + solid conduction + combustion + real EOS + geometry), then S14 (hardening + ◆C4).**- Session 27 (2026-09-01): **plan S13c — GPU RESIDENCY, the physics remainder: stable_dt + the class-D
  diffusion FORCING + the FULL resident class-D diffusion STEP + the real HDF5 TableEos projection +
  the combustion source + the blend-EOS two-branch projection — SIX legs, all validated on the RTX 4070 Ti
  SUPER.** Continues S13/S13b on
  `s13-gpu-residency`. Ben's rulings this session drove the order (physics-first; f64-only forever;
  degenerate lookups resolve by continuity not halt — PLAN §8 v1.16). Each leg = an additive doc-hidden
  CPU accessor (runs the REAL production code, changes no number) + a device kernel set + a CPU↔GPU
  cross-check bin; the CPU reference stays bit-exact, certificates byte-identical.
  **(1) stable_dt — the CFL clock** (`residency.cu` `gpu_stable_dt` + a max-tree reduction): Δt =
  cfl / max_cell σ, σ = (|u_r|+c)/dr + (|u_z|+c)/dz. The reduction is a MAX (exactly order-independent),
  so the resident march self-clocks deterministically. **CPU↔GPU rel = 0.0 (bit-identical)**, rerun
  bit-identical (added as PHASE-1b of `residency_xcheck`).
  **(2) The class-D diffusion FORCING** (`residency_diffusion.cu` `kd_lag_grads` + `kd_assemble_rates`):
  the affine viscous-stress physics that builds the CG's RHS `b` — the full τ stress tensor
  (τ_rr/τ_zz/τ_rz + the −τ_θθ/r geometric source), viscous work, Fourier conduction, species diffusion +
  its enthalpy flux — bit-for-formula from `assemble_rates`. A pure gather; **worst rel 2.1×10⁻¹²** over
  all 5 diffusion components (lag≠sol + nonzero ∂h/∂Z exercised), rerun bit-identical.
  **(3) The FULL RESIDENT class-D diffusion STEP** (`gpu_class_d_iterate` + `kd_fill_gas_rhs` +
  `run_cg_resident`): one complete SDC-inner implicit-diffusion sweep marched ENTIRELY on-device, fields
  resident throughout — assemble(sol,lag) → {Ur,Uz,Om} CG solves → re-assemble → {T,C} CG solves. Composes
  the FORCING (2) + the S13b SOLVER + the RHS builder. Vs the CPU `xcheck_class_d_iterate_dense` (the real
  assemble/fill_gas_rhs/cg_solve in SDC-inner order): **worst rel 3.2×10⁻¹⁰** (converged-solve ECT 1e-8,
  5 CG solves chained), rerun bit-identical. **The gas class-D diffusion is now fully resident.**
  **(4) The real HDF5 TableEos (p,h,Z) projection** (`residency_eos.cu` `k_project`): the actual
  equilibrium-surface EOS replacing the GammaLaw stand-in — the fixed 8-corner multilinear interp
  (interp_rule space) + the deterministic Illinois regula-falsi projection (warm hint + cold bracket),
  the two-root case resolved by CONTINUITY (warm keeps the near root, cold takes the first scan crossing —
  Ben ruling, and the CPU already behaves so). Marshaled the production `lox_lh2_v0.4.0` surface
  (`BoundColumn::marshal`, additive). Vs the CPU `TableEos` over 245 on-surface states: **worst rel
  2.1×10⁻¹⁵ (machine precision)** on p, 5.0×10⁻¹⁶ sound, 7.0×10⁻¹⁶ temperature (the interp reduction order
  is identical CPU↔GPU); rerun bit-identical. The rare near-vacuum golden-section tangency corner is a
  documented follow-on (device returns NaN; the interior fixture never reaches it).
  **(5) The COMBUSTION SOURCE (SOLV-4.4)** (`residency_combustion.cu` `k_comb_source`): the flame-
  propagation physics — the bistable-Nagumo pushed front ρ_u·K·b(1−b)(b−a) + the matched front-thickening
  diffusion ∇·(ρD_c∇b), reading the unburnt (p,h,Z) T_u/ρ_u surfaces + the ignition (p,T_u,Z) S_L surface
  (all direct interps — the cell's p is already in the prim). Production `Combustion::accumulate`
  refactored to delegate to `accumulate_inner` so the CPU accessor runs the IDENTICAL code (zero
  divergence). Vs the CPU on the design flame state (P0=1e5, H0=4.364e4, Z0=0.167) over 1379 genuinely-
  active reacting cells (max source 33 kg/m³/s): **worst rel 1.4×10⁻¹³**, rerun bit-identical.
  **(6) The BLEND EOS (p,h,Z) two-branch projection** (`residency_blend_eos.cu` `k_blend_project`): the
  shifting-equilibrium (unburnt ⊕ burnt) projection — the mass-weighted specific-volume density closure
  1/ρ = (1−b)/ρ_u(p,h_u,Z) + b/ρ_b(p,h_b+off,Z) on the partitioned sub-state enthalpies + the same Illinois
  + first-crossing scan (two-root → continuity per ruling #3; device carries no >1-crossing refusal, the
  unique-root fixture matches the CPU). The **class-R auto-ignition prerequisite** (its per-refreeze state
  query is this projection). Vs the CPU `BurnBlendEos::prim_checked` on the production unburnt v0.3.0 ⊕
  burnt v0.4.0 surfaces over 64 mid-b (both-branch) states: **worst rel 5.9×10⁻¹⁴**, rerun bit-identical.
  **CARRIED to a follow-on:** the class-R implicit auto-ignition node solve — now a SMALL step (a fixed-
  N_TAU_REFREEZE(=2) BE-at-frozen-τ loop over the blend projection (6) + the ignition τ_ign interp), but its
  faithful port needs the full `prim_checked` branch switch (pure-unburnt / mid-b / pure-burnt) since x→cap
  crosses into the pure-burnt delegation — the pure-burnt path IS the resident TableEos (4), so it is a
  branch-select wiring job. The θ-stress tensor (N_θ>1) + cut apertures + mixed-N_θ reflux (+ real
  NoSlip/Robin walls + solid-conduction CG — off the GPU-◆C3-parity path since ◆C3 is ADIABATIC free-slip,
  the CPU's own S11 deferral); the near-vacuum EOS tangency corner. The full-physics resident step on a
  coarse 3-D RL10 (the S13c acceptance) is reached once the θ/geometry generality + class-R lands. **CPU
  reference untouched but for additive doc-hidden accessors + one refactor-extract; `check.sh` green
  (gates 1/2/3/5; gate 4 is Python-only, untouched); certificates byte-identical.** **Reusable device
  pieces now in-tree:** the 8-corner multilinear interp (interp_rule space) + the Illinois projection +
  the fixed-topology CG reduction + the max-tree; table marshaling via `BoundColumn::marshal`.
  **Docs:** this entry, CLAUDE.md State, PLAN §8 v1.17. **NEXT = the S13c remainder (class-R branch-wire +
  θ/geometry), then S14 (throughput + ◆C4).**
- Session 28 (2026-09-03/04): **plan S13c — GPU RESIDENCY, two more legs: the CLASS-R implicit auto-ignition
  node solve + the class-A step on the REAL 3-D CUT GEOMETRY (the ◆C3 world); the SOLV-4 multi-root refusal
  relaxed to continuity (doc-first).** Continues `s13-gpu-residency`. The box was unreachable for the first
  ~half of the session (Tailscale showed `backhouse` offline), so both legs were written BLIND against the CPU
  code and validated on the first box trip — the S27 pattern (additive doc-hidden CPU accessors running the
  REAL production code + a device kernel set + a CPU↔GPU cross-check bin; the CPU reference bit-exact,
  certificates byte-identical, `check.sh` gates 1/2/3/5 green).
  **(0) SOLV-4 0.4.9 — mid-transition multi-root resolves by CONTINUITY (Ben ruling 2026-08-31, PLAN §8
  v1.16 (3), confirmed in the session brief).** The v0.4.8 `>1 crossing ⇒ typed refusal` is retired:
  `scan_first_crossing` breaks at the FIRST crossing from the bracket's cold end (the pre-v0.4.8 loop
  structure — identical floats for every single-root bracket, so every previously-accepted root is
  bit-identical; `MULTI_ROOT_REFUSAL` deleted; the 7 counting tests replaced by first-crossing pins). Matches
  `TableEos`'s warm-path rule and what the device projections (S27 legs 4/6) already carried. Gate 5: all five
  certificates byte-identical (the stations never enter the mid-`c` path).
  **(7) The CLASS-R implicit auto-ignition node solve** (`residency_class_r.cu` `k_class_r`, 116 regs):
  `Combustion::implicit_auto_update` bit-for-formula — BE in ρb at frozen τ_ign, `N_TAU_REFREEZE = 2`, the
  exact cap parking, the symmetric base projection onto [0, cap], both cold floors — with THE new device
  surface: the blend's full `prim_checked` **three-way branch select** on b = x/ρ (pure-unburnt TableEos /
  mid-b two-branch blend / pure-burnt TableEos at e + h_off), which a stiff node walks base → mid → cap within
  one solve. Marshals the τ_ign column (`IgnitionColumns::xcheck_marshal_delay`, additive). Fixture: the hot
  unburnt design state (T_u ≈ 1050 K, τ = 2.41e-4 s) × b₀ ∈ {0, 0.3, 0.6} × (p, h) spread + the three guard
  corners (advected ρb past the cap, the measured −0.08ρ negative quadrature base, base == cap), the node weight
  swept w/τ ∈ {1e-3, 0.3, 1, 3, 30, 1e6}: **worst rel 4.8×10⁻¹³ (ρb) / 2.0×10⁻¹² (realized rate)** over 186
  comparisons (0 CPU-refused), final states mid-b 115 / parked pure-burnt 71, **every CPU cap-parking
  reproduced EXACTLY** (61 parkings; at w/τ = 1e6 rel 0.0), rerun bit-identical. Finding while porting: the
  S27 blend-EOS kernel shipped `H_BRACKET_MARGIN = 1e-12` where `blend_eos.rs` uses **1e-9** — latent (the S27
  fixture's bracket was envelope-dominated); corrected in that kernel and the class-R inline copy.
  **(8) The class-A step on the REAL 3-D CUT GEOMETRY** (`residency_geometry.cu`): the ◆C3 world — the
  geometry-of-record RL10 contour revolved at **uniform N_θ = 8** through `build_with_geometry_theta` at the
  ◆C3 desktop dial 5 (40×119×8 = 38 080 cells, 21 336 active, **r_min = 0**), GammaLaw. θ-plane-major dense
  layout; host-precomputed geometry-time tables (the (r,z) activity map, per-sector κ + six apertures, the
  pencil RUN tables — maximal active runs per line with the ghost KIND at each end + the wall normal there —
  and the SRD member lists in the CPU's own build order via the additive `euler::xcheck_srd_neighborhoods`).
  Kernels: `k3_rate_r` / `k3_rate_theta` / `k3_rate_z` / `k3_sources` — the r-sweep area-weighted
  `(A·ap)·F` single difference over κV, the periodic θ-sweep with the canonical θ+ aperture side and the
  `A_θ/V` ratio, the z-sweep metric-ratio form, the **FND-2 §3.2 axis parity-pair gather** (θ+π plane, u_r and
  u_θ negated), interior wall run-boundary ghosts (**slip-reflect about the contour's true wall normal**, the
  ◆C3 `slip_wall_z_faces = true` setting; grid-aligned mirror available), Reflecting/Transmissive domain
  ghosts, the SOLV-1 §3.3 geometric sources + the **per-sector embedded-interface pressure closure**
  `p·W/(κV)` (`wall_closure_cell`, the θ-limb applied only when nonzero); **State Redistribution** as two
  gathers (`k3_srd_q` per small neighbourhood, `k3_srd_apply` per affected cell over the host-fixed inverse
  map, the CPU's member/accumulation order); the **3-D `stable_dt`** with the θ-arc member (max-tree); the
  **resident march with stagewise SRD** (compose → SRD after every composition, `Sdc::step_flow`'s sequence).
  Every ghost rule is on-device, so the cross-check compares **ALL active cells** (149 352 scalars): **RHS
  worst component-scaled rel 7.6×10⁻¹⁵** (the per-cell rel 2.4×10⁻⁹ is the well-balance cancellation — a
  4×10⁻⁶ net r-momentum rate against a 2×10³ component max, an FMA-order absolute 5.8×10⁻¹¹; the harness now
  reports both measures, the ECT is the component-scaled one, a 1e-6 per-cell rail stands as sanity);
  **`stable_dt` bit-identical** (rel 0.0); **SRD 3.3×10⁻¹⁶** (8160 scalars moved, 680 small sector-cells);
  **5-step march + SRD 1.1×10⁻¹⁵** (per-cell 4.2×10⁻¹²), resident-march rerun **bit-identical**, checkpoint
  `march(2)▸march(3) == march(5)` **bit-identical**. **Throughput (measured, un-tuned): 3.0×10⁷ cell-RHS/s** on
  this 21k-cell world — the 3-D r/z rate kernels compile to **254 registers + 488 B stack spill** (θ 212 regs):
  the fused 7-pencil PPM + HLLC + ghost logic per cell is exactly the register-heavy shape S13 measured (228)
  and the small world under-occupies the card; the S14 flux-buffer lever (each face once into a device flux
  buffer + a gather divergence, localized PPM temporaries) is now measured to matter on the 3-D kernels too.
  Scope as agreed with Ben (session 28): interior + reflective cut walls + Reflecting/Transmissive domain BCs;
  the injector mass-flow inflow + pressure-outflow ghosts, the general-EOS HLLC aux-slot path (blend), and
  the COUP-2 ledger reductions ride the composed-step session. Regression: the seven S13–S27 cross-check
  bins re-run on the box (see the closing note below).
  **REMAINING toward the S13c acceptance (the full-physics resident step on the coarse 3-D RL10) — now a
  COMPOSITION job, every physics piece being resident:** wire the blend EOS (S27 leg 6 + the three-way select
  of leg 7) into the 3-D `fill_prims` + the general-EOS HLLC (aux e/Γ₁ slots, `roe_sound_speed`), the
  combustion source (leg 5) on the 3-D layout (θ-faces + per-sector κ), class-R (leg 7) after each accepted
  composition, the igniter deposit, the inflow/outflow ghosts, `stable_dt`'s front-carrier member, and the
  ledger/audit + reacting-measure reductions; then the CPU↔GPU cross-check on `rl10_startup_3d`'s own march.
  Walls (NoSlip/Robin + solid conduction) stay off the ◆C3-parity path. Then **S14** (flux-buffer register
  reduction + the overnight harness) and **◆C4**.
  **CPU reference: the SOLV-4 0.4.9 continuity change (a ruling, doc-first) + additive doc-hidden accessors
  (`xcheck_marshal_delay`, `xcheck_srd_neighborhoods`); `check.sh` gates 1/2/3/5 green; certificates
  byte-identical.** **Docs:** this entry, CLAUDE.md State, PLAN §8 v1.18 + §5 S13c, SOLV-4 0.4.9.
  **Closing regression (box):** all seven S13–S27 cross-check bins re-run green after the session's changes —
  class-A 1.2×10⁻¹⁰ / marched 3.1×10⁻¹¹ (8.0×10⁷ cell-RHS/s), class-D CG 7.6×10⁻¹², forcing 2.1×10⁻¹²,
  resident diffusion step 3.2×10⁻¹⁰, TableEos 2.1×10⁻¹⁵, combustion 1.4×10⁻¹³, blend EOS **1.3×10⁻¹³** (moved
  from 5.9×10⁻¹⁴ by the 1e-9 bracket-margin correction — the scan grid shifts with the bracket; still the
  machine-precision class), every rerun bit-identical.
- Session 28, continued (2026-09-04/05): **S13c CLOSE — the COMPOSED full-physics 3-D resident ENGINE step
  validated on the real `rl10_startup_3d` assembly; S14 — the overnight harness (FND-6 checkpoint/resume,
  COUP-2 audit every step, COUP-4 verdict probes) built + smoke-tested; ◆C4 launched.** Ben's directive:
  "just get the full GPU port done", then "kick off a longer run".
  **(9) The composed engine step** (`crates/gpu/cuda/residency_engine.cu` + the host bridge
  `crates/gpu/src/engine_host.rs`, the crate's new lib half): ONE device-resident SDC step behind a persistent
  handle composing every validated leg — the blend EOS `prim_checked_hinted` (three-way branch select, warm-start
  hints on the pure limits kept RESIDENT in the prim cache, the mass-weighted blended sound speed → the e/Γ₁ aux
  slots), the general-EOS HLLC-Batten (Γ₁ slot, datum-free Roe c²), the r/θ/z sweeps on the cut geometry with the
  axis parity pair + slip walls, the COUP-7 **injector mass-flow inflow ghost** (the unburnt TableEos
  fixed-count solve, sonic-capped) and the **pressure-outflow ghost** (the pump-down ambient), the geometric +
  wall-closure sources + the sector-gated ramped **igniter deposit**, the SOLV-4.4 combustion source on the 3-D
  layout (θ-faces, per-sector κ, the loud [0,1] guard), the 2-node Lobatto IMEX-SDC composition with stagewise
  SRD and the **class-R node solve after every accepted composition** (trapezoid base, realized-rate roll-over),
  the **COUP-2 ledger** (port + source, net + gross) and the κV-weighted stored totals as fixed-topology tree
  reductions, `stable_dt` with the θ-arc + **front-carrier** members and the scale-separation guard. The audit
  CHECK runs host-side on the reduced operands through the CPU's own `Sdc` arithmetic (additive doc-hidden
  `Sdc::xcheck_audit_flow` — `audit_check` refactored into a cell-count-taking body, behaviour-preserving). The
  run SCHEDULE (`crucible_engine::run`'s injector ramp, pump-down ambient, igniter kernel, audit reference
  scales) is replicated ONCE in `engine_host::Schedule` — the cross-check scores that replication against the
  CPU `Sdc::step`, the harness drives it.
  **Cross-check (`residency_engine_xcheck`)** on the real assembly (40×119×8, tables of record, the run's own
  BCs/walls/igniter), the CPU pre-marched through the real start to **t = 3.0 ms** (2306 steps, 0.277 s/step
  on the box's 28 cores; the spark at 2 ms, reacting measure 2.5×10⁻⁸ kg/s, max b 1.0×10⁻⁴ — the reacting
  regime; the state saved so reruns skip the 10-min pre-march): **PHASE 0 prims** (cold both sides) worst
  component-scaled rel **8.5×10⁻¹²** (= the projection tolerance class); **PHASE 1 full RHS** worst
  component-scaled rel **9.5×10⁻⁹**, at the θ-momentum of an AXIS cell beside the spark kernel — diagnosed as
  the (p,h,Z) projection's own 1e-11 tolerance (FMA-distinct Illinois paths on the two sides) divided by the
  ~10⁴ well-balance cancellation of the θ-flux difference × 1/(r̄Δθ) at the axis ring: **the composed step's
  honest cross-device ECT is declared 1e-7 component-scaled** (prims stay at 1e-9; per-cell rel is diagnostic
  only — at a cancellation cell it is ECT × comp_max/|value| by construction); the ledger (port/src net+gross,
  gross-scaled per COUP-2 §3.1.1) **1.8×10⁻¹⁵**; rerun bit-identical. **PHASE 1b `stable_dt`** (θ-arc + front
  carrier) rel **1.6×10⁻¹⁶**. **PHASE 2 — 10 audited coupled steps**: the final state worst component-scaled rel
  **5.2×10⁻¹¹** (the RHS's axis discrepancy does NOT accumulate — the projection re-anchors every stage),
  **every GPU step passes the CPU's own COUP-2 identity** (worst |Δ−applied|/tol 4.2×10⁻⁵), resident-march rerun
  **bit-identical**, the checkpoint split (cons + the prim hints) **bit-identical**. **PHASE 3 throughput:
  0.137 s/step** at 3 ms (21 336 active gas cells × 8 sectors) — only ~2× the 28-core CPU: the cold mid-b blend
  scan (64 nodes × two branches, always-first by the SOLV-4 0.4.7 rule; the class-R re-projects every live cell
  twice per sweep) dominates once class-R seeds b ~ 1e-5 across the hot region — **the throughput lever for the
  composed step is a warm-started mid-b projection (a SOLV-4 amendment, accuracy-neutral) + the S14
  flux-buffer register work**, both recorded, neither taken this session.
  **(10) The overnight harness (`gpu_engine_run`, S14):** the resident engine marched under the run's schedule;
  the COUP-2 audit every step (a violation halts); the COUP-4 verdict probes at the run's `PROBE_EVERY` cadence
  on the downloaded state through the engine's OWN probe functions (`plane_mdot`, `plane_thrust`,
  `injector_end_stagnation_p`, `reacting_measure`) — ignition floor, NEVER_IGNITED, FLAMEOUT, dwell → WORKS,
  FAILED_TO_REACH; **FND-6 §3.8 checkpoints** (the dense state + the warm-start prim cache + the clock + the
  verdict trackers + the pinned-input digests: config sha256, the three table pins, the world shape), written
  atomically, `--resume` verifying the digests (a mismatch REFUSES); halt artifacts (`fields.csv` /
  `crash_fields.csv` via `run::fields_csv`, `verdict.txt`, `progress.csv`). **Smoke (box):** 600 steps in 22 s;
  a 300-step run resumed to 600 reproduces every probe of the uninterrupted run and its **checkpoint is
  byte-identical** (`cmp`) — the FND-6 continuation contract across a real process restart.
  **Box mechanics (docs/gpu-box.md updated):** the `powershell Start-Process wsl` recipe does NOT run scripts
  (measured: an inline `bash -c 'sleep 100; echo …'` survives, a script path never executes, in every quoting
  form); **`tmux` inside WSL is the long-job mechanism** — a tmux session outlives the launching ssh session
  (a 40 s CUDA job proved it), the ◆C4 run is launched that way and polled by log.
  **◆C4 LAUNCHED:** `gpu_engine_run configs/rl10_startup_3d.toml` (the ◆C3 desktop spec: dial 5, N_θ = 8, the
  bounded spark, 30 flow-throughs = 84.06 ms horizon, checkpoint every 2000 steps) under tmux on the box; the
  GPU trajectory tracks the CPU pre-march to the printed digits (R = 2.13×10⁻⁹ at step 1000 on both).
  Outcome: see the closing note.
  **(11) The mid-transition WARM START — SOLV-4 0.4.10 (doc-first), the composed step's measured EOS lever:**
  the S7 rule "the mid-`c` projection ignores the hint (front cells are a thin minority)" is measured false in
  a lit chamber — class-R seeds `c ~ 10⁻⁵–10⁻⁴` across every cell whose closures are live, and `c ∈ (10⁻⁹,
  1 − 2·10⁻³)` is the mid-`c` path, so the majority of the ◆C3 world's active cells ran the always-first 64-node
  × two-branch cold scan every projection. Amended: the mid-`c` projection takes the cell's previous projected
  pressure under EXACTLY `TableEos`'s warm-path uniqueness precondition (full-bracket straddle ⇒ Illinois on
  the tight bracket `[p/1.05, 1.05p] ∩ admissible`, or on the full bracket; else the cold first-crossing scan
  unchanged); the root-residual acceptance applies to the warm root; `HINT_SPREAD` hoisted to one owner. **An
  acceleration, never physics**: CPU test `mid_transition_warm_start_is_a_pure_acceleration` — 300 mid-b probes
  (b 1e-6..0.5, hints ±20 % around the cold root) agree with the cold root to **worst rel 3.0×10⁻¹²**;
  certificates byte-identical (gate 5 — the stations never enter mid-`c`). Device: `prim_blend`'s mid path takes
  the resident prim cache's hint under the same rule; the composed-step cross-check re-run from the saved 3 ms
  state **ALL PASS** — RHS 9.5×10⁻⁹ (unchanged: the axis cancellation), **10 audited steps: state 4.0×10⁻¹²**
  (tighter than 5.2×10⁻¹¹: with hints live on both sides the two projections re-anchor together), rerun +
  checkpoint bit-identical. Its throughput gain is measured after ◆C4 releases the GPU (the shared-card
  number is not a measurement). Build plumbing: the crate's new lib target made the per-leg bins' CUDA
  symbols drop out of the link (an rlib surrenders only the objects it references) — `crucible_gpu::
  keep_kernels()` (a function reading a table of every per-leg entry point) is called once from each bin.
  **(12) The S14 FLUX-BUFFER sweep path** (`residency_engine.cu`, `CRUCIBLE_GPU_FLUXBUF=1` /
  `DeviceEngine::set_fluxbuf`): each face's PPM + HLLC computed ONCE — a cell computes its LOW face from a 6-cell
  pencil, component by component with ~12 live scalars (the same formulas as the fused `pencil_faces`: the low
  face of cell c+1 IS the high face of cell c), run-end cells also their ghost-side high face — into a device
  flux buffer; a gather-divergence kernel forms the single well-balanced difference + the ledger ports.
  **ptxas: `k_face_rz` 88 registers (from the fused 230 + 720 B spill), `k_div_rz` 54, θ 134/40.** Cross-path
  contract measured: the two device paths are NOT bit-identical (nvcc contracts FMAs per expression context —
  the scalar-per-component PPM vs the array form; 7.6×10⁻¹⁵ component-scaled, ledger 1.4×10⁻¹⁶), so the
  contract is bit-identical WITHIN a path (same-build reruns) and ECT-class ACROSS paths — a build variant,
  like CPU↔GPU; against the CPU both paths score the identical 9.5×10⁻⁹. Its throughput is measured after
  ◆C4 releases the card (shared-card numbers are not measurements). Also added: env-gated per-kernel-group
  profile hooks (`CRUCIBLE_GPU_PROFILE=1`, cudaEvent pairs, dumped at destroy) and a `CRUCIBLE_NVCC_FMAD=0`
  build switch (IEEE per-operation arithmetic, no contraction — the cross-device-ECT experiment) for the
  post-◆C4 profiling session.
