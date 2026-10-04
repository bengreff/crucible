# Architecture Audit and Strategic Plan

Audit of all 18 design documents (~5,750 lines) for physics errors, internal contradictions,
computational infeasibility, missing specifications, and 12-month implementation risk.
Quantitative claims verified numerically. Date: 2026-07-04.

**Verdict:** The architecture is unusually coherent for its scope — the uniform gravity rule,
physics/display separation, and table-based subsystems are sound and mutually reinforcing. But
there are **7 architecture-level defects** (three of which silently break the determinism
guarantee, the single most load-bearing invariant), ~10 localized physics errors, and a 12-month
timeline oversubscribed by roughly 3×. All are fixable on paper now; several would be
catastrophic to discover in month 8.

---

## TIER 1 — Architecture-breaking defects

### 1.1 The radiation feedback loop breaks bit-determinism by design — **RESOLVED 2026-07-04**
**Fix applied:** deferred-apply protocol specced in radiation_shielding.md ("Deterministic
scheduling protocol"): trigger classes (ACUTE/GEOMETRY/DRIFT/PERIODIC) with fixed sim-time
apply delays and primary counts; jobs snapshot inputs at the trigger tick; results bind at
predetermined sim time T_apply (tick loop blocks if late — wall-clock cost only).
determinism.md and ship.md §7 updated to match. Original finding retained below.
`radiation_shielding.md:76-78` states the background thread's *timing* "may vary" and the tick
loop "always uses whatever cached map is currently available," with
`sleep(check_interval) // 1 second real-time`. But the dose map feeds `damage_dpa` →
`yield_strength(T, dpa)` → structural failure → trajectory (ship.md Phase 3b/3f). So **which
tick a new map lands on depends on wall-clock speed → hardware changes the trajectory** —
contradicting determinism.md's table ("Running on faster/slower hardware → None"). This is the
main coupling path of the flagship research tool, not an edge case.

**Fix:** schedule recomputation in *sim time*: triggers fire at deterministic ticks, runs use
fixed primary counts, and results apply at a fixed sim-time delay after the trigger (trigger at
tick N → map active at tick N+K, K deterministic). The MC still executes on a background
thread — the tick loop blocks only if it reaches N+K before the result exists. Determinism
preserved; wall-clock only affects throughput.

### 1.2 No deterministic ephemeris scheme exists for planets — including Earth
Three specs contradict each other:
- physics.md §2.1 table — planets feel "host star + sibling bodies (**local N-body** within the
  planetary system)"
- generation.md §1 — every position is "a pure function of its seed and T. No integration history"
- SPICE DE441 covers roughly ±14,000 years; the simulation window is **10 million years**.

Local N-body means integration history, which means results depend on *when the system was first
instantiated* — a player who flies past a system at T=2 Myr and returns at T=5 Myr gets different
planet positions than one arriving cold at T=5 Myr. That violates query-independence. And nothing
specifies where Earth is at T = 1 Myr — which the core narrative ("returns home, millions of
years later") requires.

**Fix:** all planetary motion analytic: Keplerian elements + precomputed secular theory
(precession rates for ω, Ω), with the host star's reflex wobble added analytically (this also
rescues the "solar system barycenter" validation target). Solar system: SPICE inside its
validity window, blended handoff to the analytic representation beyond it — document that
positions beyond ~10 kyr are self-consistent fiction (they must be: the real solar system is
chaotic with ~5 Myr Lyapunov time).

### 1.3 Time-dependent spiral modulation breaks the seed architecture
generation.md §3.2 computes per-cell star counts as N_cell = ∫n dV with spiral modulation
S(R, φ − Ω_spiral·T) — **explicitly a function of T** — while §3.1 requires "same cell always
produces the same stars" from time-invariant cell coordinates. Both cannot hold: a cell's count
would change as arms sweep past it, breaking the order-independence test in determinism.md.

**Fix:** a spiral density wave *is* stars on epicyclic orbits with correlated phases (Lin-Shu
kinematic wave). Keep N_cell time-invariant (axisymmetric density only); encode the arm as a
seed-derived correlation between a star's epicyclic phases and the arm phase at its guiding
center. The overdensity then propagates at pattern speed automatically, with zero determinism
cost. Needs a worked derivation (a small research note in itself).

### 1.4 "Self-consistent by construction" is currently false
generation.md §11.5 claims the density sourcing Φ is the same density that sets star counts. It
isn't: gravity comes from the evolved SCF particle distribution; star counts come from *analytic*
McMillan profiles + Cox & Gomez spiral. The SCF's dynamically-emerged bar/spiral will NOT
coincide with the analytic loci — stars get generated in arms displaced from the arms the ship's
accelerometer feels. **Fix:** evaluate ρ directly from the SCF basis expansion (same coefficients
give the density) and use *that* for star counts, with analytic profiles only as
IMF/population weighting.

### 1.5 The SCF plan as specified cannot represent the Milky Way
Three compounding problems (physics.md §3.1):
- **Basis mismatch:** the Hernquist-Ostriker basis is spherical-harmonic; representing a cold thin
  disk (h_z/R ≈ 0.036) needs l ≳ 40–50, i.e. tens of thousands of coefficients, not 500. With 500
  you get a smeared disk, wrong vertical forces, wrong κ/ν frequencies — which the epicyclic star
  orbits are built on. The 46 MB / 3-5 μs budgets rest on the wrong count.
- **Initial-condition problem:** free SCF evolution will *not* end at the observed MW (bar at 27°
  with Ω=38, arms at observed loci, LMC at 49.5 kpc, warp phase). Matching an N-body model to
  today's MW is made-to-measure modeling — an open research field, not a 30-minute precompute.
- **Shot noise:** 10⁶ particles → percent-level coefficient noise = spurious time-dependent force
  noise larger than the claimed accuracy targets (rotation curve ±3 km/s, >99.7% orbit accuracy).

**Fix (recommended):** Φ = analytic axisymmetric MW (McMillan 2017, matches observations by
construction) + rigidly rotating analytic bar + Cox-Gomez spiral + warp, pinned to observed
parameters at T=0 — and reserve basis-function machinery (EXP or AGAMA, which have
disk-appropriate bases) for the perturbation that genuinely can't be done analytically: the LMC
wake and halo response. Over the 10 Myr window (bar turns 22°; a star covers ~4% of an orbit),
analytic time-dependence is defensible and the accuracy claims become checkable. "One force rule
everywhere" is fully preserved — only the offline construction of Φ_mean_field changes.

### 1.6 hawking.md did not describe the paper's engine — **RESOLVED 2026-07-04**
Original finding (partially misdiagnosed): hawking.md's luminosity column was 10³–10⁶ below
both the Page formula and black_hole_paper's own numbers, and inconsistent with its own
lifetime column. On reading black_hole_paper (draft 6), the paper itself is Page-consistent
and self-consistent: P = ℏc⁶f/(G²M²) with species-summed f = 3.503×10⁻³ gives 60,200 TW at
10⁹ kg, thrust F = β·0.83·P_H/c = 55.3 MN with Usov pair exhaust at β = 0.332. The audit's
waste-heat objection was also wrong — the CFL shell's gap decoupling + three-layer firewall
limits off-axis EM leakage to < 1 kW by design; no radiators exist or are needed.

The real defect: hawking.md was a stale summary that contradicted the paper on numbers
(10³–10⁶), mechanism ("thrust reflector" / capture-fraction control — the reflection approach
the paper explicitly rejects), controls (an "on/off" switch that Hawking physics doesn't
permit), geometry (called the 6×10⁷ kg, 2 m shell "sub-voxel"), and failure modes (missed
confinement loss — lateral instability τ ≈ 2.8 s — the actual catastrophic mode).

**Fix applied:** hawking.md rewritten from the paper (design space 5×10⁷–10⁹ kg, f(M) species
sum, Usov equilibrium β per design point, punch-through/stopping-margin tables, feed-based
throttle with months-scale time constants, no off switch, confinement-loss failure mode,
drive assembly as real components). ship.md Ship 3 components/characterization/validation
updated to match. black_hole_paper is the authoritative source for all BH-drive physics.

### 1.7 The relativistic-cruise gravity formulation is invalid at its marquee regime
The scheme is Newtonian ∇Φ + EIH 1PN corrections applied to a relativistic-momentum state
(physics.md §2.5, SCOPE). EIH is a **slow-motion** (v≪c) expansion; at 0.3–0.99c its velocity
terms are not small. Known limit: ultrarelativistic deflection in a weak field is 2× the
Newtonian value (the light-bending factor) — plain −γm∇Φ or Newtonian+EIH is wrong by up to 2×.
The selling point is regime *bridging*, and this bridge is broken in the middle.

**Fix:** for a test particle in a weak field, use the geodesic of the linearized metric, exact in
velocity: dp/dt = −γm[(1+β²)∇⊥Φ + ∇∥Φ] (plus gravitomagnetic terms where needed). Barely more
code than specced; exact for photons; Newtonian at low v; reduces to the current scheme at v≪c.
This formulation itself is a publishable engineering note. Keep EIH for the slow-massive-body
regime (Mercury); use the weak-field-geodesic form for the fast-ship regime.

---

## TIER 2 — Localized physics errors (verified numerically)

1. **γ³ thrust transform is backwards** (ship.md Phase 1b: "F_coord = γ³ × F_proper").
   Longitudinal 3-force is boost-**invariant**: in the dp/dt formulation, coordinate thrust =
   proper thrust (transverse divides by γ). The γ³ belongs to the acceleration relation
   (a = F/γ³m), already handled implicitly by the momentum formulation. As written, thrust at
   γ=2 is 8× too strong.
2. **ISM drag misses momentum-flux γ²** (ship.md 1d: F = ρv²A). Relativistic momentum flux of
   stopped ISM is γ²ρ₀v²A — at 0.9c a factor 5.3.
3. **NS GR-transition radius**: physics.md §3.6.2 says 100 R_s ≈ 1200 km; actually
   R_s(1.4 M_☉) = 4.14 km → 414 km (1200 km used the 12 km NS radius, not R_s).
4. **BH drive self-gravity**: hawking.md says 10⁸ kg at 100 m gives ~10⁻⁵ m/s²; actual
   6.7×10⁻⁷ m/s² (15× off).
5. **Leak/hatch flow formulas dimensionally wrong** (atmosphere_ship.md:
   `ṁ = C_d·A·P·√(2/(ρRT))` and `flow = C_d·A·(P_a−P_b)/√(2ρ)` — neither yields kg/s). Use
   choked-orifice flow ṁ = C_d·A·P₀·√(γ/RT)·(2/(γ+1))^((γ+1)/(2(γ−1))) for vacuum leaks and
   ṁ = C_d·A·√(2ρΔP) subsonic.
6. **Chen–Kipping inverted in the Jovian regime** (generation.md §4.5): M(R) with exponent 0.01
   is degenerate — every giant gets ~0.3 M_jup, contradicting §4.3 Step 2 which correctly draws
   giant mass from dN/dm ∝ m⁻¹·¹. Rule: giants mass-first → radius; small planets
   radius-first → mass.
7. **Structural math inconsistent between docs**: ship.md 3b says "Stress = stiffness_matrix ×
   load_vector" (dimensionally wrong); structural.md says "stress = K_reduced⁻¹ × load" (that's
   displacement). Correct pipeline: u = K⁻¹F, then σ per element via stress-recovery matrices —
   which must be added to the FEA export spec or failure cannot be evaluated at all.
8. **Reactor model contradiction**: physics.md §1.2 specifies a *runtime* point-kinetics ODE;
   nuclear.md/SCOPE say runtime is a static table (thermal is "the only runtime ODE"). A static
   steady-state table cannot produce rod-insertion transients, xenon behavior, or startup
   dynamics. Recommendation: allow the point-kinetics ODE at runtime (6 tiny variables,
   deterministic, ~μs) and amend the "only ODE" claim.
9. **GR handoff contradiction**: determinism.md says at 100 r_s "both formulations agree to ~1%"
   *and* "forces are continuous… no discontinuous jump." A 1% force step is a discontinuity, and
   it kills the "energy conserved to machine precision" target. Push the handoff out to where
   mismatch is below integrator tolerance (~1000 r_s) or blend over the hysteresis band. Also:
   **total energy is not conserved in a time-dependent mean field** (rotating bar, evolving SCF)
   — that validation target is ill-posed as stated; use the Jacobi integral in the bar frame or
   restrict the test to a static-field configuration.
10. **Novelty framing error in research #1**: "the transition from Bethe-Bloch to hadronic
    cascade is uncharacterized" — 100 MeV–few GeV protons are the *best*-characterized regime in
    radiation physics (the GCR peak; HZETRN/Geant4 validated exactly there). The genuinely
    uncharacterized part is **time-resolved shield–flux co-evolution under sustained relativistic
    bombardment** (geometry/damage feedback). Reframe before a referee does.

---

## TIER 3 — Missing specifications that will block implementation

1. **Flare/event schedules don't scale**: generation.md §3.4 defines flare times as a cumulative
   sum of exponential draws — an active M-dwarf accumulates ~3.6×10⁹ draws over 10 Myr; "is it
   flaring at T?" is O(N). Needs hierarchical interval seeding (seed per time bucket → Poisson
   count per bucket → positions within bucket): O(1), still deterministic.
2. **Guiding-center query unwinding**: disk stars drift Ω·T ≈ 2.4 kpc in azimuth over 10 Myr
   (verified). The spec's drift buffer (~225 pc, dispersion only) misses this: queries must map
   φ_query → φ_guide = φ − Ω(R)·T per radius bin, then buffer by cell shear ΔΩ·T·R + epicyclic
   amplitude. Also determinism.md's `star_position` omits the φ_guide phase constant.
3. **Gaia count off ~100×**: "all stars within 500 pc (~500K)" — at the project's own validation
   density (0.1 pc⁻³) that sphere holds **~52 million** stars. Decide: full 100 pc completeness
   (~331K, GCNS-like) + sampled shells beyond; re-derive the 500 MB budget.
4. **Save/load loses octree damage**: ship.md §8 saves per-component scalars, but the transport
   engine mutates *voxel-level* state (ablation profiles, per-voxel dpa, transmutation sources).
   Either serialize octree diffs or constrain damage to per-component parameterizations (e.g.,
   ablation-depth field per surface) that regenerate the octree exactly.
5. **No ground contact model**: launches start on a pad; Falcon 9 has landing legs; no doc
   specifies surface contact forces, static friction, or touchdown. A minimal pad constraint +
   leg spring-damper spec is needed for Falcon 9 validation.
6. **No attitude-control / autopilot spec**: interface.md exposes `attitude prograde`,
   `maneuver <dv> at <time>`, but nothing defines control laws, actuator allocation (wheels vs
   RCS), or the trajectory predictor behind maneuver planning (display-path, needs its own
   integrator budget).
7. **Timewarp validity limits unstated**: attitude uses explicit Euler (invalid when ω·dt ≳ 1);
   thermal linearization at dt = 10⁶ s needs a fixed iteration policy; a star system can be
   entered *and* exited inside one 10⁷ s step (needs swept-sphere checks along the step, not
   endpoint refresh); warp in LEO is compute-bound (~10³–10⁴× at dt≈10 s and 40 μs/tick —
   "unlimited" is only true in smooth regimes). Specify per-regime dt caps and effective warp
   ceilings.
8. **NRLMSISE-00 inputs**: F10.7/Ap must be pinned to fixed deterministic values (or a
   seed-derived synthetic solar cycle) — live indices don't exist at T = +1 Myr.
9. **physics.md is a Python-era fossil**: it prescribes scipy/numpy/astropy/SpiceyPy *at runtime*
   (§5), "scipy.sparse handles the linear algebra," the nrlmsise00 pip package, and "wrap scipy's
   stepper" — contradicting SCOPE's "everything else is custom Rust, two C FFI deps." One full
   reconciliation pass needed; also fix the research-paper numbering clash (physics.md §4.2
   "paper #2" vs SCOPE's laser sail #2) and SCOPE's "~55 ns each, 10,000 per second" arithmetic.
10. **Transport-table resolution risk**: ~1 KB/entry for double-differential secondary
    distributions is very coarse against a 10%-of-Geant4 acceptance bar; budget for adaptive
    refinement after the first slab benchmark. The ablation criterion (energy-per-batch > h_vap)
    omits the thermal-relief coupling SCOPE's research #1 promises — the tool needs the
    MC→thermal→ablate macro-timestep loop, with an explicit fluence↔time mapping, in its spec.

---

## TIER 4 — Feasibility

Per-tick budgets, memory estimates, and table sizes are individually credible (arithmetic
checked; only the SCF eval cost is suspect, per 1.5). The infeasibility is in the **build plan**:

- The 12-month timeline packs: dd-arithmetic + custom DOP853 + galaxy pipeline + procedural
  generation + Kerr integrator + a *validated* MC transport crate (its own doc budgets 2 months
  for what MCNP-class teams spend years validating) + a characterization pipeline that
  auto-builds FEniCS FEA and CFD models from component graphs (auto-meshing arbitrary primitive
  assemblies is a project by itself) + 3 ships + TUI + 3–4 papers. Realistic single-developer
  estimate: 2.5–3×.
- Geant4 element tables: ~330M events, "1–10 days on a small cluster" — cluster access must be
  lined up, or a laptop-scale plan B (fewer elements/bins first).
- Research novelty needs repositioning: the laser sail (#2) has a dense Breakthrough Starshot
  literature (sail stability & thermal limits), and ACMF (#3) has ICAN-II/AIMStar heritage —
  both are *revisit-with-modern-tools* papers, not blank fields. The genuinely open lanes:
  (a) shield–flux co-evolution at relativistic speeds, (b) deterministic galaxy-scale procedural
  dynamics (kinematic-density-wave generation), (c) ship-systems characterization in curved
  spacetime.

---

# STRATEGIC PLAN (revised)

**TIMELINE RE-BASELINED 2026-07-04 — SCOPE.md "TIMELINE" is now canonical.** Planning basis
agreed with Ben: scope fixed / time flexes; ~20 h/wk + AI-assisted implementation; success
bar = working validated simulator (Gate 5); laptop-only compute (staged Geant4 tables +
LMC run start as background jobs in month 1). Result: ~15 months to Gate 5, band 13-20.
The phase month-ranges below are superseded by SCOPE.md's gate-based schedule (Phase 1:
months 1-2, Phase 2: 3-5, Phase 3: 5-9, Phase 4: 9-12, Phase 5: 12-15). The three
non-compressible constraints that set the floor: developer verification hours, debug cycles
against physical ground truth, offline compute wall-clock.

**Governing principle:** the determinism kernel and the gravity spine are the two things
everything else stands on — build them first, validate continuously, and nothing merges that
breaks their test suite. Research papers are sequenced so each produces an artifact the
simulator needs anyway.

### Phase 0 (2 weeks): Spec repair — no code
Fix on paper what this audit found: rewrite physics.md for Rust runtime; adopt the
weak-field-geodesic force law (1.7); analytic-ephemeris decision (1.2); sim-time-scheduled
radiation protocol (1.1); kinematic density-wave generation scheme (1.3/1.4); mean-field
decomposition decision (1.5); reconcile hawking.md against black_hole_paper (1.6); the Tier-2
errata; the Tier-3 specs (ground contact, autopilot, warp limits, flare seeding, save format,
Gaia sampling, NRLMSISE inputs, table resolution). Cheapest two weeks of the year.

### Phase 1 (months 1–3): Determinism kernel + gravity spine
- dd arithmetic, libm policy, fixed-order summation, seed/hash infrastructure — **with
  cross-platform + query-independence CI tests from day one** (the regression oracle for
  everything after).
- Custom DOP853 with momentum state and the weak-field-geodesic force law; analytic MW potential
  (McMillan) + bar + spiral + warp; SPICE ingest + analytic handoff.
- Gates: Mercury 43″/century, S2 precession, rotation curve 233±3 km/s, bit-identical x86/ARM
  trajectory.
- **Research thread A starts:** ISM cascade Geant4 runs (offline, independent of engine code).

### Phase 2 (months 3–6): Ship kernel + Falcon 9
- Component graph → octree → thermal network → structural (with stress recovery) →
  resources/electrical/atmosphere; ground contact + basic attitude control.
- Characterization v1 uses **analytic correlations only** (Cantera equilibrium + isentropic
  nozzle, Barrowman/handbook aero, handbook conductances). FEniCS/CFD deferred to v2 — the
  single biggest schedule lever; Falcon 9 validation targets don't need CFD fidelity.
- Gates: Falcon 9 launch-to-LEO validation table; determinism suite green.

### Phase 3 (months 5–8, overlapping): Transport engine + NTP ship
- Radiation transport crate against its benchmark ladder (slab → GCR sphere → multi-layer →
  dynamic ablation), with the sim-time deterministic protocol from day one and the
  MC→thermal→ablation macro-loop.
- NTP ship with runtime point kinetics + shadow shield.
- **Paper 1 (shield–flux co-evolution / shield lifetime at 0.3–0.9c)** from the dynamic
  benchmark + Phase-1 Geant4 tables — target month 8–9. Strongest novelty claim, already on the
  critical path.

### Phase 4 (months 7–10): Galaxy + relativistic regime
- Procedural generation with the kinematic-wave scheme; Gaia overlay (100 pc complete + sampled
  shells); query system.
- SR cruise + Kerr integrator with tolerance-matched handoff radius; BH drive ship with corrected
  Hawking tables.
- **Paper 2 (deterministic kinematic-density-wave procedural galaxy, or the velocity-complete
  weak-field force formulation)** — whichever Phase-0 derivation proved richer.

### Phase 5 (months 10–12): Integration
- Three ships through all regimes; full SCOPE validation table; multi-object/docking; interface
  polish.
- **De-scope levers, in order:** biosphere Layer 3 detail → cascading-failure explosion model →
  CFD-grade aero → laser sail/ACMF papers (move to year 2; independent of the engine).

**Cut outright from year 1:** antimatter production (paper #4 — pure literature study, no engine
dependency), IMBHs, ring-gap resonance detail, the 20+ future ships list. None touch the
critical path.

---

**Phase 0 progress (2026-07-04): ALL SEVEN TIER-1 ITEMS RESOLVED.**
- 1.1 radiation determinism → deferred-apply protocol (radiation_shielding.md, determinism.md, ship.md §7)
- 1.2 ephemerides → analytic Kepler + first-order Laplace-Lagrange everywhere; SPICE C¹-blend
  for the solar system; N-body demoted to offline validation (generation.md §4.10, physics.md §2.1). User-selected.
- 1.3/1.4 spiral generation → material arms frozen at epoch in guiding-center space
  (time-invariant cell counts, drift bound documented); query unwinding φ_g = φ − Ω(R)T specced;
  φ_guide added to the epicyclic formula (generation.md §3.1/3.2/14.3, determinism.md). User-selected.
- 1.5 mean field → hybrid analytic (McMillan + rigid bar + Cox-Gomez + warp, potential-density
  pairs, self-consistent with star generation by construction) + low-order basis expansion for
  LMC wake/halo response only (physics.md §3.1-3.2, SCOPE.md, generation.md §11.5). User-selected.
- 1.6 hawking.md → rewritten from black_hole_paper (authoritative); ship.md Ship 3 updated.
- 1.7 force law → velocity-complete weak-field geodesic form, canonical momentum
  (physics.md §2.5, SCOPE.md; photon-deflection validation target added).
- Tier-2 errata batch: γ³ thrust, γ² ISM drag, NS radius, leak/hatch flow formulas,
  Chen-Kipping giant mass-first rule, stress-recovery pipeline, runtime reactor point
  kinetics, GR handoff mismatch corrected to ~3×10⁻⁵ + integrator-restart event, energy
  validation target reformulated, runtime-ODE claim amended.

**Tier-3 spec gaps: ALL RESOLVED (2026-07-04).**
- T3.1 flare schedules → hierarchical bucket seeding, O(1) at any T (generation.md §3.4)
- T3.2 query unwinding → resolved with 1.3 (generation.md §3.1, determinism.md φ_guide)
- T3.3 Gaia → 100 pc complete + landmark/OB/exoplanet anchors + recorded-sampling-fraction
  shells; suppression divides by sampling fraction (physics.md §3.7.1)
- T3.4 save format → per-component damage parameterization (ablation depth maps, dpa
  profiles, activation inventory); transport engine writes through it; octree re-rasterizes
  exactly on load (ship.md §8, radiation_transport.md §7). User-selected.
- T3.5 ground contact → pad kinematic constraint + leg spring-damper + friction cone +
  tip-over check (ship.md §5)
- T3.6 autopilot → quaternion-PD attitude hold + actuator allocation + burn executor +
  display-path predictor + open-loop scripted guidance profiles (ship.md §9, interface.md
  pointer). User-selected (option: + scripted profiles).
- T3.7 timewarp → per-subsystem max-dt table, analytic torque-free attitude propagation,
  swept-sphere system coverage, honest compute-bound warp ceilings (ship.md §10)
- T3.8 NRLMSISE → C FFI + deterministic synthetic solar-activity model (physics.md §2.3)
- T3.9 physics.md reconciliation → runtime = Rust + 2 C FFI libs only; scipy/numpy/astropy
  demoted to offline/validation; custom Rust DOP853; research-numbering clash fixed
  (physics.md §1.3, §2.2, §2.3, §4.2, §5)
- T3.10 transport tables → staged refinement plan + laptop-scale fallback; ablation coupled
  to 1D thermal relief with explicit fluence↔time mapping (radiation_transport.md §3, §7)

**Phase 0 is COMPLETE.** All Tier-1, Tier-2, and Tier-3 findings are resolved in the design
docs. Next: Phase 1 (determinism kernel + gravity spine) per the strategic plan below.

---

## INDEPENDENT CONSISTENCY RE-READ (2026-07-05)

Two fresh-context reviewer agents (one cross-document consistency, one quantitative
recompute of ~120 numbers) audited the post-repair docs. Core result: the physics core held
(force-law limits verified analytically; BH chain, star-count derivation, mission table,
precision claims all reproduce). Findings, ALL FIXED same day:

- **C1 (critical):** the deferred-apply fix itself had introduced "auto-reduce warp on
  wall-time" — hardware-dependent dt. Rule now: granted dt NEVER depends on wall-clock;
  compute shortfall only slows wall execution (radiation_shielding.md, ship.md §7/§10).
- **Thermal update equation in ship.md** had a factor-2 conduction error (double-counted A);
  corrected to thermal.md's (C − dt·A)T_new = C·T_old + dt·q.
- **Stopping-budget threshold:** punch-through hits 1% of P_H at M ≈ 2×10⁸ kg (3000 fm),
  not 10⁷ kg — the event detector as previously specced fired ~20× too late in mass.
- characterization.md BH block rebuilt to hawking.md Tables 1-4 (capture_fraction removed);
  S_e stress-recovery matrices added to the structural deliverables.
- SCF purge completed in physics.md (5 stale hits) + generation.md §13.5 diagram; "only
  runtime ODE" claim fixed in SCOPE/thermal.md; restored generation.md's Layer-3 header
  (deleted during the §4.10 insertion); DoseRateMap struct unified; per-element table dirs;
  ~20 further minors (per-tick budget canonicalized, 1-3 μs mean-field cost, active-source
  counts ~1,000-2,500 with the 1024 cap binding in the solar neighborhood, query-radius
  ladder recomputed + progressive-search mode, MS luminosity branch boundary 20→55 M_☉,
  Mestel T-exponent, IFMR intercept 0.48, O₂ metabolic rate kg/mol confusion, krad/kGy,
  solar-degradation units, decay heat at 1 yr, binary fractions aligned, HUD epoch).

### BLACK_HOLE_PAPER ERRATA — INDEPENDENTLY VERIFIED 2026-07-05 (second agent, against the
### authoritative final draft; the paper is submitted; Ben to disposition)

**Version sync:** /Users/ben/black_hole_paper (main) and the inquiry_project snapshot are
IDENTICAL (same git HEAD a391472, same uncommitted modifications; diff -rq clean).
**Authoritative final draft:** drafts/GREFF_2026_06_28.tex — built by build_latex.py from
draft_6_submission.md and content-identical to it in every checked section. No draft_7.

Verified verdicts (recomputed from first principles AND against the paper's own scripts):

1. **CONFIRMED — §6.2 stopping-budget mass** ("M < 10⁷ kg"): kT = 240 GeV at M = 4.41×10⁷ kg;
   punch-through = 1% P_H at M ≈ 1.85×10⁸ kg (3000 fm). No script computes it (derivation
   gap, not transcription). §6.2's unfed Δv ≈ 0.18c endpoint inherits the error.
2. **CONFIRMED — §4.1 coil stored energy "~50 GJ":** Biot-Savart integral gives 848 GJ in
   the cavity sphere alone; total ≈ 17 TJ (matches an (L−M)I² inductance estimate). No
   calculation exists anywhere; 50 GJ is asserted, and bh_stability_control.py hardcodes it
   (its own "L ~ 1 H" comment is incoherent with it: 2E/I² = 5.4×10⁻¹¹ H).
3. **CONFIRMED — §4.2 readout field "30 T / 30 pm":** actual 0.0296 T and 33.8 nm.
   TRANSCRIPTION error: bh_stability_control.py COMPUTES the correct values but a hardcoded
   summary line prints "~30 pm". Control still closes (33.8 nm ≪ 100 mm offset).
4. **CONFIRMED — §2.4 sub-extremality "1.53×10⁻¹⁰":** actual g/g_ext = 4.58×10⁻²
   (g_ext = M√(4πG/μ₀) = 2.583×10⁷ A·m). ROOT CAUSE FOUND: verify_design.py has an extra
   factor of c (`g_ext = M*c*sqrt(4πG/μ₀)`); the paper's number is exactly the correct one
   divided by c. Qualitative claim survives (still sub-extremal) but "RN correction ~10⁻²⁰"
   becomes ~2×10⁻³.
5. **WITHDRAWN — §3.4 muon flux/dose:** NOT an erratum. The first reviewer divided pre-loss
   channel power by mean EXIT energy; the paper's pair (10¹⁰ m⁻²s⁻¹, ~4×10⁴ Sv/yr) is
   self-consistent — dose is number-flux (MIP) driven. Recomputed: 1.66×10¹⁰ m⁻²s⁻¹,
   1.7×10⁴ Sv/yr. Paper stands.
6. **CONFIRMED — §5.1 monopole spacing "~2.5 μm":** (V/N)^⅓ = 0.180 μm for 7.18×10¹⁴ in a
   10 mm sphere. Errs conservative (the BPS range condition gets EASIER).
7. **CONFIRMED — Table 2 top-quark S = 1.2×10⁻⁵:** Fermi integral gives 7.29×10⁻⁵; the
   table's contribution column already uses the correct value (verify_design.py agrees);
   f = 3.503×10⁻³ unaffected. Transcription.
8. **CONFIRMED — §3.6 electrosphere scale height "~2000 fm":** the paper's own formula at
   μ_e = 9 MeV gives 278.6 fm; cfl_electrosphere.py computes exactly 278.6 fm. Transcription.

**NEW findings from the verification pass:**
- **§4.3 equatorial field 227 T is a linear-gradient extrapolation invalid at r = R_c**:
  exact Biot-Savart gives 155.5 T at (ρ=2 m, z=0). Propagates to the mirror ratio and
  f_cap: 1−√(1−155.5/1000) ≈ 8.1%, not 12.1%. Net thrust is nearly unaffected (f_cap only
  sets the internal Usov recycling equilibrium; β is logarithmically robust), but Eq. 10's
  L_pair rises ~1.5× and Eq. 19's Larmor radius changes accordingly. Pole value 390 T is
  exact.
- Escaping muons exit at ~11.5 GeV (γ ≈ 110, decay length ~72 km — still free-streaming),
  not γ ~ 2,300 / 1,500 km; escaping POWER is ~4 MW, not the 84 MW above-budget channel
  power (conservative direction; dose unchanged).
- verify_design.py §2 uses a provisional g = 1.303×10⁶ A·m vs. the paper's final
  1.182×10⁶ (explains small internal mismatches, e.g., B_horizon 5.9 vs 5.4×10³⁴ T).

Script bugs to fix in black_hole_paper/calculations/: verify_design.py g_ext factor-of-c;
bh_stability_control.py hardcoded "50 GJ"/"~30 pm"/"L ~ 1 H" lines; provisional-g cleanup.

hawking.md carries corrected values with erratum notes where they differ from the draft.
