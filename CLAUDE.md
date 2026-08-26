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

## State (2026-08-25 — session 21 = plan S9 complete)

Design complete: all 29 critical-path Layer-2 docs **Reviewed 2026-08-14** (the review register was
closed 68/68-discharged and deleted 2026-08-19; findings live in doc change logs + git history).
Eighteen sessions, every one gates-green and committed; three multi-agent code reviews
(sessions 6, 10, 12) plus per-session review waves (14, 15, 16, 17, 18) with every confirmed finding fixed. `scripts/check.sh` = the 5-gate battery
(fmt, clippy, cargo test, offline pytest, certificate regen + diff). **243 Rust + 50 Python
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
refusals, thread bit-identity; a two-agent review wave (SDC/audit + an independent continuum oracle) found **no correctness bugs**, closed a real test hole (the battery was blind to the compressible dilatation terms — mutation-proven), and hardened eight latent hazards. **Certificates byte-identical** (gas diffusion is opt-in config;
stations don't schedule it — S4's ◆C1 turns it on with real transport tables). Deferrals:
skin-friction gas debit (mount-reaction wave); species-enthalpy flux + TableEos-T refresh (S4
spine); COUP-8 row + config grammar (S4); serial assembly (GPU wave); θ-diffusion (S8). Finding:
impulsive walls/drives at stiffness ring the truncated sweeps past gas positivity — declared
ramp schedules (the COUP-7 discipline) are the cure; RL10 already carries its injector ramp.

**Session 16 = plan S4: REAL PROPERTIES + ◆C1.** `crucible-solvers::transport` = the **FND-7 §3.3
spine's transport slot**, filled for the chemical regime (FND-7 0.5 / OFFL-5 0.3 §3.1a / OFFL-3
0.5 / SOLV-1 0.4.1 / COUP-3 0.4.3 / META-3 0.8 §6.10 — amended first, then coded). Backbone =
**mixture-averaged Chapman-Enskog** (the chemical sibling of Saha/QEOS: a physics model defined
everywhere, not a data corner), keyed on the **runtime local state (p, h, Z)** — the S22 rule;
(T, p, Z) is the generation coordinate only, because a runtime `T(p,h,Z) → μ(T,p,Z)` chain would
compound two `interp_error_bound`s into a quantity neither describes. **One flux, decomposed:**
`Σ h_k j_k = −ρD[(c_p,eq − c_p,fr)∇T + (∂h/∂Z)|_{p,T}∇Z]` exactly on the equilibrium manifold, so
the ∇T limb folds into `k_eff = k_fr + ρD·Δc_p` (the classical equilibrium conductivity, > 3× k_fr
at the dissociated corner) and the ∇Z limb is the resolved species-enthalpy flux — never counted
twice, and Pr stays 0.2–1.5 at both ends. **One owner:** `WallLaw` is stateless (its private
c_p/μ/Pr *was* the degenerate occupant, now stated once as `transport_constant`); `gas_diffusion`
holds no transport either — a per-cell `GasTransportField` refreshed once per Picard iterate,
faces = the two-cell arithmetic average (exact in the constant limit ⇒ certificates byte-identical;
**mutation-proven** — reading one side leaves the whole S3 battery green and fails only the new
test). Offline: CEA supplies the caloric columns (same solver/mode/coordinate as the EOS surface),
**Cantera evaluates transport on that composition and never equilibrates anything**; six columns,
14 763 nodes spanning exactly the EOS envelope, rule-space bounds ≤ 3.3% against the declared
10–20% band; two generation refusals armed (unmapped species; a Cantera↔CEA two-fit c_p
cross-check that *bounds* the > 3500 K extrapolation). **Liner ρc_p is physical** (ruling #4): the
grid-thickened ring now carries **two** declared homogenizations — κ resistance-preserving, ρc_p
**capacitance-preserving** (`ρc_p_metal·t_real/t_model`). Wall-law band corners are a declared
`band_factor` on h, replacing the `cp` proxy. **◆C1 MET:** coarse RL10 marched 15 819 steps to a
settled readout, no halt, no schedule tuning, audit clean — F +0.30%, Isp +0.25%, c\* −0.07%,
p_c −0.02%, **jacket heat −19.8%**, liner T_max 410.1 K, vs the S3 spine. The shape of that answer
is the finding: the missing forces barely move plane-integrated performance at dial 5 (the wall law
already owned the wall) and move the wall term by two orders of magnitude more. **Method note
worth keeping:** the PRE-review march reported jacket heat −3.27%, which read as harmless; the
physics reviewer called that a coarse-tier cancellation hiding a 1.9× error in the wall law's
driving potential, and fixing it turned −3.27% into −19.8%. A small delta after a large model
change is a reason to investigate, not a clean bill. **KNOWN LIMIT (recorded):** with physical areal
capacitance the liner's clock (~37 ms) exceeds the 12-flow-through budget (~11 ms) — the gas field
settles, the wall is still warming; the cure is run length (ruling #4 accepted it), not tuning.
Cost ~7× per march; station-5 scores stay S2/S3-spine of record with the measured S4 delta
recorded in the certificate — re-scoring rides the COUP-5 ensemble wave.

**Session 17 = plan S5: THE COLD/UNBURNT BRANCH (split).** The burn-progress **c = 0 branch**
shipped: `crucible_offl::FrozenReactantEngine` = the **gas-phase ideal-gas frozen reactant mixture**
(gaseous H₂+O₂ at the proportions Z sets) via the **same CEA engine and enthalpy reference** as the
equilibrium surface — a correctness requirement, since the SOLV-4 §3.6 blend `h = (1−c)h_u + c·h_b` is
a category error on two references (verified: h_u−h_b = the physical heat of reaction, +12.4→13.5 MJ/kg
across MR 4–6). Artifact `tables/chem/lox_lh2_unburnt_v0.1.0.h5` (670 KB, 9×121×15, cryo-valid to
~100 K; bounds T 1.66 K / ρ 0.81% rule-space) uses the equilibrium surface's FND-5 schema, so the
**existing `TableEos` occupant binds it with no new code** (SOLV-1 §3.4) — the whole runtime cost of the
cold branch. The **frozen↔shifting bracket** is shipped as a declared model-form band
(`KineticEfficiencyBand`: raw bracket a few %, ~3.8% at ε=61; the `jannaf-eff` ~0.8–1% is the delivered
estimate *inside* it). Docs amended first: OFFL-3 0.6, FND-7 0.5.2, SOLV-1 0.4.3, SOLV-4 0.4.1,
META-3 0.8.1. **The split (Ben):** the frozen-mode `{ρX_k}` **field-advection** widening (a `NCOMP` change
to the flat conserved state; stable Rust blocks the clean const-generic via `NPRIM=NCOMP+2`; no consumer
before S18) → **S5b**; per-species diffusion + Soret/Dufour + Stefan-Maxwell ride it. The establishment-
refusal / overshoot-ceiling retirement is **enabled** here (the branch is representable + marchable) and
**realized at S6** (the c-blend routes transient cells to it); the pinned v0.3.2 burnt surface and ◆C1
config are untouched, so the shifting stations stay **byte-identical**. **Two flagged rulings (mine):**
settle budget grows as run length at S7 (ruling #4, not paid now); station-5 re-score left to the COUP-5
wave (S18). A two-agent review wave (code/physics + doc-claims) ran before commit.

**Session 18 = plan S6: IGNITION.** The chemical regime lights: `NCOMP` 6→7 (the fixed-`+1` burn-progress
slot `ρc`, inert by default — shifting stations untouched); `BurnBlendEos` = the **energy-conserving
flamelet blend** (both branches at the cell's own (p, h, Z) on S5's shared reference — heat release is
EOS-implicit; c = 1 recovers shifting **bit-for-bit**); `Combustion` = the SOLV-4.4 **bistable-Nagumo
pushed-front** rate law (matched (D_c, K) ⇒ front speed closure-set at S_T, width Θ·Δ — FSD |∇c| and
Fisher-KPP rejected as degenerate/pulled); closures = ONE offline Cantera `S_L`/`τ_ign` (p, T_u, Z)
surface. **The S6-close envelope set (OFFL-3 0.6.2 — Ben: expand the tables, never throttle the spark):**
the build's tables were ignition-incompatible (burnt ceiling *below* unburnt = an inversion the blend trips
mid-ignition; blast compression then drives mid-transition cells past the unburnt ceiling) — cured by
**burnt v0.4.0** (h to +1.225e7; rule: burnt ≫ unburnt), **unburnt v0.2.0** (t_ceil 2900 K), **ignition
v0.2.0** (T_u ~3000 K; **contract:** its T_u envelope covers T_u at the unburnt h-ceiling), **transport
v0.2.0** (re-derived; the armed transport↔EOS refusal caught it; hot-side c_p two-fit tolerance 4% =
measured 2.44% max ×1.6) — all **strict extensions** (old nodes bit-exact, verified) ⇒ certificates
byte-identical (station 1 gains only the ρc MMS column at order 2; station 5 provenance text).
**The igniter is a literal electrical spark** (`spark-igniter-class` PINNED, META-3 0.8.3: the one cited
datum is deposited energy — H₂ MIE 0.017 mJ → exciter 0.1–20 J; COUP-7 0.4.1: placement is config), a
**bounded pulse ending ~at ignition** (the S3 impulsive-drive discipline). Blend pure-limit
thresholds **asymmetric** (SOLV-4 0.4.3): pure burnt from 1−2·`BURN_COMPLETE` (one owner — c pins
~2e-7 *below* the reaction's fixed point; equal thresholds are a knife edge), pure unburnt at 1e-9
(no cold-side attractor; the cold crossing step scales ×v_b/v_u ≈ 7.8 — review-wave finding).
**Battery 5/5:** `flame_1d` THE gate 2.0% coarse-vs-fine; `spark_box` lights; `lean_no_light` refuses;
`ignition_delay` t/τ = 0.12; **finding:** a lit closed **adiabatic** box cannot flame out (quench is a
heat-loss phenomenon) — pinned as `adiabatic_box_cannot_flame_out`; the FLAMEOUT `quench_box` rides **S7**
(cold-wall conduction via the blend↔class-D coupling S7 builds anyway; the 0.8 mm fixture is at the
quench-distance scale — H₂/O₂'s own distance is a few× under H₂/air's ~0.6 mm, so lower p or a narrower
gap may be needed). Deferrals → S7: igniter config-grammar/`run.rs` wiring (◆C2), stiff class-R
implicit auto-ignition, near-vacuum blend tangency, `turbulent-flame-speed` pin, **the cold-side/low-p
closure-envelope guard** (the ignition surface floors at T_u ≈ 230 K / p ≈ 6.8 kPa; a cryo-fill cell
refuse-halts at S7 wiring without the cold analogue of the BURN_COMPLETE guard or wider floors);
N_θ > 1 → S8.

**Session 19 = plan S7: FIRST STARTUP VERIFICATION + ◆C2.** "DOESN'T START" is a computed outcome:
every march ends in a typed **COUP-4 verdict object** (WORKS/DOESN'T-WORK + mechanism/location/time/
diagnosis + the recorded criterion; `runs/<name>/verdict.txt`). WORKS criterion = config grammar
(commanded p_c/thrust arm it; `eps_works` 0.02 default, dwell, `flowthroughs` = T_S1_HORIZON;
FAILED_TO_REACH at expiry); NEVER_IGNITED/FLAMEOUT on the **ṁ-scaled reacting-measure floor**
(`max(EPS_IGNITED, 1e-4·delivered ṁ)`); melt + **primary-stress burst** margins via SOLV-6 v1
(`structural_margins.rs`, Roark subset, FS 1.1/1.4; the ASME primary/secondary split — strain-controlled
thermal stress is a reported diagnostic, never a halt; declared shell = the RL10 tube, |p_gas−p_coolant|
load, 77 K cold allowables anchor). **The class split of record:** auto-ignition = **class-R implicit**
(fixed-structure BE-with-τ-refreeze node solve in the SDC sweeps, realized-rate quadrature + audit row,
node-0 rate capped at the realizable (cap−ρc)/Δt; marched at dt/τ = 3.04 where the S6 explicit tier
halted — that tier retired; parks bit-exactly at 1−BURN_COMPLETE). **The real quench_box:** blend↔class-D
coupling closed (blend `temperature_w` + spine `interrogation_php` through the ONE config-selected
`ChemEos` seam) — a lit 0.8 mm-gap front between 420 K walls at 0.1 atm FLAMES OUT (R 2.06e-4 → exactly
0; adiabatic control stays lit); 420 K = the gas-phase model's honest coldest wall (the products surface's envelope floor reads ~407 K at the fixture state — S15 owns colder). Cold-side validity set: **ignition v0.3.0** (strict extension, verified:
cold T_u rows 60/90/120 K → envelope floor 75 K; p top-cell midpoint 4.44 MPa → ceiling 6.67 MPa — the
stiff mini-sim caught mid-transition compression past the old inset-artifact ceiling; isolated-zero
generation refusal armed) + **unburnt v0.3.0** (declared re-gridded narrow-Z variant, Z = the burnt band:
mid-Z validity to ~96 K for the 120 K wall-cooled fill; wide-Z v0.2.0 kept as the S5b base) + the
runtime **declared non-reactive floor** (below the ignition surface's own envelope floors, or below the
reactant branch's h-floor — the cold analogue of the BURN_COMPLETE guard) + the blend's **cold-side
partition extension** (reactant sub-state pins at its floor under heat loss, products absorb the balance;
`B_PARTITION_MIN` = 0.01). `turbulent-flame-speed` **PINNED** (Zimont/Peters; first consumer = the
declared startup wrinkling — 1.0, the laminar tier: measured, wrinkling > 1 in the S_L crossover band
drives S_T toward the sound speed and trips the explicit class-A front-carrier guard; the turbulent
value lands with the S_T-CFL hardening, S8/S16). Startup shake-out cures: `EPS_ROBIN_RESID` 1e-6 → 1e-4 + absolute 1e-2 W
floor (near-vacuum exchange contraction — halt-gate only, no accepted number moves); declared gas-phase
injection state −2.5e5 J/kg (~215 K; cryo-liquid = S15); the 10–20 mbar declared altitude cell bounds the
pre-ignition cold jet; the spark is a **short exciter burst** (~50 µs ≈ kernel residence — a longer
deposit pays the ignition enthalpy once per gas replacement in the 200–400 m/s fill stream; 10 J,
innermost ring, fired at 3 ms ≈ 0.9 bar chamber). ◆C2 = `configs/rl10_startup.toml`:
thirteen shake-out attempts, every halt a cured finding (the spark's short-burst form; the
pulse-cut-at-light discipline; the class-R quadrature invariant projection; the two-sided-plus-
fallback partition; the products-window mid-transition bracket + scan-first + bisection backstop in
the ONE root finder; the establishment grace; and the headline: **a lit kernel cannot anchor in the
100–400 m/s fill stream — blowoff — so flame-holding is the ASI TORCH's job**, exactly as the real
engine does it; the igniter object gained the cited kJ-class torch tier). The certified march:
spark 3.04 ms → chamber-scale light-off ~24 ms (R ×4 decades) → stable torch-anchored flame
(R ≈ 1.6e-2, 10× the FLAMEOUT floor) → p_c plateaus ~0.90 MPa / F ~22 kN (mostly-cold chamber:
a laminar front cannot spread across the swept face; the spreaders — turbulent S_T, distributed
elements, 3-D recirculation — are S8/S16 physics) → the full 84 ms horizon marches clean (55,658
steps, audit green, R steady at 1.7e-2 for 60 ms) and ends in the typed verdict **DOESN'T WORK
(FAILED_TO_REACH; physical): p_c −71.2% / F −70.8% of commanded** — "doesn't start" is a computed
outcome whose mechanism maps onto the plan's own next phases (S8 S_T-CFL guard + turbulent
consumer, S16 distributed injection). Dial-16: cure demonstrated
in kind at dial 5; the 33×-cost dial-16 blend startup is a desktop command. Deferrals → S8+: N_θ > 1 combustion; establishment-grace
constant if a young kernel needs it; S5b species work; dataviz gate unchanged.

**Session 20 = plan S8: AZIMUTHAL CAPABILITY (Phase 3 opens).** The flow operator is genuinely 3-D:
**r = 0 axis at N_θ > 1** (FND-2 §3.2 θ↔θ+π parity-pair gather, u_r AND u_θ negated — arithmetically
identical to the mirror at N_θ = 1, certified corner unmoved) and **mixed per-brick N_θ** (FND-2 §3.4
made concrete: 2:1 ladder adjacency + ≥ NGHOST proper nesting validated loudly; maximal uniform-N_θ
pencil segments, TWO-pass sweeps — the fine side owns each jump face's flux, the coarse side applies
the exact aggregate, accumulation stays the single well-balanced difference — the push-style first cut
broke the uniform bitwise fixed point at the interface and the gates caught it at step 0). Gate battery
(solv1_azimuthal, 11 gates): axisymmetric N_θ = 8 ≡ N_θ = 1 **bitwise per plane** through axis AND
active refluxing (both fix-up orientations); mixed-world uniform bitwise fixed point; 3-D pulse crosses
the axis (conservation 1e-12, mirror symmetry 1e-10, far-side signal); transverse-flow u_θ-sign
discriminator; thread-count bit identity; Δt-vs-N_θ-profile gate; controller collapse/expand mid-march.
**Audit finding:** mirror-symmetric swirl collapses the momentum_theta tolerance (net stored cancels;
θ-sweep increments unledgered) → S[q]'s stored term = gross Σκv|q| both endpoints (COUP-2 0.3.1,
halt-gate only). **N_θ > 1 combustion** (◆C3 prerequisite): per-θ-plane re-key + θ-direction D_c faces;
a point-in-θ spark spreads to adjacent sectors first (gated). **The S_T-CFL** (SOLV-4 0.4.6/0.4.7):
σ_front joins the class-A Δt reduction (wrinkling > 1 marches honestly); the scale-separation guard is
the MODEL-FORM limit S_T > (2/3)·c (the mesh-rate form was measured resolution-dependent and cured in
review). **Review wave** (axis/reflux math + S7-carry fresh-eyes + doc-claims): every CONFIRMED finding
fixed in-session — headline: the S7 blend partition's balance↔trace fallback was discontinuous inside
its own scan bracket and B_PARTITION_MIN was not trace-bounded → superseded by the continuous window
form `h_b = clamp(balance, products window)` with the volume-weighted rule-space acceptance armed on
EVERY accepted root (the blend's true cold edge enforced there); class-R base projection symmetric,
advected content never shaved. **First turbulent consumer:** the ◆C2 rerun (`configs/rl10_startup_s8.toml`,
declared wrinkling 2.0 = the pin's band factor; raw Zimont estimate ~6 exceeds the single-constant tier)
— result in SESSION_LOG session 20. **SPLIT (improvisation rule, plan §8 v1.8):** the F_visc
θ-extension — design fixed (COUP-3 0.4.5), build → **S9** (a partial θ-stress tensor is wrong physics;
refusal stands loud); mixed-N_θ class-D + per-θ wall patches + mixed-N_θ combustion → **S11**; S9/S10
must land the uniform θ-tensor before ◆C3. Deferrals also recorded: controller 2:1-enforcement pass
(S11); class-R dt-Richardson gate, mid-transition root uniqueness guard, N_TAU_REFREEZE mid-band check
(S9 review wave).

**Session 21 = plan S9: FULL GEOMETRY (+ the S8 carries).** **The FND-3 kernel exists**
(`crucible_grid::geom3d`): the CSG analytic SDF tree (sphere/box/cylinder/cone/torus + the
revolved-profile leaf; 1-Lipschitz ⇒ sound pure-cell bound tests), STL import via the **exact**
van Oosterom–Strackee winding number (Barnes-Hut = a recorded perf deferral, FND-3 0.4),
jittered-stratified voxelization in the cylindrical measure (SplitMix64 counter keys excluding
the θ-sector ⇒ sectors are exact rotations; measured C_jitter 0.907 → declared 1.4; rate
exponent −0.723 vs the derived −2/3), and PLIC (normal from the authored geometry — the
winding number of a watertight mesh is piecewise constant, a build finding; volume-exact S-Z
offset by fixed-count bisection; area = dV/dd exact). **The S8 carry repaid (◆C3 unblocked):
F_visc carries the full θ-stress tensor at uniform N_θ** (COUP-3 0.4.6) — θ-θ implicit cores on
periodic ring stencils, curvature/cross limbs on the S3 Picard lag, θ work fluxes in the
total-energy bookkeeping; the CG dot nests per-(brick,θ-plane) partials so axisymmetric
N_θ = 2^k reduces bit-exactly like N_θ = 1 (gated bitwise per plane); θ-MMS orders 2.03–2.61 on
all seven components; the transverse-flow trap gate measures the m = 1 damping at 0.9% → 0.09%
(N_θ 8 → 32) of the partial-tensor scale. **Cut geometry is legal at uniform N_θ > 1**
(FND-2 0.5.3): six-aperture per-θ-sector `BrickGeom`, `build_with_geometry_theta` with per-sector
validation, the wall-closure θ-limb, aperture-weighted θ-faces in the sweeps, per-sector SRD with
θ-neighbor candidates — mutation-proven load-bearing; `stl_toy_chamber` marches END TO END on the
real sampled path (voxelize → ingest-as-is → 55 audited steps, drift < 1e-12). **Scope rule:**
every gas ring keeps κ > 0 in every sector (per-sector activity masks → S10/S11). **The carries
landed with a FINDING** (SOLV-4 0.4.8): the blend scan counts crossings (>1 = typed multi-root
refusal); class-R dt-Richardson order 1.80 at dt/τ ≈ 0.4–3; the N_TAU_REFREEZE witness
**falsified the 0.4.4 weak-τ-dependence premise** (reaction-driven compression: τ ×1/93 within
one node solve, lag 3.0e-1 at w/τ = 0.3) — cured by RESTATING the contract to the S3
truncated-Picard split (fixed count = structure, accuracy owned by the composed-order gate, the
witness armed as a measured-envelope pin), never by widening a tolerance. Still refused, typed
(→ S11): mixed-N_θ cut worlds, mixed-N_θ class-D, cut-θ class-D, per-θ wall patches, adaptive
N_θ on cut worlds, combustion on cut θ > 1 worlds. CSG/STL config grammar + engine 3-D assembly
ride S10/S11 (the kernel's first config consumer). Certificates byte-identical (gate 5).

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
  **NEXT = plan S10**: refinement + the AMR gate — static (r,z,θ) refinement tiles with
  conservative level interfaces (declared zones: walls, injector face, throat); the 2-D
  smeared-vs-sharp front study → the measured **go/no-go on dynamic front-tracking refinement**
  (if go: insert S10b/S10c; requires an FND-2 static-topology amendment). Candidates riding S10:
  the CSG/STL config grammar + engine assembly ingest (the kernel's first config consumer —
  or S11 with ◆C3), per-sector activity masks, the FND-3 §3.4 coarsest-reproducing-N_θ floor
  search, the voxelizer↔grid face-order unification. Read FND-2 §3.5/§3.6 + the S10 §5 entry
  before coding.
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
