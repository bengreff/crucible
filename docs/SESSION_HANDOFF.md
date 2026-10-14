# Session handoff

20 October 2026 · branch `main` (one line). The spray worker (step 5) on `main`; inquiry2 works in `.claude/worktrees/inquiry2` (RL10 contour, AMR, table generalisation, and transport cost in `core/medium.*`, `core/transport.cpp`, `core/flow.*`, `contour.*`: not touched from here). The 6 October state is in this file's git history (`git show 15d6777:docs/SESSION_HANDOFF.md`).

## Stopped at (20 October, about 19:40 CDT; session after /clear, READY TO CLEAR)

The Director's brief (after /clear, to about 20:00) had two parts.
- (1) Entropy as Table B's progress variable is approved. State its checks before code: 2(a) at every density, the resolution near equilibrium (nodes clustered toward it, the T and Y_H2O interpolation error measured there against 3 K and the composition gate), and the cost per cell. Code it next session unless the statement takes under 15 min.
- (2) The 12,176 K igniter cell is unphysical. Size the igniter to a physical one, add a guard above 6,000 K on h2o2_highT (a constant-cp continuation, labelled outside the model), rerun C1 on h2o2_highT, and report the fallback rate in the chamber and the convergent section.

Worker rules 13 to 16 were named in the brief but are not in the repository (AGENTS.md, TECHNICAL_PLAN, this file and memory searched). The brief was followed as written.

1. **Entropy's checks, stated (TABLE_B.md, *The specific entropy as progress variable: its checks*), not coded.** E0 to E5 are stated.
   - **A design consequence found while stating them (derived):** equilibrium is the entropy maximum at fixed (u, ρ, elements), so the **upper branch also rises in s** toward s_eq. The builder's upper-branch node placement and the runtime's 0-to-2 coordinate assume it falls, as Y_H2O does.
   - The branch rule is stated (engineering): a cell takes the nearer branch's curve by max |ΔY_k|, and ties go to CVODES, counted as E0.
   - The node rule: s_eq − s = span (1 − ξ)², uniform in the distance, because s is quadratic in |Y − Y_eq| near equilibrium. λ tends to 2k there, so it stays finite.
   - **Not coded, though the statement took under 15 min, so the brief's condition for coding it now was met.** The consequence above means reworking the builder's upper branch and its dcdt (ds/dt = −Σ μ_k/W_k dY_k/dt / T at fixed u and ρ). I judged that it would not fit beside the C1 runs and the handoff, and committed the statement alone (ca403ff), before any entropy code.
2. **The guard (`Species::tHigh` in `core/medium.*`, set in `adapters/reaction.cpp` when the fits reach 6,000 K).** Above it: constant cp, h linear, labelled outside the model. h2o2.yaml (fits to 3,500 K) is unchanged, so no earlier result moves. `ReactingFlow` checks tHigh between the media.
   - **inquiry2: `core/medium.*` is in your lane.** This is an additive change: a field, two helpers, and the four evaluators calling them. It will conflict textually with your pre-divided enthalpy coefficients on rebase. Keep the continuation (h linear, cp constant above tHigh) in your form.
   - Fast tests: core, thermo_cea, reaction and mixture pass.
   - The new `crucible_mixture_tests` 1b passes: cv continuous at 6,000 K to 4e-15, e linear above it to 2e-16, the 12,176 K cell inverts to 5e-12 K, and gri30 gets no continuation.
3. **The igniter, sized (TABLE_B.md, *The igniter and a guard above 6,000 K*).**
   - A spark (1 to 120 mJ per spark, sourced) is sub-grid, so the model's igniter is a torch: gas at about 3,100 K (Repas 1994, provisional). Declared power 50 to 160 kW for C1 (*inferred*: a torch flow of 1 to 3% at the stoichiometric heating value; no sourced torch flow was found).
   - The run uses 50 kW, 50 J over 1 ms from 0.2 ms. The old 300 J igniter was 300 kW.
   - **Derived, and the open part:** a volumetric deposit into the 0.007 kg/m³ early fill drives slow cells above a torch's temperature at any power. The physical representation is a hot-gas inflow, or a deposit that stops at the torch temperature. That is in `core/flow.*`, inquiry2's lane.
4. **C1 on h2o2_highT, rerun (both runs pass the 0.206 ms failure point).**
   - **The torch run (50 kW), to 4 ms:** T_max 4,029 K at 0.288 ms (3,489 K after 1 ms). The guard is never engaged.
   - **The fallback rate on h2o2_highT** (`table_b/c1_fallback_torch_highT_2026-10-20.txt`):
     - **Chamber and convergent section (columns 0 to 15): 4.146% of their updates over the whole run, 0.000% over the last ms.**
     - The 4.1% is all in the first 0.5 ms, at ρ < 0.28 kg/m³ before the chamber pressurises: 3.45% failed corners (the low-density 2(a) turn-back), 0.65% purged ambient N2 above the Z_N axis, and 0.04% η above the axis.
     - Divergent: 86.0% and 82.3%. All cells: 45.1% and 41.2%.
   - **The 300 J igniter, as the guard's check:** T_max 7,840 K at 0.216 ms, with the guard on from 0.212 to 0.226 ms.
     - Stopped at 0.828 ms with the igniter still on: my watcher read a half-written CSV line. It shows the guard carrying the old igniter through its peak, not to 4 ms.
     - Against the h2o2.yaml 300 J run, the torch removes the igniter's η overflow from the chamber (4.85% to 0.04%).
5. **Also measured: the split on the h2o2.yaml samples (e6d4378's run).** Chamber and convergent section: 9.24% of their updates over the whole run (4.85% the igniter's cells above the η axis, 3.76% failed corners, 0.62% Z_N), and 0.000% over the last ms. Divergent: 86.7% and 82.4%. `tools/fallback_rate.py` now takes the throat station as an optional 4th argument.

**Next, in order:**
1. Code entropy as builder candidate 3, with the upper branch rising in s (E2's node rule, midpoint outputs). Run E1 to E3 on the small, O/F 4 to 6 and C1-range tables, then E4, then E5 on both sample sets (`/tmp/c1fb/states_highT.csv`; `/tmp/c1g/states_torch50.csv`, already in h2o2_highT's energies). The /tmp files are not in git: if they are gone, rerun the commands in TABLE_B.md.
2. The torch as a temperature-capped deposit or a hot-gas inflow in `core/flow.*`: inquiry2's lane, or the Director's call on who does it.
3. The runtime's branch rule (E0) once s is in the builder.

## Stopped at (20 October, about 19:15 CDT; session after /clear, READY TO CLEAR)

The Director's brief (after /clear, to about 20:00): (1) an exact Jacobian for the builder's integration, check pre-registered, then the reaction tests and C1; (2) the small Table B in O/F 4 to 6 and the RL10 pressure band with its gate table; (3) the fallback rate in a C1-like run, measured, then domain against progress variable by the Director's rule (a declared domain if under about 1% of cell updates). All three done; TABLE_B.md has every result under its stated check.

1. **The analytic Jacobian (cdce5fc; TABLE_B.md, *The analytic Jacobian, stated before its code* and *Results*).** `ReactionStep` (the builder's integrator and the runtime's CVODES fallback) now uses Cantera's net-production derivatives plus T and the constant-volume energy (`Jacobian::Analytic` default; `DifferenceQuotient` kept as an option).
   - **Check 1 (FD match) passes after a flagged restatement.**
     - Found a Cantera 3.2 bug: `getNetProductionRates_ddT` is stale after a composition change at the same T and ρ. Fixed with `invalidateCache()`.
     - Cantera's one-sided `rtol-delta` is set to 1e-7, not 1e-8 (round-off scan).
     - The reference is Richardson central differences at 1e-5.
     - Worst column 3.1e-7 (h2o2) and 4.4e-7 (highT), bar 1e-6.
   - **Check 2 (the 23 cold dense-radical starts) failed with the Jacobian alone.** The mechanism was written down, then measured: CVODES's first step is 3e-7 to 2e-6 s against recombination times of 5e-14 to 2e-10 s, and its default 10 corrector failures stop short.
     - Fixed with `CVodeSetMaxConvFails(40)`: 0 integrator failures, and the 6,000 K branch integrates.
   - **Check 3 passes.** PaSR output is identical at 10 and 40 failures. Against the pre-Jacobian build, the PaSR integration errors move within their limits.
   - **Check 4 passes.** C1 4 ms: c*, Isp and thrust identical to printed digits, injector pressure 0 difference. The reaction wall time is 19% slower, on a loaded machine.
2. **O/F 4 to 6 and a declared 1 to 10 MPa band (b52679d): Y_H2O passes 2(a) and 2(b) (0.010 K).** 2(c)'s 2 failures are the builder's refusals at η 0.6, ρ 80 kg/m³ (about 200 MPa, outside the band). No RL10 measurement was read.
3. **The fallback rate in C1 (e6d4378; TABLE_B.md, *The fallback rate: measured*): 47.97% of cell updates over the run, 41.21% over the last ms. Over 1%, so by the Director's rule the progress variable is reopened.**
   - The steady part is entirely the divergent nozzle (columns 19 to 31, ρ < 0.28 kg/m³). Every failing corner there is Y_H2O's 2(a) turn-back, on the C1-range table's 83 low-density nodes (ρ ≤ 0.069). The chamber and convergent section are on the table at every update.
   - **Proposed for the Director (a design choice, derived, not tested): the specific entropy as progress variable.** At fixed (u, ρ), dS/dt = Σ (r_f − r_r) R ln(r_f/r_r) / ρ ≥ 0, since all 29 reactions are reversible with reverse rates from Kc (checked). So 2(a) passes by construction at every density. Still to measure: the N gate near equilibrium, the cost per cell, and the span on weakly reacting branches.
   - Tooling: the builder's `CRUCIBLE_TABLE_B_FH_NODES`, the chamber study's `CRUCIBLE_STATE_SAMPLES`, `tools/fallback_rate.py`, and `crucible_state_samples_eta` (re-expresses η in another mechanism's energies).
   - The C1-range table command is in TABLE_B.md. The samples are in `/tmp/c1fb/states*.csv` (not in git, 32 MB).

**Open findings:**
- **(Resolved at 19:40, see above: the guard and the torch-sized igniter.) C1 on h2o2_highT fails at igniter onset (0.206 ms), so the fallback run used h2o2.yaml (a flagged deviation).** Measured mechanism: the igniter drives a cell to 12,176 K (e 79.5 MJ/kg, mostly O and H). h2o2_highT's 1,000 to 6,000 K fits, extrapolated, turn cv negative above 10,630 K (OH), 11,340 K (O) and 12,140 K (H2O), so `Medium::temperature` cannot invert. GRI's fits extrapolate monotonically.
  - Fix options (engineering): thermo above 6,000 K (TM-4513's 6,000 to 20,000 K range, or constant cp past the fit), or a gentler igniter deposit.
  - `core/medium.*` is inquiry2's; coordinate there.
- **h2o2.yaml's η scale:** its u(6,000 K) is extrapolated 8% high, so η from h2o2.yaml states reads up to 0.079 low. Any table lookup from an h2o2.yaml run must convert.

**Next, in order:**
1. The Director's call on the progress variable (entropy proposed).
2. If entropy: state its checks in TABLE_B.md first (2(a) and the N gate on the small, O/F 4 to 6 and C1-range tables; the cost of s per cell). Add it as a builder candidate, then rerun the fallback measurement on the same samples (`python3 tools/fallback_rate.py /tmp/c1fb/states_highT.csv <nodes csv> <candidate>`).
3. The h2o2_highT thermo above 6,000 K, so that C1 runs on the table's own mechanism.

## Stopped at (20 October, about 18:30 CDT; session after /clear, READY TO CLEAR at about 100k)

The Director's brief (after /clear, session to about 20:00): (1) Table B: fix the 40 integrator failures first, then a small table in O/F 4 to 6 and the RL10's pressure band, then domain against progress variable (a declared domain with a CVODES fallback is acceptable if its rate in a C1-like run is measured and under about 1% of cell updates); (2) measure how much h2o2.yaml's 1-atm standard state moved earlier results, C1 at least; (3) declare three RL10 tier-1 modelling choices with expected effects.

1. **Table B's integrator failures (3f97e80; TABLE_B.md, *The integrator failures: two causes, measured*).** All 40 are lower branches from the frozen start, failing before the first step.
   - **(A) 17 at η = 0, fixed:** the T* bisection had no sign change, so it drifted to 6,000 K (15) or a round-off root (2). The bracket is now checked first. At or below its low end the start is the complete-combustion products. 0 failures at η = 0.
   - **(B) 23 at η 0.2 to 0.8, ρ 2.28 and 30 kg/m³, f_H 0.056 and f_st: the cold dense-radical start (the Director's suspect), open.** Flag −4 before a step, dT/dt 9e11 to 4e15 K/s. Probe (reverted): a forced 1e-18 s first step lets 22 of 23 integrate, but it breaks a 6,000 K upper branch that passed before. Inferred cause: the difference-quotient Jacobian at atol 1e-30.
   - **Next, in order:** an analytic Jacobian in `ReactionStep::follow()` (Cantera 3.2's analytic net-production derivatives, plus T and the constant-volume energy), or a scaled perturbation. Check with `ctest -R reaction` and C1, since the runtime's fallback shares the Jacobian. Rerun the small table. Then the small table in the C1 and RL10 range (f_H nodes at O/F 4 to 6; ρ over the RL10's chamber pressure band), then domain against progress variable by the Director's rule.
   - The command: `cd build && CRUCIBLE_TABLE_B_NODE_REPORT=<csv> ./crucible_build_manifold_table <out> 2 2 0 0,0.3 6 5 1e-3 30 2 ../tools/mechanisms/h2o2_highT.yaml` (11 s per candidate).
2. **h2o2.yaml's 1-atm standard state: measured and recorded** (TABLE_A.md, *What h2o2.yaml's 1-atm reading moved*). Its reacting species are 1 bar data read at 1 atm (Ar is 1 atm data; Ar and N2 are inert in C1).
   - Declared at 1 bar, C1 32x6 `eqt` moves c* +1.05e-4, thrust +4.2e-5, Isp +4.7e-5 and injector pressure +1.09e-4; ideal Tc 3,351.7 to 3,353.1 K. So every earlier C1 result moved by about 1e-4 or less.
   - The correction moves h2o2.yaml's CEA check worst from +6.13 to +7.60 K: the 1-atm reading partly cancelled OH's heat-of-formation offset. h2o2_highT corrects both, so nothing is re-run.
3. **RL10 tier-1 modelling choices (3f97e80; `docs/validation/rl10a33a/tier1_skeleton.md`, *Tier 1's modelling choices*).** Each choice has its expected direction stated before any run, all inferred:
   - **Rings as continuous annular inlets** that keep the elements' flow areas (LOX slot 0.85 mm wide at R 50 mm). More interface per unit flow area on the outer rings (R / 23.9 mm against round jets), so c* biased up and the near-face heat flux up.
   - **Swirl neglected:** larger drops and a longer flame, so c* biased low. 2-D with swirl is a tier-2 candidate if S7 shows the spray is sensitive.
   - **Transpired H2** as a uniform inflow over the porous face at about 10%, so the elements run at O/F 5.62. A lower face heat flux; c* slightly low.
   - This replaces protocol §4's uniform face for tier 1. The protocol (inquiry2's line) needs the matching amendment.
   - Still open: the per-ring counts, if the fuel-only outer row and the inner row without swirlers are not 27-element rings.

## Stopped at (20 October, about 17:50 CDT; extra session from 17:10, READY TO CLEAR at about 100k)

Ben's direction: straight line to RL10; speed must not block fidelity ("slow runs are fine for now"). inquiry2 owns `core/flow.*`, `core/amr.*` and the hot-pocket fix (lever 3): not touched from here. The Director's brief for this session, all four items committed:

1. **OH and the equilibrium check against CEA (b83e522; TABLE_B.md, *OH's heat of formation, the reference pressure, and the equilibrium check*).**
   - h2o2_highT now carries OH's heat of formation from Ruscic 2002, as CEA does (the Director's decision), and every species at TM-4513's 1 bar standard state.
   - Found: Cantera's `nasa_gas.yaml` declares no reference pressure, so the fits were read at 1 atm, putting each Kp off by 1.01325^Δn.
   - The check passes: worst |ΔT| 1.79 K over 153 cases, bar 3 K. With TM-4513's OH it would be +7.42 K.
   - **Open item:** GRI `h2o2.yaml`, the mechanism of every run before today, also declares no reference pressure. It is unchecked, and its CEA check gives +6.13 K beside.
2. **Table A to 6,000 K on h2o2_highT (f8c2c07; TABLE_A.md, *To 6,000 K*): passes.**
   - Criterion 2(a): 99th percentile 1.086 K at 200 to 4,000 K, and 1.515 K on the extension to 4,000 to 6,000 K.
   - C1 32x6, the table against direct Cantera on the same mechanism: within 0.007%; histories 1.2e-4; audit 0.825 K.
   - The mechanism change alone moves C1 by 0.03%, toward CEA (Tc 3,348.1 K against CEA's 3,347.7 K; GRI 3,351.7 K).
   - The table file is `build/table_a_highT.bin`, not in git. Its provenance holds the mechanism as the path `../tools/mechanisms/h2o2_highT.yaml`, so run from `build/` with `CRUCIBLE_MECHANISM=../tools/mechanisms/h2o2_highT.yaml`: the driver refuses a mismatched string.
3. **Table B small rebuild on h2o2_highT (df882ce; TABLE_B.md, *The small table rebuilt on h2o2_highT: the gate table*): 2(a), 2(b) and 2(c) still fail for every candidate; criterion 1 now passes.**
   - Nodes outside the thermo range fell from 146 to 0, so the review's open question about whether outside nodes enter the gates is moot for now.
   - Per-node report: `docs/evidence/table_b/small_highT_nodes_2026-10-20.csv`, from the builder's new `CRUCIBLE_TABLE_B_NODE_REPORT`.
   - **(a) 2(c): 40 CVODES failures at t = 0 (flags −3 and −4), all at ρ ≥ 0.17 kg/m³, mostly at η ≤ 0.6, at f_H 0.056 and f_st.** Not atol: 38 at 1e-20 (measured, then reverted).
     - Inferred: the lower branches from the cold frozen start, a dense radical pool at about 200 K.
     - **Next:** tag the branch in `NodeOutcome::failure`, and print the start T, Y and the first right-hand side of one failing node. That is an engineering fix; restate it in TABLE_B.md if the start or the integration changes.
   - **(b) 2(a) and the N gate fail at f_H 0.556 (O/F 0.8) and hot f_st (η ≥ 0.4).** Physics of linear-in-Y variables (inferred). **A design question is raised for the Director in TABLE_B.md:** state the table's domain and send everything outside it to CVODES, or choose a different progress variable.
   - **(c) No node of the small table is in the C1 range** (f_H 0.143 to 0.2). After (a), build a small table with f_H nodes there (for example, lean and rich counts that put nodes at O/F 4 to 6) to measure the gates where C1 and RL10 run.
4. **RL10A-3-3A tier-1 skeleton (f8b3196; `docs/validation/rl10a33a/tier1_skeleton.md`, a new name, no conflict with inquiry2's protocol files).**
   - It holds tier 1's settings to fix, and declared ranges for the missing element dimensions, all inferred from analogs (`analog_coax_elements.md`, Sonnet agent, screened before reading).
   - The ranges: LOX post ID 2.62 to 4.78 mm (central 3.54); tip wall 0.20 to 0.37 d_L; annulus gap 0.09 to 0.80 d_L; recess 0 to 1 d_L; pitch from face area 15.7 mm.
   - The agent's two incidental RL10 exposures are declared there; it recorded nothing.
   - **Open:** the ring layout is under-determined (27 per ring will not fit one pitch from the axis; `injector_element_inputs.md` also says the outer row is fuel-only and the inner row has no swirler), swirl (no swirl velocity in 2-D), and the transpired H2.
   - The protocol (`RL10A33A_PROTOCOL.md`) is still only in inquiry2's line.

## Stopped at (14 October, about 02:25 CDT; parked by the Director)

- **Cost note (b89f3fa, e18c710, `docs/evidence/COST.md`): decided by the Director at about 01:25 CDT, recorded in its *Decision* section.** The S7 step is the acoustic CFL across the 0.4 mm radial cells (measured by grid scaling), 17,500 to 22,000 steps per ms, 7 to 30 times below the convective limit (derived). Profile of blob: spray 29%, transport 40% of busy. Order: lever 3 (low's over-temperature, a correctness item) first, then lever 2 (spray subcycling, N chosen by measurement on S1 to S5); lever 1 (adaptive resolution with level subcycling) goes to inquiry2's B1; lever 4 (physical-time IMEX with implicit acoustics, about 3.5 to 23 times from guessed costs) is a note for later, no code. Profiles and `parse_sample.py` in `docs/evidence/cost/`.
- **S7 low runs above the flame temperature (lever 3; cause open).** T_max 3,634 K at 2.4 ms against 3,392 K for H2/O2 at 1 MPa (Cantera, measured; 3,367 K from 150 K reactants); pressure 1.00 MPa, budgets closed to 1e-14. Low's T_max is already 3,373 to 3,424 K from 0.8 to 1.2 ms and first passes 3,440 K between 2.3 and 2.4 ms (`/tmp/masc/low.log`); blob peaks at 3,396 to 3,414 K. Cells above Table A's 3,500 K call Cantera at about 75 ms each, so low's pace swings from about 38 to 220 ms per step.
  - **Spray coupling excluded** by reading the code (no clipping; exact deposits; liquid energy anchored at 298.15 K) and by the trace (no parcels and no deposit in the hot cells).
  - **The trace (e18c710; findings in `COST.md`, *Lever 3*, evidence in `docs/evidence/cost/low_hot_*`):** the rerun reproduces low step for step and crossed 3,450 K at 2.311 ms. Every cell above 3,450 K is in the first axial column against the injector face, at r 0 to 1.2 mm: inside the LOX post's radius and inside the 13.7 mm liquid core, where the RR mode leaves gas (it releases drops from the core's surface). The gas there is hydrogen-rich (f_H 0.14 to 0.29; H2 comes along the face across the post's lip) and 3.7 to 5.6 MJ/kg above the propellants' mixing line, against 0.5 to 0.8 MJ/kg in the flame's hottest cells elsewhere. So that zone gains enthalpy, by differential diffusion or a numerical source at the wall and axis corner. Blob injects through the post's area and stays near 3,400 K.
  - The trace run (`/tmp/masc/hot_low.log`, binary `/tmp/masc/crucible_mascotte_a10_hot`) ends by itself at 2.6 ms; its full CSVs are `/tmp/masc/hot_low_hot_{trace,cells}.csv`.
- **Wall layer ported to C++ (c7472b9, `crucible_wall_layer`)**, by a Sonnet agent. It reproduces the recorded Python a priori output on all 14 values to every printed digit (checked here), at 0.2 to 0.3 s per inversion. Eigen is optional in CMake: backhouse has none, so criterion 0(e)'s table is built on the Mac or with Eigen's headers copied into the tree there.
- **Table B: the builder and the runtime stated before code** (240bec1, e16ae9a, in `TABLE_B.md`). A Sonnet agent (aa381dda8aa9ee1e4) wrote the builder; it finished at about 01:45 CDT. Its files are uncommitted and not yet reviewed here: `core/manifold_table.{hpp,cpp}` (layout, save and load), `tools/build_manifold_table.cpp`, `ReactionStep::follow()` added in `adapters/reaction.{hpp,cpp}` (no unit test of its own), and the CMakeLists lines adding `core/manifold_table.cpp` to `crucible_core` and the `crucible_build_manifold_table` target. It builds clean, and `ctest -R reaction` passes.
  - **First small table (5 f_H x 2 Z_N x 6 eta x 5 ln rho, 300 nodes, 31 s per candidate on 2 threads), as reported by the agent:** every gate fails for all three progress-variable candidates, so no file was written, as specified. Worst turn-back 2.05e5 (candidates 0 and 1) and 2.67e6 (candidate 2); 26 to 27 N-gate failures; 2(b) worst dT 503 K; 2(c) 39 integrator failures and 2 of 11 frozen branches; criterion 1's sum check 3.64e-11 against 1e-12 (elements 2.8e-16 after projection). 146 of 300 nodes are flagged outside the thermo range, 120 degenerate, 16 zero-span.
  - **The agent's diagnosis (to verify):** the failures cluster at Z_N about 0 and eta 0 or 1, where the energy range puts equilibrium above the mechanism's 3,500 K, so the branch-start bisection runs on extrapolated thermo, CVODES fails from there, and near-zero spans blow up the turn-back ratio. Nothing was tuned. It is the same 3,500 K ceiling as low's Table A fallbacks (`COST.md`). **Next:** review the code against `TABLE_B.md`, then check whether the flagged-outside nodes are wrongly entering the gates, before any change to the eta range or the thermo ceiling (a stated decision, not a tune).
  - The agent's estimate for the full table (166,320 nodes): about 4.8 h per candidate on 2 threads.
- **The pipe test rerun is running** (Mac, slot.py, `/tmp/pipe/pipe.log`): nr 16 and nr 32 done, nr 64 running (0.5 ms at 02:13). Every check passes except 4b, the mass budget: 4.6e-12 at nr 16 and 1.9e-11 at nr 32 against a 1e-12 limit (measured; round-off size, not yet compared with the last run or judged).
- **Backhouse re-estimate (02:12 CDT, measured):** load rose to 66 after the cap (47 at 01:30), not down. high 1.2 to 1.4 ms in 1,159 s (about 5,800 s per ms; at 1.4 ms, about 22 h more to 15 ms); blob-fine 0.6 to 0.8 ms in 1,709 s (about 8,500 s per ms; about 34 h more). The Mac blob runs at about 1,100 s per ms at 5 ms, so Mac copies would finish high in about 5 to 6 h and blob-fine in about 11 to 12 h (estimate, twice the cells). Not started: the Mac's four slots are full, and lever 3 may change the code all four runs use.

- **Step 5, liquid LOX: the spray is in the engine** (`core/spray.*`, design and results in `docs/evidence/SPRAY.md`). Lagrangian parcels on the NIST liquid, Abramzon-Sirignano drops with film integrals, KH-RT, walls, Gosman-Ioannides dispersion (S8 passes), exact exchange through `Flow::exchange` inside ReactingFlow's step. S0, S1, S3, S4, S5, S8 pass; S0(b) fails as stated (NASA-7's extrapolation, kept). S6: 3.7 times the gas step on the 96-cell S3 run; `Spray::wallTime()` (7897b1d) now measures it on the S7 runs (`spray_wall_s` in the history).
  - **S2(c) judged (1c5f577): fails as stated** against Lafon et al. (2008) at 10 atm: d² slope +30.4%, lifetime −17.5%. Read against S2 (the miss is physics both our model and the S2 film lack: real-fluid film, finite liquid conductivity, H2 dissolution); nothing changed. Digitised tables in `docs/evidence/spray/lafon_2008_digitised.md`.
  - **Cantera's default transport fits fail below 300 K**; `thermo::transportFits(mechanism, 60)` fixes it. **inquiry2 should know this for the medium rework.**
- **S7, Mascotte A-10, pre-registered (8a5b8d0)** in `docs/validation/MASCOTTE_A10_PREDICTION.md`, driver `tests/mascotte_a10.cpp` (`crucible_mascotte_a10`), table `build/table_mascotte.bin` (check in `docs/evidence/mascotte/table_build_check.txt`). Blind: no measured value of the case has been read. §3 (Burick's bracket) and §6 (predicted values) are still to fill, each committed before the step after it.
  - **§3 filled (bda3ab2):** q 3.207 (Burick's Eq. 5), SMD low 223 μm, high 574 μm, core 13.7 mm.
  - **KH child-release fault found and fixed (61e8517):** children were released per 3% of the parcel's current mass (about 750 per blob; 326,061 parcels at 1.1 ms); now per 3% of its injected mass, as published (630 parcels at 1.1 ms). Spray tests rerun: as before.
  - **Running, all four from t = 0 at 00:30 CDT, CFL 0.8 (the 1 ms test held):**
    - Mac, slot.py: `blob` and `low` (`rr:223:3.207:13.7`) at 239×70; logs `/tmp/masc/{blob,low}.log`, outputs `/tmp/masc/{blob,low}_*`.
    - Backhouse, tmux `crucible_s7_high` (`rr:574:3.207:13.7`, 239×70) and `crucible_s7_blobfine` (478×70), tree `/home/greff/crucible_s7` (main at 61e8517's spray code; not `/home/greff/crucible`), logs `runs/{high,blob_fine}.log`. The grid check takes about 13 h.
    - Each run: settled if the halves differ ≤ 50 K at r ≤ 20 mm, otherwise extend by 10 ms (§2). Then `tools/mascotte_band.py <low> <high> <blob> <blob_fine> band.csv` gives §6.
  - The Mac and backhouse both carry other users' jobs, so times are loose.
  - **Pace at 00:50 CDT (measured):** Mac blob about 530 s per ms, low about 660 s per ms (both done near 04:30 to 05:30 CDT). Backhouse is at load 45 on 28 threads (other users): high about 2,100 s per ms (about 15 h), blob-fine about 3,650 s per ms (about 25 h). **Plan:** when the Mac runs finish, start a second blob-fine (and high, if backhouse's is still short of 15 ms) on the Mac through slot.py with the same arguments; use whichever copy finishes first and stop the other. Same code and inputs, so either copy is the run.
  - **S6, early (measured):** the RR low run's parcels level at about 6,500 to 7,800 from 1.2 ms (no runaway after the KH fix); the spray took 24% of its wall time to 1.9 ms (about 35% of the latest 0.2 ms); the blob run 8%. Below the gas, so S6's film-table trigger has not fired; the settled share goes in S7's record.
- **Wall functions: 2(a)'s resolved half restated like for like** (898f631, after the result, flagged): engine against tool on the same 64 rings, passes (+0.27%, −0.51%); the shared grid error is a stated limit.
- **fast ctest** (13 October): 7 of 9 pass; the failures are the known ones (`transport_verification` radial order 1.218 and 1.444, `spray_verification` S0(b)). The slow pipe test has still not been rerun.

## inquiry2's lane (worktree `.claude/worktrees/inquiry2`, 14 October, 02:15 CDT)

- **RL10A-3-3A tier 0 prediction committed** (`docs/validation/RL10A33A_PREDICTION.md`). No measured performance has been read; none may be until tier 1 is committed.
- **AMR step A done** (`docs/evidence/AMR.md`, *Results*): the explicit-geometry `Mesh` with `Mesh::coarsened`. 0(a), 0(b) and 0(c) pass; 0(d) fails as written (closure 1.9e-14 to 8.1e-14 at levels 1 to 3, from the direct formulas). Also the wall functions' y+ (`Flow::wallYPlus()`) and `chamber_study`'s end-of-run fix.
- **AMR step B1 done** (`docs/evidence/AMR.md`, *Step B1: results*; `core/amr.*`, ctest `amr_subcycling_verification`). Berger-Colella subcycling with refluxing; every level's kernels run over the whole level, and non-active cells are filled from below before each stage; inviscid only.
  - **1(a) and 1(b) pass.** Quiescent RL10 chamber on three levels: |u|/a0 at most 2.0e-13, each level at its own mesh's single-level value. 0(d)'s closure fix is not needed.
  - **3(a).** The reflection passes (3.4e-4 on 768 columns, falling). The order, 1.59 and 1.53 against 1.8, fails as written. A uniform control without AMR gives 1.63 and 1.50 on the same pulse (minmod at the peak), and the AMR's difference is a flat 0.55 of the uniform mesh's.
  - **3(b).** L1, the orders and the shock place pass. The pressure peak behind the interface (5.5e-3 to 1.4e-3 of the jump, against 1e-3) fails as written. It is a wave the shock sheds as it leaves the fine level: its peak is about 6e-3 of the jump at every resolution, while its integral falls as h. Composite conservation holds to 5e-15.
  - **Both failures stand as written.** The record gives the reading a later amendment could judge (the ratio to the uniform runs' own error). Whether to amend is the Director's call; nothing is applied.
- **AMR criterion 2 on B1, pattern fixed** (`docs/evidence/AMR.md`, *Step B1: criterion 2 results*; `tests/amr_budget.cpp`, run as `crucible_amr_budget 8e-3 4` with `CRUCIBLE_EQ_TABLE`, about 15 min). C1 `eqtt` with resolved walls, 32 x 6 plus one level, to 8 ms. The driver gained viscous registers, owned-cell boundary and source accounting, and the igniter region from the finest level. Single-level bits, B1's values and ctest (9 of 12; `spray_verification`'s (b) fails the same before and after) are unchanged.
  - **Budgets met against 1e-12:** mass 1.3e-14, energy 2.7e-15, momentum 9.3e-16.
  - **Not met as written:** ten times the single level's at 8 ms, for mass (6.0e-15) and momentum (1.1e-16). The errors are rounding walks that change sign, alike over the run (largest 1.96e-14 against 1.51e-14 for mass). A step C reading is recorded: the largest over the run against ten times the single level's largest.
  - **Not met as written:** the igniter, 1.4e-14 relative against 1e-14. The test first judged 1e-14 J and now judges the relative form. A compensated sum is planned before step C.
  - **Observed:** light-off passes Table A's 3500 K top on the AMR (coarse cells in the convergent section, up to 3524 K, at 0.34 ms), not on uniform 32 x 6 or 64 x 12. The counted Cantera fallback handles it; criterion 4 judges the solution.
- **RL10 development chamber** (`docs/evidence/RL10_DEV.md`, development only, pre-vaporized). Settled thrust in 3.3 min (uniform rings) and 5.9 min (stretched) on 4 threads; the stretched mesh binds in the throat's wall row. The implicit-time rule provisionally says no implicit step on this mesh. With AMR the target y+ decides it (300: about 1.7 h projected; 956: under 35 min), to be measured in AMR step D.
- **Table A for any pair: design and criteria** (`docs/evidence/TABLE_A_PAIRS.md`, nothing built). Diffusion moves elements up to 0.056 off the two-stream line in CH4/O2 flames, so the design is an exact element-potential solve seeded from a Bilger-f table, on `gri30_highT.yaml`.
- **AMR step B2 code committed, checks running** (`docs/evidence/AMR.md`, *Step B2: setup* and *what was built*; commit 7d28c17).
  - **What it is:** restricted loops over N_k, the listed faces, react over N_1, and the temperature-guess reset (in B1 too). `Amr::setRestricted(false)` selects B1.
  - **Smoke only so far:** 3(b) on 1 thread is cmp-identical between B2 and B1 with the reset, and B1's printed values are unchanged.
  - **The five checks are queued on backhouse:**
    - tmux `inquiry2_b2` (ours), script `/home/greff/crucible_b2/b2_checks.sh`;
    - trees `/home/greff/crucible_b2pre` (9fdab9b, the reference) and `/home/greff/crucible_b2` (B2);
    - each job waits for a 1-minute load under 24 (worker rule 16); it was 44 to 47 at 02:05;
    - progress is in `b2out/queue.log`, the comparisons go to `b2out/summary.txt`, and `b2out/done` marks the end;
    - about 1.5 h of running once the load allows.
  - **After the checks:** record *Step B2: results* and judge 3(c). The results are the cmp verdicts; check 2 against the recorded B1 values and against the 9fdab9b tree's own runs on backhouse; and the wall times, with the load.
  - **Process note:** one cmake build and two 4 s 3(b) smoke runs ran on the Mac outside slot.py on 14 October (worker rule 3).
- **For the spray worker:** inquiry2 has no pending or planned edits in `core/medium.*`, so a lean λ and D entry point there will not conflict. AMR step B2 edits `core/flow.*`, `core/transport.cpp`, `core/amr.*` and `adapters/reacting_flow.cpp` (react's cell list). The sub-300 K transport-fit finding is noted for any medium work.
- **`main` not fast-forwarded.** `worktree-inquiry2` is 18 commits ahead of `main` (e18c710) and contains it. `main` is checked out in the main checkout, and the worktree guard refuses commands aimed there. Fast-forward from the main checkout: `git merge --ff-only worktree-inquiry2`.
- **Next for inquiry2:**
  - B2's results and 3(c), then a compensated sum for the igniter's delivered energy before step C;
  - then Table A pairs step A, and the engine-spec schema for the blind-prediction database.

## Exact next steps (in order)

0. **Lever 3 (the Director's order, `COST.md` *Decision*):** the one-cell budget on a face-column cell (r 0.4 to 0.8 mm) before 2.311 ms: inviscid, conductive and species-diffusion enthalpy fluxes per face. The fluxes are in `core/flow.*` (inquiry2's), so either agree a read-only flux hook with inquiry2 or recompute them in the driver from the state through `Medium::transport`. Then the remedy: the open gas volume inside the RR core (block it or feed it from the core) if the transport is physical, or the wall and axis corner if numerical. Then decide the S7 runs: all four use this code; blob looks unaffected; low needs the fix; high (3,406 K at 1.4 ms) may reach the same zone. Then start the Mac copies of high and blob-fine (backhouse is far slower; see above). **Then lever 2**, spray subcycling: advance the parcels every N gas steps and apply their exchange as rates over the next N steps (the last step takes the remainder, so the deposit stays exact); N chosen by measurement, S1 to S5 at N = 1 and N. Also review the Table B agent's uncommitted files.
1. **S7**: when all four runs finish (check each for settling; the S6 spray share is in `spray_wall_s`): §6 (the band per §5), committed before any measured value is digitised. Then an agent digitises [N]'s CARS profiles and [R]'s OH image, and the judging goes into SPRAY.md S7. Copy the run CSVs into `docs/evidence/mascotte/`.
2. Carried from 6 October: criterion 3's reference (unsettled at 1 ms), the `turbulence_pipe_verification` rerun, criterion 0(e) (the offline wall-model table), the Table B builder, then the Back, Massier and Gier comparison (blind: never open the TR's results pages).

## For Ben

- **Hot-wall heat flux: which reference?** Standard SST, resolved to the wall, puts about 20% more shear and heat through a hot wall (gas 4.9 times the wall temperature) than the published variable-property wall laws. The DNS-backed scaling sides with the laws. Option A: a wall function that copies standard SST's hot layer (it carries SST's known error into the RL10 heat ledger). Option B: the semi-local wall model plus the published Hasan-Pecnik correction to SST in the engine and the reference, with criterion 2 restated against the corrected SST. Recommended: B. Raised with the Director 6 October (commit 329a17d, 00:43 CDT). **Answered: B** (6 October, before 00:55 CDT).

- **The omega-production item is closed.** The published SSTs form replaced the bespoke fix (the Director's instruction: published over bespoke). The cold start survives. The settled 32x6 chamber moved by 3e-5 in c* and 9e-5 in Isp (`PASR_C2.md`, *SSTs*).
- **The PaSR closure's laminar limit** (segregation-weighted kappa) is still flagged in TECHNICAL_PLAN step 7, as before.

## Where things stand

- **One project.** Post-pivot work is `main`. The docs are VISION_SCOPE.md (what), TECHNICAL_PLAN.md (how), RESEARCH.md (evidence, open model decisions, pre-pivot lessons in section 7) and this file (state). Verification records are in `docs/evidence/` ([index](evidence/README.md)), real-engine cases in `docs/validation/`, the MHD spike in `docs/parked/`. `archive/` is historical.
- **Ben, 4 October** (VISION_SCOPE *Changes since the pivot*): shifting equilibrium is the validated baseline. The chamber always evolves in time; a chemical run starts "valves open, ignite". Evaporation (liquid LOX with breakup) is built before RL10. All four regimes work by January.
- **RL10 ruling, 5 October morning** (TECHNICAL_PLAN *Compute budget*, Ruling): dual time stepping on the wall-resolved grid. Superseded at 21:45 by *Lightweight engine*: tables, wall functions as the default, block-structured AMR in index space, implicit time only by the decision rule.
- **Compute budget (measured, `docs/evidence/step_cost_2026-10-05.txt`, `chemistry_cache_2026-10-05.txt`).**
  - Before tables, chemistry was 96% of the single-thread viscous step. With Table A, the reaction is 10% to 15% of the inviscid C1 step, and the flow step (1.37 us per cell update, serial) is the rest (`TABLE_A.md`).
  - Backhouse (i7-14700K, 28 threads): 7.3e4 cell updates/s on 28 threads, about 1.7 to 2 times the Mac.
  - Local equilibrium is 2.0 times cheaper than finite rate on 1 thread, only 1.25 times on 10.
  - A chemistry cache would buy about 2 times on the step, at most 5 times; reuse is mostly a cell's own previous call.
- **Engine** (all verified, see the index): axisymmetric finite volume (HLLC, MUSCL, SSPRK2, device-thrust ledger), thermally perfect mixture, Cantera thermochemistry matching CEA, CVODES kinetics, Strang coupling, mixture-averaged transport, wall-clustered rings, SST-2003 with turbulent Chamber supplies, the PaSR closure, a subsonic partially non-reflecting outlet. C1 (chamber, "valves open, ignite") is done (`CHAMBER_C1.md`).
- **SST-2003 (`TURBULENCE_C2.md`).** Checks 1, 2, 3, 5 pass. Check 4 (pipe against `tools/sst_pipe_1d.py`) matches on nr 16, 32 and 64 (nr 64 on backhouse: within 3.2e-5 in u_b and c_f); its 1e-12 mass budget fails on all three at 2e-11 to 6e-11 (rounding of stored densities, inferred; the bound was corrected by a factor of 2). After the omega change, nr 16 and 32 still match (u_b within 2.07e-4 and 6.51e-5) but moved by more than round-off (u_b +4.4e-6 on nr 16). A reordered exact form gives the old values, so the pipe has non-zero divergence with the limiter active somewhere in the run. Where is an open item (the open ends of the nz 4 duct are a candidate, not measured).
- **PaSR (`PASR_C2.md`).** Criteria 1 to 3 pass (run 5). Criterion 4 was stated before any run (c44a933); its harness is `crucible_chamber_study frp|frt <nz> <nr> <end> <threads> <prefix>`.
  - Its smoke runs found two k-omega faults that stop any viscous turbulent cold start, in the closure run and the control alike. (1) The axis-row r^2 reconstruction was not bounded by the cell value: now clamped like the wall row (f09b077). (2) The omega-production choice above.
  - 32x6 to 8 ms (Mac): (a) to (d) pass for both. Light-off at 0.206 ms. The closure changes the light-off transient by up to 2.7% in vacuum thrust; the settled state is identical to 7 digits, because this chamber is premixed and nothing is segregated once it has burned. The closure is not yet tested on a non-premixed flame (interleaved rings do that).
  - 64x12 to 8 ms (backhouse): (a) to (d) pass for both. Light-off at 0.208 ms. The closure changes the light-off transient by up to 18.5% in vacuum thrust (at 0.244 ms); from 1 ms on, below 3.2e-5; the settled state is identical to every printed digit. c* +0.21% (was +1.54% on 32x6, falling with the grid as in C1), vacuum Isp -1.26% (does not fall; finite-rate recombination and viscous losses are the candidates, inferred).
- **Flame speed (`FLAME_C2.md`), done.** S_c against S_L 2.330470 m/s: +1.339% (40 um), +0.492% (20 um), +0.279% (10 um); S_d +0.044% at 10 um; all four criteria pass. Two caveats are on the page. The 10 um run has the open-face fix and the coarser runs do not. And the outlet's offset (an exact engine would show S_c +0.32%) is as large as the 10 um gap, so the check bounds the sum of the two.
- **Dual time (TECHNICAL_PLAN step 16, design only; parked at 21:45 on 5 October behind the *Lightweight engine* decision rule).** ESDIRK from ARK3(2)4L[2]SA against BDF2 (measured choice); local pseudo time; approximate Jacobian with exact residual; wall-normal block-tridiagonal lines, with point LU-SGS as a measured alternative; conservation from the stage residuals; criteria 1 to 5 stated in the plan. Read: Jameson's dual-time IRK paper, Yoon-Jameson LU-SGS, Kennedy and Carpenter (TM version) and Bijl et al. (AIAA 2001-2612 precursor). These added a PID step controller and a stage predictor. They also added two measured choices: ESDIRK3 against ESDIRK4 (Bijl et al. found fourth order most efficient), and the residual-weighted update against the last stage (Bijl et al. call the weighted update "potentially damaging"; it may amplify inner-iteration error on stiff wall cells, inferred). Not yet read: Jameson 1991.

## Next, in order

TECHNICAL_PLAN *Lightweight engine*, order of work:
1. ~~Thread transport~~: done 6 October (`THREAD_POOL.md`).
2. **Wall functions** (item 3): the hot-wall decision is answered (B) and in. Follow *Exact next steps* above, then set the AMR wall strip from the largest passing y+.
3. **Table B, the finite-rate manifold** (item 4): designed and probed (`TABLE_B.md`); the builder is next, after the wall-model table.
4. **AMR** (item 5), then the **RL10-like cold start** measured against the 5-minute target (item 6).
5. Optional or carried:
   - the pipe divergence with the SST limiter active (`TURBULENCE_C2.md`);
   - the C1 finite-rate and finer-grid reruns after the axis clamp;
   - PaSR criterion 4 on 64x12 with SSTs.
6. **C3, liquid LOX** (TECHNICAL_PLAN step 5; in progress, S0 to S6 run, see *Stopped at*), then **C4, RL10A-3-3A** under the blind protocol (a Sonnet agent extracts inputs only; pre-register vacuum Isp and thrust in git before any measured value is read).

## Backhouse (WSL2 Ubuntu, user greff; reach it with `ssh backhouse 'wsl bash -c "..."'` or stdin-piped scripts)

- `/home/greff/crucible`: the engine at 4e2ee5c (f10 ran from here; outputs `/home/greff/c2flame/`).
- `/home/greff/crucible_c4`: the tree committed as 03d56e3 (the criterion 4 runs; outputs `/home/greff/c4runs/`, done).
- `/home/greff/crucible_wf`: the wall-function study binaries (`wfs_hp_bin` has the corrected SST), Table A (`table_a.bin`), the law states, and the c2a runs.
- `/home/greff/crucible_c3`: the criterion 3 binaries (`study_bin`, `study_c3r`) and runs.
- `/home/greff/crucible_s7`: the S7 tree (main at 8a5b8d0 plus the spray timer 7897b1d; built with `-DCRUCIBLE_BUILD_APP=OFF`), the Mascotte table copied from the Mac (md5 identical), runs in `runs/`. Copy files in with `rsync -az -e ssh --rsync-path="wsl rsync" <files> backhouse:/home/greff/crucible_s7/<dir>/`.
- `/home/greff/crucible_b2pre` and `/home/greff/crucible_b2`: inquiry2's AMR step B2 checks (9fdab9b and the B2 code; tmux `inquiry2_b2`, outputs `crucible_b2/b2out/`). Files go in by `scp` to the Windows home (`/mnt/c/Users/ben/` in WSL), because the worktree guard refuses rsync's `--rsync-path`.
- Not ours, do not touch: gate_wg, ltaste2, battery2, flyapp-lib, lsug3, flyapp-m9r, s12spec_* (other users' jobs).
- Never write to `/home/greff/inquiry-project` (the old Rust project).

## Open items carried forward

- **Transport truncation order:** axial momentum 1.799 and radial momentum 1.223 and 1.446 fail the 1.8 criterion (`TRANSPORT_C2.md`).
- **Pipe mass budget** fails 1e-12 at 2e-11 to 6e-11 (rounding of stored cell densities, inferred).
- **C1 vacuum Isp** is not converging under refinement on 128x24; settling it needs 256x48 or a method-of-characteristics reference.
- **Low-Mach accuracy** of HLLC at chamber Mach about 0.15 (`LOW_MACH.md`).
- **Reacting step re-plans** almost every step (63759 of 68284 in the 32x6 criterion 4 run, measured).
- **Still missing:** save/load and control histories, editable geometry, retained comparison runs, probe picking, snapshot decimation, cross-platform builds and packaging.
- **Parked MHD findings:** `docs/parked/MHD_SPIKE.md`.

## Continue from here

- Numerical core: `core/flow.*`, `core/transport.cpp`, `core/medium.*`. Reacting coupling: `adapters/reacting_flow.*`. Chemistry: `adapters/reaction.*`. Thermo and the ideal rocket: `adapters/thermo.*`.
- Tables: `core/equilibrium_table.*` (Table A lookup, no Cantera), `tools/build_equilibrium_table.cpp` (build and the criterion 2(a) check), `tools/compare_histories.py` (two chamber histories at the printed times).
- Tests: `tests/`. Headless builds: `build/` (core plus Cantera, no app) and `build/native` (app).
- Cantera 3.2.0 is built from source in `~/src/cantera` (static library).
- Heavy Mac jobs go through `python3 ~/director/harness/slot.py run --label ... -- <command>`, bounded with `perl -e 'alarm shift; exec @ARGV' N`.
