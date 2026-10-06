# Session handoff

6 October 2026, about 00:50 · branch `main` (one line).

## Stopped at (5 October, evening)

- **Ben's direction, 21:45** (through the Director; TECHNICAL_PLAN *Lightweight engine*, verbatim there): a light engine, a cold start to full thrust in minutes on the Mac, heavy chemistry in tables built offline, adaptive resolution, no needless turbulence resolution. Dual time stepping (step 16) is parked behind a decision rule. Wall functions become the default wall treatment (step 17). The C2 to C4 dates were revised (T and L added; C4 still ends 8 November).
- **SSTs (8914753).** NASA TMR's SSTs production (P = mu_t S^2 in both equations) replaced the CRUCIBLE-only omega fix. Checks 1-3 and 5 pass. The pipe matches on nr 16 and 32 (within 1.0e-4 and 2.1e-5). PaSR criterion 4 on 32x6 passes: light-off 0.206 ms, c* +1.5468%, Isp -1.2393%. The 64x12 SSTs rerun was cut by the 21:45 direction.
- **Table A, tabulated equilibrium: done, all four criteria pass** (`docs/evidence/TABLE_A.md`, criteria committed first in db994a4).
  - Axes as built: f_H, T, ln rho_r, with 131 x 129 x 97 nodes; 144 MB, built in 11 s on the Mac.
  - Accuracy: 99th-percentile |dT| 1.06 K on 1e5 random states. The largest error on C1's own states is 0.83 K.
  - C1 64x12 cold start to 8 ms: 163 s, against Cantera's 2,048 s (12.6 times faster). Full thrust (4 ms) in 82 s, or 84 s on 1 thread.
  - Agreement with Cantera: c* -0.0053%, Isp equal to the printed digit, vacuum thrust -0.0006%.
  - The flow step is now 85% to 90% of the step.
  - Use it: `crucible_build_equilibrium_table build/table_a.bin 20 64 1e-4 129 97 1e-8 30 4`, then `CRUCIBLE_EQ_TABLE=build/table_a.bin crucible_chamber_study eqt ...`. The table file is not in git.
- **Threaded flow step: done** (`docs/evidence/THREAD_POOL.md`; the Director's 22:50 item).
  - One persistent pool (`core/pool`) runs the flow step and the reaction call.
  - The result is bit-identical to the serial code on 1, 4 and 6 threads: eqt 64x12 to 1 and 8 ms, frp 32x6 to 0.25 ms. Faces are computed alone and summed per cell in the serial order; the fused-multiply-add lesson is in the record.
  - C1 64x12 to full thrust (4 ms): 28 s on 4 threads (77 s before), 20 s on 6. To 8 ms: 56 s and 41 s (153 s before).
  - Transport (`core/transport.cpp`) threaded on 6 October, the same pattern: turbulent C1 (`eqtt`, SST, Table A) 64x12 to 4 ms in 79 s on 4 threads (169 s with serial transport), 57 s on 6; bit-identical to serial transport on 1, 4 and 6 threads (`THREAD_POOL.md`, *Transport on the pool*).
- **Stale notes fixed** (the Director's item 3):
  - the `archive/README.md` link now points to RESEARCH.md section 7;
  - the pre-pivot memory files are marked historical;
  - the backhouse memory now gives `/home/greff/crucible`.
- **Wall functions (6 October; `docs/evidence/WALL_FUNCTIONS.md`, commits d972db1, 64b0c9e, e10356a, b29adb2).**
  - Nichols and Nelson's law is in the engine on every no-slip face (side wall and plate rings). The first cell's k and omega are prescribed and booked; switch with `CRUCIBLE_WALL_FUNCTIONS=law`. With it off, the engine is bit-identical to before.
  - Pass: criterion 0 (a) to (d); 1(b) at y+ 99, 289 and 956 (c_f -2.5% to -3.6%, band 5%); 4 (cost +4.6% to +6.8%, band 10%).
  - **Fails: criterion 2(b)**, the hot pipe (T_axis/T_w 4.9) against standard SST. c_f is -17% to -19% and St -15% (bands 5% and 10%). A semi-local ODE wall model (Kawai and Larsson, `tools/eq_wall_model.py`) misses by as much. The reference SST's heated buffer layer sits 1.7 lower in u+ than the cold law under van Driest scaling, and DNS does not do this (Hasan, Elias, Menter and Pecnik, JFM 2025: standard SST's errors against DNS reach about 15% in velocity and 40% in temperature over 39 boundary layers; an earlier version of this line quoted 23% and 29%, which are not in the paper). **Decided 6 October by Ben: Option B**, the corrected SST in the engine and the reference, criterion 2 restated against it (`WALL_FUNCTIONS.md`, *Restated 6 October*).
  - Criterion 3's explicit reference (32 rings, y+ 1 on C1) is estimated at about 60 h; raised. 3(b), C1 to 8 ms with the law on 12 rings: c* -0.47%, Isp_vac -3.1%, 313 kW to the wall (34 kW with the no-slip wall), 218 s on 4 threads.
- Running when stopped: on the Mac, 1(b) at nr 167 (`/tmp/wf1b/law167.txt`) and 2(b) at nr 50 (`/tmp/wf2/law50.txt`), under slot.py. On backhouse, 1(a) in tmux `crucible_wf_res` and `crucible_wf_law` (`/home/greff/crucible_wf/c1a_*.txt`, about 6 to 7 h to steady state from 00:27).

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
2. **Wall functions** (item 3): built and partly passing (`WALL_FUNCTIONS.md`, *Where this leaves the criteria*). Blocked on the hot-wall decision (*For Ben*). Then: collect 1(a) and 1(b) nr 167; criterion 2(a); set the AMR wall strip from the largest passing y+.
3. **Table B, the finite-rate manifold** (item 4). Do the design and research pass before code. Open: a progress variable that keeps ignition delay, the table size, and mixing states off the manifold. Table A's lessons apply: use axes on which the diluent drops out, and cluster at kinks.
4. **AMR** (item 5), then the **RL10-like cold start** measured against the 5-minute target (item 6).
5. Optional or carried:
   - the pipe divergence with the SST limiter active (`TURBULENCE_C2.md`);
   - the C1 finite-rate and finer-grid reruns after the axis clamp;
   - PaSR criterion 4 on 64x12 with SSTs.
6. **C3, liquid LOX** (TECHNICAL_PLAN step 5), then **C4, RL10A-3-3A** under the blind protocol (a Sonnet agent extracts inputs only; pre-register vacuum Isp and thrust in git before any measured value is read).

## Backhouse (WSL2 Ubuntu, user greff; reach it with `ssh backhouse 'wsl bash -c "..."'` or stdin-piped scripts)

- `/home/greff/crucible`: the engine at 4e2ee5c (f10 ran from here; outputs `/home/greff/c2flame/`).
- `/home/greff/crucible_c4`: the tree committed as 03d56e3 (the criterion 4 runs; outputs `/home/greff/c4runs/`, done).
- `/home/greff/crucible_wf`: the wall-function study binary and the criterion 1(a) runs (tmux `crucible_wf_res`, `crucible_wf_law`, ours).
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
