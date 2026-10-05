# Session handoff

5 October 2026, 05:18 (parked at the Director's 05:30 line) · branch `main` (one line).

## Stopped at (5 October)

- PaSR criterion 4 passes on 32x6 and 64x12 for the closure and the control, after two k-omega fixes found by its smoke runs (`docs/evidence/PASR_C2.md`). The 64x12 pair ran on backhouse, 04:26 to 05:16.
- f10 (the 10 um flame) finished on backhouse at 04:46: all four FLAME_C2 criteria pass (S_c +0.279%, S_d +0.044%, drift 0.018%). Committed in 3c2c86d.
- Running on backhouse when parked: nothing of ours (checked with pgrep at 05:17; the tmux sessions crucible_f10, crucible_c4run and crucible_c4build have exited).
- Stale items noted, not edited (never touch `~/.claude` memory):
  - `archive/README.md` line 20 links to the removed `../REBUILD_NOTES.md`.
  - Several memory files predate the pivot: project_crucible, crucible_v13_architecture, crucible_pre_review_and_timeline, crucible_w2_chemical_slice, project_gpu_residency_s14_fork, project_crucible_commitments, feedback_fidelity_doctrine, feedback_resolve_everything_no_torch. MEMORY.md still says "VISION_SCOPE.md = source of truth (now v1.3)".
  - The backhouse memory (reference_gpu_box_backhouse) syncs to `/home/greff/inquiry-project`, the old Rust project, which must not be written. This engine lives at `/home/greff/crucible`.

## For Ben

- **A turbulence model choice (CRUCIBLE choice, not a DECISION NEEDED).** In omega's production only, the dilatation part of P uses omega for rho k / mu_t. It is identical wherever the SST limiter is inactive. Where the limiter is active in a strong expansion, the exact form is a sink that does not scale with omega, and it collapsed omega to 1e-138 in the divergent nozzle, stopping the cold start. The k equation keeps the exact P. The published alternative is NASA TMR's SSTs form (P = mu_t S^2 in both equations), which also drops the turbulent pressure's expansion work from k. Recorded in TECHNICAL_PLAN step 7; evidence in PASR_C2.md.
- **The PaSR closure's laminar limit** (segregation-weighted kappa) is still flagged in TECHNICAL_PLAN step 7, as before.

## Where things stand

- **One project.** Post-pivot work is `main`. The docs are VISION_SCOPE.md (what), TECHNICAL_PLAN.md (how), RESEARCH.md (evidence, open model decisions, pre-pivot lessons in section 7) and this file (state). Verification records are in `docs/evidence/` ([index](evidence/README.md)), real-engine cases in `docs/validation/`, the MHD spike in `docs/parked/`. `archive/` is historical.
- **Ben, 4 October** (VISION_SCOPE *Changes since the pivot*): shifting equilibrium is the validated baseline. The chamber always evolves in time; a chemical run starts "valves open, ignite". Evaporation (liquid LOX with breakup) is built before RL10. All four regimes work by January.
- **RL10 ruling, 5 October** (TECHNICAL_PLAN *Compute budget*, Ruling). A wall-resolved 1 s RL10 run is out of reach with the explicit step (the wall cells' acoustic limit is about 7e-11 s). Main line: dual time stepping on the wall-resolved grid (step 16), plus the cuts that cost no physics. Parallel: wall functions as a declared, switchable option, never the default, checked against the dual-time result on one short case before any use (step 17). Local time stepping is rejected. C4 slips by up to a week (recorded in the dated plan; A's buffer week is used up).
- **Compute budget (measured, `docs/evidence/step_cost_2026-10-05.txt`, `chemistry_cache_2026-10-05.txt`).**
  - Chemistry is 96% of the single-thread viscous step.
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
- **Dual time (TECHNICAL_PLAN step 16, design only).** ESDIRK from ARK3(2)4L[2]SA against BDF2 (measured choice); local pseudo time; approximate Jacobian with exact residual; wall-normal block-tridiagonal lines, with point LU-SGS as a measured alternative; conservation from the stage residuals; criteria 1 to 5 stated in the plan. Read: Jameson's dual-time IRK paper, Yoon-Jameson LU-SGS, Kennedy and Carpenter (TM version) and Bijl et al. (AIAA 2001-2612 precursor). These added a PID step controller and a stage predictor. They also added two measured choices: ESDIRK3 against ESDIRK4 (Bijl et al. found fourth order most efficient), and the residual-weighted update against the last stage (Bijl et al. call the weighted update "potentially damaging"; it may amplify inner-iteration error on stiff wall cells, inferred). Not yet read: Jameson 1991.

## Next, in order

1. **Pipe divergence (optional, small):** find where the pipe has non-zero divergence with the SST limiter active (`TURBULENCE_C2.md`, amended model note). It does not block dual time; skip it if dual time is the priority.
2. **C1 after the clamp:** the 32x6 equilibrium rerun (5 October) matches every printed result; only the round-off budgets moved (CHAMBER_C1.md). The finite-rate 32x6 run and the finer grids were not rerun; rerun them if a C1 number is quoted against post-clamp code.
3. **Dual time (step 16):** read Jameson 1991; then build it to criteria 1 to 5, starting with the open measured choices on a wall-clustered chamber.
4. **Leaner chemistry** (96% of the step): the cache measurement says about 2 times; a lean per-cell integrator is the other candidate.
5. **Wall functions (step 17)** as the declared option, after dual time gives its reference case.
6. **C2's remaining items:** turbulent inflow on the nozzle inlet; interleaved fuel and oxidizer rings (the first non-premixed test of the closure).
7. **C3, liquid LOX** (TECHNICAL_PLAN step 5), then **C4, RL10A-3-3A** under the blind protocol (Sonnet agent extracts inputs only; pre-register vacuum Isp and thrust in git before any measured value is read).

## Backhouse (WSL2 Ubuntu, user greff; reach it with `ssh backhouse 'wsl bash -c "..."'` or stdin-piped scripts)

- `/home/greff/crucible`: the engine at 4e2ee5c (f10 ran from here; outputs `/home/greff/c2flame/`).
- `/home/greff/crucible_c4`: the tree committed as 03d56e3 (the criterion 4 runs; outputs `/home/greff/c4runs/`, done).
- None of our tmux sessions remain. Not ours, do not touch: gate_wg, ltaste2, battery2, flyapp-lib, lsug3 (other users' jobs were using about 12 cores this morning).
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
- Tests: `tests/`. Headless builds: `build/` (core plus Cantera, no app) and `build/native` (app).
- Cantera 3.2.0 is built from source in `~/src/cantera` (static library).
- Heavy Mac jobs go through `python3 ~/director/harness/slot.py run --label ... -- <command>`, bounded with `perl -e 'alarm shift; exec @ARGV' N`.
