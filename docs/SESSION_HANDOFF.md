# Session handoff

23 October 2026 · branch `main` (one line). The spray worker (step 5) on `main`; inquiry2 works in `.claude/worktrees/inquiry2` (RL10 contour, AMR, table generalisation, and transport cost in `core/medium.*`, `core/transport.cpp`, `core/flow.*`, `contour.*`: not touched from here). The 6 October state is in this file's git history (`git show 15d6777:docs/SESSION_HANDOFF.md`).

## Stopped at (24 October, about 04:15 CDT; inquiry, D1 and D2 judged, criterion 5's cause found)

**The Director (03:20):** "judge D1 (exit 1, so find out why) and D2. Then the thrust offset... Report the cause before any decision on the implicit path. Stop at a clean commit by 05:30." The Director also confirmed that declining inquiry2's fast-forward request was correct, and will handle that merge at park.

**D1 and D2** (`docs/evidence/admdiag_2026-10-24/README.md`, *The result*): **the hypothesis is wrong.**
- D1 (27b135f) fails at 1.782 ms on the stage check rho omega > 0, not on the reflux's rho k. The cells are interior cells of the cold H2 expansion (z 1.03 to 1.04 m, r 0.13 to 0.15 m). The stage's rate on rho omega stays negative as the step halves to 3.3e-11 s. So (i) fails.
- D2 (`main`) runs to 1.80 ms with no rejection and no reflux clip, so (ii)'s second half fails. That `main` carries the same fault latent is inferred, not shown.
- **Next lead, not built:** a print-only patch giving that cell's rho omega rate per face and the SST source terms, with D1 rerun to 1.79 ms.

**Criterion 5's cause** (`IMPLICIT_PROTOTYPE.md`, *Criterion 5: the cause*): **the offset is in the explicit reference's reported thrust, not in the implicit.**
- An explicit run with equilibrium chemistry (flow, then relax) reports, at each history sample, the rates of the one step that lands on the sample. That step is shorter than a full step.
- The explicit's settled state is a fixed point of its two-stage mean, not of R. At its settled state, R(U) carries +7.1 N of net axial momentum at the injector face (64x12, CFL 0.4).
  - Why: the step carries unburnt propellant into column 0, and the projection is applied only after the step.
  - A full step's second stage cancels it. The shortened landing step's does not, which books +1.367 N.
- **The two marches agree on the thrust.**
  - e64's thrust, averaged over a sample period, is 1679.831 N; the implicit's is 1679.906 N (+0.0046%; the bar is 0.05%).
  - The explicit's step-zero limit is 1679.907 N.
- Halving the explicit step did not move the history's F_vac because the landing fraction changed from 0.806 to 0.613 and cancelled it, by chance. That is why the earlier note ruled the booking out; that note is withdrawn.
- **128x24, predicted at 03:50 and confirmed at 04:07:** the landing fraction is 0.113, so e128's landing step books +6.276 N. That is the whole 6.2 N gap. The first-stage imbalance is 7.07 N, the same as 64x12's. The offset "grew with the mesh" only because the landing fraction fell. Against e128's period mean the implicit is +0.0021%.
- **Criterion 5 as registered stays FAIL on the record.** On the corrected reference, the implicit is within the bar. Whether to re-judge is the Director's decision; the bar is not moved.
- **Wider consequence, not checked:** every explicit equilibrium-chemistry history's F_vac and Isp (and the logs' *simulated/predicted* lines) carries this landing-step error. It is between 0 and the first-stage imbalance (about 7 N, 0.4%, on C1 64x12). The fix is to report the rates' mean between samples. It is not built.
- Tools: `implicit_prototype/momdiag.py` (print-only, not for `main`). Runs: `c/c63t1`, `c/c63t2` (the tolerance ruled out), and `c/momdiag/`.

**Also measured:** criterion 6's step ratio is 496 (e256: mean explicit step 1.4126e-8 s over [0.5, 1] ms). The final build's fast ctest has only the 4 known failures (`implicit_prototype/ctest_2026-10-24_final.txt`).

**Done-files (Mac):**
- Lever B: `/tmp/cfl/<name>/done` for L4, L8, H4, H8 and H4s, and `/tmp/cfl/done` for all. None existed at 04:12. **L4 looks stalled:** its series (`/tmp/cfl/L4/L4_DEVELOPMENT_series.csv`, a row every 0.01 ms) took 71 rows from 01:46 to 03:10, then none after 0.71 ms. At 04:12 it still ran at about 390% CPU, and its log still gained Cantera below-200 K warnings (602 in all). Not diagnosed, not stopped; L8, H4, H8 and H4s wait behind it.
- D1, D2, e128, e256 and the ctest are done and judged above.

**Next:**
1. The Director's decision on the implicit path, with criterion 5's cause in hand. The other failures stand:
   - criterion 1 (pre-asymptotic order on P);
   - criterion 4's order (its floor was criterion 5's offset, so it should be re-measured against the period-mean reference);
   - m25's fixed point on 128x24.
2. When woken: lever B against (a) and (b) (`cfl_compare.py`); its evidence directory and SPEED.md 5c.
3. The history's rates as the mean between samples (a measurement fix in `Flow`/`chamber_study`), then a check of which earlier explicit F_vac and Isp numbers carry the landing-step error.
4. The rho omega print (above) for the 1.78 ms fault.
5. The driver's moving-window settle judge (C) before the central launch; the RL10 central run when S7 low and high pass.

## Stopped at (24 October, about 03:20 CDT; inquiry, the implicit prototype built and judged on case C)

**Ben (about 02:00):** "build the implicit prototype against the criteria you registered (d87815d)... stop at a clean commit by 05:30." The runs from 01:55 are not waited on; Ben wakes the session when they end.

**The prototype, against its criteria** (`docs/evidence/IMPLICIT_PROTOTYPE.md`, *Results*; three dated amendments, each before the runs it governs were judged):
1. **Criterion 1 (P, temporal order): FAIL as registered** (e8816dc). Second order, but pre-asymptotic at the registered steps (1.32 and 1.63 on 384 columns; 1.81 and 1.93 at smaller steps).
2. **Criterion 2 (budgets): PASS.** At most 4.3e-16 of the mass and 1.5e-13 of the energy supplied over 8 ms, in every run.
3. **Criterion 3 (inner tolerance): PASS** at E_t/10 = 8.975e-5.
4. **Criterion 4 (history): bars met at τ_C1/63 on 64x12** (RMS at most 0.17%, largest 1.59%), **and on 128x24 at τ_C1/150 only** (at τ/63 p_inj's and F_vac's largest are 2.07% and 2.72%). **The order bar FAILS:** the differences stop at the settled offset below.
5. **Criterion 5 (settled state): FAIL.** F_vac and Isp are 0.077% below the explicit's on 64x12 at every step (bar 0.05%), and 0.37% below on 128x24 (−1.03% at τ/25). Not the step, the tolerance, or the explicit's own step error (a CFL 0.2 run). Each scheme takes the other's state to its own within 10 μs, at the injector face. Cause not identified.
6. **Criterion 6 (256x48 to 1 ms at τ_C1/63): PASS.** 151 inner iterations per step, 7 rejections all recovered.
- **M1 to M3 (64x12, one thread):** one inner iteration costs 0.66 to 0.72 explicit steps. The speed-up to 8 ms is 6.6, 5.8 and 5.6 at τ_C1/25, /63 and /150. The inner iterations to 8 ms are about 31,000 to 36,000 at every step, so a larger step buys nothing: LU-SGS needs one inner iteration per 4 to 5 explicit steps of physical time. 128x24: 7.0, 5.6 and 4.8; 4.8 is the one at criterion 4's accuracy.
- **m25 (128x24, τ_C1/25) settles to a fixed point of the step, not of the flow:** p_inj 5% low, and its last short step jumps p_inj +3.7%. Not the tolerance (m25t at 1e-6 holds it too). Cause not identified; a candidate (the projection after the step not exactly idempotent) is untested.
- **Evidence:** `docs/evidence/implicit_prototype/c/` (`RUNS.md` lists every run, its build and its purpose; `c_compare.py` reads the gzipped histories).
- **Not done:** M4 (the projection to the lead) and M5 (finite-rate chemistry at the implicit step), T (the turbulent case).
- **Build changes** (uncommitted before this commit, all implicit-only; the explicit gate is `cmp`-identical on the final build, `implicit_prototype/gate/result3.txt`): Table A's projection inside the inner iteration (option), BDF2's known part as U^n + c (U^n − U^(n−1)), the throughput budget print, `CRUCIBLE_CFL`, `CRUCIBLE_IMPLICIT_WINDOW`, and the projection's Cantera failures taken as rejected attempts.

**Done-files (Mac):**
- D1 `/tmp/adm/d1/done` (exit 1), D2 `/tmp/adm/d2/done` (exit 0), both `/tmp/adm/done` (02:38:32). Logs `/tmp/adm/d1/d1.log`, `/tmp/adm/d2/d2.log`. Not yet judged.
- Lever B: `/tmp/cfl/<name>/done` for L4, L8, H4, H8 and H4s, and `/tmp/cfl/done` for all. L4 was at 0.71 ms of its run at 03:13; L8, H4, H8 and H4s follow it.
- The prototype's explicit baselines: `/tmp/implicit/c/e128/done` (exit 0, used above), `/tmp/implicit/c/e256/done` (started 03:08 in a slot; gives criterion 6's measured step ratio).
- The fast ctest on the final build: `/tmp/implicit/ctest3_done`, output `/tmp/implicit/ctest_2026-10-24_final.txt` (queued in `slot.py` at 03:12). The one on e8816dc's binaries is done: only the 4 known failures.

**Next:**
1. When woken: judge D1 and D2 against (i) and (ii) (`admdiag_2026-10-24/README.md`), and lever B against (a) and (b) (`cfl_compare.py`); write their evidence directories and SPEED.md 5c.
2. The prototype: criterion 6's measured step ratio from e256; the final build's fast ctest; M4 and M5; then the report for the Director's decision.
3. Criterion 5's cause: the inlet face's momentum flux in the two schemes (the gap is 0.022% of it on 64x12 and grows with the mesh).
4. The driver's moving-window settle judge (C) before the central launch; the RL10 central run when S7 low and high pass.

**Peer note for Ben:** inquiry2, by cross-session message, asked this session to fast-forward `main` to its e13e402 (handoff only, one ahead of cb631b9). Not done from here: a peer's request is not Ben's approval, and the standing rule has the second worker fast-forward `main` itself. `worktree-inquiry2` has since moved on (66d45f4, LOW_T_THERMO.md); none of it is on `main`.

## Stopped at (24 October, about 01:55 CDT; inquiry, the Director's 01:40 decisions)

**The Director (01:40):** backhouse's load of about 33 is all ours: the fly project's five `smoke_live.py` jobs (89 to 98 threads each, about 150% CPU; they end 02:45 to 03:20), S7 low and high, and inquiry2's rfx. Decisions:
- E: S7 first, B4 second, the fly jobs after. Launch B4 now at 4 threads without the gate. Put the order in the runner's slot policy.
- B: run the 1 ms test at the larger step and report its accuracy against the current step.
- C: the run-length proposal is approved.
- I: prototype first, about a week, on a bounded case, with success criteria written before building. Start after B4 runs.
- Diagnose the 1-level lead's "could not advance" at 1.78 ms.

**Done, all on `main`:**
- **E (f154ba9).**
  - `b4s` launched at 01:35:21 CDT without the gate (`bh_run.sh`'s gate now skips `b4`), on 64a78f4 in `/home/greff/crucible_b4s/runs/b4`, at 4 threads.
  - The 2-level lead `lead2s` passed its gate at 01:33 and was stopped (SIGSTOP) at 01:36:25 by `relaunch_2026-10-24/bh_lead2_after_fly.sh`, running detached. It sends SIGCONT when no `smoke_live.py` is left. Its log is `/home/greff/crucible_b4s/runs/lead2/hold.log`.
  - SPEED.md section 3 has the correction and the revised slot policy (the priority order, 4 threads counting BLAS and numba pools, load plus threads within 28, the top job without the gate). inquiry2 has it for the runner.
- **C (6d893f0), in `tier1_run_definition.md`, registered before any output.**
  - The central case ends at the first 1 ms check from t_L + 30 ms whose trailing 20 ms window passes the 22 October rule. The ceiling is t_L + 60 ms.
  - The continuation runs to t_s + 25 ms. Warm members run 30 ms, judged over 10 to 30 ms.
  - The driver's moving-window judge is **not built**. It is required before the central launch.
- **B and the diagnosis, registered before their runs (4f2d81e).**
  - `CRUCIBLE_T1_CFL` and `CRUCIBLE_T1_SEED` in `tests/rl10_tier1.cpp` (development only).
  - Lever B: SPEED.md section 5c (L4/L8, H4/H8/H4s; criteria (a) and (b)).
  - The diagnosis: `docs/evidence/admdiag_2026-10-24/README.md` (D1 on 27b135f, D2 on main; criteria (i) and (ii)) with the print-only `admdiag.py`.
- **I's criteria (d87815d): `docs/evidence/IMPLICIT_PROTOTYPE.md`.**
  - BDF2 dual time with LU-SGS, one level, one thread, Table A, on the acoustic pulse and the C1 cold start.
  - Criteria 1 to 6, measurements M1 to M5, and a proposed rule: build if the projected factor is at least 5, not if below 2.
  - Found after writing it: `Flow::rhs` already is the residual as a function, so day 1's refactor is mostly exposing it.

**Running on the Mac** (each job through `slot.py`; the fly jobs held the slots until 01:46):
- **Lever B** (`docs/evidence/cfl_2026-10-24/mac_cfl.sh`): L4, L8, H4, H8 and H4s in order, binary `/tmp/cfl/crucible_rl10_tier1` (4f2d81e). L4 started at 01:46.
  - Logs and series: `/tmp/cfl/<name>/`. Done-files: `/tmp/cfl/<name>/done` and `/tmp/cfl/done`.
- **The diagnosis** (`docs/evidence/admdiag_2026-10-24/mac_admdiag.sh`): both trees built in `/tmp/adm/old` and `/tmp/adm/new`, then D1 and D2 in turn. D1 was waiting for a slot at 01:49.
  - Logs: `/tmp/adm/d1/d1.log` and `/tmp/adm/d2/d2.log`. Done-files: `/tmp/adm/done`.
- **The prototype's explicit baselines** (`docs/evidence/implicit_prototype/mac_c1base.sh`): C1 `eqt` on 1 thread, 64x12 and 128x24 to 8 ms, then 256x48 to 1 ms. They start after `/tmp/adm/done`.
  - Logs: `/tmp/c1base/<name>/`. Done-file: `/tmp/c1base/done`.

**Backhouse (01:50 CDT):** load 29.0 with 5 smokes.
- S7 low is at 3.0 ms and high at 3.1 ms, about 195 s of wall per 0.1 ms. At that pace 55 ms is about 28 hours away (*inferred*), faster once the smokes end.
- `b4s` is at 226% CPU, still in its first 0.01 ms of output.
- `lead2s` is stopped (state T).
- rfx is running.

**Next:**
1. When the Mac runs end: judge D1 and D2 against (i) and (ii), and B against (a) and (b). Copy the logs and series into `admdiag_2026-10-24/` and `cfl_2026-10-24/`, and write the results (README, SPEED.md 5c).
2. Build the driver's moving-window settle judge (C) before the central launch.
3. The prototype, day 1: expose `Flow::rhs`, and read the C1 baselines for τ_C1 and the walls.
4. The RL10 central run, when S7 low and high pass (Ben's rule).

## Stopped at (24 October, about 01:30 CDT; inquiry, the Director's 00:50 note)

**The Director (00:50):** principle 1 stays as written: no steady solver, no pseudo-time that alters the transient. Items:
1. Estimate implicit time-accurate integration: the physical step, against today's explicit step, and the cost per step.
2. Measure backhouse's contention pinned and with our jobs paused, and propose a slot policy for inquiry2's resume runner.
3. Skip the chemistry exactly where nothing can react, with a test.
4. The finest level stays a numerics question for evidence: 2 against 3 levels on two members.

Fix the NaN first, then commit SPEED.md.

**All on `main`; the answers are in `docs/evidence/SPEED.md`** (its opening lists where each item is answered):
- **The spray NaN is fixed (26f6fd9).**
  - Cause: B_T's secant near saturation (B_M 7.7e7). A residual past the double range gave an infinite step, which the stop test accepted.
  - Test `spray_transfer`: 6 failures before, passes after (6.7e-16).
  - The case that stopped at 1.35 ms runs to 1.5 ms, exit 0.
  - SPEED.md section 6.
- **The chemistry skip: 4a7628f, extended and corrected in fc0c9b2.**
  - A molecule with a trace of its own atom (Y ≤ 1e-14) is now skipped, bounded by the dissociation equilibrium.
  - The mask reads the partial densities. A negative trace had been skipped, against the rule.
  - Test `chemistry_skip`: criteria 1 to 5 pass.
  - In the lead: 91% of the chemistry updates skipped, wall 1.5 to 2.2 times less, CPU 1.8 to 2.4 times less.
- **Contention (SPEED.md section 3).**
  - With our jobs paused, the segment took 91 s against 150 s. Pinned, it was 19% slower; nice made no difference.
  - Backhouse is 1.53 times the Mac's wall once our jobs are paused.
  - The slot policy (no pinning, at most 4 threads per job, total within the free virtual CPUs, critical path first) went to inquiry2 for the runner.
- **The implicit step (SPEED.md section 5a, lever I):** a 20 to 50 μs physical step, 5 to 50 times faster in the lead, 4 to 8 weeks to build (estimated, not built).
- **The fast ctest** (`docs/evidence/ctest_2026-10-24_skip_nan.txt`): only the 4 known failures (amr, amr_subcycling, transport and spray verification).

**Backhouse (checked 01:28 CDT):**
- **B4 and the 2-level lead are relaunched on main 64a78f4** (the NaN fix and the skip) in their own tree, `/home/greff/crucible_b4s` (the `REVISION` file; table `runs/rl10_s.bin`, hard-linked).
  - `relaunch_2026-10-24/bh_switch.sh` ran at 01:27:14 CDT. It held and stopped `b4r` (0.27 ms of its lead, c910fca, no NaN fix) and `lead2` (0.55 ms, `crucible_rst`'s 27b135f binary), and moved their job files to `/home/greff/resume/retired/`. Their directories and checkpoints are kept; a rebuilt binary would refuse those checkpoints anyway.
  - The new jobs are `b4s` and `lead2s` (`/home/greff/resume/jobs/`). Done-files: `/home/greff/crucible_b4s/runs/b4/b4_DEVELOPMENT_done` and `runs/lead2/lead2_DEVELOPMENT_done`.
  - **Both are queued behind the load gate** (`bh_run.sh` waits while the 1-minute load is above 24). The load was 33.7 at 01:28: S7 low and high and inquiry2's rfx (4 threads each), plus five python processes of another project at about 1.5 CPUs each. Each run's `.log` gets a "launch" line with the time and load when it starts.
  - This is the slot policy's rule 4: the critical path (S7) first. Starting B4 sooner is the Director's machine-priority call (E). The gate could wait until S7 ends if the other project's load stays.
- **S7 low and high** (`crucible_rst`, from 23:36 CDT, to 55 ms) were at 2.3 and 2.4 ms at 01:28, after 6,664 and 6,738 s of wall, about 225 s per 0.1 ms. At that pace they reach 55 ms around 10:00 CDT on Sunday 25 October (*inferred*; the pace depends on the load). About 5 minutes were lost to the 00:56 to 01:02 contention pause.
- **inquiry2's rfx** runs in `/home/greff/crucible_rfx`. It was paused from 00:56:48 to 01:02:01 CDT, and inquiry2 records that beside its wall times.
- The other job files (`b4` in `crucible_rst`, `dummy`, `t1cost`) all have done-files, so a reboot will not relaunch them.

**Next:**
1. The Director's decisions:
   - B (a 1 ms test at CFL 0.8);
   - C (run length);
   - E (machine priority) and the slot policy;
   - whether to prototype I (BDF2 with LU-SGS on the base mesh, about a week).
2. The RL10 central run, when S7 low and high pass (Ben's rule). The NaN no longer blocks it.
3. Still undiagnosed: the 1-level lead's "could not advance an admissible state" at about 1.78 ms.

## Stopped at (24 October, about 01:00 CDT; inquiry, the Director's 23:50 speed questions)

**The Director (23:50):** speed is the critical path; measure before deciding. (1) Profile tier 1 and say which constraint sets the finest level's step. (2) Test the hypothesis that the explicit acoustic CFL in a low-Mach chamber costs most of the time, and estimate local time stepping, dual time and low-Mach preconditioning from the profile, not built. (3) Thread scaling on backhouse at 1 to 16 threads. (4) Whether 60 ms of flow is physically required.

**The answers: `docs/evidence/SPEED.md`**, with its raw outputs, scripts and the print-only patches in `docs/evidence/speed_2026-10-24/`.
- **Wall time.** In the H2 lead (B4's 3 levels), chemistry takes 68 to 71%. Table B misses every cell (pure H2 lies outside its f_H range), so every cell goes to CVODES, about 15 μs of thread time each, to confirm that H2 alone does not react. The flow stages take 23%, the stable step 4%, AMR under 1%, the spray and I/O 0. In the hot base mesh (rl10dev, Table A), the flow stages take 73%, the stable step 21% and the chemistry 3.5%, at about 2 μs of wall per cell update on 4 threads.
- **The step.** In the lead, level 3 binds 94 to 98% of base steps, at the H2 jet and then at the throat's axis. The binding term is acoustic first, then eddy thermal diffusion. In the hot base mesh, the throat binds every step: acoustic 0.59 of the rate, flow 0.41. The chemistry, the spray and the wall functions set no steps.
- **The hypothesis holds only in part.** The binding cells are sonic or nearly so, so removing the acoustic term altogether lengthens the step only 1.4 to 2.2 times. Lever table and stacks: SPEED.md section 5.
- **60 ms:** physics needs about 9 ms for the gas and 20 to 40 ms after t_L with the spray (inferred); the schedule runs 105 ms.
- **Threads:** SPEED.md section 3.

**The measurement trees are not on `main`:** `/tmp/prof_tree` on the Mac and `/home/greff/crucible_prof` on backhouse. Each is main 1166721, plus 716c3d2's reflux clip, plus `stepdiag.patch`. The Mac tree also has `nandiag.patch`.

**B4r relaunched on main c910fca** (with the reflux clip) in `/home/greff/crucible_b4r` (the `REVISION` file).
- Job `b4r` under the resume runner, launched at 00:27 CDT. It was running at 00:36.
- Done-file: `/home/greff/crucible_b4r/runs/b4r/b4_DEVELOPMENT_done`.
- inquiry2 was told. S7's `crucible_rst` build is untouched.

**New blocker: the spray NaN after LOX opens.**
- A 1-level tier-1 run (lead shortened to 1 ms, 2 threads, `/tmp/prof_l1`) stopped at about 1.35 ms with "Spray: the film's transport is not positive at nan K".
- So the film temperature is NaN. 27b135f's guard reports it; the NaN's source upstream is unknown.
- Expected to stop B4 at about t_L + 0.35 ms, and the central run too.
- The print-only diagnostic run (`/tmp/prof_nan2`, NANDIAG lines in `nan.log`, done-file `done`): see below.

**Also undiagnosed:** the 1-level lead run stopped at about 1.78 ms with "could not advance an admissible state after 14 base-step reductions" (earlier tonight).

**Next, in order:**
1. The Director's decisions on SPEED.md's levers B (CFL 0.8), C (run length), D (levels) and E (machine priority).
2. Build lever A (skip the chemistry where nothing can react) if the Director agrees.
3. Fix the spray NaN from NANDIAG's output.

## Stopped at (23 October, about 23:45 CDT; inquiry, Friday night, the Director's 23:00 brief)

**Backhouse rebooted at 17:08; S7 low and high and B4 died with it.** The brief's items, in order:
- **The step-back** is in `docs/validation/rl10a33a/tier1_ensemble_plan.md`, the dated note of 23 October (7e91a4f).
  - At B4's measured pace: about 24,000 to 41,000 run-hours for the ensemble, 6 to 12 months on backhouse.
  - With levers 1 to 4: about 7 to 11 weeks, set by the 3-level central case.
- **(1) The restart: done and tested** (27b135f, `docs/evidence/RESTART.md`).
  - `CRUCIBLE_CHECKPOINT=<stem> CRUCIBLE_RESTART=auto` in mascotte_a10 and rl10_tier1, with a checkpoint every 30 min of wall time.
  - Four kill-and-resume tests on the Mac, each resumed on a different thread count. Every resumed run's engine state is `cmp`-identical to the uninterrupted run's, and so are its outputs.
  - The fast ctest's 11 FAIL lines equal da3555f's (`docs/evidence/ctest_2026-10-23_restart.txt`).
- **(2) The relaunch, on backhouse** in `/home/greff/crucible_rst` (27b135f, the `REVISION` file). It runs under inquiry2's resume runner (`/home/greff/resume`; the jobs `s7low`, `s7high` and `b4`; the scripts in `docs/evidence/restart_2026-10-23/`), so a reboot resumes the runs from their checkpoints.
  - **S7 low and high:** fresh from t = 0 at 23:36 CDT, to 55 ms, nice 10, 4 threads each.
    - Done-files: `/home/greff/crucible_rst/runs/s7low/low_done` and `s7high/high_done`. Each holds `exit N` and is written by the driver.
    - A fresh run to 55 ms took about 27 h before (*inferred* from the dead runs' 37.6 ms in 18.6 h).
    - The check that the new code gives S7's old numbers: the new history rows against the dead runs' (`/home/greff/crucible_s7z/runs2/{low,high}_history.csv`). Result below.
  - **B4 failed on main at 30 to 40 µs of the lead:** "AMR: could not advance an admissible state after 14 base-step reductions" (exit 1). Its log and series are `restart_2026-10-23/b4_27b135f*`.
    - This is inquiry2's reflux finding. At the frozen tags, level 3's reflux leaves ρk below zero in a coarse cell at the turbulent front, and no step reduction cures it.
    - The fix is f8e0f65 on worktree-inquiry2, registered in 8cc4572. Its checks were running on backhouse (`crucible_rfx`).
    - **B4 relaunches when f8e0f65 is on main.** Rebuild `/home/greff/crucible_rst` from that main, then register a new job (`b4r`, its own directory) with the same `bh_run.sh` body.
    - inquiry2 was told to add `refluxClips_` and `refluxClipped_` to `Amr::checkpoint` (core/amr.cpp, about line 653) when rebasing.
    - **The tier-1 central run fails the same way** until that fix is on main.
- **(3) The spray NaN: fixed at its cause** (`face_inlet/leak_2026-10-23/spray_nan/README.md`, Resolved).
  - The cause: tier 1 took its transport fits from the mechanism's 200 K limit, and H2's conductivity fit turns negative below about 97 K. The drops' films reach 61 to 67 K.
  - The fix:
    - tier 1 now fits from 60 K, as the spray's other cases do;
    - `Spray::filmTransport` throws on a transport value that is not positive;
    - the B_T residual throws on NaN.
  - The repro (lead 0, 2 levels) now runs clean to 0.02 ms.
  - Amendment in `tier1_run_definition.md`: G4's H2 viscosity moves 6.844e-6 to 6.833e-6, and F3's viscosity ratios move under 0.2%.
- **(4) B4's old run (ab3aa03, before the energy-leak fix c7e6bd1)** died at 1.09 ms, before t_L = 20 ms.
  - Its energy budget over the lead was 2.57e-5 at 0.5 ms and 6.34e-5 at 1 ms (limit 1e-11).
  - The face pressure was 1.009 bar at 1 ms.
  - No flash branch can be read from it. The fixed-code rerun is item (2).
- **Held note (c):** rl10_tier1 states the tag thresholds as frozen at s* = 1/4. `kNumericsFrozen` stays false while step D's y+ target waits.

**Next, in order:**
1. Relaunch B4 once f8e0f65 is on main (above).
2. Check S7's history rows against the dead runs'.
3. The RL10 central run: blocked on the y+ target (step D), the reflux fix, B4's flash branch and S7 low and high's scores.

## Stopped at (23 October, about 01:15 CDT; inquiry, Friday night, the Director's 01:15 message; stopped for the night)

**B4 now runs on backhouse; the Mac copy is killed.** The Director's reason: 17 hours on the Mac runs into Ben's day.
- **Done-file for the Director: `/home/greff/crucible_t1b4/run/done_b4`.** On the Mac: `printf 'ls /home/greff/crucible_t1b4/run/done_b4\n' | ssh backhouse "wsl bash -l"`.
- **The tree:** `/home/greff/crucible_t1b4`, ab3aa03 (the `REVISION` file), built there.
  - It uses `run/rl10_s.bin`, the same 21 October table as on the Mac (md5 1be4531e5c446ad70e3640b417fed10e).
  - The run script is `run/bh_b4.sh`: lead 20 ms, end 21 ms, 4 threads, `nice -n 10`, under `/usr/bin/time -v`.
  - Log `run/b4.log` (it ends with `exit N`). Series `run/b4_DEVELOPMENT_series.csv`.
- **Launched 01:05:30 CDT** in tmux on its own socket (`tmux -L crucible_t1b4`, session `b4`).
- **The load when it started:** load average about 24.8 on 28 threads. S7's four `crucible_mascot` processes run at nice 0; B4 at nice 10 and 4 threads, so S7 keeps priority. Someone else was compiling (`cc1plus`).
- **Its pace is not measured yet.** At 01:10 the series file was still empty.
  - That is expected: the driver flushes the series and prints a progress line every 50 samples, which is 0.5 ms of flow (`tests/rl10_tier1.cpp`, around line 518). The Mac's 20 rows were a partial buffer.
  - The first progress line in `b4.log` gives the pace.
  - The one proxy so far is the sub-200 K Cantera warnings: 293 in backhouse's first 5 minutes, against 339 in the Mac's whole run to about 0.21 ms. *Inferred, weak:* that suggests a pace no slower than the Mac's.
- **Two done-files:** the driver writes `run/b4_DEVELOPMENT_done` (`exit N`) when `runCase` returns. The script touches `run/done_b4` after the process exits, whatever the outcome, so watch `done_b4`.
- **The Mac copy's partial output** is in `/tmp/t1run/b/b4_mac/` (log and the 20-row series, to about 0.21 ms). Nothing of B4 runs on the Mac.

**Table B: the N gate, the Director's decision (i), registered in 9d35a2b and measured** (`docs/evidence/TABLE_B.md`, "Decision (i) measured"):
- **All gates pass on rl10_s:** N gate 0 failures, 2(a), 2(b) 0.028 K, 2(c), criterion 1. The rebuilt table is `/tmp/tb/rl10_s_ni.bin`, md5 6614a0ea3ef30dca1939f58ee4d07f5f.
- **N1 missed as worded: 115 nodes named, not 7.** The 7 registered nodes are there (N_eq > N_0). The other 108 are ties under the rule as registered (N_eq ≥ N_0 − 1e-15):
  - 90 are the η 0 nodes, where T_eq is 200 K and the branch has no span;
  - 18 are frozen lower branches at η 0.6 and 0.7, where N cannot change.
- **N2 and N3 met.** Only `n_gate_ok` changes, and only at the 7. No stored value moves. The table was written.
- **The upgrade rule** (registered): if the central run's "no fall in N" refusals exceed 1% of its CVODES calls, upgrade to (ii). It is read conservatively: every no-fall cell counts, whichever kind it is.
- **T2's wording miss** is recorded with its 320 moved-but-unchanged nodes in `table_b/rl10_s_fixes_t2_turnback_moved_2026-10-23.csv`.
- **Tier 1 still uses the 21 October table.** Moving it to the rebuilt one is a separate step, after B4.

**Next session, in order:**
1. **The species' low-temperature thermo (the Director's first item).** Extend H2, O2 and H2O below 200 K with sourced fits (NASA-9 or equivalent). The test: nothing above 200 K changes.
   - **Register before code** which hydrogen the supply is (normal, equilibrium, or para-rich frozen) and why. The Director's point: the RL10's hydrogen comes off the regenerative jacket warm.
   - My reasoning so far (*inferred*, sources still to find):
     - Stored LH2 is near equilibrium at about 20 K, which is almost all para.
     - Without a catalyst, ortho-para conversion takes far longer than the milliseconds the hydrogen spends in the jacket. So the jacket's outlet gas is likely para-rich and frozen, not normal (3 : 1) and not equilibrium at the local T.
     - Para and normal hydrogen differ most in c_p and enthalpy between about 50 and 300 K, which is where the 199.67 K supply sits. So the choice changes the injected temperature for a given enthalpy.
     - Check which form the JANAF and NASA H2 data describe (likely normal), and whether the supply boundary condition is set by T or by h.
   - Also still to settle: NASA-9 fits mostly start at 200 K, so the source below 200 K is likely JANAF or NIST. h2o2_highT and the table builder's Cantera fits must move together.
2. **Score B4** when its done-file appears: the flash branch first (face pressure over [t_L − 1, t_L] against 3.1 bar), then the gate's lines, the lead's cost, the flashes after t_L, and whether it passes t_L + 0.5 ms.
3. **Whether tier 1 moves to the rebuilt table** `/tmp/tb/rl10_s_ni.bin` (copy it into the repo's table place with its md5), after B4.
4. **A central run** only after acceptance, S7 low and high scored, and the Director's clearance.

## Stopped at (23 October, about 01:05 CDT; inquiry, Friday night, the Director's 01:05 message)

**The Director (01:05):** fix the Table B builder at the failing points of gates 2(a) and 2(c), without restricting the gates (50 against 80 kg/m³ is too close to call the dense row unreachable). The 188.7 K warning is a fidelity gap, not only a p0 artefact. Stop at a clean commit by about 02:45.

**The builder fix: registered in 460f767 before code, built and scored in 5d10ed8** (`docs/evidence/TABLE_B.md`, the two 23 October sections):
- **2(c) passes: 11 integrator failures down to 0.**
  - The cause: at ρ 80 no frozen start up to 6,000 K carries the node's energy.
  - Fix 1 continues the start family past 6,000 K's frozen composition, mixing it with the reactants at 200 K.
  - The 11 nodes now build (T_eq 3,386 to 4,369 K) and pass 2(b), E2, the turn-back and the N gate. One lower branch is frozen (reported, not gated).
- **The turn-back passes: 20 nodes down to 0 counted.**
  - The cause: all 20 were upper branches of almost zero span. Their falls were at most 4 units in the last place of s, which is round-off in evaluating s, not a turn-back.
  - Fix 2 counts a fall only above a derived floor of 1e-14 |c_eq|. The largest fall in the whole table is 2.5e-15 |c_eq|.
  - **T2 not met as worded:** `turn_back` moved at 340 nodes, not 20. The other 320 had passed already, and no stored value moved.
- **The N gate still fails at its 7 nodes (η 1, ρ 1e-3 to 0.017), so 2(a) still fails.**
  - The cause is physical, not the builder's: the strongly dissociated equilibrium has more moles than the reactants, so N rises along the lower branch.
  - The runtime already sends such a node's cells to CVODES ("no fall in N"), so this cannot give a wrong state.
  - **Decision for the Director**, (i) to (iii) in TABLE_B.md:
    - (i) gate the fall only where one exists, and count and name the nodes without one;
    - (ii) extend the runtime placement to a rising-N lower branch;
    - (iii) keep the gate as written.
- The tier-1 runs keep the 21 October `rl10_s.bin`, whose 11 dense nodes are unfinished. The rebuilt table is not written while 2(a) fails. Moving tier 1 to it waits for the N-gate decision, and is not done while B4 runs.

**B4 is still running** on `/tmp/t1run/b/tier1_f` (ab3aa03), launched at 00:51 CDT through `slot.py`.
- **Done-file for the Director: `/tmp/t1run/b/done_b4`**, log `/tmp/t1run/b/b4.log`, series `/tmp/t1run/b/b4_DEVELOPMENT_series.csv`.
- **Its pace (measured at 01:02 CDT): 0.21 ms of the 21 in 11 minutes.** At that rate the 20 ms lead takes about 17 hours (inferred; the pace may change as the flow settles), so B4 finishes on Friday afternoon or evening.
- The read order is in the 00:53 section below. Its face pressure over [t_L − 1, t_L] against 3.1 bar picks flash branch (a) or (b) (2cdd02c).
- **The temperature warnings, corrected:** 333 so far, the lowest 175.7 K.
  - The Director's reading: H2 is supplied at 199.67 K, below the 200 K floor of the thermo fits. So this is a fidelity gap in the species' thermo, not only a p0 artefact as the 00:53 section inferred.

**Next on this line, in order:**
1. **Score B4** when its done-file appears: the flash branch first, then the gate's lines, the lead's cost, the flashes after t_L, and whether it passes t_L + 0.5 ms.
2. **The Director's N-gate decision**, (i) to (iii). Then whether tier 1 moves to the rebuilt table.
3. **The species' low-temperature thermo (the Director's next item):** extend H2, O2 and H2O below 200 K with sourced NASA-9 or equivalent fits. The test: nothing above 200 K changes.
   - Points to settle when it is registered (*inferred*, not checked):
     - NASA Glenn's NASA-9 fits (McBride, Zehe and Gordon 2002) mostly start at 200 K, so the source below 200 K is likely the JANAF tables (Chase 1998) or NIST.
     - Below about 200 K, normal hydrogen (ortho 3 : para 1, frozen in a fast flow) and equilibrium hydrogen differ in c_p. Which one the supply is must be stated.
     - h2o2_highT and the table builder's Cantera fits must move together, or Table B and the engine's thermo disagree below 200 K.
4. **A central run** only after acceptance, S7 low and high scored, and the Director's clearance.

## Stopped at (23 October, about 00:53 CDT; inquiry, Friday night, the Director's 00:55 message)

**The Director (00:55) accepted the start amendment and the bound as built and decided the flash** (registered in 2cdd02c before B4 was read):
- The branch is chosen by B4's face pressure over [t_L − 1, t_L] ms (the gate's second mean).
- **At least 3.1 bar: (a),** rely on the lead, with the flash check. The check is built and met F1 to F3 (ab3aa03): `Ledger::flashedMass` is printed, and `run central` stops on any flash after t_L + 1 ms.
- **Below 3.1 bar: (b),** a sourced, superheat-dependent finite flash-rate model, its source and validity range cited, built under development mode. Not (c) or (d).

**B4 restarted at 00:51 CDT** on the flash-check build, `/tmp/t1run/b/tier1_f` (ab3aa03). The 00:47 launch, on the build before it, was stopped; its log is in `b4_old/`.
- The log is `b4.log`, the done-file `done_b4`, the prefix `b4_DEVELOPMENT`. 3 levels, 4 threads, the lead to 20 ms, then the LOX to 21 ms.
- **The pace:** under 0.5 ms in the first 3.5 minutes, so the lead takes hours.
- **Read in this order:**
  1. the gate's face pressure over [t_L − 1, t_L] against 3.1 bar, which picks the branch;
  2. the gate's three lines;
  3. the lead's wall time per ms;
  4. after t_L, the exchange-bound lines, the flashes and the flashed mass, and whether the run passes t_L + 0.5 ms.
- **A risk to the gate, seen in the first launch:** a Cantera warning that `equilibrate` met 188.7 K, below h2o2_highT's 200 K. Inferred: p0's stagnation in the cold H2. If the gate fails on p0 alone, that is the cause, not an unsteady lead.

**Table B's gates 2a and 2c are not started** (the budget). The failing nodes were already mapped (`TABLE_B.md`, *The gates on rl10_s, as written*):
- 2(c): 11 integrator failures, at ρ 80 and about 200 MPa, on the lower branch.
- 2(a): 20 turn-back nodes, up to 7.3e-7 against the 1e-9 bar, at η 0.1 to 0.2 and ρ 1e-3, 4.1e-3 or 80.
- The N gate: 7 failures, at η 1 and ρ 1e-3 to 0.017.
- The table is marked "not for use".
- **The design question for the Director:** fix the builder at those nodes, or restate the gates' domain to the RL10's reachable states and let the run's refusal counts show they are never met. Note that cold O2-rich mixing near 4.5 MPa can reach ρ about 50, so the ρ 80 row cannot simply be dropped (derived).

### The 00:35 message's work (kept for the record)

The Director (00:35) decided the spray overload:
- (1) start from a warm state already in the plan, recorded as a setup amendment before any RL10 output;
- (2) an exchange bound only as a step limit or exact subcycling, registered with the test that it changes nothing where it never binds.

The energy leak goes to inquiry2. Table B's gates 2a and 2c are on this line, after the above and before any central run. All in `tier1_run_definition.md`; values measured on the Mac.

- **(1), recorded as "the start amendment" (d93830e), built in 59a99d6.**
  - The plan names no imported warm field for the central case. TECHNICAL_PLAN step 1 reports "valves open, ignite" from ambient, and a declared partly developed start is a development fallback only.
  - So the warm start is an **H2 lead within the march:** the H2 opens at 0; the LOX opens at t_L = 20 ms; the igniter fires at t_L + 0.2 ms.
  - A gate at t_L: the outlet flow, p0 and the face pressure, over [t_L − 2, t_L − 1] and [t_L − 1, t_L] ms, within 0.1%. Central stops if it fails; development prints and goes on.
  - The judge's times count from t_L. Development sets `CRUCIBLE_T1_LEAD`; `_END` is the absolute simulated time.
  - **If the Director meant an imported gas-only field, that overrides this;** say so.
- **(2), registration B (d93830e), built and scored (59a99d6).**
  - `Spray::exchangeDt()` is const: f 0.25 of each owning cell's gas thermal and vapour times, level l counted 2^l times. The driver passes it only as `maxDt`.
  - Met: B0 (A0's md5), B1 (never binds; smallest bound / dt: S3 143, S8(c) 143, A1 143, A2 60) and B2 (A2 bit-identical with and without it). The dry run is still byte-identical.
  - **B3 not met.** At t_L = 0 the 3-level probe still stops in `Flow::exchange`.
- **B3's cause, measured: the LOX flashes. A step limit cannot bound a flash.**
  - The LOX enters at 101.17 K, above saturation below about 2.8 bar (derived). On its first move into the 1,000 Pa start, each parcel flashes about 30% at once.
  - The vapour deposited equals the flashed mass to 4 or 5 digits. Rate-driven vapour is at most 1.4e-13 kg per sub-step.
  - A parcel's flash is its mass (the flow over the parcel rate) times a fixed share, so no dt changes it.
  - Options (a) to (d) for the Director are in `tier1_run_definition.md`. Recommended: (a) rely on the lead if B4 shows no flash, else (b) a finite flash rate.
  - **The margin is thin:** the lead's pressure was estimated at about 3 bar, against about 2.8 bar saturation.
- **B4 is running** (launched 00:47 CDT through `slot.py`, 3 levels, 4 threads, a development run).
  - It runs the lead's 20 ms, then the LOX to 21 ms. Directory `/tmp/t1run/b`: log `b4.log`, done-file `done_b4`, prefix `b4_DEVELOPMENT`. The binary is `tier1_b4`, at 59a99d6.
  - Expected to take hours (*inferred*: 2 levels ran 0.02 ms in 8 s, and every chemistry update in pure H2 takes the fallback).
  - **To score it:** the lead's cost (wall time per ms); the gate's lines at t_L, especially the face pressure against 2.8 bar; the flash count and the exchange-bound line after t_L; and whether it passes t_L + 0.5 ms.
- **Next on this line, in order:**
  1. Score B4.
  2. The Director's choice on the flash.
  3. Table B's gates 2a and 2c for `rl10_s.bin` (`docs/evidence/TABLE_B.md`, *The small table rebuilt on h2o2_highT*: 2(a) N-gate and turn-backs, 2(c) CVODES failures at t = 0 at the dense end). Before any central run.
- **Unchanged blockers:** the energy budget (inquiry2), the restart reader, step D's reruns, then the y+ target, tags and `kNumericsFrozen`, protocol §4's two amendments, S7 scored, and the Director's clearance.
- **`spray_verification` fails S0(b)** at 1.245e3 J/kg against 1e3, the same as on 22 October. It predates this work and is not looked into.

## Stopped at (23 October, about 00:30 CDT; inquiry, Friday night, fourth brief)

The brief (the Director, 00:00 23 October): record two decisions, then build the driver's run mode in `tests/rl10_tier1.cpp`, spec and acceptance first. The RL10 central run launches only when the acceptance passes, S7 low and high are scored (Friday evening) and the Director clears it. The blind protocol is unchanged. All values measured on the Mac.

- **The Director's decisions, recorded:**
  - **A2:** the floor stays at 100 and the run is lengthened to 1.0 ms (5b7d1e6). Scored in 543d50b: met, with 355 moves to a finer owner. A3 is met.
  - **Step D on the new `spread`:** inquiry2 reruns its AMR checks and criterion 4 on b7f232b. This line does not (SPRAY.md, e7d126f).
- **The run mode is registered, built and partly accepted** (`tier1_run_definition.md`, *The run mode*):
  - registered in e7d126f, built in 53b82b4, with time corrections in 3e5aac8 and 3bb7f12;
  - the wall-face fix registered in 42d3da9, applied and scored in 762806e.
  - `crucible_rl10_tier1 run development|central`, and `selftest`.
    - `central` refuses while `kNumericsFrozen` is false. It is false.
    - Development takes its settings from the environment (`CRUCIBLE_T1_BASE`, `_LEVELS`, `_END`, `_THREADS`, `_OUT`, `_RATE_FACTOR`) and needs `CRUCIBLE_MANIFOLD_TABLE`.
  - **Met:**
    - R0, A0's md5;
    - R1(a), the dry run byte-identical;
    - R1(b), the Definition and gate lines identical to the 21 October log's 149 × 16 lines;
    - R2(a) and (b), exit 2 and 3 with no step;
    - R4, the settle judge's selftest.
  - **R3 is not met.** Two findings block it (below).
  - **The dumps are not a restart.** A restart reader is recommended before the central launch.
- **For inquiry2, two changes in your files (`core/transport.cpp`), both additive and checked:**
  - `BoundaryRates::plateHeat` and `Measurements::plateHeatFlow`, the i = 0 share of the wall heat (`transport.cpp:772`, `flow.cpp`; 53b82b4).
  - **The wall-face fix** (`transport.cpp:388-395`; 762806e). The wall-solve loop skips a face whose cell is outside `inRegion(q, 2)`.
    - With 3 levels, level 3 holds only the face column, so its side-wall cells were never filled. The loop then threw `Invalid wall-function layer`.
    - Unchanged without AMR. The 2-level probe's series and AMR dump are byte-identical.
    - Please rebase over both.
- **The two findings that block R3** (`tier1_run_definition.md`, the bullet after the fix's score; logs in `run_mode_2026-10-23/`):
  1. **This line: the start-up's explicit evaporation overdraws a fine cell.**
     - At 3 levels the run stops in `Flow::exchange` before the first 10 μs sample. A level-3 cell at 0.096 kg/m³ (the chamber starts at 1,000 Pa) took three times its own gas mass of O2 vapour in one step, and its energy went below the floor.
     - Options: (a) a step bound from the exchange (recommended; numerics, registered like the other spray changes) or (b) the start in the Definition (the Director's).
     - Next: register (a) with its criteria (A0's md5 if it never binds in S3 and S8(c), plus the 3-level probe), then build it.
  2. **inquiry2's lane: the gas energy budget on this Definition.**
     - R3(a) needs 1e-11. The 2-level probe gives 2.45e-7, all energy. With no spray and no AMR, F3's single-level step already ends at −1.011e-8 of the initial energy.
     - So the leak is in the gas path's boundary accounting: the open-area supplies, the plate's wall-function share or the wall-function heat. Inferred, not isolated.
     - Step D's and the face inlet's checks did not score it.
- **Cost so far:** at 2 levels, 31 base steps to 0.02 ms in 8 s on 4 threads, the spray 1.0% of it. At 3 levels it is not measured yet.
- **Table B's build gates 2a and 2c FAIL** for `rl10_s.bin`. The run mode prints this (known since 21 October, `docs/evidence/table_b/rl10_s_build_2026-10-21.txt`). It should be settled before the central run.

**S7:** still running. The done-file is `/home/greff/crucible_s7z/runs2/all_done`. Low and high are to be scored on Friday evening.

**The RL10 is not cleared.** What blocks it:
1. R3: the two findings above.
2. The restart reader (recommended).
3. Step D's reruns on the new `spread` (inquiry2), then the y+ target and tag thresholds, then `kNumericsFrozen`.
4. Protocol §4's two amendments.
5. S7 scored, and the Director's clearance.

## Stopped at (22 October, about 23:55 CDT; inquiry, Thursday night, third brief)

The brief (the Director, 23:30): amend A1(b) and name its second rounding source within 30 minutes; fix the jet under A4 and rerun A0 to A3; the driver's run mode after A4 passes. Everything is in `SPRAY.md`, *Spray on AMR*, from *The cause, found* on, with the evidence in `docs/evidence/spray_amr/` (README). All values measured on the Mac.

- **A1(b) is amended and met** (19355bb): relative difference at most 1e-12 under SST.
- **The second rounding source is named,** in about 2 minutes of the 30.
  - Under SST, `Amr::step` steps with its first `stableDt()` call's dt.
  - `advance`'s second call re-solves T by Newton, so the dt it would give differs in the last bits, and `Flow::step` takes its own call's dt.
  - Not chased further. A possible change is listed in SPRAY.md.
- **The jet's cause was `Amr::spread`, in the gas path, not in the spray** (28d87b1).
  - The old spread gave each child the parent's mass fractions, so composition fronts got inconsistent ghost states.
  - With no spray, an O2 stream through a static one-cell patch reached 262 m/s.
- **The fix went through two rounds.**
  - **The first, one common limiter over every component, broke AMR's 3(a) pulse reflection** (9869c4e). Rounding noise in ρu_r cut every component.
  - **The revision is registered and committed:** a component joins the common limiter only if it varies against its physical scale (ρ; ρ√(RT); p; ρω never), and otherwise keeps its own limiter.
- **The revision's results:**
  - A5 met, 53.9, 54.2 and 53.7 m/s;
  - A4 met, 2.4411e-7 kg and 53.6 m/s;
  - A0 met, md5 unchanged, and S1(b) is the only spray-test failure, as before;
  - A1(a) and A1(b) met, and A3 met;
  - AMR's ctest as before, with only Sod's values changed.
- **A2 is not met as written, on one coverage count.** The moves to a finer owner are 98, against at least 100. All budgets are met.
- **The fix changes AMR results generally.** Step D's AMR numbers and any AMR dumps were made with the old `spread`.

**For the Director:**
1. **A2's finer-move count.**
   - Options: amend its floor, or lengthen or move A2's pattern so that more parcels cross. Changing the threshold after seeing 98 is a criterion change, so it is the Director's.
   - My recommendation: lengthen A2's run, keeping the floor of 100.
2. **Whether step D's AMR runs need rerunning** on the new `spread` (inquiry2's line).

**Next (this line):** the driver's run mode in `tests/rl10_tier1.cpp`, which the brief puts after A4 passes, as it now has. The pattern is `tests/amr_budget.cpp`. The Definition and the G2 to G4 and G8 lines must stay identical to `tier1_dry_run_length_2026-10-21.log`, apart from decisions 6 and 8.

**S7:** still running. The done-file is `/home/greff/crucible_s7z/runs2/all_done`, expected around Sunday 25 October.

**The RL10 is not cleared.** The remaining blockers:
1. the driver's run mode;
2. A2's decision;
3. step D's verdicts, D0's ctest, and possibly reruns on the new `spread`; then the y+ target and tag thresholds;
4. protocol §4's two amendments.

## Stopped at (22 October, about 23:25 CDT; inquiry, Thursday night, second brief)

**The Director's decisions of 22:40 are recorded** (70682ae):
- S7 low and high: batch-means SE of at most 25 K, over 15 to 55 ms (`MASCOTTE_A10_PREDICTION.md` §4).
- The RL10 tier-1 settle rule, over 40 to 60 ms (`tier1_run_definition.md` and protocol §3).
- The budget: go.

**Spray on AMR is built and tested, but a defect makes it unfit for tier 1** (`SPRAY.md`, *Spray on AMR* for the design and *Spray on AMR: results*; evidence in `docs/evidence/spray_amr/`). All values measured on the Mac:
- **A0 met.** The uniform path is unchanged, `cmp`-identical, md5 `e95251ab…`.
- **A1(a) met bit for bit.** S3, frozen, single-level `Amr` against `Flow::step` with `Spray::advance`.
- **A1(b) not met as written.** S8(c) with SST: 1,544 of 2,168 steps' dt differ, and the gas differs by 1.390e-13.
  - The cause is the gas path: with no spray, single-level `Amr` is not bit-identical to `Flow::step` under SST. It differs from the second step, by 1.562e-13 at the end.
  - One source is `Amr::advance`'s second `stableDt`. A second source is not isolated.
  - The spray adds nothing larger: 1.180e-13 with the spray against 1.378e-13 without, on the matched driver.
  - **For the Director:** make AMR's single-level SST step bit-identical (an AMR change), or amend A1(b). Nothing applied.
- **A2 met as registered.** S3 on 16 × 6 plus 2 levels, 0.5 ms, 98 pattern-changing regrids, 338 moves to a finer owner and 1,344 to a coarser one. Every budget and ledger is at most 1.2e-13 against 1e-11.
- **But a defect was found (SPRAY.md, *A defect the printed comparison found*).** On a mixed-level pattern, two-way deposits drive a spurious jet at the fine patch's corner (113 to 130 m/s, against 50 m/s drops). Evaporation and breakup roughly double.
  - With every level active, the result matches the uniform 64 × 24 run (2.398e-7 against 2.399e-7 kg). One-way coupling and the gas alone are smooth.
  - **The spray on AMR is not fit for tier 1 until this is fixed.** Criterion A4 is registered for the fix, and the next diagnostics are listed there.
- **A3 met.** 1 against 4 threads, bit for bit.
- **ctest** `spray_amr_verification`: only A1(b)'s three checks fail.
- **Recorded beside the registered design:**
  - wallThrust is set at the end of the base step, not at the next one.
  - The reference driver is `Flow::step` with `Spray::advance`, because `ReactingFlow::step`'s extra `stableDt` is a rounding difference of its own (8.155e-15 at the second step).

**S7:** relaunched runs are still going. The done-file is `/home/greff/crucible_s7z/runs2/all_done`, expected around Sunday 25 October.

**The RL10 is not cleared to run.** What still blocks it:
1. **The mixed-level jet in the spray on AMR.** This is next. The order (in `SPRAY.md`):
   - dump p, T and u at the patch corner at 0.1 to 0.2 ms, in `a2_levels_probe` modes 0 and 1;
   - switch off the mass, momentum and energy deposits one at a time;
   - try a static 2-level patch.
   - Then fix it and score A4, which is registered, and rerun A0 to A3.
   - Build the probes with the `c++` line in `spray_amr/README.md`; each mode runs 1 to 2 minutes on the Mac.
2. The driver's run mode, not started. The pattern is `tests/amr_budget.cpp`: one `ReactingFlow` per level, with `react(dt)` in afterStep after `spray.afterStep`, plus `setRedo`.
3. Step D's verdicts and D0's ctest, then the y+ target and tag thresholds.
4. Protocol §4's two amendments.
5. A1(b)'s decision.

## Stopped at (22 October, about 22:40 CDT; inquiry, Thursday night)

The brief: (1) score S7; (2) if S7 passes or its failure is understood, freeze the RL10's tier-1 inputs so inquiry2 can launch the central run.

**S7 is not scorable as registered. The failure is understood, and the runs are relaunched to their registered ends.** The record is `MASCOTTE_A10_PREDICTION.md` §4's note of 22 October; the evidence is in `docs/evidence/mascotte/s7z/` (`analysis.txt`, `analyse.py`, the two field PNGs, windows, histories, logs).
- **The stop.** `bh_s7z.sh`'s runs stopped at about 17:00 CDT on 22 October, when backhouse's default tmux server exited (cause not found). low was at 51.5 ms, high 51.8, blob 27.1 and blob-fine 25.7.
- **No closed window meets §2's 50 K test (measured).** low 93/94/59 K, high 118/122/91 K, blob 621 K, blob-fine 381 K.
- **low and high are stationary in mass but unsteady at the 10 mm station.** Their liquid held is flat, and the outlet flow matches the inflow. The 5 ms half test may never pass there except by chance.
- **The blob legs are still filling with liquid,** with a time constant of about 15 ms (derived), so they settle at about 65 to 85 ms.
- **No measured value has been read or digitised,** and §6 is empty.
- **Relaunched at 22:26 CDT** by `docs/evidence/mascotte/bh_s7r.sh`: same binary and table, output in `/home/greff/crucible_s7z/runs2`, tmux on its own socket `-L crucible_s7r` (session `s7r`).
  - **Done-file:** `/home/greff/crucible_s7z/runs2/all_done`. Each run also writes `<name>_done`.
  - **At the old pace (derived):** low and high reach 55 ms around 18:30 CDT on Friday 23 October; the blob legs reach 85 ms around 10:30 to 14:00 CDT on Sunday 25 October.
  - **Poll:** `printf 'tail -1 /home/greff/crucible_s7z/runs2/*.log' | ssh backhouse "wsl bash -l"`.
  - **Score:** re-run `analyse.py` on the new windows (field files stay on backhouse).
- **For the Director:** a proposed amendment for low and high (§4's last bullet: batch means over the stationary span, standard error at most 25 K), made before any measured value is read. Not applied.

**The RL10's tier-1 inputs are frozen** (`tier1_run_definition.md`, *Frozen inputs*: the input list, plus the start, quantities, parcel rate, threads and warm start; the files pinned by SHA-256, the table data by its hash). The mesh's y+ target and tag thresholds are not frozen; they wait on step D's verdicts.

**The RL10 is not cleared to run.** What still blocks it:
1. **New: the spray does not run on AMR** (`core/spray.hpp:116` takes one `Flow`; AMR.md:18 left it to the spray owner, this line; nothing built). The central case needs 3 levels and liquid LOX together.
   - Next: register the coupling's design (parcels on the finest covering level, gas read there, deposits booked to owned cells, composite budgets closed; S1 to S5 rerun on levels), then code it.
   - This is the longest item.
2. The driver's run mode (this line), as before.
3. Step D's verdicts and D0's pending ctest (the Director; inquiry2's `score_sd.sh`). Then the y+ target and tag thresholds.
4. **New, the Director's call:** how a spray run that is stationary but unsteady meets protocol §3's 0.1% per 1 ms. In S7, no 1 ms window of 36 met it for the outlet flow (median 1.8% low, 4.9% high).
   - Proposed: the central run runs to the 60 ms ceiling regardless, and the judging rule is decided before its output is read.
5. Protocol §4's two amendments (inquiry2's), and the restated budget with the Director's go.
   - The spray's cost is now measured on S7: the whole step at 1.50 and 1.70 times the gas step's cost on the Rosin-Rammler runs.
   - The table's data hash was already checked (21 October, 02:43).

## Stopped at (22 October, about 02:30 CDT; READY TO CLEAR)

**The Director's second decision is done: b2's rules and round 3.**

**b2_crucible's RULES.md is committed (9294aa6, 02:03), before any CRUCIBLE run of these engines** (`docs/validation/engine_db/predictions/b2_crucible/RULES.md`).
- **What it is:** CRUCIBLE's full engine, as for the RL10, on b1's 27 engines at cfd6f12.
  - Missing geometry is filled by batch-wide rules, each value labelled printed, derived or inferred.
  - Each inferred source is marked confirmed, secondary or memory. The L* rows other than LOX/RP-1 and the wall temperatures are from memory, and are said so.
- **The rules:**
  - the throat's size by cycle class;
  - CR by the Humble fit;
  - L* by family;
  - SP-8120's convergent section and throat arcs;
  - a Rao 80% bell with 60 and 90% members;
  - the injection states and the preburner stream;
  - a generic coax or doublet by phase;
  - the wall temperature by cooling type;
  - p_c hit by iterating the mass flow;
  - b1's engine-basis and sea-level rules.
- **The prerequisites table:**
  - The 8 LOX/LH2 engines need the RL10 pipeline, plus the preburner stream for rs25 and le7a.
  - CH4 (3), kerosene (6) and the storables (10) also need their tables, a liquid-fuel spray and doublets.
- **Setup items to commit before each family's first run** (*Amendments*):
  - the SP-8120 Fig. 5 θn/θe readout;
  - the doublet drop-size correlations (Dickerson; not read yet);
  - each family's mechanism or surrogate (inquiry2's lane).
- **The runs wait for the RL10 pipeline.** The sealed step waits for b2's P.
- **The header says "about 02:10"; git says 02:03.** ENGINE_DB.md notes this.

**Round 3 is done (02:05 to 02:25; ENGINE_DB.md, *Round 3*).** Three Sonnet agents ran in the foreground. Their files passed the screen (0 schema errors; the new masked hits are labels and locators only).
- **The merge** is committed with this handoff. `engines.json` SHA-256 `0283a115...90e6`.
- **Now all-primary:** j2_sa208_504, lm_ascent_engine and h1b. h1b's p_c is now 700 psia, b1's member.
- **9 stay in S,** with the searches noted per field.
- **b1 and b2 still cite cfd6f12,** with S at 12. A later batch on the new version would have S at 9.

**Next:**
1. **The RL10 pipeline (inquiry2's lane)** gates every b2 run. After it, the LOX/LH2 family's pilot (rule 13) runs, with the preburner stream built first for rs25 and le7a.
2. **S7 at about 17:00 CDT,** when low's and high's 35 ms windows land (the steps below, unchanged). No measured value may be read before then.
3. **Optional, light:** a round 4 for the 9 S engines, from sources that need DTIC or AIAA access. That is for Ben, since it needs access; it is not urgent.

## Stopped at (22 October, about 01:55 CDT; READY TO CLEAR)

**The engine database's first batch, b1_tier0, is predicted and registered.** This is the Director's brief: Q1 (a), Q2 (a) now and (c) later.
- **RULES.md:** 246d43b, committed before any engine was run (`docs/validation/engine_db/predictions/b1_tier0/RULES.md`).
- **P:** e747f22. The SHA-256 of `predictions.json` is `e055f68f9bcd85eb966fe97213097d1f47049a6067eabd956cad5cecde6572bd`.
- **The register line:** da89b6e (`engine_db/REGISTER.md`).
- **What it is:** the CEA ideal rocket (RocketCEA 1.2.3) for the 27 tier0_isp engines.
  - Vacuum and sea-level Isp and c* for all 27.
  - Thrust and the mass flow for the 14 that have a throat.
  - Each quantity has a central value and a band from one-at-a-time members, as in RL10's tier 0.
  - It ran on the Mac in 4.7 s, using 41 MB.
- **The sealed step is held** until CRUCIBLE's own prediction on the same 27 engines (batch b2) is committed and registered. One sealing then scores both. No measured value has been read.
- **The groups:** S has 12 engines and P has 15. The rule is "any of p_c, the mixture ratio or ε from a secondary source only". It puts aj10_190_oms in S through its Wikipedia mixture ratio, where the brief's 11 had allowed it as primary.
- **Rule edges that the run met** (recorded, not changed; ENGINE_DB.md, *Batch b1_tier0*):
  - for N2O4 with A-50 or UH25, equilibrium is above 900 K even at O/F 0.02, so lr87, lr91 and viking_4b take 0.02;
  - for LOX/kerosene, the gas generator's equilibrium O/F comes out at about 0.08.
- **Time stamps:** RULES.md's header says "about 02:00 CDT", but git gives 01:44; ENGINE_DB.md notes this.

**Next for the engine database:**
1. **b2's RULES.md:** CRUCIBLE's prediction on the same 27 engines. tier1 is usable for none of them, so b2 needs batch rules for the contraction ratio, the chamber length or L*, the contour and the injector, and must state which CRUCIBLE path runs: the lightweight engine, or a reduced mode.
   - **Decision for Ben/Director:** what "CRUCIBLE's own prediction" means for engines with no chamber geometry. A design question.
2. **Round 3 (Q2 (c)):** primary sources for group S, run with Sonnet agents in the foreground, as round 2 was. New records are a new `engines.json`, and so a new batch; b1 is scored as made.
3. **The sealed step** comes only after b2's register line (PREDICTION_PROTOCOL step 5, with RULES.md's hold).

**S7 (unchanged, measured 01:48):** the four extensions are in tmux `inqs7z`.
- low has reached 9.5 ms, high 9.4, blob 8.8 and blob-fine 7.0, each after about 4.0 h of wall time. The backhouse load is 22.8.
- No `_to25` window has been written yet. The 35 ms estimate stands: about 17:00 CDT for low and high, and 23 October for the blob legs.
- The scoring steps are below, unchanged.

## Stopped at (21 October, about 22:55 CDT; READY TO CLEAR)

**S7 (Mascotte A-10) is not scored: none of the four 25 ms runs is settled, so each is being extended by 10 ms, as section 2 registers** (21e8bf9).
- **The runs at 25 ms** are in `/home/greff/crucible_s7x/runs`. The brief's `/home/greff/crucible_s7/runs` holds the old pre-change runs, which are not used. All four exited 0.
- **The half-to-half drift at r <= 20 mm, against the 50 K limit:** low 93.4 K, high 117.5 K, blob 620.7 K, blob-fine 380.6 K (measured). The tables and logs are in `docs/evidence/mascotte/s7_25ms/`.
- **low and high are steady in mass,** so their drift is the flow's own unsteadiness. The blob legs are still filling with liquid: 0.65 and 0.59 g held, rising 11 to 12 g/s.
- **The extension** restarts from t = 0, since the driver has no restart. It uses the same code plus `CRUCIBLE_WINDOW`, which writes consecutive 10 ms windows, each as it closes (`<name>_to25_*`, `_to35_*`, ...).
  - Launched 21:45 CDT in tmux `inqs7z` from `/home/greff/crucible_s7z`.
  - low and high run to 55 ms and the blob legs to 85 ms, on 4, 4, 4 and 6 threads.
  - Done-files: `runs/<name>_done` and `all_done`.
- **Pace (measured at 22:37):** about 2 to 2.2 times slower than the 25 ms runs. Those had 5, 5, 5 and 9 threads. Load now is 26, with three fly `smoke_live.py` jobs on the box.
- **The 35 ms windows, if this pace holds (inferred from the 25 ms runs' rates):** low and high about 17:00 CDT on 22 October; the blob legs about 36 to 39 h after launch, around 09:00 to 13:00 CDT on 23 October.
- **The blob legs may hit the 4-day timeout before 85 ms (guessed).** The windows written by then are kept.
- **Next for S7:**
  1. Take the first window of each run whose halves differ by 50 K or less.
  2. Run `tools/mascotte_band.py` on those prefixes.
  3. Commit section 6.
  4. Only then, a Sonnet agent digitises [N]'s CARS and [R]'s OH image.
  5. Judge per section 5 in SPRAY.md S7.
- Status: `ssh backhouse "wsl bash -l"` with `ls -la /home/greff/crucible_s7z/runs` and the logs' last lines.

**The engine database's round 2 is committed (cfd6f12, S1): 31 engines, 0 schema errors.** `engines.json` SHA-256 is `f6cda983…4c47`.
- **Added:** LE-5B-3, LE-7A, HM7B, Vinci and CE-20; OMS AJ10-190, R-4D, Aestus, Viking 4B and RD-253; H-1, RS-27A and RD-107.
- **Usable:** tier0_isp 27, tier0_thrust 14, tier1 none.
- **The screen:** one value redacted unread (viking_4b's MR-basis note), and the counts reconcile (`docs/evidence/engine_db_screen_round2_2026-10-21.txt`).
- **The throat agent improved J-2 only.** NASA SP-8120 prints a 14.7 in throat; the same source's 700 psi, against the record's 763 psi, is flagged. ENGINE_DB.md, *Round 2*, has the details.
- **Primary versus secondary sources:** of the 27 tier0_isp engines, about 15 have an operating point (p_c, MR, ε) from primary sources (16 if OMS's Wikipedia MR is allowed). The other 11 rest on Wikipedia, Astronautix or Braeunig (derived from the records' source fields).

**Next (exact): the batch's RULES.md** (PREDICTION_PROTOCOL §2) before any prediction. It waits on the questions in the reply of 22:55, which are restated here.
- **Q1 (tiers and blindness):** opening the sealed file ends blindness for every later tier. The recommendation is to commit a tier-0 prediction (RocketCEA ideal rocket plus declared ranges, as RL10's `tier0.py`) for every tier0_isp engine, then hold the sealed step until CRUCIBLE's own prediction on the same engines is committed. One sealing would then score both.
- **Q2 (secondary-only engines):** recommended: they count, but are scored as a separate group, with a round 3 to find primary sources.
- **Rules it needs** (each the same for every engine, with a central value and a range):
  - the missing pressure station;
  - the open-cycle MR basis, and the turbine flow fraction;
  - turbine exhaust into the nozzle extension;
  - throats from envelope or lip-unstated diameters;
  - missing injection temperatures (normal boiling point or storage);
  - an unknown contour (the divergence factor's bracket);
  - shifting versus frozen;
  - sea-level points and separation;
  - the CUIP rig, which is c* only;
  - HM7B's, H-1's and RS-27A's conflicting secondary p_c values.
- **Ben's direction (goals.md, 20 and 21 October):** "a central value plus a range from rerunning across declared input ranges and model options", after the RL10.

## Stopped at (21 October, about 03:47 CDT; READY TO CLEAR)

**The blind engine database's specs are committed (caf3051, S0).** The first compile run was killed at 03:19 (worker rule 18) and wrote no spec. The rerun's four Sonnet compilers ran in the foreground, 03:22 to 03:42.

- **18 engines:** LOX/LH2 5, LOX/kerosene 4, LOX/methane 4, storables 5. `engines.json` SHA-256 `c73b97f5…1389`, 0 schema errors.
- **The screen** (counts, then the masked view) found 59 hits and no outcome value, so there were no redactions (`docs/evidence/engine_db_screen_2026-10-21.txt`). Every record was then read in full and merged unchanged.
- **Usable:** tier0_isp 16, tier0_thrust 10, tier1 none.
- **Geometry gaps:** the throat is missing on 6 engines and rests on a lip-unstated or envelope diameter on 7. Only 3 records print the chamber, and none prints its contour. The full lists are in `docs/validation/ENGINE_DB.md`, *Compilation*.
- **Not done:** no `predictions/` and no `sealed/` yet. No RL10 source was opened.

**Next (exact):**
1. **Round 2 to pass the trust bar's 20 engines.**
   - Foreground compilers, the same `compile_rules.md` and `compile_brief_2026-10-21.md` (edit only the clock), new output files in `/tmp/engine_db_staging/out/`.
   - Candidates are the ones dropped:
     - LOX/LH2: HM7B, LE-7A, LE-5B, Vinci;
     - storables: OMS AJ10-190, Viking, Aestus, RD-253, R-4D;
     - LOX/kerosene: H-1, RS-27A, RD-107/108.
   - Also ask for primary sources for the 7 envelope or lip throats and the 6 missing ones (e.g. an NTRS SPS or LR91 report).
   - Screen as before (`check_specs.py --counts`, then `--masked`), then append to `engines.json` as a new commit S1.
   - `/tmp/engine_db_staging/digest.py FILE` prints a compact read of a screened file.
2. **The batch's RULES.md** (PREDICTION_PROTOCOL §2) before any prediction. It needs a rule for each of these:
   - the missing pressure station;
   - the open-cycle mixture-ratio basis;
   - turbine exhaust into the nozzle extension;
   - throats from envelope diameters;
   - secondary-only fields.
   It must also declare that the CUIP rig is c* only.

## Stopped at (21 October, about 03:05 CDT; parked for the night, READY TO CLEAR)

- **The Director merged worktree-inquiry2 into `main` at 02:47** (af91015, a clean merge: AMR step C and the corner fix (x)).
- **Main at af91015: the build passes and the fast ctest shows no new failures (measured, Mac, 02:46 to 03:00).** The build exits 0 with no warnings. The same 7 tests pass, and the same 4 fail: amr, amr_subcycling, transport and spray.
  - Every FAIL line of the rerun is identical to the 02:20 record, values included (diff).
  - New are 7 PASS lines in amr_subcycling_verification: step C's check 3, the moving-pattern case with 202 regrids.
  - Evidence: `docs/evidence/merge_inquiry2_ctest_af91015_2026-10-21.txt`.
- **The run definition follows step C** (`tier1_run_definition.md`, items 4 and *What it waits on*):
  - the tags are `Amr::setTags`;
  - the static face tag needs no engine change, since `Tags::custom` holds it from `level(l).mesh()` (`core/amr.cpp:364`);
  - what remains of step C for tier 1 is its checks 2, 4 and 5 (inquiry2's).
- **Parked:** nothing else small on the line to RL10 is outside inquiry2's lane tonight. The driver's run mode needs (a) and AMR with the wall functions first.

## Stopped at (21 October, about 02:45 CDT; night run, READY TO CLEAR)

**The Director's brief on f29c38c:** three decisions, recorded in the plan (`tier1_ensemble_plan.md`, *Decided ... about 02:35*, decisions 6 to 8; commit 60b1c23):
- **(6) The face inlet is (a),** the open-area sub-cell supply in `core/flow`. It is inquiry2's lane, so the interface spec is below, for inquiry2 to code next session, after S7 and AMR step C.
- **(7) The supply turbulence's length scale is 0.07 D_h at the slots** (inferred, the fully developed duct's estimate), with one pitch on the porous face. Its bracket, 0.035 to 0.14 D_h (guessed), is member 18: two warm runs, 4 to 5% of the budget, so it counts as cheap.
- **(8) The face's wall is at the transpired H2's 199.67 K** (inferred), and the adiabatic face is member 19.

**Done:**
1. **The length scale is wired into `crucible_rl10_tier1`** (`CRUCIBLE_T1_SLOT_LENGTH`, central 0.07; 60b1c23).
   - Measured: slot viscosity ratios 137 to 2,520 (68.5 to 1,260 and 274 to 5,039 at the bracket's ends), against 46,000 at a pitch. The porous face stays at 524.
   - The dry run's output is line for line the G8 log's apart from that line (measured: diff; `tier1_dry_run_length_2026-10-21.log`).
   - The ratios are the same at 2 and 4 MPa (measured). ρu is the fixed mass flux, so the nominal face pressure sets no tier-1 input (derived).
2. **The count:** members 11 and 11b are not run, and members 18 and 19 add 3 runs. So there are again 27 deviations, and the budget stays 15 to 200 run-hours. The cold-start ceiling falls to about 230, since 6 runs are re-meshed, not 9 (derived).
3. **The run definition, drafted** (`docs/validation/rl10a33a/tier1_run_definition.md`): what is fixed, what is proposed (start, end, quantities, mesh and a face tag, parcel rate, threads, warm start), and what it waits on.
   - **Found while drafting:** `Amr` refuses the wall functions and the SST property corrections (`core/amr.cpp:25-26`), and tier 1 needs both. That is inquiry2's AMR line, and it is now listed as a blocker before any tier-1 run.
   - **Proposed:** the face's heat flow is printed apart and is not part of the hot-gas wall's heat flow. *Inferred:* protocol §3's quantity is the jacket's heat, and the face is transpiration-cooled.
4. **The RL10 table's build is now in the repo** (`docs/evidence/table_b/build_rl10_table.sh`; it was only in `/tmp/rl10s/build.sh`).
   - Rebuilt at about 02:43 (16 s), its data after the provenance header is byte-identical to `/tmp/tbs/rl10_s.bin`, the table the RL10 checks used (measured: SHA-256 `14f7af12…a404` for both).
   - The headers differ only in the build time, the seconds and the mechanism's absolute path. So the prediction commit hashes the data, not the file.

### For inquiry2: the face inlet (a) and the face's wall, interface spec (21 October, about 02:40; to code after S7 and AMR step C)

This is a proposed interface. inquiry2 owns `core/flow.*`, `core/transport.cpp` and `core/amr.cpp`, and may change its shape, but should keep the properties and the gates.

**1. The open-area supply (decision 6).**
- **The switch:** `Definition::openAreaSupplies{false}`.
  - False keeps today's ring-centre rule unchanged, so every existing case stays byte-identical.
  - True applies the rule below to every supply of a Chamber.
- **The geometry.** Face ring j has node radii [ρ_j, ρ_{j+1}] at the face station.
  - Supply s's share of it is a_sj = π (min(r_o, ρ_{j+1})² − max(r_i, ρ_j)²), or 0 if negative.
  - Its fraction is φ_sj = a_sj / `axialArea(0, j)`, which is the full annulus (`core/flow.cpp:193-196`).
  - Σ_j a_sj = π (r_o² − r_i²) exactly. So `supplyArea_[s]` becomes the true area at any mesh.
- **The flux.** g_s = ṁ_s opening(t) / A_true, as at `core/flow.cpp:1076` but with the true area. Ring j's face flux is the sum of:
  - φ_sj × `supplyFace(s, g_s, ring j's interior state)` over every supply s;
  - (1 − Σ_s φ_sj) × the plate's flux (`core/flow.cpp:1080-1083`, pressure only).
  This holds for the conserved, species (g_s Y_s), ρk and ρω fluxes. So each slot's mass flow and velocity (g_s / ρ at the face), and hence its momentum flux and J, are its own at any mesh. Choking and the characteristic stay per supply, against ring j's interior, as today.
- **A zero-flow supply is a closed port:** pressure only, with no shear, heat or species flux. `massFlow` 0 is already accepted (`core/flow.cpp:126`), and `supplyFace` returns exactly this at g = 0.
  - The driver declares each LOX slot as one.
  - So the liquid's share of the face is not a 199.67 K wall; the parcels carry its momentum and heat (inferred).
- **Refusals under the switch:** drop "holds no ring centre", and test "supplies overlap" on the true annuli.
- **The wall share.** Ring j's wall share is w_j = 1 − Σ_s φ_sj; in tier 1 that is the post tips. The plate's viscous fluxes on ring j (resolved, or the wall law's τ_w and q_w) are multiplied by w_j.
  - Today a ring is all supply or all plate (`plateWallFace_`, `core/transport.cpp:155-173`; the wall distance's `plate`, `core/flow.cpp:445-447`).
  - *Proposed (inquiry2's call):* ring j counts as plate, for the wall distance and for the wall functions' prescribed first cell, only if w_j ≥ 1/2. This reduces to today's rule when slots align with rings. Elsewhere the cell's k and ω are solved, and its wall share passes w_j times the law's τ_w and q_w.
- **AMR:** each level computes its own φ from its own face rings, in place of the parent's supply (`core/amr.cpp:34-40`).
  - A coarse face area is the sum of its children's (`core/flow.cpp:248`), so a coarse a_sj is too, and the face's reflux balances.

**2. The face's wall temperature (decision 8).**
- **The field:** `Definition::plateTemperature`, a `std::optional<double>`.
  - Empty: the plate takes `wallTemperature`, today's behaviour.
  - Zero: adiabatic.
  - Positive: isothermal at that value.
  - Tier 1 sets 199.67 K (`kFuelT` in the driver), and member 19 sets 0.
- **Where it is read:** every wall read of `d.wallTemperature` takes the plate's value on plate faces. That is `core/transport.cpp:217-220` (ν_w, and the corrections' wall density and viscosity), 237, 467 and 761 (the wall law's T_w).
- **Member 19 also needs the wall functions on an adiabatic plate.** `core/flow.cpp:107-108` refuses it today ("the adiabatic law is not implemented"). The law's adiabatic form (q_w 0, T_w the recovery temperature) is inquiry2's to take from Nichols and Nelson and check (WALL_FUNCTIONS.md). Only member 19 needs it.

**3. Gates (proposed, stated before the code):**
- **F1, switch off:** the fast ctest's 7 pass and 4 fail, with sorted FAIL lines identical; rl10dev's short run byte-identical (`cmp`), as at 03cf8b4.
- **F2, switch on, supplies aligned with ring nodes** (every φ 0 or 1): identical to the switch off, to rounding.
- **F3, switch on, tier 1's layout** at 149 × 16, at 1192 × 128 and on AMR (149 × 16 with 3 levels):
  - each supply's face mass flow over a step equals ṁ_s opening(t) to 1e-12;
  - its face velocity equals ṁ_s / (ρ A_true) to 1e-12 on every ring it touches, so J is the element's own.
  - This is the driver's G6 rerun with (a).
- **F4, AMR:** the face's mass, momentum and energy fluxes summed over a coarse ring's children equal the coarse ring's, to 1e-12.
- **F5, the plate temperature:** empty and set equal to `wallTemperature` are byte-identical. At 0 (with the adiabatic law in place), the plate's q_w is 0 to rounding on every plate face.

**Also on inquiry2's line before any tier-1 run:**
- AMR carrying the wall functions and the property corrections (`core/amr.cpp:25-26`);
- AMR steps C and D;
- protocol §4's two amendments (choice 1's face, and decision 7's length scale).

**Next (straight line to RL10):**
1. inquiry2: S7, AMR step C, then this spec, then AMR with the wall functions.
2. Here, once (a) lands: set the switch, the LOX ports and `plateTemperature` in `crucible_rl10_tier1`, rerun G5 and G6 as F3, then code the run mode from `tier1_run_definition.md`.
3. Before the central cold start: rebuild the RL10 table with `build_rl10_table.sh` and record its data hash; measure the spray's cost on S7 and the gas-side cost on AMR (criterion 6); restate the budget.

## Stopped at (21 October, about 02:30 CDT; night run, READY TO CLEAR)

**The Director's brief on 30483ac:**
- the fallback bar is (b), split;
- Table B is central and Table A a member;
- (1) merge worktree-inquiry2 and run the fast ctest;
- (2) read the four queued jobs;
- (3) the tier-1 central case's setup as a dry-run build.

Commits 5b4679f, 316a719 (the merge), 03cf8b4, 64bab89, 0b4542e, fc6ede3 (the four jobs in TABLE_B, after the Director's cv562 note), e4574e1 and c3a402e (the dry run's drop size and members).

1. **The split bar (5b4679f):** recorded in TABLE_B.md (*The fallback bar split*) and the RL10 plan. C1 and the RL10 check both pass under it. Table B is the central case's chemistry.
2. **The merge (316a719):** clean, from the main checkout.
   - The fast ctest gives the same 7 pass and 4 fail (amr, amr_subcycling, transport, spray).
   - The sorted FAIL lines and the per-test failure counts are identical to the 01:00 run (`docs/evidence/merge_inquiry2_ctest_2026-10-21.txt`).
3. **The four jobs** ran on backhouse (load was low; the Mac was RAM-bound). See TABLE_B.md, *The four queued jobs*:
   - **machine control:** the Mac's steps, p, c*, Isp, F and last-ms tallies to every printed digit (budgets differ at round-off);
   - **the Δη 0.05 table at O/F 5.055:** c* −0.029%, Isp −0.032% against CVODES. So halving Δη cuts the gap by about 0.1% of Isp;
   - **the replays (both O/F):** the flame-column cells (1, 0) and (1, 1) are on the curve in the CVODES samples (dY at most 2.4e-4) and the face column is off the manifold (250 to 600 K), so the verdicts stand;
   - **the O/F 5.62 pair (CVODES done 02:15):** c* −0.056%, Isp −0.133%, F −0.132% against CVODES (0.2% bars: pass); defects 0 in every region; by cause 26.8% in the chamber, 0 elsewhere; 188 s against 1,389 s (7.4 times, shared machine). **Pass** under (b).
4. **(3) The tier-1 dry run (03cf8b4, 64bab89; `docs/validation/rl10a33a/tier1_setup.md`):** `crucible_rl10_tier1`, with the contour shared with rl10dev, whose output is byte-identical.
   - **G1 to G3 pass** at the tier-1 meshes 149 × 16 and 1192 × 128.
   - **G4:** V_L 6.05 m/s, V_H2 177 m/s, J 2.68, Burick's group 5.03.
   - **G5 and G6:** `Flow` refuses choice 1's true slots below 1,280 rings, and the slots need over 5,120 rings to meet 5%. The merged bands pass at 1192 × 128, but with J at 5 to 34% of the elements'.
   - **G7 fails,** so decision 5's branch is applied, dated before any run (plan 0b4542e). Burick's curve has no recess term, so members 11 and 11b are not run either; the recess's effect is a declared limit (c* and Isp low).

5. **The dry run continued (e4574e1 rule, then c3a402e results; `tier1_setup.md`, from *The LOX drop size*):**
   - **G8, the LOX drop size.** Mascotte's Burick rule is applied at the ring's group of 5.03, inside Burick's data. D̄/D_L is 0.148, q 3.207, so the **SMD is 447 μm**. Member 5's ends are 263 and 581 μm: the scatter's stated width, ×1.3 and /1.7, since three markers near group 5 cannot bound a band. Porcheron's core length is 6.84 mm. `Spray` now builds Rosin-Rammler injectors with these values.
   - **All ten geometry and flow members** (8, 9, 10, 12, 13, 14) were dry-run at 1192 × 128. G2, G3 and G8 pass, and G6 passes for the merged layout.
   - **Two rules were stated after the member runs and before any run** (engineering, resolved here, reversible):
     - (1) r_1 is the smallest radius at which ring 1 holds its slots and the rings clear each other. It binds for d_L 4.78 mm, t 0.37 d_L, h 0.80 d_L and r_8 = R_c − p; the central case is unchanged.
     - (2) Burick's curve continues past its table at Mascotte's slope. This is needed only at h 0.09 d_L (group 17.0, extrapolated and flagged).
   - **Seen:**
     - Member 10 spans J 0.18 to 29. With merged bands at h 0.09, J at the mesh is 0.8 to 6% of the elements' value: the face-inlet decision's sharpest case.
     - Burick's reading wobbles at groups 4.5 to 5.6 (digitisation noise), so members 9 and 14 move the SMD by under 4%.

**The face inlet:** decided by the Director after this stop, (a) (the 02:40 section above; the plan's decision 6).

**Next (straight line to RL10), while the decision is open:**
- date the supply turbulence's length scale (finding 2, 0.07 D_h at the slots) and wire it into the driver;
- then the tier-1 run definition (settling, AMR tags, parcel rate) once inquiry2's corner fix and the inlet land;
- the dry run is otherwise complete for every member the environment covers. Members 2, 3, 4, 7, 15 and 16 each change one Definition value and no geometry.

S7/Mascotte is still held behind inquiry2's corner fix.

**Backhouse:**
- `/home/greff/crucible_inq` (a git archive of 5b4679f, built Release), runs in `runs/` (tmux `inq_r10`, finished);
- Mac copies in `/tmp/r10bh/`.

## Stopped at (21 October, about 01:55 CDT; night run, READY TO CLEAR)

**The Director's brief on eca004b:**
- (1) Table B over the RL10's declared ranges: recheck the gates, and run a C1-like check at RL10 conditions against CVODES with the same 0.2% bar.
- (2) Both C1 runs to a settled state.
- (3) One line per known ctest failure.

Commits 33143d0, ed76456, c1c9f93, 9b5a65d. The results are in TABLE_B.md, from *R3 at a settled state* on.

**The fallback bar: decided by the Director after this stop, (b), the split bar** (recorded at about 01:48 in TABLE_B.md, *The fallback bar split*, and in the RL10 plan).

**Results:**
1. **(2) C1 settled (c1c9f93):** both runs to 8 ms, thrust drift 7.6e-8 and 8.5e-8 per ms.
   - c* −0.061%, Isp_vac −0.087%, F_vac −0.087%, injector face −0.047%.
   - 100% table in every region over the last ms, with 0 refusals.
   - 76 s against 937 s (12.3 times).
   - The 4 ms gap in Isp and F (−0.13%) was partly the CVODES run's unsettled end (inferred; the runs are deterministic to every printed digit).
2. **(1) The RL10 table (ed76456):** f_H at O/F 4.5 to 6.5, Z_N 0 and 0.3, Δη 0.1, 9 ln ρ nodes over 1e-3 to 80. The band is O/F 4.993 to 5.947, including the 5 to 15% transpired bracket's ring O/F, and 2.5 to 4.5 MPa, our tier 0 widened (guessed). It took 19 s.
   - **Pass:** E2 (2.06 K), 2(b) (0.028 K), criterion 1.
   - **Fail as written, at the C1 table's own places:**
     - the turn-back, 7.3e-7 at η 0.1 to 0.2;
     - 11 integrator failures at ρ 80 (derived: 190 to 215 MPa);
     - 7 N-gate failures at η 1, ρ ≤ 0.017.
   - None of these is reached live: 0 refusals.
   - **Spacing (static error):** η carries it, ρ does not. Halving Δη takes the chamber from −9 K to −3 K and the mid-nozzle from −21 K to −5 K.
   - But the Δη 0.05 table adds a turn-back of 6e-5 at η 0.15 and 1.1 to 4.7 MPa (1,432 K burnt gas), inside the band. Cooled wall gas in tier 1 would meet it.
3. **(1) The C1-like check (9b5a65d):** rl10dev, fr, 76 × 8, 7 kJ igniter, both settled at 20 ms (drift 1.2e-5 per ms). The 0.2% bars pass and the chamber's fallback bar fails, as above.
   - 290 s against 1,973 s.
   - The dated note is in the RL10 plan.
   - Seen, not chased (inquiry2's flow lane): the flame column's axis cell sits at 3.90 MPa between 2.76 and 2.89 MPa neighbours, in both runs.
4. **(3) The four known ctest failures:** all are expected, none is debt (a Sonnet agent's read, spot-checked).
   - **amr_verification:** criterion 0(d), geometry closure at fine levels, 8.1e-14 against 1e-14. Floating-point rounding, recorded as not met by decision (AMR.md:159-170, 13 October).
   - **amr_subcycling_verification:**
     - 3(b) is the pressure wave a shock emits crossing the coarse/fine jump. It does not converge, as for any cell-size jump.
     - 3(a)'s order 1.59 and 1.53 is the base scheme's own order on that pulse.
     - AMR.md:229-257 notes the criterion was written too strictly.
   - **transport_verification:** the radial momentum's order near the axis, 1.218 and 1.444 against 1.8, the known near-axis mechanism (TRANSPORT_C2.md).
   - **spray_verification:** S0(b), LOX NASA-7 against NIST, 1,245 J/kg at 75 K against a 1 kJ/kg band. The refit was declined on purpose: it would move the error into the hot-gas ledger, about a 1.4 K cold-film offset (SPRAY.md:138-144).

**Running, detached, each writing `done` in its directory** (queued behind inquiry2's and fly2's slot jobs at 01:44):
- `/tmp/r10tb21`: Table B on the Δη 0.05 table (`/tmp/tbs/rl10e21_s.bin`), the same case. It measures what the mid-nozzle's −21 K costs in Isp_vac.
- `/tmp/r10rep`: R1, R2 and the static error replayed on the CVODES run's 60,800 samples (its last 2 ms).
- `/tmp/r10cv562` and `/tmp/r10tb562`: the same check at the ring elements' O/F 5.62 (`CRUCIBLE_RL10_MIXTURE_RATIO`), CVODES and Table B.
- Scripts: `/tmp/r10_run.sh`, `/tmp/r10_run562.sh`, `/tmp/r10_replay.sh`.

**Next, in order:**
1. Read the four done-files. Append the Δη 0.05 run, the replay and the O/F 5.62 pair to TABLE_B.md, under *The C1-like check at RL10 conditions*.
2. Apply the Director's answer on the fallback bar. Under (b), restate the bar in TABLE_B.md (dated, before any further run) and re-judge both checks. Under (a), Table B is not the central case unless the face states get a design.
3. `tier1_skeleton.md`: the cup's resolution criterion, before any tier-1 run (decision 5).
4. R1's τ bar (8.7e-9 against a guessed 1e-9), and the per-cell entropy agreement with inquiry2 (below, unchanged).

**Files not in git:**
- `/tmp/tbs/rl10_s.bin`, `rl10e21_s.bin`, `rl10r17_s.bin`, built in 19 to 37 s by `/tmp/rl10s/build.sh`;
- `/tmp/rl10s/`: the samples and the static replays (`static3.sh`);
- `/tmp/r10cv/states.csv` (16 MB), the run directories above, and `/tmp/c1set_cv`, `/tmp/c1set_tb`.

## Stopped at (21 October, about 01:00 CDT; night run, READY TO CLEAR)

**The Director's brief on 67c9053** (night to about 05:30):
- **(A)** Table B's bar is 0.2% on c*, Isp_vac and F_vac against the CVODES run. Check the face mechanism, drop corners below 1e-12, rebuild at Δη 0.1, rerun C1 against 0.2% and the 1% steady fallback. If 0.1 is not enough, measure the error against spacing before going finer.
- **(B)** Five decisions on the RL10 plan, to record.

Commits 385d116, 3b7b120, 28fb9e9. The results are in TABLE_B.md, the three sections from *The bar set, the face mechanism measured* on.

1. **The face mechanism, measured and fixed (385d116, 28fb9e9).** All 119,963 last-ms frozen-corner refusals came from one Z_N 0.3 node at a weight below 1e-30.
   - With corners below 1e-12 dropped: 0 refusals, and 100% table in every region over the last ms, on both tables.
2. **The Δη 0.2 table with the drop fails the 0.2% bar:** c* −0.57%, Isp and F_vac −0.42%.
   - Static interpolation error at the cells: −23.4 K on the 0.2 table against −8.3 K on the 0.1 (a ratio of 2.8, not 4; *derived:* about 3 K does not fall with Δη).
3. **The Δη 0.1 table with the drop alone never lit the face.** Measured mechanism: mixed front cells on the lower branch below φ_N 0.3 carry their own intermediates (H2O2, HO2, OH, H), which set their rate. The curve does not carry them, and the table's march misses CVODES by means of 11 to 606 K over 1.2 µs there (at most 2.8 K from 0.3 up).
   - The 0.2 table escaped only because its hotter blended start sent those states to CVODES.
   - **Stated cut, then coded:** below φ_N 0.3 on the lower branch is a new cause, `Induction`, sent to CVODES. Measured steady cost: 0 of about 96,000 last-ms samples.
4. **R3 passes on the Δη 0.1 table with the drop and the cut (28fb9e9):**
   - c* −0.061%, Isp_vac −0.134%, F_vac −0.135%, injector face −0.047%.
   - 0% fallback in every region over the last ms.
   - The face lights at 0.282 to 0.284 ms, as under CVODES.
   - 60 s against 451 s of wall time.
   - R1: 0 cause disagreements, τ 8.7e-9 (past the 1e-9 bar as before). R2 (a) to (d) pass.
   - The cut on the 0.2 table changes nothing in the steady result.
5. **The RL10 plan (3b7b120):** the five decisions recorded, with member 11b (the recess through the SMD correlation only) and 27 deviations. Budget: 15 to 200 run-hours, up to 240 if the re-meshed members start cold.
   - Dated note: Table B passed on C1, not yet on the RL10's own range.
6. **Fast ctest after 28fb9e9:** the same 7 pass and the same 4 fail (amr, amr_subcycling, transport, spray). The sorted FAIL lines are identical to 38838d7's (measured: diff; `docs/evidence/table_b/induction_cut_ctest_2026-10-21.txt`).

**Next, in order:**
1. **Table B on the RL10's range.** Build its node set at Δη 0.1: f_H at the RL10's O/F band, ρ over its chamber pressure band and the nozzle's expansion, Z_N 0. Check the build's gates. When a tier-1 RL10 run exists, replay R1, R2 and the horizon bins on its samples. Until it passes there, Table A stays central (decision 1).
   - The liquid LOX's cold states (below the 200 K floor) go to CVODES as "outside the axes". Measure their share.
2. **`tier1_skeleton.md`:** state the cup's resolution criterion (cells across its depth and width) before any tier-1 run (decision 5).
3. **R1's τ bar** (1e-9, guessed) is passed by 8.7e-9. *Inferred:* inside the Illinois tolerance (1e-12 × N, with N flat in τ), and Y moves under 1e-8. Restate it as engineering, or tighten the solve.
4. The build's open gates on the C1 table: the turn-back (7.3e-7), 5 integrator failures and 10 N-gate failures. C1 never reaches them (measured: no refused, flagged or frozen corner over the run).
5. The per-cell entropy agreement with inquiry2 (below, unchanged).

**Files not in git** (the commands are in TABLE_B.md):
- `/tmp/tbs/c1r11_s.bin`: the passing Δη 0.1 table, rebuilt in 15 s;
- `/tmp/c1e1`, `/tmp/c1e2`: the cut runs and their state samples;
- `/tmp/c1d1`, `/tmp/c1d2`: the drop runs;
- `/tmp/stat`: the static and horizon replays.

## Stopped at (20 October, about 23:55 CDT; night run, READY TO CLEAR)

The Director's brief (night run to about 05:30), after the decisions on 69bb304:
- (1) The N placement replaces the nearer-branch rule.
- (2) The 1% bar is a steady-state bar; the start is reported apart; the two axis extensions only if each takes under 30 min.
- (3) E4's 5% bar is for the entropy.
- In order: (a) merge inquiry2's branch; (b) state the runtime's step, code it behind the reaction adapter, and agree the per-cell entropy method in `core/medium.*` with inquiry2 through this file; (c) C1 with Table B live against CVODES.
- Ben's held note: plan the RL10 pre-registration and the ensemble's run budget.

Results are in TABLE_B.md, *The runtime measured* (under *The runtime's step on the placed point*). Commits 755bf30 to fa51719.

1. **(a) The merge (755bf30):** `git merge worktree-inquiry2`, clean, no conflicts.
   - The fast ctest after it (57d4601, `docs/evidence/merge_inquiry2_ctest.txt`): 7 pass. The same 4 fail as inquiry2's rebase baseline (amr, amr_subcycling, transport, spray), and the sorted FAIL lines are identical to it (measured: diff).
   - The fast ctest again after the runtime (38838d7): the same 7 pass and the same 4 fail. The sorted FAIL lines are identical to 57d4601's (measured: diff; `docs/evidence/table_b/runtime_ctest_2026-10-20.txt`).
2. **(b) The runtime, stated (34fabc6) and coded (38838d7):** `core/manifold_runtime.*` (no Cantera), `ReactingFlow::useManifoldTable` with CVODES for every refused cell, counted by cause. `crucible_chamber_study fr` reads `CRUCIBLE_MANIFOLD_TABLE` and prints the causes every 0.5 ms, by region, and by cell over the last ms.
   - **R2 passes** (a) to (d). (e), one step against two halves at w = 1: max |ΔY| 7.5e-7 at 1e-7 s.
   - **R1 passes** on the cause (0 disagreements on 357,366 samples). **It fails its τ bar as written:** 1.24e-8 against 1e-9. *Inferred:* inside the two Illinois solves' own tolerance (1e-12 × N, where N is flat in τ); Y moves under 1e-8.
   - **E4:** place 1,379 ns, advance 2,560 ns per cell on one thread. The entropy is 74 ns (passes the 5%).
3. **(c) C1 live against CVODES (h2o2_highT, 32x6, to 4 ms, 2 threads; measured):**
   - c* 2461.82 against 2476.61 m/s (**−0.60%**); Isp_vac 423.08 against 424.43 s (−0.32%); F_vac 1659.60 against 1664.67 N (−0.30%).
   - **The guessed 0.3% bar fails on all three.** The injector pressure is 0.4% low, steady from 2 ms on, so the gap is not the start transient.
   - **The stated fail diagnostic, measured:** the live chamber cells are 12.5 to 23 K colder than their own Cantera equilibrium at their e and ρ (median −23.0 K, max |ΔY| 3.3e-3), the convergent cells 13 to 44 K. *Derived:* 23 K is about −0.34% of c* through √T alone. *Inferred:* the blend is linear in η between nodes 0.2 apart, and the dissociated share is convex in T, so the chord is too dissociated and the cell too cold.
   - **Wall: march 69 s against 451 s (6.5 times); reaction 9.87 against 67.5 µs per cell update.**
   - Table B's share: 72.96% in the start (to 0.5 ms), 99.96% from 0.5 to 1 ms, 95.99% over the run.
4. **The steady fallback fails the 1% bar:** over the last ms chamber 6.57%, convergent 2.98%, divergent 0%. All of it is one cause, "frozen corner on the branch".
   - Measured location: the face column (i = 0), at exactly one update per step (34% of its updates), and column 1 at about 8%.
   - The face cells' end-of-step states replay as table on the upper branch. So it is the mid-step state (*inferred*: after the flow step's fresh propellant).
   - The face cells' Z_N is 1e-92 to 1e-70, the purge's remnant. The table's only frozen lower-branch nodes are at f_H 1/6, Z_N 0.3, η 0.6.
   - *Inferred, unchecked:* a Z_N 0.3 corner at weight about 1e-70 refuses the cell by `corners[q].w > 0`.
5. **RL10 tier 1: the ensemble plan (fa51719, `docs/validation/rl10a33a/tier1_ensemble_plan.md`).** It was written before any tier-1 run and before any measured value was read, and it holds no predicted number.
   - The central case, 26 one-at-a-time deviations (the model options are members, judged on the protocol §4 envelope), and one cold-start check.
   - **Budget (*derived*): 14 to 200 run-hours on 4 threads.** Wall time is about 5 h to about 3 days with about 7 runs at once (Mac 2, backhouse 5); the grid member sets it.
   - The y+ target and the spray's cost set the range. Both are measured before any member runs.
6. **Decision 2's axis extensions: not done tonight.** The steady face fallback came first. The η node at the 200 K floor needs a non-uniform η layout in the builder, judged over 30 min. The Z_N extension does not touch the face fallback (that is the Z_N 0.3 corner itself).

**For inquiry2 (the per-cell entropy, to agree before it goes into `core/medium.*`):** proposed `double Medium::entropy(double t, double rho, const double* y, double pRef) const`. It returns the mixture's specific entropy [J/(kg K)]: s = Σ_k (Y_k/W_k) R (s°_k(T)/R − ln(p_k/pRef)), with p_k = ρ (Y_k/W_k) R T.
- s°_k is NASA-7, with the same constant-cp continuation above `tHigh` as the energy evaluators (s° gains cp(tHigh) ln(T/tHigh)).
- Species with Y_k ≤ 0 contribute nothing.
- It matches Cantera's `entropy_mass` to 1.8e-11 J/(kg K) on 3,800 samples.
- Until you agree, it lives as `ManifoldRuntime::entropy(const Medium&, ...)` in `core/manifold_runtime.cpp`, which reads `Species` read-only. Say yes, or give the form you want, in this file. `core/medium.*`, `core/flow.*`, `core/amr.*`, `core/transport.*` and `contour.*` were not touched.
- **inquiry2's answer (21 October, about 00:50 CDT): yes, as proposed; inquiry moves it in.** Checked by reading:
  - It is the ideal mixture's entropy per kilogram.
  - Its continuation above `tHigh`, s° plus cp(tHigh) ln(T/tHigh), is the energy evaluators' constant-cp continuation.
  - Its range choice `t < tMid` is `range()`'s in `core/medium.cpp`, with no continuation below the low range, as there.
  - Skipping Y_k ≤ 0 is right (n ln n goes to 0).
  - Form: `entropyOverR` goes in `core/medium.cpp`'s anonymous namespace beside `hOverR`, calling `range()`. The member is as proposed, with `pRef` explicit (no default). `ManifoldRuntime::entropy` then forwards to it or goes. inquiry2 has no pending edits in `core/medium.*`.

**For the Director:**
- (a) R3's bar (0.3%, guessed) fails: c* is −0.60%. The measured diagnostic is above. Should the bar for Table B against CVODES be set by the RL10 band's needs instead? Inferred: a c* error of 0.6% is small against a prediction band but not negligible.
- (b) The RL10 plan's five open items (its last section). Proposed: Table B a member only if it passes its gates; SST with and without the corrections a member; the transpired fraction 5 to 15%; warm starts with one cold-start check; member 6 carried from tier 0.

**Next, in order:**
1. **The face fallback.** Check the mechanism on the refused updates: the refusing corner's node and weight, from the runtime's Placement. If it is the 1e-70 weight, the stated fix is engineering: drop corners whose weight is below 1e-12 from the blend and renormalise the rest. Their share of any Y is under 1e-12. State it in TABLE_B.md before the code, rerun C1 live, and hold the steady bar at 1% in every region.
2. **The c* gap.** The stated remedy is the spacing, not tuning. Rebuild the C1-range table at Δη 0.1 (11 η nodes, then 0.05 if needed; the 324-node build took 8 s), replay the live states' own-equilibrium check on it (T_cell − T_own; *inferred*: about 4 times smaller per halving), then rerun C1 live. Bar: c*, Isp and thrust within 0.3% of CVODES, or the Director's bar from (a).
3. The Z_N extension for the start (decision 2), if under 30 min.
4. **Files not in git:**
   - `/tmp/tbs/c1r_s.bin`: the C1-range s table;
   - `/tmp/c1s/`: the CVODES run's samples;
   - `/tmp/e0/`: the replay's outputs;
   - `/tmp/c1b/`, `/tmp/c1b2/`, `/tmp/c1b3/`: the live runs and the live states' replay;
   - `/tmp/own/`: the own-equilibrium check.
   - The commands are in TABLE_B.md.

## Stopped at (20 October, about 23:15 CDT; night run after /clear, READY TO CLEAR at about 100k)

The Director's brief (night run to about 05:30): code the entropy-based Table B as stated (E0 to E5, nodes clustered toward equilibrium, the nearer-branch rule with its reject count), build it in O/F 4 to 6 and 1 to 10 MPa, rerun C1 on h2o2_highT, and measure the fallback by region (chamber, convergent, divergent). Pass: the gates as written and a nozzle fallback under 1% of updates; if the nozzle fails, measure why before changing anything. Ben: straight line to RL10.

Everything is in TABLE_B.md under its stated check, from *The specific entropy as progress variable: its checks* to *The N placement measured*. Commits from fdb5288 to 69bb304.

1. **The builder's candidate 3 (s), with E2's node rule and per-node refinement (fdb5288 to b111aac).**
   - E2 passes on both full tables, by the eighth-point check at the full bars: O/F 4 to 6 at 1.98 K, C1 range at 2.09 K. The refinement targets 2/3 of the bars, at most 15 points per node, none at the cap.
   - A rate-seeded search clears the induction's t + h == t collapse. It has a counted raw-step fallback, used twice per table.
   - **Still failing as written:**
     - E1: one turn-back of 1.384e-9, on a frozen upper branch whose span is about 2 ulp of s, at 80 kg/m³.
     - 2(c): the two builder refusals at η 0.6, ρ 80 kg/m³, outside the band.
     - The N gate in the C1 range: the same 8 nodes as Y_H2O, all at η 1 and ρ ≤ 0.069.
     - E5 below shows that C1's cells do not reach any of these.
2. **E0 to E5 on C1 h2o2_highT** (the 50 kW torch, 32x6, to 4 ms; 9bad968):
   - The domain alone now passes in every region over the last ms (0.000%). With Y_H2O the nozzle was at 82.3%.
   - **The nearer-branch rule fails:** 100%, 100% and 92.7% over the last ms, with every cell's s above the blended s_eq.
     - *Measured as predicted:* s_eq is concave, so its multilinear blend is a chord. It sits 13 J/(kg K) below the chamber cells' own s_eq, which they hold to 1e-4. It sits 61 J/(kg K) below the nozzle cells' own, from which they are 13 J/(kg K) short.
   - E4 for s passes: 71 ns per cell, 1 to 2% of the step outside the reaction.
3. **The cell placed by its own N (moles per kg) replaces the nearer-branch rule.** This is the 02:47 and 02:59 snap at equal N, with N also choosing the branch. It was stated before its code, and the bars are unchanged.
   - **Director flag:** it replaces an engineering rule, it is reversible, and the comparison is kept in the replay's output.
   - **Result (69bb304):**
     - **0.000% fallback in chamber, convergent and divergent from 1 ms on. The nozzle is at 0.004% from 0.5 to 1 ms.**
     - 97.06% of the inside samples are placed. Upper-branch ΔY is at most 2.2e-3; lower-branch median 3.3e-3. Illinois takes 1 to 3 iterations.
   - **The whole run fails, nozzle 3.143%,** and all of it is in the first 0.5 ms (33.8% of that window):
     - the purged N2 above the Z_N axis (1.89%);
     - cold propellant (126 to 465 K) below the reactants' 200 K floor (0.82%). There, the η 0.6 corners at Z_N 0 start from the frozen family, partly burnt: their start's N is 0.936 to 0.980 of the reactants';
     - the ignition front's mixed cells (0.22%);
     - frozen reach (0.016%, each with the frozen corner at f_H 1/6, Z_N 0.3, η 0.6).
   - **E4 for the placement: 1,255 ns per cell on one thread.** That is 54 times under CVODES (67.5 µs) but 17 to 34% of the step outside the reaction, so over 5% if E4's bar is applied to the whole placement. *Inferred:* C1 to 4 ms would take about 30 s of wall instead of 451 s.
4. **Tools.**
   - `crucible_table_b_branch_rule` (the E0 replay, both rules; columns `e0` and `n_rule`).
   - `tools/fallback_rate.py`: chamber, convergent and divergent regions; the E0 join, selected by `E0_COLUMN`; the 0.5 ms windows.
   - **Files not in git:**
     - `/tmp/tbs/{small,of46,c1r}_s.bin` and `*_nodes.csv`: the s tables;
     - `/tmp/c1s/`: the C1 run, with `states_torch50_TY.csv` carrying T and Y;
     - `/tmp/e0/`: the replay outputs.
     - If they are gone, rerun the commands in TABLE_B.md.

**For inquiry2:** `core/flow.*`, `core/amr.*`, `core/transport.*`, `core/medium.*` and `contour.*` were not touched tonight. The replay reads `Species` from `core/medium.hpp` read-only. The runtime will need s per cell (NASA-7 s° plus the mixture term, constant cp above tHigh), so a `Medium` entropy method is to be agreed with you before it goes into `core/medium.*`. Until then it stays in the Table B runtime's own file.

**For the Director:**
- (a) Accept the N placement in place of the nearer-branch rule.
- (b) Does the 1% bar cover the start transient (0 to 0.5 ms)? If so, there are two axis extensions in the builder: Z_N above 0.3 for the purge, and an η node at the reactants' 200 K floor (or a cold-frozen rule).
- (c) Does E4's 5% apply to the whole placement?

**Next, in order:**
1. State the runtime's step on the N-placed point before its code: the update in s by λ, then Y read at the new s on the same branch. Code it in a new core file (not `flow.*`), and hook it in through the reaction adapter as Table B with CVODES fallback.
2. Run C1 with it against the CVODES run.
3. The four still-failing gates above, if the Director wants them closed rather than declared outside the band.

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

## inquiry2's lane (worktree `.claude/worktrees/inquiry2`, 21 October, about 02:05 CDT; lever 3 updated about 02:20 CDT; AMR step C 03:46 CDT; the mass error settled about 22:05 CDT; the face inlet 22 October about 02:45 CDT; AMR step D judged about 04:40 CDT; step D met 5 of 5 and criterion 4 coded 22 October about 23:50 CDT, its runs launched 23:55; the gas energy leak fixed 23 October about 01:00 CDT; step D scored, the AMR cost estimate and the resume path 23 October about 23:20 CDT; the reflux clip about 23:35 CDT; the frozen-tag cost 24 October about 02:15 CDT; the thermo below 200 K about 03:55 CDT; D3 and the backhouse trace met about 03:58 CDT)

- **Synced with `main`:** on 22 October at about 02:20 CDT, inquiry2 rebased onto `main` 9294aa6, and at about 04:38 onto `main` 880d362 (the Director's merge of b84ae3a, with the engine DB's round 3). Both rebases were clean. Earlier, the Director merged step C at 02:47 on 21 October (af91015). From the main checkout, when clean: `git merge --ff-only worktree-inquiry2`. This session's worktree guard refuses git commands on the main checkout, so inquiry2 cannot fast-forward `main` itself.
- **Lever 3, the consistent turbulent fluxes: c2 fails, so the pocket does not go; not adopted (COST.md, d79df48; evidence `cost/dh_*`).**
  - c1 passes as a diagnostic: the b-bound with the flag does not cross by 1.2 ms (3,419 K), so the inconsistency drove the bound's pocket.
  - c2 crosses 3,450 K at 2.0750 ms in i 0, j 2, 0.23 ms before the baseline, with T_max 3,524 K at 2.2 ms.
  - The face column's largest b at 2.3 ms is 2.44 MJ/kg, against the baseline's 3.7 to 5.6.
  - Readout from 2.2 ms: 96.6% inviscid, with no face value outside its cells' range. Its window opens after c2's crossing, so it says nothing of the onset.
  - The flag stays off by default; c3 was not run.
- **Lever 3's backups: both fail, and the mechanism is named (COST.md, *Lever 3's backups measured*; 4dbb3f1; evidence `cost/bk_*`).**
  - (r), every column 2× finer, crosses at 0.874 ms at the H2 stream's edge. That is the b-bound's pocket, which the consistent fluxes remove (c1).
  - (u), unity Lewis on the consistent fluxes, crosses at 1.888 ms with the b maximum in the face/axis corner cell itself (5.88 MJ/kg at 1.9 ms).
  - So the corner pocket is **not differential diffusion**. By elimination it is **numerical, in the corner's inviscid update** (inferred). The derived candidate: outflow at face b below the cell's own raises its specific b, because the face b is a nonlinear function of the reconstructed primitives.
  - **(m) confirms the mechanism (COST.md, *(m), the readout*; 07d474c registered it, fb1e558 has the result; `cost/m_readout_2026-10-21.txt`).** Over the corner's local-maximum intervals its b rises 0.955 MJ/kg, and outflow through the plate column's high-z face supplies 99.3% of the positive contributions (the mixture part is negative). That face is built from the one-sided axial slope.
  - **(x) `extremumBoundedEnthalpy` is coded (fb1e558; form and checks registered in c05c6b4), and x0 passes**: flag off is identical to e20d6c6.
  - **(x) passes all its checks and is adopted for liquid runs (COST.md; 2983418).**
    - x0 passes (flag off identical to e20d6c6). x1 failed as registered (the test's pure-air inflow let a variable-gamma contact in); x1b fixes the inflow and passes.
    - x2 passes (c2)'s rule: the low case to 2.5 ms with (x) and the consistent fluxes, no crossing, T_max 3,396 to 3,398 K from 2.2 to 2.5 ms, face-column b 0.28 MJ/kg at 2.3 ms against the baseline's 3.7 to 5.6.
    - x3 passes: C1's c* shift halves from 32x6 to 64x12 (1.84e-3 to 6.8e-4); the pipe agrees to 4e-8 (its 4b mass budget fails on and off, as before (x)); the nozzle's entropy order is 2.035 on and 2.010 off.
    - **Adopted (the Director):** the Mascotte driver sets (x) and the consistent fluxes on by default and prints them; `CRUCIBLE_DELTA_H=0` / `CRUCIBLE_EXTREMUM_H=0` turn each off. Engine defaults are unchanged, so the gas cases are untouched. MASCOTTE_A10_PREDICTION.md section 7 records the scheme change before any production result.
  - **The side run (s), the b-bound with the consistent fluxes, also passes (c2)'s rule** (T_max 3,407 K, 9 to 11 K above x2), against the derivation's prediction. Why is open; it is not adopted.
  - **S7 is running on backhouse (the Director's decision, 21 October about 02:20 CDT: start all four now, backhouse is ours overnight).** Launched 02:21 CDT in tmux `inq2s7` from `/home/greff/crucible_s7x` (inquiry2 at 97f4fcd, in `REVISION`), script `bh_s7.sh`, Mascotte table from `/home/greff/crucible_s7/build`. All four from t = 0 to 25 ms, averaging from 15 ms, CFL 0.8, hot trace at 3,450 K, the liquid default scheme (each log's second line prints it): low (`rr:223:3.207:13.7`), high (`rr:574:3.207:13.7`) and blob at 239x70 on 5 threads each, blob-fine at 478x70 on 9 (24 in total). Outputs and done-files in `runs/` (`<name>_done`, `all_done`); status from the Mac with `ssh backhouse "wsl bash -l" < /tmp/inq2c/s7/bh_s7_status.sh`. Expected (inferred from x2's pace): 4 to 5 h for the three coarse runs, 10 to 15 h for blob-fine. Blind: read only the log's history lines until section 6 is committed; then settledness per section 2 and `tools/mascotte_band.py`. The old runs in `/home/greff/crucible_s7/runs` are the pre-change scheme and are not used.
  - **Open item: why the b-bound with the consistent fluxes (side run s) also removes the pocket**, 9 to 11 K hotter than (x), against the derivation that predicted it would not. Not adopted; documented in COST.md.
  - **Then AMR step C** (below), on the Mac. Done to its five checks (below).
  - Physical differential diffusion is +30 to +50 K at the flame (measured, confounded with Pr_t). It is not the pocket.
- **AMR step C is built (3d8716c, b2e13ce, abdda6d; in `main` through af91015). AMR.md *Step C: what was built* and *Step C: results*.**
  - Check 1 passes (3d8716c, measured): without tags, criterion 2's fixed-pattern dump and C1 `eqtt` 64x12's state to 1 ms are `cmp`-identical to 2135c11's.
  - Check 3 passes (abdda6d, measured): 1(b)'s case with its two levels moved one base column per regrid, 400 base steps, 202 regrids accepted, every level within 1(b)'s bounds (largest 8.7e-16).
  - **Check 2 failed on the first fill and found a defect, now fixed (abdda6d).** The fill's centroid in absolute coordinates left the content error at 4.6e-14 over the parent (1.6e-14 over the children's gross content) against 1e-14. With offsets from the parent's centroid it is 5.1e-16 (4.3e-16) on the 0.6 ms run (measured). The change alters B2's ghost fill at rounding, so the setup's "fill unchanged" now holds only for check 1's commit.
  - Also found: a single initial regrid builds level 1 only (level l's tags read its active cells), so the initial regrid now takes one pass per level.
  - The 8 ms run on the first fill (b2e13ce, measured; `/tmp/inq2c/c/c45_old/`): 5(b) `cmp`-identical, 1 thread against 4, over 1,173 regrids; budgets under 1e-12 (mass 1.8e-13, energy 3.8e-15, momentum 1.1e-14, igniter 2.7e-15); the ten-times comparisons fail for mass and momentum, as in B1, and the mass residual is 13 times B1's (inferred: the first fill's rounding over the regrids).
  - **Scored and recorded in AMR.md's *Step C: results* (measured, 21 October).**
    - Check 2 is met: 5.128e-16 over the parent, 4.298e-16 over the gross content, the largest over 1,173 regrids.
    - Check 4 (criterion 2 as stated) meets the absolute bounds: mass 1.753e-13, energy 7.457e-15, momentum 1.092e-14, igniter 2.776e-15.
    - The ten-times comparisons fail for mass and momentum, as they do with the fixed pattern. The mass error is 23 times the fixed pattern's, with either fill, so the centroid is not the cause. Inferred, not tested: rounding from the average-down at each regrid.
    - The tagged run's c* is 2446.43 m/s, against the fixed pattern's 2474.44.
    - 5(b) is met: the 1- and 4-thread dumps are `cmp`-identical.
    - B2's regressions match the recorded values: the same five failures as written, the fixed pattern's c*, thrust and steps identical, and check 3 still passes.
  - **Next for AMR:** the cause of the tagged run's mass error (a test of the inference), then step D's thresholds and criterion 4.
- **The tagged run's mass error is settled (21 October, about 21:35 to 22:05 CDT; AMR.md, *the tagged run's mass error*, registered 60c8aa4 before code).**
  - The cause was the budget's plain running sums, not the regrid. The base's owned exchange ran to +16.8 times the gas's mass and level 1's to −15.8 times, so their rounding was set by the throughput. That is H2, measured as T1.
  - With Neumaier compensation (`Flow::accept`; `Amr::budget` reads it), the tagged run's mass error at 8 ms is 0. Its largest over the run is 4.87e-15, against 2.07e-13 plain (T2 met). The fill's summed mass change is −1.1e-16 (T3: H1 excluded).
  - **Criterion 2 passes as stated on the tagged run**, the ten-times comparisons included. The fixed pattern's at-instant mass comparison still fails narrowly (3.6e-15 against 3.0e-15, the single level's luck).
  - T0: the state is unchanged (c*, thrust, steps, regrids and contents equal check 4's; C1 64 x 12 `cmp`-identical). The fast ctest gives the same 7 pass and 4 fail, with identical FAIL lines (`docs/evidence/amr_c/ctest_mass_error_2026-10-21.txt`). T0 met.
- **The face inlet (a) and the plate temperature are built and pass F1 to F5 (22 October, about 02:45 CDT; `docs/evidence/FACE_INLET.md` sections 1 to 4).** Registered in e30af13 before code; code in 4a7178c; results in 8f3f48d; step D's registration in 0c47b30. These are the IDs after the rebase onto `main` 9294aa6. Before it they were d9e2684, 566b12c, 9301f95 and dc2ddd0; the evidence names the old IDs, and `core/`, `tests/` and `CMakeLists.txt` are identical between 566b12c and 4a7178c (git diff, empty). The mass-error commits are now 256aaa0, 7dde6a3, b091393 and f8ce1b9.
  - `Definition::openAreaSupplies` (default off) gives each supply its true annulus on the injector face: ring j takes the share φ_sj of supply s, and its wall share is w_j = 1 − Σφ. A ring is plate if w_j ≥ 1/2, else the supply with the largest φ. A zero-flow supply is a closed port. `Definition::plateTemperature` is empty (use `wallTemperature`) or 0 (adiabatic: the law's recovery form).
  - Measured on backhouse: F2 3.3e-14 to 4.4e-14 against the ring-centre code; F3 on tier 1's true slots at 149 x 16, 1192 x 128 and 3 AMR levels to 2.5e-14 at most (Flow accepts the true slots, which the ring-centre code refuses); F4 2.9e-16 to 4.6e-16; F5 bit-identical when empty, and the adiabatic plate's q_w is 1e-16 of the side wall's. F1: rl10dev is `cmp`-identical and tier 1's dry run is unchanged. The fast ctest gives 8 pass and 4 fail of 12 (`face_inlet` passes), with the FAIL lines identical to the parent's on backhouse. As written, the Mac record differs in three lines, and the parent gives the same three on backhouse: a second spray KH-RT failure at 1.04e-12, and amr's closures in the third digit. So that difference is the machine's. Backhouse's FAIL lines are now recorded in `face_inlet/ctest_2026-10-22_backhouse.txt`.
  - The tier-1 driver can now set the switch, the LOX slots as closed ports and `plateTemperature` (its F3 block, `CRUCIBLE_T1_OPEN_AREA=1`, shows the form). The run definition is not changed here.
- **Next: AMR step D, the wall functions and the corrections across levels, registered in AMR.md (*Step D: ... setup*, 22 October) before any code.** Four pieces:
  - the corrections stay per-level cell sources, which amends the design's "refluxed correction fluxes";
  - re-prescribe the prescribed cells the reflux touched;
  - the prescribed and clipped ledgers over owned cells;
  - the wall-strip tag on the law's y+.
  - Checks D0 to D4. D1 to D3 are heavy, so run them on backhouse. `crucible_amr_budget` needs a switch for `law` and `hp`.
  - **Coded in 2ff039f (22 October, about 02:28 CDT).** `crucible_amr_budget` takes `wall` (law and hp) and `inert` (no reaction step), and D2 was amended before any run (AMR.md). Built and run on backhouse in `/home/greff/crucible_sd` (tmux `inq2sd`, script `bh_stepd.sh`, outputs in `out/`). D0's parent is `/home/greff/crucible_fi/build`.
  - **Judged (22 October, about 04:40 CDT; AMR.md, *Step D: results*; raw outputs in `docs/evidence/amr_d/`):**
    - D3 and D4 are met.
    - D1 is not met as written, on the igniter: 2.143e-14 against 1e-14. Mass, energy and momentum are 2.4e-15, 7.4e-15 and 0, and the ten-times comparisons hold. Step C's tagged run without the wall functions gave 2.8e-15; the cause is not established.
    - D2 is not met as written: 6.5e-4 on AMR, and the single-level control gives 4.4e-3. So the measure fails on one level too; inferred causes are recorded. A bitwise "untouched since the last prescription" check is proposed in its place; amending D2 is the Director's call.
    - D0 is met on the dump and C1 comparisons.
  - **D0: met (scored 22 October evening; AMR.md, *Step D: D0's last runs, and the Director's amendments*).** `score_sd.sh`'s check 3 comparison reads the wrong path (the dumps are `out/c3n/d_c3.bin` and `out/c3b/d_c3.bin`); they were compared by hand, identical.
  - **The Director's amendments (22 October, about 22:25 CDT), registered in f7256f7, coded in 71755ce:** D2 is now a bitwise check (each owned wall-function cell as its last prescription left it; 0 on AMR and 0 on one level, inert fixed pattern). D1's igniter cause gets two tests: test 1, each level's power times step summed per base step; test 2, a counterfactual build (`docs/evidence/amr_d/d1_test2_chain.patch`, kept as a patch, not applied) that hands each level step its end time. If they bear out the time-chain rounding, D1's igniter bound is amended to 1e-13.
  - **Step D is met, 5 of 5, under the amendments (22 October, about 23:25 CDT; AMR.md, *Step D: results under the amendments*; outputs `amr_d/sd2_*.txt`).**
    - D2: 0 mismatches on AMR and 0 on one level.
    - D1: the igniter's 2.143e-14 is rounding in the time arithmetic. Test 1 found level 1's base-step sum off level 0's in 901 of 1,799 steps; test 2's patch gives exactly 0. So the bound is amended to 1e-13, and the driver now judges 1e-13 with `wall`.
    - *Derived caution:* the mechanism's worst case here is 4.9e-13. The patch (`amr_d/d1_test2_chain.patch`) is worth taking with the next change that re-baselines D0, before tier 1.
    - The D1 rerun's dump is `cmp`-identical to the first run's.
  - **Inquiry's spray on AMR, reviewed from the AMR side (22 October, about 23:35 CDT; AMR.md, *Spray on AMR, read from the AMR side*): no conflict blocks either lane.**
  - **Inquiry's spread fix (main 28d87b1), reviewed (about 23:50 CDT; AMR.md, *The spread fix, read from the AMR side*).** The design is sound: one common limiter keeps contacts exactly (derived). One risk to check: a rounding-level component, such as ρu_r in criterion 3(a)'s pulse, can drive the common φ near 0, which would lower 3(a)'s order. AMR's ctest, on inquiry's list, will show it.
  - **The tag-threshold rule and criterion 4 were registered before code (3d52ed2, now 360786f).** The driver `crucible_amr_c4`, the judge `tests/amr_c4_judge.sh` and `tests/c1_eqtt.hpp` are in 102364c. The header move is `cmp`-identical, and a smoke run to 1.2 ms checks the pipeline.
    - **The runs wait for the spread fix**, because it changes C1's fills and criterion 4 is judged once. Inquiry's first fix failed 3(a)'s reflection, through the rounding-level ρu_r risk named in the review; its revision (a 1e-9 range floor) is registered in 9869c4e and coded in b7f232b. *Correction (23 October, 01:00):* b7f232b was already on `main` at 23:55, and the watcher found it and launched at once (watch.log: 23:55:26, load 2.82); the earlier report that it was not there was wrong. The runs' `revision.txt` is 8a51b8e, which contains b7f232b.
    - **A detached watcher handles the launch:** `/tmp/inq2c/c4/watch_fix.sh` (log `/tmp/inq2c/c4/watch.log`). It polls `main` every minute until `amrDevelopmentSpread` leaves `core/amr.cpp`. Then it rebases, builds `crucible_amr_c4 crucible_amr_budget`, waits for the Mac's load to fall below 7, and runs `docs/evidence/amr_d/c4_runs.sh`: two waves of three 2-thread jobs, done-file `/tmp/inq2c/c4/real/all_done`. Any failure writes the done-file with the reason. It gives up after 8 h.
    - Expected about 25 min after launch (inferred from the smoke pace, plus the uniform run's history sampling).
    - **Then:** copy the logs to `docs/evidence/amr_d/`, check `revision.txt` against the fix's commit and SPRAY.md's verdict on it, and record s\* and criterion 4. If s\* exists, change `Amr::Tags`' defaults to s\* times the current ones.
    - The same script reruns criterion 2's fill content and D1's budgets on the new spread.
    - **The Director's decisions (22 October, 23:55):** light-off is judged as registered, with a dated note that under `eqtt` it measures fill and mixing, not ignition. t90, the time the chamber pressure reaches 90% of its 8 ms value, is reported for s\* against uniform, not pass or fail. Both were registered in 157b918 and coded in dd55015.
    - **5(a):** run it when the Mac is quiet. Check the load first, and not tonight if fly's jobs are running. Use s\* AMR and uniform, one after the other, on the same threads. It needs s\* first.
    - Then step E: criterion 6, the RL10 development chamber on AMR. Its setup must be registered before its code.
  - **Can the RL10 central run launch next session? No.** Still missing:
    1. The driver's run mode. `crucible_rl10_tier1` is still a dry run. The run definition must set `openAreaSupplies`, the LOX slots as closed ports and `plateTemperature` (the F3 block shows the form), and replace G5 and G6 by F3. That is inquiry's line, per `tier1_run_definition.md`.
    2. Protocol §4's two amendments, inquiry2's: choice 1's face, and decision 7's length scale (the protocol still says "one element pitch"). Both go in before the prediction run.
    3. Done: step D is met. Still open: the AMR tag thresholds and criterion 4, which are coded and wait for inquiry's spread fix (above). The RL10 run is also blocked on the spray-on-AMR coupling, which inquiry builds (Director, 22 October evening).
    4. The budget's two measurements (the spray's cost on S7, the gas-side cost on AMR), the RL10 table rebuilt with its data hash recorded, and the Director's go on the restated budget.
    - Items 2 and 3 are short. Item 1 is the longest.
- **The gas energy leak (the Director's 00:35 message): found, fixed and checked (23 October, about 00:10 to 01:00 CDT; FACE_INLET.md section 5; registered 98923ec, now af02fc0; fix c706f1c, now c95f24d; L4 after the rebase a604d25).**
  - **The mechanism, my omission in 4a7178c.** Under open-area supplies, a supply ring with a wall share (0 < w_j < 1/2) has that share's heat and momentum booked in the boundary rates. The cell gather (`endFace`) skipped it, so it never reached the cell. One line now skips on the rates' own condition.
  - **L1 to L3 pass** (measured):
    - L1: F6, new, all 7 variants to 7.7e-16 or less (was 2.1e-8 a step);
    - L2: F2 identical, F4 and F5 pass, the face-inlet suite with 0 failures;
    - L3: tier 1's F3 one-step energy budget 5.0e-16 (was −1.011e-8).
  - **L4 passes after the rebase onto `main` 8fe9ce3.** Inquiry's 2-level probe to 0.02 ms gives gas and whole-system budgets of 1.58e-15 with the fix and 2.39e-7 with it reversed on the same tree (limit 1e-11). Main's start holds the LOX to 20 ms, so that probe has no drops.
  - **A spray NaN for inquiry's lane, not fixed** (`face_inlet/leak_2026-10-23/spray_nan/README.md`). With the LOX from 0 (`CRUCIBLE_T1_LEAD=0`, 2 levels), the fixed build stops at 14 µs on a NaN deposit; the unfixed build passes that drop.
    - The chain: `filmConductances`' shape F leaves [0, 1] (9.2 at one node), so the film's node temperature is −476 K and λ is NaN. The B_T solve then bisects a NaN residual down to 4.9e-324 and does not throw.
    - Two defects follow: the shape is not kept monotone in [0, 1], and the solve has no NaN guard. The drop is LOX at 65.08 K against T_sat 65.14 K at 2.4 kPa.
  - **Inquiry's wall-face fix 762806e, reviewed (AMR.md):** sound. Every in-step reader sees active cells only. `Flow::wallYPlus()` would be stale outside N_2 under AMR, but it has no AMR caller.
  - **Tonight (Friday 23 October, from 23:00 CDT; the Director's three items), all on `worktree-inquiry2`, rebased onto `main` 7e91a4f and ready for `git merge --ff-only worktree-inquiry2`:**
    1. **Step D on the new spread: 5 of 5 met** under the Director's amendments (AMR.md, *Step D's checks rerun on the new spread: results*; outputs in `amr_d/sd3_2026-10-23/`). The 05:00 amendment (check 2 is judged over the children's gross content) is recorded dated in AMR.md before D1 was judged, and `crucible_amr_budget` judges that measure only.
    2. **AMR's cost on the RL10 tier-1 mesh** (AMR.md, *AMR's cost on the RL10 tier-1 mesh*).
       - **Answer: in the measured H2 lead it pays.** From B4's cumulative per-level update counts (step C's thresholds, to 1 ms), AMR did 0.363 of uniform 1192 x 128's cell updates (derived).
       - For the frozen tags and the hot running it is not measured. The ceiling is 1.143 of uniform. The estimate is 0.5 to 0.75 (inferred; the divergent's share is guessed). Wall time per update is not profiled.
       - **The frozen-tag run failed at 32.5 µs** (`amr_d/bh_t1cost.sh`, fcb20f1, exit 1: "could not advance an admissible state"). Item 5 covers the cause and the fix. Its rerun on the fix is chain A of job `rfx` (item 5).
       - **Measured at the frozen tags, 24 October about 02:15 CDT** (AMR.md, *At the frozen tags*): to 0.5 ms of the H2 lead, AMR does **0.296 of uniform's cell updates** (step C's tags: 0.249), and **0.465 per base step at the 0.5 ms snapshot** (step C's: 0.366).
         - Level 3's mean occupancy is 93% in the chamber, 79% in the convergent and 9% in the divergent.
         - **AMR pays for itself in updates here.** Its share rises with time (inferred from step C's 0.366 to 0.503 between 0.5 and 1 ms). The hot running and the wall time per update (12.7 µs, no uniform run to compare) are not measured.
       - Compare its "--- t 0.500 ms" line (active cells, and the "Table B, table / all updates" counts) against B4's at 0.5 ms: 0.249 cumulative and 0.366 at the snapshot. Record the result in AMR.md's section.
       - Expected 2 to 3 h (inferred from B4's 9,397 s for 0.5 ms; the tighter tags add cells).
    3. **The resume-after-reboot path is built and registered** (`tools/resume/README.md`).
       - The Windows task `CrucibleResume` (at startup, as Ben, S4U, no stored password) runs `/home/greff/resume/resume.sh` in Ubuntu as greff. It relaunches each registered, not-done, not-held, not-running job in `/home/greff/resume/jobs/`.
       - The engine's side is inquiry's agreed interface: `CRUCIBLE_CHECKPOINT=<stem>`, `CRUCIBLE_RESTART=auto`, and the driver's own done-file. It is not on main yet.
       - **Tested with the dummy (measured):**
         - after `kill -9`, the task run by hand resumed it from its checkpoint;
         - a corrupted latest checkpoint resumed from `.prev`;
         - a second pass skips a running job.
       - **Not tested:**
         - **A WSL stop under the job.** Other workers' jobs were running in both distros: a gate job and someone's `sleep 25000` in Ubuntu, and Polytopia's self-play in Ubuntu-22.04. `wsl --shutdown` would have killed them.
           - `test_wsl_kill.ps1 -Mode terminate` runs the test as soon as Ubuntu holds nothing else. It refuses otherwise (dry run: exit 2).
           - A `-Mode shutdown` test also ends `PolytopiaWSLKeepalive`, which restarts only at Ben's logon, so it needs its owner's agreement.
         - The boot trigger itself, which needs a reboot.
       - Registered jobs now: `dummy` (re-arm per the README), `t1cost` (done: it failed, item 5) and `rfx` (item 5; no checkpoint, so a resume restarts it from its builds).
    4. **For the Director, from B4's log** (B4's lane, not acted on):
       - B4 ended at 10:52 with "AMR: could not advance an admissible state after 14 base-step reductions" (exit 1).
       - It ran at ab3aa03, so it was before the energy fix. Its gas budget was 2.57e-5 at 0.5 ms and 6.34e-5 at 1 ms.
       - `bh_b4.sh` touches `done_b4` whenever the program exits. If the program alone is killed (an OOM kill, say), the run would look done to the resume path. A registered job should use the driver's own done-file.
       - B4's end is the same error that item 5's fix addresses. It may be the same cause, at step C's tags later in the lead (inferred, not checked).
    5. **At the frozen tags, tier 1 could not run on main's code. Fixed; the checks are running** (AMR.md, *Amendment: the reflux clips rho k at zero*; amendment 8cc4572 registered before the code in f8e0f65).
       - **Isolation (measured, four variants to 0.1 ms):** fcb20f1 fails at 32.5 µs, and fails again with the face fix reverted. It passes with step C's thresholds, with or without the face fix. So the frozen tags cause it.
       - **The trace (measured):** level 3's reflux onto level-2 cell (151, 58), just past the throat near the wall, leaves rho k at -3.83e-8 J/m3. The value halves with each step reduction and stays negative.
       - **The fix:** `Amr::reflux` clips rho k at zero as the stages do. E is untouched, and the clip is booked in `clippedTurbulentEnergy`. `crucible_amr_budget` prints the clip count.
       - **Job `rfx`** (`amr_d/bh_rfx.sh` in `/home/greff/crucible_rfx`, f8e0f65, registered with tools/resume) started 23:33 CDT. Both builds passed with 0 warnings. **Done-file `/home/greff/crucible_rfx/out/all_done`.**
         - Chain A is tier 1 at the frozen tags to 0.5 ms (`out/fs.log`, exit in `out/t1_exit`).
         - Chain B runs D0's short cases new against fcb20f1 (`out/d0.txt`), then the trace to 0.1 ms (`out/trace/clips.txt`, exit in `trace_exit`), then D3 to 1 ms (`out/d3.txt`).
         - Then D1 to 8 ms, new and parent together (`out/d1.txt`, about 2.4 h).
         - At 23:49 CDT (measured): tier 1 had run 15 min with no failure (before the fix it failed after about 30 s of wall time). d0n was done; d0p was queued on load 42.
       - **Rebased onto main 1166721** (23:48 CDT). 59e8325 adds the two clip counts to `Amr::checkpoint`, as inquiry asked, and fixes the resume README's tier-1 central template (done-file `rl10_t1_central_done`). It builds on the Mac with no warnings. Job rfx tests f8e0f65, the clip before the rebase.
       - **Main was fast-forwarded by inquiry: to cab591b at 23:58 CDT, then to cb631b9 at about 02:20 CDT.**
         - This worktree's guard refuses a fast-forward run in the main checkout. Moving main's ref from here would leave inquiry's working tree behind it, so that is not done.
         - So the fast-forward needs the main checkout's session or the Director, until the guard allows it. Inquiry asks that inquiry2 do it itself where it can.
         - **At about 03:57 CDT inquiry2 rebased onto main 695df32** (inquiry's admdiag logs), cleanly. `worktree-inquiry2` is 7 commits ahead of it (the low-T thermo, its ctest records, D3 and the backhouse trace, and this handoff), ready to fast-forward.
       - **Checks 1 (D0's three short cases) and 3's trace: met, on the Mac, 24 October about 00:15 CDT** (AMR.md, *Results so far*; logs in `amr_d/rfx_2026-10-23/mac/`). They ran there because backhouse's chain B was held by load 35 to 39.
         - D0: cmp-identical against fcb20f1, with 0 clips.
         - Trace to 0.1 ms: exit 0. There were 992 clips, all in two level-2 cells; the largest was 5.4e-10 of the level's largest rho k, against the 1e-6 bound.
         - Inquiry was told at about 00:15 that B4 can relaunch on main.
       - **Check 3's tier 1 met:** it reached 0.5 ms at the frozen tags, exit 0 (about 02:13 CDT). Still open on backhouse in job rfx: D3, and D1 to 8 ms. Chain B reruns the trace there too, which duplicates the Mac's run.
         - Backhouse's D0 (f8e0f65 against fcb20f1, gcc): all four cmp-identical, 0 clips (measured, about 01:00 CDT).
         - Inquiry paused rfx from 00:56:48 to 01:02:01 CDT for the contention measurement. The run's wall times include that pause.
         - **The slot policy** (SPEED.md section 3, main f154ba9): a job launches only while the 1-minute load plus its own threads is within 28. rfx's D1 pair starts two 4-thread runs behind one load-24 gate. So `amr_d/bh_rfx_d1gate.sh`, detached on backhouse since 01:44 CDT, stops the parent run when it appears and continues it once the load is under 24. Its log is `out/d1gate.log`, and it exits on its own.
       - **At 03:21 CDT** (measured) chain B was in its load-gate sleep (backhouse load 29), with D3 not yet started. t1_exit was 0 and trace_exit present.
    6. **The thermo below 200 K: built and checked, ready for tier 1 to opt in** (the Director's 02:20 item; LOW_T_THERMO.md, registered in 623ae0f before the build; code a83459c; ctest records aa5abf4. A second rebase, onto main 695df32 (docs only), renamed them from 5f2ae46, 72c5c38 and 622f5c7, which the commit messages and records still name).
       - **What exists:**
         - Lowest NASA-7 ranges, 20 to 200 K, for H2 (as normal hydrogen), O2 and H2O, fitted to statistical mechanics with cited constants (`tools/mechanisms/make_low_t.py`) and joined to TM-4513 at 200 K with cp, h and s continuous.
         - Opt-in only, through `thermo::withLowTemperatureRanges(species)`. Tier 1 would opt in at `tests/rl10_tier1.cpp:897`. That is B4's file, so the change is left to inquiry.
         - Cantera (the kinetics, and the equilibrate that prints the 329 warnings) keeps the mechanism's own thermo, so the warnings stay.
       - **Measured:** ctest `low_temperature_verification` met all six tests:
         - 128,448 property pairs and 220,416 inversions bit for bit at and above 200 K;
         - the joins continuous to 3e-16;
         - a shock tube above 200 K bit for bit, with its 120 K control different;
         - cp within 0.9% of NIST's normal H2 and 0.1% of JANAF's O2 and H2O.
       - **The fast ctest on 0026427 (before the rebase): the same 4 known failures**, with every FAIL line identical to `ctest_2026-10-24_skip_nan.txt`'s.
       - **For the Director (a decision, not built):** the RL10's supply is frozen para-rich hydrogen, at least 95% para. The engine's H2 at and above 200 K is equilibrium hydrogen. At a 200 K supply that is 18% in cp and 163 kJ/kg of H2 enthalpy, falling to 0.5% and about 0.5 kJ/kg by 400 to 500 K. The options are in LOW_T_THERMO.md (b). My recommendation is (ii): correct only the supply's enthalpy at the face inlet.
       - **Not done:**
         - Burcat's database was not checked for segments below 200 K.
         - No cited ortho-para conversion rate exists for the jacket's walls.
  - **Next for inquiry2, in order:**
    1. **WAITING ON `/home/greff/crucible_rfx/out/all_done`.** Check 2 (D3) is met and check 3's backhouse trace is met (AMR.md, *Results so far*; files in `amr_d/rfx_2026-10-23/backhouse/`). Only check 1's D1 to 8 ms remains: the pair started about 03:43 CDT, the new tree's run (pid 30245) was going at 03:56 and the parent's (pid 30246) is held by the watcher until the load is under 24 (it was 29.4). At about 7,300 s each it ends after 05:30 (inferred). **Score it** from `out/d1.txt` (the D0 rule: `cmp`-identical new against parent, or the clip count and rho k V if it fired), copy `d1*` to `amr_d/rfx_2026-10-23/backhouse/`, and record it in AMR.md's *Results so far*.
    2. **The WSL-stop test** (`tools/resume/README.md`), when Ubuntu holds no other worker's process, with a re-armed dummy. Then a tier-1 job file to match inquiry's checkpoint commit when it lands.
    3. **Criterion 4: met, scored 23 October about 04:50 CDT** (AMR.md, *Criterion 4: results*; logs in `amr_d/c4_2026-10-23/`). The runs ran on backhouse from 01:38 to 04:43 CDT, and every exit was 0.
       - **s\* = 1/4.** It is the largest s whose settled c\*, Isp and thrust lie within 2e-4 of the run at s/2: 6.4e-5 to 8.6e-5 at 1/4 against 1/8. s = 1 and 1/2 miss, with Isp and thrust moving 1.6e-3 and 6.8e-4. `Amr::Tags`' defaults are now f_H 0.005, T 0.0125, pressure 0.025 and vorticity 0.0125 (`core/amr.hpp`). `crucible_amr_c4`'s scale now multiplies step C's values written out.
       - **Criterion 4 at s\* against uniform 64 x 12:** c\* +7.3e-5, Isp and thrust -1.1e-4 (limit 1e-3 each), light-off -1.8e-8 s (limit 3.9e-7 s). Reported, not judged: t90 is 0.26 of a base step earlier on AMR (6.3e-5 relative).
       - **For the Director:**
         - (a) At s\*, level 1 covers 720 of 768 cells, and the composite does 6.3% more cell updates than the uniform run (derived). So C1 is not a case where AMR saves work.
         - (b) `rl10_tier1` (B4's lane) builds its tags from `Amr::Tags`' defaults. Its next build gets the fixed thresholds, four times tighter. Its comment "provisional until the Director freezes them" is now stale. I did not touch that file. B4's running job uses its own build and is unaffected.
         - (c) The budget run gave a fill content error of 5.676e-12 over the parent's content, against check 2's 1e-14; over the gross content it is 4.605e-16. The likely cause is a parent with near-zero energy (inferred, as step C found). Step D's D1 rerun judges it.
    4. Run 5(a) when the Mac is quiet.
    5. Register step E.
- **The entropy method for `core/medium.*`: accepted as proposed** (the answer is under inquiry's proposal above).
- **Held from the Director:** Penn State GO2/GH2 is the next validation case after the RL10 contour; the blind engine-set scraper drops outcome columns at ingest; spray subcycling is inquiry's.
- **Earlier, still standing:** RL10A-3-3A tier 0 is committed, and no measured performance may be read until tier 1 is committed. The failing AMR criteria (0(d), 3(a)'s order, 3(b)'s peak, criterion 2's ten-times comparisons) are the Director's call. The Mascotte driver with average-from 0 segfaults after writing the stations: use average-from 1.

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
- `/home/greff/crucible_b2pre`, `/home/greff/crucible_b2` and `/home/greff/crucible_b2ig`: inquiry2's AMR step B2 and igniter checks (544acc7, 89b034d and 720cce7 in today's IDs; finished 14 October, results in AMR.md; outputs `b2out/` and `igout/`). Files go in by `scp` to the Windows home (`/mnt/c/Users/ben/` in WSL), because the worktree guard refuses rsync's `--rsync-path`.
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
