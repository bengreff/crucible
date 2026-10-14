# Session handoff

14 October 2026 · branch `main` (one line). The spray worker (step 5) on `main`; inquiry2 works in `.claude/worktrees/inquiry2` (RL10 contour, AMR, table generalisation, and transport cost in `core/medium.*`, `core/transport.cpp`, `core/flow.*`, `contour.*`: not touched from here). The 6 October state is in this file's git history (`git show 15d6777:docs/SESSION_HANDOFF.md`).

## Stopped at (14 October, about 00:35 CDT)

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
- **Wall functions: 2(a)'s resolved half restated like for like** (898f631, after the result, flagged): engine against tool on the same 64 rings, passes (+0.27%, −0.51%); the shared grid error is a stated limit.
- **fast ctest** (13 October): 7 of 9 pass; the failures are the known ones (`transport_verification` radial order 1.218 and 1.444, `spray_verification` S0(b)). The slow pipe test has still not been rerun.

## Exact next steps (in order)

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
