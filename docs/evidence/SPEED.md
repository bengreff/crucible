# Speed: where tier 1's wall time goes and what sets its step (24 October 2026)

Asked by the Director at 23:50 CDT on 23 October, after the cost estimate of 7 to 11 weeks for the central case (`docs/validation/rl10a33a/tier1_ensemble_plan.md`, *Per run*). Four items: (1) the profile and the binding constraint on the finest level, measured; (2) the levers, with the speed-up each would give and its cost to build, reported, not built; (3) thread scaling on backhouse; (4) whether 60 ms of flow is physically required. Labels: measured, derived (stated arithmetic), inferred, guessed.

**Added for the Director's 00:50 note** (principle 1 kept as written: no steady solve, no pseudo-time that alters the transient):
1. The implicit time-accurate step: the physical step the dynamics need, against today's explicit step, and the cost per implicit step. Section 5a, lever I.
2. Backhouse contention, measured with our jobs paused, pinned and unpinned, and a slot policy, revised for the Director's 01:40 decisions (the priority order; every job on backhouse is ours). Section 3.
3. The chemistry skip: built, tested and measured. Section 5b, lever A.
4. The finest level stays a numerics question for evidence: 2 against 3 levels on two members. Lever D.
5. The spray NaN: fixed at its cause, with a test. Section 6.
6. Lever B's test, CFL 0.8 against 0.4 (the Director, 01:40), registered before its runs. Section 5c.

Sections 1 to 5a measure and estimate without changing the engine. Sections 5b and 6 record the two engine changes on `main`.

Raw outputs, scripts and the instrumentation patch are in `speed_2026-10-24/`.

## The measurement tree

- `main` at 1166721, plus inquiry2's reflux clip (716c3d2's `core/amr.cpp`, `core/amr.hpp`; AMR.md, *Amendment: the reflux clips rho k at zero*), plus print-only step diagnostics (`speed_2026-10-24/stepdiag.patch`). Not committed to `core/`. Built on the Mac (`/tmp/prof_tree`) and on backhouse (`/home/greff/crucible_prof`).
- **The diagnostics** (`CRUCIBLE_STEPDIAG=1`; off, the build is the same code path). In `Flow::stableDt`, each active cell's step is cfl V / rate, with rate = Σ(|u·n| + a)A + 2ν ΣA²/V (`core/flow.cpp:971`). Per level, the diagnostics also take the smallest step with each term removed: no acoustics (|u·n| only, plus diffusion), no diffusion, and flow only. At the binding cell they record the three rates per volume and which diffusivity is largest (ν is the largest of 4/3(μ+μt)/ρ, (λ+λt)/(ρ cv) and D + μt/(ρ Sc)). `Amr::step` keeps each level's record and whether the step was cut to land on an event. The tier 1 driver adds them up over each 10 μs sample window and prints:
  - **what bound each base step:** a level (the smallest 2^l dt_l), an event (the LOX valve, the igniter), the spray's exchange bound, or the sample time or the end;
  - **the rate shares** at the binding cell, and which term is largest;
  - **the global-step speed-up** if a term were gone: steps taken over Σ dt / min(cap, 2^l dt_l without the term), the min over levels, with the same event caps;
  - **the work ratio against local time stepping (LTS):** the cell updates actually done (Σ over steps and levels of active cells times 2^l) over the cell updates needed if every owned cell took its own stable step (Σ dt Σ 1/dt_cell). This is the most LTS could save per unit of physical time, not a count of LTS iterations to a steady state.
- **The wall-time profile:** macOS `sample` of the running process, 8 to 22 s per take, attributed by the innermost engine frame on the main thread (worker 0 of each level's pool) by `sample_cat.py` (categories) and `sample_chem.py` (inside the chemistry). The main thread is one of 4 workers; its time waiting in `Pool::run` for the others measures the load balance.

## 1. Where the wall time goes

### The H2 lead (B4's configuration: 3 levels, the 20 ms lead, 4 threads; Mac, 0.05 to 0.25 ms)

Measured, `lead_profile_samples.txt`. Three takes as level 3 grew from 6,640 to 29,144 cells:

| Phase, % of the main thread's samples | 0.05 to 0.1 ms | 0.15 to 0.2 ms | 0.25 ms |
|---|---|---|---|
| **Chemistry** (`ReactingFlow::react`) | **68.4** | **69.5** | **71.1** |
| Flow stages (`Flow::firstStage`, `secondStage`) | 22.8 | 24.2 | 23.4 |
| of which: boundaries and wall functions | 6.1 | 7.1 | 6.7 |
| of which: transport properties | 3.5 | 4.0 | 3.9 |
| of which: gradients and eddy viscosity | 3.2 | 3.3 | 3.2 |
| of which: refresh (T from e, sound speed) | 2.0 | 2.1 | 2.1 |
| of which: SST sources | 1.1 | 1.1 | 1.0 |
| of which: the flux loops and the rest | 6.9 | 6.6 | 6.5 |
| Stable step (`Flow::stableDt`) | 4.1 | 4.4 | 4.4 |
| The driver's 10 μs sample (c* by Cantera equilibrate, budgets) | 2.5 | 0.7 | under 0.7 |
| AMR fill, reflux, regrid and tags, average down, save and restore | 0.8 | 0.4 | 0.1 |
| Spray | 0 | 0 | 0 |
| Checkpoint and other I/O | 0 | 0 | 0 |
| Waiting in `Pool::run` for the other 3 workers | 0.6 | 0.5 | 0.7 |

- **Inside the chemistry** (`sample_chem.py`): Cantera's rates and their derivatives 78% (`BulkKinetics::updateROP`, Eigen sparse products for the Jacobian, many small mallocs), CVODES and its dense LU 18%, the loop itself 4%, Table B's lookup 0.4%.
- **Why the chemistry costs this much:** Table B placed none of the updates. By 0.15 ms there were about 38 million chemistry updates over the 4 levels (level 3: 16.7 million in the chamber, 8.6 million in the convergent, 4.1 million in the divergent) and none hit the table (the run's `--- t 0.150 ms` line). Pure H2 has f_H = 1, outside Table B's f_H range of 0.133 to 0.182, so every cell goes to CVODES. About 15 μs of thread time per update (derived: 69% of about 213 s of wall on 4 threads, over 38 million updates). In the lead there is no oxygen in the H2, so nothing in these cells can burn; CVODES spends that time confirming no change.
- **Load balance:** the main thread waits for the others 0.5 to 0.7% of the time, so 4 threads run at about 3.9 effective.
- **AMR's own cost** (fill, reflux, regrid, average down) is under 1%. Regridding every 2 base steps does not show.

### The hot regime (the pre-vaporized base mesh: rl10dev, 150 × 16 uniform, Table A, wall functions, 4 threads; Mac)

The tier 1 run cannot reach a hot state tonight. The 3-level lead advances about 0.35 ms an hour on the Mac, and the 1-level run stopped at 1.35 ms on the spray NaN (section 6). The hot state measured here is RL10_DEV's base-mesh chamber (`crucible_chamber_study eqtt 150 16 10e-3 4`, `CRUCIBLE_CHAMBER=rl10dev`, Table A, `law`, `hp`), with the same diagnostics added to the study (`stepdiag.patch`, `tests/chamber_study.cpp`). It is hot and steady by 7 ms (thrust within 0.1%). It has no spray, no AMR and an equilibrium table, so it shows the flow's own cost and step in the hot state, not tier 1's.

Measured, `rl10dev_hot_profile.txt` (one 12 s take at about 3 ms) and `rl10dev_uniform_stepdiag.log`:

| Phase, % of the main thread's samples | Hot base mesh |
|---|---|
| Flow stages | **73.4** |
| of which: boundaries and wall functions | 21.4 |
| of which: transport properties | 12.7 |
| of which: gradients and eddy viscosity | 10.0 |
| of which: refresh | 6.2 |
| of which: SST sources | 3.1 |
| of which: the flux loops and the rest | 20.0 |
| Stable step (it computes the step's transport, gradients and eddy viscosity, which the first stage reuses) | 21.4 |
| Chemistry (Table A's lookup) | 3.5 |
| Other | 1.7 |
| Waiting in `Pool::run`, in total | 4.9 |

- **Cost:** 1.96 μs of wall per cell update on 4 threads, about 7.8 μs of thread time, of which chemistry is 0.08 μs. 25.6 s per simulated ms on 2,400 cells (measured, the run's `cost:` line).
- **The same flow work in the lead** costs about 1.5 μs of wall per cell update on 4 threads (derived: 5.0 μs per update at 0.345 ms, from 1.92 s per base step over 384,568 cell updates, less the 69% chemistry). So once the chemistry is cheap, tier 1's cost per cell update is the flow's, about 2 μs of wall on 4 threads.

## 2. What sets the step

### The H2 lead (3 levels, 4 threads; Mac)

Measured, `prof_mac_stepdiag.txt` (the run's STEPDIAG lines, tabulated by `stepdiag_table.py`):

| t, ms | Base steps | Mean dt, μs | Level 3 bound | Shares at the binding cell: flow, acoustic, diffusive | Largest term (steps) | No acoustics | No diffusion | Flow only | Work / LTS | Level 3 cells |
|---|---|---|---|---|---|---|---|---|---|---|
| 0 to 0.045 | 86 | 0.513 | 81 | 0.39, 0.54, 0.07 | acoustic 82 | 1.94 | 1.04 | 2.05 | 3.27 | 6,640 |
| to 0.095 | 118 | 0.423 | 113 | 0.26, 0.42, 0.32 | acoustic 113 | 1.64 | 1.18 | 2.01 | 2.79 | 12,012 |
| to 0.145 | 123 | 0.406 | 118 | 0.22, 0.41, 0.37 | acoustic 113, diffusive 5 | 1.62 | 1.19 | 1.94 | 2.34 | 17,416 |
| to 0.195 | 125 | 0.358 | 121 | 0.22, 0.34, 0.45 | diffusive 121 | 1.47 | 1.33 | 2.16 | 2.33 | 23,104 |
| to 0.245 | 159 | 0.314 | 154 | 0.18, 0.30, 0.52 | diffusive 154 | 1.41 | 1.51 | 2.45 | 2.42 | 29,144 |

The same run continued to 0.345 ms before it was stopped to free its slot (`lead_to_0.345ms_stepdiag_table.txt`): 250 base steps at a mean dt of 0.195 μs, level 3 binding 250 of them with the diffusive term at 0.62 of the rate, and level 3 at 41,440 cells. At 0.345 ms every level binds at the throat: level 0 at (33, 8), level 1 at (67, 17), level 2 at (134, 36), level 3 at (265, 0) on the axis, where the H2 jet's core has reached the throat. Each level's 2^l dt was then 0.258, 0.249, 0.231 and 0.171 μs. Level 3 holds 86% of the cell updates per base step (331,520 of 384,568, derived).

- **The largest diffusivity at level 3's binding cell** is the eddy thermal one, (λ + λt)/(ρ cv) (kind 4), in every sample from 0.15 ms: the turbulent H2 jet's shear layer on the finest cells.
- **No step was set by the chemistry, the spray or the wall functions.** The chemistry is split from the flow and has no step bound of its own. The spray's exchange bound set 0 steps (no liquid yet). No wall row bound a step.

### One finer level less (1 level above the base, 2 threads; Mac, 0 to 1.35 ms)

Measured, `l1_stepdiag_table.txt`. Level 1 bound 96 to 97% of base steps throughout, with the acoustic term largest (0.39 to 0.65 of the rate) except from 0.75 to 0.9 ms, where eddy thermal diffusion took over (0.43 to 0.50). Removing the acoustic term lengthens the step 1.6 to 2.3 times; removing diffusion 1.0 to 1.3 times; flow only 2.0 to 2.7 times. The work over LTS is 2.6 to 4.0. LOX opened at 1 ms; the run stopped at 1.35 ms before any liquid bound a step.

### The hot regime (the pre-vaporized base mesh, rl10dev; Mac)

Measured, `rl10dev_uniform_stepdiag.log`, 4,501 samples from 1 to 10 ms:

| Rings | Mean dt, μs | Binding cell | Shares at it: flow, acoustic, diffusive | No acoustics | No diffusion | Flow only | Work / LTS (without acoustics) |
|---|---|---|---|---|---|---|---|
| Uniform (tier 1's) | 0.194 | the throat, (33, 13), 100% of samples | 0.408, 0.591, 0.001 | 2.17 | 1.00 | 2.18 | 1.77 (1.42) |
| Stretched to the wall (1.7347) | 0.080 | the throat's wall row, (33, 15) | 0.16, 0.81, 0.03 | 3.7 | 1.04 | 3.8 | 3.8 to 4.0 (2.0 to 2.25) |

(The stretched row is from 1 to 1.5 ms; `rl10dev_stretched_stepdiag.log` holds the rest.)

- **On tier 1's uniform rings the throat sets the step in the hot state.** The flow there is sonic by design, so the acoustic and flow terms are of one size. Removing the acoustic term altogether lengthens the step 2.17 times, and the bound then moves elsewhere at almost the flow-only step (2.18). Diffusion does not matter here (0.1%).
- **On rings stretched to the wall, the throat's wall row sets it,** where the gas is slow and the acoustic term is 81% of the rate. Removing it gains 3.7. That is TECHNICAL_PLAN's "wall strip's acoustic limit". Tier 1 does not stretch its rings (`CRUCIBLE_RADIAL_STRETCHING` defaults to 0).
- **The chamber's low-Mach cells never bind.** They are larger than the throat's (the chamber's radius is about twice the throat's, with the same 16 rings), so their acoustic step is about twice the throat's.

### Answer to the Director's question 1 (the binding constraint)

Measured: **the finest level binds in the lead** (94 to 98% of base steps), at the H2 jet's front and, from about 0.3 ms, at the throat on the axis. The term there is acoustic first (0.54 of the rate), then eddy thermal diffusion in the jet's shear layer (0.52 to 0.62 from 0.15 ms). **In the hot state on the base mesh the throat binds**, acoustic 0.59 and flow 0.41. The chemistry, the spray's exchange bound and the wall functions set no steps in any run. The local Mach number at the binding cells is not small: the throat is sonic, and the lead's binding cells have an acoustic-over-flow rate of 1.3 to 2.1.

## 3. Thread scaling

Measured on one segment: tier 1 development in B4's configuration (3 levels, the 20 ms H2 lead), restarted from a 0.15 ms checkpoint and run to 0.17 ms, with level 3 at about 17,000 cells. The measurement tree has no chemistry skip, so 69% of each step is CVODES. On the Mac my other run was paused during each timed run, but the other projects' slot jobs kept the load at 3 to 7 on 10 cores (6 performance and 4 efficiency, Apple M2 Pro). On backhouse my S7 legs were paused (SIGSTOP), while inquiry2's rfx and B4r runs and the fly project's 5 `smoke_live.py` jobs at nice 0 kept the load at 28 to 38 on 28 virtual CPUs. Our run was at nice 10. Files: `mac_scaling.log`, `bh_scaling.log`, `bh_scaling_time_v.txt`.

| Threads | Mac wall (s) | Mac speed-up | backhouse wall (s) | backhouse speed-up | backhouse CPU share | backhouse user time (s) | voluntary context switches |
|---|---|---|---|---|---|---|---|
| 1 | 186.9 | 1 | 333.7 | 1 | 92% | 309 | 9 |
| 2 | 102.7 | 1.82 | 222.2 | 1.50 | 163% | 362 | 19,638 |
| 4 | 59.7 | 3.13 | 149.9 | 2.23 | 256% | 382 | 94,100 |
| 6 | 46.5 | 4.02 | | | | | |
| 8 | 42.1 | 4.44 | 149.8 | 2.23 | 278% | 413 | 336,343 |
| 16 | | | 144.5 | 2.31 | 329% | 467 | 716,621 |

- **On the Mac,** threads pay up to 6 (the performance cores) and little beyond, at 3.1 times for 4 and 4.4 times for 8.
- **On backhouse under tonight's load, scaling stops at 4 threads.** From 4 to 16 threads, CPU share rises only from 256% to 329% while the user time grows by 22%. The pool's workers wait at each barrier for a worker the scheduler has not run (voluntary switches 94,000 to 717,000). More threads than the CPUs the other jobs leave free only add waiting.
- **Backhouse against the Mac at 4 threads: 2.5 times the wall.** That splits into two factors:
  - **1.56 from the CPU share.** At nice 10 against the fly jobs at nice 0, the run got 256% of 400%.
  - **1.65 from CPU time per step.** At 1 thread the user time was 309 s against the Mac's 187 s. A fixed arithmetic probe (`probe.c`, one thread, no memory traffic) takes 0.70 to 0.78 s of CPU on every backhouse virtual CPU and 0.55 s on the Mac, so 1.3 to 1.4 of the 1.65 is the core itself. The rest is probably the shared cache and memory bandwidth under the other jobs (*inferred*).
- **No E-core signal is visible from inside WSL.** The probe ran pinned (taskset) to each of the 28 virtual CPUs in turn, twice (`bh_probe.log`). Every CPU gave 0.70 to 0.78 s, with no slow subset and no CPU slow in both passes. WSL2 exposes 14 cores with 2 threads each, all at the same nominal clock, and no `cpu_core`/`cpu_atom` split. Hyper-V schedules virtual CPUs onto whichever physical cores it chooses (*inferred* from the uniform probe). Pinning inside the guest therefore cannot select P-cores. That would need a Windows-side setting (processor affinity for the WSL VM), which applies to every WSL user on the machine. That is Ben's call, not a slot-policy item.

**Contention, with every job of ours paused** (the Director's 00:50 item 2; `bh_contention.sh`, `bh_contention.log`). The same segment at 4 threads, with all six of our processes on backhouse stopped (SIGSTOP from 00:56:48 to 01:02:01 CDT, inquiry2 agreeing for rfx and B4r). The fly project's 5 `smoke_live.py` jobs at nice 0 kept running. "Machine busy" is the share of all 28 virtual CPUs busy over the run, from `/proc/stat`, our run included.

| Run (4 threads) | Wall (s) | User time (s) | CPU share | Machine busy |
|---|---|---|---|---|
| Unpinned, nice 10, ours running (section 3's row) | 149.9 | 382 | 256% | 92.9% |
| Unpinned, nice 10, ours paused | **91.2** | 288 | 317% | 64.1% |
| Pinned to virtual CPUs 0, 2, 4, 6, nice 10, ours paused | 108.3 | 304 | 281% | 68.6% |
| Unpinned, nice 0, ours paused | 94.0 | 305 | 325% | 67.2% |

- **Our own jobs were most of the 2.5.** With them paused, the segment took 91 s against 150 s: 1.64 times faster. Backhouse against the Mac at 4 threads is then 91.2 / 59.7 = **1.53**, nearly all of it the slower core (1.3 to 1.4 from the probe).
- **Pinning hurts.** Pinned to four virtual CPUs, the run was 19% slower: the four were shared with the fly jobs' threads, and the run could not move off them.
- **Nice barely matters once ours are paused.** Nice 0 against nice 10 was 94.0 against 91.2 s, within run-to-run spread. Our jobs at nice 10 do not lose to the fly jobs when there are free CPUs; they lose to each other.
- **Oversubscription costs user time, not just wall.** With ours running, the same work took 382 s of user time against 288 s: about 33% spent in barrier waits and cache contention.
- About 15 virtual CPUs were busy with the fly jobs (64% busy less our 11%), so about **13** were free for inquiry's jobs tonight. Our six processes held 30 or more threads between them.
- **Correction (the Director, 01:40).** The 5 python jobs are not another project's. They are the fly project's `smoke_live.py` runs, also under the Director, each holding 89 to 98 threads (BLAS and numba pools) for about 150% of a CPU. The whole load of about 33 was ours: the fly smokes, S7 low and high at 4 threads (getting about 158% each), and inquiry2's rfx (13 threads, 256%). The fly workers are capping their BLAS and numba threads. Where this section says "ours", it means inquiry's and inquiry2's jobs.

**Slot policy for backhouse** (proposed at 01:00, revised for the Director's 01:40 decisions; for inquiry2's resume runner, which already sets nice). Every job on backhouse is ours (inquiry, inquiry2 and the fly workers), so it covers them all:
1. **No pinning** inside WSL. Leave placement to the scheduler.
2. **At most 4 threads per job, counting every pool.** Scaling stopped at 4 under load (section 3), and the probe shows no P-core subset to aim more threads at. For Python jobs that means capping BLAS and numba too (`OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, `NUMBA_NUM_THREADS`); an uncapped fly smoke held 89 to 98 threads.
3. **Threads in total no more than the 28 virtual CPUs.** A job launches only while the 1-minute load plus its own threads stays within 28 (the load gate: 24 for a 4-thread job). One that does not fit waits in the queue rather than slowing every running job.
4. **Priority order (the Director, 01:40): S7 low and high first, B4 second, the fly jobs after;** then the H2-only leads (the 2-level lead) and the measurements. A job above every running job in the order launches without the gate. A lower job yields to a higher one that is waiting: at its next checkpoint, or stopped (SIGSTOP) and continued when there is room. The RL10 central run's place is the Director's when it launches.
5. **Nice 10 stays** for long runs. It costs nothing measurable (section 3); the order is kept by the queue, not by nice.
6. **Windows-side processor affinity for the WSL VM** (to keep it on the P-cores) is the only way to choose physical cores. It applies to every WSL user on the machine, so it is Ben's call, not part of the runner.

**Applied on 24 October** (`relaunch_2026-10-24/`): B4 (job `b4s`) launched at 01:35 CDT without the gate, at 4 threads. The 2-level lead (`lead2s`) passed its gate at 01:33 when the load dipped to 23.9; it was stopped at 01:36 and continues when no `smoke_live.py` runs (`bh_lead2_after_fly.sh`, logged in `runs/lead2/hold.log`).

## 4. Is 60 ms of flow physically required?

**What the registered schedule runs.** The central case is the 20 ms H2 lead, then 60 ms judged over t_L + 40 to t_L + 60, then the continuation at twice the parcel rate to t_L + 85 ms: 105 ms of flow in all. A warm member runs t_L + 60 ms (`tier1_run_definition.md`, lines 25 to 40 and 253). The 60 ms was set as a ceiling (*guessed*, line 49: "rl10dev settled within 20 to 25 ms pre-vaporized, and the liquid's residence adds to that").

**The gas.** Chamber volume to the throat 0.01323 m³ and throat area 0.01237 m², so L* = 1.07 m.
- Hot: c* 2,425 m/s, γ about 1.2 (assumed), RT = (c* Γ)²/γ = 2.06e6 J/kg. The fill and emptying time constant is τ = L* c* / (RT) = **1.26 ms** (derived). The bulk flow-through time is 0.33 m at about 197 m/s, **1.7 ms** (derived: 16.6 kg/s at 1.58 kg/m³ over the 0.0533 m² chamber).
- Gas settles to 0.1% in 7τ, about **9 ms** (derived). The base mesh, pre-vaporized, measured **7.5 ms** to within 0.1% of its final vacuum thrust and **10 ms** to within 0.01% (RL10_DEV.md). The 20 to 25 ms quoted in the definition was for a stricter settledness test, not for the thrust.
- The H2 lead: pure H2 at 200 K (γ 1.4) has c* about 1,325 m/s and τ = **1.7 ms** (derived), so 7τ is about **12 ms**. A 20 ms lead is about 12 τ. The lead's own gate (the means over [t_L − 2, t_L − 1] and [t_L − 1, t_L]) judges it either way.

**The liquid, the slowest scale** (*inferred*, not measured on this engine):
- A 447 μm drop (the central SMD) by the d² law: K = 8 k ln(1 + B)/(ρ_l c_p) is about 2.4e-6 m²/s (k 0.3 W/m K, c_p 3,900 J/kg K, B about 83, ρ_l 1,120 kg/m³), so the lifetime is d²/K, about **83 ms** in still gas. With Ranz–Marshall convection at a slip Reynolds number of about 700 it is about **10 ms** while the slip lasts. Hence **10 to 30 ms**. Drops at the large end of the distribution live 4 to 9 times longer (d² at 2 to 3 times the SMD) unless they leave first.
- Drag relaxation of the same drop is about 9 ms (Schiller–Naumann at Re about 700). Crossing the 0.33 m chamber takes about **7 to 16 ms** from 6 m/s injection. A drop that has not evaporated by then leaves through the throat.
- So the spray's source field reaches its statistically steady state about **10 to 30 ms after t_L** (*inferred*). Gas settling on top of the spray's last change adds about 9 ms. Settledness is expected at roughly **t_L + 20 to t_L + 40 ms** (*inferred*).
- **On S7 (Mascotte, measured on this engine):** the low and high legs (SMD 223 and 574 μm, beside tier 1's 263 to 581 μm) were steady in mass by 25 ms: outlet 73.6 to 75.0 g/s against 73.7 g/s in (21e8bf9). They then failed the 50 K half-to-half test by staying unsteady near the injector, not by drifting.
- The blob legs filled with liquid with a time constant of about 15 ms (derived, b490437) and settle at about 65 to 85 ms. Tier 1 has no blob leg.

**What the rest of the schedule is for.** The 20 ms window (t_L + 40 to 60) and the 25 ms continuation are not physical settling.
- The window's length is set by the noise: batch means need enough integral times τ_int for an SE of 0.2%.
- The continuation tests the parcel count.
- A quiet RANS state reaches the SE bar in a shorter window; a noisy spray state needs a longer one. The S7 amendment's batch means measured the noise there. Tier 1's noise is not measured yet.

**Answer.** Physics needs about 9 ms for the gas and, on the inferred spray scales, 20 to 40 ms after t_L for the liquid and gas together. The 60 ms ceiling is about 1.5 to 3 times that. With the window and the continuation, 105 ms of flow is about 2.5 to 4 times what settling needs (*inferred*). A settle-triggered end (checked every 1 ms, as the definition had before 22 October) would end each run when it is settled rather than at a fixed time. That changes a registered rule before any tier-1 output exists, so it is the Director's decision, as is the lead's length (12 ms by the 7τ arithmetic against 20).

## 5. The levers

**The Director's hypothesis against the profile.** Half of it holds. The step is partly acoustic, but not in a low-Mach chamber: the cells that bind are at the throat (sonic) and in the H2 jet, and their acoustic share is 0.3 to 0.6. Removing the acoustic term altogether lengthens tier 1's step **1.4 to 2.2 times** (measured, every run on uniform rings). Most of the 5 to 11 h per ms is work per cell update times many cell updates:
- in the lead, 69 to 71% of the time is CVODES confirming that pure H2 does not react;
- in every phase, about 2 μs of wall per cell update on 4 threads goes to the flow itself;
- level 3 is 86% of the cell updates.

The large factor a steady method could give does not come from the sound waves. It comes from not marching 105 ms of flow at about 0.2 μs a step, which is 525,000 base steps, each with 8 substeps of level 3.

Speed-ups are on wall time per simulated ms unless stated. The basis is the ensemble plan's 7.5 to 11 h per ms (B4 on backhouse), 105 ms for the central case.

| Lever | Speed-up | Label | Cost to build | Whose decision |
|---|---|---|---|---|
| **A. Skip the chemistry where nothing can react** (**built**, 4a7628f and fc0c9b2; section 5b) | lead: **2.5 to 3.3** derived (chemistry 69 to 71% of the time); measured **1.5 to 2.2** in wall and 1.8 to 2.4 in CPU time (section 5b: 91% of the chemistry updates skipped; the derived figure assumed all of them). Hot phase: none (every hot cell holds fuel and oxidizer). | derived; measured on one segment | done | engineering |
| **B. CFL 0.4 → 0.8** (S7 already runs at 0.8 on the same RK2 and face-sum rate) | **2** if the 1 ms test holds | derived | hours: an environment hook (tier 1 hardcodes 0.4) and a 1 ms test as S7's | tier 1's run definition: the Director |
| **C. Less simulated time** (a 12 ms lead and a settle-triggered end, section 4) | **1.7 to 2.6** on the run | inferred | none | registered schedule: the Director |
| **D. Two finer levels, not three** | **7 to 9** per ms in the lead (level 3 is 86% of the cell updates; the step would rise about 1.3 times) | derived, lead only | none | a numerics question for evidence (the Director, 00:50): run at 2 levels, and check the deviation against 3 levels on two members |
| **E. Where the runs execute** | backhouse took **2.5 times** the Mac's wall on the same segment at 4 threads (section 3: 1.6 from the CPU share at nice 10 under the fly jobs' nice-0 load, 1.6 from slower CPU time per step). The 6.6 first quoted was the seed's first 10 μs window, during level 3's growth; the seed overall was 5.3. | measured | none | machine priority: the Director |
| **F. Threads** | 4 threads: **3.1** on the Mac, **2.2** on backhouse; 8 threads: 4.4 on the Mac, 2.2 on backhouse (section 3) | measured | none | engineering (the slot policy, section 3) |
| **G. The flow's per-update cost** (about 8 μs of thread time per cell update; boundaries and wall functions 21%, transport 13%, gradients 10%) | 1.3 to 2 | guessed (typical of profiling and tuning, nothing identified yet) | 1 to 2 weeks | engineering |
| **H. Low-Mach preconditioning or an IMEX acoustic step, time-accurate** | ceiling **1.4 to 2.2** (the step with the acoustic term removed, measured). Real: less, for the implicit solve's cost per step. On rings stretched to the wall, ceiling 3.7. | measured ceiling | 3 to 6 weeks | engineering, within principle 1 |
| **I. Implicit time-accurate integration** (BDF2 dual time with LU-SGS inner iterations, or ESDIRK with Newton–Krylov; inner iterations converged to a tenth of the step's temporal error, as Bijl, Carpenter and Vatsa do; section 5a) | **5 to 50** in the lead and on 3-level tier 1; **1 to 14** on the 1-level hot base mesh | derived from the measured steps; the inner-iteration counts guessed from the literature | BDF2 and LU-SGS 4 to 8 weeks; ESDIRK and Newton–Krylov 6 to 10 weeks | engineering, within principle 1; the physical step (which frequencies are kept): the Director |
| **J. A steady solve** | out: principle 1 (the Director, 00:50) | | | |
| **K. Local time stepping in the explicit march to the steady state** | at most **1.8 to 4** fewer cell updates per unit of physical time (the work / LTS ratio, measured); the iterations to converge are unknown | measured bound | 1 to 2 weeks | rejected by Ben's earlier ruling |

### 5a. The implicit time-accurate step (the Director's 00:50 item 1)

Principle 1 kept as written: the transient is marched in physical time, and each implicit step's inner iterations are converged, so no pseudo-time remains in the answer. Only the physical step changes. It is then set by the dynamics to be resolved, not by the throat's sonic CFL.

**The physical step the dynamics need** (derived from section 4's scales):

| Dynamics | Time scale | Step that resolves it |
|---|---|---|
| Chamber fill and emptying | τ = 1.26 ms hot, 1.7 ms in the H2 lead | 25 to 60 μs (20 to 50 steps per τ, second order) |
| Bulk flow-through | 1.7 ms | the same |
| LOX valve ramp | 0.5 ms | 25 μs (20 steps) |
| Settle window | 9 ms (gas) to t_L + 20 to 40 ms (spray) | any of the above |
| 1L acoustic mode (if kept) | a/(2L) with a ≈ 1,570 m/s (γ 1.2, RT 2.06e6 J/kg), L 0.33 m: **≈ 2.4 kHz** | ≈ 20 μs (20 steps per period) |
| 1R acoustic mode (if kept) | 3.83 a/(2πR), R 0.130 m: **≈ 7.4 kHz** | ≈ 7 μs |

The axisymmetric model has no tangential modes. Under the performance framing (honest bands, not instability CFD), the acoustic modes need not be resolved. Unresolved modes are damped by an L-stable scheme, and that damping is a resolution choice, like the mesh, to be shown not to move the settled observables. **Δt_phys ≈ 20 μs keeps the 1L; 50 μs keeps only the fill, the flow-through and the valve.** The igniter's kernel and the light-off are faster than either. An adaptive step (ESDIRK's embedded error estimate) would take small steps there and grow afterwards.

**Against today's explicit step** (measured, section 2; CFL 0.4):

| Where | Explicit step | Δt_phys 20 μs over it | Δt_phys 50 μs over it |
|---|---|---|---|
| Hot base mesh, 1 level, uniform rings | 0.194 μs | 103 | 258 |
| Hot base mesh, stretched rings | 0.080 μs | 250 | 625 |
| The lead, base step (3 levels) | 0.195 to 0.41 μs | 49 to 103 | 122 to 256 |
| The lead, level 3 (2^-3 of the base step) | 0.024 to 0.051 μs | 390 to 830 | 980 to 2,100 |

**The cost per implicit step.** An implicit method steps all levels together, with no subcycling. One inner iteration sweeps every cell once.
- **BDF2 dual time, LU-SGS inner iterations** (Jameson 1991; Yoon and Jameson's LU-SGS): one residual evaluation and two sweeps with a scalar diagonal. That is about 0.75 to 1 explicit RK2 step's cost per cell (*guessed*). Converging to a tenth of the temporal error takes 10 to 50 inner iterations with low-Mach preconditioning in the inner loop, and up to about 100 without (*guessed*, typical reported counts). The chamber's Mach number is about 0.13, and the inner acoustic CFL at 20 μs is tens on the base mesh and hundreds on level 3 (*derived*). Hence **7.5 to 100 explicit-step equivalents per physical step**, over all cells once.
- **ESDIRK4 with Newton–Krylov** (Bijl, Carpenter, Vatsa 2002): 5 implicit stages, each a Newton solve of 3 to 5 iterations of 10 to 30 Krylov iterations. That is about 150 to 750 residual evaluations per step (*guessed*). Bijl et al. found ESDIRK4 cheaper than BDF2 at tight error targets because its steps can be several times larger. At our looser target (resolving a 1.3 ms fill to 0.1%), BDF2 at 20 to 50 μs is the cheaper start.

**Speed-up per simulated ms** (*derived* from the measured steps; the inner counts are guessed):
- **The lead and 3-level tier 1.** Explicitly, level 3 subcycles 8 times and holds 86% of the updates, so a base step costs about 6.5 all-cell sweeps (384,568 updates over about 59,000 cells at 0.345 ms). At a base step of 0.195 μs, 20 μs is about 670 all-cell RK2 steps. Implicitly it is 7.5 to 100. **5 to 90; quoted as 5 to 50** for the implicit operator's extra memory traffic and Jacobian work per cell.
- **The 1-level hot base mesh.** 103 explicit steps per 20 μs against 7.5 to 100: **1 to 14**.
- With CFL 0.8 (lever B) the explicit side halves, so these factors halve against B.

**What it does not give for free:**
- **Convective CFL.** At 20 μs, level 3's convective CFL in the H2 jet is about 50 (1,000 m/s over cells of about 0.4 mm). The jet's front is then smeared over the steps, not resolved. For a settled RANS state this does not matter, but it has to be shown on the settled observables (the same evidence as D: two steps, two members).
- **Chemistry splitting.** One Strang split per 20 μs step puts the split error on flame zones whose reaction–diffusion time is about 10 μs to 2 ms (δ²/D with δ 0.1 to 0.5 mm). Split once per physical step, the error must be measured; inside the inner loop it costs CVODES per inner iteration.
- **The spray exchange.** The exchange bound (registration B) is an explicit-coupling limit. Small drops relax faster than 20 μs, so the exchange needs point-implicit coupling.
- **Threads.** LU-SGS's sweeps are sequential within a block. Threading needs block-Jacobi across patches or a coloured ordering, at some loss of convergence per iteration.

**Build** (*guessed*): **BDF2 and LU-SGS 4 to 8 weeks.** That covers:
- the synchronous-level step on AMR;
- the LU-SGS sweeps with point-implicit SST, wall-function and axisymmetric sources;
- preconditioning in the inner loop;
- the spray and chemistry coupling;
- verification: temporal order on an acoustic pulse and a fill transient, Bijl's inner criterion, and the lead and the hot base mesh against the explicit march.

ESDIRK with Newton–Krylov takes 6 to 10 weeks. TECHNICAL_PLAN lists implicit time-accurate integration as an open item and dual time as step 16 (parked).

**What stacks.**
- **A is built and measured** (section 5b): 1.5 to 2.2 in wall (1.8 to 2.4 in CPU time) in the lead. It does nothing after t_L, because every hot cell holds fuel and oxidizer.
- **B, C and E need no build.** B needs a hook and a 1 ms test; C and E are the Director's decisions.
  - In the lead, A × B is 3 to 4.4, if B's factor of 2 holds in its test.
  - After t_L, B alone halves the cost per ms. The hot phase's cost per ms in tier 1's configuration is not measured.
  - C shortens the run 1.7 to 2.6 times.
  - E under the slot policy (section 3) recovers the 1.6 that our own jobs cost each other on backhouse.
- **D is the largest factor with no build.** It is a resolution question for evidence: 2 against 3 levels on two members.
- **I, the implicit time-accurate step, is the largest lever within principle 1:** 5 to 50 in the lead and on 3-level tier 1. It costs 4 to 8 weeks to build and verify. It replaces the explicit step, so it does not stack with B, and its factor against B halves.
- **H is at most 1.4 to 2.2,** because the binding cells are sonic, and I subsumes it.

**Which one to test first (engineering).**
1. B: one 1 ms test on the Mac, with only a hook to build (CFL 0.8, as S7 already runs).
2. D's evidence: 2 against 3 levels on two members. It needs machine time but no build.
3. If the Director wants I, start with a prototype on the 1-level base mesh: BDF2 with LU-SGS, the fill transient against the explicit march, about a week (*guessed*). It replaces the guessed inner-iteration counts with measured ones before 4 to 8 weeks are committed.

**The Director, 01:40 24 October:** prototype first, about a week, on a bounded case, with its success criteria written before building. They are in [`IMPLICIT_PROTOTYPE.md`](IMPLICIT_PROTOTYPE.md) (the acoustic pulse and the C1 cold start; BDF2 with LU-SGS on one level and one thread; criteria 1 to 6; measurements M1 to M5; a proposed rule: build if the projected factor is at least 5, not if below 2).

### 5b. The chemistry skip, as built (lever A; the Director's 00:50 item 3)

Built on `main` (4a7628f; extended and corrected in fc0c9b2), in `ReactingFlow::react` after Table B's lookup and before the CVODES call. It acts only on the finite-rate cells Table B does not place; equilibrium cells are untouched.
- **Exact:** a cell with no reactive species present is skipped. A species is reactive if some reaction has a nonzero net coefficient for it (`Reaction::reactionSpecies`). With none present, every rate of progress is zero, so the integration returns the state unchanged. Tier 1's N2 fill is such a cell.
- **Bounded:** a cell whose only reactive species is H2 (or only O2) is skipped at or below 500 K. This is not exact, because pure H2 can dissociate.
  - At construction the mechanism is checked: every reaction direction startable from H2 alone (with inert colliders) may change only H2 and H, and likewise for O2 and O.
  - From no atoms, the composition can only move toward the dissociation equilibrium. The change at any step length is therefore below the atom's equilibrium mole fraction: x_H 2.2e-18 and x_O 3.7e-21 at 500 K and 1 Pa, smaller at higher pressure and lower temperature.
  - That is below CVODES's absolute tolerance (1e-12) by six orders. The lead's cells (H2 and N2, 200 to 300 K) are covered.
- **Trace atoms (the extension):** the molecule may carry its own atom at 0 < Y ≤ 1e-14 (`kAtomTrace`).
  - **Why it was needed.** The first build skipped only 7% of level 3's updates in the lead, against 93% expected. CVODES leaves the H2 it integrates with a positive trace of H (criterion 3 measures changes of up to 6e-36 over 1 ms), advection spreads it, and a cell with any H counted as reactive (*inferred*; with traces up to 1e-14 allowed, 91% of level 3's updates are skipped).
  - **When it is allowed.** At construction the mechanism is checked again: no direction that needs no reactive species but the molecule and the atom may change any other species. Only dissociation and recombination qualify.
  - **The bound.** The cell can then only move toward the dissociation equilibrium, so it changes by at most the larger of Y_atom (at most 1e-14) and the atom's equilibrium fraction. The heat of recombining 1e-14 of H is about 1e-10 K (derived).
- **A fault found by the extension's test.** The skip built its species mask from `Flow::massFractions`, which clips a negative partial density to 0. A cell with a negative trace was therefore skipped, and kept the negative partial density that the integrator's path clips. That breaks 4a7628f's own rule that a negative fraction counts as present. The mask now reads the partial densities.
- **Test** (`tests/chemistry_skip_tests.cpp`, ctest `chemistry_skip`; criteria written before each first run):
  - (1) The reactive and startable sets of h2o2.yaml and h2o2_highT.yaml, written independently of the engine.
  - (2) The dissociation bound below 1e-16 at 500 K and 1 Pa.
  - (3) The integrator itself changes every Y_k by under 1e-16 over 1 ms, from pure and diluted H2 and O2 at 200 to 500 K and 1e2 to 1e7 Pa.
  - (4) A duct in zones (H2 and N2 cold, N2 hot, O2 and N2 cold, H2 and N2 at 600 K, a cold H2–O2 mix). The skipped cells are exactly the first three zones and keep their state bit for bit. On against off, every partial density agrees within 1e-15 of the density, and the integrated zones bit for bit.
  - (5) Trace atoms: (a) the closure of {H2, H} and {O2, O} in both mechanisms; (b) a duct with Y_H 1e-15 and Y_O 1e-15 (skipped), Y_H 1e-10, Y_H 1e-15 at 600 K and Y_H −1e-15 (each integrated); (c) the skipped cells unchanged bit for bit, on against off within 1e-14 of the density, the integrated cells bit for bit.
  - Run history: criterion 4's run 1 failed 4(a) on the integrated-cell counter, which was fixed with the criteria unchanged; run 2 passed. Criterion 5's run 1 (`chemistry_skip_tests_c5_run1.txt`) failed 5(b) and 5(c): 40 cells skipped where 32 were expected, the 8 extra being the negative-trace cells, with a gap of 1.0e-15 of the density. That is the fault above, fixed with the criteria unchanged. Run 2 (`chemistry_skip_tests_run.txt`) passes every check.
- **Counted:** `ReactingFlow::Stats` has skippedCells and integratedCells (checkpointed). Tier 1 prints them per level beside Table B's tally. `CRUCIBLE_CHEM_SKIP=0` turns the skip off in tier 1, for timing.
- **Measured in the lead** (`skip_timing.sh`, `skip_timing.log`). The segment was 3 levels, the 20 ms lead, 4 threads, 0 to 0.1 ms on the Mac. Skip on and off were alternated twice, under the other projects' load of 8 to 12.
  - **Skipped / integrated updates, skip on:**
    - level 0: 494,938 / 20,006;
    - level 1: 674,690 / 71,674;
    - level 2: 3,764,769 / 343,767;
    - level 3: 12,383,954 / 1,288,078.
  - That is 91% of all 19.04 million updates skipped (90.6% at level 3). Skip off integrates all of them.
  - **Wall time:** on, 80.2 and 78.7 s; off, 176.5 and 122.0 s. That is **1.5 to 2.2 times** faster, 1.9 on the means.
  - **CPU time:** on, 239 and 249 s; off, 565 and 457 s. That is **1.8 to 2.4 times** less, 2.1 on the means.
  - **The first build**, before the trace extension, skipped 10% of the updates (7% at level 3). It saved 2.6% of the wall: 130.9 s against 134.4 s.
  - **Results:** the two skip-on runs are bit-identical. On against off, the 12 series rows differ by at most 3.3e-8 relative at 0.1 ms. The skip changes a cell by at most 1e-14 per step, and the flow's own sensitivity grows that (*inferred*).
  - The derived 2.5 to 3.3 assumed every lead cell skipped. The 9% still integrated carry more than a trace of H or are above 500 K.

### 5c. Lever B's test: CFL 0.8 against 0.4 (the Director, 01:40: "run the 1 ms test at the larger step and report accuracy against the current step")

**Registered before any run (24 October, about 01:43 CDT).** The hooks: `CRUCIBLE_T1_CFL` sets `Definition::cfl` and `CRUCIBLE_T1_SEED` the spray's seed, both in development only (the central case keeps CFL 0.4 and seed 1). Every run is tier 1 development on `main` at the commit that adds them, on the Mac at 4 threads through the slot harness, with `scripts` and logs in `cfl_2026-10-24/`.
- **L4 and L8:** 3 levels, the 20 ms H2 lead, 0 to 1 ms, at CFL 0.4 and 0.8. The finest level binds here (the H2 jet, then the throat's axis; section 2).
- **H4 and H8:** 1 level, the lead shortened to 1 ms (LOX at 1 ms on its 0.5 ms ramp, the igniter at 1.2 ms), 0 to 2 ms, at CFL 0.4 and 0.8. This is the hot transient after LOX opens, where B alone acts (A does nothing there).
- **H4s:** H4 with the spray's seed 2. The parcels' random draws depend on the step, so H8 against H4 holds sampling noise as well as the step's error. H4s against H4 measures that noise (the floor). Before LOX opens, H4s must equal H4 bit for bit.

**Measured:** the wall time, the base steps and the redone steps; and the series (every 10 μs) of the outlet mass flow, F_vac, p0, the face pressure and the hot-gas wall heat. For each quantity, the difference at each common row is divided by the CFL 0.4 run's mean |Q| over the window. The windows are [0.1, 1.0] ms for L and [1.0, 2.0] ms for H. Its largest value and its RMS are reported, with the difference at the window's end.

**Criteria, stated now:**
- **(a) Stable.** L8 and H8 end with exit 0, budgets within 1e-11 and no non-finite value. Their redone steps are reported beside L4's and H4's.
- **(b) Accurate.** For F_vac, the outlet flow and p0:
  - In L, the RMS difference is at most 0.2%. That is the run definition's SE target for the window means; a step error below it does not move a reported mean beyond its own noise.
  - In H, the RMS difference of H8 against H4 is at most the larger of 0.2% and the floor (H4s against H4).
- **What follows.** If (a) and (b) hold, B holds for these two regimes: about 2 times fewer steps, with the wall time as measured. Whether tier 1's registered CFL changes stays the Director's decision. If either fails, the record says which and by how much.
- RK2's error is second order, so the step's own error at CFL 0.8 is about 4/3 of the 0.8 to 0.4 difference, and at 0.4 about 1/3 (*derived*, where the difference is not noise).

## 6. Found on the way: the spray NaN after LOX opens (fixed)

The 1-level run (`/tmp/prof_l1`: the lead shortened to 1 ms, LOX opening at 1 ms on its 0.5 ms ramp, the igniter at 1.2 ms) stopped at about 1.35 ms with "Spray: the film's transport is not positive at nan K". 27b135f's guard caught the NaN film temperature, but not what made it. The same would have stopped B4 at about t_L + 0.35 ms, and the central run with it (inferred: the same code path once LOX flows).

**Where it came from** (measured; `nandiag.patch` on the same case from t = 0, `mac_nan_repro.sh`, `nan_diag.log`). The first non-finite film temperature was:
- a 0.70 mm drop at 97.427679 K, 1.3e-6 K below its saturation temperature at 2.034e5 Pa;
- so 1 − Y_s = 7.1e-9 and B_M = 7.7e7;
- in gas at 187.12 K with Y_O2 0.450 and the rest H2.

B_T had come back infinite.

**The cause.** `Spray::rates` solves ln(1 + B_T) = φ(B_T) ln(1 + B_M) by secant steps on the residual expm1(φ ln(1 + B_M)) − B_T.
- With ln(1 + B_M) at 18, an iterate below the root whose φ is above about 39 sends the residual past the double range.
- The fallback step B_T + r is then infinite.
- The step test |next − b| ≤ 1e-11 next holds for next = inf, so the solve returned B_T = inf.
- The film's shape a_T is then infinite, and its temperature inf/inf = NaN.

A second fault in the same loop showed up only in the test. With the rest N2 at a slip of 3 and 30 m/s, a huge residual was followed by a tiny secant step far from the root, which passed the same step test. B_T came back 16% off its own equation, with no NaN to show it.

**The fix** (26f6fd9, `core/spray.cpp`, the secant's stop test). An iterate counts as settled only if the step is finite, at rounding (1e-11 relative), *and* the residual is within 1e-6 of B_T. The loop also stops on an infinite iterate. Anything not settled goes to the bisection that was already there, bracketed from the last finite iterate.

**Test** (`tests/spray_transfer_tests.cpp`, ctest `spray_transfer`; criteria written before the fix, from the diagnosed drop). The case is the drop and gas above, with the rest H2 or N2, at slips of 0, 3 and 30 m/s, at T_sat − δ for δ from 3 K to 1e-12 K (121 log steps) plus the diagnosed 97.427679170131341 K.
- (a) No call throws, and every rate is finite.
- (b) Where the drop evaporates, c_pF ṁ = conductance × B_T within 1e-9 relative.
- (c) `spray_verification` prints the same lines before and after the fix.
- (d) The tier 1 case runs past 1.4 ms.

Results:
- **Before the fix:** 6 failures (`spray_transfer_before_fix.txt`). With the rest H2 at 3 and 30 m/s, 58 and 61 of 122 calls threw or gave a NaN. With the rest N2 at 3 and 30 m/s, B_T was 16% off.
- **After:** (a) and (b) pass, the largest error 6.7e-16 (`spray_transfer_after_fix.txt`). (c) holds: the two outputs differ only in their wall-time lines. (d) passes. The same case on main with the fix ran to its 1.5 ms end with exit 0, with 0 non-finite values in its 150 samples, 13,984 parcels held and 4,580 flashes by then (`nan_fix_rerun.txt`, `nan_fix_rerun_series.csv`).

B_T itself is legitimately huge here: up to 2.6e52 with the rest H2 at zero slip, where φ ln(1 + B_M) is about 120. Only the stop test was wrong.

## Limits

- **The hot regime is not tier 1's.** The hot measurements are the pre-vaporized base mesh with an equilibrium table. Tier 1's hot phase (3 levels, the spray, Table B) was not reached tonight. Its step and its share of Table B misses after t_L are not measured. B4r's progress lines will give both.
- **One machine state per number.** The Mac was shared with other projects' slot jobs during every run. The macOS `sample` takes are 8 to 22 s each.
- **The work / LTS ratios bound what local time stepping could save per unit of physical time.** They do not count the iterations a local-step march would need to converge.
- **Speed-ups for levers G, I and J are guessed,** with the arithmetic shown. Only the step ceilings (no acoustics, no diffusion, flow only) and the phase shares are measured.
- **The thread scaling and the contention runs used the lead without the skip,** where 69% of each step is CVODES and parallel. With the skip, the flow's share and its barriers dominate. Scaling with the skip on is not measured.
- **Contention was measured once per configuration,** on one evening, under the fly project's 5 jobs. The 1.6 that our own jobs cost is for tonight's mix of six processes.
- **The skip's timing is one segment** (0 to 0.1 ms of the lead, on the Mac), with on and off alternated twice under the other projects' load (section 5b).
