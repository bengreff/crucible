# Speed plan: from 46 to 75 days to under an hour per RL10 run (4 November 2026)

Written for the roadmap's AT RISK inquiry box. Ben's bar (24 October): a full RL10 run in under an hour, by a GPU port, and the RL10 is one engine of 20 or more. Principle 1 (always time-dependent) is canon, so every lever here marches the flow in physical time. SPEED.md holds the measurements this plan builds on; this file does the arithmetic and the order.

Labels: **measured** (a run or a profile), **derived** (arithmetic on measured numbers, shown), **guessed** (an estimate with its basis stated), **spec** (a manufacturer's figure).

## 1. Verdict

- No CPU lever, nor any product of CPU levers, reaches the bar. The CPU levers that exist multiply to about 40 to 110 (derived, section 3). The gap is 1,100 to 1,800 (derived, section 2).
- A GPU port of the gas step is necessary. On its own it is not enough, guessing its rate on backhouse's card at 50 to 400 times backhouse's 4 threads.
- Under the guess, the bar is reachable as the product **GPU port × D (two finer levels, not three) × B (CFL 0.8)**. That product is 700 to 7,200 (derived from the guess). Two conditions:
  - the chemistry's CVODES fallback and the spray both leave the CPU's critical path (section 4);
  - D and B pass their accuracy tests.
- The implicit step (I) is the hedge if D fails or the GPU's rate lands low, and the CPU path's largest lever.
  - On tier 1's full lead physics it measures **7.1 over the lead and 11 in its second half** (5 November, backhouse, base mesh, one thread). It passes window L at time tolerance 3e-3 (section 7.4). That beats the inviscid prototype's 5 to 6, because the driver's 10 μs samples, not the error, set the step.
  - Not yet shown: the hot window (the p0 instrument stalls there first), AMR, and robustness. One trajectory in three failed on a latent collapse of omega, now guarded (7.4).
  - Its inner solver (LU-SGS) is sequential, so on a GPU it needs a different smoother.
- Nothing here needs Ben unless the GPU's measured FP64 rate falls short. Then the choice is mixed precision (engineering) or a card with fast FP64 (money: Ben's call). Section 6 has the details.

## 2. The size of the job

- **The basis** (tier1_ensemble_plan.md): B4 on backhouse at 4 threads, before the chemistry skip, ran 7.5 to 11 h per ms of gas (11 measured, 7.5 inferred without contention). The spray multiplies that by 1.5 to 1.7 (measured on S7).
- **The central case** is 20 ms of lead plus 85 ms after it, so 1,100 to 1,800 h, or 46 to 75 days (derived there).
- **The bar** is under 1 h, so the factor needed is **1,100 to 1,800** (derived).
- **Cell updates per central run: about 2e11 to 3e11** (derived, two ways that agree).
  - (a) 105 ms at a base step of about 0.2 μs is 525,000 base steps (SPEED.md 5). At 0.345 ms there were 384,568 cell updates per base step (measured). That gives 2.0e11 at the early mesh's size.
  - (b) 790 to 1,155 h of gas time (105 ms × 7.5 to 11 h per ms) is 2.8e6 to 4.2e6 s. Divided by about 12.5 μs of wall per update, that gives 2.3e11 to 3.3e11.
    - The 12.5 μs is the Mac's 5.0 μs per update with chemistry (SPEED.md 1) times backhouse's 2.5 (lever E).
- **The rate the bar needs:** 2e11 to 3e11 updates in 3,600 s is **6e7 to 9e7 cell updates per second** (derived), with the chemistry and the spray inside the same hour.
- **Today's rates** for the flow alone, on 4 threads (derived from SPEED.md 1's 2 μs of wall per update on the Mac and lever E's 2.5):
  - the Mac: 5e5 updates per second;
  - backhouse: about 2e5.

## 3. Every lever

Speed-ups are on wall per simulated ms, against the basis above.

| Lever | Speed-up | Label | What it takes | Status tonight |
|---|---|---|---|---|
| **GPU port of the gas step** (reconstruction, HLLC, transport, SST, wall functions, the table chemistry; FP64) | **50 to 400** against backhouse's 4 threads | guessed: 0.69 TFLOPS FP64 (spec, RTX 4070 Ti SUPER, 1/64 of FP32) at 30 to 60% use over 5,000 to 20,000 FP64 flops per update gives 1e7 to 8e7 updates per second; memory at 672 GB/s (spec) and about 3 kB per update allows 2e8, so the FP64 units bind | one level first, then AMR's fill, reflux and regrid on the device; guessed 6 to 12 weeks | not started; GPU in use by another job tonight (98%, 11.4 GB of 16), so no benchmark |
| **D. Two finer levels, not three** | **7 to 9** in the lead | derived (level 3 holds 86% of the updates, and the step rises about 1.3 times) | an accuracy test: 2 levels against 3 on two members | not run |
| **B. CFL 0.4 to 0.8** | **2** | derived; the 1 ms test's outputs were lost (section 5) | rerun the test (hours on backhouse) | to rerun; moot under I, whose step no CFL sets |
| **A. The chemistry skip** | **1.5 to 2.2** in the lead's wall | measured (built) | done | done, in the basis's successors |
| **Tables for every chemistry state** (no CVODES in the march) | lead: **1.5 to 2.2** more on the CPU (the fallback is still 33 to 55% of the lead's wall after A, derived from A's 1.5 to 2.2 over a 69 to 71% share); on the GPU path it is **necessary**, not optional (section 4) | derived | widen Table B's f_H range to the lead's pure-H2 and mixing states, or a second table | not started |
| **Spray on the GPU, or fewer parcels** | removes the 1.5 to 1.7 multiplier on 85 ms; **necessary** on the GPU path (section 4) | measured multiplier | a device parcel kernel, or a parcel-count study (kParcelRate 5e5) | not started |
| **I. Implicit time-accurate step** (BDF2, dual time, LU-SGS) | **7.1** over tier 1's lead [0.5, 20] ms and **11** on [10, 20] ms (backhouse, base mesh, one thread, time tolerance 3e-3, section 7.4). The Mac's 2.6 to 3.3 (7.3) carried the old p0 instrument. The inviscid 1-level prototype gave 5 to 6 | measured | built: the port to the viscous, SST and wall-function step on tier 1's base mesh, with time-error step control and the omega guard (section 7) | window L passes (i10u); window H waits on the p0 instrument for oxygen-bearing cold gas; the step is capped by the driver's 10 μs samples; the explicit path unchanged |
| **C. Less simulated time** | 1.7 to 2.6 against the old schedule; the approved central case (lead plus 85 ms) already takes most of it | inferred | the Director's schedule | taken |
| **G. The flow's per-update cost** | 1.3 to 2 | guessed | profiling and tuning; folded into the port | not started |
| **F. Threads** | 2.2 at 4 threads on backhouse, 3.1 on the Mac; 4.4 at 8 on the Mac | measured | the slot policy | in use |
| **E. The machine** | backhouse at 4 threads runs 2.5 times the Mac's wall | measured | none | known |
| **H. IMEX acoustics** | ceiling 1.4 to 2.2 | measured ceiling | 3 to 6 weeks | dominated by I |
| **J. A steady solve** | out | principle 1 | | out |
| **K. Local time stepping** | at most 1.8 to 4 | measured bound | | rejected by Ben |

**The CPU-only product:** D × B × the tables × G × F (8 threads against 4, Mac) is 7 to 9 × 2 × 1.5 to 2.2 × 1.3 to 2 × 1.4 = **38 to 110** (derived, before I). The run time comes out at 10 to 47 h.
- With I's measured 7 to 11 in place of B, which it makes moot, the product is 7 to 9 × 1.5 to 2.2 × 1.3 to 2 × 1.4 × 7 to 11 = **130 to 610** (derived). The run time is 1.8 to 13 h.
- That assumes the gains compose, which is unmeasured: I is built on the base mesh alone, and its hot window is not run.
- Both figures assume the spray's cost falls with the gas's, which is not measured.
- The bar is not met, and this path assumes backhouse's whole CPU, which is shared.

**The GPU product:** GPU × D × B = 50 to 400 × 7 to 9 × 2 = **700 to 7,200** (derived on the guessed GPU rate), against 1,100 to 1,800 needed.
- At the low end of the guess, the bar is missed by up to 2.6 times.
- Mixed precision or I would close the low end. Mixed precision means FP32 fluxes with FP64 sums, which keeps conservation exact to the sums' rounding.

## 4. What caps the GPU path (Amdahl)

- **The chemistry fallback.** After A, CVODES is 33 to 55% of the lead's wall (derived above). If the gas went infinitely fast and CVODES stayed on the CPU, the lead would gain at most 1.8 to 3 times. So every state the march meets must be tabled, or integrated on the device.
  - In the hot phase, Table A's chemistry is already 0.08 of 1.96 μs per update (measured, SPEED.md 1), so it is cheap.
- **The spray.** It adds 0.5 to 0.7 of the gas time over 85 ms (measured multiplier). Left on the CPU, it caps the run's gain at 2.4 to 3 times (derived: 1.5 to 1.7 over 0.5 to 0.7).
- **The driver's samples** (p0 and c* by Cantera's equilibrate every 10 μs) were 0.7 to 2.5% of the AMR lead's wall (measured). On the base alone in the cold lead they grew to 84 to 100% of it from about 11 ms (measured, 7.3), until 7eb8956 took oxygen-free cold gas frozen. At a 1,000-times gas speed-up they would be most of the run, so the sample becomes a table lookup or a sparser sample.
- **AMR's bookkeeping** (fill, reflux, regrid) is under 1% on the CPU (measured). On the device it must stay there, so it moves too.

## 5. Ranked order (by factor, certainty and cost)

1. **GPU port of the gas step, one level, FP64**, with a measured updates-per-second figure on backhouse's card in its first week. That figure decides everything below. It is the only lever that can reach the bar.
2. **D, two finer levels**: the largest factor that costs only runs. Its accuracy test should run now, on the CPU, while the port is built.
3. **B, CFL 0.8**: rerun the lost test.
4. **Tables for every state**: necessary for the port.
5. **The spray on the device**, or a parcel-count study: necessary for the port.
6. **I, the implicit step**: built on tier 1's base, measured at 7.1 to 11 in the lead (7.4). It is the hedge for the GPU path's low end, and on the device it needs a parallel smoother in place of LU-SGS's sweeps. Before any more weight goes on it:
   - the hot window, after the p0 instrument takes cold oxygen-bearing gas frozen;
   - samples from BDF2's interpolant, so the step can pass 10 μs;
   - k and omega in the time error;
   - unconverged steps rejected rather than accepted.
   It ranks below the port because it does not reach the bar on the CPU (130 to 610 against 1,100 to 1,800, derived).
7. G, F and E: CPU-side and small; G is folded into the port.

## 6. What needs Ben

None tonight. One possible call later: if the measured FP64 rate of the 4070 Ti SUPER falls short and mixed precision fails its accuracy test, a card with fast FP64 is a purchase (Ben's call). The GPU is shared (someone's job held 98% of it and 11.4 GB tonight), and the 17:00 to 22:00 cap is 0.3 duty. A one-hour run therefore assumes the card is free when the run goes.

## 7. Tonight's lever: the implicit step on tier 1's physics

### 7.1 What was built (4 November, uncommitted at the time of 7.2)

- **`Flow::implicitStep` on the full gas step** (`core/flow_implicit.cpp`):
  - transport, SST and wall functions, on one level without regions;
  - the SST sources Strang-split around the transport exactly as in the explicit step: half a step with the coefficients at the start, half with those at the new state;
  - the wall functions prescribe the wall cells at every iterate and at the end; the end's prescription and the clip of rho k at 0 are booked as the explicit stages' are;
  - LU-SGS gains a thin-layer viscous term per face, l_v = nu |S|^2 / V, with nu as `stableDt` takes it;
  - the flux Jacobian's product carries rho k and rho omega (and k in the pressure).
- **The history is now the step's own increment**, not U^n - U^(n-1).
  - Why: with an operator split after the step (the chemistry, the spray, the half sources), the difference form settles where a R + (the split change)/h = 0, so at a = 2/3 the transport is weighted by 2/3 against the split change (derived). The increment form settles where R + (the split change)/h = 0.
  - Without split operators the two are the same number: case P at 192 columns is identical to 17 digits in error and both budgets on all six implicit runs (measured, against 5d4a41c's binary).
  - With Table A's projection inside the step, the increment is taken after the projection (it is part of the step's equations, keeping rho, E and the elements).
- **`Amr` takes 0 finer levels**, and `Amr::implicitStep` steps the base alone: it lands on events and lets the step at most double from one step to the next (BDF2's variable step is zero-stable below a ratio of 1 + sqrt 2).
- **Tier 1:** `CRUCIBLE_T1_IMPLICIT="dt,tolerance[,iterations[,pseudo CFL[,extrapolate]]]"`, development only, with `CRUCIBLE_T1_LEVELS=0`. The BDF2 history is not checkpointed, so a restarted run's first step is backward Euler.
- **Probe** (measured, 0.1 ms, base 149x16, one thread): the explicit stable step is 5.2e-7 s in the lead's first 0.1 ms; the implicit at 10 μs took 21 inner iterations per step, 0 rejected. Too short to judge anything; it fixed the test's lengths below.

### 7.2 The agreement test, pre-registered (written 4 November at 22:55 CST, before the runs)

**Case.** Tier 1 development, base 149x16 alone (`CRUCIBLE_T1_LEVELS=0`), central element, one thread each, Table B, from t = 0 to 30 ms:
- the H2 lead to t_L = 20 ms (transport, SST, wall functions, the chemistry hook);
- then the LOX, the igniter at t_L + 0.2 ms and 10 ms of hot flow with the spray.

**Runs.**
- e1: explicit at CFL 0.4 (the reference), spray seed 1.
- e2: explicit, spray seed 2 (the hot window's noise floor; identical to e1 before t_L, which is checked).
- i10 and i5: the implicit step at most 10 μs and 5 μs a step, seed 1. Tolerance 8.975e-5, at most 400 inner iterations, pseudo CFL 3, extrapolated first iterate (the prototype's chosen settings).

**Observables.** From the 10 μs series: outlet mdot, F_vac, p0, p_face, q_wall (q_face reported, not judged).

**Window L, the lead, [0.5, 20] ms.** For each observable:
- the RMS of (implicit - e1) over the window's samples, over e1's RMS magnitude on the window;
- the largest |implicit - e1| over the same scale.
- **Pass at a step size:** every RMS ≤ 0.5% and every largest ≤ 2%.

**Window H, the hot flow, [23, 30] ms.** The spray makes the flow chaotic, so the pointwise difference is not judged (it is reported). For each observable, the window mean's difference from e1's over e1's mean.
- **Pass:** |difference| ≤ max(1%, 2 × the e2 - e1 difference) for mdot, F_vac, p0 and p_face; ≤ max(3%, 2 × the floor) for q_wall.
- If e1 does not reach 30 ms, window H is not judged, and the reason is reported.

**Also judged.** The gas budget within the run's own 1e-11, and no non-finite value, in every run.

**Reported, not judged.**
- The speed-up: e1's wall over each implicit run's, per window, from the progress lines. The four runs share the Mac, one thread each, so this is measured under that load.
- The observed order from i10's and i5's RMS differences.
- The inner iterations, rejections and unconverged steps.

The step size used on tier 1 is the largest that passes both windows.

### 7.3 Results (4 and 5 November, Mac, one thread each; all numbers measured unless labelled)

**The pre-registered runs failed, so the pre-registered test has no passing step size.**
- i5 (5 μs fixed): stopped at 2.49 ms; all 14 step reductions were rejected.
- i10 (10 μs fixed): CVODES failed with flag -10 at 3.11 ms.
- Everything below is a post hoc redesign of the step, judged on the same windows and limits but not pre-registered.

**Two defects, found and fixed** (diagnostic runs in `~/crucible_runs/2026-11-04/i5diag/`).
1. **The stopping rule saw only the four bulk defects** (rho, the two momenta, E).
   - So rho k and rho omega could be left unconverged, and the final update (the known part plus a h R at the last iterate) took rho omega below 0. In the repro d, cell 1702 at 2.49 ms: rho omega -58 against the iterate's 111.
   - Fix: the turbulence defect enters the stopping rule, its scale taken over the cells the wall functions do not prescribe. The final rho k and rho omega are the last iterate's, which is admissible by construction (they carry no budget).
   - With this alone, f5 (5 μs) still failed in CVODES after 3 ms, and f10 and g10 (10 μs) at 0.5 ms with about 131 inner iterations per step.
2. **A fixed 5 to 10 μs step does not resolve the start-up**, where the cold H2 front runs down the nozzle.
   - The coldest base cell, checked after every step to 0.6 ms: explicit (x) 22.5 K; implicit at a fixed 10 μs (h10) 0.57 K at 0.48 ms. In g10, a cell at 0.47 K reached CVODES, whose rates threw (flag -10).
   - Fix: **time-error step control.**
     - The error compares the step's increment with BDF2's predictor (U^n plus h / h_prev times the previous increment).
     - It is taken per cell, on the start's own scales: rho; |m| against sqrt(|m|^2 + rho θ); θ, the internal energy above the floor with k. Then an RMS over the cells.
     - The next step is h × clamp(0.9 sqrt(tol / err), 0.2, 2). A step over the tolerance is redone at that step. The first step runs at the explicit stable step.
   - With it, a3 (tolerance 3e-3) and a10 (1e-2) end 0.6 ms at a coldest 22.4 K and 22.3 K, against x's 22.5 K. a10 took 351 steps to 0.6 ms at 5.0 inner iterations each; a3 took 624 steps at 3.0.

**The agreement runs after the fixes.** i10t and i5t: the step at most 10 and 5 μs, time tolerance 1e-2, otherwise 7.2's settings (binary `agree/bin_t`, md5 in `bin_t.md5`, its source as `build_t.patch` on 07a332f).

**The agreement so far** (e1 at 14.5 ms at 00:41 on 5 November; window L judged on [0.5, 14.5] ms, so the verdict holds on that span only). Each cell: RMS / largest |difference| over e1's RMS magnitude on the window, in %.

| Run | mdot | F_vac | p0 | p_face | q_wall | q_face (reported) | Window L so far |
|---|---|---|---|---|---|---|---|
| i10t | 0.26 / **3.79** | 0.025 / 0.20 | 0.012 / 0.21 | 0.012 / 0.084 | 0.011 / 0.044 | 0.22 / 1.88 | **FAIL** (mdot's largest) |
| i5t | 0.22 / **2.43** | 0.025 / 0.18 | 0.012 / 0.21 | 0.012 / 0.082 | 0.011 / 0.042 | 0.20 / 1.89 | **FAIL** (mdot's largest) |

- Every RMS is within 0.5%, and every largest but mdot's is within 2%.
- Where it fails: the outlet surge at 1.6 ms, when mdot rises from 6.7 to 10.4 kg/s in about 50 μs. To 12.8 ms, 8 of i10t's samples and 6 of i5t's were over 2%.
- Also judged, all met: gas budgets at most 2.2e-13 (e1), 1.3e-13 (i10t) and 1.2e-13 (i5t), against 1e-11; non-finite samples only at the start (3 or 4, p0 and c* before the end column flows), none in window L. e2 is identical to e1 on all 1,379 common samples before t_L.

**The time tolerance is the knob** (a sweep on the same binary, 0 to 3 ms, the step at most 10 μs; window [0.5, 3] ms, scale e1's RMS over that window, so stricter than window L's):

| Time tolerance | mdot RMS / largest | F_vac largest | p0 largest | Wall on [0.5, 3] ms | Speed-up there |
|---|---|---|---|---|---|
| 1e-2 (i10t) | 0.43 / 2.69 | 0.22 | 0.28 | 95 s | 2.0 |
| 3e-3 (u3e-3) | 0.16 / 1.11 | 0.07 | 0.09 | 119 s | 1.6 |
| 1e-3 (u1e-3) | 0.078 / 0.50 | 0.03 | 0.03 | 130 s | 1.45 |

- e1 took 189 s on the same window. The error falls about as the tolerance to the power 0.73 (derived from both pairs).
- So 3e-3 is the candidate for window L. Its cost is in the start-up only: later in the lead the step sits at its cap, where the time error is far below any of these tolerances (i10t at 14.5 ms: 3e-6).

**Speed-up** (measured on the walls; the Mac shared by four runs, one thread each; the walls include the p0 instrument, below):

| Span | e1 | i10t | i5t |
|---|---|---|---|
| [0.5, 10.5] ms | 577 s | 212 s (**2.72**) | 221 s (**2.61**) |
| [3, 10.5] ms | 388 s | 117 s (**3.32**) | 125 s (**3.10**) |
| [0.5, 3] ms | 189 s | 95 s (1.99) | 96 s (1.97) |

- To 10.5 ms, e1 took 35,247 base steps (0.30 μs mean) and i10t 1,879 (5.6 μs mean, derived): 18.8 times fewer, at 9.0 inner iterations per step. So an inner iteration costs about 0.77 explicit steps here (derived from the walls).
- Against the prototype's 5 to 6: the start-up's small steps, and the viscous and SST terms' extra iterations (inferred, not separated).
- At 14.5 ms, i10t's step sits at the 10 μs cap with a time error of 3e-6 against 1e-2. **The cap binds, not the error**, so a larger cap is the next lever. It is untested, and window L would have to be re-judged at it.

**The p0 instrument now dominates every run.**
- From about 11 ms in e1, and 12.5 to 13 ms in i5t and i10t, the wall per 0.5 ms grew from about 25 s to 856 s (e1) and from about 7 s to 638 s (i10t).
- Stack samples (macOS `sample`, 3 s each, 23:41 to 23:44): 100% of the time of e1, e2 and i10t, and 84% of i5t's, was in Cantera's `equilibrate`, called from p0's `stagnationPressure` in the driver (a TP solve, then up to 50 SP solves per cell of the end column). ChemEquil's Brinkley estimate took most of it.
- The gas there is below 200 K: Cantera warns about 600 times per 0.5 ms that T is outside its 200 to 6000 K range.
- So past 11 ms the march's speed-up cannot be read from the walls. The driver now times p0 and counts its SP solves ("p0's share" on the progress line); the probe from e2's 13.5 ms checkpoint waits on the Mac's RAM floor.
- The onset differs by about 2 ms of simulated time between the schemes, though their p0 agree to 0.2%. Inferred, not checked: a cell of the end column cooling through a threshold that the two reach at slightly different times.
- Section 4's figure (samples 0.7 to 2.5% of the lead's wall, measured on the AMR lead, whose steps cost far more) does not hold for the base alone in the cold lead.
- **Consequence:** at this pace e1 does not reach t_L by 05:30 (derived: about 1,200 s per 0.5 ms at 00:41, with 11 intervals left), and window H is out of reach tonight.

**Verdict so far.**
- The pre-registered test failed: i5 and i10 did not run.
- The redesigned step meets window L's RMS bars and every largest bar except mdot's at 1.6 ms, at time tolerance 1e-2. A tolerance of 3e-3 meets that bar on [0.5, 3] ms. It needs its own run over the whole window before it counts.
- The measured speed-up is 2.6 to 3.3 in the lead, not the prototype's 5 to 6.

### 7.4 The p0 fix, and the agreement set on backhouse (5 November; all numbers measured unless labelled)

**The p0 fix (7eb8956).** Gas with no oxygen and no H atoms below 600 K takes p0 frozen. The reasoning is in the code's comment (`tests/rl10_tier1.cpp`).
- Selftest P0 (backhouse, gcc 13): frozen against equilibrium to 1.61e-12 on seven states, limit 1e-10, each on its expected path. The whole selftest: 0 checks failed.
- **The march is unchanged.** e1o (28865e3's driver, equilibrium p0) against e1 (7eb8956) on backhouse, 1,099 common samples to 10.99 ms:
  - the 10 march columns are identical as text (17 digits);
  - p0 and c* agree to 1.54e-10 relative, with no NaN mismatches.
  - The 1.54e-10 is over P0's 1e-10, which was set on the selftest's seven states. Inferred, not checked: it is the equilibrium path's own SP solve on gas below 200 K, outside Cantera's 200 to 6000 K fits. Either way it is 3e7 times below window L's 0.5% bar (derived).
- **The cost removed.** p0's share of the wall in e1o was 62.6% at 0.5 ms, 55.0% at 1 ms, 45.2% at 2 ms and 31.0% at 11 ms. e1o took 1,290 s to 11 ms against e1's 899 s (1.43 times, derived). In every 7eb8956 run the share is 0.0 to 0.1% to t_L.
- On backhouse the old instrument was costly from the start, not only past 11 ms as on the Mac. Inferred: a different Cantera build. Not checked.

**The set.**
- Binary `t1_agree`, 7eb8956, md5 e23d649afa924f7bf2fc6a6642bacce0, Release. Table B is `rl10_s.bin`, md5 equal to the Mac's.
- Base 149x16 alone, one thread each, nice 10, all launched together at 00:58:52 at load 13. The load was 20 to 25 while they ran (these runs plus the box's other jobs).
- The runs: e1, e2 and i10t, i5t as in 7.3, plus **i10u** (step at most 10 μs, time tolerance 3e-3) and i20u and i40u (at most 20 and 40 μs, 3e-3).
- Directory `/home/greff/crucible_p0/agree/`; Mac copy `~/crucible_runs/2026-11-05/bh_agree/`. Computed by `tools/agreement_tier1.py`.

**Window L, the whole lead [0.5, 20] ms** (1,951 samples each). Each cell: RMS / largest |difference| over e1's RMS magnitude on the window, in %.

| Run | mdot | F_vac | p0 | p_face | q_wall | q_face (reported) | Window L |
|---|---|---|---|---|---|---|---|
| i10t (10 μs, 1e-2) | 0.23 / **3.87** | 0.021 / 0.19 | 0.010 / 0.20 | 0.010 / 0.082 | 0.010 / 0.042 | 0.22 / 2.11 | **FAIL** (mdot's largest) |
| i5t (5 μs, 1e-2) | 0.19 / **2.48** | 0.021 / 0.18 | 0.010 / 0.20 | 0.010 / 0.079 | 0.009 / 0.040 | 0.19 / 2.12 | **FAIL** (mdot's largest) |
| **i10u (10 μs, 3e-3)** | 0.093 / 1.60 | 0.007 / 0.058 | 0.003 / 0.066 | 0.005 / 0.027 | 0.006 / 0.016 | 0.14 / 0.63 | **PASS** |

- Also judged, all met: gas budgets at most 2.65e-13 (e1, e2), 1.16e-13 (i10t), 9.95e-14 (i5t) and 1.28e-13 (i10u), against 1e-11. Non-finite samples only at the start (3 or 4), none in window L. e2 is identical to e1 on all 1,999 samples before t_L.
- **The platform's own difference** (the same binary's source, Mac clang against backhouse gcc; computed as above):
  - i10t on 1,622 samples to 16.71 ms: mdot 0.0075 / 0.076%; every other judged observable at most 0.0085% largest.
  - e1 on 271 samples to 3.2 ms: mdot 0.036 / 0.27%.
  - So mdot's largest moves by up to about 0.3% between machines, a fifth of i10u's 1.6% (derived).

**Speed-up**, e1's wall over i10u's, from the progress lines' cumulative walls. The runs were concurrent, one thread each, and the scheduler's placement on P and E cores is not controlled. e1 and e2, doing identical work, took 1,582 and 1,577 s to 20 ms.

| Span | e1 | i10u | Speed-up | Note |
|---|---|---|---|---|
| [0.5, 20] ms | 1,548 s | 217 s | **7.13** | the window L span (i10t 7.37, i5t 4.73) |
| [0.5, 3] ms | 232 s | 66 s | 3.5 (derived) | the start-up, where the time error sets the step |
| [3, 10] ms | 558 s | 83 s | 6.7 (derived) | |
| [10, 20] ms | 758 s | 68 s | **11.1** (derived) | the step at its 10 μs cap |

- Steps to 20 ms:
  - e1: 65,642 base steps (0.30 μs mean, derived).
  - i10u: 3,599 steps (89 redone on the time error) at 6.16 inner iterations each; 9 steps were left unconverged at the 400-iteration limit.
  - i10t: 2,785 steps at 8.26 iterations; i5t: 4,444 at 7.70.
- At 20 ms i10u's time error was 2.9e-6 against its 3e-3 tolerance, with the step at 10 μs. **The cap binds, and the cap is the driver's 10 μs sample interval**, since every step ends on a sample (below).
- On the error alone, h scales as (tol / err)^(1/2), so the step could grow about 30 times there (derived). The inner iterations would grow with it, by an unmeasured amount.
- Against 7.3's Mac figures (2.6 to 3.3): those walls carried the old p0 instrument from the start.

**A larger cap is not a lever as tier 1 stands.**
- i20u and i40u are identical to each other in every byte of their series, and their steps never exceeded 10 μs (the last: 1e-5 s).
- They leave i10u's series at 3.34 ms by 1e-14 relative. Inferred: a step ending an ulp apart, at min(the cap, the distance to the next sample).
- The difference grows to 1e-5 relative by 4 ms, and back to about 3e-7 by 8.8 ms.
- To use a step over 10 μs, the samples must come from the step's own history (BDF2's interpolant) rather than end the step. That is engineering, and not built.

**A latent defect, found by i20u and i40u: rho omega's collapse.** Both failed at 8.975 ms. The error was "Could not advance an admissible gas state after 14 implicit step reductions", from the first half source at cell 2073, at every step down to 6.1e-10 s.
- A probe rerun (i20p: cell 2073's turbulence after every accepted step, a local patch, not committed) gives the chain. The cell is thin, cold nozzle gas: rho 1.33e-3 kg/m^3, T about 7 K, k about 4,200 m^2/s^2.
  1. From 8.95 ms, rho omega fell over three 10 μs steps from 670 to 288, 81 and 5.4. Two of those steps ended unconverged at 400 inner iterations. The time error does not see k or omega, so the step stayed at 10 μs.
  2. At 8.975 ms (a 5 μs step to a sample, c = 0.125), BDF2's history alone asked for rho omega = 5.38 - 0.125 × 73.45 = -3.81, the known part.
  3. The inner damping (rho omega at most halved per iteration) then halved it on each of 61 iterations. The second half source's cross-diffusion destruction, |CD| / omega (singular at omega = 0, so the ODE itself reaches 0 in finite time), took it to 1.7e-118. That state is admissible, so the step was accepted.
  4. In the next step's first half source, the MPRK22 stage underflows to omega = 0, and its ratio (d0 y0 + d1 y1) / y1 is 0/0. The NaN does not depend on the step.
- **Why only i20u:** i10u's trajectory differs from i20u's by 1e-14 at 3.34 ms and does not fail. So the failure is marginal: a rounding-level perturbation decides it, and passing i10t, i5t and i10u does not show the step is robust.
- **The explicit path does not show it.** A per-sample probe of e1's trajectory (e1p) (`t1_probe2`, a local read-only patch; its series identical to e1's to 9.1 ms, measured) gives, over 0 to 9.1 ms, a smallest omega over all cells and samples of 77.0 1/s (at 1.46 ms), and cell 2073 steady at rho k 5.5 and rho omega 29 (k 4,140 m^2/s^2, omega 2.19e4 1/s) from 8.9 to 9.0 ms.
- **The cause: one unconverged step, accepted.** i20p against e1p in cell 2073:
  - The two agree to 2% in rho k and 6% in rho omega up to 8.88 ms (measured, every sample).
  - At 8.89 ms one step ended at the 400-iteration limit, unconverged, and was taken: rho k 5.59 to 0.115 and rho omega 31.5 to 652, where e1p holds 5.51 and 29.1. T stayed at 7.0 K, so the time error, which does not see k or omega, passed it.
  - The next steps pulled it back (rho omega 667, 1,123, 1,518, 1,207, then 670 at 8.95 ms), and the fast fall from there is the chain above.
- **The fix: with the time error on, an unconverged step is redone at half the step** (the open item below, closed). Without the time error (the caller's fixed step: case P, the prototype) the step is taken as before. The run line counts the redone steps as "unconverged (redone)".
  - Measured (backhouse; `t1_r`, md5 e516d90d; i20r and i10r at i20u's and i10u's settings, to 20 ms): both reach 20 ms. i20r redid 10 unconverged steps and i10r 7; every other rejection was on the time error (92 each), so **the guard below never fired.**
  - Window L: **PASS** on both, every judged figure again equal to i10u's to the digits printed.
  - The inner iterations fell to 5.26 (i20r) and 5.53 (i10r) per step from 6.24 to 6.53, the 400-iteration steps gone (the most now 240 and 387). Speed-up on [0.5, 20] ms: **7.51** (i10r, 206 s) and 7.44 (i20r, 208 s) against e1's 1,548 s (measured; concurrent with three other runs, so the walls carry load as above).
  - Tier 1's selftest with it: 0 checks failed (backhouse). The rule acts only with the time error on, which only tier 1 sets (derived from the code), so the core and turbulence ctests do not reach it; they were rerun on the final sources (below).
- **The guard (committed with this section):** an attempt whose history alone puts rho omega at or below 0, in a cell the wall functions do not prescribe, is redone at half the step. Halving the step shrinks c = w^2 / (1 + 2w) at least twofold (derived).
  - Measured (backhouse, gcc 13; `t1_fix`, md5 78f7cde5; i20f and i10f at i20u's and i10u's settings, to 20 ms): **i20f passes 8.975 ms and reaches 20 ms.**
  - The guard redid 8 attempts in i20f and 2 in i10f (the rejections not on the time error: 97 and 91, against 89 on it in each). So it also fires where the old step recovered: i10f leaves i10u's series at 5.75 ms, i20f leaves i20u's at 8.93 ms.
  - Window L on both: **PASS**, every judged figure equal to i10u's to the digits printed (mdot RMS 0.0933%, largest 1.598%; the rest as i10u's row above).
  - Speed-up on [0.5, 20] ms: 6.76 (i10f, 229 s) and 6.53 (i20f, 237 s) against e1's 1,548 s (measured). i10f took 3,601 steps at 6.24 iterations against i10u's 3,599 at 6.16, 1.4% more inner iterations (derived); the rest of its 5.5% longer wall is the different concurrent load (inferred: i10f ran beside i20f and one explicit probe, i10u beside e1, e2 and four implicit runs).
  - With the guard: tier 1's selftest 0 checks failed, and the core and turbulence ctests passed (backhouse, 01:47).
- Still open: the time error does not include k and omega. With unconverged steps redone, it is the remaining way a wrong turbulence state can be taken.

**Window H needs the p0 instrument fixed first.**
- Past t_L the LOX brings oxygen into p0's end column. So the instrument is back on Cantera's equilibrium path for gas below 200 K (the warnings are those of 7.3). Every run's progress stopped there: none printed a line past 20.000 ms in over 20 minutes.
- Whatever the solver, the hot window needs the instrument fixed first. Taking oxygen-bearing gas frozen does not follow from the lead's argument: unburned H2 and O2 at equilibrium is mostly water, so the two paths give different p0 there.
- For cold gas the frozen value is the one a probe would read: a stagnation process lasts microseconds, and H2 with O2 does not react below 600 K in that time (inferred; no ignition delay was computed). So p0 is now frozen below 600 K whatever the composition, with its own check (7.5). That is engineering on a development instrument, not a call for Ben.
- At about 01:30 backhouse's distro restarted (`tools/resume/README.md`) and killed every run; the runner relaunched e1, e2 and i10u from 0, since none had a checkpoint yet. They were stopped at once (exit 143): past t_L they would stall in p0 again. Their series to 20 ms are on the Mac. i10t and i5t had been stopped at 20 ms (exit 143), having failed window L.

**Verdict (7.2's windows, judged post hoc on the redesigned step, as 7.3 says).**
- **Window L passes at a 10 μs cap with time tolerance 3e-3 (i10u), at a measured speed-up of 7.1 over the whole lead and 11 in its second half; the committed step (unconverged steps redone, i10r) passes at 7.5.** The tolerance 1e-2 fails on mdot's largest at the 1.6 ms surge.
- That is above the inviscid prototype's 5 to 6, because the cap binds at an error a thousand times below the tolerance (derived: 2.9e-6 against 3e-3) and the driver's samples, not the method, set the step.
- Window H is not judged yet: the p0 instrument stalled on oxygen-bearing cold gas. Its fix and the hot-window runs are 7.5.
- The step is not yet robust. One perturbed trajectory in three at the same settings failed at 8.975 ms, by omega's collapse. Its cause was an unconverged step taken. With unconverged steps redone (and the guard kept as a backstop), i10r and i20r pass 8.975 ms and window L at 7.5 times e1's speed, and the guard never fires (measured). Two trajectories are not a proof: the time error still does not see k or omega.

### 7.5 Past t_L: the p0 instrument, cold gas with oxygen, and the hot window (5 November; backhouse, gcc 13; all numbers measured unless labelled)

**p0 frozen below 600 K whatever the composition.** The rule of 7.4: the frozen path whenever the cell and the frozen T0 are both below 600 K.
- Selftest P0: frozen against equilibrium to 1.27e-12 on the lead's gas (limit 1e-10), every state on its expected path (7 frozen, 2 past 600 K on equilibrium).
- On cold H2 with O2 the frozen value lies 4.78e-4 from the equilibrium one (reported, not judged). Equilibrium there is mostly water, which a stagnation process of microseconds does not make (7.4).
- It was not enough. i10h (a 10 μs cap and time tolerance 3e-3, as i10u; `t1_h0`, this rule and the guard; from 0) passed t_L and then advanced 0.3 ms in about 8 minutes of wall, against i10k's 222 s for the whole 20 ms lead.

**The stall: ChemEquil on a trace element.**
- A timing probe on i10h's state (`t1_p0t`, a local patch) found one end-column cell taking 40.0 s over the loop's 50 SP solves. The cell: 633 K, so the equilibrium path; 736 m/s; H2 0.615 and O2 0.385 by mass; and N2 at 7.6e-70, the N2 fill washed out.
- A standalone benchmark on that cell (Cantera 3.2, ChemEquil, the "auto" solver), one TP equilibrate:

| N2 mass fraction | 7.6e-70 | 1e-30 | 1e-20 | none |
|---|---|---|---|---|
| ChemEquil | 41.3 s | 4.6 ms | 3.7 ms | 2.4 ms |

  The gibbs solver took 0.1 ms on the same cell.
- ChemEquil's SP solve also has a floor. From the fourth solve the loop's step stayed at 1.368e-9, so the loop ran all 50 solves and drifted 6.3e-8 (derived: 332,538.429 to 332,538.409 Pa). It ended 6.8e-8 below the converged value that gibbs and vcs agree on (332,538.431 Pa, in 5 solves).
- Why not use gibbs throughout: on ten hot states ChemEquil and gibbs agree to 2.1e-7, but gibbs took 21.6 ms against ChemEquil's 0.2 ms at 1,500 K.
- **The fix** (`tests/rl10_tier1.cpp`, stagnationPressure):
  - species under a mass fraction of 1e-20 are dropped before equilibrating. A species enters h and s in proportion to its mass fraction, so this moves p0 by of order 1e-19 (derived).
  - The loop stops at a step of 1e-13, or once a step under 1e-6 fails to halve.
- A second P0 check:
  - the washed-out cell against the N2-free one: 0 difference, in 0.3 ms (limit 1e-12);
  - seven hot states (1e5 to 3e6 Pa; O/F 2 to 8; burnt, half burnt, and unburnt at 900 K; 100 to 700 m/s) against gibbs: 1.57e-7 (limit 1e-6), at 8.4 SP solves a state for the two solvers together.
  - The whole selftest: 0 checks failed (`t1_g`, md5 c6defe3a).

**CVODES on cold gas with oxygen.**
- i10k (i10r's settings and binary, from 0, then `t1_g` from its 20.3 ms checkpoint) failed at 20.35 ms: CVODES flag -4 at cell 2305, T 2.40 K, rho 6.77e-4 kg/m^3, dt 1.77 μs.
- This is nozzle gas expanded far below the mechanisms' 200 K fits, which the LOX's O2 has reached. The skip integrates it because it holds both H2 and O2.
- **The cold limit** (`ReactingFlow::setColdLimit`; off by default; tier 1 sets it at kColdLimit, 200 K, and CRUCIBLE_CHEM_COLD=0 turns it off): every cell below the limit is left unchanged, whatever its species.
  - This is a model, not exact: frozen chemistry below the fits. Inferred, not computed: below 200 K the rates with an activation barrier are nil, and the barrierless recombinations have only the radicals' traces to act on.
  - It acts after the skip, so a cell it counts is one that would otherwise have been integrated.
- Check 6 in `tests/chemistry_skip_tests.cpp`, pre-registered in its header before it was built: 4's duct with columns 14 and 15 (H2 0.2, O2 0.1, N2 by moles) at 150 K, with the limit and without.
  - Counts: skipped 56, cold 8 and integrated 8 with the limit; 48, 0 and 16 without, as stated.
  - Columns 14 and 15 are unchanged bit for bit with the limit, and columns 0 to 13 are identical bit for bit with and without.
  - Both mechanisms pass.
- **Probe i10c** (`t1_c`, md5 86aa8b02, from i10k's 20.3 ms checkpoint, to 20.8 ms): it passes 20.35 ms and reaches 20.8 ms. The limit froze 333,637 cell updates on the way, against 1,010,006 integrated (both counts over the whole run).
- **The limit also acts in the lead, on the implicit path only.**
  - i10k rerun from 0 on `t1_c` counts 93,378 cold updates by 0.5 ms and about 40,000 per ms through the lead (831,233 by 18 ms).
  - e1h and e2h, restarted on `t1_c` from 16 ms, count 0 to 18.7 ms.
  - Inferred, not checked: the implicit path's cold nozzle cells carry H above the skip's 1e-14 trace, which the explicit path's do not.
  - So on the implicit path the limit changes the lead's trajectory, and window L is rejudged on it below.

**Window L on the cold limit (i10k, `t1_c`, from 0; Mac copy `~/crucible_runs/2026-11-05/bh_agree/i10k/`).**
- **PASS**, every judged figure equal to i10r's to the digits printed except mdot's RMS (0.0930% against 0.0933%) and q_wall's RMS (0.0063% against 0.0064%). The largest stay mdot 1.598%, F_vac 0.058%, p0 0.066%, p_face 0.027% and q_wall 0.016%. Gas budget at most 1.16e-13.
- Against i10r directly (the same step and settings; the flow differs only by the limit, since the p0 changes touch the instrument alone): the series part at 60 μs. The largest difference over window L, as a share of i10r's RMS, is 2.9e-4 on mdot, at most 5.3e-5 on the other judged observables and 2.0e-3 on q_face (reported).
- So the limit's effect on the lead (0.029% on mdot, derived) is below the platform's own (0.076% on i10t, 7.4) and 70 times below window L's 2% bar (derived).
- Speed-up on [0.5, 20] ms: 6.82 (i10k 227 s against e1's 1,548 s). i10k took 3,622 steps (107 redone) to 20 ms, against 3,607 (99) for its first attempt without the limit (now `agree/i10k_nc`), so the limit saves no work in the lead.
- With the limit and check 6: the core, turbulence and chemistry_skip ctests passed (3 of 3, 02:33), and tier 1's selftest failed 0 checks (`t1_c`).

### 7.6 The hot window: every run fails about 1 ms past t_L, on the spray's sub-step (5 November; backhouse, gcc 13; all numbers measured unless labelled)

**The failure.** e1h, e2h (explicit, seeds 1 and 2) and i10k (implicit) all run window H on `t1_c` from their lead checkpoints. All three break within about 1 ms of t_L:
- the liquid's energy budget, (b), is 4e-14 at 20.5 ms and 1.7e-8 to 5.0e-7 at 21.0 ms;
- e1h's last sample is at 21.03 ms and i10k's at 21.19 ms, both stopped by "Spray: the heat-transfer number's residual is NaN" (the guard added for the 23 October NaN, `face_inlet/leak_2026-10-23/spray_nan/`, reached here by another route);
- e2h was stopped by hand at 03:14 (its last sample 21.14 ms; (b) 1.7e-8 at 21.0 ms), since it was on the same path.

**Two probes, from i10k's 20.81 ms checkpoint** (print-only copies of the code; switches outside the CRUCIBLE_* environment, so the checkpoint's signature holds):
- Probe 1 prints the three liquid budgets after every base step. The energy budget is 2.7e-14 at 20.810 ms, then jumps to 5.23e-9 in the one step ending at 20.8378 ms, and climbs in steps from there; mass and momentum stay at about 5e-14.
- Probe 2 prints a sub-step whose corrector leaves the enthalpy below the LOX table. The first: a 1.54 μm drop at 93.89 K, condensation clamped (so the evaporation is zero at stage 1), in gas at 183.93 K and 579,378 Pa with Y_O2 0.8545 (T_sat 110.94 K); sub-step 0.773 μs. Its corrector puts h at -4.694e5 J/kg against the table's floor of -4.649e5 (54.361 K, the triple point), and 43% of the drop's mass goes in that one sub-step.

**The mechanism** (the drop's state measured, the reading of it from the code):
- Stage 1 sees a condensing drop, so its stiffness limit counts the conductance alone and allows a heating of up to a quarter of the gas-to-drop difference, about 22 K here (derived).
- The predictor lands just short of saturation. There the surface's O2 mole fraction p_sat/p approaches 1, and the Spalding number (Y_s - Y_inf)/(1 - Y_s) grows without bound, since the other gas is H2, so Y_s runs to 1 faster still.
- Stage 2's evaporation is then far stiffer than stage 1 allowed. The corrector takes h below the table, and `LiquidOxygen::temperature` clamps silently at 54.361 K: the parcel's energy is recomputed from the clamped T, which creates energy (the budget's jump).
- A later corrector on such a drop takes s = m^(2/3) below zero: the mass is NaN, so the Spalding residual is NaN.

**The fix** (core/spray.cpp, `Spray::moveParcel`):
- After stage 2, its own limits (a quarter of the mass time, and a quarter of the thermal time when not boiling, both at the predictor) are computed.
- A sub-step over twice them is retaken at them, at most 30 times. The factor 2 leaves every sub-step that was already smooth exactly as it was.
- A corrector that takes the mass below zero now throws by name. Bound (derived): with both limits met, s falls by at most a quarter in a sub-step, so this cannot fire unless the 30 retakes run out.

**Check F3, pre-registered** (written 5 November at about 03:15 CST and committed in 95f0ed3, before it was built or run; in `tests/spray_tests.cpp`):
- The probed drop alone, in a chamber at rest (4 x 2 cells, 4 mm by 2 mm), filled with O2 0.27 and H2 0.73 by moles (Y_O2 0.8545, derived) at 579,378.25 Pa and 183.932 K.
  - The drop: d 1.544204 μm, T 93.894757 K, at rest, count 1. One-way coupling; breakup and dispersion off.
  - At rest its Nusselt number is 2, not the run's, so the step that lands the predictor just short of saturation differs from the probe's 0.773 μs. Hence a sweep: one `advance(dt)` for each dt from 0.10 to 1.70 μs by 0.02 μs (81 runs), each run a fresh spray on the same flow.
- Pass, on every run:
  - nothing throws;
  - the drop either vaporised, or its T is finite, at least its start (it is heated, and the clamped condensation adds nothing) and above the table's floor;
  - the liquid energy budget (`balanceErrors()[1]`) is at most 1e-12 in magnitude.
  - Reported: the worst budget, the coldest end T, and the most sub-steps in a run.
- Expected: FAIL on the unfixed spray.cpp (some dt ends on the floor, or throws), PASS with the fix. The other spray checks keep their verdicts.

**The fix in tier 1, pre-registered** (written 5 November at about 03:23 CST, before any of it was seen). `t1_f` (md5 1ef89249) is `t1_c`'s source plus the fix (the two trees differ in `core/spray.cpp` and `tests/spray_tests.cpp` alone, by `diff -r`). e1s, e2s and i10s are e1h's, e2h's and i10k's settings on `t1_f`, from 0 to 30 ms, launched at 03:22 under the resume runner. Every surviving checkpoint of the failed runs is past t_L, so none is reused.
- (a) Before t_L there are no parcels, so i10s's series to 20 ms is identical to i10k's, and e1s's and e2s's to e1h's and e2h's, bit for bit, comparing columns as printed.
- (b) i10s passes 21.19 ms (i10k's last sample) and e1s passes 21.03 ms (e1h's), with the liquid budget (b) under the run's own 1e-11 at every progress line, and no throw.
- (c) Window H then goes to the agreement test of 7.2 as it stands, unchanged.
