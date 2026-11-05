# Speed plan: from 46 to 75 days to under an hour per RL10 run (4 November 2026)

Written for the roadmap's AT RISK inquiry box. Ben's bar (24 October): a full RL10 run in under an hour, by a GPU port, and the RL10 is one engine of 20 or more. Principle 1 (always time-dependent) is canon, so every lever here marches the flow in physical time. SPEED.md holds the measurements this plan builds on; this file does the arithmetic and the order.

Labels: **measured** (a run or a profile), **derived** (arithmetic on measured numbers, shown), **guessed** (an estimate with its basis stated), **spec** (a manufacturer's figure).

## 1. Verdict

- No CPU lever, nor any product of CPU levers, reaches the bar. The CPU levers that exist multiply to about 40 to 110 (derived, section 3). The gap is 1,100 to 1,800 (derived, section 2).
- A GPU port of the gas step is necessary. On its own it is not enough, guessing its rate on backhouse's card at 50 to 400 times backhouse's 4 threads.
- Under the guess, the bar is reachable as the product **GPU port × D (two finer levels, not three) × B (CFL 0.8)**. That product is 700 to 7,200 (derived from the guess). Two conditions:
  - the chemistry's CVODES fallback and the spray both leave the CPU's critical path (section 4);
  - D and B pass their accuracy tests.
- The implicit step (I) is the hedge if D fails or the GPU's rate lands low.
  - Its measured 5 to 6 on the CPU is for the inviscid prototype only.
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
| **B. CFL 0.4 to 0.8** | **2** | derived; the 1 ms test's outputs were lost (section 5) | rerun the test (hours on backhouse) | to rerun |
| **A. The chemistry skip** | **1.5 to 2.2** in the lead's wall | measured (built) | done | done, in the basis's successors |
| **Tables for every chemistry state** (no CVODES in the march) | lead: **1.5 to 2.2** more on the CPU (the fallback is still 33 to 55% of the lead's wall after A, derived from A's 1.5 to 2.2 over a 69 to 71% share); on the GPU path it is **necessary**, not optional (section 4) | derived | widen Table B's f_H range to the lead's pure-H2 and mixing states, or a second table | not started |
| **Spray on the GPU, or fewer parcels** | removes the 1.5 to 1.7 multiplier on 85 ms; **necessary** on the GPU path (section 4) | measured multiplier | a device parcel kernel, or a parcel-count study (kParcelRate 5e5) | not started |
| **I. Implicit time-accurate step** (BDF2, dual time, LU-SGS) | **5 to 6** measured on the inviscid 1-level prototype, one thread; bounded near 6.6 by LU-SGS's inner iterations | measured (IMPLICIT_PROTOTYPE.md) | tonight: the port to the viscous, SST and wall-function step on tier 1's base mesh (section 7) | in progress |
| **C. Less simulated time** | 1.7 to 2.6 against the old schedule; the approved central case (lead plus 85 ms) already takes most of it | inferred | the Director's schedule | taken |
| **G. The flow's per-update cost** | 1.3 to 2 | guessed | profiling and tuning; folded into the port | not started |
| **F. Threads** | 2.2 at 4 threads on backhouse, 3.1 on the Mac; 4.4 at 8 on the Mac | measured | the slot policy | in use |
| **E. The machine** | backhouse at 4 threads runs 2.5 times the Mac's wall | measured | none | known |
| **H. IMEX acoustics** | ceiling 1.4 to 2.2 | measured ceiling | 3 to 6 weeks | dominated by I |
| **J. A steady solve** | out | principle 1 | | out |
| **K. Local time stepping** | at most 1.8 to 4 | measured bound | | rejected by Ben |

**The CPU-only product:** D × B × the tables × G × F (8 threads against 4, Mac) is 7 to 9 × 2 × 1.5 to 2.2 × 1.3 to 2 × 1.4 = **38 to 110** (derived, before I). The run time comes out at 10 to 47 h. With I's 5 to 6 on top, if the gains compose (unmeasured), it is 1.7 to 9.5 h. Both figures assume the spray's cost falls with the gas's, which is not measured. The bar is not met, and this path assumes backhouse's whole CPU, which is shared.

**The GPU product:** GPU × D × B = 50 to 400 × 7 to 9 × 2 = **700 to 7,200** (derived on the guessed GPU rate), against 1,100 to 1,800 needed.
- At the low end of the guess, the bar is missed by up to 2.6 times.
- Mixed precision or I would close the low end. Mixed precision means FP32 fluxes with FP64 sums, which keeps conservation exact to the sums' rounding.

## 4. What caps the GPU path (Amdahl)

- **The chemistry fallback.** After A, CVODES is 33 to 55% of the lead's wall (derived above). If the gas went infinitely fast and CVODES stayed on the CPU, the lead would gain at most 1.8 to 3 times. So every state the march meets must be tabled, or integrated on the device.
  - In the hot phase, Table A's chemistry is already 0.08 of 1.96 μs per update (measured, SPEED.md 1), so it is cheap.
- **The spray.** It adds 0.5 to 0.7 of the gas time over 85 ms (measured multiplier). Left on the CPU, it caps the run's gain at 2.4 to 3 times (derived: 1.5 to 1.7 over 0.5 to 0.7).
- **The driver's samples** (c* by Cantera's equilibrate every 10 μs) were 0.7 to 2.5% of the lead's samples (measured). At a 1,000-times gas speed-up they would be most of the run, so the sample becomes a table lookup or a sparser sample.
- **AMR's bookkeeping** (fill, reflux, regrid) is under 1% on the CPU (measured). On the device it must stay there, so it moves too.

## 5. Ranked order (by factor, certainty and cost)

1. **GPU port of the gas step, one level, FP64**, with a measured updates-per-second figure on backhouse's card in its first week. That figure decides everything below. It is the only lever that can reach the bar.
2. **D, two finer levels**: the largest factor that costs only runs. Its accuracy test should run now, on the CPU, while the port is built.
3. **B, CFL 0.8**: rerun the lost test.
4. **Tables for every state**: necessary for the port.
5. **The spray on the device**, or a parcel-count study: necessary for the port.
6. **I, the implicit step**: the hedge. It is worth building only if the measured GPU rate or D's test says the product falls short. On the device it needs a parallel smoother in place of LU-SGS's sweeps.
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
