# The implicit time-accurate step: a bounded prototype, criteria before building

Written 24 October 2026, about 01:50 CDT, before any prototype code exists. Asked for by the Director at 01:40: "4-8 weeks is too big to commit on an estimate. Prototype it first on a bounded case (1D or the 2D pulse/C1 case), about a week, and measure the real speed-up and accuracy against explicit before deciding. Write the prototype's success criteria before building."

The full step's design and its guessed costs are in SPEED.md section 5a (lever I) and TECHNICAL_PLAN step 16 (parked). This prototype replaces the guesses that decide the build with measurements:
- the inner iterations per physical step;
- the cost of one inner iteration, in explicit RK2 steps;
- the accuracy of the history at the physical steps the dynamics need.

It does not decide anything itself. The decision is the Director's, on the numbers in *What is reported* below.

## What is built (about a week)

- **The spatial residual as a function.** `Flow` gets `residual(state, out)`: the existing operator (HLLC with limited reconstruction, the axisymmetric sources, the boundaries, and transport and SST where the case has them), with no change of arithmetic.
  - Gate: the explicit march through the refactored code is bit-identical to `main` on C1 64x12 to 1 ms, and the fast ctest is unchanged (4 known failures only).
- **One implicit stage solver.** It solves (U − U_known)/(a Δt) + R(U) = 0 by dual time:
  - pseudo-time steps local to each cell;
  - LU-SGS sweeps (Yoon and Jameson) on the Rusanov split of the inviscid Jacobian, with the scalar diagonal plus 1/(a Δt);
  - a damped update that keeps rho, p, T, Y_k, k and omega positive, and a rejected step retried at Δt/2, as the explicit steps do. No clipping.
- **BDF2 in physical time.** a = 2/3 and U_known = (4 U^n − U^(n−1))/3. The first step is backward Euler. Then the final update is formed from the last inner iterate's residual, so the boundary fluxes are booked with the same weights as the cells and the budgets can close to round-off (TECHNICAL_PLAN step 16, *Conservation*).
- **One level, no AMR, one thread** for the implicit march. The explicit comparison is also run on one thread, so the measured factor is the algorithm's. Threading (red-black or block-Jacobi) and AMR's synchronous levels are outside the week (see *Not measured*).
- **Chemistry: Table A only** (local equilibrium, an algebraic closure). There is no reaction split, so BDF2's second history level has no split to cross. The finite-rate question is measured apart (M5 below).
- **If the week has room** (none of these is a criterion):
  - the implicit stages of ESDIRK3 (ARK3(2)4L[2]SA, TECHNICAL_PLAN step 16), which reuse the same stage solver;
  - low-Mach preconditioning in pseudo time only.

## The cases

- **P, the acoustic pulse** (smooth, for the temporal order). The simple-wave pulse of `tests/amr_subcycling_tests.cpp` 3(a) on one level, 192 and 384 columns. The reference is the explicit march at CFL 0.1, checked against CFL 0.05.
- **C, the C1 cold start** ("valves open, ignite", inviscid, Table A, `crucible_chamber_study eqt`). It is the bounded version of the RL10's binding constraint: a sonic throat sets the explicit step.
  - Explicit step on 64x12: about 0.06 μs (*measured*, `c1/eq64_log.txt`: 17,010 steps in the first 1 ms).
  - Its fill constant τ_C1 is about 0.5 ms (*derived*: vacuum thrust within 0.1% of settled from 3.5 to 4 ms, 7 τ). The explicit 64x12 run measures it before any implicit run.
  - **The physical steps are set in steps per fill constant, to match the RL10's:** Δt = τ_C1/25, τ_C1/63 and τ_C1/150. The RL10's 50 μs and 20 μs are τ/25 and τ/63 of its 1.26 ms (SPEED.md 5a). τ_C1/150 is the third step for the convergence.
  - **Meshes:** 64x12 and 128x24 for the accuracy, plus 256x48 to 1 ms for the inner iterations at high stiffness. The ratio of the physical step to the explicit one then spans about 70 to 1,300 (*derived*; the RL10's lead spans 49 to 2,100, SPEED.md 5a).
- **T, the turbulent C1** (the PaSR criterion 4 chamber: viscous, SST, wall-clustered rings, 32x6 and 64x12). This is stage 2, only if C passes within the week. It needs point-implicit SST sources and the wall-normal viscous Jacobian. If it is not reached, the decision is made on P and C, and the turbulence's effect on the inner iterations is stated as unmeasured.

## Criteria (pass or fail, judged before the decision)

1. **Temporal order** (P). The L2 error of p against the reference, at three physical steps halving from an acoustic CFL of 4, gives an observed order between 1.8 and 2.2 at the two finest pairs.
2. **Budgets** (C, every run). Mass, energy and each element close within 1e-11 of the throughput, at every inner tolerance tried.
3. **Inner convergence settled** (C, 64x12, Δt = τ_C1/63; Bijl, Carpenter and Vatsa's test). Run the inner tolerance at 1/10, 1/100 and 1/1,000 of the temporal error. The temporal error is measured as the difference between Δt and Δt/2.
   - Pass: at some tolerance, the observables of criterion 4 differ from the 1/1,000 run's by less than 10% of the temporal error.
   - The tolerance that passes is the one used for criteria 4 and 5, and for every cost measurement.
4. **Accuracy of the history** (C, 64x12 and 128x24, against the explicit march at CFL 0.4 on the same mesh). The observables are the injector pressure (C1's p_inj), the outlet mass flow and the vacuum thrust, each normalised by its settled explicit mean.
   - Over [0.5 ms, 8 ms] at Δt = τ_C1/63, the RMS difference is at most 0.5%, and the largest is at most 2%.
   - The difference falls between τ_C1/25, τ_C1/63 and τ_C1/150, with an observed order of at least 1.
   - The first 0.5 ms (the starting shock and the igniter) is reported but not judged. Its smearing at large steps is expected (SPEED.md 5a, *Convective CFL*).
5. **The settled state** (C). c*, vacuum Isp and vacuum thrust over [6 ms, 8 ms] are within 0.05% of the explicit values at every Δt. This is the same bar as Table A's criterion 3. For a steady inviscid state it tests the converged residual, not the time scheme.
6. **Stable at the RL10's stiffness** (C on 256x48, to 1 ms, Δt = τ_C1/63). The run completes without a rejected step left unrecovered. It reports the inner iterations per step against the explicit-to-implicit step ratio.

If 1, 2, 3 and 5 pass, the prototype's numbers are trusted. If 4 fails, the steps that pass it are found (smaller Δt), and the speed-up is reported at those steps.

## What is reported (the measurements the decision needs)

- **M1:** inner iterations per physical step, against the step ratio and against the stage of the transient (the start, the fill, settled). Taken at criterion 3's tolerance.
- **M2:** the cost of one inner iteration in explicit RK2 steps. Measured by wall time on one thread on the same mesh: (implicit wall per inner iteration) / (explicit wall per step).
- **M3:** the measured speed-up per simulated ms at criterion 4's accuracy. This is the explicit wall to 8 ms over the implicit wall to 8 ms, one thread each, on 64x12 and 128x24.
- **M4:** the projection to the 3-level tier 1 and the lead. Explicit all-cell RK2 sweeps per physical step (measured: about 6.5 per base step, SPEED.md 5a) over M1 × M2 at the matching step ratio, against CFL 0.4 and against CFL 0.8 (lever B's result). It is labelled *derived*, with the parts below stated as unmeasured.
- **M5:** the finite-rate chemistry's cost at the implicit step, with no build. CVODES on the lead's states (`tests/chemistry_cost.cpp`, or the lead's own cells dumped) integrating 20 μs against 0.2 μs per call. The chemistry is 68 to 71% of the lead's wall (SPEED.md section 1), so its scaling with the step is half of the lead's answer.

**Proposed rule for the Director's decision** (a proposal, not a criterion): build the full step if M4 is at least 5 against the explicit march at the CFL it would otherwise run (0.8 if lever B passes); do not build it if M4 is below 2; between 2 and 5, the Director decides with the other levers' results (A, B, D) in hand.

## Not measured by the prototype (stated, so the projection carries them)

- **Turbulence and the walls,** unless T is reached: point-implicit SST, the wall functions' prescription inside the inner loop, and the wall-clustered aspect ratios.
- **AMR:** synchronous levels, and fill and reflux at the implicit step.
- **Threads:** LU-SGS's sweeps are sequential. A red-black or block-Jacobi version loses some convergence per iteration. The explicit march's measured 3.9 effective on 4 threads (SPEED.md section 1) is set against an assumed 3 for the implicit one, and the assumption is labelled.
- **The reaction split at 20 to 50 μs** on flame zones, and finite-rate chemistry inside the step.
- **The spray's exchange,** which needs point-implicit coupling at these steps.

## Schedule (*guessed*)

1. Day 1: the residual function and its bit-identity gate. The explicit baselines on C (τ_C1, the walls, the 256x48 step).
2. Days 2 and 3: the stage solver, BDF2, and criteria 1 and 2.
3. Days 4 and 5: criteria 3 to 6, and M1 to M5.
4. Days 6 and 7: T if C passed; otherwise the report with T stated as unmeasured.

The report goes in this file under *Results*, with the logs in `implicit_prototype/`.
