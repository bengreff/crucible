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

## Amendment, 24 October 2026, about 02:20 CDT (built; before any case C run is judged)

What the build changed from the plan above. Below that are the definitions criteria 3 and 4 need, fixed before the C runs.

- **No refactor of the residual was needed.** The implicit step calls the existing `Flow::rhs` (the explicit stages' operator) on its iterate, so the explicit march's code path is untouched.
  - Gate *(measured, 02:11 to 02:12)*: C1 64x12 to 1 ms on one thread. The new build's `CRUCIBLE_STATE_DUMP`, history and field files are `cmp`-identical to those of the binary built at 4f2d81e (`/tmp/c1base`).
  - The fast ctest result is below under *Results*.
- **Table A is a split. Line 22's "There is no reaction split" was wrong.**
  - `ReactingFlow::implicitStep` calls `Flow::implicitStep`, then Table A's projection over the step actually taken. The explicit march does the same after each of its steps. So both of BDF2's history levels are post-projection states.
  - The projection keeps ρ, ρu, E and the elements, so the budgets are unaffected.
  - Within a step, the composition is carried by the flow and does not relax with the temperature. So where the projection changes the state (the flame zones), the composite step can be first order in Δt *(derived)*.
  - P has no chemistry and is unaffected. Criterion 4's order bar (at least 1) stands.
- **BDF2 takes variable steps.** A rejection halves the step, so the coefficients follow the step ratio w = h/h_(n−1): U − (1+c) U^n + c U^(n−1) = a h R(U), with a = (1+w)/(1+2w) and c = w²/(1+2w). At w = 1 this is the plan's a = 2/3 and (4 U^n − U^(n−1))/3.
  - Case C's march goes to the next multiple of Δt. So a halved step is followed by one that lands back on the grid, and the runs at Δt and Δt/2 share their sample times.
- **The sweeps' Jacobian** is the Rusanov split with γ frozen at ρa²/p per cell. The bulk and the partial densities are swept together.
  - The residual keeps HLLC, the limited reconstruction and the full thermodynamics. The approximate Jacobian changes the convergence, not the converged answer.
  - Each cell's update is damped until the cell is admissible (density at most halved). No clipping.
- **The inner defect.** d = U − U_known − a h R(U). Its norm is the largest, over ρ, ρu_z, ρu_r and E, of the volume RMS of d over that component's scale. The scales are the volume RMS of ρ, ρa, ρa and E at the step's first iterate.
  - A step ends when the norm reaches the tolerance, or at the iteration limit. A step that ends at the limit is accepted and counted as *unconverged*, and every run reports that count.
- **Criterion 3's numbers, fixed now:**
  - **The observables** (criteria 3 and 4): p_inj, the outlet mass flow and F_vac, from the history file. Each is normalised by the explicit run's mean over [6, 8] ms on the same mesh.
  - **Comparisons** are made at the coarser run's sample times in [0.5, 8] ms. The other run's history is linearly interpolated there. The explicit history is sampled every 2 μs; the implicit history, after every step.
  - **The temporal error E_t:** the largest, over the three observables, of the RMS of the Δt run minus the Δt/2 run.
    - Both runs use tolerance 1e-6, with at most 400 inner iterations. If E_t/1,000 comes out below 1e-6, the pair is repeated at E_t/1,000.
  - **The sweep:** the inner tolerance (the defect norm above, a relative quantity like E_t) is set to E_t/10, E_t/100 and E_t/1,000 at Δt = τ_C1/63.
    - Pass: at the loosest tolerance where each observable's RMS difference from the E_t/1,000 run is below 0.1 E_t.
- **τ_C1 is measured as written above:** t_0.1% / ln 1000, from the explicit 64x12 run.
  - t_0.1% is the first time after which F_vac stays within 0.1% of its [6, 8] ms mean.
  - The three steps are τ_C1/25, τ_C1/63 and τ_C1/150 of the measured value.
- **M2's explicit wall per step** is measured by the explicit baseline on the same mesh, in the same batch, on one thread.
  - Both runs share the machine's load at the time. The load is stated with each number.

## Results

### Criterion 1 (P): FAIL as registered *(measured, 24 October, about 01:00 to 02:05 CDT, Mac, one thread)*

The error is the volume RMS of p − p_ref over the pulse's amplitude, to t = 2/3 L/a0.

| Acoustic CFL | 192 columns | 384 columns |
|---|---|---|
| 8 | 1.690e-1 | 1.453e-1 |
| 4 | 1.289e-1 | 8.812e-2 |
| 2 | 7.168e-2 | 3.519e-2 |
| 1 | 2.540e-2 | 1.135e-2 |
| 0.5 | 7.084e-3 | 3.239e-3 |
| 0.25 | | 8.491e-4 |

- **The judged orders, 4 to 2 and 2 to 1, are below the bar of [1.8, 2.2].** On 192 columns they are 0.85 and 1.50. On 384 they are 1.32 and 1.63.
- **The printed orders at smaller steps reach the bar.** On 192, 1 to 0.5 gives 1.84. On 384, 1 to 0.5 gives 1.81 and 0.5 to 0.25 gives 1.93.
- **The reference's own error is far below these errors.** CFL 0.1 against 0.05 differs by 1.2e-5 on 192 and 7.9e-6 on 384.
- **The errors are the scheme's, not the inner solve's.** At acoustic CFL 1 and above, the errors agree to 4 digits at inner tolerances 1e-8, 1e-9 and 1e-11.

**Reading.** BDF2 is second order here: the order rises toward 2 with each halving. But acoustic CFL 4 to 1 is before the asymptotic range for this pulse.
- Its half-width is 8 cells on 192 columns and 16 on 384. So at CFL 4 there are 2 and 4 steps per half-width, and the error is 13% and 9% of the pulse.
- The bar was set at steps too coarse for this pulse. It is not moved after the fact: criterion 1 fails as registered.

What this tells the decision *(derived from the table)*: BDF2 needs about 16 steps per half-width of an acoustic feature for an error near 1% over 2/3 of a domain crossing, and about 4 steps per half-width gives about 10%.

### Criterion 2 on P (printed; judged on C)

- The mass and energy budgets are at most 8.6e-14 in every implicit run, against the bar of 1e-11.

### M1 and M2 on P (one thread)

- **Inner iterations per step** at tolerance 1e-9 on 384 columns: 40, 23, 10, 6, 4 and 3 at acoustic CFL 8, 4, 2, 1, 0.5 and 0.25.
- **One inner iteration costs 0.54 to 0.81 explicit RK2 steps.** An inner iteration is one residual and two sweeps; an explicit step is two residuals.

### Case C development runs (not criteria)

These ran at Δt = 7.94 μs on 64x12, a step chosen before τ_C1 was measured. Logs: `implicit_prototype/cdev/`.
- **Budgets:** mass 6e-15 and energy 6e-14 over 0.3 ms.
- **The LU-SGS inner iteration contracts the defect by about 0.8 per iteration in the start-up** (`tr.log`, iterations 15 to 30). So a tolerance of 1e-6 needs about 50 to 100 iterations per step there.
  - With tolerance 1e-6 and a limit of 100 iterations, 35 of the first 43 steps (0.3 ms) ended unconverged, and 10 attempts were rejected.
  - With a pseudo CFL of 3 and a limit of 300, the first 0.1 ms took 150 iterations per step, with none unconverged.
- **One inner iteration costs 6.6e-4 to 9.0e-4 s,** with Table A's projection's share included. The explicit 64x12 step in the gate run costs 1.7e-3 s (29 s for 17,008 steps). That puts M2 at about 0.4 to 0.5 *(measured walls, both under a load of about 10 to 12 on 10 CPUs; not the judged M2)*.
