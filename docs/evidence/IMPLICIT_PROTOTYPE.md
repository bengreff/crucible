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

## Second amendment, 24 October 2026, about 02:45 CDT (after case C's first runs at the registered steps; the criteria are unchanged)

Four changes to the build and one to how criterion 2 is measured. They were made after the first case C runs, and every run is reported below, including the ones from before the changes.

- **Table A's projection can run inside the inner iteration.** This is an option: `ReactingFlow::implicitStep(dt, control, inside)`, which is the seventh field of `CRUCIBLE_IMPLICIT`.
  - Why: the split of the first amendment is first order at the injector face, and at τ_C1/63 it misses criterion 4 by a factor of about 30 (below).
  - With the option, each inner iterate is projected before its residual is formed. The final update is projected before it is judged admissible.
  - The final update's projection was added after the first inside runs rejected every full step. The update's composition, transported but unreacted, fell below the floor at the reacted gas's energy.
  - The budgets are unaffected: the projection keeps ρ, ρu, E and the elements.
- **The inner settings are pseudo CFL 3 with the extrapolated first iterate.** They were chosen from seven probes to 0.5 ms (table below). They are settings of the solver, not of a criterion.
- **BDF2's known part is formed as U^n + c (U^n − U^(n−1)).** The old form was (1 + c) U^n − c U^(n−1).
  - The two are equal in exact arithmetic. In floating point the old form rounds with a bias of about −6e-17 of U at c = 1/3. The bias does not vanish when U^n = U^(n−1), so the budgets drifted by it every step.
  - *Measured:* with the old form, the energy budget fell linearly to −8.3e-11 of the initial fill at τ/63, and to −1.6e-10 at τ/126, which has twice the steps. With the new form, the largest is 8.8e-12 and 6.7e-12, with no trend.
- **Criterion 2 is now measured as registered, against the throughput.**
  - The history's budgets are normalised by the initial fill. In C1 that is nitrogen at 1 kPa: about 0.25 J of energy, against 18.6 J supplied over 8 ms.
  - `crucible_chamber_study` now prints, for every run, the largest residual over the samples in kg and J. It gives each against the mass supplied, against the supplied total enthalpy at 300 K plus the igniter's energy, and against |E| at the end.
  - The enthalpy reference is the mechanism's (formation at 298.15 K), so propellant at 300 K carries little energy. That is why |E| at the end (416 J) is printed as a second scale.
  - The bar (1e-11) is unchanged.
- **`CRUCIBLE_CFL` sets the explicit march's CFL.** The default stays 0.4. It is used for criterion 5's diagnosis below.

## Third amendment, 24 October 2026, about 02:55 CDT (after criteria 2 to 5 were measured on 64x12; the criteria are unchanged)

- **`CRUCIBLE_IMPLICIT_WINDOW=from,until`** (seconds) runs the implicit step only over [from, until), and the explicit march outside it. It is for criterion 5's diagnosis: each scheme started from the other's settled state.
- **A failure inside the projection is a failed attempt.** Cantera's own exception (a non-finite temperature at an iterate) is now rethrown by the adapter as `std::runtime_error`, which the step takes as a rejected attempt and halves the step.
  - Why: m150 (128x24, τ_C1/150) stopped at 0.147 ms on `CanteraError: temperature must be positive. T = nan` from the projection of an inner iterate. The exception was not a `std::runtime_error`, so the step did not catch it.
  - The explicit march does not call this code; its gate is unaffected.

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

### Case C: the split is first order at the injector face *(measured, 64x12, 24 October, about 02:15 to 02:30 CDT)*

These runs used Table A's projection after the step (the first amendment's split), with tolerance 1e-6, at most 400 iterations, pseudo CFL 3 and the extrapolated first iterate. Logs: `implicit_prototype/c/t63/` and `t126/`.
- **p_inj sits above the explicit run's by an amount proportional to the step.** Over [0.5, 8] ms the RMS difference is 14.3% at τ_C1/63 (largest 14.4%) and 7.3% at τ_C1/126 (largest 7.4%).
  - Why: each step transports unreacted propellant into the first columns and then projects it. The injector face sees that composition for a whole step.
- **The explicit march is the same split at its own step** (0.06 μs, 117 times smaller than τ_C1/63). By the same scaling, its own p_inj carries about 0.12% *(derived)*.
- **The split's F_vac jumps on the run's last, short step.** At τ/63 it goes from 1668 to 2140 N on the 1.6 μs step that lands on 8 ms. The device thrust comes from the inner iterate's rates, taken before the projection, while the cell pressures are taken after it.
  - Not pursued: the split is not the build carried forward.
- **E_t for the split pair is 7.0%** (p_inj). Criterion 4 is not met at τ/63.
- 96 inner iterations per step at τ/63 over 8 ms, against 39 with the projection inside (below).

### Case C: the inner settings *(measured, probes to 0.5 ms at τ_C1/63, split, tolerance 1e-6, at most 400 iterations)*

| Probe | Pseudo CFL | Extrapolated first iterate | Iterations per step (most) | Unconverged of 72 |
|---|---|---|---|---|
| pa | 1e6 | off | 353 (400) | 62 |
| pc | 1e6 | on | 354 (400) | 61 |
| pe | 30 | off | 344 (400) | 57 |
| pd | 10 | off | 308 (400) | 47 |
| pg | 1.5 | off | 210 (272) | 0 |
| pb | 3 | off | 151 (217) | 0 |
| pf | 3 | on | 133 (230) | 0 |

- A large pseudo CFL stalls the LU-SGS iteration on the Rusanov-split Jacobian. Pseudo CFL 3 is the fastest that converges. The extrapolation saves 12%.

### Criterion 2 on C: PASS *(measured, with the projection inside, the new known part)*

- **In every run with the second amendment's build,** the largest residual over the samples is:
  - mass: at most 4.3e-16 of the mass supplied over 8 ms (c25 to c150, s10, s100, k63, k126, m25, m63, w1), and 4.2e-15 over v63's 1 ms;
  - energy: at most 1.5e-13 of the supplied total enthalpy and the igniter's energy over 8 ms, and 1.2e-12 over v63's 1 ms; at most 9e-15 of |E| at the end.
  - The bar is 1e-11. w2, mostly explicit, gives 5.8e-16 and 4.2e-13.
- **The explicit 64x12 run, for comparison:** its history's largest energy budget is 3.2e-11 of the initial fill. That is about 8e-12 J, or 4e-13 of the throughput *(derived from the fill's 0.245 J)*.
- **The runs with the old known part** (i63 and i126, logs kept) drifted to 8.5e-11 and 1.65e-10 of the initial fill. That is 1.1e-12 and 2.2e-12 of the throughput *(derived)*. So they would pass the registered bar too, but the drift was a defect, and it is fixed.

### Criterion 3: PASS at tolerance E_t/10 = 8.975e-5 *(measured, 24 October, about 02:38 to 02:41 CDT)*

- **The pair at τ_C1/63 and τ_C1/126, tolerance 1e-6:** E_t = 8.975e-4, set by F_vac. The other two observables' RMS differences are 1.43e-4 (p_inj) and 1.45e-4 (outlet flow).
  - E_t/1,000 came out below 1e-6, so the pair was repeated at 8.9748e-7, as registered. That gave E_t = 8.9749e-4.
- **The sweep at τ_C1/63.** The E_t/1,000 run is the repeated pair's τ/63 member, k63. The bar is 0.1 E_t = 0.0090%.

| Tolerance | p_inj RMS | outlet flow RMS | F_vac RMS | Iterations per step (most) |
|---|---|---|---|---|
| E_t/10 = 8.975e-5 | 0.00033% | 0.00062% | 0.00194% | 31.0 (102) |
| E_t/100 = 8.975e-6 | 0.00005% | 0.00013% | 0.00032% | 33.9 (214) |
| E_t/1,000 = 8.975e-7 | (reference) | | | 39.7 (243) |

- **The loosest tolerance passes, by a factor of 4.6** on F_vac. 8.975e-5 is the tolerance for criteria 4 and 5 and for M1 to M3.
- **The iterations barely depend on the tolerance:** 31 against 40 for 100 times tighter. Most steps reach the tolerance in a few iterations from the extrapolated iterate. The mean is set by the start and the fill (M1 below).

### Criterion 4 on 64x12: the bars are met at τ_C1/63; the order FAILS as registered *(measured, 24 October, about 02:40 to 02:50 CDT)*

Against the explicit run at CFL 0.4 (e64), over [0.5, 8] ms, tolerance 8.975e-5, each observable normalised by e64's [6, 8] ms mean:

| Step | p_inj RMS (largest) | Outlet flow RMS (largest) | F_vac RMS (largest) | Iterations per step (most) |
|---|---|---|---|---|
| τ_C1/25 (c25) | 0.459% (3.63%) | 0.128% (0.684%) | 0.257% (1.92%) | 69.1 (400) |
| τ_C1/63 (c63) | 0.0993% (0.330%) | 0.0166% (0.236%) | 0.169% (1.59%) | 31.0 (102) |
| τ_C1/150 (c150) | 0.0970% (0.159%) | 0.0064% (0.081%) | 0.129% (0.571%) | 12.4 (400) |

- **At τ_C1/63 the bars are met:** the largest RMS is 0.17% (bar 0.5%), and the largest single difference is 1.59% (bar 2%). At τ_C1/25 the RMS passes, but p_inj's largest (3.63%) does not.
- **The order bar (at least 1) fails for p_inj and F_vac.** The observed orders are:
  - p_inj: 1.66 from τ/25 to τ/63, then 0.03 from τ/63 to τ/150;
  - outlet flow: 2.2, then 1.10 (passes);
  - F_vac: 0.46, then 0.31.
- **Why: the differences stop at a floor that does not depend on the step.** The floor is the settled offset of criterion 5 below. The implicit's settled state is the same to 7 digits at all three steps (p_inj 3.274683 MPa, F_vac 1679.906 N), so at τ/150 the RMS is almost all offset.
- The first 0.5 ms (not judged): the outlet flow's largest difference is 11.9%, 10.6% and 7.7% at the three steps; the starting shock is smeared, as expected.

### Criterion 5: FAIL *(measured, 64x12; 128x24 below)*

Over [6, 8] ms against e64: F_vac and vacuum Isp are −0.0768% at τ_C1/25, /63 and /150, against the bar of 0.05%. c* is +0.012% (2447.80 against 2447.50 m/s; from the logs' *simulated* line, because the history's c* column is empty in these runs).

**The offset is not the time step, not the inner tolerance, and not the explicit reference's step error.** *(measured unless labelled)*
- **Not the step:** the implicit's settled F_vac is 1679.906 N at all three steps.
- **Not the tolerance:** criterion 3's sweep moves F_vac by at most 0.0019% between tolerances 8.975e-5 and 8.975e-7.
- **The explicit's own step error, from a run at CFL 0.2** (x64h, `CRUCIBLE_CFL=0.2`, settled means over [6, 8] ms):

| | e64 (CFL 0.4) | x64h (CFL 0.2) | First-order limit *(derived)* | Implicit (all three steps) |
|---|---|---|---|---|
| p_inj (MPa) | 3.277835 | 3.276252 | 3.274669 | 3.274683 |
| Wall axial force (N) | −4471.224 | −4471.564 | −4471.904 | −4471.911 |
| F_vac (N) | 1681.198 | 1681.236 | 1681.273 | 1679.906 |

- **The explicit's p_inj and wall force move with its step, and the implicit sits on their limit,** within 0.0004% and 0.007 N. That is the split at the injector face (first amendment), at the explicit's own 0.06 μs step.
- **F_vac does not:** it barely moves with the explicit's step (+0.04 N), and the implicit is 1.37 N (0.081%) below its limit.
  - With the wall force agreeing, the 1.37 N is in the inlet's momentum flux (F_vac − wall force, about 6,153 N). That is a difference of 0.022% in the inlet's gross flux, which the net thrust (3.7 times smaller) shows as 0.08% *(derived)*.
  - It is not the explicit's booking of its two stages: that term is proportional to the explicit step, and halving the step moved F_vac by 0.04 N, not by half of 1.37 N.
- **Each scheme takes the other's settled state to its own within about 10 μs** (`CRUCIBLE_IMPLICIT_WINDOW`, below):
  - w1 (implicit to 7.9 ms, then explicit): p_inj rises from 3.27468 to 3.27821 MPa and F_vac from 1679.91 to 1681.3 N within 10 μs. p at the chamber end does not change.
  - w2 (explicit to 7.9 ms, then implicit): p_inj falls from 3.27784 to 3.27421 MPa within the first two steps (7.4 μs), and F_vac is 1680.03 N at 8 ms, approaching 1679.9.
  - So the two discretisations hold different stationary states at the injector face. The difference is set within a few cells of it, and is not a failure of the implicit to settle.
- **On 128x24 the offset is larger** (below). An offset that grows as the mesh is refined is not a truncation error of either scheme.
- **Cause: not identified.** The candidates are how the inlet face's flux sees the projected composition in the two schemes, and the face's boundary state. Neither was tested tonight.

`CRUCIBLE_IMPLICIT_WINDOW=from,until` (in seconds) runs the implicit step only over [from, until) and the explicit march outside it. It was added for this diagnosis (third amendment, above).

### Criterion 6: PASS *(measured, v63: 256x48 to 1 ms at τ_C1/63, tolerance 8.975e-5, 24 October, about 02:43 to 02:48 CDT)*

- **The run completes:** 148 steps, 7 attempts rejected (418 inner iterations spent on them), each recovered at the halved step. One step ended at the limit of 400 iterations (unconverged).
- **151.3 inner iterations per step** (most 400). At the same step, 64x12 takes 52.4 per step over the first 0.5 ms (c63), and 128x24 takes 68.4 per step over 8 ms (m63).
- **The step ratio** is about 480 *(derived: the explicit step halves with each doubling of the mesh, 5.916e-8 s on 64x12 and 2.936e-8 s on 128x24, measured)*. The explicit 256x48 baseline (e256) was queued behind the other jobs and gives the measured step when it runs.
- Budgets: 4.2e-15 of the mass supplied, 1.2e-12 of the energy supplied (criterion 2).

### M1 to M3 on 64x12 *(measured, one thread, tolerance 8.975e-5; machine load about 13 to 15 on 10 CPUs throughout)*

The explicit baseline e64b (same build, same batch): 136,227 steps to 8 ms in 296.0 s, 2.173e-3 s per step.

| Step | Step ratio | M1 start [0, 0.5] ms | M1 fill [0.5, 3] ms | M1 settled [3, 8] ms | M2 (explicit steps per inner iteration) | Wall to 8 ms | M3 |
|---|---|---|---|---|---|---|---|
| τ_C1/25 | 300 | 165.7 (most 400) | 63.1 | 62.6 | 0.659 | 44.9 s | 6.59 |
| τ_C1/63 | 119 | 52.4 (most 102) | 29.7 | 29.5 | 0.662 | 51.1 s | 5.79 |
| τ_C1/150 | 50 | 22.4 (most 400) | 11.8 | 11.7 | 0.722 | 52.7 s | 5.61 |

- **The speed-up is about 6 and barely depends on the step.** The inner iterations to 8 ms are 31,366, 35,508 and 33,628 at the three steps.
  - The settled iterations per step are 0.21, 0.25 and 0.23 of the step ratio: LU-SGS on this Jacobian needs about one inner iteration per 4 to 5 explicit steps of physical time, whatever the step *(derived from the table)*.
  - So a larger step does not buy more speed-up here; the inner solver's convergence rate sets it. A better inner solver (a stronger preconditioner, or Newton-Krylov) is what would raise it.
- **M3 is at criterion 4's accuracy only at τ_C1/63 and τ_C1/150,** where criterion 4's bars are met (its order bar is not). At τ_C1/25 p_inj's largest difference (3.63%) misses the bar.
- **Criterion 5's offset applies to all of them:** the settled thrust is 0.077% below the explicit's.

### Case C on 128x24 *(measured, one thread, tolerance 8.975e-5 unless stated; 24 October, about 02:42 to 03:10 CDT)*

Against the explicit run at CFL 0.4 on the same mesh (e128, explicit step 2.936e-8 s, so the step ratios are 601, 239 and 100):

| Run | Step | Iterations per step (most) | Unconverged, rejected | p_inj RMS (largest) | Outlet flow RMS (largest) | F_vac RMS (largest) | F_vac over [6, 8] ms |
|---|---|---|---|---|---|---|---|
| m25 | τ_C1/25 | 137.7 (400) | 5, 0 | 5.60% (7.52%) | 0.187% (0.812%) | 1.09% (2.75%) | 1669.36 N (−1.03%) |
| m63 | τ_C1/63 | 68.4 (199) | 0, 2 | 0.206% (2.07%) | 0.084% (0.591%) | 0.400% (2.72%) | 1680.52 N (−0.370%) |
| m150b | τ_C1/150 | 27.0 (76) | 0, 1 | 0.104% (0.226%) | 0.011% (0.128%) | 0.353% (0.883%) | 1680.52 N (−0.370%) |

- **Criterion 4 on 128x24: not met at τ_C1/63, met at τ_C1/150.** At τ/63 the RMS values pass, but the largest differences do not: p_inj 2.07% and F_vac 2.72%, against 2%. At τ/150 every bar is met. The order bar fails as on 64x12 (F_vac's RMS: 1.09%, 0.400%, 0.353%).
- **Criterion 5 on 128x24: FAIL, and the offset is larger.** e128's settled F_vac is 1686.76 N. m63 and m150b settle to the same state to 6 digits (p_inj 3.264434 MPa, F_vac 1680.525 N), 0.370% below it.
  - p_inj is 0.10% below e128's (3.267774 MPa), as on 64x12: the explicit's own split error at the injector face.
  - The wall force differs by 0.34 N, so the inlet's momentum flux carries about 5.9 N of the 6.2 N. On 64x12 it was 1.37 N. The gap grows about four times when the mesh is refined twice *(derived)*.
- **m150 stopped at 0.147 ms** on the Cantera exception in the projection (third amendment). m150b is the same run with the fix: the one attempt that failed was rejected and recovered at the halved step.
- **m25 settles somewhere else, and that state is not a steady state of the equations.**
  - From 4 to 8 ms it holds p_inj 3.093075 MPa and F_vac 1669.612 N, constant to 7 digits. m63's are 3.264434 MPa and 1680.525 N. The histories part already in the fill: at 0.55 ms p_inj is 1.467 MPa in m25 and 1.652 MPa in m63.
  - Every settled step starts from a defect of 5.4e-2 and takes 123 iterations. At a steady state the extrapolated first iterate should already be near the answer (m63's settled steps start at 1.8e-3).
  - **The run's last, short step (4.67 μs, step ratio w = 0.26) moves p_inj by +3.7% and F_vac by −1.8%** (to 3.2075 MPa and 1640.3 N). A steady state of the semi-discrete equations does not move on any step. So m25's state is a fixed point of the BDF2 map at this step, not of the flow.
  - **Not the tolerance:** m25t, the same run at tolerance 1e-6 to 5 ms (330 iterations per step, 31 unconverged at 400), holds the same state (p_inj 3.0930 MPa, F_vac 1669.6 N at 4.4 ms). Its last, short step jumps the same way (to 3.2096 MPa and 1643.2 N).
  - **Cause: not identified.** One candidate, untested: the state is projected again after the step (`ReactingFlow::implicitStep`'s `react`), and if Table A's projection is not exactly idempotent, each step's result is kicked off the inner iteration's fixed point by the same amount. That would give a constant first defect. A step ratio of 600 amplifies such a kick most.
  - No 64x12 run shows this: c25's last steps start at 2.2e-3, and its last short step moves nothing.

### M1 to M3 on 128x24 *(measured, one thread; load about 8 to 31 on 10 CPUs, varying: the fast ctest took 4 threads from 02:57)*

The explicit baseline e128: 275,738 steps to 8 ms in 1,755 s, 6.365e-3 s per step (02:39 to 03:08).

| Step | Step ratio | M1 start | M1 fill | M1 settled | M2 | Wall to 8 ms | M3 |
|---|---|---|---|---|---|---|---|
| τ_C1/25 (m25) | 605 | 286.3 (most 400) | 138.6 | 122.7 | 0.633 | 251.9 s | 6.97 |
| τ_C1/63 (m63) | 241 | 102.8 (most 199) | 65.2 | 66.5 | 0.629 | 313.3 s | 5.60 |
| τ_C1/150 (m150b) | 101 | 38.8 (most 76) | 26.2 | 26.2 | 0.786 | 368.0 s | 4.77 |

- **At criterion 4's accuracy on 128x24 (τ_C1/150 only), M3 is 4.8.** m25's 7.0 is at the wrong settled state above.
- The settled iterations are 0.20, 0.28 and 0.26 of the step ratio, as on 64x12 (0.21 to 0.25).
- All three ran while e128 ran, as the amendment asks, but the load changed over e128's 29 minutes. m150b (M2 0.79) ran beside m25t and the ctest. So M2 on 128x24 carries the load's variation *(the walls are measured; how evenly the load fell is not)*.

### What the measurements say for the decision *(derived from the tables above; M4 and M5 not done)*

- **The measured speed-up is 5 to 6 on one thread, at the steps where the history is accurate** (64x12 at τ_C1/63: 5.8; 128x24 at τ_C1/150: 4.8).
- **It is set by LU-SGS's convergence, not by the step.** The inner iterations per unit of physical time are about the same from τ_C1/25 to τ_C1/150: about one inner iteration per 4 to 5 explicit steps, at 0.63 to 0.79 explicit steps each. The bound this gives is about 1/(0.23 × 0.66), or 6.6, whatever the step.
- **Two failures stand before trusting the numbers** (the plan's rule: 1, 2, 3 and 5 must pass): criterion 1 (pre-asymptotic order on P) and criterion 5 (a settled thrust 0.08% low on 64x12 and 0.37% on 128x24, growing with the mesh, cause not identified). The m25 fixed point is a third problem.
- **M4 (the projection to the lead) and M5 (finite-rate chemistry at the implicit step) are not done.** The proposed rule (build at 5 or more, not below 2) needs M4. On these measurements, a projection that keeps LU-SGS would start from about 5 to 6 before threads, AMR and turbulence, which the prototype does not measure.

### The gate and the fast ctest

- **The explicit march is unchanged by every build** (C1 64x12 to 1 ms, one thread: state dump, history and 10 field files `cmp`-identical to 4f2d81e's). First build 02:11 (`gate/result.txt`), second 02:38 (`result2.txt`), third 03:10 (`result3.txt`, all targets rebuilt).
- **The fast ctest on the first build's binaries (02:13, e8816dc's code):** only the 4 known failures (amr, amr_subcycling, transport and spray verification); 11 of 15 pass. Ran 02:57 to 03:07 (`ctest_2026-10-24_e8816dc.txt`).
- **The fast ctest on the third build:** queued at 03:12 behind the other slot jobs. Output `/tmp/implicit/ctest_2026-10-24_final.txt`, done-file `/tmp/implicit/ctest3_done`.
