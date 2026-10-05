# C2: PaSR turbulence-chemistry closure

Status (5 October 2026, 03:12): implemented; criteria 1, 2 and 3 pass (run 5). Criterion 2's test-power check failed twice for the fractional-segregation case because of where its samples fell. It was restated at 03:08 (53f8870), before run 5, to read a fine record of the references, and run once: it passes. Criteria were stated in the header of `tests/pasr_tests.cpp` on 4 October 2026, before the closure was written (c259fb5). The amendments are dated there. Criterion 4 (a turbulent reacting chamber with the closure, judged against the same run without it) was stated at 03:09 on 5 October (c44a933). Its harness is written (`tests/chamber_study.cpp`, modes `frp` and `frt`), but it cannot run yet: a smoke run on 32x6 stopped at 0.136 ms, before the igniter fires, in both the closure run and the control. Three faults in the k-omega numerics stop the viscous cold start one after another, none in the closure. The first (axis-row reconstruction) is fixed. The other two (off-line diffusion, omega underflow in the nozzle) are diagnosed and not yet fixed (see Criterion 4 below).

All numbers are measured unless they are marked derived or inferred.

Raw outputs, all in `pasr/`:
- `pasr_run1` to `pasr_run5` (`*_2026-10-05.txt`): the five runs of `crucible_pasr_tests`.
- `rk4_gaps_plain` and `rk4_gaps_compensated`, `heat_release` and `duct_3d_diagnosis` (`*_2026-10-05.txt`): the diagnostics behind the amendments, made by `diag.cpp` and `rise.cpp`. These programs include the test file; they are diagnostics, not tests.
- Fast ctest suite after the change: `ctest_main_2026-10-05_pasr.txt`.

## What was added

- **The closure** (TECHNICAL_PLAN step 7). In each cell and reaction substep, CVODES integrates dz/dt = kappa_eff(z) f(z).
  - kappa_eff = 1 − s(1 − kappa), with kappa = tau_c / (tau_c + tau_mix).
  - tau_c is recomputed at every right-hand-side call: the largest Y_i / |dY_i/dt| over H2, O2 and H2O, using only the species that are present and changing.
  - tau_mix = C_mix sqrt(nu_eff / epsilon), with epsilon = beta* k omega.
  - The segregation s is the largest over those species of min(1, (nu_t / Sc_t) |grad X_i|^2 / (beta* omega X_i (1 − X_i))). It is the ratio of the scalar-variance production to its dissipation, at equilibrium.
  - tau_mix and s are frozen over the substep.
- **Where it lives.**
  - `Flow::mixingInputs` computes tau_mix and s from the flow's own transport, gradients and eddy viscosity.
  - `ReactionStep::advance` takes them through a `Mixing` argument.
  - `ReactingFlow::setMixingClosure` switches the closure on. It is accepted only for FiniteRate chemistry on a turbulent flow.
  - With s = 0 the code takes the laminar path bit for bit.

## Results (run 5)

| Criterion | Measured | Limit | Result |
|---|---|---|---|
| 1(a) kappa_eff against the formula, 3 states x 3 closures | 0 | 1e-13 | pass |
| 1(a) pure N2 gives exactly 1 | exact | exact | pass |
| 1(b) tau_mix and s from the duct, k = 10 and 25 (s 0.25 to 1.0) | 2.6e-16, 2.7e-15 | 1e-12 | pass |
| 1(b) k = 0 gives s = 0 in every cell | exact | exact | pass |
| 1(c) accepted only for FiniteRate with turbulence | 4 of 4 | all | pass |
| 2, (1e-5 s, 1): T and Y against RK4 at rtol 1e-10 | 1.4e-8, 2.6e-9 | 1e-7, 1e-8 | pass |
| 2, (1e-5 s, 1): the closure changes T (test power), fine record | 0.50 | > 1e-2 | pass |
| 2, (1e-5 s, 0.9): T and Y against RK4 | 5.4e-10, 3.1e-10 | 1e-7, 1e-8 | pass |
| 2, (1e-5 s, 0.9): the closure changes T (test power), fine record | 0.46 | > 1e-2 | pass (restated check; 3.5e-3 at the 20 samples in run 4) |
| 2, laminar: T and Y against RK4 | 6.9e-10, 4.0e-10 | 1e-7, 1e-8 | pass |
| 2, every case: the error at rtol 1e-10 is below the error at rtol 1e-7 | yes | | pass |
| 3(a) s = 0 equals the laminar substep, 3 states | bitwise | bitwise | pass |
| 3(b) order in tau_mix as it goes to 0 | 1.00 | at least 0.8 | pass |
| 3(c) pure N2 unchanged | exact | exact | pass |
| 3(d)(i) inert duct: with the closure against without it | bitwise | bitwise | pass |
| 3(d)(ii) inert duct against a plain Flow step | 1.6e-16 | 1e-14 | pass |

## Runs and amendments

- **Run 1.** The test duct switched turbulence on without positive ambient k and omega, and the engine rejected it. Fixed in the test: the ambient values are now the uniform field, k 10 and omega 1e3.
- **Run 2.** Criterion 2's reference did not converge.
  - Classical RK4 on this mechanism is unstable at every step up to 2^17 per sample (h 1.1e-10 s). It is stable from 2^18 on (`rk4_gaps_plain`).
  - So the stiffest rate is about 2.5e10 /s (derived from RK4's stability bound of 2.78).
  - Between 2^18 and 2^19 per sample, T and Y change by 9e-13 and 7e-13, inside the stated 1e-10. With compensated summation the change is 1e-16.
  - The test's halving stopped at 2^17, a cap the criterion does not state. It now goes to 2^20.
- **Run 2, criterion 3(d).** It failed at 0.77, but not because of the closure.
  - ReactingFlow without the closure gives the same gaps (`duct_3d_diagnosis`).
  - The gaps are 1e-15 in axial momentum and 8e-15 in radial momentum. At rest those fields have maxima of only 1.5e-4 and 1e-14, so "relative to the field's max" turns rounding into large numbers.
  - The reaction substep rewrites Y at the 1e-16 level, even for inert gas (inferred).
  - Amended: (i) with the closure, the step must equal the step without it bit for bit; (ii) against plain Flow, momentum is scaled by sqrt(rho_max p0).
- **Run 3.** The case (1e-4 s, 0.5) passed its accuracy checks but changed T by only 1.5e-3 against laminar. It was replaced by (1e-5 s, 0.9), with the T change predicted above 1e-2 before the run.
- **Run 4.** The prediction failed: 3.5e-3.
  - The heat-release diagnosis (`heat_release`) gives the mechanism. Every case ignites near 18 us, and the closure stretches the 10-90% temperature rise:

    | Case | Rise (10-90%) | Smallest kappa_eff |
    |---|---|---|
    | laminar | 3.6 us | 1 |
    | s 0.5 | 4.5 us | 0.52 |
    | s 0.9 | 7.8 us | 0.14 |
    | s 1 | 12.9 us | 0.046 |

  - The samples are 15 us apart. The test-power check passes only if the stretched rise is still running at the 30 us sample: with s = 1 it ends at 32.4 us, but with s = 0.9 at 27 us.
  - So the check measures where the samples fall, not the closure.
  - After two attempts at the fractional case, work on it stopped (the working rule). It stays failing until the check is restated.
  - The closure's fractional-s arithmetic is covered by 1(a), at s 0.5 and 0.2.
- **Restatement (03:08, 5 October, at the Director's instruction; committed in 53f8870 before run 5).** The check now reads T from both converged references at 1280 equal times (every 1/64 of a sample, 0.23 us apart), so the rise is resolved wherever it falls. The cases, the 20 samples and the accuracy checks are unchanged. Predicted for (1e-5 s, 0.9): order 0.1.
- **Run 5, run once.** Every check passes.
  - The fractional case: the largest gap is 0.457 at 21.1 us, where T is 1648 K with the closure and 3033 K laminar. For s = 1 it is 0.498 at 22.0 us.
  - The prediction was low by about 4.6 times. It assumed both rises start together. In the record, the closure run is 21% through its rise (1200 K to 3369 K) when the laminar run is 85% through, so the closure also delays the start of the rise.
  - Every accuracy number is the same as in run 4. The fine record only adds stores of T to the reference integrations.

## Criterion 4: harness and first smoke run (5 October, 03:45)

- **Harness.** `crucible_chamber_study frp|frt <nz> <nr> <end> <threads> <prefix>` builds the C1 chamber, made viscous and turbulent as in `tests/step_cost.cpp`, with the closure on (`frp`) or off (`frt`). It judges (a) to (d) and prints the reported items.
  - The igniter is the one C1's FiniteRate run used: 300 J over 1 ms from 0.2 ms. The criterion's statement says only "igniter", so this is written down here and in the harness header before any run.
  - Light-off is the first 2 us sample at which some cell has Y_H2O > 0.5. The chamber volume for the s > 0.01 fraction is the cells upstream of the throat. Both definitions were declared with the harness, before any run.
  - `ReactingFlow::Stats::maxClippedFraction` is new. It reports the largest negative mass fraction set to zero after a reaction substep.
- **Smoke run, 32x6 to 0.5 ms (measured, [criterion4_smoke_2026-10-05.txt](pasr/criterion4_smoke_2026-10-05.txt)).** Both runs stop at 0.136273 ms with "Could not advance an admissible gas state after 14 timestep reductions". This is before the igniter starts at 0.2 ms, so the closure plays no part.
- **Where it fails (measured, with a temporary print, since reverted).**
  - The failing cell is the axis cell one column past the throat (i 18, j 0).
  - The failing check is rho*omega > 0 after the first explicit stage.
  - omega there has fallen to 1.9e-3 1/s. Its neighbours hold 2 to 5e3 1/s along the axis and 4.5e3 1/s off it. The ambient value is about 78 1/s (derived from the Spalart-Rumsey formula).
  - The stage drains rho*omega at about 6e5 kg/(m^3 s^2) (derived from two halvings), whatever the cell holds. So halving the step only helps until rho*omega is below that rate times the step.
- **Isolation, one fix, and two more faults (03:45 to 03:57; details in the evidence file).** Three faults stop the cold start one after another. Each lets omega's update remove more than the cell holds, or lets omega reach zero:
  1. **Axis-row reconstruction (fixed, kept).** The axis row reconstructs mass fractions, k and omega with an r^2 curvature that nothing bounds by the cell value. A cell near zero beside large off-axis values gets a large face value, and outward flow past the throat carries it out. The wall row already clamps its slope so faces stay within half the cell value. The axis row now gets the same clamp (`core/flow.cpp`, `radialProfiles`). With it, the run passes 0.136 ms and lights at 0.206 ms.
  2. **Off-line part of the k and omega diffusion (inferred from a test; not fixed).** With the clamp the run stops at 0.346 ms in the divergent nozzle (i 26, j 1), where omega is about 1e-32 1/s and is drained at about 3.6e3 per second. With the k and omega face gradient replaced by the two-point difference (a temporary switch), it reaches 0.5 ms. This is the term `core/flow.cpp` already documents for k.
  3. **omega decays to underflow in the nozzle (inferred; not isolated, not fixed).** With both changes the run stops at 0.524 ms: omega is exactly 0 and k about 1e-109 in two nozzle cells (i 26 j 4, i 27 j 3). The SST sources take omega down with no floor. The candidates are the dilatation part of the production and the cross-diffusion sink. The engine fills with Spalart-Rumsey ambient values but does not have the SST-2003sust sustaining terms that hold them.
  - The first test, two-point diffusion without the clamp, still stopped at 0.136 ms with the same drain rate. That is why the reconstruction, not the diffusion, is the first fault.
  - Stopped after three distinct faults, before writing fixes for 2 and 3.
- **Verification after the clamp (Mac).** The fast suite: 5 of 6 pass. `transport_verification` prints output identical, line for line, to the recorded open-face run, with the same three truncation-order failures. `turbulence_verification` passes. The slow pipe verification (TURBULENCE_C2 check 4) was then run with the clamp ([pipe_after_axis_clamp_2026-10-05.txt](pasr/pipe_after_axis_clamp_2026-10-05.txt)). nr 16 and 32 print every history line and every judged value as the recorded Mac run before the clamp, including the failing mass budgets 3.42e-11 and 1.93e-11. The nr 64 run was stopped at 0.5 ms to free the Mac; its 0.5 ms line matches the backhouse run without the clamp. So the clamp does not touch the pipe (measured).
- **Next.**
  1. Limit the off-line part of the k and omega face gradients, so a face's flux keeps the sign of the two-point difference (as OpenFOAM's limited Laplacian does). That makes their diffusion monotone under the existing step limit.
  2. Add the SST-2003sust sustaining terms (published on NASA TMR; a model variant, recorded as a choice) or a floor at the ambient values. The sustaining terms are preferred.
  3. Rerun the turbulence verification and this smoke case, then criterion 4.
  4. C1's recorded runs predate the clamp. In inviscid C1 the clamp acts on the axis-row mass fractions only where a cell holds less than about a third of the next ring's value (derived for uniform rings). A C1 rerun shows whether its numbers move.
- **Consequence.** No viscous turbulent chamber can start from cold until faults 2 and 3 are fixed. That covers criterion 4, C3 and the RL10 case. The step-cost runs did not see this: they start from a settled equilibrium field, with the valve fully open.

## Limits

- The reference shares Cantera's rates and the closure formula with the engine. It checks the integration and the wiring, not the closure's physical validity.
- The closure's constants (C_mix 1, Sc_t 0.7, beta* 0.09) and its form are model choices. Whether PaSR is right for the RL10's turbulence-chemistry interaction is a modelling question for validation, not verification.
- The fast suite after the change: 5 of 6 pass. The output of `transport_verification` is identical to the recorded open-face run: the same three truncation-order failures, with the same values.
