# C2: PaSR turbulence-chemistry closure

Status (5 October 2026, 03:12): implemented; criteria 1, 2 and 3 pass (run 5). Criterion 2's test-power check failed twice for the fractional-segregation case because of where its samples fell. It was restated at 03:08 (53f8870), before run 5, to read a fine record of the references, and run once: it passes. Criteria were stated in the header of `tests/pasr_tests.cpp` on 4 October 2026, before the closure was written (c259fb5). The amendments are dated there. Criterion 4 (a turbulent reacting chamber with the closure, judged against the same run without it) was stated at 03:09 on 5 October (c44a933); it has not been run.

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

## Limits

- The reference shares Cantera's rates and the closure formula with the engine. It checks the integration and the wiring, not the closure's physical validity.
- The closure's constants (C_mix 1, Sc_t 0.7, beta* 0.09) and its form are model choices. Whether PaSR is right for the RL10's turbulence-chemistry interaction is a modelling question for validation, not verification.
- The fast suite after the change: 5 of 6 pass. The output of `transport_verification` is identical to the recorded open-face run: the same three truncation-order failures, with the same values.
