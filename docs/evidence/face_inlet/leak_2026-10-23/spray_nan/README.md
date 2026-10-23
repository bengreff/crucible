# The spray NaN met after the leak fix (23 October 2026)

Found while checking L4 of FACE_INLET.md, section 5. It is a spray defect, not part of the gas budget. It is handed to the spray's owner and not fixed here.

**Repro** (measured on the Mac, 23 October, about 00:58 to 01:01 CDT). Use main 8fe9ce3 with the leak fix (this branch). Run `CRUCIBLE_T1_LEAD=0 CRUCIBLE_T1_LEVELS=2 CRUCIBLE_T1_END=2e-5 CRUCIBLE_T1_THREADS=4 crucible_rl10_tier1 run development` with Table B `rl10_s.bin`.
- It stops with `Flow::exchange left a cell non-admissible` (`lead0_fixed.log`).
- The same run with the fix reversed reaches 0.02 ms (`lead0_unfixed.log`, gas budget 2.45e-07, the leak).
- Under the run Definition's start (LOX at t_L = 20 ms) the probe to 0.02 ms has no drops and does not meet it.

**The chain**, from temporary prints that were reverted and not committed (`diag_prints.txt`):
1. **The drop.** It is LOX at 65.0845 K against T_sat 65.1366 K at the cell's 2404 Pa, so not boiling. Its diameter is 0.49 mm. The gas is at 227.0 K with Y_O2 0.524. In stage 1 the rates are finite: B_M 104.3, B_T 44.0.
2. **Stage 2's rates.** In `Spray::filmConductances` the shape F is carried from the previous film pass, and it leaves [0, 1]: 9.20 at the middle Gauss node (s 0.5) and 7.57 at the outer one (s 0.887). F is the interpolant's integral of φ over its whole.
3. **The node temperatures.** T = T_s + (T_∞ − T_s) expm1(a_T (1 − F)) / expm1(a_T) is then −475.7 K and −430.5 K. `filmTransport` there gives D = 0 and λ NaN.
4. **The film's λ.** c_pF ΣρD / Σφ is NaN, while ρD stays finite.
5. **The B_T solve.** Its residual is NaN on every secant step. Its bisection then moves `hi` down to 4.9e-324, because a NaN residual always sends it to the lower half. That B_T passes the `isfinite` check, so nothing throws.
6. **The result.** The rates' mass, heat and conductance are NaN. The parcel's m1, h1 and T1 are NaN, and so is the deposit. `Flow::exchange` then refuses the cell: finest level, cell (0, 21).

**Two defects, for the spray's owner to decide:**
- (a) The film pass's shape update is not kept monotone in [0, 1], so a strongly varying φ can push the node temperatures outside [T_s, T_∞]. Why φ varies so strongly at this drop is not traced. The near-saturation state at 2.4 kPa is the likely setting (inferred).
- (b) The B_T solve does not detect a NaN residual. It returns a denormal in place of throwing, as it does for a non-finite B_T.

The leak fix only changed the gas state at the face cells, by at most the w_j share of the wall's heat. That moved a drop onto this state. The unfixed run's drops did not reach it within 0.02 ms (inferred from the two logs).
