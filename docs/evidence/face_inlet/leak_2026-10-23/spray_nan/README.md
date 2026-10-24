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

## Resolved (23 October 2026, about 23:20 to 23:25 CDT, Mac; measured)

**Why F leaves [0, 1].** Temporary prints, reverted and not committed, traced it.
1. The first film pass starts from the constant-phi shape (F = s, inside [0, 1]). It puts the node nearest the surface (s 0.887) at 66.87 K, with Y_O2 0.992.
2. There `filmTransport` returned a **negative conductivity**, λ = −0.00589 W/(m K), so phi = c_pF ρD / λ came out at −1.58.
3. The sum of phi over the nodes nearly cancelled (0.216, then 0.051). The interpolant's integral divided by it gave F = 2.1, then 9.2. Everything after that follows: the node temperatures go negative, then λ goes NaN.

**The cause: the fits' range.**
- The tier-1 driver built its transport fits with `thermo::transportFits(mechanism)`. Those run from the mechanism's lower thermo limit, 200 K or more.
- Below about 97 K, H2's conductivity polynomial turns negative (`adapters/reaction.hpp`).
- The film between a drop at T_s 61 to 67 K and the gas spans exactly that range.
- Mascotte and the spray tests already fit from 60 K (`transportFits("h2o2.yaml", 60)`, SPRAY.md, Results: S2). Tier 1 did not.

**On a positive phi.** With phi > 0 at every node, F at the nodes can leave [0, 1] only by the 3-point interpolant's overshoot: about 0.08 at most, with phi all at one node (its Lagrange weights). The run showed F = −0.0017 at s 0.113, which moves that node's T about 1% above T_∞. That is bounded and is not the NaN.

**The fixes.**
- (a) Tier 1 fits from 60 K (`tests/rl10_tier1.cpp`), registered as an amendment in `tier1_run_definition.md`. G4's printed H2 viscosity at 199.67 K moves from 6.844e-06 to 6.833e-06 Pa s.
- (b) `Spray::filmTransport` throws when the conductivity is not positive, or the diffusivity is negative or not finite. The message names the temperature and the fits' range.
- (c) The B_T solve throws on a NaN residual, in place of bisecting to 4.9e-324.

**The checks.**
- The repro with (a) to (c) reaches 0.02 ms (`lead0_fits60.log`): 108 flashes, budgets (a) 1.83e-15, (b) 6.64e-16, (c) 1.82e-15, 0 non-finite values.
- With (b) and (c) but the old fits, it stops at the first base step: `Spray: the film's transport is not positive at 66.870342 K: the transport fits do not reach the film's temperature` (`lead0_loud.log`).
- The same pattern holds for the flash branch's cold drops.
