# The cost of the viscous step: profile, decision and criteria

Status: profile measured and criteria stated 13 October 2026, before code (TECHNICAL_PLAN *Lightweight engine*, candidate 1 of *What the budget says*; the Director's brief of 8 October, 11:10, put "tabulated transport properties (profile first)" ahead of the RL10 development run).

## The profile

The turbulent C1 (`crucible_chamber_study eqtt 64 12`, Table A, SST with the property corrections, wall functions `law`) on one thread, sampled for 90 s during its first 1.2 ms (`transport_cost/profile_eqtt64x12_2026-10-08.txt`; measured). It ran 1 ms in 17,048 steps and 78 s, 4.6 ms per step, 6.0 us per cell update.

| Part of the step | Share (inclusive, measured) |
|---|---|
| Mixture-averaged transport properties, `Medium::transport` | 19.8% |
| Thermodynamics: the e-to-T Newton inversion, cv, enthalpies (`Medium::properties`, `temperature`, `cv`, `enthalpy`, `speciesEnthalpies`) | about 16% |
| Gradients (`transportGradients`) | 15.1% |
| Face fluxes of the transport operator (`transportFluxes`' own work) | 13.3% |
| Inviscid faces (axial, radial, supply) | 19.3% |
| Wall functions (`wallSolve`, `prescribeWallTurbulence`) | 10.1% |
| Eddy viscosity and SST sources | 10.3% |
| Chemistry (Table A lookups) | 4.2% |

- The chain refresh, transport properties, gradients and eddy viscosity runs four times per step: in `stableDt`, in each of the two right-hand sides, and after the step for the second half of the SST source split. `stableDt` alone is 22% of the step. The chemistry step changes the state between the last of these and the next step's first, so none is a duplicate.
- Inside `Medium::transport` the time is in the divisions of the inner loops: the diffusion sums `x_j / D_jk` (28% of the function), the Wilke weights' `mu_k^(1/2) / mu_j^(1/2)` and their division by `sqrt(8 (1 + W_k/W_j))` (30%), and the 55 binary-diffusion polynomials (22%). Inside the thermodynamics the time is in the NASA-7 enthalpy polynomial, which divides by 3, 5 and T for every species at every call.
- Threading at this size scales poorly (13 October, Mac load average 5.5, to 0.3 ms): 22.7 s on 1 thread, 15.2 s on 2, 13.1 s on 4 and 19.2 s on 8 (measured). 768 cells give each of the step's parallel loops about 0.1 ms of work, comparable with the cost of handing it out. The RL10's 8,500 cells (budget guess) give ten times more per loop.

## Decision

A table of transport properties is not the lever the plan guessed. The properties are 20% of the step, so removing them entirely would buy 1.25 times. The tables considered:
- **On Table A's manifold** (mu, lambda and D_k as functions of f_H, T and rho_r). The mixing rules are not linear in the diluents (N2 fills the C1 chamber at the start), so the diluent fractions would need axes, and states off the manifold (Table B's mixed states, the finite-rate modes, the frozen nozzle) would need a fallback. Not adopted.
- **In T alone, per species and pair** (the Wilke factors and the inverse binary diffusion coefficients on a T grid, the mixing done exactly per cell). It removes the polynomials and divisions but replaces them with as many interpolations, about the same work at 10 species. Kept as the next step only if the binary polynomials remain the largest part after the change below.

Adopted instead: the same formulas with the divisions taken out of the inner loops. Reciprocals of the species' viscosity roots, of the Wilke denominators, of the molar masses and of each binary coefficient are formed once per call, and the NASA-7 enthalpy coefficients are divided by 2 to 5 once, when the medium is built. The model is unchanged; only the rounding changes. (An earlier draft of this paragraph said the NASA-7 change keeps the enthalpy bit for bit. It does not: the old polynomial's innermost term is `t*a4/5`, which rounds `t*a4` before dividing, while `t*(a4/5)` rounds the quotient first; the other quotients `a1/2` to `a3/4` are the same numbers either way. Corrected before any test ran; T1 below states the bound that replaces it.) The gradients gather their fields into one contiguous array per cell instead of calling the per-field accessor for every neighbour, with the same operations in the same order, so the result is unchanged bit for bit.

The larger levers stay where the plan put them: fewer cells (AMR) and threading that pays at the RL10's size.

## Criteria (stated before code)

- **T1 (same formulas).** On 10^5 random states (T from 200 to 6,000 K, p from 1e3 to 2e7 Pa, compositions drawn at random including near-pure gases and exact zeros), the new `Medium::transport` agrees with the previous implementation, kept in the test as the reference, within 1e-13 relative in mu, lambda and every D_km. `Medium::properties`, `enthalpy` and `speciesEnthalpies` agree within 1e-14 of cv T (cp T for the enthalpies), `cv` within 1e-14 relative, and the e-to-T inversion within 1e-12 relative (amended from "bit for bit" before the test ran; see the decision).
- **T2 (thread independence).** The turbulent C1 to 0.3 ms is the same bit for bit on 1 and 4 threads.
- **T3 (round-off only).** The turbulent C1 64x12 to 1 ms agrees with the previous build within 1e-9 relative in the injector pressure, the outlet mass flow and the vacuum thrust. The gradient change alone is checked bit for bit against the previous build.
- **T4 (cost).** The old and new property functions' cost per call (microbenchmark on the same states) and the one-thread step cost of the turbulent C1 to 1 ms are reported. No threshold: it is a measurement, and the change is kept only if the step is faster.

## Results (13 October 2026)

Raw outputs: `transport_cost/medium_cost_2026-10-13.txt` (T1, T4 microbenchmark), `transport_cost/state_comparison_2026-10-13.txt` (T2, T3), `transport_cost/step_timing_2026-10-13.txt` (T4 step).

- **T1 passes** (measured, `crucible_medium_cost`, 10^5 states, 10 species). Largest relative differences: mu 5.7e-16, lambda 4.3e-16, D_km 1.0e-15 (bound 1e-13). e 3.4e-15 of cv T, cv 0 (unchanged bit for bit), h 1.8e-15 of cp T, h_k 2.9e-15 of cp_k T (bound 1e-14). T from e 3.8e-15 relative (bound 1e-12).
- **T2 passes** (measured). The turbulent C1 64x12 to 0.3 ms gives byte-identical raw state dumps on 1 and 4 threads with the new build, as with the previous one.
- **T3 passes.**
  - Measured: at 1 ms (17,048 steps) the new state differs from the previous build's by at most 1.5e-13 of each field's largest magnitude (worst: a minor species), and 8.8e-14 at 0.3 ms. In the bulk fields the gap is at most 2.5e-14 of the field's largest magnitude (rho E). The largest pointwise relative difference is 5e-11, in rho u_r, where it is near zero.
  - Measured: every printed digit of the last history row is equal. That covers the injector pressure 2.750266e6 Pa, the outlet mass flow 0.3248790 kg/s and the vacuum thrust 1,381.944 N, as well as the wall heat flow, -2.878696109e5 W (ten digits).
  - Inferred, not printed: the three judged quantities agree within 1e-9. The history prints them to 7 digits, so this rests on the state bound above.
  - Measured: the gradient change alone (the "grad" build against the medium-only build) leaves the state byte-identical at 0.3 ms and 1 ms.
- **T4** (measured; Mac shared with other jobs, load average 7.6 to 11).
  - Microbenchmark, per call: transport 206.0 to 187.2 ns (1.10x); `properties` 16.0 to 17.6 ns (0.91x, within the run-to-run scatter of a 16 ns call); the e-to-T inversion from a 1% guess 134.6 to 110.7 ns (1.22x). An earlier run the same evening gave 1.17x for transport, so the microbenchmark moves by several percent between runs.
  - Step, one thread, C1 64x12 to 1 ms, the two builds alternated twice: previous build 80.03 and 90.78 s wall (78.78 and 84.93 s user); new build 73.03 and 72.72 s (69.81 and 70.92 s user). By mean user time the step is 1.16 times faster: 6.25 to 5.37 us per cell update. The change is kept.
- **Not measured.** No new profile was taken after the change. The table in T alone stays unadopted until a profile shows the binary-diffusion polynomials as the largest part of the step.
- **Other checks** (ctest, the same evening): core, wall function, thermo against CEA, reaction, turbulence, wall-function threads, PaSR and mixture verification pass. transport_verification fails only its two known radial-order criteria, at the values recorded before this change (1.218 and 1.444 against 1.8). The slow core and pipe tests were not run.
