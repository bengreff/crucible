# Stiff reaction integration: verification

Item 4, first piece (2026-09-30). TECHNICAL_PLAN: Cantera supplies local rate and property
evaluations, one context per thread; CRUCIBLE owns time advancement with CVODE variable-order
BDF over a consistent species/energy system; verify elemental conservation and thermodynamic
recovery.

## What exists

- `adapters/reaction.hpp/.cpp`
  - `ReactionSource`: Cantera rates for a closed, adiabatic, constant-volume parcel.
    - State z = [T, Y_1..Y_K] at fixed density.
    - dY_k/dt = w_k W_k / rho.
    - dT/dt = -sum(u_k w_k) / (rho cv), with u_k the partial molar internal energies.
    - Mass fractions are set without renormalising, so the integrator sees its own iterates.
  - `ReactionStep`: CVODES 5.3 BDF from the Cantera build, dense direct solve, difference-quotient
    Jacobian.
    - Reinitialised per step (the use pattern inside a split flow step).
    - Optional sampled output from one uninterrupted integration.
- `tests/reaction_tests.cpp` (ctest `reaction_verification`, under 1 s).

## Test

- Reference: Cantera's own `IdealGasReactor` + `ReactorNet` on an independent Solution, at rtol
  1e-12, atol 1e-20.
- CRUCIBLE integration: rtol 1e-8, atol 1e-14.
- Both are sampled on the same 200 000-point grid to t_end = 20 ms.
- Tolerances were written in the test header before the first run.

| check | tolerance | H2/O2 1000 K 1 atm (h2o2) | CH4/O2 1400 K 20 bar (gri30) |
|---|---|---|---|
| ignition delay (T0+400 K), rel. | 1e-3 | 2.6e-6 | 6.6e-8 |
| T(t_end), rel. | 1e-4 | 6.7e-10 | 4.2e-11 |
| max abs Y difference at t_end | 1e-6 | 1.3e-9 | 3.4e-11 |
| element mass fraction drift, rel. | 1e-8 | 2.0e-13 | 2.2e-12 |
| internal energy drift / (cv0 T0) | 1e-6 | 2.4e-8 | 5.5e-10 |
| T_end vs UV equilibrium [K] | 0.5 | 9e-13 | 8.1e-6 |
| max abs Y difference vs UV equilibrium | 1e-4 | 1.8e-15 | 9.8e-10 |

All values are measured.

- Ignition delays: 1.6313e-4 s (H2/O2) and 4.9755e-5 s (CH4/O2).
- End states: 3378.1 K and 3985.4 K.
- CVODES took 1093 and 1630 steps over the full 20 ms.
- Raw output: `docs/evidence/reaction_verification_2026-09-30.txt`.

## What this does not show

- Coupling to transport and the splitting error. Those need the multi-species flow state and are
  the next piece.
- Per-cell cost at flow time steps. The 20 ms single-parcel runs above are not a flow benchmark.
