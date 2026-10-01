# Multi-species core: thermally perfect mixture

Item 4, second piece (2026-09-30).

## What changed

- `core/medium.hpp/.cpp`: an ideal-gas mixture of thermally perfect species.
  - Each species carries NASA 7-coefficient polynomials.
  - The calorically perfect gas is now one species with constant cp. The core has one gas law
    with no per-material branch.
  - The temperature from internal energy uses a safeguarded Newton solve.
- `core/flow`:
  - The state carries partial densities rho*Y_k next to (rho, rho u_z, rho u_r, E).
  - Mass fractions are reconstructed like density: limited slopes, r^2 curvature on the axis row,
    and a one-sided limited slope at boundaries.
  - Face thermodynamics (e, frozen a, energy floor) come from the reconstructed (rho, p, Y).
  - HLLC is unchanged except that it takes these face values instead of a gamma.
  - Species fluxes follow the mass flux with the upwind face composition (Larrouturou, JCP 95,
    1991). This preserves positivity, and the species fluxes sum to the mass flux.
- `ReactionSource::medium()` exports a Cantera mechanism's species to the core, so the flow and the
  chemistry use the same thermo data.
- Inlet and outlet boundary models:
  - The reservoir inlet and the subsonic outlet use the frozen gamma of the local composition.
  - This is exact for a calorically perfect gas and an approximation for a thermally perfect one.
  - The TUM case needs injector mass-flow inlets, which come with the injector work.
- Build: now Release by default.
  - Before this, the build had no optimisation flags.
  - The old core ran its fast suite 12x faster once optimised (64 s to 5.6 s). Earlier wall-clock
    figures were unoptimised.

## Results (measured)

Regression check, fast core suite, single perfect gas:
- Every reference number is unchanged: Sod L1 0.00868557/0.00345083, nozzle mdot 0.879089, thrust
  ratio, entropy order, coil-force residuals.
- Only round-off-level diagnostics differ (planar-row deviations of order 1e-14).
- Cost: 12.1 s against 5.6 s for the old core, both optimised. That is 2.2x for single-gas runs.
- A bug in the temperature solve that first made this 4x was found by profiling and fixed. The
  Newton step landing exactly on the root was rejected by the bracket test.

`tests/mixture_tests.cpp`, tolerances written before the run:

| check | tolerance | result |
|---|---|---|
| core NASA-7 thermo vs Cantera (gri30; R, e, cv, a; 5 mixtures, 300-3500 K) | 1e-12 rel. | 6.3e-15 |
| temperature from Cantera's internal energy | 1e-9 K | 3.1e-10 K |
| two-gas shock tube (gamma 1.4 / 5/3) density L1, 400 vs 100 cells | ratio <= 0.5 | 0.386 (0.00811 to 0.00313) |
| star pressure, contact-shock midpoint, 400 cells | 1% | 0.49% |
| species masses sum to total mass | 1e-12 | 1.1e-16 |

The first thermo run failed at 5.7e-10. The cause was a truncated universal gas constant
(8314.462618 instead of the exact SI 8314.46261815324). It was fixed, and the tolerance was not
changed.

Reported, with no pass mark:
- Across the whole star region of the shock tube, the maximum |p/p* - 1| is 1.8%. This includes the
  contact and the smeared waves.
- A variable-gamma contact advected at uniform p and u gives spurious pressure of 4.6e-4 (100 cells)
  and 1.9e-4 (400 cells), and spurious velocity of 1.1e-3 and 4.8e-4. Both fall with resolution.
- This error is the known property of fully conservative multi-component schemes. A correction
  such as double flux (Abgrall and Karni 2001; Ma, Lv and Ihme 2017) is deferred until a reacting
  case shows it matters at the resolution used.

Evidence: `docs/evidence/mixture_verification_2026-09-30.txt`, `docs/evidence/core_suite_mixture_core_2026-09-30.txt`.
