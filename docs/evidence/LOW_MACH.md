# Low-Mach accuracy of the HLLC core (item 3, 2026-09-30)

## Question

The core venturi test shows -9.2% (40x6) and -2.6% (80x12) outlet mass flow at exit Mach 0.147 from rest
(**measured**). Combustion chambers run at Mach 0.01-0.3. Upwind dissipation in HLLC scales with the sound
speed c, while the physical pressure variations scale with rho u^2, so the relative error at low Mach grows
like 1/M (the standard Guillard-Viozat result). Should the core use a low-Mach correction?

## Candidates and why preconditioning is out

- **Thornber et al. 2008** (JCP 227:4873). Before the flux is computed, the left/right velocities at a face are
  moved toward their mean: the jump is scaled by z = min(1, max(M_L, M_R)). It has no tunable parameter.
- **HLLC-LM, Fleischmann, Adami & Adams 2020** (JCP 423:109762). HLLC is written as a central flux plus
  dissipation, and the two acoustic-wave terms are scaled by phi = sin(pi/2 min(1, M/0.1)). The 0.1 cutoff is
  the published value and was not tuned here.
- **Turkel/Weiss-Smith preconditioning: rejected without implementation.**
  - Applied with explicit time-accurate stepping, the stable time step shrinks like M^2 (Birken & Meister 2005).
  - Applied with dual time stepping, it replaces the physical march, which the startup mandate forbids.

Both implemented candidates sit behind `Definition::lowMach` (default `None`). They act on interior faces and
on the mirror wall; inlet and outlet boundary models keep their own states.

## Evidence

### Steady venturi accuracy

Source: `build/crucible_low_mach_study 0.04 40,80`. The run starts from the quasi-1D isentropic solution,
reports error at t = 0.04 s, and covers four unchoked exit Mach numbers. All values are **measured**.

| exit M | HLLC 40 / 80 | Thornber 40 / 80 | HLLC-LM 40 / 80 |
|---|---|---|---|
| 0.054 | -10.85% / -2.83% | -1.66% / -0.61% | -10.38% / -2.71% |
| 0.120 | -8.06% / -2.07% | -2.69% / -0.76% | identical to HLLC |
| 0.147 | -7.38% / -1.89% | -3.11% / -0.86% | identical to HLLC |
| 0.170 | -6.94% / -1.81% | -3.49% / -0.98% | identical to HLLC |

- **The runs are still drifting.** Between 0.03 and 0.04 s the error grows by 10-20% of its value, toward the
  from-rest numbers. So these are lower bounds on the settled error, and the ranking is what counts.
- **Thornber** cuts the error 2-6x, by most at the lowest Mach, and keeps the order near 1.8-1.9.
- **HLLC-LM** changes nothing once M > 0.1 anywhere near the throat, which is by design. Its cutoff could only
  be raised by tuning a published constant, which is an arbitrary choice made to fit a test.
- **Discarded design:** the first study, started from rest at p_b/p0 = 0.96, chokes the throat. Its -19% on both
  grids is the choked limit, not a dissipation error.

### Stability: Thornber fails

Case: planar Sod in a 5-row duct (`build/crucible_low_mach_stability`). Quantity: maximum
relative density difference between radial rows, which is exactly zero in exact arithmetic.

| scheme, nz | t = 0.05 | 0.10 | 0.15 | 0.20 | 0.30 |
|---|---|---|---|---|---|
| HLLC 400 | 4e-15 | 7e-15 | 8e-15 | 2e-14 | 3e-14 |
| Thornber 400 | 6e-12 | 6e-10 | 1.2e-8 | 4.8e-7 | 1.9e-6 |
| Thornber 800 | 4e-10 | 2.5e-7 | 2.7e-5 | 2.4e-5 | 7.2e-5 |

- The difference grows exponentially from roundoff, at about 100x per 0.05 time units, then saturates. It grows
  faster on the finer grid.
- It sits in the rarefaction fan near its head, a low-Mach region, not at the shock.
- HLLC-LM shows the same kind of row growth in the core test (5e-8 at 400 cells).

## Mechanism (written down under the two-fix rule)

- In HLLC the term ~ rho c (jump in normal velocity) is the dissipation that couples pressure and velocity
  between neighbouring cells. It is what suppresses odd-even (checkerboard) modes.
- Both fixes scale that term by M. Across radial faces in a planar flow, the jump in normal velocity is pure
  perturbation, so its damping drops from O(rho c) to O(rho u).
- The remaining coupling comes from the axisymmetric geometric source p/r and the centroid-referenced radial
  reconstruction. Without the acoustic damping, that coupling amplifies a transverse mode instead of damping it.
- This is the known pressure-velocity decoupling of low-Mach fixes (Dellacherie 2010, checkerboard modes),
  showing up in the transverse direction. It is **inferred**: I have not yet isolated which term amplifies.
- A reacting chamber has low-Mach regions everywhere (recirculation, the injector face) and is resolved on fine
  grids, which is exactly where the growth is fastest. An unstable mode there would contaminate the measurements
  the validation compares.

## Decision (engineering, mine)

1. **Default stays HLLC with no low-Mach correction.** It is stable, verified, and second-order convergent on
   the venturi.
2. **Chamber accuracy is controlled by resolution.** The validation reports a measured grid-convergence error
   band for every observable. The venturi numbers above set the expected size: at M ~ 0.15, about 2% mass-flow
   error at 80 axial cells across a 3:1 area change, falling at order ~1.9.
3. **Thornber and HLLC-LM stay in code as non-default variants.** The test suite runs them with
   `crucible_tests --low-mach=thornber|hllclm`; there, the strict row-identity check is skipped and the radial
   velocity check fails, recording the instability.
4. **Revisit only with a stable low-Mach flux.** Candidates are an all-speed scheme with proven checkerboard
   suppression (e.g. AUSM+-up, Liou 2006, whose pressure diffusion term targets exactly this coupling), or
   Rieper's LM-Roe (normal-jump-only scaling) with an explicit stability test like the one above. Either must
   pass the planar-row growth test at 800 cells to t = 0.45 before it is considered.
