# C2 step 1: molecular transport (viscosity, conduction, mixture-averaged diffusion)

Status: properties verified; operator partly verified. Criteria stated 4 October 2026 in the header of `tests/transport_tests.cpp`, before the first run. Two criteria still fail and are reported as failures below. Raw output: `transport_verification_2026-10-04.txt`; after the open-face fix, `transport_verification_openface_2026-10-04.txt` (last section).

## What was added

- **Properties** (`Medium::setTransport`, `Medium::transport`): Cantera's mixture-averaged model evaluated in the core from exported fits. Wilke viscosity, the averaged conductivity and D_km.
- **Operator** (`core/transport.cpp`): viscous stress with the hoop term, Fourier conduction plus species enthalpy flux, and mixture-averaged diffusion with the correction velocity (Soret off).
  - Gradients are least squares over the four neighbours of each cell.
  - Face gradients are the mean of the two cell gradients, with the component along the line of centres replaced by the compact difference.
  - Walls: no-slip or slip, isothermal or adiabatic, impermeable to species. Wall shear and heat enter the thrust and energy ledgers.
  - Open boundaries (nozzle inlet, outlet) use a zero normal gradient. Supply rings exchange no diffusive flux (plug inflow, declared).
- **Time step**: the explicit diffusive limit is added to the convective one.

## Properties against Cantera (reaction_verification)

Measured: the worst relative difference from Cantera's `MixTransport` over 5 mixtures, 250 to 3400 K and 1e4 to 2e7 Pa, for h2o2.yaml and gri30.yaml, is at most 8.8e-14 for viscosity, conductivity and D_km (limit 1e-12). Pass.

## Operator results (transport_verification)

| Check | Criterion | Measured | Verdict |
|---|---|---|---|
| Pipe decay (Bessel J0 profile, one e-fold), observed order nr 16 to 32 | ≥ 1.8 | 1.988 | pass |
| Pipe decay L1 error at nr 32 | < 1e-3 | 1.94e-4 | pass |
| Pipe decay momentum balance (wall shear) | < 1e-12 | 1.4e-18 | pass |
| Chamber budgets after 200 steps: mass, energy with wall heat, momentum with wall shear | < 1e-12 | 6.8e-16, 8.0e-16, 2.6e-17 | pass |
| Operator order, axial momentum (128x32 to 256x64) | ≥ 1.8 | 1.799 | **fail** (by 0.001) |
| Operator order, radial momentum | ≥ 1.8 | 1.189 | **fail** |
| Operator order, energy | ≥ 1.8 | 1.929 | pass |
| Operator order, species | ≥ 1.8 | 1.928 | pass |

The operator-order test measures truncation error: the discrete divergence of the transport flux against the exact cell average, on a contoured duct whose wall slope reaches 0.2. The pipe decay measures the solution error, which is what a simulation actually delivers.

## What was fixed on the way (in order, measured)

1. **Open ends carried no stress.** In the first pipe-decay run every grid had an L1 error of 0.54, and the decay was about 5 times too slow.
   - Mechanism: the end cells got no shear through their open faces, while interior faces carried mu du_z/dr. That imbalance drove a radial flow of d(rho u_r)/dt about 56 at the ends, against about 1e-3 inside.
   - Fix: open boundaries take a zero normal gradient. Only the stress of the tangential gradients crosses them.
2. **The first face off the axis.** The compact difference across that face estimates the slope at the mean of the two centroid radii, 1.11 dr, not at the face, dr. That is an 11% error for an even field, and it does not shrink with refinement.
   - Fix: even fields are differenced in r^2 with the cells' exact second moments, which is exact for a + b r^2.
   - In the same place, u_r/r is now averaged to the face as itself, rather than as the mean u_r divided by the face radius.
   - Result: the pipe decay converges at order 1.99.
3. **Test correction (the operator was not changed).** The first species fields varied as s^2 at the wall, so they carried a species flux through a wall that the operator makes impermeable.
   - The reference integrated that flux, and the wall-row gradient fit assumed none. Energy and species errors therefore grew toward the wall.
   - Measured with the s^2 fields: energy order 0.50, species order 0.28. Wall-row species error 18, 34, 61, 113 on the four grids.
   - The composition now varies as 1 - (1 - s^2)^2, which is flat at the wall along both directions. Energy and species then converge at 1.93.
   - The criterion and its threshold are unchanged. The header of the test records the change.

## The two failures: mechanism as far as measured

A second fix was tried and had no effect. Following the two-failed-fixes rule, I stopped there and wrote the mechanism down instead.
- The fix: the axis ring's least-squares gradient of even fields was refitted on (z, <r^2>).
- It was reverted.

**Axial momentum (1.799).**
- The error sits in the row next to the wall. Per-row relative error, rows nr-2: 1.5e-2 at 64x16 and 1.3e-2 at 128x32. It is not falling.
- The wall row itself, which is excluded from the norm, has relative errors of 1.5, 3.6, 9.1 and 8.8 on the four grids of the contoured duct. It grows over the first three grids and then holds. On a straight duct it stays constant (3.2 and 3.2 at 64x16 and 128x32).
- So the wall shear carries an error that depends on the wall slope and does not fall with refinement, and it leaks into the next row through the shared face.
- Not yet located.

**Radial momentum (1.19; 1.52 on a straight duct).**
- The error is spread over the rows near the axis. At a fixed radius it falls only about 2.5 times per halving: row 4 at 64x16 is 6.2e-3, and row 8 at 128x32 is 2.4e-3.
- The axis ring has an O(1) relative error, 0.30 then 0.41. The exact radial-momentum divergence vanishes on the axis, so relative errors there are inflated.
- The second ring has a constant 3% error in every field.
  - Derived mechanism: the face above it still uses the plain compact difference, whose centroid offset (about h^2 / 12 r) no longer cancels against the face below, which is now exact.
  - A candidate fix is the r^2 difference on every radial face, but it is not tried yet.

**What this means for use.**
- The solution-level test (pipe decay) converges at second order. Energy and species converge at second order in truncation.
- The momentum truncation error converges at first to second order near the axis and the sloped wall. That is typical for cell-centred finite volumes with least-squares gradients, where the solution can still converge faster than the truncation error. But it is not what I stated, and it stays a failure until it is fixed or explained by a solution-level test on a sloped wall.
- The FreeFlame comparison (next) exercises species and energy, not these two terms.

## Wall-clustered rings (check 4, `Definition::radialStretching`)

Added for the near-wall resolution SST needs (y+ about 1). Ring boundary j sits at tanh(b s) / tanh(b) of the local radius, s = j / nr. Zero keeps equal rings, with the old expressions taken bit for bit. Criteria were stated in the header of `tests/transport_tests.cpp` before the first run. Raw output: `transport_verification_stretch_2026-10-04.txt`.

- **Corrected before any run used it.** The first form, 1 - tanh(b (1 - s)) / tanh(b), has slope b / tanh(b) > 1 at the wall, so it clustered rings toward the axis. Caught on reading; no result came from it.
- **Unchanged on equal rings (measured).** Checks 1 to 3 give output identical to `transport_verification_2026-10-04.txt`. The core checks (`crucible_tests`), built from the code of 807701b (which also contains SST, off in these checks), give output identical to `core_verification_outlet_2026-10-04.txt`.

At b = 2 the wall ring is 0.166 and the axis ring 2.06 of the equal height (64 x 16).

| Check | Criterion | Measured |
|---|---|---|
| 4a. cell volumes against the exact duct volume | < 1e-13 | 1.8e-15 |
| 4a. axial-face areas against pi R^2, worst column | < 1e-13 | 3.3e-16 |
| 4b. operator order, axial momentum | >= 1.8 | 1.879 |
| 4b. radial momentum | >= 1.8 | 1.287, **FAIL** |
| 4b. energy | >= 1.8 | 1.885 |
| 4b. species | >= 1.8 | 1.902 |
| 4c. pipe decay order, nr 16 to 32 | >= 1.8 | 1.692, **FAIL** |
| 4c. pipe decay L1 error at nr 32 | < 1e-3 | 1.89e-4 |
| 4c. momentum balance | < 1e-12 | 1.0e-17 |

- **4b radial momentum** fails as it does on equal rings (1.189 there): the near-axis mechanism above. Stretching does not add a new failure mode in the operator.
- **4c.** The L1 errors at nr 8, 16 and 32 are 3.17e-3, 6.10e-4 and 1.89e-4 (orders 2.38, then 1.69). On equal rings they are 3.01e-3, 7.69e-4 and 1.94e-4 (1.97, 1.99), so the absolute error at nr 32 is the same. A non-monotone order suggested two error parts of opposite sign. A scratch diagnostic (old code) also measured the error against the exact cell average: 6.21e-3, 1.47e-3 and 3.67e-4 (orders 2.08, 2.00). It was stopped before nr 64 when the cause below was found. The open-face fix (next section) makes 4c pass.

## Open-face neighbour fix (1e0fd9e) and the rerun

Found by the SST fully developed pipe (`docs/evidence/TURBULENCE_C2.md`), where the bulk velocity climbed from 105 to 136 m/s in 1.5 ms on an axially uniform flow.
- **Mechanism.** In the least-squares gradient, an open or supply face (zero normal gradient) contributed a neighbour at the face midpoint carrying the cell's own value. The midpoint sits at the ring's middle radius, not at the centroid's, so that neighbour also claimed a zero radial derivative over the offset between the two. It biased the radial gradients of the end columns, and through the face-gradient average the next column too.
- **Fix.** The neighbour sits at the face on the cell's axial line (at the centroid's radius), where a zero normal gradient puts the cell's own value. The injector plate (a wall) keeps the face midpoint with the wall value.
- **Measured.** The SST pipe stays axially uniform to 13 digits and the drift is gone. On the straight slip-wall duct of the SST check 3, the laminar radial-momentum order went from 1.589 to 1.953.

Rerun of every transport check on the fixed code. Raw output: `transport_verification_openface_2026-10-04.txt`.

| Check | Before | After | Verdict after |
|---|---|---|---|
| 1. operator order, axial momentum (contoured duct, equal rings) | 1.799 | 1.799 | **fail** (by 0.001) |
| 1. radial momentum | 1.189 | 1.223 | **fail** |
| 1. energy, species | 1.929, 1.928 | 1.929, 1.928 | pass |
| 2. pipe decay order, nr 16 to 32 | 1.988 | 1.992 | pass |
| 2. pipe decay L1 at nr 32 | 1.94e-4 | 1.64e-4 | pass |
| 3. chamber budgets: mass, energy, momentum | 6.8e-16, 8.0e-16, 2.6e-17 | 1.4e-16, 3.5e-16, 2.9e-17 | pass |
| 4b. radial momentum, b = 2 | 1.287 | 1.446 | **fail** |
| 4b. axial momentum, energy, species, b = 2 | 1.879, 1.885, 1.902 | 1.881, 1.885, 1.902 | pass |
| 4c. pipe decay order, b = 2, nr 16 to 32 | 1.692 | 2.009 | pass |
| 4c. pipe decay L1 at nr 32 | 1.89e-4 | 3.70e-4 | pass |
| 4c. momentum balance | 1.0e-17 | 3.1e-17 | pass |

- 4c now converges cleanly: 6.29e-3, 1.49e-3 and 3.70e-4 at nr 8, 16 and 32 (orders 2.08, 2.01). The absolute error at nr 32 doubled. The old, smaller error came partly from the end-column bias cancelling some of the scheme's own error; that is inferred from the order becoming regular, not shown separately.
- 4c's energy balance at nr 32 is 1.05e-12 after 600 308 steps (reported, not a criterion; about 1.7e-18 per step, round-off).
- Three failures remain, down from four: the axial-momentum order on the contoured duct (unchanged, the wall-row mechanism above) and the radial-momentum order on equal and clustered rings. The radial part is now the near-axis and sloped-wall error only. The straight-duct part was largely the end-column bias.
