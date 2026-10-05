# SST-2003 pipe reference (independent 1-D solver)

Status: reference built and grid-converged. It is the target for the engine's periodic-pipe verification (TECHNICAL_PLAN step 7, verification 4). It is not a validation: no measured data is involved.

## What it is

`tools/sst_pipe_1d.py` solves steady, fully developed, compressible pipe flow under SST-2003. It shares no code with the engine.
- **Equations.** Axial momentum, k, omega and total energy (with k in the energy, the shear work and the k diffusion flux). The radial momentum gives p + 2/3 rho k uniform; p = rho R T; the mass per unit length is fixed by the fill. Constants, blending, the production limiter in both equations, mu_t with S F2, and the wall omega = 10 * 6 nu / (beta1 d1^2) all follow TMR `sst.html`, as in the TECHNICAL_PLAN design.
- **Properties.** Read from Cantera directly (`tools/n2_properties.cpp` writes `tools/n2_properties.csv`: Cantera 3.2.0, h2o2.yaml, mixture-averaged, pure N2). The engine uses its own fits of the same data.
- **Numerics.**
  - Finite volumes on rings r_j = R tanh(b j / N) / tanh(b) (the engine's wall-clustered family).
  - Coarse start: a segregated iteration on 50 rings. Each grid then gets Newton's method, with a sparse finite-difference Jacobian (5-colour stencil, analytic mass row) and pseudo-time continuation, starting from the previous grid's solution.
  - Converged to relative updates of 1e-12 or below on every grid.

## Case

N2 filled at 101325 Pa and 300 K, R = 0.4 mm, isothermal wall at 300 K, body force 2.95e5 N/m^3, b = 2, Pr_t = 0.9.
- Chosen so that Re_tau is near 180, the lowest pipe-DNS point for later validation, and the bulk Mach number is about 0.28. That keeps the engine's low-Mach HLLC dissipation out of the comparison.
- All numbers below are measured outputs of the solver unless marked.

## Results (`sst/reference_sweep_2026-10-04.txt`)

| N | y1+ | c_f | u_b (m/s) | T_axis (K) | force balance | energy balance |
|---|---|---|---|---|---|---|
| 100 | 0.136 | 0.0102696040 | 100.48379 | 304.36102 | 0 | 1.8e-13 |
| 200 | 0.068 | 0.0103798730 | 99.94863 | 304.31688 | 0 | 2.3e-13 |
| 400 | 0.034 | 0.0104342341 | 99.68793 | 304.29561 | 0 | 3.2e-13 |
| 800 | 0.017 | 0.0104607894 | 99.56132 | 304.28534 | 0 | 2.8e-12 |
| 1600 | 0.0084 | 0.0104738127 | 99.49940 | 304.28033 | -1.1e-16 | 3.3e-12 |
| 3200 | 0.0042 | 0.0104802390 | 99.46889 | 304.27787 | 0 | -1.7e-12 |

Re_tau is 182.43 throughout; c_f is 2 tau_w / (rho_b u_b^2) with the mean density.

**First-order convergence.**
- The observed order is 1.02 to 1.03 for c_f and u_b on every triple of grids, steady from 100 to 3200 rings.
- Richardson limits (derived), from 800-1600-3200: c_f 0.0104865, u_b 99.439 m/s, T_axis 304.2755 K. Successive triples agree to 2e-5 relative in c_f.

**The wall omega is the first-order term.**
- With the wall factor 5, 10 and 20, the finite-grid c_f at 100 rings is 0.0103262, 0.0102696 and 0.0102193: a 1.0% spread.
- The three Richardson limits (100-200-400) are 0.0104868, 0.0104871 and 0.0104870, and u_b agrees to 0.003%.
- Mechanism (inferred): the wall value 60 nu / (beta1 d1^2) is tied to the first-centroid distance d1, so its error scales with d1, which is proportional to 1/N here. The limit does not depend on the factor.

## Consequence for the engine test

- The engine applies the same wall omega, so it will also converge at first order. At practical ring counts (y1+ 0.1 to 0.3) the boundary-condition error in c_f is 1 to 3% (measured here at 100 rings: -2.07% against the limit).
- The engine therefore has to be compared on matched rings (the same N and b, where this term is common) and through its own extrapolated limit against this reference's limit. The test header will state both criteria before the first run.
