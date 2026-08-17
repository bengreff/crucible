# CRUCIBLE Station 2 Certificate — the De Laval Nozzle

Goal-B station 2 (SOLV-1 §6-4; SOLV-7 §6-1; META-3 `maccormack-nozzle` lineage): a stagnation reservoir (p0 = 1, ρ0 = 1) feeds an axisymmetric converging–diverging nozzle — parabolic radius R(z) = R*(1 + 0.25((z−3)/3)²), R* = 1, length 6, max wall slope 9.5° — on the one cylindrical operator at N_θ = 1, axis inside the domain. The flow **chokes at the throat on its own** and exits supersonic; nothing about the operating point is imposed (SOLV-7 §3.2). Oracle: quasi-1-D isentropic theory (area–Mach relation, choked ṁ, vacuum C_F).

Geometry enters through the grid's config-time activity seam in its binary (stair-step) degenerate form, with **slip-ghost walls**: ghost states mirror about the true contour normal (ghost-cell immersed boundary), cutting spurious wave generation from O(wall slope) to O(h·curvature). Cost of that trade: a small stair-face transpiration flux, reported below as the plane-ṁ spread, shrinking with resolution — retired when FND-3's partial apertures + cut cells land. The quasi-1-D oracle also carries its own 2-D model floor (centerline ≠ area mean); both gaps are inside the declared bands. Criteria are the named constants in `crates/solvers/src/station2_nozzle.rs`, asserted by `crates/solvers/tests/solv1_station2_nozzle.rs`.

## Choked-flow ladder (marched to t = 25)

| n_r×n_z | steps | Cd = ṁ/ṁ_ideal | p_c/p0 | steady resid | max Mach dev (z > 1) | exit M (1-D 1.906) | C_F | ideal C_F | ṁ spread |
|---|---|---|---|---|---|---|---|---|---|
| 16×96 | 3632 | 1.0060 | 0.99999 | 1.69e-2 | 0.040 | 1.821 | 1.339 | 1.399 | 0.072 |
| 24×144 | 5479 | 1.0022 | 0.99999 | 5.28e-3 | 0.036 | 1.858 | 1.413 | 1.411 | 0.052 |
| 32×192 | 7325 | 1.0020 | 0.99999 | 7.60e-4 | 0.042 | 1.849 | 1.346 | 1.399 | 0.038 |

**Criteria:** Cd in [0.985, 1.015] at every level and |Cd−1| < 0.008 at the finest (measured 0.0020 — the choked mass flow matches ideal theory to 0.2%), never growing under refinement; p_c/p0 in [0.998, 1.001]; steadiness residual < 0.03 (finest < 0.003); centerline Mach within 0.06 of the area–Mach relation past the entrance band; exit Mach supersonic and within 0.08 of the 1-D exit value; C_F within 0.06 of ideal vacuum C_F; plane-ṁ spread < 0.1 and strictly shrinking (the declared slip-wall transpiration).

## Centerline Mach profile at 32×192

| z | A/A* | M (computed) | M (quasi-1-D) | rel. dev |
|---|---|---|---|---|
| 0.02 | 1.556 | 0.4687 | 0.4109 | +0.1407 |
| 0.27 | 1.459 | 0.4769 | 0.4461 | +0.0691 |
| 0.52 | 1.372 | 0.4982 | 0.4837 | +0.0298 |
| 0.77 | 1.297 | 0.5294 | 0.5240 | +0.0103 |
| 1.02 | 1.231 | 0.5673 | 0.5668 | +0.0009 |
| 1.27 | 1.174 | 0.6096 | 0.6123 | -0.0043 |
| 1.52 | 1.126 | 0.6563 | 0.6604 | -0.0062 |
| 1.77 | 1.086 | 0.7055 | 0.7112 | -0.0081 |
| 2.02 | 1.055 | 0.7568 | 0.7647 | -0.0103 |
| 2.27 | 1.030 | 0.8122 | 0.8208 | -0.0105 |
| 2.52 | 1.013 | 0.8679 | 0.8794 | -0.0131 |
| 2.77 | 1.003 | 0.9272 | 0.9405 | -0.0142 |
| 3.02 | 1.000 | 0.9879 | 1.0040 | -0.0161 |
| 3.27 | 1.004 | 1.0509 | 1.0699 | -0.0177 |
| 3.52 | 1.015 | 1.1143 | 1.1379 | -0.0207 |
| 3.77 | 1.033 | 1.1816 | 1.2080 | -0.0218 |
| 4.02 | 1.058 | 1.2512 | 1.2801 | -0.0225 |
| 4.27 | 1.091 | 1.3209 | 1.3540 | -0.0244 |
| 4.52 | 1.132 | 1.3922 | 1.4295 | -0.0261 |
| 4.77 | 1.181 | 1.4675 | 1.5066 | -0.0260 |
| 5.02 | 1.238 | 1.5451 | 1.5852 | -0.0253 |
| 5.27 | 1.306 | 1.6243 | 1.6650 | -0.0244 |
| 5.52 | 1.382 | 1.6722 | 1.7460 | -0.0423 |
| 5.77 | 1.470 | 1.7645 | 1.8281 | -0.0347 |

The entrance band (z < 1) carries the first-order stagnation-BC adjacency + 2-D entrance turning (decays from ~+14% at the first cell to ~+1%); it is reported, excluded from the pointwise gate, and shrinks with the COUP-7 injector-object wave.

## Masked-sweep exactness & determinism

- **Uniform gas at rest in the stair-stepped cavity is a bitwise fixed point** — grid-aligned mirror: true; slip-ghost wall: true (reflecting a zero velocity about any normal is the identity, so the wall machinery adds no drift).
- **Determinism:** reruns are byte-identical (asserted in the test suite on the full conserved state).
