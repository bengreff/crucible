# CRUCIBLE Station 1 Certificate — the Bursting Diaphragm (Sod)

Goal-B station 1 (SOLV-1 §6-1; VAL-2 §3.3 `sod-shock`): compressible Euler on the unified cylindrical grid at N_θ = 1 — PPM reconstruction, HLLC flux with Batten wavespeeds, well-balanced geometric sources, one passive advected composition riding the contact. The tube is a full cylinder (r = 0 axis inside the domain), radially uniform, marched by the explicit MOL SSP-RK2 reference integrator (superseded by COUP-3's SDC-IMEX when it lands; the flux-form spatial operator carries over unchanged). Oracle: the exact Riemann solution (Toro exact solver; star state verified against Toro Table 4.2 in the unit tests). Criteria are the named constants in `crates/solvers/src/station1_sod.rs`, asserted by `crates/solvers/tests/solv1_station1_sod.rs`, and enforced against this committed file by gate 4 of `scripts/check.sh`.

## Sod tube vs exact Riemann solution (t = 0.2)

| n_z | steps | L1(ρ) global | order | L1(ρ) fan (0.3, 0.45) | order | L1(ρ) star-L (0.52, 0.66) | order |
|---|---|---|---|---|---|---|---|
| 100 | 134 | 3.7927e-3 | — | 5.0186e-4 | — | 2.0205e-4 | — |
| 200 | 244 | 1.7242e-3 | 1.137 | 2.4089e-4 | 1.059 | 4.1631e-5 | 2.279 |
| 400 | 464 | 9.2220e-4 | 0.903 | 1.2229e-4 | 0.978 | 9.1607e-6 | 2.184 |
| 800 | 902 | 5.0115e-4 | 0.880 | 6.1555e-5 | 0.990 | 1.5442e-6 | 2.569 |

**Criteria:** global L1 falls ≥1.5× per refinement and is < 1e-3 at the finest level; star-plateau window order ≥ 1.8 (formal order in smooth regions — measured 2.2–2.6); fan-interior order ≥ 0.8 (centered-rarefaction startup singularity gives the known ~1st-order fan interior; gated against further degradation, reported honestly).

## Wave structure at n_z = 800

Exact star state: p* = 0.30313, u* = 0.92745, ρ*L = 0.42632, ρ*R = 0.26557 (Toro Table 4.2).

- **Shock position:** off by 0.155 cells — criterion ≤ 1 (VAL-2: shock within 1 cell).
- **Contact position** (composition midpoint): off by 1.108 cells — criterion ≤ 2 (the Batten-restored wave).
- **Star plateaus:** max errors ρ*L 4.66e-5, ρ*R 1.72e-5, u* 7.24e-5, p* 2.37e-5 — criterion ≤ 5e-4 each.
- **Composition bounds:** C ∈ [0.00e0, 1 + 4.44e-16] over the whole tube — criterion within [−1e-12, 1 + 1e-12] (species ride the contact, no new extrema).
- **Radial uniformity:** a radially-uniform tube stays radially uniform **bitwise** through the whole shock evolution: true (one law on the one grid; the axis metric and geometric sources cancel exactly).

## Smooth formal-order study (advected tanh ramp, width 0.1, t = 0.25)

| n_z | L1(ρ) | order | L1(C) | order |
|---|---|---|---|---|
| 50 | 4.8448e-5 | — | 1.8060e-4 | — |
| 100 | 6.8586e-6 | 2.820 | 3.9584e-5 | 2.190 |
| 200 | 2.2680e-6 | 1.596 | 9.9298e-6 | 1.995 |
| 400 | 6.2780e-7 | 1.853 | 2.4952e-6 | 1.993 |

**Criteria:** mean ρ order across the ladder ≥ 1.8 (measured 2.09; per-step ≥ 1.4 — ρ couples to the acoustic families and oscillates per step); composition order per step in [1.8, 2.3] (pure contact family — textbook 2nd order). Levels [50, 100, 200, 400].

## Whole-operator MMS (SOLV-1 §6-2)

A manufactured smooth field with radial flow, swirl, and axial flow all active — every flux direction and all three geometric source terms (pressure, centrifugal ρu_θ², swirl advection ρu_ru_θ) carry nonzero operands, which the Sod run structurally cannot exercise (its u_r ≡ 0). The analytic residual enters through the operator's source intake; the computed field must recover the manufactured one at formal order in **L1 per conserved component** (the shock-capturing verification norm: the limiter's clipping at the θ-mode's smooth extrema is locally 1st-order over an O(h) measure, which L2 amplifies to a ~O(h^1.6) tail while L1 retains the formal order — measured and recorded here as the honest caveat). Levels [16, 32, 64].

### Euler MMS 2-D axisymmetric with swirl (N_θ = 1)

| n | L1(rho) | order | L1(mom_r) | order | L1(mom_theta) | order | L1(mom_z) | order | L1(rho_e) | order | L1(rho_c) | order |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 16 | 9.41e-5 | — | 1.23e-4 | — | 3.20e-5 | — | 7.51e-5 | — | 3.77e-4 | — | 8.90e-5 | — |
| 32 | 2.55e-5 | 1.88 | 2.93e-5 | 2.07 | 7.55e-6 | 2.08 | 1.82e-5 | 2.05 | 9.95e-5 | 1.92 | 2.36e-5 | 1.92 |
| 64 | 6.37e-6 | 2.00 | 7.33e-6 | 2.00 | 1.90e-6 | 1.99 | 4.44e-6 | 2.03 | 2.48e-5 | 2.00 | 5.95e-6 | 1.99 |

**Criterion:** every component's observed order in [1.8, 2.3] (formal = 2) at every refinement.

### Euler MMS 3-D with m = 2 θ-mode (N_θ = n)

| n | L1(rho) | order | L1(mom_r) | order | L1(mom_theta) | order | L1(mom_z) | order | L1(rho_e) | order | L1(rho_c) | order |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 16 | 7.00e-4 | — | 2.27e-4 | — | 1.16e-3 | — | 3.37e-4 | — | 2.78e-3 | — | 6.84e-4 | — |
| 32 | 1.53e-4 | 2.20 | 5.88e-5 | 1.95 | 2.51e-4 | 2.21 | 7.97e-5 | 2.08 | 6.07e-4 | 2.19 | 1.56e-4 | 2.13 |
| 64 | 3.64e-5 | 2.07 | 1.48e-5 | 1.99 | 5.37e-5 | 2.22 | 1.78e-5 | 2.17 | 1.43e-4 | 2.09 | 3.47e-5 | 2.17 |

**Criterion:** every component's observed order in [1.8, 2.3] (formal = 2) at every refinement.

## Well-balance, conservation, determinism

- **Uniform gas at rest is a bitwise fixed point** over 50 steps — axis config (N_θ = 1, r = 0 inside): true; annulus config (N_θ = 8): true. The geometric sources use the same A·p products as the face fluxes, so the balance is exact, not approximate.
- **Closed reflecting tube through wall reflections** (t = 0.6): mass drift 3.927e-16, energy drift 3.213e-16 relative — criterion < 1e-12 each (flux-form telescoping; measured at round-off). Momentum is exchanged with the walls — that ledger is COUP-2 §3.1.2's mount-reaction term, a later wave.
- **Determinism:** reruns are byte-identical (asserted in the test suite on the full conserved state).
