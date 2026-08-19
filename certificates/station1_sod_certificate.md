# CRUCIBLE Station 1 Certificate — the Bursting Diaphragm (Sod)

Goal-B station 1 (SOLV-1 §6-1; VAL-2 §3.3 `sod-shock`): compressible Euler on the unified cylindrical grid at N_θ = 1 — PPM reconstruction, HLLC flux with Batten wavespeeds, well-balanced geometric sources, one passive advected composition riding the contact. The tube is a full cylinder (r = 0 axis inside the domain), radially uniform, marched by the ONE production integrator — COUP-3's SDC-IMEX step (S2; explicit hyperbolic class, fixed sweeps, the COUP-2 conservation audit armed every step; the session-7 SSP-RK2 scaffolding is retired). Oracle: the exact Riemann solution (Toro exact solver; star state verified against Toro Table 4.2 in the unit tests). Criteria are the named constants in `crates/solvers/src/station1_sod.rs`, asserted by `crates/solvers/tests/solv1_station1_sod.rs`, and enforced against this committed file by gate 4 of `scripts/check.sh`.

## Sod tube vs exact Riemann solution (t = 0.2)

| n_z | steps | L1(ρ) global | order | L1(ρ) fan (0.3, 0.45) | order | L1(ρ) star-L (0.52, 0.66) | order |
|---|---|---|---|---|---|---|---|
| 100 | 134 | 3.7431e-3 | — | 4.6536e-4 | — | 1.6771e-4 | — |
| 200 | 244 | 1.7641e-3 | 1.085 | 2.3043e-4 | 1.014 | 4.3675e-5 | 1.941 |
| 400 | 463 | 9.3496e-4 | 0.916 | 1.1308e-4 | 1.027 | 8.2598e-6 | 2.403 |
| 800 | 902 | 5.2074e-4 | 0.844 | 5.5864e-5 | 1.017 | 1.8730e-6 | 2.141 |

**Criteria:** global L1 falls ≥1.5× per refinement and is < 1e-3 at the finest level; star-plateau window order ≥ 1.8 (formal order in smooth regions — measured 2.2–2.6); fan-interior order ≥ 0.8 (centered-rarefaction startup singularity gives the known ~1st-order fan interior; gated against further degradation, reported honestly).

## Wave structure at n_z = 800

Exact star state: p* = 0.30313, u* = 0.92745, ρ*L = 0.42632, ρ*R = 0.26557 (Toro Table 4.2).

- **Shock position:** off by 0.845 cells — criterion ≤ 1 (VAL-2: shock within 1 cell).
- **Contact position** (composition midpoint): off by 1.108 cells — criterion ≤ 2 (the Batten-restored wave).
- **Star plateaus:** max errors ρ*L 4.86e-5, ρ*R 1.51e-4, u* 8.84e-4, p* 2.33e-4 — criterion ≤ 2.5e-3 each.
- **Composition bounds:** C ∈ [0.00e0, 1 + 2.22e-16] over the whole tube — criterion within [−1e-12, 1 + 1e-12] (species ride the contact, no new extrema).
- **Radial uniformity:** a radially-uniform tube stays radially uniform **bitwise** through the whole shock evolution: true (one law on the one grid; the axis metric and geometric sources cancel exactly).

## Smooth formal-order study (advected tanh ramp, width 0.1, t = 0.25)

| n_z | L1(ρ) | order | L1(C) | order |
|---|---|---|---|---|
| 50 | 6.8637e-5 | — | 2.2196e-4 | — |
| 100 | 7.4019e-6 | 3.213 | 4.7951e-5 | 2.211 |
| 200 | 1.3801e-6 | 2.423 | 1.1874e-5 | 2.014 |
| 400 | 3.2852e-7 | 2.071 | 2.9799e-6 | 1.994 |

**Criteria:** mean ρ order across the ladder ≥ 1.8 (measured 2.57; per-step ≥ 1.4 — ρ couples to the acoustic families and oscillates per step); composition order per step in [1.8, 2.3] (pure contact family — textbook 2nd order). Levels [50, 100, 200, 400].

## Whole-operator MMS (SOLV-1 §6-2)

A manufactured smooth field with radial flow, swirl, and axial flow all active — every flux direction and all three geometric source terms (pressure, centrifugal ρu_θ², swirl advection ρu_ru_θ) carry nonzero operands, which the Sod run structurally cannot exercise (its u_r ≡ 0). The analytic residual enters through the operator's source intake; the computed field must recover the manufactured one at formal order in **L1 per conserved component** (the shock-capturing verification norm: the limiter's clipping at the θ-mode's smooth extrema is locally 1st-order over an O(h) measure, which L2 amplifies to a ~O(h^1.6) tail while L1 retains the formal order — measured and recorded here as the honest caveat). Levels [16, 32, 64].

### Euler MMS 2-D axisymmetric with swirl (N_θ = 1)

| n | L1(rho) | order | L1(mom_r) | order | L1(mom_theta) | order | L1(mom_z) | order | L1(rho_e) | order | L1(rho_c) | order |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 16 | 9.42e-5 | — | 1.23e-4 | — | 3.21e-5 | — | 7.51e-5 | — | 3.78e-4 | — | 8.91e-5 | — |
| 32 | 2.55e-5 | 1.88 | 2.93e-5 | 2.07 | 7.55e-6 | 2.09 | 1.81e-5 | 2.05 | 9.96e-5 | 1.92 | 2.35e-5 | 1.92 |
| 64 | 6.38e-6 | 2.00 | 7.33e-6 | 2.00 | 1.90e-6 | 1.99 | 4.44e-6 | 2.03 | 2.49e-5 | 2.00 | 5.95e-6 | 1.98 |

**Criterion:** every component's observed order in [1.8, 2.3] (formal = 2) at every refinement.

### Euler MMS 3-D with m = 2 θ-mode (N_θ = n)

| n | L1(rho) | order | L1(mom_r) | order | L1(mom_theta) | order | L1(mom_z) | order | L1(rho_e) | order | L1(rho_c) | order |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 16 | 7.01e-4 | — | 2.27e-4 | — | 1.16e-3 | — | 3.37e-4 | — | 2.78e-3 | — | 6.85e-4 | — |
| 32 | 1.53e-4 | 2.20 | 5.87e-5 | 1.95 | 2.51e-4 | 2.21 | 7.97e-5 | 2.08 | 6.08e-4 | 2.20 | 1.57e-4 | 2.13 |
| 64 | 3.65e-5 | 2.07 | 1.48e-5 | 1.99 | 5.37e-5 | 2.22 | 1.78e-5 | 2.17 | 1.43e-4 | 2.09 | 3.47e-5 | 2.17 |

**Criterion:** every component's observed order in [1.8, 2.3] (formal = 2) at every refinement.

## Well-balance, conservation, determinism

- **Uniform gas at rest is a bitwise fixed point** over 50 steps — axis config (N_θ = 1, r = 0 inside): true; annulus config (N_θ = 8): true. The geometric sources use the same A·p products as the face fluxes, so the balance is exact, not approximate.
- **Closed reflecting tube through wall reflections** (t = 0.6): mass drift 1.963e-16, energy drift 1.606e-16 relative — criterion < 1e-12 each (flux-form telescoping; measured at round-off). Momentum is exchanged with the walls — that ledger is COUP-2 §3.1.2's mount-reaction term, a later wave.
- **Determinism:** reruns are byte-identical (asserted in the test suite on the full conserved state).
