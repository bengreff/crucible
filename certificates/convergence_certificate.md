# CRUCIBLE Convergence Certificate (Goal A)

Heat conduction on the unified cylindrical world-state grid, driven config → loader → registry → grid → operator → this data. Criteria are the named constants in `crates/solvers/src/certificate.rs`, asserted by `crates/solvers/tests/goal_a_certificate.rs` and enforced against this committed file by gate 4 of `scripts/check.sh` (regenerate-and-diff). Regenerate: `cargo run --bin convergence_certificate` (deterministic — no RNG, no wall-clock; the git commit records provenance).

## MMS 2-D axisymmetric (N_θ = 1 asserted)

| n | h | L2 error | observed order |
|---|---|---|---|
| 8 | 0.12500 | 1.171240e-2 | — |
| 16 | 0.06250 | 2.944537e-3 | 1.992 |
| 32 | 0.03125 | 7.370861e-4 | 1.998 |

**Criterion:** observed order in [1.8, 2.2] (formal = 2).

## MMS 3-D with m = 2 θ-mode

| n | h | L2 error | observed order |
|---|---|---|---|
| 8 | 0.12500 | 9.122688e-3 | — |
| 16 | 0.06250 | 2.286054e-3 | 1.997 |
| 32 | 0.03125 | 5.717438e-4 | 1.999 |

**Criterion:** observed order in [1.8, 2.2] (formal = 2).

## Analytic anchors

- **Steady annulus (log profile), 32 radial cells:** max error 4.279e-4 relative to ΔT = 200 K — criterion < 2e-3.
- **Transient cylinder vs 5-term Bessel series at t̃ = 0.1 (crosses the r = 0 axis):** max error 9.466e-3 K on T₀ = 100 K — criterion < 0.5 K.
- **Closed insulated sweep, 500 steps:** total-energy drift 7.887e-16 relative — criterion < 1e-12 (flux-form telescoping).

Determinism: reruns are byte-identical (asserted in the test suite: field bits and config content hash).
