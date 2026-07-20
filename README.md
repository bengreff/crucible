# CRUCIBLE

**C**ross-**R**egime **U**nified **C**omparator **I**ntegrating **B**ounded-**L**ifetime **E**ngines —
a generalized physics sandbox for in-vacuum propulsion, from chemical to antimatter, evaluated under one
thermal / radiation / structural / uncertainty framework so regimes are cross-comparable and every output
is a distribution with a pedigree.

A 12-month inquiry project (Jul 2026 – Jun 2027). This repository currently holds the **design docs**
(Layer 1 + Layer 2) — no runtime code yet.

## Reading order

1. **[`VISION_SCOPE.md`](VISION_SCOPE.md)** — Layer 1, the source of truth (what we build and why).
   Currently at v1.3 (the unified-grid pivot).
2. **[`docs/meta/META-0_master-catalog.md`](docs/meta/META-0_master-catalog.md)** — the catalog of every
   design doc and what each owns.
3. **[`docs/meta/META-1_design-philosophy.md`](docs/meta/META-1_design-philosophy.md)** — governing
   principles (Rules 12/13 supremacy; determinism; predictive-validation doctrine).
4. Then the specific module docs under `docs/` (each header notes what to read first).

## Architecture in one paragraph

The runtime is **one 3-D multiphysics field solver** on a single world-state grid that evolves all matter
and fields over a local medium-state vector `M` (Rule 12: one law per physical quantity, no `if(material)`
seams; Rule 13: engines are pure-data configurations, no per-concept code). Reactions are source terms
feeding shared transport. Dimensionality adapts (spectral azimuthal-mode reduction) to spend 3-D cost only
where symmetry breaks. Runtime = 100% Rust; offline table generation = Python; the only seam is versioned
HDF5 tables. Every result is a p-box with a PCMM pedigree; the goal is to **resolve** untested-regime
uncertainty, not accept it.

## Layout

| Path | Contents |
|---|---|
| `VISION_SCOPE.md` | Layer-1 vision & scope (source of truth) |
| `docs/meta/` | Catalog, design philosophy, conventions, sources ledger |
| `docs/foundations/` | Spine: data/UQ model, grid+solver, geometry, config, tables, results, constitutive spine |
| `docs/coupling/` | Coupling & orchestration (time integration, audit, registry) |

All docs are **Draft** — none is Reviewed/Frozen, so nothing is cleared for implementation yet.
