# CRUCIBLE

**C**ross-**R**egime **U**nified **C**omparator **I**ntegrating **B**ounded-**L**ifetime **E**ngines —
a generalized physics sandbox for in-vacuum propulsion, from chemical to antimatter, evaluated under one
thermal / radiation / structural / uncertainty framework so regimes are cross-comparable and every output
is a distribution with a pedigree.

A 12-month inquiry project (Jul 2026 – Jun 2027). This repository holds the **design docs**
(Layer 1 + Layer 2) and the **Rust runtime workspace**: the FND-4 config loader + run manifest,
the FND-5 HDF5 table seam (statically pinned libhdf5), the FND-2 cylindrical world-state grid,
and the first two SOLV-1 operators — conduction (the Goal-A convergence certificate) and the
compressible-Euler flux operator (PPM + HLLC-Batten, the Goal-B Station-1 Sod certificate), each
certified by a committed, regenerable record in `certificates/`. Current goal: **Goal B, the
blind RL10** (see CLAUDE.md for station status).

## Build & test

Rust 1.93.1, pinned in `rust-toolchain.toml` (META-3 §2). `./scripts/check.sh` runs the VAL-3 §3.2
per-commit gate battery (fmt → clippy `-D warnings` → tests → certificate regenerate-and-diff) in
fixed order; CI mirrors it. Every session ends with the battery green and the work committed
(VISION_SCOPE §12).

## Reading order

1. **[`VISION_SCOPE.md`](VISION_SCOPE.md)** — Layer 1, the source of truth (what we build and why).
   Currently at v1.4 (§15 is the amendment log).
2. **[`docs/meta/META-0_master-catalog.md`](docs/meta/META-0_master-catalog.md)** — the catalog of every
   design doc and what each owns.
3. **[`docs/meta/META-1_design-philosophy.md`](docs/meta/META-1_design-philosophy.md)** — governing
   principles (Rules 12/13 supremacy; determinism; predictive-validation doctrine).
4. Then the specific module docs under `docs/` (each header notes what to read first).

## Architecture in one paragraph

The runtime is **one 3-D multiphysics field solver** on a single world-state grid that evolves all matter
and fields over a local medium-state vector `M` (Rule 12: one law per physical quantity, no `if(material)`
seams; Rule 13: engines are pure-data configurations, no per-concept code). Reactions are source terms
feeding shared transport. Dimensionality adapts (adaptive azimuthal resolution N_θ, FND-2 §3.4) to spend
3-D cost only where symmetry breaks. Runtime = 100% Rust; offline table generation = Python; the only seam is versioned
HDF5 tables. Every result is a p-box with a PCMM pedigree; the goal is to **resolve** untested-regime
uncertainty, not accept it.

## Layout

| Path | Contents |
|---|---|
| `VISION_SCOPE.md` | Layer-1 vision & scope (source of truth) |
| `docs/meta/` | Catalog, design philosophy, conventions, sources ledger |
| `docs/foundations/` | Spine: data/UQ model, grid+solver, geometry, config, tables, results, constitutive spine |
| `docs/coupling/` | Coupling & orchestration (time integration, audit, registry) |
| `crates/` | The Rust runtime workspace — one crate per doc area (META-2 §4); each cites the doc it implements |
| `scripts/check.sh` | The VAL-3 per-commit gate battery |

The 29 critical-path Layer-2 docs are **Reviewed (2026-08-14)** — the coding gate is open
(`REVIEW_FINDINGS.md` is the record). Still unreviewed: SOLV-5, OFFL-4 (deferred set — do not
implement); COUP-1 is a retired tombstone.
