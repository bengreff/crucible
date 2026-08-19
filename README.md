# CRUCIBLE

**C**ross-**R**egime **U**nified **C**omparator **I**ntegrating **B**ounded-**L**ifetime **E**ngines —
a generalized physics sandbox for in-vacuum propulsion, from chemical to antimatter, evaluated under one
thermal / radiation / structural / uncertainty framework so regimes are cross-comparable and every output
is a distribution with a pedigree.

A 12-month inquiry project (Jul 2026 – Jun 2027). This repository holds the **design docs**
(Layer 1 + Layer 2), the **Rust runtime workspace**, and the **Python offline pipelines**.
Certified so far (each by a committed, regenerable record in `certificates/`): the Goal-A
conduction convergence certificate and Goal-B stations 1–5 — Sod shock tube, choked De Laval
nozzle, the CEA→HDF5 flame seam, the conjugate cooled wall, and the RL10 assembly vs the
TM-107318 reference p-box (coarse tier). Current goal: **finish the chemical-regime sandbox**
per `PLAN_CHEMICAL_SANDBOX.md` — the plan of record (full-3-D spark-to-steady certification;
CLAUDE.md carries current state; SESSION_LOG.md the history).

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
| `PLAN_CHEMICAL_SANDBOX.md` | The plan of record: session-by-session build plan to the finished chemical-regime sandbox |
| `docs/meta/` | Catalog, design philosophy, conventions, sources ledger |
| `docs/foundations/` | Spine: data/UQ model, grid+solver, geometry, config, tables, results, constitutive spine |
| `docs/coupling/` | Coupling & orchestration (time integration, audit, registry) |
| `crates/` | The Rust runtime workspace — one crate per doc area (META-2 §4); each cites the doc it implements |
| `offline/` | The Python OFFL pipelines (VISION_SCOPE §6 two-language rule; `crucible_offl`, venv-pinned) — offline only, never at simulation time |
| `tables/` | Committed production HDF5 tables (FND-5 schema, digest-pinned; the one cross-language seam) |
| `certificates/` | Goal-A/Goal-B station certificates — regenerated and diffed by the gate battery |
| `data/anchors/` | Cached validation-anchor sources (archival rule, META-1 P7; sha256 in META-3) |
| `scripts/check.sh` | The VAL-3 per-commit gate battery (5 gates: fmt, clippy, cargo test, offline pytest, certificate diff) |
| `SESSION_LOG.md` | Per-session history: measured data, findings, review waves |

The 29 critical-path Layer-2 docs are **Reviewed (2026-08-14)** — the coding gate is open
(the review register was closed 68/68-discharged and removed 2026-08-19; findings live in the
docs' change logs + git history). Still unreviewed: SOLV-5, OFFL-4 (deferred set — do not
implement); COUP-1 is a retired tombstone.
