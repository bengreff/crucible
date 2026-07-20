# COUP-1 — Solver–Grid Binding & Conservative Source-Term Mapping

| Field | Value |
|---|---|
| **ID** | COUP-1 |
| **Family** | COUP (Coupling & orchestration) |
| **Status** | **Superseded (v1.3, 2026-07-19)** |
| **Depends on** | FND-1, FND-2 |
| **Version** | 0.1 |

---

> **⚠ SUPERSEDED by the v1.3 unified-grid pivot.** There is now **one discretization** — the world-state
> grid *is* the solver (FND-2 v0.4) — so there are no reduced-dimension native meshes to bind, and the
> gather/scatter solver–grid mapping this doc specified is **retired**. Its role (spending 3-D cost only
> where it's needed) is filled instead by **adaptive spectral azimuthal-mode reduction**, owned by
> **FND-2 §3.4**. The conservative-transfer ideas below (donor-cell, the `2πr` factor, kernel line
> sources, extensive/intensive operators) are kept **for history and for any future offline
> mesh-to-grid mapping**, but nothing in the runtime binds through this contract. The ID is retired and
> never reused (META-0 §6). Do not build against this doc.

## 0. Purpose

The **spatial-transfer primitive** on which the whole accountant↔physicist architecture rests: how a
reduced-dimension solver's native mesh (1-D flow path, 2-D axisymmetric r-z zone, magnetic flux tube,
1-D pellet line) **binds** to a region of the 3-D world-state grid (FND-2), and how a single quantity
moves between them **conservatively and deterministically** across non-conforming meshes. COUP-2 builds
the seven named couplers on top of this; COUP-1 owns the mapping math and the binding contract.

## 1. Scope & razor ruling

**Owns:** the region-binding descriptor; the precomputed geometric **weight maps**; the **gather**
(grid→solver) and **scatter** (solver→grid) operators; the conservation guarantee for a single quantity.
**Defers:** the catalog of *what* is exchanged and the global audit → COUP-2; the reduced solvers' native
meshes themselves → each SOLV doc; conservative *time* integration/operator splitting → COUP-3.

Razor: this is infrastructure, but it is the load-bearing wall — if the mapping is not conservative, no
downstream result can be. It obeys Rule 12 (one uniform transfer law over all bindings; no per-solver
special case) and serves Rule 13 (a new mechanism binds through this contract with no bespoke code).

## 2. Interfaces & contracts

| Interface | Consumed by | Contract |
|---|---|---|
| **Region binding** | every reduced solver, COUP-2 | a solver declares its native-mesh elements + the grid region each spans; COUP-1 returns a frozen weight map |
| **Weight map** | COUP-2 couplers | sparse `{reduced_element → [(cell, weight)]}` with the *same* weights used both directions |
| **Gather(field)** | solvers | reduced-element input = weighted integral/average of grid-cell state |
| **Scatter(source)** | solvers | deposits a solver source term onto the bound cells with the *same* weights |
| **Field kind tag** | all callers | each transferred field is `extensive` or `intensive` (§3.3), selecting the correct operator |

**The one invariant promised to everyone (§3.1):** for every extensive quantity, *what a solver emits
equals what the grid receives, to machine precision*, independent of thread count and iteration order.

## 3. Method

### 3.1 The single load-bearing idea: scatter is the transpose of gather
Bind on **one set of fixed geometric weights** `w(e, c)` (reduced element `e`, grid cell `c`). Then:
- **Gather:** `state(e) = Σ_c w(e,c) · state(c)` (weighted integral/average).
- **Scatter:** `source(c) += w(e,c) · source(e)`.
Because scatter uses the *same* weights as gather (the discrete adjoint/transpose), the deposited integral
equals the emitted integral exactly, and the result is order-independent → deterministic. This one
relationship is the mathematical core shared by MOOSE's conservative transfers, Farrell's supermesh
projection, and the Peaceman well index. [META-3: `moose-multiapp`, `supermesh`, `peaceman-well`,
`consistent-vs-conservative`]

### 3.2 The weight map — donor-cell / exact-intersection volumes
Weights are **exact geometric intersection volumes** between a reduced element's swept region and each grid
cell it overlaps (donor-cell). This is the practical minimum that is **both conservative and bounded** (no
over/undershoot) — critical when depositing heat/mass that must stay physical. Higher-order Galerkin/
supermesh projection is available where first-order smearing is unacceptable, but it is *not* bounded and
needs a positivity pass, so it is opt-in, not default. Weights are computed **once at config time** (static
topology, FND-2) and frozen. [META-3: `donor-cell`, `supermesh`]

### 3.3 Extensive vs intensive — the transpose distinction (do not mix up)
- **Extensive** (heat, mass, momentum, charge, particle count): use the **conservative** operator (weights
  sum to 1 over *sources* → global sum preserved). This is scatter.
- **Intensive** (temperature, pressure, velocity): use the **consistent** operator (weights sum to 1 over
  *targets* → a constant field is reproduced exactly). This is gather of a driving BC.
The two operators are transposes; on non-conforming meshes you cannot have both properties at once, so every
coupled field is explicitly tagged and the correct operator is selected. Using the wrong one silently breaks
either conservation or constant-reproduction. [META-3: `consistent-vs-conservative`]

### 3.4 Dimensional-reduction binding (the hard part, with a clean rule)
Every reduced solver is an integral over a swept geometry of the 3-D field; the gather/scatter pair on shared
weights handles all of them uniformly. Specifics:
- **1-D flow path:** a station represents a swept cross-section; weights = the volume of each cell within the
  station's slab. Where a station meets 3-D, use **interface coupling** (cross-section integral as the
  conservative interface). [META-3: `peaceman-well`]
- **2-D axisymmetric r-z zone:** an r-z cell is a revolved annulus. **Every weight carries the `2πr`
  revolved-volume factor** — omitting it silently loses mass. (Non-negotiable; validated by the
  `axi-2pir` result.) Azimuthal swirl momentum is carried as an extra component, transferred like any field.
- **Magnetic flux tube / pellet line (1-D in 3-D):** use a **kernel-distributed source** over a tube around
  the segment, *not* a singular Dirac line — the Dirac line makes the 3-D field singular on the centerline
  and wrecks convergence. [META-3: `koch-line-source`]

### 3.5 Determinism
Weights are frozen f64 values (FND-1/FND-2). Gather/scatter sums use the grid's canonical Morton order and
fixed-order reductions (FND-2 §3.7) → bit-identical on any thread count. No weight depends on schedule,
wall-clock, or hash order. The binding map is built at config time using position→cell lookup (a hash map
used only there), then baked into sorted sparse weight arrays.

## 4. Coupling relationships
- **FND-2** provides cell geometry (fractions, apertures, volumes) and canonical traversal; COUP-1 builds
  weights from them and never mutates topology.
- **COUP-2** consumes weight maps to implement the 7 couplers and runs the global audit (which relies on
  §3.1 holding per-coupler).
- **Every SOLV reduced solver** binds through §2; the native mesh remains the solver's, the authoritative
  matter stays on the grid (FND-2 §3.4).
- **COUP-3** wraps gather→solve→scatter in the time loop / Picard iteration.

## 5. Uncertainty & validity
COUP-1 adds **discretization/transfer error**, not model-form error: donor-cell transfer is first-order in
the mesh-size ratio (a smearing of gradients across the coarser mesh), budgeted like interpolation error and
fed to UQ (FND-1). Binding is valid where reduced-element and cell resolutions are compatible; a binding
whose weight map leaves cells or elements unclaimed is a config error and halts. Validation: infrastructure
(rung i), CI-gated (§6).

## 6. Validation plan
1. **Exact conservation:** `Σ_c scatter = solver-emitted source` to machine precision, for 1-D, 2-D-axi,
   flux-tube, and pellet-line bindings (including the `2πr` case explicitly).
2. **Constant reproduction:** the intensive operator reproduces a uniform field exactly.
3. **Determinism:** identical transfers at 1 vs N threads and across runs.
4. **Round-trip:** gather then scatter of a conserved quantity returns it unchanged (adjoint check).
5. **Order-of-accuracy:** donor-cell transfer converges first-order under refinement; supermesh higher.

## 7. References
META-3 keys: `moose-multiapp`, `supermesh`, `donor-cell`, `consistent-vs-conservative`, `peaceman-well`,
`koch-line-source`, `axi-2pir`. Depends on FND-1 (`Quantity`), FND-2 (cell geometry, traversal).

## 8. Change log
| Date | Version | Change |
|---|---|---|
| 2026-07-14 | 0.1 | Initial draft. Adjoint gather/scatter on shared frozen weights (scatter = transpose of gather) as the conservation guarantee; donor-cell exact-intersection weights (conservative + bounded); extensive/intensive operator distinction; dimensional-reduction binding with the mandatory `2πr` factor and kernel-distributed line sources; deterministic Morton-order transfer. |
