# FND-3 — Geometry & Voxelization

| Field | Value |
|---|---|
| **ID** | FND-3 |
| **Family** | FND (Foundations / spine) |
| **Status** | Draft |
| **Depends on** | FND-1, FND-2 |
| **Version** | 0.2 |

---

## 0. Purpose

How geometry is **authored** and turned into the grid's per-cell material fractions + face apertures at
**config time**. After voxelization the result is just matter in the grid, identical to everything else
(Rule 13, FND-2 §1.2) — so this doc's whole job is the config-time act of getting accurate, deterministic,
sourced geometry *into* cells. It is where the robustness work lives.

## 1. Scope & razor ruling
**Owns:** the CSG kernel; STL import + voxelization; conversion of both into per-cell material volume
fractions, face apertures, and interface surfaces; config-time refinement tagging. **Defers:** how cells
store/evolve that geometry → FND-2; conservation of the on-grid update → COUP-2 (there is no reduced-mesh
conservative transfer post-v1.3 — the former COUP-1 is retired). An in-simulator CAD/B-rep kernel is
**permanently out of scope** (VISION_SCOPE §13); voxelization fidelity is controlled by grid refinement.

**Authoring-path decision** *(confirmed by Ben, 2026-07-19):*
**config-CSG is the first-class, polished path; STL import is fully supported.** Rationale: the backbone
sweep (VISION_SCOPE §3.0) requires *parametric, sweepable* geometry, which only CSG gives; STL covers
one-off CAD parts. Both produce identical downstream cell data, so this sets *investment emphasis*, not
capability.

## 2. Interfaces & contracts
| Interface | Consumed by | Contract |
|---|---|---|
| **CSG config** (§3.1) | config author, FND-4 | primitives + booleans + revolved profiles → an analytic SDF tree |
| **STL/mesh import** (§3.2) | config author | watertight-ish surface → robust inside/outside + partial fractions |
| **Voxelize()** (§3.3) | grid construction (FND-2 §3.6) | emit per-cell material fractions + 6 face apertures + interface (centroid/normal/area) |
| **Refinement tags** (§3.4) | FND-2 | cells to refine at interfaces/gradients |

**Invariant:** emitted per-cell fractions are **partial** (accurate sub-cell coverage in [0,1]), not boolean
occupancy; and the pipeline is **deterministic** (exact predicates + fixed seeds).

## 3. Method

### 3.1 CSG path (primary) — analytic SDF tree
Geometry is an **analytic signed-distance tree**: primitive SDFs (sphere, box, cylinder, cone, torus),
booleans (**union = min(a,b), intersection = max(a,b), subtraction A−B = max(a, −b)** — always the correct *sign*
everywhere), and **revolved profiles** (2-D profile → axisymmetric solid, matching the axisymmetric engine
bias). Per voxel, the **exact partial volume** is obtained by analytic half-space/box clipping where a leaf
is a plane/box, else by **stratified point sampling** of the SDF (fixed jittered pattern). Caveat: `max`/
subtraction give a distance *bound* (sign correct, magnitude unreliable near concave joins) — fine for
occupancy, and the reason we sample rather than trust the magnitude. Because CSG is analytic and mesh-free,
it is exact, cheap, deterministic, and **parametric/sweepable**. [META-3: `gen-winding-number` (contrast),
`eb-cutcell`]

### 3.2 STL path — generalized winding number, not raw ray-parity
Raw **ray-parity fails silently and catastrophically** (whole regions flipped) on non-watertight or
degenerate meshes, so it is **not** the primary classifier. Instead use the **generalized winding number**
(Jacobson 2013): a smooth harmonic field that is integer for watertight meshes and **degrades gracefully**
for imperfect ones — robust inside/outside by **thresholding at 0.5**, no mesh repair. Evaluate with the
**Barnes-Hut fast approximation** (Barill 2018), O(log n)/query. **Partial fractions** come from averaging
the winding number over **stratified per-voxel samples**. Ray-parity is kept only as a corroborating check on
known-watertight input. Triangle-touch gating uses the Akenine-Möller SAT test. STEP→STL conversion is the
supported CAD route (FreeCAD etc.). [META-3: `gen-winding-number`]

### 3.3 Voxelize() output — fractions, apertures, surfaces
Emit, per cell: **material volume fractions** α_k (partial, Σ = 1 incl. vacuum); **six face apertures**
(open-area fraction per face); and, for cut cells, the **interface centroid, outward normal, and area** via
**PLIC** (volume-exact, closed-form Scardovelli–Zaleski offset — deterministic). PLIC areas (not Marching-
Cubes areas, which are ~8% biased) feed radiation view factors and flow-path wall area, keeping radiation
reciprocity consistent with the fractions. Where both CSG and STL describe the same part, **cross-check**
their fractions — disagreement beyond sampling tolerance flags a defective STL export or CSG eval (a free
consistency check). [META-3: `vof-plic`, `mc33-dc`]

### 3.4 Refinement tagging (config-time)
Tag cells for refinement where the **material-interface / volume-fraction gradient** is high (primary) and
where initial-field gradients are steep (Löhner). Hand tags to FND-2 §3.6, which buffers, 2:1-balances, and
freezes. The **small-cut-cell** stability issue (a tiny κ collapsing the explicit step) is handled downstream
by **State Redistribution** in the field solvers (FND-2/SOLV), designed in from the start — FND-3 just
produces the fractions. [META-3: `lohner`, `state-redistribution`]

### 3.5 Determinism
Exact geometric predicates + **fixed tie-breaking** + **fixed sample seeds** (the root cause of ray-parity/
SAT/winding-number flakiness is inconsistent FP predicates). Voxelization is a pure function of
{geometry inputs, resolution, seed} → reproducible.

## 4. Coupling relationships
- **FND-2** receives fractions/apertures/surfaces and freezes topology; it never re-derives geometry.
- **FND-4** config carries the CSG tree / STL references + resolution + seeds.
- **SOLV-2** (radiation) and wall-area accounting consume the PLIC surfaces; **SOLV-8** recession updates
  fractions FND-3 produced (same representation, no re-voxelization at runtime).

## 5. Uncertainty & validity
Voxelization introduces **geometric discretization error** (fraction/interface-position error ~ cell size),
budgeted and declared (FND-1); a feature thinner than the finest cell is flagged under-resolved, not smeared.
Infrastructure rung; CI-gated.

## 6. Validation plan
1. **Analytic volumes:** recovered fractions of sphere/cone/revolved-profile converge to exact volume at the
   expected order under refinement.
2. **Imperfect-STL robustness:** winding-number classification stays correct on meshes with holes/self-
   intersections where ray-parity fails.
3. **PLIC area/volume consistency:** extracted interface area is consistent with fractions; radiation
   reciprocity holds.
4. **CSG↔STL cross-check:** identical part via both paths agrees within sampling tolerance.
5. **Determinism:** identical voxelization across runs; fixed seeds reproduce fractions bit-for-bit.

## 7. References
META-3 keys: `gen-winding-number`, `vof-plic`, `mc33-dc`, `eb-cutcell`, `lohner`, `state-redistribution`.
Depends on FND-1, FND-2.

## 8. Change log
| Date | Version | Change |
|---|---|---|
| 2026-07-14 | 0.1 | Initial draft. CSG-first (analytic SDF, parametric/sweepable) with full STL support (generalized winding number, not raw ray-parity); partial volume fractions + face apertures + PLIC interface surfaces; config-time interface refinement tagging; deterministic (exact predicates + fixed seeds); CSG↔STL cross-check. Authoring-path emphasis marked as an open fork (FND-3-Q1). |
| 2026-07-19 | 0.2 | **CSG-first authoring confirmed by Ben** (2026-07-19); FND-3-Q1 resolved and its open-questions register removed (v1.3 convention: no in-doc question registers). |
