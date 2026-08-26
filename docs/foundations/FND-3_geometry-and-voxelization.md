# FND-3 — Geometry & Voxelization

| Field | Value |
|---|---|
| **ID** | FND-3 |
| **Family** | FND (Foundations / spine) |
| **Status** | Reviewed (2026-08-14) |
| **Depends on** | FND-1, FND-2 |
| **Version** | 0.4 (2026-08-25 plan S9: the build wave — §3.1/§3.3/§3.5 as-built specifics for the sampled path (per-cell counter-based jitter keys, six-aperture emission incl. θ-faces, PLIC normal from the authored classifier's gradient, offset by deterministic bracketed solve on the monotone S-Z volume function), §3.2 exact winding evaluation with Barnes-Hut a recorded perf deferral, §3.4 conservative S9 `N_θ^geom` floor form). 0.3 (review fix wave S10 + v1.4 cylindrical metric) |

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
| **Refinement tags** (§3.4) | FND-2 | cells to refine at interfaces/gradients + per-region `N_θ^geom` azimuthal floors |

**Invariant:** emitted per-cell fractions are **partial** (accurate sub-cell coverage in [0,1]), not boolean
occupancy; and the pipeline is **deterministic** (exact predicates + fixed seeds).

## 3. Method

### 3.1 CSG path (primary) — analytic SDF tree
Geometry is an **analytic signed-distance tree**: primitive SDFs (sphere, box, cylinder, cone, torus),
booleans (**union = min(a,b), intersection = max(a,b), subtraction A−B = max(a, −b)** — always the correct *sign*
everywhere), and **revolved profiles** (2-D profile → axisymmetric solid, matching the axisymmetric engine
bias). Per voxel, the **exact partial volume** is obtained by analytic half-space/box clipping where a leaf
is a plane/box, else by **stratified point sampling** of the SDF — a **fixed jittered pattern** of
**`N_FRAC_SAMPLES` = 8³ = 512** strata per cut cell (named constant, config-overridable;
interior/exterior cells are decided by corner/SDF-bound tests, not sampled). The **declared per-cell
fraction-error bound is derived from the pattern** (S10, 2026-08-14): jittered-stratified sampling of
an indicator with a smooth interface converges at `N^(−1/2−1/(2d))` = `N^(−2/3)` in 3-D, so
`ε_α = C_jitter · N_FRAC_SAMPLES^(−2/3)`, with `C_jitter` calibrated once against the analytic
primitives of §6.1; the derived bound is **recorded in the run manifest** (FND-6) and feeds FND-2 §5's
geometric-error budget. [META-3: `jittered-sampling`] Caveat: `max`/
subtraction give a distance *bound* (sign correct, magnitude unreliable near concave joins) — fine for
occupancy, and the reason we sample rather than trust the magnitude. Because CSG is analytic and mesh-free,
it is exact, cheap, deterministic, and **parametric/sweepable**. [META-3: `gen-winding-number` (contrast),
`eb-cutcell`]

### 3.2 STL path — generalized winding number, not raw ray-parity
Raw **ray-parity fails silently and catastrophically** (whole regions flipped) on non-watertight or
degenerate meshes, so it is **not** the primary classifier. Instead use the **generalized winding number**
(Jacobson 2013): a smooth harmonic field that is integer for watertight meshes and **degrades gracefully**
for imperfect ones — robust inside/outside by **thresholding at 0.5**, no mesh repair. **Evaluation is the
exact triangle sum at S9** (0.4, as built): the generalized winding number is evaluated by the closed-form
per-triangle solid-angle sum (van Oosterom–Strackee), O(n_tri)/query — exact, deterministic, and affordable
at config time for the mini-sim/import meshes the sessions use; the **Barnes-Hut fast approximation**
(Barill 2018, O(log n)/query) is a **recorded perf deferral** (owner: the production-STL wave, with GPU
voxelization if ever needed) — it changes cost, never the classification contract. **Partial fractions**
come from averaging the winding number over **stratified per-voxel samples**. Ray-parity is kept only as a
corroborating check on known-watertight input. Triangle-touch gating uses the Akenine-Möller SAT test.
STEP→STL conversion is the supported CAD route (FreeCAD etc.). [META-3: `gen-winding-number`]

### 3.3 Voxelize() output — fractions, apertures, surfaces
*(v1.4 — cylindrical cell metric, D-A.)* Voxelization targets the grid's **cylindrical-structured ring
cells** (FND-2 §3.2): cell volume `½(r_o²−r_i²)·Δθ·Δz` (∝ r̄), six faces `{r−, r+, θ−, θ+, z−, z+}`.
Sample points are stratified **in the cylindrical volume measure** (uniform in `r²`, θ, z — the jitter
pattern is uniform over the cell's actual volume) and evaluated in Cartesian coordinates, where the SDF /
winding number is authored (the mapping is exact); PLIC planes are fitted in the cell's local Cartesian
frame. Voxelization always runs at the config-declared finest azimuthal resolution `N_θ^max`; FND-2's
adaptive coarsening starts from these fractions.

Emit, per cell: **material volume fractions** α_k (partial, Σ = 1 incl. vacuum); **six face apertures**
(open-area fraction per face); and, for cut cells, the **interface centroid, outward normal, and area** via
**PLIC** (volume-exact, closed-form Scardovelli–Zaleski offset — deterministic). PLIC areas (not Marching-
Cubes areas, which are ~8% biased) feed radiation view factors and flow-path wall area, keeping radiation
reciprocity consistent with the fractions. Where both CSG and STL describe the same part, **cross-check**
their fractions — disagreement beyond sampling tolerance flags a defective STL export or CSG eval (a free
consistency check). [META-3: `vof-plic`, `mc33-dc`]

**As built (0.4, plan S9) — the pieces the paragraph above left open:**
- **Aperture sampling measures.** Each face aperture is the stratified sample fraction open in that face's
  own area measure: r-faces uniform in (θ, z) at fixed r_f; z-faces uniform in (r², θ); θ-faces uniform in
  (r, z) on the face's flat (r, z) rectangle (the θ-face of a ring cell is planar; its area element is
  dr·dz). One fixed 2-D jittered pattern per face (`N_FRAC_SAMPLES^(2/3)` = 64 strata), same key scheme
  as §3.5.
- **The PLIC normal comes from the authored geometry, not the fractions:** for CSG the outward normal is
  the normalized central-difference gradient of the SDF at the cell centroid, evaluated in Cartesian
  where the field is authored, then the plane is fitted in the cell's local Cartesian frame as stated.
  For STL the winding number of a watertight mesh is **piecewise constant** (0/1) — its gradient is
  analytically zero away from the surface and cannot supply a normal — so mesh cut cells take the
  **area-weighted outward facet normal** of the triangles touching the cell (fixed mesh order): exact
  for flat facets, first-order in curvature — the same honesty class as PLIC itself. (A Youngs'-style
  fraction-gradient normal needs a materialized neighborhood and is *less* accurate where the authored
  geometry is available — which at config time it always is.)
- **The volume-exact offset is solved, not case-tabled:** the plane constant d is the root of the
  monotone Scardovelli–Zaleski volume function V(d) = α, found by a deterministic fixed-structure
  bracketed bisection (fixed iteration count, bracket = the cell's diagonal support interval). This is
  numerically equivalent to the closed-form case inversion the 0.3 text named (same monotone cubic,
  root to round-off) with none of its case-table surface; the closed form remains a legal optimization.
- **The exact contour clip stays the axisymmetric exact path.** The engine's contour-of-revolution
  ingest (analytic clipping, zero sampling error) is the degenerate exact form of this section and is
  untouched — the sampled path is for genuinely 3-D authored geometry, and the two paths cross-check on
  a revolved solid exactly like CSG↔STL above (the certified stations ride the exact path unchanged).

### 3.4 Refinement tagging (config-time)
Tag cells for refinement where the **material-interface / volume-fraction gradient** is high (primary) and
where initial-field gradients are steep (Löhner). Hand tags to FND-2 §3.6, which buffers, 2:1-balances, and
freezes. **Geometry-driven N_θ floor (v1.4):** from the emitted fractions, compute per (r,z) region the
azimuthal variation of geometry (variance over θ of α_k and face apertures); a region whose *geometry* is
non-axisymmetric receives an **`N_θ^geom` floor** — the coarsest N_θ that reproduces its fractions within
the declared ε_α (§3.1) — handed to FND-2 §3.4 as a per-region azimuthal-resolution floor that holds
**independent of the flow state**. **S9 conservative form (0.4):** a brick whose sampled fractions or
apertures vary over θ beyond the declared ε_α receives the floor `N_θ^geom = N_θ^max` (the voxelized
resolution itself) — geometry-bearing θ-variation is never coarsened at S9; the coarsest-reproducing-N_θ
search this paragraph specifies rides the S10 refinement-tile wave, where the projection machinery it
needs lands anyway. A θ-uniform brick (an axisymmetric part) receives no floor. The **small-cut-cell** stability issue (a tiny κ collapsing the explicit step) is handled downstream
by **State Redistribution** in the field solvers (FND-2/SOLV), designed in from the start — FND-3 just
produces the fractions. [META-3: `lohner`, `state-redistribution`]

### 3.5 Determinism
Exact geometric predicates + **fixed tie-breaking** + **fixed sample seeds** (the root cause of ray-parity/
SAT/winding-number flakiness is inconsistent FP predicates). Voxelization is a pure function of
{geometry inputs, resolution, seed} → reproducible. **Key scheme (0.4, as built):** each stratum's jitter
comes from a counter-based hash (SplitMix64) keyed on `(seed, (r,z)-cell / canonical-face index,
face/volume tag, stratum index)` — a pure function of the addressed sample, so the pattern is independent
of evaluation order and thread schedule, and any cell can be re-sampled in isolation bit-identically.
The key deliberately excludes the θ-sector index: every sector of an (r,z) cell reuses one pattern
(sector samples are exact rotations of each other), so an axisymmetric solid voxelizes bit-identically
per sector up to the classifier's own rotational round-off (the SDF is evaluated at rotated Cartesian
points; only a sample within an ulp of the surface could classify differently — an ulp-thin set the
§6.6 zero-azimuthal-variance gate pins in practice). Shared faces key canonically from either side
(the FND-2 §3.6 bitwise face-coherence contract).

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
   expected order under refinement, and at the declared `N^(−2/3)` jittered-sampling rate under
   sample-count growth — the fit constant is the `C_jitter` of §3.1's error bound.
2. **Imperfect-STL robustness:** winding-number classification stays correct on meshes with holes/self-
   intersections where ray-parity fails.
3. **PLIC area/volume consistency:** extracted interface area is consistent with fractions; radiation
   reciprocity holds.
4. **CSG↔STL cross-check:** identical part via both paths agrees within sampling tolerance.
5. **Determinism:** identical voxelization across runs; fixed seeds reproduce fractions bit-for-bit.
6. **Cylindrical metric & floors:** an axisymmetric revolved solid yields zero azimuthal geometry
   variance (no `N_θ^geom` floor imposed); an off-axis feature yields the correct floor; ring-cell
   volumes ∝ r̄ recover analytic volumes.

## 7. References
META-3 keys: `gen-winding-number`, `vof-plic`, `mc33-dc`, `eb-cutcell`, `lohner`, `state-redistribution`,
`jittered-sampling` (v1.4). Depends on FND-1, FND-2 (§3.2 cylindrical index space).

## 8. Change log
| Date | Version | Change |
|---|---|---|
| 2026-07-14 | 0.1 | Initial draft. CSG-first (analytic SDF, parametric/sweepable) with full STL support (generalized winding number, not raw ray-parity); partial volume fractions + face apertures + PLIC interface surfaces; config-time interface refinement tagging; deterministic (exact predicates + fixed seeds); CSG↔STL cross-check. Authoring-path emphasis marked as an open fork (FND-3-Q1). |
| 2026-07-19 | 0.2 | **CSG-first authoring confirmed by Ben** (2026-07-19); FND-3-Q1 resolved and its open-questions register removed (v1.3 convention: no in-doc question registers). |
| 2026-08-14 | 0.3 | **Review fix wave (S10; ruling D-A).** §3.1: fraction sampling pinned to a **fixed jittered pattern with named `N_FRAC_SAMPLES` = 8³ = 512** strata per cut cell; **declared per-cell fraction-error bound** `ε_α = C_jitter·N^(−2/3)` derived from it (C_jitter calibrated in §6.1) and recorded in the run manifest. §3.3 restated on the **cylindrical cell metric** (ring cells, volume ∝ r̄, sampling stratified in the cylindrical volume measure, evaluated in Cartesian; voxelization at finest `N_θ^max`). §3.4: **geometry-driven `N_θ^geom` floor** computed from azimuthal fraction variance and emitted to FND-2 §3.4 (holds independent of flow state). §2/§6/§7 aligned. |
| 2026-08-25 | 0.4 | **Plan S9 (the build wave — landed with the code).** §3.2: the winding number is evaluated by the **exact per-triangle solid-angle sum** at S9 (deterministic, config-time-affordable); Barnes-Hut becomes a **recorded perf deferral** (production-STL wave) — cost only, never the classification contract. §3.3 as-built block: face-aperture sampling measures per face class (r: (θ,z); z: (r²,θ); θ: (r,z) on the planar face rectangle); the **PLIC normal comes from the authored geometry** — CSG: the SDF's central-difference gradient; STL: the area-weighted outward facet normal (a watertight mesh's winding number is piecewise constant, so its gradient cannot supply one) — never a fraction-gradient reconstruction; the **S-Z offset is the root of the monotone V(d) by deterministic fixed-structure bracketed bisection** (equivalent to the closed-form case inversion, which stays a legal optimization); the exact contour clip is recorded as the axisymmetric exact path, untouched (certificates ride it). §3.4: **S9 conservative `N_θ^geom` floor** — θ-varying geometry pins the floor at the voxelized N_θ; the coarsest-reproducing search rides S10. §3.5: the **counter-based jitter key scheme** (SplitMix64 on (seed, (r,z)-cell/canonical-face, tag, stratum) — the θ-sector deliberately excluded, so axisymmetric solids voxelize bit-identically per sector and shared faces cohere bitwise) — order- and schedule-independent by construction. Measured at the gates: `C_jitter` 0.907 → shipped 1.4 (×1.5 margin); sampling-rate exponent −0.723 vs the derived −2/3. |
