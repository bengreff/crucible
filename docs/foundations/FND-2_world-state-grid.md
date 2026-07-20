# FND-2 — World-State Grid & Unified Field Solver

| Field | Value |
|---|---|
| **ID** | FND-2 |
| **Family** | FND (Foundations / spine) |
| **Status** | Draft |
| **Depends on** | FND-1, COUP-3 (time integration), SOLV-1 (field operator) |
| **Version** | 0.4 (v1.3 unified-grid pivot) |

---

## 0. Purpose

The world-state grid is **the single representation of all matter in the simulator** *and*, since the
**v1.3 pivot, the unified 3-D field solver itself**: one grid that **evolves** all matter and fields
(flow, plasma/MHD, conduction, radiation, transport, reactions-as-sources) over the medium-state
vector `M` every tick, **and audits its own conservation**. There is no separate "accountant" grid and
"physicist" reduced-dimension solvers — that split was a Rule-12 seam (two discretizations reconciled
by transfer machinery) and is retired (META-1 §5, v1.3). This doc fixes *what a cell is*, *how the grid
is structured*, *the unified conserved state it evolves*, and *how dimensionality adapts* — at spec
level (no Rust; Layer 3 owns the `struct`s).

The governing sentence, from Ben (2026-07-14, reaffirmed under v1.3): **this design is the backbone for
everything; all matter is treated the same, except for abstracted point-source engineering components
that are fully decoupled from the reaction.** The v1.3 corollary (Ben, 2026-07-19): **the grid does not
just store the backbone — it solves it, in full 3-D, resolving as much as the physics permits.**

## 1. Scope & razor ruling

### 1.1 What the grid is
- **The authoritative store of all matter.** Every piece of matter the simulation contains — imported
  solids, procedurally-generated solids, liquids, gases, plasmas, and vacuum itself — exists as cells in
  this grid, in **one uniform representation** (§3.3). There is no second matter store.
- **The mesh for *all* field physics** — the whole unified solver runs on these cells: the field
  operator (reacting flow + two-phase + two-fluid/MHD + conduction, SOLV-1), radiation transport
  (M1 thermal + multigroup Sₙ nuclear, SOLV-2), energetic-particle transport (SOLV-3), reaction
  sources (SOLV-4), structural margins (SOLV-6), and degradation/recession (SOLV-8) are all computed
  *on these cells* (§3.4).
- **The solver *and* the auditor** (VISION_SCOPE §7.1, v1.3): the grid evolves the conserved state `U`
  (§3.4) over `M` and, each step, audits its own global conservation (COUP-2). It no longer delegates
  flow/plasma/transport to separate reduced-dimension meshes — there is one discretization.

### 1.2 Rule 12 & Rule 13 compliance (explicit — this doc must obey them)
- **Rule 13 (sandbox).** Geometry may be *authored/imported* two ways (config CSG, STL import; FND-3),
  but **import is a config-time act only**. After voxelization, an imported solid is **indistinguishable
  from any other matter** — same cell model, same fields, same physics. There is no runtime "imported
  object" type, no per-geometry code path. Defining a new engine or reactor is pure data (FND-3/FND-4);
  the grid treats the result uniformly.
- **Rule 12 (no seams).** All matter is carried through the **medium-state vector `M`** (FND-1 §3.3) and
  the material sub-states (§3.3). No physical law reads a material label or an "is-imported" flag; laws
  read `M`. "Solid," "liquid," "plasma," "vacuum" are values of the cell's state, not types — a cell is
  reclassified continuously as its state evolves (a heated solid that melts and ionizes crosses no code
  boundary).

### 1.3 The one exclusion — abstracted point-source engineering (reaction razor)
The **only** things not represented as matter in the grid are engineering components **fully removed from
the reaction**: lasers/particle beams, pulsed-power drivers, pumps/pressurization, power supplies, the
radiator interface, reactivity-control actuators, antiproton delivery optics. These are **boundary
objects** (COUP-7) — point/boundary sources with a citation, validity envelope, and error band — because
their internal physics does not carry energy release or thrust (VISION_SCOPE §4.1, Reaction Razor).

**The line is the reaction, and it is drawn tightly (Ben, 2026-07-14):** *anything the reaction touches
is simulated in full — no magic interfaces.* A laser is abstracted (its photons enter as a boundary
source); but the moment that energy reaches matter, the matter, its heating, ablation, phase change,
and any reacting mixture are **full grid physics**. Ablating throats, NSWR fuel-in-water, and the CNTR
liquid-uranium film on hydrogen are **fully simulated multi-material cells** (§3.3), not interface
closures.

### 1.4 Not owned here (deferred, owner named)
- The value/units/uncertainty types and `M`'s field list → **FND-1**.
- Voxelizing CSG/STL into cell fractions → **FND-3** (FND-2 receives the fractions).
- Operator coupling on the one grid + the every-step conservation audit → **COUP-2**; time integration
  / operator split → **COUP-3**. *(The former reduced-dimension solver–grid binding, COUP-1, is retired
  — §3.4.)*
- The **constitutive spine** (EOS, transport, opacity, stopping over `M`) → **FND-7**.

## 2. Interfaces & contracts (what downstream binds to)

| Interface | Consumed by | Contract |
|---|---|---|
| **Cell model** (§3.3) | all operators (SOLV-1…8) | the uniform per-cell state: material sub-states, apertures, `M`, per-cell fields |
| **Conserved state `U`** (§3.4) | SOLV-1/2, COUP-2/3 | the evolved conserved variables per cell (mass, species, momentum, energy, face-`B`, radiation moments); `M` is reconstructed from `U` + the spine each step |
| **Cell sweep API** (§4) | SOLV-1/2/3/4/6/8, COUP-3 | ordered, deterministic traversal of active cells + face connectivity for stencils; operators read `M`/sub-states and write `U` **in place** |
| **Conservation ledger** (§4) | COUP-2 | per-quantity stored totals + port flux accounting for the every-step audit |
| **Geometry ingest** (§3.6) | FND-3 | accept per-cell material fractions + apertures at config time; freeze topology |

**Invariant promised to everyone:** there is exactly **one authoritative discretization** — the grid cells —
and every operator reads and writes them **in place**. There are no reduced-dimension meshes and no
gather/scatter "views" (retired with COUP-1); adaptive dimensionality is the spectral azimuthal-mode
reduction of §3.4. The grid's traversal order is deterministic and independent of thread count (§3.7).

## 3. Method & governing structure

### 3.1 One representation of all matter (the backbone)
A cell does not know whether its matter was imported, generated, or condensed from a reaction. It stores
matter as a set of **material sub-states** + a shared **medium-state vector `M`**, and that is all. This
is what makes the grid the backbone: adding a new geometry, swapping a reactor core, or having a solid
ablate into a plasma are all just changes to cell contents, never new representations or code paths
(Rules 12/13, §1.2).

### 3.2 Data structure
- **A static-topology, VDB-shaped sparse "brick tree"** (NanoVDB-style): shallow-wide, fixed-depth, with
  small dense leaf bricks (e.g. 8³) under one or two coarse index levels, linearized into a contiguous
  arena with computed offsets rather than pointers. Topology is built once at config time and **frozen**
  for the run; only cell *values* change (recession updates fractions, not topology — VISION_SCOPE §7.2).
  This is why the hand-rolled implementation stays small: the hardest VDB machinery (dynamic insertion,
  rebalancing, pruning) is unneeded. [META-3: `vdb`, `nanovdb`]
- **Coordinates are f64** (Ben, 2026-07-14). Cell positions are *derived* from the integer cell index and
  the grid's f64 origin/spacing (not stored per cell); all geometry the grid carries in f64 — interface
  positions, recession depths, EB centroids/normals, vertex coordinates. **Physical state fields are f64**
  as well (META-1 §3: f64 for all physics state/accumulation). f32 is used only for optional visualization
  export, never in the physics or the audit.
- **Variable-size cells via uniform grouping** (§3.5) — the "efficiency algorithm groups uniform cells
  into larger cells."
- **Deterministic layout:** active bricks are held in a **Morton (Z-order) array**; all sweeps and audits
  iterate that array, never a hash map (§3.7). [META-3: `p4est` for Morton order only]

### 3.3 The cell state model (uniform for all matter)
Every cell — vacuum, solid, fluid, plasma, or a mixture — carries the **same** structure. Single-material
and even vacuum cells are the degenerate case (one sub-state), so there is no seam between "simple" and
"mixed" cells (Rule 12). Everything is **f64**; a cell in one ensemble member holds **concrete values only**
— no ± / no distributions (UQ is the outer loop) and **no per-scalar provenance tags** (those live on
declarations/results, FND-1 §3.6). The organizing principle: **conserved variables (mass, energy, momentum,
composition) are authoritative and stored; everything a law wants — p, per-phase T, `M`, effective
coefficients — is *derived* from them**, so no two representations can drift.

**Identity (not stored as payload):** position and size are derived from the cell's Morton index + the grid's
f64 origin/spacing + level. No per-cell xyz.

**(1) Geometry** *(static; changes only via recession)* — material volume fractions α_k (Σ = 1, **including a
vacuum material** — "in vacuum" is just α_vacuum ≈ 1); **six face apertures** (open-area fraction per face,
gating the flux stencils; a receded face → aperture 0 drops out with no remeshing); and for cut cells the
interface **centroid, outward normal, area** (from PLIC — feeds radiation view factors + flow-path wall area).
[META-3: `eb-cutcell`]

**(2) Matter — segregated material sub-states** *(dynamic; one per phase, usually just 1)* — per material k:
- `material_id` → the material's **identity / fundamental data** in FND-7 (composition {Z_k, A_k}, ρ₀, cohesive energy, ionization potentials) + genuinely-static engineering limits (strength/allowables, DPA/burnup, emissivity class). **State-dependent constitutive quantities — EOS, conductivity, viscosity, opacity, stopping — are NOT keyed by `material_id`; they are computed from `M` by the constitutive spine (FND-7), never a per-material lookup** (Rule 12; cf. §3.4);
- authoritative **conserved** variables: partial mass (α_k ρ_k) and partial internal energy (α_k ρ_k e_k) —
  T_k and partial pressure are *derived* via EOS_k, so each phase holds its **own temperature**;
- composition sub-state where the physics needs it: reacting gas → elemental/mixture composition (equilibrium
  species from OFFL-3 tables); plasma → ⟨Z⟩, free/bound-electron split, T_e vs T_i; fissile → burnup; ablating
  solid → degree-of-char.

**Full segregated multi-material** (Ben's Q1 ruling): a 3000 K graphite grain and a 20 000 K plasma in one cell
are **not averaged**. Phases reconcile to a single cell `(p, u)` for the shared laws via a **rate-based
pressure-relaxation closure** (Tipton-style) that only *partitions* under the hard invariants Σα_k = 1,
mass = Σα_kρ_k, energy = Σα_ke_kρ_k — instantaneous relaxation → homogenized limit, finite rate → segregated;
one dial, no seam. [META-3: `multimat-closure`]

**(3) Reconciled cell state** *(derived from the sub-states via the closure; cached)* — common cell pressure p
and bulk velocity **u** (the shared momentum the flow solver reconciles to).

**(4) The medium-state vector `M`** *(a derived view — the only argument any constitutive law reads, Rule 12;
cached, rebuilt when sub-states change)* — number densities n_s & Z_s, ρ, free/bound-electron densities, ⟨Z⟩,
T_e/T_i, degeneracy θ, **B** (magnetostatics = config-time setup of SOLV-1, static during a run), flow **u**, plus a derived cache (plasma
frequency, Coulomb log, Debye length). One consistent `M` per cell that every regime sees identically. (Full
field list: FND-1 §3.3.)

**(5) Field-physics state** *(dynamic; written by on-grid solvers each step)* — volumetric nuclear heating rate,
n/γ dose rates, radiosity/exchange state for surface cells, structural margin flags (hoop/thermal/burst).

**(6) Degradation / lifetime clocks** *(dynamic but slow; stage-2)* — recession depth, accumulated neutron
fluence & DPA, burnup, corrosion depth, accumulated dose, decay-heat state.

**Not held per cell:** xyz (derived), any uncertainty/distribution, any provenance tag. *(There is no
solver-binding map — post-v1.3 there are no reduced-dimension regions; every operator runs on the cell in
place, §3.4.)* Fields are grouped **hot** (M, densities, T, p, heating — every step) vs **cold** (clocks,
fluence, flags — rarely) for
cache efficiency (§3.2).

*Effective scalars* (e.g. one conductivity for a diffusion step in a mixed cell) are computed *only where a
law needs a single number*, by the physically-correct mixing rule (series/parallel/Hashin–Shtrikman) — the
thermodynamic *state* is never homogenized, only the specific transport coefficient a law asks for.
[META-3: `multimat-closure`]

### 3.4 The unified conserved state and adaptive dimensionality *(v1.3)*
There is now **one discretization** — the grid cells — and the whole solver runs on it. The
reduced-dimension native meshes and the gather/scatter binding (former COUP-1) are **retired**.

**The conserved state `U` per cell** (evolved by SOLV-1/2, coupled by COUP-3, audited by COUP-2):
`U = ( ρ, {ρX_k}_species, ρ𝐮, ρE_total, 𝐁_face [constrained-transport, staggered], E_rad, 𝐅_rad )`.
`M` and everything derived (T_e, T_i, ⟨Z⟩, p, effective coefficients) is **reconstructed from `U` +
the constitutive spine (FND-7) each step**, so no two representations drift (§3.3). Reaction rates,
opacities, cross-sections, stopping, and EOS are **M-indexed table lookups, never branches** (Rule 12).
Radiation moments (E_rad, 𝐅_rad) live *in* `U` so the same finite-volume update transports them (M1
two-moment); ∇·𝐁 is held by **constrained transport** (𝐁 on faces) — a structural invariant, not a
runtime correction. [META-3: `castro-source`, `athena-ct`, `radiation-m1`]

**Time integration** is owned by COUP-3: an **SDC-coupled IMEX** scheme — explicit hyperbolic hydro/MHD
(PPM/PLM reconstruction + HLLD/HLLC flux) with the stiff pieces (M1 radiation source, conduction,
reaction sources) as **cell-local implicit** updates, iterated to 2nd order (Strang splitting is a
documented failure mode for stiff/energetic reactions — the nuclear/antimatter end — so it is not used
there). FND-2 provides the deterministic swept cell traversal + face connectivity the update runs on
(§3.7). [META-3: `sdc-imex`, `stiff-reactions`]

**Adaptive dimensionality — spectral azimuthal-mode truncation (the compute lever, §3.8).** 3-D is the
default, but the solver collapses the azimuthal direction where geometry *and* state are axisymmetric:
- Represent azimuthal fidelity **spectrally**: mode **m=0 everywhere** (the exact axisymmetric field),
  plus modes **m=1…M(r,z)** carried **only where** an azimuthal-energy indicator
  `A(r,z) = Σ_{m≥1}|û_m|² / (|û_0|²+ε)` exceeds a threshold τ. Collapse = truncate to m=0 (≈ "2-D-axi");
  re-expand = admit more modes where symmetry breaks (instabilities, plume asymmetry, film breakup,
  tilt). This is the FBPIC/QPAD quasi-3-D mechanism made *adaptive in M per region*. [META-3:
  `spectral-azimuthal`, `symmetry-indicator`]
- It is **one uniform law projected onto fewer azimuthal DOF** — no separate solver, no stitched
  dimensional interface. Zonal 1-D/2-D/3-D stitching is **rejected**: every precedent reflects waves /
  needs iterative coupling / adds an interface multiplier — a Rule-12 seam. [META-3: `dim-hetero-coupling`]
- **Conservation & honesty:** azimuthal truncation is an *exact conservative projection* for the m=0
  conserved integrals; the discarded Σ_{m≥1} energy is **bounded by τ**. So it is **conservative always,
  controlled-error by design, lossless never.** Between regions of differing mode count, missing modes
  carry zero flux (the mode-space analog of AMR refluxing).
- **Determinism & anti-thrash:** the indicator, the azimuthal transform, and the collapse/expand
  decision use fixed reduction order and are gated to coarse time intervals; **dual thresholds with
  hysteresis** (τ_collapse < τ_expand) + a dwell time prevent thrashing (§3.7).

*Consequence for scale/cost (VISION_SCOPE §8):* cost scales with **resolved azimuthal modes × (r,z)
cells**, not the nominal 3-D voxel count. A mostly-axisymmetric engine costs ≈ 2-D (m=0), with local
3-D only where symmetry breaks — which is why the ~5000 GPU-h backbone *sweep* is affordable at adaptive
dimensionality while full-3-D is reserved for anchors/validation and UQ at full-3-D is multi-fidelity
(META-1 §9; COUP-5).

### 3.4.1 Sub-cell reaction scales: why the grid needn't resolve the reaction
Reactions occur at scales (nuclear ~fm; atomic/chemical ~Å; reaction zones ~µm) **far below any cell**. The
grid never resolves them geometrically, and it does not need to — this is scale separation, done rigorously,
not a shortcut. The reason "full simulation of the reaction" (§1.3) is compatible with mm–cm cells:

- **A cell carries a reaction *rate density*, not a reaction.** A cell holds ~10¹⁹–10²³ nuclei/atoms; what it
  represents is a continuum field — reactions·m⁻³·s⁻¹, energy·m⁻³·s⁻¹, particle-birth·m⁻³·s⁻¹ (e.g. fission
  R = Σ_f·φ; fusion R = n₁n₂⟨σv⟩; chemistry = equilibrium heat release). Averaged over that many sites the
  law of large numbers makes the rate **deterministic to ~10⁻¹⁰**; the fm-scale event locations are
  irrelevant to the aggregate the simulation actually needs.
- **The small-scale physics IS simulated — offline, at its correct scale.** OpenMC (OFFL-1) resolves the
  actual nuclear transport / cross-sections / resonances over a representative unit cell; Geant4 the
  annihilation cascade; Cantera/CEA the chemical equilibrium. They emit **homogenized, self-shielded**
  effective cross-sections, reactivity coefficients, ⟨σv⟩, and deposition kernels, tabulated vs local state.
  At runtime a cell simply **interpolates** these given its `M`. This is the standard validated two-step
  method (fine lattice calc → homogenized constants → coarse core calc); the grid is the coarse step.
- **Sub-cell heterogeneity & self-shielding live in the table, not the grid.** Flux depression inside a fuel
  grain, resonance self-shielding, etc. are captured by the offline homogenization, so a coarse cell using
  the homogenized constant reproduces the fine-scale rate — benchmarked (KRUSTY ≤300 pcm, NERVA ≤10%).
- **Cell size is set by continuum-field gradients (mm–cm), not the reaction scale.** Refinement (§3.6) goes
  where T/flux/composition/interface gradients are steep — never "where the reaction is small" (it always is).
- **Energy deposition is handled by range, not site:** short-range products (fission fragments, α; ≪ cell)
  deposit locally in the source cell; long-range (n, γ; cm–m) are transported to other cells by the
  precomputed kernels (SOLV-6). The range test (R vs cell size) decides which, per product.

So this is the **opposite** of a magic interface: a homogenized cross-section from OpenMC + SANDY covariances,
validated against hardware and carrying propagated bands, is the reaction done *right* at the scale where its
fine physics belongs. "Full simulation" = physics fidelity + the sourcing/transport/feedback loop, **not**
geometric resolution of a nucleus (which is impossible and physically pointless — the aggregate is the
observable). Where homogenization strains (steep sub-cell flux gradients, sharp reacting interfaces, kinetic
non-equilibrium), the guards are refinement, table regeneration for the geometry class, the full multi-material
cell (§3.3), declared model-form bands (pedigree-flagged), and validity-envelope enforcement that flags/halts
rather than extrapolating (FND-1 §3.7).

### 3.5 Uniform grouping (variable-size cells)
Storing and sweeping many identical fine cells in a uniform region is wasteful. The grid represents a
uniform region as a single **tile / super-cell** carrying one shared state (a VDB *tile value*), so the
field physics and the matter ledger operate on **variable-size effective cells** — fine where the action
is, coarse where the state is uniform. [META-3: `vdb`]

Reconciling this with static topology (VISION_SCOPE §7.2): the **finest resolution is fixed at config
time** (§3.6); a tile is a *storage/compute representation* of a uniform sub-region of that fixed finest
grid, **not** a coarser topology. When a tile's contents become non-uniform beyond a deterministic
tolerance during a run (e.g. a heating front enters a previously-uniform region), the tile is
**materialized** — expanded to its constituent finest cells — a value-representation change within the
frozen finest topology, triggered by a deterministic uniformity test so it never depends on scheduling.
Tiles expand monotonically within a stage; re-coarsening happens only at stage boundaries (deterministic
re-gridding). This gives the efficiency Ben asked for without breaking determinism or the static-topology
contract.

### 3.6 Config-time construction & refinement
Built once, then frozen:
1. **Ingest geometry** from FND-3: per-cell material fractions + apertures (CSG via analytic sampling; STL
   via generalized winding number — FND-3 owns the algorithms; robustness is FND-3's problem, the grid just
   receives fractions).
2. **Refine** to the chosen finest resolution using indicators evaluated on the initial geometry/fields:
   the **Löhner normalized-second-difference** (steep gradients, needs no time history) plus a
   **material-interface / volume-fraction-gradient** tag (refine at interfaces). Anticipate where runtime
   gradients will develop (reaction zones, heat sources, interfaces) and refine there, since topology is
   frozen. [META-3: `lohner`, `berger-amr`]
3. **Buffer** tags, enforce **2:1 balance** and **proper nesting**, cluster into bricks, group uniform
   regions into tiles (§3.5), and **freeze** the topology.

### 3.7 Determinism (Tier-1 bit-exact, any thread count — META-1 §2)
- **Canonical iteration:** every sweep/audit walks the **Morton-ordered active-brick array**; never a hash
  map (Rust `HashMap` iteration is randomly seeded → non-deterministic). Hash maps are used only at
  config-time construction, then baked into the sorted index array. [META-3: `gamer2-determinism`]
- **Order-independent reductions:** conservation totals use a **fixed brick-chunk decomposition + fixed-shape
  tree combine** (optionally compensated), giving bit-identical results on any core count — determinism does
  not depend on thread count. No `rayon::sum`/atomics in the physics/audit paths. [META-3: `repro-sum`,
  `fp-nonassoc`]
- **Tile materialization** (§3.5) is triggered by a deterministic per-tile uniformity test, so the storage
  representation — and therefore every downstream reduction order — is a function of the data, not the
  schedule.

### 3.8 Memory & scale (honest, at f64)
Per-cell payload at f64 with the multi-material model is heavier than the f32 figures from the raw survey:
a single-material cell is on the order of ~0.2 KB (M + fields), and a mixed cell scales with its number of
sub-states. Structural overhead (bitmasks/headers) is <0.3 %. So **stored memory ≈ (number of *distinct*
stored cells) × payload**, and the operative word is *distinct*: uniform grouping (§3.5) means the stored
cell count is the number of genuinely non-uniform/refined cells, typically **far below** the nominal finest
voxel count.

- **10⁷–10⁸ distinct f64 cells fit comfortably** on a 64–128 GB box (single-GB to tens-of-GB), and are the
  expected operating range for a per-engine model.
- **10⁹ *distinct* f64 heavy cells would exceed 128 GB** — so a "billion-cell" run means a billion-cell
  *nominal finest resolution* reached largely through uniform grouping (few distinct heavy cells), or a
  deliberately capped high-cost study, **not** a billion distinct multi-material cells. The number that
  scares the memory budget is *distinct grid cells* — and with adaptive azimuthal-mode reduction (§3.4) a
  mostly-axisymmetric engine stores ≈ its (r,z) plane (m=0), not the nominal 3-D voxel count.
- Sweeps are **bandwidth-bound** (~seconds per billion cells touched); minimizing bytes-touched-per-cell
  (hot/cold field grouping) matters more than FLOPs. [META-3: `amrex` layout, `aosoa-cabana`]

## 4. Coupling relationships
- **COUP-3 (time integration):** wraps the whole per-step update — the SDC-coupled IMEX advance of `U`
  (§3.4), the adaptive azimuthal-mode collapse/expand decisions, and the pulsed-event sequencing. FND-2
  supplies the deterministic swept traversal + face connectivity; COUP-3 owns the operator-split schedule.
- **COUP-2 (conservation audit + operator coupling):** the grid provides the per-quantity stored ledger
  and port hooks; COUP-2 runs the every-step `Δ(stored) = Σ port fluxes + sources` audit (Modelica-style
  flow connectors), halting on violation beyond tolerance. Because the field update is flux-form, the audit
  telescopes to boundary/port accounting. [META-3: `fv-telescoping`, `modelica-connector`]
- **The operators sweep the cells directly** (SOLV-1/2/3/4/6/8), reading `M` and material sub-states,
  writing to `U` and the derived/field state. Recession (fractions per §3.3; physics owned by SOLV-8) is
  the one that mutates geometry. *(COUP-1 solver–grid binding is retired — there are no reduced-dimension
  native meshes to bind; the adaptive azimuthal-mode reduction of §3.4 replaces it.)*
- **`M` is read identically by every operator** (Rule 12): the grid guarantees one consistent `M` per cell,
  reconstructed from `U` + the constitutive spine (FND-7); no operator recomputes matter classification or
  branches on a material/regime label.

## 5. Uncertainty & validity
The grid introduces **discretization/geometric error**, not physical model-form error: voxelization error in
material fractions and interface position, and field-physics truncation error. Both are **budgeted and
declared** (FND-1 uncertainty typing; the interpolation/discretization error feeds the UQ bands, META-1 §3).
Validity is a resolution statement: a run declares its finest resolution and the resulting geometric fidelity;
a feature thinner than the finest cell is flagged as under-resolved rather than silently smeared. Validation-
ladder status: infrastructure (rung i — analytic), CI-gated (§6).

## 6. Validation plan
1. **Voxelization accuracy:** recovered volume fractions of analytic CSG primitives (sphere, cone, revolved
   profile) converge to exact volumes under refinement at the expected order.
2. **Conservation on the grid:** a closed field-physics sweep (conduction with no ports) conserves energy to
   the stated tolerance; port accounting closes `Δ(stored) = Σ fluxes` to ~1e-10 relative.
3. **Determinism:** identical results (bytes) for a representative model at 1 vs N threads and across repeat
   runs; tile materialization is reproducible.
4. **Recession vs analytic:** a receding-front test reproduces an analytic Stefan/ablation solution within
   tolerance (guards the fixed-grid fraction-update approach). [META-3: `ablation-recession`]
5. **Uniform-grouping fidelity:** a tiled uniform region gives bit-identical field-physics results to the
   fully-materialized region until a gradient triggers materialization; no accuracy loss from grouping.
6. **Multi-material invariants:** Σα_k = 1, mass, and energy are conserved through a pressure-relaxation step;
   two phases retain distinct temperatures.

## 7. Open questions
| # | Kind | Question |
|---|---|---|
| FND-2-Q1 | needs-analysis | Confirm fixed-grid recession-front accuracy is adequate vs body-fitted TPS practice by validating against an analytic ablation/Stefan solution (test §6.4). Closes at SOLV-8 time; does not block the grid design. |

## 8. References
META-3 keys: `vdb`, `nanovdb`, `amrex`, `p4est`, `aosoa-cabana`, `lohner`, `berger-amr`, `repro-sum`,
`gamer2-determinism`, `fp-nonassoc`, `eb-cutcell`, `state-redistribution`, `multimat-closure`,
`ablation-recession`, `fv-telescoping`, `modelica-connector`; and (v1.3 unified solver, §6.5) `castro-source`,
`athena-ct`, `radiation-m1`, `sdc-imex`, `stiff-reactions`, `spectral-azimuthal`, `symmetry-indicator`,
`dim-hetero-coupling`. Depends on FND-1 (`Quantity`, `M`, uncertainty types), COUP-3 (time integration),
SOLV-1 (field operator), FND-7 (constitutive spine).

## 9. Change log
| Date | Version | Change |
|---|---|---|
| 2026-07-14 | 0.1 | Initial draft. Backbone/one-matter-representation principle (Rules 12/13); full segregated multi-material cells at reacting interfaces (no magic interfaces — only reaction-decoupled components abstracted); f64 coordinates & state; uniform grouping into variable-size tiles; explicit grid-cells-vs-reduced-solver-meshes resolution distinction; static-topology sparse brick tree; deterministic Morton traversal. |
| 2026-07-14 | 0.2 | Added §3.4.1 (sub-cell reaction scales: scale separation & homogenization — reactions carried as continuum rate densities, fine physics done offline and homogenized, cell size set by continuum gradients, deposition by range). Fixed a stale §3.9 recession reference. |
| 2026-07-14 | 0.3 | Expanded §3.3 into a complete categorized cell-data enumeration (geometry, segregated matter sub-states incl. composition/ionization/burnup/char, reconciled state, `M`, field-physics state, degradation clocks; conserved-authoritative/derived split; explicit "not held per cell": xyz, uncertainty, provenance, binding map). |
| 2026-07-20 | 0.4* | **Consistency-review fixes** (same version, contract-surface completion): the §2 interface table was rewritten (dropped the retired COUP-1 region-binding + gather/scatter rows and reduced-solver "views" invariant → Conserved-state `U` + in-place Cell sweep API on the one grid; SOLV IDs → operators); §3.3(2) `material_id` restricted to identity/fundamental data + static handbook limits (EOS/transport/opacity/stopping come from the spine over `M`, not a per-material lookup — closes a Rule-12 seam the review flagged); §3.3 "not held per cell" solver-binding-map clause removed. |
| 2026-07-19 | 0.4 | **v1.3 unified-grid pivot.** Retitled *World-State Grid & Unified Field Solver*. The grid now **evolves** the physics (is the solver) and audits itself; the accountant/physicist split and reduced-dimension native meshes retired (§0, §1.1). **§3.4 rewritten**: one discretization; the unified conserved state vector `U` (with radiation moments in `U`, constrained-transport `B`); time integration = SDC-coupled IMEX (COUP-3); **adaptive spectral azimuthal-mode dimensional reduction** (m=0 everywhere + adaptive m≥1 where an azimuthal-energy indicator fires; conservative + controlled-error, *not* lossless; hysteresis; zonal stitching rejected). §4 coupling updated (COUP-1 retired, COUP-3 central, SOLV IDs → operators); §3.8 scale note updated; references updated. §3.3 cell model, §3.2 sparse structure, §3.4.1 sub-cell homogenization, §3.5–3.7 all retained (still valid). |
