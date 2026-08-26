# FND-2 — World-State Grid & Unified Field Solver

| Field | Value |
|---|---|
| **ID** | FND-2 |
| **Family** | FND (Foundations / spine) |
| **Status** | Reviewed (2026-08-14) |
| **Depends on** | FND-1, COUP-3 (time integration), SOLV-1 (field operator) |
| **Version** | 0.5.4 (2026-08-26 plan S10: §3.6.1 — static (r,z) refinement level interfaces (the meridional sibling of §3.4's ring rule: fine-owns-flux aggregate, annular-metric single-difference well-balance) + the AMR gate's MEASURED NO-GO on dynamic front-tracking; the conservative primitive built/gated, the cross-pencil flux-register integration → S11). 0.5.3 (2026-08-25 plan S9: §3.4(iv) superseded — cut geometry legal at **uniform** N_θ ≥ 1 (per-θ-plane six-aperture cell geometry, wall-closure θ-limb, FND-3 S9 geometry floor); mixed-N_θ cut worlds still refuse → S11). 0.5.2 (2026-08-25 plan S8: §3.4 ring-interface exchange as built — 2:1 adjacency, fine side owns the flux, aggregate applied coarse-side, axis parity pairing as built). 0.5.1 (2026-08-19 plan-of-record note; prior: 0.5 v1.4 review fix wave) |

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
closures. The sharp liquid-surface capability this promises is designed into the cell state now as a
**dormant reserved capability** (§3.3(7), v1.4) and built with the liquid-interface wave.

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
gather/scatter "views" (retired with COUP-1); adaptive dimensionality is the adaptive azimuthal
resolution N_θ of §3.4 (v1.4). The grid's traversal order is deterministic and independent of thread count (§3.7).

## 3. Method & governing structure

### 3.1 One representation of all matter (the backbone)
A cell does not know whether its matter was imported, generated, or condensed from a reaction. It stores
matter as a set of **material sub-states** + a shared **medium-state vector `M`**, and that is all. This
is what makes the grid the backbone: adding a new geometry, swapping a reactor core, or having a solid
ablate into a plasma are all just changes to cell contents, never new representations or code paths
(Rules 12/13, §1.2).

### 3.2 Data structure *(v1.4 — natively cylindrical-structured, D-A)*
- **A natively cylindrical-structured index space:** cells are indexed `(i_r, i_θ, i_z)` about **one
  config-declared symmetry axis** per run (FND-4). Metric factors are precomputed per radial ring: cell
  volume `½(r_o²−r_i²)·Δθ·Δz` (∝ r̄), face areas likewise — the conservative flux-form update uses the
  true cylindrical metric, never a Cartesian approximation of it. At `r = 0` the axis gets a
  **reflecting treatment**: axis faces have zero area (they drop out of the flux stencil geometrically),
  and cross-axis stencil needs are met by **θ↔θ+π parity pairing** of the innermost rings (the standard
  conservative polar-axis treatment). [META-3: `cyl-axis-fv`]
- **A static-(r,z)-topology, VDB-shaped sparse "brick tree"** (NanoVDB-style): shallow-wide, fixed-depth,
  with small dense leaf bricks in the **(i_r, i_z) plane** (e.g. 8×8) under one or two coarse index
  levels, linearized into a contiguous arena with computed offsets rather than pointers. Each brick
  stores its azimuthal ring as **N_θ(brick) contiguous θ-planes** (§3.4). The **(r,z) topology is built
  once at config time and frozen** (recession updates fractions, not topology — VISION_SCOPE §7.2);
  **θ-resolution is dynamic**: an N_θ change is a deterministic value-representation change within the
  config-declared finest `N_θ^max`, exactly like tile materialization (§3.5), never a topology change.
  This is why the hand-rolled implementation stays small: the hardest VDB machinery (dynamic insertion,
  rebalancing, pruning) is unneeded. [META-3: `vdb`, `nanovdb`]
- **Coordinates are f64** (Ben, 2026-07-14). Cell positions are *derived* from `(i_r, i_θ, i_z)` and
  the grid's f64 axis frame (origin, axis direction) + spacings (not stored per cell); all geometry the
  grid carries in f64 — interface
  positions, recession depths, EB centroids/normals, vertex coordinates. **Physical state fields are f64**
  as well (META-1 §3: f64 for all physics state/accumulation). f32 is used only for optional visualization
  export, never in the physics or the audit.
- **Variable-size cells via uniform grouping** (§3.5) — the "efficiency algorithm groups uniform cells
  into larger cells."
- **Deterministic layout:** active bricks are held in a **Morton (Z-order) array over (i_r, i_z)**,
  θ-planes in fixed ascending `i_θ` order within each brick; all sweeps and audits
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
vacuum material** — "in vacuum" is just α_vacuum ≈ 1); **six face apertures** (open-area fraction per face
`{r−, r+, θ−, θ+, z−, z+}` of the ring cell, v1.4,
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

**(7) Reserved sharp-interface fields** *(dormant capability — designed now, built with the
liquid-interface wave; v1.4, D-B)* — per interface-bearing cell and tracked material pair: the interface
**unit normal n̂** and **PLIC plane constant d** (the same volume-exact plane representation FND-3 §3.3
emits at config time, so voxelization output is the initial condition when the capability activates).
**Capability spec:** a general **geometric sharp liquid-surface capture** — VOF/PLIC-class
volume-fraction advection: the tracked α_k becomes an *advected* conserved fraction (geometric,
conservative fluxes against the reconstructed interface plane) instead of static-with-recession, with
n̂/d re-reconstructed from the α_k field each step. It is a *general* free-surface / film /
immiscible-interface capability of the one solver — never tied to an engine or regime — and it is
**dormant**: the fields are reserved in the state layout now (activating it changes no layout and opens
no seam), the α_k-advection contract is stated here, and the numerics land with the future
liquid-interface SOLV wave. [META-3: `vof-plic`]

**Not held per cell:** xyz (derived), any uncertainty/distribution, any provenance tag. *(There is no
solver-binding map — post-v1.3 there are no reduced-dimension regions; every operator runs on the cell in
place, §3.4.)* Fields are grouped **hot** (M, densities, T, p, heating — every step) vs **cold** (clocks,
fluence, flags — rarely) for
cache efficiency (§3.2).

*Effective scalars* (e.g. one conductivity for a diffusion step in a mixed cell) are computed *only where a
law needs a single number*, by the physically-correct mixing rule (series/parallel/Hashin–Shtrikman) — the
thermodynamic *state* is never homogenized, only the specific transport coefficient a law asks for.
[META-3: `multimat-closure`]

### 3.4 The unified conserved state and adaptive dimensionality *(v1.3; adaptive mechanism v1.4)*
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
(PPM/PLM reconstruction + HLLD/HLLC flux), **spatially-coupled implicit diffusion** (conduction/viscous —
neighbor-coupled, solved by COUP-3 §3.1's deterministic fixed-cycle class-`D` solver; a cell-local conduction
solve does not exist), and the genuinely **cell-local implicit** stiff sources (reactions, M1 radiation
*source* coupling), iterated to 2nd order (Strang splitting is a
documented failure mode for stiff/energetic reactions — the nuclear/antimatter end — so it is not used
there). FND-2 provides the deterministic swept cell traversal + face connectivity the update runs on
(§3.7). [META-3: `sdc-imex`, `stiff-reactions`]

**Adaptive dimensionality — adaptive azimuthal resolution N_θ(r,z) (the compute lever, §3.8).** *(v1.4,
D-A — supersedes the v1.3 spectral azimuthal-mode formulation, which had no shock-capturing story for
nonlinear fluxes on spectral θ-modes; REVIEW_FINDINGS S3/S5.)* 3-D is the default; where geometry *and*
state are axisymmetric the solver **coarsens the azimuthal direction** — ordinary conservative
finite-volume **ring cells**, refined/coarsened in θ per region:

- **Mechanism.** Each (r,z) brick carries its ring at a resolution `N_θ(brick) ∈ {N_θ^max, N_θ^max/2, …,
  N_θ^guard}` (factor-2 levels of the config-declared finest `N_θ^max`). The full operator — PPM/PLM
  reconstruction, HLLC/HLLD fluxes, all source terms — runs **unchanged** at every N_θ: fewer θ-cells,
  the *same* shock-capturing numerics (what the spectral formulation could not offer). **`N_θ = 1` *is*
  the axisymmetric-with-swirl calculation** (the ring cell carries ρu_θ — "2.5-D swirl"). Between
  neighboring regions of differing N_θ, face fluxes are conservatively aggregated/subdivided
  (AMR-refluxing-style), so the audit telescopes as before. Zonal 1-D/2-D/3-D stitching remains
  **rejected** — this is one uniform law on fewer θ-cells, not a stitched dimensional interface.
  [META-3: `adaptive-theta-coarsening`, `dim-hetero-coupling`]
- **The ring-interface exchange, as built (S8, 0.5.2) — the "conservatively aggregated/subdivided"
  clause made concrete.** (i) **Adjacency is 2:1**: face-adjacent bricks may differ by at most one
  ladder factor (the θ-ladder's alignment then nests every fine sector exactly inside its coarse
  parent `j_c = ⌊j_f/2⌋`); a steeper jump refuses at operator validation — the same 2:1 balance the
  (r,z) refinement already enforces (§3.6). (ii) **The fine side owns the interface flux**: at an
  N_θ-jump face the flux is computed once per fine sub-face — the fine pencil's own reconstruction
  against ghost states **prolonged piecewise-constant** from the coarse neighbor — and the coarse
  cell receives the **area-weighted aggregate `Σ A_f·F_f`** of its children's fluxes. One computed
  number on both sides ⇒ the interface telescopes exactly and the COUP-2 audit sees an interior
  face, not a port. (iii) Ghosts for the coarse side's own reconstruction are the
  **equal-volume pair mean of the fine cells' states** (primitive operands — ghosts feed
  reconstruction only; the conservation statement is carried entirely by the flux ownership of (ii)). Piecewise-constant prolongation is
  locally first-order at the jump — declared, exactly like the one-sided lagged stencils at
  boundaries (the composed order is owned by the order gates); the *conservation* statement is
  exact regardless. (iv) **Cut-geometry worlds are legal at uniform N_θ ≥ 1 as of S9 (0.5.3 — the
  FND-3 3-D aperture wave, built):** cell geometry is stored **per θ-plane** with the **six**-face
  aperture set of §3.3(1) (grid storage order is the FaceDir index order {r−, r+, z−, z+, θ−, θ+};
  as of S10 the FND-3 voxelizer **emits in this same FaceDir order** — the single face-order owner —
  so the ingest seam is a plain identity copy, the S9 {r,θ,z}-emit permutation wart retired,
  FND-3 0.5); at N_θ = 1 every (r,z) value is bit-identical to the prior (r,z)-shaped form (the
  storage additionally carries the two θ-aperture arrays, κ-filled on the revolved path and
  provably not load-bearing there), and full-box worlds keep the no-geometry
  arithmetic-identity defaults. The discrete wall-closure identity gains its θ-limb
  `W_θ = (a_θ₊ − a_θ₋)·A_θ` (θ-face areas are θ-independent, so uncut cells still cancel
  bitwise). **Mixed-N_θ cut worlds still refuse** (the reflux of a *cut* jump face is S11 content,
  with mixed-N_θ class-D); EVERY geometry-bearing brick pins `n_theta_geom_floor` at the built N_θ (the hard S9 form —
  stricter than the θ-varying-only minimum FND-3's kernel computes; coarsen/refine/assert refuse
  on such bricks), so adaptive θ-resolution on cut worlds is wholly deferred to S11.
  (v) The **r = 0 axis at N_θ > 1** is the §3.2 parity pairing as built: the innermost
  ring's cross-axis ghosts gather from the **θ+π partner** cells with `u_r` and `u_θ` negated (the
  basis flip); at N_θ = 1 the partner is the cell itself and the gather is arithmetically identical
  to the reflecting mirror, so the axisymmetric corner is bit-unchanged.
- **Coarsening/refinement are conservative projections of `U` (S6).** Coarsen = volume-weighted
  averaging of conserved variables over merged θ-cells (ring integrals of mass, species, momentum,
  energy preserved **exactly**); refine = conservative prolongation (limited piecewise-linear in θ —
  introduces no new extrema, preserves the ring integrals exactly). **Collapse thermalizes discarded
  azimuthal kinetic energy:** projecting `ρ𝐮` onto merged cells conserves momentum but not the kinetic
  energy of the sub-ring velocity deviations; since `ρE_total` is conserved, that ΔKE reappears as
  internal energy. This is stated, not hidden: each collapse event **logs its thermalized ΔKE into the
  conservation-ledger diagnostics** (COUP-2), and the accumulated ΔKE is a term of the run's **declared
  truncation-error bound** (§5).
- **Symmetry indicator (S7) — concrete.** Per brick `b` and conserved field `q ∈ U`:
  `A_q(b) = Σ_cells V·(U_q − ⟨U_q⟩_ring)² / ( Σ_cells V·⟨U_q⟩_ring² + V_b·(U_q^floor)² )` — a
  **normalized azimuthal-variance energy norm**, with `⟨·⟩_ring` the volume-weighted ring mean and
  `U_q^floor` a per-field **absolute floor** (config-documented default: 10⁻⁶ × the field's global
  reference magnitude taken from the run's initial/inflow state), so the indicator has both a relative
  and an absolute scale and cannot blow up where `⟨U_q⟩ → 0` (stagnation, near-vacuum). The brick
  indicator is `A(b) = max_q A_q(b)` — over **ALL conserved fields**, so asymmetry in any quantity
  (composition striations, field asymmetry) is seen, not just velocity. Granularity: **per brick**.
  [META-3: `symmetry-indicator`]
- **Thresholds, cadence, dwell (named defaults with rationale).** `τ_collapse = 10⁻⁶`,
  `τ_expand = 10⁻⁴` (config-documented defaults): τ_collapse sits orders below any declared physics
  band, so a collapse discards only content already negligible against the result's error budget;
  τ_expand > τ_collapse by 10² gives **hysteresis** against thrashing while still re-resolving
  asymmetry when its energy fraction is far below the physics band. Evaluated every `N_sym = 32` steps,
  with a dwell of `N_dwell = 4` consecutive over/under-threshold evaluations before any change (both
  config-documented defaults); the evaluation uses fixed-order reductions (§3.7), so decisions are a
  function of the data, never the schedule.
- **Guard resolution (S4) — asymmetry stays detectable.** Regions eligible for *adaptive* collapse
  never drop below **`N_θ^guard` = 4** (named constant). Rationale: 4 is the coarsest ring that carries
  **both phases of m = 1** (sin θ and cos θ) — the first symmetry-breaking mode — with no blind
  orientation: at N_θ = 2 a perturbation aligned with the cell boundaries is invisible, and at N_θ = 1
  everything aliases into the ring mean and the indicator is structurally zero (the S4 defect). The
  **re-expansion trigger reads guard-ring variance growth, which is nonzero and can actually fire**, at
  a cost of only 4× the axisymmetric ring. Full `N_θ = 1` is available **only as a recorded config
  assertion** (a per-region or global axisymmetry assertion, carried into the results bundle and
  **pedigree-visible** — COUP-6/FND-6), never an adaptive decision.
- **Geometry floor.** A region whose *geometry* is non-axisymmetric receives an **`N_θ^geom` floor from
  FND-3 §3.4** (computed from the voxelized fractions' azimuthal variation), which holds **independent
  of the flow state** — flow-adaptive collapse can never under-resolve the geometry itself.

*Consequence for scale/cost (VISION_SCOPE §8):* cost scales with **Σ_bricks N_θ(brick) × (r,z) cells**,
not the nominal 3-D voxel count. A mostly-axisymmetric engine costs ≈ its (r,z) plane × N_θ^guard (× 1
only under a recorded axisymmetry assertion), with full azimuthal resolution only where symmetry breaks
— which is why the ~5000 GPU-h backbone *sweep* is affordable at adaptive resolution while full-3-D is
reserved for anchors/validation and UQ at full-3-D is multi-fidelity (META-1 §9; COUP-5).

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
  deposit locally in the source cell; long-range (n, γ; cm–m) enter the one transport operator (SOLV-2) —
  precomputed kernels are the **precomputed solution mode of that same operator** (valid within its
  tabulated geometry class), runtime Sₙ the general mode, with **exactly one mode per particle-class+band
  per run, config-declared** (E-1). The range test (R vs cell size) decides local-vs-transported, per product.

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
1. **Ingest geometry** from FND-3: per-cell material fractions + apertures at the finest `N_θ^max`, plus
   per-region `N_θ^geom` azimuthal floors (§3.4) (CSG via analytic sampling; STL
   via generalized winding number — FND-3 owns the algorithms; robustness is FND-3's problem, the grid just
   receives fractions).
2. **Refine** in the (r,z) plane to the chosen finest resolution using indicators evaluated on the initial geometry/fields:
   the **Löhner normalized-second-difference** (steep gradients, needs no time history) plus a
   **material-interface / volume-fraction-gradient** tag (refine at interfaces). Anticipate where runtime
   gradients will develop (reaction zones, heat sources, interfaces) and refine there, since topology is
   frozen. [META-3: `lohner`, `berger-amr`]
3. **Buffer** tags, enforce **2:1 balance** and **proper nesting**, cluster into bricks, group uniform
   regions into tiles (§3.5), and **freeze** the topology.

### 3.6.1 Static (r,z) refinement level interfaces *(S10 — the meridional sibling of §3.4's ring rule)*
A declared refinement zone (walls, injector face, throat — ruling #7) is a region held at a finer (r,z)
cell size than its surroundings; **tiles are the value representation, never a topology change** (§3.5) —
the finest resolution is fixed at config time and a coarse region is a super-cell of it. Where a coarse
region (cell size 2h) meets a fine region (size h) they share a **level interface**, and the conservative
exchange across it is the **exact meridional analogue of the §3.4 ring-interface rule**:

- **2:1 adjacency + proper nesting.** Face-adjacent (r,z) levels differ by at most one factor of two (the
  §3.6 balance already enforced), and a level owning an interface spans ≥ `NGHOST` uniform cells (the same
  nesting the θ-ladder carries) — a steeper jump refuses at operator validation.
- **The fine side owns the interface flux; the coarse cell applies the area-weighted aggregate.** A coarse
  face abutting two fine sub-faces takes the aggregate `Σ A_f·F_f` of its children's fluxes — **one computed
  number on both sides**, so the interface **telescopes exactly** (to the bit) and is an **interior face to
  the COUP-2 audit, never a port** (COUP-2 §3.1). The coarse side's own reconstruction reads ghosts prolonged
  piecewise-constant from across the jump (locally first-order at the interface — declared, order owned by the
  order gates; conservation exact regardless), exactly as §3.4(ii)–(iii).
- **The metric is the new content.** Unlike the θ-ladder, where every sub-sector has equal arc, an (r,z)
  **z-interface**'s fine children carry *unequal* **annular** z-face areas `½(r²_{k+1} − r²_k)·Δθ`; the
  aggregation weights are those annular ratios (an **r-interface**'s children share a radius, so they are
  equal-area — the θ-like case, exact). The **well-balanced uniform fixed point turns on the single-difference
  form**: the coarse cell reconstructs its interface-face area as the **children-sum**, not an independently
  metricked `½(r²_hi − r²_lo)·Δθ` — split accumulation / an independent coarse area breaks the fixed point
  (the trap the S8 gates caught, restated for the meridional metric). The annular reconstruction is in fact
  **bitwise-exact on a 2:1 off-axis interface** (the r²-band differences are Sterbenz-exact — measured 0 on
  the gate); the **S12 round-off class** (≤ ~1e-14 relative) is the declared *safety bound* for the general /
  near-axis case where Sterbenz can fail, reported like the S9 θ-varying well-balance, never a widened
  conservation tolerance.
- **Scope, as built (S10):** the conservative level-interface **primitive** (fine-owns-flux telescoping, the
  annular aggregation, the well-balanced fixed point) is built and gated on the real HLLC flux + the real
  cylindrical metric at N_θ = 1, class-A hydro, uncut. The build that threads it through the parallel-pencil
  `sweep_r`/`sweep_z` and the brick-arena refinement topology — a **cross-pencil flux register** (Berger–
  Colella), genuinely more than the within-pencil ring reflux — is recorded to **plan S11**, landing with its
  consumer (a refined 3-D RL10). This staging is not a shortcut but the **measured** disposition of the AMR
  gate below: mixed level × N_θ, level × cut geometry, and level × class-D all refuse there, typed.

**The AMR gate (ruling #7 — a MEASURED go/no-go on *dynamic* front-tracking refinement, not an assumption).**
The closure-set pushed front (SOLV-4 §3.6) travels at `S_T` with a width of a fixed `Θ` cells at every
resolution, so refinement sharpens *where* but not *when*. The S10 smeared-vs-sharp study marched the same
flame at coarse `h` and fine `h/2`: the front **timeline** (position vs time) agreed to **< 0.5 coarse
cells** over the march (a bounded sub-cell registration offset that *shrinks* as the march proceeds,
0.50 → 0.25 cells — the closure-set front carries the same speed on both grids; the separate `flame_1d`
gate pins the consumption speed grid-independent to ~2%), while refinement bought only a **1.93× sharper**
front at **4.0× cost per level**. **Verdict:
NO-GO** — dynamic front-tracking cannot move the timeline-driven COUP-4 verdict, so the plan's declared
fallback (**static refinement + closure-set speed = blurry front, correct timeline**; §7) stands *measured*,
and dynamic front-tracking (which would require amending the frozen-topology contract, §3.2/§0.5.1) is not
built. Static refinement's purpose is genuine geometry/wall/throat gradient resolution, never front-chasing.

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
  scares the memory budget is *distinct grid cells* — and with adaptive azimuthal resolution (§3.4) a
  mostly-axisymmetric engine stores ≈ its (r,z) plane × N_θ^guard (× 1 only under a recorded axisymmetry
  assertion), not the nominal 3-D voxel count.
- Sweeps are **bandwidth-bound** (~seconds per billion cells touched); minimizing bytes-touched-per-cell
  (hot/cold field grouping) matters more than FLOPs. [META-3: `amrex` layout, `aosoa-cabana`]

### 3.9 Execution model & data layout *(v1.4, D-H)*
- **GPU-portable SoA layout from day one.** Within a brick, per-field **structure-of-arrays** storage:
  each `U`/`M` component is a contiguous array over the brick's cells (hot/cold grouping per §3.3), and
  hot loops are index arithmetic over those arrays — **no pointer-chasing in hot loops**, no per-cell
  heap objects or dynamic dispatch. The same layout serves CPU SIMD and a GPU port unchanged; layout is
  a day-one constraint precisely so the port is a port, not a rewrite. [META-3: `amrex`, `aosoa-cabana`]
- **Sequencing: CPU fixed-order reference solver first.** The first implementation is the **CPU
  fixed-order solver — the Tier-1 correctness oracle** (byte-identical at any thread count, META-1
  §2.1). The **GPU port comes later as a *declared relaxed-reduction path***, validated against the CPU
  reference (VAL-3: negligible-tolerance + non-spiraling gate) and permitted in **non-chaotic regimes
  only** (S6; META-1 §2.1); chaotic/turbulent regimes stay on fixed-order paths on every target.
- **Which paths may later be declared relaxed-reduction:** the bulk per-cell
  reconstruction/flux/update sweeps and ensemble-member throughput — the bandwidth-bound majority
  (§3.8). **Never relaxed:** the conservation audit and ledger (COUP-2), halt logic, the
  symmetry-indicator evaluation and N_θ decisions (§3.4), and every RNG-consuming path — these remain
  fixed-order on every target, so verdicts and adaptivity decisions can never diverge between CPU and
  GPU.

## 4. Coupling relationships
- **COUP-3 (time integration):** wraps the whole per-step update — the SDC-coupled IMEX advance of `U`
  (§3.4), the adaptive-N_θ coarsen/refine decisions (§3.4), and the pulsed-event sequencing. FND-2
  supplies the deterministic swept traversal + face connectivity; COUP-3 owns the operator-split schedule.
- **COUP-2 (conservation audit + operator coupling):** the grid provides the per-quantity stored ledger
  and port hooks; COUP-2 runs the every-step `Δ(stored) = Σ port fluxes + sources` audit (Modelica-style
  flow connectors), halting on violation beyond tolerance. Because the field update is flux-form, the audit
  telescopes to boundary/port accounting; N_θ-collapse events log their thermalized ΔKE into the ledger
  diagnostics (§3.4). [META-3: `fv-telescoping`, `modelica-connector`]
- **The operators sweep the cells directly** (SOLV-1/2/3/4/6/8), reading `M` and material sub-states,
  writing to `U` and the derived/field state. Recession (fractions per §3.3; physics owned by SOLV-8) is
  the one that mutates geometry. *(COUP-1 solver–grid binding is retired — there are no reduced-dimension
  native meshes to bind; the adaptive azimuthal resolution of §3.4 replaces it.)*
- **`M` is read identically by every operator** (Rule 12): the grid guarantees one consistent `M` per cell,
  reconstructed from `U` + the constitutive spine (FND-7); no operator recomputes matter classification or
  branches on a material/regime label.

## 5. Uncertainty & validity
The grid introduces **discretization/geometric error**, not physical model-form error: voxelization error in
material fractions and interface position, and field-physics truncation error. Both are **budgeted and
declared** (FND-1 uncertainty typing; the interpolation/discretization error feeds the UQ bands, META-1 §3).
Validity is a resolution statement: a run declares its finest resolution and the resulting geometric fidelity;
a feature thinner than the finest cell is flagged as under-resolved rather than silently smeared. The
adaptive-N_θ machinery (§3.4) contributes a **declared truncation term**: azimuthal content discarded at
collapse, bounded by τ_collapse, plus the logged thermalized-ΔKE ledger — both reported with the run's
discretization budget. Validation-
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
7. **Conservative θ-coarsen/refine:** a coarsen→refine round trip preserves the ring integrals of all
   conserved variables exactly; face fluxes across an N_θ jump close the audit.
8. **Collapse thermalization accounting:** a seeded azimuthal velocity perturbation collapsed to guard
   resolution logs a thermalized ΔKE equal to the resolved-KE difference; total energy is unchanged.
9. **Guard-band re-expansion fires:** a growing seeded m=1 perturbation in an `N_θ^guard = 4` region
   drives the §3.4 indicator across τ_expand and triggers re-expansion (the S4 regression); the r=0
   axis treatment preserves a uniform free stream across the axis.

## 7. Open questions
| # | Kind | Question |
|---|---|---|
| FND-2-Q1 | needs-analysis | Confirm fixed-grid recession-front accuracy is adequate vs body-fitted TPS practice by validating against an analytic ablation/Stefan solution (test §6.4). Closes at SOLV-8 time; does not block the grid design. |

## 8. References
META-3 keys: `vdb`, `nanovdb`, `amrex`, `p4est`, `aosoa-cabana`, `lohner`, `berger-amr`, `repro-sum`,
`gamer2-determinism`, `fp-nonassoc`, `eb-cutcell`, `state-redistribution`, `multimat-closure`,
`ablation-recession`, `fv-telescoping`, `modelica-connector`; (v1.3 unified solver) `castro-source`,
`athena-ct`, `radiation-m1`, `sdc-imex`, `stiff-reactions`, `symmetry-indicator`,
`dim-hetero-coupling`; and (v1.4 cylindrical grid / adaptive N_θ / dormant interface capture)
`cyl-axis-fv`, `adaptive-theta-coarsening`, `vof-plic` — `spectral-azimuthal` retired with the
superseded formulation. Depends on FND-1 (`Quantity`, `M`, uncertainty types), COUP-3 (time integration),
SOLV-1 (field operator), FND-7 (constitutive spine).

## 9. Change log
| Date | Version | Change |
|---|---|---|
| 2026-08-26 | 0.5.4 | **Plan S10 (landed with the code): §3.6.1 — static (r,z) refinement level interfaces + the AMR gate.** The §0.5.1-staged static-refinement design is concretized as the **meridional sibling of §3.4's ring rule**: 2:1 (r,z) adjacency + proper nesting; the fine side owns the interface flux, the coarse cell applies the area-weighted aggregate (one number, telescopes to the bit, interior to the COUP-2 audit); tiles are value representation, never topology (§3.5). The **new content is the metric** — a z-interface's fine children carry unequal **annular** z-face areas, so the well-balanced uniform fixed point turns on the **single-difference form** (the coarse face area is the children-sum, not independently metricked — the S8 trap restated; the annular reconstruction residual is the S12 round-off class ≤ ~1e-14). **Built + gated (S10):** the conservative level-interface primitive on the real HLLC flux + real cylindrical metric (fine-owns-flux bitwise telescoping; the naive coarse-owns-flux leak mutation-proven load-bearing; the annular aggregation; the r-interface equal-area fixed point bitwise), N_θ = 1 / class-A / uncut. The **cross-pencil flux-register integration** through `sweep_r`/`sweep_z` + the arena refinement topology → **plan S11** (with its consumer, a refined 3-D RL10). **The AMR gate is MEASURED (ruling #7): NO-GO on dynamic front-tracking** — the smeared-vs-sharp study measured the front timeline grid-independent (< 0.5 coarse cells, non-growing) while refinement buys only 1.93× sharpness at 4.0× cost/level, so the §7 static-refinement-plus-closure-set-front fallback stands measured and dynamic front-tracking (a frozen-topology amendment) is not built. |
| 2026-08-25 | 0.5.3 | **Plan S9 (landed with the code): §3.4(iv) superseded — the FND-3 3-D aperture wave, uniform-N_θ tier.** Cell geometry is stored per θ-plane with **six** face apertures (the §3.3(1) face set, finally carried in full); `build_with_geometry`'s N_θ > 1 refusal is retired for **uniform** N_θ (mixed-N_θ cut worlds still refuse — a cut jump face's reflux is S11 content with mixed-N_θ class-D). At N_θ = 1 the storage layout and every value are bit-identical to the prior (r,z)-shaped form (gate 5's proof); full-box worlds keep the no-geometry arithmetic-identity defaults. The wall-closure identity gains `W_θ = (a_θ₊ − a_θ₋)·A_θ` (θ-face areas θ-independent ⇒ uncut cells cancel bitwise). EVERY geometry-bearing brick pins its floor at the built N_θ (the hard S9 form — stricter than the θ-varying-only minimum the FND-3 kernel computes; coarsen/refine/assert refuse), so adaptive θ-resolution on cut worlds is wholly deferred to S11. |
| 2026-08-25 | 0.5.2 | **Plan S8 (landed with the code): §3.4's ring-interface exchange made concrete.** The "conservatively aggregated/subdivided (AMR-refluxing-style)" clause is specified as built: 2:1 ladder adjacency (steeper refuses); **the fine side owns the interface flux** (computed once per fine sub-face from the fine reconstruction against piecewise-constant-prolonged coarse ghosts; the coarse cell applies the area-weighted aggregate — one number both sides, so the interface telescopes exactly and is an interior face to the COUP-2 audit, never a port); coarse-side ghosts = equal-volume pair mean of the fine states (primitive operands; conservation carried by the flux ownership); prolongation locally first-order at the jump (declared, like one-sided boundary stencils — order owned by the order gates, conservation exact regardless); cut-geometry worlds uniform-N_θ until FND-3's 3-D aperture wave. The §3.2 axis parity pairing recorded as built: cross-axis ghosts = θ+π partner with `u_r`/`u_θ` negated, arithmetically identical to the reflecting mirror at N_θ = 1. |
| 2026-08-19 | 0.5.1 | **Plan-of-record note (VISION_SCOPE v1.5, `PLAN_CHEMICAL_SANDBOX.md`).** Refinement is staged: **static declared (r,z,θ) refinement zones** (the §3.5 tile machinery — walls, injector face, throat; build wave S10) land first, inside the frozen-finest-topology contract (§3.2); **dynamic front-tracking refinement** is a measured go/no-go at S10 and, if taken, requires amending the static-topology ruling here first. GPU execution note: §3.9's "GPU as a declared relaxed-reduction path" is superseded by META-1 §2.5 — the GPU build is **bit-exact per device** (Ben 2026-08-19). No other contract change. |
| 2026-07-14 | 0.1 | Initial draft. Backbone/one-matter-representation principle (Rules 12/13); full segregated multi-material cells at reacting interfaces (no magic interfaces — only reaction-decoupled components abstracted); f64 coordinates & state; uniform grouping into variable-size tiles; explicit grid-cells-vs-reduced-solver-meshes resolution distinction; static-topology sparse brick tree; deterministic Morton traversal. |
| 2026-07-14 | 0.2 | Added §3.4.1 (sub-cell reaction scales: scale separation & homogenization — reactions carried as continuum rate densities, fine physics done offline and homogenized, cell size set by continuum gradients, deposition by range). Fixed a stale §3.9 recession reference. |
| 2026-07-14 | 0.3 | Expanded §3.3 into a complete categorized cell-data enumeration (geometry, segregated matter sub-states incl. composition/ionization/burnup/char, reconciled state, `M`, field-physics state, degradation clocks; conserved-authoritative/derived split; explicit "not held per cell": xyz, uncertainty, provenance, binding map). |
| 2026-07-20 | 0.4* | **Consistency-review fixes** (same version, contract-surface completion): the §2 interface table was rewritten (dropped the retired COUP-1 region-binding + gather/scatter rows and reduced-solver "views" invariant → Conserved-state `U` + in-place Cell sweep API on the one grid; SOLV IDs → operators); §3.3(2) `material_id` restricted to identity/fundamental data + static handbook limits (EOS/transport/opacity/stopping come from the spine over `M`, not a per-material lookup — closes a Rule-12 seam the review flagged); §3.3 "not held per cell" solver-binding-map clause removed. |
| 2026-07-19 | 0.4 | **v1.3 unified-grid pivot.** Retitled *World-State Grid & Unified Field Solver*. The grid now **evolves** the physics (is the solver) and audits itself; the accountant/physicist split and reduced-dimension native meshes retired (§0, §1.1). **§3.4 rewritten**: one discretization; the unified conserved state vector `U` (with radiation moments in `U`, constrained-transport `B`); time integration = SDC-coupled IMEX (COUP-3); **adaptive spectral azimuthal-mode dimensional reduction** (m=0 everywhere + adaptive m≥1 where an azimuthal-energy indicator fires; conservative + controlled-error, *not* lossless; hysteresis; zonal stitching rejected). §4 coupling updated (COUP-1 retired, COUP-3 central, SOLV IDs → operators); §3.8 scale note updated; references updated. §3.3 cell model, §3.2 sparse structure, §3.4.1 sub-cell homogenization, §3.5–3.7 all retained (still valid). |
| 2026-08-14 | 0.5 | **Review fix wave (S3–S9, E-1; rulings D-A/D-B/D-H).** §3.2 restated on the **natively cylindrical-structured index space** (`(i_r, i_θ, i_z)` about a config-declared axis; volumes ∝ r; reflecting r=0 axis treatment; static (r,z) topology with dynamic θ-resolution inside a config-declared `N_θ^max`) (S3). §3.4 adaptive dimensionality rewritten as **adaptive azimuthal resolution N_θ(r,z)** — conservative ring-FV coarsening/refinement running the same shock-capturing operator at every N_θ, superseding the v1.3 spectral-mode formulation (S3/S5); collapse = **conserved-variable projection with thermalized ΔKE logged per event** into the ledger diagnostics and the declared truncation bound (S6); **concrete symmetry indicator** (normalized azimuthal-variance energy norm over all conserved fields, per-brick, absolute+relative floor, named τ_collapse/τ_expand + cadence/dwell defaults with rationale) (S7); **`N_θ^guard` = 4** guard resolution so re-expansion can actually fire, full N_θ=1 only a recorded pedigree-visible config assertion (S4); **geometry-driven `N_θ^geom` floor** from FND-3. §3.3(7): **reserved sharp-interface fields** — dormant VOF/PLIC-class liquid-surface capture (normal + plane constant reserved, α_k-advection contract stated; general capability, built with the liquid-interface wave) (S9/D-B). New §3.9 **execution model**: GPU-portable SoA layout mandated day one; CPU fixed-order reference solver first (Tier-1 oracle), GPU later as a declared relaxed-reduction path on throughput sweeps only (S8/D-H). §3.4.1 deposition wording aligned to **E-1** (kernels = precomputed solution mode of the one transport operator; one mode per particle-class+band per run, config-declared). §2/§3.6/§3.8/§4/§5/§6/§8 aligned. |
