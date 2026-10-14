# Adaptive resolution (TECHNICAL_PLAN *Lightweight engine*, order of work item 5)

Status (13 October 2026, about 22:45 CDT): design and criteria, stated before any AMR code. Nothing here is verified yet.

The plan chose block-structured AMR in the logical (i, j) index space of the body-fitted mapping: refinement ratio 2 in space and time, time-accurate subcycling, conservative averaging with refluxing (Berger and Colella, J. Comput. Phys. 82, 1989), and geometry taken from the finest level. It stated four criteria. This record turns that choice into a data layout, a step, tags and numbers, and adds three criteria of its own (the geometry refactor's identity, cost and determinism, and the RL10 measurement the decision rule needs).

The physical picture, for orientation. The RL10A-3-3A gas path is about 1.5 m long. Most of it is slow, smooth gas in the chamber and nearly one-dimensional supersonic gas in the bell. The cells need to be small in only a few places: around the throat, where the area changes fastest and the gas turns sonic; in a thin strip along the wall, where the wall function needs the first cell inside the log layer; at the flame and the mixing layers near the injector; and at the shocks of the start. A uniform mesh fine enough for those places spends most of its work on gas that a coarse mesh already describes. AMR keeps a coarse base everywhere and lays finer copies of the same mapping over the places that need them. Each finer copy takes two steps of half the size per step of the one beneath, so no level is held to the smallest cell's step.

Labels: *derived* is my own reasoning, *guessed* a number not yet measured, *cited* from a source I name.

## Scope and declared limits

- **Index space, not physical space.** Refinement splits a cell into 2 x 2 children in (i, j). The mapping (the contour r(z) and the ring fractions) is evaluated once, on the finest level allowed. A coarse cell is the exact union of its children, so the wall and the throat area are as accurate on every level as on the finest (TECHNICAL_PLAN, *Adaptive resolution*).
- **At most three levels above the base** (levels 0 to 3), declared. The RL10 base of 150 x 16 at level 3 is 1,200 x 128, which holds the finest geometry in memory easily (*derived*: about 150,000 cells of node and face data).
- **Refinement in whole parent cells.** A level's active region is any union of parent cells, held as sorted j-intervals per column. This is more general than the plan's "rectangles of index space" (see *Deviations*). The kernels still run unchanged: each level's arrays span the whole domain in that level's index space, so every stencil reads its neighbours where it always has.
- **Explicit steps only.** Every level uses the engine's SSPRK2 at its own stable step. Implicit time is not part of this design; the decision rule decides it later from the measurement in criterion 6.
- **Chemistry modes.** Table A local equilibrium (`eqtt`) and frozen chemistry act cell by cell inside the flow step, so they carry over without change. The finite-rate path (ReactingFlow's split reaction call) is per cell with the level's step; it is not exercised by these criteria and is declared untested under AMR.
- **Liquid spray is out of scope.** Parcels (core/spray.*, the other worker's lane) would need a location on the finest covering level; that is left to the spray owner.

## Geometry: one mapping for every level (implementation step A)

Today `Mesh` computes face areas, reference radii and moments from `radius[]` (the wall at each axial face), `fraction[]` (the ring boundaries) and `dz`, through accessors such as `axialArea(face, j)`, `radialAreaVector(i, face)`, `radialFaceRadius(i, face)`, `ringMiddle(face, j)` and the `Cell` records. A coarse cell built from fine cells has a polyline for its sloped faces, which those formulas cannot describe.

- **The refactor.** `Mesh` stores every geometric quantity the kernels read as explicit arrays: per cell the `Cell` record; per axial face the area and the ring middle; per radial face the area vector, the reference radius r_f, the r^2 moment and the middle radius; the wall as a polyline of fine segments (for the wall distance, `walls.cpp`). The current constructor fills them with the current formulas, so a single-level run does the same arithmetic in the same order. Every call site reads the arrays.
- **Coarsening.** A second constructor builds the level-l mesh from the level-(l+1) mesh:
  - volumes, `radialPressureMeasure` (the integral of dV / r) and the r-weighted integrals behind `r` and `radialSecondMoment` are sums of the children's;
  - axial faces lie on lines of constant z, so a coarse axial face is two fine ones stacked: area summed, ring middle at the mean of its end radii;
  - a coarse radial face is the polyline of two fine radial faces: area vector summed; r_f and the r^2 moment are the means weighted by the radial area component, which is exact for anything linear in r (or r^2), the property r_f is defined by; the middle radius is the shared fine node;
  - the wall polyline is the finest one on every level.
- **Nesting holds.** The ring fractions are tanh(beta s) / tanh(beta) with s = j / nr (`Mesh::Mesh`), so the level-(l+1) fraction at 2j equals the level-l fraction at j. The axial stations nest the same way.
- **Closure.** For any closed polygon the sum of r n dl over its edges is (0, its meridional area), so the coarse cells' face sums close exactly onto their summed `radialPressureMeasure`. Quiescent gas stays quiescent on every level to round-off, as it does on one (*derived*; criterion 1 checks it).
- **The cost of the alternative.** If a coarse cell had straight chords instead, a 10 mm chord on the RL10 throat arc would put the throat area off by about 0.7% (TECHNICAL_PLAN, *derived*). With the PCHIP contour's throat radius of curvature, 0.137 r_t = 8.6 mm, against the 10 mm base spacing, the error would be larger still.

## Levels, active regions and ghost cells

- **A level is a full-domain array.** Level l holds every field of `Flow` on (nz 2^l) x (nr 2^l) cells, with three per-cell states: *active* (updated by this level), *ghost* (filled from the level below before each stage) and *off*. Memory is trivial at RL10 size. The physical boundaries (axis, wall, injector face, outlet) stay where the kernels expect them, and no patch-to-patch neighbour search is needed.
- **Active region.** Per column, a sorted list of j-intervals whose ends are even (whole parent cells). Level 0 is active everywhere.
- **Ghost width 2**, from the widest stencil:
  - MUSCL: the slope of the ghost cell next to the interface needs the ghost beyond it;
  - the viscous face flux needs the least-squares gradient of the ghost next to the interface, which reads its own neighbours.
  - Criterion 3(c) checks that width 3 changes nothing.
- **Proper nesting.** The active region of level l+1, widened by the ghost width, lies inside the active region of level l. Every ghost of level l+1 then sits over active, uncovered level-l cells. Tags are widened on the coarser level until this holds (Berger and Colella, p. 68).
- **Kernels restricted to the active region.** The cell and face loops run over interval lists instead of whole ranges. `axialSlopes(i)` and `radialProfiles(i)` take a j-range, so a fine wall strip of two rows does not pay for whole columns. Pool blocks split the lists exactly as they now split index ranges.

## The step: subcycling and refluxing (implementation step B)

Recursive, as in Berger and Colella (pp. 67 to 70):

```
advance(l, t, dt):
    save the level-l state at t (the "old" state for the level above's ghosts)
    SSPRK2 step of level l over its active cells, its ghosts at t and t + dt
        interpolated from level l-1 (space: conservative linear; time: linear between
        level l-1's old and new states); record the effective flux at every face
        between an active cell and a ghost (as a fine face) and at every face between
        an uncovered cell and a cell covered by level l+1 (as a coarse face)
    if level l+1 exists:
        advance(l+1, t, dt/2); advance(l+1, t + dt/2, dt/2)
        average level l+1 down onto the level-l cells it covers
        reflux: add (sum of fine fluxes over both substeps - coarse flux) to each
        uncovered level-l cell beside level l+1, divided by its volume
```

- **The effective flux of SSPRK2.** The step is u^(1) = u^n + dt L(u^n), u^(n+1) = (u^n + u^(1) + dt L(u^(1))) / 2, so the face flux that moved the cells is (F(u^n) + F(u^(1))) / 2 times dt. The registers accumulate that for every face array that enters `gather` and the transport gather:
  - the convective bulk flux and the transported mass fractions, k and omega;
  - the viscous bulk flux and the species, k and omega diffusion;
  - the property-correction fluxes.
  - The cell sources (the pressure source, the correction's non-divergence part, the igniter, body force, heating) are not face fluxes and are not refluxed.
- **The base step** is dt_0 = min over levels of 2^l times that level's stable step (its active cells only), ended on schedule events as now. Events then fall on base-step boundaries.
- **Rejection.** If a stage on any level is not admissible, the whole base step is redone from the saved states of every level at half the base step, which keeps the ratio 2 between levels.
- **Averaging down** is the volume-weighted mean of the children for rho, rho u, rho v, E, the partial densities, rho k and rho omega. It conserves each exactly (up to rounding), so the composite's content is level 0's content after each average.
- **Split pieces.** The Strang-split SST sources, the wall functions' prescription of the first cells' k and omega, and the clip of negative rho k act cell by cell, with the level's own step. They are applied to ghosts as well, so a ghost wall cell carries the k and omega its own distance gives.

## Boundaries and sources on fine levels

- **Supplies.** A fine injector ring takes the supply of its base parent ring. The supply area summed over fine rings equals the base one (the faces nest), so the scheduled mass flow is delivered exactly on every level.
- **Igniter.** The power density is fixed on the base: P over the base igniter cells' volume. Every fine cell over a base igniter cell takes that density. The energy delivered by the composite is then exactly P times the time, whatever is refined.
- **Body force and heating** are inherited from the base parent in the same way.
- **Walls and wall functions.** Each level builds its wall faces from its own geometry and solves the law at its own first cells; the wall distance comes from the finest wall polyline on every level. Wall forces and heat flow are tallied on the finest level covering each wall face (*Measurements*).
- **Axis and outlet.** Unchanged on every level: face boundary conditions, not ghost cells.

## Tags and regridding (implementation step C)

Tags are evaluated on level l to decide level l+1, every k_r base steps (default k_r = 2, *guessed*), and widened by a buffer of b cells of level l (default b = 2, *guessed*) so a feature does not leave its refinement between regrids. The defaults are tuned in step D, and criterion 4 is judged with the values then fixed.

- **Static tags.**
  - Throat: the columns within a declared axial half-width of the throat station, to a declared level (default 0.5 r_t to level 2, *guessed*).
  - Wall strip: refine the first rows of a wall column until the first cell's y+ (from the wall function's u_tau) is at or below the target y+. The target is an input with default 300: the largest first-cell y+ at which the corrected SST's own layer passes criterion 2 a priori (`WALL_FUNCTIONS.md`, restated criterion 2; the law passes criterion 1(b) to 956). The wall-function owner may move it; the value used is printed with every run. Before the flow exists the strip takes the throat's level.
- **Dynamic tags**, each with a threshold *guessed* now and fixed in step D:
  - flame and mixing layers: |grad f_H| times the cell size above 0.02;
  - temperature fronts (the flame in `eqtt`, where no heat-release rate exists): |grad T| times the cell size above 0.05 T;
  - shocks: a pressure jump across a face above 10% of the smaller pressure (Berger and Colella tag shocks apart from their error estimate, pp. 75 to 76);
  - shear layers: |vorticity| times the cell size above 0.05 of the local speed of sound.
- **New fine cells** are filled from the level below:
  - rho, rho u, rho v and E linearly, with each parent's limited least-squares gradient times the child's offset from the parent's centroid. That centroid is the volume-weighted mean of the children's, so the sum over the children is exactly the parent's (*derived*).
  - The mass fractions and specific k and omega are taken as the parent's, so the partial densities and rho k, rho omega are conserved through the interpolated density and still sum to rho.
  - If any child is then not admissible, that parent's children all take its value (zero slope).
- **Cells leaving a level** need nothing: the level below already holds their average.

## Measurements on the composite

- Contents (mass, energy, axial momentum) are level 0's after averaging down.
- Boundary exchanges (supply and outlet fluxes, wall force and heat flow, the igniter, body force, heating, clipped and prescribed turbulent energy, the correction energy) are summed over the uncovered cells of every level, each level with its own step. With refluxing this closes the composite's budgets to rounding, which criterion 2 checks.
- Injector pressure, exit Mach and the exit-plane thrust use the finest level covering each face.

## Threads and determinism

Each level keeps the engine's rule: every cell and face is computed alone, and every sum keeps the serial order. The registers are indexed by coarse face and summed in a fixed order (fine substep, then fine face from low to high j or i). Regridding is serial and depends only on the state. The result is the same bit for bit on any number of threads (criterion 5).

## Criteria (stated before code)

0. **Geometry refactor (step A).**
   - (a) The C1 chamber, `eqtt` 64 x 12 to 1 ms on 4 threads, prints the same history to the last digit as the build before the refactor. The turbulent C1 case of `core_tests` *Threads* has `memcmp`-identical states.
   - (b) ctest results are unchanged: 8 of 9 pass, with transport_verification's two known radial-order failures at their recorded values (1.218, 1.444).
   - (c) For a cylinder and for a cone, at uniform and at stretched rings, the coarsened level-1 mesh matches the directly built base mesh within 1e-14 relative in every array.
   - (d) On the RL10 development contour at levels 0 to 3, every coarse cell closes: its face area vectors sum to (0, its `radialPressureMeasure`) within 1e-14 relative.
1. **Free-stream preservation on a refined curved mesh.**
   - (a) Quiescent gas in the RL10 development contour with closed valves and the outlet at the same state, on the base plus two levels in a fixed pattern that crosses the throat, the wall, the axis, the injector face and the outlet. After 1,000 base steps: largest |u| below 1e-10 a, and largest |p / p0 - 1| below 1e-12, on every level. The single-level run's values are printed beside.
   - (b) Uniform axial flow in a cylinder with stretched rings and slip walls, fixed refinement through the middle of the domain. Inside the region no wave from an end can yet have reached, the same bounds hold.
2. **Budgets across coarse-fine faces under subcycling.** The C1 chamber, `eqtt` on 32 x 6 plus one level, with dynamic tags and regridding, to 8 ms:
   - mass and energy budget errors below 1e-12 of the gas's largest mass and largest |energy| in the run, and no more than 10 times the single-level 32 x 6 run's errors normalised the same way;
   - the axial momentum balance error below 1e-12 of the largest gross momentum exchange, likewise;
   - the igniter's delivered energy exact to 1e-14 relative.
3. **Waves through a coarse-fine boundary.** Inviscid, one species, a cylinder with slip walls; level 1 fixed over the middle third.
   - (a) A smooth acoustic pulse (Gaussian, amplitude 1e-3 p0, half-width 8 base cells) crossing the interface both ways, on base meshes of 64, 128 and 256 columns. The error against the uniform finest mesh converges at observed order 1.8 or better. The reflected wave, measured where the pulse came from after it has gone through, is below 1e-3 of the incident amplitude on 256 columns and falls with resolution.
   - (b) A Sod shock tube with the shock and the contact passing through the interface into the coarse level. The density's L1 error against the exact solution lies between the uniform coarse and uniform fine meshes' errors, and converges at observed order 0.8 or better over the three meshes. The shock is within one base cell of its exact place. Behind it, upstream of the interface, the pressure departs from the uniform-coarse run by at most 1e-3 of the shock's jump once the transient of crossing has passed (5 base cells behind the shock).
   - (c) Ghost width 3 gives the same active-cell states as width 2, bit for bit, in (a) and (b).
4. **C1 on base plus levels against the uniform finest mesh.** C1 `eqtt` on 32 x 6 plus one level against the uniform 64 x 12 run, both to 8 ms: settled c*, vacuum Isp and vacuum thrust within 0.1%, light-off time within one base step plus 1%. With the tags and thresholds fixed in step D before the comparison is run.
5. **Cost and determinism.**
   - (a) The time per cell update (cell updates summed over levels) is at most 25% above the uniform mesh's on the same machine and threads, as the plan budgeted.
   - (b) 1 and 4 threads give `memcmp`-identical states and budgets on the criterion 2 case.
6. **The RL10 development chamber on AMR.** The pre-vaporized RL10 development case on a 150 x 16 base with the static and dynamic tags, from a cold start to settled thrust on the Mac:
   - the measured wall time to settled thrust (vacuum thrust within 0.1% of its settled value);
   - the binding step, by level and place (the step-limit tallies);
   - the decision rule for implicit time applied (`TECHNICAL_PLAN`, *Lightweight engine*).
   - This criterion is a measurement, not a pass mark. The case stays development only, never the reported prediction.

## Order of work

- **A.** The explicit-geometry `Mesh`, its coarsening constructor and the restricted loops: criterion 0. Commit.
- **B.** Two levels in a fixed pattern with subcycling, registers and refluxing: criteria 1 and 3, then 2 with the pattern fixed. Commit.
- **C.** Tags, regridding and interpolation: criterion 2 as stated (dynamic), 5(b). Commit.
- **D.** C1 on base plus one level: thresholds fixed and recorded, then criteria 4 and 5(a). Commit.
- **E.** The RL10 development chamber on AMR: criterion 6, and the decision rule. Commit, and update `TECHNICAL_PLAN`'s budget with the measured rows.

## Deviations from the plan's text

- **Active regions are unions of parent cells, not rectangles.** The plan says patches are rectangles of index space so that each is a small structured mesh for the existing kernels. Full-domain level arrays give the kernels their structured neighbours anyway, and interval lists make clustering into rectangles unnecessary. A rectangle cover would refine whole blocks of rows to hold a two-row wall strip. Nothing in the conservation argument depends on rectangles.
- **Three new criteria (0, 5, 6).** They do not relax the plan's four (here 1 to 4); they record the refactor's safety, the cost the budget assumed, and the measurement the decision rule needs.

## Literature

- M. J. Berger and P. Colella, "Local adaptive mesh refinement for shock hydrodynamics", J. Comput. Phys. 82, 64 to 84 (1989): subcycling (p. 67), proper nesting (p. 68), refluxing (eqs. 1 to 3, pp. 69 to 70), shock tagging apart from the error estimate (pp. 75 to 76). Full text read on 5 October 2026 (TECHNICAL_PLAN).
- M. J. Berger and J. Oliger, "Adaptive mesh refinement for hyperbolic partial differential equations", J. Comput. Phys. 53, 484 to 512 (1984): known by summary only.
- Conservative linear interpolation to new fine cells with a limited gradient and a fallback to constant: standard practice in block-structured codes (AMReX's and Chombo's conservative interpolators), stated here from general knowledge, not re-read.
