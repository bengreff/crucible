# Adaptive resolution (TECHNICAL_PLAN *Lightweight engine*, order of work item 5)

Status (14 October 2026, 01:10 CDT): step A implemented, results below. 0(a), 0(b) and 0(c) are met; 0(d) is not met as written, with the decision recorded. Step B is under way in two parts (*Deviations*). B1 is implemented and judged (*Step B1: results*): 1(a) and 1(b) are met. 3(a)'s reflection is met and its order is not met as written: 1.59 and 1.53, the base scheme's own order on that pulse. 3(b) is met except its pressure peak, which is not met as written: a wave shed by the shock's crossing, whose peak is fixed by the scheme. Criterion 2 is run on B1 with the pattern fixed (*Step B1: criterion 2 results*). The mass, energy and momentum budget errors are 1.3e-14, 2.7e-15 and 9.3e-16, met against 1e-12. The comparisons with ten times the single level's are not met as written for mass and momentum: two rounding walks compared at one instant, alike over the run. The igniter is 1.4e-14 relative, not met as written against 1e-14. B2 and 3(c) come next. Steps C to E are not started. The design and criteria were committed at 22:34 CDT by author time, before any AMR code; this line said "about 22:45" until 22:36.

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
  - a coarse radial face is the polyline of two fine radial faces: area vector summed; the middle radius is the shared fine node. Its reference radius r_f is f times the column radius R_col that also places the centroids (the cell's r is its fraction moment times R_col, `Mesh::radialFaceRadius`). R_col is 1.5 times the column's first-moment factor over its volume factor, both exact integrals, so the coarse R_col is the same ratio of the summed factors and the centroids stay between their faces. The r^2 moment (f^2 times the trapezoid mean of R^2) takes the mean of the two fine columns';
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
- **Igniter.** The power density is fixed on the base: P over the base igniter cells' volume. Every fine cell over a base igniter cell takes that density. The energy delivered by the composite is then exactly P times the time, whatever is refined. (Superseded on 14 October by *Step B1: how criterion 2 is run*, (c): the region is the finest level's cells by centroid, a coarser cell's igniter volume is the sum of its children's, and each owned cell takes P times its igniter volume over the total.)
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
   - (c) For a cylinder and for a cone, at uniform and at stretched rings, the coarsened level-1 mesh matches the directly built base mesh within 1e-14 relative in every array, with one stated exception. The radial face's r^2 moment, f^2 (r0^2 + r1^2) / 2, is the trapezoid mean of R^2 along the column, not an exact integral. Its coarse value is the mean of the two fine columns' (the trapezoid over the finer nodes), which on a cone differs from the direct coarse value by exactly f^2 (r1 - r0)^2 / 8; the test checks that difference within 1e-14 relative. (Amended 13 October, 22:37 CDT, before any code: found while reading the formulas for step A.)
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

- **A.** The explicit-geometry `Mesh`, its coarsening constructor and the restricted loops: criterion 0. Commit. (The restricted loops moved to B; see *Deviations*.)
- **B.** Two levels in a fixed pattern with subcycling, registers and refluxing: criteria 1 and 3, then 2 with the pattern fixed. Commit.
- **C.** Tags, regridding and interpolation: criterion 2 as stated (dynamic), 5(b). Commit.
- **D.** C1 on base plus one level: thresholds fixed and recorded, then criteria 4 and 5(a). Commit.
- **E.** The RL10 development chamber on AMR: criterion 6, and the decision rule. Commit, and update `TECHNICAL_PLAN`'s budget with the measured rows.

## Results: step A, the explicit-geometry mesh (13 October 2026)

`Mesh` now holds every geometric quantity the kernels read as arrays: per axial face the area and the ring middle; per radial face the area vector, the middle, r_f and the r^2 moment; per column the integrals of R^2 and R^3 over z; and the side wall as the finest polyline (`wallRadius`, `wallDz`, `wallStride`, `Mesh::wallRadiusAt`). The accessors read the arrays, so no call site changed. The direct constructor fills the arrays with the former formulas in the same arithmetic. `Mesh::coarsened` builds the level below as designed above. `walls.cpp` and the side wall's neighbour in `neighbours()` read the polyline. `radius`, `fraction` and `dz` remain, holding a coarsened mesh's stations, so `core/spray.cpp` reads them unchanged. Tests: `tests/amr_tests.cpp` (ctest `amr_verification`).

- **0(a): met.** C1 `eqtt` 64 x 12 to 1 ms on 4 threads (Table A, wall functions `law`, the `hp` corrections). The history files before and after the refactor are identical byte for byte (`cmp`), and so are the final states dumped through `CRUCIBLE_STATE_DUMP` (98,304 bytes). The baseline used the build with the y+ accessor and the end-of-run fix of `chamber_study` (both below), taken before any geometry edit.
- **0(b): met in results; the criterion's count was wrong.** Every test except turbulence_pipe_verification ran on both builds (the multi-hour pipe case was stopped unrun on both). Each gave the same result before and after the refactor.
  - Before: 9 of 11 pass.
  - After: 9 of 12 pass. The extra test, amr_verification, fails on 0(d) alone.
  - transport_verification fails in both builds at its recorded radial orders, 1.218 and 1.444.
  - spray_verification fails in both on SPRAY.md's S0(b): 1,245 J/kg against a 1,000 J/kg limit, recorded with the spray commit on main.
  - The criterion's "8 of 9" was written before the rebase onto that commit brought in spray_verification. It counted the tests other than the two longest (pasr, turbulence_pipe).
  - pasr and core_slow pass in both builds. Logs are in the scratch tree, `build/prof/amr0/ctest_*.log`.
- **0(c): met.** Cylinder and cone (r from 0.05 to 0.02 m over 0.1 m), uniform rings and stretching 1.7347: 32 x 16 coarsened against 16 x 8 built directly. Every array is within 7.5e-16 of its largest magnitude. The radial face's r^2 moment differs by -f^2 (r1 - r0)^2 / 8, as amended before the code, to within 1.9e-16. The side-wall distances match within 6.4e-16.
- **0(d): not met as written.** I chose the denominator with the test code, before its first run: the largest component among the cell's four face area vectors. The largest closure errors on the RL10 development contour (base 150 x 16, stretching 1.7347) are:

  | Level | Mesh | Relative to the faces (limit 1e-14) | Relative to pi R^2 (printed) |
  |---|---|---|---|
  | 3, built directly | 1200 x 128 | 8.1e-14 | 7.6e-16 |
  | 2, coarsened | 600 x 64 | 3.9e-14 | 6.7e-16 |
  | 1, coarsened | 300 x 32 | 1.9e-14 | 5.7e-16 |
  | 0, coarsened | 150 x 16 | 8.3e-15 (met) | 5.5e-16 |

  - **The cause** lies in the direct formulas, which the refactor does not touch; it is not in the coarsening. Level 3 is built directly. The worst cell is in the bell, at R = 0.32 m, dz = 1.25 mm. The error is in the axial sum: A(i+1, j) - A(i, j) + the radial faces' z components. Each of those terms is computed from values the size of pi R^2, such as f R1 - f R0 inside -pi (r0 + r1)(r1 - r0). Their rounding is therefore a few epsilon of pi R^2, which is about epsilon R / dz of the faces. A directly built uniform 1200 x 128 mesh gives 4.4e-14, and the 150 x 16 base 8.2e-15. Coarsening sums the fine faces, so the absolute error carries over while the faces double: the ratio halves at each level, as the table shows.
  - **The consequence** for quiescent gas is a spurious force of p times the closure error. At level 3's worst cell that is about 6e-10 N on about 7e-6 kg (*derived*, at 3 MPa and 3 kg/m^3). The resulting speed after 1,000 steps of 1e-8 s is about 1e-12 a, against criterion 1's bound of 1e-10 a.
  - **Decision (13 October, 22:56 CDT).** No change to the arithmetic in step A, because 0(a) requires the same bits. 0(d) stays recorded as not met. Criterion 1 measures the consequence directly and is unchanged. If criterion 1 fails on this, the fix is already known: cancellation-free forms, -pi f^2 (R0 + R1)(R1 - R0) for the radial z component and pi R^2 (f_{j+1} - f_j)(f_{j+1} + f_j) for a stretched axial face. That fix would be a separate commit, checked against 0(d) with the same denominator.
- **Printed beside, not judged.** Level 0 coarsened from level 3 against the 150 x 16 mesh built directly: the largest cell-volume difference is 2.1% and the largest side-wall face-area difference 0.74%. These are the chords' error on the throat arc, which the polyline removes. The throat station's area is identical, because both meshes place a station on the throat.

## Step B1: how criteria 1 and 3 are run (recorded 13 October, 23:59 CDT, before any test code)

I wrote these readings before the tests, so they cannot follow the results. The tests are `tests/amr_subcycling_tests.cpp` (ctest `amr_subcycling_verification`).

- **1(a).** The setup:
  - **Case.** The RL10 development wall (as 0(d)) as a Chamber: base 150 x 16 at stretching 1.7347, with levels 1 and 2 (the finest, 600 x 64, built directly).
  - **Gas.** Perfect gas (gamma 1.4, R 287.05 J/kg K), inviscid, so every wall is a slip wall.
  - **Fill and outlet.** 1 MPa, 300 K, at rest; outlet relaxation 0.
  - **Supplies.** Two supplies, an inner core and an outer annulus, whose valves open at 1 s, after the run. The injector face therefore runs its supply-face path at zero flow beside the plate.
  - **The pattern.** Level 1 covers the base's top four rows in every column (the wall strip from the injector face to the outlet), and every row of the 20 base columns centred on the throat (the narrowest station). Level 2 covers level 1's top six rows in every column, and every row of the 36 level-1 columns centred on the same throat.
  - **Judged.** After 1,000 base steps: the largest |u| / a0 and |p / p0 - 1| over each level's active cells (on level 0, all cells).
  - **Printed beside.** The base's mesh alone and the finest mesh alone, 1,000 of their own steps each.
- **1(b).** The setup:
  - **Case.** A cylinder of radius 0.035 m and length 0.6 m: base 400 x 8 at stretching 1.7347, transmissive ends, slip wall.
  - **Flow.** Uniform axial flow at Mach 0.3, 1 bar, 300 K.
  - **The pattern.** Level 1 covers base columns 150 to 249 and base rows 2 to 5. Level 2 covers base columns 170 to 229 and base rows 3 to 4. Both levels therefore have coarse-fine faces in both directions, on stretched rings.
  - **Judged.** After 400 base steps: |u - u0| / a0 (both components) and |p / p0 - 1|, over the cells farther than (u0 + a0) t from the inlet and (a0 - u0) t from the outlet.
  - **Printed beside.** The whole domain.
- **3(a): an amendment and two readings.**
  - **The amendment.** A Gaussian with a half-width of 8 base cells cannot start inside a 64-column mesh's left third (21 columns) or end inside its right third. Read as half the width at half maximum, the half-width leaves the tail above 1e-3 of the peak over ±25 columns; even the narrowest reading (the e-folding half-width) leaves it above 1e-3 over ±21. The pulse would lie across both interfaces and both ends from the start, so neither "crossing" nor "where the pulse came from" would mean anything. The amendment keeps the pulse (Gaussian, amplitude 1e-3 p0, half-width at half maximum 8 cells of the coarsest base) and the three resolutions, and makes the cylinder three times as long. The base meshes become 192, 384 and 768 columns, with thirds of 64, 128 and 256. The pulse starts at the middle of the left third, where its tail at the third's ends is 2^-16 = 1.5e-5 of the peak.
  - **Reading 1: the width is fixed in length.** It is 8 cells of the coarsest base, and so 16 and 32 cells on the finer bases, as convergence requires.
  - **Reading 2: each run has its own uniform finest mesh.** "The uniform finest mesh" is the uniform mesh at the run's finest level's resolution: 384, 768 and 1,536 columns.
    - **The error** is the L1 norm of the density at the base's resolution: the uniform run averaged onto the base cells against the AMR run's composite, taken through its average-down.
    - **When.** The time at which the pulse's centre reaches the middle of the right third, 2L / (3 a0) after the start.
    - **Judged** on both pairs of successive meshes.
    - **Why.** A single common reference (the finest of the three) would carry its own error into the coarser comparisons and inflate the observed order. Each run's own uniform finest gives a difference of order h^2 if the coarse-fine interface is second order, and of order h if it is not, which is what the criterion tests.
  - **The pulse** is a simple wave moving toward the outlet. It is isentropic, with the incoming Riemann invariant u - 2a / (gamma - 1) exactly constant. A linear pulse would carry a left-moving part of order A^2 / p0 = 1e-3 A, the size of the bound. The initial values are cell averages by 6-point Gauss quadrature in z, on every level and on the uniform runs, so the composite starts consistent.
  - **The reflection.** The largest |p - p_N| over the left third's base cells, divided by A, where p_N is the uniform base run's pressure.
    - **When.** At 300 equal output times from the start to L / a0. By then the reflection of the second crossing has reached the middle of the left third, and the outlet's reflection of the pulse is still five sixths of the length away from it.
    - **Why it isolates the reflection.** Both runs carry the incident pulse identically there, so the difference is what came back from the interfaces.
    - **Judged** on 768 columns (below 1e-3), and as falling over both pairs.
  - **Gas and mesh.** Perfect gas at 1 bar and 300 K, at rest. Radius 0.035 m, length 1.8 m, 2 base rows, slip wall, transmissive ends. Level 1 covers the middle third, all rows.
- **3(b).** The setup:
  - **Case.** The shock tube (Sod, as core_tests: the diaphragm at the middle, length 1, 2 rows) on 64, 128 and 256 columns, level 1 over the middle third, run to t = 0.25. The shock is then at 0.938 and the contact at 0.732, both past the interface at 2/3. The rarefaction's head is at 0.204, past the other interface.
  - **The L1 error** is the density's against the exact solution at the cell centres, on the composite (the finest level covering each place). The order is judged on both pairs.
  - **The shock's place** is where the density crosses the middle of its two exact states, interpolated between centres.
  - **Reading.** "Behind it, upstream of the interface" is read as the coarse cells between the coarse-fine interface and 5 base cells behind the shock. They hold the post-shock gas the crossing could disturb. The contact lies inside this stretch, where the pressure is continuous. The other reading (the fine cells from the rarefaction's tail to the interface) is printed beside, not judged.

## Step B1: results (14 October 2026)

The driver is `core/amr.*`. `Flow` gains a constructor from a given mesh, a step split into its two stages (the driver fills between them), per-cell activity, and per-face flux integrals for the registers. Before the driver existed, the split step was checked on the C1 chamber, `eqtt` 64 x 12 to 1 ms: byte-identical history and dumps. The tests ran on the shared Mac on 4 threads under a load of 10 to 19; logs are in `docs/evidence/amrb/`.

- **1(a): met.** After 1,000 base steps (t = 3.73e-4 s, none redone, 88 s of wall), with the throat at base column 32 (R 0.0652 to 0.0627 m):

  | Level | Mesh | Active cells | Largest \|u\| / a0 (limit 1e-10) | Largest \|p / p0 - 1\| (limit 1e-12) |
  |---|---|---|---|---|
  | 0 | 150 x 16 | 2,400 (all) | 2.6e-14 | 2.2e-15 |
  | 1 | 300 x 32 | 3,360 | 8.2e-14 | 8.4e-15 |
  | 2 | 600 x 64 | 10,944 | 2.0e-13 | 2.2e-14 |

  - **Printed beside.** The 150 x 16 mesh alone, 1,000 steps: 1.8e-14 and 2.6e-15. The 600 x 64 mesh alone, 1,000 steps: 1.3e-13 and 1.7e-14. Each level sits at its own mesh's single-level values; the interfaces add nothing measurable.
  - **0(d)'s consequence.** The closure rounding predicted about 1e-12 a after 1,000 steps at level 3; at level 2 it gives 2.0e-13, five hundred times inside the bound. The cancellation-free forms are not needed and are not applied.
- **1(b): met.** Over the judged span z in (0.0603, 0.5675) m, after 400 base steps (t = 1.34e-4 s, none redone, 59 s of wall), on levels 0, 1 and 2:
  - largest |u - u0| / a0: 5.9e-16, 6.3e-16 and 5.1e-16 (limit 1e-10);
  - largest |p / p0 - 1|: 3.3e-16, 4.4e-16 and 3.3e-16 (limit 1e-12).
  Over all active cells the values are the same, so nothing grew at the ends either.
- **3(a): the reflection is met; the order is not met as written.**

  | Base columns | Base steps | Difference from its uniform 2N run (of the density amplitude) | Uniform N against uniform 2N, no AMR | Ratio | Reflection (of the pressure amplitude) |
  |---|---|---|---|---|---|
  | 192 | 1,500 | 4.951e-3 | 9.070e-3 | 0.546 | 3.54e-3 |
  | 384 | 2,700 | 1.644e-3 | 2.923e-3 | 0.562 | 1.15e-3 |
  | 768 | 4,500 | 5.694e-4 | 1.035e-3 | 0.550 | 3.39e-4 |

  - **Reflection.** 3.4e-4 of the pulse on 768 columns (limit 1e-3), falling by 0.32 and 0.29 per doubling: met.
  - **Order.** Observed 1.59 (192 to 384) and 1.53 (384 to 768), against 1.8: not met. (The log prints the limit 1.8 rounded, as "2e+00"; the format is fixed after this run.)
  - **The base scheme's own order on this pulse.** The control (`pulse_control.cpp`, no AMR) runs uniform N against uniform 2N with 3(a)'s pulse, time and norm. It gives orders 1.63 and 1.50, the same as the AMR's. The AMR's difference is a constant 0.55 of the uniform mesh's at all three resolutions: the share of the coarse level's own truncation error the pulse carries after two thirds of its path on the base. The interface adds no measurable loss of order.
  - **Why the scheme is below 2 here.** The minmod limiter clips the slope at the pulse's peak, so the scheme is first order at the extremum. On uniform N the peak has fallen to 0.829, 0.924 and 0.967 of the amplitude: 17%, 7.6% and 3.3% lost, order 1.2 at the peak.
  - **The reading was wrong, not the driver.** The reading said "of order h^2 if the coarse-fine interface is second order". That assumed the base scheme second order in L1 on this pulse, which the limiter denies. Comparing the AMR's difference with the uniform runs' own difference is what isolates the interface, and that ratio is flat. A later amendment could judge that ratio; it is not applied here.
- **3(b): met except the pressure check, which is not met as written.**

  | Base columns | L1, composite | L1, uniform N | L1, uniform 2N | Its place between them (0 fine, 1 coarse) | Shock, base cells from exact | Pressure, judged (limit 1e-3) |
  |---|---|---|---|---|---|---|
  | 64 | 8.355e-3 | 1.1456e-2 | 6.397e-3 | 0.387 | 0.36 | 5.5e-3 |
  | 128 | 4.717e-3 | 6.398e-3 | 3.566e-3 | 0.41 | 0.40 | 3.7e-3 |
  | 256 | 2.613e-3 | 3.566e-3 | 2.042e-3 | 0.37 | 0.41 | 1.4e-3 |

  - **Orders.** Observed 0.825 (64 to 128) and 0.852 (128 to 256), against 0.8: met.
  - **The thirds** are whole base cells: level 1 covers base columns lround(N/3) to lround(2N/3) - 1, so the right interface sits at 0.672, 0.664 and 0.668.
  - **The pressure departure is a wave the shock emits as it leaves the fine level.** In the scratch profiles (`sod_profile.cpp`, in the evidence folder):
    - **Shape.** A dip of 5.5e-3, 6.5e-3 and 6.4e-3 of the jump (64, 128 and 256 columns), centred just left of the interface, inside the fine level at the two finer meshes. Uniform 2N departs from uniform N there by under 1.5e-4.
    - **Source.** It is already there at t = 0.15, after the shock crosses (t = 0.095) and before the contact arrives (t = 0.18): 3.9e-3 of the jump on 64 columns and 6.0e-3 on 256, at 0.654 on 256 columns. That is where a left-going wave from the crossing would be at u - a = -0.34 in the gas between contact and shock.
    - **Scaling.** Its peak does not fall with resolution; its width does. Integrated from 0.55 to 5 base cells behind the shock, |p - p_N| / jump is 6.9e-4, 4.9e-4 and 2.1e-4 in length, falling about as h. The judged values fall only because the dip moves into the fine level as it narrows, so the judged stretch catches less of its tail.
    - **Not a conservation defect.** Mass and energy balance errors of the composite at t = 0.25 are 4.1e-15 and 5.1e-15 on 64 columns (uniform 64: 4.1e-15 and 4.1e-15).
    - **Reading.** A captured shock is a few cells wide. When it passes from cells of h/2 into cells of h, its profile widens and the conserved excess within it is shed as waves. Their strength in length scales as h and their width as h, so their peak is fixed by the scheme and does not converge. This is a property of a shock crossing any jump in cell size, not of the registers. Bounding the peak at 1e-3 of the jump asks for more than any resolution gives; a bound on the integrated departure, or on the peak relative to the uniform coarse run's own shock error, would test the interface. That is a reading for a later amendment, not applied here.
  - **Printed beside, not judged.** From the rarefaction's tail to the interface: 4.4e-2, 2.5e-2 and 1.4e-2. That is mostly the fine level resolving the rarefaction's tail better than uniform N, not reflection.
- **Housekeeping.**
  - The RL10 development wall moved from `tests/amr_tests.cpp` to `tests/rl10_wall.hpp`, for both tests. `amr_verification` gives the same values after the move: 0(c) within 7.5e-16, and 0(d) at 8.1e-14, 3.9e-14, 1.9e-14 and 8.3e-15.
  - ctest `amr_subcycling_verification` therefore reports 5 failures as recorded: 3(a)'s two orders and 3(b)'s three pressure peaks. The failures stay as written rather than being moved to a later amendment's bounds.

## Step B1: how criterion 2 is run (recorded 14 October, before any of its code)

- **Fixed pattern.** The order of work runs criterion 2 in step B "with the pattern fixed"; the dynamic tags and regridding are step C's, judged there with the same bounds.
- **Case.** C1 as `chamber_study eqtt`:
  - Table A local equilibrium, transport, SST with resolved walls (no wall functions, no property corrections);
  - one supply over the whole injector face;
  - the igniter, 0.5 J over 0.2 ms from 0.2 ms in z 5 to 15 mm, r up to 10 mm.
  - Base 32 x 6, one level (64 x 12), to 8 ms, on 4 threads. The single-level 32 x 6 run is computed beside it.
- **The pattern.** Level 1 covers base rows 0 and 1 of columns 0 to 15 (the inner third from the injector face, past the igniter, to the throat at station 16) and base row 5 in every column (the wall strip through the throat to the outlet). Three things therefore lie across a coarse-fine boundary:
  - **the injector face**, rows 0, 1 and 5 fine, rows 2 to 4 coarse;
  - **the igniter**, coarse ring 2's outer child is in the region by its centroid (9.4 mm) while ring 2 itself is coarse;
  - **the throat and outlet**, on the wall row.
- **What the driver adds, before the run.**
  - **(a) Registers.** The registers carry the viscous face fluxes with the convective ones. `transportFluxes` already stores each face's viscous flux and transported fluxes alone, in the same orientation as the convective faces, so the stage integral is area (F_c + F_v), and the moved amount plus area times the viscous transported flux.
  - **(b) Composite accounting.** Each level counts its boundary exchange only where it owns the cell next to the face. A cell is owned where it is active and not covered by the level above. This applies to inlet, outlet and wall faces, the transport fluxes' wall heat and force, and the cell sources (igniter, heating, body force). The composite budget is the base's content (which holds the composite after average-down and reflux) against the sum over levels of their integrated owned exchanges. With no pattern set, `Flow` counts every face as before, so single-level bits are unchanged.
  - **(c) The igniter's region is the finest level's.** The region is the finest level's cells by centroid, as a single-level run on that mesh would choose. A coarser cell's igniter volume is the sum of its children's. Each cell receives power times its igniter volume over the total, and the owned igniter volumes sum to the total on any pattern. On a single level the igniter volume of an igniter cell is its volume, so the arithmetic is the same.
  - **(d) Still refused:** wall functions and the property corrections. Their prescription and correction fluxes need their own handling across levels; that comes before step E, with its own check.
- **Judged** at 8 ms, with the gas's largest mass and largest |energy| taken over every base step:
  - the composite mass and energy budget errors below 1e-12 of those, and no more than 10 times the single-level 32 x 6 run's, normalised the same way;
  - the axial momentum balance error below 1e-12 of the largest gross momentum exchange, likewise;
  - the composite igniter energy (the levels' owned deliveries summed) within 1e-14 of 0.5 J.
- **Printed beside, not judged:** the single-level values, c* and vacuum thrust at 8 ms on both runs, and the wall time.

## Step B1: criterion 2 results (14 October 2026)

The driver is `tests/amr_budget.cpp` (`crucible_amr_budget 8e-3 4`). It is not a ctest because it needs `CRUCIBLE_EQ_TABLE`. It ran on the shared Mac on 4 threads under a load of about 19 to 20. Its log, with Cantera's 47,882 fallback warnings removed, is `amrb/criterion2_b1.log`.

**Before the run: nothing else moved.**
- **Single-level bits.** The C1 chamber, `eqtt` 64 x 12 to 1 ms on 4 threads, gives a `cmp`-identical history and state to the build before criterion 2's edits (0(a)'s check, repeated).
- **B1's tests.** Every printed value is unchanged from *Step B1: results*, including the 5 recorded failures (`amrb/subcycling_recheck_after_c2.log`).
- **ctest** (`ctest_2026-10-14_amr_c2_after.txt`): 9 of 12 pass, the same as after B1 (`ctest_2026-10-14_amr_b1_after.txt`), and every judged line is identical. The three failures are:
  - `amr_verification`: 0(d), as recorded;
  - `transport_verification`: 1.218 and 1.444, as recorded;
  - `spray_verification`: its (b), the medium's enthalpy against NIST, 1.245e3 J/kg against 1e3. That is the spray worker's lane, failing the same way before these edits.

**The run.** It reached 8 ms in 67,813 base steps (none redone) and 854 s of wall. The single-level 32 x 6 run took 67,624 steps and 47 s.

| Judged at 8 ms | AMR | Bound | Single level | Verdict |
|---|---|---|---|---|
| Mass budget error / largest mass | 1.33e-14 | 1e-12 | 6.0e-16 | met |
| the same, against 10 x the single level's | 1.33e-14 | 6.0e-15 | | not met as written |
| Energy budget error / largest \|energy\| | 2.67e-15 | 1e-12 | 4.0e-15 | met |
| the same, against 10 x the single level's | 2.67e-15 | 4.0e-14 | | met |
| Axial momentum balance error / gross exchange | 9.3e-16 | 1e-12 | 1.1e-17 | met |
| the same, against 10 x the single level's | 9.3e-16 | 1.1e-16 | | not met as written |
| Igniter energy, \|E / 0.5 J - 1\| | 1.42e-14 (7.1e-15 J) | 1e-14 | 1.1e-15 (5.6e-16 J) | not met as written |

- **The comparisons with the single level compare two rounding walks at one instant.**
  - **The AMR's errors change sign.** Its mass error at each half millisecond runs from 4.6e-16 down to -1.8e-14 and back up to +1.4e-14. A defect in the registers would drift one way, and its size would be a face's flux over a step: a tenth of the 0.4 kg/s through one face for one 1.2e-7 s step is 3e-5 of the gas's mass. The errors seen are nine orders of magnitude below that.
  - **The largest errors over the run are alike** (printed, not judged):

    | | AMR | Single level | Ratio |
    |---|---|---|---|
    | Mass | 1.96e-14 | 1.51e-14 | 1.3 |
    | Energy | 9.5e-15 | 1.1e-14 | 0.86 |
    | Momentum | 3.5e-15 | 1.6e-15 | 2.1 |

  - **Why the instant fails.** At 8 ms the single level's walk happens to sit near zero for mass (6.0e-16) and momentum (1.1e-17), so ten times its value at that moment is a bound on luck, not on the interface.
  - **The bound as written stays not met**, as in 3(a) and 3(b). For step C, where criterion 2 is judged as stated, a reading is recorded before that run: the largest error over the run against ten times the single level's largest over the run. It is not applied here.
- **The igniter is 1.4e-14 relative: not met as written.**
  - **The setup misread the bound.** The setup above says "within 1e-14 of 0.5 J" and the test judged 1e-14 J. Criterion 2 says "exact to 1e-14 relative", which is 5e-15 J. The test now judges the relative form; the logged run printed the absolute one, PASS at 7.1e-15 J.
  - **What the error is.** Each level adds dt times its owned power to a running sum at every stage. Over the 0.2 ms window that is about 3,400 additions on the base and 6,800 on level 1. An error of 7.1e-15 J is consistent with that rounding; it was not traced further.
  - **For step C.** A compensated sum in the delivered-energy measurement is planned before step C's run. It does not touch the state.
- **The AMR and single-level contents differ because the wall strip is resolved at 64 x 12.**
  - At 8 ms the energy content is -425.73 J on the AMR and -423.78 J alone.
  - At 50 µs (`amrb/c2_uniform_content_50us.log`, from `amrb/c2_uniform_content.cpp`):
    - every run starts at -0.2450 J;
    - the AMR is at -0.3750 J, uniform 64 x 12 at -0.3792 J and uniform 32 x 6 at -0.3069 J;
    - uniform 64 x 12 loses 7.43 kW to the 600 K wall and uniform 32 x 6 loses 3.77 kW.
  - **Reading.** The resolved-wall heat flow at these wall distances doubles when the wall cell halves, so it is not mesh converged. The AMR's wall strip is the fine mesh's, and its content follows the fine wall.
  - **Why the budgets are untouched.** Each run is judged against its own exchanges.
  - **The two uniform runs start 2e-4 apart in mass** (`c2_uniform_content_1ms.log`). The AMR's base is the 64 x 12 mesh coarsened, which follows the contour more closely than a directly built 32 x 6; the single-level run beside the AMR uses that same coarsened mesh.
- **Observation: light-off passes Table A's 3500 K top on the AMR, not on either uniform mesh.** From a scratch scan to 0.5 ms (`amrb/c2_hot_cells.cpp`, log `c2_hot_cells_0.5ms.log`):

  | Cells | Peak T | Cell-steps above 3500 K |
  |---|---|---|
  | Owned base cells, in the convergent section (base column 11, between the chamber's end at column 8 and the throat at 16), ring 3, near z 0.060 m, r 0.012 m | 3523.9 K | 3,371 |
  | Filled level-1 cells (not active) | 3529.3 K | 16,888 |
  | Owned level-1 cells | 3390.0 K | none |
  | Covered base cells | 3377.0 K | none |

  - **When.** It starts at 0.340 ms.
  - **Fallback.** These states take Table A's counted direct-Cantera fallback: 4,823 lookups on level 0 and 48,860 on level 1 over the run. B1 reacts every cell of a level, filled ones included. All of the warnings come before 1 ms.
  - **The uniform runs.** 32 x 6 and 64 x 12 to 1 ms print no fallback at all.
  - **Reading, not checked.** Here rows 0 and 1 and the wall row are fine and rings 2 to 4 coarse, so the flame crosses from fine rows into coarse ones during light-off. The filled cells are a limited linear reconstruction from those coarse cells, which can overshoot a coarse cell's temperature.
  - **Where it is judged.** It is a difference in the solution, not in the budget, and it belongs to criterion 4 (the AMR against uniform 64 x 12), where light-off time is judged.
- **Printed beside, not judged:**

  | | c* | Vacuum thrust | Outlet mdot |
  |---|---|---|---|
  | AMR | 2474.44 m/s | 1677.93 N | 0.400053 kg/s |
  | Single level | 2477.43 m/s | 1681.83 N | 0.400000 kg/s |

  - **Speed.** The AMR's wall time is 18 times the single level's. B1 runs every level's kernels over all 768 level-1 cells, 256 of them active, twice per base step: 9 times the single level's cell updates.
  - **Cost.** That is B2's business, and the cost is 5(a)'s, judged after B2.

## Step B2: setup (recorded 14 October, before any of its code)

B2 replaces B1's whole-level loops with loops over lists of cells and faces. It is checked bit for bit against B1 on the active cells, then 3(c). The regions are written in terms of N_k: the cells of a level within Chebyshev distance k (in its own index space) of an active cell, so N_0 is the active set A. The ghost width g is the depth that is filled.

- **The ghost width is 3 with SST, not 2.** *Levels, active regions and ghost cells* chose width 2 "from the widest stencil". That reading took the flux stencils only and missed the Strang-split SST source, which B1 applies to every cell before the stage.
  - An active cell beside the interface reads ghosts at distance 1 and 2 (the ghost's MUSCL slope, and its least-squares gradient for the viscous face). Those ghosts must therefore carry the first half-source, as they do in B1.
  - The half-source at a distance-2 ghost uses that ghost's gradients (strain, cross-diffusion), which read its face neighbours at distance 3.
  - So g = 3 when the level has turbulence, and g = 2 for inviscid and laminar flow. The driver refuses turbulence with g < 3. Proper nesting is checked with g.
- **Each filled cell's temperature guess becomes its parent's.**
  - `Flow::refresh` solves for T by Newton from the cell's last temperature, and stops at |T_next - T| <= 1e-13 T. The answer therefore depends on the guess at the last bits, and so on the cell's whole history of refreshes. This holds for the perfect gas too: one Newton step from any guess returns e / cv only to rounding.
  - In B1 every cell is refreshed in every refresh. In B2 a ghost at distance 3 is refreshed at the step's start but not in the stages. Its guess history then differs from B1's, and through the gradients and the half-source the active cells would differ in the last bits. No restriction of the loops can be bit-identical to B1 without this change.
  - The fill now sets each child's guess to the parent's temperature. A filled cell's temperature then depends only on its fill, in either mode. The active cells' histories are identical anyway, because every refresh covers A.
  - **B1 changes with it at rounding level.** So B1 stays selectable (`Amr::setRestricted(false)`), and B1's checks are re-run as B1 with the reset (below).
- **The regions, per kernel.** Each set follows from what the next kernel reads (MUSCL slopes and profiles, and the one-sided boundary slopes, reach distance 2 at most from a face of A; the gradients read face neighbours; the face fluxes read both cells' slopes, gradients and properties):
  - **The fill** covers the parents P_l of the cells in N_g \ A (all four children of each, since the limiter and the admissibility fallback act on the four together). Level l - 1's time-interpolated values are formed over P_l and its face neighbours.
  - **stableDt** (at the base step, and with turbulence again after the fill):
    - with turbulence: refresh and properties over N_g, gradients and eddy viscosity with the source coefficients over N_{g-1};
    - without: refresh and properties over A;
    - the step over A, as in B1.
  - **firstStage:**
    - the half-source over N_{g-1};
    - update, clip and admissibility over A;
    - the copy that keeps a non-active cell's state over N_g \ A.
  - **rhs:**
    - refresh over N_2;
    - axial slopes and radial profiles over N_1, column by column on its j-intervals;
    - the convective faces with an active cell on either side;
    - gather and finish over A;
    - boundary rates as B1 (owned cells only).
  - **transportFluxes:**
    - properties over N_2 (with the wall-face omega of the wall cells among them);
    - field values over N_2, gradients over N_1;
    - eddy viscosity over N_1;
    - the viscous faces with an active cell on either side;
    - the per-cell sum over A.
  - **secondStage:**
    - the update over A, the copy over N_g \ A;
    - with turbulence: refresh and properties over N_2, gradients and eddy viscosity over N_1, and the half-source over N_1.
  - **Why N_1, not A, for the end of the step.** Level l's ghosts next to its active cells (N_1) are face neighbours of the parents that level l + 1's fill reads, and that fill reads their state at the end of level l's step. So the second half-source, and the afterStep hook (`ReactingFlow::react`, which gains the cell list), run over N_1.
  - **The face registers** are zeroed, recorded and summed into fineSum only on the listed faces. Every coarse-fine face, and every face between a covered and an uncovered cell, has an active cell on a side.
- **Off cells** (outside N_g) are never read or written by the kernels. When the regions are set, their stage and next arrays are set to their state, so the swap at accept leaves them unchanged.
- **A single level is untouched.** With no active pattern the lists are empty and every loop runs over the whole level, as now.

**The checks, in order.**
1. **Single-level bits.** The C1 chamber, `eqtt` 64 x 12 to 1 ms on 4 threads, is `cmp`-identical to the build before B2. ctest gives the same judged lines as after criterion 2.
2. **B1 with the reset.** B1 re-runs 1(a), 1(b), 3(a), 3(b) and criterion 2's case to 8 ms. The printed values are compared with *Step B1: results* and *Step B1: criterion 2 results*; differences at rounding level are expected and recorded.
3. **B2 against B1, bit for bit.** The same five runs in both modes, each writing every level's active-cell state (conserved, partial densities, k and omega), the base step count and the composite budget. The two files must be `cmp`-identical.
4. **3(c).** B2 with width 3 against width 2 on 3(a) and 3(b): `cmp`-identical files.
5. **Printed, not judged:** width 4 against 3 on criterion 2's case to 1 ms, and B2's wall times beside B1's.

## Step B2: what was built (14 October)

The code follows the setup's regions (`Flow::setRegions`, the lists `within_`, `faces_` and `columns_`; `Amr::check` builds each level's fill parents and the cells around them). These points differ from the setup's text or were left open by it:
- **The copy covers every filled cell, not N_g \ A.** The fill writes all four children of a ghost's parent, so some filled cells lie just beyond N_g. The copy that keeps a non-active cell's state through the stages runs over the filled cells (`Flow::filled_`). The off cells are the non-active cells outside that list. Nothing reads a filled cell beyond N_g, and it holds the same values as in B1.
- **The default ghost width follows the level's turbulence.** The constructor sets 3 with turbulence and 2 without. `setGhostWidth` still overrides it, within the refusals (below 2, or below 3 with turbulence).
- **Two whole-level paths remain.**
  - `Flow::mixingInputs` gives the PaSR closure's inputs, which `ReactingFlow::react` forms first. It still refreshes and differentiates every cell.
  - `Flow::transportDerivative`, a test probe, does the same.
  - No check runs the PaSR closure under AMR, so that path is unverified in B2. Its whole-level refresh is also a cost to remove before criterion 5(a).
- **The reaction counters count the cells a level steps.** Table cells, clamped and fallbacks cover N_1 in B2 and the whole level in B1, so B2 prints smaller counts. The clipped-energy sum is the same in both modes, because a non-active cell's entry is zeroed by the copy. The dumps hold no counters.
- **The step's first `stableDt` reads stale ghosts.** With turbulence, `Amr::step` forms the step size after refreshing N_g and forming the eddy viscosity over N_{g-1}. These use the ghosts as the last step left them, and outside N_1 they differ between the modes, because B1's second half-source and `react` also covered them. The step size is taken over A, whose gradients read N_1 only, so it is the same in both modes. The ghosts are filled and recomputed before any stage reads them.

## Step B2: results (run 14 October on backhouse, collected 20 October)

The checks ran on backhouse from 03:15 to 05:01 CDT on 14 October, on 4 threads, before its reboot on 16 October; SSH was down until 20 October, so they were collected then. Trees: the build before B2 at 544acc7 (B2's setup, B1's code), and B2 at 89b034d. Evidence on backhouse in `/home/greff/crucible_b2/b2out/` (`summary.txt`, `queue.log`, the suite, budget and ctest logs, the dumps). All values measured.
- **1. Single-level bits: met.** C1 `eqtt` 64 x 12 to 1 ms: the state dump, the history, the mesh and all eleven field files are `cmp`-identical before and after B2; the logs differ only at byte 884, the wall time. ctest: the same 5 of 14 tests fail before and after (`amr_verification` 0(d), `amr_subcycling_verification`, `transport_verification`, the pipe's 1,500 s timeout, `spray_verification`), and the 17 sorted FAIL lines are identical.
- **2. B1 with the reset: rounding-level differences only, as expected.**
  - 3(a) and 3(b): every printed value is the same in the build before B2, in B1 with the reset and in the B1 results above (for example 3(a)'s orders 1.591 and 1.529, 3(b)'s judged pressures 5.458e-3, 3.657e-3 and 1.384e-3).
  - 1(b): largest |u - u0| / a0 on levels 0 to 2 is 6.460e-16, 5.595e-16 and 3.529e-16 before, 4.354e-16, 5.190e-16 and 5.422e-16 with the reset (the B1 results: 5.9e-16, 6.3e-16, 5.1e-16); |p / p0 - 1| is unchanged at 4.441e-16, 3.331e-16 and 3.331e-16. Limit 1e-10 and 1e-12: met in all three.
  - 1(a): largest |u| / a0 on levels 0 to 2 is 2.681e-14, 7.810e-14 and 2.172e-13 before, 2.639e-14, 7.809e-14 and 2.145e-13 with the reset; |p / p0 - 1| is 2.665e-15, 8.438e-15 and 2.176e-14 in both. Met.
  - Criterion 2's case to 8 ms (the build before, then B1 with the reset): mass 1.086e-14 and 9.350e-15, energy 6.008e-16 and 2.403e-15, momentum 8.559e-16 and 9.337e-16, igniter |E / 0.5 J - 1| 1.354e-14 and 1.421e-14. The verdicts are those of *Step B1: criterion 2 results*: the three absolute bounds met, the ten-times comparisons met for energy and not met as written for mass and momentum, the igniter not met as written. c* 2474.44 m/s, vacuum thrust 1677.93 N, outlet 0.400053 kg/s and 67,813 base steps are identical to the recorded values, and so is the single level (2477.43 m/s, 1681.83 N, 0.400000 kg/s, 67,624 steps).
  - The recorded B1 values (Mac) and these (backhouse, a different compiler and library) also differ in the last digits; the printed integrals agree.
- **3. B2 against B1, bit for bit: met.** The active-cell dumps (every level's conserved state, partial densities, k and omega, the base step count and the composite budget) are `cmp`-identical on 1(a), 1(b), 3(a) at 192, 384 and 768 columns, 3(b) at 64, 128 and 256 columns, and criterion 2's case to 8 ms. Every printed value in the suite and budget logs is the same but the wall times.
- **4. 3(c): met.** B2 with ghost width 3 against width 2: `cmp`-identical dumps on 3(a) at all three resolutions and 3(b) at all three.
- **5. Printed, not judged.**
  - Width 4 against width 3 on criterion 2's case to 1 ms: `cmp`-identical dumps.
  - Wall times, B1 with the reset against B2, on a shared machine (load 20 to 23 on 28 threads, with `ctest_pre` running beside the B2 runs from 03:29): 1(b) 31.6 s against 3.4 s; 1(a) 74.8 s against 25.1 s; 3(a) at 192, 384 and 768 columns 4.1, 12.4 and 176.4 s against 2.0, 6.1 and 20.6 s; criterion 2's case to 8 ms 465.6 s against 222.9 s (the single level 38.8 and 42.9 s beside them). The loads were not controlled, so these are an indication, not criterion 5(a)'s measurement. B2's criterion 2 run is still 5.2 times the single level's wall time.

## Before step C: the igniter's compensated sum (recorded 14 October, before its code)

Criterion 2 on B1 found the delivered igniter energy 1.4e-14 relative from 0.5 J, against 1e-14 (*Step B1: criterion 2 results*). The compensated sum planned there:
- **What changes.** Each level's delivered energy (`Flow::integratedHeat_`) is a running sum of dt times the step's mean igniter power. It becomes a Neumaier compensated sum.
  - The compensation is saved and restored with the other integrals when the driver redoes a step.
  - The diagnostics and the composite budget read the sum plus its compensation.
  - Nothing else moves; the state does not read the sum.
- **Why this should be enough.**
  - A step's energy is E / D times its overlap with the firing window, to a few roundings of that one term. Over N terms of about E / N each, those roundings come to about 1e-16 E / sqrt(N).
  - The overlaps telescope. Each is the difference of two adjacent step times, which is exact (Sterbenz's lemma), so their exact sum is the window.
  - What remains is the running sum's own rounding, about sqrt(N) half-ulps: about 6e-15 relative for 3,400 additions, the size seen.
  - Compensated, that falls to about one rounding of the result. What is left beside it is each level's owned share of the igniter volume over the total, a few 1e-16.
- **The checks.**
  1. Criterion 2's case to 0.5 ms (the igniter fires from 0.2 to 0.4 ms), on B2. Judged: |E / 0.5 J - 1| at most 1e-14 on the AMR. Printed: the single level's value.
  2. The state does not move. That run's dump is `cmp`-identical to the build before over everything except its last 8 bytes, the budget's igniter entry.
  3. Single-level bits: the C1 chamber, `eqtt` 64 x 12 to 1 ms, has a `cmp`-identical state dump. The history prints the igniter energy to 7 digits; any change there is listed.

## Before step C: the igniter's results (run 14 October on backhouse, collected 20 October)

Run from 05:02 to 05:04 CDT on 14 October, on 4 threads, the build before at 89b034d (B2) and after at 720cce7. Evidence on backhouse in `/home/greff/crucible_b2ig/igout/`. All values measured.
1. **Met.** Criterion 2's case to 0.5 ms on B2: |E / 0.5 J - 1| is 1.332e-15 (limit 1e-14), against 1.421e-14 before; the delivered energy is 0.499999999999999 J (0.499999999999993 J before). The single level's |E - 0.5 J| falls from 4.996e-16 J to 0. The AMR's mass, energy and momentum errors at 0.5 ms are the same before and after (6.40e-16, 2.19e-15, 2.25e-16), as are the base steps (4,053, none redone).
2. **Met.** The two dumps are `cmp`-identical over their first 57,408 bytes and differ only in the last 8, the igniter entry.
3. **Met.** The C1 `eqtt` 64 x 12 state dump, history and mesh at 1 ms are `cmp`-identical to B2's single-level run (`c1_post` in *Step B2: results*); the history's igniter column did not change.

With criterion 2's igniter now met, criterion 2's remaining failures as written are the ten-times comparisons for mass and momentum (*Step B1: criterion 2 results*); amending them is the Director's call. Step C can start.

## Step C: setup (recorded 20 October, before any of its code)

Step C adds *Tags and regridding* to the B2 driver. What it builds, and the points that section left open:
- **Interface.** `Amr::setTags(Tags)` turns regridding on; without it the driver is B2 as it stands, with the pattern set by `setActive`. `Tags` holds the regrid interval k_r (2 base steps) and the buffer b (2 cells); the throat's half-width (0.5 r_t) and level (2); the wall strip's target y+ (300); the dynamic thresholds (0.02, 0.05, 0.10, 0.05); and per-species weights for f_H (the H element's mass in each species over the species' mass, supplied by the driver, as for the budget probe; empty turns the f_H tag off). A threshold of zero or less turns its tag off.
- **When.** Once before the first base step, then after every k_r-th accepted base step. The levels are then at the same time and averaged down, so a covered cell holds its children's mean. Regridding is serial and depends only on the state (criterion 5(b)). A redone base step does not regrid.
- **Tags on level l** (0 to L - 1, over its active cells; every cell of the base), each from that level's refreshed primitives:
  - **Across a face** to a face neighbour that is active on level l: f_H when |f_H,n - f_H,c| > 0.02; T when |T_n - T_c| > 0.05 min(T_n, T_c); p when |p_n - p_c| > 0.10 min(p_n, p_c). The difference across one face stands in for |grad| times the cell size, and the smaller of the two values for "T" and for "the smaller pressure". Both cells of the face are tagged.
  - **Vorticity:** |du_r/dz - du_z/dr| Δ > 0.05 a_c, the gradients the fill's unlimited least-squares fit over the face neighbours (the axis mirror below j = 0), Δ = sqrt(V / (2π r_c)), the side of a square of the cell's meridional area.
  - **The throat:** the cells within 0.5 r_t of the throat station in z, while l + 1 is at most the throat's level.
  - **The wall strip:** a cell in the wall row is tagged while its first-cell y+ exceeds the target and l + 1 is at most the throat's level. y+ is the laminar estimate the chamber study reports (u_τ from the resolved wall shear, μ |u_t| / y over the first cell's centroid distance), because the AMR still refuses the wall functions (a deviation from the plan's "from the wall function's u_tau"). At the first regrid, before the flow exists, the whole wall row is tagged to the throat's level, as the plan says.
- **New active sets, from the top down**, so that each nests in the new set below it with ghost width g:
  - For l = L - 1 down to 0, the level-l set is the buffered tags (widened by b cells, Chebyshev, clipped to the domain), together with the parents of the level-(l + 1) cells within g of level l + 2's new active cells (the cells that nesting requires).
  - Level l + 1's new active cells are the children of that set. The nesting then holds by construction; `Amr::check` still verifies it at the next step.
  - Where no tag is set the level empties, and a level above an empty one is empty.
- **The fill of new cells, from the bottom up.** It is the fill's spatial part (the limited least-squares slopes, conservative about the children's volume-weighted centroid; the partial densities, ρk and ρω scaled with the child's density; the parent's value if a child is not admissible; the parent's temperature as the guess), applied from level l's state at the regrid time to the parents of the cells that become active. That code moves into a function both paths share, and `fill` is unchanged in what it computes. A new cell's face neighbours on level l lie within g ≥ 2 fine cells, so they are active on the new level l and hold valid values.
- **Cells leaving a level** keep nothing; the level below already holds their mean. If the pattern does not change, the regrid touches nothing (no fill, no new `check`).
- **The budget.** The composite's content is the base's, which a regrid does not touch. The next average-down replaces a newly covered base cell by its children's mean, which differs from the parent by rounding only. The owned cells, where each level books its exchanges and its igniter share, are rebuilt by `check` after the regrid.

**The checks, in order.**
1. **Nothing moves without tags.** Criterion 2's case with the fixed pattern (B2) to 1 ms: a `cmp`-identical dump to the build before; C1 `eqtt` 64 x 12 to 1 ms single level: `cmp`-identical state.
2. **The fill keeps content.** In criterion 2's case at its first regrid after light-off, for every parent whose children become active: the children's mass, momentum, energy and partial densities within 1e-14 relative of the parent's (printed: the largest).
3. **Nesting under motion.** A unit case on 1(b)'s cylinder with two levels and a tag that moves one base cell per regrid across a level-1 boundary: `check` accepts every pattern, and the free-stream bounds of 1(b) hold on every level after 400 base steps.
4. **Criterion 2 as stated.** C1 `eqtt` 32 x 6 plus one level, dynamic tags (f_H, T, p and vorticity, the throat and the wall strip) and regridding, to 8 ms: mass and energy budget errors below 1e-12 of the run's largest mass and |energy|, and at most 10 times the single-level 32 x 6 run's; the axial momentum balance likewise; the igniter within 1e-14. Printed: the regrids that changed the pattern, level 1's active cells over time, c*, vacuum thrust and wall time beside the fixed-pattern run's.
5. **5(b).** Criterion 2's case on 1 and 4 threads: `cmp`-identical dumps and budgets.

## Step C: what was built (21 October)

The code follows the setup: `Amr::setTags`, `Amr::regrid`, and `Amr::spread`, the fill's spatial part moved out of `fill`. `step` regrids before the first base step and after every k_r-th accepted one. These points were left open by the setup's text:
- **The throat station** is the wall's narrowest station strictly inside the domain (the lowest index on a tie). There is none if the narrowest station is an end. A cell is a throat cell when its centroid's z lies within 0.5 r_t of the station's z.
- **The wall strip's y+** uses the wall segment over the cell's column. u_t is the velocity's component along that segment, and y is the centroid's normal distance to it. After the first regrid, the strip refuses a level without molecular transport.
- **The nesting term** marks a level-l cell when one of its level-(l + 2) descendants lies within g (Chebyshev, clipped to the domain) of a child of a level-(l + 1) cell in the new set.
- **The initial regrid takes one pass per level above the base.** Level l's tags are read over its active cells, so a single pass before the first step builds level 1 only, and level 2 would wait for the regrid after base step k_r. Every initial pass tags the whole wall row. With one level above the base this is the setup's single pass.
- **`Tags::custom`**, a caller's own tag on level l's cell (i, j), was added for check 3.
- **`fillError()` and `fillGross()`** measure the content error over the parents filled at the last regrid that changed the pattern, for mass, both momentum components, energy and the partial densities. `fillError()` takes it relative to the parent's own content, as stated (momentum's magnitude for both components). `fillGross()` takes it relative to the children's gross content, the sum of their |content|. `regrids()` counts the regrids that changed the pattern.
- **The fill's centroid is now taken as offsets from the parent's centroid (abdda6d), a change to `fill`'s arithmetic.** In absolute coordinates the offsets' weighted sum cancels only to rounding of |z|, which is up to about 150 child spacings in C1. Check 2 found this (below). The setup's "`fill` is unchanged in what it computes" held for check 1 (3d8716c) and no longer holds bit for bit: B2's ghost fill now differs at rounding.

## Step C: results

- **Check 1 passes (measured, 21 October, Mac).**
  - The build with step C was compared with the build before it (`2135c11`), neither with tags.
  - Criterion 2's case with the fixed pattern to 1 ms: the dump is `cmp`-identical (8,429 base steps both).
  - C1 `eqtt` 64 x 12 to 1 ms, single level: the state dump is `cmp`-identical.
  - So moving the fill's spatial part into `spread` changed nothing it computes.
- **Checks 2 to 5 and the regressions after the fill change** ran on the Mac on 21 October at abdda6d, 4 threads unless stated. Evidence is in `/tmp/inq2c/c/c45/` and `/tmp/inq2c/c/regress/`, with the first fill's run in `/tmp/inq2c/c/c45_old/` (b2e13ce). All values are measured unless labelled.
- **Check 2: failed on the first fill, then passed after the fix.**
  - The first fill (b2e13ce) gave a content error of 4.6e-14 over the parent and 1.6e-14 over the children's gross content, against a limit of 1e-14.
  - The cause was the centroid in absolute coordinates. Energy carries formation enthalpies and changes sign, so near the throat a parent's energy can be close to zero, and an error taken relative to it is ill-conditioned. That is why both measures are printed.
  - With the centroid taken from the parent's (abdda6d), the driver judges the largest error over all 1,173 regrids of criterion 2's case to 8 ms, not only the first regrid after light-off. It is 5.128e-16 over the parent and 4.298e-16 over the gross content (limit 1e-14). Met.
- **Check 4 (criterion 2 as stated): met on the absolute bounds; the ten-times comparisons fail for mass and momentum, as with the fixed pattern.**
  - Mass budget error over the largest mass: 1.753e-13 (limit 1e-12). Met.
  - Ten times the single level's mass error is 6.02e-15. Not met as written.
  - Energy: 7.457e-15 (limit 1e-12). Ten times the single level's is 4.02e-14. Both met.
  - Axial momentum: 1.092e-14 (limit 1e-12). Met. Ten times the single level's is 1.14e-16. Not met as written.
  - Igniter: |E / 0.5 J - 1| = 2.776e-15 (limit 1e-14). Met.
  - The mass error is 23 times the fixed pattern's (7.540e-15). The first fill gave about the same (1.8e-13), so the centroid is not the cause.
  - Inferred, not tested: the rounding from replacing newly covered base cells by their children's mean, at each of 1,173 regrids.
  - Printed, not judged:
    - 1,173 regrids changed the pattern. Level 1 had 500 to 504 of its 768 cells active from 3.5 ms on, and ended with 500.
    - With tags: c* 2446.43 m/s, vacuum thrust 1676.79 N, 67,699 base steps, 227.1 s of wall.
    - Fixed pattern: 2474.44 m/s, 1677.93 N, 67,813 steps, 228.8 s.
    - Single level: 2477.43 m/s, 1681.83 N, 67,624 steps, 30.8 s.
    - The tagged run's c* is 1.1% below the fixed pattern's. Criterion 4 waits for step D's thresholds.
- **Check 5(b): met.** The 1-thread run (AMR only) and the 4-thread run gave `cmp`-identical dumps, which include the composite budget, over 1,173 regrids. Their printed budgets, c* and thrust also agree. The first fill had also given identical dumps.
- **Regressions after the fill change (abdda6d): rounding-level differences only, as expected.**
  - 1(b): largest |u - u0| / a0 on levels 0 to 2 is 5.936e-16, 3.577e-16 and 4.391e-16, and |p / p0 - 1| is 3.331e-16, 4.441e-16 and 3.331e-16. Met.
  - 1(a): largest |u| / a0 is 2.567e-14, 8.241e-14 and 2.015e-13, and |p / p0 - 1| is 2.665e-15, 8.438e-15 and 2.176e-14, the same as B2's. Met.
  - 3(a): the orders are 1.591 and 1.529, and the reflection on 768 columns is 3.390e-4.
  - 3(b): the judged pressures are 5.458e-3, 3.657e-3 and 1.384e-3.
  - So the five failures as written are B1's, at B1's values.
  - Check 3 still passes: 202 regrids, largest 8.708e-16.
  - Criterion 2's fixed pattern to 8 ms:
    - Budget errors: mass 7.540e-15, energy 2.670e-16, momentum 9.726e-16, igniter 1.554e-15.
    - The verdicts are B2's, except the igniter, which is now met through its compensated sum.
    - c* 2474.44 m/s, vacuum thrust 1677.93 N, outlet 0.400053 kg/s and 67,813 base steps are identical to the recorded values.
  - The same case to 1 ms: the dump differs from check 1's (3d8716c), as the fill change implies.
    - c* 2554.92 m/s, vacuum thrust 1408.02 N and 8,429 base steps are the same.
    - The budget errors change only at rounding: mass -5.503e-16 against 4.586e-16, and energy 2.361e-15 against 1.653e-15.
- **Step C is complete to its five checks.** The open items are criterion 2's two ten-times comparisons, which also fail with the fixed pattern, and the tagged run's 23-fold mass error, whose cause is inferred, not tested.
- **Settled since (below, *the tagged run's mass error*):** the cause was the budget's plain sums, not the regrid; with compensated sums, criterion 2 passes as stated on the tagged run.

## Step C: the tagged run's mass error, the test (recorded 21 October, about 21:35 CDT, before its code; 60c8aa4's text said 04:20, a misreading of the clock)

Check 4's tagged run has a mass error 23 times the fixed pattern's (1.753e-13 against 7.540e-15 at 8 ms; largest over the run 2.070e-13). Its momentum error is 11 times, and its energy error 28 times. Criterion 2 is not settled until the cause is known. There are two candidates.

- **H1, the regrid's rounding** (the inference recorded above). Each refinement fills a parent's children, which hold its content to rounding. The next average-down then replaces newly covered base cells.
  - *Arithmetic:* check 2 bounds each filled parent's error at 5.1e-16 of its content. A cell holds about 1/192 of the gas, and the pattern changes by about one parent per regrid (level 1 moves between 500 and 504 active cells). So 1,173 regrids give at most about 1e-14 of the gas's mass even with every error of one sign (derived). That is below what was measured.
- **H2, the budget's own sums.** Each level adds dt times its owned boundary rates to a plain running sum (`Flow::accept`, `core/flow.cpp:1459`). The composite residual is the base's content less the sum of every level's.
  - Suppose the tags put the whole injector face on level 1 while the base owns most of the outlet. Then level 1's sum climbs towards +ṁ t and the base's towards −ṁ t, and the residual is their difference.
  - ṁ t is 0.4 kg/s × 8 ms = 3.2e-3 kg, 18 times the gas's largest mass (1.78e-4 kg; measured).
  - So each addition rounds at about 1e-19 kg, set by the throughput, not the gas. About 2e5 additions (one per level step) walk to about 6e-17 kg, 3e-13 of the gas's mass (derived). That is the measured size.
  - With the fixed pattern, level 1 owns rows 0, 1 and 5 of the face and only row 5 of the outlet. Its net is a fraction of the throughput, so its sums stay smaller.
  - The state itself would then be conserving all along.

**The test.**
- Each level's integrated mass, energy, correction and momentum exchanges get the igniter's Neumaier compensation (`core/flow.cpp:1453-1456`).
  - The plain sums stay as they are, bit for bit. A compensation term per sum sits beside each.
  - `Amr::budget` reads the compensated values.
  - `Flow::measurements` keeps its plain balance errors, so every existing output stays identical. It adds the compensated ones beside them.
- The regrid also adds up, signed and compensated, each filled parent's mass change (the children's content less the parent's), to measure H1 directly.
- Run criterion 2's tagged case to 8 ms on 4 threads, as check 4 did, with the single-level run beside it. Then run the fixed pattern to 8 ms.

**Judged:**
- **T0, nothing else moves.**
  - The tagged run's c*, vacuum thrust, base steps, regrids and printed contents equal check 4's to every printed digit.
  - C1 `eqtt` 64 x 12 to 1 ms, single level: the state dump is `cmp`-identical to the build before.
  - The fast ctest's 7 pass and 4 fail, with sorted FAIL lines identical to af91015's record.
- **T1, H2's mechanism.** At 8 ms, H2 predicts that level 1's and the base's integrated mass exchanges each exceed 5 times the gas's largest mass, with opposite signs (printed per level).
- **T2, decisive for H2.** With compensation, the tagged run's largest mass error over the run is:
  - at most 2.07e-14, ten times below the uncompensated 2.070e-13; and
  - at most 10 times the single level's largest over the run, the single level also compensated. This is B1's recorded reading.
  - The energy and momentum errors are printed beside, with the same reading.
- **T3, H1 measured.** The summed fill mass change over the regrids, over the gas's largest mass. H1 is excluded as the cause if it is below 1e-14.
- **Criterion 2 as written** (10 times the single level's error at the 8 ms instant) is judged again and printed. It may still fail on the single level's luck at that instant (B1).

**Outcomes.**
- If T2 passes, the excess was the budget's arithmetic, not the state. The compensated sums are then kept, and criterion 2 is read on them.
- If T2 fails, H2 is rejected, and the cause is looked for in the step and the regrid with the per-step residual.

## Step C: the tagged run's mass error, results (21 October, about 21:50 CDT)

The runs are on the Mac, 4 threads, with three jobs at once under a load of about 10. The code is the registered test, uncommitted on 60c8aa4 when it ran; it is committed with this section. Evidence: `docs/evidence/amr_c/` (the logs' judged lines).

**The cause is H2, the budget's own sums. The state was conserving all along.**

- **T1 (measured).** At 8 ms the base's integrated mass exchange is +2.988e-3 kg, 16.8 times the gas's largest mass, and level 1's is −2.813e-3 kg (−15.8 times). They have opposite signs, each above 5 times. Met.
  - The direction is the reverse of the narrative in the registration. The base owns the inflow and level 1 the outflow (the wall strip and throat tags put the outlet's rows on level 1), not the face on level 1. The judged prediction did not depend on the sign.
  - With the fixed pattern the two are −1.39 and +2.37 times the gas, which is why its plain error was 23 times smaller.
- **T2 (measured).** With compensation, the tagged run's largest mass error over the run is 4.874e-15, against 2.070e-13 plain: 42 times lower. The limit was 2.07e-14. Met.
  - Against ten times the single level's largest (5.717e-15, also compensated), the limit is 5.72e-14. Met.
  - Energy 6.337e-15 against 9.094e-15 alone, and momentum 1.384e-15 against 1.636e-15. The AMR's largest errors are now below the single level's.
- **T3 (measured).** The filled parents' summed mass change over the 1,173 regrids is −1.111e-16 of the gas's largest mass. That is below 1e-14, so H1 is excluded as the cause.
- **T0.**
  - The tagged run's c* 2446.43 m/s, vacuum thrust 1676.79 N, 67,699 base steps, 1,173 regrids and both printed contents equal check 4's to every printed digit. Met.
  - C1 `eqtt` 64 x 12 to 1 ms: the state dump is `cmp`-identical to the build before. Met.
  - The fast ctest: the same 7 pass and 4 fail as af91015, and the 15 sorted FAIL lines of the rerun are identical, values included (`amr_c/ctest_mass_error_2026-10-21.txt`). Met.
- **Criterion 2 as stated, on the tagged run: every line passes, the ten-times comparisons included.**

  | At 8 ms | AMR | Bound | Verdict |
  |---|---|---|---|
  | Mass | 0.000e+00 | 1e-12; 10 x single 3.01e-15 | met; met |
  | Energy | 2.169e-15 | 1e-12; 10 x single 1.74e-14 | met; met |
  | Axial momentum | 0.000e+00 | 1e-12; 10 x single 3.64e-17 | met; met |
  | Igniter | 2.776e-15 | 1e-14 | met |
  | Check 2 (fill) | 5.128e-16 | 1e-14 | met |

- **The fixed pattern (B2's case), with the same sums, printed.**
  - At the 8 ms instant, its mass error is 3.619e-15, against ten times the single level's 3.01e-15: not met as written.
  - Its largest over the run is 6.334e-15, against the single level's 5.717e-15 (B1's recorded reading: met).
  - Energy and momentum are met at the instant (8.679e-16 against 1.74e-14, and 0 against 3.64e-17).
  - So the one remaining failure as written is B's fixed pattern at one instant, where the single level's walk sits at 3.0e-16. That is the luck B1 described, now at a fifth of its old size.
- **Kept.** `Amr::budget` reads the compensated sums. `Flow::measurements` keeps its plain errors, and prints the compensated ones beside.
- **Wall time** is 565 s here against 227 s in check 4, from the load. It says nothing about cost.

## Step D: the wall functions and the property corrections across levels, setup (recorded 22 October, about 02:20 CDT, before any of its code)

`Amr` refuses both today (`core/amr.cpp:25-26`), and tier 1 needs both. What the code does on one level, read on 22 October:

- **The corrections are not a face flux.** The corrected k and omega diffusion is formed per cell: the cell's F1, S_n/mu, 1/sqrt(rho) and rho/mu weight the sums of the face values in `axialCorrection_` and `radialCorrection_` (`core/transport.cpp:676-705`). Two cells sharing a face weigh it differently, so even on one level it is not conservative. Its difference from the conventional k sum, Phi_k, enters E and is booked as `correctionEnergy` over the owned cells (`transport.cpp:762`). The design above lists "the property-correction fluxes" among the refluxed arrays, but there is no such flux to reflux. That line is amended here: what is refluxed is the conventional k and omega diffusion the correction replaces, which the registers already carry.
- **The prescription.** After each stage, each wall cell's rho k and rho omega are set from the law at its state, with E unchanged, and their change is booked as `prescribedTurbulentEnergy` (`core/flow.cpp:1400-1454`). A prescribed cell has no k or omega equation (its turbulence derivative is set to 0). Under `Amr`, `Amr::reflux` then changes the coarse cells beside a finer level after their last prescription, in every component.
- **The ledgers.** The prescribed and the clipped turbulent energy are summed over all wall cells and all active cells (`flow.cpp:1447, 1453`), covered cells included. Both move energy between k and the internal energy with E unchanged, so neither enters the energy budget; they are reported.
- **The wall-strip tag** reads the laminar y+ (`core/amr.cpp:369-375`), where the design says the wall function's u_tau.

What step D builds:

1. **Corrections:** Phi_k and Phi_omega stay cell sources from each level's own faces, as on one level, and are not refluxed. The refusal goes for `propertyCorrections`.
2. **Prescription after the reflux:** after `Amr::reflux(l)`, every prescribed cell of level l the reflux touched is prescribed again from the law at its new state (E unchanged, the change booked), and its admissibility checked. Then at the end of every level step each owned prescribed cell holds the law's k and omega at its state, as on one level. Covered cells keep their children's mean, as every covered cell does. The refusal goes for `wallFunctions`.
3. **Ledgers over owned cells:** under `Amr` the prescribed and clipped energies are summed over owned cells only. One level is unchanged (its `owned_` is empty, so every cell counts).
4. **The tag:** with the wall functions on, the wall strip reads the law's y+ at the first cell; otherwise the laminar estimate as now.

Checks, judged as written:

- **D0, nothing else moves.** Without the wall functions and corrections, criterion 2's fixed-pattern dump and check 3's run are `cmp`-identical to 566b12c; C1 `eqtt` 64 x 12 to 1 ms with `law` and `hp` on one level is `cmp`-identical. The fast ctest gives 8 pass and 4 fail, with the sorted FAIL lines identical to the same machine's record: `face_inlet/ctest_2026-10-22_backhouse.txt` on backhouse, `amr_c/ctest_mass_error_2026-10-21.txt` on the Mac.
- **D1, the budgets with both on.** Criterion 2's tagged case (C1 `eqtt` 32 x 6 plus one level, to 8 ms) with `law` and `hp`. Mass and energy budget errors (the energy net of the summed correction energy) below 1e-12 of the gas's largest mass and largest |energy|; the axial momentum balance below 1e-12 of the largest gross exchange; the igniter exact to 1e-14. The single-level 32 x 6 run's errors with `law` and `hp` are printed beside, and the ten-times comparisons are reported.
- **D2, the prescription holds.** At the end of every base step, every owned prescribed cell's rho k and rho omega equal the law's at its state to 1e-12 relative, on every level. *Amended 22 October, about 02:25 CDT, while coding and before any run:* this is judged on criterion 2's fixed pattern with the reaction step off (`crucible_amr_budget wall inert`), not on D1's run. The reaction step after each level step and the regrid's new cells move the state after the prescription on one level too, so on D1's run the bound would fail without any defect. The single level's value is printed beside as the control, and D1's run prints its own value, not judged.
- **D3, determinism:** D1's run on 1 and 4 threads gives `memcmp`-identical states and budgets.
- **D4, the face inlet on levels with both on:** F3's AMR case (149 x 16, 3 levels) with the wall functions and corrections on passes F3's per-ring and per-supply checks to 1e-12 (`docs/evidence/FACE_INLET.md`).
- **Reported, not gated:** D1's c*, thrust and light-off time against the uniform 64 x 12 run with `law` and `hp`. Criterion 4 itself still waits for the tag thresholds.

## Step D: results (22 October, about 04:40 CDT; backhouse, 4 threads, at load about 24)

Run by `bh_stepd.sh` in `/home/greff/crucible_sd` (inquiry2 at b84ae3a, code as 2ff039f), with D0's parent built in `/home/greff/crucible_fi/build` (566b12c, code as 4a7178c). The raw outputs are in `amr_d/` (`d0.txt` to `d4.txt`). All values are measured unless labelled.

**In short** (D0, D1 and D2 are settled in the sections after this one: step D is met, 5 of 5, under the Director's amendments): D3 and D4 are met. D1 and D2 are not met as written: D1 on the igniter, and D2 fails on one level as well, by more. D0 is met on the runs finished; its fast ctest and check 3 were still running at parking (how to score them is below).

- **D0: met on the runs finished. The fast ctest and check 3 are pending.**
  - Criterion 2's fixed-pattern dump to 1 ms is `cmp`-identical to the parent's. The printed budgets, c* 2554.92 m/s, vacuum thrust 1408.02 N and 8,429 base steps are equal too.
  - C1 `eqtt` 64 x 12 to 1 ms with `law` and `hp`, one level: the state dump and the history are `cmp`-identical.
  - Check 3's run (new and parent, dumps compared) and the fast ctest were still running at about 04:40, in that order after the reported 64 x 12 run.
  - To score them: `bash /home/greff/crucible_sd/score_sd.sh` on backhouse, once `out/all_done` and `out/c3_done` exist. Met if three things hold:
    - the ctest gives 8 pass and 4 fail of 12, and the sorted rerun FAIL lines (16) are identical to `face_inlet/ctest_2026-10-22_backhouse.txt`'s;
    - `out/d0_c3.txt` reads "check 3 dump identical";
    - the two c3 logs differ only in wall time.
- **D1: not met as written, on the igniter.** The tagged case to 8 ms, `law` and `hp`, 181 regrids; level 1 ends with 636 of 768 active.
  - Mass 2.441e-15, energy 7.360e-15 and momentum 0 against 1e-12: met.
  - Igniter |E / 0.5 J − 1| = 2.143e-14 against 1e-14: not met. Step C's tagged run without the wall functions gave 2.776e-15 (check 4), and step D's inert fixed pattern with them gives 3.331e-15 (D2's run).
  - The igniter fires from 0.2 to 0.4 ms. By 0.5 ms the run had 121 regrids, and the delivered energy stays at 0.499999999999989 J from then on. The cause is not established. *Guessed:* the igniter's split between the levels while the pattern changes during the deposit.
  - Reported: the ten-times comparisons all hold (mass 2.441e-15 against 5.27e-14, energy 7.360e-15 against 1.94e-14, momentum 0 against 6.11e-18). The single level's errors: mass 5.265e-15, energy −1.938e-15, momentum 6.108e-19, igniter 0.
  - Also met: step C's check 2 (the largest fill error, 4.182e-16) and T2 (the largest mass error over the run, 5.493e-15 against 2.07e-14).
- **D2: not met as written, and the single-level control fails it by more.**
  - The inert fixed pattern to 8 ms, over 25,564 base steps: the largest relative residual is 6.458e-4 on AMR, against 1e-12. The single level, same case: 4.373e-3.
  - So D2's premise ("as on one level") is false for this measure: one level does not hold its wall cells at the law's k and ω at their state to 1e-12. AMR's residual is below one level's.
  - The wall law itself solves to 1e-13 in u+ and keeps no history between calls (`core/walls.cpp:111-145`, read), so its tolerance is not the cause.
  - *Inferred, not tested:* the two causes below.
    - The prescription takes the cell's temperature with its k before the prescription, and E stays fixed. So afterwards the temperature differs by the change in k over cv, and the law at the new state differs.
    - The law's eddy viscosity is μ_w max(ratio, 0), a difference of near-equal terms close to the sublayer. There ρk tends to 0, and a relative residual is ill-conditioned.
  - D1's run shows the extreme case: 1.985e+307 on AMR and 1.773e+307 on one level. That is a cell where the law gives ρk = 0 exactly, against a stored 0.44 (derived: 1.985e307 × 2.2e-308).
  - The inert run's budgets all pass: mass 1.937e-15, energy 3.475e-16, momentum 0, igniter 3.331e-15, the ten-times comparisons included.
  - The 2ff039f driver also printed D1's D2 value as a failed check line. By the amendment it is not judged there, so the driver now judges D2 only with `wall inert` on the fixed pattern and prints it otherwise. That change was made after these runs and is in the results commit.
  - *Proposed for next session, to register before it runs (amending D2 is the Director's call):* check what D2 is meant to show directly. Each owned prescribed cell's state, species, ρk and ρω at the end of every base step should be bitwise what its last prescription (or re-prescription after the reflux) left. Run it on AMR and on one level.
- **D3: met.** D1's run on 1 and 4 threads gives `cmp`-identical dumps, which carry the states and the composite budget. The printed budgets, c* and thrust are equal.
- **D4: met.** F3's AMR case (149 x 16, 3 levels) with the wall functions and the corrections on:
  - per ring: 5.76e-16, 1.71e-15 and 7.32e-15;
  - per supply: 1.01e-15, 3e-15 and 1.2e-15;
  - supply areas: 8.57e-16, 3.06e-15 and 1.19e-15 (limit 1e-12).
  - These equal F3's record with them off, digit for digit (`face_inlet/f3.txt`), and `Amr` now builds with them on.
  - The dry run's G5, G6 and G7 failures are tier 1's standing mesh gaps, as in that record.
- **Reported, not gated:**
  - D1 (AMR 32 x 6 plus one level, `law` and `hp`): c* 2425.21 m/s, vacuum thrust 1633.34 N. The single level 32 x 6 with both: 2459.99 m/s, 1633.47 N.
  - The uniform 64 x 12 with both to 8 ms (136,265 steps): c* 2426.18 m/s, vacuum thrust 1630.67 N. So AMR's c* is 0.04% below it and its thrust 0.16% above. The single level 32 x 6 is 1.39% and 0.17% above (derived). Its log and history are in `/home/greff/crucible_sd/out/u64/`.
  - `crucible_amr_budget` prints no light-off time, so none is reported.
  - Wall time: D1 2,352 s on AMR against 186 s on one level, at load about 24. It says nothing about cost.

## Step D: D0's last runs, and the Director's amendments to D1 and D2 (recorded 22 October, about 22:25 CDT, before the code and runs they govern)

**D0: met.** Scored on backhouse on 22 October (evening) from `/home/greff/crucible_sd/out`; the outputs are copied to `amr_d/`.

- The fast ctest gives 8 pass and 4 fail of 12. Its 16 sorted rerun FAIL lines are identical to `face_inlet/ctest_2026-10-22_backhouse.txt`'s, values included (`amr_d/ctest.txt`, `amr_d/ctest_rerun.txt`).
- Check 3's dump is `cmp`-identical, new against the parent. The driver writes it as `out/c3n/d_c3.bin` and `out/c3b/d_c3.bin`; `score_sd.sh` compared another path and printed a false "differs", so the comparison was made by hand.
- The two check 3 logs differ only in wall time (10.8 s against 11.0 s).

**D2, amended (the Director's decision, relayed 22 October evening).** The 1e-12 relative residual is ill-posed where the law's ρk goes to 0, and one level fails it too (4.373e-3). It is replaced by a direct check of what D2 is meant to show:

- At the end of every base step, every owned wall-function cell's conserved state, mass fractions, ρk and ρω are bitwise (`memcmp`) what that cell's last prescription left, on the stage or after the reflux. `Flow::prescribeCell` keeps the values it leaves beside it.
- Judged on the inert fixed pattern (`crucible_amr_budget 8e-3 4 wall inert`): 0 mismatches over every base step, on every level, and 0 on the single level as the control. If the control fails, D2 is not met, because the check would then not show what it is meant to.
- D1's run prints its count, not judged: the reaction step moves the state after the prescription there.

**D1, the igniter's cause (the Director's decision: up to an hour on it; if it is summation order or rounding, the bound is amended to 1e-13 with the evidence; if it is real, it is fixed).**

- *Inferred from the code (`core/amr.cpp:515-548, 592-615`), not yet tested:* level 1's time chain does not close on the base step's. Take a base step over [t, t + dt]. Level 1 steps over [t, m] and [m, fl(m + dt/2)], with m = fl(t + dt/2). The next base step starts at fl(t + dt), or at the event when it lands on one.
- Each level's igniter power is (E/D) × overlap / dt over its own [time, time + dt]. So level 1's deposit time in a base step is fl(m + dt/2) − t, against level 0's fl(t + dt) − t. The two can differ by one ulp of t, 5.4e-20 s near 0.3 ms (derived).
- Level 1 then delivers about (E/D) × ulp = 1.4e-16 J too much or too little in that base step, scaled by its share of the igniter (derived). Over the deposit's base steps that sums to the order of the 2.143e-14 measured, either as a random walk or as a partly systematic drift.
- **Test 1 (printed, not judged).** For each base step (its last attempt's), sum each level's igniter power times its step in long double. Print three things over the run:
  - each level's |Σ / 0.5 J − 1|;
  - the number of base steps where level 1's sum differs from level 0's by more than 1e-17 J (the products round at about 1e-19 J, derived);
  - those differences' net.
  - The hypothesis predicts level 0 at 1e-16 or below, level 1 near 1e-14, and a count on the order of the base steps in the deposit.
- **Test 2 (a counterfactual build, not committed unless adopted).** Each level's step is handed its end time: level 1's second substep ends at its parent's end, and a landing base step at the event. The igniter overlap uses that end, with dt unchanged. The hypothesis predicts D1's igniter error falls to about 1e-15 or below, with the other budgets unchanged in their leading digits.
- **Reading.** If both tests bear the hypothesis out, the cause is rounding in the time arithmetic, and D1's igniter bound is amended to 1e-13 on this evidence. Test 2's change is then reported as an option, not adopted: it moves criterion 2's and check 3's dumps by ulps, so D0 would need a new baseline. If the tests do not bear it out, the cause is looked for further within the hour, and fixed if it is real.
- D1's rerun on the new code must reproduce its first run's dump (`out/d1_t4.txt`) under `cmp`: the D2 check and Test 1 read the state and do not change it.

## Step D: results under the amendments (22 October, about 23:25 CDT; backhouse, 4 threads, at load about 22 to 26)

Run by `bh_sd2.sh` from 22:28 CDT. Build a is inquiry2 at 71755ce in `/home/greff/crucible_sd2`; build b is 71755ce plus `amr_d/d1_test2_chain.patch` in `/home/greff/crucible_sd2b`. The raw outputs are `amr_d/sd2_d1a.txt` (D1's rerun and test 1), `amr_d/sd2_d1b.txt` (test 2), `amr_d/sd2_d2a.txt` (D2) and `amr_d/sd2_d1a_cmp.txt`. All values are measured unless labelled.

**In short: step D is met, 5 of 5, under the Director's two amendments.** D1 is met only under the amended igniter bound.

| Check | Verdict | Evidence |
|---|---|---|
| D0, nothing else moves | Met | Above (*D0's last runs*) and *Step D: results*. |
| D1, the budgets with both on | Met at the amended 1e-13 | Igniter 2.143e-14. The cause is established below; mass 2.441e-15, energy 7.360e-15, momentum 0. |
| D2, the prescription holds (amended) | Met | 0 mismatches on AMR and 0 on the single level, over 25,564 and 24,923 steps. |
| D3, determinism | Met | *Step D: results*. Build a's D1 dump is `cmp`-identical to that run's. |
| D4, the face inlet on levels | Met | *Step D: results*. |

- **D2: met.** On the inert fixed pattern to 8 ms, no owned wall-function cell on either level ever differed, at the end of a base step, from what its last prescription left: 0 mismatches. The single level gives 0 too, so the control holds. The budgets equal the first D2 run's (mass 1.937e-15, energy 3.475e-16, momentum 0, igniter 3.331e-15).
- D1's tagged run prints 4,324,450 mismatches on AMR and 2,145,386 on one level, counted over the base steps. That is the reaction step and the regrid moving the state after the prescription, as the amendment says, and it is not judged.
- **D1's rerun reproduces the first run:** the dump is `cmp`-identical to `out/d1_t4.txt`, and every printed budget is equal. So the D2 check and test 1 change nothing.
- **Test 1 bears the hypothesis out.** Over the 1,799 base steps with the igniter on, each level's power × step sum gives |Σ / 0.5 J − 1| = 4.825e-17 on level 0 and 2.218e-14 on level 1. Level 1's sum differs from level 0's by more than 1e-17 J in 901 of those base steps, net −1.111e-14 J.
- That net is the composite's shortfall, scaled by level 1's share of the igniter. The composite has |E − 0.5 J| = 1.071e-14 J, against level 1's own 2.218e-14 × 0.5 J = 1.109e-14 J (derived).
- **Test 2 bears it out too.** With each level step handed its end time, the composite igniter error is 0, from 2.143e-14. Level 1's sum differs from level 0's in 0 base steps (net +1.0e-19 J).
- As predicted, mass (2.441e-15), momentum (0), the base steps (67,635), the regrids (181), c* and thrust (2425.21 m/s, 1633.34 N) are unchanged. The energy error moved from 7.360e-15 to 6.599e-15. So "unchanged in their leading digits" does not hold for it as written. It is itself a rounding-level quantity, moved by the ulp changes in the igniter's power.
- **Reading, as registered:** the cause is rounding in the time arithmetic. Level 1's two substeps end at fl(fl(t + dt/2) + dt/2), and the base step at fl(t + dt); the two differ by one ulp of t in about half the base steps. So D1's igniter bound is amended to 1e-13, and 2.143e-14 meets it.
- *Derived, a caution on the 1e-13:* the mechanism's worst case is every differing step erring the same way. That is 1,799 × (E/D) × ulp(t) = 1,799 × 2,500 W × 5.4e-20 s, about 2.4e-13 J, or 4.9e-13 relative: above 1e-13. The measured 2.1e-14 is a partly cancelling sum. The bound holds here by a factor of 4.7, but nothing guarantees it for a deposit over more base steps.
- Test 2's patch makes the igniter exact (0 here), and is about ten lines. It is not adopted, as registered, because it moves criterion 2's and check 3's dumps by ulps and D0 would need a new baseline. *Inferred:* it is worth taking with the next change that re-baselines D0 anyway, before tier 1, since the run's deposit spans more base steps.
- Wall time, at load about 22 to 26: D1 3,159 s on AMR and 93 s on one level; D2 1,117 s on AMR and 31 s on one level.

## Spray on AMR, read from the AMR side (22 October, about 23:35 CDT)

The Director asked this lane to review inquiry's spray on AMR (`SPRAY.md`, *Spray on AMR*, as built in f2f8026 and amended in 19355bb) for conflicts with the AMR code. Read: the `Amr` diff (`core/amr.cpp`, `core/amr.hpp`), `Spray::afterStep`, `Spray::redo`, `Spray::release`, `Spray::owner` and `Spray::mark`. The rebased tree builds (`crucible_amr_budget`, Mac, 23:36).

**No conflict blocks either lane.** The details, each read from the code:

- **The exchange sums in `Amr::save` and `restore` close a hole that step B left.** Without them, a redone base step would keep the failed attempt's exchanges. The gas alone never exchanges on `Amr`, so no recorded AMR result is touched.
- **The order per level is the uniform one.** `Amr::advance` calls the after-step hook (the spray, then the reaction) straight after the level's `accept`, so after the stage's prescription, and before the finer level's substeps, the average-down and the reflux (`core/amr.cpp:527-543`).
  - One asymmetry, minor and not new. A coarse wall-function cell beside a finer level that the reflux touches is prescribed again after its deposit. One the reflux does not touch keeps the prescription from before its deposit, as on one level. The reaction step already has the same split.
  - So D2's bitwise check counts every wall-function cell a deposit lands in, as it counts the reaction's. That is why D2 is judged inert and without spray; it is not a spray defect.
- **D1's test 2 patch would not disturb the spray's clock.** The clock counts each level's calls within the base step, not times (`spray.cpp:1122-1123`). Taking the patch later, at a D0 re-baseline, changes no parcel's move.
- **A1(b)'s gas-only difference under SST (1.562e-13) lies in `Amr::advance`.** Its sequence of `stableDt`, `beginStep`, the stages and `accept` differs from `Flow::step`'s somewhere not yet isolated. The Director amended A1(b) to 1e-12 rather than change `Amr`, so nothing is owed. If an exact match is wanted later, this lane can isolate it and take it with test 2's patch, so that D0 is re-baselined once.
- **`Spray::release` books [t1 − dt, t1], with t1 the base level's time after its step** (`spray.cpp:1133-1135`). The uniform path does the same (`spray.cpp:1087`). fl(fl(t + dt) − dt) need not equal t, so consecutive intervals can overlap or leave a gap of one ulp of t. That is D1's class, on both paths alike.
  - It does not open the ledger: the injected mass is what `release` booked.
  - Only the delivered total against the schedule's integral moves. *Derived:* by at most massFlow × ulp(t) per base step.
  - Carrying the last release's end as the next start would remove it. That is inquiry's call.
- **The clock uses `Amr::top()`, and `redo` handles the edges.** A failure before level 0's hook in the first base step finds no snapshot of this step and leaves the parcels alone. A second failure after a restore restores the same snapshot again. Both were read, not run.
- **The corner jet (A4) is inquiry's to find.** One fact from this side may help. Level 0's parcels move and deposit at level 0's hook, before level 1 steps. A deposit that `mark` resolves into a covered cell lands in level 1's cells before level 1's first substep. One that a level 1 parcel spreads into a cell level 1 does not own lands in the base cell after the base's step. Both are conservative (A2 shows it). This lane can print per-level exchanges near the corner if inquiry asks.

## Fixing the tag thresholds, then criterion 4 (recorded 22 October, about 23:40 CDT, before any of its code)

Step D's last item, as the order of work states it, is "thresholds fixed and recorded, then criteria 4 and 5(a)". Criterion 4's comparison was printed once already, at the guessed thresholds (*Step D: results*, reported, not gated: AMR's c* 0.04% below the uniform 64 x 12 and its thrust 0.16% above). So the thresholds are not chosen by looking at the uniform run. The rule below reads only the AMR runs and picks the value mechanically. The uniform run is judged once, after.

- **The case.** C1 `eqtt` on 32 x 6 plus one level, with `law` and `hp` (tier 1's configuration, as in D1), to 8 ms.
  - Criterion 4 was written when `Amr` refused the wall functions. It is read here with them on, because that is what tier 1 runs.
  - The resolved-wall variant is not run.
- **The runs.** The four dynamic thresholds (f_H 0.02, T 0.05, pressure 0.10, vorticity 0.05) are scaled together by s, for s = 1, 1/2, 1/4 and 1/8.
  - Kept as they are, declared: k_r = 2, b = 2, the wall y+ target 300, and the throat tag, which is capped at level 1 here.
- **The quantities,** defined the same way on every run, the uniform one included:
  - **Settled:** the mean of 21 samples at 7.00, 7.05, ..., 8.00 ms, each taken at the first step ending at or after its time.
    - The criterion fixes the end at 8 ms, so "settled" is read as the last millisecond's mean. Chamber_study's settling test is not applied.
  - **c\*:** the mass-flux-weighted equilibrium stagnation pressure on the fine chamber-end column, times A_t over ṁ_out.
    - On AMR the column is the composite's finest cells: level 1's two cells where level 1 is active over a base row, and the base cell where it is not.
    - Step D's printed c\* read level 0's averaged column, which is not the same quantity at 0.1%.
  - **Vacuum thrust:** the inlet momentum flux plus the wall and body axial forces, summed over the levels' owned faces, as now.
  - **Vacuum Isp:** F_vac / (ṁ_out g0).
  - **Light-off:** the end of the first step (on AMR, the first base step) at which the largest Y_H2O over the owned cells exceeds 0.5. That is chamber_study's threshold, checked here at every step.
- **The rule.**
  - s\* is the largest s in {1, 1/2, 1/4} that meets both of these against the run at s/2:
    - its settled c\*, Isp and thrust each lie within 0.02%, a fifth of criterion 4's tolerance;
    - its light-off lies within one base step plus 0.2%.
  - The thresholds are then fixed at s\* times the defaults, and `Amr::Tags`' defaults change to them.
  - If no s qualifies, they are not fixed. I report which quantity moves and by how much, and the Director decides.
- **Criterion 4,** judged once and as written, at s\*: against the uniform 64 x 12 run with `law` and `hp`, on the same driver and the same definitions.
  - The tolerances are unchanged: 0.1% for the settled c\*, Isp and thrust, and one base step plus 1% for light-off.
- **Printed, not judged:**
  - for each run: level 1's active cells at the end and on average, the base steps, the regrids, each quantity at 8 ms, the composite budget errors and the wall time;
  - for 5(a): the cell updates summed over levels, each level's active cells times its steps, and the wall time per update on both runs.
- **5(a) is not judged here.** Its 25% bar needs a machine under controlled load. These runs share the Mac with other workers, so their cost ratio is printed with the load beside it.

## Criterion 4's driver, and why its runs wait for the spread fix (22 October, about 23:50 CDT)

- **Built.** `crucible_amr_c4` (`tests/amr_c4.cpp`) runs one case per call and prints one `c4` line; `tests/amr_c4_judge.sh` reads the five lines, applies the rule and judges criterion 4. C1 `eqtt` and the stagnation pressure moved from `amr_budget.cpp` into `tests/c1_eqtt.hpp`, unchanged.
- **The move changes nothing (measured).** `crucible_amr_budget 4e-4 2 wall tags amr`, built before and after it, writes `cmp`-identical dumps and the same log.
- **The pipeline runs (measured; a smoke run, not criterion 4).** All five runs went to 1.2 ms on 2 threads, with the window moved to the last millisecond; the judge read them and applied the rule and the criterion. The logs are in `amr_d/c4_smoke_2026-10-22/`. Their numbers are not the registered case and are not results.
  - The fine chamber-end column is 17. It is odd, so base column 8 ends with it and the composite column is well defined.
- **Light-off, an observation from the smoke run (measured).** Under `eqtt` every cell is put at local equilibrium every step, whatever its temperature (`Chemistry::LocalEquilibrium` from Table A). Propellant therefore burns where it arrives.
  - The registered light-off (some cell's Y_H2O above 0.5) falls at 0.0102 ms on all five smoke runs. The igniter starts at 0.2 ms.
  - So for `eqtt` the item measures how fast the first cells at the injector face fill under the supply's ramp (inferred), not how a kernel grows.
  - It stays as registered: it is the criterion's text and chamber_study's definition. Whether criterion 4 wants an event that ignition sets, for an equilibrium case, is the Director's question.
- **The runs wait for the spread fix.** 28d87b1 (main) registers a change to `Amr::spread` that, in its words, "changes every AMR run whose fills meet gradients". C1's fills meet them.
  - The thresholds should be fixed for the spread tier 1 will run, and criterion 4 is judged once. So the five runs start when the fix is on `main`, from the rebased build.
  - None of them has been run to 8 ms.
- **What this lane reruns after the fix,** on the Mac: the step C and D results that inquiry's list (A5, A4, A0 to A3, AMR's ctest) leaves out.
  - Criterion 2's fill content (1e-14), from `crucible_amr_budget 8e-3 4 wall tags`.
  - D1's budgets and igniter. The s = 1 run of criterion 4 prints them on the same case.
  - Step D's recorded results stand as measured on the spread before the fix (7d16d45, step D's commit after the rebase onto 28d87b1).

## Criterion 4: the Director's additions (22 October, about 23:50 CDT, recorded before their code and before the runs)

- **Light-off is judged as registered.** Criterion 4 is a consistency check, AMR against uniform, and an event set by mixing is still a fair comparison (the Director).
  - *Note, dated 22 October 2026:* in a case that is put at equilibrium every step, as `eqtt` is, "light-off" measures how the propellant fills and mixes, not ignition. The igniter cannot influence it.
- **A second quantity, measured and reported, not pass or fail: t90, the time the chamber pressure reaches 90% of its 8 ms value,** AMR at s\* against uniform.
  - **Chamber pressure:** the quantity in c\*'s numerator, the mass-flux-weighted equilibrium stagnation pressure on the same composite chamber-end column (chamber_study's chamber pressure).
  - **The history:** that pressure is sampled at the first step ending at or after each multiple of 10 µs from 0, and at the end. A sample whose column carries no net outflow (mass flux ≤ 0) has no chamber pressure, and counts as below.
  - **t90:** after the last sample below 0.9 p(8 ms), the time where the pressure crosses 0.9 p(8 ms), linear between that sample and the next. Taking the last crossing means an early transient that touches 90% does not count.
  - **Printed:** t90 and p(8 ms) for every run. For s\* against uniform, the difference in t90, absolute, relative and in base steps.
- **5(a)'s wall time excludes the history's sampling.** The sampling's own wall time is printed beside it.

## The spread fix, read from the AMR side (22 October, about 23:50 CDT)

The fix registered in 28d87b1 (`SPRAY.md`, *The fix, registered before its code*) changes this lane's `Amr::spread`. Read against `core/amr.cpp:223-288`.

- **The defect is the AMR's, and the design removes it.** *Derived:* with one φ, every child is one combination of the parent's and its neighbours' states whose weights sum to one.
  - At a contact (common p, T and u), the molar concentrations ρ_k/W_k sum to p/(R_u T).
  - E is the sum of (ρ_k/W_k) times each species' molar energy at T, plus ½ρu². At fixed T and u that is linear in the partial densities.
  - So the child keeps p, T and u, for any thermally perfect mixture, as the registration says. Range, positivity and content hold as stated, since a smaller φ keeps each component inside the range its own φ allowed.
- **One risk to check, not a conflict: a common limiter is as small as its worst component.** Barth-Jespersen's φ is 0 at a parent where any one component has an extremum.
  - There are 4 + n_s + n_t components. Radicals peak inside a flame, so more parents will fall back to piecewise constant than now. That is safe (it is what removed the jet), but it is first order there.
  - A component at rounding level is the sharper case. In criterion 3(a)'s axial pulse (`amr_subcycling_tests.cpp`, order 1.8 or better), ρu_r is 0 or rounding. If it is rounding, its extrema are everywhere and can set φ near 0 at many ghost parents.
  - AMR's ctest is on inquiry's list and will show this. If it bites, one remedy is to leave a component out of the common φ, under its own φ, when its neighbours' range is below a floor relative to its scale. Contact preservation then holds to that floor. That is inquiry's call.
- **The children's ρ.** With each partial density linear as well, a child's sum of partial densities equals its linear ρ to rounding only. Either one can be the one kept.
- **Conflicts with this lane's code: none.** This lane's uncommitted work touched only `tests/` and `CMakeLists.txt`. The ordering is the only overlap, and it is above.

## Inquiry's wall-face fix (762806e), read from the AMR side (23 October, about 01:00 CDT)

The change: `Flow::transportFluxes` skips `wallSolve` for a wall face whose cell is outside N_2. Under AMR such a cell has no computed state, so the solve read stale primitives. The review finds it sound (derived from the code):
- **The flux reader.** `lawFlux` (transport.cpp) reads a wall solution within the same evaluation, for a cell whose derivative is used. Those cells are active, inside N_2, and freshly solved.
- **The tag reader.** The wall-strip tag (amr.cpp, `lawYPlus`) fills its array from every wall face but reads it only on the level's active cells (`on(i, j)`). Activity changes only at a regrid, so each cell it reads was active through the step before and was solved in that step's last evaluation.
- **The accessor.** On a single level every cell is in N_2, so nothing changes there; the dry run and the probes are byte-identical, as 762806e records. Under AMR, `Flow::wallYPlus()` would return stale values for faces outside N_2. It has no AMR caller today (only `face_inlet_tests` and `chamber_study`, both single-level). A future AMR caller must filter by activity.
- **No conflict** with the leak fix of FACE_INLET.md section 5. That fix is in the end-face gather and does not touch the wall solve.

## Deviations from the plan's text

- **Active regions are unions of parent cells, not rectangles.** The plan says patches are rectangles of index space so that each is a small structured mesh for the existing kernels. Full-domain level arrays give the kernels their structured neighbours anyway, and interval lists make clustering into rectangles unnecessary. A rectangle cover would refine whole blocks of rows to hold a two-row wall strip. Nothing in the conservation argument depends on rectangles.
- **The restricted loops move from step A to step B** (decided 13 October, after the step A geometry). Which loops must also cover the ghost cells (refresh, slopes, profiles, gradients, wall solves) and which cover only the active cells (gather, sources, admissibility) is fixed by step B's ghost handling. Writing the interval loops before that would fix the wrong shape. They carry step A's check: a single level with everything active must give the same bits.
- **Step B in two parts** (decided 13 October, 23:37 CDT, before any step B code).
  - **B1.** Every level's kernels run over the whole level. Before each stage the driver fills *every* cell of level l that is not active from level l-1, not only the ghosts. The active cells then see exactly what width-2 ghosts would give them, plus filled values further out that no active stencil reads. B1 is inviscid only: the driver refuses transport and turbulence until the viscous and correction registers exist. Criteria 1 and 3(a), 3(b) are judged on B1.
  - **B2.** The restricted interval loops, with real ghost and off cells. Checked bit for bit against B1 on the active cells, then 3(c), which in B1 would hold trivially (every width fills everything). B2 comes before criterion 5(a), the cost.
- **Three new criteria (0, 5, 6).** They do not relax the plan's four (here 1 to 4); they record the refactor's safety, the cost the budget assumed, and the measurement the decision rule needs.

## Literature

- M. J. Berger and P. Colella, "Local adaptive mesh refinement for shock hydrodynamics", J. Comput. Phys. 82, 64 to 84 (1989): subcycling (p. 67), proper nesting (p. 68), refluxing (eqs. 1 to 3, pp. 69 to 70), shock tagging apart from the error estimate (pp. 75 to 76). Full text read on 5 October 2026 (TECHNICAL_PLAN).
- M. J. Berger and J. Oliger, "Adaptive mesh refinement for hyperbolic partial differential equations", J. Comput. Phys. 53, 484 to 512 (1984): known by summary only.
- Conservative linear interpolation to new fine cells with a limited gradient and a fallback to constant: standard practice in block-structured codes (AMReX's and Chombo's conservative interpolators), stated here from general knowledge, not re-read.
