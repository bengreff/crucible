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
