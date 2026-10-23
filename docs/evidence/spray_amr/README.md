# Spray on AMR: evidence (docs/evidence/SPRAY.md, *Spray on AMR*)

- `a1_gas_probe.cpp`: A1's precondition, measured before the design was registered (22 October 2026, about 22:45 CDT). The gas alone (S3's chamber, no spray) on a single-level `Amr(def, 1)`, level 1 inactive, against a `Flow` on the AMR base's mesh (`Mesh::coarsened` of the finest). Result: 1,083 steps to 0.1 ms, every step's dt equal, and the final state and partial densities bit-identical (`yes`).
- `a0_dump.cpp`: A0, the uniform path unchanged. It writes the raw end state of S3's run (0.5 ms) and S8(c)'s (0.2 ms), with each parcel's fields written one by one.
  - Before the change, on the tree at b490437 (built 22 October, about 22:48 CDT): S3 took 5,415 steps and ended with 247 parcels; S8(c) took 2,168 steps and ended with 48 parcels. The file's md5 is `e95251ab17e43a961a392dcd4bd0912f`.
- Both are built against the tree's libraries:
  `c++ -std=c++20 -O2 -I. -DCRUCIBLE_SOURCE_DIR='"'$PWD'"' FILE.cpp build/libcrucible_thermo.a build/libcrucible_core.a $CANTERA_ROOT/build/lib/libcantera.a -framework Accelerate`.

## Results (22 October 2026, about 22:55 to 23:15 CDT; SPRAY.md, *Spray on AMR: results*)

All of these were run on the Mac, on the change's tree. Every value is measured.

- `a0_spray_ctest_after_2026-10-22.log`: `crucible_spray_tests` after the change. Its one failure is S1(b), which predates the change. A0's dump after the change is `cmp`-identical to the one before (md5 `e95251ab17e43a961a392dcd4bd0912f`).
- `a1_2026-10-22.log`: `crucible_spray_amr_tests a1`.
  - A1(a) is met bit for bit.
  - A1(b) is not met as written: its three FAIL lines.
- `a2_a3_2026-10-22.log`: `crucible_spray_amr_tests a2`. A2 is met on 1 thread and on 4, and A3 is met.
- `common.inc`: the test's helpers (chamber, settings, oxygen, the byte and state comparisons), shared by the probes below.
- `a1_reference_probe.cpp`: S3, frozen. It runs `ReactingFlow` (modes 0 and 1), `Flow::step` with `Spray::advance` (mode 2), and mode 2 with an extra `stableDt()` (mode 3), each against the single-level `Amr`. It prints the first step that differs.
- `a1b_sst_probe.cpp`: S8(c), SST. It runs the gas alone through `Flow::step` (mode 0), the same with a `stableDt()` before each step (mode 1), and mode 1 with the spray (mode 2), each against the single-level `Amr`.
- `a1_probes_2026-10-22.log`: the probes' output. In both probes the steps are numbered from 0, so "step 1" means the second step.
- `a2_uniform_probe.cpp` and `a2_uniform_2026-10-22.log`: S3's chamber with the spray on uniform 16 × 6 and 64 × 24 meshes, printed beside A2. Not judged.
- `a2_levels_probe.cpp` and `a2_levels_2026-10-22.log`: the defect's localisation (SPRAY.md, *A defect the printed comparison found*).
  - The modes: 0, every level active; 1, the static corner tag; 3 and 4, modes 1 and 0 with no injectors; 5, one-way; 6, no point-implicit loads.
  - Modes 0, 3 and 4 ran before the 500-step time trace was added, and modes 1, 5 and 6 with it.
  - Each mode prints the end field's u_z along finest rings 0 and 4 (the level that owns each cell, then u_z), and the parcels per 4 finest columns.
