# Spray on AMR: evidence (docs/evidence/SPRAY.md, *Spray on AMR*)

- `a1_gas_probe.cpp`: A1's precondition, measured before the design was registered (22 October 2026, about 22:45 CDT). The gas alone (S3's chamber, no spray) on a single-level `Amr(def, 1)`, level 1 inactive, against a `Flow` on the AMR base's mesh (`Mesh::coarsened` of the finest). Result: 1,083 steps to 0.1 ms, every step's dt equal, and the final state and partial densities bit-identical (`yes`).
- `a0_dump.cpp`: A0, the uniform path unchanged. It writes the raw end state of S3's run (0.5 ms) and S8(c)'s (0.2 ms), with each parcel's fields written one by one.
  - Before the change, on the tree at b490437 (built 22 October, about 22:48 CDT): S3 took 5,415 steps and ended with 247 parcels; S8(c) took 2,168 steps and ended with 48 parcels. The file's md5 is `e95251ab17e43a961a392dcd4bd0912f`.
- Both are built against the tree's libraries:
  `c++ -std=c++20 -O2 -I. -DCRUCIBLE_SOURCE_DIR='"'$PWD'"' FILE.cpp build/libcrucible_thermo.a build/libcrucible_core.a $CANTERA_ROOT/build/lib/libcantera.a -framework Accelerate`.
