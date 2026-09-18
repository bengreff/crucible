# First implementation session handoff

17 September 2026 · branch `codex/first-native-nozzle`.

Ben selected a native converging-diverging nozzle experiment: watch verified gas flow, change inlet pressure live, inspect measurements. This increment is implemented. The app starts paused; README has build/run instructions. Use `docs/IMPLEMENTATION.md` for exact evidence and scientific limitations. The archive remains unchanged (all 240 manifest entries verified).

## Continue from here

- Numerical core: `core/flow.*`; independent references and worker tests: `tests/core_tests.cpp`.
- Worker/control/snapshot ownership: `core/session.*`.
- Native UI and headless runner: `app/`; CMake builds both from the same core.
- Tested dependencies on this Mac: Qt 6.11.2, VTK 9.7.0, AppleClang 17. Homebrew installed Qt/VTK and their dependencies. Build directories are ignored, not portable artifacts.
- Verified release suite, ASan/UBSan suite, and real native rendering/control/pause smoke test. Screenshot is local at `build/native/smoke.png`. Cross-platform builds and distribution packaging remain untested.

Next recommended increment: review the numerics independently, establish full axial momentum/force accounting, and add saved experiment definitions and measurement/control histories. Then select a well-documented real propulsion dataset and integrate thermodynamics toward the reacting-chamber milestone. Do not label this ideal-gas preset a validated chemical engine or turn the project into a nozzle-only tool.

Important gaps: no editable geometry, chemical reactions, turbulence, plasma, save/restart, retained comparison runs, probe picking or scalable snapshot decimation yet. Prepared flowing initialization does not simulate startup. Outlet-plane force is not total thrust. Mesh/boundary generalization and AMReX selection remain open. The million-cell ten-step benchmark demonstrates memory/execution feasibility, not useful physical run duration in minutes.

The native shell is intentionally small and should be factored into presentation, controls, export and session interaction as those behaviors grow. Do not build a general framework before these behaviors need it. Numerical reference tests must survive any reorganization.
