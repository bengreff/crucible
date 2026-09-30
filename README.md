# CRUCIBLE

A human-facing technical sandbox for investigating how reacting gas and plasma, geometry, and fields interact to produce propulsion.

**Status: first native gas-flow prototype, 17 September 2026.** Operate a converging-diverging nozzle, change reservoir pressure live, inspect fields and measurements, and export the current field. Ideal-gas numerical checks pass; combustion and real-engine validation are still ahead.

- [Vision and scope](VISION_SCOPE.md) is the authoritative definition of the project, including the visual human workflow.
- [Research](RESEARCH.md) records sources, transport/combustion model candidates, feasibility experiments, and scientific evidence requirements.
- [Technical plan](TECHNICAL_PLAN.md) defines the native app, dependencies, state ownership, candidate algorithms, live interaction, and implementation sequence.
- [Rebuild notes](REBUILD_NOTES.md) preserve lessons, useful decisions, and specific archived code references.
- [Archive](archive/README.md) explains what was preserved, how to verify it, and what is local-only.

The intended workflow is **construct → operate → observe → change → compare**, inside one native desktop app running locally on Mac, Windows, or Linux. The combustion chamber's reacting contents belong inside the model; equipment is represented through physical interfaces. Visualization and usability develop alongside the physical model.

One shared simulation engine gains capabilities for all compatible experiments. Chemical propulsion validated against real measurements is the first milestone, followed by magnetic, pulse and energetic-product physics in that same engine. The desired development horizon is a few months; implementation and validation remain to be done.

## Build and run

Requires CMake 3.24+, a C++20 compiler, Qt 6.5+ Widgets/OpenGLWidgets, and VTK 9.3+ built with Qt support. Tested here with AppleClang 17, Qt 6.11.2 and VTK 9.7.0. On this Mac the dependencies were installed with `brew install qtbase vtk`. Other platforms need equivalent packages and may need `CMAKE_PREFIX_PATH` pointing to their installation. Packaged distribution is not ready.

```sh
cmake -S . -B build/native -DCMAKE_BUILD_TYPE=Release
cmake --build build/native --config Release
ctest --test-dir build/native -C Release --output-on-failure
open build/native/crucible.app
```

On Linux run `build/native/crucible`; on Windows run the generated `crucible.exe` (usually under `Release/`). These platforms have not been tested yet.

For the dependency-light headless core:

```sh
cmake -S . -B build/core -DCRUCIBLE_BUILD_APP=OFF -DCMAKE_BUILD_TYPE=Release
cmake --build build/core --config Release
ctest --test-dir build/core -C Release --output-on-failure
build/core/crucible_run --nz 160 --nr 24 --time 0.004
```

The app starts paused. Press **Run**, change the inlet reservoir pressure and press **Apply live**, then **Pause** to inspect. **Restart** restores the default initial state and pressure. Field choices are Mach, pressure, temperature and density. The metrics show device thrust (supply plane + wall forces - ambient), its exit-plane estimate, and mass/energy/momentum balance errors. CSV saves a field snapshot, not a restart.

A native integration check can run with `crucible --smoke-test /absolute/path/smoke.png` using the executable inside the Mac bundle. It needs a graphical desktop and exits after running, applying a control and pausing. Sanitized core builds use `-DCRUCIBLE_SANITIZERS=ON` with Clang/GCC.

Read [implementation, evidence and limitations](docs/IMPLEMENTATION.md) before interpreting results or continuing development. The inherited archive gate battery is not an active build requirement.
