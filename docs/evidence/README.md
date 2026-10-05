# Verification evidence

Each record says what was checked, against which independent reference, the measured result and the limits. Raw outputs sit beside the records. Verification shows the code solves its equations; it is not validation against a real engine (that lives in `docs/validation/`).

| Record | What it verifies | Reference | Check |
|---|---|---|---|
| [IMPLEMENTATION.md](IMPLEMENTATION.md) | Axisymmetric gas core: HLLC, axis treatment, nozzle limits, device-thrust ledger, live control | Exact Riemann, quasi-1-D nozzle theory, grid study | ctest `core_verification`, `core_slow_verification` |
| [THERMO_VERIFICATION.md](THERMO_VERIFICATION.md) | Cantera thermochemistry and the ideal rocket | NASA CEA (RocketCEA), 18 points, worst 0.131% | ctest `thermo_cea_verification` |
| [REACTION_VERIFICATION.md](REACTION_VERIFICATION.md) | Stiff reaction integration (CVODES BDF over Cantera rates) | Cantera ReactorNet, UV equilibrium | ctest `reaction_verification` |
| [MIXTURE_CORE.md](MIXTURE_CORE.md) | Thermally perfect multi-species core, species fluxes | Cantera thermo, exact two-gamma shock tube | ctest `mixture_verification` |
| [REACTING_FLOW.md](REACTING_FLOW.md) | Strang coupling of chemistry and flow | Equilibrium Hugoniot (piston-supported H2/O2/Ar detonation): front speed within 0.24% | study `crucible_detonation_study` |
| [LOW_MACH.md](LOW_MACH.md) | Low-Mach accuracy of HLLC; Thornber and HLLC-LM variants | Subsonic venturi, transverse stability | studies `crucible_low_mach_study`, `crucible_low_mach_stability` |
| [CHAMBER_C1.md](CHAMBER_C1.md) | Time-evolving chamber and nozzle: contour as data, supply faces with valve ramps, igniter, chemistry modes; "valves open, ignite" premixed H2/O2 | One-dimensional ideal rocket (CEA problem) at the simulated chamber pressure, with the Kliegel-Levine discharge coefficient and the conical divergence factor | study `crucible_chamber_study`; fast checks in `core_verification` |
| [TRANSPORT_C2.md](TRANSPORT_C2.md) | Mixture-averaged transport properties and the transport operator (stress, conduction, diffusion, walls) | Cantera MixTransport; exact cell averages of the continuum operator; Bessel pipe decay | ctest `reaction_verification`, `transport_verification` (2 order criteria fail) |
| [FLAME_C2.md](FLAME_C2.md) | Freely propagating premixed flame through the full engine (flow, transport, CVODES chemistry); the partially non-reflecting subsonic outlet | Cantera Flow1D free flame (S_L 2.330470 m/s); Poinsot-Lele reflection theory | study `crucible_flame_study`; outlet checks in `core_verification` (in progress: 40 um S_c +1.34%) |
| [SST_REFERENCE.md](SST_REFERENCE.md) | Independent 1-D SST-2003 solver for fully developed compressible pipe flow, the target of the engine's periodic-pipe check | Grid convergence to 3200 rings, Richardson limits, wall-omega factor sensitivity | `tools/sst_pipe_1d.py` (first order from the Menter wall omega; limit c_f 0.0104865) |
| [turbulence_verification_2026-10-04.txt](turbulence_verification_2026-10-04.txt) | SST machinery in the engine, step 1: exact wall distance (`core/walls.cpp`) | Straight duct d = R - r; ternary search per segment on the nozzle with plate walls | ctest `turbulence_verification` (all pass, worst 4e-16 of R) |

The ctest suite on `main` (4 October 2026, after C2 step 1): 5 of 6 pass in 163 s, Release build, this Mac (`ctest_main_2026-10-04_c2.txt`). `transport_verification` fails on two truncation-order criteria (axial and radial momentum), as recorded in TRANSPORT_C2.md. The earlier 5/5 record before transport is `ctest_main_2026-10-04.txt`.
