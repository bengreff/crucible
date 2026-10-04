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

The ctest suite on `main` (4 October 2026): 5/5 pass in 141 s, Release build, this Mac (`ctest_main_2026-10-04.txt`).
