# Session handoff

4 October 2026 · branch `main` (one line; `claude/verify-core` was fast-forwarded into it and is kept only as a pointer).

## Where things stand

- **One project.** Post-pivot work is `main`. The docs are VISION_SCOPE.md (what), TECHNICAL_PLAN.md (how), RESEARCH.md (evidence, open model decisions, pre-pivot lessons in section 7) and this file (state). Verification records are in `docs/evidence/` ([index](evidence/README.md)), real-engine cases in `docs/validation/`, and the MHD spike in `docs/parked/`. `archive/` is unchanged and historical.
- **Ben, 4 October 2026** (recorded in VISION_SCOPE *Changes since the pivot*): shifting equilibrium is the validated baseline, produced by the general reacting machinery. The chamber always evolves in time; a chemical run starts "valves open, ignite" (ambient chamber, valves open on a schedule, igniter fires, transient to steady). "Add evaporation first": liquid LOX with breakup and evaporation is built before RL10; no pre-vaporized RL10. By January all four regimes work (chemical, magnetic nozzle, ICF post-burn pulse, antimatter); the MHD spike un-parks right after chemical; the app draft ("Draw and run", "Probe and plot") follows chemical. Dated order: TECHNICAL_PLAN *Dated plan to January*.
- **Engine today** (all verified, see the evidence index): axisymmetric finite-volume gas core (HLLC, SSPRK2, device-thrust ledger), thermally perfect multi-species mixture, Cantera thermochemistry matching CEA (worst 0.131% over 18 points), CVODES kinetics matching Cantera, Strang coupling verified on a detonation (front speed within 0.24%). The native app runs the gas nozzle live. Release is the default build.
- **C1 done (measured, `docs/evidence/CHAMBER_C1.md`).** Chamber case: contour as data, supply faces with valve ramps, igniter, ambient N2 start, chemistry modes (frozen, finite rate, local equilibrium). Local equilibrium against the ideal rocket with stated 2-D corrections: c* +1.66% (32x6) and +0.40% (64x12); vacuum Isp -0.26% and -0.11%. Finite-rate light-off with the 300 J igniter lit at about 0.22 ms; settled Isp 424.31 s lies between frozen 412.77 and shifting 429.09. The energy budget fails its 1e-11 criterion only because it is divided by the fill energy (about -0.25 J); recorded, criterion unchanged. The third grid (128x24, to 5 ms) was restarted 4 October about 18:57 in `/tmp/c1/eq128/` from `build/alt/crucible_chamber_study` (do not rebuild `build/alt` while it runs; at 0.5 ms after 2194 s wall, so about 6 h); when it ends, plot with `tools/chamber_plots.py`, add the row and verdict 5 to CHAMBER_C1.md, copy the log and PNGs to `docs/evidence/c1/`.
- **C2 step 1 (measured, `docs/evidence/TRANSPORT_C2.md`).** Mixture-averaged properties match Cantera to 8.8e-14. The transport operator (stress, conduction, diffusion, walls) passes pipe decay (order 1.99) and the budgets, and energy and species truncation order (1.93). Axial momentum order 1.799 and radial momentum order 1.19 fail the 1.8 criterion; two fixes were tried, the mechanism is written down, and the untried candidate is the r^2 difference on every radial face.
- **Subsonic outlet option (measured, `docs/evidence/FLAME_C2.md`).** `Definition::outletRelaxation` (sigma) makes the outlet partially non-reflecting (Poinsot and Lele), in the characteristics W = p +- rho a u_z; zero keeps the fixed-pressure node bit for bit. Core checks against the theory impulse response and an entropy-wave check pass.
- **C2 step 2, flame speed (in progress, `docs/evidence/FLAME_C2.md`).** The fixed-pressure outlet made the duct resonate (40 um: S_c +20.4%), so the flame case now uses sigma 0.25 and a burned-gas ambient (header of `tests/flame_study.cpp`, criteria unchanged). With that outlet the 40 um flame settles near +1.4%. Runs in `/tmp/c2flame/nr/` (f40, f20; f10 next, about 16 times the 40 um cost); plot with `tools/flame_plots.py`.
- **ctest on `main`:** 5 of 6 pass; `transport_verification` fails on those two order criteria (`docs/evidence/ctest_main_2026-10-04_c2.txt`).

## Next, in order (TECHNICAL_PLAN, *Dated plan to January*)

1. **C1:** done except the 128x24 grid row (above).
2. **C2, 12-18 Oct:** finish the flame-speed grids (20 and 10 um) and judge the four criteria in FLAME_C2.md; run the full ctest. Then SST-2003 URANS, PaSR, and interleaved fuel and oxidizer rings.
3. **C3, 19-25 Oct:** liquid LOX as Lagrangian parcels with breakup and evaporation (TECHNICAL_PLAN step 5; evidence questions in RESEARCH *Liquid injection*).
4. **C4, 26 Oct-1 Nov: RL10A-3-3A.** Geometry from NASA TM-107318 Appendix E Table E1 (area against axial station, -12 in to +41.84 in; r = r_t sqrt(A/A*)). Operating inputs come from the same or primary sources. Pre-register the vacuum Isp and thrust in git before reading measured values. Known exposure, to declare in the pre-registration: the widely quoted nominal values for this engine are general knowledge. Use a Sonnet agent to extract inputs only, with an explicit instruction never to report measured performance.
5. App draft (2-8 Nov), then magnetic nozzle, ICF post-burn pulse, antimatter (dates in TECHNICAL_PLAN).

## Open items carried forward

- **Low-Mach accuracy.** HLLC dissipation scales with sound speed. At Mach about 0.15, total-pressure error is amplified about 33x in mass flow. Subsonic venturi mass flow: -9.49% (40x6) to -2.63% (80x12), order 1.85. The default stays plain HLLC; chamber accuracy comes from grid convergence with reported error bands (`docs/evidence/LOW_MACH.md`).
- **Reacting step cost.** Almost every reacting step re-plans once (the first half reaction lowers the CFL limit), about a third more reaction work. A sound-speed predictor would remove most of it.
- **Static equilibrium maximum is first order** at the wall row and at smooth axial extrema; L1 is second order.
- **Still missing:** save/load and control histories, editable geometry, retained comparison runs, probe picking, snapshot decimation, cross-platform builds and packaging.
- **Parked MHD findings** (for when it resumes): `docs/parked/MHD_SPIKE.md`. Its short version: keep the body-fitted RZ mesh with constrained transport through nodal psi; resistive MHD first, Hall next; a characteristic outlet is needed for plume domains.

## Continue from here

- Numerical core: `core/flow.*`, `core/medium.*`. Reacting coupling: `adapters/reacting_flow.*`. Chemistry: `adapters/reaction.*`. Thermo and the ideal rocket: `adapters/thermo.*`.
- Tests: `tests/`. Headless builds: `build/` (core plus Cantera, no app) and `build/native` (app).
- Cantera 3.2.0 is built from source in `~/src/cantera` (static library).
- Tested on this Mac: Qt 6.11.2, VTK 9.7.0, AppleClang 17.
- Heavy jobs go through `python3 ~/director/harness/slot.py run --label ... -- <command>`.
