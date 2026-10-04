# Session handoff

4 October 2026 · branch `main` (one line; `claude/verify-core` was fast-forwarded into it and is kept only as a pointer).

## Where things stand

- **One project.** Post-pivot work is `main`. The docs are VISION_SCOPE.md (what), TECHNICAL_PLAN.md (how), RESEARCH.md (evidence, open model decisions, pre-pivot lessons in section 7) and this file (state). Verification records are in `docs/evidence/` ([index](evidence/README.md)), real-engine cases in `docs/validation/`, and the MHD spike in `docs/parked/`. `archive/` is unchanged and historical.
- **Ben, 4 October 2026** (recorded in VISION_SCOPE *Changes since the pivot*): shifting equilibrium is the validated baseline, produced by the general reacting machinery. The chamber always evolves in time; a chemical run starts "valves open, ignite" (ambient chamber, valves open on a schedule, igniter fires, transient to steady). "Add evaporation first": liquid LOX with breakup and evaporation is built before RL10; no pre-vaporized RL10. By January all four regimes work (chemical, magnetic nozzle, ICF post-burn pulse, antimatter); the MHD spike un-parks right after chemical; the app draft ("Draw and run", "Probe and plot") follows chemical. Dated order: TECHNICAL_PLAN *Dated plan to January*.
- **Engine today** (all verified, see the evidence index): axisymmetric finite-volume gas core (HLLC, SSPRK2, device-thrust ledger), thermally perfect multi-species mixture, Cantera thermochemistry matching CEA (worst 0.131% over 18 points), CVODES kinetics matching Cantera, Strang coupling verified on a detonation (front speed within 0.24%). The native app runs the gas nozzle live. Release is the default build.
- **ctest on `main`:** 5/5 pass (`docs/evidence/ctest_main_2026-10-04.txt`).

## Next, in order (TECHNICAL_PLAN, *Dated plan to January*)

1. **C1, 5-11 Oct: time-evolving chamber and nozzle.** Wall contour from data, gas supply faces (imposed mass flux, total enthalpy and composition per radial face) with valve ramps, an igniter energy deposit, and a local-equilibrium chemistry option as the verification limit. A premixed gaseous H2/O2 "valves open, ignite" run settles, and its end state is checked against CEA.
2. **C2, 12-18 Oct:** molecular transport, SST-2003 URANS, PaSR, no-slip and thermal walls.
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
