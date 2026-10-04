# Session handoff

**MHD spike PARKED until after the chemical milestone, per Ben 2026-09-30.** Active order: Cantera thermochemistry vs CEA, validation-data survey, low-Mach accuracy, reacting chamber.

4 October 2026 · branch `claude/verify-core` (pushed to origin; not merged to main, which is pre-pivot; that merge is a separate decision).

**Stopped for the night on 4 October at the Director's request (weekly usage budget).** The last item finished was the archive fold below. Nothing is half-done. Item 4 step 1 (molecular transport) has **not** been started; pick up there.

The native converging-diverging nozzle experiment exists: verified axisymmetric gas flow, live inlet pressure, measurements. The app starts paused; README has build/run instructions. `docs/IMPLEMENTATION.md` holds the evidence and scientific limitations.

**Archive fold (4 October).** The docs from the old separate repo `bengreff/inquiry_project` (July 2026) now live in `archive/inquiry_project-2026-07/`, so CRUCIBLE is one project. They are historical only. `archive/README.md` lists both archives. `archive/pre-pivot-2026-09-17/` is unchanged.

## Chemical milestone (active, 2026-09-30)

- **Item 1 done.** Cantera 3.2.0 is built from source (`~/src/cantera`, `scons build`, static lib) and wrapped in `adapters/thermo.*`. The ideal rocket matches NASA CEA (RocketCEA) at 18 points (H2/O2, CH4/O2, LOX/LH2 at 3 O/F each, equilibrium and frozen): worst 0.131% against the 0.5% tolerance. Evidence: `docs/THERMO_VERIFICATION.md`, ctest `thermo_cea_verification` (about 40 s).
- **Item 2 done.** `docs/VALIDATION_SURVEY.md` recommends TUM GOX/GCH4 SFB/TRR40 Test Case 1 (20 bar, O/F 2.6; wall pressure, heat flux, combustion efficiency) and Penn State GO2/GH2 (wall heat flux). Next input: the Test Case 1 geometry/BC description and digitised curves.
- **Ben's decisions (via Director):** primary case is the TUM round chamber, with exact published geometry (stop if an input is missing); strict pre-registered blind prediction; second case is RL10A-3-3A. Inputs and gaps: `docs/VALIDATION_TUM_ROUND.md`. Missing: convergent nozzle contour, the TRR40 Test Case 1 document, nozzle-wall thermal condition, propellant purity. Blind exposure is declared: efficiency and Pc at O/F 2.2 were seen.
- **Item 3 done.** `docs/LOW_MACH.md`. Thornber is 2-6x more accurate on the venturi but grows a transverse odd-even mode exponentially (row deviation 1e-14 to 7e-5 at 800 cells). HLLC-LM has no effect above M 0.1. The default stays plain HLLC; chamber accuracy comes from grid convergence with reported error bands. Both variants remain behind `Definition::lowMach`.
- **Item 4 in progress (reacting chamber).** Done and pushed:
  - **Stiff reaction integration** (`adapters/reaction.*`, `docs/REACTION_VERIFICATION.md`, ctest `reaction_verification`): CVODES BDF over Cantera rates, matching Cantera ReactorNet (ignition delay to 2.6e-6) and UV equilibrium.
  - **Multi-species core** (`core/medium.*`, `docs/MIXTURE_CORE.md`, ctest `mixture_verification`): thermally perfect NASA7 mixture with Larrouturou species fluxes. Matches Cantera thermo to 6e-15 and an exact two-gamma shock tube.
  - **Release is now the default build.** All earlier timings were -O0. The single-gas suite takes 12 s against 5.6 s for the old core.
  - **Strang coupling** (`adapters/reacting_flow.*`, `docs/REACTING_FLOW.md`, study `crucible_detonation_study`): exactly symmetric, re-planning the step when heat release lowers the CFL limit.
    - Verified on a piston-supported H2/O2/Ar detonation against the equilibrium Hugoniot.
    - Front speed is within 0.24% on every grid (1% tolerance). Grid and dt contributions are below 0.03 and 0.001 percentage points.
    - The remainder is a start-up transient that decays with distance (+0.05% at 2 m).
    - The late burned plateau matches exact equilibrium to 0.01%.
  - `ctest` 5/5 pass (159 s).
- **Next in item 4, in order:**
  1. Mixture-averaged molecular transport: export Cantera's fits to the core, verify against Cantera.
  2. 1-D premixed laminar flame speed against Cantera FreeFlame (the first deflagration check, needs 1).
  3. SST URANS + PaSR.
  4. Injector mass-flow inlets.
  5. Chamber contour as data (exact TUM shape drops in; meanwhile the 30/45 deg, sharp/rounded bracket).
  6. Pre-registered TUM predictions at O/F 2.6/3.0/3.4, committed before digitising the measured curves.
  7. Grid convergence on backhouse.
- **Efficiency note:** almost every reacting step re-plans once, about a third more reaction work. A sound-speed predictor would remove most of it.

## Done on this branch

- **Verification review.**
  - Second-order boundary reconstruction.
  - Independent numerics tests (exact Riemann, HLLC properties, rarefaction, strong shock, internal normal shock, venturi, acoustic mode).
  - Nozzle grid study (`docs/evidence/`).
- **Device-thrust ledger.**
  - `deviceThrust = inlet momentum flux + wall force + bodyAxialForce − ambient` = exit-plane thrust + dP_z/dt.
  - Measured 517.88 N at 320×48 against quasi-1D ideal (C_F/ideal − 1 = +1.4e-4).
  - `bodyAxialForce` is the empty slot for a Lorentz force whose reaction acts on coils. The UI and README report thrust by component.
- **MHD architecture spike:** `docs/MHD_SPIKE.md`, throwaway code in `spike/`.
  - Recommendation: keep the body-fitted RZ mesh; use CT through nodal ψ.
  - Before any field work: fix axis reconstruction (done), then the split B0 + B1 form feeding `bodyAxialForce`, explicit magnetic boundary conditions, and a conductivity model.
  - Split B0 + B1 (`split` mode in the spike, memo case 5): static coil exact, coil reaction in `bodyAxial` (ledger closes). It does **not** fix low-β positivity. Measured mechanism: ideal flux freezing carries the throat's coil flux (about 1.2 mWb at 1 T) downstream, and its magnetic pressure evacuates the wall region. Needs resistivity and wall magnetic conditions, not numerics.
  - Conductivity: Director default (reversible by Ben) is resistive MHD first, Hall term next; recorded in `RESEARCH.md`. Outside review: `docs/ASTRA_PLASMA_MODEL_2026-09-30.md` (opinion, not canon).
  - Wall conditions and positivity (memo case 6): transparent, insulating (vacuum-matched exterior) and conducting walls; dual energy plus first-order HLL retry. The 1 T crash was pressure recovery, not the wall. With robustness, thrust depends on the wall condition at O(1) and is not mesh-converged. The ideal insulating wall was stopped under the two-fix rule (an unresolved current sheet at the wall needs resistivity). Conducting at 160 cells blows up at the outlet (backflow through an extrapolation outlet). That needs a characteristic outlet.
  - Resistive term (memo case 7): verified at second order (sinusoid, axial, Bessel with axis); Joule closure; Spitzer checked against Formulary hand values. Still first order: energy through open zero-gradient end planes.
  - Design notes, no code: `docs/VALIDITY_MONITOR.md` and `docs/VALIDATION_MAENO2013.md`. Maeno decided (Director, reversible by Ben): a prescribed plume from independent ablation data with propagated ranges; recorded in `RESEARCH.md`.
  - Resistive 1 T retry (memo case 8): with eta_m = 1 or 0.1 m^2/s (assumed scan), insulating and conducting walls both run cleanly (zero retries) and agree within 1 to 8%. They are mesh-stable from 80 to 160 cells. The transparent wall lets plasma currents cancel the coil field, so it is unphysical. The outlet did not block. A characteristic outlet is still needed for plume domains.
  - Next: the validity-monitor maps (needs a plasma EOS/ionisation state to be meaningful), then the Maeno plume inputs from published ablation data.
- **Axis reconstruction fixed in core** (centroid-referenced radial slopes, parity on the axis row, exact p/r source). Evidence in `docs/IMPLEMENTATION.md`: axis-row acoustic error first order → about third order; coil-force rest residual converges; nozzle grid study and venturi unchanged or slightly better.

## Open items

- **Low-Mach accuracy.**
  - HLLC dissipation scales with sound speed, so at Mach ~0.15 mass flow amplifies total-pressure error by about 1/(γM²) ≈ 33x.
  - Measured on the subsonic venturi: mass flow −9.49% (40×6) → −2.63% (80×12), converging at order 1.85.
  - Combustion chambers run at Mach 0.1 to 0.3, so decide on a low-Mach correction (preconditioned or all-speed flux) with chamber evidence before trusting chamber observables.
- **Save/load and control histories** were deferred for the spike (former item 4). They are still needed.
- Still missing: editable geometry, chemistry, turbulence, plasma, retained comparison runs, probe picking, and snapshot decimation.
- **Static equilibrium maximum is first order** at the wall row (one-sided limited slope) and at smooth axial extrema (minmod clipping). L1 is second order. Matters for magnetic nozzles, whose force peaks at the wall.
- Prepared flowing initialization does not simulate startup.
- Cross-platform builds and packaging are untested.

## Continue from here

- Numerical core: `core/flow.*`; tests: `tests/core_tests.cpp`, `tests/nozzle_convergence.cpp`.
- Worker/control/snapshot ownership: `core/session.*`.
- Native UI and headless runner: `app/`.
- Tested on this Mac: Qt 6.11.2, VTK 9.7.0, AppleClang 17.
- Do not build a general framework before behaviours need it. Numerical reference tests must survive any reorganization.
