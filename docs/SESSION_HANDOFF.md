# Session handoff

**MHD spike PARKED until after the chemical milestone, per Ben 2026-09-30.** Active order: Cantera thermochemistry vs CEA, validation-data survey, low-Mach accuracy, reacting chamber.

30 September 2026 · branch `claude/verify-core` (pushed to origin; not merged to main, which is pre-pivot; that merge is a separate decision).

The native converging-diverging nozzle experiment exists: verified axisymmetric gas flow, live inlet pressure, measurements. The app starts paused; README has build/run instructions. `docs/IMPLEMENTATION.md` holds the evidence and scientific limitations. The archive is unchanged.

## Chemical milestone (active, 2026-09-30)

- **Item 1 done.** Cantera 3.2.0 is built from source (`~/src/cantera`, `scons build`, static lib) and wrapped in `adapters/thermo.*`. The ideal rocket matches NASA CEA (RocketCEA) at 18 points (H2/O2, CH4/O2, LOX/LH2 at 3 O/F each, equilibrium and frozen): worst 0.131% against the 0.5% tolerance. Evidence: `docs/THERMO_VERIFICATION.md`, ctest `thermo_cea_verification` (about 40 s).
- **Item 2 done.** `docs/VALIDATION_SURVEY.md` recommends TUM GOX/GCH4 SFB/TRR40 Test Case 1 (20 bar, O/F 2.6; wall pressure, heat flux, combustion efficiency) and Penn State GO2/GH2 (wall heat flux). Next input: the Test Case 1 geometry/BC description and digitised curves.
- **Next:** item 3 (low-Mach correction or preconditioning for HLLC; venturi -9.5% mass flow at Mach 0.15 coarse), then item 4 (reacting chamber).

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
