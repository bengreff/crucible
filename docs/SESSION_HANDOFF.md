# Session handoff

30 September 2026 · branch `claude/verify-core` (pushed to origin; not merged to main, which is pre-pivot; that merge is a separate decision).

The native converging-diverging nozzle experiment exists: verified axisymmetric gas flow, live inlet pressure, measurements. The app starts paused; README has build/run instructions. `docs/IMPLEMENTATION.md` holds the evidence and scientific limitations. The archive is unchanged.

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
  - Conductivity: Director default (reversible by Ben) is resistive MHD first, Hall term next; recorded in `RESEARCH.md`. An outside check on the model class is pending.
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
