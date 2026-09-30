# MHD architecture spike (September 2026)

Question: can the current body-fitted RZ mesh and data layout carry magnetic fields (div B control, an MHD Riemann solver, axis treatment) before months are committed to them?

Answer: yes for the mesh and layout. The field transport machinery worked. What broke is (a) axis reconstruction, which is an RZ issue on any mesh and is also latent in the gas core, and (b) the physics model and magnetic boundary conditions. Neither is fixed by changing the mesh. Recommendation at the end.

Throwaway code: `spike/mhd_rz.cpp` (target `crucible_mhd_spike`), plot script `spike/plot_coil.py`. The spike is not wired into the app or the core and is not a framework. Raw tables and figures are in `docs/evidence/mhd_spike/`. Every number below is **measured** unless labelled otherwise.

## What was built

- Ideal MHD with 9 variables (ρ, ρu_z, ρu_r, ρu_θ, E, B_z, B_r, B_θ, ψ_GLM) in b = B/√μ0 units. It uses the core's mesh geometry: exact frustum volumes, integrated face vectors and the meridional measure for the p/r source.
- HLLD (Miyoshi-Kusano) in a rotated face frame (n, t1, θ). The same code serves axial faces and sloped radial faces.
- Axisymmetric geometric sources for all momentum and field components. The axisymmetric solenoidal source (B_θ and B_r terms) is included.
- Three div B options:
  - none;
  - GLM (Dedner) with parabolic damping;
  - CT through the flux function ψ = r A_θ stored at mesh nodes. dψ/dt = r (u×B)_θ, and the node EMF is the average of the adjacent Riemann fluxes. Cell B is the exact volume average of the ψ field.
  - CT needs one extra scalar per node and nothing else.
- Applied field: a single current loop (radius 5 cm at the throat plane z = 0.216 m) from AGM elliptic integrals. The field is specified by B at the coil centre.
- Same gas boundaries as `core/flow`. The wall and inlet copy B from the interior (discussed under "What broke").

## Verification cases

### 1. Brio-Wu shock tube in RZ rows (`briowu.csv`, `briowu_800.png`)

- **Axial tube:** B_θ is transverse, run in a full cylinder with rows.
- **Radial tube:** propagates in r inside a thin annulus at r ≈ 1e6 m, with B_r = 0.75 normal and B_z transverse.
- **Reference:** the same solver at 12800 cells. The compound wave has no closed-form solution, so this is a self-convergence reference, not an exact one.

| cells | L1 ρ axial | L1 ρ radial | L1 u_n | L1 B_t | HLL fallbacks |
|---|---|---|---|---|---|
| 200 | 7.825e-3 | 7.825e-3 | 1.594e-2 | 1.055e-2 | 0 |
| 400 | 4.655e-3 | 4.655e-3 | 8.763e-3 | 6.007e-3 | 0 |
| 800 | 2.658e-3 | 2.658e-3 | 4.612e-3 | 3.255e-3 | 0 |
| 1600 | 1.437e-3 | 1.437e-3 | 2.325e-3 | 1.668e-3 | 0 |

- Observed order in ρ: 0.75, 0.81, 0.89. That is the expected sub-first-order rate for a solution with discontinuities.
- Axial and radial agree to 4 digits, and mass/energy residuals are ≤ 6e-14.
- The profiles in the figure show the compound wave, the slow shock and the rarefactions in place.
- Test-design note: an earlier radial case with u_r transverse was invalid because the walls block it. It was replaced; the solver was not the problem.

### 2. Field-aligned flow (`aligned.csv`)

- Sod in a full cylinder (100/200/400 × 8) with uniform B_z along the flow. The exact solution is the gas Riemann solution, and every transverse quantity must stay zero.
- **Structure:**
  - Transverse u_r, B_r and B_θ stay at roundoff (≤ 7e-13 without cleaning, exactly 0 for B_r with CT).
  - Conservation errors are ≤ 5e-14, with no HLL fallbacks.
  - none and CT agree to all printed digits.
- **Accuracy degrades as β falls:**

| B_z (code units) | β right state | L1 ρ at 100 / 200 / 400 |
|---|---|---|
| 0 | ∞ | 6.87e-3 / 4.01e-3 / 2.23e-3 |
| 4.47 | 0.01 | 1.22e-2 / 6.61e-3 / 3.53e-3 |
| 14.14 | 0.001 | 1.96e-2 / 1.07e-2 / 5.68e-3 |

  - At β = 0.001 the error is 2.5x the gas-only error at equal cells. HLLD bounds its waves by the fast speed, so pure acoustic structures get Alfvén-speed dissipation, and the timestep shrinks with the fast speed.
  - This matters for low-β magnetic nozzles (inferred).

### 3. Magnetostatic equilibrium: gas at rest in a 1 T coil field (`static.csv`, `residual.csv`)

- The exact answer is that nothing moves, because a vacuum coil field has J = 0.
- Measured: the maximum speed after 1 ms, relative to the Alfvén speed.

| grid | none | GLM | CT |
|---|---|---|---|
| 40×6 | 79 m/s (0.096 v_A) | 26 (0.031) | 19 (0.023) |
| 80×12 | 45 (0.055) | 34 (0.041) | 16 (0.019) |
| 160×24 | fails at 0.84 ms | 36 (0.043) | 26 (0.031) |

- **This does not converge.**
- The initial force residual isolates the problem:
  - The RMS acceleration converges at about first order: 4.9e5, 2.5e5, 7.0e4, 2.9e4 m/s² (40 to 320).
  - The maximum is stuck at 1.01e6 m/s² in the axis row (j = 0) at 320×48. The scale b²/(ρL) is 1.37e7.
- Mechanism (measured with the axis probe at 640 rows):
  - RZ cell averages are r-weighted, so they belong at volume centroids.
  - The minmod reconstruction assumes cell midpoints. Near the axis the gap is O(Δr), and B_r ∝ r there.
  - The axis cell B_r is correct (volume average 1.0094). Both sides reconstruct 1.682 at the first radial face against an exact 1.514 (+11%).
  - That O(1) face error does not shrink with refinement.
- Two fixes were tried and failed: an axis parity slope (no effect), and the exact volume-average cell B from ψ (it fixed the bulk but not the axis).
- Per the two-fix rule, the third candidate (centroid-referenced slopes, or reconstructing r-weighted quantities) was not attempted.
- **The gas core uses the same midpoint assumption.** It is invisible there because no gas case in the suite has a quantity that is large and ∝ r at the axis. Swirl (u_θ ∝ r) or any field will expose it.

### 4. Coil on the nozzle (`coil_80.csv`, `coil_160.csv`, `coil_axis_bz.png`, `coil_1T_ct_160_failure.png`)

- The default core nozzle (300 kPa, 300 K, air) starts from the prepared quasi-1D state and runs for 8 ms, the same as the core thrust study.
- The gas is treated as a perfect conductor. That is not physical for air; it is a transport test of the layout.
- **B = 0 check:**
  - MHD device thrust 517.5521 N (80×12) and 517.7730 N (160×24).
  - The gas core gives 517.5420 and 517.7718 N.
  - The difference is 2e-5 and 2e-6 relative: HLLD reduces to the gas solver as it should.

| B centre | grid | none | GLM | CT |
|---|---|---|---|---|
| 0.3 T | 80×12 | fails 3.9 ms | fails 5.4 ms | runs; thrust 517.31 N, face div B 3e-17 |
| 1.0 T | 80×12 | fails 3.9 ms | fails 3.9 ms | runs; thrust 517.13 N, face div B 2e-17 |
| 0.3 T | 160×24 | fails 1.7 ms | fails 5.7 ms | runs; thrust 517.70 N |
| 1.0 T | 160×24 | fails 1.8 ms | fails 4.1 ms | **fails 0.58 ms** (50 HLL fallbacks) |

- Every failure is a pressure-positivity loss. The minimum β at failure is 1e-8 to 5e-7.
- **Where the CT runs "succeed", the field has been swept out** (see `coil_axis_bz.png`).
  - Frozen-in field lines are carried downstream within about one transit time.
  - The field that remains is whatever flux the inlet admits (about 0.01 T at 1 T applied), compressed by the area ratio.
  - So the near-zero thrust change is not a result about coils.
  - This is ideal MHD behaving correctly under these boundary conditions (inferred from the flux-freezing estimate: inlet 0.007 T × area ratio 3 ≈ 0.02 T at the throat, against about 0.017 to 0.02 T measured).
- **The 160×24, 1 T CT failure** is at the outer wall at z ≈ 0.34 m:
  - The swept flux piles up against the wall.
  - There p falls to 0.013 Pa while B ≈ 0.9 T.
  - p = E − ½ρu² − ½b² then loses all precision.
- The 80×12 run survives only because it is more diffusive.
- Cost: HLLD with 9 variables runs at 6.0 to 6.8e6 cell-stage updates/s on one core. The fast speed adds 0.6% (0.3 T) to 12% (1 T) more steps.

## What broke, and whose problem it is

1. **Axis reconstruction (mesh-independent RZ issue; also in the gas core).**
   - Midpoint-based slopes on r-weighted averages give O(1) face errors at the axis.
   - Needed before any field or swirl work: centroid-referenced reconstruction, verified by the static-coil residual converging in its maximum as well as its RMS.
2. **div B: only CT is usable.**
   - Without cleaning, the solver fails on the static case at 160 and on every coil case.
   - GLM fails on every coil case (a single axis cell reached Mach 2000 in the 80×12, 1 T run).
   - CT through nodal ψ is exact to 1e-17 and cheap on this mesh.
3. **Low-β positivity.**
   - The total-energy form cannot hold p when ½b² ≫ p.
   - Candidates:
     - the split B = B0 + B1 form, where B0 is the analytic curl-free coil field, so B0² never enters E;
     - an entropy or internal-energy fallback.
   - Not attempted.
4. **Magnetic boundary conditions are physics, not numerics.**
   - Copying B at the wall makes it neither insulating nor line-tying. Field-line footpoints slide along it.
   - An insulating wall needs matching to the vacuum field outside. A conducting wall needs line-tying.
   - The inlet needs a stated magnetization of the source.
5. **Ideal MHD is the wrong model for a coil around a flowing gas.**
   - Infinite conductivity sweeps the applied field out.
   - Applied-field devices work through finite conductivity (magnetic Reynolds number), and for EP regimes through Hall/two-fluid effects (inferred from standard theory, not measured here).
6. **Low-β accuracy.** Acoustic waves take Alfvén-speed dissipation (table 2).

## Coil force accounting

- The conservative MHD ledger closes. Momentum residual is ≤ 5e-16 in every coil run.
- It cannot say what force the coil receives, because the Lorentz force is hidden in the Maxwell stress of the fluxes.
- The spike's diagnostics ∫(J×b)_z and ∫(J×b0)_z use a cell-centred curl of the total field. They disagree with each other and between grids (for example −0.055 and −16.8 N at 1 T, 80×12, CT). **They are not trusted.**
- Recommended structure (derived):
  - Carry B = B0 + B1, with B0 the coil vacuum field (curl-free in the domain) and B1 induced by plasma currents, CT on ψ1.
  - The force on the medium is J1×(B0+B1).
  - The J1×B0 part is exactly minus the force on the coils (the Biot-Savart pair between plasma and coil currents). It enters `Measurements::bodyAxialForce`, which already carries the reaction into device thrust.
  - J1×B1 is the plasma's self-force and reduces to Maxwell stress on the domain boundary.
  - This also removes B0² from the energy (item 3).

## Recommendation

**Keep the body-fitted structured RZ mesh for fields.**

- The measured case for it:
  - exact CT through one nodal scalar;
  - a trivial axis for ψ (ψ = 0);
  - HLLD working unchanged on curved faces;
  - the gas limit reproduced to 2e-6.
- AMReX (inferred, not measured): its RZ mode uses Cartesian cells with embedded-boundary cut cells at the curved nozzle wall.
  - CT on cut cells is not a mature capability. The practical route there is GLM-type cleaning, and GLM failed here.
  - AMReX fixes none of items 1 to 6.
  - Revisit it when a result needs adaptive refinement, for example an external plume or multi-body domains.
- Nothing measured here points to a third option.
- Before committing months to fields, in order:
  1. centroid reconstruction (benefits the gas core too);
  2. the split B0 + B1 form with the coil reaction in `bodyAxialForce`;
  3. explicit wall and inlet magnetic conditions;
  4. a conductivity model.
- Items 2 and 3 are verifiable with the static-coil and aligned cases already in the spike.

DECISION NEEDED: which conductivity model the first field physics targets. The options:

- resistive MHD with a magnetic Reynolds number from the plasma state (applied-field MPD / magnetic nozzle regime);
- Hall or two-fluid (most electric-propulsion regimes);
- keep ideal MHD for fusion-like high-Rm cases.

This sets the next physics milestone. It is Ben's choice, not an engineering one.
