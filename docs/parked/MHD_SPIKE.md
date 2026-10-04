> **PARKED (Ben, 2026-09-30): magnetic/plasma work waits until the chemical milestone is done. Historical design record, not active work; see docs/SESSION_HANDOFF.md for the current state.** The uncommitted Maeno 2013 case driver (magnet field verified, plume inputs drafted, no impulse result) was committed as parked work.

# MHD architecture spike (September 2026)

Question: can the current body-fitted RZ mesh and data layout carry magnetic fields (div B control, an MHD Riemann solver, axis treatment) before months are committed to them?

Answer: yes for the mesh and layout. The field transport machinery worked. What broke is (a) axis reconstruction, which is an RZ issue on any mesh and is also latent in the gas core, and (b) the physics model and magnetic boundary conditions. Neither is fixed by changing the mesh. Recommendation at the end.

Throwaway code: `spike/mhd_rz.cpp` (target `crucible_mhd_spike`), plot script `spike/plot_coil.py`. The spike is not wired into the app or the core and is not a framework. Raw tables and figures are in `docs/parked/mhd_spike/`. Every number below is **measured** unless labelled otherwise.

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
  - p = E − ½ρu² − ½b² then loses all precision. (Corrected by case 5: precision is not the cause. The split form never subtracts b0² and still fails.)
- The 80×12 run survives only because it is more diffusive.
- Cost: HLLD with 9 variables runs at 6.0 to 6.8e6 cell-stage updates/s on one core. The fast speed adds 0.6% (0.3 T) to 12% (1 T) more steps.

### 5. Split field B = B0 + B1 (`static_split.csv`, `split_coil.csv`, `split_1T_160_failure.png`), 30 September 2026

- **Built** (`split` mode, run with `crucible_mhd_spike split 80,160 30,100 0.008 dump`):
  - CT on the induced flux ψ1 only. The coil field b0 is steady, curl-free, and taken analytically at face centres.
  - The energy slot is E1 = p/(γ−1) + ½ρu² + ½b1² (Tanaka 1994), so b0² never enters the pressure.
  - HLLD runs on the total field (b0 is identical on both sides of a face, so only b1 jumps). The momentum flux then drops T(b0), and the cross stress b0b1 moves into a cell body force J1×b0. The energy flux drops b0 · (induction flux).
  - The axial integral of J1×b0 is `bodyAxial`, and device thrust = inlet + wall + bodyAxial − ambient.
- **Static coil (1 T):** max speed 1.2e-12 and 2.9e-12 m/s at 40 and 80, against 19 and 16 m/s for CT. Exact by construction (b1 stays 0), so this verifies only the bookkeeping.
- **Coil nozzle, 8 ms:**

| B centre | grid | CT thrust / result | split thrust / result | split bodyAxial | J1×b0 diagnostic | split min β |
|---|---|---|---|---|---|---|
| 0.3 T | 80×12 | 517.31 N | 519.50 N | 38.52 N | 36.18 N | 5.4 |
| 0.3 T | 160×24 | 517.70 N | 516.08 N | 44.61 N | 41.83 N | 5.9 |
| 1.0 T | 80×12 | 517.13 N | fails 7.73 ms | | | |
| 1.0 T | 160×24 | fails 0.58 ms | fails 0.91 ms | 24.867 N | 24.863 N | |

- **Ledger:** momentum residual ≤ 5e-16 in every split run.
  - `bodyAxial` agrees with the independent cell-curl J1×b0 diagnostic to 6% at 0.3 T and to 2e-4 just before the 1 T failure.
  - The coil reaction now has a place in device thrust.
- **The split does not fix low-β positivity.** It fails at 160 a little later than CT, and it also fails at 80, where CT survived.
- **Mechanism (measured):** the failure cell is at the wall at z = 0.44 m, downstream of the coil, not at it. There:
  - p = 7.6e-4 Pa and u_z = −84 m/s;
  - the induced field is b1 = 0.61 T while b0 = 0.009 T.
- **Enclosed flux 2πψ along the wall, initial → at failure (1 T, 160):**
  - z = 0.15 m: 0.35 → 0.08 mWb;
  - z = 0.44 m: 0.027 → 1.20 mWb;
  - exit: 0.008 → 1.02 mWb.
- The throat's coil flux (about 1.2 mWb) has been carried downstream by the flow. Its magnetic pressure (about 1.4e5 Pa at 0.6 T, derived) is far above the expanded gas pressure, so it pushes gas back and off the wall, and the fluid pressure goes to zero.
- So the failure is ideal flux freezing plus the sliding-footpoint wall (items 4 and 5 below), not a discretization error. Split cannot fix it, and neither can an internal-energy fallback.
- At 0.3 T both forms also sweep the flux out (enclosed wall flux at 8 ms is −0.02 mWb for CT and −0.06 mWb for split at 160, from 0.39 mWb at the throat).
- **Inferred:**
  - CT keeps less field (max |B| 0.04 T against 0.24 T for split) because its reconstruction of the smooth b0 creates face jumps. HLLD dissipation then turns those jumps into spurious EMF, which acts as numerical resistivity. The same defect causes CT's non-converging static residual.
  - The 80×12 CT survival at 1 T is that numerical resistivity.
  - Neither form's 0.3 T thrust is converged between 80 and 160 (split moves 3.4 N), because the 8 ms state is mid-way through the sweep.
- Two-fix rule: this was the first fix tried on low-β positivity. The mechanism is in the model, so I stopped here rather than trying numerical fixes.

### 6. Wall magnetic conditions and positivity at 1 T (`wall_1T.csv`, `wall_1T_robust.csv`, `wall_1T_robust.png`, `vacuum_check.csv`), 30 September 2026

The question (from the outside review, `docs/parked/ASTRA_PLASMA_MODEL_2026-09-30.md`): was the 1 T failure caused by the wall condition? Short answer, measured: **no. It was a pressure-recovery failure in the energy equation. Changing the wall condition alone did not prevent it. Once pressure recovery was made robust, the result depended at O(1) on the wall condition.**

What was added (split form only; b1 is the induced field):

- Three wall modes for b1.
  - `transparent`: the ghost copies b, which was the old behaviour.
  - `insulating`: b1 outside the wall is the vacuum field, solved from the wall flux.
  - `conducting`: ψ1 is frozen at the wall nodes, so the wall holds its flux (line-tying), and the ghost b_n is reflected.
  - For insulating and conducting, the inlet ghost has b1 = 0, so the source is unmagnetised by plasma currents.
- Vacuum exterior. Linear-triangle FEM for div((1/r) grad ψ) = 0 on an annulus from the wall to radius Rf, with ψ = 0 at Rf and b_r = 0 on the end planes. It is precomputed once as a matrix from wall ψ to the exterior field at each wall face.
  - Checked against a current loop (r = 1 cm, z = 0.3 m, 1 kA) placed outside the nozzle. Field error at the wall faces, relative to the maximum field, rms (measured):

    | nz | 40 | 80 | 160 | 320 |
    |---|---|---|---|---|
    | rms error | 0.058 | 0.039 | 0.014 | 0.0043 |

  - Rf from 0.3 to 3 m changes this by under 10%. The maximum error is at the loop's closest point and falls from 0.22 to 0.037.
- Positivity robustness.
  - The reconstruction was already positivity-limited: primitive minmod keeps each face ρ and p between the neighbouring cell values.
  - Added a dual energy (an entropy variable ρs, advected with the mass flux). When the thermal pressure from total energy is below 1e-3 of the total energy, pressure is recovered from the entropy instead.
  - After each step the two are synchronised. Where E is used, s is reset. Where s is used, E is reset, and the change in E is booked as `sync_J` (not hidden in the energy residual).
  - Stage rejection: a stage with a non-admissible state is retried at the same dt with first-order reconstruction and HLL (the HLLE fallback), then with dt halved. These retries are counted in `robust_steps`.

Results at 1 T, coil at the throat, 8 ms target (measured; wall time on the Mac):

| nz | wall | robust | outcome | device thrust (N) | coil reaction body_axial (N) | sync (J) | robust steps |
|---|---|---|---|---|---|---|---|
| 80 | transparent | no | failed 7.73 ms | | | | |
| 80 | insulating | no | failed 1.33 ms | | | | |
| 80 | conducting | no | failed 1.69 ms | | | | |
| 160 | transparent | no | failed 0.909 ms | | | | |
| 160 | insulating | no | failed 0.248 ms | | | | |
| 160 | conducting | no | failed 0.613 ms | | | | |
| 80 | transparent | yes | 8 ms | 318.0 | 226 | 13.9 | 18 of 15409 |
| 80 | conducting | yes | 8 ms | 43.8 | -486 | 6.5 | 14 of 30458 |
| 160 | transparent | yes | 8 ms | 519.4 | 43.8 | 51.7 | 6734 of 90257 |
| 160 | conducting | yes | blew up 2.59 ms | | | | 20302 of 48858 |

No-field gas thrust at the same conditions is 517.5 N (measured, earlier runs).

What this says:

1. **Mechanism of the original failure (measured).** In every non-robust failure, the failing cells had ordinary density (0.1 to 0.7 kg/m^3) but pressure of 1e-6 to 1e-3 Pa (T about 1e-5 K). They were not evacuated. The thermal energy came out as a small difference of large magnetic and kinetic energies and lost its positivity. That is a numerical failure of pressure recovery, and it occurred with every wall condition. Dual energy fixes it in three of the four runs.
2. **The wall condition matters at O(1) (measured).** The 80-cell runs differ only in the wall condition: 318 N against 44 N of device thrust, and +226 N against -486 N of coil reaction. In the figure, the conducting wall holds about 1 T at the throat, and the flow separates into a fast core jet and a slow wall layer. The transparent wall lets the flux leave.
3. **The transparent result is not mesh-converged (measured).** Its thrust is 318 N at 80 cells and 519 N at 160. The 160 run needed 7.5% robust steps. None of these numbers is a result. They show sensitivity.
4. **Conducting at 160 is a new failure, at the outlet, not the wall (measured location, inferred cause).** The blown-up cells are in the last two columns (i = 157 to 159), mid-radius, with u_z about -14.5 km/s (inflow through the outlet) and |b| up to 4e8. The admissibility check tests positivity only, so an unbounded growth that stays positive was accepted. Inferred cause: the extrapolation outlet becomes ill-posed when the flow reverses there. The conducting wall's slow layer reaches the outlet (see the 80-cell u_z panel), and an extrapolation boundary feeds whatever comes in. Not yet fixed. It needs a characteristic outlet condition, and that belongs with the chamber or plume boundary work, not here.
5. **Insulating wall: stopped after two fixes (two-fix rule). The mechanism is written down here.**
   - Fix 1 put the exterior vacuum b1 into the wall ghost. The interior and exterior tangential b1 differed by up to 0.13 T, while b_n matched. The Riemann solver turned that jump into wall heating: ρ 0.003 kg/m^3, p 8e5 Pa, u_z -1400 m/s at the wall.
   - Fix 2 ran HLLD with a jump-free ghost and added the wall current's stress and Poynting flux as explicit corrections. It produced a 4 km/s wall jet, b1z of 1.5 T in one cell, and thrust of 1039 to 1422 N.
   - Mechanism: an insulating wall carries no current, so a tangential field jump between plasma and vacuum is a current sheet inside the plasma. In ideal MHD nothing sets that sheet's thickness. It sits in one cell layer, and its whole force (fix 2) or its whole dissipation (fix 1) lands on that layer. In a real device the sheet has a resistive thickness. So the insulating condition is not posed properly without resistivity (task 2). It will be retried once η J is in.

Answer to the review's point: the sweep-out under the earlier transparent condition does not by itself show that ideal MHD is unsuitable, and the review is right about that. A wall that holds its flux keeps the applied field in place in this run (conducting, 80 cells). But the earlier crash was not a boundary-condition failure. It was pressure recovery, and that is fixed. The physically intended condition here (an insulating nozzle with an external coil) cannot be run in ideal MHD on this mesh. That is consistent with the resistive-first default.

### 7. Resistive term (`resistive_decay.csv`), 30 September 2026

What was added:

- E_theta = eta_m J_theta at the nodes, with eta_m = eta/mu0 the magnetic diffusivity.
  - The term enters the flux function as dpsi/dt = -r E_theta, so CT and div B = 0 are unchanged.
  - J_theta = -div((1/r) grad psi) in the meridional plane, from two linear triangles per cell with a lumped mass.
- Boundary conditions:
  - The boundary term carries the tangential field: the interior value (zero-gradient) at a transparent wall and at the end planes; zero at the inlet of a nozzle whose source is unmagnetised; the matched exterior vacuum field at an insulating wall (no surface current, so b_t is continuous).
  - At a conducting wall, E = 0 (flux frozen).
- Energy:
  - Total energy keeps its conservative form. The resistive Poynting flux E_theta (b_z n_r - b_r n_z) of the induced field is added at every face.
  - Joule heating is then implicit (magnetic energy lost becomes internal energy). It is also integrated separately as a diagnostic, eta_m J^2 V, and added to the entropy variable.
- Explicit time step limit: dt <= 0.25 / (eta_m (1/dz^2 + 1/dr^2)).
- Conductivity options: constant, or transverse Spitzer.
  - Spitzer is eta_perp = 1.03e-4 Z lnL T_e^-1.5 ohm m (NRL Formulary 2019 p.34).
  - Electron-ion Coulomb log: 24 - ln(sqrt(n_e[cm^-3]) / T_e[eV]) above 10 Z^2 eV, 23 - ln(sqrt(n_e) Z T_e^-1.5) below, floored at 2. This is the same formula as in `docs/parked/VALIDITY_MONITOR.md`.
  - Spitzer assumes full ionisation at a stated Z, with n_e = Z rho / m_i.
- Not included: resistive diffusion of b_theta. It is zero in every case so far.

Spitzer check against hand values from the Formulary (measured, against derived):

| n_e = 1e20 m^-3, Z = 1 | code lnL / hand lnL | code eta (ohm m) / hand eta (ohm m) |
|---|---|---|
| T = 100 eV | 12.487 / 12.49 | 1.2862e-6 / 1.287e-6 |
| T = 1 eV | 6.882 / 6.88 | 7.0884e-4 / 7.09e-4 |

Decay cases. Setup:

- A single diffusion eigenmode in heavy gas at rest (rho = 1e4, p = 1, b = 0.01, so u x b / eta J ~ 1e-5).
- eta_m = 0.01 m^2/s, run for one field e-folding time.
- Exact answer: psi ~ exp(-eta_m lambda t). The internal energy gained must equal the magnetic energy lost.

The cases:

- radial: an annulus at r ~ 1e6 (the Cartesian limit), b_z = b cos(pi x / W). This is the Director's sinusoid B ~ exp(-(eta/mu0) k^2 t).
- axial: the same annulus, b_r = b cos(2 pi z / L).
- bessel: a full cylinder with the axis, b_z = b J0(alpha r), J1(alpha R) = 0. Run with both conducting and transparent walls.

Results (measured):

| case | n | decay rate / exact - 1 | psi L2 rel. error | Joule integral / magnetic loss - 1 | internal gain / magnetic loss - 1 |
|---|---|---|---|---|---|
| radial | 16, 32, 64, 128 | -3.2e-3, -8.0e-4, -2.0e-4, -5.0e-5 | same | -9.6e-3 ... -1.5e-4 | below 1e-7 |
| axial | 16, 32, 64, 128 | -1.3e-2, -3.2e-3, -8.0e-4, -2.0e-4 | same | -3.8e-2 ... -6.0e-4 | 7.2e-3, 3.4e-3, 1.8e-3, 1.2e-3 |
| bessel (both walls) | 16, 32, 64, 128 | -4.8e-3, -1.2e-3, -3.0e-4, -7.4e-5 | same | -9.3e-3 ... -1.5e-4 | -1.1e-7 |

The energy residual is below 1.2e-14 in every run.

- The decay rate, the field profile and the Joule integral all converge at second order.
- One fix along the way (measured): the first version took 1/r at each triangle's centroid. A uniform field (psi = r^2) then has a spurious current of order h/r, which measured as a first-order psi error growing toward the axis (16% at r = 0.06 with 32 cells). Taking 1/r at the cell centre for both triangles makes psi = r^2 exact on straight columns. The error is now uniform in r and second order.
- What is still first order (measured): the energy entering through the open, zero-gradient end planes in the axial case (0.12% of the Joule energy at 128 cells). The end-node current comes from a one-sided boundary closure. Physically that flux is zero. Field evolution is unaffected (second order). Open boundaries need a better closure before a plume run.

### 8. Resistive nozzle at 1 T: wall conditions retried (`wall_1T_eta.csv`, `wall_1T_eta.png`), 30 September 2026

Same nozzle and coil as case 6, 8 ms, robust mode on. eta_m is constant. It is an **assumed** scan parameter, not a property of the gas: the spike gas is cold air, which does not conduct, so Spitzer does not apply to it. eta_m = 1 and 0.1 m^2/s give magnetic Reynolds numbers of about 10 and 100 (derived: u L / eta_m with u = 500 m/s, L = 2 cm). The insulating wall keeps the stress and Poynting corrections from case 6 (fix 2), now with the resistive boundary condition.

Results (measured):

| eta_m (m^2/s) | wall | device thrust, 80 / 160 cells (N) | coil reaction body_axial, 80 / 160 (N) | min beta, 160 | robust steps |
|---|---|---|---|---|---|
| 1 | transparent | 395.8 / 448.3 | -53.7 / -21.2 | 3.65 | 0 |
| 1 | insulating | 316.4 / 315.5 | -158.3 / -154.0 | 0.268 | 0 |
| 1 | conducting | 312.0 / 308.7 | -157.7 / -158.4 | 0.279 | 0 |
| 0.1 | transparent | 507.1 / 515.7 | 41.2 / 30.3 | 151 | 0 |
| 0.1 | insulating | 166.7 / 160.0 | -259.6 / -284.7 | 0.200 | 0 |
| 0.1 | conducting | 151.9 / 148.3 | -321.1 / -325.5 | 0.237 | 0 |

No-field thrust is 517.5 N. All 12 runs reached 8 ms with zero retries and zero sync energy, and the energy residual is below 7e-15.

What this says:

1. **With resistivity the insulating wall works.** This answers the case 6 stop. The wall current sheet now has a resistive thickness and needs no fallback. Insulating and conducting walls agree to within 1% at Rm ~ 10 and to within 8% at Rm ~ 100. They change by 0.3 to 4% from 80 to 160 cells.
2. **The transparent wall is unphysical, now measured directly.** At Rm ~ 100 the figure's bottom row shows a total |B| of at most 0.055 T. Plasma currents have cancelled the 1 T coil field, because a zero-gradient wall lets induced field leave and return without any exterior constraint. Its thrust (516 N) is the no-field value for the wrong reason. Transparent is dropped as a physical option. It stays as a numerical reference only.
3. **The Director's point is confirmed.** With the field maintained by a proper wall condition, the coil field holds near 1 T at the throat. It turns the flow into a fast core jet with a slow wall layer downstream (figure, rows 1 to 4). This geometry is a coil on a gas nozzle, not a thruster design. The thrust numbers are spike sensitivities, not predictions: the conductivity is assumed, and the gas is not a plasma.
4. **The outlet did not block these runs.** The conducting, ideal, 160-cell failure from case 6 did not recur once resistivity was in. The extrapolating outlet is still ill-posed for backflow, and a plume domain will need a characteristic condition. That is recorded, not fixed, because nothing here was blocked.

## What broke, and whose problem it is

1. **Axis reconstruction (mesh-independent RZ issue; also in the gas core). Fixed in the gas core on 30 September 2026.**
   - Midpoint-based slopes on r-weighted averages give O(1) face errors at the axis.
   - Needed before any field or swirl work: centroid-referenced reconstruction, verified by the static-coil residual converging in its maximum as well as its RMS.
2. **div B: only CT is usable.**
   - Without cleaning, the solver fails on the static case at 160 and on every coil case.
   - GLM fails on every coil case (a single axis cell reached Mach 2000 in the 80×12, 1 T run).
   - CT through nodal ψ is exact to 1e-17 and cheap on this mesh.
3. **Low-β positivity.**
   - The split form (case 5) removes b0² from the energy and T(b0) from the fluxes. It is kept: it makes the static coil exact and puts the coil reaction in `bodyAxial`.
   - It does not prevent the failure. Ideal flux freezing sweeps the coil flux downstream, and that flux evacuates the wall region.
   - This needs the model items 4 and 5 (wall magnetic conditions; resistivity, now the Director default), and possibly a stated treatment of near-vacuum regions (`RESEARCH.md`, "Expansion into vacuum is a model decision").
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
- Implemented in the split form (case 5). There `bodyAxial` is the discrete J1×b0 integral, consistent with the momentum ledger to roundoff.
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

Conductivity, resolved 30 September 2026 as a Director default (reversible by Ben; recorded in `RESEARCH.md`): resistive MHD first, because plasma detachment, the central magnetic-nozzle question, needs field-line slippage that ideal MHD cannot represent. The Hall term is designed in as the next increment. Ideal MHD remains a limiting case (infinite conductivity), not the target.

Item 1 (centroid reconstruction) landed in the gas core on 30 September 2026; see `docs/evidence/IMPLEMENTATION.md`.
