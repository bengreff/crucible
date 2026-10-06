# C2: SST-2003 turbulence in the engine

Status: checks 1, 2, 3 and 5 pass. Check 4 (fully developed pipe against an independent solver) matches the reference on nr 16, 32 and 64, but its mass-budget criterion fails on all three, as reported below. The nr 64 column comes from a complete run on backhouse (5 October), after the Mac run was stopped at 3.5 of 4 ms when the session was parked on 4 October. Criteria were stated in the header of `tests/turbulence_tests.cpp` before the first run of each check (4 October 2026), and the amendments are dated there.

Raw outputs:
- checks 1 to 3: `turbulence_verification_2026-10-04.txt` (before the open-face fix) and `turbulence_supply_run3_2026-10-04.txt` (after it, with check 5);
- check 4: `turbulence_pipe_2026-10-04.txt` (Mac; nr 16 and 32 complete; nr 64 stopped at 3.5 ms) and `turbulence_pipe_2026-10-05_backhouse.txt` (backhouse: WSL2 Ubuntu, gcc 13, Cantera 3.2.0, commit 4e2ee5c; the pipe test is unchanged since 96422fa; nr 16, 32 and 64 complete);
- check 5: `turbulence_supply_run1_2026-10-04.txt`, `turbulence_supply_run2_2026-10-04.txt` and `turbulence_supply_run3_2026-10-04.txt`.

All numbers are measured unless they are marked derived, inferred or estimated.

## What was added

- **Model** (807701b; TECHNICAL_PLAN step 7). SST-2003 as in TMR `sst.html`.
  - rho k and rho omega are carried by the mass flux, and k is part of E.
  - mu_t = rho a1 k / max(a1 omega, S F2).
  - The turbulent stress includes −2/3 rho k. Conduction uses lambda + cp mu_t / Pr_t and diffusion D + mu_t / (rho Sc_t).
  - Sources are integrated by MPRK22 (Patankar-weighted, positive for any non-negative input) inside the Strang split, at fixed rho and E. The cross-diffusion term is split by sign.
  - Amended 5 October 2026, evening: the production is TMR's SSTs form, P = mu_t S^2 in both equations, with -2/3 rho k delta_ij kept in the stress (TECHNICAL_PLAN step 7; [NASA TMR](https://tmbwg.github.io/turbmodels/sst.html)). It replaces the exact P and the CRUCIBLE-only omega fix of that morning (the reason is in PASR_C2.md, criterion 4). Checks 1, 2, 3 and 5 (`turbulence_verification`) pass on it ([turbulence_checks_ssts_2026-10-05.txt](pasr/turbulence_checks_ssts_2026-10-05.txt)). Check 4 (the pipe) on nr 16 and 32 ([pipe_ssts_2026-10-05.txt](pasr/pipe_ssts_2026-10-05.txt)):
    | | nr 16, SSTs | nr 16, before (omega fix) | nr 32, SSTs | nr 32, before (omega fix) |
    |---|---|---|---|---|
    | u_b against the reference (limit 1e-3) | 5.03e-5 | 2.07e-4 | 1.05e-5 | 6.51e-5 |
    | c_f against the reference (limit 2e-3) | 1.01e-4 | 4.15e-4 | 2.10e-5 | 1.30e-4 |
    | T_axis - T_wall against the reference (limit 1e-2) | 5.90e-5 | 5.03e-4 | 1.08e-5 | 1.53e-4 |
    | mass budget (limit 1e-12) | 5.68e-12 FAIL | 2.38e-11 FAIL | 2.37e-11 FAIL | 2.35e-11 FAIL |
    - The agreement improves 4 to 8 times. The reference, `tools/sst_pipe_1d.py`, already used P = mu_t S^2, so the engine and the reference now solve the same model; the earlier gap included the dilatation and -2/3 mu_t div^2 parts that only the engine had (inferred, not separated). The pipe's non-zero divergence noted in the morning's amendment no longer enters the production; where it arises is still not measured.
    - The mass-budget criterion fails as before (open item, unchanged cause). nr 64 was not rerun (the Director's 21:45 instruction: nr 16 and 32 only).
  - Earlier the same day (superseded): in omega's production only, the dilatation part used omega for rho k / mu_t. Its pipe record is [pipe_after_omega_fix_2026-10-05.txt](pasr/pipe_after_omega_fix_2026-10-05.txt).
- **Walls** (ca85613, 807701b).
  - The wall distance is exact for each centroid, to the side wall and to flagged plate rings.
  - No-slip faces carry mu_t = 0 and k = 0, with omega_w = 10 · 6 nu / (beta1 d1^2).
  - Slip walls pass no k or omega flux.
- **Turbulent supplies** (dc4cbd8).
  - A Chamber supply declares an intensity I and a viscosity ratio. Its face carries k = 3/2 (I u)^2 and omega = rho k / (ratio mu).
  - The energy relation is h(T) + u^2/2 + 5/3 k = h0, and the momentum flux is g u (1 + I^2) + p.
- **rho k clip** (dc4cbd8).
  - After each RK stage, rho k = max(rho k, 0). E is not changed, so the energy budget stays exact.
  - The rho k V added is summed in `Measurements::clippedTurbulentEnergy` and judged in check 5.
  - Why it is needed: on contoured rings the tangential part of the face gradient makes k diffusion non-sign-preserving (check 5, run 1).

## Check 1: wall distance

| Case | Criterion (of R) | Measured | Verdict |
|---|---|---|---|
| Straight duct, equal rings and b = 2: d = R − r | < 1e-13 | 0 | pass |
| Nozzle, b = 1.5, against a per-segment ternary search | < 1e-13 | 4.0e-16 | pass |
| Nozzle with plate rings 0, 3, 4, 7 as walls | < 1e-13 | 2.5e-16 | pass |
| No walls: every distance infinite | yes | yes | pass |

## Check 2: decaying homogeneous turbulence (source split, energy accounting)

N2 at rest, 300 K, 1e5 Pa, slip walls, k0 5e4 m^2/s^2, omega0 1e5 1/s, marched to beta2 omega0 t = 0.8.

| Criterion | Limit | Measured | Verdict |
|---|---|---|---|
| 2a observed order, N 16 to 32 | ≥ 1.9 | 1.998 | pass |
| 2a error at N 32 | < 1e-4 | 2.17e-5 | pass |
| 2b against the MPRK22 formula evaluated in the test | < 1e-12 | 1.1e-16 | pass |
| 2c rho, E unchanged; k, omega uniform | < 1e-13 | 0 | pass |
| 2c T against e(T0) + k0 − k | < 1e-8 K | 5.7e-14 K | pass |
| 2d stiff (beta2 omega0 dt = 300): positive, finite, decreasing every step | yes | yes | pass |

The stiff case is positive but not accurate. After 10 steps the test measured k 2.93 m^2/s^2 against the exact 8.30, and omega 7.05e3 1/s against 3.33e4. The criterion asks only for positivity. The limit is recorded in TECHNICAL_PLAN: where the turbulence time scale is far below the step, the source is a damped, positive approximation.

## Check 3: turbulent transport operator (slip-wall duct, F1 = F2 = 0, mu_t about 8 mu)

The table gives the volume-weighted L1 truncation error against the exact cell average, on grids 32x8 to 256x64 (run 3, after the open-face fix).

| Quantity | 128x32 | 256x64 | Order | Limit | Verdict |
|---|---|---|---|---|---|
| Axial momentum | 2.20e-3 | 5.32e-4 | 2.046 | ≥ 1.8 | pass |
| Radial momentum | 4.65e-3 | 1.23e-3 | 1.920 | ≥ min(1.8, p_lam − 0.1) = 1.8 | pass |
| Energy | 1.83e-3 | 5.04e-4 | 1.859 | ≥ 1.8 | pass |
| Species | 2.15e-3 | 5.10e-4 | 2.073 | ≥ 1.8 | pass |
| rho k | 2.20e-3 | 5.35e-4 | 2.043 | ≥ 1.8 | pass |
| rho omega | 2.19e-3 | 5.30e-4 | 2.046 | ≥ 1.8 | pass |

- 3a: mu_t equals rho k / omega within 4.4e-16 (limit 1e-12). Pass.
- The laminar orders on the same field are: axial 1.973, radial 1.953, energy 1.982, species 1.993.
- The open-face fix (1e0fd9e) raised the laminar radial order on this straight duct from 1.589 to 1.953. The radial-momentum failure on the contoured duct in TRANSPORT_C2.md remains: this check does not clear it.

## Check 4: fully developed pipe against `tools/sst_pipe_1d.py`

Case: N2, R 0.4 mm, length 2R, no-slip wall at 300 K, axial body force 2.95e5 N/m^3 (Re_tau about 182), rings clustered with b = 2.
- Each run starts from the reference solution on the same rings and marches 4 ms.
- The reference is an independent steady solver (SST_REFERENCE.md). It shares the equations, the rings, the wall-omega rule and, in this one-dimensional flow, the centroid-difference gradients. So the comparison measures the implementation, not the grid.

How the check got to its final form:
1. **nz.** It was stated as nz 2, but the mesh needs at least 4 axial cells, and the first launch stopped before any output. It now uses nz 4; the flow is axially uniform. The header records the change.
2. **Open-face bias, found by this check (1e0fd9e).**
   - On the first attempt, started from the reference, the engine's bulk velocity rose from 105.14 to 113.16 m/s by 0.5 ms, and to 136.2 m/s by 1.5 ms.
   - Measured cause: an axial pressure gradient of about 6.5e4 Pa/m (22% of the body force) grew between the open ends. At t = 0, only the radial-momentum transport derivative varied along z, by ±98 off the mean on the axis ring at the end columns.
   - Mechanism: the zero-gradient neighbour for an open face was placed at the face midpoint, (−dz/2, ring middle − r_c). That also asserted a zero radial derivative toward a radially offset point. It biased du_z/dr in the end columns, and with it tau_zr and the k production.
   - Fix: that neighbour now sits at the face on the cell's own axial line. Afterwards the transport derivative is uniform along z to 13 digits, and the run below holds steady.
   - The laminar operator had the same bias, scaled by mu. The transport checks were rerun on the fix (TRANSPORT_C2.md).

Results:

| Criterion | Limit | nr 16 | nr 32 | nr 64 |
|---|---|---|---|---|
| 4a u_b change, last 1 ms | < 1e-5 | 2.0e-10 | 6.4e-11 | 2.1e-11 |
| 4a T_axis change, last 1 ms | < 1e-5 | 5.7e-10 | 1.3e-10 | 6.3e-12 |
| 4a wall shear force against body force | < 1e-5 | 1.4e-9 | 1.5e-10 | 1.4e-11 |
| 4b mass budget | < 1e-12 | **3.42e-11 FAIL** | **1.93e-11 FAIL** | **6.14e-11 FAIL** |
| 4b axial momentum budget | < 1e-9 | 1.8e-15 | 9.3e-15 | 1.6e-14 |
| 4b energy budget | < 1e-9 | 2.1e-11 | 7.2e-11 | 3.6e-10 |
| 4c u_b against the reference | < 1e-3 | 2.03e-4 | 6.37e-5 | 1.60e-5 |
| 4c c_f against the reference | < 2e-3 | 4.06e-4 | 1.27e-4 | 3.19e-5 |
| 4c T_axis − T_wall against the reference | < 1e-2 | 4.93e-4 | 1.50e-4 | 3.83e-5 |

The nr 16 and 32 columns are the Mac run; the nr 64 column is the backhouse run. Backhouse also reran nr 16 and 32. Its printed histories of u_b, c_f, T_axis and step count match the Mac's in every printed digit; only the force-balance residual (last digit or two) and the budgets differ, at round-off: mass 3.40e-11 and 1.68e-11, energy 2.06e-11 and 6.73e-11 (measured).

| nr | engine u_b (m/s) | reference u_b | engine c_f | reference c_f | engine T_axis (K) | reference T_axis |
|---|---|---|---|---|---|---|
| 16 | 105.16278 | 105.14145 | 0.00937609 | 0.00937989 | 304.78557 | 304.78321 |
| 32 | 102.75047 | 102.74393 | 0.00982151 | 0.00982276 | 304.55590 | 304.55521 |
| 64 | 101.10569 | 101.10408 | 0.01014366 | 0.01014398 | 304.41319 | 304.41302 |

The engine-reference difference falls by a factor of about 3.2, then 4.0, per ring halving (u_b 2.03e-4, 6.37e-5, 1.60e-5; derived). So it behaves like a second-order discretization difference between the two codes. It is not zero although the two codes share the discretization; its cause has not been isolated.

The nr 64 run took 8.91e6 steps and 2868 s wall on backhouse (measured).

The grid error is the model's, and is reported, not judged. The wall omega rule makes the solution first order (SST_REFERENCE.md). Against the reference's Richardson limits (c_f 0.0104865, u_b 99.439 m/s), the engine is:
- nr 16: u_b +5.8%, c_f −10.6%;
- nr 32: u_b +3.3%, c_f −6.3%;
- nr 64: u_b +1.7%, c_f −3.3%.

The deficit roughly halves per ring halving (c_f −10.6%, −6.3%, −3.3%; derived), as a first-order solution should.

**4b, the mass budget.** It fails on all three grids, and the threshold stays as stated.
- Mechanism (inferred, not demonstrated): rounding of the stored cell densities, accumulated over millions of steps.
  - Every step stores each cell density rounded to a double, an error of up to eps/2 of the cell's mass. In this steady flow the true change per step is far below that, so the rounding need not average out.
  - Worst-case bound (derived): (eps/2) × M × steps, with eps/2 = 1.11e-16 for doubles. That is 2.0e-10 of the initial mass on nr 16 (1.81e6 steps), 4.4e-10 on nr 32 (3.99e6 steps) and 9.9e-10 on nr 64 (8.91e6 steps). The measured 3.4e-11, 1.9e-11 and 6.1e-11 are inside it.
  - Correction (5 October): the first version of this bullet gave 1.0e-10 and 2.2e-10. Those used eps/2 = 5.6e-17, half the right value, so they were a factor of 2 too low. The conclusion (measured values inside the bound) is unchanged.
  - The ledger side is far smaller. In 4 ms, about 526 (nr 16), 514 (nr 32) and 506 (nr 64) domain masses pass the open ends (derived: u_b × 4 ms / 2R). The in-minus-out flux, rounded each step, can contribute at most about eps × 526 ≈ 6e-14 (derived).
- A scratch diagnostic (`/tmp/c2out/pipediag/budget.cpp`, nr 16, not in git) is consistent with growth by step count. At 1 ms the mass budget was 4.3e-12 at CFL 0.4 (4.5e5 steps) and 1.6e-11 at CFL 0.2 (9.1e5 steps). Halving the step at the same physical time multiplied the error by 3.7, more than a random walk (1.4) or a linear bias (2) predicts. The budget also changes sign during the run, so the scaling is not clean, and the mechanism stays inferred.
- The energy budget is the same size (2.1e-11, 7.2e-11 and 3.6e-10). It passes only because its limit, 1e-9, was stated as a round-off allowance for about 4e6 steps (nr 64 took 8.9e6 and still passes); the mass limit was not given that allowance. A mass error of 3e-11 of the contents has no physical consequence here. It is reported as the stated criterion failing, not re-judged.
- What would test the mechanism: a compensated update of the cell densities (carry each cell's rounding residual into its next step) on nr 16. It has not been done.

## Check 5: turbulent supplies (Chamber, N2)

Run 3 results:

| Criterion | Limit | Measured (worst) | Verdict |
|---|---|---|---|
| (a) energy relation h + u^2/2 + 5/3 k = h0, subsonic g 50 and choked g 2000 | < 1e-12 of cp T0 | 1.0e-16 | pass |
| (a) subsonic face p on the outgoing characteristic; u = g R T / p | < 1e-12 | 1.6e-15; 0 | pass |
| (a) choked: u = a(T) | < 1e-12 | 0 | pass |
| (a) k = 3/2 (I u)^2 | < 1e-14 | 2.2e-16 | pass |
| (a) omega = rho k / (ratio mu) | < 1e-12 | 2.2e-16 | pass |
| (a) five fluxes from the face state | < 1e-14 | 0 | pass |
| (b) I = 0 or ratio = 0 rejected with turbulence, accepted without | yes | yes | pass |
| (c) k against plug-flow decay, nz 64 and 128 | < 2e-3 | 3.51e-4, 1.37e-4 | pass |
| (c) omega against plug-flow decay, nz 64 and 128 | < 2e-3 | 3.05e-4, 1.16e-4 | pass |
| (d) mass, energy, axial momentum budgets, (c) runs and no-slip chamber | < 1e-12 | 7.8e-15 | pass |
| (d) no-slip chamber: k, omega positive and finite | yes | yes | pass |
| (d, added) clipped rho k V over the integral of rho k, every run | < 1e-12 | 0 | pass |

- (c) reported: the observed order of the larger error is 1.68 (nz 32 to 64), then 1.36 (64 to 128). The order falls as the error nears the axial-diffusion floor estimated before the run, about 1e-4 (inferred). The change over the last 0.25 ms is 7e-6 to 8e-6.
- The clip ledger is 0 in every run 3 case. In scratch runs of the no-slip chamber with the run 1 ambient, the clip added at most 6e-89 of the integral.

How it got there (amendments dated in the test header; thresholds unchanged):
1. **Run 1.**
   - (c) failed on every grid (k errors 1.1e-2 to 1.2e-2). The 1.5 ms march ended while the chamber pressure was still relaxing from the start-up overshoot (e-fold about 0.28 ms, scratch diagnostic). The march is now 4 ms.
   - (d) could not finish its 17th step: rho k reached −1.3e-198 in a quiescent cell on the contoured rings.
   - The k diffusion there is not sign-preserving, and no step reduction cures that. This led to the clip and its ledger, with the added criterion.
2. **Run 2.**
   - (d) stopped at 134 µs with omega negative on the axis ring next to the outlet.
   - Mechanism: with the ambient omega of 1 1/s, the negative cross-diffusion term −2 (1 − F1) rho sigma_w2 |∇k·∇ω| / ω drives omega to zero in about 1e-8 s once the jet's k front arrives (omega^2 falls linearly). Diffusion then carries it below zero.
   - That ambient (k 0, omega 1, so mu_t/mu = 0) also broke the plan's own rule that ambient gas carries a viscosity ratio in TMR's free-stream range, 1e-5 to 1e-2.
   - Fix: the fill and ambient of (d) now take the Spalart and Rumsey free-stream values, k = 9e-9 a^2 and omega = 1e-6 rho a^2 / mu, which give mu_t/mu 0.009 and u' about 0.03 m/s.
   - (c) keeps k 0 and omega 1: it passed and has no backflow.
3. **Run 3.** 0 failures.

These were two different mechanisms, each fixed where it arose, so the two-failed-fixes rule was not reached.

## What this does not cover

- Nozzle-inlet boundaries do not yet carry turbulent inflow; only Chamber supplies do.
- Validation (pipe DNS at Re_tau 180, the TMR axisymmetric jet) comes after verification.
- The flame runs in FLAME_C2.md predate the open-face fix.
