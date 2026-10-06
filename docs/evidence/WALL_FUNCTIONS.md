# Wall functions: compressible, heated (Nichols and Nelson)

Status: criteria stated 5 October 2026, about 23:30, before any wall-function code (TECHNICAL_PLAN *Lightweight engine*, order of work item 3). Results are appended below the line; they do not change the criteria.

## What was read

- **The method in full**, from its author: R. H. Nichols, *Turbulence Models and Their Application to Complex Flows*, Revision 4.01 (University of Alabama at Birmingham; hosted by NASA), Chapter 10, "Wall Function Boundary Conditions", pp. 10-1 to 10-20. It is the Nichols and Nelson formulation and cites the article as its Ref. 6. Read: theory, eqs. 10.1 to 10.22, the procedure, the four test cases and the application hints.
- **Not read: the journal article itself** (Nichols and Nelson, AIAA J. 42(6), 1107-1114, 2004). AIAA returns HTTP 403, and no open copy was found.
- **A published implementation, read as a cross-check:** SU2's `CNSSolver::SetTau_Wall_WF` (`SU2_CFD/src/solvers/CNSSolver.cpp`, master, 5 October 2026). Its comments cite the article's page 1108 for Gamma, beta, Q and phi, and "eq. 11" for Crocco-Busemann. These match the chapter. SU2 issue #2955 (overwriting the wall node's temperature) is a lesson taken below.
- Earlier, read in full: Menter, Carregal Ferreira, Esch and Konno, IGTC2003-TS-059. The omega blend of eqs. 10.17 to 10.19 is the one Nichols takes from Veiser, Esch and Menter.

**The method** (chapter 10):
- Six assumptions in the wall layer: analytic velocity and temperature profiles; analytic turbulence values at the first point; constant pressure; constant shear stress; constant heat flux; frozen chemistry.
- **Velocity.** Spalding's single formula over the sublayer, buffer and log layer (eq. 10.3), with its log-law term replaced by White and Christoph's compressible, heated outer law (eqs. 10.4, 10.5, 10.7). This gives eq. 10.6:
  - y+ = u+ + y+_White − e^(−κB) [1 + κu+ + (κu+)²/2 + (κu+)³/6];
  - y+_White = exp{(κ/√Γ)[asin((2Γu+ − β)/Q) − φ]} e^(−κB);
  - Γ = r u_τ² / (2 c_p T_w), the compressibility;
  - β = q_w μ_w / (ρ_w T_w k_w u_τ), the heat transfer;
  - Q = (β² + 4Γ)^(1/2), φ = asin(−β/Q), r = Pr^(1/3).
- **Temperature.** Crocco-Busemann, eq. 10.8: T = T_w (1 + βu+ − Γu+²).
- **Wall values.** The wall pressure is the first point's. The wall density follows from the equation of state at T_w.
  - An isothermal wall: eqs. 10.6 and 10.8 are solved together for u_τ and q_w.
  - An adiabatic wall: β = 0, and T_w comes from eq. 10.8.
- **Into the scheme.** The wall shear stress and heat flux replace the viscous flux at the wall face (cell-centred), along the first point's tangential velocity and the wall normal (Sondak and Pletcher; or, equivalently, a scaling of the computed flux).
- **Turbulence at the first point.** μ_t from the constant-stress relation (eqs. 10.12, 10.13). For k-omega:
  - ω = (ω_i² + ω_o²)^(1/2), with ω_i = 6 μ_w / (0.075 ρ_w y²) and ω_o = u_τ / (√C_μ κ y);
  - k = ω μ_t / ρ.

  These values are prescribed. The turbulence equations are not solved in the first cell (hint 3).
- **The author's limits.**
  - Hint 1: "The first point off the wall should not exceed y+=100. A good rule of thumb is to place the first point off the wall at y+=50."
  - The axisymmetric bump results "start to diverge for y+=200".
  - The supersonic nozzle with a cooled wall (SST): heat transfer "adequate for many applications even for an initial wall spacing of y+=100".
  - Hint 2: "Wall functions should not generally be used for calculations requiring highly accurate friction drag or heat transfer."

**Two printed equations disagree with their own derivation.**
- **Eq. 10.13** prints dy+_White/du+ = 2 y+_White (κ√Γ/Q) [1 − (2Γu+ − β)²/Q²]^(+1/2). Differentiating eq. 10.7 gives the exponent −1/2, since d asin(x)/dx = (1 − x²)^(−1/2).
  - A central finite difference of eq. 10.7 agrees with the −1/2 form to 1e-9 (measured, a scratch script, 5 October).
  - The two forms agree only where (2Γu+ − β)/Q is small. In the incompressible adiabatic limit both are exact. At Γ = 1.4e-4, β = 0 the printed form is 0.5% to 13% low over u+ 5 to 30.
  - Under a cold wall with combustion gas, the printed form is 3 to 6 times low over u+ 5 to 30 (β = 0.064, Γ = 1.4e-4, guessed inputs).
  - SU2 implements the printed (+1/2) form. The error enters the first cell's μ_t (eq. 10.12), and through it k (eq. 10.20). It does not enter u_τ.
- **Eq. 10.9** prints the adiabatic wall temperature as T = T_w (1 + (γ − 1) r u² / 2), which is not dimensionally consistent. Eq. 10.8 with β = 0 gives T_w = T_1 + r u_1² / (2 c_p), the recovery temperature. That form is used.

## Design, declared before code

1. **One formula at every y+.** Spalding's form is valid down to the wall, where u+ → y+ and τ_w → μ_w u_1 / y_1. So there is no switch between a wall function and a resolved wall: one law on every no-slip wall (the side wall and the injector plate rings), whatever the first cell's size. Criterion 1(a) tests that the formula reproduces the resolved wall where the grid resolves it.
2. **The solve (isothermal walls, which C1 and RL10 use).** For a trial u_τ, eq. 10.8 at the first cell's T_1 gives β explicitly: β = (T_1/T_w − 1 + Γu+²) / u+. The residual of eq. 10.6 at the first cell's y+ is then one scalar equation in u_τ. It is solved by a bracketed Newton-bisection to 1e-13 relative. q_w follows from β.
   - Adiabatic walls: β = 0, and T_w = T_1 + r u_1² / (2 c_p) inside the iteration.
   - The Γ → 0 limit is evaluated without cancellation (a series in √Γ).
3. **Properties.**
   - The first cell's composition is frozen through the layer (assumption 6).
   - ρ_w comes from p_1 and T_w; μ_w, k_w and Pr come from the core's mixture-averaged transport at T_w. r = Pr_w^(1/3).
   - c_p is the mean over the layer, (h(T_1) − h(T_w)) / (T_1 − T_w). Crocco-Busemann is an enthalpy relation, so this choice makes it exact at both end points for a thermally perfect mixture. Its error inside the layer is what criterion 2 measures.
4. **The derivative.** Eq. 10.13 is used with the derived exponent −1/2 (above), checked by finite difference in criterion 0.
5. **Fluxes.**
   - The wall face's viscous momentum flux is τ_w along minus the first cell's tangential velocity. Its heat flux is q_w (C1 and RL10 walls are isothermal).
   - The wall passes no species, k or omega.
   - T_w and ρ_w live only inside the wall function and are never written into a cell's state (the lesson of SU2 issue #2955).
6. **Turbulence in the first cell.** k and omega are prescribed from eqs. 10.12, 10.13 and 10.17 to 10.20 after every stage. The cell's total energy is unchanged, so its internal energy takes the change of rho k. That change is booked in the turbulent-energy ledger beside the clip (`Budgets::clippedTurbulentEnergy`), and the budgets still close.
7. **κ and B are SST's own.** They are read from the resolved reference profile of criterion 1, a least-squares line of u+ against ln y+ over 50 < y+ < 0.1 Re_τ. They are recorded here before any wall-function run. Nichols' 0.4 and 5.5 are reported beside them. The wall function must reproduce the turbulence model's own log layer, or the answer drifts with y+. These constants come from the model, not from engine data.
8. **Threads.** The wall faces are computed alone into the face arrays like every other face (`THREAD_POOL.md`). The result must not depend on the thread count.

Test infrastructure this needs (none of it changes the criteria):
- `tools/sst_pipe_1d.py`: a uniform volumetric heating option, and `tools/n2_properties.csv` regenerated from Cantera up to 3,500 K (today it stops at 600 K);
- the engine: a uniform volumetric heating hook for the pipe, beside `setBodyForce`;
- `crucible_chamber_study`: a turbulent mode with Table A chemistry (SST, no PaSR), and a wall-function switch;
- a wall-function test program for criterion 0.

## Criteria

**0. The algebra and the scheme** (a test program; exact or analytic references).
- (a) **Incompressible adiabatic limit.** With β = 0 and Γ → 0, eq. 10.6 equals Spalding's eq. 10.3, and eq. 10.12 equals eq. 10.11, to 1e-12 relative for u+ in [0, 40]. y+(u+) is continuous as Γ goes from 0 to 1e-2.
- (b) **Inversion.** For first-cell y+ in {0.5, 1, 5, 11, 30, 100, 300, 1,000, 3,000}, Γ in {0, 1e-5, 1e-4, 1e-3} and β in {−0.05, 0, 0.02, 0.05, 0.1}, the solved u_τ reproduces the input (u_1, y_1) pair to 1e-12 relative. For an isothermal wall, the solved q_w reproduces T_1 through eq. 10.8 to 1e-12 relative.
- (c) **Derivative.** dy+/du+ from the corrected eq. 10.13 agrees with a central finite difference of eq. 10.6 to 1e-7 relative over the cases of (b). The printed form's error over the same cases is reported.
- (d) **Threads and budgets.** A turbulent chamber with wall functions on both walls takes 200 steps on 1 and on 4 threads, and `memcmp` finds the states and budgets identical (as `core_tests` *Threads*). The mass and energy budgets close to CHAMBER_C1's 1e-11, with the wall heat flow and the first-cell rho k changes booked.

**1. Fully developed pipe at Re_τ about 10,000** (N2, wall 300 K, p0 1 MPa, R 5 mm; the axial force sized by the reference for Re_τ 10,000; bulk Mach about 0.2, derived). Re_τ 10,000 is about the RL10 throat's in wall units (derived: a guessed throat boundary layer of 2 to 4 mm over the viscous length of 0.24 to 0.4 um in TECHNICAL_PLAN *Lightweight engine* gives 5,000 to 17,000). The reference is `tools/sst_pipe_1d.py`, wall-resolved, at the Richardson limit over its finest grids. The engine marches from the reference profile on 4 axial columns, as check 4 of `TURBULENCE_C2.md` does, until u_b changes by less than 1e-5 per ms.
- (a) **The resolved limit.** On 64 rings with stretching 3.5 (first-cell y+ about 1), the engine with wall functions against the engine without them: u_b within 0.5% and T_axis − T_wall within 1%.
- (b) **First-cell y+ 30, 100, 300 and 1,000.** Uniform rings, nr = Re_τ / (2 y+): 167, 50, 17 and 5 rings. Pass at each y+: c_f within 5% of the reference limit (u_b within about 2.5%).
- Reported: T_axis − T_wall (viscous heating only here; criterion 2 tests heat transfer), and the first-cell y+ the engine computes from its own u_τ.

**2. Heated compressible pipe at the RL10's density ratio** (N2, wall 600 K, p0 3 MPa, R 5 mm, a uniform volumetric heating sized by the reference for T_axis about 3,000 K, so a density ratio of about 5; Re_τ in wall units about 10,000). N2 does not dissociate measurably below 3,500 K at this pressure, so the frozen-chemistry assumption is exact here. This criterion isolates compressibility, heat transfer and variable properties. The reference and the march are as in criterion 1. In this fully developed case the force fixes τ_w and the heating fixes q_w, so the model's outputs are u_b and the temperature difference.
- (a) **The resolved limit**, as 1(a): u_b within 0.5% and the Stanton number St = q_w / (ρ_b c_p,b u_b (T_b − T_w)) within 1%.
- (b) **First-cell y+ 30, 100, 300 and 1,000**, as 1(b). Pass at each y+: c_f within 5% and St within 10% of the reference limit.
- Why 5% and 10%. The compressible law's own error against experiment is 12.5% adiabatic and 18.5% heated (De Chant and Tattar, CR-191185). The wall function's departure from the resolved model it stands in for must be clearly smaller. A 5% friction error moves C1's Isp by about 0.05%, and a 10% heat error by about 0.1% (derived from a guessed 1% friction loss and a guessed 1% heat loss in Isp).
- Reported: the same runs with the printed eq. 10.13. This measures the printed form's effect, whose cause criterion 0(c) has already established.

**The largest tested y+ at which criteria 1(b) and 2(b) both pass is y+_max.** It sets the wall strip's refinement level in AMR (TECHNICAL_PLAN *Adaptive resolution*). The plan's design point is y+ 1,000 to 2,000 at the throat. The author's own guidance stops at 100, and his bump results diverge at 200. If y+_max is below 300, the wall-strip budget in *Lightweight engine* is re-derived from the measured y+_max before AMR is built.

**3. C1, turbulent, in its own geometry** (64 columns, to 8 ms; Table A chemistry with SST, no PaSR, so finite-rate light-off is not tested here; isothermal walls as in CHAMBER_C1). The reference is the resolved wall on 32 rings, with the stretching chosen so the largest first-cell y+ is 1 or less (the study's laminar estimate). The chemistry is equilibrium, so the resolved layer recombines as it cools. The frozen-layer assumption is measured here.
- (a) **The resolved limit.** Wall functions on the reference's rings against the reference: c*, vacuum Isp and vacuum thrust within 0.02%, wall heat flow and wall axial force within 2%.
- (b) **The production grid.** Wall functions on 12 uniform rings (first-cell y+ about 80 today) against the reference. Pass:
  - c*, vacuum Isp and vacuum thrust within 0.2%. That is half of the 64x12 grid error in c* (+0.40%, CHAMBER_C1), so the wall model's effect is smaller than the discretisation's.
  - The wall heat flow within 10%, and the wall axial force within 5%. Both are averages over the last 1 ms.
- Reported: the same 12-ring grid with today's resolved-wall rule (at y+ about 80), so the change the wall function makes is measured. Also the first-cell y+ distribution and the time to 8 ms.

**4. Cost.** On the 12-ring C1 of criterion 3(b), the wall-function step against today's resolved-wall step on the same grid. Pass: 10% or less extra cost per cell update. Reported: the wall time of 3(b) against the reference's, the payoff.

Every run labels its wall treatment. The wall heat flux and friction become model outputs. Their error bands in the RL10 loss ledger come from criteria 2 and 3.

### Restated 6 October 2026 (committed 00:55 CDT, 5c3b747), before any code or rerun: Ben's decision on the hot wall

*Times in this section were corrected at 01:10 CDT. The first stamps (01:05, 01:30, 02:20) were estimates, not read from the clock, and they were late. The commit times are the record.*

**The decision** (Ben, through the Director, 6 October, between 00:43 and 00:55 CDT): Option B. The RL10 is an expander engine, so wall heat is the quantity that matters, and DNS beats standard SST; the correction is published, so it is not a bespoke fix. The correction goes into the engine and into the reference alike. Criterion 2 is restated against the corrected SST before it is rerun. The wall model is tabulated offline. A measured hot-wall heat-flux dataset is added, so the truth is not only model against model. The change of SST form is recorded in TECHNICAL_PLAN. The results of criterion 2 below (against standard SST) stay as the record of why.

**The model, in the engine and the reference alike: SSTs with the corrections of Hasan, Elias, Menter and Pecnik** (J. Fluid Mech. 1019, A8, 2025; arXiv 2410.14637; the authors' solver is github.com/Fluid-Dynamics-Of-Energy-Systems-Team/RANS_Scaling2024). In their full form (their sections 3 to 5, the 2-D test of section 7 and appendix C):
- ψ = √ρ/μ, ℓ the wall distance, n the unit wall normal, S_n = 1 / (ψ + ℓ n·∇ψ); μ_k = μ + σ_k μ_t and μ_ω = μ + σ_ω μ_t, σ blended by F1 as in SSTs.
- k gains Φ_k = F1 Φ_k^in + (1 − F1) Φ_k^out, with Φ_k^in = (S_n/μ) ∇·[μ_k (S_n/μ) ∇(ρk)] − ∇·[μ_k ∇k] and Φ_k^out = (1/√ρ) ∇·[μ_k (1/√ρ) ∇(ρk)] − ∇·[μ_k ∇k]. The blend leaves the total k diffusion as F1 times the inner form plus (1 − F1) times the outer form.
- ω gains Φ_ω (the same blend) and Φ_CD, with Φ_ω^in = (ρ S_n/μ²) ∇·[μ_ω (S_n/μ) ∇(μω)] − ∇·[μ_ω ∇ω], Φ_ω^out = ∇·[μ_ω (1/√ρ) ∇(√ρ ω)] − ∇·[μ_ω ∇ω], and Φ_CD = 2(1 − F1) σ_ω2 / (√ρ ω) ∇(ρk)·∇(√ρ ω) − 2(1 − F1) ρ σ_ω2 / ω ∇k·∇ω. F1's own argument keeps the conventional CD_kω (the paper changes the transport term only).
- μ_t is multiplied by D^ic = D(R_t, M_t) / D(R_t, 0), D = [1 − exp(−R_t / (3.5 + 0.39 M_t^0.77))]², R_t = ρk/(μω), M_t = √(2k)/a.
- **Energy.** The paper's 2-D energy equation (their eq. 7.1) carries the corrected k diffusion, so Φ_k enters the total energy too. Here E includes ρk, so the correction leaves the internal energy unchanged. Φ_k is not a divergence: ∫Φ_k dV is booked in its own energy ledger, and the energy budget must still close to round-off with it.
- **Discretisation.** Face coefficients (μ_k S_n/μ, μ_k/√ρ, μ_ω S_n/μ, μ_ω/√ρ) are the two-cell means, as the conventional coefficients are; the prefactors (S_n/μ, 1/√ρ, ρS_n/μ²) are the cell's own. The corrected diffusion is computed in full and the conventional one subtracted, so the corrected terms enter as cell sources. n is the unit vector from the nearest wall point, stored with the wall distance. n·∇ψ comes from the cell gradient. In the 1-D reference, ℓ = R − r, n = −r̂ and ∇· = (1/r) d/dr (r ·).
- **One guard, chosen here (the paper gives none).** S_n is singular where ℓψ has a maximum along n. The denominator of S_n is floored at ψ/10. Every run reports the F1-weighted fraction of cells where the floor acts.
- **Not taken.** The wall dissipation ε_w (their Φ_e,2): their own 2-D test leaves it out, because with a fixed wall temperature it adds to the wall heat flux (their section 7). Pr_t stays 0.9, their choice for boundary layers and high-Mach channels.

**Checks on the reference, before it is used** (`tools/sst_pipe_1d.py`; a failure stops the reference):
- **R1, the operator.** With constant ρ and μ, every Φ is zero to 1e-12 of the conventional diffusion. On a smooth manufactured profile (ρ, μ, k, ω and μ_t analytic in r, a heated-pipe shape), the tool's discrete Φ_k, Φ_ω and Φ_CD converge to a fourth-order finite-difference evaluation of the same expressions on 100,000 points, with observed order 2 ± 0.25 over three grids.
  - *The norm, made precise 6 October about 01:00 CDT, before R1 ran.* The tool takes the wall face's gradient from the wall value and the last centroid. So the last cell's truncation error is of order one for any diffusion operator, the conventional one included. Derived for a uniform planar grid: there the discrete operator is 0.75 of the exact one. Every other cell's error is second order. The order is therefore judged on the volume-weighted L1 error over all cells but the wall cell, relative to the reference's L1 norm over the same cells. Reported beside it: the wall cell's error, the largest error, and the error in ∫Φ_k dV over the whole pipe, with the conventional operator's own wall-cell error for comparison. The reference is not a 100,000-point grid: three nested derivatives on a step of 1e-5 R would lose about 1e-6 to round-off. It is a nested five-point (fourth-order) stencil at each centroid, with step h = min(R/1000, r/10).
- **R2, the authors' code.** Their published solver, unchanged except for imports that newer SciPy renamed, on their low-Mach channel of Pecnik and Patel (2017): Re_τ 950, uniform heating 75 μ_w c_p T_w / (Pr_w h²), μ ∝ T^0.7, Pr ∝ μ, Pr_t 1, the inner corrections throughout and no Φ_CD, as in their paper. The tool runs the same configuration (planar, the same property laws, the inner form throughout, their P_k limit of 20 β*ρkω). Both are grid-converged: each value changes by less than 0.1% on doubling the grid. Pass: the centreline u+ and (T_c − T_w)/T_w within 0.5%.
  - *What the first runs showed, and R2 restated (6 October, 01:14 CDT, before the longer runs).*
    - Their solver converges at first order in the grid, as the tool does. Their u_c+ is 41.340, 40.275 and 39.760 at 2N 120, 240 and 480 (the last doubling changes it by 1.3%). The tool's is first order too, observed order 1.04 from 800 to 3,200 rings.
    - At 2N 960 their Picard loop stops at its 100,000-iteration cap with its change still above its 1e-7 tolerance. That value (39.320) is not converged and is not used.
    - So "less than 0.1% on doubling" cannot be met by either code at a grid it can afford. Restated:
      - (i) One change to their solver: the iteration cap is raised from 1e5 to 2e6. Their tolerance stays 1e-7. Their 2N 960 and 1,920 are rerun with it.
      - (ii) Each code's value is its Richardson extrapolation over its last three grids. It is accepted when the observed order is between 0.7 and 1.3 and the extrapolation moves by less than 0.1% from the previous triplet's.
      - (iii) Iteration error: their 2N 480 is rerun with tolerance 1e-8. The change is reported and must be under 0.05%; if not, every grid is rerun at the tighter tolerance.
    - The pass is unchanged: u_c+ and (T_c − T_w)/T_w within 0.5%.
    - *(iii) measured, 01:16 CDT.* Their 2N 480 at tolerance 1e-8 against 1e-7: u_c+ 39.78010 against 39.76012 (+0.050%), (T_c − T_w)/T_w 3.85251 against 3.84729 (+0.136%). That is over 0.05%, so every grid of theirs is rerun at 1e-8 (2N 240, 480, 960 and 1,920). Their 2N 480 at 1e-9 is reported as the check on 1e-8.
- **R3, constant properties.** Criterion 1's cold pipe with the corrections: the Richardson c_f within 0.5% of the uncorrected 0.00337113. The corrections vanish for constant properties. What is left is the pipe's small viscous heating and D^ic at a turbulent Mach number of about 0.02 (derived).
  - *What the first two grids showed, and R3 restated (6 October about 01:07 CDT, before the grid sequence ran).*
    - On 200 and 400 rings, the corrected cold pipe's c_f is 0.96% below the uncorrected one on the same grid. Without D^ic the gap is 0.16%. So the variable-property terms do vanish, apart from the pipe's small heating; D^ic does not.
    - D^ic lowers c_f by about 0.8% at M_τ 0.009. The 0.5% band rested on my own derivation that D^ic is negligible at M_t 0.02, and that derivation was wrong. f(M_t) = 0.39 M_t^0.77 has an infinite slope at zero. The paper's calibration (C − 5.2 = 7.18 M_τ, their appendix B) gives a log-law shift of 0.065 here, about 0.5% in c_f (derived).
    - So R3 as stated fails on D^ic, which is the model working as published, not an error in the tool. Restated:
      - (i) without D^ic, the Richardson c_f within 0.5% of 0.00337113 (the check R3 was meant to be);
      - (ii) with D^ic, the Richardson c_f and the equivalent shift in u_b+ are reported beside the paper's 7.18 M_τ, not judged;
      - (iii) D^ic itself: unity at M_t = 0, and equal to a hand evaluation at R_t 3.5, M_t 0.1 (0.978300). Measured: 1 exactly at R_t 0.5, 3.5 and 20; 0.9783002707 against 0.9783002707 (4e-16).
- Reported, not judged: the corrected heated profile's u+ under semi-local (Trettel-Larsson) scaling against the cold one at y* 30, 100 and 300. This is the collapse the correction is built for.

**κ, B and A+ are re-read from the corrected cold reference** by design item 7's rule, before any wall-model run. If they move by less than the fit's own range (κ 0.364 to 0.375), the present values stay.

**Criterion 2, restated.** The same case and the same bands. N2, wall 600 K, p0 3 MPa, R 5 mm, force 19,800 N/m³ and heating 1.95e9 W/m³. These are held fixed, so Re_τ and T_axis move with the model and are reported. The reference is the corrected tool, 200 to 6,400 rings, at its Richardson limit.
- (a) As before: u_b within 0.5% and St within 1%.
- (b) As before: c_f within 5% and St within 10% at first-cell y+ 30, 100, 300 and 1,000.
- **The wall model under test** is the Kawai-Larsson equilibrium ODE (`tools/eq_wall_model.py`: ρ, μ, λ and c_p at the local temperature, mixing length with semi-local van Driest damping, Pr_t 0.9). Nichols and Nelson's law stays in the code only as the comparison.
- **A priori first.** The ODE is given the corrected heated reference's (u, T, p) at y+ 1, 5, 30, 100, 300, 1,000 and 3,000. Pass from y+ 30 to 1,000: τ_w within 2.5% and q_w within 5%, half the bands of (b). Only then is the engine run.
- **The fallback, stated now.** If the a priori check fails, the wall model becomes the corrected SST's own inner layer. That is the 1-D constant-stress, constant-heat-flux layer with the corrected k and ω equations, as in the paper's section 4.1 solver, tabulated the same way, with the a priori check repeated. Nothing else is tried without a new statement.
- **The table, offline.** The wall model is tabulated on backhouse and the engine reads only the table. For N2 the inputs are Re_1 = ρ_w u_1 y_1/μ_w, T_1/T_w, T_w and u_1²/(c_p,w T_w); p drops out of an ideal gas. The outputs are y_1+ and B_q = q_w/(ρ_w c_p,w u_τ T_w). **Criterion 0 (e):** against direct solves at 2,000 random points inside the table's range, y_1+ and q_w within 0.5% at the 99th percentile and within 1% at most.
- For C1 the layer is a mixture that recombines as it cools. The composition axes are a Table B design item. Criterion 3 runs once that is built.

**Criterion 3, the shorter reference.** The explicit 32-ring run to 8 ms is estimated at about 60 h (below). Instead:
- **Start.** The resolved reference (64 columns, 32 rings stretched so the largest first-cell y+ is 1 or less, the corrected SST to the wall) starts from the settled 64x12 wall-function solution at 8 ms, interpolated. Each column is mapped by η = r / R_wall(z). Primitive values (ρ, u_z, u_r, T, the mass fractions, k and ω) are linear in η between the old centroids. Between the wall and the first old centroid the interpolation runs to the wall values: u = 0, T = T_w, k = 0, and ω from the wall rule.
- **The settling test.** History is written every 2 µs and grouped into 0.05 ms windows. For each window: the means of the wall heat flow and the wall axial force, and c*, vacuum Isp and vacuum thrust at the window's end. A quantity is settled when its last three window values x1, x2 and x3 contract, with ratio ρ = (x3 − x2)/(x2 − x1) between 0 and 0.9. Its geometric remainder |x3 − x2| ρ/(1 − ρ) must also be at most a quarter of its tightest band in criterion 3. That is 0.005% for c*, Isp and thrust (a quarter of 3(a)'s 0.02%) and 0.5% for the wall heat flow and force (a quarter of 2%). Without contraction it is not settled. The run stops when all five are settled. Reported values are means over the last 0.1 ms (two windows), in place of the last 1 ms.
- **Like with like.** The wall-function runs of 3(a) (on the reference's rings) and 3(b) (12 rings) use the corrected SST, the same start, the same test and the same averaging.
- **Cost first.** The first 0.1 ms of the reference is timed on backhouse. The projection takes the settling time as 1 ms (guessed: about five chamber residence times of about 0.22 ms, derived). If the projection is over about 8 h, the columns are cut to 32 for the reference and for 3(a) and 3(b) alike. The cost of that in coverage: the axial resolution is halved; it is not the production 64-column grid, so 3(b) then measures the wall model on a grid that is not the production one; the comparison stays like with like. The timing can use the present engine (standard SSTs), since the correction changes the step's cost but not dt.

**The measured check (Ben's item).** A measured hot-wall heat flux (a heated tube, or calorimetry in a rocket chamber) is being found. The search reports the inputs only. The case, its inputs, its quantity and its band are stated here before any measured value is read, as for the RL10 validation.

---

## Results

### Criterion 0 (a) to (c): pass (6 October 2026, about 00:20)

The law is `crucible::wallLaw` in `core/walls.cpp`. The test program is `tests/wall_function_tests.cpp` (ctest `wall_function_verification`, 0.2 s). Its output is in `wall_functions_criterion0_2026-10-06.txt`. κ 0.4 and B 5.5 were used; the algebra does not depend on them.

| Check | Measured | Limit |
|---|---|---|
| (a) eq. 10.6 against Spalding's eq. 10.3, u+ 0 to 40 | 2.1e-16 | 1e-12 |
| (a) eq. 10.12 against eq. 10.11, u+ 0 to 40 | 1.9e-15 | 1e-12 |
| (a) continuity: u_eq+ against an independent quadrature, 77,564 points, Γ 0 to 1e-2 | 1.7e-15 | 1e-13 |
| (a) continuity: largest decrease of y+ as Γ grows, over its rounding allowance | 0.54 | 1 |
| (b) solved u_τ against the constructed one (180 cases, at most 16 Newton iterations) | 5.9e-16 | 1e-12 |
| (b) eq. 10.6 at the solved u+, Γ, β | 2.5e-14 | 1e-12 |
| (b) T_1 back from the solved q_w through eq. 10.8 | 2.8e-14 | 1e-12 |
| (c) corrected eq. 10.13 against a five-point central difference | 4.2e-10 | 1e-7 |
| (c) reported: printed eq. 10.13 against the same difference | 1.10 (largest, at y+ 30, Γ 0, β 0.1) | none |

**How the law is evaluated.**
- u_eq+ (the integral of eq. 10.5) is evaluated in one closed form, the arctangent of an angle difference. It reduces to (2/β)(√(1 + βu+) − 1) at Γ = 0 and to u+ at Γ = β = 0. There is no series in √Γ and no separate branch at small Γ: the closed form has no cancellation there.
- The bracket of eq. 10.6, exp(κu_eq+) − (1 + x + x²/2 + x³/6) with x = κu+, is computed in whichever of two forms rounds less:
  - e^x expm1(κu_eq+ − x) plus the series of the exponential's tail;
  - or the plain difference.
- The isothermal solve uses u_eq+(u+; Γ, β) = u+ · u_eq+(1; G, θ) on its path, where G = r u_1²/(2 c_p T_w) and θ = T_1/T_w − 1 + G are fixed by the first cell. So one quadrature-free constant serves every iteration. The solve is Newton on u+ y+(u+) = Re_1 inside a bracket doubled up from u+ = 1, with bisection when a step leaves the bracket.

**What failed first, and why, stated plainly.** The criteria did not change. The test's implementation of two of them did, and the first attempt's code had two bugs.
- **(a) Continuity was first tested by three sub-checks that were wrongly posed. All three failed, and the failures were the exact law's own behaviour, not faults in the code.**
  - y+ at Γ ≤ 1e-14 against Γ = 0, to 1e-12: measured 4.3e-11. The law's real first-order sensitivity is κ e^(−κB) e^(κu+) Γ u+³/6 relative to y+. At u+ 40 that is about 4,300 Γ, so a 1e-12 bound below Γ = 1e-14 cannot hold.
  - The slope (y+(Γ = 1e-10) − y+(0)) / 1e-10 against the analytic slope, to 1e-6: measured 1.7e-4. At small u+ the change in y+ is below one unit in the last place of y+ (about u+), so the difference quotient measures rounding.
  - y+ must not decrease as Γ grows, to 1e-15 relative: measured 8.1e-15. exp(κu_eq+) amplifies a relative rounding of u_eq+ by κu_eq+, about 16 at u+ 40.

  They were replaced by two checks that test continuity directly. (1) Γ enters eq. 10.6 only through u_eq+, so u_eq+ is checked against an independent Gauss-Legendre quadrature of its integral (64 panels of 20 points) on a ladder of 2,804 values of Γ (0, 1e-300, 1e-100, then 10^(e/100) from 1e-30 to 1e-2). (2) y+ may not decrease along that ladder by more than 1e-15 (1 + κu_eq+) relative, its rounding.
- **(a) eq. 10.12 against eq. 10.11: first measured 3.6e-9.** The test's reference for eq. 10.11 subtracted the polynomial from the exponential at small x and cancelled. The reference now sums the exponential's tail as a series below x = 1: 1.9e-15.
- **(b) Two bugs in the solve, both fixed.**
  - A false root. When heating makes u_eq+ < u+, the bracket subtracted two terms of size e^(κu+). At u+ 300 every digit was lost, so the solve stopped at a spurious root at u+ 339 (u_τ wrong by 89%). Eighteen cases were wrongly listed as having no layer. Fix: the two-form evaluation above, and a bracket that grows up from u+ = 1.
  - A Newton safeguard. It required a step strictly inside the bracket, so a converged step that landed on the bracket's end was rejected, and the solve bisected away from the root (u_τ errors near 6e-14, up to 175 iterations). Fix: a step below 1e-13 u+ ends the solve.
- **(c) Printed eq. 10.13.** Its error reaches 110% (its derivative has the wrong sign). Where Γ = 0 and β ≠ 0, (2Γu+ − β)²/Q² = 1, so the printed bracket [1 − (·)²]^(+1/2) is zero and the log-law term vanishes from dy+/du+. With the derived exponent −1/2 the factor cancels: since Q² − (2Γu+ − β)² = 4ΓD², the derived form is κ y+_White / D, with D = (1 + βu+ − Γu+²)^(1/2), and that is what the code evaluates. The printed form equals 4κΓD y+_White / Q², which is zero at Γ = 0. Its smaller errors at β = 0 are in *What was read* above.

**A limit of the law found here: strong wall cooling** (derived from eq. 10.6 on the isothermal path, Γ = 0; a scratch scan, 6 October).
- Along the isothermal solve's path, y+(u+) = u+ + e^(−κB)[exp(κu+ I) − P_4(κu+)], with I = u_eq+(1; G, θ) below 1 when T_1 > T_w.
- Its slope dy+/du+ first becomes negative near T_1/T_w = 11 (I = 0.46, at u+ 20). Above this the law's velocity profile folds back, and no single y+ belongs to some u+.
- The solve's own function u+ y+(u+) stays increasing until T_1/T_w ≈ 13.5. Above that, a first-cell Reynolds number can have more than one root.
- C1's wall is 600 K against at most about 3,350 K in the gas, so T_1/T_w ≤ 5.6, where the smallest slope is 0.75 or more (derived). An RL10 hot-gas wall at 500 to 800 K gives about 4.5 to 7 (guessed wall temperatures).
- The engine will refuse (throw) when T_1/T_w > 11 at a wall cell, rather than return a doubtful root. Whether that can happen in a cold start (a cold wall under a hot flame) is a question for criterion 3's runs.

Criterion 0 (d) is under *The law in the engine* below.

### κ and B (design item 7), recorded 6 October 2026, about 00:15, before any wall-function run

**The reference for criterion 1** (`tools/sst_pipe_1d.py --radius 5e-3 --p0 1e6 --twall 300 --force 45690 --stretch 3.5`; outputs in `wall_functions/`).
- The force 45,690 N/m³ was sized for Re_τ 10,000 at the wall's properties (derived: u_τ = Re_τ μ_w/(ρ_w R), f = 2ρ_w u_τ²/R). The run gives Re_τ 9,941.
- 100 to 400 rings are not yet in the asymptotic range (observed order 0.72 to 0.82), so the run was taken to 6,400 rings. From 800 rings on, the observed order is 0.99 to 1.01, the first order that the wall-omega rule gives (`SST_REFERENCE.md`).
- Richardson limits from 1,600-3,200-6,400 rings (derived): **c_f 0.00337113, u_b 77.6781 m/s, T_axis 302.37868 K**. The 800-1,600-3,200 triple agrees to 2.6e-5 in c_f and 8e-6 in u_b. Bulk Mach about 0.22 (derived).
- These are criterion 1's reference values.

**The fit** (`tools/log_law_fit.py`, on the 6,400-ring profile; `wall_functions/log_law_fit_re_tau_1e4.txt`). By the stated rule, a least-squares line of u+ against ln y+ over 50 < y+ < 0.1 Re_τ (2,512 points):

**κ = 0.3697, B = 3.752.** These become the defaults. Nichols' 0.4 and 5.5 are reported beside them.

What the fit stands on, measured:
- **The model's log layer is not a straight line at this Re_τ.** The local κ, 1/(du+/d ln y+), rises from 0.33 at y+ 50 to 0.378 at y+ 500 and falls to 0.369 at y+ 2,000. The line misses the profile by at most 0.12 in u+ over the fitted range.
- **The fit depends on its range.** Over 30 to 0.1 Re_τ: κ 0.364, B 3.49. Over 100 to 0.1 Re_τ: κ 0.375, B 3.97. Over 50 to 0.05 and 50 to 0.2 Re_τ: κ 0.366 and 0.371. The 3,200-ring profile gives the same κ to 4 digits and B 0.01 higher.
- **This is pure k-omega with the limiter off.** F1 is 1.0000 from y+ 30 to 3,000, the cross-diffusion is zero, and a_1 ω / S is 1.009 to 1.037, so the Bradshaw limiter is off (a 1,600-ring solve). Its constant-stress log layer has κ² = √β* (β_1/β* − γ_1)/σ_ω1, so κ = 0.408 with γ_1 = 5/9 (derived). The profile approaches that slope slowly.
- **A diagnostic at Re_τ 100,000** (the same u_τ, R 5 cm; 6,400 rings, not converged, observed order 0.82; `wall_functions/log_law_fit_re_tau_1e5_diagnostic.txt`).
  - Out to y+ 500 the inner profile is the same as at 10,000 within 0.7% in u+ (u+ 20.549 against 20.553 at y+ 500; 14.301 against 14.212 at y+ 50).
  - The local κ reaches 0.395 at y+ 1,000 (y/R 0.01), then falls in the outer layer. At Re_τ 10,000, y+ 1,000 is already y/R 0.1, so the stress falling across the pipe caps the local κ at 0.378.
  - The stated rule's fit there gives κ 0.388, B 4.49. The constants therefore depend on how much of the outer layer the range includes. They describe this model's inner layer, not a universal log law.
- **Spalding's formula with the fitted constants against the model's own profile, y+ 1 to 2,000:**
  - at Re_τ 10,000, within 1.0% in u+ (largest −1.0% at y+ 5; 0.3% or less from y+ 10 to 2,000);
  - against the Re_τ 100,000 diagnostic with the same constants (a scratch check), within 1.2% up to y+ 1,000, +1.1% at y+ 2,000 and +1.5% at y+ 5,000. The diagnostic is not grid-converged (its u_b moves 0.4% from 3,200 to 6,400 rings), so its near-wall differences of about 1% are partly grid error.
  - Nichols' 0.4 and 5.5 are 2% to 5% high in u+ from y+ 10 to 500 at both Re_τ, which is 4% to 10% in c_f at a first cell there (derived, c_f ∝ 1/u+²).

So the wall function carries the SST model's own inner layer. With Nichols' constants it would carry a different one, and criterion 1(b)'s 5% band in c_f would be about used up by that alone.

### The law in the engine (6 October 2026, about 00:20 to 01:00; commit e10356a)

**As built** (`core/transport.cpp`, `core/flow.cpp`):
- One law on every no-slip wall face: the side wall (one face per column; the normal is the contour's, so a contoured wall is not radial) and the plate rings outside the supply. The corner cell takes the nearer face, the side face on a tie.
- The wall face carries the law's traction and heat flux and nothing else: no species, k or ω flux. For gradients the wall still carries u = 0 and T_w, and the cell's own k and ω (zero normal gradient).
- The first cell's k and ω are prescribed after each stage at fixed total energy: μ_t = μ_w · max(0, eq. 10.12); ω = hypot(6μ_w/(β_1 ρ_w y²), u_τ/(√β* κ y)); ρk = ω μ_t. The change in ρkV is booked in the measurement `prescribedTurbulentEnergy`, so the budget stays closed. The cell's turbulence derivative and source are zeroed. Its eddy viscosity is ρk/ω without the SST limiter, because eq. 10.12 is already the constant-stress μ_t.
- Layer properties: ρ_w = p_1/(R T_w); μ_w, k_w and Pr_w at T_w; c_p is the layer mean (h(T_1) − h(T_w))/(T_1 − T_w); the recovery factor is Pr_w^(1/3). The engine throws if T_1/T_w > 11 (the fold-back found under criterion 0).
- A known artefact, small and bounded: in gas colder than the wall, eq. 10.12 gives μ_t up to μ_w − μ_1 even as u_1 → 0.
- `Flow::setHeating` adds a uniform volumetric source for criterion 2 (booked in the energy budget). `crucible_wall_function_study` has the modes `threads`, `pipe` and `apriori`. `crucible_chamber_study` takes `CRUCIBLE_WALL_FUNCTIONS=law|printed` and `CRUCIBLE_RADIAL_STRETCHING`.

**Criterion 0 (d): pass** (`wall_functions/wall_functions_criterion0d_2026-10-06.txt`; ctest `wall_function_threads`). A turbulent N2 chamber, 16x6, contoured (5 to 3.5 mm), supply over r < 3 mm so the plate rings including the corner carry the law, wall 400 K, 200 steps:
- 1 and 4 threads: state, partial densities, ρk and ρω, the budgets, the wall heat flow and both turbulent ledgers are bit-identical (memcmp);
- mass budget 3.4e-16, energy 3.2e-16 (limit 1e-11);
- reported: wall heat flow 11.55 W with the law, 6.14 W with the no-slip rule on the same grid; wall axial force 0.6855 N against 0.6847 N; prescribed ρkV 1.5e-6 J.

**Wall functions off: unchanged.** C1 eqtt 64x12 to 4 ms on 4 threads is bit-identical (cmp of the raw state and the history) to the 04cba9b run, 80 s.

**Criterion 1 (b), the cold pipe** (`crucible_wall_function_study pipe`, 4 columns, started from the 6,400-ring reference interpolated to the ring centroids at its bulk density; steady when u_b changes less than 1e-5 per ms). Against c_f 0.00337113 and u_b 77.6781 m/s:

| Rings | First-cell y+ (engine) | Steady at | c_f | c_f error | u_b error | T_axis − T_w | Pass (5%) |
|---|---|---|---|---|---|---|---|
| 5 | 956 | 42.5 ms | 0.00324881 | −3.63% | +1.86% | 1.78 K | yes |
| 17 | 289 | 50.5 ms | 0.00327398 | −2.88% | +1.46% | 1.92 K | yes |
| 50 | 99.0 | 51.0 ms | 0.00328661 | −2.51% | +1.27% | 2.05 K | yes |
| 167 | about 30 | running (12 ms: c_f 0.00331035, −1.80%; not steady) | | | | | |

The reference's T_axis − T_w is 2.38 K. The law alone is within 1% of the reference's wall shear at these y+ (the a priori check below), so most of the 2.5% to 3.6% is the coarse grid's outer flow.

**Criterion 1 (a)** (64 rings, stretching 3.5, the law against the no-slip wall) is running on backhouse (`/home/greff/crucible_wf`, tmux `crucible_wf_res` and `crucible_wf_law`). dt is 5.9e-10 s (1 ms in 1,705,015 steps, measured), so 1 ms takes about 8 minutes on 4 threads and steady state (about 50 ms, as in 1(b)) about 6 to 7 hours (derived). At 1 ms: u_b 78.0212 (law) against 77.9839 m/s (no-slip), +0.05%. Neither is steady (force balance still 6% to 7% off at 1.5 ms).

### Criterion 2: the published law fails against SST's heated pipe (6 October 2026, about 00:30 to 01:00)

**The reference** (`tools/sst_pipe_1d.py --radius 5e-3 --p0 3e6 --twall 600 --force 19800 --heat 1.95e9 --stretch 3.5`, 200 to 6,400 rings; `wall_functions/reference_heated_2026-10-06.txt`). The force and heating were sized on 200 rings for Re_τ 10,000 and T_axis about 3,000 K. The run gives Re_τ 9,950, T_axis 2,941 K (T_axis/T_w 4.9), bulk Mach about 0.03 (derived). Observed order 0.98; Richardson limits from 1,600 to 6,400 rings (derived): **c_f 0.00725704, u_b 28.45659 m/s, T_axis 2,940.790 K, T_b 2,568.923 K, St 0.00394250**.

**A priori, without a march** (`crucible_wall_function_study apriori`; `wall_functions/apriori_nichols_nelson_2026-10-06.txt`). The law, with the engine's property choices, is given the reference profile's (u, y, T, p) at one point and compared with the reference's own wall shear (f R/2) and wall heat flux (from the energy balance):

| y+ | T_1/T_w | cold: τ_w error | heated: τ_w error | heated: q_w error |
|---|---|---|---|---|
| 1 | 1.09 | −0.01% | −3.1% | −3.7% |
| 5 | 1.39 | +1.0% | −12.8% | −14.9% |
| 30 | 2.19 | +0.27% | −24.4% | −25.1% |
| 100 | 2.75 | +0.12% | −20.6% | −19.7% |
| 300 | 3.26 | −0.55% | −17.9% | −16.0% |
| 1,000 | 3.83 | −1.04% | −15.6% | −13.1% |

(T_1/T_w is the heated case's.) The cold pipe is within about 1% everywhere. Criterion 0 verified the law's algebra to 1e-12, so this is the law's form, not the code.

**In the engine** (criterion 2(b), the heated pipe started from the reference profile; against c_f 0.00725704, St 0.00394250, u_b 28.45659 m/s):

| Rings | First-cell y+ | Steady at | c_f | c_f error (band 5%) | St | St error (band 10%) | u_b error |
|---|---|---|---|---|---|---|---|
| 5 | 989 | 98.0 ms | 0.00586155 | **−19.2%** | 0.00334083 | **−15.3%** | +11.3% |
| 17 | 298 | 98.0 ms | 0.00600353 | **−17.3%** | 0.00336784 | **−14.6%** | +9.9% |
| 50 | about 100 | running (20 ms: c_f 0.00614, St 0.00341; not steady) | | | | | |

Criterion 2(b) fails at y+ 300 and 1,000, as the a priori check predicts.

**Why: two separate effects, measured on the reference profile** (`tools/eq_wall_model.py --transforms`; `wall_functions/heated_sst_scaling_2026-10-06.txt`):
1. **The sublayer's viscosity.** The law uses μ_w throughout. At y+ 1 the gas is already 9% hotter than the wall, and by y+ 5 39% hotter, so μ is 6% to 25% higher (derived). That alone is the −3% at y+ 1.
2. **SST's heated buffer layer.** Under van Driest scaling the heated profile keeps the cold log slope (local κ 0.35 to 0.375 from y+ 100 to 1,000, against 0.357 to 0.377 cold), but sits 1.7 lower in u+ from y+ 30 on. Under semi-local (Trettel-Larsson) scaling it sits 2.0 lower, and the shift builds across y* 5 to 100. Nichols and Nelson assume the cold intercept in van Driest coordinates, so they miss by the 1.7, which is about 20% in τ_w at these u+ (derived).

**A published model that handles effect 1 does not fix effect 2.** The equilibrium ODE wall model (Kawai and Larsson, Phys. Fluids 24, 015105, 2012) integrates the constant-stress, constant-heat-flux layer with ρ(T), μ(T), λ(T) and c_p(T), and a mixing length with semi-local van Driest damping. With κ 0.3697 and A+ 14.0 fitted to the cold profile, a priori (`wall_functions/apriori_ode_model_2026-10-06.txt`):
- y+ 1: within 4e-4 in both pipes, so effect 1 is gone; heated y+ 5 and 11 are within 2.5%;
- cold, y+ 5 to 3,000: τ_w within 2.1%;
- heated, y+ 30 to 3,000: τ_w −17% to −25%, q_w −16% to −24%. No better than Nichols and Nelson.

**What the literature says.** DNS of variable-property channels collapses onto the cold law under semi-local scaling (Patel et al. 2015; Trettel and Larsson 2016). Standard SST does not follow it. Hasan, Elias, Menter and Pecnik (J. Fluid Mech. 1019, A8, 2025; arXiv 2410.14637; read in full 6 October) compare SST with DNS.
- **39 boundary layers.** The uncorrected model's errors reach about 15% in velocity and 40% in temperature. The Catris-Aupoix/Otero corrections leave them similar or slightly worse. Their corrections bring both within about 5% and 15%.
- **11 channels.** Uncorrected: up to about 25% and 30%. Corrected: within about 10% and 15%.
- **A 2-D plate at Mach 10.9.** c_f and c_h come closer to DNS than the uncorrected model's. They are shown against a ±5% margin; no number is stated.

The correction changes the k and ω diffusion so that the model follows semi-local scaling (after Pecnik and Patel, JFM 823 R1, 2017), and it adds a damping for intrinsic compressibility. So the 20% gap measured here is a known property of standard SST under strong heating. Both wall laws follow scalings that DNS supports; the reference does not.

*Corrected 6 October (committed 00:55 CDT).* This paragraph first said "Hasan and Pecnik ... errors up to 23% in velocity and 29% in temperature ... brings this to 3% and 8%". Those numbers are not in the paper; I had written them before reading it. The milestone reply sent just before carried them too.

**The question this raises, for the Director and Ben.** Criterion 2's reference is standard SST resolved to the wall. A wall function can match it only by carrying SST's own heated buffer layer (for example a 1-D SST sub-grid in the first cell). That would carry SST's known heated-wall error, about 20% here, into the hot-gas wall heat flux. The alternative is semi-local physics in both places: a variable-property wall model (the ODE model above, or Nichols and Nelson with a semi-local correction) and the Hasan-Pecnik correction in the engine's SST. Criterion 2's reference would then be the corrected SST, with its published DNS errors as the reference's own band. Nothing is changed until this is decided.

### Criterion 3: the stated reference cannot be run explicitly

Criterion 3's reference is 32 rings with the largest first-cell y+ at 1 or less. On C1 that puts the wall ring at about 0.35 µm and dt at about 1.2e-10 s (estimated tonight from the 12-ring run's first-cell y+, not measured; for scale, the 1(a) pipe's measured dt is 5.9e-10 s at a ring of about 1 µm in 300 K N2). That is about 7e7 steps to 8 ms over 2,048 cells. At C1's measured 1.5 µs per cell update on 4 threads (eqtt 64x12 to 4 ms in 79 s, 04cba9b), it would take about 60 hours (derived). Even a factor of 4 error in dt leaves it at 15 hours or more. The criterion was stated without this estimate. It needs implicit time stepping, or a smaller reference (fewer columns or a shorter run), and that choice is raised, not made here.

**Reported meanwhile: C1 eqtt 64x12 to 8 ms with the law** (12 uniform rings, wall 600 K; `/tmp/wf3`, not committed): c* 2,426.11 m/s (−0.47% against the predicted 2-D value), vacuum Isp 415.68 s (−3.12%), vacuum thrust 1,630.6 N; wall heat flow −312.8 kW (heat leaving the gas), against −33.7 kW with the no-slip rule on the same grid at 4 ms (`/tmp/wfreg/off4.txt`: c* +0.29%, Isp_vac −0.02% against the same prediction); wall axial force −4,469 N, the same as the no-slip run to 0.01% (it is mostly pressure). Laminar-estimate first-cell y+: median 60, largest 80. Drift over the last ms 5e-8. 218 s on 4 threads alongside other jobs. The prediction has no wall heat loss. With the law, 313 kW leaves through the wall, about 8% of the flow's sensible enthalpy above T_w (derived, c_p about 3.5 kJ/(kg K)), and the Isp gap opens from −0.02% to −3.1%. In the hot pipe of criterion 2 the law gives about 15% less wall heat than resolved SST.

### Criterion 4: cost, pass

C1 eqtt 64x12 to 1 ms, the law against the no-slip wall, 2 threads each, run side by side in two repeats (same load; other jobs were running on the Mac). Wall time per step from whole seconds: no-slip 46 s and 46 s for 17,088 steps (3.51 µs per cell update); law 49 s and 48 s for 17,049 steps (3.70 µs). Extra cost per cell update: **+4.6% to +6.8%** (+9.0% at worst with each time off by half a second), band 10%: **pass**. The four outputs are in `/tmp/wf4` (not committed).

Reported, the payoff: criterion 3(b) itself, C1 to 8 ms with the law on 12 rings, took 218 s on 4 threads; its stated reference is estimated at about 60 hours (above).

### The corrected SST: checks R1 and R3, the corrected references and the a priori check (6 October 2026, 01:00 to 01:20 CDT)

Outputs: `wall_functions/hp_r1_operator_2026-10-06.txt`, `hp_r3_cold_nodic_2026-10-06.txt`, `hp_r3_cold_dic_2026-10-06.txt`, `reference_hp_cold_2026-10-06.txt`, `reference_hp_heated_2026-10-06.txt` (with their 6,400-ring profiles), `log_law_fit_hp_re_tau_1e4.txt`, `apriori_ode_model_hp_2026-10-06.txt` and `heated_hp_scaling_2026-10-06.txt`.

**R1, the operator: pass.**
- With constant ρ and μ, the largest |Φ| is 1.9e-16 of the largest conventional diffusion (limit 1e-12).
- On the manufactured heated profile, the observed orders over N 100, 200 and 400 are 2.001 and 2.000 for Φ_k, 1.990 and 1.998 for Φ_ω, and 2.001 and 2.000 for Φ_CD (band 2 ± 0.25). The L1 error at N 400 is 7.9e-5 for Φ_k, 6.8e-5 for Φ_ω and 4.3e-5 for Φ_CD.
- The error in ∫Φ_k dV is 2.2e-5 at N 400, and it falls as second order.
- The wall cell's Φ_ω error is 0.09 to 0.10. The conventional ω diffusion's own wall-cell error is 0.11 to 0.12, as derived. The axis cell's Φ_k error is 5% and does not shrink with the grid, as the conventional k diffusion's (8%) does not. Its volume shrinks as h², so it does not change the order. It is reported here because the engine's axis cells use the same face rule.

**R3, constant properties: pass on (i) and (iii); (ii) reported.**
- (i) Without D^ic, the Richardson c_f is 0.00336607 (observed order 1.00), 0.150% below 0.00337113.
- (ii) With D^ic, the Richardson c_f is 0.00333874, 0.81% below (i). That is a shift of +0.0996 in u_b+, against the paper's 7.18 M_τ calibration of 0.065.
- (iii) D^ic: as recorded in the restatement above.

**κ, B and A+ from the corrected cold reference** (the R3 pipe with D^ic, 6,400 rings, Re_τ 9,941.5), by design item 7's rule:
- κ is 0.3693 and B is 3.853, against the present 0.3697 and 3.752. Both move less than the fit's own range (κ 0.363 to 0.375, B 3.49 to 4.08), so the present values stay. B's rise of 0.10 is D^ic's shift in (ii).
- A+ refits to 14.5 against the present 14.0, one step of the fit's grid. The a priori check is run with both.

**The corrected heated reference** (criterion 2's case, held fixed: N2, wall 600 K, p0 3 MPa, R 5 mm, f 19,800 N/m³, Q 1.95e9 W/m³; 200 to 6,400 rings).
- Richardson limits (1,600 to 6,400 rings, observed order 0.93): c_f 0.00565337, St 0.00310987, u_b 32.2411 m/s, T_axis 3,167.0 K and T_b 2,792.4 K. Re_τ is 10,382 on 6,400 rings.
- Against standard SST on the same case (c_f 0.00725704, St 0.00394250, T_axis 2,941 K, Re_τ 9,950), the correction lowers c_f by 22.1% and St by 21.1%.
- ∫Φ_k dV is −1.4e-6 of the supplied heat. The S_n floor acts in no cell.
- **Scaling (reported):** under Trettel-Larsson scaling the corrected heated profile lies on the cold one. u_TL+ minus the cold u+ is −0.03, +0.04 and +0.20 at y* 30, 100 and 300. Standard SST sat 2.0 below. This is the collapse the correction is built for.

**A priori, the Kawai-Larsson ODE against the corrected heated reference: fails on τ_w.** Errors, model minus reference:

| y+ | 30 | 100 | 300 | 1,000 |
|---|---|---|---|---|
| τ_w, A+ 14.0 | −1.1% | −1.8% | +1.8% | **+3.9%** |
| q_w, A+ 14.0 | −1.3% | −1.8% | +2.0% | +4.2% |
| τ_w, A+ 14.5 | **−2.6%** | **−4.4%** | −0.9% | +1.4% |
| q_w, A+ 14.5 | −2.6% | −4.3% | −0.5% | +1.9% |

- The pass needs τ_w within 2.5% and q_w within 5% from y+ 30 to 1,000. q_w passes with either A+; τ_w fails with either.
- Against standard SST, the same check was −17% to −25%. The correction removes most of that gap. What is left is the mixing-length model's own buffer layer: on the corrected cold profile with A+ 14.5, its τ_w error at y+ 30 is −2.9% too.
- So the stated fallback applies. The wall model becomes the corrected SST's own inner layer: the 1-D constant-stress, constant-heat-flux layer with the corrected k and ω equations, tabulated the same way, with the a priori check repeated. Nothing else is tried.

### The measured check, pre-registered (6 October 2026, 01:20 CDT, before any measured value is read)

The search (inputs only) found five candidates. Their measured values were not read.
- Back, Massier and Gier, JPL TR 32-415 (1965), a cooled nozzle on heated air.
- Schacht, Quentmeyer and Jones, NASA TN D-2832 (1965), an H2/O2 heat-sink nozzle.
- Schacht and Quentmeyer, NASA TN D-7207 (1973), an H2/O2 calorimeter chamber.
- Marshall, Pal, Woodward and Santoro, AIAA 2005-3572, the GO2/GH2 single element (RCM-1).
- Celano et al., EUCASS 2015, a GOX/GCH4 single element.

**The case: Back, Massier and Gier, JPL TR 32-415** (NASA-CR-57326, free on NTRS).
- **Geometry.** A water-cooled convergent-divergent nozzle with throat diameter 45.8 mm (1.803 in), contraction area ratio 7.75 and expansion area ratio 2.68. The half-angles are 30° convergent and 15° divergent. Upstream is a cooled approach section, 129 mm (5.07 in) in diameter.
- **The gas.** Air heated by burning a little methanol, mixed in a calming section before the nozzle.
- **Why this case.** The calming section separates the wall heat transfer from injector mixing and combustion. So the comparison tests the wall model and the turbulence model, not a combustion model. It has 21 axial stations, and the inlet boundary layer is measured.
- **Its limits.** It is air, not combustion products. Its highest pressures (about 1.7 MPa) are at the bottom of the rocket range.

**The run, chosen by its inputs alone:** the highest stagnation pressure at the highest stagnation temperature, with the longest approach section (the most developed inlet boundary layer). If the TR lists several runs at those settings, the first listed.

**Inputs taken from the TR:**
- the contour;
- the stagnation pressure and temperature;
- the measured wall temperature along the wall, as the boundary condition;
- the inlet boundary-layer thickness, matched by the length of the approach section in the domain;
- the gas, as air with the stated methanol products.

An agent extracts these under the same rule as the RL10 inputs: no measured heat flux is reported.

**The quantity:** the wall heat flux along the nozzle at the TR's stations, from the nozzle inlet to the exit, and its integral over the nozzle wall.

**The runs:**
- the engine as it will run the RL10: corrected SST, the tabulated wall model, and the production grid rule;
- the same with the wall resolved, so the wall model's error is separated from the turbulence model's.

**The bands, stated now:**
- the integrated heat flow within 10%;
- the local heat flux within 20% at every station.

The strongly accelerated stations near the throat (K = ν/u² du/dz above 2e-6 in the simulation) are judged like the rest. RANS without a transition model is known to be weakest there, and the RL10's throat is such a region, so excluding them would hide the error that matters. A failure is not tuned away: its size becomes the wall-heat error band in the RL10 loss ledger.

### Where this leaves the criteria (6 October, 00:45)

| Criterion | State |
|---|---|
| 0 (a) to (d) | pass |
| 1 (a) | running on backhouse; steady state about 6 to 7 h out |
| 1 (b) | pass at y+ 99, 289 and 956; y+ 30 running |
| 2 (a) | not run (the 64-ring stretched pair); expected to fail like 2(b) |
| 2 (b) | **fails**: c_f −17% to −19%, St −15% at y+ 300 and 1,000 |
| 3 | reference infeasible explicitly; 3(b) run and reported |
| 4 | pass, +4.6% to +6.8% |

The AMR wall strip is not set: it waits on the criterion 2 decision.
