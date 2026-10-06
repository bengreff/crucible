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

---

## Results
