# Liquid oxygen spray (TECHNICAL_PLAN step 5)

Status (13 October 2026): the design and criteria below were stated before any spray code ran (8 October). The code is in (`core/spray.hpp`, `core/spray.cpp`, `tests/spray_tests.cpp`), and S0 to S6 have been run (*Results*, below the criteria). S0(b) fails as stated and is kept with the reasons given. S2 fails as stated for the 1/3 rule (+24% to +78% in ṁ); the film integrals that replace it, and a diffusivity fix, are changes after that result and are labelled so; with them S2 passes. S2(c), the cross-check against Lafon et al.'s resolved real-fluid drop, fails as stated (d² slope +30%, lifetime −17.5%); it is read against S2 and nothing is changed. S7 is not yet run.

Step 5 puts liquid oxygen into the engine as Lagrangian parcels: injection with a coaxial-injector drop-size bracket, KH-RT secondary breakup, Abramzon-Sirignano evaporation, and two-way exchange with the gas inside the shared step. The plan's verification is a single drop against a reference solution, mass and energy budgets closing with parcels present, and a published subcritical LOX/GH2 spray flame. The flame is Mascotte A-10 (ONERA, 1 MPa), chosen by inputs only: `docs/validation/MASCOTTE_A10_INPUTS.md`. Its prediction is pre-registered blind in a separate record, committed before any measured value is read.

The physical picture, for orientation: a jet of liquid oxygen at about 90 K leaves a 5 mm post at a couple of metres per second, inside an annulus of hydrogen at room temperature moving more than a hundred times faster. The hydrogen strips drops off the liquid core; the drops are dragged up to the gas speed, heat up to a wet-bulb temperature a little below boiling, and evaporate into a flame that sits in the mixing layer. The model follows the drops (not the liquid core's surface) and asks the gas for its temperature, velocity and composition where each drop is.

Labels: *derived* is my own reasoning, *guessed* a number not yet measured, *cited* from a source I name. The literature behind the choices is a Sonnet agent's search on 8 October; where it confirmed or corrected a constant, the section says so.

## Scope and declared limits

- **Subcritical only.** Oxygen's critical pressure is 5.043 MPa. The liquid model stops with an error if a parcel's cell is above 5.01 MPa (the table's last subcritical row). The RL10A-3-3A chamber (about 3.2 MPa) and Mascotte A-10 (1 MPa) are inside. Near the critical point the surface tension and latent heat fall to zero, and the model follows them through the table, but two effects of high pressure are not modelled: hydrogen dissolving in the liquid, and the real-gas behaviour of the vapour film (the film is the engine's ideal-gas mixture).
- **Dilute spray.** The liquid's volume is not taken out of the gas. Near the injector, where the spray is dense, this is wrong; the largest liquid volume fraction in any cell is reported in every run so the reader can see where.
- **Drops.** Spherical for heat and mass transfer, uniform temperature inside (infinite liquid conductivity), pure O2.
- **Gas at a parcel.** The containing cell's values, not interpolated.
- **No swirl in the gas.** Parcels move in three dimensions, but the gas has no azimuthal velocity, so the azimuthal momentum the gas would receive is dropped (counted in the ledger) and its kinetic energy enters the gas as heat. Energy is kept exactly; angular momentum is not.
- **Walls.** Parcels that reach a wall rebound (normal restitution 1 by default, declared). A hot chamber wall is far above liquid oxygen's Leidenfrost temperature, so a vapour cushion rather than wetting is the expected contact (*derived*, not tested). No wall film.
- **No envelope flame around a drop.** The film around each drop is non-reacting; the flame lives in the gas cells. The film theory is that of an evaporating, not a burning, drop.

## Liquid properties

- **Source.** NIST Chemistry WebBook, oxygen saturation properties (`tools/o2_saturation_nist.tsv`, 600 rows uniform in T from 54.361 K to 154.414 K, fetched 8 October): liquid density, enthalpy, cp, viscosity, conductivity, surface tension, vapour pressure, and the saturated vapour's enthalpy. Linear interpolation in T between rows; ln p for the vapour pressure. The enthalpy's inverse is exact for the interpolant, so the drop temperature follows from its enthalpy with no inversion error.
- **One energy zero.** The engine's gas uses NASA-7 enthalpies (zero for O2 at 298.15 K). The liquid is placed on the same scale through NIST's own ideal-gas enthalpy at 298.15 K (271.2601698 kJ/kg on NIST's reference, the 298.15 K isotherm at 1e-6 MPa):
  h_l,engine(T) = h_O2,medium(298.15 K) + h_l,NIST(T) − 271.2601698 kJ/kg.
  The vapour leaving a drop at T_s carries h_l(T_s) + L(T_s), with L = h_v − h_l on saturation from the same table: the real saturated vapour's enthalpy. The ideal-gas medium then holds that energy at a temperature somewhat below T_s; the difference is the vapour's real-gas departure, about 15 kJ/kg at 1 MPa (*derived* from the table and the ideal-gas isobar). Energy is exact either way; this choice gives the drop its true latent heat.
- **The cold gas.** The medium's NASA-7 fits for O2 start at 200 K. Below it they are extrapolated, and the film around a LOX drop sits at about 100 to 250 K. Criterion S0(b) measures what the extrapolation costs against NIST's ideal-gas enthalpy.

## The drop model

- **Evaporation and heating**: Abramzon and Sirignano (Int. J. Heat Mass Transfer 32, 1989), as in the paper.
  - Spalding mass number B_M = (Y_s − Y_∞)/(1 − Y_s), with the surface vapour mass fraction from Raoult's law for a pure liquid, X_s = p_sat(T_s)/p. Where the gas already holds more oxygen than the surface (B_M < 0), B_M is set to zero: no condensation (counted).
  - Sh_0 = 1 + (1 + Re Sc)^(1/3) f(Re), Nu_0 likewise with Pr, f = max(1, Re^0.077); F(B) = (1 + B)^0.7 ln(1 + B)/B; Sh* = 2 + (Sh_0 − 2)/F(B_M), Nu* = 2 + (Nu_0 − 2)/F(B_T).
  - ṁ = π d ρ_g D Sh* ln(1 + B_M); B_T = (1 + B_M)^φ − 1 with φ = (c_pF/c_pg)(Sh*/Nu*)/Le, iterated to 1e-12; heat conducted to the surface Q = ṁ c_pF (T_∞ − T_s)/B_T.
  - Film properties by the 1/3 rule: T_r = T_s + (T_∞ − T_s)/3 and each mass fraction likewise, with the engine's mixture-averaged transport at (T_r, p, Y_r): μ, λ, c_pg and O2's mixture diffusion coefficient D. c_pF is O2's own cp at T_r. Re uses the cell's gas density and the film viscosity.
  - **Changed after S2's result** (*Results*, S2): the 1/3 rule still sets Re, Pr, Sc, c_pg and the drag, but the conductances λ and ρD in ṁ, B_T and Q come from integrals across the stagnant film (`Spray::filmConductances`, 3 Gauss points, 3 passes against B_T), and D is the Blanc's-law (mass-fraction) diffusivity. `Settings::filmPoints = 0` restores the paper's 1/3 rule.
  - The drop: dm/dt = −ṁ and m dh_l/dt = Q − ṁ L.
  - **Boiling and flashing.** If T_s would pass the saturation temperature at the cell's pressure, the drop is held there and ṁ = Q/L with the classical boiling form (B_T = c_pF (T_∞ − T_s)/L). A drop that finds itself above saturation (it moved to a lower pressure, as in the nozzle) evaporates the enthalpy above saturation at once (a flash, counted).
- **Drag.** A sphere's drag with the evaporating-drop correction. Form and constants: *Literature*, below.
- **Breakup**: KH-RT, Reitz (1987) wave model with the Rayleigh-Taylor branch of Patterson and Reitz (1998) and Beale and Reitz (1999). Equations and constants: *Literature*, below.
  - KH: each drop's radius shrinks toward B0 Λ at the rate set by B1; the mass shed accumulates on the parcel and is released as a child parcel of drops of radius B0 Λ once the parcel's shed mass reaches 3% of its injected mass, with the parent's velocity. *(Corrected 14 October: the code first compared with 3% of the parcel's current mass; see Results, KH child release. OpenFOAM's ReitzKHRT, read in its source, compares with 3% of the cloud's average injected parcel mass, `msLimit` 0.03; Reitz 1987 and Beale and Reitz 1999 were not reached.)*
  - RT: active beyond the breakup length; when its wavelength fits in the drop and its time has run, the drop splits to the RT wavelength, keeping mass.
  - Drops too small to follow (below 0.1 μm, *guessed*, far below anything that matters) become vapour at once; counted.
- **Integration.** Per parcel, in sub-steps:
  - Each sub-step is bounded by a quarter of the parcel's relaxation times (momentum, mass, heat) and by half the cell's size along its velocity.
  - Velocity: the exact exponential relaxation toward the gas velocity at the sub-step's drag.
  - Mass and enthalpy: Heun's method on s = m^(2/3) (proportional to d²) and h_l. The d² law is linear in s, so it is integrated exactly; S1(a) checks this.
  - Position: the trapezoid of the old and new velocity, in three-dimensional Cartesian coordinates with the parcel's meridional plane as the x-z plane. Afterward the parcel is rotated back to its new meridional plane. A parcel crossing the axis comes out on the other side with no special case.
- **Turbulent dispersion** is not in this first build. The eddy-interaction model (Gosman and Ioannides 1981, with the SST k and ω) is planned before the Mascotte prediction, with its own check (a tracer's dispersion against the model's diffusivity). Without it, the spray spreads only as the mean gas carries it.

## Injection and the drop-size bracket

- **Liquid streams.** Each injector is a release annulus on the injector plane with Supply's valve schedule (opening time, linear ramp). Parcels are released at a fixed rate (2 million per second by default, *guessed*), each carrying the mass delivered since the last one. Mass not yet released is held by the injector and counted in the liquid's ledger.
- **Two size modes.**
  - **Distribution (the plan's primary).** Diameters drawn from a Rosin-Rammler distribution at a Sauter mean diameter from a published coaxial-injector correlation, released from the liquid core's surface along a correlation's core length. The **bracket** is the spread of the published correlation at the case's inputs (restated under *Literature*: Burick's curve and its data scatter, as no closed-form coaxial correlation was found): the low and high diameters give two runs, and the prediction is a band between them. Correlations fitted to the validation case's own drop measurements are excluded, so the case stays blind.
  - **Blob.** Drops the size of the post (Reitz's blob method), with KH-RT doing the atomisation. It is reported beside the bracket as a third run, not used to set it.
  - In both modes KH-RT acts on every drop above its stable size.
- **Device thrust.** The injector's reaction to the liquid it delivers (ṁ v_inj, positive against the exhaust) and the walls' impulse on rebounding parcels join the gas ledger's device thrust. The Isp denominator includes the liquid flow.

## Coupling inside the shared step

- **Where.** In the shared step after the flow step, with the same dt:
  - finite rate: react(½ dt), flow(dt), spray(dt), react(½ dt);
  - equilibrium and frozen: flow(dt), spray(dt), react.
- **Exact by construction.** Each parcel's sub-step leaves a record: its cell and the change of its mass, axial and radial momentum, and enthalpy plus kinetic energy. The gas receives minus that change, as O2 mass, momentum and total energy (`Flow::exchange`, booked in the gas budgets as an external exchange). So mass, axial momentum and energy are conserved to rounding whatever the integration error of the drop model.
- **Point-implicit seen state** (*derived*, the particle-source-in-cell practice made implicit per cell). Near the injector a cell can hold many times its own gas mass in liquid. A drop model that saw the gas frozen over the step would then draw more heat than the cell holds. So before the parcels move, each cell's drag and heat conductances are summed over its parcels (Σ m_i/τ_i and Σ h_i A_i). The parcels then see the backward-Euler end state of the cell's velocity and temperature under those loads:
  - u* = (ρV u + dt Σ K_i v_i)/(ρV + dt Σ K_i);
  - T* = (ρ c_v V T + dt Σ H_i T_i)/(ρ c_v V + dt Σ H_i).
  The heat drawn in a step therefore cannot take the cell below the drops' temperature. In a dilute cell u* and T* equal the cell's own values to within the loads' small ratio. A flag turns it off for the checks.
- **Threads.** Parcels move in parallel on the flow's pool; each writes only its own record. Records are summed into cells, children appended and dead parcels removed in parcel order, serially. Each parcel draws from its own random stream (splitmix64), seeded from its injector's stream at release. The result is the same bit for bit on any number of threads.
- **Mesh.** The cell is found from z (uniform dz), the wall radius linear between stations, and a binary search on the ring fractions. This is the one place the spray knows the mesh layout; AMR (inquiry2) will replace it. On AMR's levels: *Spray on AMR* (below), stated 22 October before its code.

## Turbulent dispersion (stated 13 October, before its code)

With SST on, a parcel sees the cell's mean velocity plus a fluctuation u′, after Gosman and Ioannides (1983), the eddy-interaction model.
- **The eddy.** u′'s three components (axial, radial, azimuthal) are independent Gaussians of variance 2k/3, drawn from the parcel's own stream (Box-Muller). The parcel keeps the same u′ for the interaction time t_i, as it crosses cells. It then draws again in the cell it is in. The sub-step ends at the interaction's end. k and ω are the cell's own, ρk/ρ and ρω/ρ. Where k ≤ 0, u′ = 0 and t_i = 0, so it draws again at the next sub-step.
- **Interaction time.** t_i = min(τ_e, t_c), both drawn at the eddy's start.
  - The eddy lifetime is τ_e = l_e/√(2k/3).
  - The crossing time is t_c = −τ_p ln(1 − l_e/(τ_p |u + u′ − v|)), and is infinite when l_e ≥ τ_p |u + u′ − v|.
- **The eddy length (the one choice made here).** l_e = (3/Sc_t)·√(2/3)·√k/ω, where Sc_t is the gas's turbulent Schmidt number (0.7).
  - *Derived:* the velocity is held constant for τ_e at variance 2k/3 per component. So a tracer's long-time diffusivity is (2k/3)·τ_e/2 = k/(Sc_t ω), which is the gas's own ν_t/Sc_t (SST's k/ω away from its limiter). A drop small enough to follow the gas therefore spreads as fast as the vapour it becomes.
  - Gosman and Ioannides' own length is C_μ^{3/4} k^{3/2}/ε = C_μ^{−1/4}√k/ω, about 1.83√k/ω. That gives tracers a diffusivity 0.52 times the gas's (a Schmidt number of 1.34). This choice uses 3.50√k/ω.
  - In Lagrangian integral-time terms (T_L = τ_e/2 for a held velocity), this is T_L = 0.19 k/ε. Gosman and Ioannides' is 0.10 k/ε, and the values used in the literature run from 0.15 to 0.3 k/ε (Fluent's discrete random walk 0.15, *recalled*).
  - This sets a model constant from the gas model's own constant, not from any measurement.
- **Coupling.** The drag relaxes the parcel toward the seen velocity plus u′. The azimuthal part relaxes toward u′_θ, and its momentum goes where the azimuthal momentum already goes (the ledger). The point-implicit loads take the parcel's velocity minus its u′, so the backward-Euler seen state stays the end state under the drag the gas actually feels. The exchange stays minus the parcel's change, so the budgets stay exact. Children inherit the parent's eddy and its remaining time (same place, same eddy).
- **Declared limits.**
  - Isotropic u′ (SST carries no anisotropy).
  - In non-uniform k the model drifts tracers toward low turbulence (the well-known violation of the well-mixed condition), with no drift correction.
  - No correlation between the drag and the evaporation fluctuations beyond the slip's.
  - With SST off there is no dispersion.

## Walls, axis and outlet

- **Injector plate (z = 0).** Mirror, with the restitution.
- **Side wall.** Mirror about the local wall line, the frustum's generator in that column. The parcel's change of axial momentum is the wall's thrust impulse; with restitution below 1 the lost normal kinetic energy is booked as wall loss.
- **Axis.** No boundary (the Cartesian move).
- **Outlet.** The parcel leaves; its mass, energy and axial momentum are booked as outflow.

## Budgets

- **The liquid's ledger.** For each of mass, energy and axial momentum: held now − held at start = injected − exchanged to the gas − outflow (+ wall impulse − wall loss).
- **The gas's budgets.** Flow's as before, with the exchange added as an external term (`Measurements::exchangedMass`, `exchangedEnergy`, `exchangedMomentum`).
- **The whole system.** Gas plus liquid closes with no exchange term.

## Criteria (stated before code)

- **S0. Properties.**
  - (a) The liquid temperature recovered from the enthalpy equals the temperature put in within 1e-12 relative, at 10,000 temperatures across the table.
  - (b) The medium's O2 enthalpy (NASA-7 from Cantera's h2o2.yaml, extrapolated below 200 K) against NIST's ideal-gas enthalpy (the 1e-6 MPa isobar, 60 to 300 K, `tools/o2_lowp_nist.tsv`), with both on the engine's zero at 298.15 K: within 1 kJ/kg everywhere (about 1 K), *guessed* band. If it fails, the remedy is to fit the liquid's zero and the cold gas together, and the design comes back here.
- **S1. The drop model's algebra and integrator.**
  - (a) d² law. A drop at its own wet-bulb temperature (found by bisection on Q = ṁ L with the model's own rates) in stagnant gas, one-way (the gas is not changed): ṁ ∝ d, so d² falls linearly. The time to reach the vaporisation diameter must equal (d0² − d_min²)/K within 1e-10 relative, K from the rates at the start.
  - (b) The A-S forms (Sh*, Nu*, B_T) against the paper's formulas written out independently in the test, for Re from 0 to 400 and B_M from 1e-6 to 20: within 1e-12 relative.
  - (c) (Added 13 October, with the code it checks, after S2's result.) The film's lean transport (`Spray::filmTransport`) against Medium::transport at 10,000 random states (60 to 3,500 K, 0.1 to 5 MPa, each species absent with probability 1/2, every 100th pure O2): conductivity, O2's mixture D and the gas constant within 1e-13 relative.
- **S2. A single drop against a reference solution.**
  - **The reference.** The steady spherically symmetric evaporation of a drop in stagnant gas, solved resolved in r with variable properties: Cantera's h2o2.yaml with mixture-averaged transport and a correction velocity. Independent Python, `tools/drop_film_1d.py`, with the same surface condition (Raoult with NIST's vapour pressure) and the same latent heat. (Restated: written in C++ against Cantera's C++ API, `tools/drop_film_1d.cpp`, because the Mac has no Python Cantera. It shares only the LiquidOxygen table with the engine.) It needs no 1/3 rule, so it tests the model's film-property approximation in a hydrogen film, which is far from the hydrocarbon-in-air cases the rule was built on.
  - **Cases.** Wet-bulb drops at 1 and 3 MPa in pure H2 at 300 K and 1,000 K, and in a 1,500 K H2/H2O mixture (Y_H2O 0.6).
  - **Bands** (*guessed*): ṁ/d within 10% and the wet-bulb temperature within 2 K. If the 1/3 rule fails here, it is reported and the design comes back here.
  - (c) Cross-check against Lafon et al. (2008), the 10 atm, 1,000 K hydrogen, 100 μm, 90 K case: the drop's lifetime and the slope of d² after heat-up, digitised from their figures, within 15% (*guessed*). Their model has finite liquid conductivity and real-gas properties, so a miss here is read against S2's own result before anything is changed. Stated before the figures are read.
- **S3. Budgets with parcels.** A 16 by 6 cold chamber (H2 supply, a LOX injector, frozen chemistry) run to 0.5 ms with injection, evaporation, wall rebounds and parcels leaving:
  - (a) the gas's mass, energy and axial-momentum balance errors (with the exchange terms) and the liquid's ledgers, each within 1e-11 of its throughput;
  - (b) the whole system's, with no exchange terms, likewise.
- **S4. Threads.** The S3 run on 1, 4 and 6 threads: the final gas state and parcel list identical bit for bit (`cmp` of raw dumps).
- **S5. Breakup.**
  - (a) The KH and RT wavelength and growth-rate fits against the published formulas written out independently in the test, within 1e-12.
  - (b) A blob in a uniform fast gas stream, one-way: total drop mass conserved through every child to 1e-13 relative; each child's diameter equals 2 B0 Λ at its creation.
- **S6. Cost** (measured, not a pass/fail): the spray's share of the step on the S3 run, and parcels times sub-steps per second. If the spray costs more than the gas step, the film properties are tabulated (not built until needed).
- **S8. Turbulent dispersion** (stated 13 October, before code).
  - (a) Spreading. Drops in a frozen, uniform gas at rest, one-way, with uniform k and ω. The gas is pure O2 at 1 atm and 85 K, with the drops at 85 K, so no mass or heat passes (the clamped branch; corrected after the run: no sub-step was clamped, and the check that no mass passed holds). The drops are sized so that τ_p = τ_e/20 at zero slip, which keeps the crossing time from ever binding, and the slip Reynolds number is below about 0.01. 10,000 of them start at one point on the axis. After N = 10 eddy lifetimes, compare ⟨Δz²⟩ and ⟨r²⟩ against the exact expectation for linear drag and a held velocity, σ²[N τ_e² − 2τ_p τ_e(1 − q^N) + τ_p²(1 − q)(1 − q^{2N})/(1 + q)] with q = e^{−τ_e/τ_p} and σ² = 2k/3 (twice that for r²). Each must fall within 3 standard errors of the sample: √(2/n) relative for Δz², 1/√n for r². ⟨Δz⟩ must be within 3 standard errors of zero. *Derived* expectation. The band is statistical only.
  - (b) The interaction time against the formulas above written out independently in the test, both branches, within 1e-14 relative.
  - (c) S3's chamber with SST on (inflow intensity 5%, viscosity ratio 10) and dispersion, to 0.2 ms: the budgets of S3(a) and (b) within 1e-11, and S4's bit-for-bit identity on 1 and 4 threads.
- **S7. Mascotte A-10.** Validation, pre-registered in its own record before any measured value is read: the case's modelling choices (equivalent round chamber, igniter, wall condition, the bracket's correlations), the predicted quantities and their bands. Requires turbulent dispersion first.

## Results (13 October 2026)

`build/crucible_spray_tests` runs S0, S1, S3, S4 and S5 in about 90 s; `crucible_spray_tests cost [threads]` runs S6. All figures are *measured* unless labelled otherwise.

- **S0(a)** passes: 3.4e-16 for T(h(T)) and 2.4e-16 for T_sat(p_sat(T)).
- **S0(b) fails as stated: 1,245 J/kg at 75 K** against the 1 kJ/kg band. Here is h_medium − h_NIST on the engine's zero at each temperature:

  | T [K] | 75 | 90 | 100 | 120 | 150 | 200 | 250 | 300 |
  |---|---|---|---|---|---|---|---|---|
  | J/kg | −1245.1 | −873.7 | −673.5 | −370.2 | −106.7 | +23.5 | +10.9 | +0.2 |

  The error is the NASA-7 extrapolation below 200 K. The pre-stated remedy was to fit the liquid's zero and the cold gas together. It is not taken. Bringing 75 K inside the band needs a shift of at least 245 J/kg, and a fit over the cold range about 600 J/kg (*derived* from the table). Either moves that much error into the enthalpy the LOX carries into hot gas, where the medium is accurate and the chamber's energy budget is decided. The 298.15 K anchor keeps that delivered energy right. The cost falls only on the cold film: the vapour leaving a drop is placed at a temperature off by Δh/c_p, about 1.4 K at 75 K, 0.95 K at 90 K and 0.40 K at 120 K (*derived*, c_p of O2 gas 915 J/(kg K)). Mass and energy are still exchanged exactly (S3). The film's transport properties see the same offset, which is far inside S2's 2 K band. This is recorded as a failed criterion with its engineering resolution, not as a pass.
- **S1(a)** has two cases at 1 MPa, both evaporating with the film model (the H2 case boiled before it).
  - N2 at 300 K: the wet bulb is 108.5532 K (T_sat 119.6212 K), B_M 1.102 and K 9.716e-8 m²/s. d² against d0² − Kt at 1,000 times: 2.9e-15. The drop is alive at (1 − 1e-10) of the predicted life and vaporised by (1 + 1e-10). Pass.
  - H2 at 1,200 K: the wet bulb is 115.1556 K, B_M 51.1 and K 1.873e-6 m²/s. d² is within 7.0e-15 and the lifetime test passes. Pass.
  - **Superseded finding: "in hydrogen the A-S wet bulb is saturation, to double precision".** It held for the 1/3 rule with the mole-form D (both below, S2): its φ at 1,000 K was about 0.044, and the wet-bulb balance then needed B_M ≈ 2e17, beyond a double, so hot-H2 drops heated to saturation and boiled. The resolved film (S2) puts the wet bulb 4.9 K below saturation at 1 MPa and 1,000 K, and the film model reproduces it within 0.04 K. The 1/3 state sits in hot, hydrogen-rich gas, while the film's φ = c_pF ρD/λ is weighted by the cold, oxygen-rich layer next to the drop, where it is about 0.6 (*derived*). Even the 1/3 rule with the corrected D puts the 1,200 K wet bulb only 2.7e-5 K below saturation; that case took 17 million sub-steps and was not kept.
  - Two engine findings came out of S1(a) and still stand. (1) A boiling drop is held while heat arrives, judged by the heat's sign. Q − ṁL is zero only to within the rounding of a contracted multiply-add, and a −1e-19 W residual had released the drop from saturation. (2) **The medium's T(e) has two roots 2.6e-6 K apart at NASA-7's 1,000 K junction for H2**, and which one Medium::temperature returns depends on its guess (from 300 K: 999.99999999999989 K; from that value: 999.99999738817553 K). This is harmless for the physics and deterministic (the guesses do not depend on threads). But it moved the seen saturation point by 4.5e-8 K, so S1(a)'s hydrogen case runs at 1,200 K.
  - **Declared limit: the boiling hold does not release.** A drop at saturation stays held while heat arrives. With a carrier gas present its wet bulb lies below saturation (4 to 6 K in hydrogen, S2), so a drop that reaches saturation should cool back to it. A drop reaches saturation by a flash (a pressure drop) or in nearly pure O2. From its wet bulb in a hydrogen-rich cell, a flash needs the saturation temperature to fall 4 to 6 K, a pressure drop of order 40% at 1 MPa (*derived* from the NIST table), which happens in the nozzle, not the chamber. Flashes are counted (`Ledger::flashes`); the count on the Mascotte run decides whether the release is built.
- **S1(c)** passes: 6.9e-16 largest relative difference over the 10,000 states.
- **S1(b)** passes, run with `filmPoints = 0` (the paper's forms) and the Blanc's-law D on both sides. Across 208 evaporating cases (B_M 1e-6 to 20, Re 0 to 400) the largest difference is 6.5e-16, in B_T; the boiling cases show none. B_T is solved by secant steps from A-S's fixed-point start. That needs 5.9 residual evaluations per solve on the S3 run, against 10.0 for the plain fixed point. The film passes restart the secant from the last pass's B_T. **Solver fault found by S6 and fixed:** a residual of exactly zero inside the secant loop left the solve marked unsettled, so it fell to bisection; with the film's warm starts that was 24% of solves on the S3 run. The settle test is also now a step of 1e-11 relative (was 1e-14, at the residual's rounding floor for large B: about 4e-15 at B 50, *derived*); the superlinear step after it lands at rounding. Now 4.5 residual evaluations per solve and no bisection on the S3 run; S1(b) is 8.0e-16.
- **S2 fails as stated for the 1/3 rule, and passes with the film integrals** (changes after the result). `build/crucible_drop_film_1d [steps]` runs it. The reference integrates d ln X_c/dξ and dT/dξ by RK4 in ξ = a/r from the far state, shooting on q = ṁ/(4πa) and the wet bulb. At 2,000 steps, doubling them moves T_wb by 1.2e-11 K at most and ṁ/d by 6e-13 relative. Results (model at its defaults, 3 points and 3 passes; the 1/3 rule beside it):

  | Gas, p | T_wb ref [K] | ṁ/d ref [kg/(m s)] | film: ṁ/d, T_wb | 1/3 rule: ṁ/d, T_wb |
  |---|---|---|---|---|
  | H2 300 K, 1 MPa | 108.7334 | 3.415e-4 | −0.04%, +0.008 K | +24.4%, +9.97 K |
  | H2 300 K, 3 MPa | 127.5070 | 3.703e-4 | −0.06%, +0.007 K | +45.2%, +13.74 K |
  | H2 1,000 K, 1 MPa | 114.7248 | 1.265e-3 | −0.01%, −0.011 K | +61.2%, +4.90 K |
  | H2 1,000 K, 3 MPa | 135.7895 | 1.372e-3 | −0.18%, −0.039 K | +78.4%, +5.91 K |
  | H2/H2O 0.4/0.6 1,500 K, 1 MPa | 115.9175 | 1.660e-3 | −6.63%, +0.533 K | +35.9%, +3.69 K |
  | H2/H2O 0.4/0.6 1,500 K, 3 MPa | 137.2477 | 1.781e-3 | −6.83%, +0.708 K | +48.3%, +4.44 K |

  - Bands: ṁ/d 10%, T_wb 2 K. The 1/3 rule fails all six; its wet bulbs sit within 0.9 K of saturation.
  - **Bug found (*derived*): A-S's ρD must be the mass-fraction diffusivity.** A-S write Fick's law on Y_O2. Medium's D_km belongs with gradients of X_k. For O2 against carriers in fixed proportions (molar mass W_c) the two differ: D_Y = D_km W̄/W_c (Blanc's law), about 0.12 at the surface in hydrogen, about 0.4 at the 1/3 state and about 1.1 in nitrogen. Fixed in the model and in S1(b)'s independent forms.
  - **The film integrals (*derived*).** With η = −ln(1 − Y_O2), O2 the only species with a net flux, Fick's law on Y with the Blanc's-law D and the wet-bulb energy balance: ṁ/(4πa) = ∫ρD dη and ln(1 + B_T) = ∫φ dη over the film's (T, Y) path, φ = c_pF ρD/λ, and θ = T − T_s + L_eff/c_pF falls as e^(−∫φ dη). With constant properties these are A-S's own ρD ln(1 + B_M) and φ ln(1 + B_M). So the model keeps A-S's forms and takes ρD as the mean of ρD over η and λ as c_pF times that over the mean of φ. On s = η/ln(1 + B_M) the mass fraction is exact by construction; T follows the θ law with the shape F(s) = ∫₀ˢφ/∫₀¹φ, taken from the previous pass (the first pass uses constant φ). Gauss-Legendre in s, split at 8/ln(1 + B_M) when that is below 1. Boiling uses the pure-O2 film, where λ = 1/∫dξ/λ is exact. φ varies about 20-fold across a LOX film in hydrogen (about 0.03 far out, about 0.6 in the cold O2-rich layer), which is why one film state cannot carry it.
  - Convergence (measured): with the path's shape iterated to 1e-10, 4, 8 and 16 points agree within 0.05% in ṁ/d and 2 points are 1% to 13% low; at 3 fixed passes, 3 points are within 0.4% of 16. Against the iterated result (4 points), one pass (the constant-φ shape) is 13% to 23% low, 2 passes are within 1% and 3 within 0.1%. The defaults (3 points, 3 fixed passes) were chosen after these results, for cost; a fixed pass count keeps the rates a smooth function of the drop's state.
  - The H2/H2O residual (−6.6% to −6.8%) is *inferred* to be the model's carriers in fixed proportions: the reference lets H2 and H2O separate across the film by their own diffusivities. Pure-hydrogen cases are within 0.2%.
  - **Engine finding (measured): Cantera's default transport fits fail below 300 K.** Cantera fits each species' transport from the highest lower limit of the species' thermo ranges, 300 K for h2o2.yaml (AR and N2). Evaluated below that, as the engine does (`Medium::transport`, unclamped), H2's conductivity goes negative below about 97 K (−0.099 W/(m K) at 60 K), and O2's is +137% at 90 K and +232% at 73 K against NIST. LOX films and the cold O2 they shed live there. `thermo::transportFits(mechanism, minTemperature)` now widens every species' first range to minTemperature before Cantera fits (60 K for liquid-oxygen flows; H2O's collision integrals are tabulated down to T* = 0.1, which is 60 K). Against NIST (`tools/o2_gas_0.01MPa_nist.tsv`, `tools/h2_gas_0.1MPa_nist.tsv`, `tools/h2_gas_1MPa_nist.tsv`): O2 viscosity within 1.2%; O2 conductivity +13% at 73 K, +10% at 90 K, +4.6% at 150 K; H2 viscosity −4% to +1.5% at 0.1 MPa (−7.6% at 60 K and 1 MPa, where real-gas effects enter); H2 conductivity −7% to +12%. The price is about 2% at 300 K and above (H2 conductivity −1.9% against NIST at 300 K, where the default fit is +0.1%). The default (no minTemperature) is unchanged, so other tests are untouched; the spray tests and S2 use 60 K. The same flaw affects any engine run with gas below 300 K.
- **S8** passes (13 October, first run):
  - (b) 400 cases, 102 bound by the crossing time and 298 by the lifetime: 3.9e-16.
  - (a) Drops of 15.08 μm (τ_p 2.143e-3 s, τ_e 4.286e-2 s, slip Re at 3σ 7.8e-3), 10,000 of them, 8.1 million sub-steps; all kept, no mass passed. ⟨Δz⟩ 0.48 standard errors from zero; ⟨Δz²⟩ +0.28% (0.20 s.e.); ⟨r²⟩ −1.11% (1.11 s.e.).
  - (c) 2,168 steps to 0.2 ms: 39 parcels injected, 10 children, 48 held, all inside an eddy at the end. Gas balances 1.9e-16, 1.4e-16 and 1.7e-18; liquid ledgers 2.7e-14 to 2.9e-14; whole system 2.7e-15, 1.1e-15 and 2.7e-18. Bit-identical on 4 threads.
  - S3 and S4 (SST off) are unchanged to the last digit, so the laminar path is untouched.
- **S2(c), the model's side, committed before Lafon et al.'s figures are read** (13 October; `build/crucible_spray_tests lafon [dt]`). Their drop (100 μm, 90 K, quiescent H2 at 1,000 K and 10 atm), one-way, at the defaults:
  - lifetime **7.007 ms**; −d(d²)/dt **1.602 mm²/s**, fitted by least squares over (d/d0)² from 0.6 down to 0.1;
  - wet bulb 114.926 K (T_sat 119.847 K at 10 atm), reached to 99% of the rise by 0.557 ms; the drop swells to (d/d0)² 1.058 while it heats (the liquid's density falls).
  - The outer step moves none of these beyond the fourth figure (2, 1 and 0.5 μs).
  - To be judged: their (d/d0)² curve is digitised, the same window fitted, and their lifetime read where it reaches zero, against the 15% bands stated above. Expected sources of a difference (*inferred*, written before reading): their finite liquid conductivity lengthens the heat-up; their real-gas surface equilibrium (hydrogen dissolved in the liquid, non-ideal vapour) moves the wet bulb.
- **S2(c) fails as stated** (judged 13 October). Their curves are digitised in [spray/lafon_2008_digitised.md](spray/lafon_2008_digitised.md) (Figs. 8a and 4, by a research agent, after the model's numbers were committed).
  - **Slope.** Same window: their −d(d²)/dt is 1.228 mm²/s (*derived* from the digitised points), and ours is 1.602. Ours is +30.4% against the 15% band.
  - **Lifetime.** Theirs reaches zero at 8.49 ms by the same line (8.47 from the last six points; their Fig. 9a reads about 8.3). Ours is 7.007 ms, −17.5% (−15.6% against Fig. 9a's 8.3).
  - **Surface temperature.** Theirs plateaus at 111 to 112.3 K and is reached in under 0.02 ms at the surface. Ours is a uniform drop at 114.93 K, reached in 0.557 ms.
  - **Swelling.** Their drop does not swell. Ours swells to (d/d0)² 1.058, because the whole drop heats.
  - **Read against S2, as stated.** S2 passes against our resolved quasi-steady film reference (`tools/drop_film_1d.cpp`), which shares the model's assumptions: ideal-gas film properties, a uniform-temperature liquid, no hydrogen in the liquid, and a quasi-steady gas phase. The miss therefore lies in physics both lack, not in the Abramzon-Sirignano implementation. Nothing is changed.
  - **Candidate causes (*inferred*).** The liquid's own heating, as an effective latent heat, explains only about 4%: B_T is about 60 in hydrogen, so ṁ goes as ln(1 + B_T) and is insensitive to it. The rest lies in what their model has and ours does not: dense-fluid properties of the cold vapour layer, hydrogen dissolution, corresponding-state transport and an unsteady gas phase. These were not separated here.
  - **For S7 (*inferred*).** In hot hydrogen-rich gas the drop model may evaporate up to about 30% too fast. That would shorten the spray and the flame. The Mascotte prediction states this bias direction.
- **S3** passes (rerun with the film integrals and the 60 K fits). The run takes 5,415 steps to 0.5 ms and includes 99 parcels injected, 654 children, 4 wall hits, no flashes, 3 vaporised and 0.95 million parcel sub-steps. Of the 1.0e-6 kg of liquid injected, 2.58e-7 kg evaporated, 1.48e-7 kg left through the outlet and 5.94e-7 kg is held. Balance errors:
  - Gas (with the exchange): mass 1.3e-15, energy 6.3e-16, axial momentum 2.3e-18.
  - Liquid ledgers: 6.7e-14, 5.3e-14 and 4.0e-14.
  - Whole system (no exchange terms): 1.4e-14, 5.3e-15 and 1.1e-17.

  The wall impulse is zero here: every hit was on the straight side wall, which has no axial normal.
- **S4** passes: the final gas state and parcel list on 4 and 6 threads are identical to 1 thread, bit for bit.
- **S5** passes.
  - (a) The largest difference is 3.4e-16.
  - (b) One-way, a 1 mm blob at 90 K in H2 at 1 MPa and 300 K, streaming at 250 m/s:
    - KH alone (RT disabled), to 0.2 ms: 207 children, liquid mass through every child within 2.7e-15. The 23 children of the blob match 2 B0 Λ at their creation within 1.3e-16. *After the release fix (14 October): 17 children, all of the blob, mass within 5.3e-16, diameters exact.*
    - KH-RT, to 0.5 ms: 4 RT events of the blob, mass within 1.4e-14. All of the liquid has evaporated by 0.5 ms.
  - The diameter check covers the blob's own children. The blob takes one sub-step per step, so its state at a child's creation is known from outside; smaller drops take several.
- **S6** (S3 run, one thread, before the secant change). The gas alone took 1.71 s; with the spray, 5.19 s. **The spray is 3.48 s, 67% of the step, 2.0 times the gas step.** That is 3.4e5 parcel sub-steps per second (2.9 μs each), against 3.3 μs per cell-step for the gas. On 4 threads both were slower than on 1 (96 cells is too small to share).
  - A sampled profile of the spray's time: the B_T solve and the rest of `rates` (expm1, pow, log1p and its own body) about 57%, mixture-averaged transport 28%, NASA-7 properties 10%.
  - The stated trigger (spray costlier than the gas step) is met. **The film-property tables are still not built**, for two reasons. First, they would remove at most the 38% that is transport and properties. Second, S3's gas step is unrepresentatively cheap: 96 cells hold 2.3 parcels each, and the ratio scales with parcels per cell. This changes the trigger after the result, and is labelled so. The ratio will be measured on the Mascotte setup before its run, and the tables follow if the spray costs more than the gas there.
- **S6 with the film integrals** (S3 run, one thread, 3 points and 3 passes). Each `rates` call makes 10 transport evaluations, not 1. In order:
  - First build, load average about 13: the spray 41 s against the gas's 1.8 to 2.9 s. Profile: transport about 40%, the B_T solves about 47%.
  - Warm-starting the secant from the last pass: 26.5 to 28.0 s (same load).
  - A lean film transport (`Spray::filmTransport`: Medium's conductivity and O2 mixture-D formulas alone, absent species skipped, from the flow's own fits; Medium itself not touched; held to Medium::transport within 6.9e-16 by a new check, S1(c)): 15.7 s, load about 5.
  - The solver fault above, fixed: **the spray 6.15 s against the gas's 1.66 s, 3.7 times the gas step, 6.5 μs per parcel sub-step** (load about 5). The 1/3 rule's figure (2.0 times, 2.9 μs) was taken with the same fault in place.
  - S6's trigger stays as restated above: the ratio is measured on the Mascotte setup before its run, and the film tables follow if the spray costs more than the gas there.
  - The Mac carried other users' jobs (load average 10 to 14) for the secant timing, so only its residual counts are quoted.

- **KH child release: an implementation fault, found on S7's first production run (14 October).** The code released a child each time a parcel's shed mass reached 3% of the parcel's *current* mass. So a 5 mm blob shrinking toward B0 Λ (about 1 μm at 300 m/s slip) kept releasing children all the way down: about 750 per blob (*derived*, 3 ln(a0/a_stable)/0.03). On the Mascotte blob run at 239×70, parcels went from 767 (0.9 ms) to 28,873 (1.0 ms) to 326,061 (1.1 ms), and that 0.1 ms took 439 s, 405 s of it in the spray. The run was stopped there.
  - **Fixed** to the published rule (OpenFOAM's ReitzKHRT, read in its source): a child is released when the parcel's shed mass reaches 3% of its injected mass (`Parcel::injectedMass`, set at injection or `add`, inherited by children). That bounds the children per parcel at about 33, and a child, holding 3% of the injected mass, makes no grandchildren.
  - **After the fix** (same run, same steps: 15,954 to 1.1 ms in both): 630 parcels at 1.1 ms, and that 0.1 ms took 29 s, 1 s of it in the spray. S0 to S8 rerun: all as before except S5(b)'s child count (17, was 207), which passes; S0(b) fails as before.
  - S3, S4 and S8 have no breakup, or breakup that releases few children, so their numbers are unchanged (the rerun prints the same values). S5(b) tested the mass carried through each child and the child's size, not the release rule. No criterion tested the release rule, which is how this got through.

## Spray on AMR (stated 22 October 2026, about 22:50 CDT, before its code)

Tier 1 needs AMR's levels and the liquid LOX together (`docs/validation/rl10a33a/tier1_run_definition.md`, *What S7 shows for tier 1*, item 4). `Spray` takes one `Flow` today. This section is the coupling's design and its criteria, registered under the Director's decision of 22 October (22:40): per-level parcels, or parcels on the finest level with the sources restricted, at the worker's choice with reasons; a conservation check across regrids; and a test that a single-level AMR run reproduces the uniform-mesh spray result to a stated tolerance.

**The choice: per-level stepping.** There is one parcel list. Each parcel is owned by the finest active level at its position. It is moved in its owner's step with that level's dt, reads the gas of its owner cell, and deposits into the owned cell where each sub-step starts.

Why, against the finest level with restriction:
1. **Time.** Each level's gas receives its sources inside its own step, after the flow and before its reaction. That is before the average-down and reflux that sync it to the level below. So per level, the order is the uniform one: flow, spray, react (`adapters/reacting_flow.cpp:286-291`). On one level the method is the uniform algorithm, operation for operation, and that is what lets the reproduction test be bit for bit.
2. **Cost (derived).** A parcel is moved 2^l times per base step, like its owner's gas, not 2^L times.
   - On tier 1's 3 levels, the finest-level option moves a base-level parcel 8 times per base step, where once suffices.
   - The spray was 33% to 41% of S7's wall time on one level (measured, `MASCOTTE_A10_PREDICTION.md` §4), so that factor matters.
3. **The one-cell rule still holds.** A base step is the smallest of 2^l times each level's stable step over its active cells (`core/amr.cpp`, `Amr::step`). A level-l step therefore moves a parcel less than one of its owner's cells, its speed being below the gas's u + c (as on one level).

Against it: each parcel needs its own clock within the base step, and a redone base step needs a snapshot of the spray. Both are below. The finest-level option needs neither.

**The design.**
- **Location.** A parcel is located on the finest mesh (its dz, its wall polyline and its ring fractions). Its cell on level l is the finest cell's indices shifted right by L − l.
  - This is exact, because a coarse cell is the exact union of its children (`Mesh::coarsened`): the coarse ring fractions are the fine ones at every other boundary, and dz doubles.
  - The owner is the highest level whose cell there is active. A composite cell is (level, cell), numbered offset_l + q. Every per-cell array of the spray is composite: the gas seen, the loads, the occupied flags, the cell sizes, the volumes and the turbulence.
  - The walls (the side-wall mirror and the outlet) use the finest mesh. On one level that is the uniform mesh itself.
- **The clock.** L is the number of levels less one. Each parcel holds its phase: the number of finest sub-intervals of the current base step it has crossed, from 0 to 2^L.
  - At the k-th step of level m within the base step, the end phase is k · 2^(L − m).
  - That step moves every parcel whose owner is at most m and whose phase is below the end phase, by the level's dt, and leaves it at the end phase.
  - *Derived:* every parcel the step moves is at phase exactly (k − 1) · 2^(L − m), so all of them move by the same dt.
    - The coarser levels step first. The last step of the finest level at the previous sync point brought every parcel to that point.
    - A parcel that crosses into a coarser owner's cells mid-cycle catches up at the next step of any level at or above that owner.
    - A parcel that crosses into a finer owner's cells waits until that level's clock reaches it.
  - The code checks the derived property and throws if it fails.
- **Release.** The injectors release at the base level's step, after its moves, at phase 2^L (the uniform order: release at the end of `advance`).
- **The occupied set.** The cells whose gas is gathered for a step are every owned cell overlapping the 3 × 3 block of owner-level cells around each moved parcel's owner cell.
  - A covered neighbour is resolved into its owned children.
  - An inactive neighbour is resolved to the owned cell of the level below that holds it.
  - On one level this is the uniform set. Owned cells of finer levels are included, because a parcel can enter them during the step.
- **The gas a parcel sees** is its owner cell's, at that level's current state. A coarser owner can be ahead in time by at most its own dt. This is a first-order lag, like the uniform method's split.
- **Deposits.** Each level's deposits go through `Flow::exchange` on that level, into owned cells only. A deposit is booked at the cell where its sub-step started, on whatever level owns that cell.
  - `Amr::budget` already subtracts every level's exchange sums.
  - The average-down after the finer steps carries the fine deposits into the base's composite content.
  - A deposit into a finer level's cell enters before that level steps. A deposit into a coarser level's cell enters after that level's step and is held to its next one. Neither is counted twice.
- **The point-implicit loads** are summed per composite cell, with that cell's own volume.
- **Regrids** change no spray state. Every owner is found again at every step, and every parcel's cell is re-found then. The only per-cell value carried over is the temperature guess of the gas inversion, which is a starting point, not a value.
- **A redone base step.** On a non-admissible stage, `Amr` restores its levels and redoes the base step at half the size.
  - The spray snapshots its parcels, its ledger, the injectors' pending stores and random streams, and its id counter at the base level's step.
  - `Amr` gains a redo hook (`setRedo`), and the spray restores its snapshot through it.
  - `Amr`'s saved levels gain the exchange sums (`exchangedMass_`, `exchangedEnergy_`, `exchangedMomentum_`), which `Amr::save` does not hold today. Without them, any after-step exchange in a redone step would be counted twice.
- **Interface.**
  - `Spray(Amr&, LiquidOxygen, Settings)`.
  - `afterStep(level, dt)`, which the driver calls from `Amr::setAfterStep` before that level's reaction.
  - `redo()`, for `Amr::setRedo`.
  - The driver bounds the base step by `nextEvent`, as `ReactingFlow::step` does.
  - `Spray(Flow&)` and `advance(dt)` stay as they are, and so does their arithmetic, bit for bit.
  - On AMR, `wallThrust` is the last complete base step's wall impulse over its dt, set at the base level's next step.
  - The ledger gains the parcels' moves between levels: counts of moves to a finer owner and to a coarser one, by motion within a step. Regrid moves are not counted in these.

**Criteria (stated before the code).**
- **A0. The uniform path is unchanged.** S3's run (0.5 ms) and S8(c)'s (0.2 ms, SST and dispersion) give the same final gas state, partial densities and parcel list, `cmp`-identical, on the build before this change and after it. The whole `crucible_spray_tests` passes as before.
- **A1. One level reproduces the uniform mesh, bit for bit.** The reference is S3's chamber with the spray on a `Flow` on the AMR base's mesh (`Mesh::coarsened` of the finest). The run is the same chamber on `Amr(def, 1)` with level 1 left inactive, the spray through `afterStep` and `redo`.
  - The final gas state, the partial densities and the whole parcel list (every field, S4's comparison) must be identical to the last bit:
    - (a) S3's run, frozen, to 0.5 ms;
    - (b) S8(c)'s, SST with dispersion, to 0.2 ms.
  - The tolerance is zero because the precondition holds. Measured at about 22:45, before this record: the gas alone, with no spray, is bit-identical between a single-level `Amr` and a `Flow` on the same mesh. The probe was S3's chamber over 1,083 steps to 0.1 ms, every step's dt equal.
  - If the precondition stops holding, A1 fails. The bar is not loosened.
  - Printed beside, not judged: the largest relative difference of the conserved state from the same run on the directly built mesh.
- **A2. Conservation across regrids.**
  - **The run.** S3's chamber on a 16 × 6 base plus 2 levels (64 × 24 at the finest). A custom tag holds the injector corner (z < 2 mm, r < 3 mm). A second custom tag marks a band, 3 mm half-width on the base and 1.5 mm on level 1, that sweeps down the chamber and back within 0.5 ms. The pattern is regridded every 2 base steps; the run goes to 0.5 ms on 1 thread.
  - **(a) The composite gas budget** (`Amr::budget`, with the exchange) must close within 1e-11:
    - mass, over the gas's initial mass plus the liquid injected;
    - energy, over |the gas's initial energy| + |injected| + |outflow| + the wall loss;
    - axial momentum, over the levels' gross momentum plus |injected| + |outflow| + |the wall impulse|.
  - **(b) The liquid's ledgers** must close within 1e-11 (`balanceErrors`).
  - **(c) The whole system**, gas plus liquid with no exchange terms, must close within 1e-11 on the same scales.
  - **The run must exercise what it checks.** It needs at least 50 regrids that changed the pattern, and at least 100 parcel moves to a finer owner and 100 to a coarser one. Otherwise it is resized and rerun before it counts, and it is not a pass.
  - Redone base steps are printed. They are not forced.
- **A3. Threads.** A2's run on 1 and 4 threads gives every level's state and the parcel list bit for bit.
- **Printed, not judged:** A2's liquid evaporated, parcels held and the spray's wall share, beside the uniform 16 × 6 run's.

## Literature

From a Sonnet agent's search (8 October, notes screened before use: no Mascotte A-10 outcome appears in them). Where a primary text could not be had, the source actually read is named.

- **Abramzon and Sirignano (1989)**, read in full by the agent from the published pages: B_M, B_T, F(B), Sh*, Nu*, the iteration between B_T and B_M, and the 1/3 reference state (their citation: Hubbard, Denny and Mills 1975) are as in *The drop model* above. The agent confirmed the structure of the no-blowing correlations but not every exponent; Sh_0 = 1 + (1 + Re Sc)^(1/3) f(Re) with f = max(1, Re^0.077) (Clift, Grace and Weber's sphere correlation, which A-S use for Re up to 400) is as I recall the paper (*inferred*, to re-check against the page).
- **Drag.** The agent could not read A-S's drag form exactly. Used: C_D = (24/Re)(1 + 0.2 Re^0.63)/(1 + B_T)^0.2 (Renksizbulut and Yuen's evaporating-sphere fit, as recalled, *inferred*), with Re on the film viscosity and the gas density. Its error goes into the momentum relaxation time only; the budgets do not depend on it.
- **KH-RT.** Reitz (1987), Patterson and Reitz (1998) and Beale and Reitz (1999) are paywalled; the equations are as reproduced by arXiv:2404.03170 (which cites Reitz 1987 for KH and Taylor 1950 for RT) and arXiv:2307.03722 (the switch rule). With a the drop radius, U the slip speed, We_g = ρ_g U² a/σ, We_l = ρ_l U² a/σ, Re_l = ρ_l U a/μ_l, Z = √We_l/Re_l, T = Z √We_g:
  - Λ/a = 9.02 (1 + 0.45 Z^0.5)(1 + 0.4 T^0.7)/(1 + 0.87 We_g^1.67)^0.6;
  - Ω (ρ_l a³/σ)^0.5 = (0.34 + 0.38 We_g^1.5)/((1 + Z)(1 + 1.4 T^0.6));
  - stable radius B0 Λ, da/dt = −(a − B0 Λ)/τ_KH, τ_KH = 3.726 B1 a/(Λ Ω);
  - RT, with g_t = (3/8) C_D ρ_g U²/(ρ_l a) the drop's drag deceleration: K = √(g_t (ρ_l − ρ_g)/(3σ)), Λ_RT = 2π C_RT/K, Ω_RT² = (2/(3√(3σ))) (g_t (ρ_l − ρ_g))^1.5/(ρ_l + ρ_g) (the maximum of ω² = k g Δρ/Σρ − σk³/Σρ, *derived*, which settles the source's ambiguous √3σ), τ_RT = C_τ/Ω_RT;
  - the switch (arXiv:2307.03722): RT acts when the drop is larger than Λ_RT and τ_RT < τ_KH; otherwise KH. RT acts only beyond the breakup length C_bl √(ρ_l/ρ_g) d_0 (Beale and Reitz's practice, *cited* through the secondary sources' description, constant C_bl = 1 *guessed*).
  - Constants. B0 = 0.61 (universal in the sources). B1 is case-dependent: 1.73 to 60 in the literature, 5 and 10 in the two sources; C_τ 0.1 to 5.33. None was calibrated on a cryogenic or rocket injector. Defaults here: B0 0.61, B1 10, C_RT 0.1, C_τ 1, *guessed* within the cited ranges, and B1 is varied across its range in the Blob run's sensitivity rather than tuned.
  - Child drops: no published magnitude for a child's normal velocity was found. A child keeps its parent's velocity (momentum exact, nothing invented); the spray's spreading then comes from the injection cone and the gas.
- **Coaxial drop size.** No closed-form coaxial SMD correlation was found. NASA SP-8089 (1976) gives none for the coaxial element and points to Burick (1972, NASA CR-120936), whose correlation is a curve on a log-log plot (Fig. 38): mass-median diameter over post diameter, from about 0.3 to about 0.1, against (V_g − V_L)/(V_L MR) from 0.4 to about 12, fitted to wax/GN2 cold-flow data with Rosin-Rammler distributions (D_30 = 0.455 D̄ for those). Mayer (1994) and Varga, Lasheras and Hopfinger (2003) were not found with explicit formulas. **The bracket is therefore restated:** the Distribution runs take Burick's curve, digitised from Fig. 38 with the scatter of its data, at the case's correlating group: the scatter band's low and high diameters give the two runs. Burick's data are simulant data unconnected to any validation case. The Blob run with KH-RT stays the third. (Restated 8 October, before any spray run.)
- **Liquid core length** (via a University of Maryland thesis's compilation, originals not read): Lasheras and Hopfinger (2000) L_b/D_l ≈ 6 M^−0.5 for M < 35 (M the gas-to-liquid momentum flux ratio), and Porcheron et al. (2002) L_b/D_l ≈ 2.85 (ρ_g/ρ_l)^−0.38 Oh^0.34 M^−0.13, the one fitted to cryogenic LOX coaxial jets. The Distribution runs release along the Porcheron length; Lasheras-Hopfinger's is reported beside it.
- **A resolved single LOX drop.** Lafon, Meng, Yang and Habiballah (2008, Combust. Sci. Technol. 180) compute a 100 μm, 90 K LOX drop in quiescent hydrogen at 1,000 K and 10 atm, subcritical, with resolved transient diffusion; results are figures only (their Figs. 4, 9 and 11). Added as criterion S2(c).
- **Surface tension.** NIST's table carries it. The agent's independent check: Mulero, Cachadiña and Parra (2012), σ = 0.03843 (1 − T/154.581)^1.225 N/m, about 13.2 mN/m at 90 K.
