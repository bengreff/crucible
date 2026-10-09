# Liquid oxygen spray (TECHNICAL_PLAN step 5)

Status (8 October 2026, about 23:40 CDT): design and criteria, stated before any spray code is run, with the literature filled in. Nothing here is verified yet.

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
  - The drop: dm/dt = −ṁ and m dh_l/dt = Q − ṁ L.
  - **Boiling and flashing.** If T_s would pass the saturation temperature at the cell's pressure, the drop is held there and ṁ = Q/L with the classical boiling form (B_T = c_pF (T_∞ − T_s)/L). A drop that finds itself above saturation (it moved to a lower pressure, as in the nozzle) evaporates the enthalpy above saturation at once (a flash, counted).
- **Drag.** A sphere's drag with the evaporating-drop correction. Form and constants: *Literature*, below.
- **Breakup**: KH-RT, Reitz (1987) wave model with the Rayleigh-Taylor branch of Patterson and Reitz (1998) and Beale and Reitz (1999). Equations and constants: *Literature*, below.
  - KH: each drop's radius shrinks toward B0 Λ at the rate set by B1; the mass shed accumulates on the parcel and is released as a child parcel of drops of radius B0 Λ once it reaches 3% of the parcel (*cited*, through arXiv:2404.03170), with the parent's velocity.
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
- **Mesh.** The cell is found from z (uniform dz), the wall radius linear between stations, and a binary search on the ring fractions. This is the one place the spray knows the mesh layout; AMR (inquiry2) will replace it.

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
- **S2. A single drop against a reference solution.**
  - **The reference.** The steady spherically symmetric evaporation of a drop in stagnant gas, solved resolved in r with variable properties: Cantera's h2o2.yaml with mixture-averaged transport and a correction velocity. Independent Python, `tools/drop_film_1d.py`, with the same surface condition (Raoult with NIST's vapour pressure) and the same latent heat. It needs no 1/3 rule, so it tests the model's film-property approximation in a hydrogen film, which is far from the hydrocarbon-in-air cases the rule was built on.
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
- **S7. Mascotte A-10.** Validation, pre-registered in its own record before any measured value is read: the case's modelling choices (equivalent round chamber, igniter, wall condition, the bracket's correlations), the predicted quantities and their bands. Requires turbulent dispersion first.

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
