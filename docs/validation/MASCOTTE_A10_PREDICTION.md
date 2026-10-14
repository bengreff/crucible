# Mascotte A-10: pre-registered prediction (SPRAY.md, S7)

**Status.** Modelling choices, judging rules and the drop-size extrapolation rule: 13 October 2026, committed before any production run and before any measured value of this case is digitised. Only development runs had been made (cost and stability; no profile was read from them). The bracket's diameters (section 3) are added once Burick's Fig. 38 is digitised, by the rule stated here. The predicted values (section 6) are added from the production runs and committed before any measured curve is digitised. Changes after this commit are labelled as changes.

Inputs come from [MASCOTTE_A10_INPUTS.md](MASCOTTE_A10_INPUTS.md), with its sources [N], [R] and [F]. The driver is `tests/mascotte_a10.cpp` (`crucible_mascotte_a10`). Labels: *stated* (a source gives it), *derived* (my arithmetic), *guessed* (a choice without data), *inferred* (reasoning, not measured).

## 1. What is predicted

- **P1, temperature.** Time-averaged gas temperature against radius at x = 10, 250 and 410 mm downstream of the injector exit, for r from 0 to 25 mm. These are the stations of [N]'s CARS temperature-against-height profiles.
  - The model's r is compared with the square chamber's height through the axis.
  - Judged over |height| ≤ 20 mm; r > 20 mm is reported only (section 5).
- **P2, flame shape.** From the time-averaged OH field:
  - the axial position of the peak of the radially integrated OH per unit length;
  - the axial position where it falls to half that peak downstream;
  - the radius of the OH maximum at x = 50, 100 and 200 mm.
  - These are compared with the Abel-transformed OH* emission image of [R] (and [N]'s OH comparison) as a **proxy**: OH* chemiluminescence is not ground-state OH (*inferred*: OH* marks the heat-release zone, while OH is broader and lasts longer downstream). So P2 is judged with a wider band (section 5).
- **Not predicted:** chamber pressure. It is imposed (section 2), so the case says nothing about c*.

## 2. Modelling choices

**Chamber**
- **Shape.** The square 50 × 50 mm section becomes a round chamber of radius 28.0 mm. Its area is 2,463 mm², 1.5% below the square's 2,500 mm² (*derived*).
  - Why not the equal-area 28.21 mm: 70 rings of 0.4 mm put ring faces exactly on the annulus radii 2.8 and 6.0 mm. That keeps the hydrogen area, and so its velocity and J, exact.
  - The corners of the square, which hold about 21% of its area beyond the inscribed circle, are not represented (*derived*).
- **Length and the protrusion.** The domain runs from the injector exit plane to 478 mm (*stated*, [N], [R]: injector exit to throat).
  - The injector protrudes 10 mm into the chamber (*stated*, [N]). The 10 mm deep space between the faceplate and the injector exit plane, around the injector body, is not modelled: at z = 0 the plane outside r = 6 mm is a wall.
- **Outlet.** The convergent (its shape is unknown: a square-to-round transition, [F]) and the 15 mm throat are not modelled, because the throat's small rings would cut the time step by about three (*derived*).
  - Instead the constant-area chamber ends in a partially non-reflecting subsonic outlet at 1.00 MPa (relaxation 0.5). 1.00 MPa is the recorded operating pressure ([F] Table 1, labelled there as an operating condition).
  - So the chamber pressure is imposed, not predicted. The last station (410 mm) is 68 mm upstream of the outlet.

**Walls and window cooling**
- **Walls.** Adiabatic and no-slip on the liner and the face, as in [R]'s simulation; the real liner's cooling is unknown.
  - SST runs to the wall without wall functions. The engine's wall functions need an isothermal wall, and none is known.
  - The near-wall layer is under-resolved (*inferred*: the first ring is 0.4 mm thick).
- **Helium window film.** Omitted: helium is not in the medium.
  - It is 6 to 10 g/s (*stated*, [F] and [N]; the two campaigns differ), 8 to 14% of the propellant flow, along two of the four walls.
  - It cools and dilutes the gas near those walls (*inferred*). This is why r > 20 mm is not judged.

**Propellants**
- **LOX.** 50 g/s at 85 K (*stated*, nominal), through the 5 mm post (r < 2.5 mm). Its density is 1,166.6 kg/m³ (NIST at 85 K), so the mean velocity is 2.183 m/s (*derived*).
- **H2.** 23.7 g/s at 287 K total (*stated*), through the annulus from 2.8 to 6.0 mm (88.47 mm²). Its density is 0.8449 kg/m³ at 1 MPa, so the velocity is 317.1 m/s (*derived*).
  - The supply's inflow turbulence is intensity 0.05 with viscosity ratio 10 (*guessed*, as in the engine's other chamber runs).
- **Ratios.** MR 2.11. J is 15.28 (*derived*; the sources give 12.63, 15.4 and 12.99). The velocity ratio is 145.3 (*derived*; [R] gives 146.3).
- **The face.** Over the post (r < 2.4 mm in rings) and the 0.4 mm lip (the real lip is 0.3 mm, *derived* in the inputs record), the face is a wall for the gas. The liquid enters through it as parcels.

**Start, chemistry and turbulence**
- **Start and ignition.** No igniter.
  - The chamber starts at rest at 1 MPa, filled with the overall mixture's complete-combustion products at 2,000 K (*guessed*): Y_H2O 0.763 and Y_H2 0.237 (*derived*).
  - Both valves ramp open over 1 ms. Local equilibrium burns any hydrogen and oxygen mixture at once, so the flame needs no ignition.
  - The start-up transient is not compared.
- **Chemistry.** Local equilibrium of the h2o2 mechanism's species, from a Table A built for this case, `crucible_build_equilibrium_table build/table_mascotte.bin 20 64 1e-4 129 97 1e-3 150 4`.
  - Its range is ρ_r 1e-3 to 150 kg/m³. Cold oxygen vapour at 1 MPa reaches about 40 kg/m³ near 100 K, past Table A's 30 kg/m³ ceiling, which would send those cells to Cantera's slow path.
  - Its check against direct Cantera over ρ 0.01 to 150 kg/m³: the 99th percentile of |dT| is 1.91 K, against TABLE_A's 3 K limit ([log](../evidence/mascotte/table_build_check.txt)).
  - Mean-state equilibrium, with no turbulence-chemistry interaction (the engine's PaSR closure needs finite-rate chemistry). Flame anchoring at the lip is therefore not a prediction: equilibrium anchors it there by construction (*inferred*).
- **Turbulence.** SST with the Spalart-Rumsey ambient on the fill. Low-Mach flux option: none (the engine's default).

**Spray**
- **Model.** Lagrangian parcels, two-way coupled. Turbulent dispersion is on. Wall restitution is 0.5 (*guessed*). Parcels are released at 5e5 per second (*guessed*, for cost).
- **Bias carried in from S2(c).** Against Lafon et al.'s resolved real-fluid drop in hot hydrogen at 10 atm, the drop model evaporates about 30% faster (SPRAY.md, S2(c)). If that carries over, the spray and the flame come out short (*inferred*). This is stated, not corrected.
- **Three runs make the bracket** (SPRAY.md, *Injection and the drop-size bracket*): Distribution low, Distribution high, and Blob. The Distribution runs draw Rosin-Rammler diameters at a Sauter mean diameter from Burick's correlation (section 3), released from the liquid core's surface along Porcheron's core length. KH-RT acts on every drop in all three runs.
- **Core length** (*derived*: NIST at 85 K, ρ_l 1,166.6 kg/m³, μ_l 2.25e-4 Pa s, σ 0.01446 N/m; ρ_g 0.8449 kg/m³; D 5 mm; M = J = 15.28; Oh = 7.75e-4):
  - **Porcheron et al. (2002):** L/D = 2.85 (ρ_g/ρ_l)^−0.38 Oh^0.34 M^−0.13 = 2.73, so **13.7 mm**. Used.
  - **Lasheras and Hopfinger (2000):** 6 M^−0.5 = 1.53, so 7.7 mm. Reported beside it.
  - **Flag.** Porcheron et al. fitted cryogenic coaxial jets and may include this facility's own jets. Their correlation is for a core length, which is not a predicted quantity, so it is kept as SPRAY.md states.

**Mesh and time**
- **Production mesh.** 239 × 70: dz 2.0 mm, dr 0.4 mm, 8 rings across the annulus and 6 across the post. CFL 0.8, the engine's ceiling, provided the 1 ms stability test holds (section 7).
- **Grid check.** The Blob run is repeated at 478 × 70 (dz 1.0 mm). The radial resolution is not varied, so its error is not measured; that is a stated limit.
- **Time.** Each run goes to 25 ms. That is about four residence times of the hot gas, about 6 ms each (*derived*: 0.44 g in the chamber against 73.7 g/s).
  - The time average runs from 15 to 25 ms, sampled every 0.05 ms.
  - **Settled** means the two halves' averages (15 to 20 ms and 20 to 25 ms) differ by at most 50 K at every station point with r ≤ 20 mm. A run that is not settled is extended by 10 ms, and the average window moves with it (stated now).

## 3. The drop-size bracket from Burick (1972), and the extrapolation rule

The correlating group, (V_g − V_L)/(V_L MR), is (317.1 − 2.18)/(2.183 × 2.11) = **68.4** (*derived*). Burick's data span about 0.4 to 12 (SPRAY.md, *Literature*), so the case lies about six times beyond his range. Burick's data are simulant (wax and GN2), unconnected to this case. The rule below is stated before the digitisation is read:

- **The diameter.**
  - Take Burick's ordinate as printed, with any property correction his report gives for scaling to a real propellant. With none, use it as is.
  - The ordinate is the mass-median diameter over the post diameter (5 mm).
  - Convert it to the Sauter mean with his Rosin-Rammler exponent q: SMD = X / Γ(1 − 1/q), with X = MMD / (ln 2)^(1/q).
- **High diameter.** The upper edge of the data scatter, at the largest group with data, held constant beyond it. The drop size stops falling.
- **Low diameter.** A power law through the lower edge of the scatter, at the slope of the correlation curve over its last half decade of data (group 4 to 12), continued to 68.4.
- **Exponent.** The Rosin-Rammler exponent is Burick's. If he gives a range, use the middle of it for both runs.
- **Assessment (*inferred*).** Extrapolating a cold-flow wax correlation by six times in its group is weak ground. The band it gives is not a statement of confidence in either end. The Blob run, which does not use Burick, is the third leg.

**Applied (14 October 2026, 00:30 CDT, before any production result).** The digitisation is `docs/evidence/spray/burick_fig38.md` (a Sonnet agent, NASA CR-120936 Fig. 38, p. 67; cold-flow content only).
- **Axes.** The ordinate is D̄/D_L, the mass-median diameter over the liquid post diameter (5 mm here). The abscissa is the group above. Data markers span x ≈ 0.3 to 11; the curve is drawn to about 15.
- **Property correction.** The report names one (Ingebo's relations, its Ref. 24, wax to FLOX) but does not give it, so by the rule the ordinate is used as printed. *Inferred* direction: LOX's viscosity is about 18 times below the wax's (0.225 against 4.0 cP), which tends toward smaller real drops. This is a declared limit of the bracket.
- **Exponent.** Burick prints no number, but his Eq. (5), D30 = 0.455 D̄ for a Rosin-Rammler spray, fixes it. For a mass-based Rosin-Rammler, D30/MMD = Γ(1 − 3/q)^(−1/3) / (ln 2)^(1/q), so **q = 3.207** (*derived*). ±0.01 on 0.455 moves q by ±0.015. Then SMD = 0.851 MMD.
- **High.** The upper edge of the scatter at the largest group with data (x ≈ 10 to 11) is 0.135. MMD 675 μm, **SMD 574 μm**.
- **Low.** The curve's slope over group 4 to 12 is −0.371 (least squares on the digitised curve; −0.323 from its end points). The power law runs through the lower edge at the largest group with data (0.105 at x 10.5), the same anchor as the high rule. At 68.4 it gives 0.0523: MMD 262 μm, **SMD 223 μm**.
- **Readings not taken, so the choice is visible.** Anchoring the low line on the agent's flagged low point (x 7.1, 0.052) would give an MMD of 112 μm. The curve divided by 1.7, the typical lower scatter, would give 164 μm. The last-group anchor is used because it mirrors the high rule and that point is flagged uncertain.
- **Runs.** `rr:223:3.207:13.7` (low) and `rr:574:3.207:13.7` (high).
- **Assessment (*inferred*).** Burick's own section 10.3 warns that this cold-flow data "may not be directly applicable to LOX/GH2", whose gas velocities are far higher. KH-RT acts on the released drops at about 300 m/s slip, so the injected size sets the start of breakup, not the final drop size.

## 4. Runs

| Run | Size mode | Mesh | Purpose |
|---|---|---|---|
| low | Rosin-Rammler at the low SMD, Porcheron release | 239 × 70 | bracket |
| high | Rosin-Rammler at the high SMD, Porcheron release | 239 × 70 | bracket |
| blob | Blob (5 mm drops), KH-RT | 239 × 70 | third leg |
| blob-fine | Blob | 478 × 70 | axial grid check |

## 5. The band and the judging rules

- **The band.** At each station and radius, the band runs from the lowest to the highest of the three production runs' time-averaged T. It is widened on each side by:
  - the grid difference |T(blob-fine) − T(blob)|, and
  - the largest drift between the halves of the three runs.
- **P1 is judged per station.** Every measured point with |height| ≤ 20 mm is judged; if [N] plots both sides of the axis, both are used.
  - A point is **inside** if it lies within the band widened by the measurement's stated uncertainty. Where [N] states none: ±5% of the measured value (*guessed*).
  - A station **passes** if at least two thirds of its points are inside.
  - Also reported: the largest and the RMS difference from the band's middle; and the points with |height| > 20 mm, as reported only.
- **P2.** The measured image's axial peak passes if it lies within the range of the three runs' peak positions, widened by ±20% (*guessed*, for the OH*/OH proxy). The same rule applies to the half-maximum position, and to the OH radius where the image covers x = 50, 100 and 200 mm.
  - If the image's window does not contain the peak, P2 reports the shape inside the window, with no pass or fail.
- **Reading a miss.** A miss is read against the declared limits before any change: the S2(c) evaporation bias, equilibrium chemistry without turbulence-chemistry interaction, the omitted helium film, round for square, the imposed pressure and the drop-size extrapolation. No model constant is tuned to this case. Any rerun after the measured values are seen is labelled as a change and does not replace the blind result.

## 6. Predicted values

*Added from the production runs and committed before any measured curve is digitised.*

## 7. Development runs before this record (cost and stability only)

- **CFL 0.4 at 478 × 70, 6 threads on the Mac.** Step 3.1e-8 s at 2,000 K. 34 ms per step. That projects to about 9 h per 30 ms run, which is why the mesh, the CFL and the run length above were chosen.
- **CFL 0.8 at 478 × 70.** Step 6.2e-8 s, 42 to 46 ms per step. A 1 ms stability test is running when this is written. If it fails, the production CFL is 0.4 and this record says so.
  - **Result (14 October, 00:17 CDT): stable to 1 ms**, so production stays at CFL 0.8. 16,845 steps; the step fell from 6.2e-8 to 5.7e-8 s as the H2 jet came in; T_max 2,001 K; no table fallbacks; mass and energy budgets at most 2.9e-14 and 5.9e-14, the liquid ledger at most 2.3e-13. KH-RT breakup began at about 0.85 ms and the parcel count went from 354 (0.7 ms) to about 21,000 to 30,000 (0.9 to 1.0 ms). The test was a stability check only; the spray's cost once breakup is under way is measured on the production runs (`spray_wall_s`).
