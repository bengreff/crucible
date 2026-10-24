# Thermo below 200 K for H2, O2 and H2O: the sources, the hydrogen form, and the bit-for-bit test

Registered 24 October 2026, about 02:45 CDT, before the build and before any test has run. The code edits exist, uncommitted, on worktree-inquiry2. Asked for by the Director at 02:20: "the species' low-temperature thermo below 200 K (H2, O2, H2O; the 329 Cantera warnings near 199 K are this). Register first: the source fits (NASA-9 or equivalent, cited, with validity range), the normal vs equilibrium H2 choice for the supply (the RL10's H2 comes off the regenerative jacket warm; state which form applies and why), and the test that nothing above 200 K changes bit for bit."

## What happens below 200 K now

- **The flow (Medium).** Every species in `h2o2_highT.yaml` carries TM-4513's NASA-7 fits (McBride, Gordon and Reno 1993, NASA TM-4513, through Cantera's `nasa_gas.yaml`), valid 200 to 6000 K. Below 200 K, Medium evaluates the low range's polynomial outside its range, silently.
- **Cantera.** The same data inside Cantera warns below 200 K. The 329 warnings in tier 1's log (AMR.md, job rfx) read "ChemEquil::equilibrate: Temperature (about 199 K) outside valid range of 200 K to 6000 K". They come from a Cantera equilibrate in the tier-1 driver's chemistry, not from Medium.
- **How large the extrapolation's error is** (derived, against the sources below):
  - At 60 K: O2's cp/R is 3.636 against 3.500 (3.9% high), and H2O's 4.099 against 4.000 (2.5% high).
  - H2 below 200 K is the hydrogen-form question below, not a fit question.

## (a) The sources, with their validity

- **At and above 200 K: unchanged.** The TM-4513 fits stay as they are, valid 200 to 6000 K, with the constant-cp continuation above 6000 K as before.
- **No standard fit goes below 200 K.** NASA-9 (McBride, Zehe and Gordon 2002, NASA/TP-2002-211556) starts its lowest range at exactly 200 K for H2, O2 and H2O. The literature agent read that in the primary coefficient tables.
- **Burcat's database** (Burcat and Ruscic, ANL-05/20) allows segments below 200 K as a capability of its format. Whether H2, O2 or H2O carry one was not checked: the file did not download.
- **The equivalent used here** is the ideal gas's statistical mechanics with cited spectroscopic constants, the same physics JANAF and the NASA fits are built from. `tools/mechanisms/make_low_t.py` evaluates it, and its docstring holds the formulas:
  - **H2, as normal hydrogen** (the form is argued in (b)):
    - ortho and para frozen at 3:1, each with its own rotational partition function: even J para; odd J ortho, with nuclear-spin weight 3;
    - the v = 0 levels B0 J(J+1) − D0 J²(J+1)² + H0 J³(J+1)³, with B0 59.3392, D0 0.04599 and H0 5.2e-5 cm⁻¹ (Huber and Herzberg 1979, *Constants of Diatomic Molecules*);
    - vibration (5,987 K) neglected: under 1e-10 in cp/R at 200 K (derived).
  - **O2:**
    - a classical rotor (θ_rot 2.08 K) and a harmonic oscillator at G(1) − G(0) = 1556.385 cm⁻¹ (Huber and Herzberg 1979);
    - neglected: the ground state's spin splitting, about 3 K. In the high-temperature Schottky form it adds of order 0.2 (3 K / T)² to cp/R: about 0.005 at 20 K and under 5e-4 above 60 K (inferred).
  - **H2O:**
    - a classical asymmetric rotor and harmonic oscillators at the fundamentals 3657.053, 1594.746 and 3755.929 cm⁻¹ (Shimanouchi 1972, NSRDS-NBS 39);
    - the rotor is classical, so with rotational temperatures up to about 40 K its quantum correction grows below about 60 K (inferred: about 1% at 60 K, from the (θ/T)² term).
    - Water vapour below 200 K at the flow's pressures is supersaturated. The condensed phase is not modelled, so this range only keeps the gas-phase state defined and smooth.
- **The fits.** Each species gets a lowest NASA-7 range, 20 to 200 K:
  - cp/R is fitted by least squares on a 1 K grid from 20 to 200 K, constrained so that cp at 200 K equals TM-4513's;
  - a5 and a6 then make h and s continuous at 200 K;
  - below 20 K the polynomial is extrapolated (it only has to stay positive there: the energy floor and the inversion pass through it).
- **Fit errors** (derived, by `make_low_t.py`):

  | Species | Largest \|cp_fit / source − 1\|, 20 to 200 K | Least cp/R, 0 to 200 K | Jumps at 200 K (cp, h, s over R) |
  |---|---|---|---|
  | H2 | 1.28e-2 | 2.4845 | 0, 0, 0 |
  | O2 | 8.63e-4 | 3.4996 | 0, 0, 0 |
  | H2O | 2.43e-3 | 3.9986 | −8.9e-16, 0, 0 |

  H2's error is the largest because the fit must meet TM-4513's cp at 200 K (3.3108). That is equilibrium hydrogen's, 0.9% above normal's 3.2797 (derived; see (b)).
- **Validity claimed:** ideal gas, 20 to 200 K, to the fit errors plus the source's neglects above. That makes H2 about 1.3%, O2 about 0.1% and H2O about 1% above 60 K.
- **Independent checks** (test 6 below):
  - H2 against normal hydrogen's ideal-gas Cp0 from the NIST WebBook (Leachman et al. 2009, J. Phys. Chem. Ref. Data 38, 721; the raw tables were read by the literature agent);
  - O2 and H2O against JANAF (Chase 1998, NIST-JANAF 4th edition), tables O-029 and H-064, which I fetched as text tonight.

## (b) The hydrogen form for the supply

- **What the engine's H2 is above 200 K: equilibrium hydrogen.**
  - TM-4513's H2 at 200 K has cp/R 3.3108. Equilibrium hydrogen (the full partition function, conversion heat included) has 3.3012, normal 3.2797 and para 3.8962 (derived).
  - Δh from 200 to 298 K: the fit 1376.79 kJ/kg, equilibrium 1376.25 kJ/kg (derived).
  - JANAF's H-050 at 200 K gives cp/R 3.3011, which is equilibrium too (derived from the table). At 100 K, JANAF's 3.386 is far from normal's 2.715 (the WebBook's Cp0), so JANAF's hydrogen is the equilibrium form.
- **What the RL10's supply actually is: frozen para-rich hydrogen, neither normal nor equilibrium.**
  - The tank holds liquid hydrogen. Liquid hydrogen is at least 95.0% para by volume under MIL-PRF-27201 (Type II, Revs D and E, read directly), and equilibrium at 20 K is 99.8% para (Lee, Dominguez-Ontiveros and Williams, ICAPP 2023, OSTI 1973337).
  - The pump and the regenerative jacket warm it in milliseconds. Uncatalysed ortho-para conversion is far slower:
    - In the gas phase, the rate constant usually quoted from Milenko et al. (J. Low Temp. Phys. 107, 77, 1997; secondary, behind a paywall) gives 1e8 to 1e9 s at these densities (inferred).
    - At the walls, conversion needs paramagnetic sites (Ilisca, Prog. Surf. Sci. 41, 217, 1992; Fukutani and Sugimoto, Prog. Surf. Sci. 88, 279, 2013). No cited rate exists for the jacket's metal. That gap is stated, not closed.
  - So the gas reaches the injector warm and still at least 95% para. NASA's own hydrogen property code GASP (Hendricks, Baron and Peller 1975, TN D-7808) is a parahydrogen code for the same reason.
- **The size of the difference** (derived, para against equilibrium, by the supply's temperature):
  - cp: 18.1% higher at 199.67 K, 9.0% at 250 K, 3.7% at 300 K, 0.5% at 400 K and 0.07% at 500 K.
  - Enthalpy h_eq − h_para (common above about 1000 K): 163 kJ/kg at 199.67 K, 70 at 250 K, 27.6 at 300 K and 0.46 at 500 K.
  - The engine specifies the supply by its temperature, so with equilibrium thermo it gives the H2 that much more energy than the para-rich gas carries.
  - Once burning, the hydrogen converts quickly: H-atom exchange at flame temperatures, and collisions with paramagnetic O2 (inferred). So the error is an energy offset at the supply and a cp error in the cold mixing region, not a lasting change of the burnt gas.
- **Why it is not built tonight.**
  - The para form changes H2 between 200 K and about 500 K, so it fails the bit-for-bit test the Director asked for.
  - It needs Table B rebuilt.
  - Unless the mechanism's H2 thermo changes too, Cantera's kinetics thermo would disagree with the flow's.
  - It is the Director's decision. The two options are:
    - (i) para-rich H2 thermo over the whole range, as a frozen mixture of para and ortho species, or as one species with the para fraction's cp;
    - (ii) keep the equilibrium thermo and correct only the supply: give the face inlet h_para(T_s) instead of h_eq(T_s) for the H2 supply. That gets the energy right, and the cold region's cp stays equilibrium's.
  - My recommendation is (ii): one inlet change, no table rebuild, and it fixes the energy, which is what the chamber's performance sees.
- **Below 200 K tonight: normal hydrogen.** The gas that goes below 200 K in the flow got there from the engine's own state at or above 200 K, by expanding or cooling over milliseconds. So its form is frozen at that state's ortho-para ratio:
  - Equilibrium at 200 K is close to normal in composition. cp differs by 0.65% at 200 K (derived, 3.3012 against 3.2797).
  - **Normal** (frozen 3:1) is therefore the consistent continuation of the engine's thermo.
  - **Equilibrium** below 200 K would add conversion heat the gas cannot release in milliseconds: cp/R 4.36 at 60 K against normal's 2.52, 73% high (derived).
  - **Para** would be the supply's physical form, but it cannot join the equilibrium thermo above 200 K continuously: there is an 18% jump in cp and 163 kJ/kg in h at 200 K.
- **The band against para below 200 K** (the WebBook's Cp0, para over normal): 1.197 at 100 K, 1.287 at 150 K and 1.188 at 200 K. That band applies to supply gas which really is para-rich and goes below 200 K. It closes only with the decision above.

## (c) The tests, with tolerances, before the run

`tests/low_temperature_tests.cpp`, ctest `low_temperature_verification`. It compares two media: Medium(h2o2_highT's species), and the same species after `thermo::withLowTemperatureRanges`. The tolerances, as in the file's header:

1. **Nothing at or above 200 K changes, bit for bit.** Checked with memcmp, for:
   - the quantities: Medium's properties (r, e, cv), gas constant, internal energy, enthalpy, cv, sound speed and species enthalpies, and ManifoldRuntime::entropy;
   - the compositions: 64 of them (the 10 pure species, then random seeded ones);
   - the temperatures: 200 K exactly; 2,000 temperatures log-spaced from 200 to 20,000 K; and the 1,000 K and 6,000 K joins with the doubles on either side.
2. **The temperature inversion is bit for bit**, for every e at or above e(200 K). The guesses lie on both sides of 200 K: 20, 60, 150 and 199.99 K, the target, 0.7 and 1.3 times it, 300, 3,000 and 20,000 K, 0, and NaN.
3. **The joins at 200 K**, taken between 200 K and the double below it, per species (H2, O2, H2O): |cv jump| / cv ≤ 1e-14; |h jump| ≤ 1e-12 |h| + 1e-9 J/kg; |s jump| / s ≤ 1e-12.
4. **Below 200 K:**
   - the inversion round-trips to |T − T0| ≤ 1e-9 T0 from 20 to 199.99 K, for every composition, from guesses T0, 300 K and 20 K;
   - cv > 0 from 1 K;
   - the energy floor lies below e(20 K).
5. **An engine-level check.** A shock tube of H2, O2 and H2O (Y 0.2, 0.3, 0.5), 10 bar at 1,500 K against 1 bar at 300 K, 400 cells, run to 0.2 ms. Every state stays above 200 K, and the run with the lowest ranges must equal the run without them bit for bit: the number of steps, the conserved state and the partial densities.
   - The control puts the right side at 120 K. There the two runs must differ, which shows the check can tell them apart.
6. **The lowest ranges against independent tables** (cp/R, with R = 8.314462618 J/(mol K)):
   - H2 within 1.5% of the WebBook's normal-hydrogen Cp0 at 100, 150 and 199.9 K (the last against its 200 K entry);
   - O2 within 0.2% of JANAF O-029 at 100 K and of the WebBook's Cp0 at 150 K;
   - H2O within 0.5% of JANAF H-064 at 100 K;
   - and H2 more than 15% from JANAF's equilibrium H-050 at 100 K, which confirms the form is normal.
- **The fast ctest** must be otherwise unchanged. Nothing outside `withLowTemperatureRanges`'s callers sees a change, because Species' new fields default to tLow = 0.
- **The arithmetic that keeps it bit for bit:**
  - At t ≥ tLow every species' range selection takes the same branch as before.
  - The Newton iteration of `temperature()` evaluates an iterate below 200 K without the lowest ranges whenever the target e is at or above e(200 K). So no iterate, whatever the guess, takes a different path.
- **One engine-level difference that test 5 cannot see.** The energy floor (e at T = 0) moves for H2, O2 and H2O. So the admissibility of a Riemann star state with e between the old floor and the new one can differ. That is a state near 0 K, which no physical run reaches.

## Scope

- **Opt-in.** Nothing changes until a driver passes its species through `thermo::withLowTemperatureRanges`. Tier 1 would do it at `tests/rl10_tier1.cpp:897`: `base.species = thermo::withLowTemperatureRanges(source.medium().species());`.
- **The kinetics and Cantera's equilibrate keep the mechanism's own thermo.** So the 329 warnings remain, and Cantera's states below 200 K are still extrapolated. Moving them would need the mechanism's YAML to carry the ranges, and that changes Cantera's tables (Tables A and B).
- **Other species are unchanged.** H and AR are monatomic: their cp/R is exactly 2.5 and extrapolates exactly. The radicals and N2 keep their extrapolated low ranges below 200 K.

## Results (24 October 2026, about 03:32 CDT; Release build on the Mac, 0 warnings)

All six tests met, 0 failures (measured; `low_t/low_temperature_tests_2026-10-24.txt`). `make_low_t.py` rerun gives the 42 coefficients in `adapters/reaction.cpp` digit for digit (`low_t/make_low_t_output.txt`).

1. **At and above 200 K, bit for bit:** 128,448 composition-temperature pairs (64 compositions, 2,007 temperatures), and none differs in any property.
2. **The inversion, bit for bit:** 220,416 inversions of e at or above e(200 K), from the twelve guesses, and none differs.
3. **The joins at 200 K:**

   | Species | \|cv jump\| / cv | \|h jump\| [J/kg] | \|s jump\| / s |
   |---|---|---|---|
   | H2 | 0 | 2.3e-10 (tol 1.4e-6) | 2.9e-16 |
   | O2 | 0 | 2.9e-11 (tol 9.1e-8) | 0 |
   | H2O | 3.3e-16 | 0 (tol 1.4e-5) | 1.9e-16 |

4. **Below 200 K:**
   - the largest round-trip error is 4.0e-13 of T0;
   - the least cv from 1 to 200 K over all 64 compositions is 312 J/(kg K);
   - e(20 K) is at least 6,244 J/kg above the floor.
5. **The shock tube:**
   - At 300 K it took 919 steps both ways, with the state and partial densities bit for bit.
   - The control at 120 K took 843 steps both ways, and the states differ.
6. **Against the tables, |cp / table − 1|:**

   | Species | Point | Table | Error | Tolerance |
   |---|---|---|---|---|
   | H2 | 100 K | NIST WebBook, normal H2 | 2.0e-3 | 1.5e-2 |
   | H2 | 150 K | NIST WebBook, normal H2 | 2.6e-3 | 1.5e-2 |
   | H2 | 199.9 K (against 200 K) | NIST WebBook, normal H2 | 9.0e-3 | 1.5e-2 |
   | O2 | 100 K | JANAF O-029 | 9.7e-5 | 2e-3 |
   | O2 | 150 K | NIST WebBook | 5.0e-4 | 2e-3 |
   | H2O | 100 K | JANAF H-064 | 1.0e-3 | 5e-3 |

   - H2 lies 19.7% from JANAF's equilibrium H-050 at 100 K, so the form is normal.
   - H2's 0.9% near 200 K is the join to TM-4513's equilibrium cp. That is the 0.9% between the forms in (b), as registered.
- **The fast ctest** (`ctest -LE slow -j4`, full rebuild, 0 warnings) gives the same 4 known failures, and every FAIL line is identical to `ctest_2026-10-24_skip_nan.txt`'s as sorted sets (14 lines, values included). This held both before the rebase (0026427, `ctest_2026-10-24_low_t.txt`) and after it, onto main 96635e8 (622f5c7, `ctest_2026-10-24_low_t_rebased.txt`). `low_temperature_verification` passes in both. A second rebase, onto main 695df32 (docs only), renamed 622f5c7 to aa5abf4 and the code commit 72c5c38 to a83459c.
- **Seen, not judged:** spray_verification's known failure "(b) |h_medium − h_NIST| on the engine's zero" (1.245e3 J/kg against 1e3) evaluates Medium below 200 K for oxygen, on the extrapolated range. That test does not opt in, so it is unchanged. Whether the lowest ranges move it was not run.

## (d) The supply's hydrogen: the Director's decision, and the bound before any conversion model

**Registered 24 October 2026, about 04:04 CDT by the Mac clock, before the computation or any run.**

**The Director's decision (24 October, 04:05).** The supply is para-rich hydrogen.
- Liquid hydrogen in the tank is about 99.8% para at 20 K.
- Uncatalysed ortho-para conversion takes hours, far longer than the jacket's residence time. So the H2 leaves the jacket at about 200 K with the tank's para fraction frozen.
- In the hot zone, H-atom exchange converts it quickly.
- Before any conversion model is built, the effect is quantified, as below.
  - If it is below 0.02% (a tenth of Table B's 0.2% bar), the supply is declared para, with instantaneous equilibration above a stated temperature, and the bound is recorded.
  - If it is larger, the form needs a sourced conversion rate.
- The LH2 para fraction and the conversion time are to be cited.

**What is computed** (`tools/mechanisms/para_supply_bound.py`, CEA-style bookkeeping; light, seconds on one thread).
- **The supply's enthalpy difference.** Δh(T_s, x_p) = h_eq(T_s) − h_s(T_s, x_p), per kg of H2.
  - h_s is a frozen mixture with para fraction x_p. h_eq is ortho and para in equilibrium at T_s.
  - Both come from the rotational partition functions on make_low_t.py's levels (Huber and Herzberg), with a common zero at para's J = 0. Translation and vibration are the same in both forms and cancel.
  - Consistency check: x_p = 1 at 199.67 K must give section (b)'s 163 kJ/kg.
- **The supply states.**
  - H2 gas at T_s = 150, 200, 250 and 300 K. Its equilibrium enthalpy is the mechanism's H2 (TM-4513, equilibrium) at and above 200 K. At 150 K it is the mechanism's value at 200 K plus the equilibrium partition function's change from 200 to 150 K.
  - x_p = 0.95 (MIL-PRF-27201's minimum) and 0.998 (equilibrium at 20 K).
  - O2 in two states, to show its weight on the result:
    - liquid at 90.17 K, h = −12.979 kJ/mol (CEA's O2(L), McBride, Zehe and Gordon, NASA TP-2002-211556's data);
    - gas at 298.15 K, h = 0.
- **The heat release.** q = h_0 − h_298, per kg of propellant.
  - h_0 = Y_H2 h_eq(T_s) + Y_O2 h_O2.
  - h_298 is the enthalpy at 298.15 K of the complete-combustion products (H2O gas plus the excess H2), on the mechanism's thermo.
  - The fraction is f = Y_H2 Δh / q.
- **c\* and Isp,** with the base supply (h_0) and the para supply (h_0 − Y_H2 Δh), on h2o2_highT.yaml's ten species:
  - The chamber is an HP equilibrium at p_c (an infinite-area combustor).
  - The expansion is shifting-equilibrium and isentropic, at s_c, one-dimensional: u = sqrt(2(h_c − h(s_c, p))) and G = ρu.
  - The throat is the maximum of G over p, found by golden section to 1e-10 in p. Then c\* = p_c / G_t.
  - Vacuum Isp at area ratios 20, 60 and 100, on the supersonic branch: u_e + p_e / (ρ_e u_e).
  - Only the relative changes Δc\*/c\* and ΔIsp/Isp are reported.
- **The grid.** O/F 4, 5 and 6; p_c 2, 3 and 4 MPa; the four T_s; the two x_p; the two O2 states. That is 144 cases.

**How it is judged, by the Director's rule.**
- The judgement is made at the Director's stated supply, T_s = 200 K, over the rest of the grid (36 cases).
- If the largest of f, |Δc\*/c\*| and |ΔIsp/Isp| there is below 2e-4: para at the supply, with instantaneous equilibration above a stated temperature, and the bound recorded.
- Otherwise: a sourced conversion rate is needed before the form is fixed.
- T_s = 150, 250 and 300 K are reported beside it to show how the bound moves with the supply temperature. They are not judged.

**Results (24 October 2026, about 04:07 CDT; derived by computation, one thread, 3 s; `low_t/para_supply_bound_2026-10-24.txt`).**

- **The check passed:** x_p = 1 at 199.67 K gives 162.7 kJ/kg, against (b)'s 163.
- **The supply's enthalpy difference** h_eq − h_frozen, at x_p 0.998 (at 0.95):

  | T_s | Δh, kJ/kg of H2 |
  |---|---|
  | 150 K | 306.4 (285.8) |
  | 200 K | 161.5 (151.0) |
  | 250 K | 70.2 (65.7) |
  | 300 K | 27.5 (25.8) |

- **Judged at T_s = 200 K, 36 cases: NOT below 2e-4.**
  - f is 1.7e-3 to 2.8e-3.
  - Δc\*/c\* is −4.6e-4 to −1.0e-3.
  - ΔIsp/Isp is −5.8e-4 to −1.22e-3.
  - The largest of the three is 2.81e-3 (O/F 4, x_p 0.998, LOX). That is 14 times the bar on f, 5 times on c\* and 6 times on Isp.
  - The changes are negative: the para-rich supply carries less energy than an equilibrium supply at the same temperature.
  - All three are linear in Δh: the ratios between x_p 0.95 and 0.998 match Δh's ratio of 0.935.
  - O/F moves the result most. On c\*, p_c from 2 to 4 MPa moves it by at most 4.8%, and the O2 state by at most 4.5%.
- **Beside it, not judged:**
  - 150 K: largest 5.40e-3.
  - 250 K: 1.21e-3.
  - 300 K: 4.68e-4. Even at 300 K, f is above the bar; Δc\*/c\* there is at most 1.64e-4, and ΔIsp/Isp at most 2.03e-4.
- **So, by the Director's rule, the form needs a sourced conversion rate before it is fixed.** No conversion model is built.

**What a conversion rate can still move (derived, not judged).**
- **Para and equilibrium hydrogen have the same enthalpy once hot.** h_eq − h_frozen(x_p 0.998) is 27.5 kJ/kg at 300 K, 3.67 at 400 K, 0.458 at 500 K, 0.0079 at 700 K and 3.2e-5 at 1000 K.
  - So conversion above 500 K moves at most Y_H2 Δh(500 K) / q = 8.1e-6 of the heat release, over the whole grid; above 1000 K, 5.7e-10.
  - Book the supply's frozen enthalpy at the face. Then the chamber's conversion changes c\* and Isp in this bookkeeping by at most 8.1e-6, 25 times under the bar.
    - That holds wherever the conversion happens, and even if it has not happened by the time the gas passes 500 K.
    - Energy is conserved, and hot unconverted gas carries equilibrium's enthalpy.
- **The c\* and Isp effect is set before the injector.** It depends on how much of the supply has converted by the time it reaches the face:
  - none: frozen, the full effect above;
  - all: equilibrium, no effect.
  - So the rate the rule needs is the one in the tank, the jacket and the lines (self-conversion in the gas, and catalysis on the jacket's metal), not the hot zone's H-atom exchange.
- **The x_p choice inside the spec's range sits under the bar.** At 200 K, x_p 0.95 against 0.998 moves c\* by at most 6.5e-5 and Isp by 8.0e-5.
- **The cold stream's local state is what in-chamber conversion does move.**
  - Option (ii) gives the face inlet h_frozen(T_s) on equilibrium thermo. That is equilibration at the inlet.
  - The inlet gas then sits at 188.1 K instead of 200 K: 11.9 K colder, and 6.3% denser at the same pressure. (At 150 K: 23.3 K colder and 18.4% denser. At 250 K: 5.0 K and 2.0%.)
  - At a fixed mass flow, the H2's injection velocity falls by 6.3%. That changes the momentum-flux ratio, which is a mixing input, not part of the energy budget.
  - A para-rich H2 thermo below some T\* ≥ 500 K (option (i)) would keep the inlet at T_s, at the cost of Table B and the mechanism's thermo, as (b) says.

**The citations (a Sonnet literature agent, 24 October, 04:04 to 04:07 CDT; what each source is, and what was read).**
- **The LH2 para fraction.**
  - Equilibrium at 20.27 K is 99.800% para, and at 20.00 K 99.822% (derived here from the partition functions).
  - The standard compilation is Woolley, Scott and Brickwedde, J. Res. NBS 41, 379 (1948), RP1932. It was located, but its tables are scanned images and were not read. Secondary sources quote 99.79% at 20.4 K.
  - The delivered minimum is 95.0% (MIL-PRF-27201, read directly, (b)).
- **The conversion time.**
  - Milenko, Sibileva and Strzhemechny, "Natural ortho-para conversion rate in liquid and gaseous hydrogen", J. Low Temp. Phys. 107, 77–92 (1997), DOI 10.1007/BF02396837. It measured the liquid at 17 to 32 K and the gas at 40 to 120 K, to 0.09 g/cm3, with an interpolation formula over the fluid range. The scope was read in the abstract; **the rate itself was not read** (paywalled).
  - Larsen, Simon and Swenson, Rev. Sci. Instrum. 19, 266 (1948), on the liquid: located, but the number was not read.
  - Farkas, "Orthohydrogen, Parahydrogen and Heavy Hydrogen" (Cambridge, 1935): not accessed.
  - **So the "hours" in the decision is consistent with these sources' scope, but no conversion time has been read from a primary source.** The gap for the jacket's metal stands ((b)).
- **The hot zone.**
  - Schulz and Le Roy, J. Chem. Phys. 42, 3869 (1965), measured H + p-H2 conversion at 300 to 444 K. That range was read in Lique, Honvault and Faure, arXiv:1402.5292. No Arrhenius form was verified.
  - Ridley, Schulz and Le Roy, J. Chem. Phys. 44, 3344 (1966): located, not read.
  - Paramagnetic O2 catalyses conversion (Nielsen and Dahler, J. Chem. Phys. 46, 732, 1967, read). No rate constant was extracted.
  - By the derived point above, the hot zone's rate does not reach the bar in any case.

**For the Director.**
- By the rule as registered, the supply's form is not fixed until a conversion rate is sourced.
- The derived point above narrows which rate that is: the rate between the tank and the face (gas self-conversion over the line's residence time, and the jacket wall). Inside the chamber, the conversion's place moves c\* and Isp by under 8.1e-6.
- If the line rate confirms the frozen supply, option (ii) at x_p 0.998 gets the energy right. The cost is the cold inlet's 6.3% density offset at 200 K.
