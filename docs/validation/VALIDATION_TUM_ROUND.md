> **Status, 4 October 2026:** Ben chose RL10A-3-3A as the next engine case (VISION_SCOPE, *Changes since the pivot*). This TUM chamber is now a component case for mixing and wall heat, still blocked on the missing inputs listed below. The blind protocol below applies to every validation case.

# First validation case: TUM single-element GOX/GCH4, round chamber (inputs inventory, 2026-09-30)

Decision (Ben, 2026-09-30):
- **Primary case:** the TUM round chamber, with the published geometry used exactly. If an input is missing,
  name it and stop; do not invent an equivalent.
- **Blind protocol:** the prediction is pre-registered. Predicted wall-pressure profile, wall heat flux and
  combustion efficiency are committed before any measured curve is digitised. Test-derived inputs, such as the
  wall temperature, are declared as inputs.
- **Scope:** this is the first case, not the milestone. The second case is a flight engine (RL10A-3-3A, below).

Sources (a Sonnet agent read these in full text; I spot-checked the TRR40 report myself):

| Key | Source |
|---|---|
| [R/S] | Silvestri, Celano, Schlieben, Knab, Haidn, "Comparison of Single Element Rocket Combustion Chambers with Round and Square Cross Section", 6th EUCASS 2015 |
| [Roth17] | Roth, Silvestri, Perakis, Haidn, ISTS 2017-o-4-01 |
| [Zhukov19] | Zhukov, van Schyndel, Gomez, ISTS 2019-a-46 |
| [Aero21] | Strokach et al., Aerospace 2021, 8(11), 341 |

## Decisions on the gaps (Ben, 2026-09-30)

- **Convergent contour:** run a declared bracket of 30 and 45 deg convergent half-angle, each with a sharp
  throat and a rounded throat. The spread is carried as a declared input uncertainty in the pre-registered
  prediction.
- **Exact shape:** the geometry code takes the contour as data, so the exact shape drops in when the TRR40
  Test Case 1 document or a drawing arrives. Obtaining it is on Ben's list.
- **Blind points:** pre-register at O/F 2.6, 3.0 and 3.4. O/F 2.2 is reported as not blind for chamber
  pressure and efficiency.
- **Research agents:** they are told explicitly never to report measured results for the validation case.
- **Order:** TUM stays first. It is gaseous GOX/GCH4, which the first chemical implementation supports.
  RL10 injects liquid oxygen, so it needs breakup and evaporation models first.

## Inputs found (published, with source)

| Input | Value | Source |
|---|---|---|
| Chamber inner diameter | 12 mm | [R/S] Table 1; [Roth17] 2.1 |
| Cylindrical length | 285 mm | [R/S] Table 1 |
| Total length | 305 mm | [Roth17] 2.1 |
| Throat diameter | 7.6 mm (contraction 2.5) | [R/S] Table 1 |
| Nozzle | "conical" | [Roth17] 2.1 |
| GOX post | inner diameter 4 mm, wall 0.5 mm, length 96 mm, flush (recess 0), tip not tapered | [R/S] Table 2; [Roth17] Table 1 |
| Fuel annulus | outer diameter 6 mm. Area ratio fuel/ox 0.69 from the dimensions; a 2013 paper's "1.26" is inconsistent with its own dimensions | [R/S] Table 2 |
| Operating points (CH4) | O/F 2.2, 2.6, 3.0, 3.4 at nominal 20 bar | [Roth17] Table 2 |
| Inlet temperatures (fuel / ox) | 268/276, 269/275, 270/274, 271/273 K | [Roth17] Table 2 |
| Mass flows at O/F 2.2 | CH4 15.3 g/s, O2 33.9 g/s | [Aero21] 2.1; [Zhukov19] Table 1 |
| Mass flows at the other points | Fig. 3 only, to be digitised as **inputs** | [Roth17] Fig. 3 |
| Wall temperature profile | figure only, to be digitised as an **input** | [Aero21] Fig. 2 |
| Pressure transducers | 9, 34 mm pitch; absolute positions only in figures | [R/S] 3.2 |
| Thermocouples | 1 mm row at 17 mm pitch; inverse heat-conduction method | [R/S], [Roth17] |
| Uncertainty | pressure 4% (dominated by mass flow); heat flux 10-15% | [Aero21]; [Zhukov19] |
| Burn | 3 s evaluation window (facility practice; not restated for this dataset) | [R/S] 3.3 |

## Missing inputs (exact)

1. **Convergent nozzle contour.** The convergent half-angle and the throat radius of curvature are not
   published in any source found.
   - Not needed for the chamber-wall observables: the divergent section is downstream of the sonic throat and
     cannot influence them.
   - Needed for the convergent section: it sets the throat discharge coefficient, and therefore chamber pressure
     for a given mass flow (inferred: of order 1-2% for conical inlets), and the flow near the last transducer.
2. **The TRR40 "Test Case 1" master description** (Celano et al., Space Propulsion 2014 / AIAA 2016-2148). It is
   not open access. It is likely to contain the dimensioned drawing and the operating table.
3. **Thermal condition on the nozzle wall.** The measured wall-temperature profile covers the chamber sidewall
   only.
4. **Propellant purity.** Not stated.

## Blind-protocol exposure (declared)

While surveying sources, a research agent reported two measured results despite instructions not to:
- the measured combustion efficiencies at O/F 2.2 ([R/S] Table 5);
- the measured chamber pressure at O/F 2.2 (approximately 19 bar).

I have seen those two items. Nothing else has been seen: no pressure profile, heat-flux curve, or efficiency at
other O/F. O/F 2.2 is therefore **not blind** for chamber pressure and efficiency.

## Second case: RL10A-3-3A (flight engine)

| Source | Key content | Status |
|---|---|---|
| NASA TM-107318 (Binder, Tomsik, Veres 1997; NTRS 19970010379) | Table 2.5.1: chamber, throat, area ratio 61, Pc, O/F, c* efficiency, Isp | |
| | Appendix E Table E1: area ratio against axial station from -12 in to +41.84 in | For an axisymmetric chamber this gives the wall radius exactly, r = r_t sqrt(A/A*), so no shape assumption is needed |
| | Section 6.1: system model against measured ground and flight tests (engines P2087/P2093) | |
| | Injector: 216 coaxial elements in 8 rings, 162 with LOX swirlers | An axisymmetric model can only represent this as annular averages, a declared limitation (TECHNICAL_PLAN) |
| PWA FR-1769 (1966, RL10A-3-3) | Injector and feed detail | Its performance table is "estimated", not measured |
| Both | Instrumentation-grade Isp uncertainty | Missing. TM-107318 itself says the provenance of the c* efficiency curves is unclear |

Agent-verified, not yet read by me.
