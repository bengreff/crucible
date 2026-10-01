# Validation-data survey for the first chemical case (2026-09-30)

Purpose: choose real hardware whose geometry, feed conditions and measured performance are
published well enough that CRUCIBLE can be run blind against it (VISION_SCOPE: chemical with
real-data validation is the first substantive milestone; Ben: "well documented first").

Method: a Sonnet research subagent surveyed the open literature. I then spot-checked the
sources behind both recommendations myself (NTRS abstracts and the TRR40 2015 summer-program
report, full text). Labels: **verified** means I read it in a primary or full-text source
tonight. **Secondary** means it is quoted by another paper. **Unverified** means it came from
search snippets only.

## Table

| Rig / engine | Propellants | Geometry published | Feed / operating points | Measured (uncertainty) | Access | Missing inputs | Sources |
|---|---|---|---|---|---|---|---|
| **TUM single-element, square chamber (SFB/TRR40 Test Case 1)** | GOX/GCH4, gas/gas | Yes. Square cross-section, one coaxial element flush with the face plate, rectangular throat, contraction ratio 2.5 (**verified**). The 12x12 mm section size is **unverified** | Pc 20 bar, O/F 2.6 nominal (**verified**). Mass flows set by sonic orifices (**verified**). The 5-20 bar and O/F 2.6/3.0/3.4 sweep is **unverified** | Axial wall-pressure profile, wall heat flux reconstructed from thermocouples in the copper wall, combustion efficiency (**verified**). Uncertainty is not stated in the report | Open. TUM report PDF; EUCASS proceedings | Exact injector dimensions (post, recess, annulus), inflow turbulence, the measured values themselves (held back from the blind test), all stated in the Test Case 1 description (Haidn et al. 2015, 6th TRR40 technical report), which is not yet in hand | Roth, Haidn et al., SFB/TRR40 Summer Program Report 2015, pp. 101ff; Celano et al., EUCASS 2015 (round vs square chambers); Celano et al., EUCASS 2016 (injector characterization) |
| TUM single-element, round chamber | GOX/GCH4 (also GH2), gas/gas | Yes, 2-D axisymmetric case (**secondary**) | O/F 2.2: CH4 0.0153 kg/s at 268 K, O2 0.0339 kg/s at 276 K (**secondary**) | Pressure about 4%, heat flux about 10% (**secondary**, from Aerospace 2021 8(11):341, citing the TUM authors) | Open (MDPI, DLR elib) | Same as above | Celano et al. EUCASS 2015; Aerospace 2021, 8, 341 |
| **Penn State single-element shear-coaxial** | GO2/GH2, gas/gas | Yes, 1.5 in (38.1 mm) round chamber (**verified**, NTRS abstract) | Pc 750/600/450/300 psia. O/F 6.6 (preburner, vitiated) and 6.0 (ambient). About 0.27-0.29 lbm/s at 750 psia, scaled linearly (**verified**) | Axial wall heat flux from Gardon gauges and coaxial thermocouples (**verified**). Uncertainty is **unverified** (needs full text) | Open. NTRS 20050209932 | Injector recess and lip dimensions, preburner product composition for the 6.6 case, inflow profiles | Marshall, Pal, Woodward, Santoro, AIAA 2005-3572; Tucker, Menon, Merkle, Oefelein, Yang, AIAA 2008-5226 (NTRS 20080048109, multi-code RANS/LES comparison, **verified**) |
| TUM 7-element | GOX/GCH4, gas/gas | Yes, 30 mm bore, 7 coaxial elements (**unverified**) | Per paper, not consolidated | Wall heat flux and temperature | Open (TRR40 reports, EUCASS) | Element-to-element flow split | Perakis et al. TRR40 SS2017; Zhang 2017/2019 |
| DLR/ONERA Mascotte | LOX/GH2, LOX/GCH4 (liquid oxidizer) | Partial, varies by campaign | Per campaign | OH* imaging, flame structure | Open (DLR elib, EUCASS) | Spray atomization (droplet sizes): an extra unresolved physics layer | Yang et al., DLR elib |
| DLR P8 BKD/BKH | LOX/H2, LOX/CH4 (liquid) | Partial | Pc up to about 65 bar | Heat flux by inverse method plus calorimetry | Open papers, partial matrices | Atomization; cooling-channel details | AIAA JPP papers (**unverified**) |
| Purdue CVRC | GO2/GCH4 (varies by campaign) | Yes | Instability operating points | High-frequency pressure, OH* | Open (AIAA 2009-234) | It is an instability rig, not a steady performance benchmark | Yu et al. AIAA 2009-234 |
| RL10A-3-3A | LOX/LH2, full engine | No chamber or injector geometry | Pc about 475 psia, O/F 5 | Thrust, Isp, Pc: global only | Open (NASA TM-107318, cycle model) | Chamber and injector geometry, any local data | NASA TM-107318 |
| SSME/RS-25 | LOX/LH2, full engine | No | Engine level | Engine level | Partial | Geometry and local data | NASA-TM-103274 |
| JAXA LE-5B | LOX/LH2, full engine | No | Pc about 3.6 MPa, O/F 5 | Isp about 447 s, global | Open summaries | Geometry and local data | MHI Tech. Review 53(4) 2016; EUCASS 2019-0626 |
| Raptor / Merlin | CH4/LOX, RP-1/LOX | No | Not published | Not published | Not open | Everything | None credible: unusable |

## Decision (Ben, 2026-09-30)

- **First case:** the TUM **round** chamber (axisymmetric, published geometry used exactly). Inputs and gaps are
  in `docs/VALIDATION_TUM_ROUND.md`.
- **Second case:** a flight engine, RL10A-3-3A (NASA TM-107318).
- **Later:** the square TRR40 case, labelled as a geometry approximation.
- **Optional:** Penn State, only as a cheap heat-flux side check.
- **Blind protocol:** the prediction is pre-registered before digitising.

## Original recommendation (superseded by the decision above)

1. **TUM single-element GOX/GCH4, SFB/TRR40 Test Case 1 (20 bar, O/F 2.6).**
   - It was built as a blind CFD validation case: seven groups were scored against data they had not seen.
   - Its validation quantities include an axial wall-pressure profile and a combustion efficiency. Those are performance quantities, which match the performance-prediction framing.
   - It is gas/gas and single-element, with a nozzle at a flight-like contraction ratio of 2.5.
   - It is open.
   - Next input to obtain: the Test Case 1 description (geometry and boundary conditions), plus the experimental curves in the 2015 report figures. Those curves must be digitised, and the digitising uncertainty must be stated.
2. **Penn State GO2/GH2 single element (Marshall et al. 2005).**
   - It is the best-known gas/gas benchmark.
   - It covers four pressure levels, and NASA's multi-code comparison (Tucker et al. 2008) shows how RANS, URANS and LES actually fared against it.
   - Caveat: its measurement is wall heat flux only. There is no chamber-performance quantity, so it tests the near-wall solution more than the thrust budget.

Hydrogen and methane are both covered by the item-1 thermochemistry verification (H2/O2 and
CH4/O2 at these rigs' pressures and O/F). LOX/LH2 is verified too, for a later RL10-class
global check, but no LOX/LH2 engine has open chamber geometry, so none qualifies as a first case.

## What a CEA match does not establish

The ideal-rocket agreement in `docs/THERMO_VERIFICATION.md` only verifies the thermochemistry. It
does not validate the chamber. The TUM data exists precisely because real chambers fall short
of ideal: incomplete mixing, heat loss to the wall and finite-rate chemistry.
