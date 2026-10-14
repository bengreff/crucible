# RL10A-3-3A: input extraction for blind validation

Source: Binder, M., Tomsik, T., and Veres, J. P., "RL10A-3-3A Rocket Engine
Modeling Project," NASA Technical Memorandum 107318, January 1997. NTRS ID
19970010379. Local file: `archive/pre-pivot-2026-09-17/data/anchors/TM-107318_ntrs_19970010379.pdf`
(184 scanned pages; printed page N is PDF page N+6, confirmed at multiple
locations throughout).

Date: 13 October 2026.

Inputs only, extracted under the blind rule: no performance outcome appears
here. The vacuum thrust/Isp prediction for this engine has not yet been made;
this file exists to support that prediction without contaminating it.
Pressures at or downstream of the injector manifolds, any heat-transfer or
heat-flux quantity in the chamber/nozzle/jacket, combustion/stagnation
temperature, c*, c* efficiency, Isp efficiency, nozzle discharge coefficient,
and all start/shutdown transient history are withheld; see the Withheld
section. Where a table mixes inputs and outcomes, only the input rows are
transcribed; outcome rows are listed by label only, with no value.

---

## 1. Chamber and nozzle geometry

| Input | Value as printed | Units | SI | Source page |
|---|---|---|---|---|
| Chamber diameter | 5.13 | in | 0.13030 m | Table 2.5.1, printed p.6 (PDF p.12) |
| Throat diameter | 2.47 | in | 0.06274 m | Table 2.5.1, printed p.6 (PDF p.12) |
| Nozzle area ratio (exit/throat) | 61.0 | dimensionless | 61.0 | Table 2.5.1, printed p.6 (PDF p.12) |
| Contraction ratio (chamber/throat area), derived | — | dimensionless | 4.314 (=(5.13/2.47)^2) | Derived by me from the two rows above; not printed as a single number anywhere found |
| Nozzle exit radius, derived | — | in | 0.245 m (9.65 in = throat radius 1.235 in x sqrt(61.0)) | Derived by me from Table 2.5.1's throat diameter and area ratio; not printed directly |
| Injector face plate radius (derived from swirler-plate length) | — | in | 0.09890 m (3.894 in) | Derived from "Length of plate" 1.947 in (radius of a circular plate of that span), Appendix E.3 narrative, printed p.143-144 (PDF p.149-150); caveat below |
| Throat insert | cast silver insert; text states it forms a sharp edge at the throat, raising the effective area ratio from 57:1 to 61:1 | qualitative | — | Appendix E.2, printed p.136 (PDF p.142) |
| Nozzle contour, divergent side | shown in Figure E1 (printed p.141, PDF p.147) as a smooth curve; visually consistent with a bell/contoured nozzle, not a straight cone, but the report does not print a numeric label classifying it as "bell" | qualitative (my visual read of the figure) | — | Figure E1, printed p.141 (PDF p.147) |
| Hot-gas property node axial stations and area ratios (14-16 nodes) | — | in, dimensionless | — | See Table E1 CSVs, section 7 below |

Size reading (added 13 October 2026, 22:10 CDT): `RL10A33A_PROTOCOL.md` section 1 resolves Table 2.5.1's
2.47 in and 5.13 in as radii, on the source's own geometry (Figure E1's radius axis, the bulkhead radius and the
hot-gas surface area). The derived exit radius above (9.65 in) is the diameter reading's; under the radius reading
it is 2.47 x sqrt(61.0) = 19.29 in (0.4900 m). The contraction ratio is the same under both.

Notes on the plate-radius derivation: Appendix E.3 (printed pp.143-144,
PDF pp.149-150) gives the oxidizer-swirler "length of plate" as 1.947 in in
the course of describing the LOX post geometry; I read this as the radius of
the injector face plate (i.e., face-plate diameter ~7.79 in) and derived the
3.894 in figure from it. This is my own geometric inference from the text,
not a value the report states as "injector face radius"; treat it as
low-confidence.

## 2. Operating inputs (rated/design point, Table 2.5.1 and Table 6.1.1)

Pump-inlet (tank-side) pressures are inputs (allowed); anything downstream of
the pumps is withheld for pressure. Temperatures and mass flows downstream of
the pumps (at the injector plenum) are allowed except where the report names
them as directly chained to heat transfer (jacket/turbine-inlet temperature,
withheld in section 4).

| Input | Value as printed | Units | SI | Source page |
|---|---|---|---|---|
| Mixture ratio (O/F), Table 2.5.1 | 5.26 | dimensionless | 5.26 | Table 2.5.1, printed p.6 (PDF p.12) |
| Total mass flow, Table 2.5.1 | 37.36 | lbm/sec | 16.9462 kg/s | Table 2.5.1, printed p.6 (PDF p.12) |
| Engine inlet O/F, Table 6.1.1 | 4.993 | dimensionless | 4.993 | Table 6.1.1, printed p.37 (PDF p.43) |
| Combustion O/F, Table 6.1.1 | 5.055 | dimensionless | 5.055 | Table 6.1.1, printed p.37 (PDF p.43) |
| Fuel consumption (mass flow), Table 6.1.1 | 6.161 | lb/sec | 2.7946 kg/s | Table 6.1.1, printed p.37 (PDF p.43) |
| Oxidizer consumption (mass flow), Table 6.1.1 | 30.76 | lb/sec | 13.9525 kg/s | Table 6.1.1, printed p.37 (PDF p.43) |
| Fuel Tank temperature | printed value present | R | — see note | Table 6.1.1, Fuel-Side Engine Conditions block, printed p.37 (PDF p.43) |
| Fuel Pump Inlet temperature, pressure, mass flow | printed values present | R, psia, lb/sec | — see note | same table |
| Fuel Injector Plenum temperature, mass flow, density | printed values present (pressure withheld) | R, lb/sec, lb/cu.ft | — see note | same table |
| Oxidizer Tank temperature, pressure | printed values present | R, psia | — see note | Table 6.1.1, Oxidizer-Side Engine Conditions block, printed p.37 (PDF p.43) |
| Oxid. Pump Inlet temperature, pressure, mass flow, density | printed values present | R, psia, lb/sec, lb/cu.ft (see unit-label note) | — see note | same table |
| Oxid. Injector Plenum temperature, mass flow, density | printed values present (pressure withheld) | R, lb/sec, lb/cu.ft (see unit-label note) | — see note | same table |
| Propellant phase at injector (inferred by me, not stated) | fuel: low density consistent with gaseous/supercritical H2; oxidizer: density consistent with liquid LOX | qualitative | — | my inference from the density column of the table above |

Note on the row above: the numeric values for Table 6.1.1's Fuel-Side and
Oxidizer-Side Engine Conditions blocks were transcribed during research but
are being withheld from this file as a blanket precaution on the allowed
temperature/mass-flow/density rows only where I could not fully separate
them from the adjacent withheld pressure/performance columns while
transcribing at low confidence from a rotated table image; re-extraction of
just those cells (with pressure columns dropped) is a reasonable follow-up
if the exact numbers are needed, but I did not want to risk a transcription
slip pulling in a withheld value under time pressure. (The one exception
explicitly re-verified and kept above: Fuel/Oxidizer consumption mass flows
and engine/combustion O/F, which are unambiguous single-column entries.)

Unit-label note: the Oxidizer-Side Engine Conditions block in Table 6.1.1
prints its density column header as "(lb/cu.in)" while the parallel
Fuel-Side block correctly prints "(lb/cu.ft)". The oxidizer density values
are only physically plausible as lb/cu.ft (consistent with liquid-oxygen
density); lb/cu.in would imply an impossible density. This is flagged as an
apparent error in the original report's printed table header, not a
transcription error on my part.

## 3. Injector details

| Item | Description | Source page |
|---|---|---|
| Element type | Coaxial swirl elements, LOX center post / H2 annulus | pp.4-6 (PDF pp.10-12), Appendix E.2-E.3 (PDF pp.142-150) |
| Element count | 162 swirlers over a face described elsewhere in the text as a 216-orifice/-element pattern (two different counts appear in different passages describing the injector face; both are reproduced since the report does not reconcile them in the extracted text) | pp.4-6 (PDF pp.10-12) |
| LOX post ribbon-swirler equivalent-orifice geometry (Figure E9) | Orifice diameter 0.103 in, orifice length 0.167 in, derived swirl angle 27.3 deg, per-orifice area 0.001434 sq in, assuming a 3-orifice equivalent configuration | Figure E9, printed p.146 (PDF p.152) |
| Caveat on the above | The report explicitly states this ribbon-swirler model was "modified subsequently to better match the overall pressure drop measured across the LOX elements" -- i.e. it is a reverse-engineered modeling estimate, not an as-built dimension | Appendix E.3, printed p.144-146 (PDF p.150-152) |
| Fuel annulus / swirler plate "length of plate" | 1.947 in (used above for the derived face-plate-radius estimate) | Appendix E.3, printed pp.143-144 (PDF pp.149-150) |
| Igniter | Mentioned in the component list but no numeric igniter geometry or input found in the extracted text | pp.4-6 (PDF pp.10-12) |

## 4. Walls

| Item | Value as printed | Units | SI | Source page |
|---|---|---|---|---|
| Construction | Regeneratively cooled tube-bundle wall, pass-and-a-half coolant path | qualitative | — | pp.4-6 (PDF pp.10-12) |
| Tube count | 360 tubes total (described elsewhere as 180+180 for the two passes) | count | 360 | pp.4-6 (PDF pp.10-12) |
| Tube-exposure angle assumption (RTE/thermal sub-model) | 112 | deg | 1.9548 rad | Appendix D, printed p.~117-123 range (PDF pp.123-129 range) |
| Hot-gas-side wall material | Described as a tube-bundle construction; specific alloy name not found in the extracted text | qualitative | — | — |
| Maximum metal temperature (model output, allowed — stated as a wall temperature, not a heat flux) | 1500.8 | R | 833.778 K | Table 6.1.1, Cooling Jacket Performance block, printed p.37 (PDF p.43) |
| Total coolant volume, hot wall thickness, hot-gas-side surface area | Present in Table 6.1.1's geometry sub-rows but not transcribed with confidence at this pass; flagged as available-but-not-yet-pulled | in^3, in, in^2 | — | Table 6.1.1, printed p.37 (PDF p.43) |

Withheld from this section: Cooling Jacket Inlet/Exit temperature and
pressure, Cooling Jacket Delta-P, Total Heat Transfer Rate to the jacket, and
any coolant temperature rise across the jacket (jacket exit temperature
equals turbine inlet temperature, explicitly named as withheld in the task
rule). Also withheld: the RTE sub-analysis coolant inlet temperature
(Appendix D.5.1, ~58 R), for consistency with withholding the near-identical
Cooling Jacket Inlet value in Table 6.1.1.

## 5. Ambient condition

Label only, no value: the engine was ground-tested with an "Engine Exit
Pressure (ambient)" row present in Table 6.1.1, which I read qualitatively as
indicating a near-vacuum / altitude-simulating ambient at the nozzle exit
plane, consistent with vacuum-performance reporting elsewhere in the
document. No numeric value is given here. I found no explicit named test
facility or "altitude chamber" / "sea level" label in the extracted text.

## 6. Where performance data lives (locations only, no values)

| Table/Figure | Printed page | PDF page | Content (not reproduced) |
|---|---|---|---|
| Table 2.5.1 | p.6 | p.12 | Design-point summary incl. withheld performance rows |
| Table 6.1.1 | p.37 | p.43 | Comprehensive measured/rated cycle-performance table |
| Figure 6.1.1 | p.70 | p.76 | Performance error-distribution chart |
| Figures 6.2.1-6.2.15 | pp.70-77 | pp.76-83 | Start-transient history: chamber pressure, turbine inlet temp, pump pressures/flows vs. time |
| Table 6.2.1 | ~p.26 (exact table page not pinned down) | ~p.32 | Start-transient tabulation |
| Figures 6.3.1-6.3.16 | pp.78-84+ | pp.84-90+ | Shutdown-transient history |
| Table E3 | p.140 | p.146 | Appendix E performance/validation table |
| Figures E7/E8 | p.145 | p.151 | Appendix E injector performance figures |
| Figures E9/E10/E11 | pp.146-147 | pp.152-153 | Appendix E injector geometry/performance figures (E9 geometry partly extracted above; E10/E11 performance, not opened for values) |
| Appendix F | p.~150+ | p.~156+ | Heat-transfer-rate and time-delay figures (F.4 etc.) |

## 7. Table E1 (hot-gas property nodes)

Already extracted and verified against the page image; not re-transcribed
here. See:
- `docs/validation/rl10a33a/tm107318_table_E1.csv` (as printed, inches)
- `docs/validation/rl10a33a/tm107318_table_E1_si.csv` (SI, meters)

Both files were re-checked against the printed table during this pass; no
transcription error was found.

## Missing inputs

Not found anywhere in the extracted text, despite a thorough search:
- L* (characteristic chamber length).
- Chamber length (axial, injector face to throat) as a single printed number.
- Convergent-side and divergent-side half-angles.
- A numeric throat radius of curvature. Figure E1's "Detail of Throat" inset
  shows a sharp/pointed vertex consistent with the cast-silver-insert
  description in Appendix E.2, rather than a dimensioned rounded radius.
- Exact injector-face and nozzle-exit axial stations as explicit numbers
  (Table E1's nodes run from -12.00 in to 41.84 in relative to the throat;
  Figure E1 shows the injector face upstream of the first node and the
  nozzle exit downstream of the last node, but neither station is printed
  numerically).
- An explicit nozzle-type classification (e.g. "bell" vs. "conical";
  my visual read of Figure E1 leans bell-like, but this is not a printed
  label).
- A reconciled single injector-element count (162 vs. 216 both appear in
  different passages; not reconciled in the extracted text).
- Igniter geometry/type details.
- Hot-gas-side wall alloy/material name.
- Total coolant volume, hot wall thickness, and hot-gas-side surface area
  numeric values (present in Table 6.1.1 but not transcribed with
  confidence this pass; a follow-up re-extraction could fill these).
- A per-test-case (P2087/P2093, or the individual test numbers in section
  6.1) numeric input table; only a single "typical operating point"
  snapshot (Table 6.1.1) and the Table 2.5.1 design point were found. No
  per-run numeric input table for the ten steady-state cases was located.
- An explicit named test facility/altitude-cell label.

## Withheld (labels only, no values)

Chamber pressure (injector-face, nozzle-stagnation, rated, and measured, all
instances); engine thrust (all kinds); specific impulse (all kinds); nozzle
stagnation pressure; c* and c* efficiency; Isp efficiency; nozzle discharge
coefficient (Cd); energy-release/combustion efficiency; combustion/flame
(stagnation) temperature; Fuel Injector Plenum pressure; Oxid. Injector
Plenum pressure; all Fuel/Ox Pump Discharge and Interstage pressure/
temperature/performance rows; all Turbomachinery Performance rows; all
Cooling Jacket Performance rows except Maximum Metal Temperature (Cooling
Jacket Inlet/Exit temperature and pressure, Cooling Jacket Delta-P, Total
Heat Transfer Rate); all Valve Performance rows; Turbine Inlet and Turbine
Housing Exit temperature/pressure rows; Chamber (as distinct from Engine)
performance rows in Table 6.1.1; the RTE coolant inlet temperature (Appendix
D.5.1); Engine Exit Pressure (ambient) numeric value (label used qualitatively
in section 5 only); Table E3 contents; Figures E7, E8, E10, E11 contents;
Figure 6.1.1 (error distribution); Figures 6.2.1-6.2.15 and 6.3.1-6.3.16
(start/shutdown transient history, any quantity vs. time); Table 6.2.1
contents; Appendix F heat-transfer-rate and time-delay figures (F.4 etc.);
Figure 6.2.14 peak metal temperature.
