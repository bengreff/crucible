# RL10A-3-3A: prediction protocol (pre-registered)

Status: protocol written 13 October 2026, before any prediction run and before any measured performance of this engine was read (TECHNICAL_PLAN step 11). It fixes what is predicted, from which inputs, by which method, with which brackets, and how the comparison is read. The predicted numbers themselves are committed in a later commit, before the comparison reads any measured value; amendments to this protocol before that commit are allowed and are dated here.

Inputs: `RL10A33A_INPUTS.md` (screened 13 October: no outcome value), `rl10a33a/tm107318_table_E1*.csv` (Table E1, area against station) and `rl10a33a/tm107318_inputs_addendum.md` (requested 13 October: plenum temperatures, the injector-face and exit stations from Figure E1, wall geometry). Source: Binder, Tomsik and Veres, NASA TM-107318 (1997).

## 1. Declarations

- **Exposure on 8 October 2026.** A keyword grep of the archived pre-pivot file `archive/pre-pivot-2026-09-17/data/anchors/rl10_contour.csv` printed its comment lines. They contained a c* computed by the pre-pivot project from a chamber pressure, the throat area and the mass flow, compared with a "record" c* value, and a throat dimension of 2.47 in. The file was not read further and was never used. c* is not one of the predicted quantities below, but it is close to them: the method below has no free parameter that could absorb the exposure, and every bracket is set by a stated rule, not by a value.
- **Prior knowledge.** The author of the prediction (Claude, through the Director's workers) carries general knowledge of the RL10 family's published ratings from its training. No such value is an input or a check anywhere in the method.
- The throat dimension: TM-107318 Table 2.5.1's 2.47 in was transcribed as the throat diameter (and 5.13 in as the chamber diameter).
- **Open question on absolute size (raised 13 October 2026, 22:05 CDT, before any prediction).** It is unresolved whether these two numbers are diameters or radii. The radii in section 2 scale with the answer by a factor of 2; area ratios and stations do not change.
  - Trigger 1, internal to the inputs: the face-plate radius derived from the source's "length of plate" (3.894 in) exceeds the 2.565 in chamber radius that the diameter reading implies.
  - Trigger 2, prior knowledge (declared): under the diameter reading, the throat mass flux implies a chamber pressure several times what the author expects of an expander-cycle chamber. That expectation comes from general knowledge, so it is not evidence and will not be used to decide.
  - Resolution: the source's own geometry only (row labels verbatim, Figure E1's scale, coolant-tube dimensions, Table 6.1.1's geometry rows), requested from the extraction agent with no performance value allowed. The answer and its quotes go into the addendum, and section 2 is amended with a dated note before the prediction commit. If the source cannot settle it, both readings become a bracket.

## 2. The case

- **Operating point P1 (primary):** TM-107318 Table 6.1.1's operating point. Oxidizer 13.9525 kg/s (30.76 lb/s). Combustion O/F 5.055, so the fuel burnt in the chamber is 13.9525 / 5.055 = 2.7601 kg/s and the chamber flow is 16.7126 kg/s. The engine-inlet fuel flow is 2.7946 kg/s (6.161 lb/s, engine-inlet O/F 4.993). The 0.0345 kg/s difference is assumed to leave the engine without thrust (its path is not in the inputs; declared).
- **Operating point P2 (secondary, fast tier only unless the full tier is cheap):** Table 2.5.1's design point, O/F 5.26 and total flow 16.9462 kg/s (37.36 lbm/s), all assumed burnt.
- **Geometry.** Throat radius r_t = 1.235 in (0.031369 m). The wall is r(z) = r_t sqrt(A/A_t) from Table E1's 16 nodes, with the throat point (0, r_t) added, interpolated by a monotone cubic (PCHIP) in r(z). Upstream of the first node (-12.00 in, r = 2.552 in, against the chamber radius of 2.565 in from Table 2.5.1) the wall joins the chamber radius at the injector face: one more node there at the chamber radius, with PCHIP through it (a taper of 0.013 in over the last stretch). The exit node is (exit station, r_t sqrt(61)). These radii are under the diameter reading of section 1's open question. The exit is where A/A_t = 61. Missing stations are brackets (section 4).
- **Propellants.** The reported case (tier 1) injects LOX as liquid (step 5, `docs/evidence/SPRAY.md`) and hydrogen as gas, each at its Table 6.1.1 injector-plenum temperature (addendum). Development runs with pre-vaporized propellants are labelled development and are never the reported case (section 6).
- **Ambient.** The predictions are vacuum quantities. If the measured thrust is printed at a test-cell pressure, the comparison converts with F_amb = F_vac - p_amb A_e, A_e = 61 A_t, using the source's printed ambient pressure (a test condition, read at comparison time).

## 3. Predicted quantities and tiers

Each tier predicts, at P1 (and P2 where stated): vacuum thrust F_vac; vacuum Isp on the chamber flow and on the engine-inlet flow; the nozzle stagnation pressure and the injector-face pressure; c* = p0 A_t / mdot; and, for tier 1, the heat flow through the hot-gas wall (the expander cycle's driver, `TECHNICAL_PLAN` step 17).

- **Tier 0, ideal rocket with NASA CEA** (RocketCEA, `tools/cea_point.py`; minutes). Reactant cards: O2(L) and H2 gas at their plenum enthalpies. The chamber pressure is not an input: p_c solves mdot = Cd p_c A_t / c*_CEA(p_c) by fixed point. Cd from Kliegel and Levine at the throat bracket below. Isp_vac from CEA at eps 61, shifting and frozen, times the divergence factor lambda. lambda is bracketed between (1 + cos theta)/2 at the exit wall angle (theta from Table E1's last nodes, 5.1 to 6.5 degrees) and at the equivalent cone half-angle from throat to exit (about 10.5 degrees). Tier 0 has no boundary-layer, kinetic, mixing or atomisation loss; it is the ideal band, not the reported prediction.
- **Tier 1, the engine (reported).** The 2-D axisymmetric engine as it stands at the prediction commit: liquid LOX parcels (step 5), equilibrium or finite-rate tables (Table A, or Table B if verified by then), SST with the Hasan-Pecnik corrections and wall functions, adaptive resolution, started "valves open, ignite" (step 1) and run until settled. F_vac comes from the device-thrust ledger plus the ambient force; Isp = F_vac / (mdot g0). Settled means: F_vac, the injector-face pressure and the outlet mass flow each drift by less than 0.1% over the last 1 ms. A run that cannot afford the full start says so with its measured cost, and its result is labelled development (step 1).

## 4. Brackets for missing inputs

Each bracket is varied alone about the central case. The judged band of a quantity is the central value plus and minus the sum of the absolute one-at-a-time deviations on each side (an envelope of intervals, not a standard deviation); the quadrature sum is printed beside it, not judged.

| Missing input | Central | Bracket | Rule |
|---|---|---|---|
| Throat shape (no printed radius; Appendix E.2 describes a sharp-edged cast silver insert) | PCHIP through the nodes and (0, r_t) | straight segments in r(z) from (-0.500, 1.3012 in) to (0, 1.235 in) to (0.485, 1.4027 in), a sharp vertex | the two shapes that pass through the same printed nodes. Tier 0's Cd: Kliegel-Levine at the radius of curvature of the central shape at the throat (the mean of its two one-sided values), and at R/r_t = 0.5 for the sharp shape (the series' edge) |
| Injector face station | Figure E1's estimate (addendum) | its stated uncertainty, or -12.50 to -14.50 in if Figure E1 gives none | the chamber is the measured part: Table E1 starts at -12.00 in and Figure E1 shows the face upstream of it |
| Exit station (A/A_t = 61) | Figure E1's estimate if given, else 45.64 in | 45.29 in (linear in r from the last two nodes) to 45.98 in (parabola in r through the last three) | the two extrapolations of the printed nodes |
| Hot-gas wall temperature (uniform) | 567 K | 300 to 834 K | the upper end is Table 6.1.1's maximum metal temperature (a cycle-model output used as an input bound); the lower is a cold wall |
| LOX drop size (step 5) | `SPRAY.md`'s central SMD | `SPRAY.md`'s bracket (Burick's curve and its scatter) | fixed by step 5's design, before this case |
| Fuel and oxidizer injection temperatures | Table 6.1.1's plenum values (addendum) | plus and minus the addendum's transcription uncertainty | inputs, not tuned |
| Supply turbulence | intensity 5%, length scale one element pitch | intensity 2% to 10% | no injector turbulence data in the source |
| Face distribution | uniform O/F and mass flux over the face | none (declared) | the source gives no radial distribution; 162 or 216 elements (unreconciled in the source) do not change a uniform face |
| Fuel flow not burnt (0.0345 kg/s) | leaves without thrust | none (both Isp definitions are printed) | its path is not in the inputs |

Grid: tier 1's settled quantities on two resolutions (the finest AMR level halved); the Richardson estimate of the remaining discretisation error joins the band as one more deviation. The wall functions' measured errors (`docs/evidence/WALL_FUNCTIONS.md`, criteria 2 and 3) join the wall heat flow's band.

## 5. The comparison

- **Where the measured values are** (pages only; none read): TM-107318 Table 2.5.1 (printed p.6), Table 6.1.1 (p.37), Table E3 (p.140), Figures E7, E8, E10 and E11 (pp.145 to 147), Figure 6.1.1 (p.70). Transients (Figures 6.2.x and 6.3.x, Table 6.2.1) are outside this protocol.
- **Order.** (1) This protocol committed. (2) The predicted numbers and their bands committed, with the commit hash recorded. (3) Only then is any performance row read, by me or an agent, into `RL10A33A_COMPARISON.md`, which quotes the prediction's commit hash.
- **What counts as measured.** Each value is taken with the report's own label. A value the report labels as test data is a validation comparison. A value labelled as the authors' cycle-model output, rated or design is recorded as a model or rating comparison and is not called validation.
- **Reading.** For each quantity: the central error, and whether the measured value lies inside the judged band. No pass threshold is set. A miss outside the band is a finding, explained through the loss ledger (divergence, kinetics, boundary layer and wall heat, mixing, atomisation and evaporation) or reported as unexplained. Nothing is retuned after the reading; a model change made afterwards is a new prediction against the same data and is labelled as not blind.

## 6. Development runs

Before the reported run, development runs of this chamber with pre-vaporized propellants (gaseous O2 and H2 premixed at the face) measure the RL10-sized cost and apply the implicit-time decision rule (`TECHNICAL_PLAN`, *Lightweight engine*). They use only the inputs above. Their Isp and thrust are printed in their logs as development quantities, beside NASA CEA, and are never the reported case. They read no measured performance.
