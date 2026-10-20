> Screened 20 October 2026, about 18:15 CDT, before being read in full: count-only greps for outcome keywords, then a digit-masked view of every hit. Hits were the agent's own screen notes and label-only withheld items. One igniter-exciter pressurisation range (not an engine outcome) was removed unread. Written by a Sonnet research agent (inputs only); 27 elements per ring is the agent's inference, not sourced.

# RL10A-3-3A injector element geometry — inputs only

Blind-validation input extraction. This file reports geometry, counts, and layout
only. No chamber pressure, flow rate, pressure drop, temperature, heat flux,
thrust, Isp, c*, or efficiency value from any source was recorded.

**Screen run before writing:** searched the draft for the strings `Isp`,
`thrust`, `lbf`, `kN`, `klbf`, `c*`, `cstar`, `efficiency`, `psia`, `psi`,
`MPa`, `bar`, `Btu`, `heat flux`, `sec`, `seconds`. Hits found and their
disposition:
- `psia` appeared twice in source material pulled during research (Table 2.5.1
  combustion-pressure row in TM-107318, and a (value removed) psia exciter-lead
  pressurization figure in the 1966 report's ignition-system section) — both
  dropped, neither appears below.
- `efficiency`, `c*`, `Isp`/specific impulse, and `psia` chamber-pressure rows
  co-occurred in TM-107318 Table 2.5.1 alongside chamber/throat diameter and
  area ratio; only the geometry cells were in scope for this task and even
  those were left out below since the case already has the chamber/nozzle
  contour — the whole table was skipped.
- No remaining hit in the text below is anything other than a geometric
  dimension, a count, a dimensionless ratio/fraction of flow split, an angle,
  or a source page number.

## Sources used

1. Binder, M., Tomsik, T., and Veres, J.P., "RL10A-3-3A Rocket Engine
   Modeling Project," NASA TM-107318, 1997. NTRS 19970010379.
   https://ntrs.nasa.gov/citations/19970010379 (PDF:
   https://ntrs.nasa.gov/api/citations/19970010379/downloads/19970010379.pdf).
   Describes the **RL10A-3-3A**. Cited below as TM-107318.
2. Pratt & Whitney Aircraft, "Design Report for RL10A-3-3 Rocket Engine,"
   PWA FR-1769, Contract NAS8-15494, 28 February 1966.
   https://www.nasa.gov/wp-content/uploads/2025/06/design-report-for-rl-10-a-3-3-1966.pdf.
   Describes the **RL10A-3-3** (predecessor variant, not the "A" suffix).
   Cited below as PWA FR-1769 (1966).
3. Di Matteo, F., De Rosa, M., and Onofri, M., "Transient Simulation of the
   RL-10A-3-3A Rocket Engine," Space Propulsion 2012, paper 2355411.
   https://www.ecosimpro.com/wp-content/uploads/2015/02/SpacePropulsion2012_2355411.pdf.
   Describes the **RL10A-3-3A**. Cited below as Di Matteo et al. (2012).
4. "Rocket Propulsion Evolution: 8.21 – RL10 Engine," Engine History
   Society. https://www.enginehistory.org/Rockets/RPE08.21/RPE08.21.shtml.
   Variant not pinned to a specific suffix; treated as general RL10A-3-3-era
   narrative, cross-checked against (1) and (2). Cited below as EHS.

## 1. Element count and type

- **216 coaxial elements total**, centrally-injected LOX with an annular
  hydrogen sleeve around each (TM-107318 §2.5, p. [body text, "The injector
  has 216 coaxial elements..."], variant A-3-3A; PWA FR-1769 (1966) §I.F,
  p. I-19, "The propellant injector consists of 216 elements arranged in 8
  equally spaced concentric circles," variant A-3-3, drawing/measured — this
  is a design-report figure (FD 1554F) and body text, not a model estimate).
- **162 of the 216 elements have ribbon (not tangential-orifice) flow
  swirlers** on the LOX post; TM-107318 states this for the A-3-3A (§2.5 and
  Appendix E, p. 136: "Except for the inner row, all of the coaxial elements
  have LOX flow swirlers").
- **Reconciling 162 vs. 216** (not done explicitly in TM-107318 itself, done
  here by combining its two statements with the 1966 report's row count):
  TM-107318 Appendix E (p. 136) states the **outer** concentric row injects
  hydrogen only (no LOX, hence no swirler), and separately that **except for
  the inner row**, all coaxial elements have swirlers (i.e., the inner row
  has a LOX orifice but no swirler). PWA FR-1769 (1966) gives 8 equally
  spaced concentric rows for 216 elements. If the 216 elements split evenly
  across 8 rows (27 per row — this even split is **inferred**, not stated
  explicitly in either source), then 2 non-swirled rows (inner + outer) ×
  27 = 54, and 216 − 54 = 162, matching the swirler count exactly. Treat the
  "27 per row" step as a derived/estimated figure, not a sourced one; the
  216/162/8-row figures themselves are sourced as cited.
- **Disagreement between sources on the outer row:** PWA FR-1769 (1966),
  describing the A-3-3, says "All elements **except those in the inner and
  outer rows** incorporate swirlers which aid in the dispersion of the
  oxidizer" (p. I-19) — implying the outer row still has a LOX orifice (just
  no swirler). TM-107318 Appendix E, describing the A-3-3A, says the outer
  row "inject[s] hydrogen only" (p. 136), i.e., no LOX orifice at all in that
  row. This may be a real hardware change between the A-3-3 and A-3-3A, or an
  inconsistency between the two documents; both readings are given here
  since they disagree.
- **Igniter elements:** not part of the 216 main elements; see §5. No
  separate "igniter ring" of main-type elements is described in either
  primary source.
- Di Matteo et al. (2012), describing the A-3-3A, independently corroborates:
  "Most of the injector elements are co-axial, the hydrogen i[s] injected
  through annular orifices around each LOX element. The outer concentric row
  of elements, however, inject hydrogen only" (p. [injector-plate
  discussion section]) — matches TM-107318 Appendix E, not the 1966 report,
  for the outer-row point.

## 2. Layout

- **8 equally spaced concentric circles (rings)** of elements on the face
  (PWA FR-1769 (1966), p. I-19, drawing/measured, variant A-3-3; TM-107318
  does not re-state the ring count in body text but its Appendix E account
  of "inner row" / "outer row" behavior is consistent with a multi-ring
  pattern for the A-3-3A).
- No ring radii, pitch, or element-to-element spacing in absolute units
  (inches or mm) were found in either primary source or in open secondary
  literature searched (see Missing).
- No element-free zones on the face are described beyond the igniter housing
  occupying the geometric center (see §5); no reserved region size is given.
- EHS (variant unspecified) adds no additional quantitative layout detail
  beyond "arranged in equally-spaced concentric circles," consistent with
  (2).

## 3. Element dimensions

- **No LOX post inner diameter, outer/tip diameter, fuel-annulus inner/outer
  diameter, or gap dimension was found in open-source literature** for any
  RL10A-3 variant, despite a Figure E2 (TM-107318, p. 142) cross-sectional
  schematic of "Detail of Typical Injector Element" and a Figure 1-18 (PWA
  FR-1769 (1966), p. I-19) exploded-element schematic — both are qualitative
  line drawings without printed dimensions in the copies examined.
- **Ribbon-swirler geometry (estimate, A-3-3A):** TM-107318 Appendix E,
  Figure E9 (p. 146) gives, for the purpose of deriving an *equivalent*
  tangential-orifice swirler model for the ROCCID combustion-performance
  code, a ribbon diameter **D = 0.103 in (2.62 mm)** and one-turn pitch
  **L = 0.167 in (4.24 mm)**, from which it computes a ribbon helix angle
  θ = tan⁻¹(πD/L) = 62.7° (sweep angle from tangential α = 90° − θ = 27.3°).
  TM-107318 is explicit that this is **not itself the real hardware
  geometry** but an estimate used to build a tangential-inflow-orifice
  surrogate for a code that cannot model ribbon swirlers directly: "it was
  necessary ... to consider the geometry of the ribbon swirlers and estimate
  an equivalent configuration for tangential flow injection" (p. 136). D and
  L are drawn directly from the real ribbon-swirler cross-section in Figure
  E9's own diagram, so they are plausibly real hardware dimensions (ribbon
  diameter and lead/pitch) used as inputs to that derivation, but this is
  not stated outright, and the document does not label them as "measured."
  Treat D = 0.103 in and L = 0.167 in as **estimate/uncertain provenance**,
  not confirmed drawing dimensions.
- **LOX post recess/protrusion relative to the face:** not quantified in
  either primary source. EHS (variant unspecified) describes elements as
  "machined tubes projecting from the center plate" that "serv[e] as
  oxidizer passages and simultaneously create fuel annuli," and separately
  describes elements "protruding through" the porous Rigimesh faceplate —
  qualitatively suggesting the LOX post extends to or through the face
  plane, but with no recess/protrusion length given. TM-107318 Figure E2
  (p. 142) shows the post terminating near the Rigimesh plate in schematic
  form only, also without a dimension.
- Per-element orifice count: TM-107318 Appendix E's equivalent-swirler model
  (Figure E9) assumes **3 tangential orifices per element** for its ROCCID
  surrogate, but explicitly flags this as an arbitrary modeling choice, not
  a hardware count: "We have chosen 3 orifices here but other values could
  be chosen. The number of orifices could be chosen, for example, to better
  match the overall resistance of the actual ribbon-swirler LOX injectors"
  (p. 146). The real hardware uses a continuous helical ribbon insert, not
  discrete orifices, so this "3" does not describe the real element.

## 4. Face and transpiration

- The injector face (innermost of three conical plates, "cone 1") is a
  **porous, welded steel-mesh plate** — PWA FR-1769 (1966) calls it
  specifically "porous-welded, steel-mesh plate" (p. I-19); TM-107318 and
  Di Matteo et al. (2012) both call the material **"Rigimesh"** (TM-107318
  Appendix E p. 136; Di Matteo et al. 2012, injector-plate section). EHS
  gives a materials description: Rigimesh "a material formed by hot-rolling
  stacks of Type 347 stainless steel mesh together," noting it replaced an
  earlier sheet-metal face design that had warping problems (variant
  unspecified).
- **Fraction of hydrogen through the porous face vs. the element annuli:**
  PWA FR-1769 (1966), p. I-19 (variant A-3-3): of the gaseous hydrogen
  entering the fuel manifold, most exits through the annular orifices around
  the elements; some flows past and cools the igniter sleeve before entering
  the igniter chamber; "the rest of the hydrogen passes through cone 1 ...
  This flow provides transpiration cooling of the injector face (cone 1) and
  amounts to approximately **10% of the total hydrogen flow**." This is a
  flow-split fraction (not an absolute flow rate or pressure), stated as a
  design figure in a primary design report. TM-107318's own modeling note
  (Appendix E, p. 136) says "the faceplate transpiration cooling flow was
  not modeled for the RL10 [A-3-3A system model]" — i.e., TM-107318 does not
  give its own number for the A-3-3A and effectively defers to this type of
  figure without restating it; no separate A-3-3A-specific fraction was
  found.

## 5. Igniter

- **Type:** augmented spark igniter (TM-107318 §2.5, p. [body text], A-3-3A)
  / spark igniter with a dedicated electrode (PWA FR-1769 (1966) §I.K,
  p. I-22, A-3-3).
- **Location:** centered on the injector face. TM-107318 (A-3-3A): "An
  augmented spark ignitor is located in the center of the injector face."
  EHS (variant unspecified): "A spark igniter recessed in the propellant
  injector face initiates combustion chamber ignition." PWA FR-1769 (1966,
  A-3-3) Figure 1-18 shows the igniter sleeve on the injector centerline
  axis.
- **Propellant supply path (qualitative, no flow values):** TM-107318
  (A-3-3A): gaseous oxygen is supplied to the igniter from a point upstream
  of the injector inlet manifold; igniter fuel comes from the fuel injector
  plenum through holes in the igniter wall; fuel and oxygen mix and are
  ignited by an electric spark. PWA FR-1769 (1966, A-3-3) §I.K and Figure
  1-22 ("Igniter Assembly," drawing FD 3134A): fuel and oxidizer are fed
  into an annulus surrounding the spark igniter; oxidizer enters the annulus
  from a line connected to the igniter oxidizer supply valve; fuel enters
  the annulus from the injector through fuel metering orifices in the
  igniter plug housing; the igniter electrode has "many sharp corners" to
  reduce/stabilize the required spark-gap breakdown voltage.
- No igniter orifice/annulus dimensions were found in either source.

## 6. Distances

- **No face-to-throat axial distance, and no element-to-wall (outermost
  ring to chamber wall) distance, was found** in either primary source or
  in the open literature searched. TM-107318's cooling-jacket figures
  (Appendix D) plot axial position "with respect to throat" but these are
  heat-transfer/pressure/temperature plots, out of scope and not used here
  even for their axial-position axis, since reading geometry off a
  performance figure risked pulling plotted values along with it.

## Missing (not found in the time available)

- LOX post inner diameter, outer/tip diameter (actual hardware dimension).
- Fuel annulus inner/outer diameter or radial gap (actual hardware
  dimension).
- Ring radii or radial pitch between the 8 concentric rows, and element
  pitch/spacing within a ring.
- Confirmed element count per ring (the "27 per row" figure used above to
  reconcile 162 vs. 216 is an inferred even split, not a sourced count).
- LOX post recess or protrusion length relative to the face plane.
- Real (as opposed to ROCCID-equivalent) ribbon-swirler lead angle / number
  of ribbon starts.
- Igniter orifice sizes, annulus gap, or electrode gap dimension.
- Face-to-throat distance and outermost-ring-to-wall distance.
- A-3-3A-specific (as opposed to A-3-3, 1966) transpiration-fraction figure,
  and confirmation that the 10% figure still applies to the A-3-3A.
- Any dimensioned engineering drawing of the injector face or element
  (only non-dimensioned schematics were located: TM-107318 Figures E2, E9;
  PWA FR-1769 (1966) Figures 1-18, 1-22).

## Withheld (labels only, no values)

- TM-107318 Table 2.5.1 ("Summary of Combustion-Chamber / Nozzle
  Characteristics"): chamber pressure, combustion temperature, mixture
  ratio, mass flow, gross thrust, c*-efficiency, c*, specific impulse rows —
  withheld in full (geometry cells for chamber diameter, throat diameter,
  and nozzle area ratio also withheld here since out of scope for an
  injector-only extraction and already in hand per the case setup).
- TM-107318 Appendix E narrative sentence on using measured LOX-element
  pressure drop to adjust the ROCCID equivalent-orifice estimate — withheld
  in full (injector pressure drop).
- TM-107318 Appendix E Figures E7/E8/E10/E11 and surrounding text
  (c*-efficiency and specific-impulse maps/predictions vs. P&W data) —
  withheld.
- TM-107318 Appendix D Figures D6/D7 (cooling-jacket coolant static pressure
  and temperature predictions vs. axial position) — withheld.
- PWA FR-1769 (1966) ignition-system section: exciter-assembly internal
  pressurization figure — withheld.
- Various mass-flow, c*, and performance mentions scattered through
  TM-107318 body text (engine mass flows, turbopump performance,
  valve-actuation pressures) — withheld throughout; not itemized line by
  line since none were geometry-adjacent.
