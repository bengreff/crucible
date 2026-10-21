# RL10A-3-3A: Table 6.1.1 / Figure E1 addendum (blind extraction)

Source: Binder, M., Tomsik, T., and Veres, J. P., "RL10A-3-3A Rocket Engine
Modeling Project," NASA TM-107318, January 1997. NTRS ID 19970010379. Local
file: `archive/pre-pivot-2026-09-17/data/anchors/TM-107318_ntrs_19970010379.pdf`.
Printed page N is PDF page N+6, confirmed against printed page-footer numbers
at every location cited below.

This file supplements `docs/validation/RL10A33A_INPUTS.md`. It fills in the
Table 6.1.1 temperature/density cells that file withheld as a transcription
precaution, scales Figure E1 for the injector/exit axial stations, and
records a text search for chamber-length-type dimensions. Blind rule in
effect throughout: no chamber pressure, c*, thrust, Isp, efficiency, heat
load/flux/rate, temperature rise, or any pump/turbine performance value
appears anywhere below, even implicitly.

---

## A. Table 6.1.1 — Fuel-Side / Oxidizer-Side Engine Conditions (temperature, density only)

Table 6.1.1, "Steady-state Cycle Performance Prediction for the RL10A-3-3A,"
printed p.37 (PDF p.43). The page is a landscape scan; rotating the rendered
image -90 deg (`pdftoppm -r 300` then rotate) makes it readable. Pressure
columns in this table are withheld per the blind rule; only temperature and
density are transcribed. Row labels are exactly as printed (note the
oxidizer-side rows are "Oxygen Tank" / "Ox Pump Inlet" / "Oxid. Injector
Plenum," not literally "Oxidizer Tank" etc., but they are the rows the task
asked for).

| Row (as printed) | Temperature (total, deg R) | SI (K) | Density | Unit as printed | SI (kg/m^3) |
|---|---|---|---|---|---|
| Fuel Tank | 38.60 | 21.44 | 4.337 | lb/cu.ft | 69.47 |
| Fuel Pump Inlet | 38.6 | 21.44 | 4.337 | lb/cu.ft | 69.47 |
| Fuel Injector Plenum | 359.4 | 199.67 | 0.2742 | lb/cu.ft | 4.392 |
| Oxygen Tank | 174.7 | 97.06 | 69.07 | lb/cu.in (as printed — see caveat) | see caveat |
| Ox Pump Inlet | 174.7 | 97.06 | 69.07 | lb/cu.in (as printed — see caveat) | see caveat |
| Oxid. Injector Plenum | 182.1 | 101.17 | 68.31 | lb/cu.in (as printed — see caveat) | see caveat |

Source: Table 6.1.1, printed p.37 (PDF p.43).

**Unit-label caveat (oxidizer-side density column):** the table's printed
column header reads "density (lb/cu.in)" for the oxidizer-side block,
against "density (lb/cu.ft)" for the fuel-side block — verified by rotating
and re-reading the image at 300 dpi; the header text is legible and says
"in," not "ft." Taken literally, 69.07 lb/in^3 = 1.91e6 kg/m^3, which is
unphysical (denser than osmium) for liquid oxygen. Read instead as the same
lb/cu.ft unit used on the fuel side (almost certainly a typesetting error
in the original 1997 report, carried through uncorrected), 69.07 lb/ft^3 =
1106.4 kg/m^3, 69.02 lb/ft^3 = 1105.9 kg/m^3, 68.31 lb/ft^3 = 1094.2 kg/m^3
— all physically reasonable for saturated/subcooled liquid oxygen at the
tank/pump-inlet temperatures and pressures in this table. I transcribe the
value and unit exactly as printed and flag the correction rather than
silently applying it.

### Geometry sub-rows (hot-gas-side surface area, hot wall thickness, coolant volume)

These three items are **not** in Table 6.1.1 — I read all of Table 6.1.1
(rotated, full page) and it has no such sub-rows; it contains only the
"Overall Engine Performance," "Turbomachinery Performance," "Combustion
Chamber / Nozzle," "Fuel-Side/Oxidizer-Side Engine Conditions," "Cooling
Jacket Performance," and "Valve Performance" blocks already covered by the
companion file. The three requested geometry items are instead in a
different table, **Table 2.4.1, "Summary of Cooling Jacket
Characteristics,"** printed p.5 (PDF p.11), which is where this addendum
gets them (transcribed as geometry, not performance — the table's other
rows, Pressure Drop/Mass Flow/Temperature Rise/Total Heat Transfer Rate,
are withheld):

| Item (as printed) | Value | Units | SI |
|---|---|---|---|
| Total Coolant Volume | 966.1 | cu.in | 0.015832 m^3 |
| Typical Hot Wall Thickness | 0.013 | in | 3.302e-4 m (0.3302 mm) |
| Hot-gas-side Effective Surface Area | 7200 | in^2 | 4.6452 m^2 |

Source: Table 2.4.1, printed p.5 (PDF p.11).

---

## B. Figure E1 — axial stations of injector face and nozzle exit plane

Figure E1, "Chamber and Nozzle Profile with Property Node Locations,"
printed p.141 (PDF p.147). Rendered at 300 dpi (`pdftoppm -r 300`); the page
is portrait, not rotated.

### Method

The plot's x-axis ("Axial Position wrt Throat (inches)") carries printed
gridlines at -15, -5, 5, 15, 25, 35, 45, 55. I located each gridline's pixel
column by thresholding the image to black/white and summing dark pixels
down each column over the plot's full vertical extent (430-1185 px rows),
isolating the tall, continuous vertical gridlines from the shorter text
glyphs. Peak columns: -15 at x=458, -5 at 688, 5 at 918, 15 at 1151, 25 at
1382, 35 at 1614, 45 at 1849, 55 at 2081 (px). A linear fit through the
endpoints gives 23.186 px/in and places x=0 (the throat) at pixel 805.8;
the plot's own dash-dot "Throat" gridline peaks at pixel 804 independently
— a 1.8 px (0.08 in) self-consistency check.

I located the "Injector" and "Nozzle Exit" dash-dot reference lines the same
way (column-wise dark-pixel peak, away from the rotated text labels that
run alongside each line): Injector peaks at pixel 503, Nozzle Exit at pixel
≈1871.5.

- **(1) Injector face:** pixel 503 → x = (503 - 805.8) / 23.186 = **-13.1
  in** relative to the throat (upstream). Cross-check: Table E1's most
  upstream property node is at -12.00 in (area ratio 4.27); the figure's
  injector line sits about 1 in further upstream of that node, consistent
  with the visual (the first open-circle marker is a short distance
  downstream of the Injector line in the image).
- **(2) Nozzle exit plane:** pixel 1871.5 → x = (1871.5 - 805.8) / 23.186 =
  **+46.0 in** relative to the throat (downstream). Cross-check: Table E1's
  most downstream node is at +41.84 in (area ratio 56.1, still short of the
  design 61.0 exit ratio from Table 2.5.1); the figure's Nozzle Exit line
  sits visibly downstream of that last node, consistent with an exit
  further out than the last tabulated node.
- **Uncertainty:** the dash-dot lines are 4-5 px wide in the render (half
  of that locates the center to about ±2-3 px), and the gridline calibration
  itself has about ±1-2 px of scatter (spacing 230-235 px per 10 in rather
  than a perfectly uniform 231.9). Combined, I estimate **±0.3 in** for the
  injector face and **±0.3-0.4 in** for the nozzle exit (the exit line sits
  in a denser, more textured part of the image next to the rotated "Nozzle
  Exit" label, slightly reducing confidence there).

### (3) Chamber shape, injector face to first Table E1 node (-13.1 in to -12.0 in)

Over this short span the plotted inner-wall curve is essentially flat: the
curve value at the first node (x=-12, pixel row 922) and the next marker
toward the throat change only gradually (markers near -13 to -8 in read
within about 5% of each other — see the "Absolute size check" section below
for the actual pixel-derived values). Visually this is a **cylindrical
(constant-radius) chamber section** near the injector, not a taper or
rounded corner; the taper visibly begins and steepens further downstream,
roughly from x=-8 in toward the throat, where the curve bends over into the
convergent section. No rounded corner is drawn at the injector end in
Figure E1 — the curve simply starts (left edge of the data) at the Injector
line already at its near-flat chamber value.

### (4) Throat detail

Figure E1 includes an inset, "Detail of Throat (not to scale)," showing the
converging and diverging "Tube Profile" meeting in a V, with a "Cast Silver
Throat Insert" (shown as a filled black wedge) occupying the valley between
the two straight tube segments — i.e., the insert rounds/fills what would
otherwise be a sharp corner, rather than the tube wall itself being cast
with a sharp edge. **No numeric radius is printed on this inset**; the
caption explicitly disclaims scale ("not to scale"). Supporting narrative
text (Appendix D.2, printed p.119/PDF p.125 — quoted in full in the
"Absolute size check" section below) confirms the insert's role: it "reduces
the effective throat diameter formed by the coolant tubes and increases the
expansion area ratio from 57:1 to 61:1," i.e., without the insert the
coolant-tube contour alone would give a 57:1 ratio, and the insert narrows
the throat to reach the 61:1 design ratio from Table 2.5.1 — geometry, not
performance.

### Other printed dimensions on Figure E1

None beyond the axis scales and the two dash-dot reference lines already
used above; the main plot carries no printed numeric lengths, radii, or
angles beyond the axis tick values, and the throat-detail inset carries no
printed dimension at all (qualitative labels only).

---

## C. Chamber length / L* / nozzle length / convergent-divergent angle — text search

Searched the full extracted text (`pdftotext -layout`) of printed pp.3-8
(PDF pp.9-14, section 2) and printed pp.135-150 (PDF pp.141-156, Appendix
E/F) for: "length," "L*"/"L-star"/"characteristic length," "chamber
length," "nozzle length," "convergent angle," "divergent angle," "half
angle," "cone angle," "injector-to-throat," case-insensitively, then
re-checked every hit by eye.

**None of chamber length, injector-to-throat distance, L*, nozzle length,
or a convergent/divergent half-angle is printed anywhere in either page
range.** The only geometry-flavored "length" hits in these ranges are the
already-known Table 2.5.1 diameters (printed p.6/PDF p.12) and a single
unrelated item: Appendix F.3 ("Heat Transfer across plenum bulkhead,"
printed p.150/PDF p.156 — at the very edge of, and arguably just past, the
pp.135-150 Appendix E range given in the task) states a flow-path length
"Lplm = 1.947 inches" used in a bulkhead heat-transfer correlation, not a
chamber/nozzle axial length; see the verbatim quote in the next section.
Section 2's Table 2.5.1 (printed p.6/PDF p.12) likewise has no length/angle
row beyond the two diameters and the area ratio already captured in the
companion file.

Figure E1 itself (Part B above) is the only place an injector-to-throat-type
distance can be read at all, and only by scaling the drawing, not from a
printed number.

---

## Absolute size check (requested follow-up)

Goal: using only the report's own printed geometry, decide whether the
2.47 in / 5.13 in numbers in Table 2.5.1 are diameters or radii. No
performance value is used or computed below.

### (a) Table 2.5.1's exact row labels

Table 2.5.1, "Summary of Combustion-Chamber / Nozzle Characteristics,"
printed p.6 (PDF p.12). Verbatim (via `pdftotext -layout`, cross-checked
against the rendered page image):

> Chamber Diameter &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 5.13 in
>
> Throat Diameter &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 2.47 in
>
> Nozzle Area Ratio &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 61.0

Both rows say **"Diameter,"** verbatim, not radius or any other term.

A second, independent mention uses the same word. Appendix D.2, "Cooling
Jacket Description," printed p.119 (PDF p.125):

> "A cast silver throat insert is installed in the throat region. This
> insert reduces the effective throat diameter formed by the coolant tubes
> and increases the expansion area ratio from 57:1 to 61:1."

This again calls the dimension a **diameter**, and ties the 57:1 → 61:1
area-ratio change (Table 2.5.1's 61.0 figure) to a reduction in that
diameter by the insert — internally consistent with Table 2.5.1 on its own
terms.

**Item (a) supports: diameter**, and it is the strongest single piece of
evidence here (explicit label, stated twice, independently, in different
sections).

### (b) Other dimensions that fix absolute size

**Figure E1's radial axis and plotted curve** (printed p.141/PDF p.147).
The axis is labeled, verbatim: "Radius of Chamber/Nozzle Inner Wall
(inches)" — this label says **radius**. But the plotted curve's own
numeric values, read by the same pixel-calibration method as Part B
(y=20 at row 588.1, y=0 at row 1030.0, 22.095 px/in; both independently
cross-checked against the figure's own labeled gridlines):

- At the first Table E1 node, x=-12 in (near the injector, where the wall
  is still essentially cylindrical): curve marker centered at row 922 →
  value **4.89** (my digitization). Compare Table 2.5.1's Chamber Diameter
  5.13 in (ratio 0.95, i.e. a 5% difference) vs. the chamber *radius* it
  implies, 2.565 in (ratio 1.91, i.e. essentially double).
- At the throat, x=0 in: curve vertex centered at row ≈982 → value
  **2.17** (my digitization, lower confidence — the vertex sits where two
  closely-spaced nodes, the dash-dot throat line, and the throat-detail
  pointer lines all overlap). Compare Table 2.5.1's Throat Diameter 2.47 in
  (ratio 0.88) vs. the throat *radius* it implies, 1.235 in (ratio 1.76).

At both independent stations, the figure's own plotted number lands close
to the Table 2.5.1 **diameter** value (within ~5-12%, plausibly digitization
error on a hand-drawn 1997 plot) and far from half that value. **Taken at
face value, Figure E1's axis label ("Radius") is inconsistent with its own
plotted numbers**, which track Table 2.5.1's diameters instead. I can't
fully resolve this from the report alone — either the axis label is a
captioning error (should read "Diameter"), or my digitization of a hand-
drawn curve is systematically off by close to 2x at two independent
stations, which seems less likely given the consistency. I report this as
an open inconsistency, not a resolved fact; see the numeric detail above so
a reviewer can judge independently of my reading.

**Coolant tube count/geometry.** Appendix D.2 (printed p.119/PDF p.125) and
section 2.4 (printed p.5/PDF p.11) both state 360 total tubes (180 short +
180 long, "pass-and-a-half" design) forming the hot-gas wall, and Table
2.4.1 gives "Typical Hot Wall Thickness" = 0.013 in (already transcribed in
Part A above). No printed tube width, OD, or ID was found anywhere in the
searched text (Appendix D's heat-transfer equations use a symbolic hydraulic
diameter D_c with no numeric value given, and the referenced "Figure D-1"
tube-contour drawing, near printed p.119-121/PDF p.125-127, carries no
printed dimensions in either the extracted text or the rendered page
images I checked). This item is **inconclusive** — it doesn't fix absolute
size on its own.

**Table 6.1.1's geometry sub-rows.** As found in Part A above, these rows
are not in Table 6.1.1 at all; they are Table 2.4.1's Total Coolant Volume
(966.1 cu.in), Typical Hot Wall Thickness (0.013 in), and Hot-gas-side
Effective Surface Area (7200 in^2) — already transcribed in Part A. I did
not attempt to back out a chamber radius from the 7200 in^2 surface area
figure; doing so needs the full wall profile integrated over its true slant
length, which is beyond what a hand check can respectably do, and I did not
want to manufacture false precision for this question.

**Printed throat area, exit diameter, or exit area.** None found anywhere
in the searched text.

**Item (b) supports: the figure's plotted numbers favor "diameter" too**
(they match Table 2.5.1's diameter values, not half of them) — but this
comes with the caveat that the figure's own axis label says "radius," so
this is evidence of an internal inconsistency in the report, read as mild
additional support for "diameter" being the physically intended quantity,
not a clean confirmation.

### (c) The face-plate "length of plate" sentence, verbatim

Found in Appendix F.3, "Heat Transfer across plenum bulkhead," printed
p.150 (PDF p.156) — immediately following the printed page-150 footer, so
I read it as still belonging to that page:

> "In this case, the flowpath length along the plate is taken to be half
> the radius of the injector bulkhead (Lplm = 1.947 inches). Half the
> radius is used to account for the fact that fuel and LOX flow are both
> bled off by the injector orifices as the flows pass over the plate."

This explicitly states Lplm (1.947 in) **is half of** the injector
bulkhead's radius, i.e., it is a direct, unambiguous statement that the
injector bulkhead radius = 2 x 1.947 = **3.894 in** — this is a radius, by
the sentence's own wording, not a diameter and not some other length. It
says nothing directly about the chamber or throat dimension; it only fixes
the injector bulkhead's own size.

**The inconsistency this creates:** an injector bulkhead of radius 3.894 in
does not fit inside a chamber of radius 2.565 in (the radius implied if
Table 2.5.1's 5.13 in "Chamber Diameter" is read literally as a diameter).
It does fit inside a chamber of radius 5.13 in (the value if that number
is read as a radius instead). So **item (c), combined with the chamber
dimension, favors "radius"** for 5.13 in — unless the "injector bulkhead"
named in this sentence is some smaller internal partition rather than the
full injector face plate (the text doesn't make this distinction explicit
either way; it only ties Lplm to "the injector bulkhead," without further
describing the bulkhead's extent relative to the chamber wall).

### Summary of which interpretation each item supports

- (a) Table 2.5.1's row labels, verbatim, both places found: **diameter**
  (the most direct, most repeated, most explicit evidence).
- (b) Figure E1's plotted curve values at two stations: numerically match
  **diameter** (not radius, despite the axis itself being labeled
  "Radius" — an unresolved internal inconsistency I'm flagging, not
  resolving).
- (b) Coolant tube width/OD, printed throat/exit area: **inconclusive**,
  nothing printed.
- (c) The injector-bulkhead "length of plate" sentence: fixes the bulkhead
  itself as a **radius** (3.894 in) by its own wording, and is only
  physically consistent with the chamber it sits in if the chamber number
  is *also* read as a radius (5.13 in) rather than a diameter.

I'm not able to fully settle this from the report alone: the explicit
textual labels (two independent instances) say diameter, while the two
numeric/geometric cross-checks (Figure E1's curve values, and the
injector-bulkhead fit) both point toward radius. I've laid out every
quote and page so this can be judged directly rather than taking my own
read as final.
