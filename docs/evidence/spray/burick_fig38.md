# Burick (1972), Figure 38: coaxial-element drop-size correlation

Source: R. J. Burick, "Space Storable Propellant Performance Program: Coaxial
Injector Characterization," NASA CR-120936 (Rocketdyne report R-9220,
prepared for NASA Lewis Research Center, 1972).

PDF: https://ntrs.nasa.gov/api/citations/19720026081/downloads/19720026081.pdf
NTRS citation page: https://ntrs.nasa.gov/citations/19720026081

Page numbers below are the printed page numbers at the bottom of each page
(these equal PDF page number minus 8 for this scan). Both are given where it
matters.

This is simulant (wax / GN2) cold-flow data from Burick's report only. No
Mascotte, ONERA, A-10/A-60, or other LOX/H2 combustor measurement was read or
used for this file.

## 1. Figure 38: axes, correction factors, digitised curve and scatter

### 1.1 Identification and axis definitions

The figure is confirmed as **Figure 38, "Correlation of Cold-Flow Atomization
Data"**, printed page 67 (PDF page 75). It is the only coaxial-element
drop-size correlation plot in the report.

Caption text (Section 5.3.2, printed page 56, PDF page 64):

> "...the resultant mean drop sizes were nondimensionalized by D_L. Figure 38
> presents the parameter D-bar/D_L as a function of (V_g - V_L)/(V_L MR). As
> indicated in the figure, the parameters provide a reasonable correlation of
> the atomization data."

Axis definitions, exactly as printed on the figure and in the surrounding
text:

- Ordinate: **D-bar / D_L** (dimensionless), log scale, printed range 0.06 to
  1.0 (gridlines at 1.0, 0.6, 0.4, 0.2, 0.1, 0.06).
  - D-bar = mass median drop diameter of the spray (see Section 2 below),
    microns.
  - D_L = liquid (oxidizer) post diameter, inches (the same D_L used
    throughout the report for the coaxial element's liquid orifice).
- Abscissa: **(V_g - V_L) / (V_L * MR)** (dimensionless), log scale, printed
  range 0.1 to 20.
  - V_g = gas (annulus) velocity at the injection plane, ft/sec.
  - V_L = liquid (oxidizer) velocity at the injection plane, ft/sec.
  - MR = mixture ratio (gas flow / liquid flow) at the element.

Legend: open triangle = "Flush Post" (zero recess, R = 0); filled triangle =
"Recessed Post" (R > 0).

The derivation of this correlating group is given in Section 5.3.2 (printed
page 56): the mixing correlation used (rho_g * V_g)^2 / (MR * V_L); for
atomization "no reasonable correlation was obtained" with that group, so the
numerator was replaced by a term proportional to the gas-liquid shear rate,
(V_g - V_L), and the result nondimensionalized by D_L. No other property
(viscosity, surface tension, density) appears explicitly in either axis of
Figure 38 itself.

### 1.2 No equation is given for the curve

Section 5.3.2 (quoted above and in full in Section 1.1) describes the
correlation only in words ("the parameters provide a reasonable correlation
of the atomization data"). No closed-form equation for the dashed curve in
Figure 38 is printed anywhere in the text, in Section 5.3.2, in the figure
caption, or in the summary/conclusions sections. The curve is graphical only.

### 1.3 Property-correction factors: none on Figure 38 itself; applied downstream

Figure 38's two axes carry no built-in viscosity, surface-tension, or
density correction; they are plotted directly from the wax/GN2 cold-flow
data of Table IV. The property correction from wax/GN2 to a real propellant
is applied as a separate, later step, described in two places:

- Section 6.1.1 (printed page 53, PDF page 61): "In the analysis, the wax
  drop sizes (D-bar) were corrected for the difference in physical
  properties of wax and FLOX using the empirical relations of Ingebo, Ref.
  24 (see Section 7.0)."
- Section 7.1 (printed page 91, PDF page 100): "To predict a mean FLOX drop
  size for the hot-fire system the mean drop size data (D-bar_wax) of Fig.
  52 were corrected for the difference in physical properties of Shell-270
  wax and FLOX. The values employed for the respective physical properties
  of 82.6% FLOX and Shell-270 wax are shown in Table VIII."

Table VIII (printed page 91, PDF page 100), "PHYSICAL PROPERTIES OF FLOX AND
WAX":

| Property | 82.6% FLOX | Shell-270 Wax |
|---|---|---|
| Viscosity (mu), centipoise | 0.23 | 4.0 |
| Surface Tension (sigma), dyne/cm | 13.8 | 17.5 |
| Density (rho), lbm/ft3 | 89 | 47.1 |

Footnoted sources: FLOX properties interpolated from NASA SP-3037,
"Handling and Use of Fluorine-Oxygen Mixtures in Rocket Systems"; wax
properties from a private communication (Dannenbrink, Shell Chemical Co., to
L. Zajac, Rocketdyne). Text comment: "Shell-270 wax simulates reasonably
well the surface tension of FLOX but differs in density and absolute
viscosity."

Reference 24 is Ingebo, R. D. (NACA, drop-size/breakup correlation; listed
in the report's reference list as the crosscurrent liquid breakup dropsize
paper). I searched the pages around Table VIII (PDF pages 92-100, printed
84-91) and the pages that cite Ingebo directly (printed pages 53, PDF 61; and
printed page 85, PDF 93) for the Ingebo correction formula itself. It is
referenced by name three times (printed pages 53, 85, and 91) but its
numerical form is **not reproduced anywhere in the body text of this
report**: the report only cites Table VIII's property values and says the
correction follows "the empirical relations of Ingebo, Ref. 24." I was not
able to find the equation itself in CR-120936; it would have to be obtained
from Ingebo's own paper (Ref. 24), not from Burick's report.

### 1.4 Digitised curve

Method: the dashed correlation line was located in a 600 dpi render of
printed page 67 (PDF page 75, file rendered via `pdftoppm -png -r 600`).
Pixel calibration was obtained by directly measuring the plot's border and
gridline pixel positions in that 600 dpi image (image size 6144 x 8712 px):

- Plot borders measured at: left edge cols 1133-1141, right edge cols
  5133-5141, top edge rows 3495-3502, bottom edge rows 5805-5816.
- Gridlines checked at ordinate values 1.0, 0.6, 0.4, 0.2, 0.1, 0.06 and
  abscissa values 0.1 and 20; all matched the formulas below to within a few
  pixels.

Calibration formulas (pixel coordinates in the full 600 dpi page image):

- Ordinate (D-bar/D_L): `y_px = -1776.83 * log10(v) + 3490.99`, i.e.
  `v = 10^((3490.99 - y_px) / 1776.83)`
- Abscissa ((V_g-V_L)/(V_L*MR)): `x_px = 1737.83 * log10(v) + 2871.32`, i.e.
  `v = 10^((x_px - 2871.32) / 1737.83)`

The dashed curve was traced column by column after masking the plot's own
dashed gridlines (rows/columns with high dark-pixel fraction) and the legend
box, then converted to data values with the formulas above and interpolated
in log-log space to a set of round x values. This was cross-checked by
re-locating the curve's start and a mid-curve point independently in a
second pass; the two passes agreed to within about 1% in D-bar/D_L.

**Digitised** curve, D-bar/D_L vs. (V_g-V_L)/(V_L*MR):

| x = (V_g-V_L)/(V_L*MR) | D-bar/D_L |
|---|---|
| 0.23 | 0.277 |
| 0.28 | 0.276 |
| 0.35 | 0.272 |
| 0.45 | 0.265 |
| 0.55 | 0.260 |
| 0.65 | 0.252 |
| 0.75 | 0.247 |
| 0.85 | 0.242 |
| 1.00 | 0.241 |
| 1.30 | 0.238 |
| 1.60 | 0.248 |
| 2.00 | 0.243 |
| 2.30 | 0.225 |
| 2.60 | 0.209 |
| 2.80 | 0.206 |
| 3.20 | 0.169 |
| 3.60 | 0.165 |
| 4.00 | 0.159 |
| 4.50 | 0.154 |
| 5.00 | 0.148 |
| 5.50 | 0.154 |
| 6.00 | 0.129 |
| 6.50 | 0.133 |
| 7.00 | 0.129 |
| 7.70 | 0.125 |
| 8.50 | 0.118 |
| 9.50 | 0.113 |
| 10.50 | 0.112 |
| 11.50 | 0.113 |
| 12.50 | 0.101 |
| 13.50 | 0.097 |
| 14.50 | 0.094 |

Endpoints: the drawn curve starts at about x = 0.23-0.25, D-bar/D_L = 0.28,
and ends at about x = 14.5-15, D-bar/D_L = 0.09-0.10. Beyond about x = 15 the
trace becomes noisy (stray marks, not a continuation of the dashed line);
I stopped the table there.

Caveat: the x = 0.5-1.6 stretch of the table runs through the two densest
clusters of data-point markers (see 1.5), where triangle markers overlap the
dashed line itself. The local small rise at x = 1.60 (0.248 vs. 0.238 at
1.30) is most likely a tracing artefact from a marker sitting on the line
rather than a real local maximum in the printed curve; treat the x = 0.5-1.6
values as accurate to roughly +/-0.01-0.02 in D-bar/D_L rather than the
better precision of the clean, isolated stretches of the curve (x < 0.4 and
x > 2.5).

### 1.5 Digitised scatter (data points)

Every visible marker could not be separated reliably by pixel analysis:
Figure 38 has two dense clusters (around x = 0.7-1.0 and x = 1.5-2.0) where
3 to 6 open and filled triangles overlap each other and the dashed curve.
**Method used: upper/lower envelope of the scatter at a set of x positions**
(the task's second option), read from a combination of (a) automated
connected-component detection of isolated markers in the same masked,
600 dpi pixel field used for the curve (shape filter: roughly square
bounding box, height and width 8-45 px, fill ratio > 0.45, converted with
the same calibration formulas as Section 1.4), and (b) direct visual
reading of 2x-enlarged crops of the plot against the gridlines, for the two
dense clusters where automated separation failed. Isolated, unambiguous
points are reported as single values; clusters are reported as a y-range.

**Digitised** scatter envelope (both Flush Post and Recessed Post symbols
pooled; D-bar/D_L):

| x position | D-bar/D_L range | basis |
|---|---|---|
| ~0.30 | 0.24 - 0.28 | 2 isolated markers |
| ~0.50 | 0.26 - 0.33 | isolated markers, includes one visibly high "Flush Post" point above the curve |
| ~0.7 - 0.9 | 0.14 - 0.25 | dense cluster, 5-6 markers, envelope by eye + automated points |
| ~1.5 - 2.0 | 0.15 - 0.26 | dense cluster, about 6 markers, envelope by eye |
| ~3.0 | 0.18 - 0.20 | isolated pair |
| ~3.9 - 4.2 | 0.157 - 0.160 | 2 isolated markers |
| ~5.4 | 0.155 | 1 isolated marker |
| ~7.1 - 9.0 | 0.052 - 0.115 | 2 isolated markers; low end flagged uncertain, see below |
| ~10 - 11 | 0.105 - 0.135 | 2 isolated markers, read directly off the full-page figure |

Uncertain point: the low value at x = 7.1, D-bar/D_L = 0.052, is a
well-formed, high-confidence marker detection (compact, high fill ratio)
but it sits well below the local curve trend (curve value near x = 7.7 is
about 0.125, roughly 2.4x higher). It may be a genuine low scatter point or
a digitisation artefact near a gridline crossing; I was not able to resolve
this from the scan alone and flag it as uncertain rather than excluding it.

Overall scatter width: across the well-populated mid-range (x = 0.5 to 2),
most points fall within roughly a factor of 1.3 above and a factor of 1.7
below the curve (about +30%/-40% in D-bar/D_L). A few points, including the
x = 7.1 outlier above, lie further below the curve, by up to about a factor
of 2.4.

## 2. Drop-size distribution form

### 2.1 Distribution and its parameters

Section 5.3.1, "Drop Size Distribution Correlation" (printed page 62, PDF
page 70):

> "Figure 37 presents normalized drop size distribution data from several of
> the pressurized atomization experiments. Also shown in the figure is the
> normalized Rosin-Rammler drop size distribution function (Ref. 23). Note
> the excellent agreement with the coaxial injector data for values of
> D/D-bar > 1.0. Use of realistic drop size distribution functions for
> D/D-bar > 1.0 in combustion model programs is critical because these drop
> size ranges significantly influence the predicted eta_c*,vap."

The distribution is plotted only as a normalized cumulative mass fraction
curve against D/D-bar (Figure 37, printed page 65, PDF page 73), not given
as a printed algebraic Rosin-Rammler formula (no explicit exponent N or
spread parameter q appears anywhere in the report text that I could find).
The curve is credited to Reference 23, an internal/unpublished source (see
2.3), and is simply overlaid on the normalized data without a restated
equation.

D-bar is explicitly defined, right above Section 5.3.2 (printed page 62,
PDF page 70): "D-bar = mass median drop size." This is the quantity plotted
on the ordinate of Figure 38.

### 2.2 D30 (volume mean) relation

Equation (5), printed page 62 (PDF page 70):

> "The mass median drop size, D-bar, for a spray sample which conforms to
> the Rosin-Rammler distribution function may be converted to an equivalent
> volume mean drop size by the equation (Ref. 23):
>
> D30 = 0.455 * D-bar   (5)
>
> where
> D30 = volume mean drop size
> D-bar = mass median drop size"

No Sauter mean diameter (D32) relation is given in this section or found
elsewhere in the report.

### 2.3 Figure 37 data and the Rosin-Rammler reference curve

Figure 37, "Core-Element Normalized Drop Size Distribution Data" (printed
page 65, PDF page 73), plots cumulative mass fraction vs. D/D-bar for three
runs:

| Run No. | D-bar, micron | D_L, in. | D_g, in. | R/D_L |
|---|---|---|---|---|
| -11 | 580 | 0.136 | 0.182 | 0 |
| -31 | 348 | 0.136 | 0.182 | 1 |
| -40 | 398 | 0.070 | 0.070 | 1.25 |

A single smooth curve, labelled "Rosin-Rammler Distribution," is drawn
through the pooled, normalized data of these three runs; a vertical dashed
line at D/D-bar = 1.0 is labelled "Monodispersed" (i.e., D/D-bar = 1 would
be a single drop size). No numeric Rosin-Rammler exponent (commonly called N
or q in other sources) is printed on the figure or in the surrounding text.

Reference 23, cited for both the Rosin-Rammler curve and Equation (5), is an
internal Rocketdyne report (identified from the report's reference list as
L. Zajac, 1971, Rocketdyne report R-8455) that is not reproduced in
CR-120936; its contents (including any stated exponent value) are therefore
not available from this document.

## 3. Test conditions behind the correlation

Section 5.3, "Core Element Atomization Results" (printed page 56, PDF page
64):

> "Parametric atomization studies analogous to the mixing experiments were
> conducted employing the three element configurations (see Fig. 7). In all
> cases, the cold-flow models utilized in the mixing experiments were also
> utilized in the atomization experiments. All experiments were conducted in
> the pressurized atomization facility employing molten wax and gaseous
> nitrogen as nonreactive propellant simulants (see Appendix C)."

Simulants: molten Shell-270 wax (liquid/oxidizer simulant) and gaseous
nitrogen, GN2 (gas simulant), per Table IV's header: "Fuel--GN2,
Oxidizer--Shell-270 Wax; rho = 47.1 lbm/ft3."

Post (orifice) diameters, from Table IV (printed pages 68-69, PDF pages
76-77) and the narrative text:

- D_L (liquid/oxidizer post diameter): 0.070, 0.106, 0.112, and 0.136 inches
  appear in Table IV's D_L column. The narrative text (Section 5.3, printed
  page 56) describes "liquid jet diameters of 0.136 in., 0.108 in., and
  0.070 in." (0.108 in. in the text vs. 0.106/0.112 in. in the table; both
  are reported here as printed).
- D_g (gas annulus diameter): 0.136, 0.166, 0.182, and 0.228 inches appear in
  Table IV.
- D_p (post outer diameter): 0.080 to 0.146 inches appear in Table IV.
- R (post recess): 0 to 0.544 inches appear in Table IV (recess normalized
  by D_L reaches about 2.5-3 D_L in some runs per Section 5.3's "recessed to
  practical depths (recess/D_L = 2.5)").

Gas gap (annulus) velocity range: Figure 30 (printed page 57, PDF page 65)
shows gas gap velocity from 134 to 630 ft/sec, for h = 0.041, 0.018, and
0.005 in. gas gaps, at the design point W_L = 0.174 lbm/sec, W_g = 0.033
lbm/sec, rho_g = 1.45 lbm/ft3, D_L = 0.136 in., V_L = 36 ft/sec, MR = 5.25,
R = 0.

Liquid velocity: the Figure 30 design-point box gives V_L = 36 ft/sec for
the D_L = 0.136 in. element; other runs in Table IV span a range of liquid
and gas flow rates (W_liq from 0.0157 to 0.3480 lbm/sec; W_gas from 0.0030
to 0.0666 lbm/sec) that correspond to a range of V_L and V_g values not
individually tabulated as velocities (only as flow rates and orifice sizes)
in Table IV itself.

Mixture ratio: Section 5.3 (printed page 56) states "Figure 33 presents the
results for mixture ratios from 3 to 7.5," with design-point mixture ratios
of MR = 5.25 (gas-gap-velocity and throttling studies) and MR = 2.25 (post
recess study at design point rho_g = 1.45 lbm/ft3, W_tot = 0.207 lbm/sec).

Pressure: Table IV lists chamber/supply pressures from 0 psig (one run with
no gas flow, used as a reference case) up to 385 psig, most runs falling
between about 169 and 362 psig.

Resultant mean drop size D-bar, the quantity nondimensionalized to make
Figure 38's ordinate: Table IV lists values from 235 to 1077 microns across
all atomization runs (005-065).

Range of the correlating group (V_g - V_L)/(V_L*MR) actually covered by
plotted data in Figure 38: by direct reading and digitisation (Section 1.5),
roughly x = 0.3 to x = 11 for the data markers; the dashed correlation line
itself is drawn slightly beyond the data, out to about x = 15-16.

## 4. Extrapolation and LOX/GH2 applicability

Section 10.3, "Use of Cold-Flow Data for Other Propellants" (printed page
130, PDF page 139), quoted in full:

> "The correlated cold-flow data of this study were presented in a form (Fig
> 29 and 38) that can be employed to determine element mixing and
> atomization levels for other gas/liquid propellant combinations (i.e.,
> LOX/GH2, LOX/propane, etc.). In addition, the correlation of chamber wall
> heat flux versus chamber wall mixture ratio (Fig. 58) can be employed as a
> guide in designing elements for enhanced injector/chamber compatibility.
>
> However, caution must be exercised when applying these data to propellant
> combinations whose elements operate in ranges considerably different than
> those employed in this study. Specifically, coaxial elements employing
> LOX/GH2, in general, operate with significantly higher gas gap velocities
> (V_g = 2000 to 3000 ft/sec) than those using methane (V_g = 300 to 500
> ft/sec). In addition, LOX/GH2 elements generally have larger gas gap
> heights than those reported herein (h = 0.018 in.). Cold-flow studies
> (Ref. 2) currently underway at Rocketdyne with LOX/GH2 coaxial elements
> indicate that the cold-flow data generated in this study may not be
> directly applicable to LOX/GH2. However, injector designers concerned with
> propellants which operate with similar injection characteristics (i.e.,
> LOX/propane) should be able to utilize the results of this study
> directly."

Note on what this report's own full-scale/hot-fire propellants actually
were: this program's single-element and full-scale hot-fire tests used
FLOX (fluorine/LOX mixture, 82.6% F2 per Table VIII) as oxidizer and gaseous
methane (CH4) as fuel, not LOX/H2. Figure 38 itself is built entirely from
wax/GN2 cold-flow simulant data (Section 5.3); the report's own
applicability caution above is specifically about extending that same
wax/GN2 cold-flow correlation (not hot-fire FLOX/CH4 data) to LOX/GH2
elements, which it says "may not be directly applicable" because of the
much higher gas velocities and larger gas gaps typical of LOX/GH2 coaxial
injectors.

No other extrapolation statement specific to Figure 38's correlating-group
range (i.e., a caution about going beyond x = 0.1-20 on the abscissa) was
found in the report.
