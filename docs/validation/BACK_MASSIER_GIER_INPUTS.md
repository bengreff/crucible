# Back, Massier, and Gier: input extraction for blind validation

Source: Back, L. H., Massier, P. F., and Gier, H. L., "Convective Heat Transfer in a
Convergent-Divergent Nozzle" (also catalogued as "Convective Heat Transfer from
Turbulent Boundary Layer in Convergent-Divergent Conical Nozzle"), JPL Technical
Report No. 32-415, Jet Propulsion Laboratory, California Institute of Technology,
Pasadena, California, February 15, 1965. Also issued as NASA-CR-57326. NASA
Technical Reports Server (NTRS) ID 19650010083, https://ntrs.nasa.gov/citations/19650010083,
PDF at https://ntrs.nasa.gov/api/citations/19650010083/downloads/19650010083.pdf.

Date: 6 October 2026.

Inputs only, extracted under the blind rule: no measured heat flux or derived quantity appears here.

Page numbers below refer to the printed JPL Technical Report 32-415 page numbers (as
reproduced in the NTRS PDF, 62 pages total), not the PDF viewer page count. All axial
distances are quoted as the axial distance ratio z/L used by the report, with L = 5.925 in
(nomenclature, p. 22), together with the absolute distance I computed from it.

## Run selection

The rule: highest stagnation pressure at the highest stagnation temperature, in the
configuration with the longest approach (inlet) section.

Table B-1 ("Summary of tabulated tests," p. 29) lists, for each test, the cooled approach
length l, the stagnation pressure p_t, and the stagnation temperature T_t0. The approach
length takes four values, 0, 6, 12, and 18 in (stated in the Introduction, p. 1: "the length
of the constant-diameter cooled approach section upstream of the nozzle inlet was changed
in 6-in. lengths from 0 to 18 in."). The longest approach length is therefore 18 in. Table
B-1 marks l = 18 in explicitly against test 310, and the test block that contains it (tests
300, 301, 303, 304, 305, 306, 308, 309, 310, 311, 288, 287, 266, 281, 280, 263, 262, 269,
267, 268, 271, 273, 277, 278, 275, 293, 294, 290, 296, 298) is the only block that reaches
a stagnation temperature near 2000 R; all l = 0, 6, and 12 in tests in the table sit at
stagnation temperatures of about 1500 R or below (p. 29).

Within the l = 18 in tests, three stagnation-temperature groups appear: about 1030 R, about
1500 R, and about 2000 R (p. 29). The highest is the approximately 2000 R group, made up of
tests 293 (p_t = 30.0 psia, T_t0 = 2001 R), 294 (49.9, 2000), 290 (75.2, 1989), 296 (101.2,
2001), and 298 (125.8, 2007) (Table B-1, p. 29). The main text independently confirms this
is the ceiling for this group: "At the highest stagnation temperature, it was not possible to
obtain data above a stagnation pressure of 125 psia because of temperature limitations on
the wall-thermocouple insulating material" (Section VI, p. 10). The highest stagnation
pressure in that group is test 298, p_t = 125.8 psia (Table B-1, p. 29), restated as "p_t =
126 psia" in the Table B-3 heading for that test (p. 54); I treat 125.8 psia (Table B-1) as
the more precise of the two printed figures. No other test ties this pressure within the
2000 R, 18 in group, so no tie-break was needed.

Chosen run: Test 298, l = 18 in, p_t = 125.8 psia, T_t0 = 2007 R (Table B-1, p. 29; Table
B-3 heading, p. 54).

## 1. Nozzle contour

Stated geometry (Introduction, p. 1, and Nomenclature, p. 22-23):
- Throat diameter: 1.803 in (0.04582 m). Throat radius r* = 0.9015 in (0.02290 m).
- Contraction area ratio (approach/throat): 7.75 to 1.
- Expansion area ratio (exit/throat): 2.68 to 1.
- Convergent half-angle: 30 deg.
- Divergent half-angle: 15 deg.
- Exit Mach number: about 2.5.
- Nozzle axial length L = 5.925 in (0.1505 m), nomenclature p. 22.
- Nozzle-inlet radius R = 2.53 in (0.06426 m), nomenclature p. 22.
- Axial origin z = 0 is the nozzle inlet (nomenclature, p. 22: "z, axial distance from
  nozzle inlet").
- Throat axial location: z/L = 0.603 (Section IV, p. 6).

The report does not give a numeric throat radius of curvature, inlet radius of curvature,
or blend radii in the text; it names the symbols r_t (nozzle-throat radius of curvature)
and r* (nozzle-throat radius) in the nomenclature (p. 22) but I could not find a stated
value for r_t anywhere in the OCR text, and the dimensioned nozzle drawing (Fig. 4, p. 3,
and the flow diagram, Fig. 1, p. 1) is a photograph/line drawing, not machine-readable text,
so no numeric value could be extracted from it here.

Table B-2 ("Tabulation of static pressure data," pp. 30-32) tabulates, for 32 fixed pressure
taps, the axial distance ratio z/L and the local-to-throat area ratio A/A*. This is geometric
(same for every test) and is the nozzle contour as tabulated. Wall radius r(z) follows from
r = r* sqrt(A/A*), using r* = 0.9015 in above; I computed it:

  tap  z/L     A/A*     r (in)   z (in)   r (mm)   z (mm)
  1    0.107   7.001    2.386    0.634    60.6     16.1
  2    0.222   4.691    1.953    1.315    49.6     33.4
  3    0.3396  2.939    1.545    2.011    39.2     51.1
  4    0.3401  2.924    1.541    2.014    39.1     51.2
  5    0.395   2.255    1.354    2.340    34.4     59.4
  6    0.455   1.631    1.151    2.696    29.2     68.5
  7    0.460   1.561    1.126    2.726    28.6     69.2
  8    0.506   1.235    1.002    2.998    25.5     76.1
  9    0.529   1.135    0.960    3.134    24.4     79.6
  10   0.554   1.054    0.926    3.283    23.5     83.4
  11   0.587   1.0078   0.905    3.478    23.0     88.3
  12   0.6018  1.0012   0.902    3.566    22.9     90.6   (near throat)
  13   0.642   1.028    0.914    3.804    23.2     96.6
  14   0.691   1.176    0.978    4.094    24.8    104.0
  15   0.706   1.231    1.000    4.183    25.4    106.3
  16   0.722   1.297    1.027    4.279    26.1    108.7
  17   0.741   1.372    1.056    4.390    26.8    111.5
  18   0.753   1.423    1.075    4.462    27.3    113.3
  19   0.771   1.505    1.106    4.568    28.1    116.1
  20   0.789   1.582    1.134    4.675    28.8    118.8
  21   0.802   1.642    1.155    4.751    29.3    120.7
  22   0.819   1.716    1.181    4.852    30.0    123.3
  23   0.838   1.815    1.215    4.965    30.9    126.1
  24   0.852   1.866    1.232    5.049    31.3    128.2
  25   0.869   1.960    1.262    5.150    32.1    130.8
  26   0.887   2.049    1.291    5.255    32.8    133.5
  27   0.900   2.112    1.310    5.333    33.3    135.4
  28   0.921   2.210    1.340    5.457    34.0    138.6
  29   0.932   2.283    1.362    5.522    34.6    140.3
  30   0.949   2.379    1.390    5.623    35.3    142.8
  31   0.965   2.471    1.417    5.718    36.0    145.2
  32   0.985   2.574    1.446    5.836    36.7    148.2

(Source for z/L and A/A*: Table B-2, pp. 30-32. The r and mm columns are my conversion,
not printed in the report.)

## 2. Approach and calming section

The report describes a single "constant-diameter cooled approach section" upstream of the
nozzle inlet; it does not separately name a "calming section." Its diameter is stated once,
at the boundary-layer traverse station: "Boundary-layer traverses were made in the 5.07-in.-
diameter cooled approach section" (Section II, p. 2) -- 5.07 in = 0.1288 m. The length of
this section was varied in four configurations, 0, 6, 12, and 18 in (0, 0.1524, 0.3048, and
0.4572 m), by changing it "in 6-in. lengths from 0 to 18 in." (Introduction, p. 1). The
chosen run, test 298, uses the 18-in (0.4572 m) configuration (Table B-1, p. 29; Table B-3
heading, p. 54).

## 3. Stagnation pressure and stagnation temperature of the run

p_t = 125.8 psia (867,360 Pa, 0.8674 MPa) -- Table B-1, p. 29 (restated as "126 psia" in the
Table B-3 heading for test 298, p. 54).
T_t0 = 2007 R = 1115.0 K = 841.9 degrees C -- Table B-1, p. 29, and Table B-3 heading, p. 54.

## 4. Measured wall temperature along the wall (test 298)

Table B-3 (concluded), p. 54, tabulates, for each axial station z/L, the wall temperature
T_w in degrees F (column unit from the Table B-3 heading block, p. 33) alongside the heat-
flux and derived columns that are excluded here under the blind rule. z is my conversion to
inches/mm using z = (z/L) x 5.925 in; T_w(K) = (T_w(F) + 459.67) x 5/9.

  z/L     z (in)   z (mm)   T_w (deg F)   T_w (K)
  0.133   0.788    20.0     850           727.6
  0.204   1.209    30.7     851           728.2
  0.276   1.635    41.5     872           739.8
  0.336   1.991    50.6     953           784.8
  0.385   2.281    57.9     968           793.2
  0.429   2.542    64.6     1032          828.7
  0.469   2.779    70.6     1116          875.4
  0.512   3.034    77.0     1156          897.6
  0.541   3.206    81.4     (not reported)
  0.573   3.395    86.2     1166          903.2
  0.603   3.573    90.8     1118          876.5   (near throat, z/L = 0.603)
  0.634   3.756    95.4     1111          872.6
  0.664   3.934    99.9     1082          856.5
  0.693   4.106    104.3    1027          825.9
  0.717   4.248    108.0    1027          825.9
  0.750   4.444    112.9    986           803.2
  0.782   4.633    117.7    981           800.4
  0.825   4.888    124.1    932           773.2
  0.864   5.119    130.0    959           788.2   (plug position A)
  0.864   5.119    130.0    940           777.6   (plug position B, different circumferential
                                                    location, p. 3)
  0.905   5.362    136.2    867           737.0
  0.938   5.558    141.2    809           704.8

At z/L = 0.541 the report records no wall temperature for test 298 (shown as "--" in Table
B-3, p. 54); the station itself is listed, consistent with the other tests in the table.

## 5. Inlet boundary-layer data (inputs only)

The Table B-3 heading for test 298 (p. 54) gives, at the boundary-layer measuring station
1.25 in upstream of the nozzle inlet (Section V, p. 7, and Section II, p. 2, state this
station for all tests): momentum thickness theta = 0.052 in (1.321 mm), displacement
thickness delta* = 0.036 in (0.914 mm), and energy thickness phi = 0.040 in (1.016 mm).

I could not find a separate numeric value, for this specific test, of the velocity
boundary-layer thickness delta or the stagnation-temperature boundary-layer thickness
delta_t (both named in the nomenclature, p. 22-23); the other l = 18 in, ~2000 R tests in
the same table (293, 294, 290) report theta, delta*, and phi in the same three-value format,
never delta or delta_t (Table B-3, pp. 52-54).

Velocity- and stagnation-temperature-profile shapes (comparison with a 1/7-power law and
with the "law of the wall," u+ = 5.5 + 2.5 ln y+ for y+ > 30) are discussed and plotted only
for the 18-in approach length, T_t0 = 1500 R group (Figs. 8 and 9, p. 8-9); I found no
profile plot or shape parameter reported for the T_t0 = 2000 R group that test 298 belongs
to, so none is given here for this run specifically.

## 6. Gas

"Compressed air was heated by the internal combustion of methanol and then mixed to obtain
uniformity before it entered the nozzle... the ratio of methanol-to-air weight flow rate was
small enough, even for the highest stagnation temperature, so that the products of combustion
could be treated approximately as air" (Introduction, p. 1). No numeric fuel-air ratio or
combustion-product mole/mass fraction is stated; the treatment is qualitative ("small enough
... treated approximately as air").

Gas properties: "The variation of viscosity, specific heat, and Prandtl number with
temperature for air was obtained from Ref. 12" (Section VI, p. 10), and Reference 12 is
"Keenan, J. H., and J. Kaye, Gas Tables, Wiley and Sons, Inc., New York, 1956" (References,
p. 23). The specific-heat ratio gamma is likewise not a single fixed number: the free-stream
static temperature and Mach number were obtained "from the experimental static-to-stagnation
pressure ratio p/p_t for isentropic flow (gamma = const.) at the corresponding stagnation
temperature value of gamma for the products of combustion" (Section III, p. 4-5, and
Appendix B introduction, p. 28) -- i.e., gamma is evaluated from the same air gas tables
(Ref. 12) at the local stagnation temperature, not quoted as a single number in the text.
The adiabatic-wall recovery factor used is Pr^(1/3) at stagnation temperature (Eq. 2,
Section III, p. 4-5).

## 7. Mass flow rate

m = 3.635 lb/sec (1.6488 kg/s) for test 298 (Table B-3 heading, p. 54).

## 8. Station list

Axial origin: z = 0 at the nozzle inlet (nomenclature, p. 22).

Heat-flux/wall-temperature measurement stations (thermocouple plug axial locations,
Fig. 3 and its table, p. 3, reproduced as the z/L values used throughout Table B-3):
z/L = 0.133, 0.204, 0.276, 0.336, 0.385, 0.429, 0.469, 0.512, 0.541, 0.573, 0.603 (throat),
0.634, 0.664, 0.693, 0.717, 0.750, 0.782, 0.825, 0.864 (two plugs, different circumferential
positions), 0.905, 0.938 -- 21 distinct axial locations, 22 physical plugs total, "located
at each of twenty-one axial locations, except at z/L = 0.864, where there were two" (Section
II, p. 3).

Separately, the static-pressure taps (Table B-2, pp. 30-32) sit at the 32 axial stations
listed in Section 1 above.

## 9. Stated measurement uncertainties (inputs only)

- Stagnation temperature: the two stagnation-temperature thermocouples, located 1 in from
  the centerline and spaced 180 deg apart circumferentially, "generally read within 2% of
  each other" (Section II, p. 2).
- Gas-side wall temperature: "The gas-side wall temperatures determined from the different
  thermocouple combinations in each plug were generally within 1%" (Section III, p. 4).
- Thermocouple-junction location: "a +-1% error due to thermocouple locations determined by
  the Kelvin bridge measurements" (end of Appendix A, p. 28).
- Wall-temperature measurement: "a +-1.5% error from inaccuracies in measurement of the
  temperatures within the nozzle wall" (end of Appendix A, p. 28).
- Stagnation-minus-wall temperature difference: "a +-5% uncertainty in the difference between
  stagnation temperature and gas-side wall temperature" (end of Appendix A, p. 28).

I found no stated numeric uncertainty for the stagnation-pressure measurement (mercury
manometers, Section II, p. 2) anywhere in the text.

## What I could get and could not get

Items 1 through 9 above were all found in some form; the specific gaps, noted inline, are:
no numeric throat or inlet radius of curvature or blend radii (item 1); no numeric value of
the velocity or stagnation-temperature boundary-layer thickness (delta, delta_t) or any
profile shape parameter for this specific run (item 5); no numeric fuel-air ratio or gamma
value as a single number (item 6, both are table/chart lookups in Ref. 12, not quoted
numbers); and no stated stagnation-pressure measurement uncertainty (item 9).
