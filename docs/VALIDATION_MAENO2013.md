# First validation case, on paper: Maeno et al. 2013 magnetic thrust chamber

Status: draft, 30 September 2026. No run. Purpose: list what is known, what must be assumed,
and what an RZ model loses, before any code is pointed at this experiment.

Labels: **reported** (in the paper, with location), **derived** (arithmetic shown from reported
numbers), **inferred** (my reading, could be wrong), **assumed** (a choice we would have to make),
**not found**.

Primary source, read in full (3 pages, open access): A. Maeno, N. Yamamoto, S. Fujioka, Y. Mori,
A. Sunahara, T. Johzaki, H. Nakashima, "Analysis of Laser Wavelength and Energy Dependences of the
Impulse in a Magnetic Thrust Chamber System for a Laser Fusion Rocket", Trans. JSASS 56(3),
170-172 (2013), doi:10.2322/tjsass.56.170. Units were checked against the rendered page, not
only the text layer. The text layer drops the micro sign in two places.

Companion (abstract only, full text paywalled): Maeno et al., Appl. Phys. Lett. 99, 071501 (2011),
doi:10.1063/1.3626600. Field-strength study (not retrieved): Maeno et al., J. Propul. Power 30(1),
54-61 (2014), doi:10.2514/1.B34911.

## 1. The measurement

| Item | Value | Label |
|---|---|---|
| Measured quantity | impulse on a pendulum carrying the magnet (not the target) | reported, p.170-171, Fig. 2 |
| Headline impulses | 1.3, 5.0, 6.4 mN s at 548, 568, 550 J for 1053, 527, 351 nm | reported, p.171 text and Fig. 5 |
| Full range | 0.02 to 10 mN s over 30 to 900 J | reported, abstract |
| Detection limit | 17 uN s | reported, p.171 (micro sign confirmed on the rendered page) |
| Calibration | 240 uN s per volt of LED signal, R^2 = 0.99953, linear fit | reported, Fig. 4 |
| Calibration range | about 0.39 to 0.75 mN s (read from Fig. 4) | inferred (read off the figure) |
| Natural period | 1.3 s | reported |
| Shots per point, error bars | none given | not found |
| Field off | impulse zero without the magnetic field (2011 paper, 0.7 J shots) | reported, APL abstract |

Two consequences:

- The 1.3 to 6.4 mN s values lie above the calibration impacts, by up to 8.5 times (derived:
  6.4 / 0.75). A pendulum is linear at small angles, so this is probably fine. But the paper gives
  no direct evidence, and no error bars. Our measurement band has to be built from the scatter
  between neighbouring shots in Fig. 5. By eye that is roughly a factor of 1.5 (inferred). This
  is a weak reference band, and it should be said so.
- What we must compute is the axial force on the magnet, which is minus the force of the magnet's
  field on the plasma currents. In the split-field form that is exactly bodyAxialForce
  (the coil reaction), integrated over time. The sensitive direction of the pendulum relative to
  the magnet axis is not stated. Assumed: along the magnet axis.

## 2. Inputs that are known

| Item | Value | Label |
|---|---|---|
| Laser | GEKKO XII, single beam, 1.3 ns pulse | reported |
| Wavelengths | 1053, 527, 351 nm | reported |
| Focusing | F/3 lens, spot diameter equal to the target diameter | reported |
| Ablation area used in their model | 3.5e-6 m^2 (target surface over 2 pi sr) | reported, p.172 |
| Peak-average intensity | 550 J / (3.5e-6 m^2 x 1.3e-9 s) = 1.2e17 W/m^2 = 1.2e13 W/cm^2 | derived |
| Target | polyacetal sphere, 1.5 mm diameter, on a carbon fibre from a glass rod | reported |
| Target mass | 1.42 g/cm^3 x (pi/6)(0.15 cm)^3 = 2.5 mg | derived (density is a handbook value, not in the paper) |
| Magnet | cylindrical NdFeB, 50 mm diameter x 40 mm long, (BH)max 0.382 mJ/mm^3 | reported |
| Geometry | target on the magnet axis, 33 mm from the magnet face | reported (Fig. 3 axis, text) |
| Laser direction | 66.5 degrees from the magnet axis | reported |
| Field at target | 0.1 T, calculated, not measured | reported, Fig. 3 |
| Shield | plate "connected to a grounded vacuum chamber", 2 mm from the magnet | reported |

A check on the field (derived): for an ideal linear magnet, (BH)max = Br^2 / (4 mu0), so
Br = sqrt(4 mu0 x 3.82e5) = 1.39 T. The on-axis field of a uniformly magnetised cylinder (radius
R = 25 mm, length 40 mm) at distance z = 33 mm from its face is
(Br/2)[(z+L)/sqrt(R^2+(z+L)^2) - z/sqrt(R^2+z^2)] = 0.695 x (0.946 - 0.797) = 0.103 T. This
reproduces their 0.1 T, so the magnet can be modelled as a fixed surface current
(magnetisation M = Br/mu0 along the axis) with no free parameter.

## 3. What must be assumed

| Item | Why it matters | Proposed handling |
|---|---|---|
| Chamber pressure | background gas slows the plume and carries momentum | not found; assume <= 1e-2 Pa and show background collisions negligible at that level; band it |
| Shield plate material, size, mounting | if mounted on the pendulum, plasma hitting it adds impulse; if metal, it excludes changing field on microsecond scales | inferred: mounted on the chamber (not measured). Model as a perfectly conducting disc fixed in the lab frame; band with and without it |
| Magnet conductivity | sintered NdFeB conducts (about 7e5 S/m, handbook). Skin depth at 1 MHz is sqrt(2/(mu0 sigma omega)) = 1.5 mm, so the magnet excludes the plasma's induced field on the run time | treat magnet surface as a conductor for b1, b0 fixed |
| Plume source | the wavelength dependence comes entirely from laser absorption (critical density 1, 4, 9 e27 m^-3 for 1, 2, 3 omega, reported) and ablation | see the decision below |
| Ionisation state | their model assumes C4+ dominant | reported as their assumption; we need an ionisation model or a band over Z |
| Laser absorption fraction | their model assumes 1 | reported as their assumption; real absorption at 1 omega is lower; band |
| Pendulum axis | not stated | assumed along the magnet axis |

An energy bound that any plume we use must satisfy (derived). If the magnet receives impulse I
from plume mass m moving at v, with I between m v (stopped) and 2 m v (reflected), and kinetic
energy at most the laser energy E, then m >= I^2 / (8 E) and v <= 4 E / I. For 3 omega:
m >= (6.4e-3)^2 / (8 x 550) = 9.3 ug and v <= 4 x 550 / 6.4e-3 = 340 km/s. Only part of the
plume heads toward the magnet, so the real mass is larger. A plume model that cannot carry
about 10 ug or more toward the magnet cannot reproduce 6.4 mN s whatever the field does.

DECISION NEEDED: laser-target physics in or out of the first validation.

- **(a) Prescribed plume.** The plume (mass, velocity distribution, temperature, charge state,
  angular shape) is an input, taken with a band from an independent published radiation-hydro
  scaling, and not from Maeno's own fitted model, whose f = 0.08 coefficient would be an anchor.
  This validates only the magnetic part, plume to impulse on the magnet, which is what the
  instrument is for. Cost: small. Weakness: the wavelength trend, the most distinctive feature
  of the data, is then an input, not a prediction.
- **(b) Laser deposition and ablation in the tool.** This needs inverse-bremsstrahlung
  absorption up to the critical density, a ray trace, radiation transport and a cold-solid EOS.
  Then the 1, 2, 3 omega trend is a prediction. Cost: a large new physics area, months.

My recommendation is (a) first, because it tests the part we have. But the "only the spark and the
fuel are boundary inputs" rule from the chemical work points toward (b), and that is Ben's call.

## 4. What an RZ (axisymmetric) model loses

The magnet, its field and the target position share one axis. That part is exactly
axisymmetric. The laser is not: one beam arriving 66.5 degrees off the axis ablates one side of
the sphere, and the plume leaves roughly back along the beam (inferred: standard for single-sided
ablation; Fig. 2 appears to show the beam arriving from the magnet side). Options in RZ:

1. **Plume along the axis, toward the magnet.** This puts all the plume momentum where the
   field is strongest. It probably overestimates the impulse (inferred), because the real plume
   meets the field at an angle and part of it escapes sideways.
2. **Isotropic plume** (a spherically symmetric ablation of the same total mass and energy).
   About half heads away from the magnet. It probably underestimates the impulse (inferred).
3. **Axisymmetric cone matched to the projection.** Choose an annular-cone plume whose axial
   momentum toward the magnet equals the real plume's projection, cos(66.5 deg) = 0.40 (derived).
   This keeps the axial momentum budget but puts the mass at the wrong azimuth relative to the
   field lines.

Options 1 and 2 bracket the answer. That bracket is itself a result: if the measurement falls
outside it, something in the physics is wrong. If it falls inside, RZ cannot say more, and the
honest output is the bracket.

Also lost or at risk in RZ and in single-fluid MHD (to be flagged by the validity monitor,
`docs/VALIDITY_MONITOR.md`). Rough estimates (inferred, from assumed plume values: C4+,
v ~ 100 km/s, T_e ~ 20 eV and n_e ~ 1e21 m^-3 near the cavity edge):

- Ion gyroradius at 100 km/s in 0.1 T: m v / (q B) = (12 x 1.66e-27 x 1e5) / (4 x 1.6e-19 x 0.1)
  = 3.1 cm, comparable to the 3.3 cm standoff. The ions are not magnetised on the device scale.
  This is outside single-fluid MHD (row 4 of the monitor).
- Ion inertial length: d_i = 2.28e7 cm x sqrt(12) / (4 sqrt(2.5e14)) = 1.25 cm, so d_i / L ~ 0.4
  to 1: the Hall term matters.
- Electron mean free path: 1.44e17 x 20^2 / (1e21 x 10) m = 5.8 mm, Kn ~ 0.2 to 0.6: marginal
  for fluid closure.
- Biermann-battery field from crossed density and temperature gradients in a one-sided plume is
  inherently 3-D and generates toroidal field; RZ with an axial plume suppresses it.
- Flute (Rayleigh-Taylor) modes at the decelerating cavity edge are azimuthal. They are known in
  laser-plasma-in-field experiments and invisible in RZ.

These numbers depend on the plume assumptions, but they say something plain. **This
experiment sits at or beyond the edge of resistive single-fluid MHD.** The first comparison
should be read as a test of whether the validity monitor flags it correctly, as much as a test
of the impulse. A Hall or hybrid model (Inatomi et al. 2023, J. Evolving Space Activities 1, 4,
doi:10.57350/jesa.4, used 3-D hybrid-PIC on a related smaller experiment and still came out
3 to 7 times below the measurement) is the natural next rung.

## 5. What a run would report

Per wavelength: impulse on the magnet (time-integrated bodyAxialForce along the axis) for RZ
options 1 and 2, with the plume band, the mesh band, and the validity-monitor shares; the ratio
to the measurement with the measurement band from section 1; and the energy budget (laser energy
to plume kinetic and thermal, field, losses). Run length: plume transit over 33 mm at
100-300 km/s is 0.1 to 0.3 us (derived), so a few microseconds covers the interaction.

Pass criterion (to be fixed before the run, never after): measured value inside the bracket
of options 1 and 2 with all bands. If outside, report which way and which flags were raised.
