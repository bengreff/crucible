# Lafon et al. (2008): the 10 atm LOX drop, digitised

For SPRAY.md criterion S2(c). Digitised on 13 October 2026 by a Sonnet research agent, after the model's numbers were committed (fe86952).

**Source:** P. Lafon, H. Meng, V. Yang and M. Habiballah, "Pressure-coupled vaporization and combustion responses of liquid-oxygen (LOX) droplets in supercritical hydrogen environments", *Combust. Sci. Technol.* 180 (2008) 1 to 26. The PDF is on the author's institutional page: https://www.yang.gatech.edu/publications/Journal/CST%20(2008,%20LOX%20droplet,%20Lafon).pdf

**Case:** D0 100 μm, T0 90 K, pure H2 at 1,000 K and 10 atm, quiescent.

**Figures used:**
- Fig. 8a, (D/D0)² against time at several pressures, with the surface taken as the critical-mixing-temperature isotherm. At 10 atm the drop never reaches critical mixing, so 8a and 8b coincide.
- Fig. 4, surface temperature against time.
- Fig. 9a, lifetime against reduced pressure (log-log). Used only as a cross-check.

**Method (the agent's):**
- The raw figures were extracted with `pdfimages`, at their native resolution (Fig. 4 is 1040 × 763 px, Fig. 8 is 900 × 1488 px).
- The pixels were calibrated from the plot frame and the axis end values.
- The 10 atm curve was isolated by connected-component labelling.
- Reading uncertainty is about ±0.01 ms and ±0.005 in (D/D0)², and ±0.3 K in T_s.

## (D/D0)² against t (Fig. 8a, 10 atm; digitised)

| t (ms) | (D/D0)² | t (ms) | (D/D0)² |
|---|---|---|---|
| 0.026 | 0.999 | 4.079 | 0.542 |
| 0.424 | 0.959 | 4.477 | 0.494 |
| 0.834 | 0.918 | 4.887 | 0.442 |
| 1.232 | 0.873 | 5.298 | 0.394 |
| 1.642 | 0.827 | 5.695 | 0.345 |
| 2.053 | 0.781 | 6.106 | 0.294 |
| 2.450 | 0.735 | 6.517 | 0.242 |
| 2.861 | 0.686 | 6.914 | 0.194 |
| 3.272 | 0.637 | 7.325 | 0.142 |
| 3.669 | 0.591 | 7.735 | 0.091 |

The plotted line ends at 7.735 ms, at (D/D0)² 0.091. There is no swelling: the curve falls from the first point.

## Surface temperature against t (Fig. 4, 10 atm; digitised)

| t (ms) | T_s (K) | t (ms) | T_s (K) |
|---|---|---|---|
| 0.023 | 109.6 | 4.09 | 112.0 |
| 0.47 | 111.0 | 4.54 | 112.0 |
| 0.93 | 111.4 | 4.99 | 112.1 |
| 1.38 | 111.5 | 5.46 | 112.2 |
| 1.83 | 111.7 | 5.91 | 112.2 |
| 2.28 | 111.8 | 6.36 | 112.2 |
| 2.74 | 111.8 | 6.81 | 112.3 |
| 3.19 | 111.9 | 7.26 | 112.3 |
| 3.64 | 111.9 | 7.72 | 112.3 |

The rise from 90 K to about 110 K takes less than one pixel column (under 0.02 ms). The paper calls the plateau a pseudo-wet-bulb state, "slightly lower than the boiling temperature of oxygen because of the dissolution of hydrogen into the droplet".

## Their model (quoted or paraphrased from the paper)

- **Thermodynamics.** Real fluid: "a self-consistent and efficient method ... which properly accounts for dense-fluid corrections" (Meng et al. 2005).
- **Transport.** Ely and Hanley's extended corresponding states for viscosity and conductivity, and Takahashi's corresponding-state scheme for binary diffusivity.
- **The liquid.** The same conservation equations are solved inside the drop on a spherically symmetric finite-volume grid that moves with the surface. The liquid's conductivity is therefore finite, not infinite or lumped.
- **Hydrogen dissolution.** Included, through interfacial jump conditions.
- **Grid and step.** Not stated.

## Derived from the tables (by us)

- **Slope.** A least-squares line through the 10 points with (D/D0)² between 0.6 and 0.1 (3.669 to 7.325 ms) has slope −0.1228 per ms. That gives −d(d²)/dt = **1.228 mm²/s**.
- **Lifetime.** The same line reaches zero at **8.49 ms**; the last six points give 8.47 ms. Their Fig. 9a reads about 8.3 ms.
