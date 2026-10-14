# RL10A-3-3A: prediction (tier 0)

Status: tier 0 committed 13 October 2026, 22:15 CDT (the commit's author time, which a rebase onto main keeps; its hash, first 92e98c7, changes with each rebase until main is fast-forwarded; this line said "about 22:20" until 22:34), before any measured performance of this engine was read. Protocol: `RL10A33A_PROTOCOL.md` as amended up to and including its 13 October clarification. Tier 1 (the engine, the reported prediction) follows in a later commit. **No performance row is read, by me or an agent, until tier 1 is committed**, so tier 0's numbers here cannot inform tier 1 through a comparison. The comparison (`RL10A33A_COMPARISON.md`) will quote this commit's hash and tier 1's.

Script and record: `rl10a33a/tier0.py`, output in `rl10a33a/tier0_2026-10-13.txt` (RocketCEA 1.2.3, NASA Glenn CEA2; SciPy's PCHIP).

## Method (protocol section 3, tier 0)

- Wall: the radius reading (r_t 2.47 in, r_c 5.13 in), Table E1's nodes plus the throat, face (-13.1 in) and exit (46.0 in) nodes, by PCHIP.
- Cards: `fuel H2 H 2 wt%=100. t(k)=199.667` and `oxid LOXP O 2 wt%=100. h,kj/mol=-12.3713 t(k)=101.167` (CEA's O2(L) at 90.17 K plus the NIST saturated-liquid difference to 101.17 K). 199.67 K is 0.33 K below the lower limit of CEA's H2 polynomial (200 K); CEA extrapolated with no warning, and its printed reactant enthalpy, -665.49 cal/mol (-2.784 kJ/mol), agrees with cp times (T - 298.15 K).
- p_c from mdot = Cd p_c A_t / c*_CEA(p_c) by fixed point (no Rayleigh loss, so the injector-face pressure and the nozzle stagnation pressure both equal p_c). c* = p_c A_t / mdot = c*_CEA / Cd.
- Cd: Kliegel-Levine at R/r_t 0.5, the series' edge (the PCHIP throat's 0.137 is sharper; its extrapolated Cd is a deviation).
- Isp_vac = lambda (Isp_CEA - P) + P, P = p_e eps c*_CEA / (p_c g0), shifting and frozen; lambda central 0.98095, the mean of 0.99232 (exit wall angle 10.054 deg) and 0.96959 (cone half-angle 20.086 deg). F_vac = Isp_vac g0 mdot on the chamber flow.
- Band: central minus and plus the sums of the one-at-a-time deviations on each side (protocol section 4); quadrature printed beside, not judged.

Tier 0 has no boundary-layer, kinetic, mixing or atomisation loss. Shifting and frozen are its two chemistry limits; neither is the reported prediction.

## P1, Table 6.1.1 (O/F 5.055, chamber flow 16.7126 kg/s, engine-inlet flow 16.7471 kg/s)

| Quantity | Central | Judged band | Quadrature (not judged) |
|---|---|---|---|
| p_c = injector-face = nozzle stagnation pressure [MPa] | 3.3009 | 3.3009 to 3.3897 | 0.089 |
| c* = p_c A_t / mdot [m/s] | 2442.3 | 2442.3 to 2508.0 | 65.7 |
| Isp_vac shifting, chamber flow [s] | 461.56 | 455.07 to 467.26 | 7.43 |
| Isp_vac shifting, engine-inlet flow [s] | 460.61 | 454.13 to 466.30 | 7.41 |
| F_vac shifting [kN] | 75.647 | 74.584 to 76.582 | 1.22 |
| Isp_vac frozen, chamber flow [s] | 441.21 | 434.98 to 446.82 | 7.13 |
| Isp_vac frozen, engine-inlet flow [s] | 440.30 | 434.09 to 445.90 | 7.12 |
| F_vac frozen [kN] | 72.311 | 71.291 to 73.232 | 1.17 |

CEA at the central p_c: c*_CEA 2391.9 m/s, Tc 3301.6 K, chamber gamma 1.1473, 1-D Isp_vac 470.20 s shifting and 449.50 s frozen. Cd 0.97937.

## P2, Table 2.5.1 (O/F 5.26, 16.9462 kg/s, all burnt)

| Quantity | Central | Judged band | Quadrature (not judged) |
|---|---|---|---|
| p_c = injector-face = nozzle stagnation pressure [MPa] | 3.3244 | 3.3244 to 3.4138 | 0.089 |
| c* = p_c A_t / mdot [m/s] | 2425.8 | 2425.8 to 2491.0 | 65.2 |
| Isp_vac shifting [s] | 460.87 | 454.41 to 466.57 | 7.41 |
| F_vac shifting [kN] | 76.590 | 75.516 to 77.537 | 1.23 |
| Isp_vac frozen [s] | 438.57 | 432.39 to 444.16 | 7.09 |
| F_vac frozen [kN] | 72.884 | 71.856 to 73.813 | 1.18 |

All flow is burnt at P2, so its chamber-flow and engine-inlet Isp are equal. P2 takes P1's plenum temperatures (the inputs give none for P2; declared in the protocol's clarification). CEA at the central p_c: c*_CEA 2375.9 m/s, Tc 3341.4 K, chamber gamma 1.1438, 1-D Isp_vac 469.49 s shifting and 446.81 s frozen.

## What sets the bands (P1; P2 is within 0.04 s and 0.01 kN of the same)

| Bracket | Isp_vac shifting [s] | F_vac shifting [kN] | p_c [MPa] |
|---|---|---|---|
| lambda at the exit wall angle / at the cone | +5.16 / -5.16 | +0.845 / -0.845 | 0 |
| exit station 46.40 / 45.29 in | +0.53 / -1.33 | +0.087 / -0.217 | 0 |
| Cd at the PCHIP R/r_t 0.137 (0.954, extrapolated) | +0.02 | +0.003 | +0.089 |
| plenum temperatures plus or minus 0.1 R | under 0.002 each | under 0.0004 each | under 1.1e-5 each |
| injector face -13.4 / -12.8 in; the sharp throat | 0 (tier 0 does not see them) | 0 | 0 |

The divergence factor dominates Isp and thrust; the throat's curvature dominates p_c and c*. The wall temperature, drop size and supply turbulence brackets act only in tier 1.

## For the record, unjudged: the diameter reading (r_t 1.235 in, r_c 2.565 in)

Section 1 of the protocol resolved the size question for the radius reading; this is printed as it promised, and is not judged. Cd at R/r_t 0.5, lambda central at its own wall.

| Quantity | P1 | P2 |
|---|---|---|
| p_c [MPa] | 13.313 | 13.421 |
| c* [m/s] | 2462.5 | 2448.2 |
| Isp_vac shifting, chamber flow [s] | 468.74 | 468.19 |
| F_vac shifting [kN] | 76.823 | 77.806 |
| Isp_vac frozen, chamber flow [s] | 454.93 | 452.76 |
| F_vac frozen [kN] | 74.560 | 75.242 |

## Carried notes

- The two exposures declared in protocol section 1 (8 October: an archived c* comparison and the 2.47 in dimension; 13 October: the fuel plenum's pressure through its temperature and density). Tier 0 has no free parameter: p_c follows from mdot, Cd and CEA, and Isp and thrust from CEA and the wall's angles.
- The first run of `tier0.py` named the H2 card "GH2P", which CEA rejected (it printed only zeros; the script now stops on a zero return). The protocol's card is CEA's H2, and that is what produced every number here.
