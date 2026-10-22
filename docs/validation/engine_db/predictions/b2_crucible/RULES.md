# Batch b2_crucible: rules (the CRUCIBLE tier)

Written 22 October 2026, about 02:10 CDT. These rules were fixed before any CRUCIBLE run of any engine in the batch. The runs come only after the RL10 pipeline works (the Director, 22 October). The protocol is `../../PREDICTION_PROTOCOL.md`. Batch b1_tier0 (`../b1_tier0/RULES.md`) is frozen at its prediction commit (e747f22), so this tier's rules are a file of their own and not a section added to b1's. A change to these rules after b2's prediction commit (P) makes a new batch. A change before P is allowed under *Amendments* below.

## The batch

- **Specs:** the same `engines.json` as b1, at cfd6f12, SHA-256 `f6cda98344ed5a230442eb8a896049c554010ee01e67d8ace5503f8860b14c47`. Round 3's records, if any, go into a new `engines.json` and so a new batch. They do not change this one.
- **Engines:** b1's 27, in the same groups (S has 12, P has 15; `../b1_tier0/RULES.md`, *The groups*).
- **Tier:** CRUCIBLE's full engine model, the same as for the RL10 (`../../../rl10a33a/tier1_ensemble_plan.md`), not a reduced mode. Ben wants the highest-fidelity model shown (the Director, 22 October).
- **Quantities predicted:** the same as b1's:
  - for all 27: vacuum Isp, sea-level Isp and c*;
  - for the 14 with a printed or derived throat: vacuum thrust, sea-level thrust and the mass flow.

  Each is given on the chamber basis and on the engine basis (rule 12).
- **The sealed step** waits for b2's P commit and its register line, then scores b1 and b2 together (b1's *Tiers and blindness*). If b2 never covers an engine, that engine is still sealed only after b2's P, and the scoring says so.

## What "inferred" means here

The full engine needs inputs that most records lack: the chamber, the contour, the injector, the wall and the injection states. Rules 1 to 9 fill each one the same way for every engine, from standard design correlations, with no choice made per engine. Each fill is one of three kinds:
- **printed:** taken from the record as it stands;
- **derived:** computed from the record's inputs by a rule;
- **inferred:** from a correlation, or guessed.

Every fill is labelled with its kind in the run's setup file and in `predictions.json`. Every inferred value has a bracket. Each end of a bracket is a member, as in the RL10 plan. A source's status is stated with its value: *confirmed* (read in the cited document or a direct quotation of it), *secondary* (read in another work's citation of it), or *memory* (the predictor's recollection of it, not yet read).

## The rules (the same for every engine)

**1. The throat's size.** Throats that are printed or derived keep b1's rule 7: used as printed, or with the lip and envelope factors and their brackets.

For the 13 engines with no throat (fastrac_mc1_15to1, nk33_aj26_58, rd170, aerojet_lox_lch4_rce_altitude_nozzle, sps_aj10_137, lr87_aj_5, le5b3_qual, le7a_nominal, hm7b_rated, ce20_gg_nominal, aestus_eps, h1b, rs27a), the size is filled by cycle class from the batch's own printed throats (derived; no outcome is used):
- **pump-fed** (gas generator, staged combustion, expander): D_t is the geometric mean of the 10 pump-fed throats in `engines.json`, 0.288 m. The bracket is their minimum and maximum, 0.090 and 0.889 m.
- **pressure-fed:** D_t is the geometric mean of the 4 pressure-fed throats, 0.097 m. The bracket is 0.022 and 0.217 m.

These 13 engines are scored on Isp and c* only, as in b1. Their thrust and mass flow are printed for information and not scored.

**2. The chamber: contraction ratio.** A printed CR is used as it stands (rs25_block2a_npl, aerojet_lox_lch4_rce_altitude_nozzle, aj10_190_oms). Otherwise CR = 8.0 D_t^−0.6 + 1.25, with D_t in cm. This fit is attributed to Humble et al.; its status is *secondary* and unconfirmed. The bracket is 1 + (CR − 1)/1.6 to 1 + 1.6 (CR − 1) (guessed). That spans most of the general 2-to-10 range by throat size (*secondary*, Humble et al.). A printed chamber diameter (rs25_block2a_npl) sets CR as derived, with no member.

**3. The chamber: characteristic length L\*.** A printed chamber length, together with CR and the contour, sets the chamber volume as derived, with no L* member (rs25_block2a_npl, aerojet_lox_lch4_rce_altitude_nozzle, aj10_190_oms). Otherwise L* comes from the propellant family. The values are the Huzel & Huang / Sutton table of L*:

| Family | Central | Bracket | Status |
|---|---|---|---|
| LOX/RP-1 and kerosene | 1.14 m (45 in) | 1.02 to 1.27 m (40 to 50 in) | secondary (Huzel & Huang, via the Sun Devil AIAA paper, p. 4) |
| LOX/LH2 (H2 injected as gas, rule 6) | 0.76 m (30 in) | 0.56 to 1.02 m (22 to 40 in): the table's gas-H2 and liquid-H2 rows together | memory |
| LOX/CH4 and LNG | 1.02 m (40 in) | 0.76 to 1.27 m (30 to 50 in): the LOX/RP-1 row widened down, since methane is not in the table | memory, guessed |
| N2O4 or MON-3 with A-50, MMH, UDMH or UH25 | 0.76 m (30 in) | 0.60 to 0.89 m (24 to 35 in) | memory |

The chamber volume is V_c = L* A_t, measured up to the throat. The cylinder's length is what remains after the convergent section (rule 4). If that length comes out below 0.5 D_c at any member, the cylinder is set to 0.5 D_c and the run is flagged.

**4. The convergent section and the throat** (inferred, from NASA SP-8120 and Huzel & Huang):
- **Convergent half-angle:** 30° central, members at 20° and 45°. The 20° to 45° range is from Huzel & Huang (*secondary*).
- **Upstream throat radius:** 1.5 r_t central, members at 1.0 and 2.0 r_t. SP-8120 p. 65 asks for more than 0.6 and calls 1.0 best (*confirmed*, read by the literature agent). The 1.5 value is Rao's convention (*secondary*, Newlands).
- **The chamber-to-cone corner:** a circular arc of radius 1.0 r_t (guessed). It has no member and is a declared limit.
- **Downstream throat radius:** 0.382 r_t (*secondary*, SP-8120 via Newlands). It has no member of its own, since it is a fixed part of each bell in rule 5.

**5. The nozzle: a Rao-type 80% bell** (inferred), built from the throat and the record's ε:
- The downstream arc of 0.382 r_t ends at the angle θn.
- A quadratic Bézier runs from there to the exit at the angle θe.
- The length is L_N = 0.8 [(√ε − 1) r_t + 1.5 r_t (sec 15° − 1)] / tan 15°, which is 80% of a 15° cone.
- θn and θe are read at the record's ε from SP-8120 Fig. 5's 80% curves (*confirmed* as the source). The readout is digitised once, before the first run, and committed with the setup files. It is a design chart, not an outcome.
- **Members:** the 60% and the 90% bells, with their own θn and θe from the same figure.
- The record's `contour_type` (`bell_truncated_ideal`, `bell_thrust_optimised` or `bell_unspecified`) does not change the rule; it is printed beside the result.
- The exit is at the record's ε in every member, so a member changes the length and the wall angles only.

**6. The injection states.**
- **A printed injection temperature** is used as it stands.
- **The phase,** by rule (inferred):
  - H2 is always injected as gas.
  - CH4 and LNG are gas in expander cycles and liquid otherwise.
  - Kerosene, RP-1 and the storables are liquid.
  - In staged-combustion engines the preburner's exhaust is one stream (rule 7).
- **The temperatures,** central and bracket:
  - LOX, 90 K (88 to 95 K);
  - liquid CH4 or LNG, 112 K (110 to 120 K);
  - kerosene and the storables, 293 K (273 to 313 K);
  - gas H2 outside staged combustion, 150 K (60 to 250 K, guessed);
  - gas CH4 in an expander cycle, 250 K (180 to 350 K, guessed).

  Each bracket's two ends are one member on the fuel and one on the oxidiser. Members whose effect is below 0.05% on every quantity at the family's pilot engine (rule 13) are carried as tier 0's injection-temperature member was (*How the band is formed*).

**7. Staged combustion's preburner stream** (rs25_block2a_npl and le7a_nominal fuel-rich; nk33_aj26_58, rd170 and rd_253 oxidiser-rich). The preburner's exhaust enters the main injector as one gas stream. It carries all of one propellant and the part of the other that makes the preburner's temperature:
- fuel-rich: 900 K central, members 800 and 1050 K (b1's rule 4 bracket for T_gg);
- oxidiser-rich: 750 K central, members 650 and 850 K (guessed).

Its composition is CEA's equilibrium at that temperature and at the chamber's pressure. The rest of the other propellant enters as rule 6 says. Turbine work does not change the stream's total enthalpy, since the preburner is adiabatic, so the main chamber's energy is the propellants' own.

**8. The injector element** (inferred). A printed element type is used if CRUCIBLE models that class of element. Otherwise the family's generic element is used, and the record's printed type is shown beside the result as a declared limit.
- **One stream is gas** (LOX/LH2, expander CH4, and every staged-combustion engine): a shear coaxial element. For fuel-rich and LOX/H2 engines, the liquid oxidiser is in the post and the gas in the annulus. For oxidiser-rich staged combustion, the gas is in the post and the liquid fuel in an annulus.
- **Both streams are liquid:** an impinging doublet, like-on-like for LOX/RP-1 and kerosene, and unlike for LOX/liquid CH4 and the storables. The element types by family are from NASA-TM-110632 (Ito 1995, pp. 11 to 13; *confirmed*).
- A printed swirl coaxial (nk33_aj26_58, aestus_eps) is modelled as swirl if CRUCIBLE has swirl, otherwise as shear coax with the limit stated. A printed pintle (lm_descent_engine), triplet (aerojet_lox_lch4_rce_altitude_nozzle) or `other` element is replaced by the family's element, with the limit stated.
- **The layout:**
  - annular rings, as for the RL10 (axisymmetric);
  - a ring pitch of 17 mm for coax (members 12 and 25 mm, guessed) and 12 mm for doublets (members 8 and 20 mm, guessed);
  - a printed element count sets the pitch instead, as p = √(A_face / N), with no member.
- **The pressure drop:** the injector's Δp is 17.5% of p_c (members 15% and 20%; SP-8089 and Huzel & Huang, *secondary*). It sets the injection velocities with the record's mixture ratio.
- **Doublets** have an impingement angle of 60° and L/D 6, with no member (SP-8089, *secondary*).
- **Drop sizes:**
  - for a liquid in a shear coax: the RL10's SPRAY.md correlation (Burick), its bracket being the correlation's stated scatter;
  - for a doublet: the like- or unlike-doublet correlation (Dickerson, or the one SPRAY.md adopts), with its bracket the correlation's stated scatter.

  The doublet correlations are not yet read (a gap). They are read and committed with the setup before the first run of a family that uses them, and that is an amendment (below).

**9. The wall** (inferred; isothermal on the hot-gas side, as the RL10's tier 1):

| Cooling (the record's `chamber_cooling`) | Central | Members | Status |
|---|---|---|---|
| regenerative | 700 K | 400 and 1000 K | memory (hot walls of about 700 to 900 °F), widened |
| ablative | 1900 K | 1500 and 2300 K | memory (a char surface of about 3000 °F) |
| film | 700 K, with a fuel film of 10% of the fuel flow injected as an outer ring at the face | 5% and 20% film; 400 and 1000 K | guessed |
| combined, unknown or missing | 700 K | 400 and 1900 K | guessed: from the regenerative low end to the ablative and radiative range |

The injector face is isothermal at the injected fuel's temperature, with an adiabatic member, as RL10's member 19.

**10. The operating point: hitting the record's p_c.** CRUCIBLE imposes the mass flows and solves for p_c. So each run, the central case and every member, iterates the total mass flow at the record's mixture ratio (by secant) until the chamber pressure at the record's `pressure_station` equals the record's p_c within 0.2%:
- an injector-end or face station is matched to the face's area-mean static pressure;
- a nozzle-stagnation station is matched to the stagnation pressure at the convergent section's entrance;
- an unstated station is matched at the injector end.

The mass flow is then an output, and c* = p_c A_t / ṁ. b1's rule 4 (the mixture-ratio basis) and rule 8 (conflicting p_c, as members) apply unchanged.

**11. What CRUCIBLE resolves.** b1's divergence factor (rule 5), discharge coefficient (rule 6) and Rayleigh loss (rule 3) are not used, since the full engine resolves the contour, the throat and the chamber's flow. b1's sea-level rule 9 is applied to the vacuum run:
- F_SL = F_vac − p_a A_e;
- if the computed wall pressure falls below s p_a (s = 0.35, with members 0.25 and 0.45), the nozzle is cut where the wall pressure equals s p_a, and the engine is flagged as separated.

**12. The engine basis.** CRUCIBLE runs the thrust chamber only. For open cycles, the engine basis is b1's rule 4 exactly (the turbine flow f, T_gg and Isp_tx, with their members), applied to CRUCIBLE's chamber values. A turbine exhaust dumped into the nozzle extension is not modelled in the flow; it enters through b1's Isp_tx fraction only, as a declared limit.

**13. The chemistry, the grid and the run** follow the RL10 tier-1 plan:
- **Chemistry:** the finite-rate table is central if it passes TABLE_B.md's gates on the family's check case (c*, Isp_vac and F_vac within 0.2% of CVODES; defect fallbacks under 1% per region). Otherwise the shifting-equilibrium table is central. The other table is the member whenever both exist.
- **SST:** the Hasan-Pecnik property corrections are central, with "none" as the member.
- **Supply turbulence:** 5% intensity (members 2% and 10%) and a length scale of 0.07 D_h (members 0.035 and 0.14 D_h).
- **Grid:** the reported AMR grid, with its finest level halved as the member.
- **Starts:** warm-started members, with one cold-start check per family.
- **Settling:** the same settle rule.
- **Failures:** a failed member is "not measured", and the band is flagged incomplete.
- **Pilot engines:** each family's first run is its pilot. It is chosen by rule as the family's group-P engine with a printed throat and the fewest inferred inputs, ties broken by alphabetical id.

## Prerequisites (capabilities each family needs before its first run)

A family runs only when the RL10 pipeline works (the Director) and when its own capabilities pass their gates. Until then its engines are "not run" with the cause, and b2's P can be committed with them missing. Those engines are still sealed only after b2's P (b1's rule). The status is from the codebase survey of 22 October.

| Family | Engines | Needs, beyond the RL10 pipeline | Has now |
|---|---|---|---|
| LOX/LH2 | rs25_block2a_npl, j2_sa208_504, j2x_100pct_pl, vulcain2_rated, le5b3_qual, le7a_nominal, hm7b_rated, ce20_gg_nominal (8) | the preburner stream (rule 7) for rs25 and le7a | Tables A and B; the LOX spray; shear coax |
| LOX/CH4, LNG | le8, lm10mira_tca_ue2, aerojet_lox_lch4_rce_altitude_nozzle (3) | CH4 tables; a liquid-CH4 spray and doublets (le8, aerojet); the film ring (aerojet) | none of these |
| LOX/kerosene | f1, fastrac_mc1_15to1, nk33_aj26_58, rd170, h1b, rs27a (6) | kerosene-surrogate tables; a liquid-fuel spray; doublets (f1, fastrac, h1b, rs27a); the oxidiser-rich preburner stream (nk33, rd170) | none of these |
| Storables | sps_aj10_137, lm_ascent_engine, lm_descent_engine, lr87_aj_5, lr91_aj_5, aj10_190_oms, r4d_apollo_rcs, aestus_eps, viking_4b, rd_253 (10) | N2O4 (and MON-3) with hydrazine-family tables (equilibrium; finite rate only if a mechanism passes the gates); two liquid sprays; doublets; the oxidiser-rich preburner stream (rd_253) | none of these |

The RL10 pipeline's own items count for every family: AMR with wall functions, the open-area supply and the isothermal wall. The table mechanisms and the surrogate are design work in inquiry2's lane. Each is chosen and committed with its evidence before its family's first run.

## How the band is formed

The band is formed as in b1 and the RL10 plan:
- one-at-a-time members about the central case, one per bracket end;
- the band is the central value with the summed negative and positive deviations;
- the quadrature sum is printed beside it, not judged.

A member whose effect is below 0.05% on every predicted quantity at the family's pilot engine is run only at the pilot. Its deviation is then carried, as a fraction, to the family's other engines, as RL10's member 6 was carried from tier 0. This is the only economy allowed.

A count, derived: up to about 30 members per engine. That is roughly 600 to 800 runs for the batch, before the carry rule. The budget is planned with the RL10 pipeline's measured run times and is not a rule.

## Amendments before P

Some setup items are mechanical readouts of sources these rules name:
- SP-8120 Fig. 5's θn and θe;
- the doublet drop-size correlations;
- each family's mechanism or surrogate.

Each is committed, dated, before the first run of the family that uses it. Any other change before P is dated, reruns every member already run in the affected family, and is listed in `predictions.json`. Every change is made blind, since no outcome of these engines is open before b2's P.

## Exposures declared before any run

- **The predictor's training knowledge** of famous engines, declared for every batch (`../../ENGINE_DB.md`, *Declarations*). The memory-status values above (L* rows, wall temperatures) are generic textbook ranges. None was set by comparing with an engine.
- **The worker who wrote these rules** read every record's inputs (b1's exposure), and on 22 October a digest of the 27 engines' geometry, cooling and injector fields. It read no outcome pointer's value.
  - Rule 1's class sizes are computed from the printed throats in that digest.
  - The worker also evaluated rule 2's fit at RS-25's printed throat (2.34 against the printed CR of 2.66). The fit was not changed after that, and RS-25 uses its printed CR anyway.
- **The literature agent** (22 October) reported design correlations and their sources:
  - NASA SP-8120, SP-8089 and SP-8087;
  - NASA-TM-110632;
  - Huzel & Huang via secondary citations;
  - Newlands;
  - Humble et al. via a secondary citation;
  - titles only for Halchak's *Materials for Liquid Propulsion Systems* and for Betts & Frederick's throttling study.

  It reported no engine's outcome.
- **The codebase agent** (22 October) reported CRUCIBLE's capabilities (the prerequisites table). That report has no engine values.
- **No CRUCIBLE run of any b2 engine has been made.** The RL10 runs are of another case.
