# Batch b1_tier0: rules

Written 22 October 2026, about 02:00 CDT. These rules were fixed before any engine of the batch was run through CEA. The only CEA calls made before this file were a species check at one made-up point (1000 psia, O/F 2, ε 40), which is no engine's point. The protocol is `../../PREDICTION_PROTOCOL.md`. A change to these rules after the prediction commit (P) makes a new batch.

## The batch

- **Specs:** `engines.json` at cfd6f12 (commit S), SHA-256 `f6cda98344ed5a230442eb8a896049c554010ee01e67d8ace5503f8860b14c47`.
- **Engines:** the 27 that are tier0_isp under `check_specs.py`. The 4 left out are rs68a_rated and rd107_8d74k (no mixture ratio), vinci_a5me_full (no chamber pressure) and cuip_psu_shear_coax_case1 (no expansion ratio).
- **Tier:** tier 0, the ideal rocket with NASA CEA (RocketCEA, with the version recorded in `predictions.json`). It uses no CRUCIBLE code, so it has no CRUCIBLE commit or table hash; the script `tier0_batch.py` is in P.
- **Quantities predicted:**
  - for all 27: vacuum Isp, sea-level Isp and c*;
  - for the 14 that are also tier0_thrust: vacuum thrust, sea-level thrust and the mass flow.

## Tiers and blindness (Ben, 22 October: Q1 (a))

**This batch is the first of two predictions on the same 27 engines.** The second is CRUCIBLE's own prediction (batch b2), with its own RULES.md, written later.

**The sealed step (protocol step 5) waits for both.** It does not start after b1's register line. It starts only after b2's prediction commit and register line too. One sealed file then scores both batches together. Both stay blind, since no outcome of these engines is opened before b2 is committed.

If b2 never covers an engine of b1, that engine's outcomes are still sealed only after b2's commit, and the scoring says so.

## The groups (Ben, 22 October: Q2 (a), with (c) later)

**Group S (secondary, 12 engines):** any of the three operating-point inputs (p_c, the mixture ratio or the expansion ratio) comes only from a secondary source. The secondary sources are Wikipedia, Astronautix (Encyclopedia Astronautica) and Braeunig. This is read from each field's `source` and its citation.
- j2_sa208_504, lm_ascent_engine, lm_descent_engine, lr87_aj_5, lr91_aj_5, hm7b_rated, aestus_eps, viking_4b, rd_253, h1b, rs27a;
- aj10_190_oms, whose mixture ratio alone is Wikipedia's. The handoff's count of 11 had allowed it as primary. The rule here does not, so the group has 12.

**Group P (primary, 15 engines):**
- rs25_block2a_npl, j2x_100pct_pl, vulcain2_rated, f1, fastrac_mc1_15to1, nk33_aj26_58, rd170;
- le8, lm10mira_tca_ue2, aerojet_lox_lch4_rce_altitude_nozzle;
- sps_aj10_137, le5b3_qual, le7a_nominal, ce20_gg_nominal, r4d_apollo_rcs.

j2_sa208_504's expansion ratio is from a site that reproduces Rocketdyne's own scans. That counts as primary, but J-2 is in S anyway, through its p_c and mixture ratio.

The scoring reports P and S apart, and P + S beside them, clearly labelled. Round 3 (a search for primary sources for group S) may later give new records. Those are a new `engines.json` and so a new batch. This batch is scored as made.

## The method

**The ideal rocket.** It has an infinite-area combustor, adiabatic and in equilibrium, at a chamber pressure p0 and the chamber mixture ratio, then a one-dimensional expansion to the record's ε. The vacuum Isp is

Isp_vac = λ (Isp_CEA − P) + P,  where P = p_e ε c*_CEA / (p0 g0)

as in RL10's `tier0.py`, with λ the divergence factor on the momentum term only.

**Thrust and mass flow** (the 14 engines with a throat):
- ṁ_ch = Cd p0 A_t / c*_CEA;
- F_vac,ch = Isp_vac ṁ_ch g0;
- F_SL,ch = F_vac,ch − p_a ε A_t, with p_a = 101.325 kPa (unless separated: see rule 9).

**c\*** (scored as the source's p_c A_t / ṁ): c*_eff = c*_CEA (p_c / p0) / Cd.

**The engine basis** (open cycles, rule 4):
- ṁ_eng = ṁ_ch / (1 − f);
- F_eng = F_ch + f ṁ_eng g0 Isp_tx;
- Isp_eng = F_eng / (ṁ_eng g0).

For closed and pressure-fed cycles, f = 0.

**The band.** The central value, then one-at-a-time members, one per rule end. The band is the central value minus the sum of the members' negative deviations, and plus the sum of their positive ones, as in RL10's tier 0. The quadrature sum is reported beside it and not used.

## The rules (the same for every engine; no per-engine choice)

**1. Chemistry.**
- **Central:** shifting equilibrium to the exit.
- **Member:** frozen at the throat (RocketCEA `frozen=1, frozenAtThroat=1`).

**2. Propellant cards.** Each is RocketCEA's card at CEA's reference state:
- LOX at 90.18 K, LH2 at 20.27 K, CH4(L) at 111.66 K;
- the storables at 298.15 K.

The propellants map to these cards:
- O2(L): LOX. H2: LH2.
- CH4, and LNG (le8, lm10mira): CH4(L). LNG's heavier fractions are not modelled; this is a declared limit.
- RP-1: RP1.
- "kerosene", with the grade not stated (nk33_aj26_58, rd170): RP1 central, with Jet-A(L) as a member.
- N2O4: N2O4. Aerozine-50: A50 (50/50 UDMH/N2H4 by mass). MMH: MMH. UDMH: UDMH.
- UH25 (viking_4b): 75% UDMH and 25% N2H4(L) by mass.
- MON-3 (aj10_190_oms): RocketCEA's MON3 card, which is pure N2O4, is central. The member is 97% N2O4 and 3% NO by mass, with NO at its gas-phase enthalpy of formation, so the member overstates the NO's energy slightly.

**Temperatures.** A record's injection temperature is not used at tier 0.
- For regeneratively heated propellants, the injection state includes wall heat that the adiabatic energy balance already counts from the engine-inlet state.
- CEA's liquid storables are defined at 298.15 K only.
- The effect is estimated at no more than about 0.2% of Isp. The largest case is aerojet_lox_lch4's LOX at 113 K, 23 K above its boiling point. It is declared, not bracketed.

**3. The pressure station.** The contraction ratio enters only here. It is the record's printed CR where there is one (rs25_block2a_npl 2.66, aerojet_lox_lch4 5.85), else 3, bracketed 2 to 5.

The loss is the Rayleigh one, in closed form. A constant-area chamber is heated from rest to the Mach number M_c, where the isentropic area ratio A/A* at γ equals CR. Then

p0 / p_inj = (1 + (γ − 1) M_c² / 2)^(γ/(γ − 1)) / (1 + γ M_c²)

with γ CEA's chamber value.

The rule depends on the record's station:
- **nozzle_stagnation** (f1): p0 = p_c. No member.
- **injector_face** (aerojet_lox_lch4, rs27a): central p0 = p_c × the ratio at the central CR. Where the CR is not printed, the members are CR 2 and CR 5.
- **unspecified or missing** (every other engine): central p0 = p_c, reading the printed value as the nozzle stagnation pressure. The member reads it as the injector-face pressure, at the record's printed CR where there is one and otherwise at CR 2, the bracket's largest loss.

**4. Open cycles** (gas_generator; expander_bleed, le5b3_qual). Staged-combustion, expander_closed and pressure-fed engines burn all their flow in the chamber, so these terms vanish there.

- **The turbine flow fraction:** f = ṁ_turbine / ṁ_engine. It is 2.5% central, bracketed 1% to 4%. This is assumed: a textbook range for gas-generator cycles, not taken from any engine of the batch.
- **The gas generator's mixture ratio:** the fuel-rich ratio at which CEA's equilibrium temperature at p_c is 900 K, bracketed 800 to 1050 K.
  - It is solved by bisection between O/F 0.02 and the ratio at which the temperature first reaches 1500 K.
  - For expander_bleed, the turbine gas is pure fuel (O/F 0).
- **The mixture-ratio basis:**
  - **A basis printed as `engine`:** the printed ratio is the engine's. The chamber ratio follows from f and the gas generator's ratio. No member.
  - **A basis printed as `thrust_chamber`:** the printed ratio is the chamber's. No member.
  - **Any other basis** (inferred, assumed, unspecified or missing): the engine's central, with the chamber's as the member.
- **Isp_tx, the turbine exhaust's specific impulse:** a fraction of the chamber's vacuum Isp. There is no turbine model at tier 0, so this is assumed.
  - dumped_overboard: 0.35, bracketed 0.2 to 0.5;
  - into_nozzle_extension: 0.6, bracketed 0.4 to 0.8;
  - unknown or missing: 0.5, bracketed 0.2 to 0.8.

  At sea level the exhaust's own pressure term is in this fraction, and only the chamber's exit area takes −p_a A_e.

**5. The divergence factor.** No record prints its contour's parameters, so every engine takes the same rule, whatever its `contour_type`: λ = (1 + cos θ) / 2 with θ 10° central, bracketed 5° to 15° (15° being the classic cone). That gives λ = 0.99240 (0.99810 to 0.98296).

**6. The discharge coefficient** (ṁ, thrust and c*). Kliegel–Levine, at CEA's chamber γ, with the throat's radius of curvature R/r_t 1 central. The members are R/r_t 0.5 and Cd = 1. The Cd = 1 member also covers sources that report c* with the throat's Cd taken out.

**7. The throat** (the 14 thrust engines). The rule follows how the record's throat was got:
- **Printed, or derived from a printed throat area:** used as it stands, no member. These are j2_sa208_504, le8, lm10mira_tca_ue2, r4d_apollo_rcs, rs25_block2a_npl and f1.
- **Derived from an exit diameter with the lip unstated:** the inner exit diameter is 0.985 D central, bracketed 0.97 to 1.00 D. These are vulcain2_rated, lm_ascent_engine and lm_descent_engine.
- **Derived from an envelope, overall or maximum diameter:** the inner exit diameter is 0.925 D central, bracketed 0.85 to 1.00 D. These are j2x_100pct_pl ("overall engine size"), lr91_aj_5 ("overall"), aj10_190_oms ("Engine Dia.", read as the exit), viking_4b ("maximum diameter") and rd_253 ("overall").

The throat scales with the inner exit diameter, so A_t scales with the factor squared.

**8. A conflicting printed input.** Where a record's flag gives a second printed value for an operating-point input, that value is a member and the record's value stays central:
- j2_sa208_504: 700 psi (4.826 MPa), from NASA SP-8120;
- h1b: 700 psi (4.826 MPa);
- hm7b_rated: 3.5 MPa.

Where a flag reports a conflict but gives no value (rs27a's second p_c), the members are p_c ± 5%.

vulcain2_rated's flag compares its mixture ratio with "commonly cited" figures but gives no source. So it is not a printed alternative, and rule 4's basis members apply as for any engine.

**9. Sea level and separation.**
- The sea-level Isp and thrust are predicted for every engine at p_a = 101.325 kPa.
- If the full-flowing exit pressure is below the separation pressure p_sep = s p_a, with s = 0.35 central and bracketed 0.25 to 0.45 (which spans Summerfield's 0.4 and Schmucker's value at exit Mach 3 to 4), the nozzle is taken as cut where the flow's pressure equals p_sep. The station is found by bisection on ε with the same chemistry. Downstream of it the wall carries p_a.
- The sea-level values are then those of the cut nozzle, and the report flags the engine as separated.
- For an engine that is not separated at any end of the bracket, s has no member.

**10. Which quantity is scored.**
- Isp and thrust are given on both the chamber basis and the engine basis. For closed cycles the two are the same.
- The sealed step records the basis its source states. The scorer compares like with like, and an unstated basis is scored against the engine basis.
- A rated or nominal value is scored apart from a test or flight value (protocol step 6). A model value (for example a c* efficiency without heat loss) is not scored as validation.

## Exposures declared before P

- **The predictor's training knowledge** of famous engines' ratings, declared once for every batch (`../../ENGINE_DB.md`, *Declarations*). No such value is an input, a check or a reason for any rule above. The rules' numbers are generic textbook ranges, and none was set by comparing with an engine.
- **The worker who wrote these rules** read every record's inputs, sources and flags: the screened `engines.json`, a digest of the operating point, the geometry, the temperatures and the throats' derivations. It read no outcome pointer's value; the pointers carry locators and labels only. It also read no source document.
- **vulcain2_rated's flag** cites "~6.1-6.7" as commonly cited mixture ratios. That is an input's value, not an outcome's, and it is not used (rule 8).
- **No CEA run of any engine's point was made before this file was committed.**
