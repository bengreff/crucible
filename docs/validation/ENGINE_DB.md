# The blind engine database

Status: started 21 October 2026, about 03:10 CDT, independent of inquiry2's RL10 lane. It serves the inquiry trust bar (`~/director/goals.md`): blind predictions of Isp and thrust on 20 or more real engines across all liquid propellants, committed before their measured values are looked up, with the error spread reported. It is TECHNICAL_PLAN step 12's "engines as data", with the outcomes split off.

## Files (`docs/validation/engine_db/`)

| File | What it is |
|---|---|
| `engine_spec.schema.json` | The input schema: one record per engine, a provenance on every field. It has no outcome field, and every object is closed, so an outcome cannot enter under a new key. |
| `engines.json` | The screened specs (*Compilation*, below). |
| `check_specs.py` | Schema and consistency check, the outcome screen (count-only, masked view, redact unread) and the usability table. Stdlib only. |
| `compile_rules.md` | The rules the compiling agents were given, verbatim. |
| `compile_brief_2026-10-21.md` | The rerun's common brief (clock, save per engine), verbatim. |
| `compile_brief_2026-10-21_round2.md` | Round 2's brief (the same, with the night's clock), verbatim. |
| `compile_brief_2026-10-22_round3.md` | Round 3's brief (primary sources for group S's 12 records), verbatim. |
| `sealed_outcome.schema.json` | The outcome schema. Its files (`sealed/<batch>.json`) are written only by the sealed step, after the batch's prediction commit. None exists yet. |
| `PREDICTION_PROTOCOL.md` | How a blind prediction is committed and later scored. |
| `predictions/<batch>/` | A batch's RULES.md, its script, its run log and `predictions.json`. |
| `REGISTER.md` | One line per batch: P's full hash and `predictions.json`'s SHA-256 (protocol step 4). |

## Inputs and outcomes

**Inputs** (the spec): the propellant pair and grades; the cycle and where its turbine exhaust goes; one operating point, as chamber pressure (and its station) and mixture ratio (engine or chamber basis); the throat diameter, contraction and expansion ratios, chamber and nozzle lengths or L*, the contour type and its parameters, the cooling and any film fraction; the injector type and element count; the propellant injection and inlet temperatures.

**Outcomes** (sealed only): vacuum and sea-level Isp and thrust, c*, and the mass flow. Two of these need a reason:
- **Mass flow is an outcome** because the inputs already hold p_c and A_t, and c* = p_c A_t / ṁ. A spec carrying ṁ would hand the predictor the measured c*. So the operating point is set by p_c and the mixture ratio, and ṁ is predicted (derived). RL10A-3-3A's protocol took the opposite choice (ṁ in, p_c solved); either is blind, but not both in one record.
- **The exit pressure, the thrust coefficient, the exhaust velocity and every efficiency** are outcomes or one step from one, so none is a field.

**No input may come from an outcome.** The throat must be printed, drawn, or derived from the exit diameter and the expansion ratio, never from thrust / (C_F p_c) or from ṁ c* / p_c. Where that is the only route, the throat is `missing` and the engine can give Isp but not thrust. The mixture ratio may be the ratio of printed propellant flows, with the flows not recorded.

**Provenance** on every field: *printed* (source and locator, and the value as printed), *digitised*, *derived* (from other inputs of the record, formula given), *inferred* (from analogues, said how), *assumed*, *missing* (value null, reason given). An operating-point value also carries its basis: test, design, rated, nominal or unknown. Every inferred, assumed or missing required field also has an entry in the record's `flags`; nothing is guessed silently.

**What each engine can be used for** (`check_specs.py`'s table, computed from the provenance):
- **tier0_isp:** the propellants, p_c, the mixture ratio and the expansion ratio. The ideal-rocket Isp needs no throat.
- **tier0_thrust:** those and the throat.
- **tier1** (the time-evolving engine): those, plus the contraction ratio, a chamber length or L*, a known contour type, the injector type and both injection temperatures.

## Compilation (21 October 2026)

Four Sonnet agents, one per family (LOX/LH2, LOX/kerosene, LOX/methane, storables with at most one other), worked from `compile_rules.md` and `compile_brief_2026-10-21.md`. Each wrote a staging file outside the repo, checked it against the schema and its own masked view, and reported back with no numbers. The RL10 family is excluded: RL10A-3-3A awaits its own prediction. RL10 variants may join after that prediction's commit.

**The screen, before any record was read** (as for RL10, `rl10a33a/injector_element_inputs.md`):
1. Count-only greps: `check_specs.py --counts` prints the hits per outcome keyword, with no text, over every free-text string (notes, the printed originals, derivations, locators, citations, labels, names, flag reasons, URLs).
2. A masked view of every hit: `check_specs.py --masked` prints its path and its text with each run of digits shown as `#`.
3. A hit that carries a value is removed unread (`--redact PATH`; a number field becomes `missing`), and logged below by path only.
4. Only then is the file read in full and merged into `engines.json`.

**The first run was lost.** It was killed at 03:19 when its turn ended (worker rule 18). It wrote no spec file, only raw downloads in `/tmp/engine_db_staging/`. The rerun (03:22 to 03:42) used the same rules plus `compile_brief_2026-10-21.md`, which adds a clock, a save after every engine, and a ban on power levels quoted as a percentage of rated thrust in ids. The four agents ran in the foreground, and the Director's worker waited for them inside the turn.

**The screen's results** (`docs/evidence/engine_db_screen_2026-10-21.txt`):
- 59 hits over the four files.
- Every hit is one of: a row label or locator ('Specific impulse' row), a citation (an author's initial "N." hits the force-unit pattern), a test duration, a section number ("Sec #"), or an input (O/F, area ratio, a length).
- **No hit carries an outcome value, so nothing was redacted.**
- The merged file's counts equal the sum of the four files' counts, keyword by keyword.

The records were then read in full (a compact digest of every field's value and provenance, every source and every flag) and merged unchanged.

**The engines** (18, `engines.json` SHA-256 `c73b97f5…1389`):

| Family | N | Engines |
|---|---|---|
| LOX/LH2 | 5 | rs25_block2a_npl, j2_sa208_504, j2x_100pct_pl, rs68a_rated, vulcain2_rated |
| LOX/kerosene | 4 | f1, fastrac_mc1_15to1, nk33_aj26_58, rd170 |
| LOX/methane | 4 | le8, lm10mira_tca_ue2, aerojet_lox_lch4_rce_altitude_nozzle, cuip_psu_shear_coax_case1 |
| Storable | 5 | sps_aj10_137, lm_ascent_engine, lm_descent_engine, lr87_aj_5, lr91_aj_5 |

**Usability** (check_specs.py):
- tier0_isp: 16. The two without it are rs68a_rated (mixture ratio missing; the baseline RS-68's value was deliberately not borrowed) and cuip_psu_shear_coax_case1 (expansion ratio missing).
- tier0_thrust: 10.
- tier1: none.
- The CUIP rig's outcome is c* efficiency, so it can serve a c* prediction only. Its pointer is the source's own estimate *without heat loss*, which makes it a modelled value, scored apart from test values.

**Geometry not sourced from a nozzle dimension** (flagged in each record):
- **The throat is missing, so these give Isp only (6):** fastrac_mc1_15to1, nk33_aj26_58, rd170 (only the four-chamber envelope is printed), aerojet_lox_lch4_rce_altitude_nozzle (only derivable from an outcome, so not derived), sps_aj10_137, lr87_aj_5 (twin chambers, so the printed diameter is ambiguous).
- **The throat is derived from a diameter that may be the outer lip or the engine envelope (7):**
  - j2_sa208_504 (lip unstated);
  - j2x_100pct_pl (an "overall engine size");
  - rs68a_rated (an infobox diameter that is not variant-specific; secondary);
  - vulcain2_rated (lip unstated);
  - lm_ascent_engine and lm_descent_engine (secondary exit diameters, lip unstated);
  - lr91_aj_5 (an overall engine diameter; secondary only).
- **The throat is from a printed nozzle or chamber dimension (5):**
  - rs25_block2a_npl (the chamber diameter and the contraction ratio);
  - f1 (a printed throat area);
  - le8 (printed);
  - lm10mira_tca_ue2 (printed under a mislabelled heading, read as the throat);
  - cuip_psu_shear_coax_case1 (the chamber diameter and the contraction ratio).
- **The chamber is sourced only for rs25_block2a_npl, aerojet_lox_lch4_rce_altitude_nozzle and cuip_psu_shear_coax_case1.** Only these three print both the contraction ratio and a chamber length. Every other record lacks them.
- **No record prints its nozzle contour.** It is inferred as a bell, or missing.

**Points a batch's rules must handle** (the records flag each of them):
- Only secondary sources give these values:
  - J-2's p_c and mixture ratio;
  - RS-68A's p_c;
  - the LM engines' p_c;
  - the LMDE's mixture ratio;
  - everything in lr87_aj_5 and lr91_aj_5 (Astronautix).
- Vulcain 2's printed mixture ratio is higher than the commonly cited ones. Its basis is missing, and it may be the chamber's.
- The pressure station is unspecified for most records.
- The mixture ratios of the open cycles (J-2, J-2X, LE-8, LR87, LR91, Fastrac) may include the gas generator's flow.
- Turbine exhaust goes into the nozzle extension on J-2, J-2X, Vulcain 2, F-1 (inferred) and LE-8.
- le8 is the as-tested low-expansion nozzle, not the flight extension.
- lm10mira_tca_ue2 is a water-cooled ground article.
- nk33_aj26_58 is Aerojet's modification at Kistler inlet conditions.
- cuip burns gaseous methane.

**What was dropped, for lack of time or of primary sources this session:**
- LOX/LH2: HM7B, LE-7A, LE-5B, Vinci and CE-7.5/CE-20.
- LOX/kerosene: H-1 (Huzel & Huang's "A-1" is fictitious), RS-27A, LR79/MA-x, RD-107/108, Merlin and Rutherford.
- LOX/methane: Morpheus and RS-18 (no printed geometry with p_c and mixture ratio); Raptor, Prometheus, BE-4 and TQ-12 were not attempted.
- Storables: OMS AJ10-190, Viking, Aestus, RD-253/275 and R-4D.

No RL10 source was opened (the agents' reports).

**Exposure noted by the agents:** several sources print performance beside the geometry the records cite. The agents named Vulcain 2's AIAA 2006 Table 1, the F-1 manual's Fig. 3-1 and an AEDC lecture's SPS sentence. As declared below, only the agents saw those values.

## Round 2 (21 October 2026, 21:47 to 22:33 CDT)

Four Sonnet agents, in the foreground and two at a time, worked from `compile_rules.md` and `compile_brief_2026-10-21_round2.md` (round 1's brief with the night's clock):
- LOX/LH2, then storables, from 21:47 to 22:14.
- LOX/kerosene, then the throat agent, from 22:15 to 22:33.

The family agents took the candidates round 1 dropped. The throat agent looked for primary sources for the 13 weak or missing throats. The screen was round 1's: counts, then the masked view, then redaction unread, then the full read (`docs/evidence/engine_db_screen_round2_2026-10-21.txt`).

**One value was redacted unread:** `r2_storable.json` engines[3].operating_point.mixture_ratio_basis.note (viking_4b). Its masked view showed a mass flow with its unit. The other hits are row labels, locators, citations, names and notes saying a value was not used, as in round 1. The merged file's counts equal round 1's, minus its J-2 record, plus the four staging files', keyword by keyword.

**The engines** (31, `engines.json` SHA-256 `f6cda983…4c47`). Round 2 added 13:

| Family | N | Added in round 2 |
|---|---|---|
| LOX/LH2 | 10 | le5b3_qual, le7a_nominal, hm7b_rated, vinci_a5me_full, ce20_gg_nominal |
| LOX/kerosene | 7 | h1b, rs27a, rd107_8d74k |
| LOX/methane | 4 | none |
| Storable | 10 | aj10_190_oms, r4d_apollo_rcs, aestus_eps, viking_4b, rd_253 |

**Usability** (check_specs.py):
- tier0_isp: 27 (LOX/LH2 8, LOX/kerosene 6, LOX/methane 3, storable 10). Round 2 adds two engines without it: vinci_a5me_full (the only chamber pressure found was a subscale diffuser rig's, and it was not used) and rd107_8d74k (no mixture ratio).
- tier0_thrust: 14.
- tier1: none.

**Throats:**
- **Printed (2):**
  - r4d_apollo_rcs: a printed throat area.
  - j2_sa208_504: now NASA SP-8120's printed 14.7 in (0.373 m). It replaces round 1's 0.390 m, derived from an exit diameter whose lip was unstated. The same source gives 700 psi at area ratio 27, against the record's 763 psi. That conflict is flagged, and the operating point was not changed.
- **Derived from a diameter that may be the lip or the envelope (4):**
  - vinci_a5me_full: lip unstated.
  - aj10_190_oms: an "Engine Dia." read as the exit.
  - viking_4b: a "maximum diameter".
  - rd_253: an overall engine diameter, from a secondary source.
- **Missing (8):**
  - le5b3_qual, le7a_nominal, ce20_gg_nominal, h1b, rs27a and rd107_8d74k: no throat or exit dimension is printed (only the area ratio). A throat could be got only from an outcome, which is forbidden.
  - hm7b_rated and aestus_eps: the only diameter is an envelope that may or may not be the exit, so it was not used.

**The throat agent improved J-2 only.** For the other 12 it found no better primary dimension:
- **LM ascent and descent:** both Apollo Experience Reports, read in full, print no diameter.
- **LR87:** the only diameter is the twin-chamber maximum.
- **RS-68A:** the ULA user's guide has no nozzle dimensions.
- **Vulcain 2:** an ESA bulletin confirms the exit diameter but not its lip.
- **J-2X and LR91:** nothing better was found.
- **Fastrac, NK-33, RD-170 and SPS:** the sources were unreadable, had no geometry, or were blocked (DTIC).
- **The Aerojet LOX/LCH4 record** was already exhausted.

**Secondary sources carry more of round 2:**
- **Operating point (p_c, mixture ratio, expansion ratio) from a secondary source only** (Wikipedia, Astronautix, Braeunig): hm7b_rated, aestus_eps, viking_4b, rd_253, h1b, rs27a and rd107_8d74k. The same is true of viking_4b's and rd_253's exit diameters, and so of their throats. HM7B's p_c differs between Wikipedia's infobox and its text. H-1's and RS-27A's pages each give a second, unreconciled p_c.
- **One secondary field:** aj10_190_oms's mixture ratio (Wikipedia).
- **Fully primary:** le5b3_qual, le7a_nominal, ce20_gg_nominal and r4d_apollo_rcs.

The batch's RULES.md must decide whether secondary-only engines count in the 20.

**Not compiled:**
- RD-108: only the RD-107 (8D74K) was taken.
- CE-7.5: CE-20 was taken instead.

No RL10 source was opened (the agents' reports).

**Exposure noted by the agents:** these sources print performance beside the cited geometry:
- EUCASS 2019 (LE-5B-3), Table 1;
- MHI Technical Review (LE-7A);
- Praveen et al., IOP (CE-20);
- the Aerojet Rocketdyne data sheet and NTRS 19850008634 (OMS);
- NTRS 19690029231 (R-4D);
- the Wikipedia infoboxes (HM7B, Aestus, Viking, RD-253, H-1, RS-27);
- Astronautix and Braeunig.

As in round 1, only the agents saw those values.

## Batch b1_tier0 (22 October 2026, 01:44 to 01:47 CDT)

**The tier-0 prediction for all 27 tier0_isp engines is committed and registered.**
- RULES.md is 246d43b, committed before any engine was run. Its header says "about 02:00 CDT"; git gives 01:44.
- P is e747f22. The SHA-256 of `predictions.json` is `e055f68f…72bd`.
- The register line is da89b6e.
- The predictor is the CEA ideal rocket (RocketCEA 1.2.3, `predictions/b1_tier0/tier0_batch.py`).
- Each quantity has a central value and a band from one-at-a-time members (chemistry, divergence, station, open cycle, Cd, throat, conflicting p_c, separation), each member the same for every engine.

**Blindness (Ben, Q1 (a)):** the sealed step is held until CRUCIBLE's own prediction on the same 27 engines (batch b2) is committed and registered. One sealing then scores both. No measured value has been read.

**The groups (Ben, Q2 (a)):** S (an operating-point input from Wikipedia, Astronautix or Braeunig only) has 12 engines, and P has 15. The handoff's 11 had allowed aj10_190_oms as primary; the stated rule counts its Wikipedia mixture ratio, so it is in S. Round 3 (primary sources for S) is still to do, and its records would be a new batch.

**What the run showed about the rules:**
- **The gas generator's rule (900 K, by bisection from O/F 0.02) has no root** for N2O4 with Aerozine-50 or UH25. CEA's equilibrium temperature at O/F 0.02 is already about 1,180 to 1,200 K. lr87_aj_5, lr91_aj_5 and viking_4b take O/F 0.02, the range's end, as recorded in their diagnostics.
- **For LOX/kerosene, equilibrium puts that ratio at about 0.08.**

Neither changes the batch: the rules are as committed, and both are in `predictions.json`.

## Batch b2_crucible: rules only (22 October 2026, committed 02:03 CDT by git; the file's header says about 02:10)

**`predictions/b2_crucible/RULES.md` (9294aa6) was committed before any CRUCIBLE run of these engines.** It is a file of its own, since b1's RULES.md is frozen at its P.
- **The model:** CRUCIBLE's full engine, the same as for the RL10 (the Director, 22 October), on the same 27 engines at cfd6f12.
- **The fills:** missing geometry is filled by batch-wide rules, each value labelled printed, derived or inferred, with its source's status (confirmed, secondary or memory):
  - the throat's size by cycle class, from the batch's own printed throats;
  - CR (Humble's fit), L* by family (Huzel & Huang), the convergent section and throat arcs (SP-8120);
  - a Rao 80% bell (60 and 90% members);
  - the injection states, the preburner stream, a generic coax or doublet by phase, the wall temperature by cooling type.
- **p_c** is hit by iterating the mass flow.
- **Prerequisites:** a table per propellant family.
  - Only the 8 LOX/LH2 engines are close: the RL10 pipeline, plus the preburner stream for rs25 and le7a.
  - CH4, kerosene and the storables need new tables, a liquid-fuel spray and doublets.
- **The runs wait for the RL10 pipeline.** The sealed step waits for b2's P.

## Round 3: primary sources for group S (22 October 2026, about 02:05 to 02:25 CDT)

**Three Sonnet compilers, run in the foreground under `engine_db/compile_brief_2026-10-22_round3.md`, worked on the 12 group-S records.** Their files passed the screen:
- `check_specs.py`: 0 errors;
- `--counts`, then `--masked` on the new strings only: no outcome value, only labels and locators;
- then the full read of the changed input fields.

Their records replace the 12 in a new `engines.json`. The other 19 records are unchanged, as checked by parsing.

**Under b1's group rule, 3 of the 12 now have p_c, the mixture ratio and ε all from primary sources:**
- **j2_sa208_504:** the NASA J-2 fact sheet. It also adds the injector (coax) and its post count.
- **lm_ascent_engine:** NASA TN D-7400, Table I.
- **h1b:** Rocketdyne's data sheet for the SA-206-on block. Its p_c, 700 psia, is b1's rule-8 member, which replaces the secondary 633 psi.

**9 stay in S, each field noted "round 3: no primary source found (searched: ...)":**
- **Partly primary now:**
  - lm_descent_engine: its mixture ratio is a flight-measured value, flagged as not from the same point as its secondary p_c;
  - aj10_190_oms;
  - viking_4b: its p_c is from NASA TM-77683, the Ariane 2 generation's figure, flagged as carried over to the 4B;
  - rs27a.
- **All three still secondary:** lr87_aj_5, lr91_aj_5, hm7b_rated, aestus_eps and rd_253.

The agents note that DTIC was blocked, and that several AIAA papers they found were inaccessible.

**b1 and b2 are unaffected:** both cite cfd6f12, with S at 12. A later batch on the new version recomputes the groups (S would have 9).

**Blindness:** the compilers read sources that print performance beside the inputs: Wikipedia, Braeunig, the Arianespace manual, the AJ10-190 maker sheet and the Apollo 12 descent evaluation. They recorded none of it, and their messages carried no numbers.

## Declarations

- **The compiling agents saw outcome values** in their sources. This exposure is the agents' only: they wrote none (the screen checks it) and they make no prediction.
- **The predictor's prior knowledge.** The predictor (Claude, through the Director's workers) knows the published ratings of the best-known engines from its training. No such value is an input or a check, and every judgement in a prediction (brackets, the tier, the settings) follows a rule fixed for the whole batch before its first engine is predicted, with no per-engine knob (`PREDICTION_PROTOCOL.md`). The exposure is declared with every score.
- **Secondary sources.** Where only a secondary source (Astronautix, Wikipedia) gives a field, the field's note says so.

## Changing a record

Once a prediction cites `engines.json` by its hash, that version is fixed for that prediction. A correction (a better source, a unit fix) is a new commit of `engines.json`, dated here. Predictions made on the old version are scored as made, and a re-prediction on the new one is a new batch.
