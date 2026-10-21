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
| `sealed_outcome.schema.json` | The outcome schema. Its files (`sealed/<batch>.json`) are written only by the sealed step, after the batch's prediction commit. None exists yet. |
| `PREDICTION_PROTOCOL.md` | How a blind prediction is committed and later scored. |

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

## Declarations

- **The compiling agents saw outcome values** in their sources. This exposure is the agents' only: they wrote none (the screen checks it) and they make no prediction.
- **The predictor's prior knowledge.** The predictor (Claude, through the Director's workers) knows the published ratings of the best-known engines from its training. No such value is an input or a check, and every judgement in a prediction (brackets, the tier, the settings) follows a rule fixed for the whole batch before its first engine is predicted, with no per-engine knob (`PREDICTION_PROTOCOL.md`). The exposure is declared with every score.
- **Secondary sources.** Where only a secondary source (Astronautix, Wikipedia) gives a field, the field's note says so.

## Changing a record

Once a prediction cites `engines.json` by its hash, that version is fixed for that prediction. A correction (a better source, a unit fix) is a new commit of `engines.json`, dated here. Predictions made on the old version are scored as made, and a re-prediction on the new one is a new batch.
