# The blind engine database

Status: started 21 October 2026, about 03:10 CDT, independent of inquiry2's RL10 lane. It serves the inquiry trust bar (`~/director/goals.md`): blind predictions of Isp and thrust on 20 or more real engines across all liquid propellants, committed before their measured values are looked up, with the error spread reported. It is TECHNICAL_PLAN step 12's "engines as data", with the outcomes split off.

## Files (`docs/validation/engine_db/`)

| File | What it is |
|---|---|
| `engine_spec.schema.json` | The input schema: one record per engine, a provenance on every field. It has no outcome field, and every object is closed, so an outcome cannot enter under a new key. |
| `engines.json` | The screened specs (*Compilation*, below). |
| `check_specs.py` | Schema and consistency check, the outcome screen (count-only, masked view, redact unread) and the usability table. Stdlib only. |
| `compile_rules.md` | The rules the compiling agents were given, verbatim. |
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

Four Sonnet agents, one per family (LOX/LH2, LOX/kerosene, LOX/methane, storables with at most one other), worked from `compile_rules.md`. Each wrote a staging file outside the repo, checked it against the schema and its own masked view, and reported back with no numbers. The RL10 family is excluded: RL10A-3-3A awaits its own prediction. RL10 variants may join after that prediction's commit.

**The screen, before any record was read** (as for RL10, `rl10a33a/injector_element_inputs.md`):
1. Count-only greps: `check_specs.py --counts` prints the hits per outcome keyword, with no text, over every free-text string (notes, the printed originals, derivations, locators, citations, labels, names, flag reasons, URLs).
2. A masked view of every hit: `check_specs.py --masked` prints its path and its text with each run of digits shown as `#`.
3. A hit that carries a value is removed unread (`--redact PATH`; a number field becomes `missing`), and logged below by path only.
4. Only then is the file read in full and merged into `engines.json`.

*The screen's results and the engine list: below, once the agents report.*

## Declarations

- **The compiling agents saw outcome values** in their sources. This exposure is the agents' only: they wrote none (the screen checks it) and they make no prediction.
- **The predictor's prior knowledge.** The predictor (Claude, through the Director's workers) knows the published ratings of the best-known engines from its training. No such value is an input or a check, and every judgement in a prediction (brackets, the tier, the settings) follows a rule fixed for the whole batch before its first engine is predicted, with no per-engine knob (`PREDICTION_PROTOCOL.md`). The exposure is declared with every score.
- **Secondary sources.** Where only a secondary source (Astronautix, Wikipedia) gives a field, the field's note says so.

## Changing a record

Once a prediction cites `engines.json` by its hash, that version is fixed for that prediction. A correction (a better source, a unit fix) is a new commit of `engines.json`, dated here. Predictions made on the old version are scored as made, and a re-prediction on the new one is a new batch.
