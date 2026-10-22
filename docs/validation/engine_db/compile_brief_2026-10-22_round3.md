# Brief for an engine-spec compiler (22 October 2026, round 3: primary sources for group S)

Read these in full BEFORE starting, and follow them exactly:
1. /Users/ben/inquiry-project/docs/validation/engine_db/compile_rules.md (the rules; they govern)
2. /Users/ben/inquiry-project/docs/validation/engine_db/engine_spec.schema.json (the authority on shape)
3. Your START file (named in your task). It holds your engines' current records, which hold inputs and outcome pointers (locators only, never values).

## The job

Each engine in your START file has at least one operating-point input (chamber pressure, mixture ratio or expansion ratio) whose only source is secondary: Wikipedia, Astronautix (Encyclopedia Astronautica) or Braeunig. For each engine, in this order:

1. **The operating point.** For every field among chamber_pressure_mpa, mixture_ratio (and its basis) and expansion_ratio whose source is secondary, look for a primary source: NASA NTRS, AIAA/JPC papers, ESA/CNES/DLR papers, NASA SP or TM reports, maker data sheets or manuals, Huzel & Huang, or Sutton.
   - If you find one, replace the field's value, provenance, source and locator.
   - If the primary value differs from the secondary one, keep the secondary value only in the field's `note` (it is an input, so that is allowed), and add a flag.
   - If you find none, leave the field as it is and add a note: "round 3: no primary source found (searched: ...)".
2. **The geometry and the injector.** Where a primary source prints them, fill or improve:
   - the throat diameter or area, the contraction ratio, the chamber diameter;
   - the chamber length or L*, the contour, the cooling;
   - the injector type and element count, the injection temperatures.

   The compile rules on the throat hold: it comes from a printed dimension, a drawing, or the exit diameter and ε, never from an outcome.
3. **The operating point must stay one point.** If the primary source's point differs from the current record's (another variant, another power level), choose ONE point that has primary inputs. Make every input consistent with it, update `label` and `basis`, and update `outcome_pointers` to where that point's performance is printed (locators only).
4. **The ids stay as they are,** unless the variant changes. If it does, use a new id and say so in your final message.

You may also use the existing downloads in /tmp/engine_db_staging/ (the earlier rounds' files are named by family or engine). New downloads go there with the prefix `r3_` in the file name.

## Clock and saving

- Run `date` (CDT) before each new engine.
- Start no new engine after the start-by time in your task. Have your final, checked file and final message done by the finish-by time.
- After EACH engine, rewrite your whole OUTFILE (valid JSON, all your engines so far, the unfinished ones as in START) and run `python3 /Users/ben/inquiry-project/docs/validation/engine_db/check_specs.py OUTFILE`. Fix the errors before the next engine.
- At the end, run it again with `--masked` and remove any string that carries an outcome value.

## Limits

- Write nothing in the repository (/Users/ben/inquiry-project). Write only OUTFILE and your downloads in /tmp/engine_db_staging/.
- Never write /home/greff/inquiry-project or anything on other machines.
- The RL10 is off limits (compile rules).
- **No outcome value in the file, in notes, in scratch files or in your final message.** Outcomes are thrust, Isp, c*, C_F, exhaust velocity, mass flows, exit pressure and any efficiency. These engines' predictions are not all committed yet.

## Final message

At most 200 words, with no performance number of any kind and no input values:
- the OUTFILE path;
- for each engine: which operating-point fields are now primary, which are still secondary, and which geometry fields you added (names only);
- any change of variant or point;
- anything that may break blindness (name the source, not its numbers).
