# Brief for an engine-spec compiler (21 October 2026, rerun after the 03:19 kill)

Read in full BEFORE starting, and follow exactly:
1. /Users/ben/inquiry-project/docs/validation/engine_db/compile_rules.md  (the rules; they govern)
2. /Users/ben/inquiry-project/docs/validation/engine_db/engine_spec.schema.json  (the authority on shape)
3. /tmp/engine_db_staging/TEMPLATE.json  (a fictitious example record)

Clock and saving (the previous run was killed mid-write and lost everything):
- Run `date` (CDT) before each new engine. Start no new engine after 04:00 CDT. Have your final, checked
  file and your final message done by 04:10 CDT at the latest.
- After EACH engine is complete, rewrite your whole output file (valid JSON, all finished engines so far)
  and run `python3 /Users/ben/inquiry-project/docs/validation/engine_db/check_specs.py OUTFILE`; fix errors
  before the next engine. A file with 4 clean engines beats one with 8 half-done.
- Then, at the end, run it again with `--masked` and remove any string that carries an outcome value.

Sources and scratch:
- /tmp/engine_db_staging/ holds raw downloads from the killed run (PDFs, .txt, .html; mixed families). You
  may reuse the ones for your family's engines instead of re-downloading. New downloads go in
  /tmp/engine_db_staging/ with your family prefix in the file name (e.g. lh2_vinci_paper.pdf).
  `pdftotext` may be available; `curl -sL -o` works for NTRS (https://ntrs.nasa.gov/api/citations/<id>/downloads/<file>.pdf).
- Prefer primary sources (NTRS, AIAA/JPC, ESA/JAXA/CNES, maker data sheets, Huzel & Huang, Sutton). Secondary
  sites are allowed, cited as such in the field's note.
- Write nothing in the repository (/Users/ben/inquiry-project). Write only OUTFILE and your scratch downloads.
- Never write /home/greff/inquiry-project or anything on other machines.

Choosing engines: pick the engines of your family that have the best-documented INPUTS (chamber pressure,
mixture ratio, expansion ratio, throat or exit diameter, contraction ratio, chamber length or L*, contour,
injector, injection temperatures), each at one operating point whose performance is printed somewhere (you
record only a pointer to it). Aim for 5 to 7 engines; quality and honest provenance over count. Use engine
ids that are lowercase with underscores and specific to the variant (e.g. rs25d_109pct is NOT allowed if the
power level is a thrust percentage; use rs25d_rpl or describe the point in `label` without an outcome).
Careful: a power level stated as "% of rated thrust" is an outcome-derived label; describe the point by its
chamber pressure or as "rated power level" instead.

Throat: printed throat diameter or area, a drawing, or derived from a printed exit diameter and expansion
ratio (say whether the exit diameter is inner or outer lip). Never from thrust or flow. If it cannot be had
that way: missing, with the note "only derivable from an outcome; not derived", and flagged.

Your final message (at most 200 words, as compile_rules.md says): no number of any kind about performance, no
input values. The file path; the engine ids; for each, a word on what is missing or flagged, saying
explicitly whether its geometry (throat, contraction, chamber length/L*, contour) is sourced or not; dropped
candidates and why; anything that may break blindness (name the source, not its numbers).
