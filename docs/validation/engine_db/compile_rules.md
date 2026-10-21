# Engine-spec compiler: rules (read all before starting)

You compile INPUT specifications of real liquid rocket engines for a blind-prediction database. A
simulation will later predict each engine's performance from these inputs alone, and the prediction is
committed before anyone looks up the measured performance. So the file you write must contain NO outcome.

## Outcomes: never write them anywhere

Outcomes are: thrust (any: vacuum, sea level, rated, thrust class), specific impulse, c* (characteristic
velocity), thrust coefficient, exhaust velocity, total or per-propellant mass flow, nozzle exit pressure,
and any efficiency (c*, Isp, combustion, energy-release). You will see them in sources. Do not write any of
them in the output file, in scratch files, in notes, in citations' titles (shorten a title that carries one,
e.g. "Development of the [thrust-class] engine"), in engine names ("RL10" is fine, "the 15k-lbf engine" is
not), or in your final message. Do not say "about X kN class" anywhere.

Never derive an input from an outcome. In particular the throat diameter must come from a printed throat
dimension, a drawing, or exit diameter and expansion ratio; never from thrust / (Cf p_c), from mass flow and
c*, or anything similar. If the throat is only obtainable that way, set it missing with the note "only
derivable from an outcome; not derived", and flag it. The mixture ratio may be computed as the ratio of
printed oxidizer and fuel flows, with derivation "ratio of printed propellant flows (flows not recorded)";
the flows themselves are never recorded.

## What to write

One JSON file, at the path your task names, valid against
`/Users/ben/inquiry-project/docs/validation/engine_db/engine_spec.schema.json` (read it; it is the
authority). `/tmp/engine_db_staging/TEMPLATE.json` is a fictitious example record showing the shape.

- Units are fixed by the field name: MPa, m, K, degrees; ratios dimensionless. Convert, and put the printed
  value with its printed unit in `original` (1 psi = 0.00689476 MPa; 1 bar = 0.1 MPa; 1 atm = 0.101325 MPa;
  1 in = 0.0254 m; T[K] = T[R] × 5/9).
- One operating point per engine, chosen where a source gives the inputs and (somewhere) the measured or
  rated performance of the same point. Prefer a test-measured point; say what the point is in `label` and
  `basis`. Keep every input consistent with that one variant and point (do not mix RS-25 versions, or the
  extension-deployed and stowed nozzles).
- Every field carries provenance (printed, digitised, derived, inferred, assumed, missing), and a source id
  plus locator (table, figure, page) when printed or digitised. Prefer primary sources: NASA NTRS reports,
  AIAA/JPC papers, ESA/JAXA/CNES papers, maker data sheets, Huzel & Huang, Sutton. Secondary sites
  (Astronautix, Wikipedia, Encyclopedia Astronautica) are allowed, cited as such; say so in the field's note.
- Do not guess silently. Anything inferred or assumed gets a note saying how, AND an entry in `flags`.
  Every required field you could not source is `missing` with a reason, and flagged. A record with honest
  gaps is worth more than a complete-looking one.
- `outcome_pointers`: where the performance of your operating point is printed (source, table/page, which
  quantity, and the source's label such as 'test' or 'rated'), never the value. A later, separate step reads
  them.
- Flag anything a predictor must know: gas-generator or turbine exhaust that adds to or dumps from the
  nozzle, film cooling, an exit diameter that is the outer lip, a chamber pressure whose station is unknown,
  a mixture ratio that is the engine's rather than the chamber's.

## Checking your file (do this before you finish)

    python3 /Users/ben/inquiry-project/docs/validation/engine_db/check_specs.py FILE
    python3 /Users/ben/inquiry-project/docs/validation/engine_db/check_specs.py FILE --masked

Fix every schema error the first prints. The second lists every string that contains an outcome keyword,
digits masked: make sure each is only a label or a locator (e.g. "Table 3, row 'vacuum thrust'"), never a
value; delete any that carries one. Do not write to the repository; write only your output file.

## RL10: off limits

The RL10A-3-3A awaits its own blind prediction. Include no RL10 variant, and do not open RL10 sources. If a
source you open covers the RL10 (a survey table, say), skip that part unread and do not record anything
from it.

## Final message (to the Director)

At most 200 words, with no performance number of any kind and no input values: the file path; the engine
ids you wrote; for each, a word on what is missing or flagged (e.g. "throat missing", "contour unknown");
candidates you dropped and why; anything you saw that you think breaks blindness (e.g. a source that prints
performance beside geometry: name it, but not its numbers).
