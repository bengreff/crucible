# How a blind prediction is committed and scored

Written 21 October 2026, before any prediction and before any outcome was compiled. It applies to every batch of the engine database (`../ENGINE_DB.md`). Changes are dated here, and a change never applies to a batch whose prediction is already committed.

## The order

1. **Specs frozen.** `engines.json` is committed (commit S). The batch names its engines (a subset is allowed) and its tiers.
2. **Rules fixed, then predictions made.** Before the batch's first engine is predicted, its rules are committed in `predictions/<batch>/RULES.md`: the tier, the engine and table versions, and for every missing or flagged input the rule that fills or brackets it (for example, "a missing contraction ratio takes the batch's value of 3, bracketed 2 to 5"). A rule is the same for every engine of the batch, with no per-engine choice.
3. **The prediction commit (P).** `predictions/<batch>/predictions.json` holds, per engine: the predicted quantities (vacuum and sea-level Isp and thrust, c*, ṁ), each with its band and how the band was built; the spec's SHA-256 and S; the CRUCIBLE commit, the table hashes and the CEA version used. P's commit message quotes the SHA-256 of `predictions.json` (`shasum -a 256`).
4. **The register.** The next commit adds a line to `REGISTER.md`: the batch, P's full hash, the SHA-256, and the date. P cannot hold its own hash, so this comes after it.
5. **The sealed step.** It starts only after step 4. A separate agent, whose brief contains no predicted value and which does not open `predictions/`, compiles the batch's outcomes into `sealed/<batch>.json` (`sealed_outcome.schema.json`). Before writing, it checks that P is an ancestor of HEAD (`git merge-base --is-ancestor P HEAD`) and that `git show P:docs/validation/engine_db/predictions/<batch>/predictions.json | shasum -a 256` equals the register's value. Its selection rule:
   - values at the spec's operating point (its `label`), from the spec's `outcome_pointers` first;
   - test or flight values over rated ones, each with its source's own label;
   - any other printed value for the same point goes under `alternates` and is never scored.
   Its commit (O) quotes P.
6. **The scoring.** The scorer runs the same checks: P is an ancestor of O, `predictions.json` is unchanged since P (`git diff --quiet P -- <file>`), and the spec's hash in the prediction matches `engines.json` at S. Then, per engine and quantity, it reports the signed error (predicted − measured) / measured and whether the measured value lies inside the band. Over the batch, by family and overall, it reports: N, the mean signed error, the RMS, the median and largest absolute error, and the fraction inside the band. Test and flight comparisons are reported apart from rating and model ones; only the first count as validation (as RL10A33A_PROTOCOL §5). No pass threshold is set; the spread is the result.

## What ends blindness

- **Opening the sealed file for an engine ends blindness for that engine, for every later tier.** So a batch states its tiers before step 5. An engine meant for a later tier-1 run stays out of the tier-0 batch's sealed file, or its tier-1 result is labelled not blind.
- **Nothing is retuned after the reading.** A model change made after a score is a new batch on the same engines, labelled not blind.
- **Exposures are declared in RULES.md before P:** anything read or recalled about an outcome of the batch (the predictor's training knowledge of famous engines is declared once for every batch; `../ENGINE_DB.md`, *Declarations*).

## What the record proves, and what it does not

Git's ancestry and the content hashes fix the order: the predictions in P could not have been changed after O. Commit dates are the local clock and prove nothing on their own. An outside timestamp would make P's date checkable by a third party: for example, pushing P to `origin` before step 5, or sending P's hash somewhere dated. `main` is not pushed routinely (157 commits ahead of `origin` on 21 October), so that is Ben's call, and until he makes it the local order is the record.
