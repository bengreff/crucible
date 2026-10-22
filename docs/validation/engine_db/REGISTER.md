# Prediction register

One line per batch, written in the commit after the batch's prediction commit P (`PREDICTION_PROTOCOL.md`, step 4). The sealed step checks that P is an ancestor of HEAD and that `git show P:<file> | shasum -a 256` equals the SHA-256 here.

| Batch | P (full hash) | File | SHA-256 | Date |
|---|---|---|---|---|
| b1_tier0 | e747f22255ee2b9bfe0cbe9beebde49c3403065f | docs/validation/engine_db/predictions/b1_tier0/predictions.json | e055f68f9bcd85eb966fe97213097d1f47049a6067eabd956cad5cecde6572bd | 22 October 2026, about 02:20 CDT |

**b1_tier0's sealed step is held** until batch b2 (CRUCIBLE's own prediction on the same 27 engines) has its own line here (`predictions/b1_tier0/RULES.md`, *Tiers and blindness*).
