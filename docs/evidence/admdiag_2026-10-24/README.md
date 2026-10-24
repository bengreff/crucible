# The 1-level lead's "could not advance an admissible state" at 1.78 ms (diagnosis, 24 October 2026)

The Director (01:40): diagnose it, "a failure we don't understand".

**The failure.** Tier 1 development ran at 1 level with the 20 ms H2 lead and 4 threads, on the Mac (`/tmp/lead1_run`, launched at 23:39 CDT on 23 October). It stopped between 1.78 and 1.79 ms with "AMR: could not advance an admissible state after 14 base-step reductions" (exit 1). The binary was built at 23:39 from `main` 1166721 (27b135f's engine). The reflux clip (716c3d2, AMR.md's amendment: the reflux clips rho k at zero) reached `main` at 23:50, so that binary did not have it. B4 failed with the same error on 27b135f at 30 to 40 μs, and inquiry2 traced that failure to the reflux leaving rho k negative.

**Hypothesis.** The same fault, so the clip removes it.

**Criteria, stated before the runs.** `admdiag.py` is a print-only patch: each rejected base step prints the check that failed and the failing cells' state. Its two runs:
- **D1:** 27b135f plus admdiag, the same case to 1.79 ms.
- **D2:** `main` plus admdiag, the same case to 1.80 ms.

The hypothesis holds if:
- **(i)** in D1, the error recurs at the same step, and every one of its 14 attempts fails at the reflux's check with rho k < 0 in a coarse cell beside the coarse–fine boundary;
- **(ii)** D2 passes 1.79 ms with exit 0, and its reflux clips at least one cell in [1.78, 1.79] ms.

If (i) shows another check, or (ii) fails, the hypothesis is wrong, and the printed state is the next lead.
