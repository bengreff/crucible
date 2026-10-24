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

## The result (24 October 2026, judged about 03:20 CDT): the hypothesis is wrong

D1 ran for 1,099 s of wall time, ending at 02:29:49 CDT when D2 was queued. D2 ran from 02:29:49 to 02:38:32 (523 s). Both used 4 threads under the slot harness. Their logs, series and done-files are in `logs/`. D1's done-file: `exit 1`. D2's: `exit 0`.

**D1 exits 1 on the same error at the same place.** Its last step, at t = 1.782158629 ms, exhausts its 14 attempts: "AMR: could not advance an admissible state after 14 base-step reductions". Every attempt fails the stage check rho omega > 0 in level-1 cell 6604 (z 1.0400 m, r 0.1496 m), after the second stage. In every attempt rho k is positive and the internal energy is above its floor. As the step halves from 6.67e-8 s to 3.3e-11 s, that cell's rho omega goes from -766 to -0.37 kg/(m3 s), halving with each halving. So rho omega starts the step at about zero, and the stage's rate on it is negative (about -1.1e10 kg/(m3 s2)) and does not vanish with it. No step size makes it positive.

The earlier rejections fail the same check. From 1.772660 ms, 27 steps have rejections, 189 rejected attempts in all. Every one fails rho omega > 0 at a stage:
- level 0, first stage, 59 attempts;
- level 0, second stage, 15;
- level 1, second stage, 115.

The cells involved are level-0 cell 1637 (column 102, radial row 5 of 16, z 1.032 m, r 0.131 m) and level-1 cells 6604 and 6572 (z 1.035 m, r 0.149 m). They are interior cells of the cold H2 expansion, about 0.7 m downstream of the throat (station 33, z 0.33 m). They are not wall cells. The print does not show whether they border a coarse–fine boundary. Every earlier step recovers by halving. The reflux's check never fails, and rho k is never below zero.

**(i) fails:** the check is the stage's rho omega, not the reflux's rho k.

**D2 passes 1.80 ms with exit 0.** It has no rejected base step and no reflux clip over the whole run (0 ADMDIAG lines). The clip's print is inside `if (rho k < 0)` and is in D2's binary, so the clip never fired. **(ii)'s second half fails.**

**What D2's pass does and does not show.** D1's and D2's series agree to 1e-10 relative or better in every column to 1.77 ms; they part at 1.78 ms, after D1's rejections change its steps. The clip is not why D2 passed, because it never fired. Between 27b135f and 4f2d81e the engine also gained the chemistry skip, which skips every cold H2 cell on `main`, where 27b135f integrated them with CVODES. It also gained the spray fix, which does nothing in the lead because no LOX flows before t_L. The cells involved hold about 1e-3 kg/m3 and contribute almost nothing to the series. So a local difference in rho omega there can be large while the series agree.

That `main` carries the same fault and missed it here by a perturbation is an inference, not shown.

**The next lead (the printed state):** why rho omega falls to about zero in an interior cell of the cold expansion, and which part of the stage's rate on rho omega stays negative as it vanishes.

`Flow::clipTurbulence` clips rho k at zero and nothing clips rho omega, so a stage that takes rho omega below zero always rejects. The test is a second print-only patch. It would print that cell's rate on rho omega split by face (advective and diffusive, each face) and, from the source half-steps, rho omega before and after with the SST source terms, rerunning D1 to 1.79 ms. It was not built or run tonight.
