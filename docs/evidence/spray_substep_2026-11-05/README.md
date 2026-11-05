# The spray's sub-step defect past t_L (5 November 2026)

Evidence for SPEED_PLAN.md 7.6. All on backhouse (gcc 13), development case only.

- `probe1_lb.log`: probe 1 (`lbprobe.py`, `lb.sh`), print-only, the three liquid budgets after every base step from i10k's 20.81 ms checkpoint. The energy budget jumps to 5.2e-9 in the step ending at 20.8378 ms.
- `probe2_nan.log`: probe 2 (`nanprobe.py`, `lb2.sh`), adding the HPROBE line (a sub-step whose corrector leaves h below the LOX table) and the NaN residual's terms at the throw.
- `f3.sh`: check F3 on the unfixed spray.cpp (A) and on the fix (B), then the spray ctests and tier 1's selftest on B. `f3_run.out` is its summary.
- `spray_tests_A_unfixed.log`, `spray_tests_B_fix.log`, `spray_tests_A_against_B.diff`: F3 FAILS on A and PASSES on B; no other check line differs.
- `ctest_B_fix.log`: face_inlet, chemistry_skip and spray_transfer pass; spray_verification and spray_amr_verification fail only the lines known on this machine (S0(b) 1,245 J/kg; KH-RT 1.044e-12; the AMR spray's dt at 1.087e-12, recorded identically on 23 October in `amr_d/sd3_2026-10-23/ctest_fail_cmp.txt`).
- `tier1_selftest_B_fix.log`: 0 checks failed.
- `f3a.sh`: the AMR spray test on the unfixed spray.cpp, to compare its check lines with B's.
