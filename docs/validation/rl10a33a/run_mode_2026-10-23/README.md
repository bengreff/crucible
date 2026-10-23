# The run mode's R3 cost probes (23 October 2026, about 00:23 to 00:28 CDT, on the Mac)

Development runs (protocol §6), not the reported case. 149 × 16, 4 threads, to 0.02 ms, Table B `/private/tmp/rl10rebuild/rl10_s.bin`, through the slot limiter. Results in `tier1_run_definition.md`, *The run mode*.

- `probe_2level_0.02ms.log`: 2 levels, on the tree with the wall-face fix and the zeroed dump padding. 31 base steps, 17 regrids, 8 s.
- `probe_3level_0.02ms.log`: 3 levels, with the wall-face fix. It stops in the first step: `Flow::exchange left a cell non-admissible.`

Two lines from temporary diagnostic builds, reverted and not committed:
- The 3-level stop. `Flow::exchange`'s message, extended:
  `q 540 (i 4 j 28 of 1192 x 128) active 1 dist 0 v 2.332e-07 u {4.0322e-01 7.8258e+00 -8.4496e-01 -1.5728e+05} a mass 7.156e-08 mom 4.330e-07 -6.004e-12 en -3.231e-02 turb 2.156e-01 1.165e+02`
  (the state u after the deposit; a is the deposit: mass, axial and radial momentum, energy).
- The 2-level budget (a), by component, at 0.02 ms. The driver's `budgets()` printed each part:
  `gas mass -3.491e-16 energy 2.446e-07 momentum -1.873e-16 (residuals -1.761e-18 kg 1.079e-04 J -6.072e-18; scales 5.044e-03 4.411e+02 3.241e-02)`

## Registration B and the start amendment (23 October 2026, about 00:40 to 00:47 CDT, on the Mac)

- `b3_probe_3level_lead0_0.02ms.log`: B3, the 3-level probe with the exchange bound and t_L = 0 (`CRUCIBLE_T1_LEAD=0 CRUCIBLE_T1_END=2e-5`). It stops in `Flow::exchange` before the first 10 μs sample, as before. CVODES warnings are stripped.
- `b3_flash_diagnostic.txt`: the same probe on a temporary build, reverted and not committed. After each level's move it printed the mass flashed and the vapour exchanged. Only the 5 sub-steps with a flash are kept; the other 120 sub-steps with vapour carry at most 1.4e-13 kg.
- B4, the lead's probe (`CRUCIBLE_T1_LEAD=20e-3 CRUCIBLE_T1_END=21e-3`, 3 levels, 4 threads), launched at 00:47 CDT from `/tmp/t1run/b` (`b4.log`, done-file `done_b4`). Not scored yet.
