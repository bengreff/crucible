# Case C runs (C1 cold start, `crucible_chamber_study eqt`, one thread, 24 October 2026)

Each directory holds `<name>.log`, `<name>_history.csv.gz` and, for implicit runs, `<name>_implicit.csv` (one row per step: iterations, first and last defect, rejections so far). `probes/` holds the development runs' logs only. `c_compare.py` (one level up) reads the gzipped histories directly.

The scripts: `base.sh` (explicit baselines through `slot.py`), `e.sh` (explicit cost baseline), `x.sh` (explicit at a set CFL), `imp.sh` (implicit: name nz nr end dt tolerance max_iterations pseudo_cfl extrapolate inside), `w.sh` (criterion 5's switch runs). Their binaries were copies of `build/crucible_chamber_study` in `/tmp/implicit/c`:
- `bin` (02:16): e8816dc's build.
- `bin3` (02:30): the projection inside, with the old known part.
- `bin4` (02:38): the second amendment (the new known part, the throughput print, `CRUCIBLE_CFL`).
- `bin5` (02:51): plus `CRUCIBLE_IMPLICIT_WINDOW`.
- `bin6` (02:56): plus the projection's failures as rejected attempts (third amendment).
- `bin7` (04:16): 063294f plus the period-mean history (`PERIOD_MEAN.md`); `p.sh` runs it.
- `momdiag/` (03:30 and 03:39): print-only builds of 96635e8 with `../momdiag.py` applied, in `/tmp/momdiag`; `run.sh` and `run2.sh` there launched them. Only the `MOMDIAG` lines and the logs are kept.

| Run | Mesh | Kind | Step | Tolerance | Build | What it is for |
|---|---|---|---|---|---|---|
| e64 | 64x12 | explicit CFL 0.4 | | | bin | τ_C1, the reference for criteria 3 to 5 |
| e64b | 64x12 | explicit CFL 0.4 | | | bin4 | M2 and M3's baseline (same batch as c25 to c150) |
| x64h | 64x12 | explicit CFL 0.2 | | | bin4 | criterion 5's diagnosis: the explicit's own step error |
| e128 | 128x24 | explicit CFL 0.4 | | | bin | the reference and baseline on 128x24 |
| e256 | 256x48 | explicit CFL 0.4, to 1 ms | | | bin | criterion 6's measured step ratio (through `slot.py`, 03:08 to 03:32) |
| t63, t126 | 64x12 | implicit, split | τ/63, τ/126 | 1e-6 | bin | the split is first order at the injector face |
| i63, i126 | 64x12 | implicit, inside | τ/63, τ/126 | 1e-6 | bin3 | the old known part's budget drift |
| k63, k126 | 64x12 | implicit, inside | τ/63, τ/126 | 8.975e-7 | bin4 | E_t (criterion 3), k63 the sweep's reference |
| s10, s100 | 64x12 | implicit, inside | τ/63 | 8.975e-5, 8.975e-6 | bin4 | criterion 3's sweep |
| c25, c63, c150 | 64x12 | implicit, inside | τ/25, τ/63, τ/150 | 8.975e-5 | bin4 | criteria 2, 4, 5; M1 to M3 |
| m25, m63 | 128x24 | implicit, inside | τ/25, τ/63 | 8.975e-5 | bin4 | criteria 2, 4, 5; M1 to M3 |
| m150b | 128x24 | implicit, inside | τ/150 | 8.975e-5 | bin6 | as m63 (m150 on bin4 stopped at 0.147 ms; its log is in probes/) |
| m25t | 128x24 | implicit, inside | τ/25, to 5 ms | 1e-6 | bin6 | is m25's state the tolerance's? (no) |
| v63 | 256x48 | implicit, inside | τ/63, to 1 ms | 8.975e-5 | bin4 | criterion 6 |
| w1 | 64x12 | implicit to 7.9 ms, then explicit | τ/63 | 8.975e-5 | bin5 | criterion 5's diagnosis |
| w2 | 64x12 | explicit to 7.9 ms, then implicit | τ/63 | 8.975e-5 | bin5 | criterion 5's diagnosis |
| c63t1, c63t2 | 64x12 | implicit, inside | τ/63 | 8.975e-6, 8.975e-7 | bin4 | criterion 5's cause: is the offset the inner tolerance's? (no) |
| momdiag/c63, e64, x64h | 64x12 | as c63, e64 and x64h | | | momdiag bin | criterion 5's cause: the booked momentum rates at each sample from 7.98 ms |
| momdiag/e64s | 64x12 | as e64 | | | momdiag bin2 | criterion 5's cause: each explicit step's two stages from 7.996 ms |
| momdiag/e128s, m63 | 128x24 | as e128 and m63 | | | momdiag bin2 | criterion 5's cause on 128x24 (e128s's prediction written first, `momdiag/e128s_prediction.txt`) |
| p64 | 64x12 | explicit CFL 0.4, as e64 | | | bin7 | the period-mean reference: criteria 4 and 5 re-judged (`rejudge/`); its march `cmp`-identical to e64's (`rejudge/march_cmp.txt`) |
| p128 | 128x24 | explicit CFL 0.4, as e128 | | | bin7 | the period-mean reference on 128x24: criteria 4 and 5 re-judged (`rejudge/`); its march `cmp`-identical to e128's |

τ_C1 = 441.243 μs: τ/25 = 1.764973e-5 s, τ/63 = 7.003860e-6 s, τ/126 = 3.501930e-6 s, τ/150 = 2.941621e-6 s. Every implicit run here uses pseudo CFL 3 and the extrapolated first iterate, with at most 400 inner iterations.
