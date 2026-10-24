# History thrust as a period mean: the fix, its test, and the audit of past results

Written 24 October 2026, about 04:20 to 04:40 CDT. Asked for by the Director at 04:20: "Fix the explicit history so it reports period-mean thrust, with a test. Then audit: list every past result that used history thrust (Table B's c*/Isp comparisons, the S-series, tier-0 or anything in the blind protocol) and say for each whether the error could matter at its bar. List only tonight; fixes next session."

The defect is described in `IMPLICIT_PROTOTYPE.md` (*Criterion 5: the cause*):
- `Flow::measurements()` reports the boundary rates of the last step only.
- A history samples with `advanceTo`, so that step is shortened to land on the sample.
- With a split after the flow stages (Table A's projection: flow, then relax), the settled state is not a zero of the spatial operator. The two stages carry opposite imbalances, about ±7.1 N of axial momentum at C1's injector face at CFL 0.4.
- Full steps cancel them; the landing step books about the first stage's imbalance × (1 − its fraction of a full step). That is up to 0.42% of C1's vacuum thrust.

## The fix *(built, 24 October, about 04:12 to 04:16 CDT)*

- **`Flow` integrates its booked rates.** `integratedRates_` is the time integral since t = 0 of each booked rate (mass flows, momentum fluxes, axial forces, heat flows).
  - The explicit march adds dt × the stage mean in `Flow::accept`, the same increment the budgets use. The implicit adds dt × its step's rates.
  - It is reset with the budgets (`resetAccounting`), saved and restored with them across a redone AMR base step (`Amr::save`, `Amr::restore`), and written to the checkpoint. The checkpoint's format therefore changes: a checkpoint from an older binary does not load.
- **`Measurements::integrals`** reports them. **`periodMean(now, previous)`** returns `now` with the boundary rates replaced by their means over (previous.time, now.time], and the device and exit-plane thrusts recomputed from those means.
- **`crucible_chamber_study`** now writes its history, the settling test and the end summary from each sample's period mean.
- **The march is unchanged.** Only what the history reports changes. Checked on p64 against e64 and on p128 against e128 (*Records*, below).
- **Not fixed tonight** (the Director's "fixes next session"):
  - AMR's composite: tier 1's series (`tests/rl10_tier1.cpp`) and `tests/amr_c4.cpp` sum each level's `measurements()`, still the last step's.
  - Each level's `Flow` now carries the integrals, so the port is a sum of their changes over the sample period (inferred: each boundary face belongs to one level, as the sums assume today).

## The test *(core_verification, measured on the Mac, 04:16 CDT)*

- The small supplied chamber of `tests/core_tests.cpp`, filling. There are 20 sample periods of 13 μs from 0.2 ms, each of several steps of about 3.1e-7 s (from the run's printed step), the last one shortened to land on the sample.
- **Over each period, the period-mean net axial momentum rate times the period equals the change in the gas's P_z** within 6.6e-17 N s (bar: 1e-10 of the inlet's momentum per period, 3.5e-3 N s).
- **The last step's rates miss it by up to 5.5e-4 N s,** 16% of that scale. The test requires at least 1,000 times the period mean's miss, so it fails if the old reporting comes back.
- A period of one step gives that step's rates (within 1e-9).
- `crucible_tests` passes: "All core verification checks passed." The fast ctest on this build is in *Records* below.

## The audit: past results that used history thrust

The list was surveyed by a read-only agent and checked here where marked. Its verdicts were corrected for one point:
- **A mean of samples does not remove the bias** when the steps are steady.
- Every sample's landing step then has the same fraction, so every sample carries the same bias.
- e64 shows it: every 2 μs sample had f = 0.806 and +1.367 N, so its mean over [6, 8] ms carried the full bias *(measured)*.
- Only a time integral (a period mean) removes the bias.

**Sizes measured on C1** (Table A, inviscid, CFL 0.4; `IMPLICIT_PROTOTYPE.md`):
- **Vacuum thrust and Isp:** from 0 up to the first-stage imbalance, about 7.1 N of 1680 N (0.42%), depending on the landing fraction.
- **The wall's axial force:** its two stages differ by 0.06 N of 4471 N, so its bias is at most about 1e-5.
- **The outlet's momentum flux:** its stages differ by 0.002 N (1.2e-6). The outlet's mass flow, c*'s denominator, was not printed. Its bias is taken to be of the same order (inferred), so **c\* is not affected at any bar below**, and only F_vac and Isp are.

The bias scales with the step (x64h: half the CFL, half the imbalance), and other splits and cases are not measured.

| # | Result (file, section) | Quantity, and how it was taken | Split; mesh, CFL | Bar | Verdict |
|---|---|---|---|---|---|
| 1 | `TABLE_A.md`, criterion 3 (and its restatement to 6,000 K) | c*, vacuum Isp and F_vac at the settling windows' ends: samples | Table A against direct Cantera equilibrium, both flow then relax; C1 32x6 and 64x12, CFL 0.4 | 0.05% | **Isp and F_vac: could matter; c\*: cannot** (inferred, above). The bar is below the 0.42% ceiling. The two runs share the mesh and CFL, so their biases may nearly cancel, but this is not shown: their steps differ slightly, and so do their landing fractions. |
| 2 | `CHAMBER_C1.md`, criteria 3 and 4 | c* and vacuum Isp (from F_vac) at window ends, against the 1-D CEA bracket | Local equilibrium, flow then relax; C1 32x6, 64x12, 128x24 | c* 0.5%, Isp 1.0% | **c\*: cannot matter** (inferred, above). **Isp: could matter** at the margin (0.42% against 1.0%, with no cancellation against CEA). Criterion 4 already fails as stated, for grid convergence; part of its unexplained increments could be this bias (inferred). |
| 3 | `PASR_C2.md`, criterion 4 | settled c* and Isp, closure against control, "to every printed digit" | Finite-rate (CVODES), split after the flow; same mesh, same CFL | agreement to the printed digits | **Unknown, likely cannot matter.** A common bias cancels only if the two runs land their samples at the same fractions. The comparison with CEA is reported, not judged, so it cannot matter. |
| 4 | `WALL_FUNCTIONS.md`, criterion 3(a) and (b) | c*, vacuum Isp and thrust at window ends; wall heat flow and axial force as 25-sample window means | eqtt (Table A, SST, wall functions); the reference and wall-function meshes differ | c*/Isp/thrust 0.02% and 0.2%; heat 10%, force 2% | **Isp and thrust: could matter** (bars at or below the ceiling; different meshes give different landing fractions). **c\*, heat and force: cannot matter** (the force's bias is about 1e-5, measured on C1; the heat flow's stage difference is not measured, but 10% is far above any of these). |
| 5 | `AMR.md`, criterion 4 and the Director's additions | settled c*, Isp and thrust as 21-sample means over the last 1 ms, AMR against uniform 64x12 | eqtt, Table A; AMR levels | 0.1%, and 0.02% for convergence | **Isp and thrust: could matter, and may already be on record.** AMR.md notes that Isp at 8 ms is 0.18% to 0.41% above the last millisecond's mean on every run, the uniform one included, and attributes it to the flow still moving (inferred there). That is the size of the landing-step bias, so it may be this effect instead (inferred, not checked). The 21-sample means keep the bias where the fractions repeat. |
| 6 | `SPEED.md` 5c, lever B (L4/L8, H4/H8/H4s; criteria (a) and (b)) | F_vac, outlet flow and p0 series every 0.01 ms, RMS and end differences | Tier 1 development (AMR, spray, finite-rate manifold); CFL 0.8 against 0.4 | RMS 0.2% and others | **Could matter.** The imbalance scales with the step, so CFL 0.8 carries about twice CFL 0.4's (inferred from x64h), and the fractions differ. The runs are on 4f2d81e (last-step series). Judge the F_vac series only after the fix is ported, or on p0 alone, which is a state value. |
| 7 | `SPEED.md` 5b, the chemistry skip | skip-on against skip-off series | Same steps on both sides; the skip changes no step | 3.3e-8 relative | **Cannot matter.** The steps and landing fractions are the same on both sides. |
| 8 | `RL10_DEV.md` | development thrust beside CEA | development | none (not judged) | **Cannot matter** (no bar). `chamber_study`'s comment of 13 October records this effect in its limit: a step of about 1e-18 s booked its first stage alone and moved the development thrust by +7.4e-4. Only that step was removed then. |
| 9 | `IMPLICIT_PROTOTYPE.md`, criteria 4 and 5 | F_vac, outlet flow, p_inj against the explicit history | Table A; 64x12 and 128x24 | 0.05% (5); RMS 0.5%, max 2%, order 1 (4) | **Mattered.** It is the cause of criterion 5's FAIL as judged. Both are re-judged on the period-mean reference tonight (that file). |
| 10 | `SPRAY.md` S0 to S6; `docs/validation/` Mascotte A10 P1, P2 | temperature and OH fields; the outlet flow only as a 5% settling gate | spray, finite-rate | P1, P2 as registered | **Cannot matter** (not thrust, Isp or c*; the gate is far above the bias). |
| 11 | RL10 tier 0 (`docs/validation/rl10a33a/`) | CEA and analytic only | none | as registered | **Cannot matter** (no engine history). |
| 12 | RL10 tier 1, the reported case (`tier1_run_definition.md`, protocol) | F_vac, Isp, c* and the wall heat as batch means of 10 μs samples, summed over AMR levels | AMR, spray, finite-rate manifold, SST with wall functions | as registered (blind: no values here) | **Could matter.** Batch means of landing-step samples keep the bias where the fractions repeat. No prediction or tier-1 output exists yet, so the fix (ported to AMR's composite and the driver's series) goes in before any tier-1 output. |
| 13 | `docs/validation/ENGINE_DB.md`, `BACK_MASSIER_GIER_INPUTS.md`, `VALIDATION_SURVEY.md`, `VALIDATION_TUM_ROUND.md` | inputs and surveys | | | **Not applicable** (no engine thrust). |

Table B: none of the survey's entries is a Table B c* or Isp comparison judged from a history. Table B's own criteria (`docs/evidence/table_b/`) were not separately surveyed tonight. They go first on next session's list, with items 1, 4, 5 and 6.

**The code paths that report last-step rates** (the survey's list, read here):
- `core/flow.cpp`, `Flow::measurements()`: by design still the last step's, now beside `integrals`.
- Its consumers that sample it: `tests/rl10_tier1.cpp` (the series, about line 429), `tests/amr_c4.cpp` (about lines 190 and 302), `tests/amr_budget.cpp` and `tests/spray_amr_tests.cpp` (budgets, not affected).
- `tests/chamber_study.cpp`: fixed tonight.

## Records

- **The march is unchanged, on both meshes** (`implicit_prototype/c/rejudge/march_cmp.txt`, 04:24 and 04:41 CDT):
  - p64 against e64: 17 files `cmp`-identical (16 field files and the mesh), the same 136,227 steps. The history's state columns are identical; only inlet, outlet, F_vac, isp_vac_supply and wall_axial_force differ.
  - p128 against e128: 17 files `cmp`-identical, the same 275,738 steps.
- **The period-mean thrust is the outlet's flux, as predicted:** p64 1679.83 N and p128 1680.49 N over [6, 8] ms, the values derived from the outlet flux at 03:45 to 04:10 (1679.831 N and 1680.488 N).
- **The fast ctest on this build** (`implicit_prototype/ctest_2026-10-24_period_mean.txt`, 04:20 to 04:27 CDT): only the 4 known failures (amr, amr_subcycling, transport and spray verification); 11 of 15 pass, `core_verification` and `spray_amr_verification` among them. Whether any test redoes an AMR base step (which `Amr::save` and `restore` now cover) was not checked.
- **The re-judgement of the implicit prototype's criteria 4 and 5** on this reference is in `IMPLICIT_PROTOTYPE.md` (*The Director's decisions*).
