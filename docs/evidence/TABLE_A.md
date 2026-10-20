# Table A: tabulated shifting equilibrium

Status: criteria stated 5 October 2026 before any table code (committed 21:58 as db994a4; this line said "about 22:10" until the results were added). Results are appended below the line; they do not change the criteria. All four criteria pass (5 October, 22:40); the C1 energy budget on 32x6 misses CHAMBER_C1's 1e-11 by the normalisation recorded there, as Cantera's run does on the same grid.

## What is tested

Table A replaces the per-cell Cantera constant-(u, v) equilibrium of `Chemistry::LocalEquilibrium` with a table built offline by the same Cantera call (TECHNICAL_PLAN, *Lightweight engine*). The model is unchanged, so the reference is the direct Cantera call at the same state. The table is built for the C1 mechanism, `h2o2.yaml`.

- **Axes.**
  - f_H = Z_H / (Z_H + Z_O), with a node at the stoichiometric f of H2O. The complete-combustion products bend there.
  - Z_N, the N2 mass fraction.
  - eta = (e − e_lo) / (e_hi − e_lo), normalised per (f_H, Z_N) node. e_lo is the energy of the complete-combustion products at 200 K; e_hi is the energy of the same products, frozen, at 6000 K.
  - ln rho.
- **Values.** Each node stores the equilibrium Y_k and T, computed by Cantera's equilibrateUV started from the complete-combustion products at the node's energy.
- **Lookup.** Multilinear interpolation, with each (f_H, Z_N) corner using its own eta.
  - Z_H, Z_O and Z_N are bilinear in (f_H, Z_N), and multilinear weights reproduce bilinear functions exactly. So the returned Y carries the query's elements to rounding, and sum Y = 1.
  - The flow keeps its total energy. The cell's T comes from the core's e-to-T inversion with the new Y, so energy is exact. The tabulated T is only a consistency check.
- **Outside the table.**
  - e below e_lo: eta is clamped to 0, which returns the complete-combustion products, the low-temperature limit of equilibrium. Counted.
  - e above e_hi, rho outside the axis, or a species not made of H, O and N2 above 1e-12 (AR): the direct Cantera call is used. Counted.

## Criteria

1. **Conservation.**
   - At every node, the equilibrium's element mass fractions equal the node's to 1e-12, and |sum Y − 1| ≤ 1e-12.
   - For every lookup in criterion 2, the returned Y has the query's element mass fractions to 1e-12, and |sum Y − 1| ≤ 1e-12.
   - In the C1 runs of criterion 3, the mass and energy budgets are reported against CHAMBER_C1's 1e-11 limits, with the same normalisation caveat.
2. **Accuracy against direct Cantera** (|dT| is the core's T from the table's Y, against Cantera's equilibrium T at the same elements, e and rho).
   - (a) Random states: 1e5 states uniform in (f_H, Z_N, eta, ln rho) over the table's range. Report the median, 99th percentile and maximum of |dT|, and the maximum |dY_k| per species. Pass: 99th percentile of |dT| ≤ 3 K for states whose equilibrium T is between 200 and 4000 K.
   - (b) The run's own states: in the C1 runs of criterion 3, an audit compares every cell on every Nth reaction call with the direct Cantera call from the same pre-reaction state (N declared per run). Pass: maximum |dT| ≤ 3 K over audited cells inside the table. Report the maximum |dY_k|, and the clamp and fallback counts.
   - Why 3 K: c* goes as sqrt(T / M). A uniform 3 K error at 3350 K shifts c* by 0.045%, inside criterion 3's 0.05% (derived).
3. **The C1 cold start.** The table (`crucible_chamber_study eqt`) is compared with direct Cantera (`eq`) on 32x6 and 64x12 to 8 ms, same tree, same machine and same threads (5).
   - Pass: settled c*, vacuum Isp and vacuum thrust agree within 0.05%. 0.05% is one eighth of the 64x12 grid error in c* (+0.40%, CHAMBER_C1), so the table's error is negligible against the discretisation's.
   - Pass: injector pressure and outlet mass flow agree within 0.1% at every printed time from 0.5 ms.
   - Pass: CHAMBER_C1's settling criterion (drift below 1e-3 over the last 1 ms).
4. **Cost.** Wall time to 8 ms on the Mac for `eq` and `eqt`, same tree and threads. `eqt` is also run on 1 thread.
   - Pass: speedup of 10 or more on 64x12.
   - Pass: the 64x12 cold start to full thrust in 5 minutes or less. Full thrust means 4 ms: the recorded 64x12 vacuum thrust is within 0.1% of its settled value from 3.5 ms (`c1/eq64_log.txt`).
   - Reported: the table's own cost per cell and its share of the step.

The table's resolution (nodes per axis) is the setting that criteria 2 and 3 select; it is reported with the result. The table file is built by `crucible_build_equilibrium_table`, is not kept in git (it is reproducible), and its provenance (mechanism, axes, Cantera version, build time, node conservation error) is printed by the builder and copied into the results below.

---

## Results

Runs on this Mac (M2 Pro), 5 October 2026, 21:48 to 22:37, built from the tree committed with this record. Two single-thread jobs of another project ran from 22:18 on (1-minute load average 3.5 to 5.8 at the starts and ends of the table runs; not recorded for Cantera's 64x12 run), so the walls carry that much noise (`table_a/load_2026-10-05.txt`). Logs are in [`table_a/`](table_a/). The runs are `crucible_chamber_study eq|eqt <nz> <nr> 8e-3 <threads> <prefix>`, with `CRUCIBLE_EQ_TABLE=build/table_a.bin` and, for the audited runs, `CRUCIBLE_EQ_AUDIT=<interval>`. The fast ctest suite afterwards: 5 of 6 pass, with the known `transport_verification` failures and values unchanged (`../ctest_main_2026-10-05_table_a.txt`).

### The table as built

`crucible_build_equilibrium_table build/table_a.bin 20 64 1e-4 129 97 1e-8 30 4`, provenance as printed:

> Table A built 2026-10-05 22:18 CDT by crucible_build_equilibrium_table from h2o2.yaml, Cantera 3.2.0; f 131 nodes (20 lean, 64 rich uniform, node at f_st 0.111907, 46 clustered at it and the ends from spacing 0.0001), T 129 over [300, 3500] K, ln rho_r 97 over [1e-08, 30] kg/m3; 2 inert species; 1639203 nodes in 11 s on 4 threads; node element error 1.93e-10 from Cantera, 3.46e-16 after the projection, |sum Y - 1| 4.44e-16, failures 0

The file is 144 MB.

### How the build differs from *What is tested*, and why

The criteria are unchanged. The design they were written against failed criterion 2(a) by about 100 times, so the axes and the energy range changed.

1. **Energy range: Cantera's thermo range, 300 to 3,500 K.** Outside the NASA-7 range of `h2o2.yaml`, Cantera's equilibrium call took about 70 ms instead of 30 to 60 us (measured: its fast solver gives up and a slower one takes over). The polynomials are also extrapolated there.
   - States whose energy lies above the 3,500 K equilibrium's are *Outside* and take the direct call.
   - States below the 300 K equilibrium's are *Clamped*: the 300 K equilibrium is returned, which is complete combustion.
   - The check therefore samples eta between the products at 300 K and the products frozen at 5,000 K, not 200 K and 6,000 K.
2. **Axes: f_H, T and ln rho_r, not f_H, Z_N, eta and ln rho.** On a small trial table (f 5 lean and 8 rich nodes), the planned axes gave a 99th-percentile |dT| of 240.7 K (median 10.8 K; 4,000 samples), with the worst states near pure N2 (Z_N 0.978). Two variants did no better:
   - each corner read at the query's own eta: 298.8 K;
   - a T axis in place of eta, keeping Z_N: 179 K, worst at Z_N 0.984.
   
   Mechanism (derived): the H/O equilibrium depends on the reacting species' own partial density, rho (1 − Z_N). At fixed rho this changes by decades as Z_N approaches 1, so it is far from linear in Z_N. Since N2 and Ar are each the only carrier of their element, the diluent drops out exactly when the axis is rho_r = rho (Z_H + Z_O). The same trial grid then gave 99 K, and the rest of the error sat at the stoichiometric kink.
3. **The f axis is clustered** geometrically (ratio 1.5) at f_st and at both ends, from a spacing of 1e-4.
4. **Nodes are projected onto exact elements.** Cantera's equilibrium carries the node's elements only to 1.93e-10, short of criterion 1's 1e-12. A 2x2 solve on the H and O carriers, done twice, brings this to 3.46e-16.

Resolution study (4,000 samples, 99th-percentile |dT| over T_eq 200 to 4,000 K):

| f nodes (lean + rich, clustering) | T | ln rho_r | 99th (K) | Median (K) |
|---|---|---|---|---|
| 10 + 32, none | 129 | 49 | 11.7 | 0.305 |
| 10 + 32, none | 257 | 49 | 11.5 | 0.251 |
| 10 + 32, none | 129 | 97 | 11.4 | 0.121 |
| 20 + 64, none | 129 | 49 | 3.6 | |
| 20 + 64, at f_st from 1e-4 | 129 | 49 | 2.34 | |
| 20 + 64, at f_st and the ends from 1e-4 | 129 | 49 | 2.005 | |
| 32 + 96, at f_st and the ends from 1e-5 | 129 | 49 | 1.875 | |

The f resolution sets the error; T and rho_r beyond 129 and 49 buy little. The chosen table adds rho_r nodes for the median: 97 nodes give 0.112 K where 49 gave 0.3.

### Criterion 1: conservation. Pass

- **Nodes:** element error 3.46e-16, |sum Y − 1| 4.44e-16 (limit 1e-12).
- **Lookups,** the 89,570 compared states of criterion 2(a): element error 4.44e-16, |sum Y − 1| 5.83e-16.
- **C1 budgets** (limit 1e-11; CHAMBER_C1's normalisation divides by the initial fill's energy, about −0.25 J):

| Run | Mass | Energy |
|---|---|---|
| Table, 32x6 | 2.02e-12 | −2.20e-11 |
| Cantera, 32x6 | −4.05e-13 | −8.35e-12 |
| Table, 64x12 | 0.00e+00 | 9.28e-12 |
| Cantera, 64x12 | −3.76e-13 | −1.53e-11 |

  The table's 32x6 energy budget and Cantera's 64x12 one exceed 1e-11 by the normalisation that CHAMBER_C1 records, not by a leak. 2.2e-11 of 0.25 J is about 5e-12 J.

### Criterion 2(a): random states. Pass

`crucible_build_equilibrium_table check build/table_a.bin 100000 3` ([log](table_a/check_2026-10-05.txt)). It samples 1e5 states uniform in f_H and Z_N over [0, 1], in eta, and in ln rho over [1e-4, 30] kg/m^3, from seed 1.
- 10,430 states are *Outside* (the run makes the direct call for these), none are clamped, and Cantera failed on none.
- Over the 89,570 states with T_eq from 200 to 4,000 K, |dT| has a median of 0.112 K, a 99th percentile of **1.063 K** (limit 3 K) and a maximum of 5.439 K.

| T_eq (K) | States | Median (K) | 99th (K) | Max (K) |
|---|---|---|---|---|
| 200–1,000 | 11,258 | 0.000 | 0.000 | 0.000 |
| 1,000–2,000 | 19,777 | 0.001 | 0.244 | 2.046 |
| 2,000–2,500 | 20,015 | 0.190 | 1.173 | 4.473 |
| 2,500–3,000 | 23,877 | 0.267 | 1.278 | 5.139 |
| 3,000–3,500 | 14,639 | 0.280 | 1.584 | 5.439 |
| 3,500–4,000 | 4 | 2.003 | 2.332 | 4.539 |

- Largest |dY_k|:

| H2 | H | O | O2 | OH | H2O | HO2 | H2O2 | AR | N2 |
|---|---|---|---|---|---|---|---|---|---|
| 2.44e-4 | 5.26e-5 | 5.24e-4 | 1.61e-3 | 8.05e-4 | 1.56e-3 | 8.75e-6 | 1.31e-6 | 0 | 2.34e-10 |

- The worst state is slightly rich of stoichiometric (f_st 0.112), at a high density: f_H 0.134, Z_N 0.247, rho 8.65 kg/m^3, T_eq 3,448.6 K.
- The maximum is above 3 K, which matters for criterion 2(b). There the run's own states stayed below 0.85 K.

### Criterion 2(b): the run's own states. Pass

| Run | Audit interval (reaction calls) | Audited cells | Max \|dT\| (K) | Max \|dY\| | Clamped | Outside |
|---|---|---|---|---|---|---|
| 32x6 | 200 | 65,472 | 0.844 | 2.11e-4 | 48,780 of 13,064,064 | 0 |
| 64x12 | 1,000 | 105,216 | 0.826 | 2.11e-4 | 443,005 of 104,623,104 | 0 |

- The clamped cells (0.4%) are the ambient N2 fill (no H/O, or rho_r below 1e-8 kg/m^3) and gas colder than 300 K. The two are not counted separately. Both are inside the audit.
- The audit keeps the table's result, so the audited runs equal the unaudited ones to the printed digit (64x12: identical end state and walls).

### Criterion 3: the C1 cold start against Cantera. Pass

Same tree, same Mac, 5 threads, to 8 ms. Cantera's references were rerun on this tree: [32x6](table_a/eq32_2026-10-05.txt) and [64x12](table_a/eq64_2026-10-05.txt). They equal the recorded runs (CHAMBER_C1) to every printed digit except the 64x12 budgets, which differ at round-off (mass −3.76e-13 against −1.36e-12, energy −1.53e-11 against −1.16e-11), as the 32x6 rerun after the axis clamp did.

| | 32x6 Cantera | 32x6 table | Difference | 64x12 Cantera | 64x12 table | Difference |
|---|---|---|---|---|---|---|
| c* (m/s) | 2,478.55 | 2,478.40 | −0.0061% | 2,447.60 | 2,447.47 | −0.0053% |
| Vacuum Isp (s) | 427.99 | 427.98 | −0.0023% | 428.59 | 428.59 | 0 at the printed digit |
| Vacuum thrust (N) | 1,678.85 | 1,678.82 | −0.0018% | 1,681.21 | 1,681.20 | −0.0006% |
| Injector face (MPa) | 3.31817 | 3.31797 | −0.0060% | 3.27798 | 3.27780 | −0.0055% |
| Drift, last 1 ms (outlet, p_inj, F_vac) | 7.8e-8, 6.3e-8, 6.4e-8 | 7.9e-8, 6.3e-8, 6.5e-8 | | 6.5e-8, 5.2e-8, 5.3e-8 | 6.5e-8, 5.3e-8, 5.4e-8 | |

- **Settled values:** every difference is within 0.05%, with margin of about 8 times or more.
  - In c*, the table runs about 0.005% low. That corresponds to about 0.4 K in chamber temperature (derived from c* proportional to sqrt(T)).
  - The table does not change the grid's error: c* against the corrected ideal rocket is +0.3956% on 64x12, against Cantera's +0.4010%.
- **History** (`tools/compare_histories.py`, from the 2 us histories at each printed time from 0.5 to 8 ms): pass. The largest differences are:

| Grid | Injector pressure | Outlet mass flow | Vacuum thrust |
|---|---|---|---|
| 32x6 | 1.15e-4 | 2.39e-5 | 8.7e-5 |
| 64x12 | 1.35e-4 | 2.06e-5 | 1.59e-4 |

  The largest fall at 0.5 ms, during light-off. They are within the 1e-3 limit by about 7 times or more ([32x6](table_a/compare_32x6_2026-10-05.txt), [64x12](table_a/compare_64x12_2026-10-05.txt)).
- **Settling:** the drift over the last 1 ms is below 1e-7 on both grids (limit 1e-3).

### Criterion 4: cost. Pass

| Run | Wall to 4 ms | Wall to 8 ms | Reaction share | Reaction per cell update |
|---|---|---|---|---|
| Cantera, 32x6, 5 threads | 135 s | 248 s (197 s recorded) | 91.0% | 17.3 us |
| Table, 32x6, 5 threads | 15 s | 29 s | 19.8% | 0.44 us |
| Cantera, 64x12, 5 threads | 1,089 s | 2,048 s (2,375 s recorded) | | |
| Table, 64x12, 5 threads | 82 s | 163 s | 10.0% | 0.156 us |
| Table, 64x12, 1 thread | 84 s | 168 s | 14.5% | 0.234 us |

- **Speedup on 64x12, 5 threads: 12.6 times** to 8 ms, and 13.3 times to 4 ms (limit 10). On 32x6: 8.6 times against tonight's Cantera run and 6.8 times against the recorded one. The coarse grid gains less because the reaction call has a fixed cost. On 32x6 it takes about 84 us per call (5.7 s over 68,042 calls), of which the lookups are about 9 us (192 cells at the one-thread 0.23 us, over 5 threads). Most of the rest is starting and joining the worker threads on every call (derived, not profiled). The reaction is 19.8% of the table run on 32x6, against 10.0% on 64x12.
- **Cold start to full thrust (4 ms) on 64x12: 82 s** on 5 threads and 84 s on 1 thread (limit 5 minutes).
- **What a step costs now.** The flow step is the step: 1.37 us per cell update on one thread (derived: 168 s less the 24.4 s reaction, over 104,623,104 cell updates). That agrees with the 1.35 to 1.39 us measured for `step_cost`. TECHNICAL_PLAN's budget derived 78 s to 4 ms on one thread; the measurement is 84 s.
- The lookup costs 0.23 us per cell on one thread. Threading it gains only 1.5 times (0.156 us at 5 threads), because the reaction call starts its worker threads afresh each step for 768 cells. Threading the flow step (TECHNICAL_PLAN, order of work item 2) should use a persistent pool for both.
- Cantera's 64x12 reference ran 14% faster tonight than on 4 October (2,048 s against 2,375 s), so its wall time carries about that much noise. The speedup is against tonight's run on the same tree. The other project's two jobs began only in that run's last 4 minutes, so it had less contention than the table runs, and the speedup is if anything understated.

### Not covered

- Table A is the equilibrium model only. It burns cold premixed gas at once, as `LocalEquilibrium` does, so ignition delay and kinetic lag need Table B.
- The table covers `h2o2.yaml` with inert N2 and Ar. A mechanism in which N reacts is refused by the builder.
- States above 3,500 K or above rho_r 30 kg/m^3 take the direct call (none in C1). An RL10 chamber at about 3 MPa sits well inside the rho_r axis. Its hottest states should be checked against the 3,500 K edge before C4.

## To 6,000 K (20 October 2026, about 17:35 CDT; stated before the build)

Why: Table B's energy axis reaches past 3,500 K by design, and Table A's T axis was the mechanism's thermo range, so both stop where `h2o2.yaml`'s fits do (TABLE_B.md, *The 3,500 K ceiling*). `tools/mechanisms/h2o2_highT.yaml` carries NASA TM-4513's fits to 6,000 K at their 1 bar standard state, with OH's heat of formation from Ruscic et al. (2002); its HP equilibrium matches CEA's to 1.79 K over 3,000 to 5,000 K (TABLE_B.md, *OH's heat of formation, the reference pressure, and the equilibrium check*).

- **What changes in the code.** `crucible_build_equilibrium_table` takes the mechanism as a last, optional argument (default `h2o2.yaml`) and prints it in the provenance; its `check` takes it after the seed and refuses a table whose provenance names another. The check samples energies up to the frozen products at 5,000 K or the table's top T, whichever is higher. `crucible_chamber_study` reads the mechanism from `CRUCIBLE_MECHANISM` (default `h2o2.yaml`) for the core thermo, transport and reaction, and refuses a Table A built on another.
- **The build.** The same f and ρ_r axes as on 5 October (`20 64 1e-4`, 97 ρ_r nodes over 1e-8 to 30 kg/m³). The T axis is h2o2_highT's range, 200 to 6,000 K, at the same 25 K spacing: 233 nodes.
- **Criteria, unchanged, plus one extension.**
  - Criterion 1 as stated.
  - Criterion 2(a) as stated (99th percentile of |ΔT| ≤ 3 K on T_eq 200 to 4,000 K, 1e5 states, seed 1, against direct Cantera on h2o2_highT). **Extension:** the same bar on T_eq 4,000 to 6,000 K, judged separately.
  - **Nothing moved below 3,500 K (C1, 32x6, cheap):** criterion 3 restated on h2o2_highT, the table run (`eqt`) against direct Cantera (`eq`) on the same mechanism, 32x6 to 8 ms on 5 threads: settled c*, vacuum Isp, vacuum thrust and injector pressure within 0.05%, histories within 1e-3 (`compare_histories.py`). Criterion 2(b) on the same run (audit every 200 reaction calls): max |ΔT| ≤ 3 K. Beside, not judged: today's code with the 5 October table on `h2o2.yaml`, to show what the mechanism change alone moves.

### Results (20 October 2026, about 17:50 CDT): all pass

- **Build:** 2,960,731 nodes in 30 s on 4 threads. Criterion 1 passes at every node: element error 3.64e-16 after the projection, |Σ Y − 1| 4.44e-16, no failures ([check](table_a/check_highT_2026-10-20.txt), provenance on its first line).
- **Criterion 2(a), 1e5 samples:** 0 outside, 1 clamped, 0 Cantera failures. Lookups pass criterion 1 (4.44e-16, 5.13e-16).

  | T_eq band (K) | States | Median \|ΔT\| (K) | 99th (K) | Max (K) |
  |---|---|---|---|---|
  | 200–1,000 | 10,211 | 0.000 | 0.000 | 0.000 |
  | 1,000–2,000 | 15,874 | 0.001 | 0.248 | 1.430 |
  | 2,000–2,500 | 16,246 | 0.191 | 1.154 | 4.319 |
  | 2,500–3,000 | 24,101 | 0.279 | 1.109 | 5.368 |
  | 3,000–3,500 | 16,016 | 0.310 | 1.416 | 5.351 |
  | 3,500–4,000 | 10,058 | 0.354 | 1.683 | 5.522 |
  | 4,000–5,000 | 7,261 | 0.427 | 1.545 | 4.339 |
  | 5,000–6,000 | 233 | 0.009 | 0.703 | 0.732 |

  - **200 to 4,000 K (as stated):** 92,506 states, **99th percentile 1.086 K** (5 October's table: 1.063 K), max 5.522 K at f 0.1315, ρ_r 26.7 kg/m³, T_eq 3,564 K: **pass**.
  - **4,000 to 6,000 K (the extension):** 7,494 states, **99th percentile 1.515 K**, max 4.339 K: **pass**.
  - Few samples land above 5,000 K: the check draws energies up to the frozen products at 6,000 K, and equilibrium dissociation takes most of them below 5,000 K.
  - Largest |ΔY|: O2 1.67e-3, H2O 1.60e-3, OH 1.04e-3.
- **C1, 32x6 to 8 ms on h2o2_highT (criterion 3 restated): pass.** The table run ([eqt](table_a/eqt32_highT_2026-10-20.txt)) against direct Cantera ([eq](table_a/eq32_highT_2026-10-20.txt)):

  | | Cantera | Table | Table/Cantera − 1 |
  |---|---|---|---|
  | c* (m/s) | 2,477.89 | 2,477.73 | −0.0065% |
  | Vacuum thrust (N) | 1,678.61 | 1,678.58 | −0.0018% |
  | Vacuum Isp (s) | 427.93 | 427.92 | −0.0023% |
  | Injector pressure (MPa) | 3.31724 | 3.31702 | −0.0066% |

  - Histories: largest relative differences 1.17e-4 (injector pressure), 2.36e-5 and 1.05e-4 against the 1e-3 limit ([compare](table_a/compare_32x6_highT_2026-10-20.txt)). 5 October's table on `h2o2.yaml` gave 1.15e-4, 2.39e-5 and 8.7e-5.
  - Criterion 2(b), audit every 200 reaction calls: 65,472 cells, max |ΔT| 0.825 K, max |ΔY| 2.05e-4 (5 October: 0.844 K): **pass**.
  - Table use: 13,063,296 cell updates, 48,788 clamped at an axis end, none outside. 4.2 s of the 47 s run against Cantera's 225 s of 327 s, on a shared Mac (2 to 4 other jobs), so the times are indicative only.
- **What the mechanism change alone moves (beside, not judged).** The 5 October table on `h2o2.yaml`, run on today's code ([log](table_a/eqt32_old_2026-10-20.txt)), reproduces 5 October's run to within 3e-5 relative (c* 2,478.43 against 2,478.40 m/s; the code has moved since). Against it, the h2o2_highT table run is lower by 0.028% in c*, 0.013% in thrust, 0.012% in Isp and 0.030% in injector pressure. The cause is the thermo, not the table: the same chamber's 1-D ideal rocket gives Tc 3,348.1 K on h2o2_highT and 3,351.7 K on `h2o2.yaml`, and CEA (RocketCEA, the same reactants at 300 K, O/F 5 and p0 3.155 MPa) gives 3,347.7 K and c* 2,428.61 m/s, against 2,428.76 (h2o2_highT) and 2,429.41 m/s (`h2o2.yaml`). So the change moves C1 toward CEA.
- **A label fixed.** The chamber study's line "1-D ideal rocket at p0 (h2o2.yaml thermo)" was hard-coded and printed on the h2o2_highT runs above, though the numbers are on the run's mechanism. It now prints the mechanism; the runs above were made before the fix.
