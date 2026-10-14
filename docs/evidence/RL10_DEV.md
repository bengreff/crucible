# RL10A-3-3A development chamber: cost and the binding step

Status: two runs on 13 October 2026 (22:26 to 23:06 CDT). **Development only, never the reported case.** Pre-vaporized H2 and O2 at 300 K enter at the face. The real engine's liquid oxygen and its plenum temperatures are not modelled (`RL10A33A_PROTOCOL.md` section 6). The runs use only the protocol's inputs and read no measured performance. Their Isp and thrust are development quantities, printed beside NASA CEA.

## What the runs are for

TECHNICAL_PLAN's implicit-time decision rule: "Design an implicit step only if the measured cold start exceeds 15 minutes on the Mac (three times the target) and the binding step is the wall strip's acoustic limit." The rule is to be applied once Table A, threading, wall functions and AMR all exist (AMR.md criterion 6). These runs measure the chamber without AMR. They give the base mesh's cost and where its step binds, and they apply the rule provisionally.

## Setup

- `crucible_chamber_study eqtt 150 16 25e-3 4 <prefix>` with `CRUCIBLE_CHAMBER=rl10dev`, Table A (`build/table_a.bin`, the 5 October table), wall functions `law`, the `hp` SST corrections, on 4 threads through `slot.py`.
- Geometry is tier 0's: TM-107318's Table E1 contour by PCHIP, the face at −12.97 in (throat station 33 of 150) and the exit at 46.0 in. P1's flow (16.7126 kg/s, O/F 5.055); the igniter delivers 0.5 J over 0.2 ms from 0.2 ms.
- **Run (a): uniform rings.** Built at the `Flow::stepLimit()` commit (author time 22:34). Ran 22:26 to 22:39.
- **Run (b): rings stretched toward the wall,** `CRUCIBLE_RADIAL_STRETCHING=1.7347`. Built from the same tree plus the wall-function y+ accessor and the end-of-run fix (both committed with AMR step A), before step A's geometry edits. Ran 22:45 to 23:06.
- The Mac was shared: other slot jobs ran during both runs (C1 baselines, the ctest baselines, other projects). The 1-minute load average was 10.1 at 23:07. The walls carry that noise, unquantified.

Logs are in [`rl10dev/`](rl10dev/) (`a_uniform_2026-10-13.log`, `b_stretched_2026-10-13.log`), with the settling script `settle.py`. Field and history files stay in the scratch tree.

## Results

| | (a) uniform rings | (b) stretched 1.7347 |
|---|---|---|
| Steps to 25 ms | 137,028 | 277,606 |
| Mean step after 1 ms | 1.951e-7 s | 9.105e-8 s |
| Where the step binds after 1 ms | the throat (column 33), 100%; wall row 0% | the throat's wall row, 100% |
| Binding cell, rates per second | (33, 13): convective 2.04e6, diffusive 2.97e3 | (33, 15): convective 4.24e6, diffusive 1.29e5 |
| Wall time to 25 ms | 752 s | 1,220 s |
| Wall time per cell update (4 threads) | 2.286 us | 1.831 us |
| Wall time per simulated ms | 30.1 s | 48.8 s |
| Vacuum thrust within 0.1% of final | from 7.5 ms, at 195 s of wall | from 7.5 ms, at 353 s |
| Within 0.01% | from 10.0 ms, at 262 s | from 10.0 ms, at 472 s |
| Drift over the last 1 ms (thrust) | 7.38e-4 (the end-of-run artefact below) | 9.02e-11 |
| Budgets (mass, energy) | −1.01e-13, −1.22e-13 | −2.00e-13, −2.98e-12 |
| States outside Table A | 15 (above 3,500 K, direct call) | 0 |
| Vacuum thrust, settled | 75,554.2 N | 75,360.7 N |
| c*, vacuum Isp (final row) | 2,419.86 m/s, 461.33 s (thrust row carries the artefact) | 2,415.23 m/s, 459.81 s |
| Wall heat flow into the gas, axial force on the gas | −9.366 MW, −101.07 kN | −9.130 MW, −100.76 kN |
| First-cell y+ at the wall row, laminar estimate: median, largest | 109.6, 172.3 | 57.0, 95.5 |
| First-cell y+, wall functions: median, largest | not printed (the accessor came later) | 623.7, 3,716 |

Settling is read from the progress rows, which are 0.5 ms apart. Run (b)'s CEA column at its own p0 (3.2643 MPa, gases at 300 K) is: Tc 3,361.6 K, c* 2,424.57 m/s, Isp_vac 478.90 s shifting and 454.85 s frozen. Simulated over CEA: c* 0.9961, Isp 0.9601 of shifting.

### What the numbers say

- **Cost on the base mesh.** Settled thrust takes 3.3 minutes of wall time on uniform rings and 5.9 minutes on stretched rings, both on 4 threads of a loaded Mac. Both are under the rule's 15 minutes.
- **The binding step moves to the wall when the rings are stretched.** On uniform rings, the step is set at the throat two rings in from the wall. Stretching halves the step (9.1e-8 s against 2.0e-7 s) and moves the binding cell into the throat's wall row. Its limit there is convective (|u| + c over the thin ring), 33 times its diffusive rate. That is the condition in the rule's second clause.
- **The radial grid moves thrust by 0.26%.** From (a) to (b), vacuum thrust falls 0.26% (settled values), c* 0.19% and wall heat flow 2.5%. 150 x 16 is not a converged grid. Its throat is also unresolved: the PCHIP contour's throat radius of curvature is 0.137 r_t = 8.6 mm, below one 10.0 mm column. So the c* gap to the 2-D prediction (CEA c* over Kliegel-Levine's Cd at R/r_t 0.5: simulated −2.5%) is not a model check. AMR's throat tag exists for this.
- **The wall cells are far above the wall functions' target.** Run (b)'s first cells sit at y+ 624 (median) and 3,716 (largest) in wall-function units. AMR.md's default target is 300, and the law alone passes its criterion 1(b) up to 956. The laminar estimate (57) is 11 times smaller, the same wall-unit difference C1 shows (52 against 1,198), because the wall function evaluates rho_w and mu_w at the cold wall.

### The end-of-run artefact (run a) and its fix

Run (a)'s last recorded thrust is 75,610.04 N, 7.4e-4 above the 75,554.2 N of every progress row from 24.5 ms. The march's accumulated sample time fell an ulp short of 25 ms. A step of about 1e-18 s then followed, and its record holds the boundary rates of its first stage alone. `chamber_study`'s march now treats a target within 1e-12 s of the end as the end (committed with AMR step A). Run (b), with the fix, drifts 9.0e-11 over its last ms.

## The decision rule, provisionally

- **On this mesh the rule says no implicit step.** The cold start is 5.9 minutes, under 15. The second clause holds in run (b), where the throat's wall row binds.
- **With AMR the first clause is likely to hold as well** (derived; AMR.md criterion 6 replaces this with a measurement). Arithmetic:
  - Under Berger-Colella subcycling, level l costs its cells times 2^l per base step.
  - Take each refined level of the wall strip as 6 rows of that level: the first row tagged, plus the 2-cell buffer, doubled.
  - The base step stays at run (b)'s 9.15e-8 s, since each level halves both the cell and the step.
  - The cost per cell update is run (b)'s 1.83 us. The transient is unchanged, settling by 7.5 ms.
  - **Target y+ 300** (AMR.md's default): 2 levels along the whole wall for the median 624, 4 levels over the throat's 7 base columns for the largest 3,716, and the static throat tag to level 2. That is about 41,900 cell updates per base step, 17.5 times the base mesh's 2,400. It comes to about 840 s per simulated ms: **1.7 hours to settled thrust**, 2.2 hours with AMR.md's 25% overhead allowance.
  - **Target y+ 956** (the law's own limit): 2 levels over the throat columns only. With the static throat tag that is about 11,300 updates per base step: 28 to 35 minutes. Without the throat tag it is about 3,240: 8 minutes.
- **The target y+ decides the rule.** At 300 both clauses would hold and step 16's dual-time design would be due. Near 1,000 the explicit march stays under 15 minutes. The target is an accuracy question for the wall functions (`WALL_FUNCTIONS.md`, restated criterion 2). It should be settled by measurement on this chamber in AMR step D (wall heat flow and thrust at targets 300 and 956), before criterion 6 is run.

## Not covered

- Liquid oxygen, injection and spray (the reported case needs them; the inquiry worker's lane).
- Grid convergence: two radial grids at one axial resolution only.
- The per-column y+ distribution is not printed, only its median and largest. The AMR projection above assumes where it exceeds each target.
