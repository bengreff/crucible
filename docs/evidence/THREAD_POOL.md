# Threaded flow step on a persistent pool

Status: done 5 October 2026, about 23:25 (TECHNICAL_PLAN *Lightweight engine*, order of work item 2). The flow step and the reaction call run on one persistent pool. The result is the same bit for bit on any number of threads and equal to the serial code's. The C1 64x12 cold start reaches full thrust (4 ms) in 28 s on 4 threads, against 77 s before, and in 20 s on 6.


Transport added 6 October 2026, about 00:10: the transport properties, gradients, eddy viscosity, SST sources and transport fluxes now run on the same pool. The turbulent C1 (`eqtt`, Table A chemistry with SST) at 64x12 reaches 4 ms in 79 s on 4 threads, against 169 s with serial transport, and in 57 s on 6. The result is the same bit for bit as with serial transport (section *Transport on the pool* below).
## What was built

- **`core/pool.hpp`, `core/pool.cpp`: a persistent pool.**
  - The threads live as long as the pool. The calling thread is worker 0.
  - An idle worker spins for 100 us and then sleeps on a condition variable. 100 us is longer than the serial gaps between the loops of one flow step, so those loops do not each pay a wake-up.
  - A lost wake-up is prevented by the Dekker pattern. A sleeper counts itself before it checks for a job under the mutex, and `run()` bumps the job generation before it reads the sleeper count; both are sequentially consistent.
  - `blocks(n, block, body)` hands out index blocks from a shared counter to whichever worker is free.
  - The first exception thrown in a job is rethrown in the caller after every worker has finished.
  - The pool is not reentrant.
  - It replaces the reaction call's per-step thread start, which cost about 84 us per call (`TABLE_A.md`, criterion 4).
- **`core/flow.cpp`: the flow step on the pool.** `Flow::setThreads(n)` makes the pool, and `ReactingFlow` sets it from its own thread count and runs its reaction loops on it. The right-hand side has four parallel stages and one serial one:
  1. per column (block 1): the axial slopes, the radial profiles and the radial faces, which use only the column's own profiles;
  2. axial faces, in blocks of 4. The dozen supply faces of the injector plate, each with its own root-finder, spread over the workers;
  3. per cell, in blocks of 32: the gather of the cell's face fluxes into its derivative (and, without transport, the division by volume);
  4. serially: the boundary exchange rates (igniter, inlet, outlet, side wall, body force), summed from the stored face fluxes;
  5. with transport: the transport fluxes (still serial), then the division by volume per cell in parallel.

  `refresh`, `stableDt` (a per-worker minimum), the stage updates and the admissibility checks also run per cell on the pool.

## Why the result does not depend on the thread count

- Each face is computed alone into a face array. Each cell then gathers its faces in the serial code's order: igniter, the axial face at i (+), at i + 1 (−), the radial face at j (+, j ≥ 1), at j + 1 (−), then the body force. A cell's sum is therefore the same sequence of floating-point operations as before, whatever thread computes it.
- The boundary rates are summed serially from the face arrays, in the serial code's order.
- The rho k clipped per cell is stored per cell and summed serially. The stable step is a minimum over cells, and a minimum does not depend on the order.
- **Fused multiply-add.** Clang on arm64 contracts `x += a * b` within one expression into one fused multiply-add, which rounds once instead of twice. The first threaded version stored area-weighted fluxes and added them, and its 64x12 energy budget moved from 9.28e-13 to −0.00 (the states differed). The fix:
  - store fluxes per unit area;
  - put each area product inside the accumulating expression, as the serial code had it, so the compiler fuses the same operations;
  - keep as separate statements what the serial code had as separate statements (the species flux `area * massFlux * y[k]`, the body-force work).

## Results (Mac, M2 Pro, 6 performance and 4 efficiency cores; Release; Table A chemistry)

**Bit identity** (`thread_pool/bit_identity.txt`: `cmp` of the raw final state, written by `CRUCIBLE_STATE_DUMP`, and of the history CSV). The old binary is fef9cde plus the state dump.

| Case | Old binary | New, 1 thread | New, 4 threads | New, 6 threads |
|---|---|---|---|---|
| `eqt` 64x12 to 1 ms (inviscid, Table A) | reference | identical | identical | not run |
| `eqt` 64x12 to 8 ms | reference (4 threads, run twice: identical) | identical | identical | identical |
| `frp` 32x6 to 0.25 ms (viscous, SST, PaSR, CVODES; 966 steps, 324 replans) | reference | identical | identical | not run |

The frp runs end with the same 2 criterion-4 failures as the old binary. The 0.25 ms window is short of the settled state those two checks need, and the run was for bit identity only.

**Wall time, C1 64x12 cold start** (`crucible_chamber_study eqt 64 12 8e-3`; 136,228 steps, 768 cells; runs back to back through slot.py, 23:08 to 23:15; load averages 1.8 to 4.9, recorded in `thread_pool/timing_summary.txt`).

| Run | To full thrust (4 ms) | To 8 ms | Per cell update | Speedup on the old 4-thread run | Reaction share |
|---|---|---|---|---|---|
| Old, 4 threads (reaction threaded, flow serial) | 77 s | 153 s | 1.46 us | 1 | 9.0% |
| New, 1 thread | 85 s | 170 s | 1.63 us | 0.90 | 14.1% |
| New, 4 threads | 28 s | 56 s | 0.53 us | 2.73 | 13.9% |
| New, 6 threads | 20 s | 41 s | 0.39 us | 3.73 | 13.4% |

- Against the new code on 1 thread, 4 threads give 3.04 times (76% parallel efficiency), and 6 threads 4.15 times (69%) (derived).
- One thread costs the same as before: 85 s to 4 ms, against Table A's 84 s on 1 thread with the old code (measured on 5 October, `TABLE_A.md`).
- The per-cell-update figures are derived: the wall time divided by steps times cells.
- All four runs end with the same settled values: c* 2447.47 m/s, vacuum thrust 1681.20 N, vacuum Isp 428.59 s. The budgets are mass 0 and energy 9.28e-12 (the same normalisation as `CHAMBER_C1.md`).

**Tests.**
- `core_tests` has a new check, *Threads*: a turbulent chamber with a valve ramp, an igniter and a body force on every cell takes 200 steps on 1 and on 4 threads. `memcmp` finds the states and budgets identical. It passes (`thread_pool/core_tests.txt`).
- The fast ctest suite: 5 of 6 pass in 222 s. `transport_verification` fails on the same three truncation-order criteria with the same values as before (`ctest_main_2026-10-05_thread_pool.txt`; compare `ctest_main_2026-10-05_table_a.txt`).

## Limits and what is still serial

- Transport was serial in the first version; it is threaded since 6 October (next section).
- The serial boundary-rate sums and the body-force loop are O(nz + nr) and O(cells) of trivial work.
- `ReactingFlow::step` computes `Flow::stableDt` and `Flow::step` computes it again: a duplicate per step (not measured).
- Parallel efficiency is 69% to 76% on 768 cells. The parallel stages are short (a 64x12 step takes about 0.4 ms on 6 threads), so the spin and wake-up costs and the serial gaps weigh more than they will on larger meshes (inferred, not measured).
- The pool is not reentrant, and `ReactingFlow` refuses a change of the flow's thread count while it lives (`checkThreads`).

## Transport on the pool (6 October 2026)

**What changed** (`core/transport.cpp`, `core/flow.hpp`, `core/flow.cpp`).
- The per-cell loops (properties, gradients, eddy viscosity, SST sources) run in blocks of 32 cells. Each worker has its own scratch (`TransportScratch`: the work array of `Medium::transport`, the wall's diffusion coefficients, the face mass fractions and species enthalpies), which replaces the shared arrays `transportWork_`, `faceEnthalpy_` and `wallDiffusion_`.
- The wall-face omega is computed per wall cell. `wallCells_` lists each cell next to a no-slip wall once: the side wall's, then the injector plate's not already listed. The serial loop visited the corner cell twice and kept the second value; the second value is the same expression on the same inputs, so listing it once changes nothing.
- The transport fluxes follow the convective pattern. The axial faces (i = 0..nz, with the end faces' kind: plate, supply ring, inlet, outlet) and the radial faces (j = 1..nr, the side wall at j = nr) are computed alone, in blocks of 4, into `axialViscous_` and `radialViscous_` (bulk per unit area, then the species and k, omega fluxes per unit area). Each cell then sums its faces in the serial loops' order: axial face i, i + 1, radial face j, j + 1, then the hoop stress. The boundary rates (wall heat flow, wall axial force) are summed serially in the serial loops' order: the end faces at i = 0, at i = nz, then the side wall by column.

**Bit identity** (`thread_pool/transport_eqtt_summary.txt`: `cmp` of the raw final state and of the history CSV). The old binary is 12613a5 with the `eqtt` mode: flow threaded, transport serial.

| Case | Old, 4 threads | Old, 6 threads | New, 1 thread | New, 4 threads | New, 6 threads |
|---|---|---|---|---|---|
| `eqtt` 64x12 to 4 ms (viscous, SST, Table A; 68,366 steps) | reference | identical | identical | identical | identical |
| `eqt` 64x12 to 8 ms (inviscid; does not use transport) | | | | | identical to the 5 October reference |

**Wall time, turbulent C1 64x12 to 4 ms** (`crucible_chamber_study eqtt 64 12 4e-3`; runs back to back, 5 October 23:40 to 6 October 00:00; `thread_pool/eqtt64_4ms_transport_*.txt`).

| Run | To 4 ms | Per cell update | Speedup on old, 4 threads | Speedup on new, 1 thread | Reaction share |
|---|---|---|---|---|---|
| Old, 4 threads (transport serial) | 169 s | 3.22 us | 1 | 1.38 | 2.6% |
| Old, 6 threads | 172 s | 3.27 us | 0.99 | 1.36 | 2.0% |
| New, 1 thread | 234 s | 4.45 us | 0.72 | 1 | 5.2% |
| New, 4 threads | 79 s | 1.50 us | 2.15 | 2.97 | 5.2% |
| New, 6 threads | 57 s | 1.09 us | 2.96 | 4.09 | 4.7% |

- Parallel efficiency against 1 thread is 74% on 4 threads and 68% on 6 (derived), as for the inviscid step.
- One thread costs 4.45 us per cell update, against the 4.2 to 4.3 us measured on 5 October for the viscous turbulent step (`step_cost_2026-10-05.txt`; a different case, so the comparison is rough).
- All runs end with the same values: c* 2444.94 m/s, vacuum thrust 1682.54 N, vacuum Isp 428.98 s; budgets mass 3.76e-13 and energy -2.25e-11. The first-cell y+ at the wall row (laminar estimate) has median 57.7 and largest 75.9: this mesh's wall cells sit in the log layer, which is what the wall functions are for (`WALL_FUNCTIONS.md`).
- The per-cell-update figures are derived: the wall time divided by steps times cells (52.5 million).

**Tests.**
- `core_tests` passes, including *Threads*, whose chamber is turbulent and so now exercises threaded transport on 1 and 4 threads (`thread_pool/core_tests_transport.txt`).
- The fast ctest suite: 6 of 7 pass, with the new `wall_function_verification` (`ctest_main_2026-10-06_transport_threads.txt`). `transport_verification` fails on the same three truncation-order criteria with the same values as before (axial momentum 1.799, radial 1.223, 4b radial 1.446, limit 1.8); its output is identical to the 5 October record apart from the test numbering.

**Still serial:** the boundary-rate sums (O(nz + nr)) and the source split's bookkeeping between the loops.
