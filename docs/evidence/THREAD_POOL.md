# Threaded flow step on a persistent pool

Status: done 5 October 2026, about 23:25 (TECHNICAL_PLAN *Lightweight engine*, order of work item 2). The flow step and the reaction call run on one persistent pool. The result is the same bit for bit on any number of threads and equal to the serial code's. The C1 64x12 cold start reaches full thrust (4 ms) in 28 s on 4 threads, against 77 s before, and in 20 s on 6.

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

- **Transport (`core/transport.cpp`) is serial:** the transport properties, gradients, eddy viscosity, SST sources and transport fluxes. The inviscid C1 above does not use it. The viscous and turbulent step (4.2 to 4.3 us per cell update on one thread, measured, `step_cost_2026-10-05.txt`) is mostly transport, so threading it is the next lever for the turbulent cases. The same face-array and gather pattern applies.
- The serial boundary-rate sums and the body-force loop are O(nz + nr) and O(cells) of trivial work.
- `ReactingFlow::step` computes `Flow::stableDt` and `Flow::step` computes it again: a duplicate per step (not measured).
- Parallel efficiency is 69% to 76% on 768 cells. The parallel stages are short (a 64x12 step takes about 0.4 ms on 6 threads), so the spin and wake-up costs and the serial gaps weigh more than they will on larger meshes (inferred, not measured).
- The pool is not reentrant, and `ReactingFlow` refuses a change of the flow's thread count while it lives (`checkThreads`).
