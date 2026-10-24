# Restart from a checkpoint (23 October 2026)

Why: backhouse restarted at 17:08 on 23 October and S7 (low, high) and B4 died with it. A run that cannot survive a reboot starts again from t = 0. The Director's brief of 23:00 asked for a checkpoint write and read, with a restart bit-identical to an uninterrupted run, tested.

## The design

- **One field list per class.** Each class lists its persistent state once, in a template over the stream: `Flow::checkpoint`, `Amr::checkpoint`, `Spray::checkpoint` and `thermo::ReactingFlow::checkpoint`. It is instantiated once with `CheckpointWriter` and once with `CheckpointReader` (`core/checkpoint.hpp`), so the write order and the read order cannot drift apart.
- **What goes in.**
  - The values, and vectors of them, go in as raw bytes.
  - Structs with padding (`Parcel`, the tier-1 driver's `Settle`) go in field by field, so equal states give equal files.
  - Mesh sizes, species counts, the tags and the levels are checked with `same()`, and the reader throws on a mismatch.
  - A named mark before each object's part stops a reader that is out of step at once.
  - After the last field the reader requires the end of the stream.
- **What counts as state.** Anything a later step reads before writing, not just the conserved fields:
  - Flow's temperatures, because they are the Newton warm start.
  - The spray's per-cell gas, including its temperature guess, its snapshot for a redone base step, its random streams, its pending stores and every injector's parcel rate.
  - Amr's clock, counts, regrid state, the budget's reference and the fill's records.
  - The reaction adapter's statistics and table counts.
  - Amr's regions, ownership and parent lists are rebuilt from the active cells at the next step, because `checked_` is false after a read.
- **The file** (`tests/checkpoint_file.hpp`).
  - It holds a header (magic, version, body length, FNV-1a hash), then the body. The body opens with the run's signature: its arguments and `CRUCIBLE_*` environment.
  - The signature leaves out the thread count (`CRUCIBLE_T1_THREADS`, and mascotte's 6th argument), because the engine is bit-identical on any thread count. It also leaves out `CRUCIBLE_T1_END`, so a run may be extended.
  - The write goes to `<stem>.ckpt.tmp`, which is fsynced. The previous file is kept as `<stem>.ckpt.prev` (a hard link), the new one is renamed over `<stem>.ckpt`, and the directory is fsynced. A kill at any point leaves a readable file.
  - `CRUCIBLE_RESTART=auto` reads `.ckpt`, else `.ckpt.prev`, else starts fresh. Every file present but unusable is an error.
- **When checkpoints are written.** Every `CRUCIBLE_CHECKPOINT_EVERY` wall seconds (default 1800) and at the run's end.
  - Mascotte checkpoints after any driver step.
  - Tier 1 checkpoints at a 10 us sample, so the series stays whole.
  - The drivers' output files (mascotte's history and open hot trace, tier 1's series) are cut back to their length at the checkpoint and appended.
  - Each driver writes its own done-file on any exit (`exit N`). A killed process writes none, so no wrapper `touch`.
  - Test hooks: `CRUCIBLE_CHECKPOINT_STEPS=<n>`; `CRUCIBLE_KILL_AT_STEP=<n>`, which calls `_Exit(75)` with no flush and no done-file; and `CRUCIBLE_STATE_OUT=<file>`, which writes the engine objects' fields alone at the end, to compare two runs byte for byte.
- **Names.**

  | Driver | Checkpoint | Done-file |
  |---|---|---|
  | mascotte | `CRUCIBLE_CHECKPOINT=<stem>` | `<prefix>_done` |
  | tier 1 | `CRUCIBLE_CHECKPOINT=<stem>` | `<CRUCIBLE_T1_OUT>_DEVELOPMENT_done` (development) or `rl10_t1_central_done` (central) |

## The tests (Mac, 23 October, about 23:19 to 23:27 CDT; measured)

Each test makes run A uninterrupted. Run B is killed by `CRUCIBLE_KILL_AT_STEP` (exit 75, no done-file) and then resumed with `CRUCIBLE_RESTART=auto` on a different thread count. Both write `CRUCIBLE_STATE_OUT`.

| Case | Kill, restart from | Engine state | Outputs |
|---|---|---|---|
| mascotte `rr:223:3.207:13.7` 60 x 20 to 0.4 ms, windows of 0.1 ms, hot trace at 1998.0 K (open from step 1, still open at the kill) | step 1000, from step 900; 4 then 2 threads | `cmp` identical (607,285 bytes) | the history (all but `wall_s`, `spray_wall_s`), the hot trace (1,557 lines), the hot cells and all nine window files identical |
| the same, hot trace at 1998.25 K (opens at step 1508, after the restart, from the restored ring); `.ckpt` corrupted by 2 bytes | step 700; `.ckpt` refused ("its hash does not match"), from `.ckpt.prev` at step 300 | `cmp` identical | the same, all identical |
| tier 1 development, `CRUCIBLE_T1_LEAD=0` (drops from t = 0), 2 levels, to 0.06 ms | sample 3, from sample 2 (31 base steps, 80 parcels); 4 then 2 threads | `cmp` identical (30,059,080 bytes) | the series (6 rows), the end dumps (`_amr.bin`, `_spray.bin`) and the end line identical |
| tier 1 the same, 3 levels, to 0.03 ms, a checkpoint every sample | sample 2, from sample 1; 4 then 3 threads | `cmp` identical (120,855,140 bytes) | the series identical |

**Limits.**

- A checkpoint is read back by the same binary. It is not a portable format.
- A restart replays nothing: the wall-clock columns (`wall_s`, the spray's share) count from the process's start.
- Mascotte refuses `CRUCIBLE_BUDGET` together with checkpoints.
- The tier-1 stops at the flash (exit 5) and the lead's gate (exit 4) end without a final checkpoint, since they are verdicts.

## The relaunch (backhouse, 23 to 24 October)

S7 low and high and B4 had no checkpoints when they died, so each starts fresh from t = 0, now with checkpoints, under inquiry2's resume runner (`/home/greff/resume`: a job file per run, relaunched at boot by the Windows task CrucibleResume; `tools/resume/README.md` on worktree-inquiry2).

- Tree: this commit's, in `/home/greff/crucible_rst`, built by `restart_2026-10-23/bh_build.sh`.
- Each run: `restart_2026-10-23/bh_run.sh <name>`, from the job files beside it (`s7low.job`, `s7high.job`, `b4.job`, copied to `/home/greff/resume/jobs/`).
  - The script waits while the 1-minute load is above 24, then starts or resumes from its checkpoint.
  - The runs use nice 10 and 4 threads each.
- The inputs are as the dead runs':
  - S7 (`bh_s7z.sh`, 97f4fcd): CFL 0.8, hot trace at 3450 K, 10 ms windows from 15 ms, to 55 ms, table md5 62e730fb.
  - B4 (`bh_b4.sh`, ab3aa03): the 20 ms lead, to 21 ms, development, table md5 1be4531e.
- S7's code is newer than 97f4fcd. The check is the new runs' history rows against the dead runs' (`crucible_s7z/runs2`) over the time both cover.
