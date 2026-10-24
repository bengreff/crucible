# Resume after reboot on backhouse

Backhouse was restarted from the Start menu at 17:08 on 23 October 2026, and every WSL job died with it. This path relaunches each registered job from its latest checkpoint when the machine comes back. It was built at the Director's request (23 October, 23:00), with the engine's checkpoint interface agreed with inquiry.

## The parts

- **The Windows task `CrucibleResume`** (`register_task.ps1`).
  - One minute after boot, it runs `wsl.exe -d Ubuntu -u greff --exec /bin/bash /home/greff/resume/resume.sh`.
  - It runs as Ben with logon type S4U: it needs no logon and stores no password.
  - It changes no account, logon or security setting.
  - It retries 3 times, 5 minutes apart, if it fails, and stops after 10 minutes.
- **The runner `/home/greff/resume/resume.sh`** (this directory's `resume.sh`).
  - With no argument, it goes through `/home/greff/resume/jobs/*.job`. It launches each job that is not done, not held and not running.
  - With a name, it does the same for that one job. **This is how a job is first started,** so every registered job runs under its lock.
  - Each job runs as `setsid nohup flock -n DIR/.resume.lock nice -n NICE bash -lc "$RESUME"`, appending to `DIR/resume.log`.
  - **Running** means the lock is held. The kernel drops the lock when the process dies, so after a WSL stop the job counts as not running.
  - The runner's own log is `/home/greff/resume/resume.log`.
- **A job file** (`/home/greff/resume/jobs/<name>.job`) is sourced by bash and sets:
  - `NAME`;
  - `DIR`, the run's directory;
  - `RESUME`, the command run in `DIR`. It must restart from the latest checkpoint, or start fresh if there is none. The first launch and every resume are the same command.
  - `DONE`, the done-file. If it exists, the job is finished and is not relaunched.
  - optionally `NICE` (default 10).
- **Threads and priority:** the runner sets only `nice` and pins nothing. On backhouse, inquiry measured two things on 24 October (`docs/evidence/SPEED.md` section 3):
  - WSL2 shows 28 uniform virtual CPUs, so `taskset` cannot pick out the P-cores.
  - Under the other projects' nice-0 jobs, a 4-thread nice-10 run got 2.56 CPUs, and its wall time stopped improving beyond 4 threads.
  - **The slot policy** for every backhouse job is in SPEED.md section 3: at most 4 threads a job; a launch only while the 1-minute load plus the job's threads stays within 28; the Director's priority order. It belongs in each job's `RESUME` command (a load gate, as `bh_rfx.sh`'s `wait_load`). The runner itself does not gate.
- **To hold a job,** touch `DIR/.resume_hold`; to release it, remove the file. To retire a job, move its job file out of `jobs/`.
  - Killing a registered job without holding it means the next boot relaunches it.

## The engine's side (inquiry's interface, agreed 23 October 2026, about 23:10 CDT; on main at 27b135f, `docs/evidence/RESTART.md`)

- **`CRUCIBLE_CHECKPOINT=<stem>`** turns checkpointing on.
  - The driver writes `<stem>.ckpt` at a step boundary every `CRUCIBLE_CHECKPOINT_EVERY` wall seconds (default 1800).
  - The order is: write `.ckpt.tmp`, fsync, hard-link the old `.ckpt` to `.ckpt.prev`, rename, then fsync the directory. A kill at any point leaves a readable `.ckpt` or `.ckpt.prev`.
- **`CRUCIBLE_RESTART=auto`** reads `<stem>.ckpt`, or `.ckpt.prev` if the latest fails its read or check.
  - With neither, it starts fresh and prints `restart: no checkpoint, fresh start`.
  - A checkpoint records the run's argv, including `argv[0]` (the binary), and its `CRUCIBLE_*` environment, less `CRUCIBLE_CHECKPOINT*`, `_RESTART`, `_KILL_AT_STEP` and `_STATE_OUT`. The run refuses (exit 1) if they differ.
  - Allowed to change across a restart: `CRUCIBLE_T1_THREADS`, mascotte's 6th argument (threads), and `CRUCIBLE_T1_END`, extended only.
  - A rebuilt binary at another path, or a changed environment, refuses the restart loudly; the job then needs a fresh start, by hand.
- **The done-file is the driver's own:**
  - `crucible_rl10_tier1` writes `<CRUCIBLE_T1_OUT>_DEVELOPMENT_done` (default `rl10_t1`) for `run development`, and `rl10_t1_central_done`, a fixed name in the working directory, for `run central`;
  - `crucible_mascotte_a10` writes `<prefix>_done`.
  - A killed process writes none. A run that fails deterministically writes it, so it is not relaunched.
  - A wrapper must not touch a done-file of its own after the run (as `bh_b4.sh`'s `touch done_b4` does), or a killed run would count as done.
- **A tier-1 job, for example** (the values are illustrative; the run's own are its launch's):

  ```sh
  NAME=t1_central
  DIR=/home/greff/crucible_t1c/run
  RESUME='CRUCIBLE_MANIFOLD_TABLE=/home/greff/crucible_t1c/run/rl10_s.bin CRUCIBLE_T1_THREADS=4 CRUCIBLE_CHECKPOINT=central CRUCIBLE_RESTART=auto /home/greff/crucible_t1c/build/crucible_rl10_tier1 run central >> central.log 2>&1'
  DONE=/home/greff/crucible_t1c/run/rl10_t1_central_done
  NICE=10
  ```

## What is tested (23 October 2026, measured)

With the dummy job (`dummy_job.sh`), which counts one step a second and checkpoints with the same file protocol:

- **First launch, then a second pass:** the runner launched it, then skipped it as running (lock held).
- **Killed, then the task run by hand** (`kill -9` of the job and its lock holder at checkpoint 30, then `Start-ScheduledTask CrucibleResume`):
  - The task ran under S4U with result 0.
  - The runner relaunched the dummy, which started at 30.
  - It kept counting after the task's `wsl.exe` had exited.
- **A corrupted latest checkpoint:** after a `kill -9` at 66, `ckpt` was overwritten with garbage. The dummy started at 65, from `ckpt.prev`.
- **The login shell:** killed at 210 and resumed through the task under `bash -lc`, it started at 210.
- **The task runs:** 23:09:39, 23:10:13, 23:13:26 and 23:13:58 CDT, each with result 0. Two runs exposed faults, now fixed:
  - At 23:09:39, a quoted glob meant the runner found no jobs.
  - At 23:13:26, rsync had copied the Mac's non-executable mode, so the job's script was refused. The scripts now carry the exec bit in git.
  - So the runner logs "launched" whether or not the command then starts. A job that fails at once shows its error in `DIR/resume.log`.
- **The refusal check** (`test_wsl_kill.ps1 -DryRun`) refused both modes at 23:13, listing the other workers' processes and not the dummy's (exit 2).
- **To re-arm the dummy for a test:** remove `dummy/done`, `dummy/ckpt` and `dummy/ckpt.prev`, then run `resume.sh dummy`.

**Not yet tested:**

- **A WSL stop under the job** (`test_wsl_kill.ps1`). It was not run tonight because other workers' jobs were running in both distros: the Polytopia self-play in Ubuntu-22.04, and a gate job in Ubuntu.
  - `-Mode terminate` stops Ubuntu only. `-Mode shutdown` is `wsl --shutdown`.
  - The script refuses if any other user process runs in a distro it would stop.
  - `wsl --shutdown` also ends the `PolytopiaWSLKeepalive` task's `sleep infinity`. That task restarts only at Ben's next logon, so a shutdown test needs its owner's agreement.
- **The boot trigger itself,** which needs a Windows reboot.
