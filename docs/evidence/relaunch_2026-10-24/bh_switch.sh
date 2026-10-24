#!/bin/bash
# Retires the old runs (b4r on c910fca, lead2 on 27b135f): hold, stop the engine (a killed run writes no done-file),
# and move the job files to /home/greff/resume/retired; their directories and checkpoints are kept. Then registers
# b4s and lead2s and starts them through the runner, so each runs under its lock.
set -e
J=/home/greff/resume/jobs
mkdir -p /home/greff/resume/retired
for old in b4r lead2; do
  ( source $J/$old.job; touch $DIR/.resume_hold
    for p in $(pgrep -f crucible_rl10_tier1); do
      if [ "$(readlink /proc/$p/cwd)" = "$DIR" ]; then echo "stopping $old: pid $p"; kill -TERM $p; fi
    done )
  mv $J/$old.job /home/greff/resume/retired/
done
cp /home/greff/crucible_b4s/docs/evidence/relaunch_2026-10-24/{b4s,lead2s}.job $J/
chmod +x /home/greff/crucible_b4s/docs/evidence/relaunch_2026-10-24/*.sh
for new in b4s lead2s; do bash /home/greff/resume/resume.sh $new; done
date '+switched %F %T %Z'
