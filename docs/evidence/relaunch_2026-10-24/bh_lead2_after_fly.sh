#!/bin/bash
# The 2-level lead (lead2s) waits for the fly project's smoke_live.py jobs (the Director, 01:40 24 October: S7 first,
# B4 second, the fly jobs after; the lead ranks below B4). It passed its load gate at 01:33 CDT while the five smokes
# held about 7 CPUs and 89 to 98 threads each. Stopped (SIGSTOP) now, continued (SIGCONT) when no smoke_live.py runs.
# Its wall times then include the stop, logged here. Run detached: setsid nohup bash bh_lead2_after_fly.sh &
LOG=/home/greff/crucible_b4s/runs/lead2/hold.log
DIR=/home/greff/crucible_b4s/runs/lead2
pid=""
for p in $(pgrep -f crucible_b4s/build/crucible_rl10_tier1); do
  [ "$(readlink /proc/$p/cwd)" = "$DIR" ] && pid=$p
done
[ -n "$pid" ] || { echo "$(date '+%F %T %Z') no lead2 process" >> $LOG; exit 1; }
kill -STOP $pid && echo "$(date '+%F %T %Z') stopped pid $pid; smokes: $(pgrep -fc smoke_live.py)" >> $LOG
while pgrep -f smoke_live.py > /dev/null; do sleep 60; done
kill -CONT $pid && echo "$(date '+%F %T %Z') continued pid $pid; load $(cut -d' ' -f1 /proc/loadavg)" >> $LOG
