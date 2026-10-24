#!/bin/bash
# Job rfx's D1 pair under the backhouse slot policy (SPEED.md section 3: a job launches only while the 1-minute load
# plus its own threads stays within 28). bh_rfx.sh starts D1 new and parent together behind one load-24 gate, so this
# watcher stops the parent run (SIGSTOP) as soon as it appears, and continues it once the load, with the new run's
# 4 threads in it, is under 24. Log: out/d1gate.log. Ends when the parent is continued, or when rfx's all_done appears.
O=/home/greff/crucible_rfx/out
echo "$(date) watching" >> $O/d1gate.log
while :; do
  [ -f $O/all_done ] && { echo "$(date) all_done before the parent started" >> $O/d1gate.log; exit 0; }
  p=$(pgrep -o -f "^/home/greff/crucible_inq2t1/build/crucible_amr_budget 8e-3|^nice -n 10 /home/greff/crucible_inq2t1/build/crucible_amr_budget 8e-3") && break
  sleep 5
done
kill -STOP $p && echo "$(date) stopped parent D1 pid $p" >> $O/d1gate.log
sleep 180
while [ "$(cut -d. -f1 /proc/loadavg)" -ge 24 ]; do sleep 60; done
kill -CONT $p && echo "$(date) continued parent D1 pid $p at load $(cut -d' ' -f1 /proc/loadavg)" >> $O/d1gate.log
