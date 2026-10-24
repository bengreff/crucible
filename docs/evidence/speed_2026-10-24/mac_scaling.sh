#!/bin/bash
# Mac thread scaling beside backhouse's: tier 1 development (B4's config), restarted from the lead profile's 0.15 ms
# checkpoint (/tmp/prof_mac/prof.ckpt, level 3 about 17,000 cells) to 0.17 ms at 1, 2, 4, 6, 8 threads. Waits for
# the rl10dev diagnostics run to finish; pauses my other profile runs, hot and 1-level (SIGSTOP) during the timed runs, resumes it after.
cd /tmp/prof_scal
while [ ! -e /tmp/prof_rl10dev/done ]; do sleep 20; done   # after the rl10dev step diagnostics, so they do not share cores
HOT=$(pgrep -f '^/tmp/prof_(hot|mac)/crucible_rl10_tier1' | tr '\n' ' ')   # the hot run and the 1-level run (prof_mac's binary)
export CRUCIBLE_MANIFOLD_TABLE=/tmp/lead_tables/rl10_s.bin CRUCIBLE_T1_LEVELS=3 CRUCIBLE_T1_LEAD=20e-3 CRUCIBLE_T1_OUT=prof CRUCIBLE_STEPDIAG=1
for n in 1 2 4 6 8; do
  mkdir -p t$n && cd t$n && cp /tmp/prof_mac/prof_DEVELOPMENT_series.csv .
  [ -n "$HOT" ] && kill -STOP $HOT
  sleep 10
  echo "threads $n launch $(date '+%F %T %Z'), load $(sysctl -n vm.loadavg), paused hot: $HOT" >> ../scaling.log
  t0=$(python3 -c 'import time;print(time.time())')
  CRUCIBLE_T1_THREADS=$n CRUCIBLE_T1_END=0.17e-3 CRUCIBLE_CHECKPOINT=/tmp/prof_scal/t$n/scal CRUCIBLE_RESTART=/tmp/prof_mac/prof.ckpt \
    /usr/bin/time -l /tmp/prof_mac/crucible_rl10_tier1 run development > scal.log 2>&1
  rc=$?; t1=$(python3 -c 'import time;print(time.time())')
  [ -n "$HOT" ] && kill -CONT $HOT
  echo "threads $n exit $rc wall $(python3 -c "print(round($t1-$t0,1))") s, load after $(sysctl -n vm.loadavg)" >> ../scaling.log
  cd ..
done
echo "exit 0" > /tmp/prof_scal/done
