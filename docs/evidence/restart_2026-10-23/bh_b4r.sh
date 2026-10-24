#!/bin/bash
# B4 relaunched on main with inquiry2's reflux clip (f8e0f65, on main at cab591b; its D0 identity and trace checks
# passed, inquiry2 24 October about 00:20 CDT). B4 as bh_run.sh's b4 case, in its own tree so S7's running build in
# /home/greff/crucible_rst is not touched: the 20 ms H2 lead, to 21 ms, development, 4 threads. Started and resumed
# by /home/greff/resume's runner from b4r.job (nice 10). Waits while the 1-minute load is above 24.
R=/home/greff/crucible_b4r
cd $R/runs/b4r || exit 1
while awk '{exit !($1 > 24)}' /proc/loadavg; do sleep 60; done
echo "launch $(date '+%F %T %Z'), load $(cut -d' ' -f1 /proc/loadavg)" >> b4.log
export CRUCIBLE_MANIFOLD_TABLE=$R/runs/rl10_s.bin CRUCIBLE_T1_LEAD=20e-3 CRUCIBLE_T1_END=21e-3 CRUCIBLE_T1_OUT=b4 \
  CRUCIBLE_T1_THREADS=4 CRUCIBLE_CHECKPOINT=$PWD/b4 CRUCIBLE_RESTART=auto
exec $R/build/crucible_rl10_tier1 run development >> b4.log 2>&1
