#!/bin/bash
# B4 and the 2-level lead, relaunched fresh on main with the spray NaN fix and the chemistry skip (the old runs'
# binaries had neither; B4r would have stopped at t_L + 0.35 ms, and the skip halves the lead's CPU time). The same
# bodies as restart_2026-10-23/bh_run.sh's b4 and lead2 cases. Started and resumed by /home/greff/resume's runner
# (b4s.job, lead2s.job; nice 10). Waits while the 1-minute load is above 24.
# usage: bh_run.sh b4 | lead2
R=/home/greff/crucible_b4s
B=$R/build
cd $R/runs/$1 || exit 1
while awk '{exit !($1 > 24)}' /proc/loadavg; do sleep 60; done
echo "launch $(date '+%F %T %Z'), load $(cut -d' ' -f1 /proc/loadavg)" >> $1.log
case "$1" in
  b4)
    # The 20 ms H2 lead, to 21 ms, development, 3 levels, 4 threads.
    export CRUCIBLE_MANIFOLD_TABLE=$R/runs/rl10_s.bin CRUCIBLE_T1_LEAD=20e-3 CRUCIBLE_T1_END=21e-3 CRUCIBLE_T1_OUT=b4 \
      CRUCIBLE_T1_THREADS=4 CRUCIBLE_CHECKPOINT=$PWD/b4 CRUCIBLE_RESTART=auto
    exec $B/crucible_rl10_tier1 run development >> b4.log 2>&1 ;;
  lead2)
    # The H2 lead alone at 2 levels to just past t_L: the 2-level member's lead for the 2-against-3-levels evidence.
    export CRUCIBLE_MANIFOLD_TABLE=$R/runs/rl10_s.bin CRUCIBLE_T1_LEVELS=2 CRUCIBLE_T1_LEAD=20e-3 CRUCIBLE_T1_END=20.01e-3 \
      CRUCIBLE_T1_OUT=lead2 CRUCIBLE_T1_THREADS=4 CRUCIBLE_CHECKPOINT=$PWD/lead2 CRUCIBLE_RESTART=auto
    exec $B/crucible_rl10_tier1 run development >> lead2.log 2>&1 ;;
  *) echo "bh_run.sh: unknown run $1" >&2; exit 2 ;;
esac
