#!/bin/bash
# One relaunched run (docs/evidence/RESTART.md, The relaunch), started and restarted by /home/greff/resume's runner
# from its job file (s7low.job, s7high.job, b4.job beside this script), which applies nice 10. It waits while the
# 1-minute load is above 24, then starts, or resumes from its last checkpoint (CRUCIBLE_RESTART=auto, every 30 min
# of wall time). The driver writes its own done-file on any exit; a killed run writes none, so the runner resumes it.
# usage: bh_run.sh s7low | s7high | b4
R=/home/greff/crucible_rst
B=$R/build
cd $R/runs/$1 || exit 1
while awk '{exit !($1 > 24)}' /proc/loadavg; do sleep 60; done
case "$1" in
  s7low | s7high)
    # S7 as bh_s7z.sh (crucible_s7z, 97f4fcd): CFL 0.8, the hot trace at 3450 K, 10 ms windows from 15 ms, to 55 ms.
    P=${1#s7}
    [ "$P" = low ] && CASE=rr:223:3.207:13.7 || CASE=rr:574:3.207:13.7
    echo "launch $(date '+%F %T %Z'), load $(cut -d' ' -f1 /proc/loadavg)" >> $P.log
    export CRUCIBLE_CFL=0.8 CRUCIBLE_HOT_TRACE=3450 CRUCIBLE_WINDOW=10e-3 CRUCIBLE_CHECKPOINT=$PWD/$P CRUCIBLE_RESTART=auto
    exec $B/crucible_mascotte_a10 $CASE 239 70 55e-3 15e-3 4 $PWD/$P $R/runs/table_mascotte.bin >> $P.log 2>&1 ;;
  b4)
    # B4 as bh_b4.sh (crucible_t1b4, ab3aa03): the 20 ms H2 lead, to 21 ms, development, 4 threads.
    echo "launch $(date '+%F %T %Z'), load $(cut -d' ' -f1 /proc/loadavg)" >> b4.log
    export CRUCIBLE_MANIFOLD_TABLE=$R/runs/rl10_s.bin CRUCIBLE_T1_LEAD=20e-3 CRUCIBLE_T1_END=21e-3 CRUCIBLE_T1_OUT=b4 \
      CRUCIBLE_T1_THREADS=4 CRUCIBLE_CHECKPOINT=$PWD/b4 CRUCIBLE_RESTART=auto
    exec $B/crucible_rl10_tier1 run development >> b4.log 2>&1 ;;
  *) echo "bh_run.sh: unknown run $1" >&2; exit 2 ;;
esac
