#!/bin/bash
# Thread scaling (the Director's 23:50 item 3) on the measurement tree (main 1166721 + inquiry2's reflux clip + step
# diagnostics, not committed): tier 1 development, B4's config (3 levels, the 20 ms H2 lead). A seed run to 0.15 ms at
# 4 threads ends with a checkpoint (level 3 then about 17,000 cells, as on the Mac); then 0.15 to 0.17 ms restarted from
# it at 1, 2, 4, 8, 16 threads, each timed. No load gate (other projects hold the 1-minute load at 32 to 37 tonight):
# each timed run records the whole machine's busy CPU over its span from /proc/stat, so the numbers carry their load.
# My own S7 runs are paused (SIGSTOP) during the timed runs, resumed after.
R=/home/greff/crucible_prof
cd $R/run || exit 1
ln -f /home/greff/crucible_rst/runs/rl10_s.bin rl10_s.bin
cpu() { awk '/^cpu /{print $2+$3+$4+$7+$8+$9, $5+$6}' /proc/stat; }   # busy, idle jiffies
export CRUCIBLE_MANIFOLD_TABLE=$R/run/rl10_s.bin CRUCIBLE_T1_LEVELS=3 CRUCIBLE_T1_LEAD=20e-3 CRUCIBLE_T1_OUT=prof CRUCIBLE_STEPDIAG=1
if [ ! -f seed15.ckpt ]; then
  echo "seed launch $(date '+%F %T %Z'), load $(cat /proc/loadavg)" >> scaling.log
  CRUCIBLE_T1_THREADS=4 CRUCIBLE_T1_END=0.15e-3 CRUCIBLE_CHECKPOINT=$R/run/seed15 nice -n 10 $R/build/crucible_rl10_tier1 run development > seed.log 2>&1
  echo "seed exit $? $(date '+%F %T %Z')" >> scaling.log
fi
for n in 1 2 4 8 16; do
  [ -f scal_$n.done ] && continue
  S7=$(pgrep -f 'crucible_rst/build/crucible_mascotte_a10' | tr '\n' ' ')
  [ -n "$S7" ] && kill -STOP $S7
  sleep 20
  echo "threads $n launch $(date '+%F %T %Z'), load $(cat /proc/loadavg), paused S7: $S7" >> scaling.log
  read b0 i0 < <(cpu); t0=$(date +%s.%N)
  CRUCIBLE_T1_THREADS=$n CRUCIBLE_T1_END=0.17e-3 CRUCIBLE_CHECKPOINT=$R/run/scal_$n CRUCIBLE_RESTART=$R/run/seed15.ckpt \
    nice -n 10 /usr/bin/time -v $R/build/crucible_rl10_tier1 run development > scal_$n.log 2>&1
  rc=$?; t1=$(date +%s.%N); read b1 i1 < <(cpu)
  [ -n "$S7" ] && kill -CONT $S7
  echo "threads $n exit $rc wall $(echo "$t1 - $t0" | bc) s, machine busy $(echo "scale=1; 100*($b1-$b0)/($b1-$b0+$i1-$i0)" | bc)% of $(nproc) threads over the run, load after $(cat /proc/loadavg)" >> scaling.log
  touch scal_$n.done
done
echo "exit 0" > prof_done
