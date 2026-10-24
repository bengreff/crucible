#!/bin/bash
# The Director's 00:50 item 2: backhouse contention. (1) A per-virtual-CPU probe: one thread of fixed work pinned to
# each of the 28 virtual CPUs in turn, twice, CPU time and wall recorded (an E-core-backed CPU would show as slower
# in CPU time, if WSL maps virtual CPUs to fixed cores). (2) The thread-scaling segment (0.15 to 0.17 ms from
# seed15.ckpt, 4 threads) with every job of ours paused (SIGSTOP; resumed after): unpinned at nice 10, pinned to
# virtual CPUs 0,2,4,6 at nice 10, and unpinned at nice 0. Other projects' jobs are not touched.
R=/home/greff/crucible_prof
cd $R/run || exit 1
L=$R/run/contention.log
gcc -O2 -o probe /home/greff/crucible_prof/run/probe.c || exit 1
cpu() { awk '/^cpu /{print $2+$3+$4+$7+$8+$9, $5+$6}' /proc/stat; }
if [ "$1" != "segment" ]; then
echo "probe launch $(date '+%F %T %Z'), load $(cat /proc/loadavg)" >> $L
for pass in 1 2; do
  for c in $(seq 0 27); do
    /usr/bin/time -f "%U %e" -o probe.t taskset -c $c ./probe
    echo "probe pass $pass cpu $c user,wall $(cat probe.t)" >> $L
  done
done
echo "probe done $(date '+%F %T %Z'), load $(cat /proc/loadavg)" >> $L
fi
[ "$1" = "probe" ] && exit 0
OURS=$(pgrep -f 'crucible_(rl10_tier1|mascotte_a10)' | grep -v "^$$\$" | tr '\n' ' ')
trap 'kill -CONT $OURS 2>/dev/null' EXIT
kill -STOP $OURS
echo "paused ours: $OURS ($(date '+%F %T %Z'))" >> $L
sleep 20
export CRUCIBLE_MANIFOLD_TABLE=$R/run/rl10_s.bin CRUCIBLE_T1_LEVELS=3 CRUCIBLE_T1_LEAD=20e-3 CRUCIBLE_T1_OUT=prof CRUCIBLE_STEPDIAG=1
run() {  # tag, nice, taskset prefix
  local tag=$1 ni=$2; shift 2
  echo "$tag launch $(date '+%F %T %Z'), load $(cat /proc/loadavg)" >> $L
  read b0 i0 < <(cpu); t0=$(date +%s.%N)
  CRUCIBLE_T1_THREADS=4 CRUCIBLE_T1_END=0.17e-3 CRUCIBLE_CHECKPOINT=$R/run/cont_$tag CRUCIBLE_RESTART=$R/run/seed15.ckpt \
    nice -n $ni /usr/bin/time -v "$@" $R/build/crucible_rl10_tier1 run development > cont_$tag.log 2>&1
  local rc=$?; t1=$(date +%s.%N); read b1 i1 < <(cpu)
  echo "$tag exit $rc wall $(echo "$t1 - $t0" | bc) s, machine busy $(echo "scale=1; 100*($b1-$b0)/($b1-$b0+$i1-$i0)" | bc)%, $(grep -E 'User time|Percent of CPU' cont_$tag.log | tr -s ' \t' ' ' | tr '\n' ';')" >> $L
}
run unpinned_n10 10
run pinned0246_n10 10 taskset -c 0,2,4,6
run unpinned_n0 0
kill -CONT $OURS
echo "resumed ours ($(date '+%F %T %Z')), load $(cat /proc/loadavg)" >> $L
echo "exit 0" > cont_done
