#!/bin/bash
# Lever B's runs (SPEED.md section 5c) on main 4f2d81e (the binary copied to /tmp/cfl), one after another, each
# through the Mac's slot harness, 4 threads.
BIN=/tmp/cfl/crucible_rl10_tier1
run() {  # name levels lead end cfl seed
  mkdir -p /tmp/cfl/$1 && cd /tmp/cfl/$1
  echo "queued $(date '+%F %T %Z')" > $1.log
  CRUCIBLE_MANIFOLD_TABLE=/tmp/lead_tables/rl10_s.bin CRUCIBLE_T1_LEVELS=$2 CRUCIBLE_T1_LEAD=$3 CRUCIBLE_T1_END=$4 \
    CRUCIBLE_T1_CFL=$5 CRUCIBLE_T1_SEED=$6 CRUCIBLE_T1_OUT=$1 CRUCIBLE_T1_THREADS=4 \
    python3 ~/director/harness/slot.py run --label "inquiry: lever B $1" /usr/bin/time -p $BIN run development >> $1.log 2>&1
  echo "exit $?" > /tmp/cfl/$1/done
}
run L4 3 20e-3 1e-3 0.4 1
run L8 3 20e-3 1e-3 0.8 1
run H4 1 1e-3 2e-3 0.4 1
run H8 1 1e-3 2e-3 0.8 1
run H4s 1 1e-3 2e-3 0.4 2
echo "all $(date '+%F %T %Z')" > /tmp/cfl/done
