#!/bin/bash
# Profile segment, hot: tier 1 development, 3 levels, a 1 ms H2 lead (LOX at 1 ms, igniter at 1.2 ms), step diagnostics on.
cd /tmp/prof_hot
export CRUCIBLE_MANIFOLD_TABLE=/tmp/lead_tables/rl10_s.bin CRUCIBLE_T1_LEVELS=3 CRUCIBLE_T1_LEAD=1e-3 CRUCIBLE_T1_END=${END:-2.5e-3} \
  CRUCIBLE_T1_OUT=hot CRUCIBLE_T1_THREADS=4 CRUCIBLE_CHECKPOINT=/tmp/prof_hot/hot CRUCIBLE_RESTART=auto CRUCIBLE_STEPDIAG=1
echo "launch $(date)" >> hot.log
/tmp/prof_hot/crucible_rl10_tier1 run development >> hot.log 2>&1
echo "exit $?" > /tmp/prof_hot/done
