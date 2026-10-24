#!/bin/bash
# Profile segment: tier 1 development (B4's config: 3 levels, the 20 ms H2 lead) from 0 to 0.15 ms, step diagnostics on.
cd /tmp/prof_mac
export CRUCIBLE_MANIFOLD_TABLE=/tmp/lead_tables/rl10_s.bin CRUCIBLE_T1_LEVELS=3 CRUCIBLE_T1_LEAD=20e-3 CRUCIBLE_T1_END=${END:-0.15e-3} \
  CRUCIBLE_T1_OUT=prof CRUCIBLE_T1_THREADS=4 CRUCIBLE_CHECKPOINT=/tmp/prof_mac/prof CRUCIBLE_RESTART=auto CRUCIBLE_STEPDIAG=1
echo "launch $(date)" >> prof.log
/tmp/prof_mac/crucible_rl10_tier1 run development >> prof.log 2>&1
echo "exit $?" > /tmp/prof_mac/done
