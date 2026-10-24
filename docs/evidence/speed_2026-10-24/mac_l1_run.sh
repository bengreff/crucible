#!/bin/bash
# Hot-regime step and cost profile at 1 level (base 149 x 16 + level 1): tier 1 development, the H2 lead shortened to
# 1 ms (LOX at 1 ms, igniter at 1.2 ms), to 12 ms, step diagnostics on, 2 threads. The 3-level run cannot reach a hot
# state tonight; this one gives the hot regime's binding constraint and wall-time shares on a coarser finest level.
cd /tmp/prof_l1
export CRUCIBLE_MANIFOLD_TABLE=/tmp/lead_tables/rl10_s.bin CRUCIBLE_T1_LEVELS=1 CRUCIBLE_T1_LEAD=1e-3 CRUCIBLE_T1_END=${END:-12e-3} \
  CRUCIBLE_T1_OUT=l1 CRUCIBLE_T1_THREADS=2 CRUCIBLE_CHECKPOINT=/tmp/prof_l1/l1 CRUCIBLE_RESTART=auto CRUCIBLE_STEPDIAG=1
echo "launch $(date)" >> l1.log
/tmp/prof_mac/crucible_rl10_tier1 run development >> l1.log 2>&1
echo "exit $?" > /tmp/prof_l1/done
