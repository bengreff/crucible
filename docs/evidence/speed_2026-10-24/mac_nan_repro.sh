#!/bin/bash
# The spray NaN at about 1.35 ms (l1's error) again, with print-only NANDIAG lines at the first non-finite film
# temperature and B_T: 1 level, the H2 lead shortened to 1 ms, 2 threads, to 1.4 ms. The binary sits in /tmp/prof_hot
# so the thread-scaling script pauses it during its timed runs.
cd /tmp/prof_nan2
export CRUCIBLE_MANIFOLD_TABLE=/tmp/lead_tables/rl10_s.bin CRUCIBLE_T1_LEVELS=1 CRUCIBLE_T1_LEAD=1e-3 CRUCIBLE_T1_END=1.4e-3 \
  CRUCIBLE_T1_OUT=nan CRUCIBLE_T1_THREADS=2
echo "launch $(date)" >> nan.log
/tmp/prof_hot/crucible_rl10_tier1 run development >> nan.log 2>&1
echo "exit $?" > /tmp/prof_nan2/done
