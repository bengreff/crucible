#!/bin/bash
# The extended skip (trace atoms; the mask on the partial densities) in the lead: 3 levels, the 20 ms lead, 4 threads,
# 0 to 0.1 ms; skip on and off alternated twice, so the Mac's load from other projects falls on both.
cd /tmp/skip_time3
export CRUCIBLE_MANIFOLD_TABLE=/tmp/lead_tables/rl10_s.bin CRUCIBLE_T1_LEVELS=3 CRUCIBLE_T1_LEAD=20e-3 CRUCIBLE_T1_END=0.1e-3 CRUCIBLE_T1_THREADS=4
for r in 1 2; do for s in 1 0; do
  echo "pass $r skip $s launch $(date '+%T'), load $(sysctl -n vm.loadavg)" >> times.log
  CRUCIBLE_CHEM_SKIP=$s CRUCIBLE_T1_OUT=r${r}skip$s /usr/bin/time -p ./crucible_rl10_tier1 run development > r${r}skip$s.log 2>&1
  echo "pass $r skip $s exit $? $(grep -E '^(real|user)' r${r}skip$s.log | tr '\n' ' ') end $(date '+%T'), load $(sysctl -n vm.loadavg)" >> times.log
done; done
echo done > done
