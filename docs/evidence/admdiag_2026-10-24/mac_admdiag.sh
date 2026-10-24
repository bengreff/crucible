#!/bin/bash
# D1 and D2 (README.md): tier 1 development at 1 level, the 20 ms H2 lead, 4 threads, with the print-only admdiag
# patch; D1 on 27b135f to 1.79 ms, D2 on main 4f2d81e to 1.80 ms. Builds, then each run through the Mac's slot harness.
for t in old new; do
  cmake -S /tmp/adm/$t -B /tmp/adm/$t/build -DCMAKE_BUILD_TYPE=Release -DCANTERA_ROOT=/Users/ben/src/cantera \
    -DCRUCIBLE_WITH_CANTERA=ON -DCRUCIBLE_BUILD_APP=OFF > /tmp/adm/$t/cmake.log 2>&1 &&
  make -C /tmp/adm/$t/build crucible_rl10_tier1 -j4 > /tmp/adm/$t/build.log 2>&1 || { echo "build $t failed" > /tmp/adm/done; exit 1; }
done
run() {  # name tree end
  mkdir -p /tmp/adm/$1 && cd /tmp/adm/$1
  echo "queued $(date '+%F %T %Z')" > $1.log
  CRUCIBLE_ADMDIAG=1 CRUCIBLE_MANIFOLD_TABLE=/tmp/lead_tables/rl10_s.bin CRUCIBLE_T1_LEVELS=1 CRUCIBLE_T1_LEAD=20e-3 \
    CRUCIBLE_T1_END=$3 CRUCIBLE_T1_OUT=$1 CRUCIBLE_T1_THREADS=4 \
    python3 ~/director/harness/slot.py run --label "inquiry: 1.78 ms diagnosis $1" \
    /usr/bin/time -p /tmp/adm/$2/build/crucible_rl10_tier1 run development >> $1.log 2>&1
  echo "exit $?" > /tmp/adm/$1/done
}
run d1 old 1.79e-3
run d2 new 1.80e-3
echo "all $(date '+%F %T %Z')" > /tmp/adm/done
