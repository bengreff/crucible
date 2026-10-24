#!/bin/bash
# The implicit prototype's explicit baselines (IMPLICIT_PROTOTYPE.md, schedule day 1): the C1 cold start,
# inviscid, Table A, CFL 0.4, one thread, on 64x12 and 128x24 to 8 ms (tau_C1, the history, the wall per step)
# and 256x48 to 1 ms (the explicit step at the highest stiffness). Queued after the 1.78 ms diagnosis
# (/tmp/adm/done) so that at most two of these chains hold Mac slots; each run goes through slot.py.
B=/tmp/c1base
while [ ! -e /tmp/adm/done ]; do sleep 60; done
run() {  # name nz nr end
  mkdir -p $B/$1 && cd $B/$1 || exit 1
  CRUCIBLE_EQ_TABLE=/Users/ben/inquiry-project/build/table_a.bin python3 ~/director/harness/slot.py run --label "inquiry: C1 explicit baseline $1" \
    /usr/bin/time -p $B/crucible_chamber_study eqt $2 $3 $4 1 $1 > $1.log 2>&1
  echo "exit $?" > done
}
run e64 64 12 8e-3
run e128 128 24 8e-3
run e256 256 48 1e-3
touch $B/done
