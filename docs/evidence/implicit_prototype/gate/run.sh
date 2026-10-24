#!/bin/bash
export CRUCIBLE_EQ_TABLE=/Users/ben/inquiry-project/build/table_a.bin
cp /Users/ben/inquiry-project/build/crucible_chamber_study /tmp/implicit/gate/new/bin
cd /tmp/implicit/gate/old && CRUCIBLE_STATE_DUMP=/tmp/implicit/gate/old/state.bin nice -n 10 /usr/bin/time -p /tmp/c1base/crucible_chamber_study eqt 64 12 1e-3 1 g > g.log 2>&1 &
cd /tmp/implicit/gate/new && CRUCIBLE_STATE_DUMP=/tmp/implicit/gate/new/state.bin nice -n 10 /usr/bin/time -p /tmp/implicit/gate/new/bin eqt 64 12 1e-3 1 g > g.log 2>&1 &
wait
cmp /tmp/implicit/gate/old/state.bin /tmp/implicit/gate/new/state.bin && echo "state bit-identical" > /tmp/implicit/gate/result || echo "state DIFFERS" > /tmp/implicit/gate/result
cmp /tmp/implicit/gate/old/g_history.csv /tmp/implicit/gate/new/g_history.csv && echo "history identical" >> /tmp/implicit/gate/result || echo "history differs" >> /tmp/implicit/gate/result
echo done >> /tmp/implicit/gate/result
