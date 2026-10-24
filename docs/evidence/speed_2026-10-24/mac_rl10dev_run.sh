#!/bin/bash
# rl10dev (RL10_DEV.md's recipe) with CRUCIBLE_STEPDIAG, uniform then stretched, to 10 ms each.
cd /tmp/prof_rl10dev
export CRUCIBLE_CHAMBER=rl10dev CRUCIBLE_EQ_TABLE=/Users/ben/inquiry-project/build/table_a.bin CRUCIBLE_WALL_FUNCTIONS=law \
  CRUCIBLE_SST_CORRECTION=hp CRUCIBLE_STEPDIAG=1
./crucible_chamber_study eqtt 150 16 10e-3 4 /tmp/prof_rl10dev/u > u.log 2>&1; echo "uniform exit $?" >> status
CRUCIBLE_RADIAL_STRETCHING=1.7347 ./crucible_chamber_study eqtt 150 16 10e-3 4 /tmp/prof_rl10dev/s > s.log 2>&1; echo "stretched exit $?" >> status
touch done
