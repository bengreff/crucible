#!/bin/bash
# AMR's cost on the RL10 tier-1 mesh at the frozen tags (AMR.md, *AMR's cost on the RL10 tier-1 mesh*): B4's
# development case to 0.5 ms of the H2 lead, on inquiry2's HEAD (REVISION), whose Tags defaults are s* = 1/4 of
# step C's. B4's 0.5 ms line (step C's thresholds) is the comparison. 4 threads, nice 10, queued until the
# 1-minute load is under 24. Table B is B4's file, copied, md5 checked. Writes out/all_done last.
N=/home/greff/crucible_inq2t1; D=$N/out
mkdir -p $D && cd $D || exit 1
date > queued
while [ "$(cut -d. -f1 /proc/loadavg)" -ge 24 ]; do sleep 300; done
date > started; cat /proc/loadavg > load_at_start
fail() { echo "$1" > $D/failed; echo done > $D/all_done; exit 1; }
cp /home/greff/crucible_t1b4/run/rl10_s.bin $D/ || fail table
[ "$(md5sum < $D/rl10_s.bin | cut -d' ' -f1)" = 1be4531e5c446ad70e3640b417fed10e ] || fail "table md5"
mkdir -p $N/build && cd $N/build || fail build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCANTERA_ROOT=/home/greff/src/cantera -DBOOST_INCLUDE=/usr/include \
  -DCRUCIBLE_BUILD_APP=OFF -DCRUCIBLE_WITH_CANTERA=ON > $D/configure.log 2>&1 || fail configure
nice -n 10 cmake --build . -j 4 --target crucible_rl10_tier1 > $D/build.log 2>&1 || fail build
cd $D && cp $N/REVISION revision.txt
export CRUCIBLE_MANIFOLD_TABLE=$D/rl10_s.bin CRUCIBLE_T1_LEAD=20e-3 CRUCIBLE_T1_END=0.5e-3 CRUCIBLE_T1_OUT=fs CRUCIBLE_T1_THREADS=4
nice -n 10 /usr/bin/time -v $N/build/crucible_rl10_tier1 run development > fs.log 2>&1
echo "fs $?" > exits
date > finished
echo done > all_done
