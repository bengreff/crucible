#!/bin/bash
# S7 extension (MASCOTTE_A10_PREDICTION.md section 2: a run not settled is extended by 10 ms and the average window
# moves with it). None of the four 25 ms runs in /home/greff/crucible_s7x/runs was settled (half-to-half drift at
# r <= 20 mm: low 93.4 K, high 117.5 K, blob 620.7 K, blob-fine 380.6 K; limit 50 K). The driver has no restart, so
# each run restarts from t = 0 on the same code (crucible_s7x, inquiry2 at 97f4fcd) with only the driver's output
# option CRUCIBLE_WINDOW added (inquiry main, tests/mascotte_a10.cpp): consecutive 10 ms windows from 15 ms, each
# written as <name>_to<ms>_* when it closes. The first window that meets the 50 K test is the run's result.
# Inputs, CFL, hot trace and table as bh_s7.sh; threads 4/4/4/6 (18 of 28, room for inquiry2).
cd /home/greff/crucible_s7z || exit 1
O=/home/greff/crucible_s7z/runs
mkdir -p $O build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCANTERA_ROOT=/home/greff/src/cantera -DBOOST_INCLUDE=/usr/include \
  -DCRUCIBLE_BUILD_APP=OFF -DCRUCIBLE_WITH_CANTERA=ON > $O/configure.log 2>&1 || { echo configure > $O/failed; exit 1; }
cmake --build . -j 12 --target crucible_mascotte_a10 > $O/build.log 2>&1 || { echo build > $O/failed; exit 1; }
echo built > $O/build_done
B=/home/greff/crucible_s7z/build/crucible_mascotte_a10
TM=/home/greff/crucible_s7/build/table_mascotte.bin   # md5 62e730fb..., the Mac's table
cd $O
run() {  # name size nz threads end
  ( CRUCIBLE_CFL=0.8 CRUCIBLE_HOT_TRACE=3450 CRUCIBLE_WINDOW=10e-3 timeout 345600 $B $2 $3 70 $5 15e-3 $4 $O/$1 $TM > $1.log 2>&1
    echo $? > $1_done ) &
}
run low rr:223:3.207:13.7 239 4 55e-3
run high rr:574:3.207:13.7 239 4 55e-3
run blob blob 239 4 85e-3
run blob_fine blob 478 6 85e-3
wait
echo done > all_done
