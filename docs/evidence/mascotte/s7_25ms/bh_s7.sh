#!/bin/bash
# S7 production runs (MASCOTTE_A10_PREDICTION.md section 4), the Director's decision of 21 October about 02:20 CDT:
# all four from t = 0 to 25 ms, averaging from 15 ms, CFL 0.8, with the liquid default scheme (the consistent turbulent
# fluxes and extremumBoundedEnthalpy, on by the driver's default; printed at each run's start). 24 threads in total.
# Tree /home/greff/crucible_s7x (inquiry2 at the commit in REVISION); outputs and done-files in runs/.
cd /home/greff/crucible_s7x || exit 1
O=/home/greff/crucible_s7x/runs
mkdir -p $O build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCANTERA_ROOT=/home/greff/src/cantera -DBOOST_INCLUDE=/usr/include \
  -DCRUCIBLE_BUILD_APP=OFF -DCRUCIBLE_WITH_CANTERA=ON > $O/configure.log 2>&1 || { echo configure > $O/failed; exit 1; }
cmake --build . -j 16 --target crucible_mascotte_a10 > $O/build.log 2>&1 || { echo build > $O/failed; exit 1; }
echo built > $O/build_done
B=/home/greff/crucible_s7x/build/crucible_mascotte_a10
TM=/home/greff/crucible_s7/build/table_mascotte.bin   # md5 62e730fb..., the Mac's table
cd $O
run() {  # name size nz threads
  ( CRUCIBLE_CFL=0.8 CRUCIBLE_HOT_TRACE=3450 timeout 172800 $B $2 $3 70 25e-3 15e-3 $4 $O/$1 $TM > $1.log 2>&1
    echo $? > $1_done ) &
}
run low rr:223:3.207:13.7 239 5
run high rr:574:3.207:13.7 239 5
run blob blob 239 5
run blob_fine blob 478 9
wait
echo done > all_done
