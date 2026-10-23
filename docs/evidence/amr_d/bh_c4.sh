#!/bin/bash
# AMR criterion 4 on backhouse (AMR.md, *Amendment: the stagnation loop exits at the equilibrium solver's
# accuracy*): c4_runs.sh's runs as registered, in its two waves, 2 threads a job, at nice 10 beside S7 and B4.
# /home/greff/crucible_c4 is inquiry2's HEAD (REVISION); Table A is the Mac's file (same md5). Writes real/all_done
# last (with real/failed if the build fails).
N=/home/greff/crucible_c4; D=$N/real; TA=/home/greff/crucible_wf/table_a.bin
mkdir -p $D && cd $D || exit 1
date > started; cat /proc/loadavg > load_at_start
fail() { echo "$1" > $D/failed; echo done > $D/all_done; exit 1; }
mkdir -p $N/build && cd $N/build || fail build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCANTERA_ROOT=/home/greff/src/cantera -DBOOST_INCLUDE=/usr/include \
  -DCRUCIBLE_BUILD_APP=OFF -DCRUCIBLE_WITH_CANTERA=ON > $D/configure.log 2>&1 || fail configure
cmake --build . -j 4 --target crucible_amr_c4 crucible_amr_budget > $D/build.log 2>&1 || fail build
cd $D && cp $N/build/crucible_amr_c4 $N/build/crucible_amr_budget . && cp $N/REVISION revision.txt
export CRUCIBLE_EQ_TABLE=$TA OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
R="nice -n 10"
( $R ./crucible_amr_c4 amr scale=1 8e-3 2 > s1.log 2>&1; echo "s1 $?" >> exits ) &
( $R ./crucible_amr_c4 amr scale=0.5 8e-3 2 > s2.log 2>&1; echo "s2 $?" >> exits ) &
( $R ./crucible_amr_c4 uniform 8e-3 2 > u.log 2>&1; echo "u $?" >> exits ) &
wait
( $R ./crucible_amr_c4 amr scale=0.25 8e-3 2 > s4.log 2>&1; echo "s4 $?" >> exits ) &
( $R ./crucible_amr_c4 amr scale=0.125 8e-3 2 > s8.log 2>&1; echo "s8 $?" >> exits ) &
( $R ./crucible_amr_budget 8e-3 2 wall tags amr > budget.log 2>&1; echo "budget $?" >> exits ) &
wait
sh $N/tests/amr_c4_judge.sh s1.log s2.log s4.log s8.log u.log > judge.txt 2>&1; echo "judge $?" >> exits
date > finished
echo done > all_done
