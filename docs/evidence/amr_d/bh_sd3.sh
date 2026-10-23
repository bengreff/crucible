#!/bin/bash
# AMR step D's checks rerun on the new spread (b7f232b) and the leak fix: the Director's held note of
# 23 October. /home/greff/crucible_sd3 is inquiry2's HEAD (REVISION); /home/greff/crucible_sd3p is the same
# tree with step D's core changes (2ff039f, bca7c94) reversed, D0's parent (its budget driver keeps only the
# lines its core still has). The verdicts must hold; bit-identity with the 22 October runs is not expected.
N=/home/greff/crucible_sd3; P=/home/greff/crucible_sd3p; O=$N/out; TA=/home/greff/crucible_wf/table_a.bin
mkdir -p $O
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
wait_load() { while [ "$(awk '{print int($1)}' /proc/loadavg)" -gt 24 ]; do sleep 60; done; }
build() {  # all targets, keep going past a parent target that needs step D's API
  mkdir -p $1/build && cd $1/build || return 1
  cmake .. -DCMAKE_BUILD_TYPE=Release -DCANTERA_ROOT=/home/greff/src/cantera -DBOOST_INCLUDE=/usr/include \
    -DCRUCIBLE_BUILD_APP=OFF -DCRUCIBLE_WITH_CANTERA=ON > $O/configure_$2.log 2>&1 || return 1
  cmake --build . -j 4 -- -k > $O/build_$2.log 2>&1
  echo "build $2 exit $?" >> $O/builds.txt
}
date > $O/started
wait_load; build $N n || { echo configure_n > $O/failed; echo done > $O/all_done; exit 1; }
wait_load; build $P p || { echo configure_p > $O/failed; echo done > $O/all_done; exit 1; }
for t in crucible_amr_budget crucible_chamber_study crucible_amr_subcycling_tests; do
  [ -x $N/build/$t ] && [ -x $P/build/$t ] || { echo "missing $t" > $O/failed; echo done > $O/all_done; exit 1; }
done
echo built > $O/build_done
# D0: criterion 2's fixed pattern to 1 ms (no wall functions), C1 eqtt 64 x 12 to 1 ms with law and hp on one
# level, and step C's check 3, each new against the parent.
for v in n p; do
  X=$N/build; [ $v = p ] && X=$P/build
  mkdir -p $O/d0$v && cd $O/d0$v
  wait_load
  CRUCIBLE_EQ_TABLE=$TA $X/crucible_amr_budget 1e-3 4 amr dump=$O/d0$v/amr_dump.txt > amr.log 2>&1
  CRUCIBLE_EQ_TABLE=$TA CRUCIBLE_WALL_FUNCTIONS=law CRUCIBLE_SST_CORRECTION=hp CRUCIBLE_STATE_DUMP=$O/d0$v/c1_state.bin \
    $X/crucible_chamber_study eqtt 64 12 1e-3 4 $O/d0$v/c1 > c1.log 2>&1
  $X/crucible_amr_subcycling_tests c3 4 dump=$O/d0$v/c3 > c3.log 2>&1; echo $? > c3_exit
done
cd $O && { cmp d0n/amr_dump.txt d0p/amr_dump.txt && echo "fixed pattern dump identical"
  cmp d0n/c1_state.bin d0p/c1_state.bin && echo "C1 law hp state identical"
  cmp d0n/c1_history.csv d0p/c1_history.csv && echo "C1 law hp history identical"
  cmp d0n/c3_c3.bin d0p/c3_c3.bin && echo "check 3 dump identical"; } > $O/d0.txt 2>&1
echo done > $O/d0_done
# D3's 1-thread run in the background (1 thread), while the 4-thread runs queue.
( OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 CRUCIBLE_EQ_TABLE=$TA \
    $N/build/crucible_amr_budget 8e-3 1 amr tags wall dump=$O/d1_t1.txt > $O/d3.txt 2>&1; echo $? > $O/d3_run_done ) &
D3=$!
# D2: the fixed pattern, inert, with law and hp, to 8 ms.
wait_load; cd $O && CRUCIBLE_EQ_TABLE=$TA $N/build/crucible_amr_budget 8e-3 4 wall inert > $O/d2.txt 2>&1; echo $? > $O/d2_done
# D1: the tagged case with law and hp, 4 threads, with its dump for D3.
wait_load; CRUCIBLE_EQ_TABLE=$TA $N/build/crucible_amr_budget 8e-3 4 tags wall dump=$O/d1_t4.txt > $O/d1.txt 2>&1; echo $? > $O/d1_done
# D0's fast ctest, new and parent (the FAIL lines compared; a parent test whose target needs step D is named).
wait_load; cd $N/build && ctest -LE slow -j 2 > $O/ctest_n.txt 2>&1
wait_load; cd $P/build && ctest -LE slow -j 2 > $O/ctest_p.txt 2>&1
wait $D3
{ cmp $O/d1_t4.txt $O/d1_t1.txt && echo "D3: 1 and 4 threads identical"; } >> $O/d3.txt 2>&1
date > $O/finished
echo done > $O/all_done
