#!/bin/bash
# Step D's reruns on the new spread, resumed (AMR.md, *Step D's checks rerun on the new spread*): bh_sd3.sh from
# D0's parent runs on, the same commands, at nice 10 (the whole script is launched under nice) and with no load
# gate. bh_sd3.sh's gate (load <= 24) cannot open beside S7 and B4 (load 35 to 37, CPU 95% busy, measured
# 23 October about 01:34 CDT); criterion 4's relaunch that night runs beside them at nice 10 by the Director's
# rule. The builds and D0's new-tree runs (out/d0n) finished under bh_sd3.sh and are kept.
N=/home/greff/crucible_sd3; P=/home/greff/crucible_sd3p; O=$N/out; TA=/home/greff/crucible_wf/table_a.bin
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
date > $O/resumed
X=$P/build
mkdir -p $O/d0p && cd $O/d0p
CRUCIBLE_EQ_TABLE=$TA $X/crucible_amr_budget 1e-3 4 amr dump=$O/d0p/amr_dump.txt > amr.log 2>&1
CRUCIBLE_EQ_TABLE=$TA CRUCIBLE_WALL_FUNCTIONS=law CRUCIBLE_SST_CORRECTION=hp CRUCIBLE_STATE_DUMP=$O/d0p/c1_state.bin \
  $X/crucible_chamber_study eqtt 64 12 1e-3 4 $O/d0p/c1 > c1.log 2>&1
$X/crucible_amr_subcycling_tests c3 4 dump=$O/d0p/c3 > c3.log 2>&1; echo $? > c3_exit
cd $O && { cmp d0n/amr_dump.txt d0p/amr_dump.txt && echo "fixed pattern dump identical"
  cmp d0n/c1_state.bin d0p/c1_state.bin && echo "C1 law hp state identical"
  cmp d0n/c1_history.csv d0p/c1_history.csv && echo "C1 law hp history identical"
  cmp d0n/c3_c3.bin d0p/c3_c3.bin && echo "check 3 dump identical"; } > $O/d0.txt 2>&1
echo done > $O/d0_done
( OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 CRUCIBLE_EQ_TABLE=$TA \
    $N/build/crucible_amr_budget 8e-3 1 amr tags wall dump=$O/d1_t1.txt > $O/d3.txt 2>&1; echo $? > $O/d3_run_done ) &
D3=$!
cd $O && CRUCIBLE_EQ_TABLE=$TA $N/build/crucible_amr_budget 8e-3 4 wall inert > $O/d2.txt 2>&1; echo $? > $O/d2_done
CRUCIBLE_EQ_TABLE=$TA $N/build/crucible_amr_budget 8e-3 4 tags wall dump=$O/d1_t4.txt > $O/d1.txt 2>&1; echo $? > $O/d1_done
cd $N/build && ctest -LE slow -j 2 > $O/ctest_n.txt 2>&1
cd $P/build && ctest -LE slow -j 2 > $O/ctest_p.txt 2>&1
wait $D3
{ cmp $O/d1_t4.txt $O/d1_t1.txt && echo "D3: 1 and 4 threads identical"; } >> $O/d3.txt 2>&1
date > $O/finished
echo done > $O/all_done
