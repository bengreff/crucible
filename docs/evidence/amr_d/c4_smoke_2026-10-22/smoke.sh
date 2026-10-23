#!/bin/sh
# Identity check of the header move, and a short run of the criterion 4 pipeline (not evidence).
D=/tmp/inq2c/c4/smoke
W=/Users/ben/inquiry-project/.claude/worktrees/inquiry2
TA=/Users/ben/inquiry-project/build/table_a.bin
mkdir -p $D && cd $D || exit 1
cp $W/build/crucible_amr_budget new_budget
cp $W/build/crucible_amr_c4 .
export CRUCIBLE_EQ_TABLE=$TA OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
( /tmp/inq2c/c4/ref_budget 4e-4 2 wall tags amr dump=$D/ref.bin > ref.log 2>&1; echo "ref $?" >> exits ) &
( ./new_budget 4e-4 2 wall tags amr dump=$D/new.bin > new.log 2>&1; echo "new $?" >> exits ) &
( ./crucible_amr_c4 amr scale=1 1.2e-3 2 > s1.log 2>&1; echo "s1 $?" >> exits ) &
( ./crucible_amr_c4 amr scale=0.5 1.2e-3 2 > s2.log 2>&1; echo "s2 $?" >> exits ) &
wait
( ./crucible_amr_c4 amr scale=0.25 1.2e-3 2 > s4.log 2>&1; echo "s4 $?" >> exits ) &
( ./crucible_amr_c4 amr scale=0.125 1.2e-3 2 > s8.log 2>&1; echo "s8 $?" >> exits ) &
( ./crucible_amr_c4 uniform 1.2e-3 2 > u.log 2>&1; echo "u $?" >> exits ) &
wait
cmp ref.bin new.bin > cmp.txt 2>&1; echo "cmp $?" >> exits
sh $W/tests/amr_c4_judge.sh s1.log s2.log s4.log s8.log u.log > judge.txt 2>&1; echo "judge $?" >> exits
echo done > smoke_done
