#!/bin/sh
# AMR.md, *Fixing the tag thresholds, then criterion 4*: the five registered runs to 8 ms, then the judge; and
# criterion 2's fill content and D1's budgets on the fixed spread (*Criterion 4's driver*). Mac, 2 threads a job,
# in two waves. Launch from the worktree after the spread fix is on main and the build is rebased:
#   nohup sh docs/evidence/amr_d/c4_runs.sh > /tmp/inq2c/c4/real.out 2>&1 & disown
# Writes /tmp/inq2c/c4/real/all_done last.
D=/tmp/inq2c/c4/real
W=/Users/ben/inquiry-project/.claude/worktrees/inquiry2
mkdir -p $D && cd $D || exit 1
git -C $W log --oneline -1 > revision.txt
cp $W/build/crucible_amr_c4 $W/build/crucible_amr_budget .
export CRUCIBLE_EQ_TABLE=/Users/ben/inquiry-project/build/table_a.bin OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
( ./crucible_amr_c4 amr scale=1 8e-3 2 > s1.log 2>&1; echo "s1 $?" >> exits ) &
( ./crucible_amr_c4 amr scale=0.5 8e-3 2 > s2.log 2>&1; echo "s2 $?" >> exits ) &
( ./crucible_amr_c4 uniform 8e-3 2 > u.log 2>&1; echo "u $?" >> exits ) &
wait
( ./crucible_amr_c4 amr scale=0.25 8e-3 2 > s4.log 2>&1; echo "s4 $?" >> exits ) &
( ./crucible_amr_c4 amr scale=0.125 8e-3 2 > s8.log 2>&1; echo "s8 $?" >> exits ) &
( ./crucible_amr_budget 8e-3 2 wall tags amr > budget.log 2>&1; echo "budget $?" >> exits ) &
wait
sh $W/tests/amr_c4_judge.sh s1.log s2.log s4.log s8.log u.log > judge.txt 2>&1; echo "judge $?" >> exits
echo done > all_done
