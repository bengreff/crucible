#!/bin/sh
# Smoke run of the driver with t90 (not evidence): the five short runs, the judge, and the c4 lines' first 18 fields
# against the first smoke run (the column's refactor must leave them unchanged).
D=/tmp/inq2c/c4/smoke2
O=/tmp/inq2c/c4/smoke
W=/Users/ben/inquiry-project/.claude/worktrees/inquiry2
mkdir -p $D && cd $D || exit 1
cp $W/build/crucible_amr_c4 .
export CRUCIBLE_EQ_TABLE=/Users/ben/inquiry-project/build/table_a.bin OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
( ./crucible_amr_c4 amr scale=1 1.2e-3 2 > s1.log 2>&1; echo "s1 $?" >> exits ) &
( ./crucible_amr_c4 amr scale=0.5 1.2e-3 2 > s2.log 2>&1; echo "s2 $?" >> exits ) &
( ./crucible_amr_c4 uniform 1.2e-3 2 > u.log 2>&1; echo "u $?" >> exits ) &
wait
( ./crucible_amr_c4 amr scale=0.25 1.2e-3 2 > s4.log 2>&1; echo "s4 $?" >> exits ) &
( ./crucible_amr_c4 amr scale=0.125 1.2e-3 2 > s8.log 2>&1; echo "s8 $?" >> exits ) &
wait
for r in s1 s2 s4 s8 u; do
  grep -h '^c4 ' $r.log | cut -d' ' -f1-18 > $r.new
  grep -h '^c4 ' $O/$r.log | cut -d' ' -f1-18 > $r.old
  cmp $r.old $r.new >> same.txt 2>&1 && echo "$r same" >> same.txt
done
sh $W/tests/amr_c4_judge.sh s1.log s2.log s4.log s8.log u.log > judge.txt 2>&1; echo "judge $?" >> exits
echo done > smoke_done
