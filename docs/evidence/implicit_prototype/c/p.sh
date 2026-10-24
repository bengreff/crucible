#!/bin/bash
# The explicit C1 references with period-mean histories (bin7: the periodMean fix on 063294f): name nz nr end; one thread.
C=/tmp/implicit/c
mkdir -p $C/$1 && cd $C/$1 || exit 1
CRUCIBLE_EQ_TABLE=/Users/ben/inquiry-project/build/table_a.bin CRUCIBLE_CFL=0.4 /usr/bin/time -p $C/bin7 eqt $2 $3 $4 1 $1 > $1.log 2>&1
echo "exit $?" > done
