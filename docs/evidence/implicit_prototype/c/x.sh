#!/bin/bash
# An explicit C1 run at a given CFL (criterion 5's diagnosis): name nz nr end cfl; one thread, nice 10.
C=/tmp/implicit/c
mkdir -p $C/$1 && cd $C/$1 || exit 1
CRUCIBLE_EQ_TABLE=/Users/ben/inquiry-project/build/table_a.bin CRUCIBLE_CFL=$5 nice -n 10 /usr/bin/time -p $C/bin4 eqt $2 $3 $4 1 $1 > $1.log 2>&1
echo "exit $?" > done
