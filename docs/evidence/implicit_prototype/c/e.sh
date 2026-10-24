#!/bin/bash
# The explicit C1 baseline in the cost batch (M2, M3): name nz nr end; one thread, no nice (as the implicit runs).
C=/tmp/implicit/c
mkdir -p $C/$1 && cd $C/$1 || exit 1
CRUCIBLE_EQ_TABLE=/Users/ben/inquiry-project/build/table_a.bin /usr/bin/time -p $C/bin4 eqt $2 $3 $4 1 $1 > $1.log 2>&1
echo "exit $?" > done
