#!/bin/bash
# One implicit case C run: name nz nr end dt tolerance max_iterations pseudo_cfl extrapolate; one thread.
C=/tmp/implicit/c
mkdir -p $C/$1 && cd $C/$1 || exit 1
CRUCIBLE_EQ_TABLE=/Users/ben/inquiry-project/build/table_a.bin CRUCIBLE_IMPLICIT="$5,$6,$7,$8,$9,0,${10:-0}" \
  /usr/bin/time -p ${BIN:-$C/bin} eqt $2 $3 $4 1 $1 > $1.log 2>&1
echo "exit $?" > done
