#!/bin/bash
# momdiag runs: name nz nr end cfl implicit-spec (empty for explicit); one thread, nice 10.
R=/tmp/momdiag/runs
mkdir -p $R/$1 && cd $R/$1 || exit 1
if [ -n "$6" ]; then export CRUCIBLE_IMPLICIT="$6"; fi
CRUCIBLE_MOMDIAG=7.98e-3 CRUCIBLE_EQ_TABLE=/Users/ben/inquiry-project/build/table_a.bin CRUCIBLE_CFL=$5 nice -n 10 /usr/bin/time -p ${BIN:-/tmp/momdiag/bin} eqt $2 $3 $4 1 $1 > $1.log 2>&1
echo "exit $?" > done
