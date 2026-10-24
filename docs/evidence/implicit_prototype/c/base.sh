#!/bin/bash
# The explicit C1 baselines for the implicit prototype: name nz nr end; one thread, through slot.py.
C=/tmp/implicit/c
mkdir -p $C/$1 && cd $C/$1 || exit 1
CRUCIBLE_EQ_TABLE=/Users/ben/inquiry-project/build/table_a.bin python3 ~/director/harness/slot.py run --label "inquiry: C1 explicit baseline $1 (1 thread)" -- \
  /usr/bin/time -p $C/bin eqt $2 $3 $4 1 $1 > $1.log 2>&1
echo "exit $?" > done
