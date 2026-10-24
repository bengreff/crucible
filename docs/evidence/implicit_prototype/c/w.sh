#!/bin/bash
# Criterion 5's diagnosis: name window(from,until); C1 64x12 to 8 ms, implicit tau/63 at tolerance 8.975e-5 inside the window.
C=/tmp/implicit/c
mkdir -p $C/$1 && cd $C/$1 || exit 1
CRUCIBLE_IMPLICIT_WINDOW=$2 CRUCIBLE_EQ_TABLE=/Users/ben/inquiry-project/build/table_a.bin CRUCIBLE_IMPLICIT="7.003860e-06,8.975e-5,400,3,1,0,1" \
  /usr/bin/time -p $C/bin5 eqt 64 12 8e-3 1 $1 > $1.log 2>&1
echo "exit $?" > done
