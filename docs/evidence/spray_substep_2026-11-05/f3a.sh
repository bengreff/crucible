#!/bin/bash
# spray_amr_verification on the unfixed spray.cpp (A), to diff its check lines against the fix's (ctest.log).
D=/home/greff/crucible_f3; A=/home/greff/crucible_f3a
{
  date; uptime
  rsync -a --exclude build $D/ $A/ || exit 1
  git -C /dev/null --version > /dev/null 2>&1
  cp $D/fix/spray.cpp /dev/null
  cd $A && tar xf /home/greff/crucible_f3a_head_spray.tar && md5sum core/spray.cpp
  nice -n 10 cmake -S $A -B $A/build -DCMAKE_BUILD_TYPE=Release -DCANTERA_ROOT=/home/greff/src/cantera -DCRUCIBLE_WITH_CANTERA=ON -DCRUCIBLE_BUILD_APP=OFF > configure.log 2>&1 || { echo configure failed; exit 1; }
  nice -n 10 cmake --build build --target crucible_spray_amr_tests -j4 > build.log 2>&1 || { echo build failed; exit 1; }
  nice -n 10 build/crucible_spray_amr_tests > amrA.log 2>&1; echo "A exit $?"
  awk '/Start 17/,0' $D/ctest.log | grep -E "PASS|FAIL" > amrB.lines
  grep -E "PASS|FAIL" amrA.log > amrA.lines
  echo "lines A $(wc -l < amrA.lines), B $(wc -l < amrB.lines)"
  diff <(sed 's/^[0-9]*: //' amrA.lines) <(sed 's/^[0-9]*: //' amrB.lines) > amrAB.diff; echo "differing: $(grep -c '^<' amrAB.diff)"; head -10 amrAB.diff
  date
} > $A.out 2>&1
echo done > $A.done
