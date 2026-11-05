#!/bin/bash
# Check F3 (SPEED_PLAN.md 7.6): the new test on the unfixed spray.cpp (A), then on the fix (B), then the spray
# ctests and tier 1's selftest on B. A scratch tree from git archive of 95f0ed3; nothing else is written.
D=/home/greff/crucible_f3
{
  date; uptime
  cd $D || exit 1
  cp fix/spray_tests.cpp tests/spray_tests.cpp
  nice -n 10 cmake -S $D -B $D/build -DCMAKE_BUILD_TYPE=Release -DCANTERA_ROOT=/home/greff/src/cantera -DCRUCIBLE_WITH_CANTERA=ON -DCRUCIBLE_BUILD_APP=OFF > configure.log 2>&1 || { echo configure failed; tail -20 configure.log; exit 1; }
  nice -n 10 cmake --build build --target crucible_spray_tests -j4 > buildA.log 2>&1 || { echo build A failed; grep error buildA.log | head; exit 1; }
  md5sum core/spray.cpp
  nice -n 10 build/crucible_spray_tests > A.log 2>&1; echo "A exit $?"
  grep -v CanteraWarning A.log | grep -A6 "^F3"
  grep -E "FAIL|failure" A.log | head -8
  cp fix/spray.cpp core/spray.cpp
  md5sum core/spray.cpp
  nice -n 10 cmake --build build -j4 > buildB.log 2>&1 || { echo build B failed; grep error buildB.log | head; exit 1; }
  nice -n 10 build/crucible_spray_tests > B.log 2>&1; echo "B exit $?"
  grep -v CanteraWarning B.log | grep -A6 "^F3"
  grep -E "FAIL|failure" B.log | head -8
  diff <(grep -E "PASS|FAIL" A.log) <(grep -E "PASS|FAIL" B.log) > AB.diff; echo "A against B, check lines differing: $(grep -c '^<' AB.diff)"
  (cd build && nice -n 10 ctest -R "spray|face_inlet|chemistry_skip" --timeout 2400 --output-on-failure > ../ctest.log 2>&1; echo "ctest exit $?"; tail -12 ../ctest.log)
  CRUCIBLE_MANIFOLD_TABLE=/home/greff/crucible_p0/rl10_s.bin nice -n 10 build/crucible_rl10_tier1 selftest > selftest.log 2>&1; echo "selftest exit $?"
  grep -E "check\(s\) failed|FAIL" selftest.log | head
  md5sum build/crucible_rl10_tier1
  date
} > $D.out 2>&1
echo done > $D.done
