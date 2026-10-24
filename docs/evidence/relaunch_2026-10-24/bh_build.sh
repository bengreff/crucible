#!/bin/bash
# Builds the relaunch tree on backhouse: main with the spray NaN fix and the chemistry skip (SPEED.md, sections 5b
# and 6) in /home/greff/crucible_b4s, the tier 1 driver only. S7's running build in /home/greff/crucible_rst and
# B4r's in /home/greff/crucible_b4r are left alone. The table is a hard link to B4's (rl10_s.bin, md5 1be4531e...).
set -e
R=/home/greff/crucible_b4s
mkdir -p $R/build $R/runs/b4 $R/runs/lead2
ln -f /home/greff/crucible_rst/runs/rl10_s.bin $R/runs/rl10_s.bin
md5sum $R/runs/rl10_s.bin > $R/runs/inputs.md5
cd $R/build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCANTERA_ROOT=/home/greff/src/cantera -DBOOST_INCLUDE=/usr/include \
  -DCRUCIBLE_BUILD_APP=OFF -DCRUCIBLE_WITH_CANTERA=ON > $R/runs/configure.log 2>&1
nice -n 10 cmake --build . -j 8 --target crucible_rl10_tier1 > $R/runs/build.log 2>&1
echo built > $R/runs/build_done
