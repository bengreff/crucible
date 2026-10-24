#!/bin/bash
# Builds the relaunch tree on backhouse (docs/evidence/RESTART.md, The relaunch): this commit's tree in
# /home/greff/crucible_rst, the two drivers only. The inputs are hard links, so the runs keep them if the old
# directories go: table_mascotte.bin (md5 62e730fb..., S7's) and rl10_s.bin (md5 1be4531e..., B4's).
set -e
R=/home/greff/crucible_rst
mkdir -p $R/build $R/runs/s7low $R/runs/s7high $R/runs/b4
ln -f /home/greff/crucible_s7/build/table_mascotte.bin $R/runs/table_mascotte.bin
ln -f /home/greff/crucible_t1b4/run/rl10_s.bin $R/runs/rl10_s.bin
md5sum $R/runs/table_mascotte.bin $R/runs/rl10_s.bin > $R/runs/inputs.md5
cd $R/build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCANTERA_ROOT=/home/greff/src/cantera -DBOOST_INCLUDE=/usr/include \
  -DCRUCIBLE_BUILD_APP=OFF -DCRUCIBLE_WITH_CANTERA=ON > $R/runs/configure.log 2>&1
nice -n 10 cmake --build . -j 8 --target crucible_mascotte_a10 crucible_rl10_tier1 > $R/runs/build.log 2>&1
echo built > $R/runs/build_done
