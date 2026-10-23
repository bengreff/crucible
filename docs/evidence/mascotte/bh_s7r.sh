#!/bin/bash
# S7 relaunch (22 October 2026, started 22:26 CDT). bh_s7z.sh's four runs were stopped at about 17:00 CDT on 22 October
# (backhouse's default tmux server exited; low at 51.5 ms, high 51.8, blob 27.1, blob_fine 25.7), before their registered
# ends (MASCOTTE_A10_PREDICTION.md section 4: 55 ms for low and high, 85 ms for the blob legs). No window of any run met
# section 2's 50 K test. The driver has no restart, so each run restarts from t = 0 with bh_s7z.sh's binary (built
# 21 October from crucible_s7z: 97f4fcd plus CRUCIBLE_WINDOW), table, inputs, CFL, hot trace and threads, unchanged.
# Output goes to runs2 so the stopped runs' files stay in runs. Started in tmux on its own socket (-L crucible_s7r),
# so a kill-server on the default socket does not reach it.
O=/home/greff/crucible_s7z/runs2
mkdir -p $O
B=/home/greff/crucible_s7z/build/crucible_mascotte_a10
TM=/home/greff/crucible_s7/build/table_mascotte.bin   # md5 62e730fb..., the Mac's table
md5sum $B $TM > $O/inputs_md5.txt
cd $O
run() {  # name size nz threads end
  ( CRUCIBLE_CFL=0.8 CRUCIBLE_HOT_TRACE=3450 CRUCIBLE_WINDOW=10e-3 timeout 345600 $B $2 $3 70 $5 15e-3 $4 $O/$1 $TM > $1.log 2>&1
    echo $? > $1_done ) &
}
run low rr:223:3.207:13.7 239 4 55e-3
run high rr:574:3.207:13.7 239 4 55e-3
run blob blob 239 4 85e-3
run blob_fine blob 478 6 85e-3
wait
echo done > all_done
