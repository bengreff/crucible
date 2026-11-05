#!/bin/bash
# Print-only probe 2 (5 November): lb's per-step liquid-budget print plus the NaN residual's terms at the throw in
# Spray::rates. A copy of crucible_lb; nothing in crucible_p0 is written. From i10k's 02:38 checkpoint under t1_r's path.
S=/home/greff/crucible_p0; L=/home/greff/crucible_lb; D=/home/greff/crucible_lb2
{
  date; uptime
  rsync -a --exclude build --exclude i10k $L/ $D/ || { echo rsync failed; exit 1; }
  python3 /home/greff/nanprobe.py $D || { echo patch failed; exit 1; }
  nice -n 10 cmake -S $D -B $D/build -DCMAKE_BUILD_TYPE=Release -DCANTERA_ROOT=/home/greff/src/cantera -DCRUCIBLE_WITH_CANTERA=ON -DCRUCIBLE_BUILD_APP=OFF > $D/configure.log 2>&1 || { echo configure failed; tail -20 $D/configure.log; exit 1; }
  nice -n 10 cmake --build $D/build --target crucible_rl10_tier1 -j4 > $D/build.log 2>&1 || { echo build failed; grep -E "error" $D/build.log | head -20; exit 1; }
  echo built; date
  mkdir -p $D/i10k && cd $D/i10k || exit 1
  cp $S/agree/i10k/ck.ckpt $S/agree/i10k/i10k_DEVELOPMENT_series.csv .
  export CRUCIBLE_MANIFOLD_TABLE=$S/rl10_s.bin CRUCIBLE_T1_LEVELS=0 CRUCIBLE_T1_THREADS=1 CRUCIBLE_T1_END=21.05e-3
  export CRUCIBLE_T1_SEED=1 CRUCIBLE_T1_OUT=i10k CRUCIBLE_CHECKPOINT=$D/i10k/ck CRUCIBLE_CHECKPOINT_EVERY=100000
  export CRUCIBLE_RESTART=auto OMP_NUM_THREADS=1 CRUCIBLE_T1_IMPLICIT="10e-6,8.975e-5,400,3,1,3e-3" LBPROBE_ON=1
  echo "probe $(date '+%F %T %Z'), load $(cut -d' ' -f1 /proc/loadavg)"
  nice -n 10 bash -c "exec -a $S/t1_r $D/build/crucible_rl10_tier1 run development" > lb2.log 2>&1
  echo "probe exit $?"; grep -v CanteraWarning lb2.log | grep -v LBPROBE | tail -4 | cut -c1-300; grep -c "^LBPROBE" lb2.log; grep -c "^NANPROBE" lb2.log; grep -c "^HPROBE" lb2.log
  date
} > $D.out 2>&1
echo done > $D.done
