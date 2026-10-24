#!/bin/bash
# The reflux's rho k clip (AMR.md, Amendment: the reflux clips rho k at zero): its registered checks on backhouse.
# /home/greff/crucible_rfx is inquiry2's HEAD (REVISION); the parent is fcb20f1, /home/greff/crucible_inq2t1.
# Chain A: check 3's tier 1 at the frozen tags to 0.5 ms (AMR's cost on the tier-1 mesh). Chain B: check 1's
# three short D0 cases, new against the parent; check 3's print-only trace to 0.1 ms; check 2 (D3, D1's case
# to 1 ms, 1 and 4 threads); then, once tier 1 is done, check 1's D1 to 8 ms, new and parent together.
# 4 threads a run, nice 10, each start queued until the 1-minute load is under 24. Writes out/all_done last.
N=/home/greff/crucible_rfx; P=/home/greff/crucible_inq2t1; O=$N/out; TA=/home/greff/crucible_wf/table_a.bin
mkdir -p $O && cd $O || exit 1
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
wait_load() { while [ "$(cut -d. -f1 /proc/loadavg)" -ge 24 ]; do sleep 120; done; }
fail() { echo "$1" > $O/failed; echo done > $O/all_done; exit 1; }
date > started
TARGETS="crucible_rl10_tier1 crucible_amr_budget crucible_chamber_study crucible_amr_subcycling_tests"
build() {  # $1 tree, $2 its label, $3 extra cmake flags
  mkdir -p $1/build && cd $1/build || return 1
  cmake .. -DCMAKE_BUILD_TYPE=Release -DCANTERA_ROOT=/home/greff/src/cantera -DBOOST_INCLUDE=/usr/include \
    -DCRUCIBLE_BUILD_APP=OFF -DCRUCIBLE_WITH_CANTERA=ON > $O/configure_$2.log 2>&1 || return 1
  nice -n 10 cmake --build . -j 4 --target $TARGETS > $O/build_$2.log 2>&1
}
wait_load; build $N n || fail build_n
wait_load; build $P p || fail build_p
cp $N/REVISION $O/revision.txt
cp /home/greff/crucible_t1b4/run/rl10_s.bin $O/ || fail table
[ "$(md5sum < $O/rl10_s.bin | cut -d' ' -f1)" = 1be4531e5c446ad70e3640b417fed10e ] || fail "table md5"
T1ENV="CRUCIBLE_MANIFOLD_TABLE=$O/rl10_s.bin CRUCIBLE_T1_LEAD=20e-3 CRUCIBLE_T1_THREADS=4"
# Chain A.
( wait_load; cd $O && env $T1ENV CRUCIBLE_T1_END=0.5e-3 CRUCIBLE_T1_OUT=fs \
    nice -n 10 /usr/bin/time -v $N/build/crucible_rl10_tier1 run development > fs.log 2>&1
  echo $? > t1_exit ) &
A=$!
# Chain B. Check 1, D0's three short cases (as bh_sd3.sh), new and parent.
for v in n p; do
  X=$N/build; [ $v = p ] && X=$P/build
  mkdir -p $O/d0$v && cd $O/d0$v
  wait_load
  CRUCIBLE_EQ_TABLE=$TA nice -n 10 $X/crucible_amr_budget 1e-3 4 amr dump=$O/d0$v/amr_dump.txt > amr.log 2>&1
  wait_load
  CRUCIBLE_EQ_TABLE=$TA CRUCIBLE_WALL_FUNCTIONS=law CRUCIBLE_SST_CORRECTION=hp CRUCIBLE_STATE_DUMP=$O/d0$v/c1_state.bin \
    nice -n 10 $X/crucible_chamber_study eqtt 64 12 1e-3 4 $O/d0$v/c1 > c1.log 2>&1
  nice -n 10 $X/crucible_amr_subcycling_tests c3 4 dump=$O/d0$v/c3 > c3.log 2>&1; echo $? > c3_exit
done
cd $O && { cmp d0n/amr_dump.txt d0p/amr_dump.txt && echo "fixed pattern dump identical"
  cmp d0n/c1_state.bin d0p/c1_state.bin && echo "C1 law hp state identical"
  cmp d0n/c1_history.csv d0p/c1_history.csv && echo "C1 law hp history identical"
  cmp d0n/c3_c3.bin d0p/c3_c3.bin && echo "check 3 dump identical"; } > $O/d0.txt 2>&1
grep -h "reflux's rho k clips" d0n/amr.log >> $O/d0.txt
echo done > $O/d0_done
# Check 3's trace: the new tree with print-only lines at each reflux clip (its rho k against the level's largest).
T=$N.trace; rm -rf $T; mkdir -p $T
rsync -a --exclude build --exclude out $N/ $T/
perl -0pi -e 's/#include <cstdio>\n//; s/#include <cmath>\n/#include <cmath>\n#include <cstdio>\n/; s/(if\(nt\) for\(auto u:touched\) if\(c\.turbulence_\[u\*nt\]<0\) \{\n)/$1        { double big=0;for(std::size_t q=0;q<cm.cells.size();++q) if(c.isActive(q)) big=std::max(big,c.turbulence_[q*nt]);\n          std::fprintf(stderr,"trace clip t %.9g onto nz %d i %zu j %zu rhok %.6e level max %.6e ratio %.3e\\n",time_,cm.nz,u\/static_cast<std::size_t>(cm.nr),u%static_cast<std::size_t>(cm.nr),c.turbulence_[u*nt],big,-c.turbulence_[u*nt]\/big); }\n/' $T/core/amr.cpp
grep -c "trace clip" $T/core/amr.cpp > $O/trace_patch.txt
TARGETS=crucible_rl10_tier1; wait_load
if build $T trace; then
  mkdir -p $O/trace && cd $O/trace && wait_load
  env $T1ENV CRUCIBLE_T1_END=1e-4 CRUCIBLE_T1_OUT=tr nice -n 10 $T/build/crucible_rl10_tier1 run development > tr.log 2> clips.txt
  echo $? > trace_exit
else echo "trace build failed" > $O/trace_patch.txt; fi
# Check 2: D3 on D1's case to 1 ms, the new tree, 4 then 1 thread.
cd $O; wait_load
CRUCIBLE_EQ_TABLE=$TA nice -n 10 $N/build/crucible_amr_budget 1e-3 4 tags wall dump=$O/d3_t4.txt > $O/d3_t4.log 2>&1
wait_load
OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 CRUCIBLE_EQ_TABLE=$TA \
  nice -n 10 $N/build/crucible_amr_budget 1e-3 1 tags wall dump=$O/d3_t1.txt > $O/d3_t1.log 2>&1
{ cmp $O/d3_t4.txt $O/d3_t1.txt && echo "D3: 1 and 4 threads identical to 1 ms"; grep -h "reflux's rho k clips" $O/d3_t4.log $O/d3_t1.log; } > $O/d3.txt 2>&1
echo done > $O/d3_done
# Check 1's D1 to 8 ms, new and parent together, once tier 1 is done.
wait $A; echo done > $O/t1_done
wait_load
( CRUCIBLE_EQ_TABLE=$TA nice -n 10 $N/build/crucible_amr_budget 8e-3 4 tags wall dump=$O/d1n.txt > $O/d1n.log 2>&1; echo $? > $O/d1n_exit ) &
( CRUCIBLE_EQ_TABLE=$TA nice -n 10 $P/build/crucible_amr_budget 8e-3 4 tags wall dump=$O/d1p.txt > $O/d1p.log 2>&1; echo $? > $O/d1p_exit ) &
wait
{ cmp $O/d1n.txt $O/d1p.txt && echo "D1 tagged 8 ms dump identical"; grep -h "reflux's rho k clips" $O/d1n.log; } > $O/d1.txt 2>&1
date > $O/finished
echo done > $O/all_done
