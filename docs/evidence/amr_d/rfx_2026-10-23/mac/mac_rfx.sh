#!/bin/sh
# The reflux clip's checks 1 (D0's three short cases) and 3's trace on the Mac, while backhouse's chain B waits on
# its load gate: as amr_d/bh_rfx.sh's chain B, new = worktree-inquiry2 HEAD (the clip, rebased), parent = fcb20f1.
W=/Users/ben/inquiry-project/.claude/worktrees/inquiry2; M=/tmp/inq2c/mac_rfx; P=$M/par; T=$M/trace; O=$M/out
TA=$M/table_a.bin; TS=/tmp/rl10_s.bin
mkdir -p $O $P $T
date > $O/started
cd $W && git rev-parse HEAD > $O/revision_new.txt && git archive fcb20f1 | tar -x -C $P || { echo archive > $O/failed; exit 1; }
rsync -a --exclude build --exclude .git $W/ $T/
perl -0pi -e 's/#include <cstdio>\n//; s/#include <cmath>\n/#include <cmath>\n#include <cstdio>\n/; s/(if\(nt\) for\(auto u:touched\) if\(c\.turbulence_\[u\*nt\]<0\) \{\n)/$1        { double big=0;for(std::size_t q=0;q<cm.cells.size();++q) if(c.isActive(q)) big=std::max(big,c.turbulence_[q*nt]);\n          std::fprintf(stderr,"trace clip t %.9g onto nz %d i %zu j %zu rhok %.6e level max %.6e ratio %.3e\\n",time_,cm.nz,u\/static_cast<std::size_t>(cm.nr),u%static_cast<std::size_t>(cm.nr),c.turbulence_[u*nt],big,-c.turbulence_[u*nt]\/big); }\n/' $T/core/amr.cpp
grep -c "trace clip" $T/core/amr.cpp > $O/trace_patch.txt
build() {  # $1 tree, $2 label, then targets
  t=$1; l=$2; shift 2
  mkdir -p $t/build && cd $t/build || return 1
  cmake .. -DCMAKE_BUILD_TYPE=Release -DCANTERA_ROOT=/Users/ben/src/cantera -DBOOST_INCLUDE=/opt/homebrew/include \
    -DCRUCIBLE_BUILD_APP=OFF -DCRUCIBLE_WITH_CANTERA=ON > $O/configure_$l.log 2>&1 || return 1
  cmake --build . -j 3 --target "$@" > $O/build_$l.log 2>&1
}
D0T="crucible_amr_budget crucible_chamber_study crucible_amr_subcycling_tests"
build $W n $D0T || { echo build_n > $O/failed; exit 1; }
build $P p $D0T || { echo build_p > $O/failed; exit 1; }
build $T trace crucible_rl10_tier1 || { echo build_trace > $O/failed; exit 1; }
date > $O/built
for v in n p; do
  X=$W/build; [ $v = p ] && X=$P/build
  mkdir -p $O/d0$v && cd $O/d0$v
  CRUCIBLE_EQ_TABLE=$TA $X/crucible_amr_budget 1e-3 4 amr dump=$O/d0$v/amr_dump.txt > amr.log 2>&1; echo $? > amr_exit
  CRUCIBLE_EQ_TABLE=$TA CRUCIBLE_WALL_FUNCTIONS=law CRUCIBLE_SST_CORRECTION=hp CRUCIBLE_STATE_DUMP=$O/d0$v/c1_state.bin \
    $X/crucible_chamber_study eqtt 64 12 1e-3 4 $O/d0$v/c1 > c1.log 2>&1; echo $? > c1_exit
  $X/crucible_amr_subcycling_tests c3 4 dump=$O/d0$v/c3 > c3.log 2>&1; echo $? > c3_exit
done
cd $O && { cmp d0n/amr_dump.txt d0p/amr_dump.txt && echo "fixed pattern dump identical"
  cmp d0n/c1_state.bin d0p/c1_state.bin && echo "C1 law hp state identical"
  cmp d0n/c1_history.csv d0p/c1_history.csv && echo "C1 law hp history identical"
  cmp d0n/c3_c3.bin d0p/c3_c3.bin && echo "check 3 dump identical"; } > $O/d0.txt 2>&1
grep -h "reflux's rho k clips" d0n/amr.log >> $O/d0.txt
date > $O/d0_done
mkdir -p $O/trace && cd $O/trace
CRUCIBLE_MANIFOLD_TABLE=$TS CRUCIBLE_T1_LEAD=20e-3 CRUCIBLE_T1_THREADS=4 CRUCIBLE_T1_END=1e-4 CRUCIBLE_T1_OUT=tr \
  /usr/bin/time -l $T/build/crucible_rl10_tier1 run development > tr.log 2> clips.txt
echo $? > trace_exit
date > $O/all_done
