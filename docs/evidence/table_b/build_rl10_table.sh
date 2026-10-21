#!/bin/zsh
# Table B on the RL10's range at delta-eta 0.1 (TABLE_B.md, *Table B on the RL10's range*): the tier-1 central case's
# chemistry (docs/validation/rl10a33a/tier1_run_definition.md). The arguments and environment are /tmp/rl10s/build.sh's
# first build, verbatim but for the paths (built 21 October 2026, about 01:00 CDT, as /tmp/tbs/rl10_s.bin).
# Usage: build_rl10_table.sh <build dir> <output .bin> [node report .csv]
set -e
B=${1:?build dir}; OUT=${2:?output .bin}
SRC=$(cd "$(dirname "$0")/../../.." && pwd)
export CRUCIBLE_TABLE_B_CANDIDATES=3 CRUCIBLE_TABLE_B_FH_NODES=0.133333,0.142857,0.153846,0.166667,0.181818 CRUCIBLE_TABLE_B_WRITE_ANYWAY=1
[[ -n "$3" ]] && export CRUCIBLE_TABLE_B_NODE_REPORT="$3"
"$B/crucible_build_manifold_table" "$OUT" 0 0 0 0,0.3 11 9 1e-3 80 3 "$SRC/tools/mechanisms/h2o2_highT.yaml"
