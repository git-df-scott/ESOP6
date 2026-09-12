#!/bin/bash
# strata/measure_v3.sh -- timing harness behind the numbers in strata/V3_LOG.md.
# Usage: measure_v3.sh [threads] [what]      what = split | scaling | all
set -u
cd "$(dirname "$0")"
T=${1:-4}
WHAT=${2:-all}
export OMP_NUM_THREADS=$T
run(){ local lab="$1"; shift
    echo "### $lab"
    echo "### ./v3engine $*"
    ./v3engine "$@" 2>&1 | grep -E "^(BAND|vtable|memory|masks)"
}
echo "== threads = $T =="
if [ "$WHAT" = all ] || [ "$WHAT" = split ]; then
  # cost split, measured on a narrower band so the three variants are cheap.
  # full - (masks+prime loop) = residual-table cost
  # (masks+prime loop) - masks = prime-loop cost
  # masks = leaf enumeration + the two composite 2-sum masks
  run "SPLIT full"              4300000 4320000 --k7 3,4
  run "SPLIT masks only"        4300000 4320000 --k7 3,4 --nostep1 --notable
  run "SPLIT masks+prime loop"  4300000 4320000 --k7 3,4 --notable
  run "SPLIT full, direct residues instead of the T_p tables" 4300000 4320000 --k7 3,4 --nores
fi
if [ "$WHAT" = all ] || [ "$WHAT" = scaling ]; then
  run "SCALE 6.000-6.010M"  6000000 6010000 --k7 3,4
  run "SCALE 7.000-7.008M"  7000000 7008000 --k7 3,4
  run "SCALE 8.000-8.005M"  8000000 8005000 --k7 3,4
  run "SCALE 10.000-10.003M" 10000000 10003000 --k7 3,4
fi
