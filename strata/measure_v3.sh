#!/bin/bash
# strata/measure_v3.sh -- timing harness for strata/V3_LOG.md.
# Usage: measure_v3.sh <threads>
set -u
cd "$(dirname "$0")"
T=${1:-4}
export OMP_NUM_THREADS=$T
run(){ # $1 = label, rest = args
    local lab="$1"; shift
    echo "### $lab : ./v3engine $*"
    ./v3engine "$@" 2>&1 | grep -E "^(BAND|vtable|memory|masks)"
}
echo "== threads = $T =="
run "timing band, full engine"        4300000 4400000 --k7 3,4
run "timing band, masks only"         4300000 4400000 --k7 3,4 --nostep1 --notable
run "timing band, masks + prime loop" 4300000 4400000 --k7 3,4 --notable
run "timing band, direct residues"    4300000 4400000 --k7 3,4 --nores
for b in "6000000 6010000" "7000000 7010000" "8000000 8010000" "10000000 10010000"; do
    run "scaling band $b" $b --k7 3,4
done
