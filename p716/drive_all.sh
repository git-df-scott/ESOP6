#!/bin/bash
# Cover (1500,3000] in 250-wide chunks, then (3000,5000] in 500-wide chunks.
# Chunk width trades checkpoint granularity against the cost of rebuilding the
# 3-sum table per chunk (the build is ~TB^3/6 inserts however many passes it is
# split into, so it is amortised over the chunk's leaf work).
cd "$(dirname "$0")" || exit 1
./run716.sh 1500 3000 250 7 || exit $?
./run716.sh 3000 5000 500 7 || exit $?
echo "ALL RANGES DONE"
