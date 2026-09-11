#!/bin/bash
# run716.sh -- production driver for the (7,1,6) search.
#
#   ./run716.sh FMIN FMAX CHUNK Q
#
# Coverage is advanced one CHUNK of f at a time.  Each chunk is a SEPARATE
# engine invocation over (lo,hi], which internally runs Q residue passes
# (pass r keeps only 3-sums == r (mod Q) and only the c with R3 == r (mod Q)),
# so the chunk is fully covered when the invocation exits 0.  Only then is the
# chunk recorded in coverage.txt.  This is the checkpoint: a container restart
# loses at most the chunk in flight, and re-running the driver resumes from
# the last recorded chunk.
#
# Q must be chosen so the per-pass Bloom fits comfortably in RAM; the engine
# prints its size on stderr.  Q=7 is exactly uniform because x^7 == x (mod 7).
set -u
cd "$(dirname "$0")" || exit 1
FMIN=${1:?FMIN}; FMAX=${2:?FMAX}; CHUNK=${3:-500}; Q=${4:-7}
COV=runs/coverage.txt
LOG=runs/prod.log
mkdir -p runs
touch "$COV"

# resume: start above the highest hi already recorded, if it is inside our range
LAST=$(awk '$1=="COVERED"{print $3}' "$COV" | sort -n | tail -1)
if [ -n "${LAST:-}" ] && [ "$LAST" -gt "$FMIN" ]; then
  echo "resuming: coverage.txt already records up to f=$LAST" | tee -a "$LOG"
  FMIN=$LAST
fi

lo=$FMIN
while [ "$lo" -lt "$FMAX" ]; do
  hi=$((lo+CHUNK)); [ "$hi" -gt "$FMAX" ] && hi=$FMAX
  echo "=== chunk ($lo,$hi]  Q=$Q  $(date -u +%FT%TZ)" | tee -a "$LOG"
  t0=$(date +%s)
  ./engine716 "$lo" "$hi" --q "$Q" --threads 4 >> "$LOG" 2>> "$LOG"
  rc=$?
  t1=$(date +%s)
  if [ $rc -eq 3 ]; then
    echo "!!! SOLUTION -- engine exited 3 on chunk ($lo,$hi]" | tee -a "$LOG"
    grep '^SOLUTION' "$LOG" | tail -1 > SOLUTION.txt
    cat SOLUTION.txt
    exit 3
  fi
  if [ $rc -ne 0 ]; then
    echo "ERROR rc=$rc on chunk ($lo,$hi] -- stopping" | tee -a "$LOG"
    exit $rc
  fi
  # last BAND line of this invocation carries the chunk totals
  band=$(grep '^BAND' "$LOG" | tail -1)
  lv=$(echo "$band" | awk '{print $4}')
  sl=$(echo "$band" | awk '{print $9}')
  echo "COVERED $lo $hi leaves=$lv solutions=$sl wall=$((t1-t0))s Q=$Q $(date -u +%FT%TZ)" | tee -a "$COV" | tee -a "$LOG"
  python3 ./covlog.py || echo "covlog.py failed" | tee -a "$LOG"
  lo=$hi
done
echo "DONE ($FMIN,$FMAX]" | tee -a "$LOG"
