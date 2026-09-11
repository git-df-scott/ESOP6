#!/bin/bash
# test716.sh -- every test mandated by strata/SPEC_716.md for the (7,1,6) engine.
# Each test prints PASS or FAIL; exit status is the number of failures.
cd "$(dirname "$0")" || exit 1
E=./engine716
FAIL=0
pass(){ echo "PASS  $*"; }
fail(){ echo "FAIL  $*"; FAIL=$((FAIL+1)); }

echo "== build =="
gcc -O3 -march=native -fopenmp -std=gnu11 -o engine716 engine716.c -lm || { fail "build"; exit 1; }
pass "build"

echo
echo "== T1: 20 planted solutions (random bases <= 300, incl. equal bases and g=1) =="
python3 - > /tmp/p716_plants.txt <<'PY'
import random
random.seed(71620260911)
plants=[]
# 2 forced special shapes first, then random
b=sorted([random.randint(1,300) for _ in range(5)],reverse=True)
plants.append(b+[1])                                   # g = 1
v=random.randint(1,300); w=random.randint(1,300)
plants.append(sorted([v,v,w,w,random.randint(1,300),random.randint(1,300)],reverse=True))  # equal bases
plants.append(sorted([200,200,200,7,7,1],reverse=True))# many equal bases and g=1
while len(plants)<20:
    plants.append(sorted([random.randint(1,300) for _ in range(6)],reverse=True))
for p in plants:
    print(sum(x**7 for x in p), " ".join(map(str,p)))
PY
n=0
while read -r S BASES; do
  n=$((n+1))
  out=$($E 2 300 --plant "$S" 2>/dev/null)
  rc=$?
  line=$(echo "$out" | grep '^SOLUTION' | head -1)
  if [ $rc -ne 3 ] || [ -z "$line" ]; then
    fail "plant $n (S=$S, planted: $BASES) not found (rc=$rc)"; continue
  fi
  outq=$($E 2 300 --q 16 --plant "$S" 2>/dev/null); rcq=$?
  lineq=$(echo "$outq" | grep '^SOLUTION' | head -1)
  if [ $rcq -ne 3 ] || [ "$lineq" != "$line" ]; then
    fail "plant $n: residue bucketing Q=16 did not reproduce it ('$lineq' vs '$line')"; continue
  fi
  set -- $line          # SOLUTION f a b c d e g
  got="$3 $4 $5 $6 $7 $8"
  if ! python3 ./verify7.py --sum "$S" $got > /tmp/p716_v.txt 2>&1 || ! grep -q "EXACT True" /tmp/p716_v.txt; then
    fail "plant $n: reported bases $got do NOT sum to S=$S"; continue
  fi
  if [ "$got" != "$BASES" ]; then
    echo "      note: plant $n returned a different (also exact) decomposition: $got vs $BASES"
  fi
  # near miss S+1 must be silent
  out2=$($E 2 300 --plant "$((S+1))" 2>/dev/null); rc2=$?
  if [ $rc2 -eq 3 ] || echo "$out2" | grep -q '^SOLUTION'; then
    fail "plant $n: near-miss S+1 produced a solution"; continue
  fi
  out2q=$($E 2 300 --q 16 --plant "$((S+1))" 2>/dev/null)
  if echo "$out2q" | grep -q '^SOLUTION'; then fail "plant $n: near-miss S+1 found under Q=16"; continue; fi
  pass "plant $n: found $got (Q=1 and Q=16), S+1 silent"
done < /tmp/p716_plants.txt
[ "$n" -eq 20 ] || fail "expected 20 plants, got $n"

echo
echo "== T2: (7,1,7) control 568^7 = 525^7+439^7+430^7+413^7+266^7+258^7+127^7 =="
if python3 ./verify7.py 568 525 439 430 413 266 258 127 | grep -q "EXACT True"; then
  pass "control identity verified in Python big ints"
else
  fail "control identity does not hold in Python"
fi
V=$(python3 -c 'print(266**7+258**7+127**7)')
q=$($E 2 300 --query3 "$V" 2>/dev/null | grep '^QUERY3')
echo "      $q"
if echo "$q" | grep -q "bloom=1 exact=1 266 258 127"; then
  pass "3-sum 266^7+258^7+127^7 present in the table and exactly recovered"
else
  fail "3-sum 266^7+258^7+127^7 not found in table ($q)"
fi

echo
echo "== T3: differential on (2,300] vs independent Python 4+3 brute force =="
$E 2 300 --threads 4 2>/dev/null | tail -1 > /tmp/p716_eng300.txt
# BAND fmin fmax leaves masked bloom_queries positives verified solutions elapsed
ELEAVES=$(awk '{print $4}' /tmp/p716_eng300.txt); ESOLS=$(awk '{print $9}' /tmp/p716_eng300.txt)
if [ ! -s /tmp/p716_ref300.txt ]; then python3 ./ref716.py 2 300 > /tmp/p716_ref300.txt; fi
RLEAVES=$(grep '^REF' /tmp/p716_ref300.txt | awk '{print $4}')
RSOLS=$(grep '^REF' /tmp/p716_ref300.txt | awk '{print $5}')
echo "      engine leaves=$ELEAVES solutions=$ESOLS ; reference leaves=$RLEAVES solutions=$RSOLS"
if [ "$ELEAVES" = "$RLEAVES" ] && [ "$ESOLS" = "$RSOLS" ] && [ "$ESOLS" = "0" ]; then
  pass "identical leaf counts ($ELEAVES) and identical zero verdicts"
else
  fail "differential mismatch (engine $ELEAVES/$ESOLS vs reference $RLEAVES/$RSOLS)"
fi

echo
echo "== T4: NB=1 vs NB=4 on (2,600] must agree exactly =="
$E 2 600 --nb 1 --threads 4 2>/dev/null | tail -1 > /tmp/p716_nb1.txt
$E 2 600 --nb 4 --threads 4 2>/dev/null | tail -1 > /tmp/p716_nb4.txt
# fields: 1 BAND 2 fmin 3 fmax 4 leaves 5 masked 6 bloom_queries 7 positives 8 verified 9 solutions
A=$(awk '{print $4,$5,$6,$9}' /tmp/p716_nb1.txt); B=$(awk '{print $4,$5,$6,$9}' /tmp/p716_nb4.txt)
PA=$(awk '{print $7}' /tmp/p716_nb1.txt);        PB=$(awk '{print $7}' /tmp/p716_nb4.txt)
VA=$(awk '{print $8}' /tmp/p716_nb1.txt);        VB=$(awk '{print $8}' /tmp/p716_nb4.txt)
echo "      nb=1: leaves/masked/queries/solutions = $A   positives=$PA verified=$VA"
echo "      nb=4: leaves/masked/queries/solutions = $B   positives=$PB verified=$VB"
if [ "$A" = "$B" ]; then
  pass "leaves, masked, bloom_queries and solutions identical under NB=1 and NB=4"
else
  fail "NB=1 and NB=4 disagree on a deterministic count"
fi
# positives/verified are NOT deterministic across NB: the per-pass filter is a
# different random object (NB times fewer entries in NB times fewer bits), so the
# SET of Bloom false positives differs while its RATE is the same.  Every positive
# is settled by the exact verifier, so no verdict depends on this.  Require the two
# positive counts to agree within 5 sigma of the Poisson fluctuation.
if python3 -c "
import sys,math
a,b=$PA,$PB
s=math.sqrt(a+b) if a+b else 1.0
sys.exit(0 if abs(a-b)<=5*s else 1)"; then
  pass "positives agree within Poisson noise ($PA vs $PB): same false-positive rate, different filter geometry"
else
  fail "positives differ by more than Poisson noise ($PA vs $PB)"
fi

echo
echo "== T5: Bloom false-negative test, 10^6 random table entries =="
bt=$($E 2 600 --bloomtest 1000000 --threads 4 2>/dev/null | grep '^BLOOMTEST')
echo "      $bt"
if echo "$bt" | grep -q "false_negatives=0"; then pass "no false negatives"
else fail "$bt"; fi
bs=$($E 2 600 --bloomstat 2000000 --threads 4 2>/dev/null | grep '^BLOOMSTAT')
echo "      $bs (informational: measured false-positive rate)"

echo
echo "== T6: residue bucketing Q=16 -- 16 passes must sum to the Q=1 counts on (2,600] =="
$E 2 600 --q 16 --threads 4 2>/dev/null | tail -1 > /tmp/p716_q16.txt
A=$(awk '{print $4,$5,$6,$9}' /tmp/p716_nb1.txt); B=$(awk '{print $4,$5,$6,$9}' /tmp/p716_q16.txt)
PA=$(awk '{print $7}' /tmp/p716_nb1.txt); PB=$(awk '{print $7}' /tmp/p716_q16.txt)
VA=$(awk '{print $8}' /tmp/p716_nb1.txt); VB=$(awk '{print $8}' /tmp/p716_q16.txt)
echo "      Q=1 : leaves/masked/queries/solutions = $A   positives=$PA verified=$VA"
echo "      Q=16: leaves/masked/queries/solutions = $B   positives=$PB verified=$VB"
if [ "$A" = "$B" ]; then
  pass "leaves, masked, bloom_queries and solutions identical under Q=1 and Q=16"
else
  fail "Q=16 pass sum differs from Q=1"
fi
# positives/verified are not comparable across bucketings: the Q=16 per-pass filter is
# sized for the LARGEST residue class, so the other 15 passes are underloaded and their
# false-positive rate is lower.  Every positive is settled exactly, so verdicts are
# unaffected; require only that Q=16 not exceed the Q=1 positive count materially.
if python3 -c "
import sys,math
a,b=$PA,$PB
sys.exit(0 if b<=a+5*math.sqrt(a+b+1) else 1)"; then
  pass "positives $PB (Q=16) <= $PA (Q=1) as expected from the lower filter load"
else
  fail "Q=16 produced more positives than Q=1 ($PB vs $PA)"
fi

echo
echo "== T7: residue bucketing, EXACT pass-sum identity on all five counters =="
# T6 compares leaves/masked/queries/solutions only, because `positives` and
# `verified` depend on the Bloom's geometry and the Q=16 per-pass filter is a
# different random object from the Q=1 filter (sized for the largest residue
# class, and holding 1/Q of the entries).  --nobloom removes that dependence
# entirely: every masked leaf is passed to the exact verifier, so all five
# counters become exact functions of the leaf set alone.  The pass sums must
# then match Q=1 digit for digit, for every Q.
REFN=$($E 2 600 --q 1 --nobloom --threads 4 2>/dev/null | tail -1 | awk '{print $4,$5,$6,$7,$8,$9}')
echo "      Q=1 (exact): leaves masked queries positives verified solutions = $REFN"
for q in 2 7 16 42; do
  GOTN=$($E 2 600 --q $q --nobloom --threads 4 2>/dev/null | tail -1 | awk '{print $4,$5,$6,$7,$8,$9}')
  if [ "$GOTN" = "$REFN" ]; then
    pass "Q=$q pass sums equal Q=1 exactly on all five counters"
  else
    fail "Q=$q exact pass-sum mismatch: $GOTN vs $REFN"
  fi
done

echo
if [ $FAIL -eq 0 ]; then echo "ALL TESTS PASS"; else echo "$FAIL TEST(S) FAILED"; fi
exit $FAIL
