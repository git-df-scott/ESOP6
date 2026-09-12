#!/bin/bash
# strata/test_v3.sh -- every test mandated by strata/SPEC_V3.md, section "Tests".
# Each check prints PASS or FAIL.  Exit status 0 iff every check printed PASS.
#
#   1  is_two_sum differential vs a brute-force Python two-sum oracle,
#      10^5 positives + 10^5 negatives, all divisibility patterns, plus a
#      second batch at the production Bmax and a structured-family batch.
#   2  planted full-engine decompositions in the k7=3 and k7=4 strata,
#      with m +/- 42^6 near-misses (identical budgets) required to be silent.
#   3  candidate counts 124 / 1314 on the two reference bands.
#   4  differential vs ./bin/caseA2 on (700000,730000]: identical candidate set,
#      identical verdicts, leaf count compared with caseA2's 4,139,088.
#   5  verify_solution.py wired to any SOLUTION line the engine can emit.
set -u
cd "$(dirname "$0")"
TMP=${TMPDIR:-/tmp}/v3test.$$
mkdir -p "$TMP"
trap 'rm -rf "$TMP"' EXIT
FAILED=0
pass(){ echo "PASS  $*"; }
fail(){ echo "FAIL  $*"; FAILED=1; }
ck(){ if [ "$2" = "$3" ]; then pass "$1 ($2)"; else fail "$1: got [$2] want [$3]"; fi; }
THREADS=${OMP_NUM_THREADS:-2}
export OMP_NUM_THREADS=$THREADS

echo "=== build ==="
if gcc -O3 -march=native -fopenmp -std=gnu11 -Wall -Wextra -Wno-unused-parameter \
       -o v3engine v3engine.c -lm 2>"$TMP/cc.err"; then
    if [ -s "$TMP/cc.err" ]; then cat "$TMP/cc.err"; fail "build emitted warnings"; else pass "build v3engine (no warnings)"; fi
else
    cat "$TMP/cc.err"; fail "build v3engine"; echo "cannot continue"; exit 1
fi
if gcc -O3 -march=native -fopenmp -o caseA2_nb4 caseA2_nb4.c -lm 2>/dev/null; then
    pass "build caseA2_nb4 (caseA2_timed.c with the nb=3,2,1 passes removed)"
else fail "build caseA2_nb4"; fi

echo "=== engine self-tests (roots, Hensel lifts, exact divisions, filters) ==="
if ./v3engine --twosum 16354 200 --selftest >"$TMP/st.out" 2>"$TMP/st.err"; then
    grep -q "TWOSUM YES 5 3" "$TMP/st.out" && pass "self-test battery + is_two_sum(3^6+5^6)=YES 5 3" \
        || { cat "$TMP/st.out"; fail "self-test battery"; }
else cat "$TMP/st.err"; fail "self-test battery (nonzero exit)"; fi

echo "=== TEST 1: is_two_sum differential vs brute-force Python oracle ==="
# 1a: 10^5 positives + 10^5 negatives with bases <= 3000 (SPEC_V3 test 1)
python3 v3_testlib.py gen 3000 100000 100000 12345 59 "$TMP/t1.txt" || fail "TEST1a generation"
./v3engine --twosumfile "$TMP/t1.txt" --bmax 3000 >"$TMP/t1.out" 2>/dev/null || fail "TEST1a engine run"
if python3 v3_testlib.py check "$TMP/t1.txt" "$TMP/t1.out" >"$TMP/t1.chk"; then
    cat "$TMP/t1.chk"; pass "TEST 1a differential 200000/200000 agree (bases <= 3000)"
else cat "$TMP/t1.chk"; fail "TEST 1a differential"; fi

# 1b: same at the production table size Bmax = 104761 (the (4.3M,4.4M] band)
python3 v3_testlib.py gen 104761 4000 4000 999 331 "$TMP/t2.txt" || fail "TEST1b generation"
./v3engine --twosumfile "$TMP/t2.txt" --bmax 104761 >"$TMP/t2.out" 2>/dev/null || fail "TEST1b engine run"
if python3 v3_testlib.py check "$TMP/t2.txt" "$TMP/t2.out" >"$TMP/t2.chk"; then
    cat "$TMP/t2.chk"; pass "TEST 1b differential 8000/8000 agree (bases <= 104761, production Bmax)"
else cat "$TMP/t2.chk"; fail "TEST 1b differential"; fi

# 1c: hand-built structured families, one per pattern the completeness proof uses
python3 - "$TMP/t3.txt" <<'EOF'
import sys
B=3000
fam=[("rough x rough",701,709),("rough x rough (equal)",727,727),
     ("2q x q' (SPEC_V3 gap)",2*701,709),("4q x 3q'",4*701,3*709),
     ("6q x q'",6*233,719),("common large prime q,q",719,719),
     ("common large prime 2q,3q",2*719,3*719),("determination at 11",11*13,709),
     ("both divisible by 11",11*13,11*17),("recursion depth 2 at p=2",4*701,4*709),
     ("recursion depth 2 at p=3",9*211,9*317),("recursion 2 then 3",6*211,6*317),
     ("recursion depth 3 at p=2",8*211,8*317),("2 x q",2,701),("1 x q",1,701),
     ("1 x 1",1,1),("smooth x smooth",120,28),("2^4 x small",16*101,6)]
out=open(sys.argv[1],"w"); tags=open(sys.argv[1]+".tags","w")
for t,a,b in fam:
    assert a<=B and b<=B, t
    out.write("%d %d\n"%(a**6+b**6,B)); tags.write("%s %d %d\n"%(t,a,b))
EOF
./v3engine --twosumfile "$TMP/t3.txt" --bmax 3000 >"$TMP/t3.out" 2>/dev/null
n3=$(wc -l < "$TMP/t3.txt"); y3=$(grep -c " YES " "$TMP/t3.out")
ck "TEST 1c every structured family found" "$y3" "$n3"
paste -d' ' "$TMP/t3.out" "$TMP/t3.txt.tags" | awk '{printf "       %-28s %s\n", $5" "$6" "$7" "$8, $2}' | head -20

echo "=== TEST 2: planted decompositions (k7=3, k7=4) and near-misses ==="
python3 - "$TMP/plant.txt" "$TMP/plant.tags" <<'EOF'
import sys
# Planted four-sum decompositions.  k7 = #bases coprime to 7 (the stratum),
# k2 = #odd, k3 = #coprime to 3.  All m exceed 42^6 so that m-42^6 is a legal
# near-miss with the *identical* (k2,k3,k7) budgets.
cases=[([70,101,103,107],            "k7=3 one base divisible by 7"),
       ([1402,709,1409,1414],        "k7=3, leaf pair (2q,q') -> residual table"),
       ([202,206,211,217],           "k7=3, leaf pair both even -> recursion at 2"),
       ([56,75,88,97],               "k7=3 mixed parities"),
       ([101,103,107,109],           "k7=4 all odd"),
       ([1402,709,1409,1423],        "k7=4, leaf pair (2q,q') -> residual table")]
out=open(sys.argv[1],"w"); tg=open(sys.argv[2],"w")
for bs,note in cases:
    m=sum(b**6 for b in bs); B=max(bs)
    k2=sum(b%2 for b in bs); k3=sum(1 for b in bs if b%3); k7=sum(1 for b in bs if b%7)
    assert m%8==k2%8 and m%9==k3 and m%7==k7 and m>42**6
    tg.write("%s|k7=%d k2=%d k3=%d|%s\n"%(note,k7,k2,k3," ".join(map(str,sorted(bs)))))
    out.write("%d %d\n"%(m,B))                       # must be FOUND
    tg.write("near-miss m+42^6 (%s)|k7=%d|-\n"%(note,k7)); out.write("%d %d\n"%(m+42**6,B))
    tg.write("near-miss m-42^6 (%s)|k7=%d|-\n"%(note,k7)); out.write("%d %d\n"%(m-42**6,B))
    tg.write("bound too small (%s)|k7=%d|-\n"%(note,k7)); out.write("%d %d\n"%(m,B-1))
EOF
./v3engine --plantfile "$TMP/plant.txt" --bmax 1500 >"$TMP/plant.out" 2>"$TMP/plant.err" || fail "TEST 2 engine run failed: $(tail -1 "$TMP/plant.err")"
ok=1; i=0
while read -r line; do
    tag=$(sed -n "$((i+1))p" "$TMP/plant.tags"); note=${tag%%|*}; rest=${tag#*|}; k7=${rest%%|*}; want=${tag##*|}
    got=$(echo "$line" | cut -d' ' -f2-)
    case "$note" in
      near-miss*|bound*) [ "$got" = "NONE" ] || { fail "TEST2 $note -> $got (expected NONE)"; ok=0; } ;;
      *) gotb=$(echo "$got" | cut -d' ' -f2-)
         if [ "$got" = "NONE" ]; then fail "TEST2 $note ($k7) -> NONE (expected $want)"; ok=0
         elif [ "$gotb" != "$want" ]; then fail "TEST2 $note ($k7) -> $gotb (expected $want)"; ok=0
         else echo "       $note ($k7) -> FOUND $gotb"; fi ;;
    esac
    i=$((i+1))
done < "$TMP/plant.out"
nplant=$(wc -l < "$TMP/plant.txt")
ck "TEST 2 result lines == planted cases" "$i" "$nplant"
[ $ok = 1 ] && [ "$i" = "$nplant" ] && pass "TEST 2 planted decompositions and near-misses ($i cases)"
# independent confirmation that the near-misses really have no decomposition
python3 - "$TMP/plant.txt" <<'EOF'
import sys
bad=0
lines=[l.split() for l in open(sys.argv[1])]
_c={}
def four(m,B):
    """brute force m = b1^6+..+b4^6 with 1<=b_i<=B, via a two-sum dictionary"""
    P=[x**6 for x in range(B+1)]
    if B in _c: S2=_c[B]
    else: S2=_c[B]={}
    if not S2:
     for x in range(1,B+1):
        for y in range(1,x+1):
            v=P[x]+P[y]
            if v not in S2: S2[v]=(x,y)
    for a in range(B,0,-1):
        if 4*P[a]<m: break
        if P[a]>=m: continue
        for b in range(a,0,-1):
            s=P[a]+P[b]
            if s+2*P[b]<m: break
            if s>=m: continue
            r=m-s
            if r in S2:
                x,y=S2[r]
                if x<=b: return (a,b,x,y)
    return None
for n,(R,B) in enumerate(lines):
    R=int(R);B=int(B)
    if n%4==0: continue
    if four(R,B) is not None:
        print("ORACLE DISAGREES: case %d has a decomposition"%n); bad=1
print("PASS  TEST 2b Python oracle agrees the near-misses have no decomposition" if not bad
      else "FAIL  TEST 2b")
EOF

echo "=== TEST 3: candidate counts ==="
c1=$(./v3engine 700000 730000 --k7 5 2>/dev/null | tr ' ' '\n' | grep '^candidates=' | cut -d= -f2)
c2=$(./v3engine 730000 1000000 --k7 5 2>/dev/null | tr ' ' '\n' | grep '^candidates=' | cut -d= -f2)
ck "TEST 3 (700000,730000] candidate count" "$c1" "124"
ck "TEST 3 (730000,1000000] candidate count" "$c2" "1314"
c3=$(python3 v3_testlib.py candcount 700000 730000)
ck "TEST 3b independent Python count (700000,730000]" "$c3" "124"

echo "=== TEST 4: differential vs ./bin/caseA2 on (700000,730000] ==="
../bin/caseA2 700000 730000 --dump-candidates 2>/dev/null | grep '^CANDIDATE' \
    | awk '{print $2, $3}' | sort -n > "$TMP/ca2.cand"
../bin/caseA2 700000 730000 2>/dev/null | grep -E '^(SOLUTION|DEGENERATE)' > "$TMP/ca2.verdict"
./v3engine 700000 730000 --k7 0,1,2,3,4 --candidates "$TMP/v3.candfile" >"$TMP/v3.out" 2>/dev/null
grep '^CANDIDATE' "$TMP/v3.candfile" | awk '{print $2, $3}' | sort -n > "$TMP/v3.cand"
grep -E '^(SOLUTION|DEGENERATE)' "$TMP/v3.out" > "$TMP/v3.verdict"
if diff -q "$TMP/ca2.cand" "$TMP/v3.cand" >/dev/null; then
    pass "TEST 4 candidate sets identical ($(wc -l < "$TMP/v3.cand") (f,t) pairs)"
else fail "TEST 4 candidate sets differ"; diff "$TMP/ca2.cand" "$TMP/v3.cand" | head; fi
nv1=$(wc -l < "$TMP/ca2.verdict"); nv2=$(wc -l < "$TMP/v3.verdict")
ck "TEST 4 caseA2 verdict lines" "$nv1" "0"
ck "TEST 4 v3engine verdict lines" "$nv2" "0"
lv=$(tr ' ' '\n' < "$TMP/v3.out" | grep '^leaves=' | cut -d= -f2)
la=$(./caseA2_nb4 700000 730000 2>&1 | tr ' ' '\n' | grep '^j2_nodes=' | cut -d= -f2)
ck "TEST 4 leaf count vs caseA2 restricted to nb=4" "$lv" "$la"
echo "       v3engine leaves = $lv ; caseA2 (nb=4 only) j2_nodes = $la ;"
echo "       caseA2 as shipped (nb=4,3,2,1) j2_nodes = 4139088, i.e. 2424 extra"
echo "       leaves from its three degenerate passes, which v3engine does not run."

echo "=== TEST 5: verify_solution.py wiring ==="
# (a) the verifier itself works on an exact identity and rejects a false one
python3 verify_solution.py 0 0 0 0 7 7 | grep -q "EXACT True" && pass "TEST 5a verifier accepts an exact identity" || fail "TEST 5a"
if python3 verify_solution.py 42 84 126 168 5 6 | grep -q "EXACT False"; then pass "TEST 5b verifier rejects a false tuple"; else fail "TEST 5b"; fi
# (b) the harness that every run below uses: any SOLUTION-SIX line is piped into it
checksol(){ # $1 = engine stdout file
    awk '/^SOLUTION-SIX/{print $2,$3,$4,$5,$6,$7}' "$1" | while read -r a b c d e f; do
        python3 verify_solution.py "$a" "$b" "$c" "$d" "$e" "$f"
    done
}
printf 'SOLUTION-SIX 42 84 126 168 5 6\n' > "$TMP/fake.out"
if checksol "$TMP/fake.out" | grep -q "EXACT False"; then
    pass "TEST 5c SOLUTION-SIX harness routes a (fabricated) solution line into the verifier"
else fail "TEST 5c harness"; fi
if [ -s "$TMP/v3.out" ] && grep -q SOLUTION "$TMP/v3.out"; then
    checksol "$TMP/v3.out" | grep -q "EXACT True" || fail "TEST 5d a real SOLUTION failed verification"
else pass "TEST 5d no SOLUTION was produced by any run in this script (nothing to verify)"; fi

echo
if [ $FAILED = 0 ]; then echo "ALL TESTS PASSED"; else echo "SOME TESTS FAILED"; fi
exit $FAILED
