#!/usr/bin/env bash
# test_k7.sh -- every test in SPEC_K7.md section "Tests", plus the extra checks
# requested for the k7 engine.  Each test prints PASS or FAIL.  Exit code is the
# number of failures.
#
#   1. planted decompositions for k7 = 0, 1, 2 (positive) and near-misses (negative)
#   2. candidate-count agreement with src/caseA2.c  (124 and 1314)
#   3. control-band end-to-end, k7 <= 2
#   4. 2000-random-m differential per stratum vs an independent Python 4-sum finder
#   5. isqrt6 unit test vs Python on 1e5 values
#   6. first-50 / last-50 candidate m-reconstruction in Python
#   7. Bloom filter: zero false negatives, measured false-positive rate vs ideal
set -u
cd "$(dirname "$0")"
ENG=./k7engine
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
FAIL=0
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}

ok(){ echo "PASS  $*"; }
no(){ echo "FAIL  $*"; FAIL=$((FAIL+1)); }

echo "=== build ==="
gcc -O3 -march=native -fopenmp -std=gnu11 -Wall -Wextra -Wno-unused-parameter \
    -o "$ENG" k7engine.c -lm 2>"$TMP/build.err"
if [ $? -ne 0 ]; then cat "$TMP/build.err"; no "build"; exit 1; fi
[ -s "$TMP/build.err" ] && { echo "--- compiler warnings ---"; cat "$TMP/build.err"; }
ok "build (gcc -O3 -march=native -fopenmp)"

########################################################################
echo "=== TEST 1: planted decompositions and near-misses ==="
# SPEC_K7.md: k7=2 bases 7*11,7*13,5,19 ; k7=1 bases 7*3,7*5,7*8,11 ;
#             k7=0 bases 7*2,7*3,7*5,7*9
python3 - "$TMP" <<'PY'
import sys
T=sys.argv[1]
sets={0:[14,21,35,63],1:[21,35,56,11],2:[77,91,5,19]}
lines=[];meta=[]
M42=42**6
for k,bs in sets.items():
    m=sum(b**6 for b in bs)
    assert m%7==k
    lines.append("%d 300"%m);            meta.append(("k7=%d planted"%k,"FOUND"," ".join(map(str,sorted(bs)))))
    lines.append("%d 300"%(m+1));        meta.append(("k7=%d near-miss m+1"%k,"NONE",""))
    lines.append("%d 300"%(m+M42));      meta.append(("k7=%d near-miss m+42^6 (same k2,k3,k7)"%k,"NONE",""))
    lines.append("%d 300"%(m-M42));      meta.append(("k7=%d near-miss m-42^6 (same k2,k3,k7)"%k,"NONE",""))
    # base above the bound B must not be found
    Bs=max(bs)-1
    lines.append("%d %d"%(m,Bs));        meta.append(("k7=%d planted with B=%d (below max base)"%(k,Bs),"NONE",""))
open(T+"/plant.in","w").write("\n".join(lines)+"\n")
open(T+"/plant.meta","w").write("\n".join("%s\t%s\t%s"%x for x in meta)+"\n")
PY
$ENG 0 0 --plantfile "$TMP/plant.in" 2>/dev/null > "$TMP/plant.out"
python3 - "$TMP" <<'PY' > "$TMP/plant.res"
import sys
T=sys.argv[1]
meta=[l.rstrip("\n").split("\t") for l in open(T+"/plant.meta")]
out={}
for l in open(T+"/plant.out"):
    p=l.split()
    out[int(p[0])]=(p[1]," ".join(p[2:]))
bad=0
for i,(name,want,bases) in enumerate(meta):
    got=out.get(i,("MISSING",""))
    okk = got[0]==want and (want!="FOUND" or got[1]==bases)
    print(("PASS  " if okk else "FAIL  ")+"TEST1 %s -> %s %s"%(name,got[0],got[1]))
    if not okk: bad+=1
sys.exit(1 if bad else 0)
PY
cat "$TMP/plant.res"
if grep -q FAIL "$TMP/plant.res"; then no "TEST 1 planted/near-miss"; else ok "TEST 1 planted/near-miss (all 15 cases)"; fi

# near-misses must be verified genuinely un-decomposable by the Python oracle too
python3 - <<'PY'
M42=42**6
sets={0:[14,21,35,63],1:[21,35,56,11],2:[77,91,5,19]}
def iroot6(n):
    if n<=0: return 0
    r=int(round(n**(1/6)))
    while r>0 and r**6>n: r-=1
    while (r+1)**6<=n: r+=1
    return r
def four(m,B):
    two={}
    for x in range(1,B+1):
        x6=x**6
        if x6>m: break
        for y in range(1,x+1):
            s=x6+y**6
            if s>m: break
            if s not in two or two[s]>x: two[s]=x
    hi=min(B,iroot6(m))
    for b4 in range(hi,0,-1):
        r=m-b4**6
        if 4*b4**6<m: break
        if r<3: continue
        h3=min(b4,iroot6(r))
        for b3 in range(h3,0,-1):
            r2=r-b3**6
            if 3*b3**6<r: break
            if r2<2: continue
            v=two.get(r2)
            if v is not None and v<=b3: return True
    return False
bad=0
for k,bs in sets.items():
    m=sum(b**6 for b in bs)
    for lab,mm,B in [("m",m,300),("m+1",m+1,300),("m+42^6",m+M42,300),("m-42^6",m-M42,300),("m,B=max-1",m,max(bs)-1)]:
        want = (lab in ("m",))
        got=four(mm,B)
        if got!=want:
            print("FAIL  TEST1-oracle k7=%d %s oracle=%s want=%s"%(k,lab,got,want)); bad+=1
print("PASS  TEST 1b Python oracle agrees the near-misses really have no decomposition" if not bad else "FAIL  TEST 1b")
PY

########################################################################
echo "=== TEST 2: candidate-count agreement with src/caseA2.c ==="
echo "  convention: caseA2 scans the CLOSED interval [fmin,fmax]; k7engine scans"
echo "  the HALF-OPEN interval (FMIN,FMAX].  They coincide whenever gcd(FMIN,42)>1,"
echo "  which holds for 700000, 730000 and 1000000 (all even), so the counts must match."
c1=$($ENG 700000 730000 2>/dev/null | tr ' ' '\n' | sed -n 's/^candidates=//p')
c2=$($ENG 730000 1000000 2>/dev/null | tr ' ' '\n' | sed -n 's/^candidates=//p')
[ "$c1" = "124" ]  && ok "TEST 2 (700000,730000]  candidates=$c1 (expected 124)"  || no "TEST 2 (700000,730000]  candidates=$c1 expected 124"
[ "$c2" = "1314" ] && ok "TEST 2 (730000,1000000] candidates=$c2 (expected 1314)" || no "TEST 2 (730000,1000000] candidates=$c2 expected 1314"
# and a third, independent Python re-derivation of the 144 roots and the counts
python3 - <<'PY'
M=42**6
def pow6m(x,m): return pow(x,6,m)
r2=[x for x in range(1,64) if x%2 and pow6m(x,64)==1]
r3=[x for x in range(1,729) if x%3 and pow6m(x,729)==1]
r7=[x for x in range(1,117649) if x%7 and pow6m(x,117649)==1]
assert (len(r2),len(r3),len(r7))==(4,6,6)
roots=[]
for a in r2:
    for b in r3:
        x12=next(t for t in range(a,64*729,64) if t%729==b)
        for c in r7:
            x=next(t for t in range(x12,64*729*117649,64*729) if t%117649==c)
            roots.append(x)
assert len(roots)==144
def count(lo,hi):   # half-open (lo,hi]
    n=0
    for f in range(lo+1,hi+1):
        if f%2==0 or f%3==0 or f%7==0: continue
        for u in roots:
            t=u*f%M
            if 0<t<f: n+=1
    return n
a=count(700000,730000); b=count(730000,1000000)
print(("PASS  " if a==124 else "FAIL  ")+"TEST 2b independent Python (700000,730000]=%d"%a)
print(("PASS  " if b==1314 else "FAIL  ")+"TEST 2b independent Python (730000,1000000]=%d"%b)
PY

########################################################################
echo "=== TEST 3: control band 700000..730000 end to end, k7 <= 2 ==="
$ENG 700000 730000 --candidates "$TMP/cand_control.txt" 2>/dev/null | tee "$TMP/band.out"
rc=${PIPESTATUS[0]}
nsol=$(tr ' ' '\n' < "$TMP/band.out" | sed -n 's/^solutions=//p')
if [ "$nsol" = "0" ] && [ "$rc" = "0" ]; then ok "TEST 3 control band: 0 solutions, exit code 0"
else no "TEST 3 control band: solutions=$nsol exit=$rc"; fi

########################################################################
echo "=== TEST 4: 2000-random-m differential per stratum vs Python brute force ==="
python3 - "$TMP" <<'PY'
import random,sys
T=sys.argv[1]
random.seed(20260911)
B=300
M42=42**6
cases=[]   # (m, B, label)
for k7 in (0,1,2):
    for _ in range(2000):
        # exactly k7 bases coprime to 7, the rest divisible by 7
        co=[random.randrange(1,B+1) for _ in range(k7)]
        co=[c+1 if c%7==0 else c for c in co]
        co=[c if c<=B else c-7 for c in co]
        div=[7*random.randrange(1,B//7+1) for _ in range(4-k7)]
        bs=co+div
        assert sum(1 for b in bs if b%7) == k7, (bs,k7)
        m=sum(b**6 for b in bs)
        cases.append((m,B,"k7=%d built"%k7))
        # a near-miss in the SAME stratum: shifting by 42^6 preserves m mod 8,9,7
        cases.append((m+M42*random.randrange(1,50),B,"k7=%d perturbed"%k7))
open(T+"/diff.in","w").write("\n".join("%d %d"%(m,b) for m,b,_ in cases)+"\n")
open(T+"/diff.lab","w").write("\n".join(l for _,_,l in cases)+"\n")
PY
$ENG 0 0 --plantfile "$TMP/diff.in" 2>/dev/null > "$TMP/diff.out"
python3 - "$TMP" <<'PY'
import sys
T=sys.argv[1]
B=300
def iroot6(n):
    if n<=0: return 0
    r=int(n**(1/6))
    while r>0 and r**6>n: r-=1
    while (r+1)**6<=n: r+=1
    return r
# independent brute-force 4-sum oracle: b1<=b2<=b3<=b4<=B, all >=1
TWO={}
for x in range(1,B+1):
    x6=x**6
    for y in range(1,x+1):
        s=x6+y**6
        if s not in TWO or TWO[s]>x: TWO[s]=x
def four(m):
    hi=min(B,iroot6(m))
    for b4 in range(hi,0,-1):
        if 4*b4**6<m: break
        r=m-b4**6
        if r<3: continue
        h3=min(b4,iroot6(r))
        for b3 in range(h3,0,-1):
            if 3*b3**6<r: break
            r2=r-b3**6
            if r2<2: continue
            v=TWO.get(r2)
            if v is not None and v<=b3: return True
    return False
cases=[l.split() for l in open(T+"/diff.in")]
labs=[l.strip() for l in open(T+"/diff.lab")]
eng={}
for l in open(T+"/diff.out"):
    p=l.split(); eng[int(p[0])]=(p[1]=="FOUND", [int(z) for z in p[2:]])
bad=0; npos=0; nneg=0
from collections import Counter
cnt=Counter()
for i,(ms,bs) in enumerate(cases):
    m=int(ms)
    e_found,e_bases=eng[i]
    if e_found:
        assert sum(b**6 for b in e_bases)==m, "engine returned a wrong decomposition at %d"%i
        assert all(1<=b<=B for b in e_bases), "engine base out of bound at %d"%i
    p_found=four(m)
    if e_found!=p_found:
        bad+=1
        if bad<10: print("FAIL  TEST4 mismatch i=%d %s m=%d engine=%s python=%s"%(i,labs[i],m,e_found,p_found))
    cnt[(labs[i],p_found)]+=1
    if p_found: npos+=1
    else: nneg+=1
for k in sorted(cnt): print("       %-16s decomposable=%-5s : %d"%(k[0],k[1],cnt[k]))
print("       positives=%d negatives=%d total=%d"%(npos,nneg,len(cases)))
print(("PASS  " if bad==0 else "FAIL  ")+"TEST 4 differential: %d/%d agree"%(len(cases)-bad,len(cases)))
PY

########################################################################
echo "=== TEST 5: isqrt6 unit test vs Python on 1e5 values ==="
python3 - "$TMP" <<'PY'
import random,sys
T=sys.argv[1]; random.seed(7)
vals=[]
for _ in range(40000): vals.append(random.randrange(0,1<<114))
for _ in range(20000):
    r=random.randrange(1,2642245); vals.append(r**6)
for _ in range(20000):
    r=random.randrange(2,2642245); vals.append(r**6-1)
for _ in range(20000):
    r=random.randrange(1,2642245); vals.append(r**6+1)
vals += [0,1,2,63,64,65,(1<<128)-1,2642245**6]
open(T+"/roots.in","w").write("\n".join(map(str,vals))+"\n")
PY
$ENG 0 0 --isqrt6 "$TMP/roots.in" > "$TMP/roots.out" 2>/dev/null
python3 - "$TMP" <<'PY'
import sys
T=sys.argv[1]
def iroot6(n):
    if n<=0: return 0
    r=int(n**(1/6))+2
    while r>0 and r**6>n: r-=1
    while (r+1)**6<=n: r+=1
    return r
vs=[int(x) for x in open(T+"/roots.in")]
gs=[int(x) for x in open(T+"/roots.out")]
assert len(vs)==len(gs), "line count mismatch %d %d"%(len(vs),len(gs))
bad=0
for v,g in zip(vs,gs):
    w=iroot6(v)
    if w!=g:
        bad+=1
        if bad<6: print("FAIL  TEST5 isqrt6(%d) engine=%d python=%d"%(v,g,w))
print(("PASS  " if bad==0 else "FAIL  ")+"TEST 5 isqrt6: %d/%d exact"%(len(vs)-bad,len(vs)))
PY

########################################################################
echo "=== TEST 6: first-50 / last-50 candidate m-reconstruction in Python ==="
$ENG 730000 1000000 --candidates "$TMP/cand2.txt" >/dev/null 2>&1
python3 - "$TMP" <<'PY'
import sys
T=sys.argv[1]
bad=0; tot=0
for fn in ("cand_control.txt","cand2.txt"):
    rows=[l.split() for l in open(T+"/"+fn) if l.startswith("CANDIDATE")]
    sel=rows[:50]+rows[-50:]
    for r in sel:
        f,t,k2,k3,k7,m=int(r[1]),int(r[2]),int(r[3]),int(r[4]),int(r[5]),int(r[6])
        tot+=1
        if m*42**6 + t**6 != f**6: print("FAIL  TEST6 m*42^6+t^6 != f^6 for f=%d t=%d"%(f,t)); bad+=1; continue
        if not (0<t<f): print("FAIL  TEST6 t out of range f=%d t=%d"%(f,t)); bad+=1; continue
        from math import gcd
        if gcd(f,42)!=1 or gcd(t,42)!=1: print("FAIL  TEST6 gcd f=%d t=%d"%(f,t)); bad+=1; continue
        if (m%8,m%9,m%7)!=(k2,k3,k7): print("FAIL  TEST6 budgets f=%d"%f); bad+=1; continue
        if pow(t*pow(f,-1,42**6),6,42**6)!=1: print("FAIL  TEST6 t/f not a 6th root of unity f=%d"%f); bad+=1
print(("PASS  " if bad==0 else "FAIL  ")+"TEST 6 m-reconstruction: %d/%d candidates exact"%(tot-bad,tot))
PY

########################################################################
echo "=== TEST 7: Bloom filter health (no false negatives, FP rate vs theory) ==="
echo "  200,010,000 pair sums inserted at 16 bits/pair, 20,000,000 non-member queries."
echo "  Reference is the BLOCKED-Bloom expectation for 512-bit lines, k=8, 16 bpp"
echo "  (expectation over Poisson(32) items per line of (1-(1-1/512)^(8n))^8 = 8.7281e-4),"
echo "  NOT the non-blocked ideal 5.7450e-4: the 512-bit blocking alone costs 1.52x."
$ENG 0 0 --bloomstat 20000 20000000 2>/dev/null > "$TMP/bs.out"
cat "$TMP/bs.out"
python3 - "$TMP" <<'BLOOMPY'
import sys
T=sys.argv[1]
d=dict(kv.split("=",1) for kv in open(T+"/bs.out").read().split() if "=" in kv)
fn=int(d["false_negatives"]); rate=float(d["fp_rate"])
blocked=8.728125e-4
print("       false negatives = %d"%fn)
print("       measured FP rate = %.6g   blocked-Bloom theory = %.6g   ratio = %.3fx"%(rate,blocked,rate/blocked))
okk = (fn==0) and (rate < 1.10*blocked)
print(("PASS  " if okk else "FAIL  ")+"TEST 7 Bloom: zero false negatives and FP rate within 1.10x of blocked-Bloom theory")
BLOOMPY

########################################################################
echo
if [ "$FAIL" = "0" ]; then echo "SUMMARY: no shell-level failures; check every PASS/FAIL line above."
else echo "SUMMARY: $FAIL shell-level failures"; fi
exit $FAIL
