# Stratified class-1 engine: the "k7 <= 2" stratum

## Setting (all facts verified independently by the auditor)

Primitive solution a^6+b^6+c^6+d^6+e^6 = f^6, class 1 (Meyrignac): one term t
carries all three exemptions (odd, coprime to 3, coprime to 7); the other four
are 42*b_i.  Then

    f^6 - t^6 = 42^6 * m,    m = b1^6 + b2^6 + b3^6 + b4^6,   1 <= b_i <= B := (f-1)/42.

t = zeta * f mod 42^6 for one of the 144 sixth roots of unity zeta mod 42^6,
with 0 < t < f, gcd(f,42)=1.

Budgets (exact, from x^6 in {0,1} mod 8, 9, 7):
    k2 := m mod 8 = #{i : b_i odd}         must be <= 4
    k3 := m mod 9 = #{i : 3 does not divide b_i}   must be <= 4
    k7 := m mod 7 = #{i : 7 does not divide b_i}   must be <= 4

## The stratum

We search exactly the candidates with k7 in {0, 1, 2}.  Heuristic mass of this
stratum within class 1: (1 + 24 + 216)/2401 ~= 10%.  Everything with k7 in {3,4}
is NOT covered and must be reported as not covered.

Why it is cheap: 7^6 = 117649 and B < 7^6 for all f < 4.94e9, so a base that is
forced to satisfy a congruence mod 7^6 is determined up to at most
6 * ceil(B/7^6) values instead of being enumerated.

## Cases

Let R7 = 7^6.  Precompute ROOTS[u] = sorted list of x in [1, R7) with x^6 == u (mod R7),
for every unit residue u (each list has 0 or 6 entries).  Also precompute
P6[x] = x^6 as unsigned __int128 for x <= B_max, and small residue tables.

### k7 = 0  (all four bases divisible by 7)
Require 7^6 | m (valuation lemma; otherwise no solution, skip).  m' = m / 7^6,
B' = B / 7.  Solve the GENERIC four-sum m' = sum b'_i^6, b'_i <= B', with the
recursive budgets recomputed from m' (k2' = m' mod 8, k3' = m' mod 9,
k7' = m' mod 7; each <= 4 else skip).  Use a 2+2 split: DFS over (b'_4 >= b'_3)
with residue masks, pair table query for the remainder.  These candidates are
astronomically rare (7^6 | m); cost is irrelevant, correctness is not.

### k7 = 1  (exactly one base coprime to 7; call it b4)
m mod R7 must be a unit sixth-power residue; otherwise skip.
For each r in ROOTS[m mod R7], for each b4 = r + j*R7 <= B (j >= 0):
    R = m - b4^6  (must be > 0 and divisible by R7 -- guaranteed by construction, assert it)
    R' = R / R7, B' = B/7.
    Need R' = b'_1^6 + b'_2^6 + b'_3^6, b'_i <= B'.
    Budgets for the three: k2' = R' mod 8 (<=3), k3' = R' mod 9 (<=3), k7' = R' mod 7 (<=3).
    Enumerate b'_3 from ceil((R'/3)^(1/6)) to min(floor(R'^(1/6)), B') (largest part),
    apply residue masks for 3-sums then 2-sums, query pair table for R' - b'_3^6,
    exact-verify any hit.

### k7 = 2  (exactly two bases coprime to 7; call them b3 <= b4)
For each b3 in [1, B] with gcd(b3,7)=1 (apply stride from k2/k3 only when the
budget is 0 or full, i.e. all-even / all-odd / all-div-3 / all-coprime-3):
    u = (m - b3^6) mod R7; if u == 0 or ROOTS[u] empty: continue
    for each r in ROOTS[u], for each b4 = r + j*R7 with b3 <= b4 <= B:
        R = m - b3^6 - b4^6; if R <= 0: continue
        assert R % R7 == 0; R' = R / R7; B' = B/7
        need R' = b'_1^6 + b'_2^6 with b'_i <= B'  (so 2 <= R' <= 2*B'^6)
        residue masks for 2-sums (mod 64, 27, 49, 13, 19, 31, 37, 43), then pair
        table query, then exact verify.

## Pair table
Sums x^6 + y^6 for 1 <= y <= x <= B'_max = floor((F_max-1)/42/7).
At F_max = 2.0e7: B'_max = 68027, ~2.31e9 pairs.  Blocked Bloom filter at 16 bits
per pair = 4.6 GB.  Any Bloom positive MUST be settled by an exact search
(loop x from ceil((R'/2)^(1/6)) to floor(R'^(1/6)), test R' - x^6 is a perfect
sixth power by integer sixth root with exact 128-bit confirmation).
Fix the known hash defect: derive the 8 in-line bit positions from independent
9-bit slices of a well-mixed 64-bit hash (e.g. two different 64-bit mixes of the
128-bit key, 9 bits each from non-overlapping positions).

## Arithmetic
All values as unsigned __int128.  m < (F/42)^6; at F = 2e7 this is < 2^114. Fine.
Compute m by exact 128-bit subtraction f^6 - t^6 then division by 42^6, and ASSERT
the remainder is zero (f^6 < 2^128 requires f < 2^21.3 = 2.6e6 -- NOT ENOUGH for
f = 2e7!).  So compute f^6 - t^6 via the factorisation
(f-t)(f+t)(f^2-ft+t^2)(f^2+ft+t^2) and strip 42^6 by gcd from the factors as the
existing caseA2.c does (caseA2.c:282-294), OR use the identity
m = ((f^6 - t^6)/42^6) computed with a small multiprecision routine.  Either way,
add a self-check: m * 42^6 + t^6 == f^6 verified in Python for the first 50 and
last 50 candidates of every run (dump them to a file).

## Integer sixth root
isqrt6(n) for 128-bit n: floating-point estimate then correct with exact
128-bit multiplications until r^6 <= n < (r+1)^6.  Unit-test against Python on
10^5 random values including exact sixth powers and sixth powers minus one.

## Output protocol
Print `CANDIDATE f t k2 k3 k7` per candidate (to a file, not stdout, when the
count is large), and `SOLUTION f t b1 b2 b3 b4` with the ORIGINAL bases
(multiply the primed ones back by 7) whenever an exact decomposition is found.
On any SOLUTION, immediately also print the six numbers
(42*b1, 42*b2, 42*b3, 42*b4, t, f) and exit non-zero so the wrapper runs
strata/verify_solution.py.

## Tests (must all pass before any production run)
1. Planted decompositions (test mode `--plant m B`): for each case k7=0,1,2 build
   m from chosen bases (e.g. k7=2: bases 7*11, 7*13, 5, 19 ; k7=1: 7*3,7*5,7*8, 11;
   k7=0: 7*2,7*3,7*5,7*9) and confirm the engine prints the decomposition.
   Also plant near-misses (m+1) and confirm silence.
2. Candidate enumeration equals the repo's for 700000..730000 (124) and
   730000..1000000 (1314) -- reuse the repo's counting convention (f <= fmax,
   f > fmin ... check caseA2.c and match it exactly, then state the convention).
3. End-to-end on the control band 700000..730000 restricted to k7<=2: report the
   count of candidates per k7 value, zero solutions, wall time.
4. Differential vs caseA2 on a synthetic 4-sum oracle: for 2000 random m built
   from random bases <= 300 in each stratum, both the new decomposer and an
   independent brute-force Python 4-sum finder must agree.

## Coverage statement to record after the run
"Class 1, k7 in {0,1,2}, F0 < f <= F1: N candidates, 0 solutions, wall time, host."
Never state more than that.
