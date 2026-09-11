# (7,1,6) engine: f^7 = a^7+b^7+c^7+d^7+e^7+g^7

Auditor's design, 2026-09-11.  A (7,1,6) solution is a counterexample to
Euler's sum-of-powers conjecture for k = 7, exactly as (5,1,4) was for k = 5.
Heuristic density per e-fold of f (unordered solutions) ~ 1.8e-3, versus
~2.6e-5 for (6,1,5): about 70x more promising per unit of height.  Status in
the literature: no (7,1,6) solution is known (MathWorld, Wikipedia); the
EulerNet search bound could not be retrieved (euler.free.fr is dead), so the
run must record its own coverage and nothing more.

## Algorithm: 4+3 meet in the middle with nested windows
Order the six bases a >= b >= c >= d >= e >= g.
Table: all sums d^7+e^7+g^7 with F >= d >= e >= g >= 1 (F^3/6 entries), in a
blocked Bloom filter (16 bits per entry; repaired hash from k7engine.c), plus an
EXACT verifier for any positive: given R, loop d from floor(R^(1/7)) down to
ceil((R/3)^(1/7)), then e from floor((R-d^7)^(1/7)) down to
ceil(((R-d^7)/2)^(1/7)), test R-d^7-e^7 is a seventh power by integer root.
Enumeration: for f in (FMIN, FMAX]:
  a in [ceil((f^7/6)^(1/7)), f-1]
  R1 = f^7 - a^7;  b in [ceil((R1/5)^(1/7)), min(a, floor(R1^(1/7)))]
  R2 = R1 - b^7;   c in [ceil((R2/4)^(1/7)), min(b, floor(R2^(1/7)))]
  R3 = R2 - c^7;   require 3 <= R3 <= 3*c^7 (so that d <= c is possible);
  masks: R3 mod 49, 29, 43 must be 3-sum residues (tables; carry residues
  incrementally, no 128-bit division); then Bloom; then exact verify.
Arithmetic: f^7 < 2^128 needs f < 2^18.28 = 318,000 -- fine for any reachable
F.  Use unsigned __int128 throughout; precompute P7[x] = x^7.
Seventh-power residues: mod 49 the units give {1,18,19,30,31,48}; mod 29 the
units give 4 values; mod 43 give 6 values.  Build the 3-sum residue bitmaps
by brute force at start-up and print their pass rates (expected 0.755, 0.862,
0.860).

## Symmetry and completeness
Every unordered solution is enumerated exactly once (a >= b >= c on the
enumerated side, d >= e >= g in the table, and the boundary c >= d is enforced
by R3 <= 3 c^7 plus the verifier's d <= ... hmm: the verifier must ALSO enforce
d <= c; pass c into the verifier and cap d at c).  Equal bases are allowed.

## Cost model
Leaves ~ kappa * F^4 with kappa measured on a narrow band; expect ~0.01-0.02.
Bloom query ~300-500 ns (cache miss); masks pass ~56%.  Table: F=2000 ->
1.33e9 entries (2.7 GB); F=2500 -> 2.6e9 (5.2 GB); F=3000 -> 4.5e9 (9 GB).
Build the table ONCE for FMAX and process f in increasing bands so partial
coverage is recorded as the run proceeds (print a BAND line every 50 f).

## Tests (mandatory)
1. Planted: pick random a..g <= 300, compute S = sum of the six seventh
   powers, and run the enumeration with f^7 replaced by S (test mode
   `--plant S`); must report the six bases.  Near-miss S+1 silent.
   20 random plants including cases with equal bases and with g = 1.
2. Known (7,1,7) solution as a control of the arithmetic:
   568^7 = 525^7+439^7+430^7+413^7+266^7+258^7+127^7 -- verify in Python and
   verify the engine's 3-sum table contains 266^7+258^7+127^7.
3. Differential on (2, 300]: run the engine and an independent brute-force
   Python 4+3 search with a dict; both must report 0 solutions and the same
   leaf count (or explain the difference exactly).
4. Bloom false-negative test: 10^6 random table entries must all query positive.
5. Any SOLUTION line goes through strata/verify_solution.py generalised to
   exponent 7 (add an optional exponent argument to that script).

## Output
`BAND fmin fmax leaves bloom_queries positives verified solutions elapsed` per
band; `SOLUTION f a b c d e g` on success, exit code 3.
