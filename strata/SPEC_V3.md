# Engine v3: table-free class-1 search with complete coverage

Auditor's design, 2026-09-11. Complements SPEC_K7.md (which stays the cheapest
route for k7 <= 2). v3 is for k7 in {3,4} (80% of candidates by count, 90% of
class-1 mass) and, more generally, for ANY candidate.

## Outer structure (unchanged from caseA2)
Class-1 candidates (f, t, m), budgets (k2,k3,k7). DFS over the two largest
bases b4 >= b3 with the nested "largest part" windows
   b4 in [ceil((m/4)^(1/6)), floor(m^(1/6))],  b3 in [ceil((R4/3)^(1/6)), min(floor(R4^(1/6)), b4)]
and the exact stride/skip rules from the budgets, exactly as caseA2.c:dfs.
Leaf remainder R = m - b4^6 - b3^6, unknowns b1 <= b2 <= b3, budgets (o2,o3,o7)
for the pair (each 0..2).

## Leaf: is_two_sum(R, bound, o2, o3, o7)  -- NO PAIR TABLE
Step 0. Residue masks for 2-sums, computed INCREMENTALLY (never a 128-bit
   division in the hot path): keep m mod Q, and tables b^6 mod Q for b <= Bmax,
   for Q1 = 64*27*49 = 84672 and Q2 = 13*19*31*37*43 = 12190001. Bitmaps
   TWO_SUM_Q1[84672] and TWO_SUM_Q2[12190001] (1.5 MB) of residues that are
   x^6+y^6 mod Q. Two lookups per leaf. Expected pass rate ~1%.
Step 1. Let P = all primes p <= 700 (125 primes), each with exponent k_p chosen
   so that p^k_p > Bmax (e.g. 2^19, 3^12, 5^9, 7^7, 11^6, 13^5, ..., p^2 for
   p > sqrt(Bmax)). For each p in P in increasing order:
     r = R mod p^k_p  (incremental: m mod p^k - T_p[b4] - T_p[b3])
     (i)  if r != 0 and r is a unit sixth-power residue mod p^k_p (bitmap
          H_p): for each root rho (at most 6; 2 when p = 2 mod 3; 4 when p=2),
          for each b2 = rho + j*p^k_p <= bound: S = R - b2^6 must be > 0 and
          divisible by p^6 (exact-multiply test as in k7engine.c:div_r7);
          S/p^6 must be a perfect sixth power b1'^6 with p*b1' <= b2. Check
          budgets. Return SOLUTION if so.
     (ii) if R is divisible by p^6 (test exactly): recurse
          is_two_sum(R/p^6, bound/p, budgets adjusted for the factor p).
          (Budgets: dividing both bases by p changes o2 only if p=2, o3 only if
          p=3, o7 only if p=7 -- recompute o's from the new R directly as
          o2 = R' mod 8, o3 = R' mod 9, o7 = R' mod 7, each must be <= 2.)
   NEVER break out of the p loop on a failed hypothesis; only on success.
Step 2. Rough case: both bases are 1 or primes > 700 (composites of primes
   > 700 exceed 709^2 = 502681 > Bmax for FMAX <= 2.1e7 -- assert this).
   ROUGH = {1} u {primes q : 700 < q <= Bmax}. Precompute the sorted array of
   64-bit fingerprints of x^6+y^6 for x >= y in ROUGH (about pi(Bmax)^2/2
   entries; 8e8 at Bmax=4.76e5 -> 6.4 GB as u64, or a 16-bpp blocked Bloom at
   1.6 GB with exact verification restricted to ROUGH x ROUGH). Query it.
   (The same table serves every recursion level since bound/p <= Bmax.)

## Completeness proof (to be included as a comment)
Let (b1,b2) be a true representation. Let p be the least prime <= 700 dividing
b1*b2. If p | exactly one of them, step 1(i) at p finds it (the unit one is b2
or b1; handle both orders: the root search must allow the determined base to
be either the larger or the smaller -- simplest: require only b2 <= bound and
p*b1' <= bound, no ordering). If p | both, step 1(ii) at p recurses on a
strictly smaller instance, complete by induction. If no such p exists, both are
700-rough and <= Bmax, hence in ROUGH, and step 2 finds it.

## Cost model
Leaves per candidate ~ B^2/1300 (measured for caseA2 at F=4M: 7e6 at
B=9.5e4). Step 0 ~10 ns; step 1 runs on ~1% of leaves at ~2 us; total ~30 ns
per leaf. Cumulative leaves to F (from the fitted candidate density
D(F)=6.45e-9 F per f-unit): L(F) ~ 7e-16 F^4. Restricted to k7 in {3,4}
(x0.4): F=8e6 -> 1.1e12 leaves ~ 3.5e4 core-s; F=1e7 -> 2.8e12 ~ 8.4e4 core-s.
Memory: T_p tables 125 x Bmax x 4 B = 240 MB at Bmax 4.76e5; H_p bitmaps and
root lists ~ sum 4*p^k_p bytes ~ 0.5 GB; ROUGH table <= 6.4 GB. No pair Bloom.

## Tests (mandatory, same discipline as SPEC_K7.md)
1. is_two_sum differential: 10^5 random R built from random (b1,b2) <= 3000
   in every divisibility pattern (one divisible by p, both, rough-rough,
   rough-prime) must return TRUE with the correct pair; 10^5 random non-sums
   (R = b1^6+b2^6+1, and random 128-bit values) must return FALSE, checked
   against a brute-force Python two-sum oracle.
2. Planted candidates: build m from four bases in each k7 stratum and confirm
   the full engine reports the decomposition; near-misses must be silent.
3. Candidate counts: 124 / 1314 on the two reference bands.
4. Differential vs caseA2 on 700000..730000: identical CANDIDATE set, identical
   verdicts, and (optionally) identical leaf counts j2_nodes = 4,139,088 when
   the same nb=4-only traversal order is used.
5. Every SOLUTION line is run through strata/verify_solution.py.
