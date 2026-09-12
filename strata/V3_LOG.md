# V3_LOG.md — the table-free-leaf engine (`strata/v3engine.c`)

Machine: 4 cores (Intel Xeon @ 2.80 GHz), 15 GB RAM, Linux 6.18.44-fc-v24,
gcc 13.3.0, Python 3.11.15.  Repo `/home/user/ESOP6`, branch
`claude/practical-bohr-ihlohk`.  **Everything written by this session is under
`strata/`; nothing under `src/` or `bin/` was modified and nothing was
committed.**

New files: `v3engine.c`, `test_v3.sh`, `v3_testlib.py`, `measure_v3.sh`,
`caseA2_nb4.c` (a two-line edit of the pre-existing `strata/caseA2_timed.c`),
this log.

The `k7engine` production run that owned the machine at the start of this
session finished on its own while the tests were being written
(`strata/runs/k7_2_20M.out`, `EXIT=0`, 1,248,387 candidates, 21,826,977,286
leaves, **0 solutions**, class 1 / `k7 <= 2` / `2 < f <= 2e7`, 1925 s).  Nothing
under `strata/runs/` was touched.  All timings below therefore ran on the idle
machine with `OMP_NUM_THREADS=4`; the differential tests were also run at 2
threads with identical results.

**No `SOLUTION` was produced by anything in this session.**

---

## 0. Conventions

Band `(FMIN, FMAX]`, half-open, exactly as `k7engine` (see K7_LOG.md §0); all
endpoints used here are even, so the counts are directly comparable with
`src/caseA2.c`'s closed-interval `[fmin,fmax]`.

`v3engine` searches for **exactly four positive bases** (`nb = 4`).  `caseA2`
additionally runs `nb = 3, 2, 1` and prints `DEGENERATE`; `v3engine` does not,
exactly as `k7engine` does not (K7_LOG.md §8.5).  This is the only semantic
difference between the two and it is visible in the leaf counts (§5).

Default stratum: `--k7 3,4` (what `k7engine` does not cover).  `--k7 0,1,2,3,4`
processes every candidate that survives the budgets.

---

## 1. Build

```
$ cd /home/user/ESOP6/strata
$ gcc -O3 -march=native -fopenmp -std=gnu11 -Wall -Wextra -Wno-unused-parameter \
      -o v3engine v3engine.c -lm
```
Zero warnings.  `test_v3.sh` fails the build check if any warning appears.

---

## 2. Two errors in SPEC_V3.md, and what was done about them

### 2.1 `Q2` is not the product it is said to be

SPEC_V3.md step 0 asks for `Q2 = 13*19*31*37*43 = 12190001`.  That product is
**12,182,287**.  Building the bitmap for 12,190,001 (which factors as
7 · 1741429) produced 1,283,160 distinct sixth powers instead of the expected
`3*4*6*7*8 = 4032` and a useless mask.  The engine uses 12,182,287.

### 2.2 The completeness proof of step 1 is false for small primes

SPEC_V3.md step 1(i) determines the base `b2` not divisible by `p` from
`r = R mod p^k_p` with `p^k_p > Bmax`.  That is valid only when `p^k_p` divides
`b1^6`, i.e. only when `6 * v_p(b1) >= k_p`, i.e. only when

        p^(6 * v_p(b1))  >  Bmax .

For `p >= 11` (and `Bmax < 11^6 = 1.77e6`, i.e. `FMAX < 7.4e7`) `k_p <= 6` and
`v_p >= 1` already suffices, so every prime factor in `[11, 700]` really does
determine.  For `p in {2,3,5,7}` it does not: if `b1 = 2*q` and `b2 = q'` then
`R mod 2^19` says nothing about `b2^6 mod 2^19`, because `R - b2^6 = 64*q^6` has
2-adic valuation 6, not 19.  The spec's claim *"if `p` divides exactly one of
them, step 1(i) at `p` finds it"* is therefore false, and its residual class
`ROUGH = {1} u {primes > 700}` is too small: the pair `(2*701, 709)` is missed by
both step 1 and step 2 as specified.  (`test_v3.sh` TEST 1c contains this exact
pair; with the residual table disabled the engine answers `NO` for it, which is
the spec's algorithm, and `YES` with the table, which is the truth.)

**The fix** (implemented, and proved as a comment at the head of `v3engine.c`).
Let `lim = floor(Bmax^(1/6))` — 6 at `Bmax = 104761`, 7 at 238095, 8 at 476190 —
and

    E = { e : every prime power p^v || e satisfies p^v <= lim }   (12 / 24 / 32 elements)
    V = { e*q <= Bmax : e in E, q = 1 or q prime > PCUT }

with `PCUT` the least prime whose square exceeds `Bmax` (so that no `b <= Bmax`
has two prime factors above `PCUT`).  Step 1 determination is applied for every
prime `p <= PCUT`; step 1 recursion for every such `p` with `p^6 | R`; step 2
queries a pair filter over `V x V` instead of `ROUGH x ROUGH`.  The induction is:

* `p <= PCUT` divides both bases  ->  `p^6 | R`, recurse on `(R/p^6, bound/p)`;
* a base has a prime power `p^v` with `p <= PCUT` and `p^(6v) > Bmax`  ->  step 1
  determines the other base exactly (it is `< p^k_p`, hence the unique
  representative of its residue class);
* otherwise each base is `e*q` with `e in E` and `q = 1` or one prime `> PCUT`
  (two such primes would exceed `PCUT^2 > Bmax`), i.e. both lie in `V`.

A common prime factor `q > PCUT` is covered by the third case: then the bases are
`q*x`, `q*y` with `x, y <= Bmax/PCUT < PCUT` and `gcd(x,y) = 1`, so `x, y in E`
and both bases are in `V`.

`PCUT` also replaces the spec's fixed 700: `PCUT^2 > Bmax` is exactly what the
"at most one rough prime factor" argument needs, and it is 137 at `FMAX = 730000`,
331 at 4.4e6, 419 at 1e7 — i.e. 33 / 67 / 81 primes in the loop instead of 125,
for free.

**Cost of the fix.**  `|V| ~ 2.8 * pi(Bmax)` instead of `pi(Bmax)`, so ~7.8x the
pairs of the spec's ROUGH table.  A sorted array of 64-bit fingerprints (the
spec's preference) would be 3.7 GB at `FMAX = 4.4e6` and 21 GB at `FMAX = 1e7`,
so **RAM forces a blocked Bloom filter**: 8 bits per pair, the repaired
`k7engine` hash family (three independent mixes, 8 independent 9-bit slices,
512-bit lines).  Measured false-positive rate 2.9 % (blocked-Bloom theory at
8 bits/pair is 2.5 %); **every** positive is settled exactly, and the exact
search is restricted to `V x V` by binary-searching the sorted `V` and testing
membership of the complementary base in a bitmap over `[0, Bmax]`.  A false
positive costs ~1600 scanned `V` elements, each rejected by one sixth-power
residue lookup mod `Q1`; measured total cost of the residual table is §4.

---

## 3. What the engine does

Outer structure: `src/caseA2.c:dfs()` replicated line for line — same nested
windows over the two largest bases, same `d2*d3*d7` stride, same skip rules, same
five residue masks (64, 27, 49, 13, 43) with the same incremental residues, same
traversal order — with `nb = 4` only.  Leaves are counted at exactly the point
`caseA2_timed.c` counts `j2_nodes`.

The `j = 2` leaf is the only thing that changes:

1. two (now three) composite 2-sum residue masks, applied from residues carried
   incrementally down the DFS (`PQ1[b]`, `PQ2[b]`, `PQ3[b]` tables; no division);
2. the prime loop (determination + recursion), never broken out of on a failed
   hypothesis;
3. the residual `V x V` pair filter.

Determination roots: for `p >= 5` a table of the sixth roots mod `p` (total
`sum p` entries) plus a Hensel lift to `p^k_p` (`(6 x^5)^-1 mod p` is tabulated,
so the lift is `k_p - 1` steps of four Barrett-reduced 32-bit multiplies).  For
`p = 2` and `p = 3` Hensel does not apply (`p | 6`), so one root per residue is
tabulated directly mod `2^k_2` and `3^k_3` and the remaining roots come from the
4 (resp. 6) sixth roots of unity.

Masks: SPEC_V3's `Q1 = 64*27*49` and `Q2 = 13*19*31*37*43`, plus a third,
`Q3 = 61*67*73*79`.  Measured residue densities of `{x^6+y^6 mod q}`:

| q | 2^6 | 3^3 | 3^6 | 7^2 | 13 | 19 | 31 | 37 | 43 | 61 | 67 | 73 | 79 | 97,103,127,151 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| density | .2656 | .2593 | .2236 | .3061 | .3846 | .5263 | .5161 | .5135 | .6744 | .8361 | .8358 | .8356 | .8354 | 1.0000 |

so `3^6` replaces `3^3` in `Q1` for free and `Q3` is the last useful modulus that
fits in 32 bits.  Combined residue density 0.00767 → 0.00032.  Measured effect on
the control band: leaves reaching the prime loop fell from 456,238 to 196,806
(2.3x), and on `(4.3e6, 4.4e6]` from 16.5 % of leaves to 7.0 %.

Everything else (144 roots, gcd-peeled `m` with the `mod 2^61-1` check, `iroot6`,
exact division by `p^6` via the 2-adic inverse, blocked Bloom with the repaired
hash, candidate file format, `SOLUTION` / `SOLUTION-SIX` protocol) is taken from
`k7engine.c` unchanged in substance.

### Memory (reported by the engine on every run)

| `FMAX` | `Bmax` | `PCUT` | primes | `|E|` | `|V|` | pairs | `T_p` | pair filter | total |
|---|---|---|---|---|---|---|---|---|---|
| 730,000 | 17,380 | 137 | 33 | 12 | 5,960 | 1.8e7 | 2.3 MB | 17.8 MB | 0.02 GB |
| 4,400,000 | 104,761 | 331 | 67 | 12 | 30,255 | 4.6e8 | 28.1 MB | 457.7 MB | 0.49 GB |
| 10,003,000 | 238,166 | 491 | 93 | 24 | (§4) | (§4) | | | |

For comparison, `caseA2`'s pair Bloom over all bases is 227 GB at `FMAX = 2e7`
and `k7engine`'s (over primed bases) is 4.63 GB; `v3engine`'s residual filter is
0.46 GB at 4.4e6.

---

## 4. Runs and measurements

(filled in below)

---

## 5. Tests — `./test_v3.sh`

(filled in below)

---

## 6. Coverage statements

(filled in below)
