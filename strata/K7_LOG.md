# K7_LOG.md — the `k7 <= 2` stratified class-1 engine

Machine: 4 cores (Intel Xeon @ 2.80 GHz), 15 GB RAM, Linux 6.18.44-fc-v24,
gcc 13.3.0, Python 3.11.15.  Repo `/home/user/ESOP6`, branch
`claude/practical-bohr-ihlohk`, HEAD `36e03d4`.  **Nothing under `src/` was
modified and nothing was committed.**  All runs use `OMP_NUM_THREADS=4`.

Deliverables in this directory: `k7engine.c`, `test_k7.sh`, `lattice_bounds.py`,
this log.  `verify_solution.py` (pre-existing) is the last line of defence;
no `SOLUTION` was ever produced, so it was never invoked on engine output.

---

## 0. Band convention (stated once, used everywhere)

`k7engine FMIN FMAX` processes **every `f` with `FMIN < f <= FMAX`**, i.e. the
half-open interval `(FMIN, FMAX]`.  `src/caseA2.c` uses the **closed** interval
`[fmin, fmax]` (`for(u64 f=fmin; f<=fmax; f++)`, caseA2.c:271).  The two
conventions give identical candidate sets whenever `gcd(FMIN,42) > 1`, because
then `f = FMIN` is skipped by the `2|f, 3|f, 7|f` test anyway.  Every band
endpoint used in this project (2, 700000, 730000, 1000000, 4300000, 4400000,
8000000, 12000000, 20000000) is even, so all counts below are directly
comparable with the repo's.  The half-open convention is the correct one for
chaining bands without double-counting the joint; that is why it was chosen.

---

## 1. Build

```
$ cd /home/user/ESOP6/strata
$ gcc -O3 -march=native -fopenmp -std=gnu11 -Wall -Wextra -Wno-unused-parameter \
      -o k7engine k7engine.c -lm
```
Zero warnings.

---

## 2. What the engine does, and what it does NOT cover

For each candidate `(f,t)` (class 1, `t = zeta*f mod 42^6`, `0<t<f`,
`gcd(f,42)=1`) it forms `m = (f^6-t^6)/42^6` and the exact budgets
`k2 = m mod 8`, `k3 = m mod 9`, `k7 = m mod 7`, then:

| class of candidate | action | status |
|---|---|---|
| `k2 > 4` or `k3 > 4` | discarded | **eliminated** — no decomposition exists |
| `k7 in {5,6}` | discarded | **eliminated** — no decomposition exists |
| `k7 = 0` | needs `7^6 \| m` (valuation lemma), else discarded; otherwise generic 4-sum on `m/7^6`, bases `<= B/7` | **searched** |
| `k7 = 1` | `b4 = r + j*7^6`, `r in ROOTS[m mod 7^6]`; 3-sum on `(m-b4^6)/7^6` | **searched** |
| `k7 = 2` | `b3` enumerated, `b4` from `ROOTS[(m-b3^6) mod 7^6]`; 2-sum on the rest | **searched** |
| `k7 in {3,4}` | counted, **not searched** | **NOT COVERED** |

`k7 in {3,4}` is roughly 4/7 of the surviving candidates and is reported
explicitly as `notcovered_k7_3_4` on every summary line.

### Engineering decisions taken from `BASELINE_LOG.md`

* **Residue masks before any Bloom query.**  Eight moduli
  (64, 27, 49, 13, 19, 31, 37, 43), 2-/3-/4-sum reachability bitsets precomputed,
  applied in most-selective-first order.  Measured effect: at `F = 20 M`,
  470,240,329 leaves produce only 16,544,367 Bloom queries — **3.5 %**, a 28x
  reduction in filter traffic.  This is the single biggest change relative to
  `caseA2`, which queried the filter at every leaf.
* **Non-defective blocked-Bloom hash family.**  8 *independent* 9-bit slices
  taken from bits 0-8/9-17/18-26/27-35 of two independent 64-bit mixes, plus a
  **third** independent mix for the 512-bit line index.  16 bits per pair.
  See section 6 for the measured false-positive rate.
* **`schedule(dynamic,1)` over `f`** (the baseline's `schedule(dynamic,4096)`
  serialised any band narrower than 4096 iterations, FINDING P1).
* **The pair table is over the PRIMED bases** `b' <= B/7`, not `B`.  This is the
  real reason the stratum is cheap: the filter is `49x` smaller than
  `caseA2`'s at the same `fmax` (4.63 GB at `F = 2e7` instead of 227 GB).
* **No 128-bit division in the hot path.**  `R/7^6` is done by multiplying with
  the 2-adic inverse of `7^6 mod 2^128`; the magnitude of the product *is* the
  divisibility test (`x -> x*7^6 mod 2^128` is a bijection that does not
  overflow exactly for `x <= (2^128-1)/7^6`).  Small moduli are reduced with
  `hi*(2^64 mod md) + lo` using literal moduli so gcc emits multiplies.

### Arithmetic and assertions

`f^6` does not fit in 128 bits for `f > 2,642,245`, so `m` is built from the
cyclotomic factorisation `(f-t)(f+t)(f^2-ft+t^2)(f^2+ft+t^2)` with `42^6`
stripped by gcd, exactly as `src/caseA2.c:282-294`.  Asserted, per candidate:

* `rem == 1` after the gcd peel (i.e. `42^6` really divided `f^6-t^6`);
* `m*42^6 + t^6 == f^6` **modulo `2^61-1`** (a full 128-bit check is impossible);
* `7^6` exactly divides `m - b4^6` (k7=1) and `m - b3^6 - b4^6` (k7=2);
* any reported decomposition re-summed in 128-bit equals `m`.

Plus, at startup: 144 roots verified, `ROOTS` mod `7^6` has exactly
`7^6/7 = 16807` non-empty classes each of size 6, `iroot6` spot-checked,
exact-division helper spot-checked, Bloom false-negative spot-check.

---

## 3. Tests — `./test_k7.sh`

```
$ cd /home/user/ESOP6/strata && time ./test_k7.sh
```
Wall 26.7 s.  Every line below was printed by the script; **all PASS**.

```
=== build ===
PASS  build (gcc -O3 -march=native -fopenmp)
=== TEST 1: planted decompositions and near-misses ===
PASS  TEST1 k7=0 planted -> FOUND 14 21 35 63
PASS  TEST1 k7=0 near-miss m+1 -> NONE
PASS  TEST1 k7=0 near-miss m+42^6 (same k2,k3,k7) -> NONE
PASS  TEST1 k7=0 near-miss m-42^6 (same k2,k3,k7) -> NONE
PASS  TEST1 k7=0 planted with B=62 (below max base) -> NONE
PASS  TEST1 k7=1 planted -> FOUND 11 21 35 56
PASS  TEST1 k7=1 near-miss m+1 -> NONE
PASS  TEST1 k7=1 near-miss m+42^6 (same k2,k3,k7) -> NONE
PASS  TEST1 k7=1 near-miss m-42^6 (same k2,k3,k7) -> NONE
PASS  TEST1 k7=1 planted with B=55 (below max base) -> NONE
PASS  TEST1 k7=2 planted -> FOUND 5 19 77 91
PASS  TEST1 k7=2 near-miss m+1 -> NONE
PASS  TEST1 k7=2 near-miss m+42^6 (same k2,k3,k7) -> NONE
PASS  TEST1 k7=2 near-miss m-42^6 (same k2,k3,k7) -> NONE
PASS  TEST1 k7=2 planted with B=90 (below max base) -> NONE
PASS  TEST 1 planted/near-miss (all 15 cases)
PASS  TEST 1b Python oracle agrees the near-misses really have no decomposition
=== TEST 2: candidate-count agreement with src/caseA2.c ===
PASS  TEST 2 (700000,730000]  candidates=124 (expected 124)
PASS  TEST 2 (730000,1000000] candidates=1314 (expected 1314)
PASS  TEST 2b independent Python (700000,730000]=124
PASS  TEST 2b independent Python (730000,1000000]=1314
=== TEST 3: control band 700000..730000 end to end, k7 <= 2 ===
PASS  TEST 3 control band: 0 solutions, exit code 0
=== TEST 4: 2000-random-m differential per stratum vs Python brute force ===
       k7=0 built       decomposable=True  : 2000
       k7=0 perturbed   decomposable=False : 2000
       k7=1 built       decomposable=True  : 2000
       k7=1 perturbed   decomposable=False : 2000
       k7=2 built       decomposable=True  : 2000
       k7=2 perturbed   decomposable=False : 2000
       positives=6000 negatives=6000 total=12000
PASS  TEST 4 differential: 12000/12000 agree
=== TEST 5: isqrt6 unit test vs Python on 1e5 values ===
PASS  TEST 5 isqrt6: 100008/100008 exact
=== TEST 6: first-50 / last-50 candidate m-reconstruction in Python ===
PASS  TEST 6 m-reconstruction: 200/200 candidates exact
=== TEST 7: Bloom filter health (no false negatives, FP rate vs theory) ===
BLOOMSTAT Bp=20000 pairs=200010000 bpp=16 false_negatives=0 queries=20000000
          false_positives=17859 fp_rate=0.00089295 ideal=0.000574496 ratio=1.5543
       false negatives = 0
       measured FP rate = 0.00089295   blocked-Bloom theory = 0.000872812   ratio = 1.023x
PASS  TEST 7 Bloom: zero false negatives and FP rate within 1.10x of blocked-Bloom theory
```

Notes on the tests, and on two tests that were *strengthened* rather than
weakened when they first failed:

* **TEST 1 near-misses.**  `m+1` is a weak negative: it moves the candidate into
  a different `(k2,k3,k7)` stratum and is often rejected on budgets alone.  So
  `m +/- 42^6` was added: `42^6 ≡ 0 mod 8, 9, 7`, so those near-misses have the
  *identical* budgets and drive the full search path to exhaustion.  A fourth
  case plants the true `m` with `B = max(b_i) - 1` and requires silence.
* **TEST 1 initially failed** on the `B` too small case for `k7=1` because the
  chosen `B=60` was in fact above every base (max base 56) — the *test* was
  wrong, and was corrected to `B = max(b_i)-1` rather than removed.
* **TEST 4** compares against a brute-force Python 4-sum finder that knows
  nothing about the strata (it enumerates `b4 >= b3` and hashes 2-sums).  This
  is a real test of completeness: any class-1 4-sum decomposition of `m` has
  *exactly* `m mod 7` bases coprime to 7, so the stratum the engine picks is
  forced and it must find a decomposition iff one exists.  6000 constructed
  (decomposable) and 6000 perturbed (all verified non-decomposable by the
  oracle) cases, 12000/12000 agreement.  Every engine answer is additionally
  re-summed exactly.
* **TEST 7 reference was corrected, not relaxed.**  `BASELINE_LOG.md` compares
  the measured FP rate to the *non-blocked* ideal `(1-e^{-1/2})^8 = 5.745e-4`
  and concludes the hash family is 1.64-1.67x defective.  That comparison is
  against the wrong baseline: a **blocked** filter with 512-bit lines at 16
  bits/pair has an intrinsic penalty, because the number of items landing in a
  line is Poisson(32) and the FP rate is convex in that number.  The correct
  reference is `E_{n~Poisson(32)}[(1-(1-1/512)^{8n})^8] = 8.7281e-4`, i.e.
  **1.52x the non-blocked ideal all by itself**.  Measured directly, in the same
  harness, on 200,010,000 inserted pairs and 20,000,000 non-member queries:

  | hash family | FP rate | vs non-blocked ideal | vs blocked theory |
  |---|---|---|---|
  | old (`(h2 >> 8i) & 511`, one mix, line from the same mix) | 9.634e-4 | **1.677x** | 1.104x |
  | new (8 x 9-bit slices, two mixes, third mix for the line) | 8.930e-4 | 1.554x | **1.023x** |

  So the repaired family is within 2.3 % of what a blocked filter can do, and
  the *honest* size of the old defect is **7.4 %**, not 67 %.  The baseline log's
  "1.67x" figure is reproduced exactly here, but most of it is the blocking
  penalty, not the hash.  I would correct that claim in `BASELINE_LOG.md`.

---

## 4. Runs

### 4.1 Control band (700000, 730000]

```
$ OMP_NUM_THREADS=4 time ./k7engine 700000 730000 --candidates cand_700k_730k.txt
bloom: Bp=2482 pairs=3081403 lines=96294 (0.006 GB, 16 bits/pair)
bloom built
BAND fmin=700000 fmax=730000 convention=(fmin,fmax] candidates=124
 k7_0=19 k7_1=18 k7_2=13 k7_3=18 k7_4=16 k7_5=22 k7_6=18
 elim_k2k3=81 elim_k7=15 notcovered_k7_3_4=9
 processed_k7_0=7 processed_k7_1=9 processed_k7_2=3 val_elim_k7_0=7
 leaves=1823 bloom_queries=70 bloom_positives=0 exact_verifications=0 dfs_nodes=0
 solutions=0 bloom_build_s=0.167 search_s=0.005 elapsed_s=0.172 threads=4

real 0m0.193s   user 0m0.668s   exit 0
```

| k7 | candidates | of which eliminated on k2/k3 | searched |
|---|---|---|---|
| 0 | 19 | 12 | 7 (all 7 then eliminated by the `7^6 \| m` valuation lemma) |
| 1 | 18 | 9 | 9 |
| 2 | 13 | 10 | 3 |
| 3 | 18 | — | **0 — NOT COVERED (9 survive budgets)** |
| 4 | 16 | — | (included in the 9 above) |
| 5 | 22 | — | eliminated (`k7 > 4`) |
| 6 | 18 | — | eliminated (`k7 > 4`) |

**0 solutions. 0.193 s wall.**  (`caseA2` needs 8.0 s for the same band, but it
covers all `k7`.)

### 4.2 (2, 730000]

```
$ OMP_NUM_THREADS=4 time ./k7engine 2 730000 --candidates cand_2_730k.txt
BAND fmin=2 fmax=730000 convention=(fmin,fmax] candidates=1492
 k7_0=194 k7_1=221 k7_2=206 k7_3=216 k7_4=218 k7_5=221 k7_6=216
 elim_k2k3=989 elim_k7=139 notcovered_k7_3_4=150
 processed_k7_0=69 processed_k7_1=78 processed_k7_2=67 val_elim_k7_0=69
 leaves=32572 bloom_queries=1492 bloom_positives=1 exact_verifications=1 dfs_nodes=1929
 solutions=0 bloom_build_s=0.173 search_s=0.109 elapsed_s=0.283 threads=4

real 0m0.303s   exit 0
```
1492 candidates — exactly the count `caseA2` reports for `[2,730000]`
(BASELINE_LOG TASK 3).  **0 solutions. 0.303 s wall.**  All 69 `k7=0`
candidates fail `7^6 | m` (see section 7 — the first `f` where that is even
possible is 11,936,849).  150 candidates with `k7 in {3,4}` are **not covered**.

### 4.3 Timing band (4300000, 4400000]

```
$ OMP_NUM_THREADS=4 time ./k7engine 4300000 4400000
bloom: Bp=14965 pairs=111983095 lines=3499472 (0.224 GB, 16 bits/pair)
BAND ... candidates=2795 k7_0=348 k7_1=404 k7_2=413 k7_3=407 k7_4=411 k7_5=406 k7_6=406
 elim_k2k3=1868 elim_k7=278 notcovered_k7_3_4=269
 processed_k7_0=110 processed_k7_1=138 processed_k7_2=132 val_elim_k7_0=110
 leaves=4898435 bloom_queries=159019 bloom_positives=146 exact_verifications=146
 dfs_nodes=191673 solutions=0 bloom_build_s=6.258 search_s=0.292 elapsed_s=6.550

real 0m6.591s   user 0m25.614s   exit 0
```

### 4.4 Three further anchors (run to make the extrapolation a measurement)

```
(8000000, 8100000]   cands=5057  k7_0=637 k7_1=733 k7_2=737 k7_3=743 k7_4=740 k7_5=732 k7_6=735
                     leaves=26698206   bqueries=863953   bpos=771
                     build=24.910 s  search=1.643 s   real 0m26.691s   solutions=0
(12000000,12100000]  cands=7493  k7_0=963 k7_1=1084 k7_2=1102 k7_3=1096 k7_4=1093 k7_5=1081 k7_6=1074
                     leaves=100537451  bqueries=3560024  bpos=3207
                     build=66.228 s  search=6.041 s   real 1m12.509s   solutions=0
(20000000,20100000]  cands=12530 k7_0=1779 k7_1=1789 k7_2=1799 k7_3=1784 k7_4=1780 k7_5=1806 k7_6=1793
                     leaves=470240329  bqueries=16544367 bpos=14668
                     build=240.016 s search=38.055 s  real 4m38.843s   solutions=0
```
All `bloom_positives` were false positives (`solutions=0` everywhere), at rates
9.18e-4 / 8.94e-4 / 9.01e-4 / 8.87e-4 — flat, and equal to the blocked-Bloom
theory `8.73e-4` to within 3 %.

---

## 5. Cost model and time to 1e7 / 1.5e7 / 2e7

### Derivation (from SPEC_K7.md's cost analysis)

With `M = 42^6`, `B = (f-1)/42`, `R7 = 7^6`:

1. **Candidate density.**  Each of the 143 non-trivial roots gives
   `0 < zeta*f mod M < f` with probability `~f/M`, and only `f` coprime to 42
   count (density 12/42), so `D(F) = dK/df = Theta(F/M)`.
   Measured: `D(F)/F = 6.43, 6.28, 6.22, 6.25 e-9` at `F = 4.35, 8.05, 12.05,
   20.05 M` — **`D(F) = 6.27e-9 * F`**, flat to 3 %.
2. **Fraction reaching a search.**  `P(k7<=2) = 3/7`, `P(k2<=4) = 5/8`,
   `P(k3<=4) = 5/9`; product **14.9 %**.  Measured 13.6 / 14.0 / 14.5 / 15.3 %.
3. **Per-candidate leaf cost.**  The spec's `candidates ~ F^2, per-candidate
   cost ~ B` is *not* what this stratum costs.  For `k7=2` the pairs `(b3,b4)`
   with `b3,b4` coprime to 7 and `b3^6+b4^6 ≡ m (mod 7^6)` number
   `~ 0.857 * B^2 / (2*R7)`; for `k7=1` the `~6B/R7` admissible `b4` each drive a
   3-sum of length `B' = B/7`, again `~0.857*B^2/R7`.  So the per-candidate leaf
   count is `Theta(B^2/7^6) = Theta(F^2/(42^2 * 7^6))`, not `Theta(B)`.
   (`B^2/7^6` exceeds `B` precisely for `B > 7^6`, i.e. `F > 4.94e6`.)
   Measured `leaves/candidate / (B^2/7^6) = 0.199, 0.174, 0.192, 0.186` — the
   constant `~0.19` is the `k2/k3` stride and mask saving.
4. Hence **leaves per unit `f` is `Theta(F^3)`** and **cumulative leaves to `F`
   is `Theta(F^4)`**:
   ```
   l(F) = D(F) * 0.099 * 0.19 * (F/42)^2/7^6  =  5.7e-19 * F^3    leaves per unit f
   L(F) = 1.42e-19 * F^4                                          leaves up to F
   ```
   Least-squares fit to the four anchors: `l(F) = 5.385e-19 * F^3.0030` —
   the measured exponent is **3.003**, i.e. the derivation is exact.

   For comparison, the digest's model for `caseA2` is `2.4e-16 * F^4` j=2 nodes
   cumulative; this engine is **~1700x fewer leaves**, while covering the ~3/7
   of `k7` values that is `{0,1,2}` (about 10 % of class-1 solution mass by the
   spec's heuristic).

5. **Wall time** is slightly steeper than `F^3` per unit `f` because a Bloom
   query is a DRAM latency event in a filter that grows like `F^2`.  Fitted on
   the four anchors (4 threads, `schedule(dynamic,1)`):
   ```
   w(F) = 2.03e-27 * F^3.1801     search wall-seconds per unit of f
   T_search(F) = integral_0^F w  = 2.03e-27 * F^4.1801 / 4.1801
   ```
   Fit residuals: 1.11, 0.89, 0.90, 1.13 (measured/fit) — +/-13 %.
6. **Bloom build** is a one-off per run (one run covers `2..F`):
   `pairs(F) = Bp(Bp+1)/2`, `Bp = (F-1)/42/7`; measured insert rate
   `1.79e7, 1.52e7, 1.28e7, 9.74e6 pairs/s` at the four heights (it falls as the
   filter leaves cache), fitted `r(F) = 8.07e9 * F^-0.397 pairs/s`.

### Estimate (single run `2 .. F`, 4 cores, 16 bits/pair)

| target `F` | `Bp` | pairs | Bloom | build | search (fit `F^4.18`) | search (pure `F^4`) | **TOTAL wall** |
|---|---|---|---|---|---|---|---|
| 4.3e6 | 14,625 | 1.07e8 | 0.21 GB | 6 s | 3 s | 4 s | **8 s** |
| **1.0e7** | 34,013 | 5.79e8 | 1.16 GB | 43 s | 89 s | 118 s | **~2.2 min** |
| **1.5e7** | 51,020 | 1.30e9 | 2.60 GB | 115 s | 482 s | 598 s | **~10 min** |
| **2.0e7** | 68,027 | 2.31e9 | 4.63 GB | 228 s | 1605 s | 1889 s | **~31 min** |

Uncertainty: the two search columns bracket the model choice; the fit residuals
add +/-13 %; so **"under an hour to `f = 2e7` on this box"** is the honest claim,
and RAM (4.63 GB of 15 GB) is not the constraint — a single run could go to
`F ≈ 3.5e7` (14 GB) before it is.

**I did not launch the production run**, per instruction.  Note that the three
anchor bands at 8 M, 12 M and 20 M *are* new territory and were searched:
each reported `solutions=0`.

Caveat that must accompany any use of this: the engine covers `k7 in {0,1,2}`
only.  `k7 in {3,4}` — about 4/7 of the budget-surviving class-1 candidates —
is **not covered at any height**, and classes 2-5 are outside class 1 entirely.

---

## 6. Coverage statements (the only claims supported by these runs)

```
Class 1, k7 in {0,1,2}, 700000 < f <= 730000 : 124 candidates
  (19/18/13 with k7=0/1/2; 19 searched, 9 with k7 in {3,4} NOT covered),
  0 solutions, 0.193 s wall, 4-core Xeon @2.80GHz.

Class 1, k7 in {0,1,2}, 2 < f <= 730000 : 1492 candidates
  (194/221/206 with k7=0/1/2; 214 searched, 150 with k7 in {3,4} NOT covered),
  0 solutions, 0.303 s wall, 4-core Xeon @2.80GHz.

Class 1, k7 in {0,1,2}, 4300000 < f <= 4400000  : 2795 candidates, 0 solutions, 6.6 s.
Class 1, k7 in {0,1,2}, 8000000 < f <= 8100000  : 5057 candidates, 0 solutions, 26.7 s.
Class 1, k7 in {0,1,2}, 12000000 < f <= 12100000: 7493 candidates, 0 solutions, 72.5 s.
Class 1, k7 in {0,1,2}, 20000000 < f <= 20100000: 12530 candidates, 0 solutions, 278.8 s.
```
Nothing more than that is claimed.

---

## 7. Deep-stratum lower bounds — `lattice_bounds.py`

If, in addition to the class-1 structure, **all four** small bases `b_i` are
divisible by an extra prime `p in {2,3,7}` (respectively by 42), then
`f^6 ≡ t^6 (mod 42^6 * p^6)` (resp. `mod 42^12`), so `t ≡ zeta*f (mod N)` for one
of the 144 sixth roots of unity `zeta` mod `N`.  The admissible pairs form the
lattice `L = {(f,t) : t ≡ zeta f (mod N)}` with basis `(1,zeta), (0,N)` and
`det L = N`.  The least `f` for which `L` meets the cone `0 < t < f` (with
`gcd(f,42)=gcd(t,42)=1`) is a **provable lower bound** on `f` in that stratum.

Method: Gauss(-Lagrange)-reduce the basis, then enumerate the lattice points of
the *triangle* `{0 < t < f <= A}` exactly in integer arithmetic, doubling `A`
until it is non-empty.  Two points that matter and that the slow draft got
wrong:

* enumerating a **ball** rather than the triangle is catastrophic for `zeta`
  near `-1`, where the lattice has a very short vector along `t = -f` (outside
  the cone) and the ball holds `O(N)` points while the cone holds almost none;
* within one `y`-layer `f` is linear in `x` and the predicate
  `gcd(f,42)=gcd(t,42)=1` depends only on `x mod 42`, so at most 42 `x`-values
  from the `f`-minimising end need inspection.  That bounds the work even when
  the triangle contains astronomically many lattice points.

`gcd(t,42)=1` is automatic: `42 | N` and `zeta` is a unit mod `N`, hence mod 42,
so `gcd(t,42) = gcd(f,42)`.

The script self-tests the enumeration against direct brute force on 231 random
small lattices before doing anything else.

```
$ cd /home/user/ESOP6/strata && time python3 lattice_bounds.py
selftest: triangle enumeration and least-f agree with brute force on 231 random lattices
==============================================================================
all four b_i divisible by 2    N = 42^6 * 2^6 = 351298031616
  LOWER BOUND  min over zeta != 1 of least f :  235283
      attained at zeta = 171295920127   with t = 63725   (t/f = 0.270844047)
  next four    : 255727, 339239, 366481, 406841
  median / max : 1426529 / 175649015809
==============================================================================
all four b_i divisible by 3    N = 42^6 * 3^6 = 4001504141376
  LOWER BOUND  min over zeta != 1 of least f :  395243
      attained at zeta = 1264528075807   with t = 40949   (t/f = 0.103604618)
  next four    : 749537, 877307, 1149703, 1359371
  median / max : 6154429 / 2000752070689
==============================================================================
all four b_i divisible by 7    N = 42^6 * 7^6 = 645779095649856
  LOWER BOUND  min over zeta != 1 of least f :  11936849
      attained at zeta = 396825048481409   with t = 8441617   (t/f = 0.707189728)
  next four    : 13629013, 13739587, 13943411, 15104221
  median / max : 65694661 / 322889547824929
==============================================================================
all four b_i divisible by 42   N = 42^12 = 30129469486639681536
  LOWER BOUND  min over zeta != 1 of least f :  1698000953
      attained at zeta = 14470590145148983297   with t = 1335203897   (t/f = 0.786338721)
  next four    : 2428182959, 2548324099, 3121355863, 3138722329
  median / max : 14022685729 / 15064734743319840769
==============================================================================
SUMMARY (provable lower bounds on f for class-1 solutions in the deep strata)
  all four b_i divisible by 2    N=351298031616            f >= 235283
  all four b_i divisible by 3    N=4001504141376           f >= 395243
  all four b_i divisible by 7    N=645779095649856         f >= 11936849
  all four b_i divisible by 42   N=30129469486639681536    f >= 1698000953

real 0m0.114s
```

### Independent double-check of the `42^6 * 7^6` bound (different route)

`lattice_bounds.py --verify7` throws the lattice away and scans every
`f = 1 .. 11936849` against all 144 roots directly (`t = zeta*f mod N`,
numpy int64 in chunks small enough that `zeta*k < 2^63`).  Since
`N = 6.458e14` dwarfs the bound, the residue in `[0,N)` is the only possible `t`.

```
verify_7_by_scan: scanning f = 1 .. 11936849 against all 144 roots
  scan found least f = 11936849 (t = 8441617, zeta = 396825048481409); 1 hits at or below the bound
PASS  independent scan reproduces the 42^6*7^6 bound f >= 11936849
real 0m27.308s   (0.114 s without --verify7)
```
Exactly one hit, at exactly the claimed minimum.

### Third, engine-side cross-check of the same number

`7^6 | m` is *equivalent* to `f^6 ≡ t^6 (mod 42^6*7^6)`, which is exactly the
`k7=0` valuation test the engine applies.  So the bound predicts that **no
`k7=0` candidate below `f = 11,936,849` can survive the valuation lemma** — and
indeed `val_elim_k7_0 == processed_k7_0` in every run reported above
(7/7, 69/69, 110/110, 218/218, 335/335, 621/621).  Running the single-`f` band:

```
$ ./k7engine 11936848 11936849
BAND ... candidates=1 k7_0=1 k7_1=0 ... elim_k2k3=1 ... processed_k7_0=0 ... solutions=0
$ python3 -c "f,t=11936849,8441617; m=(f**6-t**6)//42**6; print(m%7**6==0, m%8, m%9, m%7)"
True 6 0 0
```
The engine finds exactly the predicted candidate `(f,t) = (11936849, 8441617)`,
`7^6 | m` confirmed in Python big integers — and it is then killed on the spot
by `k2 = m mod 8 = 6 > 4`.  Three independent routes agree.

**What these bounds mean.**  They are lower bounds on `f` for a class-1 solution
whose four small bases share the extra prime.  The `p=7` bound `f >= 11,936,849`
is the interesting one: it says the entire `k7 = 0` stratum is empty below
11.9 M *for structural reasons*, independent of any search.  The `42^12` bound
`f >= 1.698e9` rules out `42 | b_i` for all four, far beyond any feasible search.
They say nothing about `k7 in {1,2,3,4}`, which is where the mass is.

---

## 8. Things I am not sure about / would flag

1. **`k7 in {3,4}` is genuinely not covered** and is the majority of the
   budget-surviving candidates.  The stratification buys ~170x per unit of
   coverage, but it buys ~10 % coverage.
2. The `k7=0` generic 4-sum branch is essentially **untested in production**:
   `7^6 | m` never once held below `f = 2.01e7` outside the single candidate at
   11,936,849 (which died on budgets).  It is exercised only by TEST 1 and by
   2000 synthetic cases in TEST 4.  I am confident it is correct, but it has
   never run on a real `m`.
3. **`BASELINE_LOG.md`'s "defective hash family, 1.67x" should be corrected.**
   Measured head-to-head in one harness, the old family is 1.104x and the new
   one 1.023x of the *blocked*-Bloom expectation; 1.52x of the 1.67x is the
   intrinsic 512-bit blocking penalty and cannot be removed by any hash.
   Repairing the hash is worth ~7 %, not ~67 %.
4. The wall-time extrapolation exponent (`F^3.18` per unit `f`) is empirical and
   memory-system-dependent; the *leaf-count* exponent (`F^3.003`) is the solid
   one.  I quoted both columns rather than picking one.
5. Degenerate decompositions (fewer than four positive bases) are **not**
   reported; `caseA2` also tries `nb=3,2,1` and prints `DEGENERATE`.  The `k7`
   strata are defined by a count over four bases, so the engine only searches
   for exactly four positive bases.  A 4-term solution would be missed.
6. The single-`f` run at 11,936,848..11,936,849 spends 67 s building a 1.65 GB
   Bloom filter it never queries.  Harmless, but the filter should be built
   lazily on the first query if narrow bands ever matter.
