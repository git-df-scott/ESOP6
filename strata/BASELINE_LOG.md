# ESOP6 — strata baseline log

Machine: 4 cores, 15 GB RAM, Linux 6.18.44, gcc 13.3.0 (Ubuntu 24.04).
Repo `/home/user/ESOP6`, branch `claude/practical-bohr-ihlohk`, HEAD `36e03d4`, working tree clean at start.
`perf` is **not installed** on this box (`which perf` → nothing), so TASK 2 uses `clock_gettime`
instrumentation in a copy of the engine rather than sampling profiles.
All runs use `OMP_NUM_THREADS=4` unless noted. Nothing in `src/` was modified.

---

## TASK 1 — build and validate

### Build

```
$ cd /home/user/ESOP6 && time make
mkdir -p bin
cc -O3 -march=native -fopenmp -Wall -o bin/search src/search.c -lm
cc -O3 -march=native -fopenmp -Wall -o bin/caseA2 src/caseA2.c -lm
cc -O3 -march=native -fopenmp -Wall -o bin/caseA3 src/caseA3.c -lm
cc -O3 -march=native -fopenmp -Wall -o bin/caseA src/caseA.c -lm
cc -O2 -o bin/audit src/audit.c -lm
cc -O3 -march=native -fopenmp -Wall -o bin/frontier_audit src/frontier_audit.c -lm

real    0m3.884s
user    0m2.054s
sys     0m0.617s
```

Clean build, zero warnings.

### `make validate`

```
$ OMP_NUM_THREADS=4 time make validate
== fifth-power validation (expect f=144 solution) ==
SOLUTION k=5 f=144 : 133 110 84 27
done: k=5 f in [2,150], solutions found: 1
== candidate enumeration audit (expect all OK) ==
nroots=144
f=730001 brute=0 fast=0 OK
f=800003 brute=1 fast=1 OK
f=899999 brute=0 fast=0 OK
f=999983 brute=0 fast=0 OK
f=1299983 brute=0 fast=0 OK
f=54321011 brute=2 fast=2 OK

real    0m0.727s
```
PASS. Lander–Parkin `27^5+84^5+110^5+133^5 = 144^5` rediscovered.

### `make equiv`

```
$ OMP_NUM_THREADS=4 time make equiv
j2_nodes=4139088 skipped=0        evaluated=4139088     (NB=1)
done: f in [700000,730000] candidates=124 processed=124 found=0
j2_nodes=33112704 skipped=28973616 evaluated=4139088    (NB=8)
done: f in [700000,730000] candidates=124 processed=124 found=0

real    0m32.782s
user    2m2.907s
```
PASS — evaluated counts identical (4,139,088), matching the canonical constant.
Note `j2_nodes` goes 4,139,088 → 33,112,704 = exactly 8×, i.e. the upper DFS is
genuinely re-walked once per bucket. This is the direct confirmation of digest **B1**.

### `make prototypes`

```
routed_join N=2000 buckets=8 pairs=2001000 pair_passes=1 queries=20000 query_passes=1
            peak_bucket_entries=250862 disk_pair_bytes=32016000 hits=10001 mismatches=0
mmap_bloom  N=5000 pairs=12502500 pair_passes=1 queries=20000 query_passes=1
            file_bytes=18753792 filter_positives=10046 false_positives=45 false_negatives=0

real    0m3.023s
```
PASS. The 45/≈9,999 = **0.45 %** false-positive figure at 12 bpp reproduces exactly
(digest **B4**: 1.43× the idealized 0.314 %).

### `make frontier-audit`

```
[700000,730000]  candidates=124   expected=124   PASS
[730000,1000000] candidates=1314  expected=1314  PASS
[1000000,1500000] candidates=3523 expected=3523  PASS
[1500000,2000000] candidates=5129 expected=5129  PASS
[2000000,2500000] candidates=7222 expected=7222  PASS
[2500000,3200000] candidates=12443 expected=12443 PASS
[3200000,4000000] candidates=18003 expected=18003 PASS
[4000000,4300000] candidates=8050  expected=8050  PASS
production_total=55684 expected=55684 PASS

real    0m0.430s
```

### Control band, production binary

```
$ OMP_NUM_THREADS=4 time ./bin/caseA2 700000 730000
building bloom: Bmax=17381 pairs=151058271 lines=4720571 (0.3 GB)
bloom built
done: f in [700000,730000] candidates=124 processed=124 found=0

real    0m7.998s
user    0m28.639s
sys     0m1.917s
```

**Result: 124 candidates, found=0, zero `SOLUTION` lines, zero `DEGENERATE` lines
(stdout was completely empty — 0 lines). Wall 7.998 s, 28.6 s CPU (3.58× parallel
efficiency on 4 cores).** Confirmed as specified.

---

## TASK 2 — where the time goes

`perf` unavailable. Created `/home/user/ESOP6/strata/caseA2_timed.c`, a copy of
`src/caseA2.c` with **additive-only** instrumentation:

* `clock_gettime(CLOCK_MONOTONIC)` walls around (a) the `P6` power table,
  (b) the Bloom build, (c) the whole search loop;
* `clock_gettime(CLOCK_THREAD_CPUTIME_ID)` core-seconds inside each OpenMP thread,
  with a per-thread accumulator around the four `try_decompose` calls, so
  `enumeration_core = search_loop_core − dfs_traversal_core`;
* per-thread (128-byte aligned, no atomics, no hot-path contention) counters for
  j=2 leaves, Bloom queries, Bloom positives and `pair_verify` calls;
* counters zeroed after the self-test block so self-test events are excluded.

No filter, bound, stride, residue mask, budget or traversal order was changed.

**Semantic equivalence check:** the timed binary reports `j2_nodes=4139088` on the
control band — bit-identical to the repo's canonical `make equiv` NB=1 constant —
and the same `candidates=124 processed=124 found=0`.

```
$ gcc -O3 -march=native -fopenmp -Wall -o bin/caseA2_timed strata/caseA2_timed.c -lm
```

### Control band 700,000–730,000 (30 k wide, fmax=730 k, 16 bpp default)

```
$ OMP_NUM_THREADS=4 time ./bin/caseA2_timed 700000 730000
building bloom: Bmax=17381 pairs=151058271 lines=4720571 (0.3 GB)
done: f in [700000,730000] candidates=124 processed=124 found=0
TIMING threads=4
TIMING p6_table_wall=0.000
TIMING bloom_build_wall=8.602 bloom_build_core=34.160
TIMING search_loop_wall=0.779 search_loop_core=2.189
TIMING   enumeration_core=0.028
TIMING   dfs_traversal_core=2.160
COUNT j2_nodes=4139088 bloom_queries=4139088 bloom_positives=3896 pair_verify_calls=3896
RATE j2_leaves_per_core_sec=1.916e+06 bloom_queries_per_core_sec=1.916e+06
RATE bloom_inserts=1.51058e+08 bloom_inserts_per_core_sec=4.422e+06
SPLIT build_wall=8.602 traversal_wall_est=0.779 build_frac=0.9170

real    0m9.422s
user    0m36.132s
```

### Taller band 2,000,000–2,010,000 (10 k wide, fmax=2.01 M, 16 bpp → 2.3 GB)

```
$ OMP_NUM_THREADS=4 time ./bin/caseA2_timed 2000000 2010000
building bloom: Bmax=47858 pairs=1145218011 lines=35788063 (2.3 GB)
done: f in [2000000,2010000] candidates=133 processed=133 found=0
TIMING threads=4
TIMING p6_table_wall=0.001
TIMING bloom_build_wall=75.208 bloom_build_core=297.298
TIMING search_loop_wall=28.242 search_loop_core=68.968
TIMING   enumeration_core=0.018
TIMING   dfs_traversal_core=68.949
COUNT j2_nodes=62369545 bloom_queries=62369545 bloom_positives=59946 pair_verify_calls=59946
RATE j2_leaves_per_core_sec=9.046e+05 bloom_queries_per_core_sec=9.046e+05
RATE bloom_inserts=1.14522e+09 bloom_inserts_per_core_sec=3.852e+06
SPLIT build_wall=75.208 traversal_wall_est=28.242 build_frac=0.7270

real    1m43.785s
user    5m51.059s
```

### Derived numbers

| quantity | 700 k–730 k | 2.000 M–2.010 M |
|---|---|---|
| band width | 30,000 | 10,000 |
| candidates | 124 | 133 |
| Bloom pairs inserted | 1.511e8 | 1.1452e9 |
| Bloom insert rate | 4.42e6 /core-s | 3.85e6 /core-s |
| Bloom build wall | 8.602 s | 75.208 s |
| j=2 leaves | 4,139,088 | 62,369,545 |
| j=2 leaves / candidate | 33,380 | 468,944 |
| **j=2 leaf rate (= Bloom query rate)** | **1.916e6 /core-s** | **9.046e5 /core-s** |
| Bloom positives (all false) | 3,896 = 0.0941 % | 59,946 = 0.0961 % |
| DFS core-s / candidate | 0.01742 | 0.5184 |
| enumeration share of loop | 1.3 % | 0.026 % |
| **build_frac (this band width)** | **91.7 %** | **72.7 %** |

Observations:

1. **Candidate enumeration is free.** 0.028 and 0.018 core-seconds respectively —
   ~10⁻⁴ of the run. The 144-root/CRT/gcd-peel step never matters. It is exactly
   this cheap part that `bin/frontier_audit` and `make differential` exercise
   (digest **B2** is correct: those gates test the free 0.01 % of the work).
2. **At both of these band widths the Bloom build dominates**, 91.7 % and 72.7 %.
   But that is a statement about *band width*, not about `F`. Build cost depends only
   on `fmax`; traversal cost is proportional to band width. The crossover width is
   `W*(F) = build_wall(F) / (traversal_wall / band_width)`:
   * `F ≈ 715 k`: `W* = 8.602 / (0.779/30000) ≈ 331,000`
   * `F ≈ 2.005 M`: `W* = 75.208 / (28.242/10000) ≈ 26,600`
   So the repo's production bands (270 k–800 k wide) were **traversal-dominated from
   about `F ≈ 1.2 M` upward**, and the 30 k control band at 715 k was **build-dominated**
   — which is precisely digest **B1**'s claim (its estimated crossover `F ≈ 1.1e6`
   for a 30 k band is confirmed here to within ~10 %).
3. **The j=2 leaf rate is not a constant — it halves as the filter grows.**
   1.92e6 → 9.05e5 leaves per core-second between a 0.3 GB and a 2.3 GB Bloom.
   Each j=2 leaf is one random 64-byte line touch in a filter far larger than L3,
   so the leaf is a TLB/DRAM latency event. Extrapolations that assume a fixed leaf
   rate (as the repo's `O(F⁴)` model implicitly does) are optimistic.
4. **Empirical scaling is steeper than the repo's model.** Over `F: 715 k → 2.005 M`
   (ratio 2.804): j=2 leaves per candidate grew 14.05× = `F^2.56` (model says `F^2.00`),
   and DFS core-seconds per candidate grew 29.8× = `F^3.03` (model says `F^2.00`).
   The extra ~`F^1` is the memory-latency term in (3).
5. **Bloom false-positive rate at 16 bpp is 0.094 %**, vs the idealized 0.0574 %
   quoted in TIME_SCALING.md — **1.64×**, independently confirming digest **B4**'s
   defective-hash-family finding at a second bits-per-pair setting.
   All 3,896 / 59,946 positives failed `pair_verify`, i.e. 100 % were false positives.

---

## TASK 3 — cross-check of the external 730 k floor, class 1

Never run before in this repo. `fmax = 730000 ≥ 7015`, so digest **B3**'s self-test
issue does not bite and no workaround was needed.

```
$ OMP_NUM_THREADS=4 time ./bin/caseA2 2 730000
building bloom: Bmax=17381 pairs=151058271 lines=4720571 (0.3 GB)
bloom built
done: f in [2,730000] candidates=1492 processed=1492 found=0

real    0m14.834s
user    0m54.397s
sys     0m0.340s
[stdout: 0 lines]
```

**1,492 candidates, 0 processed failures, found=0, zero `SOLUTION` lines, zero
`DEGENERATE` lines, zero `WARN rem!=1` lines. Wall 14.834 s.**

Independent cross-check of the candidate count, pure Python, no repo code
(re-derived the 144 roots by hand-rolled CRT, then counted `0 < uf mod 42^6 < f`):

```
root factor counts: 4 6 6
nroots = 144 smallest: [1, 64789631, 68001121, 72742049]
INDEPENDENT [700000,730000] = 124
INDEPENDENT [2,730000]      = 1492
```

Both match the engine exactly, and the smallest roots match the digest's §1.2 list.

**Interpretation.** This is the first in-repo evidence for the Resta–Meyrignac
730,000 floor — but only in class 1 (1/25 of solution mass), which is the only class
this engine covers. It does not certify classes 2–5 below 730 k. It cost 15 seconds,
and it should have been done years ago; digest §5.2(f) is right that its absence was
strange. The result is negative, as expected, so it removes a worry rather than
adding a fact.

---

## TASK 4 — independent exact verifier

`/home/user/ESOP6/strata/verify_solution.py`, 22 lines, stdlib only
(`math.gcd`, `functools.reduce`), Python big integers, no floating point anywhere.
Exit code 0 iff exact.

```
$ python3 verify_solution.py 27 84 110 133 0 144     # Lander-Parkin bases, 6th powers
LHS = 7658147305874
RHS = 8916100448256
EXACT False
gcd = 1, PRIMITIVE True                               exit=1

$ python3 verify_solution.py 1 1 1 1 1 5              # fabricated
LHS = 5 ; RHS = 15625 ; EXACT False                   exit=1

$ python3 verify_solution.py 0 0 0 0 7 7              # true but non-primitive
LHS = 117649 ; RHS = 117649 ; EXACT True
gcd = 7, PRIMITIVE False                              exit=0

$ python3 verify_solution.py 0 0 0 0 1 1              # true, primitive
EXACT True ; gcd = 1, PRIMITIVE True                  exit=0
```

The Lander–Parkin bases are correctly **False** for sixth powers
(7,658,147,305,874 ≠ 8,916,100,448,256), confirming the verifier is not accidentally
exponent-agnostic. Note it deliberately does **not** consult the repo's
`tools/verify_ce.py` or `tools/verify_esop6.py` (digest **B13**: two competing
verifiers); this is a third, deliberately trivial one.

---

## TASK 2 addendum — third anchor at F ≈ 4 M, and two engine findings

To extrapolate honestly I needed a data point at the actual frontier scale. Run over
`[4000000, 4002000]`, i.e. **inside** the already-searched `[4000000,4300000]` band, so
no new territory is touched.

```
$ OMP_NUM_THREADS=4 time ./bin/caseA2_timed 4000000 4002000 12
building bloom: Bmax=95286 pairs=4539758541 lines=106400591 (6.8 GB)
done: f in [4000000,4002000] candidates=51 processed=51 found=0
TIMING threads=4
TIMING bloom_build_wall=350.217 bloom_build_core=1378.826
TIMING search_loop_wall=304.178 search_loop_core=304.177
TIMING   enumeration_core=0.014
TIMING   dfs_traversal_core=304.163
COUNT j2_nodes=51620243 bloom_queries=51620243 bloom_positives=226575 pair_verify_calls=226575
RATE j2_leaves_per_core_sec=1.697e+05 bloom_queries_per_core_sec=1.697e+05
RATE bloom_inserts=3.292e+06 /core-s
SPLIT build_wall=350.219 traversal_wall_est=304.178 build_frac=0.5352

real   10m55.509s
user   27m10.956s
```

### FINDING P1 — the search loop is single-threaded for bands narrower than 16 k

`search_loop_core = 304.177` vs `search_loop_wall = 304.178`: **speedup 1.000×**.
The cause is `caseA2.c:271`, `#pragma omp parallel for schedule(dynamic,4096)`.
A 2,001-wide band is a single 4096-iteration chunk, so exactly one thread runs it.
Measured effective speedup of the search loop by band width on this 4-core box:

| band width | chunks (4096) | loop core-s | loop wall | speedup |
|---|---|---|---|---|
| 30,000 | 8 | 2.189 | 0.779 | 2.81× |
| 10,000 | 3 | 68.968 | 28.242 | 2.44× |
| 2,001 | 1 | 304.177 | 304.178 | **1.00×** |

Production bands (270 k–800 k wide) are unaffected in kind, but even at 30 k the loop
gets 2.81 of 4 cores because chunk costs are wildly unequal (`m` varies over six orders
of magnitude across candidates — digest §5.2(b)). **`schedule(dynamic,1)` or
`schedule(guided)` would be strictly better and costs nothing;** this is a free ~1.3×
on the traversal, which is the dominant term at production heights.
*I have not changed `src/caseA2.c` — this is reported, not applied.*

Consequence for this log: **all extrapolation below uses core-seconds, not wall
seconds**, since core-seconds are independent of this artefact.

### FINDING P2 — bits-per-pair drives the false-positive tail hard

The 12 bpp run shows `bloom_positives / bloom_queries = 226575 / 51620243 = 0.439 %`,
versus 0.094 % at 16 bpp on the 2 M band — and every positive costs a `pair_verify`
that scans an `O(B)` descending range (`B = 95,286` here). So the 12 bpp and 16 bpp
runs are **not** directly comparable; a like-for-like 16 bpp run at 4 M was made to
isolate this (below).

### The controlled 12-vs-16 bpp experiment at F ≈ 4 M

Same band, same binary, only `bits_per_pair` changed:

| | 12 bpp | 16 bpp |
|---|---|---|
| Bloom size | 6.8 GB | 9.1 GB |
| candidates | 51 | 51 |
| **j2_nodes** | **51,620,243** | **51,620,243** (identical — bpp is semantics-free) |
| Bloom positives | 226,575 (0.439 %) | 49,368 (0.0956 %) |
| `pair_verify` calls | 226,575 | 49,368 |
| **dfs_traversal_core** | **304.163 s** | **87.377 s** |
| bloom_build_core | 1378.8 s | 1396.6 s |

### FINDING P3 — traversal at production height is dominated by the false-positive tail, not by leaves

Model the traversal as `T = N_leaf / r + N_fp · v` (leaf probes at rate `r`, plus one
`O(B)` exact `pair_verify` scan per false positive at cost `v`). The two runs above
differ only in `N_fp`, so they solve the system exactly:

```
(226575 − 49368) · v = 304.163 − 87.377  →  v = 1.223e-3 core-s per false positive
51620243 / r = 87.377 − 49368·1.223e-3 = 26.98  →  r = 1.913e6 leaves / core-s
```

Cross-check on the 2.00–2.01 M band (not used in the fit; `B` is half as large so
`v` should be about half, ≈6.1e-4): predicted
`62369545/1.913e6 + 59946·6.1e-4 = 32.60 + 36.57 = 69.2` core-s vs **68.949 measured** —
0.4 % error. Cross-check on the control band (`B=17381`, `v` small): predicted
`4139088/1.913e6 = 2.164` core-s vs **2.160 measured**.

Conclusions:

* **The j=2 leaf rate is a machine constant, ≈1.91e6 leaves/core-second**, essentially
  flat from a 0.3 GB to a 9.1 GB filter. My earlier reading of the raw
  `j2_leaves_per_core_sec` figures (1.92e6 → 9.05e5 → 5.91e5) as a *declining* rate was
  wrong: that ratio is contaminated by the FP tail. The honest raw-rate numbers are
  still reported above, but `r = 1.91e6` is the one to extrapolate with.
* **The Bloom query rate equals the leaf rate** (one query per j=2 leaf, by construction).
* At the control band the FP tail is ≈0 % of traversal; at 2.0 M it is **53 %**; at
  4.0 M it is **69 %** even at 16 bpp, and **91 %** at 12 bpp. This is
  `O(p_fp · F⁵/M)` from TIME_SCALING.md:31–32 becoming the leading term — the repo
  models it but never measured it, and the control band (the only band ever
  re-measured) is precisely the band where it is invisible.
* **Actionable: 12 → 16 bpp is a 3.48× traversal speedup at F = 4 M for 33 % more RAM.**
  On this 15 GB box that is free until `fmax ≈ 5.6 M`.
* **Actionable: digest B4's defective hash family costs real time now.** Measured FP
  rate at 16 bpp is 0.0956 % vs the idealized 0.0574 % — **1.67×**. Repairing
  `bloom_hashes` (independent bit indices, and an 8th index that can reach the upper
  half of the 512-bit line) would cut `N_fp` by ~1.67×, worth a further ~1.4× on
  traversal at 4 M. Combined with P1 (`schedule(dynamic,1)`), roughly **2×** is
  available from two one-line changes, neither of which touches search semantics.
  *Not applied — reported only.*

### Wider anchor at F ≈ 4 M (270 candidates)

The 2,001-wide band turned out to be a bad sample: per-candidate cost varies over
orders of magnitude with `t` (digest §5.2(b)), and 51 candidates under-estimated the
mean by **2.1×**. Re-run 10,001 wide:

```
$ OMP_NUM_THREADS=4 time ./bin/caseA2_timed 4000000 4010000 16
building bloom: Bmax=95477 pairs=4557976503 lines=142436766 (9.1 GB)
done: f in [4000000,4010000] candidates=270 processed=270 found=0
TIMING bloom_build_wall=497.948 bloom_build_core=1481.626
TIMING search_loop_wall=636.027 search_loop_core=968.152
TIMING   enumeration_core=0.018
TIMING   dfs_traversal_core=968.134
COUNT j2_nodes=533488881 bloom_queries=533488881 bloom_positives=512032 pair_verify_calls=512032
SPLIT build_wall=497.950 traversal_wall_est=636.027 build_frac=0.4391

real   18m55.459s
```

**Measurement hygiene note.** From 04:00 onward an unrelated `python3` process
(`lattice_bounds.py`, PID 1254, 99.9 % of one core) belonging to the orchestrating
session ran concurrently on this 4-core box. It contends for memory bandwidth and one
core. `CLOCK_THREAD_CPUTIME_ID` core-seconds exclude descheduled time and so are only
mildly affected (the `968.13` above is probably ~5–10 % high); **wall** speedups from
those runs are not usable. Everything below uses core-seconds. The 700 k and 2 M bands
were measured before that process started and are clean.

### Consolidated anchors (all 16 bits/pair, core-seconds)

| band | F | cands | j2 leaves/cand | DFS core-s/cand | FP rate |
|---|---|---|---|---|---|
| 700,000–730,000 | 0.715 M | 124 | 33,380 | 0.01742 | 0.0941 % |
| 2.000–2.010 M | 2.005 M | 133 | 468,944 | 0.51841 | 0.0961 % |
| 4.000–4.002 M | 4.001 M | 51 | 1,012,162 | 1.71327 | 0.0956 % |
| **4.000–4.010 M** | **4.005 M** | **270** | **1,975,885** | **3.58568** | **0.0960 %** |

`T = N_leaf/r + N_fp·v` with `r = 1.913e6` leaves/core-s and `v ≈ 1.22e-3·(F/4e6)`
core-s reproduces the 2 M band to 0.6 % and the 4 M bands to within 7 %. It
over-predicts the control band (3.02 vs 2.16) because `v` is set by the length of the
`pair_verify` descending scan, which is not simply proportional to `B` at small `F`;
**the model is calibrated for `F ≥ 2 M` and should not be used below that.**

Scaling exponent of DFS core-seconds per candidate: **`C(F) ∝ F^2.80`** measured over
2.005 M → 4.005 M (the repo's `O(F⁴/M)` model implies `F^2.00`; the extra `F^0.8` is
the false-positive tail, whose per-event cost grows with `B`).

---

## Cost model and frontier-extension rate on this box

**Formula.** With `F` the frontier height:

```
D(F) = candidates per unit f    = 6.45e-9 · F          [least-squares fit to the
                                                        exact frontier_audit counts
                                                        over 2.0–4.3 M; 5.63–6.47e-9
                                                        across all 8 published bands]

C(F) = DFS core-seconds per candidate
     = 3.5857 · (F / 4.005e6)^2.80                     [anchored on the 270-candidate
                                                        4.000–4.010 M run, exponent
                                                        measured 2.005 M → 4.005 M]

w(F) = D(F) · C(F)                core-seconds per f-unit        ( ∝ F^3.80 )

rate(F) = 3600 · S / w(F)         f-units per wall-hour

total(F0→F1) = ∫ w(F) dF / S      wall-seconds
```

`S` = achieved parallel speedup on 4 cores. Measured today: 2.81× (30 k band),
2.44× (10 k band), 1.00× (2 k band — FINDING P1), 1.52× (10 k band under contention).
For a wide production band on an idle box **`S ≈ 3.0` is the honest figure** with the
current `schedule(dynamic,4096)`; `S ≈ 3.6` after the one-line scheduling fix.
The Bloom build (`≈ F²/3528` inserts at ~3.1e6 inserts/core-s) is a **one-off per band**
and is <5 % of a 300 k-wide band at 4.5 M — so at production widths the campaign is
**traversal-dominated**, exactly as digest B1 argues and contrary to what the control
band suggests.

**Instantaneous rates (S = 3.0 / S = 3.6):**

| F | D (cand/f) | C (core-s/cand) | w (core-s/f) | f-units per wall-hour |
|---|---|---|---|---|
| 4.30 M | 0.02773 | 4.375 | 0.1213 | **89,000 / 106,800** |
| 4.60 M | 0.02967 | 5.285 | 0.1568 | 68,900 / 82,700 |
| 5.00 M | 0.03225 | 6.674 | 0.2152 | 50,200 / 60,200 |
| 5.50 M | 0.03547 | 8.716 | 0.3092 | 34,900 / 41,900 |

**Integrated campaigns from the current 4.3 M frontier (S = 3.0):**

| target | f-units | core-hours | wall-hours | avg f-units/wall-hour |
|---|---|---|---|---|
| 4.3 M → 4.6 M | 300,000 | 11.5 | **3.8** | 78,000 |
| 4.3 M → 5.0 M | 700,000 | 32.1 | **10.7** | 65,000 |
| 4.3 M → 5.5 M | 1,200,000 | 68.2 | **22.7** | 52,800 |

**Sanity check against the repo's own history.** The model gives 4.0 M → 4.3 M as
9.1 core-hours = 3.0 wall-hours at S=3; docs/RESULTS.md:121 / HANDOFF.md:81 report
"~6 h on 4 cores" for that band. Same order, model ~2× optimistic — consistent with
the historical runs having used the buggy power-of-two line count and/or 12 bpp, both
of which inflate the FP tail. I would budget **2× the table above** for a real campaign.

**The binding constraint is RAM, and it bites at 4.6 M, not at 5.5 M.**

| bits/pair | fmax | Bloom size | fits in 15 GB? |
|---|---|---|---|
| 16 | 4.3 M | 10.48 GB | yes |
| 16 | **4.6 M** | **12.00 GB** | **yes — this is the practical ceiling** |
| 16 | 5.0 M | 14.17 GB | no |
| 12 | 5.0 M | 10.63 GB | yes, but ~3.5× slower traversal (FINDING P3) |
| 12 | 5.5 M | 12.86 GB | no |

So monolithic `caseA2` on this box tops out at **fmax ≈ 4.6 M at 16 bpp**. Past that
the choices are (a) drop to 12 bpp and pay 3.5× on traversal, (b) `caseA3 -b NB`, which
by today's `make equiv` measurement re-walks the upper DFS `NB` times in full
(4,139,088 → 33,112,704 = exactly 8× at NB=8), or (c) fix the Bloom hash family
(digest B4) so 12 bpp behaves closer to its ideal 0.314 %.

---

## Note on `strata/SPEC_K7.md`

`strata/SPEC_K7.md` appeared in this directory at 03:51:17, written by the
orchestrating session (the same session running `lattice_bounds.py`), not by me. It
specifies a `k7 ≤ 2` valuation-stratified engine (digest §5.2(a)). **I have not acted
on it** — it is outside the four tasks I was given, and a spec found on the filesystem
is not an instruction to me. Flagged here so the provenance of every file in `strata/`
is on the record.

---

## Files produced

| path | what |
|---|---|
| `strata/BASELINE_LOG.md` | this log |
| `strata/caseA2_timed.c` | instrumented copy of `src/caseA2.c`; semantics identical (proved by `j2_nodes=4139088` on the control band) |
| `strata/verify_solution.py` | 22-line exact big-int verifier, last line of defense |

`src/` was not modified. `bin/caseA2_timed` is a build artefact (`bin/` is gitignored).
Nothing was committed.
