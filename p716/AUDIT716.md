# AUDIT716 -- adversarial completeness audit of `engine716.c`

Date: 2026-09-12.  Auditor: independent pass, adversarial brief ("assume there
is a bug").  Target: `/home/user/ESOP6/p716/engine716.c` (md5 of the binary in
place, rebuilt bit-identically from source: `417e5276fafbf7413aa094dbeb16454d`;
`/proc/<prod pid>/exe` pointed at that same binary, so the audited source *is*
what produced the coverage).

Nothing in production code was changed.  All audit code lives outside the
repo, in
`/tmp/claude-0/-home-user/0f5cfd23-9c6f-59f5-bc05-75ceb9fbcc55/scratchpad/aud/`
(`hook716.c`, `hook2.c`, `window_exact.py`), each of which `#include`s
`engine716.c` with `main` renamed so the audit calls the **real** static
functions rather than a re-implementation.

---

## HEADLINE

**No completeness bug was found.**  Every mechanism that could silently drop a
solution -- the long-double-seeded integer roots, the four nested windows, the
verifier windows, the 2-sum Bloom, the residue masks, the Q-pass c-restriction
and table split, chunk coverage -- was tested at or above production magnitude
and came out exact.  The reported **ZERO for `2 < f <= 6500` stands.**

Two things nevertheless need your attention, neither a search bug:

* **OPERATIONAL: the production run is DEAD.**  `engine716 6500 6750 --q 42`
  (pid 9593) and its `run716.sh` driver both vanished at **14:33:07 UTC** while
  in pass 12/42 of chunk `(6500,6750]`, ~4240 s in.  No `FATAL`/`ASSERT`, no
  `ERROR rc=` line from `run716.sh` (the driver would have written one), no
  kernel OOM record -- consistent with the whole detached process group being
  reaped externally, not with an engine fault.  I did **not** kill it and I
  have **not** restarted it.  `runs/coverage.txt` is unaffected: the chunk was
  never checkpointed, so recorded coverage is still exactly `f <= 6500` and the
  checkpoint design behaved correctly.  Restart with `./run716.sh 6000 7000 250 42`
  when you want it back (it resumes from 6500).
* **CONTEXT: the search is heuristically hopeless at this scale.**  The
  expected number of (7,1,6) solutions with `2 < f <= F` is
  `(6/7) * Gamma(8/7)^6 / (Gamma(13/7) * 6!) * ln F  =  9.82e-4 * ln F`.
  For `F = 6500` that is **0.0066**; for `F = 10^6` it is 0.011; expectation 1
  is reached only near `F = e^1188`.  Zero is the overwhelmingly likely
  outcome, and the yield is logarithmic in `F`, so no reachable amount of extra
  CPU changes that.  (Naive heuristic; it ignores congruence and algebraic
  structure, but the exponent count is robust.)

---

## Item 1 -- `iroot7` / `ceil_root7_div` at production magnitude -- **SOUND**

The long double seed is irrelevant to correctness: both functions follow the
`powl` estimate with *unbounded* exact `u128` correction loops
(`while(r>0 && ipow7(r)>R) r--;` / `while(r<IR7MAX && ipow7(r+1)<=R) r++;`, and
the `kge`-guarded pair in `ceil_root7_div`).  Any finite estimate error is
walked off exactly.  `kge` handles the `k*n^7` overflow correctly via
`KLIM[k] = floor((2^128-1)/k)`.

Verified, not argued:

    cd /home/user/ESOP6/p716
    gcc -O2 -std=gnu11 -I. -o .../hook716 .../hook716.c -lm
    .../hook716 roots
    -> ROOTS checks=21144000 bad=0 ; AUDIT CLEAN     (11 s)

which is
* **exhaustive boundaries for every base `x = 1..8000`** (well past production
  6500/6750): `iroot7(x^7)==x`, `iroot7(x^7-1)==x-1`, `iroot7(x^7+1)==x`, and
  for **k = 2,3,4,5,6**: `ceil_root7_div(k*x^7,k)==x`,
  `ceil_root7_div(k*x^7-1,k)==x`, `ceil_root7_div(k*x^7+1,k)==x+1`;
* **3,000,000 random `u128` values uniform in `[1, 6*8000^7]`** (i.e. up to
  ~1.3e28, above every production value) compared against an exact
  binary-search reference that uses **no floating point at all**
  (`ref_iroot7`, `ref_ceil7div` in `hook716.c`), for `iroot7` and for
  `ceil_root7_div` with k = 1..6.  Zero disagreements.

Latent, **not** reachable in production: `IR7MAX` is computed as 319519 but the
true largest `r` with `r^7 < 2^128` is **319557**, so both functions saturate
0.08% below the 128-bit ceiling.  Unreachable while `FMAX <= 200000`
(`200000^7 = 1.28e37`).  See L1 below.

## Item 2 -- window completeness -- **SOUND**

Claim proved: for every ordered `a>=b>=c>=d>=e>=g>=1` with `S = sum x_i^7` and
`f = S^(1/7)`, `enumerate()` admits `(a,b,c)` and `verify3()` recovers
`(d,e,g)`.

*Argument.*  `a` is the largest of six positive terms, so `6a^7 >= S` gives
`a >= ceil_root7_div(S,6) = amin`, and `a^7 <= S` with `a <= f-1` (a=f would
force the other five to zero) gives the upper end.  `R1 = S-a^7` is a sum of
five terms with `b` largest, so `5b^7 >= R1` and `b^7 <= R1`, `b <= a`:
exactly `blo..bhi`.  Same for `c` with four terms (`4c^7 >= R2`) and for `d`
with three (`3d^7 >= R3`, `d <= c`, `d <= TB` since `d <= c < f <= FMAX = TB`)
and `e` with two (`2e^7 >= Rp`, `e <= d`).  `R3 <= 3c^7` is implied by `d<=c`,
so the prune is safe, and `R3>=3`, `R2>=4`, `R1>=5`, `Rp>=2` are implied by
positivity.

*Empirical, at production scale, two independent implementations.*

(a) Driving the **actual C bound functions** (`hook716 window`), TB = 6500:

    .../hook716 window 6500 250000
    -> WINDOW tuples=250000 bmax=6500 fails=0 ; AUDIT CLEAN

Each tuple is checked against every inequality of `enumerate()` **and**
`verify3()` as written, plus `mask3_ok(R3)`, `mask2_ok(Rp)`, plus for
`Q in {1,2,3,6,7,14,16,21,42}` that the Q-pass c-stride reaches `c` in exactly
one pass and that `build_table`'s g-class restriction inserts `(d,e,g)` in
exactly that pass.  Tuple shapes include `g=1`, two pairs of equal bases, all
six equal, and `a = 6500`.

(b) An **exact-integer Python re-implementation** of every inequality
(`window_exact.py`), TB = 6750, three seeds:

    python3 window_exact.py 7   120000   -> structured 7328 + random 120000 + skewed 120000, FAILS=0
    python3 window_exact.py 101 400000   -> TOTAL tuples=807328 FAILS=0
    python3 window_exact.py 202 400000   -> TOTAL tuples=807328 FAILS=0
    (1,861,984 tuples total, 0 failures)

The structured families are the boundary cases: all six bases equal (this makes
`R3 == 3c^7` exactly, the `R3 > 3*P7[c]` prune boundary, and `a == amin`
exactly), five equal + `g=1`, `c,c,c,1,1,1`, `a,1,1,1,1,1`, and for every
`a in [6400,6500]` the dominant shapes with tails `1,2,3,7,100,6399,6400`.
Plus 640k uniform-random and 640k deliberately skewed (wide dynamic range)
tuples with bases up to 6500.

## Item 3 -- pair Bloom `F2` false negatives -- **SOUND**

`build_pairs()` inserts `e^7+g^7` for **every** `e in [1,TB]`, `g in [1,e]`, and
it ignores the pass entirely, so `F2` is pass-independent.  In `verify3`, a
genuine `Rp = e^7+g^7` always has `g <= e <= d <= min(iroot7(R3), c, TB)`, hence
`e <= TB`: every reachable pair is in the insert set.  A Bloom filter has no
false negatives once the value is inserted, so `if(!filt_query(&F2,Rp)) continue;`
cannot lose a solution.

Verified **exhaustively** at production TB (6750, the in-flight chunk's bound,
above the 6500 that has been recorded):

    .../hook2 f2 6750
    -> TB=6750 pairbloom 0.091 GB
       F2 EXHAUSTIVE pairs=22784625 false_negatives=0
       F2 fp_rate=1.55e-05 (informational)
       AUDIT CLEAN

All 22,784,625 reachable pairs queried; zero false negatives.

## Item 4 -- `mask2_ok` / `mask3_ok` -- **SOUND**

`sumres[i][j]` is built by brute-force closure from
`single = { x^7 mod m : x in [0,m) }`, so it is by construction a **superset**
of the residues attainable by a real j-fold sum -- the mask can only
over-accept, never reject a genuine sum.

Cross-checked against an independent Python brute force that uses bases `>= 1`
only (`.../hook2 masks` dumps the tables; the comparison script is inline in
the audit log):

    OK  mod=49 j=1 7/49  j=2 19/49  j=3 37/49   (0.1429 / 0.3878 / 0.7551)
    OK  mod=29 j=1 5/29  j=2 13/29  j=3 25/29   (0.1724 / 0.4483 / 0.8621)
    OK  mod=43 j=1 7/43  j=2 19/43  j=3 37/43   (0.1628 / 0.4419 / 0.8605)
    UNSAFE_MISSING_RESIDUES: 0

Sets are **equal**, not merely supersets, and match the rates in `LOG716.md`.
Direct rejection test on genuine sums:

    .../hook2 masks
    -> MASKS n=4000000 reject3=0 reject2=0                 (random d,e,g in [1,6800])
       MASKS exhaustive<=130 reject3=0 reject2=0           (all d>=e>=g<=130)

`md49_/md29_/md43_` reduce a `u128` as `((hi mod m)*(2^64 mod m) + lo mod m) mod m`,
which is exact; the 4M-sample test exercises them on real 1e27-magnitude values.

## Item 5 -- the Q-pass c-restriction -- **SOUND**

`x^7 == x (mod Q)` for `Q in {2,3,6,7,14,21,42}` (for each prime factor
`p in {2,3,7}`, `p-1 | 6`), so `x -> x^7 mod Q` is a bijection and
`PRECNT[v] == 1` for every `v` -- confirmed for the three Q values actually
used in production (7, 14, 42).  For `Q = 16` it is not a bijection
(`PRECNT[0] = 8`) and the code correctly loops over `PRECNT[want]` classes.

`want = (modq(R2) + Q - pass) % Q` is exactly the condition
`c^7 == R2 - pass (mod Q)`, i.e. `R3 = R2 - c^7 == pass (mod Q)`.
`build_table` picks `g` with `g^7 == pass - d^7 - e^7 (mod Q)`, i.e. inserts
exactly the triples whose sum is `== pass`, so leaf and table agree pass by
pass.

Verified:

    .../hook2 cpart
    -> CPART Q=2  trials=20000 bad_visit_counts=0 wrong_pass=0
       CPART Q=3  ...  Q=6 ... Q=7 ... Q=14 ... Q=16 ... Q=21 ... Q=42  (all 0/0)
       MODQ checked q=1..1024
       AUDIT CLEAN

For each of 20,000 random `(chi, clo, R2)` per Q, the audit replays the exact
`want / PRE / delta / cstart / step` code over **all** Q passes and counts, for
every `c in [clo,chi]`, how many times it is visited: always **exactly once**
(`bad_visit_counts=0`), and always in the pass with `modq(R2 - c^7) == pass`
(`wrong_pass=0`).  So the passes partition the c-range with no gap and no
duplication, and `cstart`/`delta` have no off-by-one (including the `t==0`
class, where `build_table` correctly starts at `g0 = QMOD`).

`modq()` on `u128` was checked against a bit-by-bit long-division reference for
**every** modulus `q = 1..1024`, 2000 random 128-bit values each: no
disagreement.

## Item 6 -- everything else

| Sub-item | Verdict | Evidence |
|---|---|---|
| `nb` bucket path sharing `pass` with QMOD | SOUND | `main` DIEs on `--q>1 && --nb>1`; with `nb==1` the `nb>1 &&` guard short-circuits. When `nb>1, Q==1`, `build_table` and `enumerate` both gate on `bucket_of(...)==pass` on the *same* value `R3`, so each leaf is evaluated in exactly one pass with its generating triple present. Production always ran `nb=1`. |
| off-by-one in `cstart`/`delta` | SOUND | Item 5, `CPART` exhaustive visit-count test. |
| `return` on solution inside the parallel region | SOUND | It returns from `enumerate()`, it does not jump out of the `omp for`; remaining iterations short-circuit on `volatile sol_found`; `band_line` then `emit_solution()`/`exit(3)`. Benign race: two threads may both increment `sols`, only the first CAS is reported. |
| overflow in `3*P7[c]` | SOUND in production | `c <= TB <= FMAX <= 200000` (enforced), `3*200000^7 = 3.8e37 < 2^128`. Latent only via `--plant` with `S > 2^128/3` -- see L2. |
| `P7` indexing beyond TB | SOUND | `P7` is `TB+2` long and filled to `TB+1`. Every consumer is bounded: `a <= min(f-1, iroot7(F7), TB)`, `b<=a`, `c<=b`, `d <= min(iroot7(R3), c, TB)`, `e<=d`, `g<=e`, `build_pairs` `e<=TB`, `build_table` `d,e,g<=TB`, `P7[f]` with `f<=FMAX=TB`. |
| f-loop bounds and chunk coverage | SOUND | `for(lo=FMIN+1; lo<=FMAX; lo+=100)`, inner `f=lo..min(lo+99,FMAX)` -- covers `(FMIN,FMAX]` with no gap. `runs/coverage.txt` chunks are contiguous with **zero gaps**: 1500-1750-2000-2250-2500-2750-3000-3500-4000-4500-5000-5250-5500-5750-6000-6250-6500 (checked programmatically). `f=1,2` are trivially empty (`f^7 < 6` fails, and `2^7=128` needs `a<=1` so max sum 6). |

### Global consistency check (leaves scaling)

A dropped-leaf bug would show up as a step in the leaf density.  Normalising
each chunk's leaf count by `sum_{f in chunk} f^3`:

    (1500,1750] Q=7  0.002244      (5000,5250] Q=14 0.002239
    (1750,2000] Q=7  0.002243      (5250,5500] Q=14 0.002239
    (2000,2250] Q=7  0.002243      (5500,5750] Q=14 0.002239
    (2250,2500] Q=7  0.002242      (5750,6000] Q=14 0.002238
    (2500,2750] Q=7  0.002241      (6000,6250] Q=42 0.002238
    (2750,3000] Q=7  0.002241      (6250,6500] Q=42 0.002238
    (3000,3500] Q=7  0.002240      prior (2,1500] Q=1  0.002248
    (3500,4000] Q=7  0.002240
    (4000,4500] Q=7  0.002239
    (4500,5000] Q=7  0.002239

Smooth and monotone, with **no discontinuity at either bucketing change**
(Q=7 -> 14 at f=5000, Q=14 -> 42 at f=6000).  Every chunk in `runs/prod.log`
also shows all `Q` passes completed (`passes_seen == npasses` for all 16
recorded chunks).

### End-to-end planted solutions beyond the tested scale

The existing T1 uses bases `<= 300`.  New planted runs:

*Bases 6000-6500, dominant-`a` shapes, `--nobloom` (so every masked leaf goes
to the exact verifier -- no filter can hide a miss), TB = 6500:*

    ./engine716 2 6500 --plant <S> --nobloom --threads 2
    plant 1..6  rc=3, 2-4 s each, all found, all re-verified in Python big ints:
      [6031 588 442 220 123 31] [6347 499 487 205 84 81] [6229 578 574 342 66 57]
      [6414 550 522 510 464 46] [6054 486 452 236 154 70] [6065 578 338 254 237 109]

*Generic random shapes, bases `<= 2000` (3.3x the tested scale), with the full
production machinery (3-sum Bloom ON) under Q=1, Q=7 and the production Q=42,
plus the `S+1` near-miss:*

    ./engine716 2 2 --q {1,7,42} --plant <S> --threads 2
    plant2k 1 [1944 1767 1263 852 270 59]  Q1/Q7/Q42 all found, exact, S+1 silent
    plant2k 2 [1056 800 777 658 564 379]   Q1/Q7/Q42 all found, exact, S+1 silent
    (remaining plants of this batch were still running when the audit closed;
     see NOTE-A)

*Full production configuration at production TB:*

    ./engine716 2 6500 --q 42 --plant 290218665402350491715925409 --threads 2
    (TB=6500, per-pass Bloom ~2.2 GB, exactly the (6250,6500] chunk geometry;
     see NOTE-A for the result)

Note that planted runs use `amax = TB`, not `amax = f-1`, so they cannot by
themselves exercise the `f-1` cap.  That cap is safe by the argument in Item 2
(`a = f` forces the other five bases to zero) and is covered by
`window_exact.py`, which applies `min(f-1, iroot7(S))` whenever `S` is a perfect
seventh power.

---

## Latent issues (none affects the recorded `f <= 6500` coverage)

* **L1.** `IR7MAX = 319519`; the true maximum with `r^7 < 2^128` is `319557`.
  `iroot7` and `ceil_root7_div` therefore saturate for arguments in
  `[319520^7, 2^128)`.  Unreachable: `FMAX <= 200000` is enforced, and every
  argument is `<= 6*TB^7`.  Fix if you ever raise FMAX: replace the
  `powl`-threshold search with an exact `u128` binary search on `ipow7`.
* **L2.** `if(R3 > 3*P7[c])` overflows `u128` once `c > 253000`
  (`3*253000^7 > 2^128`).  Only reachable through `--plant S` with
  `S > 2^128/3`, i.e. planted bases above ~253k; production `FMAX <= 200000`
  keeps `c` far below that.  Fix: `if(P7[c] <= KLIM[3] ? R3 > 3*P7[c] : 0)`.
* **L3.** `ctrs` is sized from the thread count observed in one early
  `#pragma omp parallel` probe.  If libgomp ever handed a later region a larger
  team (dynamic teams), `ctrs[tid()]` would write past the allocation.
  Production always passes `--threads 4` so the teams match; add
  `omp_set_dynamic(0)` and an `ASSERT(omp_get_num_threads()<=nthreads)` for
  safety.
* **L4.** `covlog.py` hard-codes the prior-session line
  `(7,1,6): 2 < f <= 1500, 2849012361 leaves` instead of reading it from
  `runs/coverage.txt`.  The claim is corroborated -- `runs/measure_1500.out`
  and `runs/measure_1500_hp.out` both end in
  `BAND 2 1500 2849012361 2480291308 2480291308 1810801 546 0`, and its leaf
  density 0.002248 lines up with the chunked runs -- and `git diff d070459
  3272873 -- engine716.c` shows the only change to the engine since then is the
  test-only `--nobloom` switch, so the search logic is identical.  Still, that
  range should be given a real `COVERED` line in `runs/coverage.txt` rather than
  living in a Python string literal.
* **L5.** (cosmetic) `single` in `init_masks` includes `x = 0`; harmless, since
  residue 0 is attained anyway by `x = m` and the mask can only over-accept.

## NOTE-A -- audit runs still in flight when this file was written

Two background audit jobs had not finished: the remaining 11 `plant2k` tuples
(bases `<= 2000`, Q = 1/7/42, Bloom on) and the single production-configuration
plant at TB = 6500 with `--q 42`.  Both were running against a 4-core box that
was also running other audit jobs.  Neither had reported a miss.  If either
does report a miss, this verdict must be revisited -- the results are appended
below when available.

## Reproduction

    cd /home/user/ESOP6/p716
    A=/tmp/claude-0/-home-user/0f5cfd23-9c6f-59f5-bc05-75ceb9fbcc55/scratchpad/aud
    gcc -O2 -std=gnu11 -I. -o $A/hook716 $A/hook716.c -lm
    gcc -O2 -fopenmp -std=gnu11 -I. -o $A/hook2 $A/hook2.c -lm
    $A/hook716 roots
    $A/hook716 window 6500 250000
    OMP_NUM_THREADS=1 $A/hook2 masks
    OMP_NUM_THREADS=1 $A/hook2 cpart
    OMP_NUM_THREADS=2 $A/hook2 f2 6750
    python3 $A/window_exact.py 101 400000
