# Raw run log

Verbatim output from the search programs, in chronological order. Machine:
4 cores, 15 GB RAM, gcc 13.3 with `-O3 -march=native -fopenmp`.

Each `done:` line reports the range, the number of candidate (f,t) pairs that
survived the modular sieve, the number actually processed (equal, always — a
mismatch would indicate a dropped candidate), and the solutions found.

## Validation: fifth-power mode

The same DFS engine, run at k = 5, must rediscover Lander & Parkin's 1966
counterexample. It does, in ~30 ms:

```
$ ./search 5 2 150
SOLUTION k=5 f=144 : 133 110 84 27
done: k=5 f in [2,150], solutions found: 1
```

Check: 27⁵ + 84⁵ + 110⁵ + 133⁵ = 14348907 + 4182119424 + 16105100000 +
41615795893 = 61917364224 = 144⁵. ✓

## Validation: completeness audit

Brute-force count of t < f with t⁶ ≡ f⁶ (mod 42⁶), compared against the
144-root enumeration used by the sweep:

```
$ ./audit
nroots=144
f=730001 brute=0 fast=0 OK
f=800003 brute=1 fast=1 OK
f=899999 brute=0 fast=0 OK
f=999983 brute=0 fast=0 OK
f=1299983 brute=0 fast=0 OK
f=54321011 brute=2 fast=2 OK
```

No candidates are missed.

## Exhaustive search, all cases

```
$ ./search 6 2 2000
done: k=6 f in [2,2000], solutions found: 0

$ ./search 6 2000 4000
done: k=6 f in [2000,4000], solutions found: 0

$ ./search 6 2 4000        # after optimization pass
done: k=6 f in [2,4000], solutions found: 0

$ ./search 6 4000 5000
done: k=6 f in [4000,5000], solutions found: 0
```

Runtime grew from 3.4 s (f ≤ 2000) to 85 s for the 4000–5000 band alone,
confirming the ~f^3.3 scaling that makes exhaustive search impractical past
the low tens of thousands on this hardware.

## Concentrated case: control run

Inside the published exhaustively-searched region, so it must find nothing:

```
$ ./caseA2 700000 730000
building bloom: Bmax=17381 pairs=151058271 lines=8388608 (0.5 GB)
bloom built
done: f in [700000,730000] candidates=124 processed=124 found=0

real	0m7.736s
```

## Concentrated case: past the frontier

```
$ ./caseA2 730000 1000000
building bloom: Bmax=23810 pairs=283469955 lines=16777216 (1.1 GB)
bloom built
done: f in [730000,1000000] candidates=1314 processed=1314 found=0

$ ./caseA2 1000000 1500000
building bloom: Bmax=35715 pairs=637798470 lines=33554432 (2.1 GB)
bloom built
done: f in [1000000,1500000] candidates=3523 processed=3523 found=0

$ ./caseA2 1500000 2000000
building bloom: Bmax=47620 pairs=1133856010 lines=67108864 (4.3 GB)
bloom built
done: f in [1500000,2000000] candidates=5129 processed=5129 found=0

$ ./caseA2 2000000 2500000
building bloom: Bmax=59524 pairs=1771583050 lines=67108864 (4.3 GB)
bloom built
done: f in [2000000,2500000] candidates=7222 processed=7222 found=0

$ ./caseA2 2500000 3200000
building bloom: Bmax=76191 pairs=2902572336 lines=134217728 (8.6 GB)
bloom built
done: f in [2500000,3200000] candidates=12443 processed=12443 found=0

$ ./caseA2 3200000 4000000 12
building bloom: Bmax=95239 pairs=4535281180 lines=134217728 (8.6 GB)
bloom built
done: f in [3200000,4000000] candidates=18003 processed=18003 found=0

$ ./caseA2 4000000 4300000 12
building bloom: Bmax=102381 pairs=5240985771 lines=134217728 (8.6 GB)
bloom built
done: f in [4000000,4300000] candidates=8050 processed=8050 found=0
```

## Totals

| | |
|---|---|
| Candidates examined past the published frontier (730k–4.3M) | **55,684** |
| Solutions found | **0** |
| Pair sums indexed at peak | 5.24 billion |
| Wall time for the campaign | ~18 hours on 4 cores |

## Incidents

- The 4.0–4.3M chunk died silently mid-run on its first attempt — the process
  vanished without a `done:` line, coinciding with a container disruption. No
  data was lost (completed ranges were already committed) and the chunk was
  re-run from scratch to completion. This is the motivating case for the
  checkpointing upgrade described in HANDOFF.md §5.
- No `WARN` lines were emitted at any point. The `rem != 1` guard in the
  factorization path (which would indicate 42⁶ failing to divide f⁶ − t⁶,
  i.e. a broken candidate) never fired.

---

# Addendum, 2026-09: independent reproduction and two fixes

Duncan cloned the work to a Windows/WSL box and reproduced both sanity checks
independently, on different hardware and a different toolchain:

```
./search 5 2 150        -> 27^5+84^5+110^5+133^5 = 144^5 found in 18 ms
./caseA2 700000 730000  -> 124 candidates, 0 solutions, 5.2 s
```

Both match the recorded values exactly.

## Fix 1 — Bloom filter was silently over-allocating

`caseA2` rounded the filter's line count up to a power of two, so the
requested bits-per-pair was only a lower bound. At fmax = 5×10⁶ with 12
bits/pair requested:

```
requested 12 bpp        -> 10.6 GB of bits
power-of-two rounding   -> 268435456 lines = 17.2 GB  (effective 19.4 bpp)
```

Confirmed by direct computation. Sizing is now exact, using Lemire range
reduction (`(h * nlines) >> 64`) in place of power-of-two masking. The control
run is unchanged afterwards:

```
$ ./bin/caseA2 700000 730000
building bloom: Bmax=17381 pairs=151058271 lines=4720571 (0.3 GB)
done: f in [700000,730000] candidates=124 processed=124 found=0
```

## Fix 2 — the memory wall is gone (`src/caseA3.c`)

The pair table is partitioned by a hash of the *sum* into NB buckets. Pass k
inserts only the pairs whose sum lands in bucket k and answers only the
queries whose sum lands in bucket k, so each query is answered in exactly one
pass and the filter holds 1/NB of the table.

Equivalence is verified by instrumented node counts rather than assumed —
`make equiv`:

```
=== NB=1 (monolithic) ===
Bmax=17381 pairs=151058271 buckets=1 filter=0.30 GB (0.3 GB monolithic)
j2_nodes=4139088 skipped=0 evaluated=4139088
done: f in [700000,730000] candidates=124 processed=124 found=0

=== NB=8 (bucketed) ===
Bmax=17381 pairs=151058271 buckets=8 filter=0.04 GB (0.3 GB monolithic)
j2_nodes=33112704 skipped=28973616 evaluated=4139088
done: f in [700000,730000] candidates=124 processed=124 found=0
```

33,112,704 − 28,973,616 = 4,139,088 evaluated, identical to the monolithic
count to the digit, and both report 0 solutions. Filter memory drops 8×.
These figures reproduce Duncan's local prototype exactly.

Resulting reach, at 12 bits/pair:

| f | caseA2 | caseA3 ×8 | caseA3 ×16 |
|---|---|---|---|
| 3.2M | 4.4 GB | 0.54 GB | 0.27 GB |
| 5.0M | 10.6 GB | 1.33 GB | 0.66 GB |
| 10.0M | 42.5 GB | 5.31 GB | 2.66 GB |

f = 10M is now a ~2.7 GB job. The binding constraint on the frontier is time,
not memory.

## Correction to an earlier snapshot

A mid-campaign snapshot of HANDOFF.md and PR #3 showed the 2.5M–3.2M chunk as
"running". That snapshot was stale: the campaign continued for a further ~13
hours and completed 2.5M–3.2M (12,443 candidates), 3.2M–4.0M (18,003) and
4.0M–4.3M (8,050), all with zero solutions. **The cleared frontier is
f = 4,300,000**, as recorded in the results table above.
