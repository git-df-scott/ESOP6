# Engine audit

## Scope and environment

Audited commits: ESOP6 `26b4f83`, historical branch through `33690dc`, and
historical production blob `0a7d9e6`.  Replay host: x86-64 Intel Xeon Platinum
8573C, 9 logical CPUs, GCC with OpenMP, 21 GiB RAM, no swap.

## Control results

| Control | Expected | Result |
|---|---|---|
| `make validate` | Lander–Parkin at 144; root audit all OK | PASS, 0.242 s |
| `make control` | 124 candidates, 0 solutions | PASS, 6.498 s |
| `make equiv` | 4,139,088 effective j2 nodes for NB=1 and NB=8 | PASS |
| `make frontier-audit` | all eight counts; total 55,684 | PASS, 0.553 s |
| `make differential` | identical candidate sets across random bands and NB 1,2,3,4,7,8 | PASS, 217.5 s |
| `make prototypes` | exact routed join and mmap Bloom have no misses | PASS |
| UBSan | caseA3 700k–701k, NB=3 | PASS |
| ASan/UBSan | tiny-fmax self-test, leak detection disabled due ptrace host | PASS |

Randomized differential bands used fixed seed `0xE50F6` and contained 35, 27,
and 31 candidates respectively.  The comparison is of actual `(f,t)` sets,
not counts.

## Bloom sizing

For \(B=\lceil f_{max}/42\rceil\) and \(P=B(B+1)/2\) pairs, exact 12-bpp
storage is `ceil(12P/512)*64` bytes.

| fmax | pairs | fixed caseA2 | old rounded caseA2 | old effective bpp | caseA3 ×8 | caseA3 ×16 |
|---:|---:|---:|---:|---:|---:|---:|
| 4.3M | 5,240,985,771 | 7.861 GB | 8.590 GB | 13.112 | 0.983 GB | 0.491 GB |
| 5.0M | 7,086,272,676 | 10.629 GB | 17.180 GB | 19.395 | 1.329 GB | 0.664 GB |
| 5.5M | 8,574,409,581 | 12.862 GB | 17.180 GB | 16.029 | 1.608 GB | 0.804 GB |
| 10M | 28,344,971,656 | 42.517 GB | 68.719 GB | 19.395 | 5.315 GB | 2.657 GB |

The sizing bug cannot create a false negative: it allocated more bits than
requested.  It distorted memory planning and reduced false positives.

## caseA3 timing

700k–730k, 16 bpp, 8 threads:

| Buckets | Filter | Wall | Relative to NB=1 | Reached | Skipped | Bloom-evaluated |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.30 GB | 6.18 s | 1.00× | 4,139,088 | 0 | 4,139,088 |
| 2 | 0.15 GB | 7.73 s* | 1.25×* | 8,278,176 | 4,139,088 | 4,139,088 |
| 4 | 0.08 GB | 7.74 s* | 1.25×* | 16,556,352 | 12,417,264 | 4,139,088 |
| 8 | 0.04 GB | 10.57 s | 1.71× | 33,112,704 | 28,973,616 | 4,139,088 |
| 16 | 0.02 GB | 17.61 s* | 2.85×* | 66,225,408 | 62,086,320 | 4,139,088 |

`*` single timing; NB=1 and NB=8 are two-run means after final compilation.
The penalty is sublinear in NB on this control because wrong-bucket leaves
exit before Bloom and exact verification.  Larger ranges still need
calibration; it is unsafe to extrapolate a constant 1.71× to 10M.

## Correctness findings and fixes

1. Added missing `<string.h>` in caseA3.
2. Self-tests formerly indexed `P6[200]` even when `fmax/42 < 200`.  Both
   engines now use a dedicated 0..200 power table during self-test.
3. A reduced decomposition with fewer than four terms was printed as
   `SOLUTION`.  That corresponds to padding the original equation with zeros,
   forbidden by the positive-integer target.  Such diagnostics are now
   `DEGENERATE`; only four reduced parts increment `found`.
4. Misleading-indentation warnings were removed.
5. The exact CE verifier is `tools/verify_ce.py`; it rejects nonpositive input
   and uses Python big integers only.

These fixes do not invalidate any historical zero-result: no historical run
printed a shorter decomposition or a solution.

## Integer and overflow audit

- `f <= 1e8` implies reduced base `Bmax <= 2,380,953`; its sixth power and
  `(f^6-t^6)/42^6` fit unsigned 128-bit, though not signed 128-bit near the
  ceiling.
- Pair counts, byte counts, and root products fit unsigned 64/128-bit in the
  enforced range.
- `iroot6` starts from floating approximation but corrects both directions
  with exact sixth powers; floating point affects speed, not the result.
- OpenMP Bloom writes use atomic OR.  All queries occur after the build-loop
  barrier.
- Instrumentation counters are signed 64-bit and are a diagnostic ceiling for
  giant/high-NB jobs.  They are not part of the search decision.

## Single-pass prototypes

| Prototype | Control | Result | Scale verdict |
|---|---|---|---|
| exact external routed join | N=2000, NB=8, 2,001,000 pairs, 20,000 queries | one pair pass, one query pass, 32.016 MB pair disk, peak 250,862 entries, 0 mismatches | exact but disk is 16 bytes/pair plus query records |
| mmap Bloom | N=5000, 12 bpp, 12,502,500 pairs, 20,000 queries | 18.754 MB file, 45 false positives, 0 false negatives | compact, but random faults can collapse beyond RAM |

The best practical CPU choice on this host remains caseA3 with the smallest
bucket count that fits.  The best way to truly remove repeat traversal is a
single-pass routed/offline join on storage or GPU, but only if bandwidth and
query buffering are engineered explicitly.
