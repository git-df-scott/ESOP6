# ESOP6 status — 2026-09-04

## Verdict: YELLOW

Only bounded search and mathematically specified structural lanes are justified.
The RAM wall is removed and the leaf sieve is stronger, but this audit did not
prove an exponent reduction below the concentrated engine's effective
\(F^4\) traversal.  A blind run to 10M would spend time, not resolve the real
bottleneck.

## Canonical result

- Target: positive integers \(a,b,c,d,e,f\) with
  \(a^6+b^6+c^6+d^6+e^6=f^6\).
- Concentrated/Meyrignac-class-1 frontier: **historically searched through
  \(f=4,300,000\)**, with 55,684 candidate \((f,t)\) pairs and no reported
  solution.
- This audit independently regenerated every candidate count and two set
  digests, exactly totaling **55,684**.  It did not rerun the expensive
  decompositions.  The historical repository retained a curated transcript,
  not original log files or checksums, so the zero-result is
  **historically attested, control-reproduced, but not end-to-end independently
  reproduced**.
- Cheap controls pass: Lander–Parkin, 700k–730k (124 candidates, 0 solutions),
  candidate-root audit, caseA3 4,139,088-node equivalence, randomized
  candidate-set differentials, and exact oracle prototypes.
- No counterexample was found.

## New results

1. The exact Bloom-sizing correction is confirmed: at 5M and 12 bpp,
   10.629409 GB replaces the old 17.179869 GB allocation (19.395 effective
   bpp).
2. caseA3 is complete.  At NB=8 it uses 1/8 the Bloom memory and cost 10.57 s
   versus 6.18 s for NB=1 on the 700k–730k control: **1.71× measured**, not 8×.
3. Incrementally tracked mod 19 and mod 31 filters reject 3,590,968 of
   4,139,088 leaf queries (86.76%) before Bloom lookup.  They are available as
   `--extra-sieve`.
4. A proved valuation filter for a sum of at most four sixth powers is
   available as `--valuation-prune`:
   \(v_2\bmod6\in\{0,1,2\}\), \(v_3\bmod6\in\{0,1\}\), and
   \(v_7\equiv0\pmod6\).
5. Two single-pass low-memory prototypes pass exact controls: external exact
   bucket routing and a disk-backed mmap Bloom.  Both exchange repeated CPU
   work for \(O(F^2)\) disk traffic; neither dominates caseA3 on this host.
6. The production reporter now labels fewer-than-four reduced parts
   `DEGENERATE`, not `SOLUTION`; only a four-part reduced decomposition can
   increment the CE count.

## Immediate recommendation

Do not run 4.3M→5.5M in this session.  Start Astra with the full regression
gate in [ASTRA_HANDOFF.md](ASTRA_HANDOFF.md), then attack the explicit
positive rational-point problem in [GEOMETRIC_ATTACK.md](GEOMETRIC_ATTACK.md)
or a class-complete fused GPU implementation.  Search farther only as a
bounded calibration after those gates pass.
