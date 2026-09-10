# ESOP6 status — 2026-09-10

## 2026-09-10 current continuation

No counterexample. The normalized fibre pilot has **573 certified exclusions
and 54 unresolved fibres**; exclusions are global on their specified fibres,
not an unrestricted integer-height frontier. See
[REPORT.md](research/twisted_fibres_2026_09_10/REPORT.md),
[the complete ledger](research/twisted_fibres_2026_09_10/FIBRE_LEDGER.md), and
[the next arithmetic task](research/twisted_fibres_2026_09_10/QC_HANDOFF.md).
The final shortlist has **38** genus-two rank-(1,1) fibres (39 was an interim
count before the last gate closed one). All new jobs completed. Historical
remote jobs were not controlled or restarted. Read the explicit corrections
before importing the night report's conjectural or numerical closure claims.

## Latest: second strike after 80a2d9f

No ESOP6 solution or positive rational surface point was found. The complete
normalized degree-six boundary-contact polynomial ansatz is impossible over
Q, including asymmetric A,B and arbitrary rational coefficient denominators.
A 3-adic reduction and independent Python/C++ coefficient certificates give
the obstruction modulo 729. The formal branch in Q[[t]] exists, but is not
an exact rational-function curve. No elliptic residual was obtained.

All inherited gates passed; the new `make boundary-checks` gate verifies the
symbolic calculations and both finite certificates. No integer search was
repeated or extended. The historical frontier qualification is unchanged.

Read [ASTRA_SECOND_STRIKE.md](ASTRA_SECOND_STRIKE.md) and
[BOUNDARY_CONTACT_6.md](BOUNDARY_CONTACT_6.md). The single next construction
is the degree-eight extension with a positive quadratic denominator in
[RATIONAL_CURVE_ATTEMPT.md](RATIONAL_CURVE_ATTEMPT.md). This closes only the
specified equal-slope degree-six family, not all rational boundary curves.

## Latest: direct solution strike after 485390a

No solution, positive rational surface point, or rational curve on the
surface was found. No elliptic curve with a rigorous positive lifting route
was obtained. A bounded targeted integer search was completed.

- The cubic quotient conic pencil lifts generically to genus 9; its rational
  exceptional member at lambda=2 has no nonzero real lift.
- Degrees 1,2,3 cannot furnish a rational parametrization. General degree 4
  remains open here; degree <=5 is excluded through the two obvious boundary
  points by order-six contact.
- New coupled constraints on lambda: v2=4r2+1, v3=4r3-1, and
  v7=4r7 or -2r7, with each r>=1.
- Four valuation-selected conics: 19,465,408 parameter occurrences,
  7,942,986 positive quotient images, no square lift in the specified boxes.
- A new table-free factor/divisor oracle tested 10,382 sparse surface targets.
  All targets and 4,090 prime certificates were independently replayed with
  Python's standard library. This is not an unrestricted ESOP6 height bound.
- Every requested inherited control gate passed; two new standalone integer
  verifiers also agree on 18 controls.

See [ASTRA_DIRECT_STRIKE.md](ASTRA_DIRECT_STRIKE.md) for exact domains,
evidence, and the single next direct construction: independent degree-six
forms with order-six contact at `[0:0:1:1]`. The historical 4.3M frontier
and its evidence grade are unchanged.

## September 4 preparation record

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
