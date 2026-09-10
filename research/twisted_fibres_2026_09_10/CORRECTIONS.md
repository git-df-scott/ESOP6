# Corrections governing continuation

These corrections supersede the corresponding broad statements in the
night report. They preserve historical logs and do not restart completed
numerical batches.

1. **LPS is a conjecture in the cited argument.** Its proposed exclusion of
   nontrivial rational (6,1,4) boundary points cannot be used as a proved
   premise. A rational linear factor gives a rational boundary contact;
   an irreducible factor of degree k gives a degree-k contact, not necessarily
   a rational one. Reducibility does not imply a rational root.
2. **Degree six is a proved lower bound in the specified boundary-contact
   setting**, not a universal irreducibility or curve-existence theorem.
   The restricted real odd-degree/contact arguments remain separate valid
   arguments when their basepoint and degree hypotheses are satisfied.
3. **Random-start solvers are not certified exhaustive eliminations.**
   Running every seed in a finite list is not complete algebraic curve
   enumeration. No-rational-hit logs keep that exact evidential scope.
4. **Some end-of-night action items are stale.** The morning correction
   reopens the cubic-cover V3' lane. Its degree-six factor does not have
   a forced real root. This continuation does not reinstate the retracted
   closure or rerun its recorded lotteries.
5. **A partial floating-point seed scan is not an exact completeness
   certificate.** The retained N500 fast generator uses floating arithmetic
   before exact validation. Its points can be checked exactly, but its
   lack of missed candidates requires a separate error/completeness proof.
   The retained N500 list contains **31 exact equation failures** among
   2,924 entries, from comparing against `int(float(x6**6))`. The verified
   replacement has 2,893 distinct valid seeds. `seed_validation.json`
   retains each bad entry and its nonzero integer residual. The original
   lists remain unchanged; the generator's final check and the consumer's
   initial guard now use the actual exact equation. The N40 and N70 lists
   also contain duplicate representations: respectively 19 and 68 unique
   `(sorted fourtuple,x6)` keys, rather than 50 and 187 distinct points.
6. **Local insolvability cannot reject the twisted fibres.** The smooth
   rational boundary has nearby all-nonzero points at every local place.
   A missing unit-denominator chart is not a missing Q_p point.
7. **The simultaneous square cover is genus ten.** Genus four belongs to
   each intermediate double cover. Ordinary elliptic two-descent does not
   directly solve these branched square-class conditions.

No positive ESOP6 counterexample is claimed. Numerical curves, exact
algebraic seeds, rational quotient points and finite-height misses retain
their distinct meanings. See ARITHMETIC.md, INDEPENDENT_AUDIT.md, and
LITERATURE.md for the supporting arguments and source-access limits.
