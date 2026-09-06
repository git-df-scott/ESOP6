# Astra prompt — rational and elliptic curves on the sextic fourfold (compute split)

User handoff received September 5, 2026. The remote did not contain this file
at checkout. This preserves the supplied construction assignment; corrections
and execution results are in `ASTRA_CURVES_REPORT.md`. The credential and shell
wrapper are intentionally not part of the mathematical handoff.

## Target

Repository: `git-df-scott/ESOP6`.

    a^6+b^6+c^6+d^6+e^6=f^6, all positive integers.

Read README.md, STATUS.md, ASTRA_HANDOFF.md, and this file first.

## Supplied background and constraints

The user reports the concentrated sweep at f=4,310,000; do not extend it.
The supplied heuristic is c*log(F), with c of order 1e-5. Construction on the
full sextic fourfold is the assigned lane. Do not use the repeated-coordinate
surface slices. Do not rerun integer sweeps or touch the frontier.

X is the smooth sextic fourfold in P5. Rational curves have expected dimension
one, and genus-one curves expected dimension zero. The supplied prompt asserted
that positivity makes every coordinate definite; the execution report corrects
this to the necessary condition on the sixth coordinate. Expected dimensions
are not assertions of existence or definition over Q.

The other session owns full-fourfold numerical CONICS, including the
normalization a1=1, p6=g(1+t^2), b5=0, and symmetric ansatz searches. Do not
duplicate those numerical searches.

## Assigned lanes, in priority order

1. Genus-one curves of low degree.
   - Plane cubics: restrict the sextic to a plane and impose cubic times cubic.
     Set up the polynomial system, solve numerically if appropriate, refine,
     and test rationality.
   - Elliptic quartics: intersections of two quadrics in a P3 contained in X.
     Use the corresponding polynomial incidence system.
   - For a genus-one curve over Q, compute its Jacobian, rank and generators
     where possible, and search rational points with all coordinates nonzero.
     A Jacobian point must be lifted to the actual curve before counting it.
2. Degree-four rational curves under symmetry.
   - (12)(34) with t -> -t, p5 and p6 even.
   - Add (13)(24) with s <-> t.
   - Normalize the overall scalar, account for PGL2, count unknowns and
     equations, run bounded solves, refine, and test rationality/definiteness.
3. Exact conic Fano-scheme cross-check.
   - Restrict to a plane and impose conic times quartic.
   - Use symmetry before large elimination.
   - Investigate components over Q and their real/rational loci; distinguish
     the general-sextic Fano curve from the special Fermat scheme.

## Hard rules

Never announce a counterexample from floating point. A rational candidate
curve must satisfy the identity exactly over Q. A candidate sextuple must pass
both repository verifiers and be committed with their outputs before an
announcement:

    python3 tools/verify_esop6.py a b c d e f
    node tools/verify_esop6.mjs a b c d e f

For every solved system record unknowns, equations, Jacobian ranks, complex
and real counts, definite counts and rational counts. When a count is not
complete or certified, label it as such. Negative bounded results are a
permitted deliverable if no curve appears.

Commit to an `astra/curves-<lane>` branch and write retained results under
`results/astra_curves_<date>/`.

## Other-session numerical log supplied by the user

B: 1400 reported distinct solutions; real 0; real/all-definite 0;
unknowns 8; rank histogram {7:395,3:373,4:346,5:70,2:154,6:62}.

D: 1315 reported distinct solutions; real 0; real/all-definite 0;
unknowns 14; rank histogram
{3:828,6:205,9:3,5:90,13:43,8:20,7:88,12:7,4:17,11:9,10:5}.

These logs were supplied, not independently reproduced in this run. The
exact real D control now appears in `ASTRA_CONIC_CONTROL.md`.
