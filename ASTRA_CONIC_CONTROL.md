# Exact real conic control for the other session

This is a **number-field control, not a Q-rational curve or ESOP6 candidate**.
It is supplied for the existing full-conic D solver; no duplicate numerical
conic campaign was run.

Put alpha=sqrt(2)>0, delta=11^(1/6)>0, g=delta/alpha. In the D normalization
p1(0)=1, p6=g(1+t^2), b5=0, take

    p1 = 1-t^2
    p2 = 2t
    p3 = (1+2t-t^2)/alpha
    p4 = (1-2t-t^2)/alpha
    p5 = (1+t^2)/alpha
    p6 = g*(1+t^2).

The homogenizations have no common zero, span all binary quadratics, and
parametrize a genuine conic. At t=1/4 all six coordinates are positive.
The sixth coordinate is definite; several other coordinates have real zeros.

Exact verification needs only

    8U^6+8V^6+(U+V)^6+(U-V)^6 = 10(U^2+V^2)^3,
    U=s^2-t^2, V=2st, U^2+V^2=(s^2+t^2)^2.

Adding (s^2+t^2)^6 on the left gives 11(s^2+t^2)^6. Divide all
coordinates by alpha to obtain the displayed normalization. The identity
is checked symbolically modulo alpha^2-2 and delta^6-11.

The original 13-by-18 coefficient Jacobian has exact rank **11**. The
13-by-14 D-normalized coefficient Jacobian also has exact rank **11**.
Thus rank deficiency by itself is not a certificate that a point represents
a constant map or common-factor artifact.

This conic cannot provide a rational real point: p6/p5=delta is irrational
and p5 never vanishes on real P1. It must never be submitted as a rational
candidate. Nevertheless, a numerical assertion that D has *no real solutions*
is false; a bounded run that *found none* is compatible with this control.

Reproduce:

    python3 tools/astra_curves/exact.py

See `results/astra_curves_2026_09_05/exact_results.json` for exact ranks,
field relations, the plane factorization, and the coordinate arrays.
