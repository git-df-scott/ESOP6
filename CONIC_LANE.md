# Conics on the Fermat sextic fourfold — exact reduction and status

Date: 2026-09-06.  All code is in `tools/conics/`, all logs in
`results/conics_2026_09_06/`.

## 0. Why conics, and why the earlier framing was wrong

The calibrated Hardy–Littlewood heuristic (`tools/heuristic_count.py`) gives
the expected number of generic solutions of `a^6+b^6+c^6+d^6+e^6=f^6` with
`f <= F` as about `5.0e-4 * log F` (singular series 0.0990, real density
5.07e-3). For k=5 the same computation predicts 2.7 solutions below 10^6 and
three are known, so the calibration is sound. For k=6 it predicts 0.007
solutions below 10^6 and 0.12 below 10^100. A counterexample, if one exists,
is therefore structured: it lies on a rational curve or an elliptic curve of
the fourfold

    X : x1^6 + x2^6 + x3^6 + x4^6 + x5^6 = x6^6   in P^5.

Two corrections to the previous documents:

* Signs are irrelevant: `(-a)^6 = a^6`. A rational curve on X only needs
  to miss the coordinate hyperplanes at the chosen parameter. There is no
  "positive chamber" condition on the curve.
* Slices such as `2X^6+2Y^6+Z^6=W^6` are surfaces of general type, where
  rational curves have negative expected dimension. On the Calabi–Yau
  fourfold X the expected dimension of rational curves is +1 in every
  degree. Curve searches belong on X, not on slices.

## 1. Proved: real rational curves have even degree and definite p6

Let `p = (p1,...,p6)` be a real rational curve on X in reduced form. If `p6`
has a real zero `t0`, then `sum p_i(t0)^6 = 0` forces every `p_i(t0) = 0`, so
`t - t0` divides all six coordinates, contradicting reducedness. Hence `p6`
has no real zero: its degree is even and it is a definite form. Conics are
the first case.

## 2. Proved: conics are trigonometric sixth-power identities

Normalise by a real Möbius map and scaling so that `p6 = g (s^2 + t^2)`.
Every binary quadratic is a combination of `A = s^2 - t^2`, `B = 2st`,
`C = s^2 + t^2`, which satisfy `A^2 + B^2 = C^2`. Writing
`p_i = x_i A + y_i B + z_i C` and `A = C cos(th)`, `B = C sin(th)`, the conic
condition is exactly

    sum_{i=1}^{5} ( z_i + x_i cos(th) + y_i sin(th) )^6  =  g^6   for all th.   (*)

With `zeta_i = x_i + i y_i` and `N_i = |zeta_i|^2`, the vanishing of Fourier
modes 1..6 of (*) is (conjugated so that zeta, not its conjugate, appears):

    M6:  sum zeta_i^6                                   = 0
    M5:  sum z_i zeta_i^5                               = 0
    M4:  sum (10 z_i^2 + N_i) zeta_i^4                  = 0
    M3:  sum z_i (8 z_i^2 + 3 N_i) zeta_i^3             = 0
    M2:  sum (16 z_i^4 + 16 z_i^2 N_i + N_i^2) zeta_i^2 = 0
    M1:  sum z_i (8 z_i^4 + 20 z_i^2 N_i + 5 N_i^2) zeta_i = 0
    g^6 = sum [ z_i^6 + (15/2) z_i^4 N_i + (45/8) z_i^2 N_i^2 + (5/16) N_i^3 ].

A conic is defined over Q iff, after a rational rotation and scaling, all
`x_i, y_i, z_i` are rational; then `g^6` must be a rational sixth power.
Rational rotations are multiplication of all `zeta_i` by an element of
`Q(i)^*`.

Two components appear.

### 2a. The design stratum `z_i = 0` (proved reduction)

M5, M3, M1 vanish identically. With `W_i = N_i^3` and `omega_i = zeta_i /
conj(zeta_i)` (a point on the unit circle), M6, M4, M2 read

    sum_i W_i omega_i^k = 0   for k = 1, 2, 3,

i.e. the atoms `omega_i` with weights `W_i` form a weighted 3-design on
the circle, and `g^6 = (5/16) sum N_i^3`. This is a Hilbert identity
`(A^2+B^2)^3 = sum c_i l_i(A,B)^6`; the classical four-term identity
`8A^6 + 8B^6 + (A+B)^6 + (A-B)^6 = 10 (A^2+B^2)^3` is the square
configuration and was recovered numerically (`g^6 = 10` exactly).

By the Jones–Njåstad–Thron / Szegő characterisation of positive quadrature
rules on the circle, five distinct nodes carry a nonzero weight vector iff
`e2(omega_1,...,omega_5) = 0`; the weights are then unique up to scale, and
they are all positive iff `|e1(omega)| < 1`. The positivity criterion was
checked on 38,301 random `e2 = 0` configurations with no exception
(`tools/conics/szego_check.py`).

A rational conic in this stratum therefore requires:

1. five rational points `omega_i = gamma_i / conj(gamma_i)` on the circle
   (`gamma_i` Gaussian integers) with `e2 = 0` and `|e1| < 1`;
2. the kernel weights `w_i` (automatically rational) with all ratios
   `w_i / |gamma_i|^6` in one class of `Q^* / Q^{*6}`;
3. `(5/16) sum N_i^3` a cube of a rational sum of two squares.

Exhaustive search over Gaussian integers of norm up to `B^2`
(`tools/conics/design_search.py`): see §4 for counts. Rational positive
designs exist: the first ones appear at norm 256, for example the nodes
`1, (20-99i)/101, (-15+8i)/17, (11+60i)/61, (-132-85i)/157`, with
`|e1|^2 = 2029309/16443709`. So the moduli surface of the design stratum has
rational points; what fails, in every case found, is the sixth-power class
condition on the weights (not even the weaker cube condition holds).

### 2b. The general component `z_i != 0` (structure observed, then proved exact)

Every non-degenerate general conic found numerically has, after rotation,
the reflection-symmetric shape

    l1 = z + r1 cos(th),   l5 = z - r1 cos(th),   l2 = r2 sin(th),
    l3 = r3 cos(th - alpha),  l4 = r3 cos(th + alpha).

For this shape (*) is equivalent, with `r1 = 1`, `c = cos(2 alpha)`,
`m = r3^6`, `Z = z^2`, `A = r2^6`, to the exact system

    F6:  A = 2 + 2 m (4c^3 - 3c)
    F4:  Z = -( 2 + m (c+1)(4c^2 - 2c - 1) ) / 10
    F2:  4 Z^2 + 4 Z = m c (c^2 - 1).

Eliminating Z gives the plane curve

    Gamma:  m^2 P(c)^2 - m ( 6 P(c) + 25 c (c^2 - 1) ) - 16 = 0,
            P(c) = (c+1)(4c^2 - 2c - 1),

verified against a numerically found conic to 25 digits
(`tools/conics/gamma_curve.py`). Real conics occupy the arc from the design
end `c = 0` (alpha = 45 degrees, `m = 2`, the square identity) to a boundary
end near `c = -0.31` where `m -> infinity`.

A rational conic of this shape needs a rational point of Gamma, hence the
discriminant a rational square. The discriminant factors:

    D(c) = 25 (c+1)^2 ( 137 c^4 - 186 c^3 + 21 c^2 + 28 c + 4 ),

so the condition is a rational point on the genus-one curve

    E :  y^2 = 137 c^4 - 186 c^3 + 21 c^2 + 28 c + 4,

with Jacobian `Y^2 = X^3 - 611307 X + 183887334` (invariants I = 22641,
J = -6810642). Additional conditions for the conic to be rational: `m` a
rational cube `n^3` with `n(1+c)/2` and `n(1-c)/2` rational squares, `A` a
rational sixth power, `Z` a positive rational sum of two squares, and
`g^6 Z^3` a rational sixth power.

PARI 2.15.4 (`results/conics_2026_09_06/elliptic_curve_E.gp` and its log):
the Jacobian has conductor 4070 = 2·5·11·37, torsion group Z/4, `ellrank`
returns lower bound 0 and upper bound 0, and the analytic rank is 0
(L(E,1) = 2.3676...). Hence E(Q) is finite of order 4. The quartic has no
rational points at infinity (137 is not a square), so its rational points
are exactly `(c, y) = (0, ±2), (1, ±2)`. At `c = 0`, `m = 2` is not a cube
(and `Z = 0`, the degenerate design end); at `c = 1`, `Z = -1 < 0`, not a
real conic. The point `c = -1` is excluded since `P(-1) = 0` makes Gamma
read `-16 = 0`.

**Theorem (conditional only on the reflection-symmetric shape).** The
general component of conics on the Fermat sextic fourfold contains no
conic defined over Q. Consequently a rational conic on X, if any, lies in
the design stratum of §2a.

## 3. Method validation

The Gauss–Newton solver (`tools/conics/conic_general.py`) finds
Ramanujan-type real conics on the Fermat cubic surface (41 distinct on 400
starts) and genuine real conics for (6,1,6) and (6,1,7). For (6,1,5) it finds
genuine real conics at about 5% of converged starts; all are design,
boundary, or reflection-symmetric general type
(`tools/conics/pattern_check.py`).

## 4. Search record

| search | domain | result |
|---|---|---|
| rational 5-node designs, Gaussian norm <= 64 | 60 nodes, 223 e2=0 configs | 0 with positive weights |
| rational 5-node designs, Gaussian norm <= 144 | 132 nodes, 491 configs | 0 with positive weights |
| rational 5-node designs, Gaussian norm <= 256 | 240 nodes, 919 configs | 4 with positive weights (2 up to conjugation); 0 sixth-class matches, 0 cube matches |
| rational 5-node designs, Gaussian norm <= 400 | 384 nodes, 1493 configs | 30 with positive weights; 0 sixth-class matches, 0 cube matches |
| rational points on E, height <= 1500 | | c = 0, 1 only |
| real conics on X, full 14-unknown system | 3000 real starts | 72 genuine, none rational |

## 5. What this does and does not establish

* Established exactly: the reduction of conics on X to (*), the
  design/general dichotomy, the design-stratum criterion, and the curve
  Gamma with its discriminant curve E.
* Established by exhaustive small-height search: no rational design conic
  at Gaussian norm up to 400 and no rational general conic of height up to
  1500 on E.
* Established with PARI: rank(E) = 0, E(Q) = Z/4, so no rational conic in
  the reflection-symmetric general family.
* Not established: that every general conic is reflection-symmetric (all
  37 numerical general examples are, and the dimension count agrees).
* Open: a rational positive 5-node design whose Szegő weights are sixth
  powers times `|gamma_i|^6` (equivalently a rational point on the
  sixth-power cover of the design surface); rational curves of
  degree 4 (trigonometric polynomials of degree 2 whose sixth powers sum to
  a constant); genus-one curves on X.

## 6. Reproduction

    pip install numpy scipy mpmath sympy
    cd tools/conics
    python3 conic_general.py 6 5 600 1        # real conics on X
    python3 gamma_curve.py                    # verify Gamma exactly
    python3 gamma_rational.py                 # rational points on E, height 300
    python3 ec_weierstrass.py                 # Weierstrass form, height 1500
    python3 design_search.py 12               # rational designs
    python3 szego_check.py                    # positivity criterion
