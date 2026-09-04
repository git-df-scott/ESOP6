# Formal branch at [0:0:1:1]

**A branch with rational formal coefficients exists. An exact rational-function
curve was not found.** These are different statements. Here “FORMAL RATIONAL
BRANCH: YES” means a branch in Q[[t]], not a parametrization in Q(t).

## Exact recursive branch

For any A,B in Q[[t]] with constant coefficients 1, set

```text
D=(A^6+B^6)/2,
F(M,t)=M^5+(10/27)t^12 M^3+(1/81)t^24 M-D.
```

At t=0,M=1, `F=0` and `dF/dM=5`. Therefore there is a unique formal
M=1+O(t), with each successive coefficient given by equation (2) in
[BOUNDARY_CONTACT_6.md](BOUNDARY_CONTACT_6.md). Set

```text
X=tA, Y=tB, Z=M-t^6/3, W=M+t^6/3.
```

This is an exact identity of formal series. There are two arbitrary input
series A,B; M is determined, so this is an infinite-dimensional formal
choice, not a finite-dimensional rational-curve family. For polynomial or
convergent rational inputs, the implicit-function theorem gives a convergent
real analytic branch near zero. Arbitrary formal inputs need not converge.

Up to order eleven, M is simply the formal fifth root of D. The first
forced contact correction enters at order twelve. Allowing the missing
coefficients m7,m8,... always repairs a finite jet; requiring all of them
to vanish is the global polynomial restriction that fails in the proof.

## Two exact jets, and a rational reconstruction that fails

For A=B=1, put s=t^12. Exact arithmetic through t^72 gives

\[
M=1-\frac2{27}t^{12}+\frac{11}{3645}t^{24}
 -\frac{341}{13286025}t^{48}+\frac{644}{358722675}t^{60}
 +\frac{3397}{16142520375}t^{72}+O(t^{84}).
\]

The t^36 coefficient is zero. It is a gap, not truncation. Only primes
3 and 5 occur in these denominators; the recurrence proves that all
coefficients of this particular branch lie in Z[1/15]. This does not imply
a rational generating function.

The [1/1] Padé expression in s that matches the first two nonconstant
coefficients is

```text
M_Pade=(1-s/30)/(1+11s/270).
```

Its exact substitution into `M^5+(10/27)s M^3+s^2 M/81-1` is

```text
-s^3(14641s^4+1292280s^3-246600s^2-951588000s+7938810000)
 / (9(11s+270)^5).
```

The first error is `-121 t^36/196830`. This expression is not a curve on
the surface. No approximate contact was promoted to an identity.

As an asymmetric control, A=1+t and B=1-t give

```text
M=1+3t^2-15t^4+(631/5)t^6-(6177/5)t^8+(65517/5)t^10+...
```

The full exact jet through t^30 is in `symbolic_checks.json`. Its polynomial
degree-six truncation first fails at order eight. Both retained jets were
substituted using Fraction arithmetic through every recorded coefficient;
the regression test also asserts that no floating-point coefficient enters.

## Why a denominator only in Z cannot help

With A and B polynomial, the equation for M is monic over Q[t]. More
elementarily, if M=P/Q in lowest terms, multiply the equation by Q^5:

```text
P^5+(10/27)t^12 P^3 Q^2+(1/81)t^24 P Q^4-D Q^5=0.
```

It follows that Q divides P^5, so Q is constant. Any rational-function
midpoint would already be polynomial. When deg A,deg B<=5, the degree
argument forces deg M<=6, and the full obstruction applies. Thus increasing
the Padé denominator of Z alone cannot evade this failure.

For A=B=1 the impossibility is even visible directly: a polynomial M of
degree below six would leave the unique term of degree `24+deg M`; one
of degree above six would leave the unique leading term M^5. At degree
six its leading coefficient r would require
`r(r^4+(10/27)r^2+1/81)=0`, which has no nonzero real root.

There is also a useful explanation of the algebraic, non-rational branch.
For A=B=1 the equation is quadratic in s=t^12. With

```text
v=(M/3)(s+15M^2),
```

it becomes `v^2=16M^6+9M`, with base point `(M,v)=(1,5)` and
`s=-15M^2+3v/M`. The sextic is squarefree, so this intermediate
hyperelliptic curve has genus two. The additional condition s=t^12 remains.
This is not an elliptic residual or a positive rational-point construction.

## Local geometry and positivity

In the affine chart Z=1 the surface is the formal/analytic graph

```text
W=(1+2x^6+2y^6)^(1/6)
 =1+(x^6+y^6)/3-5(x^6+y^6)^2/18+...
```

The tangent direction selected by the normalized construction is (1,1,0)
in the (x,y,W-1) coordinates. The tangent plane is W=1, and its intersection
with the real surface has only the boundary point. The path therefore has
order-six contact with that plane while the surface itself is smooth.

This graph is not a rational formula in x,y: the polynomial
`1+2x^6+2y^6` is squarefree as a polynomial in x over Q(y), so it cannot be
a sixth power in Q(x,y). Smoothness guarantees the formal branch, not
rationality of this sixth root. No rational local parametrization of the
surface has been constructed here.

For the analytic A=B=1 branch, `0<t<1/2` implies `1/2<M<1`: substitute
M=1/2 and M=1 into its strictly increasing positive quintic. Hence
`Z=M-t^6/3>1/2-1/192>0`, W>0 and X=Y=t>0. These are positive **real**
points; their existence does not supply rational coordinates. The inherited
diagonal obstruction excludes rational specializations with t nonzero.

For a future exact rational curve, positivity can be certified without a
height search. Normalize all numerator and denominator constant terms to
1. If C is the maximum sum of absolute values of their nonconstant
coefficients, choose `epsilon=min(1,1/(2C))` (epsilon=1 if C=0). Every
such polynomial is positive for `0<t<epsilon`. Choose any rational t in
that interval, clear denominators, and invoke both standalone sextuple
verifiers. No such exact curve was available in this strike.
