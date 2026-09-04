# Exact geometric strike on the repeated-coordinate surface

Target: `S: 2X^6+2Y^6+Z^6=W^6` over Q with all coordinates nonzero.
No rational surface point or rational curve on S was found. The symbolic
checks are in `tools/astra_geometry.py`. Quotient points below are always
distinguished from square lifts to S.

## 1. Geometry that constrains the attack

S is smooth: its four partial derivatives vanish simultaneously only at the
excluded zero vector. Adjunction gives `K_S=O_S(2)`, which is ample. Thus S
is of general type; the ambient six-variable fourfold has a different
canonical bundle and must not be conflated with this surface.

In particular, S cannot admit a genus-one fibration, even after resolving a
rational pencil. On a resolution `pi:S'->S`, a general fibre F satisfies
`F^2=0`, whereas

```text
K_S' . F = pi*(2H) . F + (effective exceptional divisor) . F > 0.
```

Adjunction then gives `2g(F)-2>0`. This does not exclude isolated elliptic
curves on S, or genus-one constructions on a quotient or on the full fourfold.

The inherited Bremner–Choudhry–Ulas control

```text
(1-t-t^2)^3 + (1+t-t^2)^3 = 2 - 2t^6
```

was checked by exact polynomial expansion. Their square/cube lifting method
motivates the quotient attack; the displayed identity does not have the
required positive `2,2,1,-1` sixth-power signature. Real nonzero sixth-power
coordinate rescaling cannot reverse a coefficient's sign. See their
[paper, Section 5 and concluding constructions](https://arxiv.org/pdf/1402.4583).

## 2. Low parameter degrees and boundary contact

A reduced map `P^1 -> S` is represented by four homogeneous forms of the
same degree d with no common zero. At a real zero of W the surface equation
forces X=Y=Z=0. Thus W has no real projective zero; its degree must be even.
This excludes d=1 and d=3 (indeed every odd d).

There are no conics on a smooth diagonal sextic over an algebraically closed
field of characteristic zero: Salberger's Theorem 9.4 gives degree at least
`(6+1)/3` for any curve other than the standard lines. A degree-two map has
image a line or conic; a real line was already excluded. Hence d=2 is also
excluded. This uses a published theorem, not a new classification result:
[Salberger, Section 9](https://londmathsoc.onlinelibrary.wiley.com/doi/full/10.1112/plms.12508).

The elementary pairwise-coprime Wronskian check is consistent with that
bound: fourth powers of the four quadratic forms would force degree 32 in
a Wronskian whose maximum degree is 30. We use the published theorem to
cover all shared-root and dependent cases, rather than assuming coprimality.

Now consider either obvious rational boundary point
`B_epsilon=[0:0:1:epsilon]`, `epsilon=+1 or -1`. If the map passes through
it at a parameter point, set `m=min(ord X,ord Y)>=1`. In

```text
(W-epsilon Z) * product of the other five linear factors of W^6-Z^6
    = 2(X^6+Y^6),
```

the second factor is nonzero at that point. Therefore
`ord(W-epsilon Z)>=6m`. A nonzero section of `O(d)` has at most d zeroes
counted with multiplicity, so d>=6. If `W-epsilon Z` vanishes identically,
the real equation forces X=Y=0 and the map is constant.

Equivalently, in the chart Z=1, W=1+h,

```text
2(X^6+Y^6) = 6h+15h^2+20h^3+15h^4+6h^5+h^6.
```

The tangent plane contains only the boundary point over R. Tangent lines,
conics, cubics, and quartics through it cannot produce the desired move into
the positive interior.

**Scope:** d=1,2,3 are excluded globally for rational parametrizations.
For d=4, only the boundary-contact locus is excluded here. General rational
quartics avoiding those points were not classified. Other boundary rational
points were not exhaustively classified either.

## 3. The diagonal restriction is closed

Suppose `4X^6+Z^6=W^6` has a positive rational solution. Clear denominators
and divide the gcd. The usual parity argument makes W,Z odd. Any common
prime factor of W,Z also divides X, so gcd(W,Z)=1. Thus

```text
A=(W^3+Z^3)/2, B=(W^3-Z^3)/2
```

are coprime positive integers and `AB=X^6`. Consequently A=u^6, B=v^6.
But then `(u^2)^3+(v^2)^3=W^3`, contradicting Fermat's theorem for cubes.
Signs do not evade this argument. This closes `X=+Y` and `X=-Y`, but not
near-diagonal points with unequal absolute values.

Fixing any other constant ratio `X/Y=r` gives a smooth diagonal plane sextic
of genus 10. Its lower-genus quotients still require exact lifting; treating
one such quotient as an elliptic curve on S would lose conditions.

## 4. An explicit conic pencil on the cubic quotient

Set

```text
U=X^2, V=Y^2, R=Z^2, T=W^2,
A=U+V, B=U-V, C=T-R, D=T+R.
```

The cubic quotient is `2U^3+2V^3+R^3=T^3`. Its rational line
`U+V=0, T-R=0` gives the plane pencil `C=lambda A`. Substitution gives

\[
(2-\lambda^3)A^2+6B^2=3\lambda D^2. \tag{1}
\]

The exact identity checked by the script is

```text
4*(2U^3+2V^3+R^3-T^3)
  = A*((2-lambda^3)A^2+6B^2-3lambda D^2).
```

At a positive lift, `A>0`, `|B|<A`, `D>lambda A`, and
`0<lambda<cuberoot(2)`. These inequalities follow directly from positive
U,V,R,T and their cubic equation.

To lift a rational projective quotient point, multiply it to a **positive
primitive integer vector** `(U,V,R,T)`. All four entries must be perfect
squares. This criterion is necessary and sufficient: for any prime, the
common scalar's valuation is even because some primitive coordinate is a
unit, and all unit square classes must agree. An actual rational square
lift therefore has square primitive coordinates. Over R there is no further
sign difficulty once all four coordinates are nonzero.

### Why the generic lift is genus 9

The conic determinant is `18lambda(lambda^3-2)`. For rational
`lambda != 0,2`, the conic is smooth and meets each of U,V,R,T=0 in two
distinct points, with no overlaps. The tangency factors are
`8-lambda^3` and `2-4lambda^3`; overlap occurs only on a singular conic.

The projective square map is a degree-eight `(Z/2)^3` cover. There are eight
branch points of inertia order two. Riemann–Hurwitz gives

```text
2g-2 = 8*(-2) + 8*4 = 16, hence g=9.
```

In the original quadric section, four complex base lines are also present;
the genus-nine computation concerns the residual degree-eight curve, not
the entire degree-twelve intersection.

At the only nonzero rational tangency value `lambda=2`, (1) becomes

```text
(W^2+Z^2)^2 + 4X^2Y^2 = 0.
```

It has no nonzero real lift. Over Q(i), the residual member splits using
`W^2+Z^2=+/-2iXY`. Each component is the intersection with
`W^2-Z^2=2(X^2+Y^2)` of a second quadric. The determinant of its quadric
pencil is `(4alpha^2+beta^2)(beta^2-alpha^2)`, with distinct roots, so these
are smooth genus-one curves over the extension. Their lack of real points
closes this degeneration as a positive-point attack.

Thus the natural rational conic pencil supplies no genus-one fibre with a
positive rational lifting route. This does not exclude other isolated curves.

## 5. Coupled valuation restrictions on the pencil parameter

In a primitive solution of the repeated surface, reduction modulo 8,9,7
forces `42|X,Y` and `gcd(WZ,42)=1`. Put
`r_p=min(v_p(X),v_p(Y))>=1`, and write

```text
D0=W^2-Z^2, Q0=W^4+W^2 Z^2+Z^4,
D0*Q0=2(X^6+Y^6), lambda=D0/(X^2+Y^2).
```

At 2, let epsilon be 1 if both reduced bases are odd and 0 otherwise.
Then `v2(X^6+Y^6)=6r2+epsilon`,
`v2(X^2+Y^2)=2r2+epsilon`, and Q0 is odd. Hence

```text
v2(lambda)=4r2+1 >= 5.
```

At 3, `v3(X^6+Y^6)=6r3`, `v3(X^2+Y^2)=2r3`, and `v3(Q0)=1`.
The last equality follows on writing `W^2=Z^2+3k`. Therefore

```text
v3(lambda)=4r3-1 >= 3.
```

At 7, -1 is not a square, so the corresponding sum valuations are 6r7
and 2r7. If `W^2=Z^2 mod 7`, Q0 is a unit and
`v7(lambda)=4r7`. Otherwise D0 is a unit and
`v7(lambda)=-2r7`; then Q0 carries the required valuation. Both branches
must be retained. Discarding the negative-valuation branch would be wrong.

In reduced form the numerator of lambda is therefore divisible by 864.
The easily parametrized `lambda=1/2` member fails these necessary conditions
and is closed at every rational parameter. Its completed small-height test
is recorded, but should never be enlarged.

The next experiments selected
`lambda=864/(49d)` for `d=17,19,23,25,29,31,37,41`.
Four explicit conic base points were returned and verified. For each,
finite local square-lift tests had surviving residue classes; this does
**not** assert existence of a point over every local field. The four ensuing
bounded rational searches are recorded in the ledger.

## 6. Exact factorization attack on the same surface

Set X=42x, Y=42y. Then

```text
N=(W^6-Z^6)/(2*42^6)=x^6+y^6.
```

For any prime `p=3 mod 4` dividing N, reduction of
`(x^3)^2+(y^3)^2` forces p to divide both x and y. Repeating proves
`v_p(N)=0 mod 6`. This is a two-term condition; it is not a valid replacement
for the inherited valuation law for four terms.

Let `h=x^2+y^2`. The identities

\[
N=h(h^2-3x^2y^2),\qquad
(x^2-y^2)^2=\frac{4N/h-h^2}{3}
\]

give a complete oracle:

1. Factor N using `(W-Z)(W+Z)(W^2+WZ+Z^2)(W^2-WZ+Z^2)` and remove
   exactly `2^7*3^6*7^6`.
2. Enumerate all divisors h with `N<h^3<=4N`.
3. Demand that `(4N/h-h^2)/3` be an integer square `d^2`.
4. Demand that `(h+d)/2` and `(h-d)/2` be positive integer squares.
5. Recompute the original sixth-power equation before reporting a hit.

The inequalities are equivalent to a positive product `x^2*y^2` and its
nonnegative discriminant. Thus the enumeration is complete for each target.
Each factor prime is proved using recursive Lucas certificates, and the
retained search is independently replayable without a factoring package.

Work per target is factorization plus divisor enumeration, not a loop through
every possible x or a large pair table. There is no claimed general
polynomial factorization bound and no exponent reduction for full ESOP6.

## 7. Minimal boundary-contact construction and one exact dead subfamily

Degree six is the first parameter degree that can leave B+ with order-six
contact. Normalize a rational preimage to t=0 and choose

```text
X=t*A5(t), Y=t*B5(t), Z(0)=A5(0)=B5(0)=1,
W-Z=c*t^6.
```

The leading surface coefficient forces `c=2/3`. Dividing the identity by
`t^6` removes the known contact before any coefficient solving. Independent
A5,B5 of degree at most five and Z of degree at most six are the next
specified construction. Positivity near t=0 would be automatic for a
successful exact identity.

A further two-boundary-point parity restriction was tested exactly:

```text
W=1+t^6/3, Z=1-t^6/3,
X=t*A4(t), Y=t*A4(-t), A4(0)=1.
```

It requires

```text
A4(t)^6+A4(-t)^6 = 2+(20/27)t^12+(2/81)t^24.
```

Its leading coefficient would satisfy `a4^6=1/81`, impossible over Q
because the 3-adic valuation is -4. Lower degree also fails to supply the
nonzero t^24 coefficient. This closes exactly that normalized symmetric
subfamily, and explains why the next construction must retain independent
A5,B5 and a free Z.
