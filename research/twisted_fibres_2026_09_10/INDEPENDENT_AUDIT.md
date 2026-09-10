# Independent arithmetic audit of the twisted fibres

Date: 2026-09-10. Scope: exact mathematical audit; this file does not report
an executed search or a counterexample. The existing integer verifiers remain
the acceptance gate.

## 1. Objects and exact equivalences

Fix a positive rational number `c`. For the intended application initially
`c=a1^6+a2^6+a3^6+a4^6`, with all four `ai` positive integers. Define

```
C_c : A^6 + c S^6 = B^6,
E_c : X^3 + c Y^3 = Z^3.
```

These are smooth projective curves over Q. The morphism

```
C_c -> E_c,   [A:S:B] |-> [A^2:S^2:B^2]
```

has degree four. Its image on rational points consists exactly of points
having `X/Z` and `Y/Z` rational squares when `Z != 0`. To obtain an actual
positive ESOP6 solution, both squares must be nonzero. Then `X/Z>0`,
`Y/Z>0`, and the equation implies `X/Z<1`. Taking positive square roots
and clearing denominators produces

```
(a1 S, a2 S, a3 S, a4 S, A; B).
```

Signs can be chosen positive because all exponents are even. The square test
must be exact in Q, followed by denominator clearing and both independent
integer verifiers. The projective point `[1:0:1]` is excluded: it yields
four zero entries.

On `Y != 0`, put

```
u = Z/Y,   v = X/Y.
```

Then `u^3-v^3=c`. Conversely any rational pair with this cube difference
gives `[X:Y:Z]=[v:1:u]`. This is an exact equivalence, with no integrality
assumption. When `c != 0`, `u-v != 0` automatically. The additional lift
criterion is precisely that **both `u` and `v` are nonzero rational
squares**, because `X/Z=v/u` and `Y/Z=1/u`.

Rational cube differences without the square tests are not sextic points.
The cases `u=0` or `v=0` are allowed in the cube-difference equivalence but
are excluded from the positive sextic problem. There is exactly one rational
point on `E_c` with `Y=0`, namely `[1:0:1]`.

## 2. Mordell model, with inverse checked algebraically

Use the Weierstrass curve

```
M_c : eta^2 = xi^3 - 432 c^2.
```

The mutually inverse formulas away from the point at infinity are

```
xi  = 12 c/(u-v),
eta = 36 c (u+v)/(u-v),

u = (eta+36 c)/(6 xi),
v = (eta-36 c)/(6 xi).
```

For example, writing `h=u-v` gives

```
xi^3-432c^2
 = 432c^2 [4(u^2+uv+v^2)-h^2]/h^2
 = 1296c^2 (u+v)^2/h^2 = eta^2.
```

The inverse satisfies

```
u^3-v^3 = c (eta^2+432c^2)/xi^3 = c.
```

For positive `c`, no affine rational point has `xi=0`, since its square
would be negative. Thus every affine rational Mordell point is represented
by these formulas. The point at infinity corresponds to `[1:0:1]` on the
cubic. Every real affine Mordell point has `xi>0`. A positive sextic lift
necessarily has `eta>36c`; this inequality is a cheap sign filter, followed
by the exact square tests on `u,v`.

An exact rational cube-difference decision can therefore be reduced to
whether `M_c(Q)` contains any point other than its identity. A rank-zero
certificate can decide it, subject to the torsion exceptions below.

## 3. The cover genera and ramification

Let `D_X`, `D_Y`, `D_Z` denote the degree-three reduced geometric divisors
cut out by the three coordinate lines on `E_c`. They are pairwise disjoint
for `c != 0`. Each line meets the cubic in three distinct points. For

```
x=X/Z,   y=Y/Z,
```

the divisors are

```
div(x)=D_X-D_Z,    div(y)=D_Y-D_Z.
```

Consequently adjoining `sqrt(x)` gives a connected degree-two cover
branched at six geometric points, those of `D_X+D_Z`. Riemann--Hurwitz
gives `2g-2=6`, hence **genus four**. Adjoining only `sqrt(y)` likewise
gives genus four. The third intermediate double cover, adjoining
`sqrt(xy)`, branches at `D_X+D_Y`, and also has genus four.

The two square classes in the function field are independent: at a point
of `D_X`, `x` has odd valuation and `y` even valuation; at `D_Y` their
roles reverse. The composite has group `(Z/2Z)^2` and degree four. Over
each point of `D_X`, `D_Y`, or `D_Z`, its inertia has order two. In
particular, over `D_Z` both square roots change sign together; inertia is
the diagonal order-two subgroup, **not** a group of order four.

There are nine geometric branch points. Each has two points above it,
each of ramification index two, so the total ramification contribution is
`9*2=18`. Thus `2g(C_c)-2=18`, giving **genus ten**, as independently
confirmed by the smooth plane-sextic genus formula `(6-1)(6-2)/2=10`.

These are branched covers of an elliptic curve. They are not the usual
unramified covers arising directly from elliptic two-descent. Evaluation
`P |-> (x(P),y(P)) mod Q*^2` is not automatically a group homomorphism or
a finite Selmer map: the displayed divisors have odd multiplicities.
An ordinary 2-Selmer rank bound can help compute `E_c(Q)`, but does not
solve the two branched square conditions. A 3-isogeny or cubic-twist
descent may be a better rank computation for the j=0 elliptic curve;
it is still distinct from proving rational-point completeness on `C_c`.

## 4. Crucial limit on local filters: no whole fibre is locally empty

**Proposition.** For every positive rational `c`, `C_c` has real points
and Q_p-points with all three coordinates nonzero, for every prime `p`.
The same statement holds for `E_c` with `Y != 0`.

Proof. The rational boundary point `[1:0:1]` on `C_c` is smooth: in the
chart `A=1`, the derivative with respect to `B` of
`1+cS^6-B^6` equals `-6` there, which is nonzero in every Q_p and in R.
The local implicit-function theorem therefore supplies points with small
nonzero `S` and `B` close to 1. More explicitly, take `S` sufficiently
p-divisible that `1+cS^6` lies in a neighborhood of 1 consisting of sixth
powers, and take its sixth root. Over R, any sufficiently small positive
`S` works. Squaring the coordinates gives the corresponding cubic points.

Thus **neither local solubility nor unrestricted cube-difference
congruences can reject an entire `c` fibre**. This includes the full
two-square cover; it too has a rational smooth boundary point and nearby
nonzero points at every place. Local computation is useful for rejecting
residue classes, exact valuation charts, or finite-height search boxes.
It must not label an entire fibre `LOCALLY_INSOLUBLE` merely because
classes with `S` a unit are absent. Requiring primitive integer coordinates
does not repair this issue: the construction has `A=1` locally.

## 5. Exact normalization and reuse of elliptic arithmetic

The sextic coefficient is naturally considered modulo rational sixth
powers. If `c_raw=d^6*c0`, then

```
A^6+c0*S^6=B^6
```

corresponds to the raw fibre with `s=S/d`. Store that scale and the
original representations explicitly. In particular, when `d` divides
`c_raw` through sixth powers but does not divide every `ai`, the normalized
representation `ai/d` is rational, not an integer four-tuple.

For positive integer `c_raw`, the canonical positive sixth-power-free
integer is the product of `p^(v_p(c_raw) mod 6)`. Deduplicate by that
coefficient, preserving all original representations and scales. This
normalizes the stated scaling equivalence; it does not claim to classify
every abstract isomorphism of the curves.

For additional elliptic reuse write `c0=k*t^3`, with `k` positive
cube-free and `t` squarefree. Then

```
E_c0 -> E_k,   [X:Y:Z] |-> [X:tY:Z]
```

is an isomorphism over Q. Its square conditions become

```
X/Z is a square,   (Y_k/Z)/t is a square.
```

Equivalently, on cube-difference coordinates for `E_k`, both `t*u_k`
and `t*v_k` must be squares. The coefficient can therefore be indexed
by an elliptic cube class and a retained quadratic twist label. **Do not
discard that label when reusing a rank computation.** On Mordell models
the scaling is `xi_c0=t^2*xi_k`, `eta_c0=t^3*eta_k`.

## 6. Congruence and valuation charts forced by a primitive representation

First divide a raw four-tuple by its common gcd. For `p=7`, count its
entries not divisible by 7; call this count `r7`. For `p=3` do the same,
and for `p=2` count odd entries. Primitivity ensures each count is in
`{1,2,3,4}`. Sixth-power residues give

```
c = r7 (mod 7),  c = r3 (mod 9),  c = r2 (mod 8).
```

For each of these primes, a primitive fibre point has `B` a unit. If
the corresponding count is not 1, then `p|S` and `A` is a unit. If the
count is 1, exactly two charts are possible: `p|S` with `A` a unit, or
`p|A` with `S` a unit. This follows immediately by enumerating the
zero-or-one sixth-power residues in the equation. In particular,
`gcd(B,42)=1`.

For a primitive raw four-tuple,

```
v2(c) is 0, 1, or 2;
v3(c) is 0 or 1;
v7(c) is 0.
```

Sixth-power normalization preserves these valuations and residue counts,
since the extracted scale is a unit at 2, 3, and 7 and its sixth power is
1 in each displayed residue ring.

In the `p|S` chart put `k=v_p(S)>=1`. Since `A,B` are units, LTE gives

```
v2(B-epsilon*A) = v2(c)+6k-1,
v3(B-epsilon*A) = v3(c)+6k-1,
```

where for `p=2` choose `epsilon` so `4|(B-epsilon*A)`, and for `p=3`
choose it so `3|(B-epsilon*A)`. At 7, let `zeta` be the unique sixth
root of unity in Z_7 congruent to `B/A` modulo 7. Then

```
v7(B-zeta*A)=6k.
```

These exact congruences can accelerate bounded height searches. They do
not obstruct unbounded `k`, so do not close a fibre. The alternative
unit-`S` chart when the count is 1 is locally soluble: the representation
itself shows `c` to be a local sixth power, since one sixth-power term
is a unit and the other three are sufficiently p-divisible.

## 7. A genuine finite factor-allocation descent

Suppose `c` is a positive sixth-power-free integer and `(A,S,B)` is a
primitive integer point on `C_c`. Then the three integers are pairwise
coprime. If a prime divided `A,B` but not `S`, the equation would force
`v_p(c)>=6`, a contradiction. A prime dividing either `A,S` or `B,S`
would divide all three, also impossible.

Factor

```
c*S^6 = (B-A)(B+A)(B^2+AB+A^2)(B^2-AB+A^2).
```

The four factors are pairwise coprime outside 2 and 3. For a proof,
reduce modulo a common prime `p>3`: because `A,B` are units there,
`B/A` would have to be a root shared by two distinct cyclotomic factors
of `T^6-1`; that polynomial has distinct roots in characteristic `p>3`.

Therefore every prime `p|c`, `p>3`, is assigned to exactly one of the
four factors, with exponent congruent to `v_p(c)` modulo 6. Any other
prime greater than 3 occurs to an exponent divisible by 6. Thus, apart
from explicitly handled powers of 2 and 3, the factors have the form
`d_i*w_i^6`, where the `d_i` are pairwise coprime sixth-power-free
divisors whose product is the prime-to-6 part of `c`.

There are finitely many assignments: at most `4^omega(c_prime_to_6)`
before valuation and residue filters. Not every assignment is soluble;
each retains the compatibility identities between all four factors.
This is a legitimate finite descent organization for the sextic lift,
but enumerating assignments alone is not a proof of insolubility.
The primes 2 and 3 must be treated with the valuation formulas above,
not silently distributed as coprime factors. This branch can provide
certified exclusions if each resulting system is solved or locally
obstructed; the universally soluble boundary chart means that a naive
local test on the unseparated projective fibre cannot do so.

## 8. Certified rank-zero rejection and explicit torsion handling

A proved rank upper bound zero, together with an exact determination of
rational torsion, lists all of `M_c(Q)`. Neither an analytic-rank estimate
nor failure to find generators is by itself such a certificate.

One inexpensive per-fibre torsion certificate is to count `M_c(F_p)`
at several good primes `p>3`, with `p` not dividing `c`. Rational torsion
injects into good reduction at these primes, so its order divides the
gcd of the counts. A gcd of 1 proves torsion trivial. A small gcd can
be resolved by exact division-polynomial calculations.

For this particular model:

* Rational 2-torsion has `eta=0`, hence `xi^3=432c^2`. It exists exactly
  when `c=2*q^3` for rational `q>0`. Its cube-difference coordinates are
  `(u,v)=(q,-q)`, so it has no positive square lift.
* The 3-division polynomial is `3*xi*(xi^3-1728c^2)`. The root `xi=0`
  cannot give a real point. Other rational 3-torsion exists exactly
  when `c=q^3`, with `xi=12q^2`, `eta=+/-36c`. These give `v=0` or
  `u=0`, so are boundary points for the target problem.

These formulas, together with a reduction-order bound supported only
at 2 and 3, suffice to certify the torsion list without importing an
unproved global torsion classification. If a reduction-order bound
retains another prime, use more good primes or an actual torsion routine.
Do not assume a partial generator list is complete.

Final ledger replay independently proves complete torsion for all 1,984
certified elliptic models. Good-reduction order gcds suffice for 1,864.
For the other 120 the order bound is 3, whereas the listed torsion is
trivial. The exact polynomial `psi_3=3*x*(x^3+4*b)` on `y^2=x^3+b`
has no rational point over either possible x branch: `b` is not a square
and `-4*b` is not a cube. Hence no element of order three exists, and
the order bound proves trivial torsion. In particular all 573 rejecting
fibres, using 570 distinct elliptic models, now have an independent
complete torsion proof. The retained PARI rank upper bounds remain the
rank certificates; they were checked against their raw transcripts,
not recomputed by the independent verifier.

In positive rank, finite elliptic point enumeration is only a bounded
search. Finiteness of `C_c(Q)` follows from its genus, but does not give
an effective bound automatically. Potential global tools include
explicit descent on the factor assignments, the genus-four intermediate
covers and their Jacobians, or elliptic Chabauty after a correct cover
construction. Each needs hypotheses and certified arithmetic, rather
than a blanket `2-descent` label.

## 9. Genus-two intermediate curves and additional elliptic gates

The elliptic cube-difference quotient is not the only inexpensive gate.
First, if `c=h^3` is a positive rational cube, the original equation is
`(A^2)^3+(h*S^2)^3=(B^2)^3`. Fermat's theorem for exponent three excludes
nonzero rational solutions after denominator clearing. This is an exact
global coefficient gate, requiring no rank computation.
Define

```
H_c : w^2=t^6+c,
t=A/S,   w=(B/S)^3.
```

This is a smooth genus-two curve. The map `C_c -> H_c` has degree three,
obtained by adjoining a cube root of `w`. The six geometric zeros of `w`
are simple, so each ramifies with index three. The two points at infinity
have pole order three for `w` and are unramified. Riemann--Hurwitz gives
`2g(C_c)-2=3*(2g(H_c)-2)+6*(3-1)=18`, independently reproducing genus ten.

Two degree-two elliptic quotients of `H_c` are

```
E_plus : y^2=x^3+c,     (x,y)=(t^2,w),
E_sq   : y^2=x^3+c^2,   (x,y)=(c/t^2,c*w/t^3).
```

The involutions are respectively `(t,w)->(-t,w)` and
`(t,w)->(-t,-w)`. The former fixes the two points over `t=0`; the latter
fixes the two points at infinity. Their product is the hyperelliptic
involution. Each degree-two map has two ramification points, consistently
with Riemann--Hurwitz. The induced product of elliptic quotients gives
the usual isogeny decomposition of the genus-two Jacobian; in particular
its rank is the sum of their ranks. Consequently classical genus-two
Chabauty is not automatically available once both elliptic ranks are
positive.

Applying the same construction with the roles of the sextic coordinates
permuted gives all the following rational elliptic quotient maps. Every
map in this table has degree six from `C_c`. Conditions in the last column
mean **nonzero rational squares and nonzero rational cubes**, with signs
chosen so that the reconstructed `A,S,B` are positive.

| Elliptic equation | Coordinates of image | Exact lift conditions |
|---|---|---|
| `y^2=x^3+c` | `x=(A/S)^2`, `y=(B/S)^3` | `x` square, `y` positive cube |
| `y^2=x^3-c` | `x=(B/S)^2`, `y=(A/S)^3` | `x` square, `y` positive cube |
| `y^2=x^3+c^2` | `x=c(S/A)^2`, `y=c(B/A)^3` | `x/c` square, `y/c` positive cube |
| `y^2=x^3-c^3` | `x=c(B/A)^2`, `y=c^2(S/A)^3` | `x/c` square, `y/c^2` positive cube |
| `y^2=x^3+c^3` | `x=-c(A/B)^2`, `y=c^2(S/B)^3` | `-x/c` square, `y/c^2` positive cube |

Every entry follows by substitution. Conversely the stated exact lift
conditions reconstruct ratios of `A,S,B`, and the elliptic equation
recovers the sextic identity. In particular these are sufficient as well
as necessary lift tests. Points at infinity or zero ratios are boundary
points and do not satisfy the positive target.

Independent exact Laurent-polynomial checks in this audit, using Python
integer exponent tuples and `fractions.Fraction` coefficients, gave the following
residuals, writing `F=B^6-A^6-c*S^6`: in table order,
`F/S^6`, `-F/S^6`, `c^2*F/A^6`, `-c^3*F/A^6`, and
`-c^3*F/B^6`. The inverse cube-difference substitution gave
`c*(eta^2-xi^3+432c^2)/xi^3`. These checks verify the maps algebraically;
they do not assert that any lifted positive rational point was found.

For example the fourth row satisfies

```
y^2=c^4 S^6/A^6
   =c^3(B^6/A^6-1)=x^3-c^3,
```

and the fifth satisfies

```
y^2=c^4 S^6/B^6
   =c^3(1-A^6/B^6)=x^3+c^3.
```

Thus a certified rank-zero and full torsion calculation on **any one**
of these curves, followed by failure of all its torsion points to meet
the row's exact lift test, rejects the entire fibre. Requiring positive
rank on only the cube-difference model would miss these additional cheap
global obstructions.

The `+c^2` curve always has rational 3-torsion `(0,+/-c)`. It is boundary,
so torsion nontriviality is not a reason to retain a fibre. The `-c^3`
curve always has 2-torsion `(c,0)`, and the `+c^3` curve has 2-torsion
`(-c,0)`; these are likewise boundary. Every torsion point should pass
through the explicit lift test rather than a rank-only heuristic.

The coefficients `+/-c^3` depend up to rational Weierstrass scaling only
on the **square class** of `c`. If `c=d*q^2`, use
`x=q^2*x_d`, `y=q^3*y_d` to identify the models with coefficient
`+/-d^3`. These rank calculations may therefore be cached much more
aggressively than the full sixth-power-free coefficient. The lift tests
must still use the retained scale. In particular, square coefficients
reduce to `y^2=x^3+/-1`; once their rank-zero/torsion certificates are
recorded, they provide a uniform global gate for all square `c`.

The `+c^2` gate is not the same elliptic isogeny class as the original
cube-difference gate for `c` in general: its standard 3-isogenous curve
has coefficient `-27c^2`, which becomes `M_(2c)` after the rational
scaling `xi=4x`, `eta=8y`. This may give another safe cache reuse, but
the explicit quotient map is sufficient for rejection without deriving
any isogeny formulas.

## 10. Three product quotients, full Jacobian splitting, and the next global gate

Three further quotient maps are useful. Unlike the five monomial maps in
Section 9, these have degree **twelve** from `C_c` to an elliptic curve.

| Elliptic coefficient in `y^2=x^3+b` | x | y |
|---|---|---|
| `4c` | `(A^2/(B*S))^2` | `(B^6+c*S^6)/(B^3*S^3)` |
| `-4c` | `(B^2/(A*S))^2` | `(c*S^6-A^6)/(A^3*S^3)` |
| `4c^4` | `(c*S^2/(A*B))^2` | `c^2*(B^6+A^6)/(A^3*B^3)` |

For the first and third map, positive coordinates give `y>0`. In every
case `x` must be a nonzero rational square. Write `r=sqrt(x)>0`. The
following are exact inverse tests, including all remaining lift conditions:

* `b=4c`: require `y>0`; set `t=(y+r^3)/2`, `B=cuberoot(t)`, `S=1`,
  `A=sqrt(r*B)`.
* `b=-4c`: set `t=(r^3-y)/2`, `A=cuberoot(t)`, `S=1`, `B=sqrt(r*A)`.
  Both signs of `y` must be considered separately.
* `b=4c^4`: require `y>0`; set `t=(y+r^3)/(2c^2)`, `B=cuberoot(t)`,
  `A=1`, `S=sqrt(r*B/c)`.

All roots must exist in Q and be positive. Direct substitution then gives
`A^6+c*S^6=B^6`. The independent verifier uses sixth-root tests on the
reconstructed powers instead of the final square-root formulas, and agrees
with production on the synthetic sign/scale controls. A rank-zero fibre
gate can therefore check every torsion point completely, including those
with positive x; retaining every positive-x point would be unnecessarily
weak.

The genus-two origin of the third map is particularly transparent:

```
u=AB/S^2,  v=(B^6+A^6)/S^6,
H_AB : v^2=4u^6+c^2.
```

This map `C_c -> H_AB` has degree six. A preimage requires a sixth root
of `(v-c)/2`, after setting `S=1`; the other coordinate is recovered
from the product u. The six points of `C_c` with `S=0` each have
ramification index two, and there is no further ramification:
Riemann--Hurwitz reads `18=6*(2*2-2)+6`. Each of the two elliptic
quotients of `H_AB` has degree two. Its first elliptic quotient is
`y^2=x^3+16c^2`, using `x=4u^2`, `y=4v`; its second is the `4c^4`
row above. The `16c^2` elliptic curve is 3-isogenous over Q to the
original cubic-difference model `y^2=x^3-432c^2`.

The other product constructions are

```
u=AS/B^2, v=(c*S^6-A^6)/B^6, v^2=1-4c*u^6;
u=BS/A^2, v=(B^6+c*S^6)/A^6, v^2=1+4c*u^6.
```

Inverting u and dividing v by `u^3` gives the monic genus-two models
`w^2=t^6-4c` and `w^2=t^6+4c`. Their elliptic factors have coefficients
`-4c,16c^2` and `4c,16c^2`, respectively. These maps also have degrees
six then two. For example, the first new elliptic x-function has pole
degree 24 on `C_c`, whereas x has pole degree two on an elliptic curve;
this independently checks degree twelve and prevents confusing it with
the earlier degree-six maps.

### A complete decomposition over Q, proved by differentials

Write `E_b:y^2=x^3+b`. There is a Q-isogeny

```
Jac(C_c) ~ E_c * E_(-c) * E_(c^2)^2 * E_(-c^3) * E_(c^3)
             * E_(-432c^2) * E_(4c) * E_(-4c) * E_(4c^4).
```

Here the subscripts denote Weierstrass coefficients, so `E_(-432c^2)`
is the earlier Mordell model; this notation does not reuse the projective
cubic symbol `E_c` from Section 1. The statement has ten elliptic factors
and nine coefficient formulas. At special c some factors can be further
isogenous or isomorphic, but the displayed multiplicities still hold.

To prove it directly, work in the affine chart `S=1`, where
`B^6=A^6+c` and `dB=A^5/B^5*dA`. A basis of regular differentials is

```
omega_(a,b) = A^(a-1)*B^(b-6)*dA,
a,b >= 1, a+b <= 5.
```

This is the standard plane-sextic basis of monomials of total degree at
most three times `dA/B^5`, and contains ten forms. The pullbacks of
`dx/y` under the explicit elliptic maps are:

| Map | Pullback |
|---|---|
| coefficient `c` | `2*omega_(2,3)` |
| coefficient `-c` | `2*omega_(3,2)` |
| coefficient `c^2`, first map | `-2*omega_(1,3)` |
| coefficient `-c^3` | `-2*omega_(1,2)` |
| coefficient `c^3` | `-2*omega_(2,1)` |
| coefficient `-432c^2` | `(2/3)*omega_(2,2)` |
| coefficient `4c` | `2*omega_(4,1)` |
| coefficient `-4c` | `-2*omega_(1,4)` |
| coefficient `4c^4` | `-2*omega_(1,1)` |
| coefficient `c^2`, second map | `2*omega_(3,1)` |

The second map to `E_(c^2)` is `x=-c/B^2`, `y=c*A^3/B^3`.
Every table entry was independently checked by exact Laurent-polynomial
arithmetic modulo `B^6-A^6-c`; the verifier checks that the ten indices
are precisely the full genus-ten basis. The nonmonomial denominator in
the cubic-difference map is checked separately. Its derivative identity
is `dx/y=-D'/(3DN)`, where `D=B^2-A^2`, `N=B^2+A^2`, yielding the
displayed factor `2/3`.

All maps are defined over Q and extend from rational maps to morphisms
between the smooth projective curves. Translating their images of the
rational boundary basepoint to the elliptic origins gives homomorphisms
from `Jac(C_c)` to the product. The displayed independent pullbacks show
that their product has invertible differential between abelian varieties
of dimension ten. Its kernel is finite, proving the Q-isogeny.

Consequently the exact Mordell--Weil rank is the sum of the nine elliptic
ranks with the `c^2` rank counted twice. Certified individual upper and
lower bounds can be summed in exactly the same way. If all nine ranks
are proved positive, the Jacobian rank is at least ten, so the usual
rank-less-than-genus Chabauty criterion fails. Failure of that criterion
is not a proof that every more refined p-adic method fails.

### Quadratic Chabauty applicability and its limitations

There are at least eleven independent rational Neron--Severi classes:
eight singleton factor classes, and the two axes plus diagonal on the
repeated `E_(c^2)^2`. Pullback through the isogeny preserves their
independence over Q. Thus the rational Neron--Severi rank satisfies
`rho_Q >= 11`, without using extra geometric CM endomorphisms that might
not be defined over Q.

Balakrishnan--Dogra, *Quadratic Chabauty and rational points I*,
[Lemma 3.2](https://arxiv.org/pdf/1601.00388), proves finiteness of the
level-two Chabauty--Kim set when `r < g+rho_Q-1`, for a smooth projective
curve over Q with a rational basepoint and a prime of good reduction.
Our boundary supplies the basepoint, and any prime not dividing `6c`
is good. Therefore a certified Jacobian rank upper bound at most 19
satisfies this criterion. This is a finiteness/applicability result,
not an executed point-classification algorithm. Their more explicit
Theorem 1.2 has additional assumptions, including rank equal to genus
and a full-dimensional p-adic closure; it cannot simply be invoked for
all ranks between eleven and nineteen.

A more practical next target is a **genus-two quotient with two rank-one
elliptic factors**. The same paper's Theorem 1.4 explicitly treats monic
bielliptic curves `w^2=t^6+a4*t^4+a2*t^2+a0` when each elliptic factor
has rank one, using non-torsion points, p-adic logarithms and heights,
and the finite set of bad-prime local-height differences. This produces
a finite p-adic superset; rational reconstruction and a Mordell--Weil
sieve remain part of a complete rational-point determination. No such
height or integral calculation was executed in this audit.

Five immediately available monic quotient models and their rank pairs
are shown below. They all have rational points at infinity because the
leading coefficient is one.

| Monic genus-two equation | Elliptic rank pair |
|---|---|
| `w^2=t^6+c` | `(r_c,r_(c^2))` |
| `w^2=t^6-c` | `(r_(-c),r_(c^2))` |
| `w^2=t^6+4c` | `(r_(4c),r_(-432c^2))` |
| `w^2=t^6-4c` | `(r_(-4c),r_(-432c^2))` |
| `w^2=t^6+c^2/4` | `(r_(-432c^2),r_(4c^4))` |

The last three pairings use only the explicit rational model scalings
and the standard rational 3-isogeny already identified above. A `(1,1)`
pair certified in the ledger is therefore a concrete lower-genus
quadratic-Chabauty target. Point-classification on such a quotient must
still retain its exact cover condition back to `C_c`.

For instance, `c=67` has representation `(1,1,1,2)` and the retained
certified ranks of `y^2=x^3+67` and `y^2=x^3+4489` are both one.
On `H_67:w^2=t^6+67`, the exact rational point

```
t=7/6,  w=1801/216
```

satisfies `1801^2=7^6+67*6^6`. It fails the required cube condition on w
because 1801 is not a cube. In this original quotient, every positive
sextic point maps to rational `t=A/S`, `w=(B/S)^3`; thus classifying
`H_67(Q)` and checking which w are cubes would rigorously decide that
entire fibre. The point at infinity and any zero coordinate remain
boundary. The known point is a useful control, not a counterexample.

## 11. Corrections to inherited geometric language

The current night report is useful historical evidence but contains
claims that the new arithmetic pipeline must not treat as proved:

1. Section 1.2 invokes the Lander--Parkin--Selfridge **conjecture** to
   rule out rational `(6,1,4)` points. Without a separate proof this is
   conditional, and cannot certify a closed lane.
2. A homogeneous coordinate form having no rational root need not be
   irreducible over Q: a product of irreducible quadratic forms is an
   immediate counterexample to that implication. A rational linear
   factor forces a rational boundary contact; more general factors
   force contacts over their residue fields.
3. An odd-degree real factor forces a real zero and hence, for the
   stated sum of four real sixth powers, a real boundary contact.
   Rationality of that parameter requires additional hypotheses.
   The restricted linear-form argument can provide them; an unrestricted
   claim for all rational curves cannot simply inherit them.
4. A finite numerical lottery through every bounded seed is not a
   certified complete curve search merely because many random starts
   were used. Failed numerical recognition is not exact elimination.

These qualifications do not rerun or invalidate independently certified
finite calculations. They prevent conjectural or numerical premises from
being propagated into new proof labels.

## Audit conclusion

The proposed elliptic quotient is correct, but its simultaneous square
lift is genus ten and locally soluble for every coefficient. The strongest
cheap whole-fibre gate is therefore **certified elliptic rank/torsion
rejection across all nine displayed elliptic coefficient models**, with rank
computations reused by the appropriate cube or square class and lift
scales retained. Local arithmetic remains useful inside exact
valuation charts and finite descent assignments. All positive survivors
must return through the exact rational-square and integer-verifier gate.
