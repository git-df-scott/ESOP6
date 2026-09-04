# Degree-six boundary contact: an impossibility certificate

The normalized polynomial identity proposed in Strike #1 is **impossible
over Q**. The proof includes a finite coefficient certificate modulo 729,
reproduced by separate Python and C++ implementations. It allows independent,
asymmetric A and B and arbitrary rational coefficients, including denominators
divisible by 3. It is not just a test of integral coefficients.

The exact scope is

```text
X=t A(t), Y=t B(t), W=Z+(2/3)t^6,
A(0)=B(0)=Z(0)=1,
deg A,deg B<=5, deg Z<=6.
```

This does not exclude all rational curves on the surface, higher degrees,
or degree-six curves with a different ratio of the two leading tangent
coordinates. The equality of those leading coordinates is a restriction,
not a projective normalization available for every curve.

## 1. Contact, midpoint, and minimal degrees

At B=[0:0:1:1] the gradient of
`2X^6+2Y^6+Z^6-W^6` is `(0,0,6,-6)`. The surface is smooth; its tangent
plane is W=Z. A path with X=t+O(t^2), Y=t+O(t^2), Z=1+O(t) has

```text
2X^6+2Y^6 = 4t^6+O(t^7),
W^6-Z^6 = 6(W-Z)+O(t*(W-Z),(W-Z)^2).
```

Thus the contact order is exactly six and the leading coefficient of W-Z
is 2/3. For degree-six homogeneous coordinate forms, a section with this
order of vanishing is necessarily `(2/3)t^6` in the chosen affine chart.
There are no additional correction coefficients at this degree.

Put

```text
M=Z+t^6/3, U=(A+B)/2, V=(A-B)/2,
D=(A^6+B^6)/2=U^6+15U^4V^2+15U^2V^4+V^6.
```

The exact divided equation is

\[
D=M^5+\frac{10}{27}t^{12}M^3+\frac1{81}t^{24}M
 =\frac{M(3M^2+t^{12})(27M^2+t^{12})}{81}.                 \tag{1}
\]

In the original coordinates, before the midpoint substitution,

\[
A^6+B^6=2Z^5+\frac{10}{3}t^6Z^4+\frac{80}{27}t^{12}Z^3
 +\frac{40}{27}t^{18}Z^2+\frac{32}{81}t^{24}Z+\frac{32}{729}t^{30}.
\]

Degree six is necessary in this polynomial construction. More precisely,
write r=max(deg A,deg B) and h=deg M. The real leading coefficients in
A^6+B^6 cannot cancel. If h<6, the right side of (1) has degree 24+h,
so `6r=24+h` forces h=0 and r=4. Then M=1 and the leading equation is
`a4^6+b4^6=2/81`, impossible: the 3-adic valuation of a sum of two rational
sixth powers is a multiple of six, whereas `v3(2/81)=-4`.
If h>6, the term M^5 alone has degree greater than 30. Consequently any
putative solution needs h=6 and r=5. This proves that lowering both A and B
below degree five loses every possibility; it does not require each to
have degree five separately.

The leading coefficient m6 of M must be positive, since
`m6^5+(10/27)m6^3+m6/81` has its sign. If Z had degree below six, then
`m6=1/3`, giving `a5^6+b5^6=32/729`. Its 2-adic valuation is 5, while such
a two-term sum has valuation 0 or 1 modulo six. Thus Z also needs degree six.
These degree conditions are necessary, not sufficient.

## 2. Freedom count and triangular equations

Projective scaling sets Z(0)=1; rescaling the parameter sets the first
coefficient of X to 1. Requiring the first coefficient of Y to be 1 selects
the equal-slope direction. The leading equation then fixes c=2/3.

There are 16 remaining coefficients: five in each of A,B and six in Z.
After division by `4t^6`, the constant equation is automatic and there are
30 coefficient equations, indexed 1 through 30. This is an overdetermined
system; equation count alone is not a nonexistence proof.

One continuous redundancy remains. The transformation

```text
t -> t/(1+k t),
(X,Y,Z,W) -> (1+k t)^6 (X,Y,Z,W)(t/(1+k t))
```

preserves all the specified leading constants and W-Z. It changes
`u1 -> u1+5k` and `m1 -> m1+6k`. Choose `k=-u1/5`. The first coefficient
equation `5m1=6u1` then gives u1=m1=0. Exchange of X and Y sends V to -V;
it is discrete and does not remove a continuous parameter. Adding a common
polynomial to A and B is not another symmetry of (1).

Hence write

```text
U=1+u2 t^2+u3 t^3+u4 t^4+u5 t^5,
V=v1 t+v2 t^2+v3 t^3+v4 t^4+v5 t^5.
```

For each n>=1, let `M_<n=1+sum_{j<n} m_j t^j`. The exact increasing-order
coefficient equation is

\[
5m_n=[t^n]\left(D-M_{<n}^5-\frac{10}{27}t^{12}M_{<n}^3
                         -\frac1{81}t^{24}M_{<n}\right).                \tag{2}
\]

For n=1,...,6 this solves a coefficient of Z linearly. For n=7,...,30,
set m_n=0; the right sides must vanish. These are the exact 24 residual
equations in **nine** effective parameters. In particular:

```text
m1=0
m2=3(2u2+5v1^2)/5
m3=6(u3+5v1v2)/5
m4=3(u2^2-20u2v1^2+10u4-125v1^4+50v1v3+25v2^2)/25
m5=6(u2u3-20u2v1v2-10u3v1^2+5u5-250v1^3v2+25v1v4+25v2v3)/25
m6=-(4u2^3-270u2^2v1^2-30u2u4-5250u2v1^4+600u2v1v3
     +300u2v2^2-15u3^2+600u3v1v2+300u4v1^2-15775v1^6
     +7500v1^3v3+11250v1^2v2^2-750v1v5-750v2v4-375v3^2)/125
```

Recover `Z=M-t^6/3`. None of these equations forces V=0.
The complete factored coefficients at orders 7 and 8 are retained in
`results/astra_second_2026_09_04/symbolic_checks.json`. They are linear in
u5,v5, with coefficient matrix

\[
\begin{pmatrix}
6(u_2-10v_1^2)/5&30v_2\\
6(u_3-20v_1v_2)/5&-6(4u_2v_1+50v_1^3-5v_3)
\end{pmatrix}.
\]

Its determinant is `-36/5` times

```text
4u2^2v1+10u2v1^3-5u2v3+5u3v2-500v1^5+50v1^2v3-100v1v2^2.
```

When nonzero it eliminates u5,v5; its zero set would require a separate
branch. Rather than discard that branch or perform a generic elimination,
the following arithmetic argument handles both at once.

## 3. Invert the parameter and control every denominator

Set

```text
x(s)=s^5 A(1/s), y(s)=s^5 B(1/s), m(s)=s^6 M(1/s).
```

Now x,y are monic quintics, m is monic of degree six, and

\[
x^6+y^6=2m^5+\frac{20}{27}m^3+\frac2{81}m.               \tag{3}
\]

The normalization above says that the coefficient of s^4 in x+y is zero,
and the coefficient of s^5 in m is zero.

Use the 3-adic Gauss valuation: the minimum valuation of a polynomial's
coefficients. A sum of two sixth powers has valuation exactly six times
the minimum of the two polynomial valuations. Indeed, after reduction,
cancellation would give a square root of -1 in F3(s), which does not exist.
The same fact holds for weighted Gauss valuations
`w_r(P)=min_j(v3(P_j)+j r)` with rational r: pass to a totally ramified
extension to realize the weight, whose residue field is still F3.

Let `b=w_0(m)<=0` and `a=min(w_0(x),w_0(y))`.

* b=0 would give `6a=-4`, impossible.
* b=-1 gives a=-1, and is excluded below.
* b<=-2 makes the term 2m^5 uniquely dominant in (3). Thus
  `6a=5b`, so `b=-6k, a=-5k` for an integer k>=1.

### The exceptional valuation b=-1 is impossible

Put `a0=3x, b0=3y, d=3m`, which are integral polynomials with at least
one of a0,b0 a Gauss unit. Equation (3) becomes

```text
a0^6+b0^6=6d^5+20d^3+6d.
```

Modulo 3, `d=-(a0^2+b0^2)`. Substitute this modulo 9; a change of d by
3h makes no difference modulo 9, since the derivative of the right side
is divisible by 3. With `S=a0^2+b0^2`, division by 3 and reduction gives

```text
S * (S^4-(a0^2-b0^2)^2+1)=0  in F3[s].
```

S is nonzero by anisotropy of the sum of two squares. If either a0 or b0
has positive degree d0, the three terms in the second factor have degrees
8d0, at most 4d0, and zero. They cannot sum to zero. For constants, one
nonzero entry makes that factor 1 and two nonzero entries make it 2.
This closes b=-1 without assuming any coefficient is initially integral.

### The remaining valuations can all be rescaled to one integral problem

For r<=0, `w_r(m)<=-6k`, so the same unique dominant term gives
`6 min(w_r(x),w_r(y))=5w_r(m)`. On each linear segment of these Newton
polygons the slopes are integers. The slope of w_r(m) must be a multiple
of six. Its only possibilities are 6 and 0. Consequently

```text
w_r(m)=min(6r,-6k),
min(w_r(x),w_r(y))=min(5r,-5k).
```

It follows that

```text
x0(s)=3^(5k) x(3^(-k)s),
y0(s)=3^(5k) y(3^(-k)s),
m0(s)=3^(6k) m(3^(-k)s)
```

are integral monic polynomials of degrees 5,5,6. The constant coefficient
of m0 is a unit, and the coefficient of s^4 in x0+y0 remains zero. Their
equation is

\[
\frac{x_0^6+y_0^6}{2}=m_0^5
 +10\,3^{12k-3}m_0^3+3^{24k-4}m_0.                    \tag{4}
\]

Both correction terms vanish modulo 3^9 for every k>=1. It suffices to
exclude the unperturbed identity `(x0^6+y0^6)/2=m0^5` modulo 3^6.

## 4. Modulo 9, every surviving polynomial has one residue shape

In this section bars mean reduction modulo 3. Set `S=x0^2+y0^2`.
The identity modulo 3 is `S^3=2m0^5`. Unique factorization implies

```text
bar(m0)=P^3, bar(S)=2P^5,
```

for a monic quadratic P in F3[s]. Choose an integral lift of P and write
`m0=P^3+3h`, `S=2P^5+3j`. Use
`x0^6+y0^6=S^3-3x0^2y0^2S` modulo 9. After division by 3 one obtains

```text
P^5*(2P^10-2 bar(x0)^2 bar(y0)^2-P^7 bar(h))=0,
```

and hence `P^7` divides `bar(x0)^2 bar(y0)^2`.

There are only three factorization types for P:

1. Two distinct linear factors are impossible: the valuation of a sum of
   two squares at a rational linear factor is even, whereas `2P^5` has
   odd valuation there.
2. If `P=(s-r)^2`, both monic quintics must equal `(s-r)^5`. The normalized
   s^4 coefficient in their sum forces r=0. This contradicts the unit
   constant coefficient of m0.
3. If P is irreducible, let its valuations in the two quintics be e and f.
   The odd valuation 5 in their square sum forces e=f. The divisibility
   above gives `4e>=7`, and the degrees give e<=2. Thus e=f=2. Write
   `bar(x0)=P^2(s+r)`, `bar(y0)=P^2(s+q)`. Their square sum determines
   `P=((s+r)^2+(s+q)^2)/2`. The normalized s^4 coefficient forces q=-r.
   Irreducibility then gives r=+/-1 and `P=s^2+1`.

Exchanging x0,y0 therefore leaves just one residue shape. Reversing s
back to the coefficient variable t, its nine coefficients are

```text
U=1+2t^2+t^4, V=t+2t^3+t^5  (mod 3),
(u2,u3,u4,u5,v1,v2,v3,v4,v5)=(2,0,1,0,1,0,2,0,1).
```

No parity constraint is imposed on its lifts: all nine coefficients acquire
independent subsequent 3-adic digits.

## 5. Exhaustive finite obstruction modulo 729

For a coefficient vector q, compute D=(A^6+B^6)/2 and the unique monic
degree-six midpoint jet through coefficient six using (2), omitting the
correction terms already zero modulo 3^9. Let R(q) be coefficients 7,...,30
of `D-M_6^5`. These are polynomials with denominators prime to 3.

Every first derivative of R is divisible by 3 as a polynomial: D has this
property, and the recursive fifth-root equations preserve it. Therefore,
for n>=1,

```text
R(q+3^n h) = R(q)+3^(n+1) L h  (mod 3^(n+2)),
L=(Jacobian R)/3 evaluated at q mod 3.
```

The second-order term has an extra factor 3 and all remaining terms vanish
at this precision. Thus this is an exact linear lifting rule, not a
first-order approximation used outside its precision. L is constant on
the one residue shape. Its rank is four, with pivot columns u2,u3,u4,u5;
each solvable step has exactly `3^5=243` next digits. The full 24-by-9
matrix and kernel basis are in the certificate JSON.

| Equation modulus | Coefficient prefixes checked | Prefixes admitting the next digit | Children |
|---:|---:|---:|---:|
| 27 | 1 | 1 | 243 |
| 81 | 243 | 81 | 19,683 |
| 243 | 19,683 | 162 | 39,366 |
| **729** | **39,366** | **0** | **0** |

The first unavoidable obstruction in this exhaustive arithmetic reduction
is at modulus 729. This is not a claim that one fixed bare t-coefficient
contradicts every unnormalized input. The lower coefficient equations have
already been absorbed, and every remaining coefficient class is rejected
by the simultaneous residual equations at this precision.

In total, 59,293 coefficient prefixes were checked, including the initial
one. No rational-height or integer-point search occurred. Since (4) would
have to hold even modulo 3^9, failure modulo 3^6 is a contradiction.
This proves the stated polynomial impossibility.

## 6. Independent reproduction and precise limit

```bash
python3 tools/boundary_contact.py
c++ -O2 -std=c++17 tools/verify_boundary_contact.cpp -o /tmp/verify-boundary
/tmp/verify-boundary
make boundary-checks
```

Python uses A,B sixth powers and a recursive fifth root. C++ uses the U,V
norm and the explicitly displayed midpoint coefficients. The linear systems
are also solved differently. All four prefix counts and checksums agree.
The common assumptions are the mathematical reduction proved above, not
an assumed integral-coefficient restriction or an assumed parity branch.

The formal branch survives; its polynomial truncation does not. The smallest
surviving rational-function extension in this fixed tangent direction is
specified in [RATIONAL_CURVE_ATTEMPT.md](RATIONAL_CURVE_ATTEMPT.md).
