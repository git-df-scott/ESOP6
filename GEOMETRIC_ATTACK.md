# Geometric / Elkies-style attack

> **Direct-strike update:** this file is the preparation proposal. The
> executed attack and exact scope are in [GEOMETRIC_STRIKE.md](GEOMETRIC_STRIKE.md).
> The proposed global genus-one fibration on the smooth sextic surface is
> obstructed by its ample canonical bundle; quotient curves must be lifted.
> The natural conic pencil gives genus 9. Degrees 1–3 are excluded, while
> general degree 4 is still unclassified here. The next degree-six ansatz
> advances only the boundary-contact branch, which is proved empty through
> degree five; it does not assert that all degree-four curves were excluded.

## The variety

The projective solution variety is

\[
X:\;x_1^6+x_2^6+x_3^6+x_4^6+x_5^6-x_6^6=0\subset\mathbf P^5.
\]

It is smooth in characteristic zero: all six partial derivatives vanish only
when every homogeneous coordinate is zero. It has dimension four, and
adjunction gives \(K_X=\mathcal O_X(6-6)=\mathcal O_X\). Together with the
standard intermediate-cohomology vanishing for a smooth hypersurface, this
makes it a Calabi–Yau fourfold.

Over \(\mathbf Q\), visible symmetries include permutation of the first five
coordinates and independent sign changes modulo projective common sign. Over
a field containing the relevant roots of unity, it becomes a diagonal Fermat
form with a larger monomial symmetry group.

The desired point is stronger than an arbitrary rational point: every
coordinate must be nonzero and, after signs, lie in the positive real chamber.
Boundary points such as `[0:0:0:0:1:1]` are useless. Their tangent section is
positive-definite in the four newly activated coordinates, explaining why the
naive tangent-line construction at an obvious point does not enter the
positive chamber.

## What Elkies actually supplied in degree four

Elkies did not merely search a quartic box. He reduced
\(A^4+B^4+C^4=D^4\) to rational points on an elliptic curve, found a
positive-rank point, and used the group law to generate infinitely many
solutions. See [Elkies's publication list](https://people.math.harvard.edu/~elkies/math_pubs.html)
and the cited 1988 *Mathematics of Computation* paper.

The reusable ingredients are a low-dimensional slice carrying a genus-one
fibration, a rational section or positive-rank specialization, an explicit
map back to the diagonal equation, and control of signs/nonzero coordinates.
What does not transfer automatically is the K3 surface geometry: `X` is a
fourfold, and an arbitrary surface section of a sextic is usually of general
type.

## Closest constructive literature

Bremner–Choudhry–Ulas systematically search low-degree curves on diagonal
quartic and sextic surfaces and reduce square/cube conditions to genus-one
curves. Their [2014 paper](https://arxiv.org/abs/1402.4583) gives explicit
sixth-power identities and the control identity

\[
(1-T-T^2)^3+(1+T-T^2)^3=2-2T^6.
\]

It does not produce a `(6,1,5)` point: after demanding the two cubic bases be
squares it has a balanced 4-versus-2 signature. It does provide a concrete
symbolic method and regression identity for an Astra geometry pipeline.

## One bounded construction problem

Start with the repeated-coordinate diagonal sextic surface

\[
S:\;2X^6+2Y^6+Z^6=W^6\subset\mathbf P^3.
\]

Any rational point with \(XYZW\ne0\) in the positive real chamber gives the
ESOP6 solution \((a,b,c,d,e,f)=(X,X,Y,Y,Z,W)\) after clearing denominators.
Repeated positive terms are allowed.

The Astra task is exact and bounded:

> For degrees \(d=1,2,3,4\), classify rational maps
> \([s:t]\mapsto[X(s,t):Y(s,t):Z(s,t):W(s,t)]\) with homogeneous binary
> forms of degree \(d\) satisfying `2X^6+2Y^6+Z^6-W^6=0`. Quotient by common
> scale and PGL2 normalization. For every rational component, determine
> whether it has a rational point mapping into `X,Y,Z,W>0`. If the direct
> curve ansatz is empty, search genus-one pencils on `S` obtained by the
> Bremner–Choudhry–Ulas square/cube substitution pattern, and compute section
> rank exactly.

Every outcome is certifiable: an explicit identity and point, or a Gröbner
basis/resultant certificate that the normalized ansatz is empty. Degree is
increased only after certifying the previous degree.

## Secondary full-fourfold lane

If the repeated surface dies through degree four, do not enumerate arbitrary
larger coefficients. Search rational surfaces on `X` using paired binary
forms, with the BCU identity as a control. The required deliverable is an
explicit elliptic fibration over \(\mathbf Q\) with a non-torsion section and
a specialization in the positive chamber—not a numerical point cloud.

No rational curve, surface, or elliptic fibration solving ESOP6 was found in
this audit. The strongest lead is the degree-1-through-4 classification on
`S`, because it is finite, exact, positivity-aware, and directly converts one
rational point into a counterexample.
