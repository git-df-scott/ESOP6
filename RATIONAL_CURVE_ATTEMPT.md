# Exact curve attempts and the smallest surviving extension

The complete prescribed polynomial family is impossible, by
[BOUNDARY_CONTACT_6.md](BOUNDARY_CONTACT_6.md). All first-strike artifacts
are retained. There is no new surface point, sextuple, rational curve, or
elliptic residual. Construction work stopped at condition **D**, after the
polynomial obstruction and a precise minimal rational extension were obtained.

## Bounded attempts and first obstructions

The residual in the following table is
`D-M^5-(10/27)t^12 M^3-t^24 M/81`, after the coefficients of M through
degree six have been solved.

| Restriction | First obstruction | Outcome |
|---|---|---|
| U=1, V=qt | `-6177 q^8 t^8` | q=0 is forced; then the t^12 coefficient is -10/27 |
| U=1, V=qt^2 | `-75 q^4 t^8` | same terminal q=0 obstruction |
| U=1, V=qt^3 | `-5(405q^4+2)t^12/27` | no real q |
| U=1, V=qt^4 | `15q^2 t^8` | q=0, then -10/27 at t^12 |
| U=1, V=qt^5 | `15q^2 t^10` | q=0, then -10/27 at t^12 |
| A(t)=B(-t), with Z free | leading point at infinity lies on X=-Y if a5!=0 | inherited diagonal obstruction forces a5=0; the remaining degree argument forces M=1 and `a4^6+b4^6=2/81`, impossible |
| M=1 with independent A,B | degree 24 requires `a4^6+b4^6=2/81` | 3-adic valuation -4 is not divisible by six |
| Full independent A5,B5,Z6 | no coefficient lift modulo 729 after the exhaustive denominator reduction | **impossible over Q**, including all determinant-zero branches |
| Rational Z only, polynomial A5,B5 | a reduced denominator Q must divide its coprime numerator P^5 | Q is constant; already covered |
| [1/1] midpoint reconstruction for A=B=1, in s=t^12 | exact residual begins `-121t^36/196830` | approximate contact only |

For the reflection row, uniqueness of the formal midpoint implies M is even.
If its top A coefficient is nonzero, the point at infinity has X=-Y!=0.
If that coefficient is zero, the degree-six identity at infinity forces
the leading midpoint coefficient to vanish. Then deg M<=4 and the degree
argument above gives M=1. This closes the larger reflection family with
free Z, not merely the previously closed constant-midpoint subfamily.

The five sparse asymmetric attempts, one Padé control, and unrestricted
3-adic reduction are the bounded coefficient work actually performed. No
generic rational-height search, Gröbner basis, external campaign, integer
frontier extension, or repetition of the 10,382-target search was run.

## Smaller mutations do not evade the obstruction

1. **One additional coefficient in Z:** with deg Z=7 and W-Z=(2/3)t^6,
   the highest nonzero term of W^6-Z^6 has degree 41, while the other side
   has degree at most 36. The extra coefficient must be zero. The same
   leading-degree argument excludes increasing Z alone further.
2. **A t^7 correction in W-Z:** keeping Z of degree at most six gives a
   unique t^42 term in W^6, so that correction must vanish. Raising both
   W and Z to seven cannot give a new reduced degree-seven real map.
3. **A common linear denominator with degree bounds at most seven after
   clearing it:** a reduced map of odd degree is impossible. A real zero
   of W would force X=Y=Z=0, contrary to having no common base point.
   Cancelling any common factor leaves degree at most six, already excluded
   in this tangent direction. If all four coordinates are merely divided
   by the same rational function, their projective curve is unchanged.

Thus a linear denominator does not create a missing degree-six or
degree-seven curve. This statement is about the stated degree bounds;
it does not exclude arbitrary high-degree numerators over a linear
denominator.

## The single strongest third direct strike

Use the smallest genuine extension that preserves the selected tangent
direction: reduced projective degree eight, with an irreducible positive
quadratic contact factor. After a Möbius normalization, write

\[
\begin{aligned}
Q(t)&=1+q t^2,\qquad q\in\mathbb Q_{>0},\\
X&=tP_7(t)/Q(t),&Y&=tR_7(t)/Q(t),\\
Z&=T_8(t)/Q(t),&W&=Z+\tfrac23t^6,\\
P_7(0)&=R_7(0)=T_8(0)=1.
\end{aligned}
\]

This is a precise **LIVE** construction, not an identity asserted to work.
The degrees are bounds; P7 and R7 remain independent. Clear Q and put
`N=T8+t^6 Q/3`. The exact residual problem is

\[
\frac{P_7^6+R_7^6}{2}
 = QN^5+\frac{10}{27}t^{12}Q^3N^3+\frac1{81}t^{24}Q^5N.       \tag{5}
\]

Its constant derivative in N is again 5. The same increasing-order
recurrence eliminates the eight coefficients of N first. There are 23
coefficients after centering Q, and 42 nonconstant equations; eliminating
N leaves 15 coefficients and 34 residual equations. No claim is made that
these equations are independent or that this count predicts a solution.
Attack (5) structurally before any coefficient enumeration.

Why this is the minimal extension here:

* Reduced degrees below six cannot leave B with the required contact;
  degree six in this normalized direction has now been excluded. Every
  odd reduced degree is impossible over R, so eight is next.
* In degree eight, W-Z is a section with its prescribed order-six zero;
  its remaining factor is a quadratic Q with Q(0)=1.
* Any further **real** zero of W-Z would force X=Y=0 and contact order at
  least six there. Only two remaining zeros are available, so Q has no
  real projective zero. In particular it is genuinely quadratic and
  positive after normalization.
* If initially `Q=1+l t+h t^2`, the admissible parameter change
  `t -> t/(1-l t/2)` replaces it by `1+(h-l^2/4)t^2`. Thus q>0 is a
  normalization, not an extra parity restriction on P7,R7 or N. The one
  Möbius freedom is used here, so do not additionally set the linear
  coefficient of (P7+R7)/2 to zero without justification.

For an exact identity in (5), Q is positive on the entire real line and
the three other numerator constant terms are 1. The explicit coefficient
bound in [FORMAL_BRANCH.md](FORMAL_BRANCH.md) then supplies a rational
positive interval and an immediate integer sextuple specialization.

“Minimal” here means minimal reduced projective degree and a common
denominator of least positive degree in the **fixed equal-slope contact
direction**. Degree-six constructions with unequal leading slopes and
rational curves avoiding B are outside this proof. Their status has not
been changed. The rational-function architecture as a whole remains open.
