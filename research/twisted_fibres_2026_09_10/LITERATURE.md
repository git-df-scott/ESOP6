# Sixth-power literature audit and transfer obstructions

Date: 2026-09-10. No positive-integer counterexample is established here.

## Source access and limits

The two requested works are identified precisely:

| Work | Primary bibliographic link | Access during this session |
|---|---|---|
| Andrew Bremner, *A geometric approach to equal sums of sixth powers*, Proceedings of the London Mathematical Society (3) **43** (1981), 544–581 | [DOI 10.1112/plms/s3-43.3.544](https://doi.org/10.1112/plms/s3-43.3.544) | Crossref metadata obtained; publisher full text returned HTTP 403; public TDM endpoint requires an access token. **Original paper not obtained or read.** |
| Masato Kuwata, *Equal sums of sixth powers and quadratic line complexes*, Rocky Mountain Journal of Mathematics **37** (2007), 497–517 | [DOI 10.1216/rmjm/1181068763](https://doi.org/10.1216/rmjm/1181068763); [author bibliography](https://www.kuwata.r.chuo-u.ac.jp/papers_en.html) | Metadata and author bibliography obtained. Project Euclid full-text/PDF URLs returned an HTML bot challenge instead of the article. **Original paper not obtained or read.** |
| Andrew Bremner, Ajai Choudhry and Maciej Ulas, *Constructions of diagonal quartic and sextic surfaces with infinitely many rational points* (2014 preprint) | [arXiv:1402.4583](https://arxiv.org/abs/1402.4583); [PDF](https://arxiv.org/pdf/1402.4583) | PDF obtained; sixth-power sections 5–6, pp. 13–19, read in full. |
| L. J. Lander, T. R. Parkin and J. L. Selfridge, *A survey of equal sums of like powers*, Mathematics of Computation **21** (1967), 446–459 | [DOI 10.1090/S0025-5718-1967-0222008-0](https://doi.org/10.1090/S0025-5718-1967-0222008-0) | Crossref metadata obtained; AMS PDF returned HTTP 403. **Original paper not obtained or read.** |

The bibliography alone does not establish any claimed Néron–Severi lattice,
Picard rank, involution, elliptic pencil or descent theorem in either requested
original. Those details remain a literature-access dependency. No such numerical
lattice data are imported into the search as a theorem. It is also unsupported
to call these the *only* published elliptic-curve attacks on sixth powers.

## What the available primary paper actually supplies

BCU §5 searches weighted polynomial identities on surfaces of types
`P X^6 + Q Y^6 + R Z^3 + S W^2 = 0` and
`P X^6 + Q Y^6 + R Z^3 + S W^3 = 0`; extra square or cube conditions lift these
to other exponents. In §6 the construction first obtains a quartic polynomial
`F(t)` and then chooses a coefficient from `F(t0)` so that the auxiliary quartic
has the rational point `(t,W)=(t0,1)`. This does not produce a point for an
arbitrary prescribed coefficient. The displayed sixth-power example has mixed
signs, `2X^6 - 2Y^6 + Z^6 = F(t0)^3 W^6`. Matching a chosen ESOP6 slice requires
a separate coefficient match and a rational lifting argument. No universal
square-cover descent follows from this construction.

These statements summarize §§5–6 of the [primary preprint](https://arxiv.org/pdf/1402.4583),
not the inaccessible 1981 or 2007 originals. Existing repo computations based
on that control should be retained; this is an audit of applicability, not a
request to rerun them.

## An elementary identity that can be checked without either missing paper

The following derivation is independent algebra. It is the appropriate way to
verify a proposed balanced `(2,6)` auxiliary construction before importing its
geometry.

Put

\[
q_+(a,d)=a^2+ad-d^2,\qquad q_-(a,d)=a^2-ad-d^2.
\]

Then

\[
q_++q_-=2(a^2-d^2),
\]

and

\[
q_+^3+q_-^3
=2(a^2-d^2)\big((a^2-d^2)^2+3a^2d^2\big)
=2(a^6-d^6).
\]

Consequently the three quadratic equations

\[
q_+(a,d)+q_-(b,e)=0,
\quad q_+(b,e)+q_-(c,f)=0,
\quad q_+(c,f)+q_-(a,d)=0
\]

imply both

\[
a^2+b^2+c^2=d^2+e^2+f^2,
\qquad a^6+b^6+c^6=d^6+e^6+f^6.
\]

The implication is one way: equality of the two power sums does not by itself
prove that these particular three quadrics vanish. Nor does the elementary
identity establish any claimed surface isomorphism, its field of definition,
or an effective rational-point algorithm.

## Certified obstruction to the obvious Gaussian transfer

Let `M` be the simultaneous-square/sixth-power locus in the previous display.
Over `Q(i)`, make the diagonal change

\[
(a,b,c,d,e,f)=(x_1,x_2,x_3,i x_4,i x_5,x_6).
\]

Since `i^2=i^6=-1`, the two equations become

\[
\sum_{j=1}^{5}x_j^2=x_6^2,
\qquad \sum_{j=1}^{5}x_j^6=x_6^6.
\tag{G}
\]

**The real locus of (G) consists only of the ten projective points with one
nonzero left coordinate equal to `±x6`. In particular it contains no positive
ESOP6 point.**

Proof: set `t_j=x_j^2 >= 0`. Equations (G) give

\[
0=\left(\sum_j t_j\right)^3-\sum_j t_j^3
=3\sum_{i\ne j}t_i^2t_j+
6\sum_{i<j<k}t_it_jt_k.
\]

All summands are nonnegative. If two of the `t_j` are positive, at least one
term in the first sum is positive, a contradiction. The projective all-zero
case is excluded. This proves the assertion without LPS or any unproved
arithmetic conjecture.

Thus diagonal Gaussian transfer of any family confined to `M` cannot solve
ESOP6. A different twist or construction must be written explicitly and its
real locus checked anew; geometric equivalence over an extension does not
establish rational equivalence or preserve the sought real component.

## The Gaussian seed lies exactly on this trap

The campaign's Gaussian seed satisfies the additional quadratic equation:

\[
(7+11i)^2+(7-11i)^2+8^2+12^2+15^2
=-144+64+144+225=289=17^2.
\]

So its attractive `8,15,17` relation places it on (G), whose real locus was
just classified. **Any proposed family through this seed that preserves the
same quadratic equation is incapable of yielding an all-nonzero rational
point.** A useful family through the seed must leave that equation.

In swap-quotient coordinates `e1=x1+x2`, `e2=x1*x2`, the seed has
`(e1,e2)=(14,170)` and discriminant `e1^2-4e2=-484`. A rational discriminant
square would lift to real coordinates. On (G) such a lift is forced to be one
of the trivial points; the discriminant is not free to become a square at an
admissible specialization. This rejects an entire family of quotient
constructions, rather than merely a bounded parameter sweep.

The boundary `Q(sqrt(-249))` seed does not satisfy this extra quadratic
equation: its left square sum is `2(81-249)+196+324=184`, whereas `22^2=484`.
The present obstruction therefore does not eliminate all families through
that seed.

## Correction of the inherited slice-contact statements

`research/NIGHT_REPORT_2026_09_10.md` §1.2 explicitly calls LPS a conjecture,
then reasons as though its boundary nonexistence prediction were a theorem.
The following distinctions are necessary for certified pruning:

1. A rational root of a left coordinate form of a morphism `P^1 -> X`
   produces a rational point on that coordinate slice. Unless it is trivial,
   this would be an additional sixth-power problem with at most four nonzero
   left terms. Declaring it impossible on the basis of LPS is **conditional**.
2. A reducible homogeneous polynomial need not have a rational root. For
   example, `(s^2+t^2)(s^2+2t^2)` is reducible over `Q`, with no rational or
   real zero in `P^1`. Reducibility alone does not invoke rational boundary
   contact. A factor of degree `r` usually gives an algebraic contact point
   of degree up to `r`, and LPS over `Q` says nothing about that point.
3. Therefore neither the asserted irreducibility of all coordinate forms nor
   an unconditional universal minimum degree six follows from the cited
   slice-contact argument. Each excluded family must retain its actual
   independent identity, congruence certificate or real-contact proof.
4. The low-degree **real contact** argument below *is* unconditional and
   should not be conflated with the conditional arithmetic assertion.

For completeness, here is the real-contact proof in the form useful for the
odd-factor lane. Suppose binary forms `F1,...,F6` of common degree `d` have
no common projective zero and satisfy the sextic identity. At a real
parameter where `F6=±F5`, positivity gives `F1=...=F4=0`. Necessarily
`F5*F6 != 0`, since otherwise all six vanish. If
`m=min(ord(F1),...,ord(F4))`, then the sum of their sixth powers has order
exactly `6m`: the leading coefficients are real sixth powers and cannot
cancel. The factorization of `F6^6-F5^6` shows that `F6∓F5` has order `6m`,
because every other factor is nonzero there. Unless `F6∓F5` is identically
zero, `d>=6m>=6`. If it is identically zero, positivity makes the whole
map trivial. This proves impossibility of nontrivial real contact for
`d<6`.

In particular an odd-degree nonzero real binary form `F6-F5` has a real
projective root, so nontrivial odd-degree parameterizations of degree
`1,3,5` are excluded. For even degrees `2,4`, that form can be definite;
the argument alone does not exclude such maps. No completed finite search
is undone by correcting a theorem claim beyond its proved scope.

## Néron–Severi and involution transfer: required evidence

The requested direct transfer has not been established from the inaccessible
original papers. Before treating a claimed elliptic pencil as an ESOP6
algorithm, provide all of the following concrete mathematical data:

- An explicit surface in the target or quotient, equations for the map to
  the base, and verification that the generic fibre is geometrically
  integral of genus one.
- A rational divisor class/pencil and its descent to `Q`, not merely an
  invariant lattice computation over an algebraic closure. A pencil over an
  extension may descend to a nonsplit base, and an invariant divisor class
  alone does not supply the desired rational section.
- Either a `Q(t)` section or a torsor description with a demonstrated
  rational specialization. A genus-one curve is not automatically an
  elliptic curve with a chosen rational origin.
- The coordinate-lifting and discriminant conditions on that fibre, with
  their ramification and real solubility checked. An elliptic point on a
  quotient need not lift to the sextic.

These are mathematical dependencies, not new computational search ranges.
The strongest new certified conclusion of this literature lane is the
complete rejection of the simultaneous-quadratic Gaussian family above.
