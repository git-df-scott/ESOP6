# Exact quotient descent and the Pythagorean seed lane

2026-09-10. No rational point on the original sextic and no counterexample
were found. The results here are exact identities, six entire rational
coordinate slices excluded by local arithmetic, a complete local-solubility
criterion for one Pythagorean subfamily, and a bounded enumeration of that
subfamily. No numerical lottery was repeated.

Reproduce the retained certificate with standard-library Python:

```bash
python research/quotient_descent_2026_09_10/exact_quotient_audit.py
python research/quotient_descent_2026_09_10/independent_finite_verify.py
```

The run takes about one second. `exact_results.json` retains every prime
witness, residue set, rejected seed slice, and the seven intermediate
Pythagorean survivors. It records a disposition digest for all 202,861
normalized parameters. Five symbolic identities are checked coefficient by
coefficient with exact sparse-polynomial arithmetic. Both seed equations
are independently evaluated in their quadratic fields using pairs of integers.
The separate finite verifier imports no generating code; it directly
enumerates finite-field coordinates, checks the sixth-power unit images
modulo 256 and 729, and rechecks every retained local gate. It passed.

## 1. Quotient arithmetic and the genuine cover conditions

Put `e1=x1+x2`, `e2=x1*x2`. Newton's identities give

\[
x_1^6+x_2^6=e_1^6-6e_1^4e_2+9e_1^2e_2^2-2e_2^3.
\]

Write `A=e1/2`, `Omega=A^2-e2`. Then

\[
x_{1,2}=A\pm\sqrt\Omega,\qquad
x_1^6+x_2^6=2\Phi(A,\Omega),
\]

\[
\Phi=A^6+15A^4\Omega+15A^2\Omega^2+\Omega^3.
\]

A quotient point lifts over Q precisely when `Omega` is a rational square
(equivalently `e1^2-4e2=4*Omega` is a rational square). For a usable lift,
both resulting coordinates must be nonzero; their signs can then be
removed because the exponent is even. A negative discriminant has no real
lift. Positive discriminant is insufficient without an exact rational square.

The Gaussian seed has `(e1,e2)=(14,170)`, discriminant `-484`. The boundary
seed has `(18,330)`, discriminant `-996`. Neither is a rational sextic point.

There is an inexpensive additional elliptic compression, but it retains a
cube condition. Fix `A,x3,x5,x6`, let `H=x4`, and set

\[
C=x_6^6-x_3^6-x_5^6,\quad T=\Omega+5A^2,\quad x=-2T,\quad y=2H^3.
\]

Since `Phi=T^3-60A^4*T+176A^6`, the exact elliptic equation is

\[
\boxed{y^2=x^3-240A^4x+4C-1408A^6.}
\]

The inverse requires **both**

\[
y/2\in\mathbb Q^3,\qquad -x/2-5A^2\in\mathbb Q^2.
\]

The signs in the constant term matter. For the Gaussian Pythagorean
compression, `A=7`, `C=3*2040^2`, and the elliptic point is
`(x,y)=(-248,3456)`, giving `H=12` and `Omega=-121`.

Assume the displayed cubic in x has distinct roots, and `C-2A^6 != 0`.
Then `y/2` on the elliptic curve has three simple zeros and a pole of order
three at infinity. The cubic cover `H^3=y/2` ramifies at the three zeros
and not at infinity, so Riemann–Hurwitz gives genus four:
`2g-2=3*(2*1-2)+3*(3-1)=6`. Imposing the second condition gives the
double cover `w^2=Omega`; its six simple branch points occur at
`Omega=0`, `H^6=C-2A^6`. There are three points over infinity on the
genus-four curve, each with an even order-two pole of `Omega`; these do
not branch. Consequently the full cover has genus ten:
`2g-2=2*(2*4-2)+6=18`.

Thus elliptic arithmetic alone does not produce a rational fibration of
the sextic. Singular parameters and coincident branch points require
separate normalization and are not covered by these generic genus counts.

## 2. Six entire coordinate slices through the seeds are impossible

Keep `x6` fixed and keep any two of the three rational tail coordinates
`x3,x4,x5` fixed at their seed values. Allow the other tail coordinate and
both members of the conjugate pair to vary arbitrarily over Q. Every one
of these six slices is locally insoluble after asking for a rational lift.

Each would require a sum of three rational sixth powers equal to the
constant in this table.

| Seed | Free tail coordinate | Fixed tail values | Required sum | Obstruction |
|---|---:|---|---:|---|
| Gaussian, f=17 | x3 | 12,15 | 9,760,960 | unit 6 mod 7 |
| Gaussian, f=17 | x4 | 8,15 | 12,484,800 | unit 6 mod 7 |
| Gaussian, f=17 | x5 | 8,12 | 20,889,441 | unit 6 mod 7 |
| Boundary, f=22 | x3 | 18,0 | 79,367,680 | valuation 1 at 7 |
| Boundary, f=22 | x4 | 14,0 | 105,850,368 | valuation 3 at 3 |
| Boundary, f=22 | x5 | 14,18 | 71,838,144 | valuation 1 at 7 |

The constants and certificates are independently generated in
`exact_results.json`; the proof does not depend on numerical approximation.

For three rational numbers, put `r=min v_p(x_i)`. At p=7, the primitive
unit sixth powers are all 1 modulo 7. Their nonzero count is 1,2, or 3,
so no cancellation is possible:

\[
v_7(\textstyle\sum x_i^6)=6r,
\quad 7^{-6r}\textstyle\sum x_i^6\pmod7\in\{1,2,3\}.
\]

At p=3 the same argument modulo 9 shows that the valuation is `6r` if
one or two summands have minimum valuation, and `6r+1` if all three do.

For the Gaussian seed, all three required constants are 6 modulo 7,
which is impossible. For the boundary seed, the slices fixing `(18,0)`
and `(14,18)` have valuation one at 7. The remaining slice, fixing
`(14,0)`, has valuation three at 3. All are impossible even with arbitrary
rational denominators: the minimum-valuation argument accounts for them.

**Consequence:** a potentially useful quotient family through either seed
must vary at least two of the three normalized ratios
`x3/x6,x4/x6,x5/x6`. A search varying only `A,Omega` and one tail ratio
cannot work, at any height or degree. This is a certified pruning rule,
not a numerical absence claim about all conics through the seeds.

The literature worker also observed a separate real obstruction. The
Gaussian seed satisfies `sum_{i<=5} xi^2=x6^2`: the conjugate pair has
square sum `-144`, cancelled by `12^2`, and `8^2+15^2=17^2`. A quotient
family preserving this entire quadratic relation cannot have a real
lift with at least two nonzero left-hand coordinates. Indeed, for
`u_i=xi^2>=0`, `(sum u_i)^3>sum u_i^3` when at least two are positive.
Only the trivial one-term boundary survives equality. Preserving just
the Pythagorean relation on three tail coordinates is a different,
less restrictive construction considered next.

## 3. Pythagorean subfamily and complete local arithmetic

Set

\[
x_3=2mn,\quad x_5=m^2-n^2,\quad x_6=m^2+n^2,
\quad K=x_3x_5x_6=2mn(m^4-n^4).
\]

Then

\[
x_6^6-x_3^6-x_5^6=3K^2,
\quad x_1^6+x_2^6+x_4^6=3K^2.\tag{P}
\]

The Gaussian seed is `(m,n)=(4,1)`, `K=2040`. Its factorization
`7+11i=(1+3i)(4-i)` explains an available Gaussian parametrization,
but does not supply a square discriminant.

Here is a complete local-solubility test for (P), for any nonzero
rational K. It concerns the entire three-variable affine equation;
it does not assert that its particular fixed-A elliptic subcover has points.

* At **2** and **3**, (P) is soluble precisely when `v_p(K)` is a
  multiple of three.
* Let `S={7,31,67,79,139,223}`. At each p in S, write
  `K=p^(3r)*k` with k a p-adic unit. Solubility is equivalent to the
  existence of such an integer r and
  `3*k^2 in R_p+R_p+R_p`, where
  `R_p={a^6 mod p:a in F_p}`. The exact finite sets are retained.
* At every other finite prime, (P) is automatically soluble.
* At the real place, (P) has positive real points.

Proof at 2: after removing the common sixth power, one or three odd
summands give valuation zero; two give valuation one. Since the right
side has even valuation, the latter is impossible. Its normalized unit
is 3 modulo 8, so precisely three summands are odd, and `v_2(K)=3r`.
Conversely put two normalized variables equal to 1 and take a sixth root
of `3*k^2-2`, which is 1 modulo 8. The sixth-power units of Q2 are exactly
`1+8 Z2` (square units, followed by the invertible cube map).

At 3, normalized valuations are zero or one as above. The right-side
valuation is odd, so all three summands have the same minimum valuation
and `6r+1=1+2v_3(K)`. Conversely `3*k^2-2` is 1 modulo 9 and has a
sixth root in Q3. The sixth-power units are exactly `1+9 Z3`; for example,
apply Hensel to `((1+3t)^6-1)/9`, whose derivative is 2 modulo 3.

For the remaining primes consider the smooth projective plane sextic

\[
D:\ u^6+v^6+w^6=0.
\]

Its genus is ten in characteristic other than 2 and 3. The retained
exact enumeration of all primes through 397 shows that the only
anisotropic primes are S. Normalizing a nonzero coordinate to 1 loses
no projective point, because the equation is symmetric. For p>=401,
the Weil lower bound `#D(F_p)>=p+1-20 sqrt(p)>0` proves isotropy.
The bound is the standard one recalled and used, for example, in
Elkies–Howe–Kresch–Poonen–Wetherell–Zieve,
[*Curves of every genus with many points, II*, §2.1](https://arxiv.org/pdf/math/0208060).

At an anisotropic prime the minimum-valuation argument forces
`v_p(3K^2)=6r`, and the primitive residue must be in the displayed
sumset. These conditions are also sufficient: a nonzero residue
solution has a unit variable, so one-variable Hensel applies.

At an isotropic prime choose a common variable valuation r sufficiently
negative that `p^(-6r)*3K^2` is integral and divisible by p. Lift a
nonzero point of D by Hensel to a solution with this right side, and
rescale the variables by `p^r`. This use of denominators is essential.
It proves that no omitted prime can obstruct the affine rational problem.
If a chosen witness has a zero coordinate, a small perturbation of the
other free variable makes every coordinate nonzero before the same
Hensel step. Thus nonzero local points can be used throughout.

For a **primitive integral Pythagorean triple**, K is divisible by 2 and
3. The conditions therefore force both valuations to be at least three.
Also K must be divisible by 7: if all three Pythagorean coordinates
are units modulo 7, their squared product is 2 modulo 7, so `3K^2=6`,
which is excluded. Consequently `v_7(K)>=3` as well. This gives the
particularly cheap prerequisite `42^3 | K` with each of its 2,3,7
valuations a multiple of three. Exact valuations and the normalized
unit residue at 7 must still be checked.

## 4. Exhausted normalized parameter range

The domain was every coprime pair `1<=n<m<=1000` with opposite parity,
each giving one primitive Pythagorean triple. Swapping its two legs
does not change (P). This is a bounded arithmetic search on distinct
Pythagorean ratios, not an extension of the unrestricted integer frontier.

| Stage | Survivors |
|---|---:|
| Primitive normalized parameter pairs | 202,861 |
| 2-adic test | 58,019 |
| 3-adic test | 4,682 |
| 7-adic test | 7 |
| All places | 4 |

The four everywhere locally soluble slices are:

| m | n | x3 | x5 | x6 |
|---:|---:|---:|---:|---:|
| 428 | 85 | 72,760 | 175,959 | 190,409 |
| 883 | 540 | 953,640 | 488,089 | 1,071,289 |
| 964 | 343 | 661,304 | 811,647 | 1,046,945 |
| 995 | 652 | 1,297,480 | 564,921 | 1,415,129 |

The other three survivors of 2,3,7 are `(316,27)`, killed at 79;
`(540,197)`, killed at 31 (also 67); and `(999,716)`, killed at 31.
No global rational point search was performed on these four retained
slices, and local solubility does not imply such a point. In particular,
the first row being below an inherited integer frontier is no reason
to repeat that integer search. Rational denominators would change the
eventual integer height.

## 5. Scope and next construction

No rational family on Z through the two seeds was obtained. The exact
pruning results show why certain tempting low-dimensional families
cannot lift and supply a much smaller, rigorously specified alternative
Pythagorean slice list. The next useful quotient construction must move
at least two normalized tail ratios, avoid preserving the Gaussian
full-square identity, and retain every cube and square lift condition.
An elliptic point on the compression above is only a proposal until
both lifts are exact and all six reconstructed integers are positive.

This report does not use Lander–Parkin–Selfridge as a theorem. The night
report's implication from 'no rational root' to 'irreducible coordinate
form' is invalid without more hypotheses; a reducible quartic can have
two irreducible quadratic factors. Nor do numerical Newton starts or
failed recognition certify that all rational curves through a seed are
absent. Those historical overclaims are not assumptions of any result here.
