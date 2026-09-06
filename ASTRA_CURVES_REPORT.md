# ESOP6 full-fourfold curve strike — September 5–6, 2026

**No ESOP6 counterexample. No curve identity over Q was obtained.**

The plane-cubic lane is excluded over R by an exact odd-degree argument.
Five real algebraic coefficient solutions of the rational-quartic systems
were certified by exact rational contraction inequalities. They give genuine
quartic embeddings, but their coefficients have **not** been proved rational.
The elliptic-quartic and global conic-classification lanes remain open.

Branch: `astra/curves-elliptic`.
Base: `d3c30ce`, the latest inherited multiplicity/square/conic branch, including
`60d0243` and the two original strikes. Main was still at `36e03d4` at checkout.
The result directory retains the September 5 start date; work crossed into
September 6 UTC. `ASTRA_PROMPT_CURVES.md` was absent from the remote and was
created from the supplied assignment. No supplied credential was saved.

No integer sweeps, frontier extension, repeated-coordinate surface search,
or duplicate numerical conic search was run. The user-reported 4,310,000
frontier was not audited or changed. Existing source engines were untouched.

## 1. Results and exact scope

| System | Unknowns / equations | Work completed | Certified result |
|---|---:|---|---|
| Plane sextic = cubic times cubic | 28 / 28 in one affine chart | Exact ideal exported; global real obstruction proved | No real or Q-defined plane cubic on X; no complex enumeration |
| Elliptic quartic: F = Q1 A4 + Q2 B4 in a P3 | 84 / 84, with 10 cofactor syzygies removed | 32 starts: 24 real, 8 complex; 250 evaluations per start; refinement and rational reconstruction | No exact Q curve; full classification not completed |
| Single-involution rational quartic | 14 / 13 after two normalizations | 48 starts: 36 real, 12 complex; 400 evaluations per start | Two real algebraic roots certified after rational transverse sections; exact Jacobian rank 13 |
| Double-involution rational quartic | 8 / 7 after scalar normalization | Same 48-start allocation | Three real algebraic roots certified after rational transverse sections; exact Jacobian rank 7 |
| Full conic plane-section ideal | 29 / 28 in a scale chart | Exact ideal; exact number-field control; standard-plane strata over Q | No complete component decomposition or component genera |
| Degree-four coefficient identity over F7 | 30 coefficients / 25 equations | Complete support reduction, 2,801 projective forms, 35 complement-pair checks | Every tuple is an axis tuple or zero; independent squarefree-factor replay passes |

The **rational candidate count is zero** in the sense of exact verified
identities. This does not assert that every refined algebraic root is
irrational, nor exclude rational reconstructions beyond those attempted.
There is no sextuple to send through the two ESOP6 verifiers. No Jacobian,
rank, or generator computation was triggered: no eligible genus-one curve
over Q was established.

## 2. Two necessary corrections to the construction brief

### Only the sixth coordinate must be definite

For a morphism P1_R -> X represented by basepoint-free homogeneous forms,
a real zero of P6 forces a common real zero of all coordinates. Thus P6 is
definite and the degree is even. A zero of a *left* coordinate does not force
the other coordinates to vanish. It is enough that none of the six forms
vanishes at the specialization used to produce a candidate.

The exact conic in `ASTRA_CONIC_CONTROL.md` has this property, and so do
all five certified real quartic roots: P6 is definite, but at least one left
coordinate changes sign. Requiring every coordinate to be definite removes
these valid real curves. Signs at a nonzero specialization can be removed
by taking absolute values because the exponent is even.

### Expected dimension is not arithmetic existence

The virtual count is (1-g)(4-3) for this Calabi–Yau fourfold. It is one for
genus zero and zero for genus one, consistent with
[Klemm–Pandharipande](https://arxiv.org/abs/math/0702189).
It neither guarantees real points nor definition over Q; an isolated curve
may be defined over a number field. The logarithmic counting heuristic also
does not prove that any counterexample must lie on a rational or elliptic
curve. Construction is the assigned strategy, not a proved necessity.

## 3. Lane 1a: all real plane cubics are excluded

On the hyperplane H: x6=0 the real defining equation becomes

    x1^6+x2^6+x3^6+x4^6+x5^6=0.

It has no real projective point. Suppose C is a real projective curve of
odd degree contained in X. If C is not contained in H, its intersection
with H is a conjugation-invariant zero-cycle of odd degree. Nonreal points
come in conjugate pairs, so at least one intersection point is real, a
contradiction. If C is contained in H, an odd-degree real linear section
of C likewise supplies a real point, again a contradiction.

In particular **no real plane cubic, smooth or otherwise of pure degree
three, lies on X**. This excludes plane elliptic cubics over Q. It does not
exclude complex cubics or even-degree genus-one curves.

The requested incidence ideal is still supplied for independent use.
On the plane chart

    x1=u, x2=v, x6=w,
    (x3,x4,x5)=(L3(u,v,w),L4(u,v,w),L5(u,v,w)),

there are nine plane parameters. Give the first cubic nine parameters by
fixing its w^3 coefficient to one; the second cubic has ten coefficients.
Equating all 28 ternary sextic coefficients gives **28 equations in 28
unknowns**. This is `plane_cubic_ideal.json`. No many-start search was spent
on a real lane already excluded by proof. Its complex solution set was not
enumerated; planes contained in X already produce excess complex families.

## 4. Lane 1b: elliptic quartic system and bounded run

Let z=(z0,z1,z2,z3). The tested P3 chart is

    (x1,x2,x3,x4,x5,x6)=(z0,z1,z2,z3,L(z),M(z)),

with eight plane parameters. Use the open chart of the quadric-pencil
Grassmannian where

    Q1=z0^2 + eight free quadratic coefficients,
    Q2=z1^2 + eight free quadratic coefficients,

and neither free part contains z0^2 or z1^2. This contributes sixteen
parameters. Impose

    z0^6+z1^6+z2^6+z3^6+L^6-M^6 = Q1*A4 + Q2*B4.

A4 has 35 coefficients. B4 has 25: its monomials divisible by z0^2 are
removed by the syzygy (A4,B4) -> (A4+Q2*C2,B4-Q1*C2). Because Q1 is monic
with leading monomial z0^2, division eliminates those ten coefficients.
Thus there are **8+16+35+25=84 unknowns and 84 coefficient equations**.
On the regular complete-intersection locus this removes exactly the ten
cofactor syzygies. The geometric count is 8+16-24=0, since h0(E,O_E(6))=24.

This is a full ambient P3 chart, not a repeated-coordinate surface slice.
Other P3 and pencil charts were not run. Smoothness, real points, definition
over Q, and a rational point on the actual genus-one curve remain separate
gates. In particular a positive Jacobian rank would not by itself provide a
point on a genus-one torsor.

With seed 20260905, 32 starts reached six distinct *numerical coefficient
endpoints* below the weighted residual threshold 1e-9. Four used real
coefficients, two complex coefficients. Numerical Jacobian ranks at these
six endpoints, using singular-value cutoff 1e-8*max(sigma_max,1), were:

    {84: 2, 83: 1, 81: 2, 77: 1}.

These are not six certified elliptic curves. One complex endpoint exhausted
its budget; several are singular or lie in numerically singular incidence
loci. A real coefficient tuple need not give any real curve points.

Five endpoints were refined with 80-digit arithmetic. Four reached residuals
below 1e-65; real endpoint 20 stalled near 1.5e-33. These are diagnostics only.
For four real endpoints, coefficient-wise continued-fraction reconstruction
at denominator caps 100, 10,000 and 1,000,000 gave twelve rational tuples;
**all twelve failed an independent exact 84-coefficient identity check**.
The three reconstructions of the nonreal endpoint were ineligible as Q
coefficients. No curve was passed to a Jacobian/rank routine. Sage, PARI,
Singular and Macaulay2 were unavailable on this host, but that was not the
arithmetic bottleneck: no Q curve had survived the identity gate.

The full ideal, raw endpoints, ranks, refinement histories and exact nonzero
reconstruction residuals are retained. Completeness of the complex roots,
all real elliptic quartics, and the arithmetic lane is **not claimed**.

## 5. Lane 2: exact real quartic root certificates

Use coefficient order s^4,s^3*t,s^2*t^2,s*t^3,t^4.

For the single symmetry, take P2(s,t)=P1(s,-t), P4(s,t)=P3(s,-t), with
P5 and P6 even. There are sixteen coefficients before normalization. Set
P1[s^4]=1 using common scale, and P1[s^3*t]=1 using t scaling on the chart
where both coefficients are nonzero. Fourteen parameters remain, with
thirteen even coefficients of the degree-24 identity to vanish.

For the double symmetry, also set P3(s,t)=P1(t,s), and take P5 and P6 even
and reciprocal. There are nine coefficients before common scale; setting
P1[s^4]=1 leaves eight. The residual is even and reciprocal, so coefficients
at degrees 0,2,4,6,8,10,12 give seven equations. No further continuous PGL2
normalization was assumed. Special omitted normalization charts remain open.

The bounded search endpoint accounting is:

| System | Small-residual endpoints | Real coefficient endpoints | Numerical ranks at accepted endpoints | Complete complex root count |
|---|---:|---:|---|---|
| Single | 2 | 2 | rank 13: 2 | Unknown; positive-dimensional system |
| Double | 8 | 3 | rank 7: 8 | Unknown; positive-dimensional system |

To turn five of these into exact **real** statements, one free coefficient
was fixed to a nearby rational value, and each resulting square system was
refined. The actual fixed values and root certificates are:

| System / start | Fixed parameter | Certified Jacobian rank | P6 definite | All six definite | Q coefficients verified |
|---|---|---:|---|---|---|
| Single / 18 | x9=-1173/904 | 13 | Yes | No | No |
| Single / 25 | x5=-129/158 | 13 | Yes | No | No |
| Double / 18 | x5=-853/921 | 7 | Yes | No | No |
| Double / 27 | x4=-1246/705 | 7 | Yes | No | No |
| Double / 34 | x5=-3/5 | 7 | Yes | No | No |

Variable indexing is defined by the exported ideals. Distinct coefficient
boxes are not a claim of distinct geometric components or inequivalence
under parameter changes and target symmetries.

### Exact proof, not a floating-point solution announcement

Each certificate contains a rational center c, a rational approximate inverse
B of the square Jacobian, the exact fixed coefficient, and radius r=10^-25.
For the *unweighted exact polynomial equations* F, the verifier proves

    eta = ||B F(c)||_infinity <= r/2,
    K = ||I-B J(c)||_infinity + ||B||_infinity*n^2*H*r < 1/2.

Here H=720*(5M)^4 bounds every second partial throughout the cube, where
M bounds the absolute value of every coordinate-polynomial coefficient.
Indeed each second derivative has the form 30*sum(sign_i*P_i^4*dP_i*dP_i),
and each parameter derivative has polynomial coefficient l1 norm at most
two. All bounds and comparisons are computed with exact rational arithmetic.

Then T(x)=x-BF(x) maps the closed cube into itself and is a contraction.
It has a unique real fixed point, which is a zero of the exact equations
because B is invertible. Thus these are isolated real algebraic roots of
the respective rational transverse systems. The largest certified K is
less than 7.14e-14; the largest eta is less than 3.65e-55. The certificates,
not these rounded summaries, carry the proof.

The verifier also establishes throughout each cube:

- a 5-by-5 coefficient minor is invertible, so the six forms span all of
  H0(P1,O(4)); the resulting map is a basepoint-free quartic embedding;
- all three even coefficients of P6 are positive, proving definiteness;
- every coordinate is nonzero at (s,t)=(1,0);
- a left coordinate has opposite signs at two retained real parameter
  values, proving the all-definite filter would reject it.

A second symbolic implementation independently reconstructs and differentiates
the equations and agrees exactly on eta and the inverse-Jacobian error for
all five certificates. The Jacobian submatrix stays invertible, proving
ranks 13 and 7, respectively. The original systems consequently have local
one-dimensional real solution loci at these points.

These algebraic-root certificates are **not identities with coefficients
in Q[t]** and are not rational candidates. Each of the five refined roots
was reconstructed at four denominator caps: 100, 10,000, 1,000,000 and
1,000,000,000. All twenty rational reconstructions failed exact coefficient
substitution. No claim of irrationality of the exact roots follows.

## 6. A complete degree-four coefficient restriction at 7

Let six homogeneous quartics over F7 satisfy the sextic identity. At any
point of P1(F7), a sixth power is zero or one. Since there are five left
summands, at most one can be nonzero at each of the eight points.

A nonzero quartic has at most four zeros, so its support has at least four
points. There can therefore be at most two nonzero left quartics. If there
are two, each has exactly four support points and the supports are
complementary. Each prescribed set of four roots determines its quartic
up to scalar; nonzero scalar sixth powers are one. This leaves exactly
35 unordered complementary pairs among 70 four-point supports.

The certificate enumerates all 2,801 nonzero projective quartics and their
sixth powers. None of the 35 sums is a sixth power of a quartic. An independent
SymPy squarefree factorization of each sum finds a factor multiplicity not
divisible by six, reproducing every rejection without the lookup table.

Only one left polynomial can survive, and equality of its sixth power to
the right polynomial forces proportionality by a nonzero F7 constant.
If the right polynomial vanishes identically, all left polynomials do too.

Hence **every primitive integral coefficient model of a rational quartic
reduces to a constant axis map at 7**. This is a necessary bad-reduction
condition on the model. It does not exclude rational quartics: arbitrarily
bad reduction and denominators are allowed. In the two symmetry charts with
P1[s^4]=1, an all-7-integral coefficient solution is impossible, since
symmetry makes at least two left polynomials nonzero. Any rational solution
in these normalizations must have a denominator divisible by 7.

## 7. Lane 3: exact conic ideal, real control, and standard-plane strata

On x1=u,x2=v,x6=w, with the other three coordinates arbitrary linear forms,
write

    F|P = q2*r4,   coefficient(q2,w^2)=-1.

There are nine plane parameters, five conic parameters and fifteen quartic
cofactor parameters: **29 unknowns and 28 equations**. This is the exported
`conic_plane_ideal.json`; a smooth conic also needs det(q2)!=0. A global
scheme computation must handle the other coefficient charts and saturation.

The exact real control is

    [sqrt(2)*u : sqrt(2)*v : u+v : u-v : w : 11^(1/6)*w],
    u^2+v^2=w^2.

Its restricted sextic factors exactly as

    10*(u^2+v^2-w^2)*((u^2+v^2)^2+(u^2+v^2)*w^2+w^4).

Its parameter Jacobian has rank 11 of 13 rows; the plane-factorization
Jacobian has rank 26 of 28 rows. These ranks are exact over the specified
number field. After division of all coordinates by sqrt(2), it meets the
other session's D normalization and still has exact rank 11. The explicit
seed and its lack of rational points are in `ASTRA_CONIC_CONTROL.md`.

The full Fermat conic scheme is **not** the smooth curve that arises for a
general sextic fourfold; the general statement in
[Ciliberto–Zaidenberg](https://arxiv.org/abs/1910.11423) does not transfer to
this special hypersurface. To see the obstruction directly, pair the six
coordinates into three pairs. In each positive-positive pair the ratio is
a sixth root of -1; in the pair containing x6 it is a sixth root of +1.
There are 15 pairings and 6^3 choices, giving **3,240 distinct standard
complex planes**. Every such plane contains a P5 of conics, with the smooth
ones forming an open subset. This alone disproves a global dimension-one
premise. It does not assert these P5 strata are all irreducible components.

For one pairing the exact plane-parameter ideal is

    (a^6+1, b^6+1, c^6-1).

Its reduced Q algebra has four degree-two closed points with residue field
Q(i), and 52 degree-four closed points with residue field Q(zeta_12).
This was computed by the simultaneous action of {1,5,7,11} on the triples
of exponents (odd,odd,even) modulo 12. Across all pairings this gives
**60 degree-two and 780 degree-four standard plane points over Q**, accounting
for all 3,240 geometric standard planes. None is real. The corresponding
smooth conics have their unique span in that nonreal plane and are not
real-defined conics. This is an exact classification of the stated standard
strata, not all planes or all conics.

A full Groebner decomposition of the remaining conic incidence scheme,
its Q-defined components and any one-dimensional component genera was **not
completed**. Asking for a genus for every component before separating these
higher-dimensional strata is not well-posed.

## 8. Reproduction and remaining construction

Dependencies used: Python, NumPy, SciPy, SymPy and mpmath. Exact certificates
use rational arithmetic; numerical packages supply search seeds only.
All commands are from the repository root:

    python3 tools/astra_curves/models.py
    python3 tools/astra_curves/exact.py
    OPENBLAS_NUM_THREADS=1 python3 tools/astra_curves/search.py elliptic --starts 32 --real-starts 24 --max-nfev 250
    OPENBLAS_NUM_THREADS=1 python3 tools/astra_curves/search.py single --starts 48 --real-starts 36 --max-nfev 400
    OPENBLAS_NUM_THREADS=1 python3 tools/astra_curves/search.py double --starts 48 --real-starts 36 --max-nfev 400
    OPENBLAS_NUM_THREADS=1 python3 tools/astra_curves/refine.py
    OPENBLAS_NUM_THREADS=1 python3 tools/astra_curves/certify_real.py
    OPENBLAS_NUM_THREADS=1 python3 tools/astra_curves/verify.py

Rerunning the numerical searches can change endpoints across library versions;
the retained rational certificates replay deterministically without rerunning
those searches. File hashes and the environment are in `manifest.json`.

The concrete arithmetic frontier is now the **Q-point problem on the quartic
parameter loci**, with bad reduction at 7 required. The exact rational slice
x5=-3/5 in the double-symmetry system supplies one especially simple real
root isolation to use in elimination. It does not prove the root is rational
or make that rational section privileged for finding Q points. Arithmetic
elimination/descent and other sections are still required.

Still unfinished: a Q-coefficient curve or sextuple; classification of all
elliptic quartics, including missing charts; a Q genus-one curve with rational
points; rationality or field degrees of the certified quartic roots; and the
remaining global Fano-scheme components and genera. No background search is
claimed to continue after this work block.
