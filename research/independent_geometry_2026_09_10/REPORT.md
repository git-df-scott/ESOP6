# Full-fourfold centrally symmetric conic obstruction

Date: 2026-09-10. Status: **No ESOP6 counterexample found.**

This independent lane investigates a boundary-free construction on the full
six-coordinate variety, rather than the repeated-coordinate surface or the
degree-six boundary ansatz. It establishes an exact obstruction for one
natural family. No claim of novelty relative to all literature is made.

## Construction under examination

Let q(u,v) be a nondegenerate homogeneous binary quadratic over Q, and seek
five rational linear forms L_i(u,v) satisfying the polynomial identity

    sum_{i=1}^5 L_i(u,v)^6 = q(u,v)^3.                 (1)

If the conic q(u,v)=w^2 has a rational point, it can be parametrized over Q.
Substitution into (1) would give a rational curve on the ESOP6 fourfold;
a specialization with all coordinates nonzero would yield a counterexample.
This includes rational isometric-circle constructions using five linear
coordinate functions on a conic, without prescribing a known boundary seed.

**Theorem. No identity (1) exists for nondegenerate q.**
Indeed, the obstruction already holds over Q_7.

## Local lemma

Suppose a_1,...,a_n,z are in Q_7, 1<=n<=5, not all a_i zero, and

    a_1^6 + ... + a_n^6 = z^3.

Then z is a nonzero square in Q_7.

Proof. Let r=min_i v_7(a_i). Write a_i=7^r b_i. At least one b_i is a unit.
Every sixth power of a unit is 1 modulo 7, so the sum of the b_i^6 reduces
to k modulo 7, where 1<=k<=n<=5 counts the units. Consequently

    v_7(sum a_i^6)=6r, hence v_7(z)=2r.

Write z=7^(2r) z_0 with z_0 a unit. Reduction gives z_0^3=k modulo 7.
Nonzero cubes modulo 7 are 1 and 6. Since k is between 1 and 5, k=1.
Thus z_0 modulo 7 is one of 1,2,4, exactly the nonzero square residues.
Hensel's lemma, with derivative 2s a unit, lifts a square root of z_0.
Together with the even valuation this proves the claim.

There is also no nonzero tuple with sum a_i^6=0 over Q_7 for n<=5:
the same minimum-valuation reduction would give k=0 modulo 7, impossible.

## Proof of the theorem

A polynomial identity over Q remains valid at every point of Q_7^2. If
q(u,v) is nonzero, (1) ensures not all the L_i(u,v) vanish, so the lemma
forces q(u,v) to be a Q_7 square. Thus q would take only square or zero
values on Q_7^2.

Over a field of characteristic different from 2, a nondegenerate binary
quadratic diagonalizes by an invertible linear change of variables:

    q = A U^2 + B V^2, with A,B nonzero in Q_7.

The values at (1,0) and (0,1) force A and B to be squares. Write
A=alpha^2, B=beta^2. Taking U=1/alpha and V=2/beta gives q=5. But 5 is a
nonsquare unit modulo 7, hence a nonsquare in Q_7. Contradiction.

This proof neither assumes repeated coordinates nor bounds coefficients.
It excludes the entire specified identity family, not a finite search box.

## Scope and what survives

The proof excludes centrally symmetric constructions in which all five
summand coordinates depend linearly on u,v and have no w component, while
the denominator coordinate is w on q(u,v)=w^2.

It does **not** exclude general conics in the full fourfold. A general
construction can have

    x_i = a_i u + b_i v + c_i w,    q(u,v)=w^2.

After replacing w^2 by q, the identity has coupled even and odd parts in w.
The local argument above does not reduce that general system to (1).
In particular, this report does not claim that all full-fourfold conics
are absent, or that all rational curves are boundary-attached.

It also does not exclude point constructions on higher-genus curves,
rational curves of higher degree, or isolated rational ESOP6 points.
No existence claim is made for any surviving family.

## An additional necessary low-degree reduction constraint

Let homogeneous forms X_1,...,X_5,W of degree d<=3 over Z_7 satisfy
sum X_i^6=W^6, with at least one coefficient a 7-adic unit. Their reductions
modulo 7 satisfy the same identity.

At each of the 8 points of P^1(F_7), either all six evaluations vanish, or
W is nonzero and exactly one X_i is nonzero: the number of nonzero summands
is between 0 and 5, and must equal W^6, which is 0 or 1.

Every nonzero degree-d binary form has at most d distinct projective roots,
so its nonzero-evaluation support has size at least 8-d. Two nonzero X_i
would have disjoint supports of size at least 8-d each; this is impossible
for d<=3. W cannot be the only nonzero reduced form, since its evaluations
would all vanish. Thus exactly one X_j remains nonzero as a polynomial.
The polynomial identity gives W^6=X_j^6, hence W=zeta X_j for some
zeta in F_7^* (factor in the integral domain F_7[u,v]).

Therefore any full-fourfold conic represented by primitive integral binary
forms must reduce at 7 to a single boundary coordinate, with the other
four summand forms zero. This is a necessary reduction constraint only;
it is not a proof of nonexistence or a license to assume a rational
boundary preimage in characteristic zero.

## Verification

Run:

    python research/independent_geometry_2026_09_10/verify.py

The standard-library-only script checks all 16,807 ordered five-tuples
modulo 7, the residue square/cube conditions, and all 36 nonzero diagonal
coefficient pairs modulo 7. Results are in verification.json. The report's
valuation and diagonalization arguments supply the proof beyond these
finite regression checks. No integer ESOP6 search is represented as run.

## Explicit real conic: real existence is not the missing step

The parent supplied the following real curve, checked here by exact SymPy
expansion in verify_real_conic.py. Put

    U=s^2-t^2, V=2st, W=s^2+t^2, k=(63/80)^(1/6).

Then

    x=(kU, kV, k(U+V)/sqrt(2), k(U-V)/sqrt(2), W/2)

satisfies sum x_i^6=W^6 identically, since

    U^6+V^6+((U+V)/sqrt(2))^6+((U-V)/sqrt(2))^6
       = (5/4)(U^2+V^2)^3,
    U^2+V^2=W^2,
    (63/80)(5/4)+(1/2)^6=1.

Generic real parameters make all coordinates nonzero; independent signs
can make the summands positive. This curve has **no rational projective
point**, not merely no known rational parameter specialization. Indeed its
coordinates satisfy

    x_3+x_4=sqrt(2) x_1,  x_3-x_4=sqrt(2) x_2.

If all homogeneous coordinates were rational after a common scaling,
these relations force x_1=x_2=x_3=x_4=0. Then x_5=W/2 and the sextic equation
force W=0, so there is no projective point. The real construction therefore
cannot be converted to ESOP6 by rational approximation or a lucky parameter.

## Executed full-conic necessary tangent sieve

Let a primitive quadratic binary-form map have unique surviving summand Z
modulo 7, as proved above. Absorb the sixth root of unity relating it to W
into Z over Q_7. Let a>=1 be the minimum coefficient valuation of the other
four forms, and write those forms as 7^a U_i. At least one U_i has unit
coefficient. The sum of their sixth powers is not the zero polynomial
modulo 7: otherwise evaluation at all eight points of P^1(F_7) would force
each quadratic U_i to vanish identically. Hence this sum has Gauss valuation
zero. Factoring W^6-Z^6 shows W-Z has Gauss valuation exactly 6a; the other
five factors have Gauss valuation zero because their reductions are
nonzero constant multiples of Z. Write W-Z=7^(6a)H. Expansion gives

    sum U_i^6 = 6 Z^5 H + 15*7^(6a) Z^4 H^2 + ... .

Consequently the necessary equation modulo 7, 49, and 343 is

    sum_{i=1}^4 U_i^6 = 6 Z^5 H.                    (2)

The executable calls this surviving form W rather than Z.

The exhaustive modulo-7 sieve uses all 57 nonzero projective coefficient
classes of binary quadratics, plus zero. Multiplying a nonzero form by an
F_7 unit does not change its sixth power. For Q_7 lifting this normalization
is legitimate because each unit residue lifts to a sixth root of unity;
it is not an assertion of independent rational scaling symmetry.

All 1,711 unordered pairs were joined against each other (1,464,616
unordered pair-pair presentations). The targets use every nonzero H in
F_7^3 and all 57 normalized W, giving 19,494 distinct right sides.
Duplicates from different pair decompositions were removed.

Results:

- 228 proportional solutions, in which every nonzero U_i is proportional W.
- **84 nonproportional solutions**, all with four nonzero U_i.
- Every nonproportional solution has H proportional W and irreducible W;
  the discriminants 3,5,6 each occur 28 times.
- The coefficient Jacobian of (2), with 13 equations and 18 variables, has
  rank 11 over F_7 for all 84 seeds.
- **All 84 seeds lift to modulo 49.** This is checked by rank of the
  augmented linear system and by direct integer expansion of a saved lift.
- Taking free correction variables to zero picks one particular modulo-49
  lift per seed. **63 of those selected branches lift to modulo 343**;
  21 selected branches fail. These failures do not exclude any full seed
  family, because the alternative modulo-49 lifts were not exhausted.

One seed is

    U_1=1+t^2, U_2=1+2t^2,
    U_3=1+t+4t^2, U_4=1-t+4t^2,
    W=1+4t^2, H=3+5t^2=3W (mod 7).

All 84 seeds and the selected exact coefficient lifts are saved in
conic_tangent_seeds.json. The complete repeatable calculation is
conic_tangent_sieve.py, using NumPy only for vectorized modular matching.
The output summary is conic_tangent_sieve.json.

**Interpretation:** this is a concrete surviving local search space for
general full-fourfold conics. It is not a Q_7 existence proof, rational
reconstruction, or counterexample. The higher terms omitted from (2) first
enter at 7^(6a), so even indefinitely lifting (2) alone would not verify the
original conic system. The proportional solutions also remain relevant
as possible deeper degenerations and have not been discarded globally.

## Canonical continuation using the original conic equation

A further bounded run continued the 63 saved modulo-343 branches using the
**original** equation, with Xi=7U_i and the dominant summand W-7^6 H:

    F = sum U_i^6 - 6HW^5 + 15*7^6 H^2 W^4
        -20*7^12 H^3 W^3 +15*7^18 H^4 W^2
        -6*7^24 H^5 W +7^30 H^6 = 0.

This corrects the truncation before it affects the next digits. At each
step the free correction variables were fixed to zero in the specified
row-reduced linear system. Surviving branch counts were:

| Modulus | Canonical surviving branches |
|---|---:|
| 7^3 | 63 |
| 7^4 | 12 |
| 7^5 | 10 |
| 7^6 through 7^12 | 10 at every step |

Every retained lift was checked by exact integer coefficient expansion of F
at its modulus. The final modulus is 13,841,287,201.

Bounded rational reconstruction was then attempted coefficient by
coefficient with denominator at most 10,000 and numerator absolute value
at most floor((M-1)/20,000), so the standard product bound guarantees
uniqueness within that rectangle. Reconstructed coefficients, when complete,
were substituted into the original polynomial identity using exact rational
arithmetic. **No exact rational identity was obtained.** Larger-height
reconstructions and alternative choices of free lifting variables remain
untested. The ten modular branches are not rational-point certificates.

Executable: continue_original_conic.py. Full exact residues, failed-branch
levels, and reconstruction attempts: original_conic_lifts.json.
