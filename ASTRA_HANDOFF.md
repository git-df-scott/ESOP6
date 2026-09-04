# Astra handoff — after the September 4 direct strike

## Binding current state

- No positive ESOP6 solution was found in the direct strike.
- No rational curve or positive rational point on `2X^6+2Y^6+Z^6=W^6`
  was found. Rational points on its cubic quotient are not solutions.
- Every requested preparation gate passed from base commit `485390a`.
- The historical 4.3M class-1 frontier remains historically attested,
  with 55,684 candidate pairs independently regenerated; lost production
  decompositions were not rerun. Its qualification is unchanged.
- caseA3, its memory correction, NB=8 control timing, and the inherited
  four-term residue/valuation filters retain their audited scope.

Read the new work first:

```text
ASTRA_DIRECT_STRIKE.md
GEOMETRIC_STRIKE.md
RATIONAL_POINT_LEDGER.md
results/astra_direct_2026_09_04/manifest.json
```

The preceding preparation documents remain available for provenance.

## What the direct attack resolved

The natural cubic-quotient conic pencil generically has genus-nine square
lifts. The rational exceptional member lambda=2 has no nonzero real lift.
A global genus-one fibration on the smooth sextic surface is ruled out by
adjunction; an isolated curve or a quotient construction is a different
question and must preserve all lift conditions.

Reduced parametrization degrees 1,2,3 are excluded. General degree 4 was not
classified. Every nonconstant curve through `[0:0:1:+/-1]` needs parameter
degree at least six. The exactly diagonal X=+/-Y restriction is closed by
Fermat's theorem for cubes.

For a primitive positive point and
`lambda=(W^2-Z^2)/(X^2+Y^2)`, with `r_p=min(v_p(X),v_p(Y))>=1`:

```text
v2(lambda)=4r2+1, v3(lambda)=4r3-1,
v7(lambda)=4r7 OR -2r7.
```

The negative 7-adic branch is essential. lambda=1/2 is globally closed by
these necessary restrictions. Four selected admissible conics were tested
over explicitly bounded parameter boxes, without a lift; they are not
globally closed.

A table-free exact divisor oracle tested 10,382 sparse CRT targets, with
factor primes proved and all dispositions independently replayed. This
does not extend any continuous height frontier or improve the unrestricted
six-variable exponent.

## Single next direct construction

Construct an exact identity in the family

```text
X=t*A5(t), Y=t*B5(t), W=Z+(2/3)*t^6,
A5(0)=B5(0)=Z(0)=1,
deg(A5),deg(B5)<=5, deg(Z)<=6.
```

Keep A5 and B5 independent. Divide the defining identity by t^6 before
coefficient work; exploit order-six contact and projective/parameter
normalization. A successful identity has positive rational specializations
for sufficiently small rational t>0. This degree increase is only in the
boundary branch proved empty through degree five; it is not a claim that
general quartics have been classified.

Do not impose the additionally symmetric family
`W=1+t^6/3, Z=1-t^6/3, X=t*A4(t), Y=t*A4(-t)`.
Its leading coefficient requires `a4^6=1/81`, impossible over Q.

Do not substitute larger empty boxes or a larger historical frontier for
this construction. No generic high-degree Groebner sweep is justified.

## Reproduction and solution certificate

Cheap new checks:

```bash
make direct-checks
```

The inherited gate, already completed in the recorded strike, is:

```bash
make clean && make -j"$(nproc)" && make validate control equiv frontier-audit prototypes && make differential
```

For an actual six-positive-integer candidate, stop searching immediately:

```bash
python3 tools/verify_esop6.py a b c d e f
node tools/verify_esop6.mjs a b c d e f
```

Both use arbitrary-precision integers and print full powers, exact sides,
positivity, gcd, and the recomputed normalized equality. Commit the exact
candidate and both outputs before making a solution claim.
