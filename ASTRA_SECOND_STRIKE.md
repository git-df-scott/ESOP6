# Astra second strike — boundary-contact-6

```text
ESOP6 SOLUTION FOUND: NO
POSITIVE RATIONAL SURFACE POINT: NO
EXACT BOUNDARY-CONTACT CURVE: NO
FORMAL RATIONAL BRANCH: YES
POLYNOMIAL ANSATZ: IMPOSSIBLE
RATIONAL-FUNCTION RELAXATION: LIVE
ELLIPTIC RESIDUAL: NO
```

“Formal rational branch” means coefficients in Q[[t]], not a rational
function in Q(t). No six-positive-integer candidate was produced.

## Result

The complete first-strike polynomial proposal

```text
X=t A5(t), Y=t B5(t), W=Z+(2/3)t^6,
A5(0)=B5(0)=Z(0)=1,
deg A5,deg B5<=5, deg Z<=6
```

is **impossible over Q**, including independent asymmetric A5,B5 and
arbitrary rational denominators in their coefficients. The obstruction
does not close all rational boundary-contact curves, unequal leading
tangent directions, or the positive rational-point problem.

The proof is in [BOUNDARY_CONTACT_6.md](BOUNDARY_CONTACT_6.md). Its decisive
steps are:

1. Center W and Z. The first six coefficient equations solve Z triangularly;
   after a legitimate Möbius normalization, nine effective coefficients and
   24 residual equations remain.
2. Invert the parameter. A 3-adic Gauss/Newton-polygon argument handles
   every possible coefficient denominator, reducing all surviving cases to
   integral monic quintics and a sextic with a unit constant coefficient.
3. Modulo 9, their only possible residue shape, up to exchange, is
   `U=1+2t^2+t^4, V=t+2t^3+t^5 (mod 3)`. Its lifts retain all nine
   independent coefficients.
4. Exact linear lifting has rank four. The surviving-prefix counts are
   `1 -> 81 -> 162 -> 0` at equation moduli `27,81,243,729`. The final
   stage rejects every one of 39,366 prefixes. Altogether 59,293 prefixes
   are checked, including the initial one.

The original rescaled equation would agree with the unperturbed equation
modulo at least 3^9. Its failure already modulo 3^6 is therefore a rigorous
contradiction. Python and a standalone C++ verifier independently reproduce
the complete finite certificate, with matching counts and checksums at all
four levels. This is a computer-assisted impossibility proof for the stated
ansatz, not a finite search over rational heights.

## Formal and rational-function work

[FORMAL_BRANCH.md](FORMAL_BRANCH.md) gives the unique exact formal branch,
two verified jets through orders 72 and 30, and the exact first error of a
Padé reconstruction. A rational Z alone cannot evade the obstruction:
its monic equation forces its denominator to cancel.

[RATIONAL_CURVE_ATTEMPT.md](RATIONAL_CURVE_ATTEMPT.md) records every bounded
subfamily, its first obstruction, and the smallest genuine extension in
the fixed tangent direction. A single additional polynomial coefficient
or a linear denominator cannot produce a new reduced degree-seven curve.

## Inheritance and validation

Base main commit: `80a2d9f5bef3c7ee4b622a8ec96254ac05c9bd0f`. It was read
and matched to the remote. The first-strike reports, scripts, exact
certificates, logs, and manifest were read before this attack; all 33
manifest entries matched at the start. First-strike result files, attack
scripts, reports, and ledger remain unchanged.

All inherited gates passed:

```bash
make clean
make -j"$(nproc)"
make validate control equiv frontier-audit prototypes
make differential
make direct-checks
```

The inherited `frontier-audit` checks the retained candidate counts and
digests. It does not redo historical production decompositions. The
`direct-checks` target replays the existing rejection certificates; it
does not relaunch the targeted search. No new integer targets, rational
height boxes, literature searches, or frontier extensions were attempted.

New gate:

```bash
make boundary-checks
```

This checks exact symbolic identities, rational formal coefficients,
24 linear-lifting controls, the full Python obstruction, and the independent
C++ certificate. Evidence and hashes are in
`results/astra_second_2026_09_04/`.

## Stop condition and single third strike

Stop condition **D**: the proposed polynomial ansatz is rigorously impossible
and a precise minimal rational relaxation is obtained. Construction work
stopped; subsequent work verified and recorded that result.

Attack the degree-eight extension with a positive quadratic denominator:

```text
Q=1+q t^2, q in Q, q>0,
X=t P7/Q, Y=t R7/Q, Z=T8/Q, W=Z+(2/3)t^6,
P7(0)=R7(0)=T8(0)=1.
```

With `N=T8+t^6Q/3`, solve the exact identity

```text
(P7^6+R7^6)/2 = Q N^5+(10/27)t^12 Q^3 N^3+(1/81)t^24 Q^5 N.
```

The quadratic is forced to have no real zero by the contact divisor, and
centering it uses the available Möbius freedom. Keep P7,R7 independent.
An identity would immediately give positive rational specializations near
zero. This is a live next construction, not a claim that a curve exists.
