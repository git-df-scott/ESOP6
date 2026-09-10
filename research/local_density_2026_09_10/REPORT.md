# ESOP6 primitive local-density audit

2026-09-10. Exact finite congruence counts. No counterexample found or claimed.

For F=x1^6+...+x5^6-x6^6 and q=p^k, let N(p,k) count all solutions modulo q
and P(p,k) count only solutions with at least one coordinate a p-adic unit.
The finite primitive density is delta(p,k)=P(p,k)/p^(5k).

## Actual computation

`local_density.py` computes sixth-power residue multiplicities, convolves
five copies, and dots the result with the sixth-power histogram. It separately
repeats the calculation with every coordinate divisible by p and subtracts
this count; singular all-divisible tuples are not silently retained.
All counts use integer arithmetic; fractions are exact until display.

Computed:

* Every good prime 5<=p<=199 (44 primes).
* p=2, exponents k=1,...,10.
* p=3, exponents k=1,...,6.
* Independent direct six-variable brute-force checks at q=2,3,4.
* Direct computations modulo p^2 for p=5,7,11,13, verifying smooth lifting.
* Closed quadratic-form point count at every tested p=5 modulo 6:
  N(p,1)=p^5+(p-1)p^2.

All controls passed. Complete integers and fractions are in results.json.

| Prime | Primitive density | Status |
|---|---:|---|
| 2 | 5/8 | Exact for every k>=3 |
| 3 | 20/81 | Exact for every k>=2 |
| 5 | 3224/3125 | Exact stable good-prime density |
| 7 | 180/16807 | Exact stable good-prime density |
| 11 | 162260/161051 | Exact stable good-prime density |
| 13 | 1050840/371293 | Exact stable good-prime density |
| 19 | 2113020/2476099 | Exact stable good-prime density |
| 31 | 21162420/28629151 | Exact stable good-prime density |

The finite product, including the stabilized factors at 2 and 3 and all
good primes through 199, is approximately **0.006348102641304657**.
This is a finite local product, not an inferred discovery probability.

## Why the bad-prime densities stabilize

Modulo 8, odd sixth powers are 1 and even sixth powers are 0. A primitive
solution therefore has exactly one odd left coordinate and an odd right
coordinate. For q=2^k, k>=3, sixth powers of units are precisely the residue
classes 1 modulo 8, with exactly four unit preimages per such class.
There are 5*(q/2)^5 admissible left tuples and four right preimages each.
Thus P(2,k)=20*(q/2)^5 and delta(2,k)=5/8 exactly.

Modulo 9, unit sixth powers are 1 and nonunit sixth powers are 0. A primitive
solution again has exactly one unit left coordinate and a unit right
coordinate. For q=3^k, k>=2, unit sixth powers are precisely 1 modulo 9,
with six preimages each. Thus

    P(3,k)=5*(2q/3)*(q/3)^4*6=(20/81)*q^5.

For p not dividing 6, every primitive solution modulo p is smooth, since
some derivative +/-6*x_i^5 is a unit. Each lifts with a factor p^5 per
precision level, proving delta(p,k)=(N(p,1)-1)/p^5 for every k>=1.

At p=7, nonzero sixth powers equal 1. There must be exactly one nonzero left
coordinate and a nonzero right coordinate, so N(7,1)-1=5*6^2=180.

## The critical-dimension trap

Define the unrestricted density alpha(p,k)=N(p,k)/p^(5k).
For k>=7, writing every coordinate as p*y gives exactly

    N_nonprimitive(p,k)=p^30*N(p,k-6).

The extra p^30 arises from the five unused p-adic digits of each of six y
coordinates. Consequently

    alpha(p,k)=delta(p,k)+alpha(p,k-6).

Since the primitive densities are eventually positive and constant here,
the unrestricted density diverges linearly along each k modulo 6.
This is the degree-equals-number-of-variables issue in concrete form.
Multiplying unrestricted limiting densities as if they were finite local
factors would be an error. The report uses primitive factors throughout.

## What this does and does not establish

These local tests show no local obstruction at the tested primes. The
trivial primitive point (1,0,0,0,0,1) already gives a reason to expect local
solubility. Local solubility does not establish a rational point with all
five left coordinates nonzero, much less a positive integer counterexample.

The strong loss at 7 is partly offset by enrichment at primes such as 13.
Neither a small nor a positive finite product measures distance to a
counterexample. A rational-point count or first-height forecast would
additionally require archimedean normalization, treatment of projective
scaling and positivity, and a justified global heuristic/asymptotic.
No such forecast is made. No convergence claim for the infinite product
is inferred from truncation at 199.
