# ESOP6: executed full-fourfold search

2026-09-10. Base repository: `git-df-scott/ESOP6`, main commit
`36e03d4ccc3a2f1a2b7d8d77ea3652cfb758a3cd`.

**No counterexample found.** Two research agents and the coordinator executed
the work below. The latest authorized account budget was14%, replacing5%.
No live account meter was available, so actual percentage consumption cannot
be certified. All computational jobs and agent tasks in this pass have ended.
No generic integer-height sweep was extended, and no remote repository push
was performed. This package preserves the new work separately.

## Main outcome: ten exact modular branches, no rational reconstruction

The full-fourfold conic lane produced actual surviving coefficient data.
It does not require two repeated pairs.

Primitive quadratic coordinate forms reduce modulo7 to a single surviving
summand X_j and W=zeta*X_j. The proof is in the independent-geometry report.
This is a reduction statement, not a rational boundary attachment assumption.

In the first nontrivial valuation chart put the other four coordinates equal
to7*U_i and X_j=W-7^6*H, with U_i,W,H binary quadratics. The exact divided
original equation is

    sum U_i^6 -6HW^5 +15*7^6 H^2 W^4 -20*7^12 H^3 W^3
      +15*7^18 H^4 W^2 -6*7^24 H^5 W +7^30 H^6 = 0.

Its first reduction is sum U_i^6=6HW^5 modulo7. Exhaustive finite matching
over all58 projective quadratic classes including zero found228 proportional
and84 nonproportional solutions. The84 have irreducible W, four nonzero U_i,
and H proportional W modulo7. They are normalized modular seeds, not84
integer counterexample candidates.

| Modulus | Surviving selected branches |
|---|---:|
| 7 | 84 nonproportional seeds |
| 49 | 84 |
| 343 | 63 |
| 7^4 | 12 |
| 7^5 | 10 |
| Each of7^6 through7^12 | 10 |

After modulo49, one canonical correction branch per surviving seed was
followed, with free correction variables set to zero. Alternative corrections
were not exhausted. Failures therefore exclude selected continuations only.
The228 proportional solutions and other valuation charts remain outside this
continuation. All displayed lifts were checked by direct integer coefficient
expansion. The original higher terms were included when their precision
became relevant; this is not merely a lift of the tangent equation.

At modulus7^12=13,841,287,201, coefficientwise rational reconstruction used
denominator bound10,000 and numerator bound floor((M-1)/20,000), a uniqueness
bound. None of the ten branches yielded a complete reconstructed tuple.
Consequently no rational conic reached final exact specialization or integer
verification. These are not proved Q7 curves and not evidence of Q-points.

The next targeted operation would vary the seven free lifting directions
systematically rather than treat canonical-digit branches as representative
of the rational locus. It needs a declared parameter/height budget; this
package does not leave that computation running.

## Numerical construction attempts actually run

The repeated-pair degree-eight pilot completed24 starts, plus8 conditioning
replays in rescaled coordinates. No residual reached its rational
reconstruction threshold. This does not exclude that family.

The full-fourfold conic pilot completed24 starts:12 near an explicit algebraic
real conic and12 independent Gaussian starts. It used W=1+q*t^2 for q=1,2,3
and fixed a5=1/2,b5=0. These are stated restrictions, not an exhaustive conic
normalization. Minimum floating coefficient residual was approximately
3.55e-15. All60 coefficientwise rational reconstructions at denominator bounds
10,100,1000,10000 failed exact polynomial verification. This reconstruction
procedure is not exhaustive over rational coefficient boxes.

All56 optimization starts have retained coefficients, residuals and evaluation
counts. No complex random-start or genuine quartic search is claimed. Cubic
reduced real maps were skipped for the odd-degree obstruction: W would have
a real projective zero, forcing all coordinates to vanish there.

See NUMERICAL_RUNS_2026_09_10.md and the degree8/full_conics subdirectories.

## Exact constraints and a realness shortcut

For the four-independent-coordinate boundary direction (1,2,4,3), c=815,
the agent established that deg M must equal6, some A_i must have degree5,
and at least one M coefficient has an even denominator. If b is the ordinary
2-adic coefficient Gauss valuation of M, b=0,-2,-3,-4 are excluded; b=-1
remains open. The proof allows arbitrary rational A_i coefficients and does
not extend its Gauss-valuation lemma to ramified coefficient fields.

A separate proof excludes every identity sum_{i<=5} L_i(u,v)^6=q(u,v)^3
with rational linear L_i and nondegenerate binary quadratic q, even over Q7.
This closes centrally symmetric conic constructions, not translated conics.

Real conics themselves are easy to exhibit. Let U=s^2-t^2, V=2st,
W=s^2+t^2, and k=(63/80)^(1/6)>0. The six coordinates

    kU, kV, k(U+V)/sqrt(2), k(U-V)/sqrt(2), W/2, W

satisfy the sextic identity exactly and are all nonzero at generic real
parameters. However x3+x4=sqrt(2)*x1 and x3-x4=sqrt(2)*x2 prove this conic
has no rational projective points. Thus discovering real conics numerically
does not resolve the arithmetic obstacle. The identity and obstruction are
retained and exactly checked.

## Local densities and the limits of a first-height estimate

The exact primitive local calculation found

    delta2=5/8, delta3=20/81, delta7=180/16807.

The product of primitive factors over primes through199 is approximately
0.006348102641304657. Bad-prime stabilization was proved, and44 good primes
were counted exactly. This is a finite product, not the entire singular
series and not a discovery probability.

For reference, the formal positive-real radial integral is J*log N, with

    J=Gamma(1/6)^5 / (6^5 Gamma(5/6))
     =Gamma(7/6)^5 / Gamma(5/6)
     ≈0.6087949289329304.

This follows by substituting y_i=x_i^6 in the five-coordinate delta-function
integral at f=1, and then integrating df/f. Combining J with the truncated
primitive product and dividing by5! for generically distinct unordered
summands gives the purely formal coefficient≈3.220577246976681e-5 of log N.
No first-height estimate is inferred. There is no justified Poisson model,
no supplied global counting theorem, and rational curves can create
accumulating families that invalidate the generic logarithmic picture.

Unrestricted local factors would be wrong here: their finite densities obey
alpha(p,k)=delta(p,k)+alpha(p,k-6), so their limits diverge. Primitive
normalization is essential. Full exact counts and the recurrence proof are
in local_density_2026_09_10/.

## Interpretation of the user's proposed pivot

Releasing repeated pairs is justified. Their surface exclusions do not
establish corresponding fourfold exclusions. Nevertheless an N^-2 box/dyadic
heuristic for the restricted surface is not a decreasing cumulative count,
and does not supply a calibrated probability that there are no rational
points. Expected curve dimensions likewise are virtual dimensions, not
existence or nonexistence theorems for this special diagonal variety.

Signs can be removed after clearing denominators because sixth powers are
even. Rationality and nonzero coordinates remain indispensable. A rational
curve construction must be nonconstant, and its parameter curve must supply
rational points; a complex curve or an unparametrized Q-conic without a
Q-point is not enough. Complex sign-twist equivalence does not automatically
descend a curve to Q on the desired twist.

The strongest retained object from this execution is therefore the ten
full-equation modular conic branches, with every limitation above preserved.
No claim is made that they are close to a rational counterexample.

## Package and reproduction

All paths below are relative to research/ in the original workspace or to
the research/ directory in this package. Python3 scripts use standard
library arithmetic and, where stated, NumPy/SciPy/SymPy. These dependencies
must be available when rerunning.

- independent_geometry_2026_09_10/REPORT.md: complete conic proofs, sieve,
  original-equation lift provenance and real-conic verification.
- independent_geometry_2026_09_10/original_conic_lifts.json: ten final branches.
- four_coordinate_2026_09_10/REPORT.md: denominator obstruction and proof.
- local_density_2026_09_10/REPORT.md: exact primitive densities.
- NUMERICAL_RUNS_2026_09_10.md: precise numerical domains and replay commands.

Proof regressions and saved modular lifts were verified during execution;
the included reports distinguish those checks from counterexample searches.
No sixteen-million or million-scale lemma/sieve count is presented as that
many tested positive-integer ESOP6 candidates.
