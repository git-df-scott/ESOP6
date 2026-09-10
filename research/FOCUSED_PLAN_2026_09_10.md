# ESOP6: focused counterexample plan

2026-09-10. Repository inspected at `36e03d4ccc3a2f1a2b7d8d77ea3652cfb758a3cd` on `main`.

**No counterexample found in this pass.** This is a bounded planning and exact-algebra pass, not an executed large search. The user authorized up to 5% usage now. No live account meter was available, so actual percentage consumption cannot be certified. No long-running jobs were launched. The larger search stages below are proposed work, not completed results or background tasks.

## 1. Target and decision

Find positive integers, with repetitions allowed,

\[
a^6+b^6+c^6+d^6+e^6=f^6.
\]

Equivalently, find a rational point with all coordinates nonzero on
\(\mathcal X:x_1^6+x_2^6+x_3^6+x_4^6+z^6=w^6\).
Clear denominators, take absolute values, and divide the common gcd. Negative rational coordinates are harmless because the exponent is even; zero coordinates fail this project's five-positive-summand target.

**Recommended order:** one bounded structural attempt at the remaining degree-eight repeated-pair family; then prioritize a degree-six boundary construction on the full fourfold with four independent small coordinates. Retain a separate bounded arithmetic lane for single repeated pairs and spread residue roles. Do not resume the generic height sweep.

The rationale is specific: the old two-variable 3-adic anisotropy argument does not extend to four independent sixth powers. This removes a known obstruction, but does not predict a solution.

## 2. What the repository actually establishes

| Result | Scope and planning consequence |
|---|---|
| Class-1 frontier through 4,300,000 | Historical zero-result; 55,684 candidate pairs independently regenerated. Lost decomposition logs prevent calling this an independent full replay. It is not an all-class bound. |
| Four selected quotient conics, no lift | Finite parameter boxes only. The generic square-lift curve has genus nine, so a rational quotient point is far from an ESOP6 point. |
| 10,382 retained sparse surface targets, no solution | Selected targets, not an interval exclusion to their largest height. |
| Equal-slope degree-six repeated-pair boundary family | Closed over Q, including rational denominators, by the recorded reduction and modulo-729 certificate. Do not rerun it as a discovery lane. |
| Degree-eight quadratic-factor extension | Proposed in main; still a construction problem, not an identity. The reduction below narrows it. |
| Four independent small coordinates | Outside the repeated-pair proof. The specific construction below has not been solved in this pass. |

Recovered conversation context also reports a later 82,915-target sample and a square-coefficient restriction for degree eight. Those artifacts are not in this inspected main snapshot. Do not add the target counts together or claim this pass independently reproduced that later campaign. The degree-eight restriction is rederived below; it is not claimed as a first discovery.

Sources: [STATUS.md](https://github.com/git-df-scott/ESOP6/blob/36e03d4ccc3a2f1a2b7d8d77ea3652cfb758a3cd/STATUS.md), [CANONICAL_FRONTIER.md](https://github.com/git-df-scott/ESOP6/blob/36e03d4ccc3a2f1a2b7d8d77ea3652cfb758a3cd/CANONICAL_FRONTIER.md), [BOUNDARY_CONTACT_6.md](https://github.com/git-df-scott/ESOP6/blob/36e03d4ccc3a2f1a2b7d8d77ea3652cfb758a3cd/BOUNDARY_CONTACT_6.md).

## 3. Corrections that determine the search

1. A rational point need not be known to lie on a rational curve. Finding a curve is a sufficient construction strategy, not a necessary condition for a counterexample.
2. The repeated-pair surface \(2X^6+2Y^6+Z^6=W^6\) is a restrictive slice. Its curve exclusions do not automatically apply to the full fourfold.
3. General-type geometry does not exclude isolated rational or elliptic curves. Conversely, a cubic quotient or formal power series does not provide an exact rational lift.
4. No numerical residual, local p-adic solution, or finite empty coefficient box settles the rational-point problem.
5. A proposed parameter normalization must preserve the family. In particular, fixing the two initial slopes equal is a restriction, not a universal symmetry.

## 4. First bounded attempt: reduce degree eight before solving

Use the repository's exact family, with independent polynomials P,R of degree at most seven and N of degree at most eight:

\[
\frac{P^6+R^6}{2}=QN^5+\frac{10}{27}t^{12}Q^3N^3+\frac1{81}t^{24}Q^5N,
\quad Q=1+qt^2,\quad q>0,
\]

with P(0)=R(0)=N(0)=1. The cleared projective coordinates are
\([tP:tR:N-t^6Q/3:N+t^6Q/3]\).

**Necessary condition for a reduced degree-eight map: q is a rational square.**

Proof: Q is irreducible over Q. In the quadratic field K=Q[t]/(Q), the identity gives P^6+R^6=0. If one is zero, both are zero. Then Q divides P and R, and comparison of Q-adic orders in the identity forces Q to divide N. All four projective coordinates have the common factor Q. Cancelling it returns the excluded equal-slope degree-six family.

Otherwise P/R in K has sixth power -1. Such a root has order four or twelve. An order-twelve root has minimal polynomial \(u^4-u^2+1\), of degree four, so cannot lie in a quadratic field. Thus K contains i and equals Q(i). Since K=Q(sqrt(-q)), q must be a rational square.

Write q=s^2 with s>0 rational. For remainders
P mod Q=p0+p1*t and R mod Q=r0+r1*t, choose a sign epsilon in {+1,-1}. Evaluation at t=i/s gives

\[
p_0=-\epsilon r_1/s,\qquad p_1=\epsilon s r_0.
\]

These are two linear relations. Keep both signs unless their equivalence under a permitted symmetry is proved. For a reduced solution, N is nonzero modulo Q and the order of P^6+R^6 along Q is exactly one. Enforce these conditions to discard cancelled maps.

**Implementation sequence:**

1. Parameterize the two remainder relations exactly and retain rational s as a variable.
2. Eliminate the eight nonconstant N coefficients recursively: the derivative of Q*N^5 at t=0 is five. Retain all remaining coefficient equations.
3. Saturate away common projective factors and any denominators introduced by elimination. Split determinant-zero cases explicitly.
4. For bounded rational s, first test coefficient systems over small good primes (initially 11, 13, 19), then lift surviving branches. Exclusion of all Q-points requires denominator/valuation control; an integral residue search alone does not supply it.
5. Try s=u/v, gcd(u,v)=1, 1<=u,v<=8 as an initial finite diagnostic. Cap a difficult specialization at 30 seconds and the entire pilot at ten CPU-minutes. These are future computational caps, not a conversion of account usage.

**Stop:** an exact identity; a proof excluding the full specified family; or expiration of the pilot with a retained list of survivors, excluded specializations, and timeouts. A timeout is unresolved. No automatic coefficient-height or degree escalation follows a null pilot.

Source family: [RATIONAL_CURVE_ATTEMPT.md](https://github.com/git-df-scott/ESOP6/blob/36e03d4ccc3a2f1a2b7d8d77ea3652cfb758a3cd/RATIONAL_CURVE_ATTEMPT.md).

## 5. Main new direction: four independent small coordinates

At the boundary point [0:0:0:0:1:1], set

\[
x_i=tA_i(t)\ (i=1,2,3,4),\quad
w=z+c t^6,\quad M=z+\frac c2t^6,
\]

where deg A_i<=5, deg M<=6, M(0)=1, and A_i(0)=alpha_i. Choose positive rational alpha_i, normalize alpha_1=1, and put

\[
c=\frac{\alpha_1^6+\alpha_2^6+\alpha_3^6+\alpha_4^6}{6}.
\]

The exact divided identity is

\[
\boxed{\sum_{i=1}^4 A_i^6
=6cM^5+5c^3t^{12}M^3+\frac38c^5t^{24}M.}
\]

This equation is sufficient: any exact solution has all six coordinates positive for sufficiently small positive rational t. A formal solution or truncation is insufficient.

**Why this differs from the closed family.** A sum of two rational sixth powers has 3-adic valuation divisible by six. Four terms can violate that pattern. For example, the slope vector (1,2,4,3) has sixth-power sum 4890, with v3=1, and hence c=815. This is a tangent-direction example, not a point on the sextic. The four A_i remain independent throughout; repeated constant slopes would also not require repeated coordinate polynomials.

**Implementation sequence:**

1. Derive the coefficient recurrence with exact rationals. Its linear coefficient for m_n is 30c, nonzero. Eliminate m_1 through m_6 and retain residual coefficients 7 through 30.
2. First run the explicitly displayed slope vector. Next use normalized triples alpha_2,alpha_3,alpha_4=u/v with coprime 1<=u,v<=4, quotienting by coordinate permutations only when justified. Keep a manifest of the finite domain.
3. Prioritize slopes with three 3-adic unit entries and one divisible by three; compare with four-unit and other valuation patterns. These are exploratory strata, not an exhaustive rational slope classification.
4. At bad primes 2,3,5,7, use the undivided integral equations and explicit valuation charts. Do not invert 30c modulo a prime dividing it. Use good primes for cheap coefficient screening.
5. Compute the Jacobian rank and kernel on surviving modular branches. Treat rank-deficient branches separately. Lift with exact arithmetic, attempt rational reconstruction only with a stated height bound, then verify the entire polynomial identity.
6. If symbolic growth is prohibitive, introduce a declared sparse coefficient support and record it as a restricted family. Do not impose reflection or repeated-pair symmetry merely to make the solver finish.

**Pilot cap:** ten slope strata or ten CPU-minutes, whichever comes first, with a 30-second per-hard-solve timeout. Continue only if an exact branch, a useful denominator theorem, or demonstrable reduction in the remaining system justifies it. Numerical small residuals alone do not justify further spending.

**Failure meaning:** this could close or fail to resolve selected degree-six boundary families. Curves avoiding the boundary, higher degrees, and isolated rational points remain outside it.

## 6. Independent fallback: one repeated pair, with spread roles

In any primitive ESOP6 tuple, exactly one summand is a unit at each of 2,3,7. Consequently every value appearing at least twice is divisible by 42. Two repeated pairs leave only one unpaired summand to carry all three unit roles, forcing class 1. A single repeated pair leaves three unpaired terms, allowing the spread cases the main historical search did not cover deeply.

Use
\[
2(42r)^6+b^6+c^6+d^6=f^6,
\]
with the three unit roles explicitly allocated among b,c,d. Also retain multiplicity patterns 3+1+1 and 4+1 in the domain specification; do not confuse a distinct-only engine with coverage of repeated summands.

Before production, benchmark an exact hash or sorted join on a declared sparse domain of (f,r,b), finishing with c^6+d^6=f^6-2(42r)^6-b^6. Pair tables must allow c=d. Recover planted two-, three-, and four-term repeated-power targets and compare all answers on a small domain with independent exhaustive enumeration.

There is no demonstrated exponent improvement here. Advance this lane only if the bounded pilot identifies genuinely new coverage and an affordable cost. Do not repeat the old two-pair target sample or extend f merely because memory is available. A missing repeated-summand guarantee in an external solver is a reason to audit its scope, not evidence that a solution is nearby.

## 7. Acceptance, evidence, and spending discipline

For an identity: verify every coefficient with exact arithmetic; remove the homogeneous gcd; prove the map is nonconstant and each desired coordinate is nonzero; choose an explicit rational specialization avoiding all zeros and poles.

For a sextuple: stop searching and run both existing verifiers:

```sh
python3 tools/verify_esop6.py a b c d e f
node tools/verify_esop6.mjs a b c d e f
```

Retain the full integers, both sums, gcd, and both outputs. The normalized exact difference must be zero and all six integers positive. No confidence score replaces this certificate.

Each future run needs source revision, exact parameters, seed, arithmetic type, bounds, elapsed time, output, and a disposition of FOUND / EXCLUDED IN STATED DOMAIN / UNRESOLVED. Record timeouts separately. An exclusion of a restricted family is not a near-counterexample.

If account metering becomes available, record the baseline and stop before a five-percentage-point increase, allowing room for verification and reporting. CPU time and token counts do not reliably measure that percentage. The pilot caps above do not override the user's account budget.

## 8. Work completed in this pass

- Cloned and read the current main snapshot, including status, frontier, handoff, geometry, and the degree-six proof.
- Rederived the degree-eight square-q restriction and its linear remainder relations.
- Derived the full-four-coordinate midpoint equation. Standard-library Fraction/binomial checks confirmed coefficients 6, 5, 3/8.
- Checked the sixth-sum factorization and 100 exact rational controls of the quadratic remainder relation.
- Verified the slope sum 4890 and v3=1.
- Ran the repository's independent verifier controls: all 18 passed using Python integers and JavaScript BigInt. There is no known positive ESOP6 fixture in those controls.
- Did not rerun the modulo-729 exhaustive certificate, extend a height frontier, solve the new coefficient systems, or claim an exact curve.

The first symbolic-check attempt encountered an unavailable SymPy dependency. The completed checks used Python's standard library instead; no dependency installation was needed.

For background, [Bremner–Choudhry–Ulas](https://arxiv.org/abs/1402.4583) construct families of diagonal quartic and sextic surfaces with rational points. Their abstract does not supply the fixed ESOP6 coefficients or an automatic lifting theorem. This plan derives its actual equations directly rather than treating that literature as a solution.

**Research judgment:** the next useful result is an exact rational branch in the four-coordinate system, or a precise obstruction showing why it fails. The degree-eight repeated-pair family deserves one reduced, capped attempt. Neither route is guaranteed to yield a counterexample.
