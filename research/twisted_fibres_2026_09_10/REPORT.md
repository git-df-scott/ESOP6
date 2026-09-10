# ESOP6: nine elliptic gates and a rigorous quadratic-Chabauty shortlist

2026-09-10. Branch `claude/inspiring-pasteur-lz5nh3`, PR #2.
Starting commit `41b0583ae18d6d27a7dd4b8beca7ff4d34a11e17`.

Before final publication, concurrent branch commit `3090aa86b9005d51b15651088e8deba162f5621c`
was integrated. Its completed N300 log snapshot and new moduli/S3 scripts
were preserved. The inherited inventory records the earlier starting snapshot;
the later logs do not upgrade numerical coverage to certified curve elimination.

**No positive-integer counterexample was found.** The main outcome is
**573 entire rational fibres excluded, out of 627**, leaving **54 unresolved**.
These exclusions have no height bound on the remaining fibre coordinates.
The new structural calculation supplies **38 surviving fibres with a
certified genus-two rank-(1,1) route**, and 35 with full genus-ten rank
bounds satisfying a quadratic-Chabauty finiteness criterion.

The arithmetic proof and implementation checkpoint was published as
`888434d801bdd2971c8e510c0d4367a94eebe3de`. This report records the completed
execution that followed. No historical integer sweep or numerical curve
lottery was restarted. No research process from this execution is left running.

## What was derived before the fibre search

For `C_c:A^6+cS^6=B^6`, the square map to `X^3+cY^3=Z^3` has degree four.
On `Y!=0`, `u=Z/Y`, `v=X/Y` give `c=u^3-v^3`; a positive lift requires
both u and v to be nonzero rational squares. The exact Mordell model and
inverse are

```
eta^2=xi^3-432c^2,
u=(eta+36c)/(6xi), v=(eta-36c)/(6xi).
```

Each individual square condition gives genus four; their composite is
genus ten. The coordinate divisors give nine branch points with order-two
inertia, including diagonal inertia at Z=0. These ramified covers are not
ordinary unramified 2-Selmer covers.

**Every c-fibre is locally soluble with nonzero coordinates at every place.**
The smooth rational boundary `[1:0:1]` has nearby nonzero local points.
Consequently an unrestricted local or cube-difference congruence filter
cannot reject c. The implementation retains the correct 2,3,7 valuation
charts and uses them only as chart restrictions.

Five additional degree-six elliptic maps arise from genus-two quotients;
three further maps from product quotients have degree twelve. Together
with the degree-four cubic quotient they give nine Weierstrass coefficients:

```
-c^3, c^3, c, -c, c^2, -432c^2, 4c, -4c, 4c^4.
```

Every map has an explicit exact inverse lift test. A certified rank upper
bound zero plus failure of **every** rational torsion point to lift excludes
the whole fibre. Exact formulas and proofs are in
[ARITHMETIC.md](ARITHMETIC.md) and
[INDEPENDENT_AUDIT.md](INDEPENDENT_AUDIT.md).

## The exhausted coefficient domain

Enumerated all 715 sorted fourtuples with entries in `1..10`; removing the
89 nonprimitive tuples leaves 626 distinct primitive rays. Added the two
retained exact Y seeds, without searching for them again. There are 628
representations but 627 normalized sixth-power classes.

The one duplicate coefficient is

```
2^6+2^6+9^6+9^6 = 3^6+5^6+6^6+10^6 = 1063010.
```

Representations and scaling factors are retained. Elliptic arithmetic is
cached after independent sixth-power normalization of its Weierstrass
coefficient, preserving all square/cube lift classes. There are 5,624
possible distinct normalized elliptic models for the full nine-test list;
early rejection required only 2,079 model tasks.

| Gate coefficient | Additional whole fibres excluded | Remaining |
|---|---:|---:|
| `-c^3` | 267 | 360 |
| `c^3` | 58 | 302 |
| `c` | 62 | 240 |
| `-c` | 35 | 205 |
| `c^2` | 62 | 143 |
| `-432c^2` | 61 | 82 |
| `4c` | 17 | 65 |
| `-4c` | 5 | 60 |
| `4c^4` | 6 | 54 |

Thus **91.4% of this coefficient pilot is globally excluded**. This is
not an unrestricted six-variable height frontier or an estimate of the
fraction of all possible ESOP6 solutions. The 54 remaining fibres are not
counterexample candidates with known positive points.

The Y seed with tail `(48,175,228,350)`, c=`2007479003912993`, is excluded
by the `-c^3` rank-zero gate. The other retained Y seed, with tail
`(75,119,273,279)`, c=`888650267843116`, remains unresolved.

## Certification and bounded point search

PARI/GP 2.17.4 completed 1,984 model tasks with certified rank intervals;
95 tasks timed out and remain explicitly inconclusive. The first six gates
used 8 seconds per model, the final three 12 seconds, with three independent
workers. Their stage wall times were approximately 668 and 338 seconds.
No failed model was treated as rank zero. The version-pinned cubic class
groups were explicitly passed through `bnfcertify` before their rank
bounds were accepted. Raw transcripts and program hashes are retained.

An independent standard-library replay checked all 627 fibre rows, 628
representations, 4,398 exact model points and 5,844 exact lift attempts.
It independently proved full torsion for all 1,984 certified models:
1,864 by finite-field order bounds, and the other 120 by those bounds plus
the exact 3-division polynomial. All 573 exclusions passed. The rank
computations themselves remain PARI descent results; the replay is not
a second implementation of elliptic descent.

The final rational-point search exhausted 2,850 coefficient–torsion tuples
over 48 tasks: 46 selected survivor-model subgroups and both retained Y
seed orbits. Rank-one coefficients lay in `[-12,12]`; two or three selected
generators used `[-2,2]` in each coordinate, with full torsion translates.
The selected points were not claimed to form full Mordell–Weil bases.

Necessary finite-field power-class tests rejected 2,796 tuples. All 54
materialized exact elliptic points were boundary: 48 identities and six
zero-coordinate lifts. **No nonboundary sextic point reached the integer
certificate gate.** In total 47 of the 54 unresolved fibres received a
selected subgroup search; seven lacked an available subgroup. These finite
misses did not change the fibre exclusion counts.

The point search was validated against independent PARI group arithmetic,
exact torsion controls and 81 constructed generic positive fibre points
across all nine models. Those generic controls are not ESOP6 solutions.
The existing Python/JavaScript integer checkers agreed on all 18 controls.

## Full Jacobian decomposition and the next serious attack

Ten explicit elliptic maps pull back to the ten independent regular
differentials of the plane sextic. This proves over Q that

```
Jac(C_c) ~ E_c * E_(-c) * E_(c^2)^2 * E_(-c^3) * E_(c^3)
             * E_(-432c^2) * E_(4c) * E_(-4c) * E_(4c^4),
```

where `E_b` means `y^2=x^3+b`. Nine coefficient formulas account for ten
elliptic factors because c^2 occurs twice. The independent verifier checks
all ten differential identities, rather than inferring the decomposition
from a dimension heuristic.

The rational Neron–Severi rank is at least eleven. Therefore a certified
Jacobian rank upper bound at most nineteen satisfies the
Balakrishnan–Dogra quadratic-Chabauty finiteness inequality. All 35 remaining
fibres with complete nine-model rank intervals meet this bound. This does
not mean their rational points have been computed.

More immediately, **38 remaining fibres have at least one genus-two quotient
with both elliptic ranks exactly one**. These fit the more explicit
bielliptic framework. Start with

```
c=67=1^6+1^6+1^6+2^6,
H_67: w^2=t^6+67.
```

Its two elliptic ranks are one. The exact control `(t,w)=(7/6,1801/216)`
fails the necessary cube test on w. A complete determination of `H_67(Q)`,
followed by exact cube testing and reconstruction, would settle that whole
fibre. Its genus-ten Jacobian has rank exactly twelve, so ordinary
rank-less-than-genus Chabauty is not the right tool.

See [QC_HANDOFF.md](QC_HANDOFF.md) and `survivors.json`. No p-adic height
calculation, Coleman integration, quadratic-Chabauty point classification,
or proof of completeness of a genus-two point list was executed here.

## Secondary quotient results and inherited-data repairs

The [quotient report](../quotient_descent_2026_09_10/REPORT.md) proves that
all six coordinate slices obtained by fixing two normalized tail coordinates
at either algebraic seed are locally impossible to lift. A viable family
must vary at least two tail ratios. The Gaussian seed also lies on the
extra quadratic locus `sum(x_i^2)=x6^2`; preserving that locus permits only
trivial real points. No rational quotient family was found.

An exact complete local test for the Pythagorean subfamily reduced **202,861**
primitive parameter pairs `1<=n<m<=1000` to **four everywhere locally soluble
slices**. Their parameters are `(428,85)`, `(883,540)`, `(964,343)`, and
`(995,652)`. This is an all-place local statement on those slices, not a
global rational-point result.

The inherited N500 Y2 list contained **31 invalid entries** among 2,924:
the fast generator compared to a rounded binary64 sixth power. For example
`[93,191,273,425], f=457, S=52473038` has exact residual `-1`.
The verified replacement contains 2,893 distinct valid seeds; corrected
>300 chunks are supplied. Original lists remain unchanged as evidence.
The generator's final comparison and the consumer's initial guard now use
the actual integer equation. The floating proposal stage still has no
completeness proof. The original N40/N70 counts also included duplicates;
their unique sorted-coordinate keys number 19/68, respectively.

All 1,239 retained N300 seeds passed exact equation checking. Historical
lottery logs were inventoried and not rerun. No corresponding process was
visible in this workspace; the other session's remote process state could
not be inspected or controlled, so no claim is made that its jobs stopped.

The night report's LPS-based boundary exclusion is conjectural; reducibility
does not imply a rational root; random numerical starts do not certify
complete curve enumeration. [CORRECTIONS.md](CORRECTIONS.md) explicitly
supersedes these broad claims while preserving valid restricted proofs.

## Literature access and evidence

Bremner 1981 and Kuwata 2007 were identified precisely, but their publisher
full texts could not be obtained. **Neither original was read.** No
Neron–Severi or involution result was attributed to them without access.
BCU 2014 sections 5–6 and the relevant Balakrishnan–Dogra primary results
were obtained and read. [LITERATURE.md](LITERATURE.md) records the gap and
the exact transferable statements; the Jacobian decomposition here is
proved directly by explicit maps.

Use [FIBRE_LEDGER.md](FIBRE_LEDGER.md) for all coefficients and dispositions,
`survivors.json` for the next targets, and `search_summary.json` for counts.
`arithmetic_evidence.zip` preserves the full ledger, every model result,
the six-gate checkpoint and logs with SHA256 manifest. The scripts and
[REPRODUCE.md](REPRODUCE.md) restore and independently audit the evidence.
