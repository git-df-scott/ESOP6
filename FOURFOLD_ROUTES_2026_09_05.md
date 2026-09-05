# Full-fourfold routes, 2026-09-05

**ESOP6 counterexample: NONE. Exact rational curve: NONE.**

Base: `36e03d4ccc3a2f1a2b7d8d77ea3652cfb758a3cd`.
This work moves the main construction to

\[
\mathcal X:\quad x_1^6+x_2^6+x_3^6+x_4^6+x_5^6=x_6^6.
\]

The target is a rational point with all six coordinates nonzero. Taking
absolute values and clearing denominators then gives positive integers.
There is no separate sign-chamber problem on this sextic. On its cubic
quotient, however, simultaneous square lifting requires a projective point
with a representative whose six coordinates are positive.

The results below are campaign derivations and reproducible computations,
not externally reviewed theorems. Finite congruences are explicitly
distinguished from characteristic-zero identities.

## What ran

| Route | Execution | Outcome and scope |
|---|---|---|
| Full fourfold, degrees 1 and 3 | Real-place argument | Every reduced odd-degree parametrization is excluded, in all degrees. |
| Full fourfold, conics | Complete F7 support argument; transverse residue census; Hensel lifting; rational reconstruction | Four nontrivial transverse patterns; each chosen full-conic path survives modulo 7^41 in the divided equation. No exact rational reconstruction. Conics over Q remain open. |
| Cubic quotient | Exact tangent construction; 84,612,096 normalized bounded directions | 1,926,206 positive cubic outputs, no square lifts. This bounded run is a control below the inherited height frontier. |
| Larger cubic points | 1,000,000 reproducibly sampled directions | 21,159 positive outputs, 21,128 with sixth-power height above 4.3M if square-liftable; no square lifts. No interval was exhausted. |
| Divisibility-based cubic seeds | A separate 1,000,000 sampled directions | 20,687 positive outputs, all above the 100M height threshold if square-liftable; no square lifts. No height-completeness claim. |
| Two square ratios by construction | Two explicit conic families; 587,692 rotations in total | 317,378 positive outputs, no complete lifts. The common equal-pair-sum restriction is now excluded exactly by the 3-adic argument below. |
| Degree-eight repeated surface | Quadratic-field restriction; sparse exact elimination; Groebner bases over F7, F11, F13 and Q | q must be a rational square for a reduced map; one explicitly stated sparse family has unit ideal over Q. The general degree-eight family remains open. |
| Full fourfold, boundary contact | Derivation with four independent transverse polynomials | An exact residual system, not solved. |
| Integer sweep 4.3M–4.4M | Existing caseA3, two threads, 90-second cap | Built the first of eight filters; timed out without a completed search result. No frontier advance, no claim that the interval is empty. |

Counts are outputs or directions, not distinct rational points: duplication
is possible across seeds and parametrizations. The H=10 control lies inside
the H=24 run and must not be added to it.

## 1. Degree constraints on the full fourfold

Represent a morphism P1_R -> X by homogeneous binary forms P1,...,P6 of
the same degree d with no common zero. At a real projective zero of P6,
the sextic identity forces all five other forms to vanish. Therefore P6
has no real projective zero. A real binary form of odd degree has such a
zero; consequently **d is even**.

This excludes lines and degree-three parametrizations over Q, including
odd-degree covers. It does not exclude conics or even-degree curves. A
genus-zero curve defined over Q but without a Q-point is not a
parametrization by P1_Q and would not itself supply a counterexample.

At an axis point, for example [0:0:0:0:1:1], a real nonconstant branch
with m the least positive transverse vanishing order satisfies

\[
\operatorname{ord}(P_6-P_5)=6m.
\]

This follows by factoring P6^6-P5^6; the other factor is nonzero there,
and real leading sixth powers cannot cancel. Hence any such curve needs
d >= 6. Normalizing a proposed conic's constant point to a known axis
point would therefore remove every possible nonconstant real conic.

The virtual dimension 6(d+1)-(6d+1)-4=1 is a deformation count, not an
existence theorem over Q. It does not imply that a finite-field component
has a rational characteristic-zero point. Negative expected dimension on
the repeated-coordinate surface likewise does not prove emptiness.

### Complete F7 coefficient constraint

For d <= 3, every nonzero homogeneous binary form has at least 8-d
nonzero values on P1(F7). At each of its eight points, a sixth power is
0 or 1. Five such terms can equal the sixth power on the right only when
zero or one of the left terms is nonzero. Two nonzero left forms would
have overlapping support, because 2(8-d)>8. Thus at most one left form
is nonzero as a polynomial, and it is proportional to P6.

This is a complete statement about **F7-rational coefficient tuples**,
not about geometric tuples over the algebraic closure. Any rational conic
must have bad reduction in this presentation at 7; bad reduction is not
a nonexistence proof. The support counts were independently enumerated
for all 48, 342 and 2,400 nonzero forms of degrees 1, 2 and 3.

### Four transverse patterns and genuine finite lifts

One precise bad-reduction chart is

\[
P_i=7R_i\ (1\leq i\leq4),\quad P_5=A,\quad
P_6=A+7^6B,
\]

where every displayed polynomial has degree at most two and A is monic.
Put delta=7^6. After dividing by 7^6, the exact equation is

\[
\sum_{i=1}^4R_i^6=
\sum_{j=1}^6\binom6j\delta^{j-1}A^{6-j}B^j.
\]

At A=t^2+1 modulo 7 its first condition is sum R_i^6=6A^5 B.
There are 58 projective quadratic forms including zero. An exhaustive
pair-sum join of 1,711 pairs finds nine four-term multisets divisible by
A^5. Five have every form divisible by A; the four others, up to individual
nonzero F7 scalars and permutation, are:

| Pattern | Four R_i modulo 7 | B modulo 7 |
|---|---|---|
| 1 | t; t^2-1; t^2+2t-1; t^2-2t-1 | 4(t^2+1) |
| 2 | t^2+3t+1; t^2-3t+1; t^2+2; t^2+4 | 3(t^2+1) |
| 3 | t^2+2t+3; t^2-2t+3; t^2+3t+5; t^2-3t+5 | 3(t^2+1) |
| 4 | t^2+t-1; t^2+3t-1; t^2-3t-1; t^2-t-1 | 3(t^2+1) |

The 13-by-17 Jacobian has rank 11 at each pattern. A selected path,
choosing zero for free new digits, lifts each full equation through 7^41.
This includes the correction terms beginning at 7^6, not merely its
unperturbed limit. Node BigInt independently checks the resulting original
sextic identities modulo 7^47 and confirms that none is an exact identity.

Bounded rational reconstruction with numerator and denominator at most
floor(sqrt(7^41/2)) reconstructs respectively 14,16,13,11 of the 17
coefficients. None reconstructs an exact curve. No exhaustion of other
free-digit paths, other bad-reduction charts, or all conics is claimed.
The stored prefixes are not a proof of an infinite 7-adic lift.

There is an additional warning about the chosen paths: patterns 2–4 keep
B=3A and hence P6=352948 P5. This ratio cannot hold at a rational ESOP6
point because P5/P6 is not 2-adically integral, whereas every coordinate
divided by P6 must be 2-adically integral. Changing free digits can change
this ratio; this does not close their residue patterns.

## 2. Cubic quotient: avoid the standard plane trap

Let F(U)=sum_{i<=5} U_i^3-U6^3. The standard paired planes, including
(a,-a,b,-b,c,c), lie inside sum_{i<=5} U_i=U6. Lines joining their
points stay there. A point with all Ui positive cannot satisfy both this
linear relation and F(U)=0: the cube of a positive sum exceeds the sum
of the individual cubes when at least two terms are nonzero.

Instead choose B=(a,-a,b,-b,c,c) and a tangent direction D satisfying
sum eps_i B_i^2 D_i=0. Define

\[
f=F(D),\quad k=\sum\epsilon_i B_iD_i^2,\qquad C=fB-3kD.
\]

Since F(lambda B+mu D)=3lambda mu^2 k+mu^3 f, F(C)=0 exactly.
Unequal a,b,c allow the construction to leave the standard hyperplane.
The search uses

\[
D=(c^2x,c^2y,c^2z,c^2r,0,a^2(x+y)+b^2(z+r)).
\]

The bounded C++ run uses 1<=a<=b<=4, 1<=c<=4, gcd(a,b,c)=1,
excluding a=b=c. Its direction vectors have entries in [-24,24], gcd=1
and first nonzero entry positive. Every positive cubic output is checked
exactly before square sieving. A documented bound keeps all intermediate
integers inside signed 128-bit arithmetic. Python uses arbitrary-precision
integers for the larger sampled runs. Sample outputs are independently
replayed by Node BigInt.

A primitive positive integer cubic vector lifts projectively to six
rational squares iff every coordinate is an integer square. Indeed the
common projective scalar must have even valuation at each prime, since
some coordinate of the primitive vector is a unit. A cubic point alone
therefore never counts as an ESOP6 candidate.

The large random run uses seed 20260905, a,b,c in [1,100], and four
direction parameters in [-10000,10000]. The additional divisibility-based
run uses a=1764i, b=1764j, c=1764k+1 with 1<=i,j,k<=8. This rule was
chosen to explore the concentrated divisibility pattern; it is not a
proved sufficient condition for square lifting.

## 3. Two squares by construction: exact closure of my equal-pair route

I imposed D1+D2=D3+D4=-1 and C6/C5=r^2. For fixed paired B the latter
condition is a conic in the two pair differences. Rational norm rotations
parametrize it. Thus the projectively normalized fifth and sixth cubic
coordinates are 1 and r^2; four square conditions remain.

Two instances were implemented:

* B=(1,-1,1,-1,2,2): with r=a/b and D=a^2-b^2, the conic is
  X^2+Y^2=2(-a^4+3a^2b^2-b^4).
* B=(1,-1,2,-2,3,3): its corresponding norm equation is
  X^2+Y^2=4500a^2b^2-2874(a^2-b^2)^2.

These were not just random cubic searches: both square ratios were checked
exactly for every positive output. However, the common pair-sum restriction
is fatal, as the following proof explains.

**Proposition.** No all-nonzero rational square lift on a tangent section
at (a,-a,b,-b,c,c), c!=0, can also satisfy C1+C2=C3+C4.

**Proof.** Suppose a lift has sixth-power coordinates x1,...,x4,z,w and
set S=x1^2+x2^2=x3^2+x4^2. The tangent-plane condition gives

\[
w^2-z^2=\frac{a^2+b^2}{c^2}S. \tag{1}
\]

Clear denominators in the sextuple and make it primitive. Modulo 9, exactly
one of the five left bases is coprime to 3, and w is coprime to 3. The
equal pair sums modulo 3 show that this one base must be z. Thus all four
xi are divisible by 3. The valuation of a sum of two rational squares at
3 is twice the smaller input valuation. Equal pair sums imply that the
minimum valuation k>=1 is attained in each pair. After dividing the equal
pair sums by 3^(2k), their residues force the same number of unit terms
in each pair. The number of xi attaining k is therefore two or four, so

\[
v_3\left(\sum x_i^6\right)=6k.
\]

For 3-adic units z,w,
v3(w^4+w^2z^2+z^4)=1. Factoring w^6-z^6 gives
v3(w^2-z^2)=6k-1, which is odd. But the right side of (1) has even
3-adic valuation: both a^2+b^2 and S are sums of two squares, and c^2
has even valuation. If a=b=0, equation (1) already gives the impossible
equality w^2=z^2. This covers the exceptional zero case as well.
Contradiction.

This closes the **entire equal-pair tangent restriction**, not all cubic
tangent constructions and not all ESOP6 points with equal pair sums.
The unpaired directions in section 2 remain outside this obstruction.

## 4. Degree-eight surface: a necessary square and a sparse certificate

The inherited residual is

\[
(P^6+R^6)/2=QN^5+(10/27)t^{12}Q^3N^3+(1/81)t^{24}Q^5N,
\quad Q=1+qt^2,
\]

with q>0 rational, degrees P,R<=7 and N<=8, constants P(0)=R(0)=N(0)=1.

**Necessary condition for a reduced degree-eight map: q is a rational square.**
At a root alpha of Q, either R(alpha)!=0 and (P/R)(alpha)^6=-1, or Q
divides both P and R. In the second case valuation of the residual forces
Q|N as well, so every cleared projective coordinate has the common factor
Q. In the first case the quadratic field Q(alpha) contains a sixth root
of -1. Such a root has order 4 or 12; order 12 has degree four over Q.
Therefore the field is Q(i), which is equivalent to q being a rational
square. This argument reduces q; it does not exclude q=s^2.

For q=s^2, Gaussian factorization motivates

\[
P=A-stB,\quad R=B+stA,\qquad P^2+R^2=Q(A^2+B^2).
\]

The executed sparse family was A=1+u t^6, B=1+v t^6, s!=0.
Eliminate N1,...,N8 recursively. Residual 9 is 732s^3(u-v)/5, so u=v.
Residual 10 then gives u=v=(3867794/5375)s^6. Substitution into residual
12 and residual 14 divided by s^2 gives respectively

\[
E=\frac{2(243601295510403s^{12}-722265625)}{3900234375},
\]
\[
G=\frac{2(3599809104248613s^{12}-2744609375)}{1300078125}.
\]

If A_E and A_G denote their s^12 coefficients, then

\[
A_G E-A_E G=-\frac{220092525210624}{144453125}\ne0.
\]

This is a characteristic-zero contradiction, including arbitrary rational
denominators (indeed it works over C). It also has unit Groebner basis over
F7, F11 and F13. It closes exactly the stated sparse family, not arbitrary
degree-six A,B or the original degree-eight problem.

## 5. Additional full-fourfold boundary system

Use four independent transverse polynomials Ai, rather than two repeated
coordinates, and write xi=tAi, z=M-k t^6/2, w=M+k t^6/2. Exact expansion
gives

\[
\sum_{i=1}^4 A_i^6=6kM^5+5k^3t^{12}M^3+(3/8)k^5t^{24}M.
\]

At Ai(0)=M(0)=1, k=2/3. This system preserves four independent tangent
directions. The identity is verified, but its coefficient equations were
not solved here; no boundary curve is claimed.

## 6. References, reproduction, and limits

The relevant square/cube construction source was checked directly:
[Bremner–Choudhry–Ulas, Sections 5–6](https://arxiv.org/html/1402.4583v1).
Its mixed-power surfaces and coefficients do not automatically have the
required five-versus-one signature. The control
(1-t-t^2)^3+(1+t-t^2)^3=2-2t^6 was verified exactly.
The historical claims about 200 CPU-years and every known rational point
were not independently audited in this work.

Reproduce the mathematical checks and finite conic calculations:

```sh
python3 tools/fourfold_symbolic_routes.py
python3 tools/fourfold_conic_F7.py
python3 tools/fourfold_conic_lift.py
python3 tools/fourfold_reconstruct.py
python3 tools/fourfold_degree8_sparse.py
python3 tools/fourfold_proof_checks.py
```

Reproduce searches (the first is the bounded census; the next two are
sampled and do not certify any height interval):

```sh
g++ -O3 -std=c++17 tools/fourfold_tangent_search.cpp -o /tmp/fourfold-tangent
/tmp/fourfold-tangent 24 4
python3 tools/fourfold_large_tangent.py --samples 1000000
python3 tools/fourfold_large_tangent.py --samples 1000000 --crt-seeds
node tools/verify_fourfold_routes.mjs
```

The two-square scripts remain for replay, not for extending a route now
excluded by the proposition. The integer timeout is a partial attempt;
restarting it would repeat its uncheckpointed filter construction.

**Still unresolved:** a rational conic on the full fourfold; global lifting
of the finite conic prefixes; a genus-one construction with all necessary
square conditions; unrestricted cubic tangent square lifts; the general
degree-eight surface identity; the independent four-direction boundary
system; and the ESOP6 counterexample itself. No background search is claimed
to continue after this work block.
