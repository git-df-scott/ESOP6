# Multiplicity, square-branch, and conic strike — 2026-09-05

**No ESOP6 counterexample was found. No rational curve was constructed.**

Base: `60d02436` (the full base SHA is recorded in the manifest), on
`codex/fourfold-routes-2026-09-05`, beyond `main`. This branch preserves the
inherited work and adds the following results. The requested ultimate target
remains unmet. In particular, the unrestricted degree-eight square branch
and the global rational-conic problem are still open.

## Results and exact coverage

Here the multiplicities describe five **positive** left-hand bases. A family
search allows further coincidences between its displayed variables; it does
not exclude repeated values in the final pair. Counts can therefore overlap
between families. The partition list has seven entries: the proposed list
omitted `5`.

| Pattern | Equation | Result in this work |
|---|---|---|
| 1+1+1+1+1 | five distinct bases | Not searched in this strike |
| 2+1+1+1 | 2a^6+b^6+c^6+t^6=f^6 | Complete bounded exclusion for primitive solutions with f <= 200,000 |
| 2+2+1 | 2a^6+2b^6+t^6=f^6 | Complete bounded exclusion for primitive solutions with f <= 50,000,000 |
| 3+1+1 | 3a^6+b^6+t^6=f^6 | Complete bounded exclusion for primitive solutions with f <= 500,000 |
| 3+2 | 3a^6+2b^6=f^6 | Impossible at every height, modulo 7 after primitive reduction |
| 4+1 | 4a^6+b^6=f^6 | Impossible at every height, by reduction to Fermat's theorem for cubes |
| 5 | 5a^6=f^6 | Impossible at every height, modulo 7 after primitive reduction |

Primitive reduction preserves multiplicities and decreases f. Thus the
bounded exclusions also exclude nonprimitive positive solutions within the
same bounds. They do **not** exclude larger solutions. Taken together, this
is a complete repeated-summand exclusion only through **f = 200,000**.
The deeper 50-million result applies only to `2+2+1`. No historical
all-class or distinct-only coverage claim is used in these conclusions.

These are computational results supported by a completeness argument and
positive/differential controls, not externally reviewed theorems. The full
production negative queries have not all been replayed by a second search
implementation.

## 1. Exact multiplicity reductions

Modulo 8, 9, and 7, a sixth power is 0 for a base divisible by 2, 3, or 7,
respectively, and 1 otherwise. In a primitive solution, f is a unit at each
of these primes: if f were divisible by one, counting at most five unit
summands would force all five bases to share that prime.

Exactly one left summand is therefore a unit at each of 2, 3, and 7.
A value occurring at least twice cannot carry any of those roles. Hence
every repeated base is divisible by 42. The roles must belong to singleton
slots. The numbers of labelled role assignments for the three surviving
patterns are 27, 1, and 8, respectively. Assigning all roles to one singleton
without proof would miss spread cases of `2+1+1+1` and `3+1+1`.

Patterns `3+2` and `5` have no singleton slot, so they are impossible.

For `4+1`, make (a,b,f) primitive. The role count makes b and f odd. Any
prime dividing both b and f also divides a, so gcd(b,f)=1. Consequently

    U=(f^3+b^3)/2,  V=(f^3-b^3)/2

are positive coprime integers, and UV=a^6. Unique factorization of integers
gives U=r^6 and V=s^6, with r,s positive. But then

    f^3=r^6+s^6=(r^2)^3+(s^2)^3,

contrary to Fermat's theorem for exponent 3. This is an all-height reduction,
not a finite congruence experiment. An independent proof of the standard
cubic theorem is available in [Joseph Lipman's Purdue notes](https://www.math.purdue.edu/~jlipman/MA598/x%5E3%2By%5E3%2Bz%5E3.pdf).

## 2. Search algorithm, completeness, and retained counts

Implementation: `tools/repeated_campaign.cpp`. It uses exact integers,
small residue masks, sorted tables or a monotone two-pointer join, and exact
replay. No floating-point root estimates or probabilistic membership filters
are used. Repeated pairs are explicitly retained.

### The `2+2+1` branch

The unique singleton t carries all three roles. Write a=42A and b=42B.
Then t and f are units modulo M=42^6=5,489,031,744 and

    T=(f^6-t^6)/(2M)=A^6+B^6.

There are 144 sixth roots of unity modulo M, independently verified by CRT
from 4 roots modulo 64, 6 modulo 729, and 6 modulo 117649. For every t<H
and each root r, enumerate every f=r*t+kM with t<f<=H. Test divisibility by
2M, sieve T, and solve A^6+B^6=T with 1<=A<=B<=floor(H/42).

The two-pointer solver starts at (1,floor(T^(1/6))), with the upper bound
clipped to H/42. A sum below T discards the current A; a sum above T discards
the current B. On equality it records the pair and advances both ends.
Strict monotonicity proves that every possible pair is considered or
excluded, including A=B. The sixth root is computed by integer bisection.

### The spread-capable `3+1+1` and `2+1+1+1` branches

Designate t as the unique singleton not divisible by 7. Every other base is
divisible by 7, and the repeated base a is divisible by 42. Enumerate all
six rays f=r*t+k*7^6, with t<f<=H and f coprime to 42. For each, enumerate
positive a=42j until the exact residual is too small.

For `3+1+1`, test whether f^6-t^6-3a^6 is a positive sixth power of a
multiple of 7. For `2+1+1+1`, query the exact sorted table b^6+c^6 for
positive multiples of 7 with b<=c<H. Its four tables distinguish the
remaining numbers of units at 2 and 3. The required budgets are

    1-[2 does not divide t],  1-[3 does not divide t].

This covers every assignment of the 2 and 3 roles. The pair-table cutoff
b^6+c^6+2*42^6+1<=H^6 is necessary, since a>=42 and t>=1. It cannot remove
a solution. All exact equal-key table entries are checked.

Before the terminal oracle, exact sumsets are used modulo
64, 729, 13, 19, 31, 37, and 43. These are necessary filters, not an
assertion of local sufficiency.

| Final run | (f,t) pairs | Residual queries | Passed final local masks | Exact hits |
|---|---:|---:|---:|---:|
| 2+2+1, f<=50,000,000 | 8,572,249 | 4,284,458 | 14,521 | 0 |
| 3+1+1, f<=500,000 | 1,752,088 | 11,312,407,443 | 150,404 | 0 |
| 2+1+1+1, f<=200,000 | 264,875 | 740,604,313 | 6,567,739 | 0 |

The largest pair table contains 163,899,685 entries across the four budgets.
Recorded elapsed times were about 60, 129, and 54 seconds respectively;
the jobs overlapped and these timings are not benchmark claims. Smaller
pilot runs are retained separately and must not be added to the final runs.

### Arithmetic safety and verification

At f=50 million, f^6 does not fit in 128 bits. The concentrated engine never
forms it. It factors the difference into

    (f-t)(f+t)(f^2+ft+t^2)(f^2-ft+t^2),

cancels 2M with gcds **before multiplication**, and only forms the quotient.
The CLI caps this family at f=100 million. Then T<10^48/(2M)<10^38,
2T<2^128, and every intermediate partial product is no larger than the final
positive product. Each individual quadratic factor fits in 64 bits. The
unscaled engines are capped at f=1 million, where even five sixth powers
fit in 128 bits. CRT products fit in 64 bits at all supported bounds.

Controls include 2,850 injected positive pair representations, all 150,000
targets in a small differential domain, 20,000 positive sieve controls,
independent CRT enumeration, 2,000 Python big-integer quotient checks up to
100 million, large positive pair controls through base 2 million, and an
independently enumerated small pair table. The first development root-list
self-test caught an overflowing raw modular-power calculation; it was fixed
before any successful `2+2+1` production run. No failed run is counted as
coverage. A proposed wide-integer dependency was unavailable; the final
implementation uses the proved factored quotient and no external C++ library.

## 3. Degree eight: the square branch is narrowed, not finished

Retain the exact inherited identity, constants P(0)=R(0)=N(0)=1, and degrees
P,R<=7 and N<=8:

    (P^6+R^6)/2 = Q N^5 + (10/27)t^12 Q^3 N^3 + (1/81)t^24 Q^5 N,
    Q=1+q*t^2.

For a reduced rational map q must be a positive rational square. Writing
q=s^2 permits the full Gaussian factorization P=A-stB, R=B+stA, with
independent degree-at-most-six A,B after choosing the orientation. No sparse
or equal-pair constraint is justified by this factorization.

**New restricted obstruction:** the displayed normalized identity has no
solution with **all coefficients P,R,N,q integral at 7**. This is true even
without assuming q is a square.

To prove it, reduce the cleared sextic identity
2(tP)^6+2(tR)^6+(N-t^6Q/3)^6=(N+t^6Q/3)^6 modulo 7. At every nonzero
t in F7, the repeated-slot count forces P(t)=R(t)=0. Therefore

    P=(1-t^6)(1+a*t), R=(1-t^6)(1+b*t) over F7.

For each q,a,b in F7, coefficients 1 through 8 uniquely determine N:
the multiplier of its next coefficient is 5, invertible modulo 7. Check
every remaining coefficient. All **343** cases fail: 246 first fail at
degree 9, 96 at degree 10, and one at degree 12. The square-q subset has
196 cases and is included. A retained certificate gives N and a nonzero
residual coefficient for every case; SymPy independently substitutes them.

Thus a hypothetical rational solution in this normalization must have at
least one coefficient with a denominator divisible by 7. Clearing such
denominators destroys the unit-constant normalization used in the argument.
It would be invalid to call this a proof for arbitrary rational coefficients.
That bad-reduction branch still needs a Gauss-valuation/Newton-polygon
analysis. Cancellation and the unrestricted rational-function route are
not closed by this new calculation.

## 4. Conics: simultaneous local checks and a globalization control

For six degree-at-most-two polynomials, we enumerated the distinct sixth
power polynomials at the indicated coefficient and equation moduli. A
pair/triple exact join solves the coefficient identity, with a primitive coefficient tuple at the relevant prime. At 2, 3, and 7,
evaluation and the degree bound force the right polynomial to be primitive
as well; at 5 the census also includes zero right-hand reduction.

| Prime | Coefficients modulo | Equation modulo | Distinct sixth polynomials | Solution multisets |
|---|---:|---:|---:|---:|
| 2 | 4 | 8 | 29 | 28 |
| 3 | 3 | 9 | 14 | 13 |
| 5 | 5 | 5 | 63 | 80,259 |
| 7 | 7 | 7 | 58 | 57 |

The mod-8 sixth polynomial depends only on coefficients modulo 4, and the
mod-9 polynomial only on coefficients modulo 3, by the binomial expansion.
At 2, 3, and 7, the solutions have just one nonzero sixth polynomial on the
left, matching the right. Hence every primitive integral conic tuple reduces
to an axis at each of those primes; the surviving coordinate may differ
between primes. Direct five-term enumeration independently reproduces the
mod-8 and mod-9 counts. At 5, ranks 1, 2, and 3 occur: that prime does not
force the same collapse.

This **does not** prove absence over Q: bad reduction is allowed. It also
does not make a Hasse principle available for the parameter space of conics.
The word "conic" describes the proposed curve on the fourfold, not the
high-degree coefficient system being solved.

For each of the four retained 7-adic prefixes, CRT combines its 17
coefficients modulo 7^41 with an exact axis solution modulo 8, 9, and 5.
The resulting integer quadratic tuples satisfy the original sextic identity
modulo **8, 9, 5, and 7^47 simultaneously**, but none is an exact identity.
Node BigInt independently verifies all four tuples. Bounded rational
reconstruction modulo 360*7^41 reconstructs respectively 9, 10, 5, and 4
of the 17 coefficients; none gives a curve. This is a concrete demonstration
that adding several primes can still leave false global candidates. It
neither excludes other free-digit paths nor proves a rational point exists.

## 5. Unequal-pair square-cube construction and its exact veto

Take the cubic tangent seed B=(1,-1,2,-2,3,3), with independent pair sums
u,v. Put

    D=((u+x)/2,(u-x)/2,(v+y)/2,(v-y)/2,0,h), h=(u+4v)/9,
    F=sum eps_i*D_i^3, k=sum eps_i*B_i*D_i^2, C=F*B-3k*D.

Then C is exactly on the cubic quotient. Imposing C6/C5=r^2 gives a conic
in x,y; it forces two projective square ratios, C5/C5=1 and C6/C5=r^2,
while leaving four ratios to lift. Unlike the earlier family, this allows
C1+C2 and C3+C4 to differ. Their ratio is u/v.

For the diagnostic choice u=1,v=182,h=81, r=a/b and d=a^2-b^2, the conic is

    X^2+182Y^2 = 9[d^2(4h^3-u^3-v^3)/3+4a^2b^2h^3],
    X=3d*x+2h*b^2, Y=3d*y+4h*b^2.

The 24 ratios a/b were first required to satisfy a^6=b^6 modulo 42^6.
Integer norm representations and rational norm rotations produced 167,040
directions and 108 positive cubic outputs. The two forced ratios and cubic
identity were checked exactly; no full square lift passed. These are
diagnostic outputs, not a height census or new locally viable frontier.

**A subsequent exact norm test closes this chosen family entirely.** If its
four remaining ratios were rational squares, u/v would be a ratio of two
nonzero sums of two rational squares. At every prime p=3 mod 4, such a ratio
has even valuation: a sum of two squares has valuation twice the minimum
input valuation, since -1 is not a square modulo p. But v7(1/182)=-1.
Contradiction. This closes the unequal ratio 1:182, not all unequal pairs.

The executable now applies this all-height gate before any production
search; `--replay-rejected-control` reproduces the diagnostic run. Future
unequal-pair selections must first pass this norm-ratio test, in addition
to the RHS-ratio congruences. RHS admissibility alone was insufficient.

## Reproduction and remaining work

From this repository root:

```sh
g++ -O3 -std=c++17 -Wall -Wextra tools/repeated_campaign.cpp -o /tmp/repeated_campaign_final
/tmp/repeated_campaign_final --self-test
/tmp/repeated_campaign_final 221 50000000
/tmp/repeated_campaign_final 311 500000
/tmp/repeated_campaign_final 2111 200000
python3 tools/multiplicity_algebra.py
python3 tools/conic_crt_reconstruction.py
/tmp/repeated_campaign_final --r-seeds > results/multiplicity_square_2026_09_05/square_cube_ratio_seeds.csv
python3 tools/unequal_square_cube.py
# Optional historical diagnostic, already proved incapable of a lift:
python3 tools/unequal_square_cube.py --replay-rejected-control
python3 tools/verify_multiplicity_campaign.py /tmp/repeated_campaign_final
node tools/verify_multiplicity_bigint.mjs
```

The 200k `2+1+1+1` table uses several GB of RAM. Search
counts, certificates, prefixes, samples, verification records and hashes
are in `results/multiplicity_square_2026_09_05/`. Existing campaign results
and source engines were not modified. No background search remains active.

What remains unresolved: a counterexample; repeated families above the
stated bounds; the 7-adic denominator branch of degree eight; a global
rational conic; and an unequal-pair square-cube family passing all the
derived local tests and actually lifting. No finite search result here
settles those questions.
