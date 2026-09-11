# Astra third strike — the degree-eight family

```text
ESOP6 SOLUTION FOUND: NO
POSITIVE RATIONAL SURFACE POINT: NO
DEGREE-EIGHT IDENTITY: NOT FOUND
DEGREE-EIGHT OBSTRUCTION: NOT OBTAINED
DEGREE-EIGHT TOP EQUATION (deg N=8): EQUIVALENT TO ESOP6
DEGREE-EIGHT 3-INTEGRAL CASE: IMPOSSIBLE
DEGREE-EIGHT REMAINING BRANCH (deg N=2): LIVE
UNEQUAL-SLOPE DEGREE SIX, deg M=0: IMPOSSIBLE
UNEQUAL-SLOPE DEGREE SIX, deg M=6: OPEN
```

No six-positive-integer candidate was produced, so neither verifier was run
on a candidate. Nothing here promotes a near miss, a formal series, or a
bounded empty search to a solution.

## 1. Result

The designated single next construction — RATIONAL_CURVE_ATTEMPT.md equation
(5), the degree-eight family with a centered positive quadratic denominator —
was set up exactly and attacked. It was neither solved nor closed. Two things
were settled, and they change what the construction is worth.

1. **The construction is not a reduction.** Its main branch
   (`deg N = 8`, `max(deg P7,deg R7) = 7`) has as its top-degree coefficient
   equation

   ```text
   (p7^6+r7^6)/2 = n8 q (27 n8^2+q^2)(3 n8^2+q^2)/81,
   ```

   and the substitution `X=p7, Y=r7, Z=n8-q/3, W=n8+q/3` turns that equation
   into `2X^6+2Y^6+Z^6=W^6` itself. A rational solution of the *single leading
   equation* of the family is therefore already a sixth-power counterexample
   (five terms when `p7 r7 != 0`, the stronger three-term relation
   `2X^6+Z^6=W^6` when exactly one of them vanishes). Solving that branch is
   at least as hard as the original problem, and any obstruction to it would
   be an obstruction to ESOP6 itself. This is a structural fact about the
   construction, verified symbolically, not an opinion about difficulty.

2. **The degree-six certificate strategy cannot be transplanted.** The
   3-adic reduction of BOUNDARY_CONTACT_6.md rescales the reversed data and
   compares with an unperturbed equation modulo a high power of 3. Here the
   corresponding unperturbed equation `(x^6+y^6)/2 = n^5 q~` has the exact
   rational solution family

   ```text
   (x, y, n, q~) = (s^5 q~, s^5 q~, s^6 q~, q~)   for every Q0,
   ```

   so no amount of lifting of that reduction can ever produce a
   contradiction. In degree six the analogous family was killed by two
   normalizations that are not available here (one Möbius freedom is already
   spent centering `Q`). Building the missing non-degeneracy input is exactly
   where this strike stalls; see §5.

The remaining branch, `deg N = 2`, is not self-referential and is the new
single next construction (§7).

## 2. Exact scope

```text
Q = 1+q t^2, q in Q, q>0,
X = t P7/Q, Y = t R7/Q, Z = T8/Q, W = Z+(2/3) t^6,
P7(0)=R7(0)=T8(0)=1, deg P7,deg R7 <= 7, deg T8 <= 8,
N = T8+t^6 Q/3,
(P7^6+R7^6)/2 = Q N^5+(10/27) t^12 Q^3 N^3+(1/81) t^24 Q^5 N.        (5)
```

The linear coefficient of `(P7+R7)/2` was **not** set to zero: the single
Möbius freedom is used to center `Q`, exactly as RATIONAL_CURVE_ATTEMPT.md
requires. P7 and R7 are independent throughout. Nothing below closes general
rational curves, other tangent directions, larger degrees, or the positive
rational-point problem.

## 3. The system and its reversed normal form

Equation (5) was re-derived from `2X^6+2Y^6+Z^6-W^6` in exact rational
arithmetic; the surface identity is

```text
2X^6+2Y^6+Z^6-W^6 = 4 t^6 (LHS(5)-RHS(5)) / Q^6.
```

Inverting the parameter (`x = s^7 P7(1/s)`, `y = s^7 R7(1/s)`,
`n = s^8 N(1/s)`, `q~ = s^2 Q(1/s) = s^2+q`) gives an exact equivalence,
coefficient by coefficient (`81*[t^d] = [s^{42-d}]`, all 43 orders):

```text
(81/2)(x^6+y^6) = n q~ (3n^2+q~^2)(27n^2+q~^2),                       (R)
x, y monic of degree 7, n monic of degree 8, q~ = s^2+q, q>0.
```

This is the degree-eight analogue of equation (3) of BOUNDARY_CONTACT_6.md,
with `t^12` replaced by `q~^2`. It is the form used for every arithmetic
argument below.

**Shape of the residual system.**

| form | unknowns | equations | degrees |
|---|---:|---:|---|
| (5), coefficients of `t^1..t^42` | 23 (`q`, `p1..p7`, `r1..r7`, `n1..n8`) | 42 | order-`d` equation has total degree `<= d` |
| after eliminating `n1..n8` by the increasing-order recurrence (constant derivative in N is 5) | 15 (`q`, `p1..p7`, `r1..r7`) | 34 (orders 9..42) | as above |
| (R), coefficients of `s^0..s^41` | 23 | 42 | every equation has total degree exactly `<= 6` |

The `s^42` coefficient of (R) and the `t^0` coefficient of (5) are automatic.
For orders 1..11 the equation is weighted homogeneous of degree `d` under
`deg p_i = deg r_i = deg n_i = i`, `deg q = 2`; from order 12 on it splits
into weighted-homogeneous pieces of degrees `d`, `d-12`, `d-24`, because the
explicit `t^12` and `t^24` break the grading.

Full symbolic elimination of `n1..n8` in the 15 remaining unknowns was
attempted in a sparse rational polynomial ring. It solved `n1..n6`
(`n6` in 80 s) and did not finish `n7` within 20 minutes; it was killed, as
the time budget requires. Nothing below depends on it: the reversed form (R)
carries the same information without any elimination.

## 4. Degree analysis and the two branches

In (5) the left side has degree `6*rho` with `rho = max(deg P7,deg R7)`; the
leading coefficients of `P7^6` and `R7^6` are real sixth powers and cannot
cancel. On the right, the three terms have degrees `2+5h`, `18+3h`, `34+h`
with `h = deg N`, and their leading coefficients all carry the sign of the
leading coefficient of `N` (because `q>0`), so they cannot cancel either.
Hence `6 rho = max(2+5h, 18+3h, 34+h)` with `rho <= 7`, `h <= 8`:

| `h = deg N` | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|---|
| degree of the right side | 34 | 35 | 36 | 37 | 38 | 39 | 40 | 41 | 42 |
| admissible | | | **rho=6** | | | | | | **rho=7** |

Only `(h,rho) = (8,7)` and `(h,rho) = (2,6)` survive.

* `(8,7)`: top equation `(p7^6+r7^6)/2 = n8 q (27n8^2+q^2)(3n8^2+q^2)/81`,
  which is the surface equation (§1, item 1). It is unobstructed at 2 and at 3.
  Writing `d = v3(q)-v3(n8)` and `e = v2(q)-v2(n8)`, the admissible classes
  are `d = 1`, `d == 2 (mod 6)` for `d<=0`, `d == 0 (mod 6)` for `d>=2`; and
  `e == 0 or 5 (mod 6)` for `e>0`, `e == 0 or 1 (mod 6)` for `e<0`, with
  `e = 0` excluded. A bounded integer search over
  `1 <= n,q <= 150`, `0 <= p,r <= 250` found no point, and a wider manual run
  over `1 <= n,q <= 300`, `0 <= p,r <= 400` found none either. That is a
  bounded empty box, not a nonexistence claim, and it is of course also a
  (tiny) ESOP6 search.
* `(2,6)`: top equation `(p6^6+r6^6)/2 = q^5 n2/81`, which is solvable over Q
  for any `p6,r6,q` by choosing `n2`. This branch is the live one.

## 5. The 3-adic work, and exactly where it stalls

Write `w` for the 3-adic Gauss valuation (minimum of the coefficient
valuations), `a = min(w(x),w(y))`, `b = w(n)`, `c = w(q~) = min(0,v3(q))`.
A sum of two sixth powers has valuation exactly `6*min` (otherwise `-1` would
be a square in `F3(s)`), and the two quadratic factors of (R) have exact
valuations because `1+2b` and `3+2b` are odd while `2c` is even. So (R)
forces

```text
4+6a = b+c+min(1+2b,2c)+min(3+2b,2c).                                 (V)
```

**Theorem (proved).** (R) has no solution with 3-integral `x,y,n,q~`. They
are monic, so `a=b=c=0`, and (V) reads `4 = 0`. Equivalently, modulo 3 the
identity (R) reads `0 = n q~^5` with both factors monic. Every solution must
therefore have coefficients with denominators divisible by 3, and in
particular `v3(q) < 0` or `n` has a 3-denominator. The first obstruction
modulus for the integral case is **3**.

Substituting `s -> 3^-k s` and scaling by `3^(42k)` turns (R) into

```text
81(x0^6+y0^6)/2 = 81 n0^5 q~0 + 30*3^(12k) n0^3 q~0^3 + 3^(24k) n0 q~0^5,
q~0 = s^2+3^(2k) q,
```

and for `k` large enough every member of the rescaled data is integral and
monic, so every solution satisfies the **unperturbed** congruence

```text
(x0^6+y0^6)/2 = n0^5 q~0   (mod 3^(12k-4)),  12k-4 >= 8.               (U)
```

Two exhaustive finite certificates for (U) were computed (standard library
only, exact integer arithmetic):

| stage | enumerated | surviving |
|---|---:|---:|
| residue shapes modulo 3 | 19,683 `(n,q~)` pairs, of which 27 are cube-admissible; 2,187 monic `x` tested for each of those (59,049 tests) | **95** shapes |
| exact linear lifting to modulus 9 | 95 shapes, rank 8 in the 9 digit unknowns of `(n,Q0)` | **27** shapes, 81 states |

Modulo 3, (U) is `2 S^3 = n^5 q~` with `S = x^2+y^2`, so `n^5 q~` must be a
cube and `S` its cube root; anisotropy of `x^2+y^2` at odd-degree primes
prunes the rest. Modulo 9 the left side still depends only on `x,y` mod 3,
so the step is a 42-by-9 linear system over `F3`; its rank is 8 for every
one of the 95 shapes.

**This cannot be pushed to a contradiction, and the reason is structural.**
(U) has the exact rational solutions `(s^5 q~, s^5 q~, s^6 q~, q~)` for every
`Q0` — three of the 27 survivors are exactly their residues — so the lifting
tree never empties. In BOUNDARY_CONTACT_6.md the analogous exact solutions
`(s^5, s^5, s^6)` were removed by two normalizations (the vanishing `s^4`
coefficient of `x+y`, and the unit constant coefficient of `m0` coming from
a tight Newton-polygon rescaling). Here the first is unavailable (the Möbius
freedom is spent on centering `Q`) and the second was not established: the
weighted Gauss valuations `w_r` for the degree-eight family admit three
regimes of (V)

```text
R1: c <= b,      4+6a = b+5c
R2: c = b+1,     a = b
R3: c >= b+2,    6a = 5b+c
```

and I did not reduce them to a single tight rescaling as the degree-six proof
does. **That is the precise stall point.** Every exact solution of (U)
exhibited here has `x = y`, i.e. `P7 = R7`, i.e. `X = Y`, which is separately
closed by the Fermat-cubes argument of GEOMETRIC_STRIKE.md §3. Whether the
solution set of (U) is exactly that locus is unknown to me; if it were, the
lifting certificate would close. That is a conjecture here, not a result.

A further honest limit: even if the shapes were pruned, exhaustive lifting
past modulus 9 is out of reach in this parametrization. The exact Jacobian
modulo 3 of the 42 equations in the 23 digit unknowns has rank 8 or 12 on the
27 survivors, so a single further level would branch by `3^15` or `3^11` per
state. No such enumeration was run and none is claimed.

2-adic and 7-adic analogues were not built for the full system. The 2-adic
analysis of the `deg N = 8` top equation (§4) is recorded and is
inconclusive, as expected for an equation equivalent to ESOP6.

## 6. The unequal-slope degree-six sibling

Scope: `X = t A5(t)`, `Y = rho t B5(t)`, `rho in Q`, `rho>0`, `rho != 1`,
`A5(0)=B5(0)=1`, `deg A5,deg B5 <= 5`, `Z(0)=1`, `deg Z <= 6`,
`W = Z + c t^6`.

* The leading contact constant is no longer `2/3`: it is `c = (1+rho^6)/3`.
* With `M = (W+Z)/2`, `g = A5`, `f = rho B5` the exact divided identity is

  ```text
  g^6+f^6 = (c/16) M (12M^2+c^2 t^12)(4M^2+3c^2 t^12).                 (6)
  ```

  For `rho=1`, `c=2/3`, this is `2/81 * M (3M^2+t^12)(27M^2+t^12)`, i.e.
  equation (1) of BOUNDARY_CONTACT_6.md.
* The degree argument survives unchanged: with `h = deg M`,
  `r = max(deg g,deg f)`, `6r = max(5h, 12+3h, 24+h)` leaves exactly
  `(h,r) = (0,4)` and `(h,r) = (6,5)`.
* **Key lemma.** `v3(c) = min(0, 6 v3(rho)) - 1`, which is `5 mod 6` for
  every nonzero rational `rho`.
* **Theorem (new).** The branch `deg M = 0` is impossible for **every**
  rational `rho>0`. The valuation balance is
  `6 min(w(g),w(f)) = nu + min(2nu,1) + min(1+2nu,0) = 5 nu + 1` with
  `nu = v3(c) <= -1`, and `5 nu+1 == 2 (mod 6)` is never divisible by 6.
  For `rho=1` this specializes to the known `a4^6+b4^6 = 2/81`,
  `v3 = -4` obstruction, so the equal-slope result is recovered as a case.
* In the branch `deg M = 6`, reversing the parameter gives
  `x^6+y^6 = (c/16) m (12m^2+c^2)(4m^2+3c^2)` with `x` monic of degree 5,
  `y` of degree 5 with leading coefficient `rho`, `m` monic of degree 6. The
  case `b = w(m) = 0` (3-integral reversed data) is impossible for every
  rational `rho>0`, by the same `5 nu+1` computation. The surviving
  valuations in a bounded window are `b = -1` and `b == 0 (mod 6)`, matching
  the equal-slope case exactly.
* **Not done.** No analogue of the modulus-729 certificate for `b<0` was
  built, so the unequal-slope degree-six family is **not** closed. The
  correct normalization to use for `rho != 1` (the Möbius freedom is still
  available there) was not fixed either.

## 7. Single next construction

Attack the `deg N = 2` branch of (5), which is the only branch of the
designated family that is not equivalent to ESOP6 itself. In the reversed
variables it loses its forced `s^6` factor and becomes, exactly,

```text
(81/2)(x^6+y^6) = mu q~ (3 s^12 mu^2+q~^2)(27 s^12 mu^2+q~^2),         (8)
x, y monic of degree 6,  mu = s^2+n1 s+n2,  q~ = s^2+q,  q in Q_{>0}.
```

15 unknowns, 36 coefficient equations (degrees 0..35; degree 36 is
automatic). Its top equation `(p6^6+r6^6)/2 = q^5 n2/81` is solvable over Q,
so unlike the `deg N = 8` branch it is not self-referential. Its 3-integral
case is impossible by the theorem of §5, which applies verbatim.

Two further specific things are worth doing before any coefficient
enumeration, and each is smaller than a full certificate:

1. finish the weighted Newton-polygon analysis of (V) for (8) — three
   regimes, one free constant `q` — to obtain a *tight* rescaling, which is
   the exact ingredient missing in §5;
2. decide whether `x = y` (equivalently `X=Y`, already closed by Fermat for
   cubes) is the only solution locus of the unperturbed equation (U); if so,
   the lifting certificate closes.

Do **not** attack the `deg N = 8` branch as a construction: §1.1 shows that
finding its top coefficients is already finding a counterexample.

## 8. Reproduction

```bash
make third-checks
```

which runs

```bash
python3 tools/degree_eight_strike.py
python3 tools/unequal_slope_six.py
python3 tests/third_strike.py
```

Retained evidence and hashes are in `results/astra_third_2026_09_11/`
(`manifest.json` lists the sha256 of every file; `third_gates.json` and `inherited_gates.log` record the gate runs). The wider bounded search of
§4 is

```bash
python3 tools/degree_eight_strike.py --nq-bound 300 --pr-bound 400
```

Inherited gates re-run on this tree at the start of the strike, both passing:

```bash
make boundary-checks
make direct-checks
```

No integer sweep was run or extended (`caseA2`/`caseA3` beyond the existing
700k–730k control were explicitly declined). No Gröbner basis was attempted.
No earlier document, result directory, or `STATUS.md` was modified.

## 9. What is not claimed

* No ESOP6 counterexample, no positive rational surface point, no rational
  curve, no elliptic residual.
* The degree-eight family (5) is **not** proved impossible. Only its
  3-integral case is, plus the structural facts of §1.
* The unequal-slope degree-six family is **not** closed; only its
  `deg M = 0` branch and the 3-integral case of its `deg M = 6` branch are.
* The 95 mod-3 shapes and the 27 mod-9 survivors are a certificate about the
  unperturbed congruence (U), not about (5).
* Bounded empty searches are recorded as bounded empty searches.
