# The heuristic constant, and where a solution could live

Two questions the repository never answered quantitatively:

1. Under the standard random model, how many solutions of
   `a^6+b^6+c^6+d^6+e^6=f^6` should exist below a given height?
2. If a solution exists at reachable height, it must come from algebraic
   structure the random model does not see. Where can that structure live,
   and where is it already ruled out?

Both answers are unfavourable to every lane currently in the repository,
and they say precisely which lane is not yet closed.

## 1. The Hardy–Littlewood constant

For `a_1^k+...+a_n^k=f^k` the random model predicts that the number of
**primitive** solutions with `f<=F` is

```text
N(F) ~ (rho/n!) * prod_p sigma_p * F^(n-k+1)/(n-k+1)     (n-k+1 > 0)
N(F) ~ (rho/n!) * prod_p sigma_p * log F                  (n = k-1)
```

with real density `rho=(n/k)*Gamma(1+1/k)^n/Gamma(1+n/k)` and p-adic
densities `sigma_p = lim_j N*_j(p)/p^(nj)`, where `N*_j` counts solutions
modulo `p^j` that are not all divisible by `p`. Primitivity is essential:
one primitive solution at height `f` contributes `F/f` imprimitive
multiples, so the log law is a statement about primitive tuples only. For
`p` not dividing `k` the congruence is smooth on primitive tuples and
`j=1` is exact; for `p | k` the ratio is computed until stationary.

`tools/singular_series.py` computes every `sigma_p` with exact integer
convolution (the FFT branch for large primes was checked against the exact
branch on all primes below 100 and agrees to every digit) and records the
stabilisation sequence for `p=2,3`. Results for `p<=20000` are in
`results/heuristic_2026_09_11/`.

### The local factors for five sixth powers

| p | sigma_p | why |
|---:|---:|---|
| 2 | 0.625 = 5/8 | stationary from j=3 (checked to j=11) |
| 3 | 0.2469 = 20/81 | stationary from j=2 (checked to j=7) |
| **7** | **0.01071** | sixth powers are 0 or 1 mod 7; exactly one term may be a unit |
| 13 | 2.830 | sixth powers are 0, ±1 mod 13; favourable |
| 19 | 0.853 | |
| 31 | 0.739 | |
| 37 | 1.641 | |
| 43 | 0.824 | |
| others | within 16% of 1, converging like p^-2 | |

Partial products of `sigma_p` over `p<=997, 4999, 19997` are
`0.007153, 0.007314, 0.007341`. The Weil bound
`|sigma_p-1| <= (k-1)^(n+1)/k * p^-2` for `p` not dividing 6 bounds the
remaining tail factor by `exp(±0.013)`.

```text
rho   = 0.60879
prod  = 0.0073412
C     = rho * prod / 5!  =  3.72e-5
```

The catastrophe is `sigma_7`. The mod-7 rigidity that the repository uses
as a sieve is the same fact that makes solutions almost never exist: the
sieve is cheap because the haystack is empty.

### Calibration on equations with known solutions

The same program, same code path, on the cases where the answer is known:

| Equation | C | Predicted | Observed |
|---|---:|---|---|
| four fifth powers = fifth power (Lander–Parkin) | 4.95e-2 | 0.25 by 144, 0.56 by 85,359 | 3 known below 85,359 |
| seven sixth powers = sixth power | 1.78e-6, growth F^2 | **1.16 by 1,141** | smallest known solution is 1,141 |
| six sixth powers = sixth power | 5.77e-8, growth F^1 | 0.06 by 10^6, 1 at 1.7e7 | none known |

The seven-term prediction lands exactly on the known minimum. The
five-power case is within a factor of about five (favourable). The six-term
case predicts the first solution near height `2*10^7`; it is a cheap test
of this model that the repository could run on a GPU and that nobody has
to trust me on.

### What it says about five sixth powers

Expected primitive solutions below `F`:

| F | expected |
|---:|---:|
| 730,000 (published all-class frontier) | 5.0e-4 |
| 4,300,000 (this repository, one class only, times 1/25) | +7.6e-7 over the frontier |
| 10^8 (hard limit of the engine) | 6.9e-4 |
| 10^12 | 1.0e-3 |
| 10^100 | 8.6e-3 |

The expectation reaches 1 at `F = exp(1/C) ≈ 10^11,660`. Under the random
model the probability that any counterexample exists below `10^100` is
under one percent, and extending the concentrated sweep from 4.3M to the
engine's hard limit of 10^8 buys an expected `3.7e-5 * ln(23) / 25 ≈ 5e-6`
solutions. No amount of CPU, GPU, or bucket engineering changes this
exponent; every integer-search lane in ATTACK_MATRIX.md is chasing a
one-in-200,000 event.

Conclusion 1: **a counterexample at any reachable height would have to be
structural, not sporadic.** The Elkies solution to the quartic case was
structural: it came from an elliptic fibration on a K3 surface, not from a
random hit. The only question that matters is whether such a structure
exists here and where.

## 2. Where structure can and cannot live

Write `X` for the sextic fourfold `x_1^6+...+x_5^6=x_6^6` in P^5
(smooth, Calabi–Yau, `K_X=0`), and `S` for the repository's
repeated-coordinate surface `2X^6+2Y^6+Z^6=W^6` in P^3 (smooth,
`K_S=O(2)` ample, general type).

### 2a. The surface lane contradicts Bombieri–Lang

Every construction in GEOMETRIC_STRIKE.md, BOUNDARY_CONTACT_6.md and
RATIONAL_CURVE_ATTEMPT.md seeks a rational curve on `S`. The Bombieri–Lang
conjecture says rational points on a surface of general type are not
Zariski dense; Lang's stronger form says the rational and elliptic curves
on it are finitely many and contain all but finitely many rational points.
On a Fermat-type sextic the expected exceptional curves are the 108 lines,
none of which is defined over Q or meets the positive chamber. The
degree-six impossibility proof and the degree ≤3 exclusions are exactly
what this conjecture predicts, and it predicts the same outcome at degree
eight, ten, and every degree after that. Continuing the boundary-contact
programme is a bet against Bombieri–Lang. That bet can be made, but it
should be made knowingly; the repository currently presents this lane as
the most promising one.

### 2b. The two classical descents, and why the Elkies analogue is closed

A sixth power is a square of a cube and a cube of a square, giving two
maps from `X`:

```text
squaring:  X -> C,  [x] -> [x_i^2],   C: z_1^3+...+z_5^3 = z_6^3  (Fermat cubic fourfold)
cubing:    X -> Q,  [x] -> [x_i^3],   Q: y_1^2+...+y_5^2 = y_6^2  (smooth quadric fourfold)
```

Both targets are rational and have dense rational points. A solution is a
rational point of `C` whose six coordinates are squares up to a common
scalar, or of `Q` whose coordinates are cubes. The natural Elkies-style
move is to take a rational subvariety `V` of the target and lift it: the
sixth-power points above `V` form an abelian cover of `V`, and one wants
that cover to have Kodaira dimension at most one so that it can carry
infinitely many rational points.

**Planes in the cubic.** Above a plane `Λ ⊂ C` the sixth-power points form
`{x : x_i^2 = ℓ_i(s)}`, the intersection of three diagonal quadrics in P^5:
degree 8, `K = O(2*3-6) = 0`, a K3 surface. This is the exact analogue of
Elkies' K3, and it is the only lift of this type with Kodaira dimension
zero: above a rational surface of degree `e>=2` in `C` the cover has
canonical class `π*(K_V + 3h)` with `h` the hyperplane class, which is
positive (for a quadric surface `(1,1)`, for a cubic scroll `4`), so those
covers are of general type. Above a line the lift is the intersection of
four diagonal quadrics, a canonical curve of genus 5, with finitely many
rational points by Faltings and none known.

So the Elkies analogue lives on planes of `C`, and `C` is projectively
equivalent over Q to the Fermat cubic (`z_6 -> -z_6`). Degtyarev, Itenberg
and Ottem ([Planes in cubic fourfolds](https://arxiv.org/abs/2105.13951))
prove that a smooth complex cubic fourfold contains at most 405 planes and
that the Fermat cubic is the unique cubic attaining this bound. Its 405
planes are the standard ones: choose a pairing of the six coordinates into
three pairs and, on each pair, a relation `z_a = -ω z_b` (or `z_a = ω z_6`
for the pair containing index 6) with `ω^3=1`. At most one pair contains
index 6. A positive real point of `X` maps to a real point of `C` with all
`z_i > 0`; on a pair `(a,b)` inside `{1,...,5}` it would satisfy
`z_a = -ω z_b` with `ω = -z_a/z_b` real and negative, impossible for a
cube root of unity. Hence

> **no positive point of `X` lies above any plane of `C`.** The K3 lane
> of the squaring descent is empty, for every plane, over every field.

**Planes in the quadric.** `Q` contains two three-parameter families of
planes. Above a plane the cube conditions give
`{x : x_i^3 = ℓ_i(s)}`, the intersection of three diagonal cubics in P^5,
`K = O(3)`, general type; above a line, a curve of genus 10. Nothing of
Kodaira dimension ≤1 arises.

### 2c. The fourfold itself

On `S`, rational curves of degree `d` have `4(d+1)-4 = 4d` effective
coefficients against `6d+1` equations: negative expected dimension, which
is why every degree has been empty. On `X` a rational curve of degree `d`
has `6(d+1)-4 = 6d+2` coefficients against `6d+1` equations: expected
dimension one, matching the virtual dimension `dim X - 3 = 1` for rational
curves on a Calabi–Yau fourfold. Rational curves on `X` are therefore not
forbidden by any conjecture, and a rational curve over Q through the
positive chamber would give infinitely many solutions.

But a positive rational curve on `X` is not a route to a first solution:
its value at any rational parameter is already a solution. The curve
multiplies solutions; it does not create the first one. That first point is
exactly what the random model says does not exist at reachable height.

## 3. What this leaves

| Lane | Status after this note |
|---|---|
| Integer sweeps, any class, any engine | expected yield ≤ 10^-5 per e-fold of height; not worth running |
| Rational curves on `S` (boundary contact, any degree) | contradicts Bombieri–Lang; expected empty at every degree |
| K3 lifts of planes in the cubic descent | **closed**: no positive point lies above any of the 405 planes |
| Lifts of other rational surfaces or lines in either descent | general type or genus ≥5; no infinite family possible |
| Rational curves on `X` directly | not forbidden, but presuppose a first point |

The single lane that is neither closed nor conjecturally empty is:

> Find a surface `V ⊂ X` defined over Q, of Kodaira dimension at most one,
> not contained in any of the standard loci above, whose real locus meets
> the positive chamber.

That is a precise question about the Fermat sextic fourfold, and it is
where a structural counterexample would have to come from. It is also
where a proof of nonexistence of such surfaces would say, together with
Bombieri–Lang, that Euler's conjecture at `k=6` is "true in practice":
false only at heights no computation will ever reach.

## Reproduction

```bash
python3 tools/singular_series.py --k 6 --n 5 --pmax 20000 --output out.json
python3 tools/singular_series.py --k 5 --n 4 --pmax 20000 --output out5.json
python3 tools/singular_series.py --k 6 --n 7 --pmax 20000 --output out7.json
python3 tools/singular_series.py --k 6 --n 6 --pmax 20000 --output out6.json
```

Each run takes a few minutes on one core. The 2-adic and 3-adic
stabilisation sequences, every `sigma_p`, the partial products, and the
predicted counts are in the JSON.
