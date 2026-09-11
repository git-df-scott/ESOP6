# Where a counterexample would live, and what a CPU-hour buys

Auditor's analysis, 2026-09-11. Everything here is heuristic except the two
lattice bounds at the end.

## 1. The heuristic count

For a^6+b^6+c^6+d^6+e^6 = f^6, model the left side as a random integer near
f^6.  Solutions with f in [F, 2F] are then expected to number about

    N([F,2F]) ~ sigma * V * (5/6) * log 2,

where V = Gamma(7/6)^5 / Gamma(11/6) = 0.7306 is the volume factor and sigma
is the singular series (product of local densities).  Computed exactly by
convolution of the sixth-power residue distributions:

| modulus | local density |
|---|---|
| 8 | 0.750 |
| 9 | 0.259 |
| 7 | 0.0108 |
| 13 | 2.830 |
| 19 | 0.853 |
| 31 | 0.739 |
| 37 | 1.641 |
| 43 | 0.824 |
| 61..127 (p = 1 mod 6) | 0.97..1.15 |
| p = 5 mod 6 | 1 + O(1/p^2) |

Product over the listed moduli: sigma ~ 0.0070.  The 7-adic factor alone
(0.0108) is why the conjecture has held so long at k = 6: exactly one of the
five terms may be coprime to 7.

Consequences:
- expected solutions with f <= 4.3e6 (all classes): about 0.065
- expected solutions with 4.3e6 < f <= 1e8 (all classes): about 0.013
- per e-fold of height, all classes: 0.0042
- the count is LOG-distributed in f: doubling the frontier is worth the same
  everywhere.  There is no "sweet spot" in height.
- Compare k = 5: known (5,1,4) solutions at f = 144, 14132, 85359 are spaced
  roughly log-uniformly, as this model predicts.

## 2. Mass of the sub-searches

Class 1 (all three exemptions on one term) is 1/25 of the mass.  Inside class 1,
the count k7 = m mod 7 of small bases coprime to 7 splits it further:

| k7 | mass within class 1 | cost per candidate |
|---|---|---|
| 0 | 1/2401 | trivial (7^6 divides m; provably f >= 11,936,849) |
| 1 | 24/2401 ~ 1% | O(B^2/7^6) |
| 2 | 216/2401 ~ 9% | O(B^2/7^6) |
| 3 | 864/2401 ~ 36% | O(B^2/1300) leaves, table-free (SPEC_V3) |
| 4 | 1296/2401 ~ 54% | O(B^2/1300) leaves, table-free (SPEC_V3) |

(B = (f-1)/42.  Masses use independent divisibility; the exact binomial
weights are C(4,k)(6/7)^k(1/7)^(4-k).)

## 3. Expected value per plan, from the frontier of 4.3e6

EV = 0.0042 * (mass fraction) * ln(F_new / F_old).

| plan | mass | reach | EV | cost on this box |
|---|---|---|---|---|
| class 1, k7 <= 2 (k7engine) | 0.04 * 0.10 | 2.0e7 | 6e-6 | 0.5 h |
| class 1, all k7, old engine to 4.6e6 | 0.04 | 4.6e6 | 1.1e-5 | 4-8 h, RAM-capped |
| class 1, all k7, table-free v3 to 8e6 | 0.04 | 8e6 | 1.0e-4 | ~6 h (model) |
| class 1, all k7, v3 to 1e7 | 0.04 | 1e7 | 1.4e-4 | ~14 h (model) |
| class 4 from 3e6 to 3.5e6 | 0.16 | 3.5e6 | 1.0e-4 | ~4 days |
| class 5 beyond 5e6 | 0.48 | -- | -- | GPU territory |

Reading: every plan available in a session is worth about one chance in ten
thousand.  The table-free engine is the best of them because it removes the RAM
wall and is ~10x cheaper per leaf, not because it changes the odds by much.

## 4. Why "cheap deep slivers" do not help

A stratum where all four small bases share an extra prime p has mass 1/p^4 but
reach only grows by a factor p (the congruence modulus grows by p^6, the bound
shrinks by p).  Mass shrinks faster than log(reach) grows, so the EV of such a
sliver is ~ ln(p)/p^4 -- negligible.  The lattice bounds below make the extreme
case rigorous: all four bases divisible by 42 forces f >= 1,698,000,953.

## 5. Why there is no algebraic shortcut

Any rational curve on x1^6+...+x5^6 = x6^6 yields rational points, so finding
one is at least as hard as finding a point: specialising the curve at a
parameter value already gives a solution (for a conic, the u^6 coefficient
equation is itself Sum a_i^6 = e^6).  Linear parametrisations are impossible
outright (at a real zero of the sixth-power side every term vanishes).  The
repeated-coordinate surface 2X^6+2Y^6+Z^6 = W^6 is of general type with no
genus-one fibration (repo, GEOMETRIC_STRIKE.md).  Sub-loci not previously
considered: 3X^6+Y^6+Z^6 = W^6 (general type, locally solvable) and the curve
W^6-E^6 = 5A^6 (genus 10, finitely many points).  3X^6+2Y^6 = W^6 has no
primitive solution (mod 7).

## 6. Rigorous lower bounds (strata/lattice_bounds.py, verified three ways)

| extra prime dividing all four small bases | f >= |
|---|---|
| 2 | 235,283 |
| 3 | 395,243 |
| 7 | 11,936,849 |
| 42 | 1,698,000,953 |
