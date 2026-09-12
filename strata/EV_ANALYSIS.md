# Expected value of a counterexample search — CORRECTED 2026-09-12

**This file replaces an earlier version whose numbers were wrong. The errors and
their sizes are documented in section 5, because the earlier numbers were quoted
in `SESSION_2026_09_11.md` and in PR #3.**

## 1. Method (the part that was previously wrong)

For `x1^k + ... + xn^k = y^k` the Hardy--Littlewood heuristic gives, for the
number of solutions with `y` in a range,

    N  ~  (1/n!) * S * A * log(range ratio),        A = Gamma(1+1/k)^n / Gamma(n/k),

where `S = prod_p sigma_p` is the singular series and the `1/n!` converts
ordered tuples to distinct (unordered) solutions.

**The trap.** The naive local density

    sigma_p(e) = #{(x,y) mod p^e : sum xi^k = y^k} / p^(e*n)

**does not converge.** Every solution has all its scalar multiples in the
solution set, and that family grows with `e`. Measured directly:

| k=6, p=7 | p^1 | p^2 | p^3 | p^4 |
|---|---|---|---|---|
| naive | 0.01077 | 0.01113 | 0.01363 | 0.03112 |
| **primitive** | **0.01071** | **0.01071** | **0.01071** | **0.01071** |

The same divergence appears at `p=2` and `p=3` for k=6 (naive 2-adic runs
0.75, 0.75, 0.875, 1.125, 1.625, ...). The quantity that converges is the
density of **primitive** solutions — those not having every coordinate
divisible by `p`:

    sigma_p = lim_e [ N_total(p^e) - N_all-divisible(p^e) ] / p^(e*n).

All numbers below use that. Convergence was verified per prime by computing
successive `e` until the value repeated to 12 decimals.

## 2. Converged local densities

k=6, n=5 (stable exponents 2^3, 3^2, 7^1, rest p^1):

| p | 2 | 3 | **7** | 5 | 11 | 13 | 19 | 31 | 37 | 43 |
|---|---|---|---|---|---|---|---|---|---|---|
| sigma_p | 0.6250 | 0.2469 | **0.01071** | 1.0317 | 1.0075 | 2.8302 | 0.8534 | 0.7392 | 1.6408 | 0.8242 |

k=7, n=6 (stable exponents 7^2, rest p^1):

| p | 2 | 3 | **7** | 5 | 29 | 43 | 71 | 113 |
|---|---|---|---|---|---|---|---|---|
| sigma_p | 0.9844 | 0.9986 | **2.1191** | 0.9999 | 0.6752 | 1.6058 | 0.9610 | 0.9448 |

The single factor `sigma_7 = 0.0107` for k=6 versus `2.119` for k=7 is the whole
story: at k=6 exactly one of the five left terms may be coprime to 7, and that
costs a factor of ~100. At k=7 the prime 7 *helps*.

## 3. Converged totals (stable over primes <= 200, 500, 1000)

| | S (primes<=1000) | A | **distinct solutions per e-fold of f** |
|---|---|---|---|
| (6,1,5) | 0.00715 | 0.60879 | **3.63e-05** |
| (7,1,6) | 2.06837 | 0.60593 | **1.74e-03** |

**k=7 advantage: 48x** (successive prime bounds give 54, 49.1, 48.0 — converged).
The count is log-uniform in `f`: no height is a sweet spot.

## 4. What the searches were actually worth

| search | expected distinct solutions |
|---|---|
| (6,1,5), all classes, f <= 4.3e6 (the pre-existing frontier) | 5.5e-04 |
| (6,1,5), all classes, 4.3e6 < f <= 2e7 | 5.6e-05 |
| **this session's k=6 sliver** (class 1 ~1/25, k7<=2 ~10%, to 2e7) | **2.2e-07** |
| **this session's k=7 run**, 1500 < f <= 6500 | **2.6e-03** |
| **(7,1,6) cumulative, 2 < f <= 6500** | **1.5e-02** |

So the k=6 sliver was worth about **one chance in 4.5 million**, and the k=7
work about **one in 65** cumulatively (one in 390 for the increment added this
session). The pivot to k=7 was correct — more strongly than the original,
badly-derived argument claimed — but no run here had a good chance of success.

## 5. Errors in the previous version, and their sizes

1. **Divergent method.** The old local densities were the naive ones. They only
   looked stable because they were evaluated at low exponents (8, 9, 7, 13, ...)
   where naive and primitive nearly agree. At higher exponents they diverge.
   Corrected S: 0.0070 -> 0.00715 (k=6), 2.19 -> 2.068 (k=7). By luck the old
   totals were close; the *method* was unsound and could not be relied on.
2. **Mismatched prime sets.** The old "70x" compared a k=6 product over
   {8,9,7,13,19,31,37,43} with a k=7 product over {8,9,49,29,43,71,113,127,197,211}.
   Corrected, converged ratio: **48x**, not 70x.
3. **Ordered vs unordered mixed.** The old EV table used the **ordered** k=6
   density (0.0042 per e-fold) without dividing by 5! = 120, while the k=7
   figure had been divided by 6!. Every k=6 EV in that table, and the
   "one chance in ten thousand" summary, is therefore **overstated by ~120x**.
   The true figure for the k=6 plans is about **one in a million**.

## 6. Unaffected by all of this

The search results themselves. Coverage, leaf counts and zero-solution verdicts
are exact machine facts; the four lattice lower bounds in `lattice_bounds.py`
are proved arithmetic. Only the *motivating estimates* were wrong.
