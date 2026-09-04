# ESOP6 — Euler's Sum of Powers, sixth-power case

A search for a counterexample to Euler's sum of powers conjecture at k = 6:

```
a⁶ + b⁶ + c⁶ + d⁶ + e⁶ = f⁶      in positive integers
```

Euler conjectured in 1769 that no such solution exists (more generally, that
a k-th power needs at least k k-th powers to be written as a sum). The
conjecture held for ~200 years, then fell twice:

| k | year | who | counterexample |
|---|---|---|---|
| 5 | 1966 | Lander & Parkin, CDC 6600 brute force | 27⁵ + 84⁵ + 110⁵ + 133⁵ = 144⁵ |
| 4 | 1988 | Elkies (elliptic fibrations), minimal case by Frye | 95800⁴ + 217519⁴ + 414560⁴ = 422481⁴ |
| 6 | — | **open** | none known |

For k = 6 the question is wide open. Published exhaustive searches
(Meyrignac's EulerNet ecosystem and successors, ca. 2002) established that
**no solution exists with f ≤ 730,000**. Distributed searches have hammered
the problem since the 1990s, so raw compute is not a winning strategy — the
cost of an exhaustive search grows like ~f^3.3, while the expected number of
solutions grows only like log f.

**Status of this project: no counterexample found.** What it did produce is a
new, mathematically-structured exclusion region and a validated toolchain,
described below.

---

## Headline result

> **The concentrated-exemption case of a⁶+b⁶+c⁶+d⁶+e⁶ = f⁶ has no solution
> for 730,000 < f ≤ 4,300,000.**
>
> 55,684 candidate pairs survived the modular sieve in that range; every one
> was tested to exhaustion and rejected. This is a slice of the search space
> extending 5.9× beyond the published exhaustive frontier, and no prior search
> is known to have covered it.

Euler's conjecture survives, but it had to defend ground nobody had attacked.

---

## The mathematics

### 1. Structure theorem (elementary, and it carries the whole project)

Sixth powers are extremely rigid modulo small numbers, because 6 = φ(7) =
φ(9) and 6 ≥ 3:

| modulus | possible values of x⁶ |
|---|---|
| 7 | 0 (if 7\|x), else 1 |
| 9 | 0 (if 3\|x), else 1 |
| 8 | 0 (if x even), else 1 |

Reduce the equation mod 7. The left side contributes 1 for each term coprime
to 7 and 0 otherwise, so it equals (number of terms coprime to 7) mod 7. The
right side is 0 or 1. With five terms, the count must be 0 or 1 exactly. But a
count of 0 means 7 divides all six numbers — the solution is not primitive.
The same argument runs mod 9 and mod 8. Therefore:

> **In any primitive solution, exactly one of a…e is odd, exactly one is
> coprime to 3, exactly one is coprime to 7 — and f is coprime to 42.**

Every solution is an integer multiple of a primitive one with smaller f, so
restricting to primitives loses nothing.

This immediately says something striking: **four of the five terms on the left
are divisible by 2, four by 3, four by 7.** The solution set is far thinner
than a naive count suggests.

### 2. The concentrated case, and why it can be swept deep

The three "exemptions" (odd / coprime-to-3 / coprime-to-7) may land on one,
two, or three distinct terms. The **concentrated case** is when all three land
on the *same* term t. Then the other four terms are each divisible by 2, 3 and
7 — that is, by 42 — and reducing the equation modulo 42⁶ kills them entirely:

```
f⁶ ≡ t⁶  (mod 42⁶),        42⁶ = 5,489,031,744
```

with both f and t units mod 42⁶. So t/f is a **sixth root of unity modulo
42⁶**. By CRT over 2⁶ · 3⁶ · 7⁶ there are exactly

```
4 × 6 × 6 = 144
```

such roots. For each f, only 144 residues of t are admissible out of ~5.5
billion, and only those with t < f matter. The expected number of candidate
pairs (f, t) with f ≤ F is

```
      2      F²
N ≈  ─── · 144 · ─────  ,
      7        2 · 42⁶
```

which stays in the tens of thousands even at F = 4.3 million. **This is the
leverage**: a sub-case that generic enumeration cannot reach becomes sweepable
millions deep, because the modular structure collapses the haystack rather
than the needle.

Each surviving candidate still has to be finished off: writing

```
m = (f⁶ − t⁶) / 42⁶
```

the remaining four terms are 42·(bases), and m must be a sum of four sixth
powers with bases ≤ (f−1)/42. In that reduced equation the same rigidity
returns as *exact* counting constraints:

- number of odd bases = m mod 8
- number of bases coprime to 3 = m mod 9
- number of bases coprime to 7 = m mod 7

each of which must be ≤ 4 or the candidate is instantly infeasible. When a
budget hits 0, every remaining base is forced divisible by that prime (the
search strides by up to 42); when a budget equals the number of remaining
bases, all of them are forced coprime.

### 3. Honest accounting of coverage

The complementary cases — exemptions spread over two or three terms — leave
only three or two terms divisible by 42. Their candidate density scales like
F³/42⁶ or worse, so the modular collapse does **not** thin them below generic
cost; they remain covered only to the classical f ≤ 730,000 bound. Under a
naive equidistribution model the concentrated case carries roughly 1/25 of the
heuristic solution mass.

So this project cleared one sub-case very deep rather than all cases a little
deeper. That was a deliberate trade: it is the only slice where modest
hardware reaches genuinely unsearched territory.

---

## Results

All runs on 4 cores / 15 GB RAM.

### Exhaustive search (all cases)

| f range | solutions |
|---|---|
| 2 – 5,000 | 0 |

(Consistency check against the literature; the exhaustive searcher's ~f^3.3
scaling makes larger ranges impractical on this hardware, which is what
motivated the structured approach.)

### Concentrated case

| f range | candidate pairs | solutions |
|---|---|---|
| 700,000 – 730,000 | 124 | 0 — *control run inside the known-clear region* |
| 730,000 – 1,000,000 | 1,314 | **0** |
| 1,000,000 – 1,500,000 | 3,523 | **0** |
| 1,500,000 – 2,000,000 | 5,129 | **0** |
| 2,000,000 – 2,500,000 | 7,222 | **0** |
| 2,500,000 – 3,200,000 | 12,443 | **0** |
| 3,200,000 – 4,000,000 | 18,003 | **0** |
| 4,000,000 – 4,300,000 | 8,050 | **0** |
| **730,000 – 4,300,000 total** | **55,684** | **0** |

The ceiling at 4.3M was the memory wall of the machine that ran the sweep:
the pair-sum filter grows as (f/42)². That wall has since been removed by the
bucketed variant (`src/caseA3.c`), so the frontier is now limited by time
rather than RAM. See [HANDOFF.md](HANDOFF.md) for how to push higher.

---

## Validation — why these zeros are trustworthy

A null result is only worth as much as the evidence that the machinery would
have reported a hit. Four independent checks:

1. **It finds a real counterexample.** Run in fifth-power mode, the same DFS
   engine rediscovers Lander & Parkin's 1966 result — 27⁵+84⁵+110⁵+133⁵ =
   144⁵ — in about 30 ms (`make validate`).
2. **Candidate enumeration is provably complete.** `src/audit.c` counts, by
   brute force over all t < f, the solutions of t⁶ ≡ f⁶ (mod 42⁶), and
   compares against the 144-root enumeration. They agree exactly on every
   sampled f, including f = 54,321,011. Nothing is being skipped.
3. **The control run agrees with the literature.** The band 700,000–730,000
   lies inside the published exhaustively-searched region; the sweep finds
   124 candidates and zero solutions there, as it must.
4. **Planted solutions are recovered; infeasible targets are rejected.** The
   decomposition routine has built-in self-tests that reconstruct known
   four-sixth-power sums (both mixed-class and all-divisible-by-42) and
   correctly reject a residue-infeasible target. Bloom-filter positives are
   always verified exactly, so a false positive costs time, never correctness.

---

## Code

```
src/search.c    exhaustive (6,1,5) DFS — class-constraint strides, incremental
                residue sieves mod 13/43, exact per-level analytic bounds.
                Doubles as the k=5 validation harness.

src/caseA2.c    the workhorse: concentrated-case sweep. 144-root candidate
                enumeration, exact class budgets with forced strides up to 42,
                residue masks mod 64/27/49/13/43 maintained incrementally, and
                a blocked Bloom filter over all pair sums x⁶+y⁶ (8 hash bits
                confined to one 64-byte cache line) with exact verification.
                128-bit arithmetic throughout.

src/caseA3.c    low-memory variant of the sweep. Partitions the pair table by
                a hash of the sum into NB buckets and makes NB passes, so the
                filter holds 1/NB of the table; every query is still answered
                in exactly one pass. Removes the RAM wall at the cost of
                repeating pair enumeration and DFS traversal per pass.
                Approach contributed by Duncan.

src/audit.c     completeness audit of the candidate enumeration.

src/caseA.c     earlier, slower iteration of the concentrated-case search.
                Kept for provenance; superseded by caseA2.c.
```

### Build and run

```sh
make                          # builds everything
make validate                 # k=5 sanity check: must print the 144^5 solution
make control                  # 700k-730k control run: 124 candidates, 0 found
make equiv                    # bucketed sweep evaluates the same nodes as monolithic

./bin/caseA2 4300000 5000000 12       # continue the sweep (~10.6 GB)
./bin/caseA3 4300000 5000000 12 -b 8  # same search, ~1.3 GB, 8 passes
./bin/search 6 2 20000             # exhaustive search, all cases
./bin/audit                        # completeness audit
```

Requires gcc with OpenMP. Hard limit: f ≤ 10⁸ (128-bit overflow guard).
Memory for `caseA2` is (f/42)²/2 × bits-per-pair / 8 bytes; `caseA3` divides
that by the bucket count:

| f | caseA2 @12 bpp | caseA3, 8 buckets | caseA3, 16 buckets |
|---|---|---|---|
| 3.2M | 4.4 GB | 0.54 GB | 0.27 GB |
| 5.0M | 10.6 GB | 1.33 GB | 0.66 GB |
| 10.0M | 42.5 GB | 5.31 GB | 2.66 GB |

Time scales roughly linearly in the bucket count, so use the smallest NB that
fits comfortably in free RAM.

**If a `SOLUTION` line ever appears: do not announce it.** Verify it first with
independent exact arithmetic (Python big integers are ideal — the whole claim
is checkable in one line). The protocol is in [HANDOFF.md](HANDOFF.md) §6.

---

## Where an actual breakthrough would come from

Not from more CPU. The expected number of solutions below F grows like
C·log F while the cost of finding them grows polynomially, and C is evidently
small for k = 6 — nothing below 730,000, whereas the k = 5 counterexample sits
at f = 144. Each e-fold of additional height costs roughly F⁴ work. That
treadmill is what stopped every brute-force campaign before this one.

The promising direction is algebraic, and the geometry is encouraging. The
variety

```
X:  x⁶ + y⁶ + z⁶ + u⁶ + v⁶ = w⁶     in P⁵
```

is a smooth sextic fourfold with trivial canonical bundle — a **Calabi–Yau
fourfold**. That is the same structural class as the quartic surface Elkies
cracked in 1988 (a K3, the Calabi–Yau surface case): not of general type, so
not obstructed by the Bombieri–Lang philosophy that predicts scarcity of
rational points. Elkies' method was to find elliptic fibrations with sections
of infinite order and search along them.

No analogous fibration argument over ℚ is known for the sextic fourfold.
Finding one — or any ℚ-rational curve on X with all coordinates nonzero —
would beat any amount of enumeration. Note that Shioda's inductive structure
(Fermat varieties dominated by products of Fermat curves) only operates over
cyclotomic extensions, and the Fermat sextic *curve* has no nontrivial
rational points by Fermat's Last Theorem at n = 6, so nothing descends that
way. This is the open door.

---

## Continuing the work

[HANDOFF.md](HANDOFF.md) is written so that a new session or a new machine can
resume with zero context loss: proven results, code map, hardware limits,
next actions, verification protocol.

The natural next steps, in ascending order of expected value:

1. **More height in the concentrated case.** Anything above f = 4,300,000 is
   unclaimed. With `caseA3` the memory wall is gone — f = 10M fits in under
   3 GB at 16 buckets — so the binding constraint is now time, which scales
   as ~F⁴ overall.
2. **The spread-exemption cases**, still covered only to 730,000. These need a
   different algorithmic idea — the modular collapse does not apply.
3. **The algebraic geometry.** See above. This is where the problem is
   actually likely to fall.
