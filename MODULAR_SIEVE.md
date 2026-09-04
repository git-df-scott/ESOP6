# Modular sieve audit

## Why 42 is special

For sixth powers,

```text
x^6 mod 8, 9, or 7 is 0 when the base is divisible by 2, 3, or 7,
and 1 otherwise.
```

In a primitive five-term solution, `f` cannot be divisible by 2, 3, or 7;
otherwise a zero right-hand residue forces every left term to share that
prime.  Therefore exactly one left term is odd, exactly one is coprime to 3,
and exactly one is coprime to 7.  When those roles coincide at `t`, the other
four terms are divisible by `lcm(2,3,7)=42`, giving

\[
f^6\equiv t^6\pmod{42^6},\qquad 42^6=2^6 3^6 7^6.
\]

Thus 42 was not selected by empirical tuning.  It is the largest cheap
role-counting collapse supplied simultaneously by the primes whose relevant
unit groups make sixth powers 1 modulo 8, 9, and 7.

## Existing cascade

1. root enumeration modulo \(42^6\) selects candidate `(f,t)`;
2. exact class budgets use `m mod 8,9,7`;
3. DFS sumset masks use mod 64, 27, 49, 13, and 43;
4. caseA3 bucket gate;
5. Bloom membership;
6. exact two-sixth-power verification.

The power moduli 64/27/49 capture lifting information beyond the top-level
count.  Mod 13 is the strongest small unrestricted two-sum filter.  Mod 43 is
weaker at j=2 and saturates by j=3, but was already in the historical engine.

## Systematic small-prime screen

`python3 tools/modulus_rank.py --limit 127` enumerates sixth-power residues and
j-fold sumsets from first principles.

| p | sixth-power residues | 2-sum pass | 3-sum pass | 4-sum pass |
|---:|---:|---:|---:|---:|
| 13 | 3 | 0.384615 | 0.538462 | 0.692308 |
| 7 | 2 | 0.428571 | 0.571429 | 0.714286 |
| 37 | 7 | 0.513514 | 1.000000 | 1.000000 |
| 31 | 6 | 0.516129 | 0.838710 | 1.000000 |
| 19 | 4 | 0.526316 | 0.842105 | 1.000000 |
| 43 | 8 | 0.674419 | 1.000000 | 1.000000 |

Primes above 43 rapidly saturate for three and four terms.  Mod 37 is useful
only at the final pair, while 19 and 31 also prune the j=3 stage and fit in
64-bit masks.  They are the best additions per table byte.

## Implemented mod-19/mod-31 leaf cascade

`caseA3 --extra-sieve` tracks remainders incrementally through the DFS and,
after the bucket gate, tests the pair-sum masks before Bloom lookup.

700k–730k, 16 bpp, 8 threads:

| Mode | Reached | Bucket-skipped | Extra-sieve rejected | Bloom evaluated | Wall |
|---|---:|---:|---:|---:|---:|
| NB=1 baseline | 4,139,088 | 0 | 0 | 4,139,088 | 6.177 s |
| NB=1 + extra | 4,139,088 | 0 | 3,590,968 | 548,120 | 5.572 s |
| NB=8 baseline | 33,112,704 | 28,973,616 | 0 | 4,139,088 | 10.568 s |
| NB=8 + extra | 33,112,704 | 28,973,616 | 3,590,968 | 548,120 | 10.402 s |

The rejection is 86.76%.  Wall improvement is about 9.8% for NB=1 and 1.6%
for NB=8 on this control.  It is enabled explicitly rather than silently
changing the canonical equivalence count.

## Recommendation

Use:

```text
class/valuation checks -> existing 64/27/49/13/43 masks
-> caseA3 bucket gate -> tracked 19/31 pair masks
-> Bloom -> exact verifier
```

Do not add large residue tables.  Before production, benchmark mod 37 as a
third pair-only mask and test whether removing mod 43 improves net time.  The
decision metric is wall time and Bloom/exact calls, not rejection percentage
alone.
