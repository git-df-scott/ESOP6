# caseA3 completeness proof

## Definitions

Let \(S(x,y)=x^6+y^6\), stored for every ordered-by-size pair
\(1\le y\le x\le B_{max}\).  For a 128-bit sum `s`, caseA3 computes

```text
h(s) = mix(lo(s) xor C*(hi(s)+1))
bucket(s) = high64(h(s) * NB)
```

where multiplication is unsigned 128-bit and `high64` is Lemire range
reduction into `[0,NB)`.  `mix` is a fixed deterministic 64-bit avalanche
function.  Bloom line selection uses an independent fixed hash and the same
range-reduction construction into `[0,nlines)`.

## Exactly-one-bucket theorem

For every j=2 query remainder \(R\), `bucket(R)` is a total deterministic
function with one value \(k\in\{0,\ldots,NB-1\}\).

During pass \(k\):

1. every pair sum \(S(x,y)\) satisfying `bucket(S(x,y)) == k` is inserted;
2. every query satisfying `bucket(R) == k` is evaluated;
3. queries assigned to other buckets return at the bucket gate.

Therefore every query is evaluated in exactly one pass.  If
\(R=S(x,y)\), equality of the full 128-bit values implies identical inputs to
the hash, hence identical buckets.  The generating pair is present in the
filter in the unique pass that evaluates the query.

## Collisions and Bloom behavior

- **Bucket collision:** unrelated sums may share a bucket.  This increases
  load only; it cannot remove the true sum.
- **Bloom hash collision:** unrelated sums may set the same bits.  This causes
  a false positive only.
- **Bloom false negatives:** impossible under the algorithmic model because
  insertion sets all eight queried bits and no operation clears the filter
  during a pass.  The filter is cleared only between completed passes.
- **Final decision:** a Bloom positive calls `pair_verify`, which checks
  \(R=x^6+y^6\) with exact unsigned-128 arithmetic and exact class budgets.
  Thus Bloom false positives cannot become reported solutions.

The per-bucket allocation uses the expected `P/NB` entry count.  Hash
imbalance can raise a bucket's false-positive rate but cannot overflow the
fixed bit array or create a false negative.

## Interaction with DFS

All j=4 and j=3 traversal is repeated in every pass.  A wrong-bucket j=2 leaf
returns failure for that pass, but the same deterministic path is revisited in
its assigned pass.  If a valid four-part decomposition exists, its final pair
remainder is evaluated in that pass and exact verification succeeds.

Early return after a found solution affects only instrumentation after the
solution; it does not suppress the solution.  On a zero-solution equivalence
control, complete node totals must match exactly.

## Width assumptions

- Pair sums and remainders are unsigned 128-bit.
- `f <= 1e8` keeps the reduced target below \(2^{128}\).
- Bucket and Bloom range reduction multiply two 64-bit values into unsigned
  128-bit and take the high half.
- `NB >= 1`; command parsing rejects a missing `-b` argument and normalizes
  zero to one.
- OpenMP atomic OR protects concurrent Bloom insertion; the implicit loop
  barrier precedes queries.

## Certificates

`make equiv` reproduces:

```text
NB=1: reached 4,139,088; skipped 0; evaluated 4,139,088
NB=8: reached 33,112,704; skipped 28,973,616; evaluated 4,139,088
```

and both give 124 candidates, 124 processed, 0 found.  Additionally,
`make differential` compares actual `(f,t)` sets on three fixed-random bands
against caseA2 for NB 1, 2, 3, 4, 7, and 8; all 18 comparisons pass.

## Verdict

**Correct for the concentrated search, subject to ordinary hardware/compiler
correctness.**  Bucketing changes time and memory, not the searched set.
