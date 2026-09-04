# Astra attack matrix

Scores are 1–5. For leverage, asymptotic improvement, distance, informative
failure, discovery probability, and certificate quality, higher is better.
For CPU and RAM, 5 means cheaper. Probability scores are comparative
judgments, not calibrated probabilities.

| Rank | Attack | Math leverage | Asymptotic | CPU | RAM | Distance | Failure value | CE chance | Certificate | Total |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | Degree ≤4 rational curves / elliptic pencils on `2X^6+2Y^6+Z^6=W^6` | 5 | 5 | 3 | 4 | 5 | 5 | 3 | 5 | 35 |
| 2 | Structured offline `(S+S)∩T` algorithm exploiting the root-polynomial target family | 5 | 5 | 2 | 3 | 5 | 5 | 3 | 4 | 32 |
| 3 | Independently certify and port all five Meyrignac classes with fused GPU peel/probe | 4 | 2 | 2 | 2 | 4 | 5 | 4 | 5 | 28 |
| 4 | p-adic valuation strata and genuine descent search | 5 | 4 | 4 | 5 | 4 | 4 | 2 | 4 | 32 |
| 5 | Close class 2–4 coverage from 3.0M toward 4.3M | 3 | 1 | 2 | 2 | 4 | 3 | 3 | 5 | 23 |
| 6 | Optimize tracked mod 19/31/37 cascade; remove weak filters by measurement | 2 | 1 | 5 | 5 | 2 | 3 | 2 | 5 | 25 |
| 7 | Scale external routed join / line-sorted mmap query batches | 3 | 2 | 2 | 4 | 3 | 4 | 2 | 4 | 24 |
| 8 | Concentrated 4.3M→5.5M caseA3 calibration | 1 | 1 | 2 | 4 | 2 | 2 | 2 | 5 | 19 |

Rank 4 ties rank 2 numerically but remains fourth because no descent invariant
is currently identified; rank 2 has a crisp complexity target. Rank 3 is
ahead of its raw score because it repairs the most important evidence gap:
spread coverage is externally reported without retained production logs in
the audited clone.

## Top three operating lanes

1. **Geometry:** solve or certify empty the bounded repeated-coordinate
   surface ansatz.
2. **Asymptotics:** either exploit the structured target family to beat the
   offline `|S||T|` join, or record a rigorous barrier for each attempted
   decomposition.
3. **Class-complete computation:** durable logs, checksums, planted controls,
   and independent CPU slices for all five classes; then close the weaker
   class frontiers before pushing class 1.

“Search farther” is deliberately not a top-three attack.
