# Concentrated and spread search regimes

## Mathematical partition

In a primitive solution there is one odd term, one term not divisible by 3,
and one term not divisible by 7.  Distribute these three labelled roles among
the five summands.  Up to permutation, the operational Meyrignac partition is:

| Class | Shape | Guaranteed scaled core | Description |
|---:|---|---|---|
| 1 | all three roles coincide | four terms `42*c` | concentrated case; ESOP6 caseA2/A3 |
| 2 | a two-role unit plus one free role | three `42*c`, one `14*d` | spread, free 3-role |
| 3 | a two-role unit plus one free role | three `42*c`, one `21*d` | spread, free parity role |
| 4 | a two-role unit plus one free role | three `42*c`, one `7*d` | spread, combined 2/3 behavior |
| 5 | three roles on three terms | two `42*c`, frees `21*d` and `14*e` | fully spread |

The exact residue seed tables differ by class; a complete search must run all
five for every eligible `f`.  ESOP6's 4.3M headline covers class 1 only.

## Verified and reported frontiers

| Domain | Frontier | Evidence here |
|---|---:|---|
| all classes, published historical | 730,000 | published floor |
| class 1, ESOP6 | 4,300,000 | candidate sets replayed; zero-result historical grade B |
| all classes, `cavedave/six-one-five` | 2,353,973 | external report, grade C |
| classes 2–4, same external campaign | 3,000,000 | external report, grade C |
| class 5, same external campaign | 5,000,000 | external report, grade C |

The external campaign changes planning: “all spread territory stops at 730k”
is no longer the strongest public computational claim.  Because its raw
production logs are not committed and its GPU run was not reproduced here,
the numbers must retain the `reported` qualifier.

## Information per CPU-hour

| Rank | Lane | Value assessment |
|---:|---|---|
| 1 | structural/geometric construction | only lane capable of bypassing the F^4 frontier; failure can prove a bounded ansatz empty |
| 2 | independently certify external all-class 730k–2.354M record | converts strategically important grade-C spread coverage into durable certificates |
| 3 | port fused class-5-style peel/probe and class-complete residue tables | large throughput gain where the leaf becomes a two-sum; preserves completeness |
| 4 | close classes 2–4 from 3.0M toward 4.3M | balances the frontier; more informative than extending already-deep class 1 |
| 5 | class 1 4.3M→5.5M | easy with caseA3 and <2 GB at NB=8, but low marginal log-height and no new structure |

Option A (class 1 4.3→5.5M) is therefore below weaker spread territory and the
structural lane.  This session deliberately did not run it.
