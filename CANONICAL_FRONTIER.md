# Canonical frontier and provenance

## Later direct strike: separate sparse domain

The direct strike based on `485390a` completed 10,382 explicitly selected
repeated-coordinate surface targets and retained a standard-library-only
independent replay. Its selected heights reach about `9.66e11`, but this is
**not** a continuous exclusion through that value. Four bounded conic
parameter searches also completed without a square lift. Neither result
extends the historical concentrated frontier or the all-class frontier.
See [ASTRA_DIRECT_STRIKE.md](ASTRA_DIRECT_STRIKE.md) and
[RATIONAL_POINT_LEDGER.md](RATIONAL_POINT_LEDGER.md).

## Evidence grades

- **A:** independently reproduced exact result with a retained command/output.
- **B:** committed contemporary result plus independent cheap consistency
  checks, but no retained original machine log/checksum.
- **C:** external campaign report not reproduced here.

The historical branch explicitly says its `/tmp` logs did not survive.  The
current `docs/RESULTS.md` is a later curated transcript.  Therefore no row is
silently promoted to grade A merely because a `done:` line was copied into
Markdown.

## Historical concentrated campaign

Every production row used Git blob
`0a7d9e6299bdea43edb7c85cf65e86ab66048d62` (`euler6/caseA2.c`).
That version used the power-of-two Bloom allocation; the over-allocation
changes memory and false-positive rate, not completeness.

| Range | Candidates | Solutions | Completion evidence | Original log | Audit status |
|---|---:|---:|---|---|---|
| 700,000–730,000 | 124 | 0 | pre-handoff control; later curated transcript | unavailable | **A**, replayed 2026-09-04 |
| 730,000–1,000,000 | 1,314 | 0 | commit `cdee5da`, no later than 2026-07-31 19:52:53Z | unavailable | **B** |
| 1,000,000–1,500,000 | 3,523 | 0 | commit `cdee5da`, no later than 2026-07-31 19:52:53Z | unavailable | **B** |
| 1,500,000–2,000,000 | 5,129 | 0 | commit `cdee5da`, no later than 2026-07-31 19:52:53Z | unavailable | **B** |
| 2,000,000–2,500,000 | 7,222 | 0 | commit `cdee5da`, no later than 2026-07-31 19:52:53Z | unavailable | **B** |
| 2,500,000–3,200,000 | 12,443 | 0 | commit `a4815d4`, 2026-07-31 20:15:09Z | unavailable | **B** |
| 3,200,000–4,000,000 | 18,003 | 0 | commit `b9a2dae`, 2026-08-01 00:44:09Z | unavailable | **B** |
| 4,000,000–4,300,000 | 8,050 | 0 | commit `4d565a4`, 2026-08-01 08:40:40Z | unavailable | **B** |
| **730,000–4,300,000** | **55,684** | **0** | chronological commits above | unavailable | **B** |

The arithmetic is exact:

```text
1,314 + 3,523 + 5,129 + 7,222 + 12,443 + 18,003 + 8,050 = 55,684.
```

## Independent candidate-set replay

`make frontier-audit` regenerated all candidate sets on an Intel Xeon Platinum
8573C using nine threads.  Wall time: 0.553 s.

| Range | Count | XOR digest | Sum digest |
|---|---:|---|---|
| 700,000–730,000 | 124 | `f3d2af1c8ee17283` | `4e79e22fe6d29995` |
| 730,000–1,000,000 | 1,314 | `db60584a67eb16a1` | `d1eda2046563a691` |
| 1,000,000–1,500,000 | 3,523 | `15c415743bf94217` | `9f8d51cd8c14e65f` |
| 1,500,000–2,000,000 | 5,129 | `1c85a7f1beb1b106` | `71771e92c9e4fc02` |
| 2,000,000–2,500,000 | 7,222 | `076b84a2a77bd61b` | `50000433d496e245` |
| 2,500,000–3,200,000 | 12,443 | `9952a747b0b87612` | `1d60ddc1e0504258` |
| 3,200,000–4,000,000 | 18,003 | `a6bfc7c71a32dfe9` | `6c92cdb8482c71f1` |
| 4,000,000–4,300,000 | 8,050 | `1824fcfe8195fe7b` | `b8507d887eab4ad5` |

This certifies the candidate enumeration and reconciliation.  It does not
replace the lost exact-decomposition logs.

## Broader all-class frontier

The independent repository
[`cavedave/six-one-five`](https://github.com/cavedave/six-one-five) reports:

| Regime | Reported frontier | Grade here |
|---|---:|---|
| all five Meyrignac classes | 2,353,973 | C |
| classes 2–4 | 3,000,000 | C |
| class 5 | 5,000,000 | C |
| class 1 | 3,680,000 | C |

Its code and commit chronology are public and its README describes planted
tests and a CPU/GPU cross-check, but production logs/checksums are not tracked
in that clone and no suitable GPU is available here.  These figures therefore
inform attack selection without replacing ESOP6's canonical grade-B class-1
frontier of 4.3M.

## Classical 730k floor

Resta–Meyrignac's 2003 paper on sixth-power equal sums reports the historical
search and is the source cited for the 730,000 bound.  It remains the
independent published floor for all cases.  The exact domain must be stated as
the positive primitive `(6,1,5)` equation; results for `(6,2,5)` are not a
substitute.
