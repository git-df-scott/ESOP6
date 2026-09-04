# HANDOFF — continuing the ESOP6 search

Everything needed to resume this project on a new machine or in a new session,
with no context loss. Read [README.md](README.md) first for the mathematics;
this file is the operational companion.

## 1. What is already established — do not re-derive

**Structure theorem.** In any primitive solution of a⁶+b⁶+c⁶+d⁶+e⁶ = f⁶,
exactly one of a…e is odd, exactly one is coprime to 3, exactly one is coprime
to 7, and f is coprime to 42. (Proof: x⁶ ≡ 0 or 1 mod each of 7, 9, 8; count
terms; a count of 0 contradicts primitivity.) Searching primitives is fully
general.

**Concentrated-case reduction.** If all three exemptions sit on one term t,
the other four terms are divisible by 42, so f⁶ ≡ t⁶ (mod 42⁶) with
42⁶ = 5,489,031,744. The admissible t are exactly u·f mod 42⁶ over the 144
sixth roots of unity mod 42⁶ (4 mod 2⁶ × 6 mod 3⁶ × 6 mod 7⁶, by CRT —
verified in code by direct powering, and audited by brute force).

**Second-stage constraints.** With m = (f⁶−t⁶)/42⁶, the four remaining bases
satisfy: (#odd) = m mod 8, (#coprime to 3) = m mod 9, (#coprime to 7) = m mod
7, each ≤ 4 or the candidate is infeasible.

**Coverage caveat.** Spread-exemption cases (2 or 3 distinct exempt terms)
have candidate density ~F³/42⁶ or worse and are **not** thinned by this
method. They remain covered only to the classical f ≤ 730,000. The
concentrated case carries ~1/25 of heuristic solution mass.

**Heuristic expectation.** Expected solutions below F grow like C·log F with C
small for k=6. The 730k→4.3M sweep carried perhaps a 1–5% chance of a hit.
That is the honest number; the null result is the expected outcome, not a
surprise.

## 2. Results table (authoritative)

| f range | case | candidates | solutions |
|---|---|---|---|
| ≤ 5,000 | all (exhaustive) | — | 0 |
| 700,000–730,000 | concentrated | 124 | 0 (control) |
| 730,000–1,000,000 | concentrated | 1,314 | 0 |
| 1,000,000–1,500,000 | concentrated | 3,523 | 0 |
| 1,500,000–2,000,000 | concentrated | 5,129 | 0 |
| 2,000,000–2,500,000 | concentrated | 7,222 | 0 |
| 2,500,000–3,200,000 | concentrated | 12,443 | 0 |
| 3,200,000–4,000,000 | concentrated | 18,003 | 0 |
| 4,000,000–4,300,000 | concentrated | 8,050 | 0 |
| **730,000–4,300,000** | **concentrated** | **55,684** | **0** |

Ranges are idempotent and independent — re-running one is always safe, and
chunks can be distributed across machines without coordination beyond claiming
the range.

## 3. Hard limits

- `fmax ≤ 1e8` — enforced in code (128-bit overflow guard).
- `caseA2` Bloom RAM = `(fmax/42)² / 2 × bits_per_pair / 8` bytes exactly.
  `caseA3 -b NB` divides that by NB:

  | f | caseA2 @12 bpp | caseA3 ×8 | caseA3 ×16 |
  |---|---|---|---|
  | 3.2M | 4.4 GB | 0.54 GB | 0.27 GB |
  | 5.0M | 10.6 GB | 1.33 GB | 0.66 GB |
  | 10.0M | 42.5 GB | 5.31 GB | 2.66 GB |

  Runtime grows roughly linearly in NB (pair enumeration and DFS traversal
  repeat per pass; exact verification does not). Use the smallest NB that fits
  comfortably in *free* RAM, not total RAM.
- Do not go below ~12 bits/pair: false-positive rate rises and exact
  verification cost explodes.
- **Fixed 2026-09:** `caseA2` used to round the filter's line count up to a
  power of two, so a 12 bpp request at fmax=5e6 allocated 17.2 GB rather than
  10.6 GB (an effective 19.4 bpp). Sizing is now exact, via Lemire range
  reduction instead of power-of-two masking. Bug found by Duncan.
- Runtime scales ~F⁴ overall (candidates ~F², cost each ~F²). Budget
  accordingly: the 4.0–4.3M chunk took ~6 h on 4 cores.

## 4. Immediate next actions

1. `make && make validate && make control` — three checks, under two minutes
   total. `validate` must print the 144⁵ solution; `control` must report 124
   candidates and 0 found. If either fails, stop and debug before trusting any
   sweep.
2. Claim a range above **f = 4,300,000** (everything below is done) and run
   either `./bin/caseA2 <fmin> <fmax> [bpp]` (one pass, needs the full filter
   in RAM) or `./bin/caseA3 <fmin> <fmax> [bpp] -b NB` (NB passes, 1/NB the
   memory). `make equiv` proves the two evaluate identical oracle nodes.
3. Record the `done:` line in [docs/RESULTS.md](docs/RESULTS.md), update the
   tables in README.md and this file, and commit. **The repo is the only
   durable memory** — a cloud container's /tmp does not survive.
4. Run long sweeps detached (`nohup`, `screen`, or a background task) and
   check on them periodically; do not block a session on them.

## 5. Suggested upgrades for a long campaign

Neither is needed for short chunks; both are worth a few hours of work before
a multi-week run:

- **Checkpointing.** Periodically write the highest completed f so an
  interrupted chunk resumes instead of restarting. (A chunk dying silently
  mid-run has already happened once — container hiccup — and cost a full
  re-run.)
- ~~Disk-backed sorted pair-sum table~~ — **done differently, and better**:
  `caseA3` partitions the table in memory by a hash of the sum, which removes
  the RAM ceiling with no disk I/O at all. f = 10M now fits in 2.7 GB at 16
  buckets. Contributed by Duncan.

A GPU port is *not* recommended: this is integer, branch-heavy,
memory-latency-bound work, not float throughput.

## 6. Protocol if a SOLUTION line ever prints

1. **Do not announce it.** Not to anyone, not yet.
2. Verify with independent exact arithmetic, in a different language than the
   searcher. Python:
   ```python
   a,b,c,d,e,f = ...          # parts printed as parts42 are already ×42
   assert a**6+b**6+c**6+d**6+e**6 == f**6
   ```
3. Sanity-check the shape: every term < f, all positive, and compute the gcd
   (a primitive solution should have gcd 1; a non-primitive one is still valid
   but implies a smaller solution exists — divide it out and re-verify).
4. Commit the tuple to this repo immediately. A commit is a timestamped,
   public priority claim and costs nothing.
5. Then tell the world. It would be the first known (6,1,5) solution and would
   settle a question open since 1769 — the arithmetic is one line, so anyone
   can confirm it in seconds, which is exactly why the verification must come
   first.

## 7. Context and precedent

- k=5 fell to brute force in 1966 (Lander & Parkin); k=4 to algebraic
  construction in 1988 (Elkies, via elliptic fibrations on a K3 surface, with
  Frye finding the minimal case). Construction has beaten enumeration whenever
  both were available.
- July 2026: a counterexample to the Jacobian conjecture (open since 1939) was
  announced by Levent Alpoge with LLM assistance — for n ≥ 3; the n = 2 case
  remains open, and it has not yet been journal peer-reviewed. Again a
  *construction* problem rather than a search problem. Worth keeping in view
  when deciding where to spend effort on this one.
