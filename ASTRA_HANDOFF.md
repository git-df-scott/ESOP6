# Astra handoff

## Binding state

- Verdict: **YELLOW**.
- Concentrated class 1: historical zero-result through 4,300,000; 55,684
  candidate sets independently regenerated, expensive decomposition logs not
  independently replayed.
- caseA3: complete; NB=8 memory ÷8, measured control penalty 1.71×.
- Best added leaf filter: tracked mod 19 + mod 31, 86.76% fewer Bloom queries.
- No proved improvement below the F^4 traversal.
- External all-class claims exist through 2,353,973, classes 2–4 through 3M,
  and class 5 through 5M; treat them as reported until logs/checksums are
  independently preserved.
- Do not run toward 10M.

## Exact first Astra command

```bash
make clean && make -j"$(nproc)" && make validate control equiv frontier-audit prototypes && make differential
```

Every target must pass. `make differential` takes about 3.6 minutes on the
audit host because it runs 18 full-engine comparisons.

## First Astra research task

Read `GEOMETRIC_ATTACK.md`, then classify degree 1–4 rational curves on

```text
2 X^6 + 2 Y^6 + Z^6 = W^6
```

with exact symbolic certificates and positive-real-component testing. Use the
Bremner–Choudhry–Ulas identity in that file as a symbolic control. Do not
advance degree after a timeout; split the normalized ansatz and certify each
component.

## First computational task, only if geometry is paused

Reconstruct all five Meyrignac class tables independently, reproduce a small
published-overlap band and planted solutions, and emit immutable JSON
manifests containing command, Git commit, compiler, hardware, range, class,
candidate/probe/hit counts, wall time, and SHA-256 of log and binary. Close
classes 2–4 before extending class 1.

For a bounded class-1 calibration after all gates:

```bash
OMP_NUM_THREADS="$(nproc)" ./bin/caseA3 4300001 4400000 12 -b 8 --extra-sieve --valuation-prune
```

This is a calibration command, not authorization for 4.4M→10M.

## Candidate protocol

For any four-part `SOLUTION` line, reconstruct the five original summands and
run:

```bash
python3 tools/verify_ce.py a b c d e f
```

Require positive inputs, exact equality, and an independently repeated check
before any announcement. `DEGENERATE` lines are not ESOP6 counterexamples.

Stop and document if any candidate count, digest, node equality, planted hit,
or CPU/GPU cross-check differs. Stop a scale-up if projected RAM exceeds the
declared ceiling, logs cannot be made durable, or timing shows the
false-positive verifier approaching the F^5 tail.
