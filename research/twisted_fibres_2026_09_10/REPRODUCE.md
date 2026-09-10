# Reproduction and resume

Do not regenerate the retained ledger to resume arithmetic. The enumeration
program refuses to overwrite an existing ledger unless explicitly requested.
No historical integer sweep or curve lottery is launched by these commands.

## Exact checks using standard-library Python

From repository root:

```bash
python research/twisted_fibres_2026_09_10/restore_evidence.py
python research/twisted_fibres_2026_09_10/fibres.py --selftest
python tests/direct_verifiers.py
python research/twisted_fibres_2026_09_10/audit_verify.py
python research/twisted_fibres_2026_09_10/verify_seed_lists.py
python research/quotient_descent_2026_09_10/independent_finite_verify.py
```

The independent arithmetic replay checks all representations, prime factors,
model scalings, exact point equations, inverse maps, rejection entries, and
the full torsion groups using finite-field group orders and independent
rational group arithmetic. It checks the retained PARI rank transcripts and
their use; it does not constitute an independent implementation of 2-descent.

`arithmetic_evidence.zip` contains the full ledger, every model transcript,
the six-gate checkpoint and execution logs. `evidence_manifest.json` gives
SHA256 hashes. The restore step verifies every byte and refuses to overwrite
changed local data. `FIBRE_LEDGER.md` and `survivors.json` are the compact
reviewable views; they do not replace the full evidence.

## PARI/GP dependency

Use official PARI/GP **2.17.4**, pinned because the class-group certification
checks inspect that version's rank-initialization structure. Its source is
at [the official release URL](https://pari.math.u-bordeaux.fr/pub/pari/unix/pari-2.17.4.tar.gz).
SHA256:

```
02651d99c391007d384b3fadbc20abc6916b77036f9e496c99e9ce8688ca4b53
```

The recorded unprivileged native-kernel build used `./Configure` with
`--without-gmp --without-readline --graphic=none --static`, then `make -j4 gp`.
The bundled ellrank, elltors, ellratpoints regressions passed. Build details
and executable provenance are in `runtime_manifest.json`.

The initial run used 3 independent workers, 8 wall seconds per model and
effort 0, through the first six quotient gates. The final three gates used
12 seconds per model. Raw stage logs and the six-gate checkpoint are retained.
Every process receives the fixed random seed 20260910. Inconclusive models
are cached as such; a future retry must be explicit:

```bash
python research/twisted_fibres_2026_09_10/run_arithmetic.py --gp /absolute/path/to/gp --seconds 30 --workers 3 --retry-failed
```

This resumes only unresolved fibres, reuses successful cached models, and
explicitly retries failed/timed-out models. It does not infer a complete
Mordell–Weil basis from a matching pair of rank bounds. The optional bounded
point hunt has its separate result file and exact declared coefficient boxes.

## New enumeration ranges

The retained pilot covers every sorted primitive integer fourtuple with
entries in 1..10, plus the two retained Y seeds. The resulting normalized
coefficients describe whole rational fibres, not a six-variable height box.
For a genuinely new pilot, choose a new output directory and carry across
cached model files as appropriate. Never count a failed rational point box
as a globally excluded fibre.

The verified N500 seed replacement and four >300 chunks are provided for
future curve jobs. They were not queued here. Original seed lists and logs
are retained unchanged; 31 invalid N500 entries are documented exactly.
