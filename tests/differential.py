#!/usr/bin/env python3
"""Randomized candidate-set differential tests for caseA2 and caseA3."""

from __future__ import annotations

import os
import random
import re
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = re.compile(r"^CANDIDATE (\d+) (\d+)$", re.MULTILINE)


def run(program: str, lo: int, hi: int, buckets: int | None) -> set[tuple[int, int]]:
    command = [str(ROOT / "bin" / program), str(lo), str(hi), "12"]
    if buckets is not None:
        command += ["-b", str(buckets)]
    command.append("--dump-candidates")
    env = dict(os.environ, OMP_NUM_THREADS="4")
    proc = subprocess.run(
        command,
        cwd=ROOT,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
        timeout=180,
    )
    return {(int(f), int(t)) for f, t in CANDIDATE.findall(proc.stdout)}


def main() -> int:
    rng = random.Random(0xE50F6)
    ranges = []
    for _ in range(3):
        lo = rng.randrange(700_000, 721_000)
        ranges.append((lo, lo + rng.randrange(6_000, 9_001)))

    for lo, hi in ranges:
        expected = run("caseA2", lo, hi, None)
        for buckets in (1, 2, 3, 4, 7, 8):
            actual = run("caseA3", lo, hi, buckets)
            if actual != expected:
                missing = sorted(expected - actual)[:5]
                extra = sorted(actual - expected)[:5]
                raise SystemExit(
                    f"FAIL range=[{lo},{hi}] buckets={buckets} "
                    f"missing={missing} extra={extra}"
                )
        print(f"PASS range=[{lo},{hi}] candidates={len(expected)} buckets=1,2,3,4,7,8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
