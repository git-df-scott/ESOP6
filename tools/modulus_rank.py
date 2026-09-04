#!/usr/bin/env python3
"""Rank modular sixth-power sumset filters by rejection strength."""

from __future__ import annotations

import argparse
import math


def primes(limit: int) -> list[int]:
    answer = []
    for n in range(2, limit + 1):
        if all(n % p for p in range(2, math.isqrt(n) + 1)):
            answer.append(n)
    return answer


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=127)
    args = parser.parse_args()
    rows = []
    for modulus in primes(args.limit):
        residues = {pow(x, 6, modulus) for x in range(modulus)}
        sums = {0}
        rates = []
        for _ in range(4):
            sums = {(x + y) % modulus for x in sums for y in residues}
            rates.append(len(sums) / modulus)
        score = -math.log2(rates[1]) / math.log2(modulus)
        rows.append((rates[1], -score, modulus, len(residues), rates))

    print("| modulus | |R6| | 2-sum pass | 3-sum pass | 4-sum pass | 2-sum bits rejected |")
    print("|---:|---:|---:|---:|---:|---:|")
    for _, _, modulus, size, rates in sorted(rows)[:20]:
        print(
            f"| {modulus} | {size} | {rates[1]:.6f} | {rates[2]:.6f} | "
            f"{rates[3]:.6f} | {-math.log2(rates[1]):.3f} |"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
