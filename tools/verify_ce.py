#!/usr/bin/env python3
"""Independent exact verifier for an ESOP6 sextuple."""

from __future__ import annotations

import argparse
import math


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Check a^6+b^6+c^6+d^6+e^6 == f^6 using Python integers"
    )
    parser.add_argument("values", metavar="N", type=int, nargs=6)
    args = parser.parse_args()
    a, b, c, d, e, f = args.values

    if min(args.values) <= 0:
        print("INVALID: all six integers must be positive")
        return 2
    lhs = sum(x**6 for x in (a, b, c, d, e))
    rhs = f**6
    if lhs != rhs:
        print(f"FAIL: lhs={lhs} rhs={rhs} delta={lhs-rhs}")
        return 1

    primitive = math.gcd(math.gcd(math.gcd(math.gcd(math.gcd(a, b), c), d), e), f) == 1
    print(f"PASS: exact equality; primitive={'yes' if primitive else 'no'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
