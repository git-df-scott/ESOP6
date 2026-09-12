#!/usr/bin/env python3
"""Exact independent verifier for a^6+b^6+c^6+d^6+e^6 == f^6 (Python big ints)."""
import sys
from math import gcd
from functools import reduce

def main(argv):
    if len(argv) != 7:
        sys.exit("usage: verify_solution.py a b c d e f")
    a, b, c, d, e, f = (int(x) for x in argv[1:])
    lhs = a**6 + b**6 + c**6 + d**6 + e**6
    rhs = f**6
    g = reduce(gcd, (a, b, c, d, e, f))
    print("a b c d e f = %d %d %d %d %d %d" % (a, b, c, d, e, f))
    print("LHS = %d" % lhs)
    print("RHS = %d" % rhs)
    print("EXACT %s" % (lhs == rhs))
    print("gcd = %d, PRIMITIVE %s" % (g, g == 1))
    return 0 if lhs == rhs else 1

if __name__ == "__main__":
    sys.exit(main(sys.argv))
