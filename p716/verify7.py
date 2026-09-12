#!/usr/bin/env python3
"""Exact independent verifier for seventh-power (k,1,n) identities.

  verify7.py F T1 T2 ... Tn          check F^7 == sum Ti^7
  verify7.py --sum S T1 ... Tn       check S    == sum Ti^7   (planted targets)
  verify7.py --exp K ...             use exponent K instead of 7

Generalises strata/verify_solution.py (which is fixed to exponent 6, 5 terms)
to exponent 7 and an arbitrary number of terms.  Pure Python big integers.
"""
import sys
from math import gcd
from functools import reduce


def main(argv):
    argv = list(argv[1:])
    K = 7
    if "--exp" in argv:
        i = argv.index("--exp")
        K = int(argv[i + 1])
        del argv[i:i + 2]
    as_sum = False
    if argv and argv[0] == "--sum":
        as_sum = True
        del argv[0]
    if len(argv) < 2:
        sys.exit("usage: verify7.py [--exp K] [--sum] TARGET T1 T2 ... Tn")
    vals = [int(x) for x in argv]
    target, terms = vals[0], vals[1:]
    lhs = sum(t ** K for t in terms)
    rhs = target if as_sum else target ** K
    ok = (lhs == rhs)
    g = reduce(gcd, terms + ([] if as_sum else [target]))
    print("exponent = %d, terms = %d" % (K, len(terms)))
    print("target %s = %d" % ("S" if as_sum else "F^%d" % K, rhs))
    print("terms  = %s" % " ".join(str(t) for t in terms))
    print("LHS = %d" % lhs)
    print("RHS = %d" % rhs)
    print("EXACT %s" % ok)
    print("gcd = %d, PRIMITIVE %s" % (g, g == 1))
    if any(t <= 0 for t in terms):
        print("NONPOSITIVE TERM -> INVALID")
        return 1
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
