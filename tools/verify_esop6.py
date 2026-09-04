#!/usr/bin/env python3
"""Standalone ESOP6 certificate printer. Python standard library only."""
import json
import math
import sys


def main():
    if len(sys.argv) != 7:
        raise SystemExit('Usage: verify_esop6.py a b c d e f')
    values = [int(t,10) for t in sys.argv[1:]]
    powers = [n**6 for n in values]
    g = math.gcd(*values)
    normalized = sorted(values[:5]) + [values[5]]
    if g:
        normalized = [n//g for n in normalized]
    lhs, rhs = sum(powers[:5]),powers[5]
    np = [n**6 for n in normalized]
    result = {'integers':[str(n) for n in values],
              'positivity':all(n>0 for n in values),
              'sixth_powers':[str(n) for n in powers],
              'lhs':str(lhs),'rhs':str(rhs),'equality':lhs==rhs,
              'gcd':str(g),'sorted_normalized_tuple':[str(n) for n in normalized],
              'normalized_sixth_powers':[str(n) for n in np],
              'normalized_lhs':str(sum(np[:5])), 'normalized_rhs':str(np[5]),
              'normalized_equality':sum(np[:5])==np[5],
              'solution':all(n>0 for n in values) and lhs==rhs}
    print(json.dumps(result,indent=2))
    return 0 if result['solution'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
