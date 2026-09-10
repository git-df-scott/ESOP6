#!/usr/bin/env python3
"""Search for points of  A^6+B^6+C^6 = F^6  (and A^6+B^6+C^6+D^6=F^6)
over real quadratic fields Q(sqrt D), D>0 non-square.

This is a NECESSARY condition for a twisted rational curve with |T|=2 (resp. |T|=1):
the two (resp. one) anti-invariant coordinates vanish at the fixed points of sigma,
leaving the remaining coordinates as a point of the surface over Q(sqrt Delta).

Method: (u+v sqrt D)^6 = A + B sqrt D with A,B in Z.  Meet in the middle on
A^6+B^6 = F^6 - C^6 as integer 2-vectors.
"""
import sys
from math import isqrt, gcd


def sixth(u, v, D):
    a, b = 1, 0
    for _ in range(6):
        a, b = a*u + b*v*D, a*v + b*u
    return a, b


def gen(H, D):
    """Representatives of (u+v sqrt D) up to sign, |u|,|v| <= H, not both 0."""
    out = []
    for v in range(0, H+1):
        for u in range(-H, H+1):
            if v == 0 and u <= 0:
                continue
            if gcd(u, v) != 1:      # projective: keep primitive only
                continue
            out.append((u, v))
    return out


def search3(D, H, verbose=True):
    reps = gen(H, D)
    pw = {r: sixth(r[0], r[1], D) for r in reps}
    two = {}
    for i, r1 in enumerate(reps):
        a1, b1 = pw[r1]
        for r2 in reps[i:]:
            a2, b2 = pw[r2]
            two.setdefault((a1+a2, b1+b2), []).append((r1, r2))
    hits = []
    for rf in reps:
        af, bf = pw[rf]
        for rc in reps:
            ac, bc = pw[rc]
            key = (af-ac, bf-bc)
            if key in two:
                for (r1, r2) in two[key]:
                    hits.append((r1, r2, rc, rf))
    return hits, len(reps), len(two)


def classify(h, D):
    """Drop degenerate hits: rational points, zero coordinates, duplicates of LHS=RHS."""
    r1, r2, rc, rf = h
    vals = [r1, r2, rc, rf]
    if all(v[1] == 0 for v in vals):
        return 'rational'
    # a coordinate equal (as a sixth power) to the RHS makes the rest sum to 0
    s6 = [sixth(*v, D) for v in vals]
    if s6[3] in s6[:3]:
        return 'degenerate(one term = RHS)'
    return 'NONTRIVIAL'


if __name__ == '__main__':
    H = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    Ds = [d for d in range(2, 40) if isqrt(d)**2 != d]
    print(f'(6,3,1) over Q(sqrt D):  A^6+B^6+C^6 = F^6,  coefficient bound H={H}')
    print(f'{"D":>4} {"#reps":>7} {"#2sums":>9}  result')
    total = 0
    for D in Ds:
        hits, nr, n2 = search3(D, H)
        nt = [h for h in hits if classify(h, D) == 'NONTRIVIAL']
        total += len(nt)
        note = 'none' if not nt else f'*** {len(nt)} NONTRIVIAL ***'
        print(f'{D:>4} {nr:>7} {n2:>9}  raw={len(hits):<5} {note}')
        for h in nt[:5]:
            print('      ', h)
    print('TOTAL nontrivial:', total)
