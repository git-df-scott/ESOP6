#!/usr/bin/env python3
"""Local analysis of the NECESSARY quadratic point.

An anti-invariant form vanishes at both fixed points of sigma (see sigma_test.py, P2).
So for a twisted curve with twist set T (|T| = k of the six coordinates) the
remaining 6-k coordinates satisfy, at a fixed point theta in Q(sqrt(Delta)),
the corresponding (6, 6-k-1, 1) relation over the real quadratic field K = Q(sqrt Delta).

  |T| = 2 (source (3,3))   ->  A^6+B^6+C^6 = F^6          over K
  |T| = 1 (source (6,2,4)) ->  A^6+B^6+C^6+D^6 = F^6      over K

This script tests solvability of those relations in the residue field of K at a
rational prime p:  F_{p^2} if p is inert in K, F_p x F_p if p splits.
A non-solvable residue field kills every Delta for which p is inert.
"""
import sys
from itertools import product


def field_p2(p, D):
    """Elements of F_p[x]/(x^2-D) as pairs; returns list of sixth powers (set)."""
    sixth = set()
    for u in range(p):
        for v in range(p):
            a, b = 1, 0
            for _ in range(6):
                a, b = (a*u + b*v*D) % p, (a*v + b*u) % p
            sixth.add((a, b))
    return sixth


def inert(p, D):
    if D % p == 0:
        return None  # ramified
    return pow(D % p, (p-1)//2, p) == p-1


def solvable(sixth, p, nterms):
    """Is  s_1+...+s_n = s_0  solvable with s_i sixth powers, not all zero?"""
    S = sorted(sixth)
    # sums of nterms sixth powers
    cur = {(0, 0)}
    for _ in range(nterms):
        nxt = set()
        for (a, b) in cur:
            for (c, d) in S:
                nxt.add(((a+c) % p, (b+d) % p))
        cur = nxt
    hits = []
    for w in S:
        if w in cur and w != (0, 0):
            hits.append(w)
    # also the all-zero-target case with non-zero summands
    return hits


def main():
    primes = [7, 13, 19, 31, 37, 43, 61, 67, 73, 79]
    for nterms, label in ((3, '(6,3,1): A^6+B^6+C^6=F^6'), (4, '(6,4,1): A^6+B^6+C^6+D^6=F^6')):
        print('=' * 74)
        print(label, ' over F_{p^2} (p inert) -- nontrivial solvability')
        for p in primes:
            if (p - 1) % 6:
                continue
            sx = field_p2(p, 0)  # placeholder
            # need an actual non-residue D mod p to build F_{p^2}
            D = next(d for d in range(2, p) if pow(d, (p-1)//2, p) == p-1)
            sixth = field_p2(p, D)
            hits = solvable(sixth, p, nterms)
            print(f'  p={p:3d}  |sixth powers in F_(p^2)| = {len(sixth):4d}  '
                  f'nontrivial solutions: {"YES" if hits else "NONE"}  '
                  f'({len(hits)} admissible RHS values)')
    print('=' * 74)
    print('Rational comparison (F_p, p=7): sixth powers mod 7 = {0,1}')
    print('  a^6+b^6+c^6=f^6 mod 7 forces (#nonzero among a,b,c) in {0,1} mod 7 -> 7-adic collapse')


if __name__ == '__main__':
    main()
