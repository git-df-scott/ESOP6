#!/usr/bin/env python3
"""Exact finite certificate for the normalized degree-six contact obstruction.

Standard library only. This is coefficient lifting, not an integer-point
search. The reduction from arbitrary rational coefficients to BASE is proved
in BOUNDARY_CONTACT_6.md. A formal series is never reported as a rational
function or as an ESOP6 solution.
"""
import argparse
import hashlib
import itertools
import json
from fractions import Fraction
from pathlib import Path
import time

VARIABLES = ['u2', 'u3', 'u4', 'u5', 'v1', 'v2', 'v3', 'v4', 'v5']
BASE = [2, 0, 1, 0, 1, 0, 2, 0, 1]


def mul(a, b, modulus=None, limit=None):
    size = len(a) + len(b) - 1
    if limit is not None:
        size = min(size, limit)
    out = [0] * size
    for i, x in enumerate(a):
        if not x or i >= size:
            continue
        for j, y in enumerate(b[:size-i]):
            if y:
                out[i+j] += x*y
    return out if modulus is None else [x % modulus for x in out]


def power(a, exponent, modulus=None, limit=None):
    out = [1]
    for _ in range(exponent):
        out = mul(out, a, modulus, limit)
    if limit is not None:
        out += [0] * (limit-len(out))
    return out


def midpoint_jet(d, degree, modulus=None):
    """Unique M=1+... with M^5=D through the requested degree."""
    m = [1]
    for k in range(1, degree+1):
        known = power(m, 5, modulus, k+1)[k]
        value = (d[k] if k < len(d) else 0) - known
        m.append(Fraction(value, 5) if modulus is None
                 else value * pow(5, -1, modulus) % modulus)
    return m


def residual(q, modulus):
    """Coefficients 7..30 of (A^6+B^6)/2-M_6^5 modulo modulus.

    The 3-adic rescaling makes the omitted midpoint correction terms
    divisible by 3^9. We require only a necessary identity modulo 3^6.
    """
    u, v = [1, 0, *q[:4]], [0, *q[4:]]
    a, b = [x+y for x, y in zip(u, v)], [x-y for x, y in zip(u, v)]
    aa, bb = power(a, 6, modulus), power(b, 6, modulus)
    d = [(x+y)*pow(2, -1, modulus) % modulus for x, y in zip(aa, bb)]
    m = midpoint_jet(d, 6, modulus)
    return [(x-y) % modulus for x, y in zip(d, power(m, 5, modulus))][7:]


class LiftSpace:
    """Solve Lh=b over F_3, retaining every free digit."""
    def __init__(self):
        r = residual(BASE, 27)
        assert all(x % 9 == 0 for x in r)
        cols = []
        for j in range(9):
            q = BASE.copy()
            q[j] += 3
            delta = [(x-y) % 27 for x, y in zip(residual(q, 27), r)]
            assert all(x % 9 == 0 for x in delta)
            cols.append([x//9 for x in delta])
        self.matrix = [[cols[j][i] for j in range(9)] for i in range(24)]
        rows = [row + [int(i == j) for j in range(24)]
                for i, row in enumerate(self.matrix)]
        self.pivots = []
        rank = 0
        for j in range(9):
            p = next((k for k in range(rank, 24) if rows[k][j]), None)
            if p is None:
                continue
            rows[rank], rows[p] = rows[p], rows[rank]
            inv = pow(rows[rank][j], -1, 3)
            rows[rank] = [x*inv % 3 for x in rows[rank]]
            for k in range(24):
                if k != rank:
                    c = rows[k][j]
                    rows[k] = [(x-c*y) % 3 for x, y in zip(rows[k], rows[rank])]
            self.pivots.append(j)
            rank += 1
        self.rank, self.rows = rank, rows
        self.free = [j for j in range(9) if j not in self.pivots]
        self.kernel = []
        for j in self.free:
            vec = [int(i == j) for i in range(9)]
            for row, p in zip(rows, self.pivots):
                vec[p] = -row[j] % 3
            self.kernel.append(vec)

    def particular(self, rhs):
        b = [sum(x*y for x, y in zip(row[9:], rhs)) % 3 for row in self.rows]
        if any(b[self.rank:]):
            return None
        out = [0]*9
        for i, j in enumerate(self.pivots):
            out[j] = b[i]
        return out

    def digits(self, particular):
        for free in itertools.product(range(3), repeat=len(self.free)):
            yield [(x + sum(c*k[j] for c, k in zip(free, self.kernel))) % 3
                   for j, x in enumerate(particular)]


def fnv1a(data):
    h = 14695981039346656037
    for byte in data:
        h = ((h ^ byte)*1099511628211) & ((1 << 64)-1)
    return f'{h:016x}'


def certificate():
    start = time.monotonic()
    space = LiftSpace()
    assert space.rank == 4 and space.pivots == [0, 1, 2, 3]
    prefixes, levels = [BASE], []
    for place in [3, 9, 27, 81]:
        modulus = 9*place
        prefixes.sort()
        data = ''.join(','.join(map(str, q))+'\n' for q in prefixes).encode()
        children, liftable = [], 0
        for q in prefixes:
            r = residual(q, modulus)
            assert all(x % (3*place) == 0 for x in r)
            rhs = [(-x//(3*place)) % 3 for x in r]
            particular = space.particular(rhs)
            if particular is None:
                continue
            liftable += 1
            for h in space.digits(particular):
                children.append([x+place*y for x, y in zip(q, h)])
        levels.append({'modulus': modulus, 'coefficient_digit_place': place,
                       'prefixes': len(prefixes), 'liftable': liftable,
                       'children': len(children),
                       'prefix_sha256': hashlib.sha256(data).hexdigest(),
                       'prefix_fnv1a64': fnv1a(data)})
        prefixes = children
    assert [(x['prefixes'], x['liftable'], x['children']) for x in levels] == [
        (1, 1, 243), (243, 81, 19683), (19683, 162, 39366), (39366, 0, 0)]
    return {'result': 'NO NORMALIZED LIFT MODULO 729', 'variables': VARIABLES,
            'base_residue_mod_3': BASE, 'jacobian_divided_by_3_mod_3': space.matrix,
            'rank': space.rank, 'pivots': space.pivots,
            'kernel_basis': space.kernel, 'levels': levels,
            'prefix_checks_including_initial': sum(x['prefixes'] for x in levels),
            'unperturbed_identity_valid_modulus_at_least': 3**9,
            'obstruction_modulus': 3**6,
            'scope': 'normalized degree-six polynomial boundary-contact ansatz over Q',
            'seconds': round(time.monotonic()-start, 6)}


def formal_midpoint(a, b, degree):
    """Exact Q[[t]] jet for the actual, perturbed midpoint equation."""
    d = [(x+y)*Fraction(1, 2) for x, y in zip(power(list(map(Fraction, a)), 6, limit=degree+1),
                                power(list(map(Fraction, b)), 6, limit=degree+1))]
    m = [Fraction(1)]
    for k in range(1, degree+1):
        known = power(m, 5, limit=k+1)[k]
        if k >= 12:
            known += Fraction(10, 27)*power(m, 3, limit=k-11)[k-12]
        if k >= 24 and k-24 < len(m):
            known += Fraction(1, 81)*m[k-24]
        m.append((d[k]-known)/5)
    return m


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = certificate()
    encoded = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded)
    print(encoded, end='')


if __name__ == '__main__':
    main()
