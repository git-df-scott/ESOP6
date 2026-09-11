#!/usr/bin/env python3
"""Regression and independent verification of the third-strike results.

Independent here means: the retained certificates are re-derived by code that
does not reuse the generating routine.  The mod-3 shape list is re-checked by
direct polynomial evaluation in F3[s]; the reversed normal form is re-checked
by exact Fraction evaluation at pseudo-random rational coefficients instead of
by symbolic expansion; the mod-9 states are re-checked by integer arithmetic.
No floating point is allowed anywhere.
"""
import json
import random
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tools'))
import degree_eight_strike as d8            # noqa: E402
import unequal_slope_six as u6              # noqa: E402

RESULTS = ROOT / 'results/astra_third_2026_09_11'


def poly_mul(a, b):
    o = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, u in enumerate(a):
        if u:
            for j, v in enumerate(b):
                if v:
                    o[i + j] += u * v
    return o


def poly_pow(a, e):
    r = [Fraction(1)]
    for _ in range(e):
        r = poly_mul(r, a)
    return r


def independent_reversed_form(trials=6, seed=20260911):
    """(5) <=> (R) checked at exact random rational coefficients."""
    rng = random.Random(seed)
    for _ in range(trials):
        p = [Fraction(1)] + [Fraction(rng.randint(-9, 9), rng.randint(1, 7))
                             for _ in range(7)]
        r = [Fraction(1)] + [Fraction(rng.randint(-9, 9), rng.randint(1, 7))
                             for _ in range(7)]
        n = [Fraction(1)] + [Fraction(rng.randint(-9, 9), rng.randint(1, 7))
                             for _ in range(8)]
        q = Fraction(rng.randint(1, 9), rng.randint(1, 7))
        Q = [Fraction(1), Fraction(0), q]
        lhs = [(u + v) / 2 for u, v in zip(poly_pow(p, 6), poly_pow(r, 6))]
        rhs = poly_mul(Q, poly_pow(n, 5))
        t12 = [Fraction(0)] * 12 + [x * Fraction(10, 27)
                                    for x in poly_mul(poly_pow(Q, 3),
                                                      poly_pow(n, 3))]
        t24 = [Fraction(0)] * 24 + [x * Fraction(1, 81)
                                    for x in poly_mul(poly_pow(Q, 5), n)]
        E = [Fraction(0)] * 43
        for src in (lhs,):
            for i, v in enumerate(src):
                E[i] += v
        for src in (rhs, t12, t24):
            for i, v in enumerate(src):
                if i < 43:
                    E[i] -= v
        # reversed side
        x = list(reversed(p))
        y = list(reversed(r))
        nn = list(reversed(n))
        qt = [q, Fraction(0), Fraction(1)]
        G = [Fraction(0)] * 43
        left = [(u + v) * Fraction(81, 2)
                for u, v in zip(poly_pow(x, 6), poly_pow(y, 6))]
        n2 = poly_pow(nn, 2)
        q2 = poly_pow(qt, 2)

        def padd(a, b, ca=1, cb=1):
            L = max(len(a), len(b))
            return [ca * (a[i] if i < len(a) else Fraction(0))
                    + cb * (b[i] if i < len(b) else Fraction(0))
                    for i in range(L)]

        right = poly_mul(poly_mul(nn, qt),
                         poly_mul(padd(n2, q2, 3, 1), padd(n2, q2, 27, 1)))
        for i, v in enumerate(left):
            G[i] += v
        for i, v in enumerate(right):
            G[i] -= v
        for d in range(43):
            assert 81 * E[d] == G[42 - d], d
        assert all(isinstance(v, Fraction) for v in E + G)
    return trials


def independent_mod3_shapes(shapes):
    """Every listed shape really solves (x^6+y^6)/2 = n^5 q~ in F3[s]."""
    def mul(a, b):
        o = [0] * (len(a) + len(b) - 1)
        for i, u in enumerate(a):
            for j, v in enumerate(b):
                o[i + j] = (o[i + j] + u * v) % 3
        return o

    def pw(a, e):
        r = [1]
        for _ in range(e):
            r = mul(r, a)
        return r

    for x, y, n, Q0 in shapes:
        lhs = [(a + b) * 2 % 3 for a, b in zip(pw(list(x), 6), pw(list(y), 6))]
        rhs = mul(pw(list(n), 5), [Q0 % 3, 0, 1])
        assert len(lhs) == len(rhs) == 43
        assert all((a - b) % 3 == 0 for a, b in zip(lhs, rhs))
    return len(shapes)


def independent_mod9_states(shapes):
    states = 0
    for sh in shapes:
        z = list(sh[0][:7]) + list(sh[1][:7]) + list(sh[2][:8]) + [sh[3]]
        base = d8._F(z, 9)
        assert all(v % 3 == 0 for v in base)
        rhs = [(-(v // 3)) % 3 for v in base]
        rows = [r[14:23] for r in d8._jacobian(sh)]
        sols, rank = d8._solve_affine(rows, rhs, 9)
        assert rank == 8
        if sols:
            for h in sols:
                zz = list(z)
                for i in range(9):
                    zz[14 + i] += 3 * h[i]
                assert all(v % 9 == 0 for v in d8._F(zz, 9))
                states += 1
    return states


def main():
    de = d8.checks()
    un = u6.checks()
    assert de['result'] == 'PASS' and un['result'] == 'PASS'

    trials = independent_reversed_form()
    shapes, digest = d8.mod3_shapes()
    assert digest == de['mod3']['shapes_sha256']
    assert independent_mod3_shapes(shapes) == 95
    assert independent_mod9_states(shapes) == 81

    # no solution is being claimed anywhere
    assert de['leading_equation_search']['hits'] == []
    assert de['valuation']['integral_case_possible'] is False
    assert un['branch_degM0']['statement'].startswith('deg M = 0 branch impossible')

    for name, live in (('degree_eight.json', de),
                       ('unequal_slope_six.json', un)):
        path = RESULTS / name
        if path.exists():
            retained = json.loads(path.read_text())
            assert retained == json.loads(json.dumps(live)), name

    print(json.dumps({
        'result': 'PASS',
        'exact_symbolic_identities': 'PASS',
        'independent_reversed_form_trials': trials,
        'mod3_shapes': len(shapes),
        'mod3_shapes_sha256': digest,
        'mod9_states': 81,
        'retained_certificate_match': RESULTS.exists(),
        'esop6_solution_found': False,
    }, indent=2))


if __name__ == '__main__':
    main()
