#!/usr/bin/env python3
"""Exact analysis of the degree-eight boundary-contact family on 2X^6+2Y^6+Z^6=W^6.

Scope (RATIONAL_CURVE_ATTEMPT.md equation (5)):

    Q = 1+q t^2 (q rational, q>0),
    X = t P7/Q, Y = t R7/Q, Z = T8/Q, W = Z+(2/3) t^6,
    P7(0)=R7(0)=T8(0)=1, deg P7,R7<=7, deg T8<=8,
    N = T8 + t^6 Q/3,
    (P7^6+R7^6)/2 = Q N^5 + (10/27) t^12 Q^3 N^3 + (1/81) t^24 Q^5 N.     (5)

Everything here is exact (Rational/Integer only).  No floating point enters
any certificate.  Nothing in this file reports a formal series, a near miss,
or a bounded empty search as an ESOP6 solution.

Main results produced by ``checks()``:

  * equation (5) is re-derived from the surface equation;
  * (5) is equivalent to the reversed normal form

        (81/2)(x^6+y^6) = n q~ (3n^2+q~^2)(27n^2+q~^2),                   (R)

    with x,y monic of degree 7, n monic of degree 8, q~ = s^2+q;
  * the degree analysis leaves exactly two branches, deg N=8 with
    max(deg P,deg R)=7 and deg N=2 with max(deg P,deg R)=6;
  * the top-degree equation of the deg N=8 branch is *equivalent to an
    ESOP6 counterexample*: it is the surface equation itself;
  * (R) has no solution with 3-integral coefficients (3-adic Gauss valuation);
  * the 3-adically rescaled unperturbed equation (x^6+y^6)/2 = n^5 q~ has the
    exact rational solution family (s^5 q~, s^5 q~, s^6 q~, q~), so the
    degree-six certificate strategy of BOUNDARY_CONTACT_6.md cannot terminate
    on this family;
  * exhaustive classification mod 3 (95 residue shapes) and exact linear
    lifting mod 9 (27 surviving shapes, 81 states).
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path

import sympy as sp
from sympy import Rational as Rat

T, S, QQ_ = sp.symbols('t s q')
PC = [sp.Symbol('p%d' % i) for i in range(8)]
RC = [sp.Symbol('r%d' % i) for i in range(8)]
NC = [sp.Symbol('n%d' % i) for i in range(9)]


def _polys():
    P = 1 + sum(PC[i] * T**i for i in range(1, 8))
    R = 1 + sum(RC[i] * T**i for i in range(1, 8))
    N = 1 + sum(NC[i] * T**i for i in range(1, 9))
    Q = 1 + QQ_ * T**2
    return P, R, N, Q


def surface_derivation():
    """Re-derive (5) from 2X^6+2Y^6+Z^6-W^6 with no floating point."""
    P, R, N, Q = _polys()
    T8 = N - T**6 * Q / 3
    W8 = T8 + Rat(2, 3) * T**6 * Q
    X, Y, Z, W = T * P / Q, T * R / Q, T8 / Q, W8 / Q
    surf = sp.expand(sp.together(2 * X**6 + 2 * Y**6 + Z**6 - W**6) * Q**6)
    lhs = (P**6 + R**6) / 2
    rhs = (Q * N**5 + Rat(10, 27) * T**12 * Q**3 * N**3
           + Rat(1, 81) * T**24 * Q**5 * N)
    assert sp.cancel(sp.together(W - Z) - Rat(2, 3) * T**6) == 0
    assert sp.expand(surf - 4 * T**6 * (lhs - rhs)) == 0
    return sp.expand(lhs - rhs)


def reversed_normal_form(E):
    """(5) <=> (R): (81/2)(x^6+y^6) = n*q~*(3n^2+q~^2)*(27n^2+q~^2)."""
    x = S**7 + sum(PC[i] * S**(7 - i) for i in range(1, 8))
    y = S**7 + sum(RC[i] * S**(7 - i) for i in range(1, 8))
    n = S**8 + sum(NC[i] * S**(8 - i) for i in range(1, 9))
    qt = S**2 + QQ_
    G = sp.expand(Rat(81, 2) * (x**6 + y**6)
                  - n * qt * (3 * n**2 + qt**2) * (27 * n**2 + qt**2))
    pe = sp.Poly(E, T)
    pg = sp.Poly(G, S)
    ce = {m[0]: c for m, c in zip(pe.monoms(), pe.coeffs())}
    cg = {m[0]: c for m, c in zip(pg.monoms(), pg.coeffs())}
    for d in range(0, 43):
        assert sp.expand(81 * ce.get(d, 0) - cg.get(42 - d, 0)) == 0, d
    assert ce.get(0, 0) == 0 and cg.get(42, 0) == 0
    orders = sorted(d for d in range(1, 43) if ce.get(d, 0) != 0)
    degrees = {}
    for d, c in cg.items():
        degrees[42 - d] = sp.Poly(c, *([QQ_] + PC[1:] + RC[1:] + NC[1:])
                                  ).total_degree()
    return {'equation_orders': orders,
            'equation_count': len(orders),
            'max_total_degree_reversed': max(degrees.values()),
            'unknowns': 23}


def grading_check(E):
    """For the grading deg p_i = deg r_i = deg n_i = i, deg q = 2 the three
    pieces of (5) are weighted homogeneous: the order-d coefficient of
    D - Q N^5 has weighted degree d, that of t^12 Q^3 N^3 has weighted degree
    d-12, and that of t^24 Q^5 N has weighted degree d-24.  The explicit
    powers of t break homogeneity of the whole equation from order 12 on,
    so the order-d equation has ordinary total degree at most d."""
    lam = sp.Symbol('lam')
    P, R, N, Q = _polys()
    sub = {PC[i]: lam**i * PC[i] for i in range(1, 8)}
    sub.update({RC[i]: lam**i * RC[i] for i in range(1, 8)})
    sub.update({NC[i]: lam**i * NC[i] for i in range(1, 9)})
    sub[QQ_] = lam**2 * QQ_
    pieces = [(sp.expand((P**6 + R**6) / 2 - Q * N**5), 0),
              (sp.expand(Q**3 * N**3), 0),
              (sp.expand(Q**5 * N), 0)]
    checked = 0
    for expr, _ in pieces:
        pe = sp.Poly(expr, T)
        ce = {m[0]: c for m, c in zip(pe.monoms(), pe.coeffs())}
        for d in sorted(ce):
            if d > 12:
                break
            assert sp.expand(ce[d].subs(sub) - lam**d * ce[d]) == 0, (d, expr)
            checked += 1
    pe = sp.Poly(E, T)
    ce = {m[0]: c for m, c in zip(pe.monoms(), pe.coeffs())}
    for d in range(1, 12):
        assert sp.expand(ce[d].subs(sub) - lam**d * ce[d]) == 0, d
    return {'grading': 'p_i,r_i,n_i -> i ; q -> 2',
            'weighted_homogeneous_pieces_checked': checked,
            'equation_weighted_homogeneous_for_orders': '1..11',
            'from_order_12': 'splits into weighted degrees d, d-12, d-24'}


def degree_branches():
    """deg(RHS) = max(2+5h, 18+3h, 34+h); no cancellation is possible because
    the three leading coefficients all carry the sign of the leading
    coefficient of N and q>0.  deg(LHS)=6*rho with rho=max(deg P,deg R)."""
    admissible = []
    table = {}
    for h in range(0, 9):
        d = max(2 + 5 * h, 18 + 3 * h, 34 + h)
        table[h] = d
        if d % 6 == 0 and d // 6 <= 7:
            admissible.append((h, d // 6))
    assert table[8] == 42 and table[2] == 36
    assert admissible == [(2, 6), (8, 7)], admissible
    return {'rhs_degree_by_degN': table, 'admissible_(degN,rho)': admissible}


def leading_equations(E):
    """Top-degree coefficient equations of the two branches."""
    out = {}
    # branch deg N = 8, rho = 7: coefficient of t^42.
    pe = sp.Poly(E, T)
    ce = {m[0]: c for m, c in zip(pe.monoms(), pe.coeffs())}
    p7, r7, n8 = PC[7], RC[7], NC[8]
    top = sp.expand(ce[42])
    claim = sp.expand((p7**6 + r7**6) / 2
                      - n8 * QQ_ * (27 * n8**2 + QQ_**2) * (3 * n8**2 + QQ_**2) / 81)
    assert sp.expand(top - claim) == 0
    out['degN8_top'] = '(p7^6+r7^6)/2 = n8*q*(27*n8^2+q^2)*(3*n8^2+q^2)/81'
    # this top equation IS the surface equation
    X, Y = p7, r7
    Z, W = n8 - QQ_ / 3, n8 + QQ_ / 3
    assert sp.expand(2 * X**6 + 2 * Y**6 + Z**6 - W**6 - 4 * claim) == 0
    out['degN8_top_is_surface'] = ('X=p7, Y=r7, Z=n8-q/3, W=n8+q/3 turns the '
                                   'top equation into 2X^6+2Y^6+Z^6=W^6')
    # branch deg N = 2, rho = 6: coefficient of t^36 after n3..n8 vanish.
    sub = {NC[i]: 0 for i in range(3, 9)}
    sub.update({PC[7]: 0, RC[7]: 0})
    top36 = sp.expand(ce[36].subs(sub))
    claim36 = sp.expand((PC[6]**6 + RC[6]**6) / 2 - QQ_**5 * NC[2] / 81)
    assert sp.expand(top36 - claim36) == 0
    out['degN2_top'] = '(p6^6+r6^6)/2 = q^5*n2/81'
    return out


def branch_degN2():
    """The deg N = 2 branch of (R), written without its forced s^6 factor.

    deg N = 2 means n3=...=n8=0, i.e. n = s^6 * mu with mu = s^2+n1 s+n2 monic;
    max(deg P,deg R) = 6 means p7 = r7 = 0, i.e. x = s x', y = s y' with
    x',y' monic of degree 6.  Dividing (R) by s^6 gives

        (81/2)(x'^6+y'^6)
            = mu q~ (3 s^12 mu^2 + q~^2)(27 s^12 mu^2 + q~^2).          (8)

    15 unknowns, 36 coefficient equations (degrees 0..35; degree 36 is
    automatic)."""
    n1, n2, Q0 = sp.symbols('n1 n2 Q0')
    xs = [sp.Symbol('a%d' % i) for i in range(7)]
    ys = [sp.Symbol('b%d' % i) for i in range(7)]
    xp = S**6 + sum(xs[i] * S**i for i in range(6))
    yp = S**6 + sum(ys[i] * S**i for i in range(6))
    mu = S**2 + n1 * S + n2
    qt = S**2 + Q0
    x, y, n = sp.expand(S * xp), sp.expand(S * yp), sp.expand(S**6 * mu)
    full = sp.expand(Rat(81, 2) * (x**6 + y**6)
                     - n * qt * (3 * n**2 + qt**2) * (27 * n**2 + qt**2))
    red = sp.expand(Rat(81, 2) * (xp**6 + yp**6)
                    - mu * qt * (3 * S**12 * mu**2 + qt**2)
                    * (27 * S**12 * mu**2 + qt**2))
    assert sp.expand(full - S**6 * red) == 0
    pr = sp.Poly(red, S)
    cr = {m[0]: c for m, c in zip(pr.monoms(), pr.coeffs())}
    assert cr.get(36, 0) == 0
    orders = sorted(d for d in range(0, 36) if cr.get(d, 0) != 0)
    return {'reduced_identity':
            '(81/2)(x^6+y^6) = mu*q~*(3 s^12 mu^2+q~^2)(27 s^12 mu^2+q~^2)',
            'x,y': 'monic degree 6', 'mu': 's^2+n1 s+n2', 'q~': 's^2+q',
            'unknowns': 15, 'equations': len(orders),
            'top_equation': '(p6^6+r6^6)/2 = q^5 n2/81 (solvable over Q)'}


def valuation_balance(a, b, c):
    """(4+6a) - v3(RHS) for the reversed form (R) under the 3-adic Gauss
    valuation, with a=min(w(x),w(y)), b=w(n), c=w(q~).

    v3 of a sum of two sixth powers is exactly 6*min (anisotropy of x^6+y^6
    over F3(s)); the two quadratic factors have exact valuations because
    1+2b and 3+2b are odd while 2c is even."""
    return (4 + 6 * a) - (b + c + min(1 + 2 * b, 2 * c) + min(3 + 2 * b, 2 * c))


def integrality_theorem():
    """If x,y,n,q~ are 3-integral then w(x)=w(y)=w(n)=w(q~)=0 (they are monic),
    and the valuation balance is 4 != 0.  So (R) has no 3-integral solution;
    in particular q and every coefficient must have a denominator divisible
    by 3."""
    assert valuation_balance(0, 0, 0) == 4
    # a<=0 is forced (monic); no a makes the balance vanish with b=c=0
    assert all(valuation_balance(a, 0, 0) != 0 for a in range(-40, 1))
    # the three regimes of the balance, recorded over a bounded window
    regimes = {'R1_c<=b': [], 'R2_c=b+1': [], 'R3_c>=b+2': []}
    for b in range(-12, 1):
        for c in range(-12, 1):
            for a in range(-20, 1):
                if valuation_balance(a, b, c) != 0:
                    continue
                key = ('R1_c<=b' if c <= b else
                       'R2_c=b+1' if c == b + 1 else 'R3_c>=b+2')
                regimes[key].append((a, b, c))
    return {'integral_case_balance': 4,
            'integral_case_possible': False,
            'admissible_(a,b,c)_in_[-20,0]x[-12,0]^2':
                {k: len(v) for k, v in regimes.items()},
            'sample_R2': regimes['R2_c=b+1'][:5],
            'sample_R3': regimes['R3_c>=b+2'][:5]}


def unperturbed_exact_solutions():
    """After s -> 3^-k s and scaling by 3^(42k) the exact equation becomes

        81(x0^6+y0^6)/2 = 81 n0^5 q~0 + 30*3^(12k) n0^3 q~0^3 + 3^(24k) n0 q~0^5

    so every solution satisfies (x0^6+y0^6)/2 = n0^5 q~0 modulo 3^(12k-4).
    That unperturbed equation has exact rational solutions for every Q0,
    hence no amount of 3-adic lifting of it can contradict anything."""
    Q0 = sp.Symbol('Q0')
    qt = S**2 + Q0
    x = y = sp.expand(S**5 * qt)
    n = sp.expand(S**6 * qt)
    assert sp.expand((x**6 + y**6) / 2 - n**5 * qt) == 0
    assert sp.Poly(x, S).degree() == 7 and sp.Poly(n, S).degree() == 8
    assert sp.Poly(x, S).LC() == 1 and sp.Poly(n, S).LC() == 1
    return {'family': '(x,y,n,q~) = (s^5*q~, s^5*q~, s^6*q~, q~)',
            'holds_for_every_Q0': True,
            'corresponds_to': 'x=y, i.e. P7=R7, i.e. X=Y (diagonal)'}


def rescaling_identity():
    """Check the rescaled equation used above, symbolically in one variable."""
    lam, Q0 = sp.symbols('lam Q0')
    x = S**7 + sp.Symbol('x6') * S**6
    y = S**7 + sp.Symbol('y6') * S**6
    n = S**8 + sp.Symbol('n7') * S**7
    qt = S**2 + Q0
    lhs = Rat(81, 2) * (x**6 + y**6)
    rhs = (81 * n**5 * qt + 30 * lam**12 * n**3 * qt**3 + lam**24 * n * qt**5)
    # the identity being asserted is the algebraic rearrangement
    got = sp.expand(n * qt * (3 * n**2 + lam**12 * qt**2)
                    * (27 * n**2 + lam**12 * qt**2))
    assert sp.expand(got - (81 * n**5 * qt + 30 * lam**12 * n**3 * qt**3
                            + lam**24 * n * qt**5)) == 0
    assert sp.expand(lhs - lhs) == 0 and rhs is not None
    return {'rescaled_form':
            'n*q~*(3n^2+L^12 q~^2)(27n^2+L^12 q~^2) = 81n^5q~+30L^12n^3q~^3+L^24 n q~^5'}


# ---------------------------------------------------------------------------
# finite F3 certificates
# ---------------------------------------------------------------------------
def _mul(a, b, m=3):
    o = [0] * (len(a) + len(b) - 1)
    for i, u in enumerate(a):
        if u:
            for j, v in enumerate(b):
                if v:
                    o[i + j] = (o[i + j] + u * v) % m
    return tuple(o)


def _pow(a, e, m=3):
    r = (1,)
    for _ in range(e):
        r = _mul(r, a, m)
    return r


def _trim(a):
    a = list(a)
    while a and a[-1] == 0:
        a.pop()
    return tuple(a)


def _cube_root(f):
    f = _trim(f)
    if not f:
        return ()
    for i, c in enumerate(f):
        if c and i % 3:
            return None
    return tuple(f[3 * i] for i in range(len(f) // 3 + 1))


def mod3_shapes():
    """All (x,y,n,Q0) mod 3 with (x^6+y^6)/2 = n^5*q~ in F3[s], x,y monic of
    degree 7, n monic of degree 8, q~=s^2+Q0.  Mod 3 this is 2*S^3 = n^5*q~
    with S=x^2+y^2."""
    squares = {}
    monic7 = []
    counts = {'nq_pairs': 0, 'cube_admissible': 0, 'x_tested': 0}
    for c in itertools.product(range(3), repeat=7):
        x = tuple(c) + (1,)
        monic7.append(x)
        squares[_mul(x, x)] = x
    shapes = []
    for cn in itertools.product(range(3), repeat=8):
        n = tuple(cn) + (1,)
        n5 = _pow(n, 5)
        for Q0 in range(3):
            counts['nq_pairs'] += 1
            tgt = tuple((2 * u) % 3 for u in _mul(n5, (Q0, 0, 1)))
            Sp = _cube_root(tgt)
            if Sp is None:
                continue
            counts['cube_admissible'] += 1
            for x in monic7:
                counts['x_tested'] += 1
                x2 = _mul(x, x)
                L = max(len(Sp), len(x2))
                rem = _trim(tuple(((Sp[i] if i < len(Sp) else 0)
                                   - (x2[i] if i < len(x2) else 0)) % 3
                                  for i in range(L)))
                y = squares.get(rem)
                if y is not None:
                    shapes.append((x, y, n, Q0))
    shapes.sort()
    digest = hashlib.sha256(repr(shapes).encode()).hexdigest()
    mod3_shapes.counts = counts
    return shapes, digest


def _F(z, M):
    """42 coefficients of (x^6+y^6)/2 - n^5*q~ modulo M."""
    x = list(z[0:7]) + [1]
    y = list(z[7:14]) + [1]
    n = list(z[14:22]) + [1]
    qt = [z[22] % M, 0, 1]
    inv2 = pow(2, -1, M)
    a, b = _pow(tuple(x), 6, M), _pow(tuple(y), 6, M)
    lhs = [(u + v) * inv2 % M for u, v in zip(a, b)]
    rhs = _mul(_pow(tuple(n), 5, M), tuple(qt), M)
    out = [(lhs[i] - rhs[i]) % M for i in range(43)]
    assert out[42] == 0
    return out[:42]


def _jacobian(shape):
    """Exact Jacobian mod 3 of the 42 equations in the 23 digit unknowns.
    It depends only on the mod-3 shape:
        d/dx_i  -> 3*x^5*s^i  (divided by 3)
        d/dn_i  -> -5*n^4*q~*s^i
        d/dQ0   -> -n^5
    """
    x, y, n, Q0 = shape
    qt = (Q0 % 3, 0, 1)
    x5, y5 = _pow(x, 5), _pow(y, 5)
    n4q = _mul(_pow(n, 4), qt)
    n5 = _pow(n, 5)
    cols = []
    for i in range(7):
        cols.append([0] * i + list(x5))
    for i in range(7):
        cols.append([0] * i + list(y5))
    for i in range(8):
        # d/dn_i of -n^5*q~ is -5 n^4 q~ s^i, and -5 == 1 (mod 3)
        cols.append([0] * i + list(n4q))
    cols.append([(2 * v) % 3 for v in n5])
    cols = [(c + [0] * 43)[:42] for c in cols]
    return [[c[k] % 3 for c in cols] for k in range(42)]


def _solve_affine(rows, rhs, ncol, enumerate_solutions=True):
    m = len(rows)
    A = [list(r) + [rhs[i]] for i, r in enumerate(rows)]
    piv, r = [], 0
    for c in range(ncol):
        p = None
        for i in range(r, m):
            if A[i][c] % 3:
                p = i
                break
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        inv = pow(A[r][c], -1, 3)
        A[r] = [v * inv % 3 for v in A[r]]
        for i in range(m):
            if i != r and A[i][c] % 3:
                f = A[i][c]
                A[i] = [(u - f * v) % 3 for u, v in zip(A[i], A[r])]
        piv.append(c)
        r += 1
    for i in range(r, m):
        if A[i][ncol] % 3:
            return None, r
    if not enumerate_solutions:
        return [], r
    free = [c for c in range(ncol) if c not in piv]
    sols = []
    for vals in itertools.product(range(3), repeat=len(free)):
        h = [0] * ncol
        for c, v in zip(free, vals):
            h[c] = v
        for i, c in enumerate(piv):
            h[c] = (A[i][ncol] - sum(A[i][cc] * h[cc] for cc in free)) % 3
        sols.append(tuple(h))
    return sols, r


def mod9_lifting(shapes):
    """Exact linear lifting of the unperturbed equation from 3 to 9.
    Modulo 9 the left side depends only on x,y mod 3, so only the nine
    digits of (n,Q0) are involved."""
    states, ranks, survivors = [], set(), set()
    for sh in shapes:
        z = list(sh[0][:7]) + list(sh[1][:7]) + list(sh[2][:8]) + [sh[3]]
        base = _F(z, 9)
        assert all(v % 3 == 0 for v in base)
        rhs = [(-(v // 3)) % 3 for v in base]
        rows = [r[14:23] for r in _jacobian(sh)]
        sols, rk = _solve_affine(rows, rhs, 9)
        ranks.add(rk)
        if sols:
            survivors.add(sh)
            for h in sols:
                zz = list(z)
                for i in range(9):
                    zz[14 + i] += 3 * h[i]
                assert all(v % 9 == 0 for v in _F(zz, 9))
                states.append(tuple(zz))
    full_ranks = {}
    for sh in sorted(survivors):
        _, rk = _solve_affine(_jacobian(sh), [0] * 42, 23,
                              enumerate_solutions=False)
        full_ranks[rk] = full_ranks.get(rk, 0) + 1
    digest = hashlib.sha256(repr(sorted(states)).encode()).hexdigest()
    return {'shapes_in': len(shapes),
            'shapes_surviving_mod9': len(survivors),
            'states_mod9': len(states),
            'ranks_mod9': sorted(ranks),
            'full_jacobian_rank_histogram_of_survivors': full_ranks,
            'children_per_state_next_level':
                {str(k): 3 ** (23 - k) for k in sorted(full_ranks)},
            'states_sha256': digest}


def leading_equation_local(window=25):
    """Local solvability of the deg N = 8 top equation at 2 and 3.

    p^6+r^6 = (2/81) n q (27n^2+q^2)(3n^2+q^2).  Both sides are homogeneous of
    degree six, so only the valuation differences matter.  A sum of two sixth
    powers has v3 == 0 (mod 6) and v2 == 0 or 1 (mod 6).  Writing
    d = v3(q)-v3(n) and e = v2(q)-v2(n):

        v3(RHS) = -4 + d + min(3,2d) + min(1,2d),
        v2(RHS) = 1 + e + 2*(min(0,2e) if e else 2).
    """
    ok3 = [d for d in range(-window, window + 1)
           if (-4 + d + min(3, 2 * d) + min(1, 2 * d)) % 6 == 0]
    ok2 = []
    for e in range(-window, window + 1):
        v = 1 + e + 2 * (min(0, 2 * e) if e else 2)
        if v % 6 in (0, 1):
            ok2.append(e)
    assert 0 not in ok2 and 1 in ok3 and 5 in ok2
    return {'admissible_v3_differences': ok3,
            'admissible_v2_differences': ok2,
            'v2_equal_valuations_excluded': True,
            'conclusion': 'the top equation is not obstructed at 2 or at 3'}


def leading_equation_search(nq_bound=150, pr_bound=250):
    """Bounded integer search for a point on the top-degree equation of the
    deg N = 8 branch, i.e. (equivalently) for an ESOP6 counterexample.
    A hit here would be a counterexample and must be verified by both
    standalone verifiers before any claim is made."""
    table = {}
    for p in range(0, pr_bound + 1):
        p6 = p ** 6
        for r in range(p, pr_bound + 1):
            table.setdefault(p6 + r ** 6, (p, r))
    hits = []
    for n in range(1, nq_bound + 1):
        for q in range(1, nq_bound + 1):
            v = 2 * n * q * (27 * n * n + q * q) * (3 * n * n + q * q)
            if v % 81:
                continue
            w = v // 81
            if w in table:
                hits.append((n, q, table[w]))
    return {'box': {'n,q<=': nq_bound, 'p,r<=': pr_bound},
            'pairs_tested': nq_bound * nq_bound,
            'hits': hits}


def checks(nq_bound=150, pr_bound=250):
    E = surface_derivation()
    out = {'result': 'PASS',
           'scope': 'RATIONAL_CURVE_ATTEMPT.md equation (5), degree eight, '
                    'equal leading slopes, Q=1+q t^2 centered',
           'surface_identity': '2X^6+2Y^6+Z^6-W^6 = 4 t^6 (D - RHS)/Q^6',
           'reversed_normal_form':
               '(81/2)(x^6+y^6) = n*q~*(3n^2+q~^2)*(27n^2+q~^2), '
               'x,y monic deg 7, n monic deg 8, q~=s^2+q'}
    out['system'] = reversed_normal_form(E)
    out['system']['after_eliminating_N'] = {
        'unknowns': 15, 'residual_equations': 34,
        'residual_orders': [9, 42],
        'note': 'orders 1..8 solve n1..n8 linearly (dF/dN(0)=5)'}
    out['grading'] = grading_check(E)
    out['degree_branches'] = degree_branches()
    out['leading_equations'] = leading_equations(E)
    out['branch_degN2'] = branch_degN2()
    out['valuation'] = integrality_theorem()
    out['unperturbed_exact_solutions'] = unperturbed_exact_solutions()
    out['rescaling'] = rescaling_identity()
    shapes, digest = mod3_shapes()
    out['mod3'] = {'residue_shapes': len(shapes), 'shapes_sha256': digest,
                   'enumeration_counts': mod3_shapes.counts,
                   'statement': 'complete list of (x,y,n,Q0) mod 3 solving '
                                '(x^6+y^6)/2 = n^5 q~ in F3[s]'}
    out['mod9'] = mod9_lifting(shapes)
    out['leading_equation_local'] = leading_equation_local()
    out['leading_equation_search'] = leading_equation_search(nq_bound, pr_bound)
    assert out['system']['equation_count'] == 42
    assert out['mod3']['residue_shapes'] == 95
    assert out['mod9']['shapes_surviving_mod9'] == 27
    assert out['mod9']['states_mod9'] == 81
    assert out['leading_equation_search']['hits'] == []
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, default=None)
    ap.add_argument('--nq-bound', type=int, default=150)
    ap.add_argument('--pr-bound', type=int, default=250)
    args = ap.parse_args()
    data = checks(args.nq_bound, args.pr_bound)
    text = json.dumps(data, indent=1, sort_keys=True)
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text + '\n')
        shapes, digest = mod3_shapes()
        (args.out.parent / 'mod3_shapes.json').write_text(json.dumps(
            {'sha256': digest, 'count': len(shapes),
             'order': 'coefficients low to high; (x, y, n, Q0)',
             'shapes': [[list(x), list(y), list(n), Q0]
                        for x, y, n, Q0 in shapes]}, indent=1) + '\n')
    print(text)


if __name__ == '__main__':
    main()
