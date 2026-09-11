#!/usr/bin/env python3
"""Degree-six boundary contact with UNEQUAL leading tangent slopes.

Scope:

    X = t A5(t),  Y = rho t B5(t),  rho in Q, rho>0, rho != 1,
    Z with Z(0)=1, W = Z + c t^6,
    A5(0)=B5(0)=1, deg A5, deg B5 <= 5, deg Z <= 6.

This is the sibling left explicitly open by BOUNDARY_CONTACT_6.md, whose
proof covers only rho=1.  Everything here is exact.

Results produced by ``checks()``:

  * the leading contact constant is c = (1+rho^6)/3 (c=2/3 only for rho=1);
  * with M = (W+Z)/2, g = A5, f = rho B5 the exact divided identity is

        g^6+f^6 = (c/16) M (12M^2+c^2 t^12)(4M^2+3c^2 t^12);            (6)

  * the degree argument leaves exactly two branches, deg M = 0 with
    max(deg g,deg f)=4 and deg M = 6 with max(deg g,deg f)=5;
  * KEY LEMMA: v3(c) = min(0,6*v3(rho)) - 1 == 5 (mod 6) for every nonzero
    rational rho;
  * BRANCH deg M = 0 IS IMPOSSIBLE for every rational rho>0 (3-adic Gauss
    valuation).  This generalizes the rho=1 obstruction a4^6+b4^6=2/81;
  * in the branch deg M = 6, the reversed equation has no solution whose
    reversed data is 3-integral, again for every rational rho>0.

Not proved here: the branch deg M = 6 with 3-denominators, i.e. the full
analogue of the modulus-729 certificate of BOUNDARY_CONTACT_6.md.
"""
import argparse
import json
from fractions import Fraction
from pathlib import Path

import sympy as sp
from sympy import Rational as Rat

T, S, RHO, C, M, LAM = sp.symbols('t s rho c M lam')


def contact_identity():
    """W=M+c t^6/2, Z=M-c t^6/2 and the exact divided equation."""
    W = M + C * T**6 / 2
    Z = M - C * T**6 / 2
    div = sp.expand((W**6 - Z**6) / (2 * T**6))
    expanded = sp.expand(3 * C * M**5 + Rat(5, 2) * C**3 * T**12 * M**3
                         + Rat(3, 16) * C**5 * T**24 * M)
    assert sp.expand(div - expanded) == 0
    factored = sp.expand(C / 16 * M * (12 * M**2 + C**2 * T**12)
                         * (4 * M**2 + 3 * C**2 * T**12))
    assert sp.expand(div - factored) == 0
    # leading contact constant: at t=0, M=1 the right side is 3c and the
    # left side of the curve equation is g(0)^6+f(0)^6 = 1+rho^6
    assert sp.expand(div.subs({M: 1, T: 0}) - 3 * C) == 0
    cc = (1 + RHO**6) / 3
    assert sp.expand((3 * cc) - (1 + RHO**6)) == 0
    # consistency with the equal-slope normalization of BOUNDARY_CONTACT_6.md
    assert sp.expand(factored.subs(C, Rat(2, 3))
                     - Rat(2, 81) * M * (3 * M**2 + T**12)
                     * (27 * M**2 + T**12)) == 0
    return {'c': '(1+rho^6)/3',
            'identity': 'g^6+f^6 = (c/16) M (12M^2+c^2 t^12)(4M^2+3c^2 t^12)',
            'rho=1_reduces_to': '2/81 * M (3M^2+t^12)(27M^2+t^12)'}


def degree_branches():
    """deg(RHS of (6)) = max(5h, 12+3h, 24+h), h = deg M; no cancellation
    because all three leading coefficients carry the sign of the leading
    coefficient of M and c>0.  deg(LHS) = 6r, r = max(deg g,deg f) <= 5."""
    table, admissible = {}, []
    for h in range(0, 7):
        d = max(5 * h, 12 + 3 * h, 24 + h)
        table[h] = d
        if d % 6 == 0 and d // 6 <= 5:
            admissible.append((h, d // 6))
    assert admissible == [(0, 4), (6, 5)], admissible
    return {'rhs_degree_by_degM': table, 'admissible_(degM,r)': admissible}


def c_valuation(rho):
    """v3(c) for c=(1+rho^6)/3."""
    rho = Fraction(rho)
    def v3(x):
        x = Fraction(x)
        if x == 0:
            raise ValueError
        n, d, v = x.numerator, x.denominator, 0
        while n % 3 == 0:
            n //= 3
            v += 1
        while d % 3 == 0:
            d //= 3
            v -= 1
        return v
    c = (1 + rho**6) / 3
    return v3(c), min(0, 6 * v3(rho)) - 1


def key_lemma():
    """v3(c) = min(0,6*v3(rho))-1 = 5 (mod 6) for every nonzero rational rho."""
    samples = ['1', '2', '3', '1/3', '5/2', '9/4', '4/9', '27', '1/27',
               '7/5', '243/2', '2/243', '10/3', '3/10', '11/9', '81/16']
    checked = []
    for r in samples:
        got, pred = c_valuation(r)
        assert got == pred, (r, got, pred)
        assert got % 6 == 5, (r, got)
        checked.append([r, got])
    return {'lemma': 'v3(c) = min(0,6*v3(rho))-1 == 5 (mod 6)',
            'samples_checked': checked}


def branch_degM0():
    """deg M = 0 (M=1, r=4): g^6+f^6 = (c/16)(12+c^2 t^12)(4+3c^2 t^12).
    3-adic Gauss valuation: 6*min(w(g),w(f)) = nu + min(2nu,1) + min(1+2nu,0)
    with nu = v3(c) <= -1, which equals 5nu+1 == 2 (mod 6).  No solution."""
    h0 = sp.expand(C / 16 * (12 + C**2 * T**12) * (4 + 3 * C**2 * T**12))
    ref = sp.expand(3 * C + Rat(5, 2) * C**3 * T**12 + Rat(3, 16) * C**5 * T**24)
    assert sp.expand(h0 - ref) == 0

    def balance(nu):
        return nu + min(2 * nu, 1) + min(1 + 2 * nu, 0)

    bad = []
    for k in range(0, 40):
        nu = -1 - 6 * k                      # every admissible v3(c)
        val = balance(nu)
        assert val == 5 * nu + 1
        assert val % 6 != 0
        bad.append([nu, val])
    # sanity: rho=1 gives nu=-1 and the known a4^6+b4^6=2/81 obstruction
    assert balance(-1) == -4
    return {'statement': 'deg M = 0 branch impossible for every rational rho>0',
            'balance_is': '5*nu+1, nu=v3(c) == 5 (mod 6)',
            'nu_values_checked': len(bad),
            'first_five': bad[:5],
            'rho=1_specialization': 'v3(2/81) = -4, not divisible by six'}


def branch_degM6():
    """deg M = 6 (r=5).  Reverse the parameter: x=s^5 g(1/s) (monic, deg 5),
    y=s^5 f(1/s) (leading coefficient rho), m=s^6 M(1/s) (monic, deg 6):

        x^6+y^6 = (c/16) m (12m^2+c^2)(4m^2+3c^2).                      (7)

    3-adic Gauss valuation with a=min(w(x),w(y)), b=w(m)<=0:

        6a = nu + b + min(1+2b,2nu) + min(2b,1+2nu).

    b=0 (3-integral reversed data) forces 6a = 5nu+1 == 2 (mod 6): impossible.
    """
    m, x, y = sp.symbols('m x y')
    lhs = sp.expand(C / 16 * M * (12 * M**2 + C**2 * T**12)
                    * (4 * M**2 + 3 * C**2 * T**12))
    rev = sp.expand(lhs.subs({M: m / S**6, T: 1 / S}) * S**30)
    assert sp.simplify(rev - C / 16 * m * (12 * m**2 + C**2)
                       * (4 * m**2 + 3 * C**2)) == 0

    def balance(a, b, nu):
        return 6 * a - (nu + b + min(1 + 2 * b, 2 * nu) + min(2 * b, 1 + 2 * nu))

    for k in range(0, 20):
        nu = -1 - 6 * k
        assert all(balance(a, 0, nu) != 0 for a in range(-60, 1))
    # what survives with b<0 (not closed here), recorded over a bounded window
    surviving = []
    for k in range(0, 3):
        nu = -1 - 6 * k
        for b in range(-12, 0):
            for a in range(-40, 1):
                if balance(a, b, nu) == 0:
                    surviving.append([nu, b, a])
    return {'reversed_identity':
            'x^6+y^6 = (c/16) m (12m^2+c^2)(4m^2+3c^2)',
            'b=0_excluded_for_every_rho': True,
            'open': 'b<0, i.e. 3-denominators in the reversed data',
            'surviving_(nu,b,a)_in_bounded_window': surviving[:12],
            'surviving_count_in_window': len(surviving)}


def checks():
    out = {'result': 'PASS',
           'scope': 'degree six, X=t A5, Y=rho t B5, rho in Q_{>0}, rho!=1'}
    out['contact'] = contact_identity()
    out['degree_branches'] = degree_branches()
    out['key_lemma'] = key_lemma()
    out['branch_degM0'] = branch_degM0()
    out['branch_degM6'] = branch_degM6()
    out['not_proved'] = ('the deg M = 6 branch with 3-denominators; no '
                         'analogue of the modulus-729 certificate was built')
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, default=None)
    args = ap.parse_args()
    data = checks()
    text = json.dumps(data, indent=1, sort_keys=True)
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text + '\n')
    print(text)


if __name__ == '__main__':
    main()
