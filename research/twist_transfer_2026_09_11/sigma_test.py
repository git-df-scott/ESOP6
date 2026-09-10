#!/usr/bin/env python3
"""Test a parametric (3,3) family of sixth powers for the twist symmetry.

Setup.  sigma in PGL2(Q) an involution with irrational fixed points is
    M = [[alpha, beta], [gamma, -alpha]],   M^2 = Delta*I,  Delta = alpha^2+beta*gamma,
Delta > 0 and not a rational square (<=> real irrational fixed points).
If F is a binary form of degree N over Q with F o M = mu F then
    F = (1/mu^2) F o M^2 = (Delta^N/mu^2) F   =>  mu^2 = Delta^N.
Hence mu in Q forces N even and mu = +- Delta^(N/2).  Consequences:

  (P1)  NO nonzero rational binary form of ODD degree is a sigma-eigenvector.
        Any family whose coordinate forms have odd degree is immediately dead.
  (P2)  An anti-invariant form (mu = -Delta^(N/2)) vanishes at BOTH fixed points
        of sigma, so it is divisible by the fixed-point quadratic
        Phi = gamma*s^2 - 2*alpha*s*t - beta*t^2   (disc 4*Delta, irreducible/Q).
        Therefore the two anti-invariant coordinates share the common quadratic
        factor Phi.  -> cheap necessary filter: gcd of the two same-side forms
        must contain an irreducible quadratic with positive non-square disc.
"""
import itertools
import sympy as sp

s, t, al, be = sp.symbols('s t alpha beta')


def homogenize(polys, var):
    N = max(sp.Poly(p, var).degree() for p in polys)
    out = []
    for p in polys:
        P = sp.Poly(p, var)
        out.append(sp.expand(sum(c[0] and 0 or 0 for c in []) +
                             sum(co * s**m * t**(N - m)
                                 for (m,), co in P.terms())))
    return out, N


def gcd_filter(l, r, var=None):
    """Pairwise gcd of same-side forms; report quadratic factors and discriminants."""
    res = []
    for side, trip in (('L', l), ('R', r)):
        for i, j in itertools.combinations(range(3), 2):
            g = sp.gcd(sp.Poly(trip[i], s, t), sp.Poly(trip[j], s, t)).as_expr()
            res.append((side, i, j, sp.factor(g)))
    return res


def sigma_systems(forms, N, signs):
    """Return the polynomial system in (alpha,beta) for gamma=1 and given signs."""
    Delta = al**2 + be
    eqs = []
    for F, eps in zip(forms, signs):
        lhs = sp.expand(F.subs({s: al*s + be*t, t: s - al*t}, simultaneous=True))
        rhs = sp.expand(eps * Delta**sp.Rational(N, 2) * F)
        d = sp.Poly(sp.expand(lhs - rhs), s, t)
        eqs.extend(d.coeffs())
    return [sp.expand(e) for e in eqs]


def run(name, l, r, var):
    print('=' * 70)
    print(name)
    forms, N = homogenize(list(l) + list(r), var) if var not in ((s, t),) else (list(l) + list(r), max(sp.Poly(p, s, t).total_degree() for p in list(l)+list(r)))
    print('  coordinate form degree N =', N, ('(ODD -> dead by P1)' if N % 2 else ''))
    L, R = forms[:3], forms[3:]
    ident = sp.expand(sum(u**6 for u in L) - sum(u**6 for u in R))
    print('  homogenised identity holds:', ident == 0)
    print('  --- gcd filter (P2): same-side pairwise gcds ---')
    for side, i, j, g in gcd_filter(L, R):
        print(f'    {side}[{i},{j}] gcd = {g}')
    if N % 2:
        print('  -> skipping sigma search (odd degree, impossible over Q)')
        return
    print('  --- full sigma search (gamma=1), all 6 same-side anti-invariant pairs ---')
    for side in (0, 1):
        for pair in itertools.combinations(range(3), 2):
            signs = [1] * 6
            for k in pair:
                signs[3 * side + k] = -1
            eqs = sigma_systems(forms, N, signs)
            gb = sp.groebner(eqs, al, be, order='lex')
            tag = f'anti = {"LR"[side]}{pair}'
            if list(gb.exprs) == [sp.Integer(1)]:
                print(f'    {tag}: NO SOLUTION (Groebner basis = <1>)')
            else:
                print(f'    {tag}: GB = {list(gb.exprs)}')
                sols = sp.solve(eqs, [al, be], dict=True)
                print('        solutions:', sols)


if __name__ == '__main__':
    import families as F
    l, r, v = F.brudno_delorme_deg4(); run('Brudno-Delorme deg 4', l, r, v)
    l, r, v = F.delorme_deg5();        run('Delorme deg 5', l, r, v)
    l, r, v = F.pell_piezas(sub=True); run('Piezas/Pell (x^2-6y^2=1), deg 4 in (s,t)', l, r, v)
