#!/usr/bin/env python3
"""Known parametric families for equal sums of sixth powers, exactly transcribed.

Each family returns (lhs, rhs, var) where lhs, rhs are 3-tuples (or (2,4)) of
sympy expressions in `var` with  sum(l**6) == sum(r**6)  identically.
Sources recorded in FAMILIES.md.
"""
import sympy as sp

n, s, t, X, Y, Z = sp.symbols('n s t X Y Z')


def brudno_delorme_deg4():
    """S. Brudno (1976) / J. Delorme, as given on Piezas 023 (deg 4, k=2 and 6)."""
    a = -n**4 - n**3 - 5*n**2 + 8*n + 8
    b = (n**3 + 7*n - 2)*(n + 2)
    c = 3*(3*n**2 + 2*n + 4)
    d = (n**2 - n + 3)*(n + 2)**2
    e = -4*n**3 - 5*n**2 - 8*n + 8
    f = -n**4 + n**2 + 14*n + 4
    return (a, b, c), (d, e, f), n


def delorme_deg5():
    """J. Delorme (Math. Comp. 59 (1992)), deg 5 identity as given on Piezas 023."""
    x1 = 3*n**5 + 8*n**4 + 9*n**3 - 4*n**2 - 9*n - 2
    x2 = -2*n**5 - n**4 + 12*n**3 + 13*n**2 + 4*n - 1
    x3 = -n**5 - 9*n**4 - 13*n**3 - 7*n**2 - 7*n - 3
    y1 = 2*n**5 + 9*n**4 + 4*n**3 - 9*n**2 - 8*n - 3
    y2 = -3*n**5 - 7*n**4 - 7*n**3 - 13*n**2 - 9*n - 1
    y3 = -n**5 + 4*n**4 + 13*n**3 + 12*n**2 - n - 2
    return (x1, x2, x3), (y1, y2, y3), n


def pell_piezas(sub=True):
    """Piezas:  (ad+b)^k+(c-2d)^k+(de+f)^k = (de-f)^k+(c+2d)^k+(ad-b)^k, k=2,6
    with {a,b,c,d,e,f} = {-4x+11y, y, 2(x-3y)(x-2y)+y^2, x-y, 2x-3y, 2x-5y},
    x^2-6y^2 = 1.  Homogenised as X^2-6Y^2 = Z^2 with
    (X,Y,Z) = (s^2+6t^2, 2st, s^2-6t^2); each coordinate is multiplied by Z^2."""
    A = (-4*X + 11*Y)
    B = Y
    C = 2*(X - 3*Y)*(X - 2*Y) + Y**2
    D = (X - Y)
    E = (2*X - 3*Y)
    F = (2*X - 5*Y)
    # a*d+b  ->  A*D + B*Z   (weight 2);  c-2d -> C - 2*D*Z ; d*e+f -> D*E + F*Z
    l = (A*D + B*Z, C - 2*D*Z, D*E + F*Z)
    r = (D*E - F*Z, C + 2*D*Z, A*D - B*Z)
    if sub:
        rep = {X: s**2 + 6*t**2, Y: 2*s*t, Z: s**2 - 6*t**2}
        l = tuple(sp.expand(u.subs(rep)) for u in l)
        r = tuple(sp.expand(u.subs(rep)) for u in r)
        return l, r, (s, t)
    return l, r, (X, Y, Z)


def brudno_kaplansky():
    """S. Brudno (1970); Brudno-Kaplansky (1974):
    a^k+b^k+c^k = d^k+(b+c)^k+(b-c)^k, k=2,6, with
    {c,d} = {a*x/(a^2+9b^2), 3*b*x/(a^2+9b^2)} and (a^2-b^2)(a^2+9b^2) = x^2.
    NOT rational (elliptic); returned in raw form for checking only."""
    a, b, x = sp.symbols('a b x')
    den = a**2 + 9*b**2
    c = a*x/den
    d = 3*b*x/den
    return (a, b, c), (d, b + c, b - c), (a, b, x), sp.Eq((a**2 - b**2)*den, x**2)


def piezas_24_quadgrade():
    """Piezas:  (a+d)^k+(-a+d)^k+(b+c)^k+(b-c)^k = 2(c+d)^k+2(-c+d)^k, k=2,6
    with 5a^2+2b^2 = 7c^2 and 2a^2+5b^2 = 7d^2.  A (6,4,4); kept for record."""
    a, b, c, d = sp.symbols('a b c d')
    return (a + d, -a + d, b + c, b - c), (c + d, -c + d, c + d, -c + d), (a, b, c, d)
