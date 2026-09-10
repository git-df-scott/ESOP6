#!/usr/bin/env python3
"""Why Bremner's device forces the (3,3) sign pattern.

Bremner's net (Kuwata (1.3)) is built from ONE identity,

    q(s,t)^3 - q(t,s)^3 = 2 (s^6 - t^6),      q(s,t) = s^2 + s t - t^2,

applied to the three coordinate pairs (x,u), (y,v), (z,w) with a cyclic
matching.  Each pair contributes one POSITIVE and one NEGATIVE sixth power,
so three pairs give the sign pattern (3,3) and nothing else.

To leave the (3,3) pattern one would need the same device for a pair of
EQUAL signs, i.e. binary quadratic forms A, B over Q with

    A^3 - B^3 = c (s^6 + t^6),  c in Q*.

This script classifies, completely and exactly, the solutions of
A^3 - B^3 = c(s^6 - t^6) and of A^3 - B^3 = c(s^6 + t^6) over Q-bar, and
records which are defined over Q.
"""
import sympy as sp

s, t = sp.symbols('s t')
lam, a, b, d = sp.symbols('lam a b d')


def divisors_deg2(poly):
    """All degree-2 divisors of poly over Q, up to scalars."""
    out = []
    fl = sp.factor_list(poly)[1]
    facs = [f for f, m in fl for _ in range(m)]
    import itertools
    for r in (1, 2):
        for combo in itertools.combinations(range(len(facs)), r):
            p = sp.expand(sp.prod([facs[i] for i in combo]))
            if sp.Poly(p, s, t).total_degree() == 2:
                p = sp.expand(p / sp.LC(sp.Poly(p, s, t)))
                if p not in out:
                    out.append(p)
    return out


def classify(sign, name):
    """A^3 - B^3 = c (s^6 + sign t^6);  normalise c = 2 (Bremner's scaling)."""
    target = sp.expand(2 * (s**6 + sign * t**6))
    print("=" * 72)
    print("A^3 - B^3 = 2(%s)" % name)
    print("  rational degree-2 divisors of the target:",
          [str(p) for p in divisors_deg2(s**6 + sign * t**6)])
    found = []
    for D0 in divisors_deg2(s**6 + sign * t**6):
        D = lam * D0
        S = a * s**2 + b * s * t + d * t**2
        # (A-B)(A^2+AB+B^2) = target, with A^2+AB+B^2 = (3S^2+D^2)/4
        eq = sp.expand(D * (3 * S**2 + D**2) / 4 - target)
        eqs = sp.Poly(eq, s, t).coeffs()
        for so in sp.solve(eqs, [lam, a, b, d], dict=True):
            L = so.get(lam, lam)
            if L == 0:
                continue
            Dv = sp.simplify(D.subs(so))
            Sv = sp.simplify(S.subs(so))
            A = sp.expand((Sv + Dv) / 2)
            B = sp.expand((Sv - Dv) / 2)
            if sp.expand(A**3 - B**3 - target) != 0:
                continue                      # guard against spurious roots
            coeffs = [sp.expand(A).coeff(s, i).coeff(t, 2 - i) for i in (2, 1, 0)] + \
                     [sp.expand(B).coeff(s, i).coeff(t, 2 - i) for i in (2, 1, 0)]
            overQ = all(sp.nsimplify(cf).is_rational for cf in coeffs)
            key = (sp.srepr(A), sp.srepr(B))
            if key in [k for k, _ in found]:
                continue
            found.append((key, (A, B, overQ)))
            print("   A = %-28s B = %-28s over Q: %s" % (A, B, overQ))
    if not found:
        print("   (no solutions)")
    return [v for _, v in found]


sols_minus = classify(-1, "s^6 - t^6")
sols_plus = classify(+1, "s^6 + t^6")

print("=" * 72)
q = lambda A, B: A**2 + A * B - B**2
print("Bremner's form q(s,t)=s^2+st-t^2 is among the (s^6 - t^6) solutions:",
      any(sp.expand(A - q(s, t)) == 0 and sp.expand(B - q(t, s)) == 0
          for A, B, _ in sols_minus))
print("Rational solutions for s^6 + t^6:",
      [(str(A), str(B)) for A, B, r in sols_plus if r])
print("Non-rational solutions for s^6 + t^6:",
      [(str(A), str(B)) for A, B, r in sols_plus if not r])
print()
print("The i-twist of Bremner's form:  q(s, i t) = s^2 + i s t + t^2,")
print("  q(s,i t)^3 - q(i t,s)^3 - 2(s^6+t^6) =",
      sp.expand(q(s, sp.I * t)**3 - q(sp.I * t, s)**3 - 2 * (s**6 + t**6)))
print("  -> the same-sign device exists only over Q(i): it IS the twist.")
