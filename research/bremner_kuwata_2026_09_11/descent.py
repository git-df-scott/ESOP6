#!/usr/bin/env python3
"""(3,3) -> (5,1) descent analysis for Bremner's K3 K_B and Choudhry's K3.

All arithmetic exact (sympy / python ints).  Run:  python3 descent.py
"""
import sympy as sp
from itertools import product, permutations

x, y, z, u, v, w = sp.symbols('x y z u v w')
vars6 = (x, y, z, u, v, w)
ok = lambda name, cond: print(("PASS " if cond else "FAIL ") + name)

# Bremner's net, Kuwata (1.3)
E1 = (x**2 + x*u - u**2) - (w**2 + w*z - z**2)
E2 = (y**2 + y*v - v**2) - (u**2 + u*x - x**2)
E3 = (z**2 + z*w - w**2) - (v**2 + v*y - y**2)
net = [E1, E2, E3]
F33 = x**6 + y**6 + z**6 - u**6 - v**6 - w**6
Q2 = x**2 + y**2 + z**2 - u**2 - v**2 - w**2          # Kuwata (1.2), = V_3

# ---------------------------------------------------------------- A
# Which diagonal sign involutions preserve the net, and do any of them
# realise the (3,3) -> (5,1) i-twist?
mons = [p*q for i, p in enumerate(vars6) for q in vars6[i:]]
cv = lambda f: [sp.Poly(sp.expand(f), *vars6).coeff_monomial(m) for m in mons]
N = sp.Matrix([cv(e) for e in net]).T
in_net = lambda f: sp.Matrix.hstack(N, sp.Matrix(cv(f))).rank() == N.rank()
aut = []
for eps in product([1, -1], repeat=6):
    if eps[0] == -1:
        continue                                       # projective: fix global sign
    sub = {g: e*g for e, g in zip(eps, vars6)}
    if all(in_net(sp.expand(e.subs(sub, simultaneous=True))) for e in net):
        aut.append(eps)
pat = [1, 1, 1, -1, -1, -1]
gives51 = lambda eps: abs(sum(p*e for p, e in zip(pat, eps))) == 4
print("--- A. diagonal sign automorphisms of Bremner's net (mod global sign)")
for e in aut:
    T = [c for c, g in zip('xyzuvw', e) if g == -1]
    print("     flip %-22s twisted sextic pattern: %s"
          % (T or ['nothing'], "(5,1)" if gives51(e) else "(3,3)"))
ok("no sign automorphism of K_B realises the (5,1) twist",
   not any(gives51(e) for e in aut))
need = sorted({tuple('xyzuvw'[i] for i in range(6) if T[i])
               for T in product([0, 1], repeat=6)
               if abs(sum(pat[i]*(-1 if T[i] else 1) for i in range(6))) == 4},
              key=lambda s_: (len(s_), s_))
print("     flip-sets that WOULD give (5,1) (2-element ones):",
      [t for t in need if len(t) == 2])

# ---------------------------------------------------------------- B
# The sigma-real locus of K_B (sigma = conjugation twisted by the i-scaling).
X, Y, Z, U, V, W = sp.symbols('X Y Z U V W', real=True)
print("--- B. sigma-real locus of K_B")
for label, sub in [("twist (u,v)", {x: X, y: Y, z: Z, u: sp.I*U, v: sp.I*V, w: W}),
                   ("twist (x,y)", {x: sp.I*X, y: sp.I*Y, z: Z, u: U, v: V, w: W})]:
    print("   ", label)
    for k, e in enumerate(net, 1):
        re, im = sp.expand(e.subs(sub)).as_real_imag()
        print("      E%d  Re: %-34s Im: %s"
              % (k, sp.expand(re), sp.expand(im)))
re2 = sp.expand(E2.subs({x: X, y: Y, z: Z, u: sp.I*U, v: sp.I*V, w: W})).as_real_imag()[0]
ok("Re(E2) is positive definite in (X,Y,U,V)  =>  X=Y=U=V=0",
   sp.expand(re2 - (X**2 + Y**2 + U**2 + V**2)) == 0)

# ---------------------------------------------------------------- C
# V_3 under the twist is the (2,1,5) multigrade, which has no positive points.
tw = sp.expand(Q2.subs({x: X, y: Y, z: Z, u: sp.I*U, v: sp.I*V, w: W}))
print("--- C. V_3 (Sum x^2 = Sum u^2) twisted:", tw, "= 0")
ok("twist of V_3 is X^2+Y^2+Z^2+U^2+V^2 = W^2 (the (2,1,5) multigrade)",
   sp.expand(tw - (X**2+Y**2+Z**2+U**2+V**2-W**2)) == 0)
# positivity lemma, numeric sanity on the strict inequality
import random
bad = 0
for _ in range(20000):
    a = [random.randint(1, 10**4) for _ in range(5)]
    f2 = sum(t*t for t in a)
    if sum(t**6 for t in a) >= f2**3:
        bad += 1
ok("Sum a_i^6 < (Sum a_i^2)^3 for 20000 random positive 5-tuples", bad == 0)

# ---------------------------------------------------------------- D
# Census: which (6,3,3) solutions lie on K_B, up to signed relabelling?
SOLS = [((3,19,22),(10,15,23)), ((15,52,65),(36,37,67)), ((23,54,73),(33,47,74)),
        ((3,55,80),(32,43,81)), ((11,65,78),(37,50,81)), ((40,125,129),(51,113,136)),
        ((25,62,138),(82,92,135)), ((1,132,133),(71,92,147)), ((26,169,225),(111,121,230)),
        ((14,163,243),(75,142,245)), ((11,188,243),(103,148,249)), ((29,197,261),(131,139,267)),
        ((36,179,275),(65,169,276)), ((113,241,282),(173,186,293)), ((1,173,294),(75,154,295)),
        ((68,249,289),(125,211,300)), ((83,211,300),(124,185,303)), ((37,199,309),(99,173,311)),
        ((23,282,311),(107,243,326)), ((1,55,330),(159,268,311)), ((92,311,317),(124,277,337)),
        ((27,317,356),(127,271,372)), ((75,226,398),(102,149,400)), ((43,371,372),(140,307,405)),
        ((113,192,410),(303,340,368)), ((120,367,431),(152,345,439)), ((23,432,479),(127,393,496)),
        ((148,299,507),(177,281,508)), ((34,229,518),(421,424,448)), ((53,176,524),(230,395,506)),
        ((93,409,512),(271,293,528)), ((281,419,534),(317,400,537)), ((1,500,515),(197,409,556)),
        ((178,461,543),(271,387,562))]
def on_net(P):                       # plain integer arithmetic, no sympy
    X_,Y_,Z_,U_,V_,W_ = P
    return (X_*X_+X_*U_-U_*U_ == W_*W_+W_*Z_-Z_*Z_ and
            Y_*Y_+Y_*V_-V_*V_ == U_*U_+U_*X_-X_*X_ and
            Z_*Z_+Z_*W_-W_*W_ == V_*V_+V_*Y_-Y_*Y_)
def reachable(A, B):
    base = list(A) + list(B)
    for perm in permutations(range(6)):
        for sg in product([1, -1], repeat=6):
            if on_net([sg[i]*base[perm[i]] for i in range(6)]):
                return [sg[i]*base[perm[i]] for i in range(6)]
    return None
print("--- D. census of (6,3,3) solutions with max coordinate <= 600")
n12 = n13 = 0
for A, B in SOLS:
    q = sum(i*i for i in A) - sum(i*i for i in B)
    r = reachable(A, B) if q == 0 else None
    n12 += (q == 0); n13 += (r is not None)
    print("   %-22s %-22s (1.2):%s  (1.3):%s" %
          (A, B, "yes" if q == 0 else "NO ", r if r else "no"))
print("   totals: %d solutions, %d satisfy (1.2), %d lie on K_B" % (len(SOLS), n12, n13))

# ---------------------------------------------------------------- E
# Choudhry's surface (r = 1, 2, 6) under the twist.
print("--- E. Choudhry's K3 (r=1,2,6) under the twist (u,v) -> (iU, iV)")
L = x + y + z - u - v - w
re, im = sp.expand(L.subs({x: X, y: Y, z: Z, u: sp.I*U, v: sp.I*V, w: W})).as_real_imag()
print("     linear condition:  Re:", sp.expand(re), "  Im:", sp.expand(im))
print("     so V = -U, and the r=2 condition becomes X^2+Y^2+Z^2+2U^2 = W^2,")
print("     i.e. the (2,1,5) multigrade again -> no positive points.")
