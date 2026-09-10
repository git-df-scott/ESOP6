#!/usr/bin/env python3
"""Exact verification of the Bremner(1981)/Kuwata(2007) machinery and of the
(3,3) -> (5,1) descent analysis.  sympy, exact rationals only."""
import sympy as sp
from itertools import product

x,y,z,u,v,w = sp.symbols('x y z u v w')
s,t = sp.symbols('s t')
ok = lambda name,cond: print(("PASS " if cond else "FAIL ")+name)

# ---------- 1. Bremner's net (Kuwata eq (1.3)) ----------
E1 = (x**2+x*u-u**2) - (w**2+w*z-z**2)
E2 = (y**2+y*v-v**2) - (u**2+u*x-x**2)
E3 = (z**2+z*w-w**2) - (v**2+v*y-y**2)
F33 = x**6+y**6+z**6-u**6-v**6-w**6
Q2  = x**2+y**2+z**2-u**2-v**2-w**2          # eq (1.2)

ok("(1.2) is the sum of the three equations of (1.3)", sp.expand(E1+E2+E3-2*Q2)==0)

# the engine identity
q = lambda a,b: a**2+a*b-b**2
ok("q(s,t)^3 - q(t,s)^3 = 2(s^6-t^6)", sp.expand(q(s,t)**3-q(t,s)**3-2*(s**6-t**6))==0)

# (1.3) => (1.1): cube and add
A1,A2,A3 = q(x,u),q(y,v),q(z,w)
B1,B2,B3 = q(u,x),q(v,y),q(w,z)
ok("A_i - B_{i-1} are exactly E1,E2,E3 (cyclic matching)",
   sp.expand(A1-B3-E1)==0 and sp.expand(A2-B1-E2)==0 and sp.expand(A3-B2-E3)==0)
ok("sum A_i^3 - sum B_i^3 = 2*F33",
   sp.expand(A1**3+A2**3+A3**3-B1**3-B2**3-B3**3-2*F33)==0)
# F33 lies in the ideal of the net: explicit cofactors
cof = sp.groebner([E1,E2,E3],x,y,z,u,v,w,order='grevlex')
ok("F33 reduces to 0 modulo the net ideal", sp.simplify(cof.reduce(F33)[1])==0)

# ---------- 2. Kuwata's change of basis (3.7) ----------
q1 = -y**2+v**2+x*u-z*w
q2 = x**2+y**2-z**2-u**2-v**2+w**2-2*x*u+2*y*v
q3 = -4*x**2+4*u**2-4*y*v+4*z*w
ok("q1 = (E1-E2-E3)/2",  sp.expand(q1-(E1-E2-E3)/2)==0)
ok("q2 = (-E1+3E2-E3)/2",sp.expand(q2-(-E1+3*E2-E3)/2)==0)
ok("q3 = -2E1-2E2+2E3",  sp.expand(q3-(-2*E1-2*E2+2*E3))==0)
x0,x1,x2,x3,x4,x5 = sp.symbols('x0 x1 x2 x3 x4 x5')
sub = {x:x0, y:(x1+x4)/2, v:(x1-x4)/2, z:x2, w:-x3, u:x5}   # inverse of (3.7)
ok("q1 -> Plucker x0x5-x1x4+x2x3", sp.expand(q1.subs(sub)-(x0*x5-x1*x4+x2*x3))==0)
ok("q2 -> F (the quadratic line complex)",
   sp.expand(q2.subs(sub)-(x0**2+sp.Rational(1,2)*x1**2-x2**2+x3**2
                           -sp.Rational(1,2)*x4**2-x5**2-2*x0*x5+x1*x4))==0)
ok("q3 -> -4x0^2-x1^2+x4^2+4x5^2-4x2x3",
   sp.expand(q3.subs(sub)-(-4*x0**2-x1**2+x4**2+4*x5**2-4*x2*x3))==0)

# ---------- 3. Kuwata Table 3: the nine smallest solutions satisfy (1.3) ----------
# NB: Table 3 row 3 is printed as (33:47:74:-73:54:23); that fails (1.3).
# The sign-corrected (33:47:74:-73:54:-23) is on K_B.
T3 = [(3,19,22,-23,10,-15),(36,37,67,-65,-15,-52),(33,47,74,-73,54,-23),
      (32,43,81,-3,80,55),(37,50,81,-78,-11,-65),(51,113,136,-125,-40,-129),
      (71,92,147,-132,133,-1),(111,121,230,26,225,169),(75,142,245,14,243,163)]
allok=True
for P in T3:
    d = dict(zip((x,y,z,u,v,w),P))
    vals=[E1.subs(d),E2.subs(d),E3.subs(d),F33.subs(d),Q2.subs(d)]
    if any(V!=0 for V in vals): allok=False; print("  offending",P,vals)
ok("all 9 Table-3 points satisfy (1.3), (1.2) and (1.1)", allok)

# ---------- 4. sign involutions preserving the net ----------
net = [E1,E2,E3]
mons = [a*b for i,a in enumerate((x,y,z,u,v,w)) for b in (x,y,z,u,v,w)[i:]]
def coeffvec(f):
    f=sp.expand(f); return [sp.Poly(f,x,y,z,u,v,w).coeff_monomial(m) for m in mons]
N = sp.Matrix([coeffvec(e) for e in net]).T           # 21 x 3, the net
def in_net(f):
    return sp.Matrix.hstack(N, sp.Matrix(coeffvec(f))).rank()==N.rank()
keep=[]
for eps in product([1,-1],repeat=6):
    if eps[0]==-1: continue                            # projective: fix global sign
    d = dict(zip((x,y,z,u,v,w),[e*g for e,g in zip(eps,(x,y,z,u,v,w))]))
    if all(in_net(sp.expand(e.subs(d,simultaneous=True))) for e in net):
        keep.append(eps)
print("sign involutions preserving Bremner's net (up to global sign):")
for e in keep: print("   flip",[c for c,s_ in zip('xyzuvw',e) if s_==-1] or ['(none)'])

# which 2-element flip sets turn the (3,3) sign pattern into (5,1)?
pat = [1,1,1,-1,-1,-1]
def is51(p): return abs(sum(p)) == 4
need=[]
for T in product([0,1],repeat=6):
    p=[pat[i]*(-1 if T[i] else 1) for i in range(6)]
    if is51(p): need.append(tuple('xyzuvw'[i] for i in range(6) if T[i]))
print("flip-sets T realising the (5,1) twist:", sorted(set(need),key=len)[:12],"...")
print("  (T is always 2 coords from one triple, or its 4-element complement)")
