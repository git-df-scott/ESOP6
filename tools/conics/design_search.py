"""
Exact search for RATIONAL 5-node weighted 3-designs on the circle (the z=0 conic stratum of ESOP6).
Nodes omega_i = gamma_i/conj(gamma_i) for primitive Gaussian integers gamma_i (rational points on S^1).
Condition for a nonzero kernel of weights: e2(omega_1..5) = 0.  Fix omega_1 = 1 by rotation.
For each rational config: exact kernel weights w (fractions), check positivity, then the arithmetic
condition: w_i / |gamma_i|^6 all in the same class of Q*/Q*^6  (then a rational conic exists).
Also report the weaker cube condition (w_i/w_j in Q^3) and square condition separately.
"""
import sys, itertools
from fractions import Fraction as Fr
from math import gcd, isqrt
B = int(sys.argv[1]) if len(sys.argv) > 1 else 60
# Gaussian rationals as pairs of Fractions
def gmul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def gadd(a,b): return (a[0]+b[0], a[1]+b[1])
def gconj(a): return (a[0], -a[1])
def gdiv(a,b):
    n = b[0]*b[0]+b[1]*b[1]; c = gmul(a, gconj(b)); return (c[0]/n, c[1]/n)
# nodes: gamma = a+bi primitive, upper half plane (b>0) or (b==0,a==1); omega = gamma^2/N
nodes = []
for a in range(-B, B+1):
    for b in range(0, B+1):
        if a*a+b*b > B*B or (a,b)==(0,0): continue
        if b == 0 and a != 1: continue
        if gcd(abs(a), b) != 1: continue
        N = a*a+b*b
        om = (Fr(a*a-b*b, N), Fr(2*a*b, N))
        nodes.append((om, (a,b), N))
print("rational circle points (nodes):", len(nodes), file=sys.stderr)
one = (Fr(1),Fr(0))
def sixth_class_key(q):
    """q positive Fraction: return the class of q in Q*/Q*^6 as a reduced (num,den) with exponents mod 6."""
    from sympy import factorint
    num, den = q.numerator, q.denominator
    f = {}
    for p,e in factorint(num).items(): f[p] = f.get(p,0)+e
    for p,e in factorint(den).items(): f[p] = f.get(p,0)-e
    return tuple(sorted((p, e % 6) for p,e in f.items() if e % 6))
def kernel_weights(oms):
    # solve sum W_i om_i^k = 0, k=1,2,3 (real and imag) : 6 x 5 rational system; return kernel vector or None
    rows = []
    for k in (1,2,3):
        pw = [one]*5
        for i in range(5):
            p = one
            for _ in range(k): p = gmul(p, oms[i])
            pw[i] = p
        rows.append([pw[i][0] for i in range(5)]); rows.append([pw[i][1] for i in range(5)])
    # gaussian elimination for kernel
    M = [r[:] for r in rows]; n=5; piv=[]; r=0
    for c in range(n):
        pr = next((i for i in range(r,len(M)) if M[i][c]!=0), None)
        if pr is None: continue
        M[r],M[pr]=M[pr],M[r]; pv=M[r][c]; M[r]=[x/pv for x in M[r]]
        for i in range(len(M)):
            if i!=r and M[i][c]!=0:
                f=M[i][c]; M[i]=[x-f*y for x,y in zip(M[i],M[r])]
        piv.append(c); r+=1
    free = [c for c in range(n) if c not in piv]
    if len(free) != 1: return None if not free else "multi"
    fc = free[0]; w=[Fr(0)]*n; w[fc]=Fr(1)
    for i,c in enumerate(piv): w[c] = -M[i][fc]
    return w
found = 0; configs = 0; stats = {}
nl = len(nodes)
for j in range(nl):
    for k in range(j+1, nl):
        for l in range(k+1, nl):
            o = [one, nodes[j][0], nodes[k][0], nodes[l][0]]
            e1 = o[0]; e2 = (Fr(0),Fr(0))
            for i in range(1,4):
                e2 = gadd(e2, gmul(e1, o[i])); e1 = gadd(e1, o[i])
            if e1 == (0,0): continue
            o5 = gdiv((-e2[0], -e2[1]), e1)
            if o5[0]*o5[0]+o5[1]*o5[1] != 1: continue
            if o5 in o: continue
            configs += 1
            oms = o + [o5]
            e1f = gadd(e1, o5); e1abs2 = e1f[0]*e1f[0]+e1f[1]*e1f[1]
            stats.setdefault("e1lt1", 0)
            if e1abs2 < 1: stats["e1lt1"] += 1; print("  |e1|<1 config:", [str(x[0])+"+"+str(x[1])+"i" for x in oms], "|e1|^2=", e1abs2, flush=True)
            w = kernel_weights(oms)
            if w is None or w == "multi": continue
            if any(x == 0 for x in w): continue
            if all(x > 0 for x in w) or all(x < 0 for x in w):
                w = [abs(x) for x in w]
                # |gamma_i|^6 : for o5 need its gamma: o5 = gamma^2/N with gamma primitive: find from o5 = (p/q, r/q)
                gam = [(1,0)] + [nodes[t][1] for t in (j,k,l)]
                # o5 = (c/d) + (s/d) i with c^2+s^2 = d^2 -> gamma = sqrt of (c + s i)*d... use: omega = gamma/conj(gamma) => gamma = 1 + omega up to real scaling (if omega != -1)
                if o5 == (Fr(-1),Fr(0)): gam5 = (0,1)
                else:
                    g = (1+o5[0], o5[1]); den = g[0].denominator*g[1].denominator//gcd(g[0].denominator, g[1].denominator)
                    ga, gb = int(g[0]*den), int(g[1]*den); d = gcd(ga,gb); gam5 = (ga//d, gb//d)
                gam.append(gam5)
                Ns = [a*a+b*b for a,b in gam]
                ratios = [w[i]/Fr(Ns[i])**3 for i in range(5)]
                keys = [sixth_class_key(r) for r in ratios]
                cube_ok = len(set(tuple((p,e%3) for p,e in kk) for kk in keys)) == 1
                sixth_ok = len(set(keys)) == 1
                found += 1
                tag = "SIXTH-CLASS-MATCH" if sixth_ok else ("cube-match" if cube_ok else "")
                if sixth_ok or cube_ok or found <= 15:
                    print(f"config gammas={gam} weights={[str(x) for x in w]} ratios-classes={keys} {tag}", flush=True)
print(f"B={B}: rational e2=0 configurations={configs}, with |e1|<1: {stats.get('e1lt1',0)}, with positive kernel weights={found}", flush=True)
