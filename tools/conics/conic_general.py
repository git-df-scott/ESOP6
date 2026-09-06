"""
Conics on  p_1^n + ... + p_m^n = p_{m+1}^n  (p_i real quadratics in t).
Normalisation: a_1 = 1 (overall scale), p_{m+1} = g (1 + t^2) (Mobius shift+scale, p_{m+1} definite),
b_m = 0 (rotation stabiliser of 1+t^2).  Unknowns 3m-1, equations 2n+1.
Real Gauss-Newton from random real starts; strict post-filter removes boundary/constant solutions.
usage: conic_general.py n m nstarts seed
"""
import numpy as np, sys
from numpy.polynomial import polynomial as P
n, m, N, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
rng = np.random.default_rng(seed)
L = 2*n+1

def unpack(x):
    v = np.zeros((m+1, 3))
    v[0,0] = 1.0; v[0,1] = x[0]; v[0,2] = x[1]; k = 2
    for i in range(1, m):
        for j in range(3):
            if i == m-1 and j == 1: continue
            v[i,j] = x[k]; k += 1
    g = x[k]; v[m] = [g, 0, g]
    return v

def resjac(x):
    v = unpack(x)
    R = np.zeros(L); Jv = np.zeros((L, m+1, 3))
    for i in range(m+1):
        p = v[i]; pn1 = P.polypow(p, n-1); pn = P.polymul(pn1, p)
        s = 1 if i < m else -1
        R += s*pn
        d = n*pn1
        for j in range(3): Jv[j:j+L-2, i, j] += s*d
    # chain rule to x
    cols = [(0,1),(0,2)] + [(i,j) for i in range(1,m) for j in range(3) if not (i==m-1 and j==1)]
    J = np.zeros((L, len(cols)+1))
    for c,(i,j) in enumerate(cols): J[:,c] = Jv[:,i,j]
    J[:,-1] = Jv[:,m,0] + Jv[:,m,2]
    return R, J, v

def solve(x0, iters=120):
    x = x0.copy()
    for it in range(iters):
        R, J, v = resjac(x); nr = np.linalg.norm(R)
        if nr < 1e-14: return x, nr
        dx = np.linalg.lstsq(J, -R, rcond=None)[0]
        lam = 1.0
        while lam > 1e-5:
            xn = x + lam*dx; Rn,_,_ = resjac(xn)
            if np.linalg.norm(Rn) < nr: x = xn; break
            lam *= 0.5
        else: x = x + 1e-4*dx
        if np.linalg.norm(x) > 1e5: return None, np.inf
    R,_,_ = resjac(x); return x, np.linalg.norm(R)

def strict(v):
    R,_,_ = resjac_v(v)
    scale = max(np.linalg.norm(P.polypow(v[i],n)) for i in range(m+1))
    rel = np.linalg.norm(R)/scale
    sv = np.linalg.svd(v, compute_uv=False)
    span = sv[2]/sv[0] if len(sv) > 2 else 0
    discs = [v[i,1]**2 - 4*v[i,0]*v[i,2] for i in range(m+1)]
    reldef = min(-discs[i]/np.sum(v[i]**2) for i in range(m+1))
    tiny = min(np.linalg.norm(v[i]) for i in range(m+1)) / max(np.linalg.norm(v[i]) for i in range(m+1))
    return rel, span, reldef, tiny

def resjac_v(v):
    R = np.zeros(L)
    for i in range(m+1):
        R += (1 if i < m else -1)*P.polypow(v[i], n)
    return R, None, v

nunk = 3*m - 1
found = []
stats = dict(conv=0, genuine=0, definite=0)
for k in range(N):
    x0 = rng.normal(0, 1.2, nunk)
    x, nr = solve(x0)
    if x is None or nr > 1e-9: continue
    stats["conv"] += 1
    v = unpack(x)
    rel, span, reldef, tiny = strict(v)
    if rel < 1e-11 and span > 1e-3 and tiny > 1e-2:
        stats["genuine"] += 1
        if reldef > 1e-3: stats["definite"] += 1
        if not any(np.linalg.norm(np.abs(v)-np.abs(f)) < 1e-6 for f in found):
            found.append(v)
            print(f"GENUINE conic ({'DEFINITE' if reldef>1e-3 else 'indefinite'}) rel={rel:.1e} span={span:.2e} tiny={tiny:.2e}:")
            for i in range(m+1): print("    p%d = %s" % (i+1, np.round(v[i], 8).tolist()))
print(f"(n={n}, m={m}) starts={N} converged={stats['conv']} genuine={stats['genuine']} definite={stats['definite']} distinct_genuine={len(found)}  unknowns={nunk} eqs={L} expected_dim={nunk-1-L}")
np.save(f"gen_{n}_{m}_{seed}.npy", np.array(found))
