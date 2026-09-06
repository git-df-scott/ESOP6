"""
Direct search for REAL DEFINITE conics on X, parametrised so definiteness is automatic:
  p_i(t) = lam_i^2 * ((t-u_i)^2 + v_i^2),  i=1..5   (lam_i real, v_i != 0)
  p_6(t) = 1 + t^2                                   (normalised by real Mobius + scale)
Identity: sum_i p_i^6 = (1+t^2)^6, 13 real equations, 15 unknowns minus SO(2) rotation.
scipy least_squares from many random starts; a zero residual = a real definite conic.
"""
import numpy as np, sys
from scipy.optimize import least_squares
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv)>1 else 0)
N = int(sys.argv[2]) if len(sys.argv)>2 else 2000
target = np.polynomial.polynomial.polypow([1,0,1], 6)

def coeffs(x):
    lam, u, v = x[0:5], x[5:10], x[10:15]
    tot = np.zeros(13)
    for i in range(5):
        p = np.array([u[i]**2 + v[i]**2, -2*u[i], 1.0]) * lam[i]**2
        tot += np.polynomial.polynomial.polypow(p, 6)
    return tot - target

best = []
hits = 0
for k in range(N):
    x0 = np.concatenate([rng.uniform(0.3,1.2,5), rng.normal(0,1.5,5), rng.uniform(0.2,2.0,5)])
    lb = np.concatenate([np.full(5,0.3), np.full(5,-6.0), np.full(5,0.05)]); ub = np.concatenate([np.full(5,1.0), np.full(5,6.0), np.full(5,6.0)])
    x0 = np.clip(x0, lb+1e-3, ub-1e-3)
    r = least_squares(coeffs, x0, bounds=(lb,ub), xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000)
    c = np.linalg.norm(r.fun)
    best.append(c)
    if c < 1e-11:
        hits += 1
        print("HIT residual", c, "x=", np.round(r.x, 8).tolist(), flush=True)
best.sort()
print(f"starts={N} hits={hits}  best residual norms: {[f'{b:.3e}' for b in best[:12]]}")
