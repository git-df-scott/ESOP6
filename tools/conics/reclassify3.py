import numpy as np, glob
import conic_solver as cs
from numpy.polynomial import polynomial as P
for f in sorted(glob.glob("sols_*_real.npy")):
    V = np.load(f); found=[]
    for v in V:
        if not np.all(np.abs(v.imag)<1e-8): continue
        vr = v.real.reshape(6,3)
        R,_ = cs.residual_and_jac(v)
        scale = max(np.linalg.norm(P.polypow(vr[i],6)) for i in range(6))
        rel = np.linalg.norm(R)/scale
        sv = np.linalg.svd(vr, compute_uv=False); span = sv[2]/sv[0]
        norms = [np.linalg.norm(vr[i]) for i in range(6)]; tiny = min(norms)/max(norms)
        if rel < 1e-11 and span > 1e-3 and tiny > 1e-2:
            found.append(vr); print(f, "GENUINE real conic rel=%.1e span=%.2e tiny=%.2e"%(rel,span,tiny)); print(np.round(vr,6).tolist())
    print(f, "genuine real conics (p6 definite, any signs):", len(found))
