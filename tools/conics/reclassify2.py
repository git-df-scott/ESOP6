import numpy as np, sys, glob
import conic_solver as cs
def strict(v):
    vr = v.real
    R,_ = cs.residual_and_jac(vr.astype(complex))
    scale = max(np.linalg.norm(np.polynomial.polynomial.polypow(vr[3*i:3*i+3],6)) for i in range(6))
    rel = np.linalg.norm(R)/scale
    M = vr.reshape(6,3); sv = np.linalg.svd(M, compute_uv=False)
    discs = [vr[3*i+1]**2-4*vr[3*i]*vr[3*i+2] for i in range(6)]
    # relative definiteness: -disc / (|a|^2+|b|^2+|c|^2)
    reldef = min(-d/np.sum(vr[3*i:3*i+3]**2) for i,d in enumerate(discs))
    return rel, sv[2]/sv[0], reldef
for f in sorted(glob.glob("sols_*_real.npy")):
    V = np.load(f); cnt=0
    for v in V:
        if not np.all(np.abs(v.imag)<1e-8): continue
        rel, spanratio, reldef = strict(v)
        if rel < 1e-11 and spanratio > 1e-3 and reldef > 1e-3:
            cnt+=1; print(f, "CANDIDATE rel=%.1e span=%.2e reldef=%.2e"%(rel,spanratio,reldef), np.round(v.real,6).tolist())
    print(f, "strict candidates:", cnt, "of", len(V))
