import numpy as np, mpmath as mp
import conic_solver as cs
from numpy.polynomial import polynomial as P
V = np.load("sols_D_real.npy"); gen=[]
for v in V:
    if not np.all(np.abs(v.imag)<1e-8): continue
    vr = v.real.reshape(6,3); R,_ = cs.residual_and_jac(v)
    scale = max(np.linalg.norm(P.polypow(vr[i],6)) for i in range(6))
    sv = np.linalg.svd(vr, compute_uv=False); norms=[np.linalg.norm(vr[i]) for i in range(6)]
    if np.linalg.norm(R)/scale < 1e-11 and sv[2]/sv[0] > 1e-3 and min(norms)/max(norms) > 1e-2: gen.append(vr)
for idx in range(min(8,len(gen))):
    vr = gen[idx]
    x = (vr[:,0]-vr[:,2])/2; y = vr[:,1]/2; z = (vr[:,0]+vr[:,2])/2; g = z[5]
    zeta = (x[:5]+1j*y[:5])/g; zz = z[:5]/g
    print(f"--- conic {idx}: g={g:.6f}")
    for i in range(5): print(f"   zeta{i+1} = {zeta[i].real:+.6f}{zeta[i].imag:+.6f}i  |zeta|={abs(zeta[i]):.6f} arg={np.degrees(np.angle(zeta[i])):+8.3f}  z={zz[i]:+.6f}")
