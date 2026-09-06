import numpy as np, glob
import conic_solver as cs
from numpy.polynomial import polynomial as P
V = np.load("sols_D_real.npy")
gen = []
for v in V:
    if not np.all(np.abs(v.imag)<1e-8): continue
    vr = v.real.reshape(6,3); R,_ = cs.residual_and_jac(v)
    scale = max(np.linalg.norm(P.polypow(vr[i],6)) for i in range(6))
    sv = np.linalg.svd(vr, compute_uv=False); norms=[np.linalg.norm(vr[i]) for i in range(6)]
    if np.linalg.norm(R)/scale < 1e-11 and sv[2]/sv[0] > 1e-3 and min(norms)/max(norms) > 1e-2: gen.append(vr)
print("genuine real conics:", len(gen))
ndesign = 0
for vr in gen:
    # p = a + b t + c t^2 = x A + y B + z C with A=(1,0,-1), B=(0,2,0), C=(1,0,1)
    x = (vr[:,0]-vr[:,2])/2; y = vr[:,1]/2; z = (vr[:,0]+vr[:,2])/2
    zrel = np.max(np.abs(z[:5]))/np.max(np.abs(vr[:5]))
    zeta = x[:5] + 1j*y[:5]
    N = np.abs(zeta)**2; om = zeta/np.conj(zeta); W = N**3
    mom = [abs(np.sum(W*om**m)) for m in (1,2,3)]
    e2 = sum(om[i]*om[j] for i in range(5) for j in range(i+1,5))
    g6 = (5/16)*np.sum(N**3)
    if zrel < 1e-6:
        ndesign += 1
        print(f"DESIGN-type conic: |z|rel={zrel:.1e} moments(1,2,3)={np.round(mom,10).tolist()} |e2|={abs(e2):.2e} g^6 check: {g6:.8f} vs {z[5]**6:.8f}")
    else:
        print(f"general conic (z!=0): |z|rel={zrel:.2e}")
print("design-type:", ndesign, "of", len(gen))
