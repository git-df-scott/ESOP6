import numpy as np
import conic_solver as cs
from numpy.polynomial import polynomial as P
V = np.load("sols_D_real.npy"); gen=[]
for v in V:
    if not np.all(np.abs(v.imag)<1e-8): continue
    vr = v.real.reshape(6,3); R,_ = cs.residual_and_jac(v)
    scale = max(np.linalg.norm(P.polypow(vr[i],6)) for i in range(6))
    sv = np.linalg.svd(vr, compute_uv=False); norms=[np.linalg.norm(vr[i]) for i in range(6)]
    if np.linalg.norm(R)/scale < 1e-11 and sv[2]/sv[0] > 1e-3 and min(norms)/max(norms) > 1e-2: gen.append(vr)
cnt = {"design(z=0)":0, "boundary(4-point+const)":0, "symmetric general":0, "OTHER":0}
for vr in gen:
    x = (vr[:,0]-vr[:,2])/2; y = vr[:,1]/2; z = (vr[:,0]+vr[:,2])/2; g = z[5]
    zeta = (x[:5]+1j*y[:5])/g; zz = z[:5]/g
    nz = [i for i in range(5) if abs(zz[i]) > 1e-6]
    if len(nz) == 0: cnt["design(z=0)"] += 1; continue
    if len(nz) == 1 and abs(zeta[nz[0]]) < 1e-5: cnt["boundary(4-point+const)"] += 1; continue
    ok = False
    if len(nz) == 2:
        i,j = nz
        # pattern: zeta_j = -zeta_i (or equal with z_j = -z_i), and the other three: one at 90deg to zeta_i, two at +-alpha with equal modulus
        same = (abs(zeta[i]+zeta[j])<1e-6 and abs(zz[i]-zz[j])<1e-6) or (abs(zeta[i]-zeta[j])<1e-6 and abs(zz[i]+zz[j])<1e-6)
        if same:
            others = [k for k in range(5) if k not in nz]
            rel = [zeta[k]/(zeta[i]/abs(zeta[i])) for k in others]   # rotate so zeta_i is real
            # find the one perpendicular
            perp = [k for k,r in zip(others,rel) if abs(r.real) < 1e-6]
            rest = [r for k,r in zip(others,rel) if abs(r.real) >= 1e-6]
            if len(perp) == 1 and len(rest) == 2 and abs(abs(rest[0])-abs(rest[1]))<1e-6 and (abs(rest[0]-np.conj(rest[1]))<1e-6 or abs(rest[0]+np.conj(rest[1]))<1e-6):
                ok = True
    cnt["symmetric general" if ok else "OTHER"] += 1
    if not ok: print("OTHER:", np.round(zeta,4).tolist(), np.round(zz,4).tolist())
print(cnt)
