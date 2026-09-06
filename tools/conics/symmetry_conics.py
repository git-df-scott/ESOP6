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
nsym=0
for vr in gen:
    x = (vr[:,0]-vr[:,2])/2; y = vr[:,1]/2; z = (vr[:,0]+vr[:,2])/2; g = z[5]
    zeta = (x[:5]+1j*y[:5])/g; zz = z[:5]/g
    # test reflection symmetry: exists axis angle phi s.t. {(e^{-2i phi} conj(zeta_i), z_i)} == {(zeta_i, z_i)}
    best = 1e9
    for i in range(5):
        for j in range(i,5):
            if abs(abs(zeta[i])-abs(zeta[j]))>1e-6 or abs(zz[i]-zz[j])>1e-6: continue
            phi2 = np.angle(zeta[i])+np.angle(zeta[j])   # reflection maps zeta_i <-> zeta_j when 2phi = arg i + arg j
            refl = np.exp(1j*phi2)*np.conj(zeta)
            # match sets
            err = 0
            for a in range(5):
                err += min(abs(refl[a]-zeta[b])+abs(zz[a]-zz[b]) for b in range(5))
            best = min(best, err)
    sym = best < 1e-6
    nsym += sym
    print(f"conic: reflection-symmetric={sym} (err {best:.1e}); |zeta|={np.round(np.abs(zeta),4).tolist()} z={np.round(zz,4).tolist()} nz={np.sum(np.abs(zz)>1e-6)}")
print("reflection-symmetric:", nsym, "of", len(gen))
