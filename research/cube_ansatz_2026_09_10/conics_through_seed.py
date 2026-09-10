"""Family I conics through the seed quadratic point P=(9+sqrt(-249), 9-sqrt(-249), 14, 18, 0, 22) on x5=0.
Normalize the base point where x5=0 to E2=0 and x5=E2.  Then
  A = 9 E1 + c E2,  Omega = -249 E1^2 + w1 E1E2 + w2 E2^2,  L3 = 14E1 + m3E2, L4 = 18E1 + m4E2, L5 = E2, L6 = 22E1 + m6E2
and the identity 2Phi(A,Omega) + L3^6 + L4^6 + L5^6 - L6^6 = 0 has its E1^6 coefficient already satisfied.
Six remaining coefficient equations in six unknowns (c,w1,w2,m3,m4,m6): zero-dimensional.  Solve numerically
from many random complex starts, cluster, and test rationality/algebraicity."""
import numpy as np, sys, mpmath as mp
def pm(a,b): return np.convolve(a,b)
def pw(a,k):
    r=np.array([1.0+0j])
    for _ in range(k): r=pm(r,a)
    return r
def pad(a,n=7):
    r=np.zeros(n,dtype=complex); r[:len(a)]=a; return r
def full(v):
    c,w1,w2,m3,m4,m6=v
    A=np.array([9,c]); Om=np.array([-249,w1,w2])
    Phi=pad(pw(A,6))+15*pad(pm(pw(A,4),Om))+15*pad(pm(pw(A,2),pw(Om,2)))+pad(pw(Om,3))
    r=2*Phi+pad(pw(np.array([14,m3]),6))+pad(pw(np.array([18,m4]),6))+pad(pw(np.array([0,1]),6))-pad(pw(np.array([22,m6]),6))
    return r
def res(v): return full(v)[1:]
def jac(v,h=1e-6):
    r0=res(v); J=np.zeros((6,6),dtype=complex)
    for k in range(6):
        dv=np.zeros(6,dtype=complex); dv[k]=h; J[:,k]=(res(v+dv)-r0)/h
    return J
def newton(v,iters=200):
    for it in range(iters):
        r=res(v); sc=1+np.abs(v).max()**6
        if np.linalg.norm(r)/sc<1e-14: return v,np.linalg.norm(r)/sc
        try: dv=np.linalg.solve(jac(v),-r)
        except np.linalg.LinAlgError: return v,1
        # damping
        lam=1.0
        while lam>1e-4:
            vn=v+lam*dv
            if np.linalg.norm(res(vn))<np.linalg.norm(r): break
            lam/=2
        v=v+lam*dv
        if np.abs(v).max()>1e6: return v,1
    return v,np.linalg.norm(res(v))/(1+np.abs(v).max()**6)
rng=np.random.default_rng(int(sys.argv[1]) if len(sys.argv)>1 else 0)
sols=[]
starts=int(sys.argv[2]) if len(sys.argv)>2 else 3000
for s in range(starts):
    scale=10**rng.uniform(0,2.5)
    v=scale*(rng.standard_normal(6)+1j*rng.standard_normal(6))
    v,nr=newton(v)
    if nr<1e-12:
        if not any(np.linalg.norm(v-u)<1e-6*(1+np.linalg.norm(u)) for u in sols):
            sols.append(v)
print("distinct solutions found:",len(sols))
assert all(abs(full(v)[0])<1e-6 for v in sols)  # E1^6 coefficient is identically satisfied
np.save("/tmp/claude-0/-home-user/15becd99-5cb2-5637-beeb-15d2ba0dfed3/scratchpad/seed_conics.npy",np.array(sols))
real=[v for v in sols if np.abs(v.imag).max()<1e-8]
print("real solutions:",len(real))
for v in sorted(sols,key=lambda v:np.abs(v.imag).max())[:40]:
    print(np.round(v,6))
