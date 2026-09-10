"""Z/2-symmetric sextic rational curves on X (x1<->x2 swapped by an involution of P^1).
On the quotient line (E1,E2) with the conic w^2 = omega(E1,E2):
   x1,x2 = A ± R w   (A cubic, R quadratic),   x3,x4,x5,x6 = L_j cubic forms in (E1,E2).
Identity: 2 Phi(A, R^2 omega) + L3^6 + L4^6 + L5^6 - L6^6 = 0, a binary form of degree 18 (19 equations).
Unknowns 4+3+3+16 = 26; symmetries GL2 (4) and (R,omega)->(lam R, omega/lam^2) (1): expected dimension 2.
Slice principle: irreducible L_j give degree-6 slice points; rational roots of L_j give quadratic slice points."""
import numpy as np, sys
from collections import Counter
def pm(a,b): return np.convolve(a,b)
def pw(a,k):
    r=np.array([1.0+0j])
    for _ in range(k): r=pm(r,a)
    return r
def pad(a,n=19):
    r=np.zeros(n,dtype=complex); r[:len(a)]=a; return r
def unpack(v):
    A=v[0:4]; R=v[4:7]; om=v[7:10]; L=[v[10+4*j:14+4*j] for j in range(4)]
    return A,R,om,L
def res(v):
    A,R,om,L=unpack(v)
    Om=pm(pw(R,2),om)
    Phi=pad(pw(A,6))+15*pad(pm(pw(A,4),Om))+15*pad(pm(pw(A,2),pw(Om,2)))+pad(pw(Om,3))
    return 2*Phi+pad(pw(L[0],6))+pad(pw(L[1],6))+pad(pw(L[2],6))-pad(pw(L[3],6))
def jac(v,h=1e-7):
    r0=res(v); J=np.zeros((19,26),dtype=complex)
    for k in range(26):
        dv=np.zeros(26,dtype=complex); dv[k]=h; J[:,k]=(res(v+dv)-r0)/h
    return J
def norm(v):
    v=v.copy(); s=np.abs(np.r_[v[0:4],v[10:26]]).max(); v[0:4]/=s; v[10:26]/=s
    # R*sqrt(om) must scale like A: scale R by 1/s, keep om
    v[4:7]/=s
    return v
def newton(v,real,iters=120):
    for it in range(iters):
        v=norm(v); r=res(v); sc=1+abs(v).max()**6
        if np.linalg.norm(r)<1e-13*sc: return v,True
        v=v+np.linalg.lstsq(jac(v),-r,rcond=None)[0]
        if real: v=v.real.astype(complex)
    return v,np.linalg.norm(res(v))<1e-11*(1+abs(v).max()**6)
if __name__=="__main__":
    n=int(sys.argv[1]) if len(sys.argv)>1 else 200
    for real in (False,True):
        rng=np.random.default_rng(3); ok=0; nd=0; dims=[]; sols=[]
        for s in range(n):
            v=rng.standard_normal(26)+(0 if real else 1j*rng.standard_normal(26))
            v,good=newton(v.astype(complex),real)
            if good:
                ok+=1; A,R,om,L=unpack(v)
                deg=min([np.linalg.norm(A),np.linalg.norm(R),np.linalg.norm(om)]+[np.linalg.norm(l) for l in L])
                disc=om[1]**2-4*om[0]*om[2]
                if deg>1e-3 and abs(disc)>1e-6:
                    nd+=1; J=jac(v); sv=np.linalg.svd(J,compute_uv=False); dims.append(int(26-(sv>1e-8*sv[0]).sum())); sols.append(v)
        print("real" if real else "complex",": converged",ok," nondegenerate",nd," local dims",Counter(dims))
        np.save("/tmp/claude-0/-home-user/15becd99-5cb2-5637-beeb-15d2ba0dfed3/scratchpad/symsext_%s.npy"%("real" if real else "cplx"),np.array(sols))
