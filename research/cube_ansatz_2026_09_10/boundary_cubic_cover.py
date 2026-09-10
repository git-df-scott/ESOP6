"""Boundary-contact genus-1 cubic covers on X (the only real members of family V3').

Odd-factor principle: N = G6^3 - K M^3 has odd degree 9, so it has a real root t1; there sum G_i^6 = N*(G6^3+KM^3) = 0
forces G_1..G_4 to vanish at t1, i.e. the curve passes through B = [0:0:0:0:1:1]. For a Q-curve t1 must be rational
(quadratic t1 would force all G_i proportional).  Put t1 = 0:
   G_i = t (a_i + b_i t)  (i=1..4),  G6 = 1 + c1 t + c2 t^2,  K = 1 + k1 t + k2 t^2 + k3 t^3,  M = 1 + m1 t,
   sum G_i^6 + K^2 M^6 = G6^6,  curve on X: x5 = u*M with u^3 = K  (genus 1 if K has distinct roots).
Normalizations: G6(0)=1 (scale), M(0)=1 and K(0)=1 ((K,M)-scaling, sign of K absorbed in u), a_1 = 1 (t-scaling),
b_1 = 0 (shear fixing 0).  Unknowns a2,a3,a4,b2,b3,b4,c1,c2,k1,k2,k3,m1 (12); equations t^1..t^12 (12): zero-dim."""
import numpy as np, sys, json
from fractions import Fraction
def pm(a,b): return np.convolve(a,b)
def pw(a,k):
    r=np.array([1.0+0j])
    for _ in range(k): r=pm(r,a)
    return r
def pad(a,n=13):
    r=np.zeros(n,dtype=complex); r[:len(a)]=a; return r
def forms(v):
    a2,a3,a4,b2,b3,b4,c1,c2,k1,k2,k3,m1=v
    G=[np.array([0,1.0,0]),np.array([0,a2,b2]),np.array([0,a3,b3]),np.array([0,a4,b4])]
    G6=np.array([1.0,c1,c2]); K=np.array([1.0,k1,k2,k3]); M=np.array([1.0,m1])
    return G,G6,K,M
def res(v):
    G,G6,K,M=forms(v)
    r=sum(pad(pw(g,6)) for g in G)+pad(pm(pw(K,2),pw(M,6)))-pad(pw(G6,6))
    return r[1:]
def jac(v,h=1e-7):
    r0=res(v); J=np.zeros((12,12),dtype=complex)
    for i in range(12):
        dv=np.zeros(12,dtype=complex); dv[i]=h; J[:,i]=(res(v+dv)-r0)/h
    return J
def newton(v,iters=80):
    for it in range(iters):
        r=res(v); sc=1+np.abs(v).max()**6
        if np.linalg.norm(r)<1e-12*sc: return v,np.linalg.norm(r)/sc
        try: dv=np.linalg.solve(jac(v),-r)
        except np.linalg.LinAlgError: return v,1.0
        lam=1.0
        while lam>1e-3 and np.linalg.norm(res(v+lam*dv))>np.linalg.norm(r): lam/=2
        v=v+lam*dv
        if np.abs(v).max()>1e5: return v,1.0
    return v,np.linalg.norm(res(v))/(1+np.abs(v).max()**6)
def ratrec(z,maxden=10**6):
    if abs(z.imag)>1e-7*(1+abs(z)): return None
    f=Fraction(z.real).limit_denominator(maxden)
    return f if abs(float(f)-z.real)<1e-8*(1+abs(z.real)) else None
if __name__=="__main__":
    starts=int(sys.argv[1]) if len(sys.argv)>1 else 2000; real='--real' in sys.argv
    rng=np.random.default_rng(0); sols=[]
    for k in range(starts):
        v=10**rng.uniform(-1,1.5)*(rng.standard_normal(12)+(0 if real else 1j*rng.standard_normal(12)))
        v,nr=newton(v.astype(complex))
        if nr<1e-13:
            J=jac(v); sv=np.linalg.svd(J,compute_uv=False)
            if sv[-1]/sv[0]<1e-10: continue
            G,G6,K,M=forms(v)
            if min(np.linalg.norm(g) for g in G)<1e-4 or np.linalg.norm(K[1:])<1e-6: continue   # degenerate
            if not any(np.linalg.norm(v-u)<1e-6*(1+np.linalg.norm(u)) for u in sols): sols.append(v)
    print("distinct isolated nondegenerate solutions:",len(sols)," real:",sum(1 for v in sols if np.abs(v.imag).max()<1e-7))
    for v in sols:
        fr=[ratrec(z) for z in v]
        if all(f is not None for f in fr): print("RATIONAL CANDIDATE:",[str(f) for f in fr])
    for v in sorted(sols,key=lambda v:np.abs(v.imag).max())[:6]:
        print(np.round(v,5))
    np.save("/tmp/claude-0/-home-user/15becd99-5cb2-5637-beeb-15d2ba0dfed3/scratchpad/bcc_%s.npy"%("real" if real else "cplx"),np.array(sols))
