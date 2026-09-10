"""Genus-1 curves on X symmetric under x1<->x2, through a rational point of the Fano quotient Z = X/tau.
Z: P(e1,e2) + x3^6+x4^6+x5^6 = x6^6 with P(e1,e2) = e1^6-6e1^4e2+9e1^2e2^2-2e2^3 (= x1^6+x2^6, e1=x1+x2, e2=x1x2).
Weighted conic on Z through z0=(e1,e2,x3,x4,x5,x6) at t=0: e1,x3,x4,x5,x6 quadratic in t, e2 quartic.
Pull-back to X: x1,x2 = (e1 ± w)/2 with w^2 = e1^2 - 4 e2 (a quartic in t): genus 1.
Normalize reparametrizations fixing t=0 by (x6 linear coeff)=0 and (x3 linear coeff)=1.  12 unknowns, 12 eqs."""
import numpy as np, sys, json
from fractions import Fraction
def pm(a,b): return np.convolve(a,b)
def pw(a,k):
    r=np.array([1.0+0j])
    for _ in range(k): r=pm(r,a)
    return r
def pad(a,n=13):
    r=np.zeros(n,dtype=complex); r[:len(a)]=a; return r
def Pform(e1,e2):
    return pad(pw(e1,6))-6*pad(pm(pw(e1,4),e2))+9*pad(pm(pw(e1,2),pw(e2,2)))-2*pad(pw(e2,3))
class Sys:
    def __init__(s,z0):
        s.e1,s.e2,s.x3,s.x4,s.x5,s.x6=[complex(v) for v in z0]
        chk=Pform(np.array([s.e1]),np.array([s.e2]))[0]+s.x3**6+s.x4**6+s.x5**6-s.x6**6
        assert abs(chk)<1e-6*abs(s.x6)**6
    def forms(s,v):
        a=v
        E1=np.array([s.e1,a[0],a[1]]); X3=np.array([s.x3,1.0,a[2]]); X4=np.array([s.x4,a[3],a[4]])
        X5=np.array([s.x5,a[5],a[6]]); X6=np.array([s.x6,0.0,a[7]]); E2=np.array([s.e2,a[8],a[9],a[10],a[11]])
        return E1,E2,X3,X4,X5,X6
    def res(s,v):
        E1,E2,X3,X4,X5,X6=s.forms(v)
        r=Pform(E1,E2)+pad(pw(X3,6))+pad(pw(X4,6))+pad(pw(X5,6))-pad(pw(X6,6))
        return r[1:]
    def jac(s,v,h=1e-7):
        r0=s.res(v); J=np.zeros((12,12),dtype=complex)
        for i in range(12):
            dv=np.zeros(12,dtype=complex); dv[i]=h; J[:,i]=(s.res(v+dv)-r0)/h
        return J
    def scale(s,v): return 1+max(abs(s.x6),np.abs(v).max())**6
    def newton(s,v,iters=60):
        for it in range(iters):
            r=s.res(v)
            if np.linalg.norm(r)<1e-11*s.scale(v): return v,np.linalg.norm(r)/s.scale(v)
            try: dv=np.linalg.solve(s.jac(v),-r)
            except np.linalg.LinAlgError: return v,1.0
            lam=1.0
            while lam>1e-3 and np.linalg.norm(s.res(v+lam*dv))>np.linalg.norm(r): lam/=2
            v=v+lam*dv
            if np.abs(v).max()>1e4*max(1,abs(s.x6)): return v,1.0
        return v,np.linalg.norm(s.res(v))/s.scale(v)
    def solve_all(s,starts,seed=0,real=True):
        rng=np.random.default_rng(seed); sols=[]; base=max(abs(s.x6),1.0)
        for k in range(starts):
            v=base*10**rng.uniform(-1.5,1.5)*(rng.standard_normal(12)+(0 if real else 1j*rng.standard_normal(12)))
            v=v.astype(complex)
            v,nr=s.newton(v)
            if nr<1e-13:
                J=s.jac(v); sv=np.linalg.svd(J,compute_uv=False)
                if sv[-1]/sv[0]<1e-10: continue
                if not any(np.linalg.norm(v-u)<1e-6*(1+np.linalg.norm(u)) for u in sols): sols.append(v)
        return sols
def ratrec_float(z,maxden=3000):
    if abs(z.imag)>1e-7*(1+abs(z)): return None
    f=Fraction(z.real).limit_denominator(maxden)
    return f if abs(float(f)-z.real)<1e-10*(1+abs(z.real)) else None
if __name__=="__main__":
    pts=json.load(open(sys.argv[1])); starts=int(sys.argv[2]) if len(sys.argv)>2 else 800
    for p in pts:
        A,Om=p['A'],p['Om']; r3,r4,r5=p['r']; r6=p['r6']
        z0=(2*A, A*A-Om, r3,r4,r5,r6)
        S=Sys(z0); sols=S.solve_all(starts)
        rat=[v for v in sols if all(ratrec_float(z) is not None for z in v)]
        print("Z-point",z0,": distinct genus-1 curves",len(sols)," near-real",sum(1 for v in sols if np.abs(v.imag).max()<1e-6)," rational-looking",len(rat),flush=True)
        for v in rat: print("   CANDIDATE (verify exactly):",[str(ratrec_float(z)) for z in v],flush=True)
