"""Genus-1 tau-symmetric quartic curves on Z = X/tau through a rational point of Z.
e1,x3..x6 quartic in (E1,E2); e2 = (e1^2 - R^2 Q4)/4, R quadratic, Q4 quartic; x1,x2 = (e1 ± R w)/2, w^2 = Q4.
Point z0 at E=[1:0]. Normalizations: R(1,0)=1, Q4(1,0)=e1^2-4e2 at z0, x6's E1^3E2 coeff = 0, x5's E1^3E2 coeff = 1.
24 unknowns, 24 equations."""
import numpy as np, sys, json
from fractions import Fraction
def pm(a,b): return np.convolve(a,b)
def pw(a,k):
    r=np.array([1.0+0j])
    for _ in range(k): r=pm(r,a)
    return r
def pad(a,n=25):
    r=np.zeros(n,dtype=complex); r[:len(a)]=a; return r
def Pform(e1,e2):
    return pad(pw(e1,6))-6*pad(pm(pw(e1,4),e2))+9*pad(pm(pw(e1,2),pw(e2,2)))-2*pad(pw(e2,3))
class Sys:
    def __init__(s,z0):
        s.e1,s.e2,s.x3,s.x4,s.x5,s.x6=[complex(v) for v in z0]; s.q0=s.e1**2-4*s.e2
    def forms(s,v):
        a=v
        E1=np.array([s.e1,a[0],a[1],a[2],a[3]]); X3=np.array([s.x3,a[4],a[5],a[6],a[7]]); X4=np.array([s.x4,a[8],a[9],a[10],a[11]])
        X5=np.array([s.x5,1.0,a[12],a[13],a[14]]); X6=np.array([s.x6,0.0,a[15],a[16],a[17]])
        R=np.array([1.0,a[18],a[19]]); Q=np.array([s.q0,a[20],a[21],a[22],a[23]])
        E2=(pm(E1,E1)-pm(pm(R,R),Q))/4
        return E1,E2,X3,X4,X5,X6,R,Q
    def res(s,v):
        E1,E2,X3,X4,X5,X6,R,Q=s.forms(v)
        return (Pform(E1,E2)+pad(pw(X3,6))+pad(pw(X4,6))+pad(pw(X5,6))-pad(pw(X6,6)))[1:]
    def jac(s,v,h=1e-7):
        r0=s.res(v); J=np.zeros((24,24),dtype=complex)
        for i in range(24):
            dv=np.zeros(24,dtype=complex); dv[i]=h; J[:,i]=(s.res(v+dv)-r0)/h
        return J
    def scale(s,v): return 1+max(abs(s.x6),np.abs(v).max())**6
    def newton(s,v,iters=80):
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
    def solve_all(s,starts,seed=0):
        rng=np.random.default_rng(seed); sols=[]; base=max(abs(s.x6),1.0)
        for k in range(starts):
            v=(base*10**rng.uniform(-1.5,1.5)*rng.standard_normal(24)).astype(complex)
            v,nr=s.newton(v)
            if nr<1e-13:
                J=s.jac(v); sv=np.linalg.svd(J,compute_uv=False)
                if sv[-1]/sv[0]<1e-10: continue
                if not any(np.linalg.norm(v-u)<1e-6*(1+np.linalg.norm(u)) for u in sols): sols.append(v)
        return sols
def ratrec(z,maxden=3000):
    if abs(z.imag)>1e-7*(1+abs(z)): return None
    f=Fraction(z.real).limit_denominator(maxden)
    return f if abs(float(f)-z.real)<1e-10*(1+abs(z.real)) else None
if __name__=="__main__":
    pts=json.load(open(sys.argv[1])); starts=int(sys.argv[2]) if len(sys.argv)>2 else 1000
    for p in pts:
        A,Om=p['A'],p['Om']; r3,r4,r5=p['r']; r6=p['r6']
        z0=(2*A,A*A-Om,r3,r4,r5,r6); S=Sys(z0); sols=S.solve_all(starts)
        rat=[v for v in sols if all(ratrec(z) is not None for z in v)]
        print("Z-point",z0,": genus-1 quartic curves found",len(sols),"real",sum(1 for v in sols if np.abs(v.imag).max()<1e-8),"rational-looking",len(rat),flush=True)
        for v in rat: print("   !!!! CANDIDATE (verify exactly):",[str(ratrec(z)) for z in v],flush=True)
