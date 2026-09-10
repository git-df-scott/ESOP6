"""Conics on the Fano fourfold  Y: x1^6+x2^6+x3^6+x4^6 + T^3 = x6^6  (weights 1,1,1,1,1,2; X -> Y is T = x5^2)
through a given point y0=(x1,x2,x3,x4,x6,T).  Parametrize  G_i(t) = x_i + g_i1 t + g_i2 t^2 (i=1,2,3,4,6),
F(t) = T + f1 t + f2 t^2 + f3 t^3 + f4 t^4.  Normalize the reparametrizations fixing t=0 by g11 = 1, g61 = 0.
12 unknowns, 12 equations (t^1..t^12 coefficients):  zero-dimensional.  A rational solution is a genus-1 curve
w^2 = F(t) over Q on X; rational points on it with w != 0 are ESOP6 solutions."""
import numpy as np, sys
def pm(a,b): return np.convolve(a,b)
def pw(a,k):
    r=np.array([1.0+0j])
    for _ in range(k): r=pm(r,a)
    return r
def pad(a,n=13):
    r=np.zeros(n,dtype=complex); r[:len(a)]=a; return r
class System:
    def __init__(s,y0):
        s.x=np.array(y0[:5],dtype=complex); s.T=complex(y0[5])
        assert abs(sum(s.x[:4]**6)+s.T**3-s.x[4]**6)<1e-9*max(1,abs(s.x).max()**6)
    def forms(s,v):
        g=v[:8]; f=v[8:12]
        G=[np.array([s.x[0],1.0,g[0]]), np.array([s.x[1],g[1],g[2]]), np.array([s.x[2],g[3],g[4]]),
           np.array([s.x[3],g[5],g[6]]), np.array([s.x[4],0.0,g[7]])]
        F=np.array([s.T,f[0],f[1],f[2],f[3]])
        return G,F
    def res(s,v):
        G,F=s.forms(v)
        r=sum(pad(pw(G[i],6)) for i in range(4))+pad(pw(F,3))-pad(pw(G[4],6))
        return r[1:]
    def jac(s,v,h=1e-7):
        r0=s.res(v); J=np.zeros((12,12),dtype=complex)
        for k in range(12):
            dv=np.zeros(12,dtype=complex); dv[k]=h; J[:,k]=(s.res(v+dv)-r0)/h
        return J
    def newton(s,v,iters=150):
        for it in range(iters):
            r=s.res(v); sc=1+np.abs(v).max()**6+abs(s.T)**3
            if np.linalg.norm(r)/sc<1e-14: return v,np.linalg.norm(r)/sc
            try: dv=np.linalg.solve(s.jac(v),-r)
            except np.linalg.LinAlgError: return v,1.0
            lam=1.0
            while lam>1e-3 and np.linalg.norm(s.res(v+lam*dv))>np.linalg.norm(r): lam/=2
            v=v+lam*dv
            if np.abs(v).max()>1e8: return v,1.0
        return v,np.linalg.norm(s.res(v))/(1+np.abs(v).max()**6)
    def solve_all(s,starts=2000,seed=0,scale_range=(-1,2)):
        rng=np.random.default_rng(seed); sols=[]
        for k in range(starts):
            v=10**rng.uniform(*scale_range)*(rng.standard_normal(12)+1j*rng.standard_normal(12))
            v,nr=s.newton(v)
            if nr<1e-12:
                J=s.jac(v); sv=np.linalg.svd(J,compute_uv=False)
                if sv[-1]/sv[0]<1e-9: continue   # skip singular/degenerate
                if not any(np.linalg.norm(v-u)<1e-7*(1+np.linalg.norm(u)) for u in sols): sols.append(v)
        return sols
if __name__=="__main__":
    rng=np.random.default_rng(5)
    x=rng.standard_normal(5)+1j*rng.standard_normal(5)
    T=(x[4]**6-sum(x[:4]**6))**(1/3)
    S=System(list(x)+[T])
    sols=S.solve_all(starts=int(sys.argv[1]) if len(sys.argv)>1 else 1500)
    print("random complex point: distinct nonsingular conics through it:",len(sols))
    # also count how many are 'degenerate' (F perfect square => plain conic on X)
    sq=0
    for v in sols:
        G,F=S.forms(v); rts=np.roots(F[::-1])
        if min(abs(rts[i]-rts[j]) for i in range(4) for j in range(i+1,4))<1e-5: sq+=1
    print("of which F has a repeated root:",sq)
