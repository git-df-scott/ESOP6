"""Genus-1 cubic-cover curves on X through a rational point of Y2.
Y2: x1^6+x2^6+x3^6+x4^6 + S^2 = x6^6  (S = x5^3).  Ansatz ('Y2-conic with a triple root'):
   G_i(t) = x_i + g_i1 t + g_i2 t^2   (i=1,2,3,4,6),   S(t) = K(t) M(t)^3,  K cubic, M linear,
   sum G_i^6 + K^2 M^6 = G6^6.   Curve on X:  x5^3 = K(t) M(t)^3, i.e. u^3 = K(t) (u = x5/M): genus 1 if K has
   distinct roots; K irreducible over Q avoids rational (6,1,4) points at the zeros of x5.
Normalizations: M(0)=1, K(0)=s0, g11=1, g61=0.  Unknowns: g12,g21,g22,g31,g32,g41,g42,g62, k1,k2,k3, m1  (12);
equations: t^1..t^12 coefficients (12).  Zero-dimensional.  Complex Newton from random starts, dedupe, polish in
mpmath, rational reconstruction, exact verification."""
import numpy as np, sys, json, mpmath as mp
from fractions import Fraction
def pm(a,b): return np.convolve(a,b)
def pw(a,k):
    r=np.array([1.0+0j])
    for _ in range(k): r=pm(r,a)
    return r
def pad(a,n=13):
    r=np.zeros(n,dtype=complex); r[:len(a)]=a; return r
class Sys:
    def __init__(s,x,x6,S):
        s.norm=float(x6)
        s.x=[complex(v)/s.norm for v in x]; s.x6=complex(1.0); s.s0=complex(S)/s.norm**3
        assert abs(sum(v**6 for v in s.x)+s.s0**2-s.x6**6)<1e-9*abs(s.x6)**6
    def forms(s,v):
        g=v[:8]; k=v[8:11]; m1=v[11]
        G=[np.array([s.x[0],1.0,g[0]]),np.array([s.x[1],g[1],g[2]]),np.array([s.x[2],g[3],g[4]]),np.array([s.x[3],g[5],g[6]])]
        G6=np.array([s.x6,0.0,g[7]]); K=np.array([s.s0,k[0],k[1],k[2]]); M=np.array([1.0,m1])
        return G,G6,K,M
    def res(s,v):
        G,G6,K,M=s.forms(v)
        r=sum(pad(pw(Gi,6)) for Gi in G)+pad(pm(pw(K,2),pw(M,6)))-pad(pw(G6,6))
        return r[1:]
    def jac(s,v,h=1e-7):
        r0=s.res(v); J=np.zeros((12,12),dtype=complex)
        for i in range(12):
            dv=np.zeros(12,dtype=complex); dv[i]=h; J[:,i]=(s.res(v+dv)-r0)/h
        return J
    def scale(s,v): return 1+max(abs(s.x6),abs(s.s0)**(1/3),np.abs(v).max())**6
    def newton(s,v,iters=120):
        for it in range(iters):
            r=s.res(v)
            if np.linalg.norm(r)<1e-11*s.scale(v): return v,np.linalg.norm(r)/s.scale(v)
            try: dv=np.linalg.solve(s.jac(v),-r)
            except np.linalg.LinAlgError: return v,1.0
            lam=1.0
            while lam>1e-3 and np.linalg.norm(s.res(v+lam*dv))>np.linalg.norm(r): lam/=2
            v=v+lam*dv
            if np.abs(v).max()>1e7: return v,1.0
        return v,np.linalg.norm(s.res(v))/s.scale(v)
    def solve_all(s,starts,seed=0,real=False):
        rng=np.random.default_rng(seed); sols=[]
        base=1.0
        for k in range(starts):
            v=base*10**rng.uniform(-1.5,1.5)*(rng.standard_normal(12)+(0 if real else 1j*rng.standard_normal(12)))
            v=v.astype(complex)
            v,nr=s.newton(v)
            if nr<1e-13:
                J=s.jac(v); sv=np.linalg.svd(J,compute_uv=False)
                if sv[-1]/sv[0]<1e-10: continue
                if not any(np.linalg.norm(v-u)<1e-6*(1+np.linalg.norm(u)) for u in sols): sols.append(v)
        return sols
def ratrec_float(z,maxden=10**6):
    if abs(z.imag)>1e-7*(1+abs(z)): return None
    f=Fraction(z.real).limit_denominator(maxden)
    return f if abs(float(f)-z.real)<1e-10*(1+abs(z.real)) else None
def exact_check(x,x6,S,fr):
    import sympy as sp
    t=sp.symbols('t')
    g=fr[:8]; k=fr[8:11]; m1=fr[11]
    G=[x[0]+t+g[0]*t**2, x[1]+g[1]*t+g[2]*t**2, x[2]+g[3]*t+g[4]*t**2, x[3]+g[5]*t+g[6]*t**2]
    G6=x6+g[7]*t**2; K=S+k[0]*t+k[1]*t**2+k[2]*t**3; M=1+m1*t
    expr=sp.expand(sum(Gi**6 for Gi in G)+K**2*M**6-G6**6)
    return expr==0, [str(Gi) for Gi in G], str(G6), str(K), str(M)
if __name__=="__main__":
    pts=json.load(open(sys.argv[1])); starts=int(sys.argv[2]) if len(sys.argv)>2 else 600; tag=sys.argv[3] if len(sys.argv)>3 and not sys.argv[3].startswith('--') else 'x'
    seen=set(); out=[]
    for p in pts:
        if sum(int(v)**6 for v in p['x'])+int(p['S'])**2 != int(p['x6'])**6:
            print('INVALID EXACT Y2 SEED; skipped',p,flush=True)
            continue
        key=(tuple(sorted(p['x'])),p['x6'])
        if key in seen: continue
        seen.add(key)
        S=Sys(p['x'],p['x6'],p['S'])
        sols=S.solve_all(starts,real=('--real' in sys.argv))
        rat=0
        for v in sols:
            fr=[ratrec_float(z) for z in v]
            if all(f is not None for f in fr):
                xn=[Fraction(v,p['x6']) for v in p['x']]; ok,G,G6,K,M=exact_check(xn,Fraction(1),Fraction(p['S'],p['x6']**3),fr)
                if ok:
                    rat+=1
                    print("!!!! RATIONAL GENUS-1 CURVE through",p, "G=",G,"G6=",G6,"K=",K,"M=",M,flush=True)
                    out.append(dict(point=p,G=G,G6=G6,K=K,M=M))
        nearreal=sum(1 for v in sols if np.abs(v.imag).max()<1e-6)
        print("seed",key,"S=",p['S'],": distinct curves found",len(sols)," real",nearreal," rational",rat,flush=True)
    json.dump(out,open("research/cube_ansatz_2026_09_10/cubic_cover_rational_hits_%s.json"%tag,"w"),indent=1)
