"""Weighted conics on Y2 through two known Y2-points y0 (t=0) and y1 (t=inf):
  G_i(t) = x_i + g_i t + mu*y_i t^2  (i=1,2,3,4,6),  S(t) sextic with S(0)=s0, lead = mu^3 s1.
Condition: D(t) := G6^6 - sum G_i^6 (degree 12, coefficients polynomial in g_1..g_4,g_6,mu) is a perfect square.
Square test: take the polynomial square root by power series from the constant term s0^2 (nonzero), the defect is the
last 6 coefficients -> 6 equations in 6 unknowns. Newton from random real starts; rational recognition; exact check."""
import numpy as np, json, sys, itertools, math, time
from fractions import Fraction
def pw(a,k):
    r=np.array([1.0+0j])
    for _ in range(k): r=np.convolve(r,a)
    return r
def pad(a,n=13):
    r=np.zeros(n,dtype=complex); r[:len(a)]=a; return r
def sqrt_defect(D,s0):
    # find s with s^2 = D up to degree 6 (7 coeffs), starting from s[0]=s0; return remaining 6 residuals
    s=np.zeros(7,dtype=complex); s[0]=s0
    for k in range(1,7):
        acc=sum(s[j]*s[k-j] for j in range(1,k))
        s[k]=(D[k]-acc)/(2*s0)
    sq=pw(s,2)
    return (D-pad(sq))[7:13], s
class Pair:
    def __init__(s,y0,y1):
        s.x=np.array(y0['x']+[y0['x6']],dtype=float); s.y=np.array(y1['x']+[y1['x6']],dtype=float)
        s.s0=float(y0['S']); s.s1=float(y1['S'])
    def D(s,v):
        g=v[:5]; mu=v[5]
        G=[np.array([s.x[i],g[i],mu*s.y[i]]) for i in range(5)]
        return pad(pw(G[4],6))-sum(pad(pw(G[i],6)) for i in range(4))
    def res(s,v):
        d,_=sqrt_defect(s.D(v),s.s0); return d
    def jac(s,v,h=1e-7):
        r0=s.res(v); J=np.zeros((6,6),dtype=complex)
        for i in range(6):
            dv=np.zeros(6,dtype=complex); dv[i]=h; J[:,i]=(s.res(v+dv)-r0)/h
        return J
    def solve(s,starts,rng):
        sols=[]; base=max(abs(s.x).max(),1.0)
        for k in range(starts):
            v=(base*10**rng.uniform(-1,1)*rng.standard_normal(6)).astype(complex)
            for it in range(60):
                r=s.res(v); sc=1+np.abs(s.D(v)).max()
                if np.linalg.norm(r)<1e-12*sc: break
                try: dv=np.linalg.solve(s.jac(v),-r)
                except np.linalg.LinAlgError: break
                lam=1.0
                while lam>1e-3 and np.linalg.norm(s.res(v+lam*dv))>np.linalg.norm(r): lam/=2
                v=v+lam*dv
                if np.abs(v).max()>1e6*base: break
            r=s.res(v); sc=1+np.abs(s.D(v)).max()
            if np.linalg.norm(r)<1e-11*sc and abs(v[5])>1e-9:
                # also the leading coefficient must match: S(inf) = ±mu^3 s1  (automatic if D's top coeff = mu^6 s1^2)
                if not any(np.linalg.norm(v-u)<1e-6*(1+np.linalg.norm(u)) for u in sols): sols.append(v)
        return sols
def ratrec(z,maxden=5000):
    if abs(z.imag)>1e-8*(1+abs(z)): return None
    f=Fraction(z.real).limit_denominator(maxden)
    return f if abs(float(f)-z.real)<1e-9*(1+abs(z.real)) else None
if __name__=="__main__":
    pts=json.load(open(sys.argv[1])); H=int(sys.argv[2]) if len(sys.argv)>2 else 150; starts=int(sys.argv[3]) if len(sys.argv)>3 else 6
    pts=[p for p in pts if p['x6']<=H]
    print("points",len(pts),"pairs",len(pts)*(len(pts)-1)//2,flush=True)
    rng=np.random.default_rng(0); t0=time.time(); n=0; found=0
    for a,b in itertools.combinations(range(len(pts)),2):
        for (p,q) in ((pts[a],pts[b]),(pts[b],pts[a])):
            # try all 24 assignments of q's four coordinates relative to p's (order matters for a conic)
            for perm in itertools.permutations(range(4)):
                q2=dict(x=[q['x'][i] for i in perm],x6=q['x6'],S=q['S'])
                P=Pair(p,q2); sols=P.solve(starts,rng); n+=1
                for v in sols:
                    rec=[ratrec(z) for z in v]
                    if all(r is not None for r in rec):
                        found+=1; print("!!!! RATIONAL Y2-CONIC CANDIDATE through",p,q2,[str(r) for r in rec],flush=True)
        if (a*len(pts)+b)%2000==0: print("progress pair",a,b,"instances",n,"elapsed %.0fs"%(time.time()-t0),flush=True)
    print("done instances",n,"candidates",found)
