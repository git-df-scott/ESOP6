"""BCU-type ansatz with five sixth powers:  G1^6+G2^6+G3^6+G4^6 + F^3 = G6^6  in C[t],
deg G_i <= 2, deg F <= 4.  A rational point + a rational point on w^2=F(t) gives ESOP6 points
(x5 = w).  Numerical existence / dimension / reality test."""
import numpy as np, sys
from collections import Counter
def pm(a,b): return np.convolve(a,b)
def pw(a,k):
    r=np.array([1.0+0j])
    for _ in range(k): r=pm(r,a)
    return r
def pad(a,n=13):
    r=np.zeros(n,dtype=complex); r[:len(a)]=a; return r
def unpack(v):
    G=[v[3*i:3*i+3] for i in range(5)]   # G1..G4, G6
    F=v[15:20]
    return G,F
def res(v):
    G,F=unpack(v)
    return sum(pad(pw(G[i],6)) for i in range(4))+pad(pw(F,3))-pad(pw(G[4],6))
def jac(v,h=1e-7):
    r0=res(v); J=np.zeros((13,20),dtype=complex)
    for k in range(20):
        dv=np.zeros(20,dtype=complex); dv[k]=h; J[:,k]=(res(v+dv)-r0)/h
    return J
def normalize(v):
    G,F=unpack(v)
    s=max(np.abs(v[:15]).max(),1e-300); v=v.copy(); v[:15]/=s; v[15:]/=s*s
    return v
def newton(v,iters=100):
    for it in range(iters):
        v=normalize(v); r=res(v)
        if np.linalg.norm(r)<1e-13: return v,np.linalg.norm(r)
        v=v+np.linalg.lstsq(jac(v),-r,rcond=None)[0]
    v=normalize(v); return v,np.linalg.norm(res(v))
real='--real' in sys.argv
rng=np.random.default_rng(11)
sols=[]
for s in range(300):
    v=rng.standard_normal(20)+(0 if real else 1j*rng.standard_normal(20))
    v,nr=newton(v)
    if nr<1e-11:
        G,F=unpack(v)
        degen=min([np.linalg.norm(g) for g in G]+[np.linalg.norm(F)])
        # F a perfect square? (then it's just a conic) : check discriminant-like: roots of F in pairs
        rts=np.roots(F[::-1]) if abs(F[4])>1e-8 else np.roots(F[:4][::-1])
        Fsq = len(rts)>=2 and min(abs(rts[i]-rts[j]) for i in range(len(rts)) for j in range(i+1,len(rts)))<1e-5
        J=jac(v); sv=np.linalg.svd(J,compute_uv=False); dim=20-(sv>1e-7*sv[0]).sum()
        sols.append((dim,degen>1e-3,Fsq,v))
print("converged",len(sols)," nondegenerate:",sum(1 for d,g,q,v in sols if g),
      " of which F perfect square (plain conic):",sum(1 for d,g,q,v in sols if g and q))
print("local dims (nondeg, F not square):",Counter(d for d,g,q,v in sols if g and not q))
np.save("/tmp/claude-0/-home-user/15becd99-5cb2-5637-beeb-15d2ba0dfed3/scratchpad/v3_%s.npy"%("real" if real else "cplx"),np.array([v for d,g,q,v in sols if g and not q]))
for d,g,q,v in sols[:3]:
    if g and not q:
        G,F=unpack(v); print("dim",d); print(" G:",np.round(np.array(G),3)); print(" F:",np.round(F,3))
