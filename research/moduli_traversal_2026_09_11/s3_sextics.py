"""S3-symmetric sextic rational curves on X.  M:(s,t)->(t,t-s) (order 3), iota:(s,t)->(t,s); iota M iota = M^-1.
x2 = x1∘M, x3 = x1∘M^2, x1 iota-(anti)invariant so iota swaps x2,x3; x4,x5,x6 in an S3-character space.
Only gauge: overall scaling => rational normal form (monic in a rational basis)."""
import numpy as np, sys, itertools
from collections import Counter
from math import comb
d=6
def actM(c):
    out=np.zeros(d+1,dtype=complex)
    for k,ck in enumerate(c):
        for j in range(k+1): out[d-k+j]+=ck*comb(k,j)*(-1)**(k-j)
    return out
def actI(c): return c[::-1].copy()
I7=np.eye(d+1)
AM=np.array([actM(I7[i]) for i in range(d+1)]).T
AI=np.array([actI(I7[i]) for i in range(d+1)]).T
def kernel(A,tol=1e-9):
    u,s,vt=np.linalg.svd(A); r=(s>tol).sum(); return vt[r:].T
def space(chiI,chiM=1):
    K=np.vstack([AM-chiM*np.eye(d+1),AI-chiI*np.eye(d+1)]) if chiM is not None else (AI-chiI*np.eye(d+1))
    return kernel(K)
B1={+1:space(+1,None),-1:space(-1,None)}      # x1: iota-(anti)invariant
B4={+1:space(+1,1),-1:space(-1,1)}            # S3-invariant / sign
def pw(a,k):
    r=np.array([1.0+0j])
    for _ in range(k): r=np.convolve(r,a)
    return r
def pad(a,n=6*d+1):
    r=np.zeros(n,dtype=complex); r[:len(a)]=a; return r
def run(c1,c4,c5,c6,starts,real):
    b1=B1[c1]; b4=B4[c4]; b5=B4[c5]; b6=B4[c6]
    dims=[b1.shape[1],b4.shape[1],b5.shape[1],b6.shape[1]]; n=sum(dims)
    def forms(v):
        i=0; x1=b1@v[i:i+dims[0]]; i+=dims[0]; x4=b4@v[i:i+dims[1]]; i+=dims[1]; x5=b5@v[i:i+dims[2]]; i+=dims[2]; x6=b6@v[i:i+dims[3]]
        return [x1,actM(x1),actM(actM(x1)),x4,x5,x6]
    def res(v):
        X=forms(v); return sum(pad(pw(X[i],6)) for i in range(5))-pad(pw(X[5],6))
    def jac(v,h=1e-7):
        r0=res(v); J=np.zeros((37,n),dtype=complex)
        for k in range(n):
            dv=np.zeros(n,dtype=complex); dv[k]=h; J[:,k]=(res(v+dv)-r0)/h
        return J
    rng=np.random.default_rng(1); found=[]; ld=[]
    for s in range(starts):
        v=(rng.standard_normal(n)+(0 if real else 1j*rng.standard_normal(n))).astype(complex)
        for it in range(150):
            v=v/np.abs(v).max(); r=res(v)
            if np.linalg.norm(r)<1e-13: break
            v=v+np.linalg.lstsq(jac(v),-r,rcond=None)[0]
            if real: v=v.real.astype(complex)
        if np.linalg.norm(res(v))<1e-11:
            X=forms(v)
            if min(np.linalg.norm(x) for x in X)<1e-3: continue
            P6=[pad(pw(X[i],6))*(1 if i<5 else -1) for i in range(6)]; sc=max(np.linalg.norm(p) for p in P6)
            if any(np.linalg.norm(sum(P6[i] for i in S))<1e-7*sc for k in range(1,6) for S in itertools.combinations(range(6),k)): continue
            J=jac(v); sv=np.linalg.svd(J,compute_uv=False); ld.append(int(n-(sv>1e-8*sv[0]).sum())); found.append(v)
    return dims,found,ld
if __name__=="__main__":
    starts=int(sys.argv[1]) if len(sys.argv)>1 else 60
    for c1 in (1,-1):
        for c4,c5,c6 in itertools.product((1,-1),repeat=3):
            for real in (False,True):
                dims,found,ld=run(c1,c4,c5,c6,starts,real)
                if found or real: print("chars x1=%+d x4..6=%+d%+d%+d dims %s %s: nondegenerate %d local dims %s"%(c1,c4,c5,c6,dims,"real" if real else "cplx",len(found),Counter(ld)),flush=True)
