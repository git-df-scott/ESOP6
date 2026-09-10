"""Z/3-symmetric sextic rational curves on X.
sigma: (s,t) -> (t, t-s) has order 3 in PGL2 (M^3 = -I).  Ansatz: x2 = x1∘M, x3 = x1∘M^2, x4,x5,x6 M-invariant sextics.
Count: x1 (7) + 3*(dim of invariant sextics) unknowns, ~13 independent equations, gauge = centralizer torus (1) + scale (1)."""
import numpy as np, sys
from collections import Counter
from numpy.polynomial import polynomial as P
d=6
# represent a binary form by coefficient vector c[k] of s^(d-k) t^k ; act by M: F(s,t) -> F(t, t-s)
def act(c):
    # F(t, t-s): substitute s->t, t->t-s ; compute via polynomial in (s,t): expand sum c_k t^(d-k) (t-s)^k
    out=np.zeros(d+1,dtype=complex)  # index = power of t ... we keep form as coefficients of s^(d-j) t^j
    for k,ck in enumerate(c):
        # t^(d-k) * (t-s)^k = t^(d-k) * sum_j C(k,j) t^j (-s)^(k-j)
        for j in range(k+1):
            from math import comb
            out[d-k+j]+= ck*comb(k,j)*(-1)**(k-j)   # s^(k-j) t^(d-k+j): power of t is d-k+j
    return out
# invariant subspace: kernel of (act - I)
A=np.array([act(np.eye(d+1)[i]) for i in range(d+1)]).T   # columns = act(e_i)
w,V=np.linalg.eig(A)
inv=V[:,np.abs(w-1)<1e-8].real
inv=np.linalg.qr(inv)[0]   # orthonormal basis of invariant sextics
m=inv.shape[1]; print("dim invariant sextics:",m)
def pw(a,k):
    r=np.array([1.0+0j])
    for _ in range(k): r=np.convolve(r,a)
    return r
def pad(a,n=6*d+1):
    r=np.zeros(n,dtype=complex); r[:len(a)]=a; return r
def forms(v):
    x1=v[:7]; x2=act(x1); x3=act(x2)
    x4=inv@v[7:7+m]; x5=inv@v[7+m:7+2*m]; x6=inv@v[7+2*m:7+3*m]
    return [x1,x2,x3,x4,x5,x6]
def res(v):
    X=forms(v)
    return sum(pad(pw(X[i],6)) for i in range(5))-pad(pw(X[5],6))
n=7+3*m
def jac(v,h=1e-7):
    r0=res(v); J=np.zeros((37,n),dtype=complex)
    for k in range(n):
        dv=np.zeros(n,dtype=complex); dv[k]=h; J[:,k]=(res(v+dv)-r0)/h
    return J
def newton(v,real):
    for it in range(150):
        v=v/np.abs(v).max(); r=res(v)
        if np.linalg.norm(r)<1e-13: return v,True
        v=v+np.linalg.lstsq(jac(v),-r,rcond=None)[0]
        if real: v=v.real.astype(complex)
    return v,np.linalg.norm(res(v))<1e-11
for real in (False,True):
    rng=np.random.default_rng(9); ok=0; nd=0; dims=[]; sols=[]
    for s in range(int(sys.argv[1]) if len(sys.argv)>1 else 200):
        v=rng.standard_normal(n)+(0 if real else 1j*rng.standard_normal(n))
        v,good=newton(v.astype(complex),real)
        if good:
            ok+=1; X=forms(v)
            if min(np.linalg.norm(x) for x in X)>1e-3:
                # exclude vanishing subsums among the six
                import itertools
                P6=[pad(pw(X[i],6))*(1 if i<5 else -1) for i in range(6)]; sc=max(np.linalg.norm(p) for p in P6); degen=False
                for k in range(1,6):
                    for S in itertools.combinations(range(6),k):
                        if np.linalg.norm(sum(P6[i] for i in S))<1e-7*sc: degen=True
                if not degen:
                    nd+=1; J=jac(v); sv=np.linalg.svd(J,compute_uv=False); dims.append(int(n-(sv>1e-8*sv[0]).sum())); sols.append(v)
    print("real" if real else "complex",": converged",ok," nondegenerate",nd," local dims",Counter(dims))
    np.save("/tmp/claude-0/-home-user/15becd99-5cb2-5637-beeb-15d2ba0dfed3/scratchpad/z3sext_%s.npy"%("real" if real else "cplx"),np.array(sols))
