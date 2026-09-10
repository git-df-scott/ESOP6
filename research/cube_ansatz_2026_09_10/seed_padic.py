"""Rational solutions of the zero-dimensional seed system via mod-p enumeration + Hensel + rational reconstruction.
System: Family-I conic through P=(9±sqrt(-249),14,18,0,22); unknowns (c,w1,w2,m3,m4,m6)."""
import sympy as S, itertools, sys, json
from fractions import Fraction
c,w1,w2,m3,m4,m6,E1,E2=S.symbols('c w1 w2 m3 m4 m6 E1 E2')
A=9*E1+c*E2; Om=-249*E1**2+w1*E1*E2+w2*E2**2
Phi=A**6+15*A**4*Om+15*A**2*Om**2+Om**3
ident=S.expand(2*Phi+(14*E1+m3*E2)**6+(18*E1+m4*E2)**6+E2**6-(22*E1+m6*E2)**6)
P=S.Poly(ident,E1,E2)
eqs=[S.expand(P.coeff_monomial(E1**(6-k)*E2**k)) for k in range(7)]
assert eqs[0]==0
eqs=eqs[1:]
vars_=[c,w1,w2,m3,m4,m6]
# exact rational verification helper
def verify(vals):
    sub=dict(zip(vars_,vals))
    return all(S.simplify(e.subs(sub))==0 for e in eqs)
# fast evaluation mod p
funcs=[S.lambdify(vars_,e,'math') for e in eqs]
polys=[S.Poly(e,*vars_) for e in eqs]
def eval_mod(vals,p):
    out=[]
    for pl in polys:
        s=0
        for mon,coef in pl.terms():
            t=int(coef)%p
            for v,ex in zip(vals,mon):
                if ex: t=t*pow(v,ex,p)%p
            s=(s+t)%p
        out.append(s)
    return out
# the first equation is linear: solve for m6 mod p
lin=polys[0]
print("linear eq:",eqs[0])
import numpy as np
terms=[[(tuple(int(e) for e in mon),int(cf)) for mon,cf in pl.terms()] for pl in polys]
def solve_mod_p(p):
    d=lin.as_dict()
    coefs={v:0 for v in vars_}; const=0
    for mon,cf in d.items():
        if sum(mon)==0: const=int(cf)
        else:
            i=mon.index(1); coefs[vars_[i]]=int(cf)
    cm6=coefs[m6]%p
    assert cm6!=0, "choose a prime not dividing the m6 coefficient"
    inv=pow(cm6,-1,p)
    g=np.indices((p,)*5).reshape(5,-1).astype(np.int64)
    cc,ww1,ww2,mm3,mm4=g
    mm6=(-(const+coefs[c]*cc+coefs[w1]*ww1+coefs[w2]*ww2+coefs[m3]*mm3+coefs[m4]*mm4))%p*inv%p
    V=[cc,ww1,ww2,mm3,mm4,mm6]
    pows=[[None]*7 for _ in range(6)]
    for i in range(6):
        pows[i][0]=np.ones_like(V[i]); 
        for e in range(1,7): pows[i][e]=pows[i][e-1]*V[i]%p
    mask=np.ones(len(cc),dtype=bool)
    for k in range(1,6):
        acc=np.zeros(len(cc),dtype=np.int64)
        for mon,cf in terms[k]:
            t=np.full(len(cc),cf%p,dtype=np.int64)
            for i,e in enumerate(mon):
                if e: t=t*pows[i][e]%p
            acc=(acc+t)%p
        mask&=(acc==0)
    idx=np.nonzero(mask)[0]
    return [tuple(int(V[i][j]) for i in range(6)) for j in idx]
# Jacobian for Hensel
Jsym=S.Matrix(eqs).jacobian(vars_)
Jpolys=[[S.Poly(Jsym[i,j],*vars_) for j in range(6)] for i in range(6)]
def evalpoly_mod(pl,vals,M):
    s=0
    for mon,coef in pl.terms():
        t=int(coef)%M
        for v,ex in zip(vals,mon):
            if ex: t=t*pow(v,ex,M)%M
        s=(s+t)%M
    return s
def solve_linear_mod_p(Jm,rhs,p):
    # gaussian elimination mod p, Jm 6x6 list, rhs list ; returns solution or None
    n=6; Aug=[row[:]+[rhs[i]] for i,row in enumerate(Jm)]
    r=0; piv=[]
    for col in range(n):
        pr=next((i for i in range(r,n) if Aug[i][col]%p),None)
        if pr is None: continue
        Aug[r],Aug[pr]=Aug[pr],Aug[r]
        inv=pow(Aug[r][col],-1,p); Aug[r]=[x*inv%p for x in Aug[r]]
        for i in range(n):
            if i!=r and Aug[i][col]%p:
                f=Aug[i][col]; Aug[i]=[(x-f*y)%p for x,y in zip(Aug[i],Aug[r])]
        piv.append(col); r+=1
    if r<n: return None
    sol=[0]*n
    for i,col in enumerate(piv): sol[col]=Aug[i][n]
    return sol
def hensel(vals,p,k):
    M=p
    vals=list(vals)
    for step in range(1,k):
        M2=M*p
        F=[evalpoly_mod(pl,vals,M2) for pl in polys]
        assert all(f%M==0 for f in F)
        Jm=[[evalpoly_mod(Jpolys[i][j],vals,p) for j in range(6)] for i in range(6)]
        rhs=[(-(f//M))%p for f in F]
        d=solve_linear_mod_p(Jm,rhs,p)
        if d is None: return None
        vals=[(v+M*dd)%M2 for v,dd in zip(vals,d)]
        M=M2
    return vals,M
def ratrec(a,M):
    # rational reconstruction: find n/d = a mod M with |n|,d <= sqrt(M/2)
    import math
    B=math.isqrt(M//2)
    r0,r1=M,a%M; s0,s1=0,1
    while r1>B:
        q=r0//r1; r0,r1=r1,r0-q*r1; s0,s1=s1,s0-q*s1
    if s1==0 or abs(s1)>B: return None
    return Fraction(r1,s1) if s1>0 else Fraction(-r1,-s1)
import time
for p in [int(x) for x in sys.argv[1:]] or [11,13]:
    t=time.time(); sols=solve_mod_p(p)
    print("p=%d: %d solutions mod p (%.0fs)"%(p,len(sols),time.time()-t),flush=True)
    found=[]
    for s in sols:
        h=hensel(s,p,14)
        if h is None: continue
        vals,M=h
        rec=[ratrec(v,M) for v in vals]
        if all(r is not None for r in rec):
            if verify(rec):
                print("RATIONAL SOLUTION:",rec,flush=True); found.append(rec)
    print("  rational solutions reconstructed:",len(found),flush=True)
