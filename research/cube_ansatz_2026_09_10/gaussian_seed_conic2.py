"""Gauge-fixed Family-I conics through the Gaussian point: fix c = 1 (E2-scaling) and m5 = 0 (shear).
Unknowns w1,w2,m3,m4,m6 (5); equations 6 (overdetermined by one: a generic point has none)."""
import sympy as sp, numpy as np, sys, math
from fractions import Fraction
w1,w2,m3,m4,m6,E1,E2=sp.symbols('w1 w2 m3 m4 m6 E1 E2')
A=7*E1+E2; Om=-121*E1**2+w1*E1*E2+w2*E2**2
Phi=A**6+15*A**4*Om+15*A**2*Om**2+Om**3
ident=sp.expand(2*Phi+(8*E1+m3*E2)**6+(12*E1+m4*E2)**6+(15*E1)**6-(17*E1+m6*E2)**6)
P=sp.Poly(ident,E1,E2); eqs=[sp.expand(P.coeff_monomial(E1**(6-k)*E2**k)) for k in range(7)][1:]
VARS=[w1,w2,m3,m4,m6]
terms=[[(tuple(int(x) for x in mon),int(cf)) for mon,cf in sp.Poly(e,*VARS).terms()] for e in eqs]
J=sp.Matrix(eqs).jacobian(sp.Matrix(VARS))
JT=[[[(tuple(int(x) for x in mon),int(cf)) for mon,cf in sp.Poly(J[i,j],*VARS).terms()] for j in range(5)] for i in range(6)]
def evalp(T,vals,M):
    acc=0
    for mon,cf in T:
        tt=cf%M
        for vv,e in zip(vals,mon):
            if e: tt=tt*pow(vv,e,M)%M
        acc=(acc+tt)%M
    return acc
def solve_mod_p(p):
    g=np.indices((p,)*5).reshape(5,-1).astype(np.int64)
    pows=[[np.ones_like(g[0])] for _ in range(5)]
    for i in range(5):
        for e in range(1,7): pows[i].append(pows[i][-1]*g[i]%p)
    mask=np.ones(g.shape[1],dtype=bool)
    for T in terms:
        acc=np.zeros(g.shape[1],dtype=np.int64)
        for mon,cf in T:
            tt=np.full(g.shape[1],cf%p,dtype=np.int64)
            for i,e in enumerate(mon):
                if e: tt=tt*pows[i][e]%p
            acc=(acc+tt)%p
        mask&=(acc==0)
    sols=[tuple(int(g[i][j]) for i in range(5)) for j in np.nonzero(mask)[0]]
    return [s for s in sols if (s[0]*s[0]+484*s[1])%p!=0]   # nondegenerate Omega
def hensel(sol,p,k):
    vals=list(sol); M=p
    for _ in range(k-1):
        M2=M*p; Fv=[evalp(T,vals,M2) for T in terms]
        if any(f%M for f in Fv): return None
        Aug=[[evalp(JT[i][j],vals,p) for j in range(5)]+[(-(Fv[i]//M))%p] for i in range(6)]
        r=0; piv=[]
        for col in range(5):
            pr=next((i for i in range(r,6) if Aug[i][col]%p),None)
            if pr is None: continue
            Aug[r],Aug[pr]=Aug[pr],Aug[r]; iv=pow(Aug[r][col],-1,p); Aug[r]=[x*iv%p for x in Aug[r]]
            for i in range(6):
                if i!=r and Aug[i][col]%p:
                    f=Aug[i][col]; Aug[i]=[(x-f*y)%p for x,y in zip(Aug[i],Aug[r])]
            piv.append(col); r+=1
        if any(all(x==0 for x in row[:5]) and row[5] for row in Aug): return None   # inconsistent
        d=[0]*5
        for i,col in enumerate(piv): d[col]=Aug[i][5]
        vals=[(v+M*dd)%M2 for v,dd in zip(vals,d)]; M=M2
    return vals,M
def ratrec(a,M):
    B=math.isqrt(M//2); r0,r1=M,a%M; s0,s1=0,1
    while r1>B:
        q=r0//r1; r0,r1=r1,r0-q*r1; s0,s1=s1,s0-q*s1
    if s1==0 or abs(s1)>B: return None
    return Fraction(r1,s1) if s1>0 else Fraction(-r1,-s1)
for p in [int(x) for x in sys.argv[1:]] or [13,17,19,23]:
    s=solve_mod_p(p); print("p=%d: nondegenerate solutions mod p: %d"%(p,len(s)),flush=True)
    lifted=0; rat=0
    for sol in s:
        h=hensel(sol,p,12)
        if h is None: continue
        lifted+=1; vals,M=h; rec=[ratrec(v,M) for v in vals]
        if all(r is not None for r in rec):
            sub=dict(zip(VARS,rec))
            if all(sp.simplify(e.subs(sub))==0 for e in eqs):
                rat+=1; print("   !!!! RATIONAL CONIC THROUGH THE GAUSSIAN POINT:",[str(r) for r in rec],flush=True)
    print("   lift to p^12: %d ; rational: %d"%(lifted,rat),flush=True)
