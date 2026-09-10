"""Plane cubic elliptic curves on X through a rational point of Y2.

Y2: x1^6+x2^6+x3^6+x4^6 + S^2 = x6^6  (S = x5^3).  A weighted line on Y2 through y0 = (x_i, x6; s0):
   L_i = x_i p0 + m_i p1 (i=1..4), L6 = x6 p0 + p1 (m6 = 1 by scaling), m1 = 0 (by shearing p0 -> p0 + b p1),
   S = s0 p0^3 + s1 p0^2 p1 + s2 p0 p1^2 + s3 p1^3,
   identity  L6^6 - L1^6 - L2^6 - L3^6 - L4^6 = S^2.
The plane spanned by the five linear forms meets X in the two cubics u^3 = ±S; the curve (L1,..,L4, u, L6) with
u^3 = S(p0,p1) is a plane cubic (genus 1) over Q on X.  Any rational point with u != 0 is an ESOP6 solution.
Unknowns m2,m3,m4,s1,s2,s3; six coefficient equations (p0^5p1 .. p1^6); the first is linear in (m2,m3,m4,s1).
Exact pipeline: enumerate solutions mod p (vectorised), Hensel-lift, rationally reconstruct, verify exactly."""
import sympy as sp, numpy as np, sys, json, math, time
from fractions import Fraction
m2,m3,m4,s1,s2,s3,p0,p1=sp.symbols('m2 m3 m4 s1 s2 s3 p0 p1')
VARS=[m2,m3,m4,s1,s2,s3]
def build(x,x6,s0):
    L=[x[0]*p0, x[1]*p0+m2*p1, x[2]*p0+m3*p1, x[3]*p0+m4*p1]
    L6=x6*p0+p1; S=s0*p0**3+s1*p0**2*p1+s2*p0*p1**2+s3*p1**3
    ident=sp.expand(L6**6-sum(l**6 for l in L)-S**2)
    P=sp.Poly(ident,p0,p1)
    eqs=[sp.expand(P.coeff_monomial(p0**(6-k)*p1**k)) for k in range(7)]
    assert eqs[0]==0, "seed is not on Y2"
    return eqs[1:]
def terms_of(e): return [(tuple(int(v) for v in mon),int(c)) for mon,c in sp.Poly(e,*VARS).terms()]
def solve_mod_p(eqs,p):
    T=[terms_of(e) for e in eqs]
    lin=T[0]
    coef=[0]*6; const=0
    for mon,c in lin:
        if sum(mon)==0: const=c
        else: coef[mon.index(1)]=c
    # solve the linear equation for s1 (index 3) if its coefficient is a unit mod p, else for another variable
    piv=next((i for i in (3,0,1,2) if coef[i]%p),None)
    if piv is None: return None
    others=[i for i in range(4) if i!=piv]   # variables m2,m3,m4,s1 except pivot ; s2,s3 free
    free=others+[4,5]
    g=np.indices((p,)*5).reshape(5,-1).astype(np.int64)
    V=[None]*6
    for k,i in enumerate(free): V[i]=g[k]
    inv=pow(coef[piv]%p,-1,p)
    V[piv]=(-(const+sum(coef[i]*V[i] for i in others)))%p*inv%p
    pows=[[np.ones_like(V[0])] for _ in range(6)]
    for i in range(6):
        for e in range(1,7): pows[i].append(pows[i][-1]*V[i]%p)
    mask=np.ones(len(V[0]),dtype=bool)
    for k in range(1,6):
        acc=np.zeros(len(V[0]),dtype=np.int64)
        for mon,c in T[k]:
            t=np.full(len(V[0]),c%p,dtype=np.int64)
            for i,e in enumerate(mon):
                if e: t=t*pows[i][e]%p
            acc=(acc+t)%p
        mask&=(acc==0)
    idx=np.nonzero(mask)[0]
    return [tuple(int(V[i][j]) for i in range(6)) for j in idx]
def evalp(e_terms,vals,M):
    s=0
    for mon,c in e_terms:
        t=c%M
        for v,ex in zip(vals,mon):
            if ex: t=t*pow(v,ex,M)%M
        s=(s+t)%M
    return s
def hensel(eqs,vals,p,k):
    T=[terms_of(e) for e in eqs]
    J=sp.Matrix(eqs).jacobian(VARS); JT=[[terms_of(J[i,j]) for j in range(6)] for i in range(6)]
    M=p; vals=list(vals)
    for _ in range(k-1):
        M2=M*p
        F=[evalp(t,vals,M2) for t in T]
        Jm=[[evalp(JT[i][j],vals,p) for j in range(6)] for i in range(6)]
        rhs=[(-(f//M))%p for f in F]
        # solve Jm d = rhs mod p
        A=[row[:]+[rhs[i]] for i,row in enumerate(Jm)]; r=0; piv=[]
        for col in range(6):
            pr=next((i for i in range(r,6) if A[i][col]%p),None)
            if pr is None: continue
            A[r],A[pr]=A[pr],A[r]; iv=pow(A[r][col],-1,p); A[r]=[v*iv%p for v in A[r]]
            for i in range(6):
                if i!=r and A[i][col]%p:
                    f=A[i][col]; A[i]=[(v-f*w)%p for v,w in zip(A[i],A[r])]
            piv.append(col); r+=1
        if r<6: return None
        d=[0]*6
        for i,col in enumerate(piv): d[col]=A[i][6]
        vals=[(v+M*dd)%M2 for v,dd in zip(vals,d)]; M=M2
    return vals,M
def ratrec(a,M):
    B=math.isqrt(M//2); r0,r1=M,a%M; s0,s1=0,1
    while r1>B:
        q=r0//r1; r0,r1=r1,r0-q*r1; s0,s1=s1,s0-q*s1
    if s1==0 or abs(s1)>B: return None
    return Fraction(r1,s1) if s1>0 else Fraction(-r1,-s1)
def run_seed(x,x6,s0,primes=(17,19,23),digits=12):
    eqs=build(x,x6,s0)
    report={}; rational=[]
    for p in primes:
        sols=solve_mod_p(eqs,p)
        if sols is None: report[p]='pivot-fail'; continue
        report[p]=len(sols)
        for s in sols:
            h=hensel(eqs,s,p,digits)
            if h is None: continue
            vals,M=h
            rec=[ratrec(v,M) for v in vals]
            if all(r is not None for r in rec):
                sub=dict(zip(VARS,rec))
                if all(sp.simplify(e.subs(sub))==0 for e in eqs):
                    rational.append(rec)
    # dedupe
    uniq=[]
    for r in rational:
        if r not in uniq: uniq.append(r)
    return report,uniq
if __name__=="__main__":
    pts=json.load(open(sys.argv[1]))
    seen=set(); hits=[]
    for pt in pts:
        key=(tuple(sorted(pt['x'])),pt['x6'])
        if key in seen: continue
        seen.add(key)
        # try each choice of which coordinate carries m=0 ... we fix m1=0 on x[0]; permute to cover cases with x1 changed
        t=time.time()
        report,rat=run_seed(pt['x'],pt['x6'],pt['S'])
        print("seed",key,"S=",pt['S'],"mod-p counts",report,"rational lines:",len(rat),"(%.0fs)"%(time.time()-t),flush=True)
        for r in rat:
            print("   !!!! RATIONAL PLANE CUBIC on X:", "x=",pt['x'],"x6=",pt['x6'],"s0=",pt['S'],"(m2,m3,m4,s1,s2,s3)=",[str(v) for v in r],flush=True)
            hits.append(dict(point=pt,sol=[str(v) for v in r]))
    json.dump(hits,open("research/cube_ansatz_2026_09_10/plane_cubic_hits.json","w"),indent=1)
