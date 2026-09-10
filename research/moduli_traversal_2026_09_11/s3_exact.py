"""Exact rational-point search on the S3-symmetric sextic family (chars x1:+1 under iota; x4:+1, x5:sign, x6:+1).
Rational bases of the character spaces; unknowns u1..u4 (x1), v1,v2 (x4), w (x5), z1,z2 (x6) with z1 = 1 (scaling).
Fix (u1,u2) on a rational grid; enumerate the remaining 6 unknowns mod p; Hensel; reconstruct; verify exactly."""
import sympy as sp, numpy as np, itertools, sys, math, time, json
from fractions import Fraction
from math import comb
d=6
s,t=sp.symbols('s t')
def form(c): return sum(sp.Integer(1)*c[k]*s**(d-k)*t**k for k in range(d+1))
def actM_vec(c):
    out=[0]*(d+1)
    for k,ck in enumerate(c):
        for j in range(k+1): out[d-k+j]+=ck*comb(k,j)*(-1)**(k-j)
    return out
def actI_vec(c): return list(reversed(c))
E=sp.eye(d+1)
AM=sp.Matrix([actM_vec(list(E.col(i))) for i in range(d+1)]).T
AI=sp.Matrix([actI_vec(list(E.col(i))) for i in range(d+1)]).T
def space(chiI,chiM):
    K=sp.Matrix.vstack(AM-chiM*sp.eye(d+1),AI-chiI*sp.eye(d+1)) if chiM is not None else AI-chiI*sp.eye(d+1)
    return K.nullspace()
B1=space(+1,None); B4=space(+1,1); B5=space(-1,1); B6=space(+1,1)
print("basis dims",len(B1),len(B4),len(B5),len(B6))
u=sp.symbols('u1:5'); v=sp.symbols('v1:3'); w=sp.symbols('w'); z2=sp.symbols('z2')
def comb_vec(B,coefs): return [sp.expand(sum(c*b[k] for c,b in zip(coefs,B))) for k in range(d+1)]
x1=comb_vec(B1,u); x4=comb_vec(B4,v); x5=comb_vec(B5,[w]); x6=comb_vec(B6,[1,z2])
x2=actM_vec(x1); x3=actM_vec(x2)
X=[x1,x2,x3,x4,x5,x6]
F=[form(c) for c in X]
ident=sp.expand(sum(f**6 for f in F[:5])-F[5]**6)
P=sp.Poly(ident,s,t)
eqs=[sp.expand(P.coeff_monomial(s**(36-k)*t**k)) for k in range(37)]
eqs=[e for e in eqs if e!=0]
VARS=list(u)+list(v)+[w,z2]
J=sp.Matrix(eqs).jacobian(sp.Matrix(VARS))
print("J shape",J.shape,flush=True)
print("nonzero coefficient equations:",len(eqs))
# generic rank
sub={vv:sp.Rational(int(np.random.randint(1,50)),int(np.random.randint(1,50))) for vv in VARS}
print("generic Jacobian rank:",J.subs(sub).rank())
polys=[sp.Poly(e,*VARS) for e in eqs]
terms=[[(tuple(int(x) for x in mon),int(c)) for mon,c in pl.terms()] for pl in polys]
json.dump({'eqs':[str(e) for e in eqs],'vars':[str(v) for v in VARS]},open('research/moduli_traversal_2026_09_11/s3_system.json','w'))
def solve_mod_p(fixed,p):
    """fixed: dict var-> Fraction for u1,u2 ; enumerate u3,u4,v1,v2,w,z2 mod p"""
    free=[VARS[i] for i in range(8) if VARS[i] not in fixed]
    assert len(free)==6
    g=np.indices((p,)*6).reshape(6,-1).astype(np.int64)
    V={}
    for k,var in enumerate(free): V[var]=g[k]
    for var,val in fixed.items():
        V[var]=np.full(g.shape[1],(val.numerator*pow(val.denominator,-1,p))%p,dtype=np.int64)
    vals=[V[var] for var in VARS]
    pows=[[np.ones_like(vals[0])] for _ in range(8)]
    for i in range(8):
        for e in range(1,7): pows[i].append(pows[i][-1]*vals[i]%p)
    mask=np.ones(g.shape[1],dtype=bool)
    for T in terms:
        acc=np.zeros(g.shape[1],dtype=np.int64)
        for mon,c in T:
            tt=np.full(g.shape[1],c%p,dtype=np.int64)
            for i,e in enumerate(mon):
                if e: tt=tt*pows[i][e]%p
            acc=(acc+tt)%p
        mask&=(acc==0)
        if not mask.any(): break
    idx=np.nonzero(mask)[0]
    return [ {var:int(V[var][j]) for var in VARS} for j in idx]
def evalp(T,vals,M):
    acc=0
    for mon,c in T:
        tt=c%M
        for vv,e in zip(vals,mon):
            if e: tt=tt*pow(vv,e,M)%M
        acc=(acc+tt)%M
    return acc
JT=[[ [(tuple(int(x) for x in mon),int(c)) for mon,c in sp.Poly(J[i,j],*VARS).terms()] for j in range(8)] for i in range(len(eqs))]
def hensel(sol,fixed,p,k):
    free=[i for i in range(8) if VARS[i] not in fixed]
    vals=[sol[var] for var in VARS]; M=p
    for _ in range(k-1):
        M2=M*p
        # fixed vars: exact rational mod M2
        for i,var in enumerate(VARS):
            if var in fixed: vals[i]=(fixed[var].numerator*pow(fixed[var].denominator,-1,M2))%M2
        Fv=[evalp(T,vals,M2) for T in terms]
        if any(f%M for f in Fv): return None
        A=[[evalp(JT[i][j],vals,p) for j in free]+[(-(Fv[i]//M))%p] for i in range(len(eqs))]
        # solve over F_p (overdetermined, consistent)
        r=0; piv=[]; n=6
        for col in range(n):
            pr=next((i for i in range(r,len(A)) if A[i][col]%p),None)
            if pr is None: continue
            A[r],A[pr]=A[pr],A[r]; iv=pow(A[r][col],-1,p); A[r]=[x*iv%p for x in A[r]]
            for i in range(len(A)):
                if i!=r and A[i][col]%p:
                    f=A[i][col]; A[i]=[(x-f*y)%p for x,y in zip(A[i],A[r])]
            piv.append(col); r+=1
        if r<n: return None
        if any(all(x==0 for x in row[:n]) and row[n] for row in A): return None
        dvec=[0]*n
        for i,col in enumerate(piv): dvec[col]=A[i][n]
        for k2,i in enumerate(free): vals[i]=(vals[i]+M*dvec[k2])%M2
        M=M2
    return vals,M
def ratrec(a,M):
    B=math.isqrt(M//2); r0,r1=M,a%M; s0,s1=0,1
    while r1>B:
        q=r0//r1; r0,r1=r1,r0-q*r1; s0,s1=s1,s0-q*s1
    if s1==0 or abs(s1)>B: return None
    return Fraction(r1,s1) if s1>0 else Fraction(-r1,-s1)
def verify(vals):
    sub=dict(zip(VARS,vals))
    return all(sp.simplify(e.subs(sub))==0 for e in eqs)
if __name__=="__main__":
    H=int(sys.argv[1]) if len(sys.argv)>1 else 4; p=int(sys.argv[2]) if len(sys.argv)>2 else 11
    grid=sorted({Fraction(a,b) for b in range(1,H+1) for a in range(-H,H+1)})
    print("grid size",len(grid)**2,flush=True)
    hits=[]; t0=time.time(); n=0
    for a in grid:
        for b in grid:
            fixed={VARS[0]:a,VARS[1]:b}
            sols=solve_mod_p(fixed,p); n+=1
            for sol in sols:
                h=hensel(sol,fixed,p,10)
                if h is None: continue
                vals,M=h; rec=[ratrec(x,M) for x in vals]
                if all(r is not None for r in rec) and verify(rec):
                    X6=[float(x) for x in comb_vec(B6,[1,rec[8]])]
                    print("!!!! RATIONAL S3-SYMMETRIC SEXTIC:",[str(r) for r in rec],flush=True); hits.append([str(r) for r in rec])
            if n%50==0: print("grid points done",n,"elapsed %.0fs"%(time.time()-t0),flush=True)
    json.dump(hits,open('research/moduli_traversal_2026_09_11/s3_hits_H%d_p%d.json'%(H,p),'w'))
    print("done; hits",len(hits))
