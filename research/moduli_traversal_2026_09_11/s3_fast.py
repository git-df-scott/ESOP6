"""Fast exact search on the S3-symmetric sextic family (setup as in s3_exact.py).
x5^6 = w^6 * B5^6 and x6 = B6a + z2*B6b.  For fixed (u1,u2) and enumerated (u3,u4,v1,v2) mod p,
G := x1^6+x2^6+x3^6+x4^6 is known; for each z2 the 37-vector R=(B6a+z2 B6b)^6 - G must equal W*B5^6 with W a sixth power."""
import sympy as sp, numpy as np, sys, math, time, json
from fractions import Fraction
from math import comb
d=6
def actM_vec(c):
    out=[0]*(d+1)
    for k,ck in enumerate(c):
        for j in range(k+1): out[d-k+j]+=ck*comb(k,j)*(-1)**(k-j)
    return out
E=sp.eye(d+1)
AM=sp.Matrix([actM_vec(list(E.col(i))) for i in range(d+1)]).T
AI=sp.Matrix([list(reversed(list(E.col(i)))) for i in range(d+1)]).T
def space(chiI,chiM):
    K=sp.Matrix.vstack(AM-chiM*sp.eye(d+1),AI-chiI*sp.eye(d+1)) if chiM is not None else AI-chiI*sp.eye(d+1)
    out=[]
    for v in K.nullspace():
        den=sp.ilcm(*[sp.fraction(x)[1] for x in v]); out.append([int(x*den) for x in v])
    return out
B1=space(+1,None); B4=space(+1,1); B5=space(-1,1); B6=space(+1,1)
print("bases",B1,B4,B5,B6,flush=True)
MM=np.array([actM_vec(list(np.eye(d+1,dtype=np.int64)[i])) for i in range(d+1)],dtype=np.int64).T
def polymul_mod(A,B,p):
    N=A.shape[0]; a=A.shape[1]; B=np.broadcast_to(B,(N,B.shape[-1])) if B.ndim==1 else B; b=B.shape[1]
    out=np.zeros((N,a+b-1),dtype=np.int64)
    for i in range(a):
        for j in range(b): out[:,i+j]=(out[:,i+j]+A[:,i]*B[:,j])%p
    return out
def polypow_mod(A,k,p):
    out=np.ones((A.shape[0],1),dtype=np.int64)
    for _ in range(k): out=polymul_mod(out,A,p)
    return out
def search_mod_p(u1,u2,p):
    inv=lambda x: pow(int(x)%p,-1,p)
    U1=(u1.numerator*inv(u1.denominator))%p; U2=(u2.numerator*inv(u2.denominator))%p
    g=np.indices((p,)*4).reshape(4,-1).astype(np.int64); u3,u4,v1,v2=g
    B1a=np.array(B1,dtype=np.int64)
    x1=(U1*B1a[0][None,:]+U2*B1a[1][None,:]+u3[:,None]*B1a[2][None,:]+u4[:,None]*B1a[3][None,:])%p
    x2=(x1@MM.T)%p; x3=(x2@MM.T)%p
    B4a=np.array(B4,dtype=np.int64); x4=(v1[:,None]*B4a[0][None,:]+v2[:,None]*B4a[1][None,:])%p
    G=(polypow_mod(x1,6,p)+polypow_mod(x2,6,p)+polypow_mod(x3,6,p)+polypow_mod(x4,6,p))%p
    B5v=np.array(B5[0],dtype=np.int64)%p; B56=polypow_mod(B5v[None,:],6,p)[0]
    k0=int(np.nonzero(B56)[0][0]); invk=inv(B56[k0])
    B6a=np.array(B6[0],dtype=np.int64)%p; B6b=np.array(B6[1],dtype=np.int64)%p
    S6=np.array(sorted({pow(x,6,p) for x in range(1,p)}),dtype=np.int64); sols=[]
    for z2 in range(p):
        X66=polypow_mod(((B6a+z2*B6b)%p)[None,:],6,p)[0]
        R=(X66[None,:]-G)%p; W=(R[:,k0]*invk)%p
        ok=np.all((R-(W[:,None]*B56[None,:])%p)%p==0,axis=1)&np.isin(W,S6)
        for j in np.nonzero(ok)[0]:
            for wv in range(1,p):
                if pow(wv,6,p)==int(W[j]): sols.append(dict(u1=U1,u2=U2,u3=int(u3[j]),u4=int(u4[j]),v1=int(v1[j]),v2=int(v2[j]),w=wv,z2=z2))
    return sols
s,t=sp.symbols('s t'); u=sp.symbols('u1:5'); v=sp.symbols('v1:3'); w=sp.symbols('w'); z2s=sp.symbols('z2')
VARS=list(u)+list(v)+[w,z2s]
def form(c): return sum(c[k]*s**(d-k)*t**k for k in range(d+1))
def comb_vec(B,coefs): return [sp.expand(sum(c*b[k] for c,b in zip(coefs,B))) for k in range(d+1)]
x1s=comb_vec(B1,u); x4s=comb_vec(B4,v); x5s=comb_vec(B5,[w]); x6s=comb_vec(B6,[1,z2s])
x2s=actM_vec(x1s); x3s=actM_vec(x2s)
Fs=[form(c) for c in [x1s,x2s,x3s,x4s,x5s,x6s]]
ident=sp.expand(sum(f**6 for f in Fs[:5])-Fs[5]**6); Pl=sp.Poly(ident,s,t)
eqs=[e for e in (sp.expand(Pl.coeff_monomial(s**(36-k)*t**k)) for k in range(37)) if e!=0]
terms=[[(tuple(int(x) for x in mon),int(c)) for mon,c in sp.Poly(e,*VARS).terms()] for e in eqs]
Jm=sp.Matrix(eqs).jacobian(sp.Matrix(VARS))
JT=[[[(tuple(int(x) for x in mon),int(c)) for mon,c in sp.Poly(Jm[i,j],*VARS).terms()] for j in range(8)] for i in range(len(eqs))]
def evalp(T,vals,M):
    acc=0
    for mon,c in T:
        tt=c%M
        for vv,e in zip(vals,mon):
            if e: tt=tt*pow(vv,e,M)%M
        acc=(acc+tt)%M
    return acc
def hensel(sol,u1,u2,p,k):
    vals=[sol[str(vv)] for vv in VARS]; M=p; free=[2,3,4,5,6,7]
    for _ in range(k-1):
        M2=M*p
        vals[0]=(u1.numerator*pow(u1.denominator,-1,M2))%M2; vals[1]=(u2.numerator*pow(u2.denominator,-1,M2))%M2
        Fv=[evalp(T,vals,M2) for T in terms]
        if any(f%M for f in Fv): return None
        A=[[evalp(JT[i][j],vals,p) for j in free]+[(-(Fv[i]//M))%p] for i in range(len(eqs))]
        r=0; piv=[]; n=6
        for col in range(n):
            pr=next((i for i in range(r,len(A)) if A[i][col]%p),None)
            if pr is None: continue
            A[r],A[pr]=A[pr],A[r]; iv=pow(A[r][col],-1,p); A[r]=[x*iv%p for x in A[r]]
            for i in range(len(A)):
                if i!=r and A[i][col]%p:
                    f=A[i][col]; A[i]=[(x-f*y)%p for x,y in zip(A[i],A[r])]
            piv.append(col); r+=1
        if r<n or any(all(x==0 for x in row[:n]) and row[n] for row in A): return None
        dv=[0]*n
        for i,col in enumerate(piv): dv[col]=A[i][n]
        for kk,i in enumerate(free): vals[i]=(vals[i]+M*dv[kk])%M2
        M=M2
    return vals,M
def ratrec(a,M):
    B=math.isqrt(M//2); r0,r1=M,a%M; s0,s1=0,1
    while r1>B:
        q=r0//r1; r0,r1=r1,r0-q*r1; s0,s1=s1,s0-q*s1
    if s1==0 or abs(s1)>B: return None
    return Fraction(r1,s1) if s1>0 else Fraction(-r1,-s1)
def verify(vals):
    sub=dict(zip(VARS,vals)); return all(sp.simplify(e.subs(sub))==0 for e in eqs)
if __name__=="__main__":
    H=int(sys.argv[1]) if len(sys.argv)>1 else 3; p=int(sys.argv[2]) if len(sys.argv)>2 else 11
    grid=sorted({Fraction(a,b) for b in range(1,H+1) for a in range(-H,H+1)})
    print("grid",len(grid)**2,"points; p =",p,flush=True); t0=time.time(); hits=[]; n=0; nsol=0
    for a in grid:
        for b in grid:
            sols=search_mod_p(a,b,p); n+=1; nsol+=len(sols)
            for sol in sols:
                h=hensel(sol,a,b,p,9)
                if h is None: continue
                vals,M=h; rec=[ratrec(x,M) for x in vals]
                if all(r is not None for r in rec) and verify(rec):
                    print("!!!! RATIONAL S3-SYMMETRIC SEXTIC ON X:",[str(r) for r in rec],flush=True); hits.append([str(r) for r in rec])
            if n%20==0: print("done",n,"grid pts, mod-p solutions so far",nsol,"elapsed %.0fs"%(time.time()-t0),flush=True)
    json.dump(hits,open('research/moduli_traversal_2026_09_11/s3_hits_H%d_p%d.json'%(H,p),'w')); print("finished; hits",len(hits),"mod-p sols",nsol)
