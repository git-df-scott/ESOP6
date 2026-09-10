"""Symmetric conic families on x1^6+..+x5^6=x6^6, written on the quotient line.

Coordinates (E1,E2) on the quotient P^1 are chosen as (x6, x5).  The conic is the
double cover O^2 = Omega(E1,E2), Omega a binary quadratic.
 Family I  (one swapped pair):  x1=A+O, x2=A-O, x3=L3, x4=L4, x5=E2, x6=E1
 Family II (two swapped pairs): x1=A+O, x2=A-O, x3=D+kO, x4=D-kO, x5=E2, x6=E1
with A,D,L3,L4 linear forms in (E1,E2).  Writing F=k^2 and
 Phi(A,B,Om) = A^6 + 15 A^4 B Om + 15 A^2 B^2 Om^2 + B^3 Om^3  (= ((A+sqrt(B)O)^6+(A-sqrt(B)O)^6)/2),
the sextic identities are
 I : 2 Phi(A,1,Om) + L3^6 + L4^6 + E2^6 - E1^6 = 0        (9 unknowns, 7 eqs)
 II: 2 Phi(A,1,Om) + 2 Phi(D,F,Om) + E2^6 - E1^6 = 0      (8 unknowns, 7 eqs)
"""
import numpy as np
def pm(a,b): return np.convolve(a,b)
def pw(a,k):
    r=np.array([1.0+0j])
    for _ in range(k): r=pm(r,a)
    return r
def pad(a,n=7):
    r=np.zeros(n,dtype=complex); r[:len(a)]=a; return r
E1=np.array([1.0,0]); E2=np.array([0,1.0])
def Phi(A,B,Om):
    return pad(pw(A,6))+15*B*pad(pm(pw(A,4),Om))+15*B*B*pad(pm(pw(A,2),pw(Om,2)))+B**3*pad(pw(Om,3))
def resI(v):
    a,c,w0,w1,w2,l3,m3,l4,m4=v
    A=np.array([a,c]); Om=np.array([w0,w1,w2])
    return 2*Phi(A,1,Om)+pad(pw(np.array([l3,m3]),6))+pad(pw(np.array([l4,m4]),6))+pad(pw(E2,6))-pad(pw(E1,6))
def resII(v):
    a,c,d,f,F,w0,w1,w2=v
    A=np.array([a,c]); D=np.array([d,f]); Om=np.array([w0,w1,w2])
    return 2*Phi(A,1,Om)+2*Phi(D,F,Om)+pad(pw(E2,6))-pad(pw(E1,6))
def numjac(res,v,h=1e-7):
    r0=res(v); J=np.zeros((len(r0),len(v)),dtype=complex)
    for k in range(len(v)):
        dv=np.zeros(len(v),dtype=complex); dv[k]=h; J[:,k]=(res(v+dv)-r0)/h
    return J
def newton(res,v,iters=100):
    for it in range(iters):
        r=res(v)
        if np.linalg.norm(r)<1e-13*(1+np.linalg.norm(v)**6): return v,np.linalg.norm(r)
        v=v+np.linalg.lstsq(numjac(res,v),-r,rcond=None)[0]
        if np.linalg.norm(v)>1e3: break
    return v,np.linalg.norm(res(v))
if __name__=="__main__":
    import sys
    from collections import Counter
    for name,res,n in [("I",resI,9),("II",resII,8)]:
        for real in (False,True):
            rng=np.random.default_rng(7)
            sols=[]
            for s in range(400):
                v=rng.standard_normal(n)+(0 if real else 1j*rng.standard_normal(n))
                v,nr=newton(res,v)
                if nr<1e-10*(1+np.linalg.norm(v)**6):
                    J=numjac(res,v); sv=np.linalg.svd(J,compute_uv=False); dim=n-(sv>1e-7*sv[0]).sum()
                    if name=="I":
                        a,c,w0,w1,w2,l3,m3,l4,m4=v; degen=min(abs(w0)+abs(w1)+abs(w2),abs(l3)+abs(m3),abs(l4)+abs(m4),abs(a)+abs(c))
                    else:
                        a,c,d,f,F,w0,w1,w2=v; degen=min(abs(w0)+abs(w1)+abs(w2),abs(d)+abs(f),abs(a)+abs(c),abs(F))
                    sols.append((dim,degen>1e-3,v))
            print("family",name,"real" if real else "complex",": converged",len(sols),
                  " nondegenerate:",sum(1 for d,g,v in sols if g),
                  " local dims of nondegenerate:",Counter(d for d,g,v in sols if g))
            np.save("/tmp/claude-0/-home-user/15becd99-5cb2-5637-beeb-15d2ba0dfed3/scratchpad/fam%s_%s.npy"%(name,"real" if real else "cplx"),np.array([v for d,g,v in sols if g]))
