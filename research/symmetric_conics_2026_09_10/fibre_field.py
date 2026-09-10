"""Fix a,c rational in Family I; Newton-solve remaining 7 unknowns in high precision; identify minimal polynomials."""
import numpy as np, mpmath as mp, sys
sys.path.insert(0,'research/symmetric_conics_2026_09_10')
from families import resI, numjac, newton
mp.mp.dps=60
X=np.load('/tmp/claude-0/-home-user/15becd99-5cb2-5637-beeb-15d2ba0dfed3/scratchpad/famI_real.npy')
def res_mp(v7,a,c):
    w0,w1,w2,l3,m3,l4,m4=v7
    def pw(P,k):
        r=[mp.mpf(1)]
        for _ in range(k):
            r=[sum(r[i]*P[j] for i in range(len(r)) for j in range(len(P)) if i+j==n) for n in range(len(r)+len(P)-1)]
        return r
    def mul(P,Q): return [sum(P[i]*Q[j] for i in range(len(P)) for j in range(len(Q)) if i+j==n) for n in range(len(P)+len(Q)-1)]
    def pad(P,n=7): return P+[mp.mpf(0)]*(n-len(P))
    A=[mp.mpf(a),mp.mpf(c)]; Om=[w0,w1,w2]
    Phi=[x+15*y+15*z+t for x,y,z,t in zip(pad(pw(A,6)),pad(mul(pw(A,4),Om)),pad(mul(pw(A,2),pw(Om,2))),pad(pw(Om,3)))]
    r=[2*p for p in Phi]
    for k,x in enumerate(pw([l3,m3],6)): r[k]+=x
    for k,x in enumerate(pw([l4,m4],6)): r[k]+=x
    r[6]+=1; r[0]-=1
    return r
def newton_mp(v7,a,c,iters=60):
    v=[mp.mpf(x) for x in v7]
    for it in range(iters):
        r=res_mp(v,a,c)
        nr=max(abs(x) for x in r)
        if nr<mp.mpf(10)**(-55): return v,nr
        J=mp.matrix(7,7); h=mp.mpf(10)**(-30)
        for k in range(7):
            vv=v[:]; vv[k]+=h
            rk=res_mp(vv,a,c)
            for i in range(7): J[i,k]=(rk[i]-r[i])/h
        try: dv=mp.lu_solve(J,mp.matrix([-x for x in r]))
        except ZeroDivisionError: return v,nr
        v=[v[k]+dv[k] for k in range(7)]
    return v,max(abs(x) for x in res_mp(v,a,c))
import random
random.seed(1)
done=0
for x in X[:40]:
    a,c=x[0].real,x[1].real
    # snap a,c to nearby rationals with small denominators
    fa=mp.mpf(round(a*4)/4); fc=mp.mpf(round(c*4)/4)
    if fa==0 and fc==0: continue
    v,nr=newton_mp(x[2:].real,fa,fc)
    if nr>mp.mpf(10)**(-40): continue
    w0,w1,w2,l3,m3,l4,m4=v
    if abs(w1*w1-4*w0*w2)<1e-8: continue
    degs=[]
    for val in v:
        pol=mp.findpoly(val,8,maxcoeff=10**6)
        degs.append(len(pol)-1 if pol else None)
    print("a=%s c=%s  minpoly degrees of (w0,w1,w2,l3,m3,l4,m4):"%(fa,fc),degs, " w=",[mp.nstr(t,8) for t in v[:3]])
    done+=1
    if done>=12: break
