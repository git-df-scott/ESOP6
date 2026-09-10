"""Quadratic points on the (6,1,4) slice with a conjugate pair:
   x1,x2 = A +- sqrt(Om),  2*Phi(A,Om) + r3^6 + r4^6 = r6^6,  Phi=A^6+15A^4Om+15A^2Om^2+Om^3.
Om a perfect square would be a rational (6,1,4) point (conjecturally nonexistent)."""
import sys, math
N=int(sys.argv[1]) if len(sys.argv)>1 else 60
S6=[r**6 for r in range(N+1)]
def int_root_cubic(A,K):
    x=K**(1/3)-5*A*A
    for _ in range(80):
        f=x**3+15*A*A*x*x+15*A**4*x+A**6-K; df=3*x*x+30*A*A*x+15*A**4
        if df==0: break
        x-=f/df
    r=round(x)
    for om in (r-1,r,r+1):
        if om**3+15*A*A*om*om+15*A**4*om+A**6==K: return om
    return None
hits=[]; cnt=0
for r6 in range(1,N+1):
    for r3 in range(0,r6):
        for r4 in range(r3,r6):
            M=S6[r6]-S6[r3]-S6[r4]
            if M<=0 or M%2: continue
            K=M//2; Amax=int(K**(1/6))+1
            for A in range(0,Amax+1):
                cnt+=1; om=int_root_cubic(A,K)
                if om is not None and om!=0:
                    s=math.isqrt(abs(om)); sq=(om>0 and s*s==om)
                    hits.append((A,om,r3,r4,r6,"SQUARE" if sq else "nonsquare"))
print("checked",cnt,"hits",len(hits))
for h in hits: print(h)
