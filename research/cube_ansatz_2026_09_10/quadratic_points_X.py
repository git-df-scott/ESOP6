"""Quadratic points on X with a conjugate pair:  x1,x2 = A ± sqrt(Om),  2*Phi(A,Om) + r3^6 + r4^6 + r5^6 = r6^6.
Om not a square (else it's an ESOP6 solution outright).  For each such point, the secant line through P and its
conjugate meets X in four more points: the residual quartic g(nu) in  x1 = A+nu, x2 = A-nu  is
 2*Phi(A, nu^2) + r3^6+r4^6+r5^6 - r6^6 = 0  divided by (nu^2 - Om); a rational root nu gives an ESOP6 solution."""
import sys, math, json
N=int(sys.argv[1]) if len(sys.argv)>1 else 40
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
hits=[]
for r6 in range(1,N+1):
    for r3 in range(1,r6):
        for r4 in range(r3,r6):
            for r5 in range(r4,r6):
                M=S6[r6]-S6[r3]-S6[r4]-S6[r5]
                if M<=0 or M%2: continue
                K=M//2; Amax=int(K**(1/6))+1
                for A in range(0,Amax+1):
                    om=int_root_cubic(A,K)
                    if om is None or om==0: continue
                    s=math.isqrt(abs(om)); sq=(om>0 and s*s==om)
                    h=dict(A=A,Om=om,r=[r3,r4,r5],r6=r6,square=sq)
                    # residual quartic in W=nu^2: (Phi(A,W) - Phi(A,Om))/(W-Om) = W^2 + (15A^2+Om) W + (15A^4+15A^2 Om+Om^2)
                    b=15*A*A+om; c=15*A**4+15*A*A*om+om*om
                    disc=b*b-4*c
                    roots=[]
                    if disc>=0 and math.isqrt(disc)**2==disc:
                        for W in ((-b+math.isqrt(disc))//2,(-b-math.isqrt(disc))//2):
                            if W>0 and math.isqrt(W)**2==W: roots.append(math.isqrt(W))
                    h['rational_residual_nu']=roots
                    hits.append(h); print(h,flush=True)
print("total",len(hits))
json.dump(hits,open("research/cube_ansatz_2026_09_10/quadratic_points_X_N%d.json"%N,"w"),indent=0)
