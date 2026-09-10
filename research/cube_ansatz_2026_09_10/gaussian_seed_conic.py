"""Family I conics through the Gaussian point P=(7±11i,8,12,15,17), P at base point E2=0.
A = 7E1 + cE2, Omega = -121E1^2 + w1E1E2 + w2E2^2, L3=8E1+m3E2, L4=12E1+m4E2, L5=15E1+m5E2, L6=17E1+m6E2.
Identity 2Phi(A,Om)+L3^6+L4^6+L5^6-L6^6=0: E1^6 coefficient holds; 6 remaining equations in 7 unknowns.
Generic point: no solution. Count solutions mod several primes as a specialness signal, Hensel+reconstruct any."""
import sympy as sp, numpy as np, sys, math
from fractions import Fraction
c,w1,w2,m3,m4,m5,m6,E1,E2=sp.symbols('c w1 w2 m3 m4 m5 m6 E1 E2')
A=7*E1+c*E2; Om=-121*E1**2+w1*E1*E2+w2*E2**2
Phi=A**6+15*A**4*Om+15*A**2*Om**2+Om**3
ident=sp.expand(2*Phi+(8*E1+m3*E2)**6+(12*E1+m4*E2)**6+(15*E1+m5*E2)**6-(17*E1+m6*E2)**6)
P=sp.Poly(ident,E1,E2); eqs=[sp.expand(P.coeff_monomial(E1**(6-k)*E2**k)) for k in range(7)]
assert eqs[0]==0; eqs=eqs[1:]; VARS=[c,w1,w2,m3,m4,m5,m6]
print("linear eq:",eqs[0])
terms=[[(tuple(int(x) for x in mon),int(cf)) for mon,cf in sp.Poly(e,*VARS).terms()] for e in eqs]
def solve_mod_p(p):
    lin=terms[0]; coef=[0]*7; const=0
    for mon,cf in lin:
        if sum(mon)==0: const=cf
        else: coef[mon.index(1)]=cf
    piv=next((i for i in range(7) if coef[i]%p),None)
    if piv is None: return None
    others=[i for i in range(7) if i!=piv]
    g=np.indices((p,)*6).reshape(6,-1).astype(np.int64); V=[None]*7
    for k,i in enumerate(others): V[i]=g[k]
    inv=pow(coef[piv]%p,-1,p); V[piv]=(-(const+sum(coef[i]*V[i] for i in others)))%p*inv%p
    pows=[[np.ones_like(V[0])] for _ in range(7)]
    for i in range(7):
        for e in range(1,7): pows[i].append(pows[i][-1]*V[i]%p)
    mask=np.ones(len(V[0]),dtype=bool)
    for T in terms[1:]:
        acc=np.zeros(len(V[0]),dtype=np.int64)
        for mon,cf in T:
            tt=np.full(len(V[0]),cf%p,dtype=np.int64)
            for i,e in enumerate(mon):
                if e: tt=tt*pows[i][e]%p
            acc=(acc+tt)%p
        mask&=(acc==0)
    return [tuple(int(V[i][j]) for i in range(7)) for j in np.nonzero(mask)[0]]
for p in [int(x) for x in sys.argv[1:] if x.isdigit()]:
    s=solve_mod_p(p); print("p=%d: %s solutions mod p"%(p,"pivot-fail" if s is None else len(s)),flush=True)
    if s and len(s)<=5: print("   ",s)

def classify(sols,p):
    from collections import Counter
    kinds=Counter()
    for (cc,ww1,ww2,mm3,mm4,mm5,mm6) in sols:
        # proportional (constant map): (c,m3,m4,m5,m6) = lam*(7,8,12,15,17), w1=-242 lam, w2=-121 lam^2
        lam=(cc*pow(7,-1,p))%p
        prop = all(((x-lam*y)%p==0) for x,y in ((mm3,8),(mm4,12),(mm5,15),(mm6,17))) and (ww1+242*lam)%p==0 and (ww2+121*lam*lam)%p==0
        disc=(ww1*ww1+4*121*ww2)%p   # disc of Omega = w1^2 - 4*(-121)*w2
        kinds[("constant" if prop else "nonconstant", "disc0" if disc==0 else "disc!=0")]+=1
    return kinds
if __name__=="__main__" and '--classify' in sys.argv:
    for p in (11,13):
        s=solve_mod_p(p); print(p,classify(s,p))
        nc=[x for x in s if not (all(((y-(x[0]*pow(7,-1,p))%p*z)%p==0) for y,z in ((x[3],8),(x[4],12),(x[5],15),(x[6],17))))]
        print("  sample nonconstant:",nc[:4])
