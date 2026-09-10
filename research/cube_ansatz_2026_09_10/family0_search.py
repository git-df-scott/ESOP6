"""Family 0 (weighted lines on the Fano threefold  x1^6+x2^6+x3^6+x4^6 + T^3 = x6^6 ; conics on X via x5^2=T):
   (d E1)^6 - (d E2)^6 - sum_{j=2..4} (n_j E1 + m_j E2)^6  =  C(E1,E2)^3   with C a binary quadratic over Q.
Here x6 = d E1, x1 = d E2, x_j = n_j E1 + m_j E2 (j=2,3,4), and x5 = w with w^2 = C/d^2 ... (F = C/d^2 after scaling).
Exact integer search with endpoint cube filters."""
import itertools, sys, json
from math import gcd
H=int(sys.argv[1]) if len(sys.argv)>1 else 8
DMAX=int(sys.argv[2]) if len(sys.argv)>2 else 6
def icbrt(n):
    if n==0: return 0
    s=-1 if n<0 else 1; n=abs(n)
    r=round(n**(1/3))
    for c in (r-1,r,r+1):
        if c**3==n: return s*c
    return None
def binom6(): return [1,6,15,20,15,6,1]
B=binom6()
hits=[]
for d in range(1,DMAX+1):
    d6=d**6
    # filter A: (n2,n3,n4) unordered multiset with d^6 - sum n^6 a cube
    A=[]
    for n2 in range(0,H+1):
        for n3 in range(n2,H+1):
            for n4 in range(n3,H+1):
                c=icbrt(d6-n2**6-n3**6-n4**6)
                if c is not None: A.append(((n2,n3,n4),c))
    # filter Bm: (m2,m3,m4) with -(d^6 + sum m^6) a cube  <=> d^6+sum m^6 is a cube
    Bm=[]
    for m2 in range(0,H+1):
        for m3 in range(0,H+1):
            for m4 in range(0,H+1):
                c=icbrt(d6+m2**6+m3**6+m4**6)
                if c is not None: Bm.append(((m2,m3,m4),-c))
    print("d=%d: A-filter %d  B-filter %d"%(d,len(A),len(Bm)),flush=True)
    for (ns,alpha) in A:
        # ns are absolute values (signs irrelevant for 6th powers but matter for cross terms); assign signs & pairings later
        for perm in set(itertools.permutations(ns)):
            for (ms,gamma) in Bm:
                for sg in itertools.product((1,-1),repeat=3):
                    m=[s*x for s,x in zip(sg,ms)]
                    n=list(perm)
                    # sextic coefficients: coef[k] of E1^(6-k) E2^k
                    coef=[0]*7
                    coef[0]+=d6; coef[6]-=d6
                    for nj,mj in zip(n,m):
                        for k in range(7):
                            coef[k]-=B[k]*nj**(6-k)*mj**k
                    # perfect cube of alpha E1^2 + beta E1E2 + gamma E2^2 ?  coef[1] = 3 alpha^2 beta
                    if alpha==0 or gamma==0: continue
                    if coef[1] % (3*alpha*alpha): continue
                    beta=coef[1]//(3*alpha*alpha)
                    cube=[alpha**3, 3*alpha**2*beta, 3*alpha*beta**2+3*alpha**2*gamma, beta**3+6*alpha*beta*gamma, 3*beta**2*gamma+3*alpha*gamma**2, 3*beta*gamma**2, gamma**3]
                    if cube==coef:
                        hits.append(dict(d=d,n=n,m=m,C=[alpha,beta,gamma]))
                        print("HIT",hits[-1],flush=True)
print("total hits",len(hits))
json.dump(hits,open("/tmp/claude-0/-home-user/15becd99-5cb2-5637-beeb-15d2ba0dfed3/scratchpad/family0_hits.json","w"))
