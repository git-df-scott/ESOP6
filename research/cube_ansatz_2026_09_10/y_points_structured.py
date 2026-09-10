"""Rational points on Y: x1^6+x2^6+x3^6+x4^6 + T^3 = x6^6 (T != 0 any sign), exploiting that at most two of
x1..x4 are coprime to 7 (sixth powers of units are 1 mod 7, cubes are 0,±1).  Meet-in-the-middle: the two
7-divisible coordinates (c,d) are enumerated, the free pair (a,b) looked up in a sorted table of pair sums."""
import numpy as np, sys, math, json
N=int(sys.argv[1]) if len(sys.argv)>1 else 120
six=np.array([x**6 for x in range(N+1)],dtype=np.float64)  # use float for speed, verify exactly after
ia,ib=np.triu_indices(N+1)
pairsum=six[ia]+six[ib]
order=np.argsort(pairsum); ps=pairsum[order]; pa=ia[order]; pb=ib[order]
mult=np.arange(7,N+1,7)
ic,idd=np.triu_indices(len(mult))
c=mult[ic]; d=mult[idd]; cd=six[c]+six[d]
pts=[]
for x6 in range(2,N+1):
    X6=x6**6
    Tm=int(round(X6**(1/3)))+1
    T=np.arange(-Tm,Tm+1); T=T[T!=0]
    targets=X6-T.astype(np.float64)**3
    need=targets[:,None]-cd[None,:]          # (nT, npairs)
    flat=need.ravel()
    pos=np.searchsorted(ps,flat)
    pos=np.clip(pos,0,len(ps)-1)
    hit=np.abs(ps[pos]-flat)<0.5
    for k in np.nonzero(hit)[0]:
        ti,pi=divmod(k,len(cd))
        a,b=int(pa[pos[k]]),int(pb[pos[k]]); cc,dd=int(c[pi]),int(d[pi]); t=int(T[ti])
        if a==0 or cc==0: continue
        if a**6+b**6+cc**6+dd**6+t**3!=X6: continue
        g=math.gcd(math.gcd(a,b),math.gcd(math.gcd(cc,dd),x6))
        if g!=1: continue
        sq=t>0 and math.isqrt(t)**2==t
        pts.append(dict(x=sorted([a,b,cc,dd]),x6=x6,T=t,T_is_square=sq))
        print(pts[-1],flush=True)
print("total",len(pts))
json.dump(pts,open("research/cube_ansatz_2026_09_10/y_points_structured_N%d.json"%N,"w"),indent=0)
