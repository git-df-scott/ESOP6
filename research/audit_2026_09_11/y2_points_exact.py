"""Exact (int64) Y2 seed search: x1^6+x2^6+x3^6+x4^6 + S^2 = x6^6, 1<=xi<x6, S>0, primitive.

Audit replacement for y2_points_fast.py, whose float64 pair sums silently lose
solutions once x6^6 exceeds 2^53 (x6 > 456).  Usage: y2_points_exact.py LO HI OUT.json
"""
import numpy as np, sys, math, json, time

lo=int(sys.argv[1]); hi=int(sys.argv[2]); out=sys.argv[3]
six=np.array([x**6 for x in range(hi+1)],dtype=np.int64)
pts=[]; seen=set()
for x6 in range(lo,hi+1):
    t0=time.time()
    ia,ib=np.triu_indices(x6)
    m=ia>=1; ia=ia[m].astype(np.int64); ib=ib[m].astype(np.int64)
    ps=six[ia]+six[ib]
    X6=int(x6)**6
    for k in range(len(ps)):
        rv=(X6-int(ps[k]))-ps[k:]
        pos=rv>0
        r=np.rint(np.sqrt(np.where(pos,rv,0).astype(np.float64))).astype(np.int64)
        hit=np.nonzero(pos&(r*r==rv))[0]
        for j in hit:
            x=[int(ia[k]),int(ib[k]),int(ia[k+j]),int(ib[k+j])]
            xs=sum(v**6 for v in x)
            target=int(x6)**6
            if target<=xs: continue
            S=math.isqrt(target-xs)
            if xs+S*S!=target or S==0: continue
            g=math.gcd(math.gcd(x[0],x[1]),math.gcd(math.gcd(x[2],x[3]),x6))
            if g!=1: continue
            key=(tuple(sorted(x)),x6)
            if key in seen: continue
            seen.add(key)
            pts.append(dict(x=sorted(x),x6=x6,S=S,S_is_cube=(round(S**(1/3))**3==S)))
    print("x6=%d cumulative %d  (%.1fs)"%(x6,len(pts),time.time()-t0),flush=True)
    json.dump(pts,open(out,"w"),indent=0)
print("total",len(pts))
