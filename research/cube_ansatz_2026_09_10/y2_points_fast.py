"""Vectorised Y2 point search: x1^6+x2^6+x3^6+x4^6 + S^2 = x6^6, x_i>=1, primitive."""
import numpy as np, sys, math, json
N=int(sys.argv[1]) if len(sys.argv)>1 else 150
six=np.array([x**6 for x in range(N+1)],dtype=np.float64)
pts=[]; seen=set()
for x6 in range(2,N+1):
    ia,ib=np.triu_indices(x6); ia=ia[ia>0]; ib=ib[-len(ia):] if False else None
    ia,ib=np.triu_indices(x6)
    m=ia>=1; ia=ia[m]; ib=ib[m]
    ps=six[ia]+six[ib]
    X6=six[x6]
    for k in range(len(ps)):
        a=ps[k]
        rest=X6-a-ps[k:]            # require pair index >= k to avoid double counting
        ok=rest>0
        r=np.floor(np.sqrt(rest[ok])+0.5)
        hit=np.nonzero(r*r==rest[ok])[0]
        if len(hit)==0: continue
        idx=np.nonzero(ok)[0][hit]
        for j in idx:
            x=[int(ia[k]),int(ib[k]),int(ia[k+j]),int(ib[k+j])]
            S=int(r[np.nonzero(ok)[0]==k+j][0]) if False else int(round(math.sqrt(int(X6)-int(a)-int(ps[k+j]))))
            xs=sum(v**6 for v in x)
            if xs+S*S!=int(X6) or S==0: continue
            g=math.gcd(math.gcd(x[0],x[1]),math.gcd(math.gcd(x[2],x[3]),x6))
            if g!=1: continue
            key=(tuple(sorted(x)),x6)
            if key in seen: continue
            seen.add(key)
            pts.append(dict(x=sorted(x),x6=x6,S=S,S_is_cube=(round(S**(1/3))**3==S)))
    if x6%25==0: print("x6=%d points so far %d"%(x6,len(pts)),flush=True)
print("total",len(pts))
json.dump(pts,open("research/cube_ansatz_2026_09_10/y2_points_N%d.json"%N,"w"),indent=0)
