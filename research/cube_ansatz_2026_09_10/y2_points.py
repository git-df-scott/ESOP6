"""Rational points on the Fano fourfold Y2: x1^6+x2^6+x3^6+x4^6 + S^2 = x6^6 (S>0 integer; S = x5^3 would be ESOP6).
Primitive, all x_i>0.  Meet-in-the-middle on pair sums."""
import numpy as np, sys, math, json
N=int(sys.argv[1]) if len(sys.argv)>1 else 60
S6=[x**6 for x in range(N+1)]
pts=[]
for x6 in range(2,N+1):
    target=S6[x6]
    # pairs (x1<=x2) and (x3<=x4) with all < x6 ; need target - (a+b) = square >0
    pairs=[(S6[a]+S6[b],a,b) for a in range(1,x6) for b in range(a,x6)]
    for (sa,a,b) in pairs:
        for (sb,c,d) in pairs:
            if c<a or (c==a and d<b): continue
            v=target-sa-sb
            if v<=0: continue
            r=math.isqrt(v)
            if r*r==v:
                g=math.gcd(math.gcd(a,b),math.gcd(math.gcd(c,d),x6))
                if g!=1: continue
                cube = round(r**(1/3))**3==r
                pts.append(dict(x=[a,b,c,d],x6=x6,S=r,S_is_cube=cube))
                print(pts[-1],flush=True)
print("total",len(pts))
json.dump(pts,open("research/cube_ansatz_2026_09_10/y2_points_N%d.json"%N,"w"),indent=0)
