"""Rational points on the Fano fourfold Y: x1^6+x2^6+x3^6+x4^6 + T^3 = x6^6 (T any sign, integer),
primitive, all x_i>0, T != 0 and T not a perfect square (a square T would itself be an ESOP6 solution!)."""
import numpy as np, sys, math, json
N=int(sys.argv[1]) if len(sys.argv)>1 else 40
S6=np.array([x**6 for x in range(N+1)],dtype=object)
Tmax=int(round((N**6)**(1/3)))+2
cubes={t**3:t for t in range(-Tmax,Tmax+1)}
pts=[]
for x6 in range(1,N+1):
    for x1 in range(1,x6):
        for x2 in range(x1,x6):
            R=S6[x6]-S6[x1]-S6[x2]
            for x3 in range(x2,x6):
                R3=R-S6[x3]
                for x4 in range(x3,x6):
                    v=R3-S6[x4]
                    if v==0: continue
                    t=cubes.get(v)
                    if t is not None:
                        g=math.gcd(math.gcd(x1,x2),math.gcd(math.gcd(x3,x4),x6))
                        if g!=1: continue
                        sq = t>0 and math.isqrt(t)**2==t
                        pts.append(dict(x=[x1,x2,x3,x4],x6=x6,T=t,T_is_square=sq))
                        print(pts[-1],flush=True)
print("total",len(pts))
json.dump(pts,open("research/cube_ansatz_2026_09_10/y_points_N%d.json"%N,"w"),indent=0)
