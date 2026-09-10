"""V3 lottery: genus-1 quadratic-cover curves (w^2 = F(t), F quartic) through rational points of Y."""
import json, sys, numpy as np
from fractions import Fraction
sys.path.insert(0,'research/cube_ansatz_2026_09_10')
from conics_through_point import System
def ratrec(z,maxden=10**6):
    if abs(z.imag)>1e-7*(1+abs(z)): return None
    f=Fraction(z.real).limit_denominator(maxden)
    return f if abs(float(f)-z.real)<1e-8*(1+abs(z.real)) else None
pts=json.load(open(sys.argv[1])); starts=int(sys.argv[2]) if len(sys.argv)>2 else 1500
for p in pts:
    y0=list(p['x'])+[p['x6'],p['T']]
    S=System(y0); sols=S.solve_all(starts=starts,scale_range=(0,3))
    nreal=sum(1 for v in sols if np.abs(v.imag).max()<1e-6)
    rat=[v for v in sols if all(ratrec(z) is not None for z in v)]
    print("Y-point",y0,": curves found",len(sols),"near-real",nreal,"rational-looking",len(rat),flush=True)
    for v in rat: print("  !!!! CANDIDATE:",[str(ratrec(z)) for z in v],flush=True)
    np.save("research/cube_ansatz_2026_09_10/v3_sols_%d.npy"%p['x6'],np.array(sols))
