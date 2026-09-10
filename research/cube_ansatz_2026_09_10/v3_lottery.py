"""V3 lottery: genus-1 quadratic-cover curves (w^2 = F(t), F quartic) through rational points of Y."""
import json, sys, numpy as np
from fractions import Fraction
sys.path.insert(0,'research/cube_ansatz_2026_09_10')
from conics_through_point import System
def ratrec(z,maxden=3000):
    if abs(z.imag)>1e-7*(1+abs(z)): return None
    f=Fraction(z.real).limit_denominator(maxden)
    return f if abs(float(f)-z.real)<1e-10*(1+abs(z.real)) else None
import sympy as sp
def exact_check(y0,fr):
    t=sp.symbols('t'); x=y0[:4]; x6=y0[4]; T=y0[5]; g=fr[:8]; f=fr[8:]
    G=[x[0]+t+g[0]*t**2, x[1]+g[1]*t+g[2]*t**2, x[2]+g[3]*t+g[4]*t**2, x[3]+g[5]*t+g[6]*t**2]
    G6=x6+g[7]*t**2; F=T+f[0]*t+f[1]*t**2+f[2]*t**3+f[3]*t**4
    return sp.expand(sum(Gi**6 for Gi in G)+F**3-G6**6)==0
pts=json.load(open(sys.argv[1])); starts=int(sys.argv[2]) if len(sys.argv)>2 else 1500
for p in pts:
    y0=list(p['x'])+[p['x6'],p['T']]
    S=System(y0); sols=S.solve_all(starts=starts,scale_range=(0,3),real=True)
    nreal=sum(1 for v in sols if np.abs(v.imag).max()<1e-6)
    rat=[v for v in sols if all(ratrec(z) is not None for z in v)]
    print("Y-point",y0,": curves found",len(sols),"near-real",nreal,"rational-looking",len(rat),flush=True)
    for v in rat:
        fr=[ratrec(z) for z in v]
        if exact_check(y0,fr): print("  !!!! EXACT RATIONAL GENUS-1 CURVE:",[str(q) for q in fr],flush=True)
    np.save("research/cube_ansatz_2026_09_10/v3_sols_%d.npy"%p['x6'],np.array(sols))
