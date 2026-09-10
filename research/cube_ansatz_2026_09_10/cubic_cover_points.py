"""Last mile for a rational cubic-cover curve through a Y2 seed:
   G_i(t) = x_i + g_i1 t + g_i2 t^2 (g11=1, g61=0), K(t) = s0 + k1 t + k2 t^2 + k3 t^3, M(t) = 1 + m1 t,
   x5 = u*M(t) with u^3 = K(t).  Search t = p/q with K_hom(p,q) = q^3 K(p/q) a perfect cube; clear denominators;
   verify the sextuple exactly."""
import sys, json, math
from fractions import Fraction as Fr
def icbrt(n):
    if n==0: return 0
    s=-1 if n<0 else 1; n=abs(n); r=round(n**(1/3))
    for c in (r-1,r,r+1):
        if c**3==n: return s*c
    return None
def search(x,x6,s0,sol,H):
    g=[Fr(v) for v in sol[:8]]; k=[Fr(v) for v in sol[8:11]]; m1=Fr(sol[11])
    den=1
    for v in g+k+[m1]: den=den*v.denominator//math.gcd(den,v.denominator)
    hits=[]
    for q in range(1,H+1):
        for p in range(-H,H+1):
            if math.gcd(p,q)!=1: continue
            t=Fr(p,q)
            Kt=s0+k[0]*t+k[1]*t**2+k[2]*t**3
            # u^3 = Kt ; write Kt = a/b, need a*b^2 a cube
            a,b=Kt.numerator,Kt.denominator
            c=icbrt(a*b*b)
            if c is None or c==0: continue
            u=Fr(c,b)
            G=[x[0]+t+g[0]*t**2, x[1]+g[1]*t+g[2]*t**2, x[2]+g[3]*t+g[4]*t**2, x[3]+g[5]*t+g[6]*t**2]
            G6=x6+g[7]*t**2; X5=u*(1+m1*t)
            vals=G+[X5,G6]
            D=1
            for v in vals: D=D*v.denominator//math.gcd(D,v.denominator)
            ints=[abs(int(v*D)) for v in vals]
            if any(v==0 for v in ints): continue
            assert sum(v**6 for v in ints[:5])==ints[5]**6
            hits.append(ints)
    return hits
if __name__=="__main__":
    hits=json.load(open(sys.argv[1])); H=int(sys.argv[2]) if len(sys.argv)>2 else 300
    for h in hits:
        pt=h['point']
        # parse the printed forms back? simpler: expect 'sol' list of 12 rationals if present
        if 'sol' not in h: print("no sol field; parse G/K/M from strings manually"); continue
        res=search(pt['x'],pt['x6'],pt['S'],h['sol'],H)
        print("curve through",pt,": ESOP6 solutions:",res[:5])
