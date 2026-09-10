"""Last mile: given a rational plane cubic on X from plane_cubics_through_seed.py, i.e. linear forms
L_i = x_i p0 + m_i p1 (i=1..4), L6 = x6 p0 + p1 and S = s0 p0^3 + s1 p0^2 p1 + s2 p0 p1^2 + s3 p1^3 with
L6^6 - sum L_i^6 = S^2, search coprime integer (p0,p1) with S(p0,p1) a perfect cube u^3, u != 0.
Then (L1,L2,L3,L4,u,L6) * common denominator is an ESOP6 solution; verify with tools/verify_esop6.py."""
import sys, json, math
from fractions import Fraction
def icbrt(n):
    if n==0: return 0
    s=-1 if n<0 else 1; n=abs(n); r=round(n**(1/3))
    for c in (r-1,r,r+1):
        if c**3==n: return s*c
    return None
def search(x,x6,s0,m,s,H):
    m2,m3,m4=m; s1,s2,s3=s
    den=1
    for v in list(m)+list(s): den=den*v.denominator//math.gcd(den,v.denominator)
    hits=[]
    for p0 in range(-H,H+1):
        for p1 in range(1,H+1):
            if math.gcd(p0,p1)!=1: continue
            S=s0*p0**3+s1*p0**2*p1+s2*p0*p1**2+s3*p1**3   # Fraction
            Sd=S*den**3
            assert Sd.denominator==1
            u=icbrt(int(Sd))
            if u is None or u==0: continue
            L=[Fraction(x[0]*p0)*den, (x[1]*p0+m2*p1)*den, (x[2]*p0+m3*p1)*den, (x[3]*p0+m4*p1)*den]
            L6=Fraction(x6*p0+p1)*den
            vals=[int(v) for v in L]+[u,int(L6)]
            if any(v==0 for v in vals): continue
            vals=[abs(v) for v in vals]
            assert sum(v**6 for v in vals[:5])==vals[5]**6
            hits.append(vals)
    return hits
if __name__=="__main__":
    hits=json.load(open(sys.argv[1])); H=int(sys.argv[2]) if len(sys.argv)>2 else 200
    for h in hits:
        pt=h['point']; sol=[Fraction(v) for v in h['sol']]
        res=search(pt['x'],pt['x6'],pt['S'],sol[:3],sol[3:],H)
        print("curve through",pt,": ESOP6 solutions found:",res[:5])
