"""Do two Y2 points lie on a common Y2-line (hence on a common plane cubic of X)?
Line: L = y0 p0 + mu y1 p1 (weight-1 coords), S = S0 p0^3 + s1 p0^2 p1 + s2 p0 p1^2 + mu^3 S1 p1^3.
Identity L6^6 - sum L_i^6 = S^2.  With A_k = x6^(6-k) y6^k - sum x_i^(6-k) y_i^k, the mu-dependence cancels and
the pair is collinear iff  15 A2 = 6 S0 A5/S1 + 9 A1^2/S0^2,  20 A3 = 2 S0 S1 + 18 A1 A5/(S0 S1),
15 A4 = 6 A1 S1/S0 + 9 A5^2/S1^2  (exact rational arithmetic).  Signs of S0,S1 (±) both tried."""
import json, sys, itertools
from fractions import Fraction as Fr
pts=json.load(open(sys.argv[1]))
def A(x,x6,y,y6,k): return x6**(6-k)*y6**k-sum(a**(6-k)*b**k for a,b in zip(x,y))
hits=0; tested=0
for P,Q in itertools.combinations(pts,2):
    x,x6,S0=P['x'],P['x6'],P['S']; y,y6,S1=Q['x'],Q['x6'],Q['S']
    for perm in itertools.permutations(range(4)):
        yy=[y[i] for i in perm]
        A1,A2,A3,A4,A5=[A(x,x6,yy,y6,k) for k in (1,2,3,4,5)]
        for e0 in (1,-1):
            for e1 in (1,-1):
                s0=e0*S0; s1_=e1*S1
                tested+=1
                c1= Fr(15*A2)==Fr(6*s0*A5,s1_)+Fr(9*A1*A1,s0*s0)
                if not c1: continue
                c2= Fr(20*A3)==Fr(2*s0*s1_)+Fr(18*A1*A5,s0*s1_)
                c3= Fr(15*A4)==Fr(6*A1*s1_,s0)+Fr(9*A5*A5,s1_*s1_)
                if c2 and c3:
                    hits+=1; print("COLLINEAR PAIR:",P,Q,"perm",perm,"signs",e0,e1,flush=True)
print("pairs tested (with perms/signs):",tested,"collinear:",hits)
