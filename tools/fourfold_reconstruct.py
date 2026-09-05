"""Bounded rational reconstruction is accepted only after exact substitution."""
import json,math
from fractions import Fraction
from pathlib import Path
from fourfold_conic_lift import residual
out=Path('results/fourfold_routes_2026_09_05')
M=7**41;H=math.isqrt(M//2)
def rr(a):
 a%=M
 if a==0:return Fraction(0)
 r0,r1=M,a;s0,s1=0,1
 while r1>H:
  q=r0//r1;r0,r1=r1,r0-q*r1;s0,s1=s1,s0-q*s1
 if not s1 or abs(s1)>H or math.gcd(r1,s1)!=1 or math.gcd(s1,M)!=1:return None
 f=Fraction(r1,s1)
 return f if (f.numerator-a*f.denominator)%M==0 else None
records=[]
for rec in json.loads((out/'conic_lifts.json').read_text()):
 fs=[rr(int(x)) for x in rec['full_conic_final_vector']]
 exact=all(x is not None for x in fs) and all(x==0 for x in residual(fs,7**6))
 records.append({'pattern':rec['pattern'],'modulus':str(M),'numerator_denominator_bound':str(H),
  'reconstructed_coefficients':sum(x is not None for x in fs),'coefficients':[str(x) if x is not None else None for x in fs],
  'exact_curve':exact})
(out/'rational_reconstruction.json').write_text(json.dumps(records,indent=2)+'\n')
print(json.dumps(records,indent=2))
