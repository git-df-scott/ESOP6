"""Exact algebra supporting the written proofs, not a substitute for them."""
from collections import defaultdict
from pathlib import Path
import json,sympy as S
out=Path('results/fourfold_routes_2026_09_05')
def v3(n):
 k=0
 while n%3==0:n//=3;k+=1
 return k
pairs=defaultdict(list)
for a in range(1,31):
 for b in range(a,31):pairs[a*a+b*b].append((a,b))
checks=0
for norm,ps in pairs.items():
 for a,b in ps:
  for c,d in ps:
   k=min(v3(x) for x in (a,b,c,d))
   assert v3(norm)==2*k
   assert v3(a**6+b**6+c**6+d**6)==6*k
   checks+=1
units=0
for a in range(1,100):
 for b in range(1,100):
  if a%3 and b%3:
   assert v3(a**4+a*a*b*b+b**4)==1;units+=1
s=S.symbols('s');E=2*(243601295510403*s**12-722265625)/S.Integer(3900234375)
G=2*(3599809104248613*s**12-2744609375)/S.Integer(1300078125)
AE=S.expand(E).coeff(s,12);AG=S.expand(G).coeff(s,12)
contradiction=S.factor(AG*E-AE*G)
assert contradiction==S.Rational(-220092525210624,144453125)
r={'equal_pair_valuation_controls':checks,'unit_factor_valuation_controls':units,
 'sparse_exact_nonzero_combination':str(contradiction),'status':'PASS'}
(out/'proof_checks.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
