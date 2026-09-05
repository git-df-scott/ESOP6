"""Exhaustive first transverse residue calculation at an irreducible F7 base point."""
import itertools,json
from collections import defaultdict
from pathlib import Path
import sympy as S
t=S.symbols('t');modulus=(t*t+1)**5
forms=[(0,0,0)]
for d in range(3):
    forms.extend(tuple(co)+(1,)+(0,)*(2-d) for co in itertools.product(range(7),repeat=d))
assert len(forms)==58
vec=[]
for co in forms:
    f=sum(co[i]*t**i for i in range(3))
    rem=S.Poly(S.rem(f**6,modulus,t),t,modulus=7)
    vec.append(tuple(int(rem.nth(i))%7 for i in range(10)))
table=defaultdict(list)
for i in range(58):
    for j in range(i,58):
        table[tuple((a+b)%7 for a,b in zip(vec[i],vec[j]))].append((i,j))
sol=set()
for key,pairs in table.items():
    for i,j in pairs:
        for k,l in table.get(tuple(-v%7 for v in key),[]):
            sol.add(tuple(sorted((i,j,k,l))))
qidx=forms.index((1,0,1))
nontrivial=[x for x in sorted(sol) if any(i not in (0,qidx) for i in x)]
report={'field':7,'divisor':'(t^2+1)^5','projective_quadratic_forms_including_zero':len(forms),
 'pair_entries':sum(map(len,table.values())),'solution_multisets':len(sol),
 'not_all_divisible_by_t2_plus_1':len(nontrivial),
 'nontrivial_solutions':[[forms[i] for i in x] for x in nontrivial],
 'scope':'complete residue problem sum of four quadratic sixth powers divisible by Q^5 over F7; not a Q-conic classification'}
Path('results/fourfold_routes_2026_09_05/conic_F7_transverse.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
