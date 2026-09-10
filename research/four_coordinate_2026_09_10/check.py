"""Exact checks accompanying the four-coordinate denominator obstruction."""
import itertools
import json
from pathlib import Path
import sympy as S

t, M, c, x, y, z, v, u, v1, v2, v3, v4 = S.symbols('t M c x y z v u v1 v2 v3 v4')
rhs = 6*c*M**5 + 5*c**3*t**12*M**3 + S.Rational(3,8)*c**5*t**24*M
assert S.expand((M+c*t**6/2)**6-(M-c*t**6/2)**6-t**6*rhs) == 0
assert sum(a**6 for a in (1,2,4,3)) == 6*815

# First anisotropic quadratic in the mod-4 reduction.
q = x*x+y*y+z*z+x*y+x*z+y*z
norm = (x+z)**2+(x+z)*(y+z)+(y+z)**2
assert all(int(a)%2 == 0 for a in S.Poly(q-norm,x,y,z).coeffs())
# After four entries agree mod 2, the mod-8 obstruction is the same norm.
vv = [v1,v2,v3,v4]
quotient = S.expand(sum((u+2*a)**2 for a in vv)/4)
V = sum(vv)
assert all(int(a)%2 == 0 for a in S.Poly(quotient-(u*u+u*V+V*V),u,*vv).coeffs())

# Exhaustive cross-check: all four degree-at-most-one polynomials modulo 8.
# Meet-in-middle tests all 64^4 ordered quadruples without enumerating them.
def mul(a,b):
    out=[0]*(len(a)+len(b)-1)
    for i,aa in enumerate(a):
        for j,bb in enumerate(b): out[i+j]=(out[i+j]+aa*bb)%8
    return tuple(out)
polys=[]
for aa,bb in itertools.product(range(8),repeat=2):
    p=(aa,bb); sixth=(1,)
    for _ in range(6): sixth=mul(sixth,p)
    polys.append((sixth, bool(aa%2 or bb%2)))
pairs={}
for a,ap in polys:
    for b,bp in polys:
        key=tuple((aa+bb)%8 for aa,bb in zip(a,b))
        counts=pairs.setdefault(key,[0,0])
        counts[int(ap or bp)]+=1
primitive=0; nonprimitive=0
for key,counts in pairs.items():
    opposite=pairs.get(tuple((-a)%8 for a in key),[0,0])
    nonprimitive+=counts[0]*opposite[0]
    primitive+=counts[1]*sum(opposite)+counts[0]*opposite[1]
assert primitive == 0
assert nonprimitive == 16**4
result={
    'counterexample_found':False,
    'midpoint_identity_verified':True,
    'initial_sixth_sum':4890,
    'c':815,
    'four_square_reduction_identities_verified':True,
    'mod8_degree_at_most_one_ordered_quadruples_covered':64**4,
    'primitive_sixth_power_zero_sums':primitive,
    'nonprimitive_sixth_power_zero_sums':nonprimitive,
    'scope':'The enumeration cross-checks a lemma; the accompanying algebra proves it for arbitrary degree.'
}
Path(__file__).with_name('checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
