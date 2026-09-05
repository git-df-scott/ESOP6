"""Combine existing 7-adic prefixes with explicit local axis points elsewhere.

This is a control against mistaking simultaneous congruences for rational curves.
No import of the old lifting script (it has mutating top-level code).
"""
import json
import math
from fractions import Fraction
from pathlib import Path
from multiplicity_algebra import mul, power

OUT=Path('results/multiplicity_square_2026_09_05')

def residual(x):
    polys=[x[3*i:3*i+3] for i in range(4)]
    A=[x[12],x[13],1];B=x[14:17]
    P6=[a+7**6*b for a,b in zip(A,B)]
    left=[sum(power(P,6)[j]*7**6 for P in polys)+power(A,6)[j] for j in range(13)]
    return [a-b for a,b in zip(left,power(P6,6))]

def reconstruct(a,M):
    H=math.isqrt(M//2);a%=M
    if not a:return Fraction(0)
    r0,r1,s0,s1=M,a,0,1
    while r1>H:
        q=r0//r1;r0,r1=r1,r0-q*r1;s0,s1=s1,s0-q*s1
    if not s1 or abs(s1)>H or math.gcd(r1,s1)!=1 or math.gcd(s1,M)!=1:return None
    f=Fraction(r1,s1)
    return f if (f.numerator-a*f.denominator)%M==0 else None

rows=[];M7=7**41;other=360;M=M7*other
axis=[0]*12+[1,0]+[0]*3
for old in json.loads(Path('results/fourfold_routes_2026_09_05/conic_lifts.json').read_text()):
    y=[int(v) for v in old['full_conic_final_vector']]
    x=[a+M7*((b-a)*pow(M7,-1,other)%other) for a,b in zip(y,axis)]
    rr=residual(x)
    for m in (8,9,5,7**47):assert all(v%m==0 for v in rr)
    assert any(rr)
    fs=[reconstruct(a,M) for a in x]
    exact=all(f is not None for f in fs) and not any(residual(fs))
    assert not exact
    rows.append({'pattern':old['pattern'],'coefficient_modulus':str(M),
                 'original_equation_moduli':[8,9,5,str(7**47)],
                 'combined_coefficients':[str(a) for a in x],
                 'reconstructed_coefficients':sum(f is not None for f in fs),
                 'rational_reconstruction':[str(f) if f is not None else None for f in fs],
                 'exact_curve':False})
report={'rows':rows,'scope':'Four explicit simultaneous local prefixes. Uses trivial axis reductions at 2,3,5. Neither excludes other rational reconstructions nor establishes a global curve.'}
(OUT/'conic_crt.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({r['pattern']:r['reconstructed_coefficients'] for r in rows}))
