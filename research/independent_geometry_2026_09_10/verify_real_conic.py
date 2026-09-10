#!/usr/bin/env python3
import sympy as s
from pathlib import Path
import json
u,v=s.symbols('u v')
base=s.expand(u**6+v**6+((u+v)/s.sqrt(2))**6+((u-v)/s.sqrt(2))**6-s.Rational(5,4)*(u*u+v*v)**3)
assert base==0
r,t=s.symbols('r t'); U=r*r-t*t;V=2*r*t;W=r*r+t*t
assert s.expand(U*U+V*V-W*W)==0
assert s.Rational(63,80)*s.Rational(5,4)+s.Rational(1,2)**6==1
result={'four_form_identity_residual':str(base),'pythagorean_identity_verified':True,'normalization_verified':True,'rational_counterexample':False,'rationality_obstruction':'x3+x4=sqrt(2)*x1 and x3-x4=sqrt(2)*x2 force first four rational homogeneous coordinates zero; x5=W/2 then forces W=0.'}
Path(__file__).with_name('real_conic_verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
