#!/usr/bin/env python3
"""Reproduce the eight bounded conic seed requests used in this strike.

A returned seed is verified by exact substitution. A solver report of no
seed is not promoted to an independently certified nonexistence result.
"""
import json
import time
from sympy import symbols
from sympy.solvers.diophantine.diophantine import diop_ternary_quadratic

A,B,D=symbols('A B D',integer=True)
rows=[]
for den in [17,19,23,25,29,31,37,41]:
    p,q=864,49*den
    eq=(2*q**3-p**3)*A*A+6*q**3*B*B-3*p*q*q*D*D
    start=time.monotonic()
    point=diop_ternary_quadratic(eq)
    if any(v is None for v in point):
        row={'lambda':[p,q], 'solver_result':
             'no rational conic point reported; not an independent obstruction certificate'}
    else:
        a,b,d=map(int,point)
        assert (2*q**3-p**3)*a*a+6*q**3*b*b-3*p*q*q*d*d==0
        row={'lambda':[p,q],'point':[a,b,d]}
    row['seconds']=time.monotonic()-start
    rows.append(row)
print(json.dumps(rows,indent=2))
