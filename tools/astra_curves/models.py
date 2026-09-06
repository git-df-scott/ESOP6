#!/usr/bin/env python3
"""Export exact plane-cubic and elliptic-quartic incidence ideals."""
import json,itertools
from pathlib import Path
import sympy as S
OUT=Path('results/astra_curves_2026_09_05')
def mons(n,d):
 if n==1:return [(d,)]
 return [(i,)+a for i in range(d,-1,-1) for a in mons(n-1,d-i)]
def mon(z,e):return S.prod(v**k for v,k in zip(z,e))
def plane_cubics():
 z=S.symbols('u v w');x=S.symbols('x0:28');ls=[sum(x[3*i+k]*z[k] for k in range(3)) for i in range(3)];m3=mons(3,3);pivot=m3.index((0,0,3));free=[i for i in range(10) if i!=pivot]
 c1=z[2]**3+sum(x[9+j]*mon(z,m3[k]) for j,k in enumerate(free));c2=sum(x[18+j]*mon(z,e) for j,e in enumerate(m3))
 pol=S.Poly(z[0]**6+z[1]**6+sum(l**6 for l in ls)-z[2]**6-c1*c2,*z)
 eq=[str(pol.coeff_monomial(mon(z,e))) for e in mons(3,6)]
 return {'variables':[str(v) for v in x],'unknowns':28,'equation_count':28,'coordinate_forms':[str(z[0]),str(z[1])]+[str(l) for l in ls]+[str(z[2])],'cubic_1':str(c1),'cubic_2':str(c2),'equations':eq,'coverage':'one plane chart and one cubic scale chart; all real plane-cubic curves excluded by the separate odd-degree proof','complex_solution_enumeration_performed':False}
def elliptic():
 z=S.symbols('z0:4');x=S.symbols('x0:84');l=sum(x[i]*z[i] for i in range(4));m=sum(x[4+i]*z[i] for i in range(4));m2=mons(4,2);m4=mons(4,4);free=[i for i,e in enumerate(m2) if e not in [(2,0,0,0),(0,2,0,0)]];bf=[i for i,e in enumerate(m4) if e[0]<2]
 q1=z[0]**2+sum(x[8+j]*mon(z,m2[k]) for j,k in enumerate(free));q2=z[1]**2+sum(x[16+j]*mon(z,m2[k]) for j,k in enumerate(free));a=sum(x[24+j]*mon(z,e) for j,e in enumerate(m4));b=sum(x[59+j]*mon(z,m4[k]) for j,k in enumerate(bf))
 pol=S.Poly(sum(v**6 for v in z)+l**6-m**6-q1*a-q2*b,*z)
 eq=[str(pol.coeff_monomial(mon(z,e))) for e in mons(4,6)]
 return {'variables':[str(v) for v in x],'unknowns':84,'equation_count':84,'coordinate_forms':[str(v) for v in z]+[str(l),str(m)],'quadric_1':str(q1),'quadric_2':str(q2),'cofactor_A':str(a),'cofactor_B':str(b),'equations':eq,'gauge':'quadric pencil RREF at z0^2,z1^2; B has no monomial divisible by z0^2, removing the 10-dimensional degree-two cofactor syzygy','conditions':'quadrics must meet in a smooth curve (pencil discriminant nonzero); all Grassmann and pencil charts needed for global completeness','parameter_dimension_count':{'P3':8,'pencil':16,'cofactors_after_syzygy':60}}
if __name__=='__main__':
 for name,fn in [('plane_cubic',plane_cubics),('elliptic',elliptic)]:
  d=fn();assert len(d['equations'])==d['equation_count'];(OUT/(name+'_ideal.json')).write_text(json.dumps(d,indent=2)+'\n');print(name,d['unknowns'],d['equation_count'])
