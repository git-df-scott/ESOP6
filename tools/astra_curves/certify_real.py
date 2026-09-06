#!/usr/bin/env python3
"""Exact rational contraction certificates for REAL ALGEBRAIC quartic roots.
These are not Q[t] identities and are not ESOP6 counterexample certificates.
"""
from fractions import Fraction as F
from pathlib import Path
import json,itertools
import sympy as S
import mpmath as mp
import numpy as np
from search import RationalQuartic
OUT=Path('results/astra_curves_2026_09_05');mp.mp.dps=100

def conv(p,q):
 r=[F(0)]*(len(p)+len(q)-1)
 for i,a in enumerate(p):
  for j,b in enumerate(q):r[i+j]+=a*b
 return r
def power(p,k):
 r=[F(1)]
 for _ in range(k):r=conv(r,p)
 return r
def ps_and_j(sys,x):
 A=sys.A.astype(int);base=sys.base.astype(int)
 ps=[[F(int(base[i,k]))+sum((int(A[i,k,j])*x[j] for j in range(sys.n)),F(0)) for k in range(5)] for i in range(6)]
 val=[F(0)]*25;jac=[[F(0)]*sys.n for _ in range(25)]
 for i,p in enumerate(ps):
  pp=power(p,5);six=conv(pp,p);sgn=1 if i<5 else -1
  for k,c in enumerate(six):val[k]+=sgn*c
  for k in range(5):
   for j in range(sys.n):
    aa=6*sgn*int(A[i,k,j])
    if aa:
     for h,c in enumerate(pp):jac[k+h][j]+=aa*c
 return ps,[val[k] for k in sys.rows],[[jac[k][j] for j in range(sys.n)] for k in sys.rows]
def mul(A,B):return [[sum((a*b for a,b in zip(row,col)),F(0)) for col in zip(*B)] for row in A]
def inf(A):return max(sum(abs(v) for v in row) for row in A)
def matvec(A,b):return [sum((u*v for u,v in zip(row,b)),F(0)) for row in A]
def frac(v):return F(str(v))
def enc(A):return [[str(z) for z in row] for row in A]
def approx_inv(J,digits=55):
 M=mp.matrix([[mp.mpf(z.numerator)/z.denominator for z in row] for row in J]);B=M**-1
 return [[F(mp.nstr(B[i,j],digits)) for j in range(B.cols)] for i in range(B.rows)]
def verify_certificate(c):
 sys=RationalQuartic(c['lane']=='double');x=[F(v) for v in c['center']];active=c['active'];B=[[F(v) for v in row] for row in c['inverse']];r=F(c['radius']);n=len(active)
 ps,val,JJ=ps_and_j(sys,x);J=[[row[j] for j in active] for row in JJ];BJ=mul(B,J)
 error=inf([[(F(int(i==j))-BJ[i][j]) for j in range(n)] for i in range(n)])
 eta=max(abs(z) for z in matvec(B,val));bnorm=inf(B)
 # Each coordinate-polynomial coefficient has at most one active parameter
 # in these charts, but the deliberately looser n*r bound is used.
 M=max(abs(z) for p in ps for z in p)+n*r
 H=720*(5*M)**4
 K=error+bnorm*n*n*H*r
 assert eta<=r/2 and K< F(1,2)
 # p6 is even with three strictly positive coefficients throughout the box.
 assert all(ps[5][k]>n*r for k in (0,2,4))
 # Rank-five polynomial coefficient matrix implies basepoint-free quartic embedding.
 rows=c['rank_minor_rows'];C=[[ps[i][k] for k in range(5)] for i in rows]
 invC=S.Matrix(C).inv();cnorm=max(sum(abs(F(z)) for z in invC.row(i)) for i in range(5))
 assert cnorm*5*n*r<1
 # Every coordinate nonzero at (s,t)=(1,0), robustly.
 assert all(abs(p[0])>n*r for p in ps)
 # Show at least one coordinate changes sign at retained rational projective points.
 change=False
 for p in ps[:5]:
  samples=[(p[0],n*r),(p[-1],n*r),(sum(p),5*n*r),(sum((-1)**k*z for k,z in enumerate(p)),5*n*r)]
  if any(v>err for v,err in samples) and any(v< -err for v,err in samples):change=True
 assert change
 return {'contraction_verified_exactly':True,'unique_real_root_in_fixed_parameter_box':True,'rhs_definite_verified':True,'degree_four_embedding_verified':True,'all_coordinates_nonzero_at_1_0':True,'all_definite':False,'rational_coefficients_verified':False,'eta':str(eta),'K':str(K),'eta_decimal':float(eta),'K_decimal':float(K),'hessian_entry_bound':str(H),'inverse_jacobian_norm':str(bnorm),'center_jacobian_inverse_error':str(error)}

def main():
 data=json.loads((OUT/'refinements.json').read_text());cs=[]
 for rr in data['records']:
  if rr['lane']=='elliptic':continue
  sys=RationalQuartic(rr['lane']=='double');fixed=rr['fixed_parameter_indices'];active=[i for i in range(sys.n) if i not in fixed]
  # 55 significant decimal digits, with fixed parameter retained EXACTLY as the rational used in refinement.
  x=[F(mp.nstr(mp.mpf(z[0]),55)) for z in rr['refined_parameters']]
  for i,z in zip(fixed,rr['fixed_parameter_values']):x[i]=F(z[0]).limit_denominator(1000)
  ps,val,JJ=ps_and_j(sys,x);J=[[row[j] for j in active] for row in JJ];B=approx_inv(J)
  rows=next(list(rows) for rows in itertools.combinations(range(6),5) if S.Matrix([ps[i] for i in rows]).det()!=0)
  c={'lane':rr['lane'],'source_start':rr['start'],'center':[str(z) for z in x],'active':active,'fixed':fixed,'inverse':enc(B),'radius':str(F(1,10**25)),'rank_minor_rows':rows,'coefficient_field':'isolated real algebraic root of the exact rational coefficient equations; NOT established to be Q'}
  c['verification']=verify_certificate(c);cs.append(c)
  print(rr['lane'],rr['start'],{k:v for k,v in c['verification'].items() if k.endswith('_decimal') or isinstance(v,bool)},flush=True)
 (OUT/'real_quartic_certificates.json').write_text(json.dumps({'proof':'T(x)=x-BF(x) on the closed infinity-norm radius-r cube. eta<=r/2 and K<1/2 give a self-map and strict contraction; Banach gives a unique real zero. H bounds every second partial throughout the cube; K=||I-BJ(center)||inf+||B||inf*n^2*H*r. All quantities used in the inequalities are rational.','records':cs},indent=2)+'\n')
if __name__=='__main__':main()
