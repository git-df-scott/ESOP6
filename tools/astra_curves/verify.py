#!/usr/bin/env python3
"""Independent symbolic checks and exact candidate rejection ledger."""
import json,itertools,math,sys
from fractions import Fraction as F
from pathlib import Path
import sympy as S
OUT=Path('results/astra_curves_2026_09_05')
u,v,w,t=S.symbols('u v w t')

def monomials(n,d):
 return [e for e in itertools.product(range(d+1),repeat=n) if sum(e)==d][::-1]
def term(z,e):return math.prod(a**b for a,b in zip(z,e))
def rational_elliptic_check(parameters):
 x=[F(z[0]) for z in parameters]
 if any(F(z[1]) for z in parameters):return {'over_Q_reconstruction':False,'exact_identity':False,'reason':'nonreal reconstructed coefficients'}
 m2=monomials(4,2);m4=monomials(4,4);m6=monomials(4,6)
 free=[i for i,e in enumerate(m2) if e not in [(2,0,0,0),(0,2,0,0)]];bfree=[i for i,e in enumerate(m4) if e[0]<2]
 q=[[F(0)]*10 for _ in range(2)];q[0][m2.index((2,0,0,0))]=1;q[1][m2.index((0,2,0,0))]=1
 for j,i in enumerate(free):q[0][i]=x[8+j];q[1][i]=x[16+j]
 a=x[24:59];b=[F(0)]*35
 for j,i in enumerate(bfree):b[i]=x[59+j]
 residual={e:F(0) for e in m6}
 for e in m6:
  coeff=math.factorial(6)//math.prod(math.factorial(k) for k in e)
  residual[e]=coeff*(term(x[:4],e)-term(x[4:8],e))+(1 if 6 in e else 0)
 for qq,h in [(q[0],a),(q[1],b)]:
  for e,c in zip(m2,qq):
   for ee,cc in zip(m4,h):residual[tuple(i+j for i,j in zip(e,ee))]-=c*cc
 nonzero=[(e,c) for e,c in residual.items() if c]
 return {'over_Q_reconstruction':True,'exact_identity':not nonzero,'nonzero_coefficient_count':len(nonzero),'first_nonzero_monomial':nonzero[0][0] if nonzero else None,'first_nonzero_coefficient':str(nonzero[0][1]) if nonzero else None}

def check_F7():
 data=json.loads((OUT/'quartic_F7_certificate.json').read_text());out=[]
 for r in data['records']:
  pp=sum(int(c)*t**i for i,c in enumerate(r['p']));qq=sum(int(c)*t**i for i,c in enumerate(r['q']))
  poly=S.Poly(pp**6+qq**6,t,modulus=7);const,factors=S.sqf_list(poly)
  exponents=[int(k) for p,k in factors];assert any(k%6 for k in exponents)
  out.append({'support':r['support'],'squarefree_multiplicities':exponents,'cannot_be_a_sixth_power':True})
 assert len(out)==35;return out

def symbolic_quartic(lane):
 if lane=='single':
  x=S.symbols('x0:14');p1=1+t+x[0]*t*t+x[1]*t**3+x[2]*t**4;p3=sum(x[3+i]*t**i for i in range(5));p5=x[8]+x[9]*t*t+x[10]*t**4;p6=x[11]+x[12]*t*t+x[13]*t**4
  ps=[p1,p1.subs(t,-t),p3,p3.subs(t,-t),p5,p6];rows=range(0,25,2)
 else:
  x=S.symbols('x0:8');p1=1+sum(x[i-1]*t**i for i in range(1,5));p3=t**4+sum(x[i-1]*t**(4-i) for i in range(1,5));p5=x[4]*(1+t**4)+x[5]*t*t;p6=x[6]*(1+t**4)+x[7]*t*t
  ps=[p1,p1.subs(t,-t),p3,p3.subs(t,-t),p5,p6];rows=range(0,13,2)
 poly=S.Poly(sum(z**6 for z in ps[:5])-ps[5]**6,t);eq=[poly.nth(i) for i in rows]
 assert all(poly.nth(i)==0 for i in range(1,25,2))
 if lane=='double':assert all(S.expand(poly.nth(i)-poly.nth(24-i))==0 for i in range(25))
 (OUT/(lane+'_ideal.json')).write_text(json.dumps({'variables':[str(z) for z in x],'coordinates':[str(p) for p in ps],'equations':[str(e) for e in eq],'unknowns':len(x),'equation_count':len(eq),'chart_conditions':'p1(s^4)=1; single additionally p1(s^3*t)=1, excluding the odd-leading-coefficient-zero subchart'},indent=2)+'\n')
 return x,ps,eq

def verify_real_symbols():
 from certify_real import verify_certificate
 data=json.loads((OUT/'real_quartic_certificates.json').read_text());models={l:symbolic_quartic(l) for l in ['single','double']};out=[]
 for c in data['records']:
  result=verify_certificate(c);x,ps,eq=models[c['lane']];center=list(map(S.Rational,c['center']));sub=dict(zip(x,center));active=c['active'];B=S.Matrix([[S.Rational(z) for z in row] for row in c['inverse']]);n=len(active)
  # SymPy builds/differentiates expressions independently of the coefficient-convolution harness.
  vals=S.Matrix([S.Poly(e,*x).eval(sub) for e in eq]);J=S.Matrix([[S.Poly(S.diff(e,x[j]),*x).eval(sub) for j in active] for e in eq]);E=S.eye(n)-B*J
  eta=max(abs(v) for v in B*vals);err=max(sum(abs(E[i,j]) for j in range(n)) for i in range(n))
  assert F(eta)==F(result['eta']);assert F(err)==F(result['center_jacobian_inverse_error'])
  out.append({'lane':c['lane'],'source_start':c['source_start'],'independent_symbolic_eta_and_inverse_error_agree':True,'exact_contraction_pass':True})
 return out

if __name__=='__main__':
 f=check_F7();r=verify_real_symbols();ell=[]
 for rec in json.loads((OUT/'refinements.json').read_text())['records']:
  if rec['lane']!='elliptic':continue
  for cand in rec['rational_reconstructions']:
   v=rational_elliptic_check(cand['parameters']);ell.append({'start':rec['start'],'cap':cand['cap'],**v});assert not v['exact_identity']
 data={'F7_independent_squarefree_replay':f,'real_root_symbolic_replay':r,'elliptic_exact_reconstruction_checks':ell,'rational_curve_found':False}
 (OUT/'independent_verification.json').write_text(json.dumps(data,indent=2)+'\n');print(json.dumps({'F7_pairs_verified':len(f),'real_certificates_verified':len(r),'elliptic_reconstructions_checked':len(ell),'exact_elliptic_identities':0},indent=2))
