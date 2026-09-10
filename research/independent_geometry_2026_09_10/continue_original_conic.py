#!/usr/bin/env python3
"""Canonical branches only; this is not exhaustive Hensel lifting."""
import json, math
from fractions import Fraction
from pathlib import Path
root=Path(__file__).parent
seeds=json.loads((root/'conic_tangent_seeds.json').read_text())
def mul(a,b):
 out=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):out[i+j]+=x*y
 return out
def power(a,n):
 out=[1]
 for _ in range(n):out=mul(out,a)
 return out
def residual(forms):
 U=forms[:4];W,H=forms[4:]
 out=[sum(power(u,6)[k] for u in U) for k in range(13)]
 for j in range(1,7):
  c=math.comb(6,j)*(-1)**j*7**(6*(j-1))
  term=mul(power(W,6-j),power(H,j))
  out=[a+c*b for a,b in zip(out,term)]
 return out
def jacobian(forms):
 U=forms[:4];W,H=forms[4:]
 bases=[[6*v%7 for v in power(u,5)] for u in U]
 bases += [[-30*v%7 for v in mul(power(W,4),H)],[-6*v%7 for v in power(W,5)]]
 cols=[]
 for b in bases:
  for j in range(3):cols.append([0]*j+b+[0]*(2-j))
 return list(zip(*cols))
def solve(jac,rhs):
 a=[[int(v)%7 for v in row]+[int(b)%7] for row,b in zip(jac,rhs)];pivs=[];r=0
 for c in range(18):
  p=next((i for i in range(r,13) if a[i][c]),None)
  if p is None:continue
  a[r],a[p]=a[p],a[r]; inv=pow(a[r][c],-1,7);a[r]=[v*inv%7 for v in a[r]]
  for i in range(13):
   if i!=r:
    z=a[i][c];a[i]=[(v-z*w)%7 for v,w in zip(a[i],a[r])]
  pivs.append(c);r+=1
 if any(not any(row[:18]) and row[18] for row in a):return None
 out=[0]*18
 for i,c in enumerate(pivs):out[c]=a[i][18]
 return out
def reconstruct(a,m,B=10000):
 a%=m
 if not a:return Fraction(0)
 A=(m-1)//(2*B);r0,r1=m,a;s0,s1=0,1
 while abs(r1)>A:
  q=r0//r1;r0,r1=r1,r0-q*r1;s0,s1=s1,s0-q*s1
 if not s1 or abs(s1)>B or math.gcd(r1,s1)!=1:return None
 f=Fraction(r1,s1)
 return f if (f.numerator-a*f.denominator)%m==0 else None
branches=[{'seed_index':i,'forms':s['canonical_mod343_lift']} for i,s in enumerate(seeds) if s['canonical_mod343_lift'] is not None]
counts={'3':len(branches)};failed=[];m=343
for exponent in range(4,13):
 keep=[]
 for b in branches:
  forms=b['forms'];F=residual(forms);assert all(v%m==0 for v in F)
  d=solve(jacobian(forms),[-v//m for v in F])
  if d is None:
   failed.append({'seed_index':b['seed_index'],'first_failed_exponent':exponent});continue
  new=[[x+m*d[3*i+j] for j,x in enumerate(form)] for i,form in enumerate(forms)]
  assert all(v%(m*7)==0 for v in residual(new))
  keep.append({'seed_index':b['seed_index'],'forms':new})
 branches=keep;m*=7;counts[str(exponent)]=len(branches)
reconstructions=[]
for b in branches:
 rr=[[reconstruct(x,m) for x in f] for f in b['forms']]
 complete=all(x is not None for f in rr for x in f)
 exact=complete and all(x==0 for x in residual(rr))
 reconstructions.append({'seed_index':b['seed_index'],'all_coefficients_reconstructed':complete,'exact_original_identity':exact,'forms':[[str(x) if x is not None else None for x in f] for f in rr]})
result={'modulus':m,'surviving_canonical_branches_by_exponent':counts,'failed_canonical_branches':failed,'final_branches':branches,'reconstructions':reconstructions,'exact_identity_count':sum(r['exact_original_identity'] for r in reconstructions),'scope':'Only canonical free-zero corrections; branch failures do not exclude seed families. Original divided equation used, including terms from 7^6 onward.'}
(root/'original_conic_lifts.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ('failed_canonical_branches','final_branches','reconstructions')},indent=2))
