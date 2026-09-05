"""First exact Hensel obstruction for the four transverse conic patterns."""
import json,itertools,hashlib
from pathlib import Path
import sympy as S
out=Path('results/fourfold_routes_2026_09_05')
patterns=json.loads((out/'conic_F7_transverse.json').read_text())['nontrivial_solutions']
t=S.symbols('t')
def mul(a,b):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return c
def power(a,k):
 c=[1]
 for _ in range(k):c=mul(c,a)
 return c
def add_at(c,a,fac=1):
 for i,v in enumerate(a):c[i]+=fac*v
def residual(x,delta=0):
 c=[0]*13
 for i in range(4):add_at(c,power(x[3*i:3*i+3],6))
 A=[x[12],x[13],1];B=x[14:17]
 add_at(c,mul(power(A,5),B),-6)
 if delta:
  for j in range(2,7):
   import math
   add_at(c,mul(power(A,6-j),power(B,j)),-math.comb(6,j)*delta**(j-1))
 return c
def jac(x):
 cols=[]
 for i in range(4):
  P=x[3*i:3*i+3];pow5=power(P,5)
  for j in range(3):cols.append([0]*j+[6*v for v in pow5]+[0]*(2-j))
 A=[x[12],x[13],1];B=x[14:17]
 for j in range(2):cols.append([0]*j+[-30*v for v in mul(power(A,4),B)]+[0]*(2-j))
 for j in range(3):cols.append([0]*j+[-6*v for v in power(A,5)]+[0]*(2-j))
 return [[col[i]%7 for col in cols] for i in range(13)]
def solve(A,b):
 M=[row[:]+[v%7] for row,v in zip(A,b)];r=0;piv=[]
 for j in range(17):
  k=next((k for k in range(r,13) if M[k][j]%7),None)
  if k is None:continue
  M[r],M[k]=M[k],M[r];inv=pow(M[r][j],-1,7);M[r]=[v*inv%7 for v in M[r]]
  for k in range(13):
   if k!=r:
    f=M[k][j];M[k]=[(a-f*b)%7 for a,b in zip(M[k],M[r])]
  piv.append(j);r+=1
  if r==13:break
 if any(not any(row[:17]) and row[17] for row in M):return None,r,17-r,[]
 x=[0]*17
 for k,j in enumerate(piv):x[j]=M[k][17]
 free=[j for j in range(17) if j not in piv];ker=[]
 for j in free:
  v=[0]*17;v[j]=1
  for k,i in enumerate(piv):v[i]=-M[k][j]%7
  ker.append(v)
 return x,r,len(free),ker
records=[]
for idx,Ps in enumerate(patterns):
 sum_poly=sum(sum(P[j]*t**j for j in range(3))**6 for P in Ps)
 quo,rem=S.div(S.Poly(sum_poly,t,modulus=7),S.Poly((1+t*t)**5,t,modulus=7))
 assert rem.is_zero
 B=[int(quo.nth(i))*pow(6,-1,7)%7 for i in range(3)]
 x=sum(Ps,[])+[1,0]+B
 J=jac(x);res=residual(x);assert all(v%7==0 for v in res)
 dig,rank,nullity,ker=solve(J,[-v//7 for v in res])
 rec={'pattern':idx+1,'initial_vector':x,'jacobian_rank':rank,'nullity':nullity,'lift_to_49':dig is not None}
 if dig is not None:
  y=[a+7*b for a,b in zip(x,dig)];assert all(v%49==0 for v in residual(y))
  rec['one_lift_mod49']=y
  # One chosen path is evidence of existence, never completeness of branching.
  steps=[]
  for n in range(2,7):
   mod=7**n;rr=residual(y);assert all(v%mod==0 for v in rr)
   dd,_,_,_=solve(J,[-v//mod for v in rr])
   steps.append({'modulus':7**(n+1),'chosen_path_lifts':dd is not None})
   if dd is None:break
   y=[a+mod*b for a,b in zip(y,dd)]
  rec['chosen_path']=steps
  rec['unperturbed_final_vector']=y
  # Full conic: first four coordinates 7*R_i, fifth A, sixth A+7^6 B.
  y=x[:];full_steps=[]
  for n in range(1,41):
   mod=7**n;rr=residual(y,7**6);assert all(v%mod==0 for v in rr)
   dd,_,_,_=solve(J,[-v//mod for v in rr])
   full_steps.append({'modulus':str(7**(n+1)),'chosen_path_lifts':dd is not None})
   if dd is None:break
   y=[a+mod*b for a,b in zip(y,dd)]
  rec['full_conic_chosen_path']=full_steps
  rec['full_conic_final_vector']=[str(v) for v in y]
  rec['full_conic_residual_zero']=all(v==0 for v in residual(y,7**6))
 records.append(rec)
print(json.dumps(records,indent=2))
(out/'conic_lifts.json').write_text(json.dumps(records,indent=2)+'\n')
