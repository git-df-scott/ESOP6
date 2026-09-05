"""An exact, explicitly restricted degree-eight Gaussian-factor attack."""
import json,time
from pathlib import Path
import sympy as S
t,s,u,v,z=S.symbols('t s u v z')
out=Path('results/fourfold_routes_2026_09_05')
P=1+u*t**6-s*t*(1+v*t**6)
R=1+v*t**6+s*t*(1+u*t**6)
Q=1+s*s*t*t
D=S.Poly(S.expand((P**6+R**6)/2),t)
def mul(a,b,n):
 c=[S.S(0)]*(n+1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):
   if i+j<=n:c[i+j]+=x*y
 return [S.expand(x) for x in c]
def power(a,k,n):
 b=[S.S(1)]+[S.S(0)]*n
 for _ in range(k):b=mul(b,a,n)
 return b
N=[S.S(1)]+[S.S(0)]*8
for n in range(1,9):
 rhs=mul([1,0,s*s],power(N,5,n),n)
 N[n]=S.expand((D.nth(n)-rhs[n])/5)
 print('eliminated',n,flush=True)
rhs=mul([1,0,s*s],power(N,5,16),16)
# The t^12 correction enters the residual but not the first eight equations.
corr=mul(power([1,0,s*s],3,4),power(N,3,4),4)
res=[]
for n in range(9,17):
 r=S.factor(D.nth(n)-rhs[n]-(S.Rational(10,27)*corr[n-12] if n>=12 else 0))
 res.append(r)
 print('residual',n,str(r),flush=True)
record={'family':'P=1+u*t^6-s*t*(1+v*t^6), R=1+v*t^6+s*t*(1+u*t^6), Q=1+s^2*t^2',
 'N_coefficients':[str(x) for x in N], 'residual_9_to_16':[str(x) for x in res], 'groebner':{}}
eq=[S.fraction(x)[0] for x in res]+[s*z-1]
for p in (7,11,13):
 start=time.monotonic();gb=S.groebner(eq,z,u,v,s,modulus=p)
 record['groebner'][str(p)]={'basis':[str(x.as_expr()) for x in gb.polys],'seconds':time.monotonic()-start,
 'scope':'necessary residual equations in this sparse family; finite-characteristic result only'}
 print('mod',p,record['groebner'][str(p)],flush=True)
start=time.monotonic();gb=S.groebner(eq,z,u,v,s)
record['groebner']['QQ']={'basis':[str(x.as_expr()) for x in gb.polys],'seconds':time.monotonic()-start}
(out/'degree8_sparse.json').write_text(json.dumps(record,indent=2)+'\n')
print('QQ',record['groebner']['QQ'],flush=True)
