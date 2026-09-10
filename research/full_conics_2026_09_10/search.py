"""Full fourfold real conic pilot. Rational reconstruction is verified exactly."""
import json,time,sys
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from numpy.polynomial.polynomial import polypow,polymul
import sympy as S

out=Path(__file__).parent
def pad(a,shift=0):
 r=np.zeros(13);r[shift:shift+len(a)]=a;return r
def main():
 rng=np.random.default_rng(6112026); rows=[];start=time.monotonic()
 for q in [1,2,3]:
  W=np.array([1.,0.,float(q)])
  k=(63/80)**(1/6); root=np.sqrt(q)
  seed=np.array([[k,0,-k*q],[0,2*k*root,0],[k/np.sqrt(2),np.sqrt(2)*k*root,-k*q/np.sqrt(2)],[k/np.sqrt(2),-np.sqrt(2)*k*root,-k*q/np.sqrt(2)],[.5,0,.5*q]])
  assert np.max(np.abs(sum(pad(polypow(a,6)) for a in seed)-pad(polypow(W,6))))<1e-10
  # Fix a5=1/2,b5=0, leaving 13 variables and 13 coefficient equations.
  def unpack(v):return np.r_[v[:12],.5,0,v[12]].reshape(5,3)
  def fun(v):return sum(pad(polypow(a,6)) for a in unpack(v))-pad(polypow(W,6))
  def jac(v):
   a=unpack(v);cols=[]
   for i in range(4):
    for j in range(3):cols.append(6*pad(polypow(a[i],5),j))
   cols.append(6*pad(polypow(a[4],5),2));return np.array(cols).T
  for n in range(4):
   v=(rng.normal(0,.4,13) if "--random" in sys.argv else np.r_[seed[:4].ravel(),seed[4,2]]+rng.normal(0,.03,13))
   direction=rng.normal(size=13);h=1e-6
   assert np.linalg.norm((fun(v+h*direction)-fun(v-h*direction))/(2*h)-jac(v)@direction)<1e-3
   fit=least_squares(fun,v,jac=jac,max_nfev=1000,ftol=1e-13,xtol=1e-13,gtol=1e-13)
   coeff=unpack(fit.x);res=float(max(abs(fun(fit.x))));checks=[]
   if res<1e-8:
    t=S.Symbol('t')
    for bound in [10,100,1000,10000]:
     polys=[sum(S.Rational(str(x)).limit_denominator(bound)*t**j for j,x in enumerate(a)) for a in coeff]
     rem=S.Poly(sum(a**6 for a in polys)-(1+q*t*t)**6,t)
     checks.append(dict(denominator_bound=bound,exact=rem.is_zero))
     if rem.is_zero:
      (out/'exact_candidate.json').write_text(json.dumps(dict(q=q,polys=list(map(str,polys)))))
      raise SystemExit('EXACT IDENTITY: SPECIALIZE AND VERIFY')
   row=dict(q=q,start=n,nfev=fit.nfev,status=int(fit.status),residual_max=res,coefficients=coeff.tolist(),rational_checks=checks)
   rows.append(row);print(json.dumps({key:row[key] for key in ['q','start','nfev','residual_max']}),flush=True)
 (out/('random_results.json' if '--random' in sys.argv else 'results.json')).write_text(json.dumps(dict(domain='W=1+q*t^2 q=1,2,3; a5=1/2,b5=0;four starts each;seed6112026',initialization=('normal(0,.4)' if '--random' in sys.argv else 'known real conic plus normal(0,.03)'),rows=rows,elapsed=time.monotonic()-start,counterexample_found=False),indent=2))
if __name__=='__main__':main()
