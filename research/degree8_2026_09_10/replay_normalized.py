"""Replay the s=2 conditioning check in u=2t coordinates."""
import json,time
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from search import setup,exact_check
rng=np.random.default_rng(6102027);rows=[];start=time.monotonic()
for eps in [-1,1]:
 unpack,fun,jac=setup(2.,eps,normalized=True)
 for seed in range(4):
  v=rng.normal(0,.08,20) if seed else np.zeros(20);v[0]=eps;v[5]=-eps
  direction=rng.normal(size=20);h=1e-6
  assert np.linalg.norm((fun(v+h*direction)-fun(v-h*direction))/(2*h)-jac(v)@direction)<1e-3
  z=least_squares(fun,v,jac=jac,max_nfev=600,ftol=1e-12,xtol=1e-12,gtol=1e-12)
  a,b,n=unpack(z.x);res=fun(z.x)
  original=[poly*2.**np.arange(len(poly)) for poly in (a,b,n)]
  exact,detail=exact_check(*original,2.) if max(abs(res))<1e-8 else (False,'above rational reconstruction threshold')
  row=dict(s=2,epsilon=eps,seed=seed,nfev=z.nfev,residual_max=float(max(abs(res))),exact_identity=exact,P=original[0].tolist(),R=original[1].tolist(),N=original[2].tolist());rows.append(row)
  print(json.dumps({k:row[k] for k in ['epsilon','seed','nfev','residual_max','exact_identity']}),flush=True)
  if exact:raise SystemExit('EXACT IDENTITY REQUIRES SPECIALIZATION')
Path(__file__).with_name('normalized_results.json').write_text(json.dumps(dict(coordinates='u=2t; original P,R,N coefficients retained',elapsed=time.monotonic()-start,rows=rows,counterexample_found=False),indent=2))
