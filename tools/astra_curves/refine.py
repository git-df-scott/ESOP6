#!/usr/bin/env python3
"""High precision refinement and exact rational rejection. No float is certified."""
import json,time
from fractions import Fraction
from pathlib import Path
import numpy as np, mpmath as mp
from scipy.linalg import qr
from search import EllipticQuartic,RationalQuartic
OUT=Path('results/astra_curves_2026_09_05')
mp.mp.dps=80

def mpvec(vals):return np.array([mp.mpc(str(z[0]),str(z[1])) if z[1] else mp.mpf(str(z[0])) for z in vals],dtype=object)
def norm(v):return max(abs(z) for z in v)
def fmt(z):return [mp.nstr(mp.re(z),78),mp.nstr(mp.im(z),78)]
def rational_parts(z,b):
    return [Fraction(str(mp.re(z))).limit_denominator(b),Fraction(str(mp.im(z))).limit_denominator(b)]
def exact_quartic(sys,x,cap):
    # Q-coefficients only. Nonreal coefficients have no Q reconstruction here.
    if max(abs(mp.im(z)) for z in x)>mp.mpf('1e-60'):return {'cap':cap,'eligible_real':False}
    xx=[Fraction(str(mp.re(z))).limit_denominator(cap) for z in x]
    base=sys.base.astype(int);AA=sys.A.astype(int)
    ps=[]
    for i in range(6):
        ps.append([Fraction(int(base[i,k]))+sum((int(AA[i,k,j])*xx[j] for j in range(sys.n)),Fraction()) for k in range(5)])
    def conv(p,q):
        r=[Fraction()]*(len(p)+len(q)-1)
        for i,a in enumerate(p):
            for j,b in enumerate(q):r[i+j]+=a*b
        return r
    res=[Fraction()]*25
    for i,p in enumerate(ps):
        h=[Fraction(1)]
        for _ in range(6):h=conv(h,p)
        res=[a+(1 if i<5 else -1)*b for a,b in zip(res,h)]
    nz=next((i for i,z in enumerate(res) if z),None)
    return {'cap':cap,'eligible_real':True,'exact_identity':nz is None,'first_nonzero_coefficient':nz,'residual_numerator':str(res[nz].numerator) if nz is not None else None,'residual_denominator':str(res[nz].denominator) if nz is not None else None,'parameters':[str(z) for z in xx]}

def refine(lane,r):
    sys=EllipticQuartic() if lane=='elliptic' else RationalQuartic(lane=='double');x=mpvec(r['parameters']);result={'lane':lane,'start':r['start'],'precision_digits':80}
    nn=np.array([complex(z) for z in x]);_,j=sys.both(nn)
    _,_,piv=qr(j,pivoting=True);cols=piv[:sys.m] if sys.n>sys.m else np.arange(sys.n)
    if sys.n>sys.m:
        fixed=[int(i) for i in range(sys.n) if i not in cols]
        for k in fixed:
            pair=rational_parts(x[k],1000);x[k]=mp.mpc(mp.mpf(pair[0].numerator)/pair[0].denominator,mp.mpf(pair[1].numerator)/pair[1].denominator)
        result['fixed_parameter_indices']=fixed;result['fixed_parameter_values']=[fmt(x[k]) for k in fixed]
    iterations=[]
    for it in range(14):
        v,j=sys.both(x);nv=norm(v);iterations.append(mp.nstr(nv,8))
        if nv<mp.mpf('1e-65'):break
        try:
            delta=mp.lu_solve(mp.matrix(j[:,cols].tolist()),mp.matrix((-v).tolist()))
        except (ZeroDivisionError,ValueError) as err:
            result['newton_error']=str(err);break
        damping=mp.mpf(1)
        for trial in range(12):
            y=x.copy()
            for k,d in zip(cols,delta):y[k]+=damping*d
            if norm(sys.both(y)[0])<nv:break
            damping/=2
        x=y
    result['residual_history']=iterations;result['final_residual_inf']=mp.nstr(norm(sys.both(x)[0]),12);result['refined_parameters']=[fmt(z) for z in x]
    result['exact_curve_certificate']=False
    if lane!='elliptic':result['rational_reconstructions']=[exact_quartic(sys,x,b) for b in (100,10000,1000000,1000000000)]
    else:
        # Exact rational checking uses the full polynomial expression later in verify.py.
        result['rational_reconstructions']=[{'cap':b,'parameters':[[str(y) for y in rational_parts(z,b)] for z in x]} for b in (100,10000,1000000)]
    return result

if __name__=='__main__':
    records=[];begin=time.time()
    for lane in ['single','double','elliptic']:
        data=json.loads((OUT/(lane+'.json')).read_text())
        for r in data['records']:
            if not r['small_residual_endpoint']:continue
            # Refine every real quartic endpoint and every elliptic endpoint except unfinished 31.
            if lane!='elliptic' and r['complex_start']:continue
            if lane=='elliptic' and r['start']==31:continue
            result=refine(lane,r);records.append(result)
            print(lane,r['start'],result['residual_history'],flush=True)
            (OUT/'refinements.json').write_text(json.dumps({'elapsed_s':time.time()-begin,'records':records},indent=2)+'\n')
