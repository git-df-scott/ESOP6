"""Bounded numerical discovery; only exact rational polynomial equality is success."""
import json, time
from fractions import Fraction
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from numpy.polynomial.polynomial import polymul as mul, polypow as power

OUT=Path(__file__).parent
def padded(a, shift=0):
    b=np.zeros(43); b[shift:shift+len(a)]=a; return b

def setup(s,eps,normalized=False):
    physical_s=s
    if normalized: s=1.
    k=physical_s**-12 if normalized else 1.
    q=np.array([1.,0.,s*s]); q3=power(q,3); q5=power(q,5)
    # 12 independent P/R coefficients; enforce both remainders at t=i/s.
    def unpack(v):
        p=np.r_[1.,v[:5],0.,0.]; r=np.r_[1.,v[5:12]]
        p[6]=s**6-s**4*p[2]+s*s*p[4]+eps*(s**5*r[1]-s**3*r[3]+s*r[5]-r[7]/s)
        p[7]=s**6*p[1]-s**4*p[3]+s*s*p[5]-eps*(s**7-s**5*r[2]+s**3*r[4]-s*r[6])
        n=np.r_[1.,v[12:20]]
        return p,r,n
    base=unpack(np.zeros(20))
    dirs=[tuple(x-y for x,y in zip(unpack(np.eye(20)[j]),base)) for j in range(20)]
    def fun(v):
        p,r,n=unpack(v)
        return (padded((power(p,6)+power(r,6))/2)-padded(mul(q,power(n,5)))
          -10/27*k*padded(mul(q3,power(n,3)),12)-1/81*k*k*padded(mul(q5,n),24))[1:]
    def jac(v):
        p,r,n=unpack(v); p5=power(p,5);r5=power(r,5);n4=power(n,4);n2=power(n,2)
        return np.array([(3*padded(mul(p5,dp))+3*padded(mul(r5,dr))
            -5*padded(mul(mul(q,n4),dn))-10/9*k*padded(mul(mul(q3,n2),dn),12)
            -1/81*k*k*padded(mul(q5,dn),24))[1:] for dp,dr,dn in dirs]).T
    return unpack,fun,jac

def exact_check(p,r,n,s):
    import sympy as S
    t=S.Symbol('t'); q=1+S.Rational(str(s))**2*t*t
    convert=lambda a:sum(S.Rational(str(x)).limit_denominator(1000)*t**i for i,x in enumerate(a))
    P,R,N=map(convert,(p,r,n))
    residual=S.Poly((P**6+R**6)/2-q*N**5-S.Rational(10,27)*t**12*q**3*N**3-S.Rational(1,81)*t**24*q**5*N,t)
    return residual.is_zero, str(S.factor(residual.as_expr())) if residual.is_zero else str(residual.degree())

def main():
    rng=np.random.default_rng(6102026); rows=[]; start=time.monotonic()
    for s in [0.5,1.,2.]:
      for eps in [-1,1]:
       unpack,fun,jac=setup(s,eps)
       for seed in range(4):
        if time.monotonic()-start>180: break
        v=rng.normal(0,0.08,20) if seed else np.zeros(20)
        v[0]=eps*s;v[5]=-eps*s
        check_v=v+rng.normal(0,.001,20)
        direction=rng.normal(0,1,20);h=1e-6
        err=np.linalg.norm((fun(check_v+h*direction)-fun(check_v-h*direction))/(2*h)-jac(check_v)@direction)
        assert err<1e-4*(1+np.linalg.norm(jac(check_v)@direction)),err
        sol=least_squares(fun,v,jac=jac,max_nfev=600,ftol=1e-12,xtol=1e-12,gtol=1e-12)
        p,r,n=unpack(sol.x); residual=fun(sol.x)
        exact,detail=exact_check(p,r,n,s) if np.max(np.abs(residual))<1e-8 else (False,'not attempted: residual above threshold')
        row=dict(s=s,epsilon=eps,seed=seed,nfev=sol.nfev,status=int(sol.status),residual_max=float(np.max(np.abs(residual))),residual_norm=float(np.linalg.norm(residual)),exact_identity=exact,exact_check=detail,P=p.tolist(),R=r.tolist(),N=n.tolist())
        rows.append(row); (OUT/'results.json').write_text(json.dumps(dict(domain='s=1/2,1,2; epsilon=+-1; four starts each; seed6102026; max600nfev/start;180s between-start cap',rows=rows,elapsed=time.monotonic()-start,counterexample_found=False),indent=2))
        print(json.dumps({k:row[k] for k in ['s','epsilon','seed','nfev','residual_max','exact_identity']}),flush=True)
        if exact: raise SystemExit('EXACT IDENTITY: SPECIALIZE AND VERIFY BEFORE CLAIM')
    print('DONE',len(rows),'No exact identity found.',flush=True)
if __name__=='__main__':main()
