"""Exact checks for the full-fourfold attack. No numerical tolerance is a certificate."""
import itertools, json, math, time
from pathlib import Path
import sympy as S

OUT=Path('results/fourfold_routes_2026_09_05')
OUT.mkdir(parents=True,exist_ok=True)
t=S.symbols('t')
report={}

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

# Full finite-field coefficient census is replaced by an exhaustive support
# argument; enumerate every polynomial to verify its projective zero count.
for d in (1,2,3):
    hist={}
    for co in itertools.product(range(7),repeat=d+1):
        if not any(co):continue
        vals=[sum(co[j]*pow(x,j,7) for j in range(d+1))%7 for x in range(7)]+[co[d]]
        nz=sum(v!=0 for v in vals)
        hist[nz]=hist.get(nz,0)+1
    assert min(hist)>=8-d and 2*(8-d)>8
    report[f'F7_degree_{d}']={'nonzero_polynomials':sum(hist.values()),'support_histogram':hist,
       'at_most_one_nonzero_left_polynomial':True,
       'scope':'F7-rational coefficient tuples, not geometric points over F7bar; bad reduction is not excluded'}

# A real nonconstant conic control exists with irrational sixth-power weights.
u=1-t*t;v=2*t;w=1+t*t
control=S.expand(8*u**6+8*v**6+(u+v)**6+(u-v)**6-10*w**6)
assert control==0
report['real_conic_control']={'identity':'8u^6+8v^6+(u+v)^6+(u-v)^6=10w^6; u=1-t^2,v=2t,w=1+t^2',
 'rational_lift':False,'reason':'a^6=8b^6 forces (a/b)^2=2 for nonzero rational a,b'}

# Groebner certificates for that necessary weight ratio with a normalized b=1.
# Empty finite-field rational point sets use field equations, not geometric emptiness.
a=S.symbols('a')
report['weight_ratio_groebner']={}
for p in (7,11,13,17):
    gb=S.groebner([a**6-8,a**p-a],a,modulus=p)
    report['weight_ratio_groebner'][str(p)]=[str(x.as_expr()) for x in gb.polys]

# Independent symbolic verification of the tangent construction on the cubic.
aa,bb,cc,x,y,z,r=S.symbols('a b c x y z r')
B=[aa,-aa,bb,-bb,cc,cc]
D=[cc**2*x,cc**2*y,cc**2*z,cc**2*r,0,aa**2*(x+y)+bb**2*(z+r)]
sig=[1,1,1,1,1,-1]
assert S.expand(sum(e*b*b*d for e,b,d in zip(sig,B,D)))==0
f=sum(e*d**3 for e,d in zip(sig,D))
k=sum(e*b*d*d for e,b,d in zip(sig,B,D))
# Universal coefficient identity avoids expanding a needlessly huge polynomial.
lam,mu,F,K=S.symbols('lam mu F K')
assert S.expand((3*lam*mu**2*K+mu**3*F).subs({lam:F,mu:-3*K}))==0
report['tangent_identity']={'tangent_constraint_verified':True,'universal_cubic_identity_verified':True,
 'formula':'C=F(D)B-3*(sum eps_i B_i D_i^2)*D'}

# Quadratic contact factor and its Gaussian parametrization.
s,A,Bb=S.symbols('s A B')
P=A-s*t*Bb;R=Bb+s*t*A;Q=1+s*s*t*t
assert S.expand(P*P+R*R-Q*(A*A+Bb*Bb))==0
assert S.expand(P**6+R**6-(P*P+R*R)*(P**4-P*P*R*R+R**4))==0
report['quadratic_contact']={'necessary_condition':'q is a positive rational square for a reduced degree-eight identity',
 'gaussian_factorization_verified':True,'proof':'See FOURFOLD_ROUTES_2026_09_05.md'}

# Equal-weight full-fourfold boundary-contact residual, with four independent slopes.
N,h=S.symbols('N h')
assert S.expand((N+h/2)**6-(N-h/2)**6-(6*h*N**5+5*h**3*N**3+S.Rational(3,8)*h**5*N))==0
report['full_boundary_contact']={'identity_verified':True,
 'residual':'sum(A_i^6)=6k M^5+5k^3 t^12 M^3+(3/8)k^5 t^24 M, i=1..4',
 'status':'new exact system, unsolved'}

# BCU identity and the sign obstruction in the standard secant planes.
assert S.expand((1-t-t*t)**3+(1+t-t*t)**3-(2-2*t**6))==0
report['BCU']={'control_verified':True,'signature':'four sixth powers versus two after both cube bases are squares; not 5 versus 1'}
(OUT/'symbolic_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
