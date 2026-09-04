#!/usr/bin/env python3
"""Small exact identities and formal jets for the second strike (SymPy)."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import sympy as s
from boundary_contact import formal_midpoint, mul, power


def checks():
    t, m, c = s.symbols('t M c')
    u2,u3,u4,u5,v1,v2,v3,v4,v5 = s.symbols('u2 u3 u4 u5 v1 v2 v3 v4 v5')
    midpoint = s.cancel(((m+t**6/3)**6-(m-t**6/3)**6)/(4*t**6))
    assert s.expand(midpoint-m**5-s.Rational(10,27)*t**12*m**3-s.Rational(1,81)*t**24*m) == 0
    assert s.factor(midpoint) == m*(3*m*m+t**12)*(27*m*m+t**12)/81
    assert s.expand(4*t**6-((1+c*t**6)**6-1)).coeff(t,6) == 4-6*c

    u = [s.Integer(1),s.Integer(0),u2,u3,u4,u5]
    v = [s.Integer(0),v1,v2,v3,v4,v5]
    def product(a,b,n=9):
        return [s.expand(x) for x in mul(a,b,limit=n)]
    def pp(a,k,n=9):
        out=[s.Integer(1)]
        for _ in range(k):out=product(out,a,n)
        return out+[s.Integer(0)]*(n-len(out))
    d = [s.expand(x) for x in map(sum,zip(pp(u,6),
         [15*x for x in product(pp(u,4),pp(v,2))],
         [15*x for x in product(pp(u,2),pp(v,4))],pp(v,6)))]
    mm=[s.Integer(1)]
    for k in range(1,7):mm.append(s.expand((d[k]-pp(mm,5,k+1)[k])/5))
    m5=pp(mm,5)
    residuals=[s.factor(d[k]-m5[k]) for k in [7,8]]
    matrix=s.Matrix([[s.diff(r,z) for z in [u5,v5]] for r in residuals])
    assert all(s.diff(r,z,2)==0 for r in residuals for z in [u5,v5])
    constant=s.symbols('q')
    sparse=[]
    for k in range(1,6):
        dd=1+15*constant**2*t**(2*k)+15*constant**4*t**(4*k)+constant**6*t**(6*k)
        mt=s.series(dd**s.Rational(1,5),t,0,7).removeO()
        e=s.Poly(s.expand(dd-mt**5-s.Rational(10,27)*t**12*mt**3-s.Rational(1,81)*t**24*mt),t)
        first=min(i[0] for i,a in e.terms() if a!=0)
        sparse.append({'U':'1','V':f'q*t^{k}','midpoint_truncation':str(mt),
                       'first_residual_order':first,
                       'first_residual_coefficient':str(s.factor(e.coeff_monomial(t**first)))})

    formal=[]
    for a,b,degree,label in [([1],[1],72,'A=B=1'),([1,1],[1,-1],30,'A=1+t, B=1-t')]:
        mmq=formal_midpoint(a,b,degree)
        dd=[(x+y)*F(1,2) for x,y in zip(power(list(map(F,a)),6,limit=degree+1),
                                  power(list(map(F,b)),6,limit=degree+1))]
        fifth=power(mmq,5,limit=degree+1);third=power(mmq,3,limit=degree+1)
        for k in range(degree+1):
            rhs=fifth[k]+(F(10,27)*third[k-12] if k>=12 else 0)+(F(1,81)*mmq[k-24] if k>=24 else 0)
            assert rhs==dd[k],(label,k)
        formal.append({'input':label,'verified_through_order':degree,
                       'nonzero_midpoint_coefficients':{str(i):str(x) for i,x in enumerate(mmq) if x}})
    # One rational reconstruction is recorded, then rejected by its exact
    # numerator. Formal matching is not an identity certificate.
    z=s.symbols('z');pade=(1-z/30)/(1+11*z/270)
    error=s.factor(s.together(pade**5+s.Rational(10,27)*z*pade**3+s.Rational(1,81)*z*z*pade-1))
    assert s.series(error,z,0,4).removeO()==-s.Rational(121,196830)*z**3

    # Verify the algebra used to exclude the exceptional Gauss valuation -1.
    a,b=s.symbols('a b');ss=a*a+b*b
    defect=s.expand(a**6+b**6+6*ss**5+20*ss**3+6*ss)
    reduced=s.Poly(defect/3-ss*((a*a-b*b)**2-ss**4-1),a,b,modulus=3)
    assert reduced.is_zero
    p,h,j,xy=s.symbols('P h j xy')
    small=s.Poly(s.expand((2*p**5+3*j)**3-3*xy*(2*p**5+3*j)-2*(p**3+3*h)**5),p,h,j,xy)
    expected=3*p**5*(2*p**10-2*xy-p**7*h)
    assert all(int(cc)%9==0 for cc in s.Poly(small.as_expr()-expected,p,h,j,xy).coeffs())

    return {'result':'PASS','leading_contact_equation':'4-6*c=0',
            'midpoint_factorization':str(s.factor(midpoint)),
            'midpoint_coefficients_1_to_6':{str(k):str(s.factor(mm[k])) for k in range(1,7)},
            'first_residuals':{str(k):str(r) for k,r in zip([7,8],residuals)},
            'linear_matrix_for_u5_v5':[[str(s.factor(x)) for x in row] for row in matrix.tolist()],
            'determinant':str(s.factor(matrix.det())),
            'sparse_asymmetric_attempts':sparse,'formal_jets':formal,
            'pade_variable':'z=t^12','pade_1_1':str(pade),
            'pade_exact_residual':str(error),
            'pade_first_residual':'-121*t^36/196830',
            'valuation_minus_one_reduction':'PASS', 'modulo_nine_factor_reduction':'PASS'}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path)
    args=parser.parse_args();out=checks();encoded=json.dumps(out,indent=2)+'\n'
    if args.output:args.output.write_text(encoded)
    print(encoded,end='')


if __name__=='__main__':main()
