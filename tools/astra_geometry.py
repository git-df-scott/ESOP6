#!/usr/bin/env python3
"""Exact symbolic certificates for the September 4 direct surface strike.

Requires SymPy. These are identities, not a certificate of nonexistence of
rational points on the surface. See GEOMETRIC_STRIKE.md for the proofs and
the precise limits of the low-degree argument.
"""
import json
import sympy as s


def main():
    t, q, A, B, D, lam = s.symbols('t q A B D lam')
    X, Y, Z, W = s.symbols('X Y Z W')
    F = 2*X**6 + 2*Y**6 + Z**6 - W**6
    bcu = (1-t-t*t)**3 + (1+t-t*t)**3 - 2 + 2*t**6
    assert s.expand(bcu) == 0

    # Cubic quotient coordinates U=X^2,V=Y^2,R=Z^2,T=W^2.
    # A=U+V, B=U-V, C=T-R=lam*A, D=T+R.
    U, V = (A+B)/2, (A-B)/2
    R, T = (D-lam*A)/2, (D+lam*A)/2
    conic = (2-lam**3)*A*A + 6*B*B - 3*lam*D*D
    assert s.expand(4*(2*U**3+2*V**3+R**3-T**3)-A*conic) == 0

    # Exact rational parametrization of the conic at lam=1/2.
    Q = [t*t+8*t-5, -t*t+8*t+5,
         2*(t*t-2*t+5), 2*(t*t+2*t+5)]
    assert s.expand(2*Q[0]**3+2*Q[1]**3+Q[2]**3-Q[3]**3) == 0
    assert s.expand(2*(Q[3]-Q[2])-(Q[0]+Q[1])) == 0

    # Rational exceptional member lam=2: no nonzero real square lift.
    lifted = s.expand(conic.subs({A:X*X+Y*Y, B:X*X-Y*Y,
                                 D:W*W+Z*Z, lam:2}))
    assert s.expand(lifted+6*((W*W+Z*Z)**2+4*X*X*Y*Y)) == 0

    # Boundary expansion in the chart Z=1, W=1+q.
    contact = s.expand((1+q)**6-1)
    assert s.expand(contact-q*(6+15*q+20*q*q+15*q**3+6*q**4+q**5)) == 0

    # S=x^6+y^6, h=x^2+y^2: divisor/discriminant reconstruction.
    h = X*X+Y*Y
    total = X**6+Y**6
    assert s.expand(total-h*(h*h-3*X*X*Y*Y)) == 0
    assert s.cancel((4*total/h-h*h)/3-(X*X-Y*Y)**2) == 0

    # Boundary, diagonal, and conic-discriminant reductions.
    assert s.expand((W**3-Z**3)*(W**3+Z**3)-(W**6-Z**6)) == 0
    # Minimal-contact boundary normalization with A5(0)=B5(0)=1.
    contact_c=s.Rational(2,3)
    parity_target=s.cancel(((1+contact_c*t**6/2)**6-
                            (1-contact_c*t**6/2)**6)/(2*t**6))
    assert s.expand(parity_target-(2+s.Rational(20,27)*t**12+
                                  s.Rational(2,81)*t**24)) == 0
    alpha,beta=s.symbols('alpha beta')
    # Matrix pencil for the two quadrics on the complex lambda=2 fibre.
    matrix=s.Matrix([[-2*alpha,s.I*beta,0,0],
                     [s.I*beta,-2*alpha,0,0],
                     [0,0,-alpha+beta,0],[0,0,0,alpha+beta]])
    assert s.factor(matrix.det()-(4*alpha**2+beta**2)*(beta**2-alpha**2)) == 0
    out = {
        'bcu_identity': 'PASS', 'cubic_quotient_identity': 'PASS',
        'lambda_half_parametrization': 'PASS',
        'lambda_two_real_obstruction': 'PASS', 'boundary_contact': 'PASS',
        'divisor_discriminant_identity': 'PASS',
        'conic_determinant': str(s.factor(s.prod([2-lam**3,6,-3*lam]))),
        'coordinate_tangency_factors': ['8-lambda^3','2-4*lambda^3'],
        'generic_branch_points': 8, 'generic_cover_degree': 8,
        'generic_lift_genus': 1 + (8*(-2)+8*4)//2,
        'degree_two_wronskian_required_degree': 32,
        'degree_two_wronskian_maximum_degree': 30,
        'degree_four_classification': 'OPEN; boundary-contact locus excluded',
        'lambda_valuations': {'2':'4*r2+1, r2>=1','3':'4*r3-1, r3>=1',
                             '7':'4*r7 or -2*r7, r7>=1'},
        'symmetric_degree_six_target':str(parity_target),
        'symmetric_degree_six_obstruction':'leading coefficient a4^6=1/81 has no rational solution',
        'complex_lambda_two_quadric_determinant':str(s.factor(matrix.det())),
    }
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
