from fractions import Fraction as Fr
from math import isqrt, gcd
import sympy as sp
c = sp.symbols('c')
P = (c+1)*(4*c**2-2*c-1)
U = 6*P + 25*c*(c**2-1); V = 8*P
D = sp.expand(U**2 + V**2)
print("D(c) =", sp.factor(D))
print("degree", sp.degree(D, c), " squarefree:", sp.gcd(D, sp.diff(D,c)) == 1)
def is_sq(fr):
    if fr < 0: return False
    n, d = fr.numerator, fr.denominator
    return isqrt(n)**2 == n and isqrt(d)**2 == d
Dp = sp.Poly(D, c); coeffs = [Fr(int(x)) for x in Dp.all_coeffs()]
def Dval(x):
    r = Fr(0)
    for co in coeffs: r = r*x + co
    return r
hits = []
H = 300
for b in range(1, H+1):
    for a in range(-b, b+1):
        if gcd(abs(a), b) != 1: continue
        x = Fr(a, b)
        if is_sq(Dval(x)): hits.append(x)
print(f"rational c=a/b, |a|<=b<={H}, with D(c) a rational square: {len(hits)}")
for x in hits:
    Pv = (x+1)*(4*x*x-2*x-1); Uv = 6*Pv + 25*x*(x*x-1)
    y = Fr(isqrt(Dval(x).numerator), isqrt(Dval(x).denominator))
    ms = [(Uv + s*y)/(2*Pv*Pv) for s in (1,-1)] if Pv != 0 else []
    print(f"  c={x}  D={Dval(x)}  m candidates={[str(m) for m in ms]}  1-c^2 square? {is_sq(1-x*x)}")
