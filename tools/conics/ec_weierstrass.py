# Convert y^2 = 137c^4 - 186c^3 + 21c^2 + 28c + 4 (point c=0,y=2) to Weierstrass form and search points.
import sympy as sp
c, y = sp.symbols('c y')
Q = 137*c**4 - 186*c**3 + 21*c**2 + 28*c + 4
# quartic y^2 = a c^4 + b c^3 + cc c^2 + d c + e^2 with e=2 : standard reduction (Mordell / Cassels):
# For y^2 = a x^4 + b x^3 + c x^2 + d x + q^2, set ... use invariants I, J of the binary quartic -> E: Y^2 = X^3 - 27 I X - 27 J
a,b,cc,d,e = 137,-186,21,28,4
I = 12*a*e - 3*b*d + cc**2
J = 72*a*cc*e + 9*b*cc*d - 27*a*d**2 - 27*e*b**2 - 2*cc**3
print("I =", I, " J =", J)
X, Y = sp.symbols('X Y')
E = Y**2 - (X**3 - 27*I*X - 27*J)
print("Jacobian curve: Y^2 = X^3 - 27 I X - 27 J  =>", sp.expand(X**3 - 27*I*X - 27*J))
disc = -16*(4*(-27*I)**3 + 27*(-27*J)**2)
print("discriminant:", disc, sp.factorint(abs(disc)))
# minimal-ish: also print the quartic's discriminant
print("quartic disc factors:", sp.factorint(sp.discriminant(Q, c)))
# naive search for rational points on the quartic with larger height using numerators/denominators
from math import isqrt, gcd
def is_sq(n): return n >= 0 and isqrt(n)**2 == n
hits=[]
H=1500
for q in range(1,H+1):
    for p in range(-H, H+1):
        if gcd(abs(p),q)!=1: continue
        val = 137*p**4 - 186*p**3*q + 21*p**2*q**2 + 28*p*q**3 + 4*q**4
        if is_sq(val): hits.append((p,q))
print("rational points on quartic with height<=%d:"%H, hits)
