import mpmath as mp, numpy as np
mp.mp.dps = 30
# Refine conic 0 in the symmetric family: unknowns r2, r3, z, alpha with r1 = 1
# F6: 2 + 2 r3^6 cos6a - r2^6 = 0 ; F4: 20 z^2 + 2 + r2^6 + 2 r3^6 cos4a = 0 ; F2: 32 z^4 + 32 z^2 + 2 - r2^6 + 2 r3^6 cos2a = 0
def F(v):
    r2, r3, z, a = v
    return [2 + 2*r3**6*mp.cos(6*a) - r2**6, 20*z**2 + 2 + r2**6 + 2*r3**6*mp.cos(4*a), 32*z**4 + 32*z**2 + 2 - r2**6 + 2*r3**6*mp.cos(2*a)]
v0 = [0.818952/0.167417, 0.758476/0.167417, 0.816756/0.167417, mp.radians(53.731)]
# 3 equations, 4 unknowns: fix alpha, solve for r2,r3,z
def solve_alpha(a):
    f = lambda r2, r3, z: F([r2, r3, z, a])
    return mp.findroot(f, (v0[0], v0[1], v0[2]))
sol = solve_alpha(v0[3]); r2, r3, z = sol
print("refined conic 0:", r2, r3, z, "alpha deg", mp.degrees(v0[3]))
c = mp.cos(2*v0[3]); n = r3**2; m = n**3
P = (c+1)*(4*c**2-2*c-1)
Z = -(2 + m*P)/10
print("check z^2 =", z**2, " vs Z(c,m) =", Z)
print("check Gamma: m^2 P^2 - m(6P+25c(c^2-1)) - 16 =", m**2*P**2 - m*(6*P + 25*c*(c**2-1)) - 16)
print("check A = r2^6 =", r2**6, " vs 2+2m(4c^3-3c) =", 2 + 2*m*(4*c**3-3*c))
# full identity check: sum ell_i^6 constant over theta
def total(th):
    L1 = mp.cos(th); L2 = r2*mp.sin(th); L3 = r3*mp.cos(th - v0[3]); L4 = r3*mp.cos(th + v0[3])
    return (z+L1)**6 + (z-L1)**6 + L2**6 + L3**6 + L4**6
vals = [total(mp.mpf(k)/7) for k in range(7)]
print("sum ell^6 at 7 angles:", [mp.nstr(x, 12) for x in vals])
# trace the family: alpha sweep, report c and whether m>0, Z>0, A>0
print("family trace (alpha deg, c=cos2a, n=r3^2, m=n^3, Z, A):")
for adeg in [45, 48, 50, 52, 53, 53.5, 54, 55, 56, 58, 60, 62, 65]:
    try:
        s = solve_alpha(mp.radians(adeg)); r2,r3,z = s
        if abs(mp.im(r2))>1e-20: continue
        c = mp.cos(2*mp.radians(adeg))
        print(f"  {adeg:5}  c={mp.nstr(c,8)}  n={mp.nstr(r3**2,8)}  m={mp.nstr(r3**6,8)}  Z={mp.nstr(z**2,8)}  A={mp.nstr(r2**6,8)}")
    except Exception as e:
        print("  ", adeg, "no solution:", str(e)[:40])
