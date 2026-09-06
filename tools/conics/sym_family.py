"""
Z/2-symmetric general conics on X:  ell_1 = z + L1, ell_5 = z - L1, ell_2 = L2, ell_4 = -L2, ell_3 = L3,
L_i = Re(conj(zeta_i) u), u = e^{i theta}.  Identity  sum ell_i^6 = g^6  <=>  E6, E4, E2 (complex):
  E6: 2 z1^6 + 2 z2^6 + z3^6 = 0
  E4: (20 z^2 + 2 N1) z1^4 + 2 N2 z2^4 + N3 z3^4 = 0
  E2: (32 z^4 + 32 z^2 N1 + 2 N1^2) z1^2 + 2 N2^2 z2^2 + N3^2 z3^2 = 0
Normalise zeta3 = 1 (rotation + scale). Unknowns: zeta1, zeta2 (complex), z (real) -> 5 real. 6 real equations.
"""
import mpmath as mp, random, sys
mp.mp.dps = 40
def eqs(v):
    a,b,c,d,z = v
    z1 = mp.mpc(a,b); z2 = mp.mpc(c,d); z3 = mp.mpc(1,0)
    N1 = abs(z1)**2; N2 = abs(z2)**2; N3 = mp.mpf(1)
    E6 = 2*z1**6 + 2*z2**6 + z3**6
    E4 = (20*z**2 + 2*N1)*z1**4 + 2*N2*z2**4 + N3*z3**4
    E2 = (32*z**4 + 32*z**2*N1 + 2*N1**2)*z1**2 + 2*N2**2*z2**2 + N3**2*z3**2
    return [E6.real,E6.imag,E4.real,E4.imag,E2.real,E2.imag]
def jac(v, h=mp.mpf('1e-20')):
    J = mp.matrix(6,5); f0 = eqs(v)
    for j in range(5):
        vv = list(v); vv[j] += h; f1 = eqs(vv)
        for i in range(6): J[i,j] = (f1[i]-f0[i])/h
    return J
def newton(v, it=60):
    for _ in range(it):
        f = mp.matrix(eqs(v)); J = jac(v)
        try: dx = mp.lu_solve(J.T*J, -J.T*f)
        except ZeroDivisionError: return None
        v = [v[i]+dx[i] for i in range(5)]
        if mp.norm(f) < mp.mpf('1e-30'): return v
        if mp.norm(mp.matrix(v)) > 1e4: return None
    return v if mp.norm(mp.matrix(eqs(v))) < mp.mpf('1e-25') else None
random.seed(int(sys.argv[1]) if len(sys.argv)>1 else 0)
sols = []
for k in range(int(sys.argv[2]) if len(sys.argv)>2 else 300):
    v0 = [mp.mpf(random.gauss(0,1)) for _ in range(5)]
    v = newton(v0)
    if v is None: continue
    a,b,c,d,z = v
    if abs(mp.mpc(a,b)) < 1e-6 or abs(mp.mpc(c,d)) < 1e-6 or abs(z) < 1e-8: continue  # degenerate
    key = tuple(round(float(abs(x)),6) for x in (mp.mpc(a,b), mp.mpc(c,d), z))
    if any(sum(abs(key[i]-s[0][i]) for i in range(3))<1e-5 for s in sols): continue
    J = jac(v); sv = mp.svd_r(J, compute_uv=False)
    rank = sum(1 for s in sv if s > 1e-12*sv[0])
    sols.append((key, v, rank))
    z1 = mp.mpc(a,b); z2 = mp.mpc(c,d)
    N1 = abs(z1)**2; N2 = abs(z2)**2
    g6 = 2*z**6 + 15*z**4*N1 + mp.mpf(45)/4*z**2*N1**2 + mp.mpf(5)/8*N1**3 + mp.mpf(5)/8*N2**3 + mp.mpf(5)/16
    print(f"sol rank={rank}: zeta1={mp.nstr(z1,15)} zeta2={mp.nstr(z2,15)} z={mp.nstr(z,15)} N1={mp.nstr(N1,12)} N2={mp.nstr(N2,12)} g^6={mp.nstr(g6,12)}")
    for name,val in [("N1",N1),("N2",N2),("z^2",z**2),("g^6",g6),("Re z1^2",(z1**2).real),("Im z1^2",(z1**2).imag)]:
        p = mp.findpoly(val, 8, maxcoeff=10**6)
        print(f"    {name}: minpoly={p}")
print("distinct nondegenerate solutions:", len(sols))
