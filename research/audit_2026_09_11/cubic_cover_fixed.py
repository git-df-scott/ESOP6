"""Audit replacement for cube_ansatz_2026_09_10/cubic_cover_through_point.py.

BUG FIXED (gauge chart).  The reparametrisation group fixing t = 0 acts on the
linear coefficients of the five weight-1 forms by

        g_i1  ->  lam * g_i1 + 2 c x_i           (i = 1,2,3,4,6)

so the quantity  W_jk = x_j g_k1 - x_k g_j1  satisfies  W -> lam W : it is a
GAUGE INVARIANT.  The old chart g11 = 1, g61 = 0 is reachable only when
W_16 = x1 g61 - x6 g11 != 0, hence every curve through the seed with W_16 = 0 --
an intrinsic codimension-1 condition, not a degeneracy -- was invisible to the
old solver.  Here the gauge is fixed by two RANDOM rational functionals
alpha.g = 1, beta.g = 0, solved by substitution so the system stays 12 x 12 with
the same conditioning.  The only curves lost are those with g_1 parallel to x,
i.e. the doubly-traversed-line parametrisations, which are degenerate.

Also: mpmath polish to 40 digits + continued-fraction rational recognition with a
denominator bound up to 10^12, so a rational curve with denominators in the
10^3..10^12 range cannot be silently rejected by a 1e-10 double-precision test.
"""
import numpy as np, sys, json, math, mpmath as mp
from fractions import Fraction

NEQ = 13

def gauge_chart(x4, x6, alpha, beta):
    """Return (p, Q) with g_1 = p + Q.u  (u in Q^3) the solutions of
    alpha.g = 1, beta.g = 0, as Fractions.  x = (x1..x4,x6) only used for report."""
    import sympy as sp
    a = sp.Matrix([[sp.Rational(v) for v in alpha], [sp.Rational(v) for v in beta]])
    b = sp.Matrix([1, 0])
    sol = a.solve_least_squares(b) if False else None
    # particular solution + nullspace
    aug = a.pinv()*b
    p = [sp.nsimplify(v) for v in aug]
    ns = a.nullspace()
    assert len(ns) == 3, len(ns)
    Q = [[sp.Rational(ns[j][i]) for j in range(3)] for i in range(5)]
    p = [sp.Rational(v) for v in p]
    # sanity
    for row, rhs in ((alpha, 1), (beta, 0)):
        assert sum(sp.Rational(row[i])*p[i] for i in range(5)) == rhs
        for j in range(3):
            assert sum(sp.Rational(row[i])*Q[i][j] for i in range(5)) == 0
    return [float(v) for v in p], [[float(Q[i][j]) for j in range(3)] for i in range(5)], p, Q

class CC:
    def __init__(s, x, x6, S, alpha, beta):
        s.norm = float(x6)
        s.xe = [Fraction(int(v), int(x6)) for v in x]
        s.s0e = Fraction(int(S), int(x6)**3)
        s.x = [float(v) for v in s.xe]
        s.x6 = 1.0
        s.s0 = float(s.s0e)
        s.pf, s.Qf, s.pe, s.Qe = gauge_chart(x, x6, alpha, beta)
        s.alpha, s.beta = alpha, beta
    def lin(s, u):
        return [s.pf[i] + s.Qf[i][0]*u[0] + s.Qf[i][1]*u[1] + s.Qf[i][2]*u[2] for i in range(5)]
    def res(s, v):
        def pw(a, k):
            r = np.array([1.0+0j])
            for _ in range(k): r = np.convolve(r, a)
            return r
        def pad(a, n=NEQ):
            r = np.zeros(n, dtype=complex); m = min(n, len(a)); r[:m] = a[:m]; return r
        u = v[0:3]; q = v[3:8]; k = v[8:11]; m1 = v[11]
        L = [s.pf[i] + s.Qf[i][0]*u[0] + s.Qf[i][1]*u[1] + s.Qf[i][2]*u[2] for i in range(5)]
        G = [np.array([s.x[i], L[i], q[i]]) for i in range(4)]
        G6 = np.array([1.0+0j, L[4], q[4]])
        K = np.array([s.s0+0j, k[0], k[1], k[2]])
        M = np.array([1.0+0j, m1])
        r = sum(pad(pw(Gi, 6)) for Gi in G) + pad(np.convolve(pw(K, 2), pw(M, 6))) - pad(pw(G6, 6))
        return r[1:]
    def jac(s, v, h=1e-7):
        r0 = s.res(v); J = np.zeros((12, 12), dtype=complex)
        for i in range(12):
            dv = np.zeros(12, dtype=complex); dv[i] = h
            J[:, i] = (s.res(v+dv)-r0)/h
        return J
    def scale(s, v): return 1+max(1.0, abs(s.s0)**(1/3.), np.abs(v).max())**6
    def newton(s, v, iters=150):
        for _ in range(iters):
            r = s.res(v)
            if np.linalg.norm(r) < 1e-13*s.scale(v): return v, np.linalg.norm(r)/s.scale(v)
            try: dv = np.linalg.solve(s.jac(v), -r)
            except np.linalg.LinAlgError: return v, 1.0
            lam = 1.0
            while lam > 1e-3 and np.linalg.norm(s.res(v+lam*dv)) > np.linalg.norm(r): lam /= 2
            v = v+lam*dv
            if np.abs(v).max() > 1e8: return v, 1.0
        return v, np.linalg.norm(s.res(v))/s.scale(v)
    def solve_all(s, starts, seed=0, real=False):
        rng = np.random.default_rng(seed); sols = []
        for _ in range(starts):
            a = rng.standard_normal(12)
            if not real: a = a+1j*rng.standard_normal(12)
            v = (10**rng.uniform(-1.5, 1.5))*a.astype(complex)
            v, nr = s.newton(v)
            if nr < 1e-13:
                J = s.jac(v); sv = np.linalg.svd(J, compute_uv=False)
                if sv[-1]/sv[0] < 1e-10: continue
                if not any(np.linalg.norm(v-w) < 1e-6*(1+np.linalg.norm(w)) for w in sols): sols.append(v)
        return sols
    # ---- high precision ----
    def res_mp(s, v):
        def pw(a, k):
            r = [mp.mpc(1)]
            for _ in range(k):
                o = [mp.mpc(0)]*(len(r)+len(a)-1)
                for i, ai in enumerate(r):
                    for j, bj in enumerate(a): o[i+j] += ai*bj
                r = o[:NEQ]
            return r
        def mul(a, b):
            o = [mp.mpc(0)]*min(len(a)+len(b)-1, NEQ)
            for i, ai in enumerate(a):
                for j, bj in enumerate(b):
                    if i+j < len(o): o[i+j] += ai*bj
            return o
        u = v[0:3]; q = v[3:8]; k = v[8:11]; m1 = v[11]
        P = [mp.mpf(s.pe[i].p)/s.pe[i].q if hasattr(s.pe[i], 'p') else mp.mpf(s.pe[i]) for i in range(5)]
        P = [mp.mpf(str(s.pe[i])) for i in range(5)]
        Qm = [[mp.mpf(str(s.Qe[i][j])) for j in range(3)] for i in range(5)]
        L = [P[i]+Qm[i][0]*u[0]+Qm[i][1]*u[1]+Qm[i][2]*u[2] for i in range(5)]
        X = [mp.mpf(s.xe[i].numerator)/s.xe[i].denominator for i in range(4)]
        s0 = mp.mpf(s.s0e.numerator)/s.s0e.denominator
        tot = [mp.mpc(0)]*NEQ
        for i in range(4):
            for j, c in enumerate(pw([X[i], L[i], q[i]], 6)): tot[j] += c
        for j, c in enumerate(mul(pw([s0, k[0], k[1], k[2]], 2), pw([mp.mpc(1), m1], 6))): tot[j] += c
        for j, c in enumerate(pw([mp.mpc(1), L[4], q[4]], 6)): tot[j] -= c
        return tot[1:]
    def polish(s, v, dps=40):
        mp.mp.dps = dps
        w = [mp.mpc(z) for z in v]
        h = mp.mpf(10)**(-dps//2)
        for _ in range(50):
            r = s.res_mp(w); nr = max(abs(z) for z in r)
            J = mp.matrix(12, 12)
            for i in range(12):
                uu = list(w); uu[i] = uu[i]+h
                r1 = s.res_mp(uu)
                for j in range(12): J[j, i] = (r1[j]-r[j])/h
            try: dv = mp.lu_solve(J, mp.matrix([-z for z in r]))
            except Exception: return w, nr
            w = [a+dv[i] for i, a in enumerate(w)]
            if nr < mp.mpf(10)**(-dps+10): break
        return w, max(abs(z) for z in s.res_mp(w))

def ratrec(z, maxden=10**12, tol=None):
    if abs(mp.im(z)) > mp.mpf(10)**(-18)*(1+abs(z)): return None
    x = mp.re(z)
    if tol is None: tol = mp.mpf(10)**(-22)
    num0, den0 = 1, 0
    ai = int(mp.floor(x)); num1, den1 = ai, 1
    frac = x-ai
    for _ in range(80):
        if den1 and den1 <= maxden:
            if abs(mp.mpf(num1)/den1-x) < tol*(1+abs(x)): return Fraction(num1, den1)
        if frac == 0: break
        a = 1/frac; ai = int(mp.floor(a)); frac = a-ai
        num0, den0, num1, den1 = num1, den1, ai*num1+num0, ai*den1+den0
        if den1 > maxden: break
    return None

def reconstruct_exact(s, fr):
    """fr = 12 Fractions in the chart; return exact forms and verify the identity."""
    import sympy as sp
    t = sp.symbols('t')
    u = fr[0:3]; q = fr[3:8]; k = fr[8:11]; m1 = fr[11]
    L = [s.pe[i]+s.Qe[i][0]*u[0]+s.Qe[i][1]*u[1]+s.Qe[i][2]*u[2] for i in range(5)]
    X = [sp.Rational(s.xe[i].numerator, s.xe[i].denominator) for i in range(4)]
    s0 = sp.Rational(s.s0e.numerator, s.s0e.denominator)
    G = [X[i]+sp.Rational(L[i])*t+sp.Rational(q[i])*t**2 for i in range(4)]
    G6 = 1+sp.Rational(L[4])*t+sp.Rational(q[4])*t**2
    K = s0+sp.Rational(k[0])*t+sp.Rational(k[1])*t**2+sp.Rational(k[2])*t**3
    M = 1+sp.Rational(m1)*t
    expr = sp.expand(sum(Gi**6 for Gi in G)+K**2*M**6-G6**6)
    return expr == 0, [str(g) for g in G], str(G6), str(K), str(M)

def run(seedfile, starts, tag, real, limit, skip, rng_seed=12345):
    seeds = json.load(open(seedfile))
    rng = np.random.default_rng(rng_seed)
    out = []; seen = set(); done = 0
    for p in seeds:
        x = [int(v) for v in p['x']]; x6 = int(p['x6']); S = int(p['S'])
        if sum(v**6 for v in x)+S*S != x6**6:
            print('INVALID EXACT Y2 SEED; skipped', p, flush=True); continue
        key = (tuple(sorted(x)), x6)
        if key in seen: continue
        seen.add(key)
        if skip > 0: skip -= 1; continue
        while True:
            alpha = [int(v) for v in rng.integers(-9, 10, 5)]
            beta = [int(v) for v in rng.integers(-9, 10, 5)]
            import sympy as sp
            if sp.Matrix([alpha, beta]).rank() == 2 and any(alpha): break
        s = CC(x, x6, S, alpha, beta)
        sols = s.solve_all(starts, real=real)
        nreal = sum(1 for v in sols if np.abs(v.imag).max() < 1e-6)
        # intrinsic invariant W16 = x1*g61 - x6*g11 of each found curve
        W = []
        for v in sols:
            L = s.lin(v[0:3].real if np.abs(v.imag).max() < 1e-9 else v[0:3])
            W.append(abs(s.x[0]*L[4]-1.0*L[0]))
        rat = 0
        for v in sols:
            w, nr = s.polish(v)
            if nr > mp.mpf(10)**(-25): continue
            fr = [ratrec(z) for z in w]
            if all(f is not None for f in fr):
                ok, G, G6, K, M = reconstruct_exact(s, fr)
                if ok:
                    rat += 1
                    print("!!!! RATIONAL GENUS-1 CURVE through", p, "G=", G, "G6=", G6, "K=", K, "M=", M, flush=True)
                    out.append(dict(point=p, G=G, G6=G6, K=K, M=M))
        print("seed", key, "S=", S, "alpha", alpha, "beta", beta, ": curves", len(sols),
              "real", nreal, "minW %.3g" % (min(W) if W else -1), "rational", rat, flush=True)
        done += 1
        if limit and done >= limit: break
    json.dump(out, open("/home/user/ESOP6/research/audit_2026_09_11/cubic_cover_fixed_hits_%s.json" % tag, "w"), indent=1)

if __name__ == "__main__":
    args = sys.argv[1:]
    real = '--real' in args
    def opt(n, d): return int(args[args.index(n)+1]) if n in args else d
    run(args[0], int(args[1]), args[2], real, opt('--limit', 0), opt('--skip', 0), opt('--rngseed', 12345))
