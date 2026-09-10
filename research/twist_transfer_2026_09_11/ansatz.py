#!/usr/bin/env python3
"""Direct twisted-curve ansatz on the ESOP6 fourfold.

MASTER REDUCTION (derived in REPORT.md).  Let sigma be a rational involution of P^1
with irrational real fixed points; in the quotient coordinates (x,y) = (E1,E2) the
sigma-invariant forms of degree 2n in (s,t) are the degree-n forms in (x,y), and the
anti-invariant ones are O*(degree n-1), where O^2 = Q(x,y) = x^2 - lambda*y^2,
lambda = 4/D > 0 non-square.  A twisted curve with twist set T gives a counterexample
iff the following holds identically in Q[x,y], with R = -Q = lambda*y^2 - x^2:

  |T|=2 (source (3,3)):   U1^6+U2^6+U3^6 - U4^6 = Q^3 ( P1^6 + P2^6 )
  |T|=1 (source (6,2,4)): U1^6+U2^6+U3^6+U4^6 - U5^6 = Q^3 ( P1^6 )

with U_i of degree n and P_j of degree n-1; and then at ANY rational (x0,y0) with
R(x0,y0) = w^2 a rational square, (U1,U2,U3, w*P1, w*P2 ; U4) is a rational point of
  x1^6+x2^6+x3^6+x4^6+x5^6 = x6^6.
(The conic w^2+x^2 = lambda*y^2 has infinitely many rational points whenever it has one,
e.g. lambda=2: (w,x,y)=(1,1,1).)

Over R, lambda is pure gauge: Q is indefinite so Q ~ x*y after a real GL2 change.
We therefore solve the REAL problem   sum eps_i U_i^6 = x^3 y^3 * sum P_j^6 .
"""
import numpy as np

def pw6(c):
    """c = coeffs of a binary form (dehomogenised); return coeffs of c^6."""
    r = np.array([1.0])
    for _ in range(6):
        r = np.convolve(r, c)
    return r

def residual(z, n, k):
    """k = |T| (number of twisted coords).  4-k+... see below.
    Layout: (6-k) forms U of degree n, then k forms P of degree n-1."""
    nU, nP = 6 - k, k
    U = [z[i*(n+1):(i+1)*(n+1)] for i in range(nU)]
    off = nU*(n+1)
    P = [z[off+i*n:off+(i+1)*n] for i in range(nP)]
    tot = np.zeros(6*n+1)
    for i, u in enumerate(U):
        s = -1.0 if i == nU-1 else 1.0      # last U is the lone x6 coordinate
        tot = tot + s*pw6(u)
    rhs = np.zeros(6*(n-1)+1)
    for p in P:
        rhs = rhs + pw6(p)
    # multiply rhs by x^3 y^3 : in dehomogenised coords (t = x/y) that is t^3 * (...)
    # degrees: tot has degree 6n in t; x^3y^3 * deg-6(n-1) form -> t^3 * poly(deg 6n-6)
    full = np.zeros(6*n+1)
    full[3:3+len(rhs)] += rhs
    res = tot - full
    # gauge: fix overall scale by normalising the P block
    res = np.append(res, np.sum(np.concatenate(P)**2) - 1.0)
    return res

def jac(z, n, k, h=1e-7):
    f0 = residual(z, n, k)
    J = np.zeros((len(f0), len(z)))
    for i in range(len(z)):
        zz = z.copy(); zz[i] += h
        J[:, i] = (residual(zz, n, k) - f0)/h
    return J

def solve(n, k, tries=4000, seed=0, tol=1e-11):
    rng = np.random.default_rng(seed)
    nvar = (6-k)*(n+1) + k*n
    found = []
    for it in range(tries):
        z = rng.normal(size=nvar)*rng.choice([0.3,1.0,3.0])
        for _ in range(400):
            f = residual(z, n, k)
            if not np.all(np.isfinite(f)) or np.max(np.abs(f)) > 1e14:
                break
            J = jac(z, n, k)
            try:
                dz = np.linalg.lstsq(J, -f, rcond=None)[0]
            except np.linalg.LinAlgError:
                break
            nrm = np.linalg.norm(dz)
            if nrm > 2.0:
                dz *= 2.0/nrm
            z = z + dz
            if np.max(np.abs(f)) < tol:
                break
        f = residual(z, n, k)
        if np.all(np.isfinite(f)) and np.max(np.abs(f)) < 1e-9:
            found.append(z.copy())
    return found

def report(n, k, sols):
    print(f'--- n={n}, |T|={k}:  {len(sols)} converged real solutions out of trials')
    nU, nP = 6-k, k
    good = []
    for z in sols:
        U = [z[i*(n+1):(i+1)*(n+1)] for i in range(nU)]
        off = nU*(n+1)
        P = [z[off+i*n:off+(i+1)*n] for i in range(nP)]
        if min(np.max(np.abs(p)) for p in P) < 1e-6:
            continue                                  # some P identically 0
        if min(np.max(np.abs(u)) for u in U) < 1e-6:
            continue                                  # some U identically 0 (lower stratum)
        # degenerate: two U equal up to sign
        deg = False
        for i in range(nU):
            for j in range(i+1, nU):
                if min(np.max(np.abs(U[i]-U[j])), np.max(np.abs(U[i]+U[j]))) < 1e-6:
                    deg = True
        if deg:
            continue
        good.append(z)
    print(f'    non-degenerate (all P,U nonzero, no repeated U): {len(good)}')
    for z in good[:4]:
        print('      ', np.round(z, 6))
    return good

if __name__ == '__main__':
    import sys
    for (n, k) in [(1,2),(1,1),(2,2),(2,1),(3,2),(3,1)]:
        sols = solve(n, k, tries=int(sys.argv[1]) if len(sys.argv)>1 else 1500, seed=n*10+k)
        report(n, k, sols)
