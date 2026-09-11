#!/usr/bin/env python3
"""
lattice_bounds.py -- rigorous lower bounds on f for the DEEP class-1 strata.

Setting.  In a class-1 primitive solution a^6+...+e^6 = f^6 the four non-exempt
terms are 42*b_i and  f^6 - t^6 = 42^6 * (b1^6+b2^6+b3^6+b4^6).  If in addition
all four b_i are divisible by an extra prime p in {2,3,7} (resp. by 42), then

        f^6 == t^6   (mod 42^6 * p^6)      (resp.  mod 42^12),

so t == zeta*f (mod N) for one of the 144 sixth roots of unity zeta mod N,
with 0 < t < f and gcd(f,42) = gcd(t,42) = 1.

For a fixed zeta the admissible (f,t) form the lattice

        L = {(f,t) in Z^2 : t == zeta*f (mod N)},   basis (1,zeta), (0,N),
        det L = N.

The least f for which L meets the open cone 0 < t < f is therefore a provable
LOWER BOUND on f in that stratum, and the minimum over zeta != 1 is the bound
for the stratum as a whole.  (zeta = 1 forces t == f mod N, i.e. f > N.)

Method.  Gauss(-Lagrange)-reduce the basis, then enumerate the lattice points of
the TRIANGLE {0 < t < f <= A} exactly (integer arithmetic only), doubling A until
the triangle contains an admissible point.  Enumerating the triangle rather than a
ball matters: for zeta close to -1 the lattice has a very short vector along
t = -f, which lies outside the cone, so a ball enumeration would walk O(N) points
while the triangle enumeration walks O(A^2/N + #y-layers).

Runtime: a few seconds for all 4*144 lattices.
"""
import sys, math
from math import gcd, isqrt

# ----------------------------------------------------------------- roots of unity
def sixth_roots_prime_power(p, e):
    """All z in (Z/p^e)^* with z^6 == 1."""
    M = p**e
    if p == 2:
        assert e >= 3
        # the unit group is a 2-group, so z^6 = 1 <=> z^2 = 1: exactly 4 solutions
        rs = sorted({1, M - 1, M // 2 - 1, M // 2 + 1})
    else:
        g = {3: 2, 7: 3}[p]          # 2 is a primitive root mod 3^e, 3 mod 7^e
        order = (p - 1) * p**(e - 1)
        assert pow(g, order // (p - 1), M) != 1 and pow(g, order // p, M) != 1
        step = order // 6
        rs = sorted({pow(g, j * step, M) for j in range(6)})
    for z in rs:
        assert pow(z, 6, M) == 1 and gcd(z, p) == 1
    assert len(rs) == (4 if p == 2 else 6)
    return rs

def crt_pair(a, m, b, n):
    return (a + m * ((b - a) * pow(m, -1, n) % n)) % (m * n)

def all_sixth_roots(factors):
    """factors: [(p,e), ...]  ->  (sorted list of the 144 roots, N)"""
    zs, M = [0], 1
    for p, e in factors:
        pe = p**e
        rs = sixth_roots_prime_power(p, e)
        zs = [crt_pair(z, M, r, pe) for z in zs for r in rs]
        M *= pe
    zs.sort()
    assert len(zs) == 144, len(zs)
    for z in zs:
        assert pow(z, 6, M) == 1
    return zs, M

# ----------------------------------------------------------------- lattice tools
def gauss_reduce(b1, b2):
    """Lagrange reduction of a 2-D integer basis; exact."""
    n = lambda v: v[0] * v[0] + v[1] * v[1]
    d = lambda u, v: u[0] * v[0] + u[1] * v[1]
    while True:
        if n(b1) > n(b2):
            b1, b2 = b2, b1
        n1 = n(b1)
        if n1 == 0:
            return b1, b2
        # mu = nearest integer to <b1,b2>/n1
        num, den = d(b1, b2), n1
        mu = (2 * num + den) // (2 * den)
        if mu == 0:
            return b1, b2
        b2 = (b2[0] - mu * b1[0], b2[1] - mu * b1[1])

def floor_div(a, b):
    """floor(a/b) for any nonzero integer b (Python's // already floors)."""
    return a // b

def ceil_div(a, b):
    return -((-a) // b)

def _x_interval(b1, b2, A, y):
    """Integer x-interval for which x*b1 + y*b2 lies in {0 < t < f <= A}."""
    cons = ((b1[1],         b2[1],         1),      # t >= 1
            (b1[0] - b1[1], b2[0] - b2[1], 1),      # f - t >= 1
            (-b1[0],        -b2[0],        -A))     # f <= A
    lo = hi = None
    for (c, d, e) in cons:
        rhs = e - d * y
        if c == 0:
            if rhs > 0:
                return None
        elif c > 0:
            v = ceil_div(rhs, c)
            lo = v if lo is None else max(lo, v)
        else:
            v = floor_div(rhs, c)
            hi = v if hi is None else min(hi, v)
    if lo is None or hi is None or lo > hi:
        return None
    return lo, hi

def _y_range(b1, b2, A):
    det = b1[0] * b2[1] - b1[1] * b2[0]
    assert det != 0
    ys = [(t * b1[0] - f * b1[1]) / det for (f, t) in ((0, 0), (A, 0), (A, A))]
    ymin, ymax = math.floor(min(ys)) - 1, math.ceil(max(ys)) + 1
    assert ymax - ymin < 10**6, "y range too wide (%d)" % (ymax - ymin)
    return ymin, ymax

def triangle_points(b1, b2, A, cap=10**6):
    """ALL lattice points (f,t) = x*b1 + y*b2 with 0 < t < f <= A.  Exact.
    Used by the self-test; the production search uses least_f_in_triangle."""
    out = []
    ymin, ymax = _y_range(b1, b2, A)
    for y in range(ymin, ymax + 1):
        iv = _x_interval(b1, b2, A, y)
        if iv is None:
            continue
        lo, hi = iv
        assert len(out) + (hi - lo + 1) <= cap, "too many points"
        for x in range(lo, hi + 1):
            f = x * b1[0] + y * b2[0]
            t = x * b1[1] + y * b2[1]
            assert 0 < t < f <= A
            out.append((f, t))
    return out

def least_f_in_triangle(b1, b2, A):
    """Least f with 0 < t < f <= A, gcd(f,42)=gcd(t,42)=1, over the lattice.

    Within one y-layer f is linear in x and the coprimality predicate depends
    only on x mod 42, so at most 42 x-values from the f-minimising end of the
    interval have to be inspected.  That keeps the cost O(42 * #layers) even
    when the triangle holds astronomically many lattice points."""
    best = None
    ymin, ymax = _y_range(b1, b2, A)
    for y in range(ymin, ymax + 1):
        iv = _x_interval(b1, b2, A, y)
        if iv is None:
            continue
        lo, hi = iv
        step = 1 if b1[0] >= 0 else -1
        start = lo if step == 1 else hi
        n = min(42, hi - lo + 1)
        for k in range(n):
            x = start + step * k
            f = x * b1[0] + y * b2[0]
            t = x * b1[1] + y * b2[1]
            assert 0 < t < f <= A
            if gcd(f, 42) == 1 and gcd(t, 42) == 1:
                if best is None or f < best[0]:
                    best = (f, t)
                break
    return best

def least_f(z, N):
    """Least f > 0 with some 0 < t < f, t == z*f (mod N), gcd(f,42)=gcd(t,42)=1."""
    b1, b2 = gauss_reduce((1, z), (0, N))
    A = isqrt(N) + 1
    while True:
        got = least_f_in_triangle(b1, b2, A)
        if got is not None:
            f, t = got
            assert (t - z * f) % N == 0 and 0 < t < f
            return f, t
        A *= 2
        assert A < 1 << 90, "no admissible point found"

def selftest():
    """triangle_points and least_f_in_triangle vs direct brute force."""
    import random
    random.seed(1)
    n = 0
    for _ in range(400):
        N = random.randrange(50, 4000)
        z = random.randrange(1, N)
        if gcd(z, N) != 1:
            continue
        A = random.randrange(5, 900)
        bf = []
        for f in range(1, A + 1):
            t = (z * f) % N
            if t == 0:
                t = N
            while t < f:
                bf.append((f, t))
                t += N
        bf.sort()
        b1, b2 = gauss_reduce((1, z), (0, N))
        assert bf == sorted(triangle_points(b1, b2, A)), (N, z, A)
        cop = [(f, t) for (f, t) in bf if gcd(f, 42) == 1 and gcd(t, 42) == 1]
        want = min(cop) if cop else None
        assert least_f_in_triangle(b1, b2, A) == want, (N, z, A)
        n += 1
    print("selftest: triangle enumeration and least-f agree with brute force on %d random lattices" % n)

# ----------------------------------------------------------------------- driver
CASES = [
    ("all four b_i divisible by 2 ", "42^6 * 2^6", [(2, 12), (3, 6), (7, 6)]),
    ("all four b_i divisible by 3 ", "42^6 * 3^6", [(2, 6), (3, 12), (7, 6)]),
    ("all four b_i divisible by 7 ", "42^6 * 7^6", [(2, 6), (3, 6), (7, 12)]),
    ("all four b_i divisible by 42", "42^12",      [(2, 12), (3, 12), (7, 12)]),
]

def main():
    selftest()
    results = {}
    for label, nname, factors in CASES:
        zs, N = all_sixth_roots(factors)
        rows = []
        for z in zs:
            if z == 1:
                continue
            f, t = least_f(z, N)
            rows.append((f, t, z))
        rows.sort()
        results[nname] = (N, rows)
        fmin, tmin, zmin = rows[0]
        print("=" * 78)
        print("%s   N = %s = %d" % (label, nname, N))
        print("  sixth roots of unity mod N : %d   (zeta = 1 excluded: it forces f > N)" % len(zs))
        print("  LOWER BOUND  min over zeta != 1 of least f :  %d" % fmin)
        print("      attained at zeta = %d   with t = %d   (t/f = %.9f)" % (zmin, tmin, tmin / fmin))
        print("  next four    : %s" % ", ".join(str(r[0]) for r in rows[1:5]))
        print("  median / max : %d / %d" % (rows[len(rows) // 2][0], rows[-1][0]))
        sys.stdout.flush()
    print("=" * 78)
    print("SUMMARY (provable lower bounds on f for class-1 solutions in the deep strata)")
    for label, nname, _ in CASES:
        N, rows = results[nname]
        print("  %-30s N=%-22d  f >= %d" % (label, N, rows[0][0]))
    return results

def verify_7_by_scan(claimed=None):
    """Independent double-check of the N = 42^6 * 7^6 bound by a DIFFERENT route:
    scan every f up to the claimed bound and test all 144 roots directly.

    No lattice, no reduction, no enumeration -- just t = zeta*f mod N.  Since
    N = 6.458e14 is far larger than the claimed bound, the residue in [0,N) is the
    only t that can possibly satisfy t < f.  gcd(t,42) = gcd(f,42) automatically
    (42 | N and zeta is a unit mod 42), so only gcd(f,42)=1 has to be imposed.
    Vectorised with numpy in chunks small enough that zeta*k stays below 2^63."""
    import numpy as np
    zs, N = all_sixth_roots([(2, 6), (3, 6), (7, 12)])
    if claimed is None:
        claimed = min(least_f(z, N)[0] for z in zs if z != 1)
    print("verify_7_by_scan: scanning f = 1 .. %d against all %d roots" % (claimed, len(zs)))
    C = 10000
    assert C * N < (1 << 63), "chunk too large for int64"
    hits = []
    f0 = 1
    while f0 <= claimed:
        c = min(C, claimed - f0 + 1)
        fs = np.arange(f0, f0 + c, dtype=np.int64)
        good_f = (fs % 2 != 0) & (fs % 3 != 0) & (fs % 7 != 0)
        if good_f.any():
            k = np.arange(c, dtype=np.int64)
            for z in zs:
                if z == 1:
                    continue
                t0 = (z % N) * (f0 % N) % N
                t = (t0 + k * (z % N)) % N
                m = good_f & (t > 0) & (t < fs)
                if m.any():
                    for idx in np.nonzero(m)[0]:
                        hits.append((int(fs[idx]), int(t[idx]), z))
        f0 += c
    hits.sort()
    if hits:
        f, t, z = hits[0]
        print("  scan found least f = %d (t = %d, zeta = %d); %d hits at or below the bound"
              % (f, t, z, len(hits)))
        okk = (f == claimed)
    else:
        print("  scan found NO admissible (f,t) with f <= %d" % claimed)
        okk = False
    print(("PASS  " if okk else "FAIL  ") +
          "independent scan reproduces the 42^6*7^6 bound f >= %d" % claimed)
    return okk

if __name__ == "__main__":
    res = main()
    if "--verify7" in sys.argv:
        print()
        verify_7_by_scan()
