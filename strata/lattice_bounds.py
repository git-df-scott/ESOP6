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

def triangle_points(b1, b2, A):
    """All lattice points (f,t) = x*b1 + y*b2 with 0 < t < f <= A.  Exact."""
    det = b1[0] * b2[1] - b1[1] * b2[0]
    assert det != 0
    # y as a rational function of (f,t):  y = (t*b1[0] - f*b1[1]) / det
    ys = []
    for (f, t) in ((0, 0), (A, 0), (A, A)):
        ys.append((t * b1[0] - f * b1[1]) / det)
    ymin, ymax = math.floor(min(ys)) - 1, math.ceil(max(ys)) + 1
    assert ymax - ymin < 10**7, "y range too wide (%d); basis not reduced?" % (ymax - ymin)
    out = []
    # three half-planes, each of the form  c*x + d*y >= e
    cons = ((b1[1],              b2[1],              1),                 # t >= 1
            (b1[0] - b1[1],      b2[0] - b2[1],      1),                 # f - t >= 1
            (-b1[0],             -b2[0],             -A))                # f <= A
    for y in range(ymin, ymax + 1):
        lo, hi = None, None
        feasible = True
        for (c, d, e) in cons:
            rhs = e - d * y
            if c == 0:
                if 0 < rhs:
                    feasible = False
                    break
            elif c > 0:
                v = ceil_div(rhs, c)
                lo = v if lo is None else max(lo, v)
            else:
                v = floor_div(rhs, c)
                hi = v if hi is None else min(hi, v)
        if not feasible or lo is None or hi is None or lo > hi:
            continue
        assert hi - lo < 10**7, "x range too wide"
        for x in range(lo, hi + 1):
            f = x * b1[0] + y * b2[0]
            t = x * b1[1] + y * b2[1]
            assert 0 < t < f <= A
            out.append((f, t))
    return out

def least_f(z, N, need_coprime=True):
    """Least f > 0 with some 0 < t < f, t == z*f (mod N), gcd(f,42)=gcd(t,42)=1."""
    b1, b2 = gauss_reduce((1, z), (0, N))
    A = isqrt(N) + 1
    while True:
        pts = triangle_points(b1, b2, A)
        good = [(f, t) for (f, t) in pts
                if (not need_coprime) or (gcd(f, 42) == 1 and gcd(t, 42) == 1)]
        if good:
            f, t = min(good)
            assert (t - z * f) % N == 0
            return f, t
        A *= 2
        assert A < 1 << 80, "no admissible point found"

# ----------------------------------------------------------------------- driver
CASES = [
    ("all four b_i divisible by 2 ", "42^6 * 2^6", [(2, 12), (3, 6), (7, 6)]),
    ("all four b_i divisible by 3 ", "42^6 * 3^6", [(2, 6), (3, 12), (7, 6)]),
    ("all four b_i divisible by 7 ", "42^6 * 7^6", [(2, 6), (3, 6), (7, 12)]),
    ("all four b_i divisible by 42", "42^12",      [(2, 12), (3, 12), (7, 12)]),
]

def main():
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

if __name__ == "__main__":
    main()
