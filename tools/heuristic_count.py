"""
Hardy-Littlewood heuristic for  a^6+b^6+c^6+d^6+e^6 = f^6  (positive, unordered).
Expected count of solutions with f <= F  ~  S * rho_inf * log F, where
  rho_inf = (5/6) * Gamma(7/6)^5 / Gamma(11/6) / 5!   (real density per f, unordered)
  S       = prod_p sigma_p,  sigma_p = p^k * #{(a..e,f) mod p^k : sum = 0} / p^{6k}
computed exactly at prime powers large enough to stabilise (2^7, 3^7, 7^3, p^1 else).
Same for k=5 (a^5+b^5+c^5+d^5=e^5) as calibration: three solutions known below 1e5.
"""
from math import gamma, factorial, log
from collections import Counter
import sys

def sigma(p, k, npow, nterms):
    m = p**k
    pw = Counter(pow(x, npow, m) for x in range(m))     # distribution of x^n mod m
    # convolve nterms times
    dist = Counter({0: 1})
    for _ in range(nterms):
        nd = Counter()
        for s, c in dist.items():
            for v, w in pw.items():
                nd[(s+v) % m] += c*w
        dist = nd
    # count (sum, f) with sum == f^n
    cnt = sum(dist[v]*w for v, w in pw.items())
    return cnt * m / m**(nterms+1)

def heuristic(npow, nterms, F):
    primes = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113]
    S = 1.0
    detail = {}
    for p in primes:
        if p == 2: k = 7
        elif p == 3: k = 6
        elif p == 5 and npow == 5: k = 3
        elif p == 7 and npow == 6: k = 3
        else: k = 1
        s = sigma(p, k, npow, nterms)
        detail[p] = s
        S *= s
    rho = (nterms/npow) * gamma(1+1/npow)**nterms / gamma(1+nterms/npow) / factorial(nterms)
    # per-f density is rho * f^(nterms-npow) = rho / f ; sum over f<=F ~ rho*log F
    return S, rho, detail, S*rho*log(F)

for (n, t) in [(5,4),(6,5)]:
    S, rho, detail, E = heuristic(n, t, 1e6)
    print(f"k={n}: singular series S={S:.4e}  real density rho={rho:.4e}  expected count(f<=1e6)={E:.3e}  count(f<=1e100)={S*rho*log(1e100):.3e}")
    print("   local factors:", {p: round(v,4) for p, v in detail.items() if abs(v-1) > 0.05})
