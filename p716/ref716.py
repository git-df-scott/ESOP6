#!/usr/bin/env python3
"""Independent brute-force 4+3 reference for (7,1,6), pure Python big ints.

  ref716.py FMIN FMAX      exhaustive search over FMIN < f <= FMAX
Prints:  REF fmin fmax leaves solutions   and one SOLUTION line per hit.

Written from the algorithm description only (nested largest-part windows plus
a dict/set of 3-sums); it shares no code with engine716.c.
"""
import sys


def iroot(n, k=7):
    if n <= 0:
        return 0
    x = 1 << ((n.bit_length() + k - 1) // k)
    while True:
        y = ((k - 1) * x + n // x ** (k - 1)) // k
        if y >= x:
            return x
        x = y


def lo_root(X, k):
    """smallest n >= 1 with k*n^7 >= X"""
    if X <= k:
        return 1
    t = iroot((X + k - 1) // k)
    while k * t ** 7 < X:
        t += 1
    while t > 1 and k * (t - 1) ** 7 >= X:
        t -= 1
    return max(1, t)


def main():
    fmin, fmax = int(sys.argv[1]), int(sys.argv[2])
    P = [x ** 7 for x in range(fmax + 2)]
    # table of all d^7+e^7+g^7 with fmax >= d >= e >= g >= 1
    tab = set()
    for d in range(1, fmax + 1):
        pd = P[d]
        for e in range(1, d + 1):
            pde = pd + P[e]
            for g in range(1, e + 1):
                tab.add(pde + P[g])
    leaves = 0
    sols = 0
    for f in range(fmin + 1, fmax + 1):
        F7 = P[f]
        amin = lo_root(F7, 6)
        amax = min(f - 1, iroot(F7))
        for a in range(amax, amin - 1, -1):
            R1 = F7 - P[a]
            if R1 < 5:
                continue
            bhi = min(a, iroot(R1))
            blo = lo_root(R1, 5)
            for b in range(bhi, blo - 1, -1):
                R2 = R1 - P[b]
                if R2 < 4:
                    continue
                chi = min(b, iroot(R2))
                clo = lo_root(R2, 4)
                for c in range(chi, clo - 1, -1):
                    R3 = R2 - P[c]
                    if R3 < 3 or R3 > 3 * P[c]:
                        continue
                    leaves += 1
                    if R3 not in tab:
                        continue
                    # exact: find d<=c
                    dhi = min(c, iroot(R3))
                    dlo = lo_root(R3, 3)
                    for d in range(dhi, dlo - 1, -1):
                        Rp = R3 - P[d]
                        if Rp < 2:
                            continue
                        ehi = min(d, iroot(Rp))
                        elo = lo_root(Rp, 2)
                        for e in range(ehi, elo - 1, -1):
                            rest = Rp - P[e]
                            if rest < 1:
                                continue
                            g = iroot(rest)
                            if 1 <= g <= e and P[g] == rest:
                                sols += 1
                                print("SOLUTION %d %d %d %d %d %d %d" % (f, a, b, c, d, e, g))
    print("REF %d %d %d %d" % (fmin, fmax, leaves, sols))


if __name__ == "__main__":
    main()
