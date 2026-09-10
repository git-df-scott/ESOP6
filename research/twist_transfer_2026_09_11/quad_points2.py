#!/usr/bin/env python3
"""Meet-in-the-middle search for A^6+B^6+C^6 = F^6 over Q(sqrt D), D>0 non-square.

Coordinates u+v*sqrt(D) with |u|,|v| <= H (NO per-coordinate primitivity filter:
projective scaling is global, not per coordinate).  (u+v r)^6 = A + B r exactly in Z.
Hash meet-in-the-middle with exact verification of every candidate.
"""
import sys
from math import isqrt
import numpy as np

MASK = (1 << 21) - 1


def sixth_arrays(H, D):
    us, vs = [], []
    for v in range(0, H + 1):
        for u in range(-H, H + 1):
            if v == 0 and u <= 0:
                continue
            us.append(u); vs.append(v)
    u = np.array(us, dtype=object); v = np.array(vs, dtype=object)
    a = np.ones(len(us), dtype=object); b = np.zeros(len(us), dtype=object)
    for _ in range(6):
        a, b = a * u + b * v * D, a * v + b * u
    A = np.array([int(x) for x in a], dtype=np.int64)
    B = np.array([int(x) for x in b], dtype=np.int64)
    return np.array(us), np.array(vs), A, B


def key(A, B):
    return ((A & MASK) << 22) ^ (B & MASK)


def search(D, H):
    us, vs, A, B = sixth_arrays(H, D)
    N = len(A)
    # build 2-sums  A_i+A_j, i<=j
    ii, jj = np.triu_indices(N)
    SA = A[ii] + A[jj]
    SB = B[ii] + B[jj]
    K = key(SA, SB)
    order = np.argsort(K, kind='stable')
    Ks = K[order]
    hits = []
    # queries: F^6 - C^6
    CH = 4000
    for start in range(0, N, CH):
        f = np.arange(start, min(start + CH, N))
        QA = (A[f][:, None] - A[None, :]).ravel()
        QB = (B[f][:, None] - B[None, :]).ravel()
        QK = key(QA, QB)
        lo = np.searchsorted(Ks, QK, 'left')
        hi = np.searchsorted(Ks, QK, 'right')
        cand = np.nonzero(hi > lo)[0]
        for q in cand:
            for pos in range(lo[q], hi[q]):
                p = order[pos]
                if SA[p] == QA[q] and SB[p] == QB[q]:
                    fi = f[q // N]; ci = q % N
                    hits.append((ii[p], jj[p], ci, fi))
    out = []
    for (i, j, c, f) in hits:
        quad = [(us[i], vs[i]), (us[j], vs[j]), (us[c], vs[c]), (us[f], vs[f])]
        s6 = [(A[i], B[i]), (A[j], B[j]), (A[c], B[c]), (A[f], B[f])]
        assert s6[0][0]+s6[1][0]+s6[2][0] == s6[3][0]
        assert s6[0][1]+s6[1][1]+s6[2][1] == s6[3][1]
        if s6[3] in s6[:3]:
            continue                      # one term equals RHS -> other two sum to 0
        if all(q[1] == 0 for q in quad):
            continue                      # rational point
        out.append(quad)
    return out, N


if __name__ == '__main__':
    H = int(sys.argv[1]); Dmax = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    Ds = [d for d in range(2, Dmax + 1) if isqrt(d) ** 2 != d]
    print(f'A^6+B^6+C^6 = F^6 over Q(sqrt D): |u|,|v| <= {H}')
    tot = 0
    for D in Ds:
        out, N = search(D, H)
        tot += len(out)
        print(f'  D={D:<3} reps={N:<6} nontrivial={len(out)}' +
              ('  *** ' + str(out[:3]) if out else ''), flush=True)
    print('TOTAL nontrivial quadratic points found:', tot)
