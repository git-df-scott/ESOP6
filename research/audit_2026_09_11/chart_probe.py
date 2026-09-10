"""Chart-coverage and Newton-completeness probe for cubic_cover_through_point.py.

The old solver fixes the gauge by g11 = 1, g61 = 0, reachable only if the GAUGE
INVARIANT W16 = x1 g61 - x6 g11 is nonzero.  Since the equation is symmetric in
x1..x4, running the SAME solver on the 4 cyclic rotations of (x1,x2,x3,x4) gives
the four charts (j,6), j = 1..4; a curve is invisible to all four only if
W_j6 = 0 for every j, which forces g_1 parallel to x (a doubly traversed line).

Each found solution is canonicalised into a gauge-invariant fingerprint
(c chosen so that sum x_i g_i1 = 0, lambda so that max|g_i1| = 1), so solutions
from different charts can be identified.  Output: per-chart counts, union count,
and how many union members each chart missed (= combined chart hole + Newton
miss rate).
"""
import sys, json, itertools
import numpy as np
sys.path.insert(0, '/home/user/ESOP6/research/cube_ansatz_2026_09_10')
import cubic_cover_through_point as cc

def canon(s, v):
    """v: solution of the old chart.  Return hashable fingerprint, gauge-normalised."""
    x = [s.x[0], s.x[1], s.x[2], s.x[3], s.x6]
    g1 = np.array([1.0+0j, v[1], v[3], v[5], 0.0])
    g2 = np.array([v[0], v[2], v[4], v[6], v[7]])
    K = np.array([s.s0, v[8], v[9], v[10]])
    M = np.array([1.0+0j, v[11]])
    X = np.array(x)
    sx2 = (X*X).sum()
    lam = 1.0
    c = -(X*g1).sum()/(2*sx2)
    # apply (t0 -> t0 + c t1, t1 -> lam t1) with lam = 1 first
    G1n = 2*c*X + g1
    G2n = c*c*X + c*g1 + g2
    Kn = np.array([K[0], 3*c*K[0]+K[1], 3*c*c*K[0]+2*c*K[1]+K[2], c**3*K[0]+c*c*K[1]+c*K[2]+K[3]])
    Mn = np.array([M[0], c*M[0]+M[1]])
    # now scale lam so that max|G1n| = 1
    mx = np.abs(G1n).max()
    if mx < 1e-12: return None
    lam = 1.0/mx
    # pick a deterministic phase: make the largest-modulus entry positive real
    i0 = int(np.argmax(np.abs(G1n)))
    ph = G1n[i0]/abs(G1n[i0])
    lam = lam/ph
    G1n = lam*G1n
    G2n = lam*lam*G2n
    Kn = np.array([Kn[0], lam*Kn[1], lam**2*Kn[2], lam**3*Kn[3]])
    Mn = np.array([Mn[0], lam*Mn[1]])
    # coordinates 1..4 are interchangeable: sort their (x, g1, g2) triples
    trip = sorted([(np.round(X[i].real, 6), np.round(X[i].imag, 6), np.round(G1n[i].real, 5),
                    np.round(G1n[i].imag, 5), np.round(G2n[i].real, 5), np.round(G2n[i].imag, 5))
                   for i in range(4)])
    tail = tuple(np.round([G1n[4].real, G1n[4].imag, G2n[4].real, G2n[4].imag,
                           Kn[1].real, Kn[1].imag, Kn[2].real, Kn[2].imag,
                           Kn[3].real, Kn[3].imag, Mn[1].real, Mn[1].imag], 4))
    return tuple(itertools.chain(*trip)) + tail

def match(a, b, tol=1e-3):
    return max(abs(np.array(a)-np.array(b))) < tol

def probe(seed, starts, rotations=4, real=True, rngseed=0):
    x = [int(v) for v in seed['x']]; x6 = int(seed['x6']); S = int(seed['S'])
    allc = []      # list of (fingerprint, set of charts)
    percharts = []
    for r in range(rotations):
        xr = x[r:]+x[:r]
        s = cc.Sys(xr, x6, S)
        sols = s.solve_all(starts, seed=rngseed+1000*r, real=real)
        fps = [canon(s, v) for v in sols]
        fps = [f for f in fps if f is not None]
        percharts.append(len(fps))
        for f in fps:
            for e in allc:
                if match(e[0], f): e[1].add(r); break
            else:
                allc.append([f, {r}])
    return percharts, allc

if __name__ == "__main__":
    pts = json.load(open(sys.argv[1]))
    starts = int(sys.argv[2]); nseeds = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    skip = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    for p in pts[skip:skip+nseeds]:
        pc, allc = probe(p, starts)
        n = len(allc)
        only = [sum(1 for e in allc if e[1] == {r}) for r in range(4)]
        print("seed", p['x'], p['x6'], "per-chart", pc, "union", n,
              "found-by-all-4", sum(1 for e in allc if len(e[1]) == 4),
              "found-by-exactly-one", only, flush=True)
