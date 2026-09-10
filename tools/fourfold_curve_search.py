#!/usr/bin/env python3
"""Numerical search for rational curves P^1 -> {F1^6+...+F5^6 = F6^6}.

Unknowns: six polynomials F_i of degree <= d (complex or real coefficients).
Equations: the 6d+1 coefficients of  sum_{i<=5} F_i^6 - F_6^6  vanish.
Gauss-Newton with the pseudo-inverse from random starts. A converged point is
classified as degenerate when some proper sub-sum of the signed sixth powers
vanishes identically (e.g. F_e^6 == F_f^6), which never yields an ESOP6 point.
"""
import numpy as np, sys, argparse, itertools, json

SIGNS = np.array([1, 1, 1, 1, 1, -1], dtype=float)

def powers(F, k):
    P = np.array([1.0 + 0j]) if np.iscomplexobj(F) else np.array([1.0])
    for _ in range(k):
        P = np.convolve(P, F)
    return P

def residual_and_jac(x, d, real):
    n = d + 1
    F = x.reshape(6, n)
    deg = 6 * d + 1
    r = np.zeros(deg, dtype=x.dtype)
    J = np.zeros((deg, 6 * n), dtype=x.dtype)
    for i in range(6):
        F5 = powers(F[i], 5)
        F6 = np.convolve(F5, F[i])
        r += SIGNS[i] * F6
        for j in range(n):
            # d/d c_{ij} of F_i^6 = 6 F_i^5 t^j
            J[j:j + len(F5), i * n + j] += 6 * SIGNS[i] * F5
    return r, J

def normalize(x, d):
    """Use the two scaling symmetries (overall scale, t -> lam*t) to balance x."""
    n = d + 1
    F = x.reshape(6, n).copy()
    a0 = np.linalg.norm(F[:, 0]); ad = np.linalg.norm(F[:, d])
    if a0 > 0 and ad > 0:
        lam = (a0 / ad) ** (1.0 / d)
        F = F * (lam ** np.arange(n))
    F /= np.abs(F).max()
    return F.ravel()

def rel_residual(x, d, real):
    r, J = residual_and_jac(x, d, real)
    n = d + 1
    F = x.reshape(6, n)
    scale = sum(np.linalg.norm(powers(F[i], 6)) for i in range(6))
    return np.linalg.norm(r) / scale, r, J

def newton(x, d, real, iters=60, tol=1e-13):
    for it in range(iters):
        x = normalize(x, d)
        nr, r, J = rel_residual(x, d, real)
        if nr < tol:
            return x, nr, it
        dx = np.linalg.lstsq(J, -r, rcond=None)[0]
        x = x + dx
    x = normalize(x, d)
    nr, r, J = rel_residual(x, d, real)
    return x, nr, iters

def classify(x, d, eps=1e-8):
    n = d + 1
    F = x.reshape(6, n)
    norms = np.array([np.linalg.norm(f) for f in F])
    if (norms < eps * norms.max()).any():
        return "zero-coordinate"
    P6 = [SIGNS[i] * powers(F[i], 6) for i in range(6)]
    scale = max(np.linalg.norm(p) for p in P6)
    for k in range(1, 6):
        for S in itertools.combinations(range(6), k):
            s = sum(P6[i] for i in S)
            if np.linalg.norm(s) < 1e-7 * scale:
                return "vanishing-subsum:" + "".join("abcdef"[i] for i in S)
    return "NONDEGENERATE"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--degree", type=int, default=2)
    ap.add_argument("--starts", type=int, default=200)
    ap.add_argument("--real", action="store_true")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    n = a.degree + 1
    tally = {}
    found = []
    for s in range(a.starts):
        if a.real:
            x0 = rng.standard_normal(6 * n)
        else:
            x0 = rng.standard_normal(6 * n) + 1j * rng.standard_normal(6 * n)
        x, nr, it = newton(x0, a.degree, a.real)
        if nr < 1e-11:
            c = classify(x, a.degree)
            tally[c] = tally.get(c, 0) + 1
            if c == "NONDEGENERATE":
                _, _, J = rel_residual(x, a.degree, a.real)
                sv = np.linalg.svd(J, compute_uv=False)
                rank = int((sv > 1e-9 * sv[0]).sum())
                tally["localdim=%d" % (6 * n - rank)] = tally.get("localdim=%d" % (6 * n - rank), 0) + 1
                found.append(x)
        else:
            tally["no-convergence"] = tally.get("no-convergence", 0) + 1
    print(json.dumps(tally, indent=1))
    if found:
        print("NONDEGENERATE curves found:", len(found))
        for x in found[:5]:
            print(np.round(x.reshape(6, n), 6))
        if a.out:
            np.save(a.out, np.array(found))

if __name__ == "__main__":
    main()
