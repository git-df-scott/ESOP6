#!/usr/bin/env python3
"""Fit the v3engine cumulative cost model from the bands recorded by
measure_v3.sh and print the projected wall time to reach given heights.

Usage: fit_v3.py meas.log [extra BAND lines on stdin]
"""
import sys, re, math

def parse(path):
    out = []
    for line in open(path):
        if not line.startswith("BAND"):
            continue
        d = dict(kv.split("=", 1) for kv in line.split() if "=" in kv)
        out.append(d)
    return out

def main():
    rows = [d for d in parse(sys.argv[1])]
    pts = []
    print("%-24s %10s %8s %14s %12s %10s %9s" %
          ("band", "width", "cands", "leaves", "search_s(4c)", "leaves/cs", "step1/leaf"))
    for d in rows:
        f0, f1 = int(d["fmin"]), int(d["fmax"])
        w = f1 - f0
        F = (f0 + f1) / 2.0
        leaves = int(d["leaves"]); s = float(d["search_s"]); th = int(d["threads"])
        if leaves == 0 or s == 0:
            continue
        lps = leaves / (s * th)
        pts.append((F, leaves / w, s * th / w, lps))
        print("%-24s %10d %8s %14d %12.1f %10.3g %9.4f" %
              ("(%g,%g]" % (f0, f1), w, d["processed"], leaves, s,
               lps, int(d["step1"]) / leaves))
    if len(pts) < 2:
        return
    # power-law fits  l(F) = c F^a   and   w(F) = d F^b   (core-seconds per unit f)
    def fit(xs, ys):
        n = len(xs)
        lx = [math.log(x) for x in xs]; ly = [math.log(y) for y in ys]
        mx = sum(lx) / n; my = sum(ly) / n
        b = sum((a - mx) * (c - my) for a, c in zip(lx, ly)) / sum((a - mx) ** 2 for a in lx)
        a = my - b * mx
        return math.exp(a), b
    F = [p[0] for p in pts]; L = [p[1] for p in pts]; W = [p[2] for p in pts]
    c, al = fit(F, L)
    d, be = fit(F, W)
    print("\nleaves per unit f   l(F) = %.4g * F^%.4f" % (c, al))
    print("core-s  per unit f  w(F) = %.4g * F^%.4f" % (d, be))
    print("residuals (measured/fit):",
          " ".join("%.2f" % (w / (d * f ** be)) for f, w in zip(F, W)))
    print("\ncumulative from f=2, k7 in {3,4}, 4 cores:")
    print("%10s %14s %16s %14s" % ("F", "leaves", "core-seconds", "wall (4 cores)"))
    for T in (4.4e6, 6e6, 7e6, 8e6, 1e7):
        Lc = c * T ** (al + 1) / (al + 1)
        Cc = d * T ** (be + 1) / (be + 1)
        h = Cc / 4 / 3600
        print("%10.3g %14.4g %16.4g %10.1f h" % (T, Lc, Cc, h))

main()
