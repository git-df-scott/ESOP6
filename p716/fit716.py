#!/usr/bin/env python3
"""Fit leaves = kappa*F^4 to BAND lines and extrapolate wall time on 4 cores."""
import sys

bands = []
for ln in open(sys.argv[1]):
    if ln.startswith("BAND"):
        p = ln.split()
        bands.append((int(p[2]), int(p[3]), int(p[6]), int(p[7]), float(p[9])))
# dedupe (last line repeats)
seen = {}
for b in bands:
    seen[b[0]] = b
bands = [seen[k] for k in sorted(seen)]
F, L, Q, P, T = bands[-1]
print("final: F=%d leaves=%.4g queries=%.4g positives=%.4g elapsed=%.1f s" % (F, L, Q, P, T))
kappa = L / F ** 4
print("kappa = leaves/F^4 = %.4g   (mask pass %.4f, bloom fp %.3g)" % (kappa, Q / L, P / Q))
# search-only time (subtract table build = elapsed of first band minus its own work)
t0 = bands[0][4] - 0.0
build = min(b[4] for b in bands) - 0.0
# steady-state leaf rate from the last half of the run
i = len(bands) // 2
dL = bands[-1][1] - bands[i][1]
dT = bands[-1][4] - bands[i][4]
rate = dL / dT
print("startup (pair bloom + table build) ~ %.1f s" % bands[0][4])
print("steady-state leaf rate  = %.4g leaves/s on 4 cores = %.4g leaves/s/core"
      % (rate, rate / 4))
print()
print("extrapolation, leaves = kappa*F^4, at the measured steady-state rate:")
print("%8s %12s %12s %10s" % ("F", "leaves", "core-hours", "wall h (4c)"))
for f in (1500, 2500, 3000, 3500, 4000, 5000):
    lv = kappa * f ** 4
    print("%8d %12.4g %12.2f %10.2f" % (f, lv, lv / (rate / 4) / 3600, lv / rate / 3600))
