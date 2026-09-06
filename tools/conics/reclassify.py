import numpy as np, sys, glob
import conic_solver as cs
for f in sorted(glob.glob("sols_*.npy")):
    V = np.load(f)
    if V.size == 0: print(f, "empty"); continue
    cl = [cs.classify(v) for v in V]
    nreal = sum(c["real"] for c in cl); ngen = sum(c.get("genuine",False) for c in cl)
    ranks = {}
    for c in cl:
        if c["real"]: ranks[c["span_rank"]] = ranks.get(c["span_rank"],0)+1
    print(f"{f}: n={len(V)} real={nreal} genuine_definite_conics={ngen} real span-rank histogram={ranks}")
    for v,c in zip(V,cl):
        if c.get("genuine"): print("   GENUINE:", np.round(v.real,8).tolist())
