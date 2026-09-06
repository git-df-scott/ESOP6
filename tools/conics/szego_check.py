# Random real e2=0 configurations of 5 points on the circle: check "all weights positive  <=>  |e1| < 1"
import numpy as np
rng = np.random.default_rng(5)
agree = 0; total = 0; pos_e1 = []; 
for trial in range(4000):
    th = rng.uniform(0, 2*np.pi, 4); om = np.exp(1j*th)
    e1 = om.sum(); e2 = sum(om[i]*om[j] for i in range(4) for j in range(i+1,4))
    if abs(e1) < 1e-9: continue
    o5 = -e2/e1
    if abs(abs(o5)-1) > 1e-9: continue   # need 5th on circle: generic 4 points won't work; instead solve properly below
for trial in range(20000):
    # choose 3 points and solve for 2 more on the circle with e2=0 via numeric root finding of a self-inversive polynomial
    th = rng.uniform(0, 2*np.pi, 3); om3 = np.exp(1j*th)
    # unknown pair (a,b) on circle: e2 = e2(3) + e1(3)(a+b) + ab = 0 -> choose a random on circle, b = -(e2_3 + e1_3 a)/(e1_3 + a): need |b|=1 -> 1 real condition; solve for arg a by bisection scanning
    e1_3 = om3.sum(); e2_3 = om3[0]*om3[1]+om3[0]*om3[2]+om3[1]*om3[2]
    phis = np.linspace(0, 2*np.pi, 2000, endpoint=False)
    a = np.exp(1j*phis); b = -(e2_3 + e1_3*a)/(e1_3 + a)
    f = np.abs(b) - 1
    for k in range(len(phis)-1):
        if f[k]*f[k+1] < 0:
            # refine by bisection
            lo, hi = phis[k], phis[k+1]
            for _ in range(50):
                mid = (lo+hi)/2; am = np.exp(1j*mid); bm = -(e2_3 + e1_3*am)/(e1_3 + am)
                if (np.abs(bm)-1)*f[k] > 0: lo = mid
                else: hi = mid
            am = np.exp(1j*lo); bm = -(e2_3 + e1_3*am)/(e1_3 + am)
            om = np.concatenate([om3, [am, bm]])
            if min(abs(om[i]-om[j]) for i in range(5) for j in range(i+1,5)) < 1e-6: continue
            M = np.array([[om[i]**kk for i in range(5)] for kk in (1,2,3)])
            Mr = np.vstack([M.real, M.imag])
            u,s,vt = np.linalg.svd(Mr); w = vt[-1]
            if s[-1] > 1e-8: continue
            allpos = np.all(w > 0) or np.all(w < 0)
            e1 = abs(om.sum()); total += 1
            agree += (allpos == (e1 < 1))
            if allpos: pos_e1.append(e1)
print(f"configurations tested={total}; 'all weights positive' agrees with '|e1|<1' in {agree} cases; positive cases |e1| range: {min(pos_e1) if pos_e1 else None:.3f}..{max(pos_e1) if pos_e1 else None:.3f}")
