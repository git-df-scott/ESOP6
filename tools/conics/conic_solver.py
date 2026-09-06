"""
Numerical search for conics on the Fermat sextic fourfold
    X: p1^6 + p2^6 + p3^6 + p4^6 + p5^6 = p6^6,
where each p_i(t) = a_i + b_i t + c_i t^2 (dehomogenised, s=1).
The identity gives 13 polynomial equations in the 18 coefficients.
An "ansatz" is an affine map x -> v = A x + v0 from free unknowns to
the 18 coefficients, encoding symmetry restrictions and normalisations.
Gauss-Newton from random complex starts; solutions are collected,
deduplicated, and classified (real? definite? candidate rational?).
"""
import numpy as np, sys, json, itertools, math

rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 0)

def polymul(a, b):
    return np.convolve(a, b)

def poly_pow(p, n):
    r = np.array([1.0+0j])
    for _ in range(n):
        r = polymul(r, p)
    return r

def residual_and_jac(v):
    """v: 18 complex coefficients (a1,b1,c1,...,a6,b6,c6). Returns R (13,) and J (13,18)."""
    R = np.zeros(13, dtype=complex)
    J = np.zeros((13, 18), dtype=complex)
    for i in range(6):
        a, b, c = v[3*i], v[3*i+1], v[3*i+2]
        p = np.array([a, b, c])
        p5 = poly_pow(p, 5)
        p6 = polymul(p5, p)
        sgn = 1 if i < 5 else -1
        R += sgn * p6
        # d(p^6)/da = 6 p^5, /db = 6 p^5 t, /dc = 6 p^5 t^2
        d = 6 * p5
        J[0:11, 3*i]   += sgn * d
        J[1:12, 3*i+1] += sgn * d
        J[2:13, 3*i+2] += sgn * d
    return R, J

def solve(A, v0, x0, iters=100, tol=1e-13):
    x = x0.copy()
    for it in range(iters):
        v = A @ x + v0
        R, J = residual_and_jac(v)
        nr = np.linalg.norm(R)
        if nr < tol:
            return x, nr, it
        Jx = J @ A
        dx, *_ = np.linalg.lstsq(Jx, -R, rcond=None)
        # damping
        lam = 1.0
        while lam > 1e-4:
            xn = x + lam*dx
            Rn, _ = residual_and_jac(A @ xn + v0)
            if np.linalg.norm(Rn) < nr:
                x = xn
                break
            lam *= 0.5
        else:
            x = x + 1e-3*dx
        if np.linalg.norm(x) > 1e6:
            return None, np.inf, it
    v = A @ x + v0
    R, _ = residual_and_jac(v)
    return x, np.linalg.norm(R), iters

def classify(v):
    """v: 18 coefficients. Returns dict with realness, definiteness."""
    real = np.all(np.abs(v.imag) < 1e-8 * max(1, np.abs(v).max()))
    out = {"real": bool(real)}
    if real:
        vr = v.real
        discs = []
        for i in range(6):
            a, b, c = vr[3*i], vr[3*i+1], vr[3*i+2]
            discs.append(b*b - 4*a*c)
        out["disc"] = discs
        out["all_definite"] = all(d < -1e-9 for d in discs)
        M = vr.reshape(6,3); sv = np.linalg.svd(M, compute_uv=False)
        out["span_rank"] = int(np.sum(sv > 1e-6*sv[0]))
        out["genuine"] = bool(out["all_definite"] and out["span_rank"] == 3)
        out["degenerate_zero"] = any(abs(vr[3*i])+abs(vr[3*i+1])+abs(vr[3*i+2]) < 1e-8 for i in range(6))
    return out

def run(A, v0, nstarts, name, real_starts=False, scale=2.0):
    n = A.shape[1]
    sols = []
    for k in range(nstarts):
        if real_starts:
            x0 = rng.normal(size=n)*scale + 0j
        else:
            x0 = (rng.normal(size=n) + 1j*rng.normal(size=n))*scale
        x, nr, it = solve(A, v0, x0)
        if x is None or nr > 1e-9:
            continue
        v = A @ x + v0
        # dedup on coefficient vector (rounded), also modulo sign flips of p1..p5 and permutations of p1..p5 handled loosely by sorting keys
        key = tuple(np.round(np.abs(v), 6))
        if any(np.linalg.norm(np.abs(v)-np.abs(s["v"])) < 1e-5 for s in sols):
            continue
        cl = classify(v)
        sols.append({"v": v, "x": x, "cls": cl, "res": nr})
    return sols

def report(sols, name, xnames):
    print(f"=== {name}: {len(sols)} distinct solutions ===")
    nreal = sum(1 for s in sols if s["cls"]["real"])
    ndef = sum(1 for s in sols if s["cls"].get("all_definite"))
    ngen = sum(1 for s in sols if s["cls"].get("genuine"))
    ndeg = sum(1 for s in sols if s["cls"].get("degenerate_zero"))
    print(f"real: {nreal}   real & all-definite: {ndef}   GENUINE (definite, spans a plane): {ngen}   real with a vanishing p_i: {ndeg}")
    for s in sols:
        if s["cls"]["real"] and s["cls"]["genuine"]:
            xs = ", ".join(f"{n}={z.real:.10f}" for n, z in zip(xnames, s["x"]))
            print(f"  REAL definite={s['cls']['all_definite']} degen={s['cls']['degenerate_zero']} disc={np.round(s['cls']['disc'],6).tolist()}")
            print(f"       {xs}")
    return sols

def ansatz_matrix(rows, n):
    """rows: list of 18 dicts {var_index: coeff} plus constant key 'c'."""
    A = np.zeros((18, n), dtype=complex); v0 = np.zeros(18, dtype=complex)
    for r, spec in enumerate(rows):
        for k, c in spec.items():
            if k == 'c': v0[r] = c
            else: A[r, k] = c
    return A, v0


def rank_at(A, v0, x):
    v = A @ x + v0
    R, J = residual_and_jac(v)
    sv = np.linalg.svd(J @ A, compute_uv=False)
    return int(np.sum(sv > 1e-8 * sv[0])), sv

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "A"
    N = int(sys.argv[3]) if len(sys.argv) > 3 else 400
    if which == "A":
        # p1=1+bt+ct^2, p2=1-bt+ct^2, p3=c+bt+t^2, p4=c-bt+t^2, p5=e(1+t^2), p6=g(1+t^2)
        b,c,e,g = 0,1,2,3
        rows = [{'c':1},{b:1},{c:1},  {'c':1},{b:-1},{c:1},  {c:1},{b:1},{'c':1},  {c:1},{b:-1},{'c':1},
                {e:1},{},{e:1},  {g:1},{},{g:1}]
        A, v0 = ansatz_matrix(rows, 4); names = ["b","c","e","g"]
        title = "A: (12)(34) t->-t and (13)(24) s<->t; p1 normalised"
    elif which == "B":
        # p1=1+b1t+c1t^2, p2=1-b1t+c1t^2, p3=a3+b3t+c3t^2, p4=a3-b3t+c3t^2, p5=a5+c5t^2, p6=g(1+t^2)
        b1,c1,a3,b3,c3,a5,c5,g = range(8)
        rows = [{'c':1},{b1:1},{c1:1},  {'c':1},{b1:-1},{c1:1},  {a3:1},{b3:1},{c3:1},  {a3:1},{b3:-1},{c3:1},
                {a5:1},{},{c5:1},  {g:1},{},{g:1}]
        A, v0 = ansatz_matrix(rows, 8); names = ["b1","c1","a3","b3","c3","a5","c5","g"]
        title = "B: (12)(34) t->-t, p5 p6 even (expected 1-dim)"
    elif which == "D":
        # full: a1=1; p6=g(1+t^2); b5=0 fixes rotation
        rows = [{'c':1}]; k = 0
        for i in range(5):
            for j in range(3):
                if i == 0 and j == 0: continue
                if i == 4 and j == 1: rows.append({}); continue
                rows.append({k:1}); k += 1
        g = k; rows += [{g:1},{},{g:1}]
        A, v0 = ansatz_matrix(rows, k+1)
        names = [f"{n}{i+1}" for i in range(5) for n in "abc" if not (i==4 and n=='b') and not (i==0 and n=='a')] + ["g"]
        title = "D: full system, a1=1, p6=g(1+t^2), b5=0 (expected 1-dim)"
    elif which == "S3":
        # surface slice 2X^6+2Y^6+Z^6=W^6 : p1=p2, p3=p4 ; a1=1, p6=g(1+t^2), rotation fixed by b5=0
        b1,c1,a3,b3,c3,a5,c5,g = range(8)
        rows = [{'c':1},{b1:1},{c1:1}, {'c':1},{b1:1},{c1:1}, {a3:1},{b3:1},{c3:1}, {a3:1},{b3:1},{c3:1}, {a5:1},{},{c5:1}, {g:1},{},{g:1}]
        A, v0 = ansatz_matrix(rows, 8); names=["b1","c1","a3","b3","c3","a5","c5","g"]; title="S3 slice 2X^6+2Y^6+Z^6=W^6 conics (8 unknowns, 13 eqs)"
    elif which == "S2":
        # surface slice 3X^6+Y^6+Z^6=W^6 : p1=p2=p3
        b1,c1,a4,b4,c4,a5,c5,g = range(8)
        rows = [{'c':1},{b1:1},{c1:1}]*3 + [{a4:1},{b4:1},{c4:1}, {a5:1},{},{c5:1}, {g:1},{},{g:1}]
        A, v0 = ansatz_matrix(rows, 8); names=["b1","c1","a4","b4","c4","a5","c5","g"]; title="S2 slice 3X^6+Y^6+Z^6=W^6 conics"
    elif which == "T4":
        # threefold slice 2X^6+Y^6+Z^6+U^6=W^6 : p1=p2 (11 unknowns, 13 eqs)
        rows = [{'c':1},{0:1},{1:1}, {'c':1},{0:1},{1:1}]; k=2
        for i in range(3):
            for j in range(3):
                if i==2 and j==1: rows.append({}); continue
                rows.append({k:1}); k+=1
        rows += [{k:1},{},{k:1}]
        A, v0 = ansatz_matrix(rows, k+1); names=[f"x{i}" for i in range(k+1)]; title="T4 slice 2X^6+Y^6+Z^6+U^6=W^6 conics (11 unknowns, 13 eqs)"
    REAL = "--real" in sys.argv

    sols = run(A, v0, N, which, real_starts=REAL, scale=(1.5 if REAL else 2.0))
    report(sols, title, names)
    ranks = {}
    for s in sols:
        r, sv = rank_at(A, v0, s["x"]); ranks[r] = ranks.get(r, 0) + 1
    print("Jacobian rank histogram over solutions:", ranks, " (unknowns:", A.shape[1], ")")
    # save
    np.save(f"sols_{which}{'_real' if REAL else ''}.npy", np.array([s["v"] for s in sols]))
