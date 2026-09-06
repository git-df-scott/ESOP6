import numpy as np, sys
sys.argv=[sys.argv[0],"A","1","40"]
import conic_solver as cs
a,b,c,e = 0,1,2,3
rows = [{a:1},{b:1},{c:1},{a:1},{b:-1},{c:1},{c:1},{b:1},{a:1},{c:1},{b:-1},{a:1},{e:1},{},{e:1},{'c':1},{},{'c':1}]
A, v0 = cs.ansatz_matrix(rows, 4)
sols = cs.run(A, v0, 40, "A")
for s in sols[:12]:
    v = A@s["x"]+v0
    R,J = cs.residual_and_jac(v)
    Jx = J@A
    sv = np.linalg.svd(Jx, compute_uv=False)
    print(np.round(s["x"],5), "rank-svals:", np.round(sv,4), "res", s["res"])
