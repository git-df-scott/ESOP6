# Next attack: classify the rank-(1,1) genus-two quotients

No counterexample exists in the retained output. Do not repeat the exhausted
coefficient enumeration or treat its subgroup misses as global results.
Start with REPORT.md, INDEPENDENT_AUDIT.md section 10, and survivors.json.

## First fully specified target

```
c=67, representation (1,1,1,2),
H: w^2=t^6+67,
E1: y^2=x^3+67,
E2: y^2=x^3+4489.
```

Both elliptic ranks are certified exactly one. Retained independent
non-torsion points include

```
P1=(49/36,1801/216),
P2=(-63/4,193/8).
```

These are not asserted to be saturated Mordell–Weil generators. The two
maps from H are `(x,y)=(t^2,w)` and `(67/t^2,67w/t^3)`. The rational
points at infinity provide basepoints. Prime 13 is a concrete good ordinary
prime for these j=0 elliptic factors; required p-adic logarithm and regulator
conditions still need to be checked in the actual implementation.

The exact control `t=7/6,w=1801/216` lies on H and fails the cube condition.
A desired sextic lift requires nonzero t and w, with w a rational cube.
Writing `B/S=cuberoot(w)` and `A/S=t` reconstructs C_67; clearing denominators
then gives `(S,S,S,2S,A;B)` in integers, followed by both independent checkers.
Signs are removable because the final powers are even.

Use the rank-(1,1) bielliptic construction in Balakrishnan–Dogra,
[*Quadratic Chabauty and rational points I*, Theorem 1.4](https://arxiv.org/pdf/1601.00388).
The next work is to compute the finite bad-prime height sets and the p-adic
height/logarithm functions, isolate the resulting finite p-adic set, and
prove the rational reconstruction/Mordell–Weil sieve complete. Merely finding
several H-points is not sufficient. No such calculation was run in this pass.

## Prioritized extension

The 38 eligible fibres are listed with their exact rank-(1,1) quotient pairs
in survivors.json. The five available monic genus-two forms are

```
w^2=t^6+c,        w^2=t^6-c,
w^2=t^6+4c,       w^2=t^6-4c,
w^2=t^6+c^2/4.
```

Use the simplest rank-one pair and best arithmetic data per fibre. Retain
the correct map back to C_c: the first two need a cube condition on w;
the product quotients require the additional exact tests given in the audit.
All five are genus-two curves defined over Q; none is a rational parametrization.

For 35 survivors all nine elliptic rank intervals are certified. Their
genus-ten Jacobian rank upper bounds satisfy `r<=19` and the rational
Neron–Severi lower bound is 11, so the general level-two finiteness criterion
also applies. This higher-genus route is available but is less direct than
the explicit genus-two rank-(1,1) construction. If using the more explicit
general theorem, verify its extra hypotheses rather than relying on the
inequality alone.

Seven unresolved fibres have no available selected non-torsion subgroup in
the retained bounded point search; nineteen do not yet have a complete set
of nine certified rank intervals. Their open statuses are preserved. Failed
8/12-second tasks may be retried explicitly, without restarting successful
model work or any historical Y2 lottery.

The secondary lane's four locally soluble Pythagorean slices are separately
available. They are not rational points on X. Any quotient family through
the algebraic seeds must move at least two normalized tail ratios and leave
the Gaussian extra-quadratic trap. General degree-six moduli remains below
these now-specified arithmetic tasks in priority.
