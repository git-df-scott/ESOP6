# ESOP6 night report — 2026-09-10

Everything below is scaffolding and discovery for the counterexample hunt. Read this first;
tomorrow's plan is at the end. All scripts referenced live under `research/`.

## 1. Discoveries

### 1.1 A small algebraic point on the (6,1,4) slice

On the slice `x5 = 0` of the fourfold X: x1^6+x2^6+x3^6+x4^6+x5^6 = x6^6, there is a quadratic point

```
(9 + sqrt(-249))^6 + (9 - sqrt(-249))^6 + 14^6 + 18^6 = 22^6
```

Exactly: 2*Phi(9,-249) = 71,838,144 with Phi(A,Om) = A^6+15A^4 Om+15A^2 Om^2+Om^3, and
71,838,144 + 7,529,536 + 34,012,224 = 113,379,904 = 22^6. It is the unique primitive solution of
`2 Phi(A,Om) + r3^6 + r4^6 = r6^6` with Om not a square and r6 <= 80 (`cube_ansatz_2026_09_10/slice_quadratic_points.py`;
an N = 220 run is in `slice_quadratic_points_N220.txt`). This is the first concrete "seed": a conjugate pair
of coordinates plus three rational ones, which is exactly the shape a symmetric conic produces where an even
coordinate vanishes. LPS predicts no *rational* (6,1,4) point; this is the nearest algebraic thing.

### 1.2 The slice-contact principle (why every conic lane failed)

A rational curve C over Q on X meets every coordinate hyperplane. If the coordinate form x_j(s,t) has a
rational root, C passes through a rational point of the Calabi–Yau threefold slice X ∩ {x_j = 0}, i.e. a
(6,1,4) solution, which the Lander–Parkin–Selfridge conjecture forbids (the only rational points on the slices
are the trivial boundary points [0:0:0:0:1:±1]). Consequences:

* every coordinate form of a Q-curve on X must be irreducible over Q, or vanish only at a trivial boundary
  point (which for conics is impossible by the order-6 contact argument);
* for conics, each x_j contributes a *quadratic* point on a slice (six sparse events at once); the repo's
  exhaustive small searches are consistent with this;
* the first degree at which contact points are "free" is six (a Q-line meets a sextic threefold in a degree-6
  point): sextic rational curves with all six coordinate forms irreducible are the natural material;
* genus-0 cubic covers u^3 = L1 L2^2 M^3 and Family 0 (below) die immediately: x5 vanishes at a rational root.

### 1.3 Quotient-line formulation of symmetric conics

A conic stable under a coordinate swap carries an involution of P^1; on the quotient line with coordinates
(E1,E2) the identity becomes a binary sextic identity. With O^2 = Omega(E1,E2) the conic,

* Family I (x1<->x2):   (A+O)^6+(A-O)^6 + L3^6+L4^6+L5^6 = L6^6, A,L_j linear in E — a 2-dimensional family;
* Family II (two swaps): 2Phi(A,Om)+2Phi_F(D,Om)+L5^6 = L6^6 — 1-dimensional;
* Family 0 (x5 on a double cover): L1^6+..+L4^6 + Om^3 = L6^6 — 2-dimensional.

All three are nonempty over R (`symmetric_conics_2026_09_10/families.py`). Exact small sweeps found no rational
points (`catalecticant_sweep.py`: 263,520 nondegenerate cases, 0 hits). The reason is structural: Family I maps
to the seed threefold W = {2Phi(A,Om)+r3^6+r4^6=r6^6} with image a surface, so a generic seed (including 1.1) has
no conic through it even over C. Verified: the 6×6 system "Family I through the seed" has no solutions mod
17, 19, 23, 31 (`cube_ansatz_2026_09_10/seed_padic.py`), and double-precision Newton "solutions" were artefacts
(50-digit polish diverges). Lesson recorded: count seed dimension against family dimension before solving.

### 1.4 The cube-of-a-square reformulation and two Fano fourfolds

Astra's theorem: every Q-conic collapses mod 7 to one coordinate, because sixth powers of units are 1. The way
out is to stop writing x5^6 as a sixth power:

* Y : x1^6+x2^6+x3^6+x4^6 + T^3 = x6^6 in P(1,1,1,1,1,2), T = x5^2  (K = O(-1), Fano);
* Y2: x1^6+x2^6+x3^6+x4^6 + S^2 = x6^6 in P(1,1,1,1,1,3), S = x5^3  (K = O(-2), Fano).

X is the double cover of Y branched over T = 0 and the cyclic cubic cover of Y2. Rational curves of low degree
on Y / Y2 pull back to genus ≤ 1 curves on X when the branch divisor restricts with the right multiplicities.
This is precisely the Elkies mechanism (r^4+s^4+t^2 = 1 then t -> t^2) and the Bremner–Choudhry–Ulaş
polynomial method (P G1^6+Q G2^6+R G3^6+S F^3 = 0), run with five sixth powers and unit coefficients.

Families found and measured numerically (`cube_ansatz_2026_09_10/v3_numeric.py`):

| family | ansatz | curve on X | dim over C | real points |
|---|---|---|---|---|
| V3  | G_i quadratic, F quartic: ΣG_i^6+F^3=G6^6 | genus 1, w^2=F(t) | 3 | yes |
| V3' | G_i quadratic, S = K·M^3, K cubic: ΣG_i^6+K^2M^6=G6^6 | genus 1, u^3=K(t) | 3 | (solver built) |
| V3'' | S = L1 L2^2 M^3 | genus 0, sextic | 2 | dead by 1.2 |

Through a point of the Fano fourfold Y there are finitely many V3 curves (≥ 981 distinct found from 1500 starts,
`conics_through_point.py`); the same holds for V3' through a point of Y2. So rational points of Y, Y2 are
lottery tickets: each gives a zero-dimensional Galois set of genus-1 curves; a rational member is an elliptic
curve over Q on X, and any rational point on it with w ≠ 0 is a counterexample.

Seeds: Y(Q) is sparse at small height (none with x6 ≤ 45, `y_points.py`; heuristically c·N with c ~ 1e-4 because
of the 7-adic density). Y2(Q) is rich: 50 primitive points with x6 ≤ 40 (`y2_points_N40.json`), e.g.
6^6+15^6+12^6+18^6 + 9980^2 = 23^6.

### 1.5 Heuristic counts (why search was never going to work)

Primitive local densities (Astra): δ2 = 5/8, δ3 = 20/81, δ7 = 180/16807; product through 199 ≈ 0.00635.
Archimedean constant J = Γ(7/6)^5/Γ(5/6) ≈ 0.6088. Expected unordered primitive solutions with f ≤ N ≈
3.2e-5·log N: about 4e-4 below the EulerNet frontier 730,000, and the first "expected" solution at
log N ≈ 31,000. Only structure (a rational or positive-rank elliptic curve over Q) can produce a findable
counterexample. Every lane tonight is a structure lane.


### 1.6 Plane cubics on X (lines on Y2) — the cleanest material

A weighted line on Y2 (x1..x4, x6 linear in (p0,p1), S cubic) exists iff the binary sextic

```
L6^6 - L1^6 - L2^6 - L3^6 - L4^6  =  S(p0,p1)^2
```

is a perfect square. The plane Π spanned by the five linear forms then meets X in two plane cubics u^3 = ±S,
each a genus-1 curve over Q on X (u = x5). Slice contacts are degree-3 points (allowed). The family is
3-dimensional over C; through a rational point of Y2 there are finitely many such lines, and the system is
6 unknowns / 6 equations with one linear equation, so the exact pipeline (mod p → Hensel → rational
reconstruction → sympy check) decides it completely per seed. Observed mod-p counts are 0–4, i.e. the Galois sets
are tiny. `plane_cubics_through_seed.py` runs it over every Y2 seed; `plane_cubic_points.py` turns a rational
line into integer ESOP6 solutions by searching rational points on u^3 = S. Klein-symmetric sub-families are dead
(they force 2X^6+2Y^6 = Z^6). Overnight log: `plane_cubics_lottery_N70.log`.

Parallel lottery: `cubic_cover_lottery_N40.log` (genus-1 cubic covers u^3 = K(t), K cubic, through Y2 seeds,
numerical + rational recognition; 40–85 curves per seed found, none rational so far).

### 1.7 A Gaussian-integer solution of Euler's equation, and the quotient Z = X/τ

```
(7+11i)^6 + (7-11i)^6 + 8^6 + 12^6 + 15^6 = 17^6        (2Φ(7,-121) = 9,498,816; 17^6 = 24,137,569)
```

All six coordinates nonzero: a point of X over Q(i) of height 17 (`quadratic_points_X.py`; the only such point
with a conjugate pair and r6 ≤ 45, besides its multiples). Its secant with the conjugate point has no rational
residual. It is a rational point of the Fano fourfold Z = X/τ (τ: x1<->x2; coordinates e1 = x1+x2, e2 = x1x2,
K_Z = O(-1)), at (e1,e2,x3,x4,x5,x6) = (14,170,8,12,15,17); X → Z is the double cover branched along x1 = x2,
and X(Q) = {z ∈ Z(Q): e1^2 - 4e2 is a square}. Family I = weighted lines on Z; weighted conics on Z (e2 quartic)
pull back to genus-1 curves w^2 = e1^2 - 4e2 on X, a 3-dimensional family, finitely many through each Z-point
(`z_conics_through_point.py`: 14 and 16 found through the two known points, none rational).

### 1.8 Periodic table of genus ≤ 1 families (via u^6 = f on P^4)

Take x1..x4, x6 of degree d in t, f = x6^6 - Σ x_i^6 (degree 6d), and the μ6-cover u^6 = f (u = x5). Root
multiplicities of f decide the genus; the slice principle decides which are admissible over Q.

| d | pattern of f | cover | genus | status |
|---|---|---|---|---|
| 1 | f = S^2, S cubic | u^3 = S | 1 | plane cubics on X; exact per-seed solver; no rational line with two points ≤ 150 |
| 1 | f = T^3, T quadratic | u^2 = T | 0 | Astra's design conics, dead 7-adically |
| 2 | f = T^3, T irreducible quartic | w^2 = T | 1 | V3, 3-dim; needs Y(Q) seeds (none ≤ 150) |
| 2 | f = K^2 M^6, K irreducible cubic | u^3 = K | 1 | V3', 3-dim; numerical lottery running |
| 2 | f = (L1 L2^2 M^3)^2 | genus 0 sextic | 0 | dead: x5 = 0 at a rational root |
| 3 | f = (Q^2 R)^3, Q,R quadratics | w^2 = R | 0 | sextic rational curves, 2-dim; two quadratic slice points each |

Collinearity test (`y2_collinear_pairs.py`): no two of the 298 Y2 points with x6 ≤ 150 lie on a common plane
cubic of X (4.2M pair/permutation/sign tests, exact).

### 1.9 The odd-factor principle (a Kepler moment, late in the night)

Along any real curve on X write the identity as a product: x6^6 - x5^6 = (x6^k - x5^k)(...) with the relevant real
factor N. If N has odd degree in the parameter it has a real root t1, and there Σ_{i≤4} x_i^6 = 0 forces
x1 = x2 = x3 = x4 = 0: the curve passes through the boundary point B = [0:0:0:0:1:1]. For a Q-curve t1 is then
rational (a quadratic t1 would make all four coordinate forms proportional). Consequences, each checked numerically:

* Plane cubics u^3 = ±S (lines on Y2): N = L6^3 - S is a real cubic; linear forms cannot all vanish at t1 unless
  proportional. **No real plane cubic of this type exists.** The 68+31 seeds tested (all mod-p counts, no rational
  lines) were complex-only Galois sets. Lane closed; the solver remains as a template.
* Cubic covers V3' (u^3 = K M^3): N = G6^3 - K M^3 has degree 9, so real members pass through B at a rational
  parameter, G_i = t·(linear). That boundary-contact system is zero-dimensional and has **no solutions even over C**
  (`boundary_cubic_cover.py`, 800 complex and 1500 real starts). V3' has no real members; the per-seed counts
  "near-real 0" were this theorem showing up as data. Lane closed.
* Quadratic covers V3 (w^2 = F): N = G6^2 - F is a quartic, no forced root; real members exist. Live, needs Y(Q)
  seeds (structured search to x6 ≤ 360 running: `y_points_structured_N360.log`).
* Z-quotient conics (τ-symmetric genus-1 curves, w^2 = e1^2 - 4e2 quartic): 3-dimensional with real members
  (local dim 7); live; tickets are Z(Q) points = quadratic points of X with a conjugate pair. Known: (7±11i,8,12,15,17)
  and (9±√-249,14,18,0,22), neither carrying a rational curve of this family (14 and 16 complex curves each).
* Sextic rational curves: N = x6 - x5 has even degree; real nondegenerate sextics exist (two with generic local
  dimension 5). Live; no seed mechanism (through a point is overdetermined).

## 2. Lanes closed tonight

* Family 0 and genus-0 cubic covers: rational root of x5 (1.2).
* Family I through any seed of W not on the special image surface; the seed in 1.1 in particular.
* Klein-4 symmetric conics: the E2^6 coefficient forces 4ν^6+8μ^6 to be a cube, i.e. rational points on the
  genus-4 curve y^3 = x^6+32 beyond the trivial ones.
* "Design" conics (four coordinates linear in two variables): Astra's 7-adic theorem (unchanged).
* The repeated-coordinate surface 2X^6+2Y^6+Z^6=W^6 as a counterexample source: general type, expected count N^-2.

## 3. Scaffolding delivered

* `tools/fourfold_curve_search.py`: complex/real Gauss–Newton search for degree-d rational curves on the full
  fourfold with degeneracy classification and local-dimension readout.
* `symmetric_conics_2026_09_10/families.py`: quotient-line formulations of Families I, II.
* `cube_ansatz_2026_09_10/seed_padic.py`: the exact adapter — mod-p enumeration → Hensel lift → rational
  reconstruction → exact verification — for zero-dimensional systems.
* `cube_ansatz_2026_09_10/conics_through_point.py`, `cubic_cover_through_point.py`: zero-dimensional solvers for
  genus-1 curves through a rational seed on Y / Y2 with rational recognition and exact check.
* `cube_ansatz_2026_09_10/y_points.py`, `y2_points.py`: seed enumerators.
* Overnight job: `cubic_cover_lottery_N40.log` (all Y2 seeds with x6 ≤ 40, 3000 starts each).

## 4. Tomorrow

0. The plane-cubic and cubic-cover lotteries are closed by the odd-factor principle (1.9); do not rerun them.
1. Read `cubic_cover_lottery_N40.log` and `cubic_cover_rational_hits.json`. Any hit: find rational points on the
   elliptic curve u^3 = K(t) (search t, then 2-descent if needed), specialize, run both verifiers.
2. Extend Y2 seeds to x6 ≤ 70 (`y2_points_N70.log`) and rerun; raise starts until the per-seed count saturates.
3. V3 on Y: Y needs a Meyrignac-scale search (targets x6^6 − T^3 with two of x1..x4 divisible by 21). Port the
   repo's C pair-sum engine to targets of that shape.
4. Sextic rational curves: run `fourfold_curve_search.py --degree 6` (complex, then real), classify
   components, and test irreducibility of coordinate forms. This is the "natural material" of 1.2.
5. Integrate the 7-adic verdicts on Families I/II when `symmetric_conics_2026_09_10/LOCAL_OBSTRUCTIONS.md` lands.
6. Literature to obtain: Bremner 1981 (PLMS 43), Kuwata 2007 (RMJM 37), Letac 1942 — the only published
   elliptic-curve attacks on sixth powers.
