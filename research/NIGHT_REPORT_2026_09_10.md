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
