# Moduli traversal lane (Claude, 2026-09-11)

Complementary to Astra's twisted-cubic fibration lane. Target: rational points on the moduli of rational curves
of degree d on X: x1^6+..+x5^6 = x6^6, since such a point IS a rational curve on X and hence infinitely many
counterexamples (any real point with nonzero coordinates; signs are irrelevant).

Normalization (rational-compatible, so rational moduli points are rational solutions):
  x6 = t^d + 0·t^{d-1} + ... with a_{d-1} = 0 (translation), a_1 = 0 (other shear), a_{d-2} = a_{d-3} (scaling),
  leading coefficient 1 (Q^* scaling). Remaining unknowns 6(d+1) - 4, equations 6d+1: one-dimensional.

Method: pseudo-arclength continuation along real components; at every step a common-denominator rational test over
all coordinates (false-positive probability ~ (2qh)^{6d-3}); candidates go through exact sympy identity check and
tools/verify_esop6.py. Negative results are recorded as: component id, arc length covered, box, step h, denominator
bound D => no rational curve on that component with denominator ≤ D inside the box.

Order: conics (d=2, validation; expected none), then sextics (d=6, generic components with local dimension 5).
Odd d impossible over R. d=4 also considered (quartics: 30 unknowns).

Note on the fibration lane: X/Z and Y/Z square conditions on E_c give a (Z/2)^2 cover (genus 4 then 10) with no
isogeny-type descent, and per-fibre expected point counts are dominated by small heights (~c^{-2/3}); rank-0 fibres
are certified dead, survivors are as hard as the original problem on that fibre.

## Results (2026-09-11)

* Generic moduli: no rational-compatible normal form exists (the two unipotent gauges do not commute; invariant
  coordinates have degree ≥ 12 in the coefficients), so traversal-with-detection is restricted to symmetric
  sub-moduli with rational anchors.
* Z/3-symmetric sextics (t ↦ 1/(1−t)): real nondegenerate members exist (`z3_sextics.py`, local dim 3 = 1 + torus +
  scale), but the residual gauge is the anisotropic torus Q(ω)^*/Q^* acting with degree 6, again with no rational
  normal form.
* S3-symmetric sextics (`s3_sextics.py`): the only character combination with real members is (x1 ι-even; x4, x6
  S3-invariant; x5 sign character), a 2-dimensional family with scaling as the only gauge. Its exact search
  (`s3_fast.py`, structured p^5 enumeration + Hensel + reconstruction) finds nothing, and the reason is structural:
  the sign-character sextic is st(s−t)·(cubic), so x5 vanishes at the rational S3-orbit {0,1,∞}; by the slice-contact
  principle every Q-curve in the family carries rational (6,1,4) points. Lane closed.
* Height-300 cubic-cover search complete (1239 seeds, no rational curve); height-500 running.
