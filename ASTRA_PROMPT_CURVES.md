# Astra prompt — rational and elliptic curves on the sextic fourfold (compute split)

Copy everything below the line into Astra's session.

---

You are working on ESOP6, the search for a counterexample to Euler's sum of
powers conjecture at k = 6:

    a^6 + b^6 + c^6 + d^6 + e^6 = f^6,   all positive integers.

Repository: `git-df-scott/ESOP6`. Read `README.md`, `STATUS.md`,
`ASTRA_HANDOFF.md`, and `ASTRA_PROMPT_CURVES.md` (this file) first.

## What is settled and must not be repeated

1. Brute force is dead. The concentrated-case sweep stands at f = 4,310,000
   and each further 10k band costs ~18 CPU-minutes on 4 cores, growing like
   F^3 per band. Do not extend it.
2. The heuristic count of "generic" solutions below F is c * log F with
   c on the order of 1e-5 after the forced 42-divisibility of four terms.
   A counterexample, if it exists, is structured: it lies on a rational or
   elliptic curve of the fourfold. Construction is the only lane.
3. Every algebraic attempt so far was on 2-dimensional slices such as
   2X^6 + 2Y^6 + Z^6 = W^6. Those are surfaces of general type where rational
   curves have negative expected dimension. Do not work on slices.

## The target variety and the key dimension count

X : x1^6 + x2^6 + x3^6 + x4^6 + x5^6 = x6^6 in P^5 is a smooth Calabi-Yau
fourfold. On a CY fourfold the expected dimension of the family of rational
curves is +1 in every degree, and 0 for genus-one curves. Conics on a general
sextic fourfold form a 1-parameter family (14 parameters, 13 conditions).
The Fermat fourfold has a large automorphism group, so it may carry more.

Two facts to use throughout:

- Positivity forces even degree. On a real rational curve p(t) in the
  positive chamber, p6(t) has no real zero (a real zero would force all
  six coordinates to vanish there, making the parametrisation reducible).
  So every coordinate is a definite binary form and the degree is even.
- The sign group (Z/2)^5 and S5 act on the curve families. Fixed loci of
  subgroups give small polynomial systems.

## Split of work

The other session has finished the CONIC reduction; read `CONIC_LANE.md`.
Summary: conics on X are identities `sum_i (z_i + x_i cos t + y_i sin t)^6
= g^6`. The general component has no rational conic (elliptic curve of
conductor 4070, rank 0). What remains open in degree 2 is the design
stratum: five rational points `omega_i` on the unit circle with `e2 = 0`
and `|e1| < 1`, whose Szegő weights `w_i` satisfy `w_i / |gamma_i|^6` all
in one class of `Q^*/Q^{*6}`. Rational positive designs exist (30 up to
Gaussian norm 400, listed in `results/conics_2026_09_06/design_20.log`).

You own the following, in priority order.

### Lane 1: genus-one curves of low degree on X (Elkies mechanism)

Expected dimension 0, so they are isolated and potentially defined over Q.

1a. Plane cubics. A plane P in P^5 with X ∩ P = C1 ∪ C2, two plane cubics.
    Parametrise planes through the coordinates, restrict the Fermat sextic,
    and demand the plane sextic factor as cubic × cubic. Set this up as a
    polynomial system in the Plücker/affine plane coordinates and solve
    numerically (homotopy or many-start Newton), then refine and test for
    rationality. A rational plane cubic with a rational point and positive
    rank yields infinitely many rational points; then check the positive
    chamber.
1b. Elliptic quartics: intersections of two quadrics in a P^3 ⊂ P^5 lying
    on X. Same recipe, larger system.

For any genus-one curve found over Q: compute its Jacobian, rank, and
generators (Sage or PARI if available; otherwise mwrank via pip if the
environment allows), then enumerate points and test positivity.

### Lane 2: degree-4 rational curves under symmetry

Coordinates p_i(t) of degree 4, definite. Normalise a1 = 1 and PGL2 as in
the conic case. Impose symmetry: (12)(34) with t -> -t, p5 p6 even; then
(12)(34) with t -> -t plus (13)(24) with s <-> t. Count unknowns versus the
25 coefficient equations and solve the 0- or 1-dimensional systems the same
way. Report all real definite solutions and any rational ones.

### Lane 3: the design stratum, exactly

3a. Prove or refute that every general conic on X is reflection-symmetric
    (the one unproved step in CONIC_LANE.md §2b). Approach: show that the
    Fourier system M1..M6 with some `z_i != 0` forces the pairing
    `(zeta, z), (-zeta, z)`; or exhibit a numerical counterexample with the
    solver in `tools/conics/conic_general.py`.
3b. Describe the moduli surface S of rational 5-node designs (`e2 = 0`,
    `|e1| < 1`) and its sixth-power cover. Find a rational curve on S if
    one exists, restrict the class conditions to it, and decide whether the
    cover has rational points. Extend `design_search.py` beyond norm 400
    with a faster integer implementation.
3c. Six nodes: the same construction with six terms gives conics on the
    (6,1,6) fivefold, also open. There the kernel condition is a single
    determinant and the moduli is larger; search rational 6-node positive
    designs with sixth-power weights.

## Hard rules

- Never report a solution from floating point. A candidate curve must be
  verified as an exact identity in Q[t] (sympy `expand` or fraction-free
  arithmetic). A candidate sextuple must pass
  `python3 tools/verify_esop6.py a b c d e f` and
  `node tools/verify_esop6.mjs a b c d e f` in the repo, and be committed
  with both outputs before any announcement.
- Record every solved system: unknown count, equation count, Jacobian rank
  at solutions, number of distinct complex solutions, real solutions,
  definite solutions, rational solutions. Negative results with these
  counts are the deliverable if no curve appears.
- Do not run integer sweeps, do not touch the 4.31M frontier, do not work on
  the surface slices.
- Commit to a branch named `astra/curves-<lane>` and write results under
  `results/astra_curves_<date>/`.
