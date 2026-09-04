# ESOP6 direct solution strike — 2026-09-04

```text
ESOP6 SOLUTION FOUND: NO
POSITIVE RATIONAL SURFACE POINT: NO
NEW RATIONAL CURVE: NO
NEW ELLIPTIC ATTACK: NO
TARGETED INTEGER SEARCH RUN: YES
```

These statuses concern this session. A curve or point on a quotient is not
counted as a curve or point on the target surface.

Base commit: `485390a17a1ae618f7f00cc547c94fa154144ebd`, fetched and matched to
`origin/main` before work. The three repository commits were read in order,
followed by the September 4 handoff, arithmetic, geometry, and engine notes.

## What was attacked

The primary target was exactly

\[
2X^6+2Y^6+Z^6=W^6,\qquad XYZW\ne0,
\]

with rational coordinates. Three concrete attacks were completed:

1. The cubic quotient and its rational conic pencil, including degenerations,
   low-degree curves, boundary contact, and exact square lifting.
2. Four conics selected using newly derived coupled valuation restrictions,
   followed by finite local tests and an explicitly bounded rational-parameter
   search.
3. A sparse direct integer search using the factorization of `W^6-Z^6` and
   an exact divisor/discriminant solver for `x^6+y^6=N`.

The last method avoids a pair table and a scan over the possible pair bases
for each factored surface target. It is a useful new bounded architecture
for the repeated-coordinate surface. It is **not** an exponent improvement
for the unrestricted six-variable problem; factorization cost is retained.
This is the scope of stop condition **E** used here. The selected batches
were completed; no exhaustive classification of all promising attacks is
claimed.

## Exact control gates

Executed, in the requested order:

```bash
make clean
make -j"$(nproc)"
make validate control equiv frontier-audit prototypes
make differential
```

All passed. The retained log is
`results/astra_direct_2026_09_04/control_gates.log`.

- Lander–Parkin fifth-power solution recovered exactly.
- 700k–730k control: 124 candidates, zero reported solutions.
- Effective oracle nodes: 4,139,088 in both the monolithic and NB=8 runs.
- All eight candidate counts and both digests matched; production total 55,684.
- Both prototypes passed their positive controls and exact checks.
- All 18 differentials passed for NB=1,2,3,4,7,8 on the three fixed bands.

The historical 4.3M decomposition campaign was not rerun, and its evidence
grade was not changed.

## Strongest mathematical results

Details and proofs are in [GEOMETRIC_STRIKE.md](GEOMETRIC_STRIKE.md).

- The natural conic pencil on the cubic quotient generically lifts to **genus
  9**, rather than genus 1. Its rational exceptional parameter `lambda=2`
  gives a real sum-of-squares obstruction. The complex elliptic components
  there do not furnish a positive solution.
- Reduced rational parametrizations of degrees **1, 2, and 3** cannot supply
  a nonconstant curve on this real surface. Degree 2 is excluded by the known
  degree bound for curves on diagonal sextic surfaces. General degree 4
  remains **open in this session**.
- Any nonconstant curve through either obvious boundary point
  `[0:0:1:1]`, `[0:0:1:-1]` needs parameter degree at least **6**.
- The exactly diagonal restriction `X=Y` has no positive rational solution,
  by reduction to Fermat's theorem for cubes.
- For a primitive positive point, put
  `lambda=(W^2-Z^2)/(X^2+Y^2)` and `r_p=min(v_p(X),v_p(Y))`.
  Then `r_2,r_3,r_7 >= 1` and

  ```text
  v2(lambda) = 4*r2 + 1
  v3(lambda) = 4*r3 - 1
  v7(lambda) = 4*r7 OR -2*r7
  ```

  These are necessary conditions, not a global obstruction. In particular,
  the initially tested `lambda=1/2` conic is entirely excluded, independently
  of its finite parameter search.

## Bounded searches actually completed

| Domain | Exact amount tested | Result |
|---|---:|---|
| `lambda=1/2`, coprime `1<=u,v<=2000` | 2,433,175 parameter pairs; 1,582,669 positive quotient images | no square lift; subsequently excluded for all parameters by valuations |
| `lambda=864/931, 864/1225, 864/1519, 864/1813`; coprime `-2000<=u<=2000`, `1<=v<=2000`, plus infinity | 19,465,408 parameter occurrences; 7,942,986 positive quotient images | no square lift |
| Specified first/last points on CRT rays for 96 selected values of `W` | 10,382 distinct `(W,Z)` targets | no positive surface point |

Parameter occurrences are not counts of distinct projective quotient points.
The conic searches use one implementation of the parameter box. Exact
symbolic identities and a separately implemented exhaustive local-control
oracle check their arithmetic ingredients. Their limits are parameter-height
limits, not bounds on all points of the surface.

The sparse integer batch selected 16 unit values of `W` from each of six
declared bands, ending at `10^12`, using seed `660604`. For each of the 144
sixth roots of unity modulo `2*42^6`, only the first and last positive `Z<W`
on that ray were included. Duplicate pairs were removed. Actual tested
`W` ranged from 7,793,543 to 965,959,379,131. This is **not** an empty frontier
through that upper value.

| Disposition | Targets |
|---|---:|
| two-sixth-power valuation rejection | 5,820 |
| pair-residue rejection | 4,278 |
| certified prime `p=3 mod 4` with exponent not divisible by 6 | 283 |
| complete divisor/discriminant test | 1 |
| **Total** | **10,382** |

The divisor oracle passed 849 comparisons with an independent monotone
search, including 465 planted positive two-sixth-power targets. A separate
standard-library-only replayer verified all 10,382 dispositions, the exact
domain, and 4,090 prime certificates, including proof dependencies and
control primes. Full target rows and certificates are retained.

No unrestricted four-part offline join, GPU search, or historical-frontier
extension was run. The integer batch retained the highest-priority repeated
surface, where the four-part core becomes two sixth powers exactly.

## Verification and reproduction

```bash
python3 tools/astra_geometry.py
python3 tools/discover_conic_seeds.py
python3 tools/surface_lift_search.py --height 2000
python3 tools/admissible_conic_strike.py --height 2000
python3 tools/surface_divisor_strike.py --per-band 16
python3 tools/replay_surface_divisor.py results/astra_direct_2026_09_04
python3 tests/direct_verifiers.py
```

The factor search used `python-flint`; SymPy is a supported fallback. The
independent replay requires only Python's standard library. The geometry and
seed-discovery scripts require SymPy. Runtime details and file hashes are in
`results/astra_direct_2026_09_04/manifest.json`.

Two standalone outsider verifiers now accept six decimal integers:

```bash
python3 tools/verify_esop6.py a b c d e f
node tools/verify_esop6.mjs a b c d e f
```

They print every sixth power, positivity, the exact two sides, equality,
gcd, the sorted primitive tuple, and its recomputed equality. Both reject
zero-coordinate arithmetic equalities as solutions. All 18 cross-language
controls passed. No known positive ESOP6 solution was available as a fixture.

## The single strongest next direct strike

Construct a **degree-six curve with order-six contact at `[0:0:1:1]`** in the
explicit normalized family

\[
X=tA_5(t),\quad Y=tB_5(t),\quad
W=Z+\tfrac23t^6,
\qquad A_5(0)=B_5(0)=Z(0)=1,
\]

where `deg A5,deg B5<=5`, `deg Z<=6`, and `A5,B5` are independent.
Impose the surface identity after dividing out `t^6`; use the contact
conditions before solving coefficients. A successful identity immediately
has positive rational specializations for sufficiently small rational
`t>0`. Do not impose the extra two-boundary-point parity symmetry: the
corresponding subfamily was excluded by a sixth-power coefficient obstruction.

This is a concrete next construction, not a claim that it exists. It is
chosen because the contact obstruction proves why the lower-degree boundary
attempts cannot work, while a successful normalized identity would directly
produce the six required positive integers.
