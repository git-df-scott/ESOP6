# Rational-point ledger — direct strike, 2026-09-04

No positive rational point on `2X^6+2Y^6+Z^6=W^6` was found. No row below
is an unverified solution candidate. Coordinate types are explicit so that
a cubic-quotient point cannot be mistaken for a sixth-power point.

## Exact points and identities

| Object | Coordinates or parameter | Exact status |
|---|---|---|
| Surface boundary B+ | `[X:Y:Z:W]=[0:0:1:1]` | valid, but X=Y=0; not a solution |
| Surface boundary B- | `[0:0:1:-1]` | valid boundary point; sign change does not remove zeroes |
| Cubic quotient conic | `lambda=1/2`, `[U:V:R:T]=[t^2+8t-5:-t^2+8t+5:2(t^2-2t+5):2(t^2+2t+5)]` | cubic identity verified; no rational square lift by the valuation obstruction |
| Degenerate residual fibre | `lambda=2` | complex genus-one components; no nonzero real lift |
| BCU polynomial identity | `(1-t-t^2)^3+(1+t-t^2)^3=2-2t^6` | symbolic control passes; not an ESOP6 point |

The conic seed triples below use `(A,B,D)`, where `A=U+V`, `B=U-V`,
`D=T+R`, and `T-R=lambda*A`. They solve
`(2-lambda^3)A^2+6B^2=3lambda D^2`. They are quotient data, not square roots.

| lambda | Exact conic seed `(A,B,D)` | Positive quotient images in tested box | Square lifts |
|---|---|---:|---:|
| 864/931 | (938448,607860,1084489) | 715,087 | 0 |
| 864/1225 | (1543500,373248,1500625) | 2,350,621 | 0 |
| 864/1519 | (109368,13428,115601) | 2,621,648 | 0 |
| 864/1813 | (456876,257400,744485) | 2,255,630 | 0 |
| **Total** | | **7,942,986** | **0** |

Each box is coprime `-2000<=u<=2000`, `1<=v<=2000`, plus `[u:v]=[1:0]`:
4,866,352 parameter occurrences per conic. Repeated images are possible.
No claim is made beyond these boxes. Finite local checks used moduli
`2^8,2^12,3^5,3^7,7^3,7^4,13^2,19^2,31^2` and left survivors in every case.
The local predicate was independently checked on all 13,058 four-coordinate
residue vectors modulo 8,9,7 using direct enumeration of unit scalars.

Seed requests at `864/833,864/1127,864/1421,864/2009` returned no point from
the conic solver. These are **solver reports**, not independently certified
empty conics; no further box search was run on them.

The earlier `lambda=1/2` box contained 2,433,175 coprime positive parameter
pairs and 1,582,669 positive quotient images, with zero square lifts. The
subsequently proved valuation obstruction supersedes that finite negative.

## Sparse integer target ledger

The complete domain and all 10,382 target dispositions are retained in:

```text
results/astra_direct_2026_09_04/divisor_domain.json
results/astra_direct_2026_09_04/divisor_targets.jsonl
results/astra_direct_2026_09_04/prime_certificates.json
results/astra_direct_2026_09_04/divisor_independent_replay.json
```

Domain SHA-256 (sorted lines `W,Z\n`):

```text
5fa7e950eb948d8f48ecf2c9360e271c8d69c22f4c7bfb34180c85f3726a3e14
```

Only one target survived the preliminary and certified inert-prime filters:

```text
W = 539374987829
Z = 3279554165
N = (W^6-Z^6)/(2*42^6)
  = 2242947361610075429337818084750779759088179669366859867748217
```

It had one divisor in `N<h^3<=4N`; the exact reconstruction produced no
positive pair `x^6+y^6=N`. This is a rejected target, not a near-solution.

No sextuple reached the solution-verification stage. The standalone Python
and JavaScript verifiers are ready and were cross-checked on 18 controls,
including boundary equalities that correctly fail positivity.

## Open items preserved honestly

- General rational quartics avoiding B+ and B- were not classified.
- Other boundary rational points were not exhaustively classified.
- The four admissible conic members are not globally closed.
- No positive-point elliptic curve or Mordell–Weil route was obtained.
- The sparse batch does not give a continuous height exclusion.
- The historical 4.3M frontier retains its existing grade-B qualification.

The next direct construction is the independent degree-six boundary-contact
ansatz stated in [ASTRA_DIRECT_STRIKE.md](ASTRA_DIRECT_STRIKE.md).
