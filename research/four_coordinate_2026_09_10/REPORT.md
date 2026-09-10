# Four independent coordinates: exact denominator obstruction

2026-09-10. No counterexample found. This is a structural computation and proof,
not an integer-height search and not an exclusion of the complete rational family.

We study the proposed identity

    sum(A_i(t)^6, i=1..4)
      = 6 c M(t)^5 + 5 c^3 t^12 M(t)^3 + (3/8)c^5 t^24 M(t),

with rational coefficients, deg A_i <= 5, deg M <= 6, M(0)=1,
(A_1(0),...,A_4(0))=(1,2,4,3), and c=815. The constant equation is
1+64+4096+729 = 4890 = 6*815. Reconstruction is
x_i=t A_i, z=M-c t^6/2, w=M+c t^6/2.

## Result 1: only the full degree branch survives

Let k=max deg A_i and d=deg M. For d<6, the RHS has degree 24+d:
the last term strictly dominates the other two terms. The LHS has degree 6k,
because real sixth-power leading coefficients cannot cancel. Thus d=0,k=4.
For d=6, its leading coefficient is

    c*m6 * (6*m6^4 + 5*c^2*m6^2 + 3*c^4/8),

which is nonzero. Hence k=5. It is positive only when m6>0.

In the first case M=1. Its leading coefficient equation is

    sum(a_i4^6) = 3*c^5/8.

A nonzero sum of at most four rational sixth powers has 2-adic valuation
6a+e with e in {0,1,2}: multiply by a sixth power to make all entries
2-integral with at least one odd; the sum modulo 8 equals the number of odd
entries, which is 1,2,3,or 4. But v2(3*815^5/8)=-3 is 3 modulo 6.
Contradiction. Therefore every solution needs deg M=6, max deg A_i=5,
and positive leading coefficient of M.

## Result 2: the midpoint must have a coefficient with an even denominator

Let v_G be the ordinary coefficient Gauss valuation at 2 on Q[s].
For any nonzero sum of four polynomial sixth powers,

    v_G(sum(P_i^6)) = 6 min_i v_G(P_i) + e,  e in {0,1,2}.

Here is a complete proof of the needed polynomial version. Normalize the
P_i to be in R=Z_(2)[s] with at least one unit coefficient. Put U_i=P_i^3;
at least one U_i remains nonzero modulo 2. Suppose sum U_i^2 vanishes
modulo 8. Modulo 2, U_4=U_1+U_2+U_3. Write this as an exact equality plus
2V in R. Dividing the square sum by 2 and reducing modulo 2 gives

    U_1^2+U_2^2+U_3^2+U_1 U_2+U_1 U_3+U_2 U_3 = 0.

This is X^2+XY+Y^2=0 with X=U_1+U_3, Y=U_2+U_3.
The quadratic is anisotropic over F_2(s): a nonzero solution would put a
root of T^2+T+1 in F_2(s), whereas the algebraic constants of F_2(s) are
exactly F_2. Consequently all four U_i have the same reduction U.
Write U_i=U+2V_i using any common lift. Divide their square sum by 4 and
reduce modulo 2. With V=sum V_i, the equation is now

    U^2+UV+V^2=0.

The same anisotropy forces U=0, contradicting primitivity. Thus the square
sum has Gauss valuation at most 2, proving the lemma for sixth powers.
This argument concerns ordinary Gauss valuation over Q; it does not assert
the same statement after arbitrary ramified field extensions.

Invert the parameter:

    P_i(s)=s^5 A_i(1/s),  m(s)=s^6 M(1/s).

Then m is monic of degree 6, and

    sum P_i^6 = 6c m^5 +5c^3 m^3 +(3/8)c^5 m.

If all coefficients of M are 2-integral, v_G(m)=0. The three RHS summands
have valuations 1,0,-3 respectively, so its valuation is exactly -3.
The lemma excludes this value on the LHS. Therefore M necessarily has a
coefficient of negative 2-adic valuation. This result allows *arbitrary*
rational coefficients in the A_i; it is not merely an integral-input test.

## Remaining denominators: not excluded

Write b=v_G(m)<=0 and a=min_i v_G(P_i). The result excludes b=0.
For b<=-2, the first RHS term is uniquely dominant, giving

    6a+e=5b+1,  e in {0,1,2}.

Thus b modulo 6 must be 0,1,or 5. In particular b=-2,-3,-4 are excluded;
b=-5,-6,-7 can satisfy this necessary valuation condition.

At b=-1 the first and last RHS terms tie, so the unique-dominance argument
does not apply. Indeed, putting n=2m gives the exact factorization

    RHS = (c/16) n (3n^2+c^2)(n^2+3c^2).

This branch remains open. No claim of rational nonexistence follows from
the denominator obstruction alone.

## Executed checks

Run `python3 research/four_coordinate_2026_09_10/check.py`.

The executable verifies the midpoint identity and both mod-2 reduction
identities symbolically. An exact meet-in-the-middle enumeration covers
all 64^4 = 16,777,216 ordered quadruples of degree-at-most-one polynomials
modulo 8. There are zero primitive sixth-power zero sums; the 65,536
nonprimitive zero sums are exactly those with every coefficient even.
This is a finite regression for the algebraic lemma, not a search through
16 million potential ESOP6 counterexamples. The general-degree conclusion
comes from the proof above. Results are retained in checks.json.

## Search implication

The independent-coordinate relaxation escapes the earlier two-coordinate
3-adic obstruction but acquires a demonstrable 2-adic denominator condition.
Any new exact coefficient solve should retain b=-1 explicitly and must not
restrict the midpoint to coefficients with odd denominators. The full
rational family and the ESOP6 counterexample target remain unresolved.
