# Structural arithmetic

Let \(g=\gcd(f,t)\), \(f=gx\), and \(t=gy\), so \(\gcd(x,y)=1\).  Since
`f` and `t` are units modulo 42, `x` and `y` are odd and coprime to 3 and 7.

## Factor gcds

Write

\[
A=x-y,\quad B=x+y,\quad C=x^2+xy+y^2,\quad D=x^2-xy+y^2.
\]

Then \(x^6-y^6=ABCD\), and elementary substitution gives:

- \(\gcd(A,B)\mid2\);
- \(\gcd(A,C)\mid3\), while \(\gcd(A,D)=1\);
- \(\gcd(B,D)\mid3\), while \(\gcd(B,C)=1\);
- \(\gcd(C,D)=1\) for coprime odd `x,y`.

Restoring `g`, common factors outside powers of `g` can occur only through 2
and 3 in the stated adjacent pairs.  The quadratic factors are Eisenstein
norm forms.  If a prime \(p\equiv2\pmod3\) divides one of them and
\(\gcd(f,t)=1\), its exponent is even; otherwise it must divide both `f` and
`t`.  These facts make factorization certificates compact, but by themselves
do not constrain an additive four-sixth-power representation enough to give
a descent.

## Valuation lemma for the reduced target

**Lemma.** If a positive integer \(S\) is a sum of at most four positive
sixth powers, then

\[
v_2(S)\bmod6\in\{0,1,2\},\qquad
v_3(S)\bmod6\in\{0,1\},\qquad
v_7(S)\equiv0\pmod6.
\]

**Proof.** Sixth powers are 0 or 1 modulo 8, 9, and 7 according as the base is
divisible by 2, 3, and 7.  With at most four terms, if the sum is 0 modulo
one of those moduli, the number of unit terms must be zero.  Hence every base
is divisible by the corresponding prime and the whole sum is divisible by
its sixth power.  Divide by that sixth power and repeat.  For 2, residual
valuation 1 or 2 can arise from two or four odd terms; residual 3–5 cannot.
For 3, residual 1 can arise from three unit terms; residual 2–5 cannot.  For
7, no positive residual below 6 can arise.  ∎

This is implemented by `--valuation-prune`.  On the 700k–730k NB=1 control it
rejects 184 decomposition calls and reduces j=2 nodes from 4,139,088 to
4,124,811 (0.34%).  It is mathematically exact but not a major speedup on that
band.  Its best use is early certification and on valuation-stratified ranges.

## Strongest p-adic pruning already available

For the four small terms in the concentrated case, let the minimum base
valuations be \(r_2,r_3,r_7\).  Then

\[
v_p(f^6-t^6)\ge6r_p.
\]

LTE/root-group analysis pins `f` to a thin p-adic ray through `t`; in
particular, extra factors 49 in all four small terms push the 7-adic modulus
to \(7^{12}\), far beyond the present frontier.  This is stronger for search
partitioning than trial-factorizing every `f^6-t^6`.

## What factorization does not prove

Every factor of `f^6-t^6` divides the additive target, but a divisor of a sum
of four sixth powers need not divide any summand.  Assigning the four
cyclotomic factors to `a,b,c,d` is therefore invalid without a new theorem.
Likewise, local sixth-power residue conditions exist at every fixed modulus;
no finite CRT cascade proves global nonexistence.

## Recommended arithmetic lane

Stratify candidates by the exact tuple
`(v2(f-t), v2(f+t), v3(cyclotomic branch), v7(root branch), gcd(f,t))`, apply
the valuation lemma to the reduced target, and measure whether rare high-
valuation strata admit a rigorous lower bound on `f` or a smaller-solution
descent.  A descent—not another residue table—is the desired outcome.
