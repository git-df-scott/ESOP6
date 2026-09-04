# Time scaling and the exponent wall

Let \(F=f_{max}\), \(M=42^6\), and \(B=\lceil F/42\rceil\).

## Current algorithm from first principles

1. There are \(\Theta(F)\) possible `f` values and 144 sixth roots modulo
   \(M\).
2. For \(F<M\), a root-generated residue `t` is below `f` with probability
   \(\Theta(f/M)\).  Summing over `f` gives
   \(K(F)=\Theta(F^2/M)\) concentrated candidates.  The observed 55,684 count
   is consistent with this quadratic law.
3. The pair Bloom contains \(P=B(B+1)/2=\Theta(F^2)\) pair sums and takes
   \(\Theta(F^2)\) build time.
4. For one four-part target, DFS chooses the first part in \(O(B)\), the
   second in \(O(B)\), and sends the remaining two parts to the oracle.
   Residue and class filters reduce constants but not this worst-case
   \(O(B^2)\) traversal.
5. Total traversal is therefore
   \(K(F)O(B^2)=O(F^4/M)\), customarily abbreviated \(O(F^4)\).

With NB buckets, pair construction and upper DFS traversal are repeated NB
times: worst-case \(O(NB\,F^4/M)\).  Only \(1/NB\) of leaf queries reaches the
Bloom filter, which explains the much smaller measured penalty on the control.

## False-positive tail

At fixed Bloom bits-per-pair, false positives occur with fixed probability.
Each exact `pair_verify` can scan \(O(B)\) candidate first parts.  Thus the
strict asymptotic upper model contains an \(O(p_{fp}F^5/M)\) term.  At 12 bpp
and eight hashes, the idealized Bloom rate is about 0.00314; at 16 bpp it is
about 0.000574.  The campaign's observed regime is traversal-dominated, but a
10M projection must measure exact-verifier time rather than assume pure
\(F^4\).

## Why the standard alternatives do not remove an exponent

Write \(S=\{x^6+y^6:1\le y\le x\le B\}\), with \(|S|=\Theta(F^2)\), and let
\(T\) be the \(\Theta(F^2/M)\) reduced targets.  The question is whether
\(s_1+s_2=t\) for some \(s_1,s_2\in S,t\in T\).  This is an offline 3SUM-type
join.

| Method | RAM | Time | I/O | Exactness | Exponent result |
|---|---:|---:|---:|---|---|
| current Bloom + DFS | \(O(F^2)\) bits | \(O(F^4/M)\) traversal | none | exact after verification | baseline |
| caseA3 NB passes | \(O(F^2/NB)\) bits | worst \(O(NB F^4/M)\) | none | exact | memory only |
| exact hash pair table | \(O(F^2)\) words | \(O(F^4/M)\) probes | none | exact | constant only |
| sorted pair table/two-pointer per target | \(O(F^2)\) words | \(O(F^4/M)\) | optional | exact | no |
| external routed join | \(O(F^2/NB)\) words | sort/join plus query generation | \(O(F^2+Q)\) | exact | avoids repeat passes; `Q` remains \(F^4\) |
| mmap Bloom | page-cache dependent | one build + one traversal | random \(O(P+Q)\) page traffic | exact after verification | no |
| modular partition | divided by residue density | same power of F | low | exact if all classes covered | constant only |
| GPU fused probes | device-table dependent | same candidate/probe exponent | PCIe/build cost | exact with host recheck | throughput only |
| dense FFT/NTT convolution | \(O(F^6)\) address range | at least \(O(F^6\log F)\) | enormous | exact with CRT | worse |

The routed prototype proves that pair and query streams can each be assigned
to one bucket in one pass.  It does not solve the dominant problem that there
are \(Q=\Theta(F^4/M)\) leaf queries to route.  Claims of an exponent drop must
address this `Q`, not just pair-table construction.

## Strongest algorithmic attack on F^4

Treat the whole band as one structured offline join rather than independent
targets.  The exact research problem is:

> Exploit the algebraic target family
> \(T=\{(f^6-t^6)/42^6:t\equiv\omega f\pmod{42^6}\}\) to solve
> \((S+S)\cap T\) in \(o(|S||T|)\) time and \(O(|S|)\) or less space.

Generic sorting/hashing does not achieve this.  A genuine improvement must
use the polynomial/root structure of `T`, a new additive-energy bound for
sixth-power pair sums, or a sparse-convolution algorithm with a proved bound.
Until that is derived, GPU/routing work is an excellent constant-factor lane,
not an asymptotic breakthrough.

## Calibration policy

Use the smallest NB that keeps peak memory safely below the machine limit;
benchmark 10k–50k-wide overlapping bands at the intended `F`.  Record pair
build, upper DFS, Bloom queries, exact-verifier calls, and wall time separately.
Do not extrapolate from 700k solely by \((F/700k)^4\).
