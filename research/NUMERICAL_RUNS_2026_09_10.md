# Executed numerical construction searches

No exact rational identity or ESOP6 counterexample was found.
These are bounded local optimizations, not exhaustive coefficient searches.

## Degree eight, repeated pairs

The exact residual is

    (P^6+R^6)/2-Q*N^5-(10/27)t^12 Q^3*N^3-(1/81)t^24 Q^5*N,
    Q=1+s^2 t^2, P(0)=R(0)=N(0)=1.

The two necessary linear remainder relations at t=i/s eliminate P6,P7.
Twenty real variables remain; all 42 nonconstant residual coefficients
are fitted using an analytic Jacobian, checked against directional finite
differences at every start.

`degree8_2026_09_10/search.py`: s=1/2,1,2, both remainder signs,
four starts per pair, RNG seed6102026, maximum600 evaluations per start.
All24 starts completed. Minimum maximum absolute residual coefficient
was 0.004753005309566127. Two starts reached their evaluation cap.

The s=2 coordinates were poorly conditioned. A separate replay in u=2t
used eight starts and seed6102027. All eight reached the600-evaluation
cap; the minimum maximum coefficient residual in these rescaled coordinates
was3.567777906859995e-6. This number is not directly comparable with the
unscaled residual. None reached the1e-8 rational reconstruction threshold.
The replay is retained as `degree8_2026_09_10/replay_normalized.py`.

## Full-fourfold conics

`full_conics_2026_09_10/search.py` fits every coefficient in

    sum_i X_i(t)^6=(1+q*t^2)^6,
    X_i=a_i+b_i*t+c_i*t^2,

for q=1,2,3, with a5=1/2,b5=0 fixed. Thus13 real variables satisfy13
coefficient equations. Fixing these coefficients is a restricted family,
not a normalization proved to cover every rational conic. Only quadratic
coordinate forms are attempted; no genuine quartic or cubic search is claimed.

There were12 perturbations of an explicitly known algebraic real conic,
plus12 independent Gaussian starts (`--random`), each capped at1000
evaluations. The analytic Jacobian was directionally checked at every start.
The smallest maximum residual was3.552713678800501e-15; these small floating
residuals establish neither rationality nor exactness.

For runs with maximum residual below1e-8, each coefficient was individually
approximated by a rational with denominator at most10,100,1000,10000.
All60 reconstructed polynomial tuples failed exact SymPy coefficient
verification. This reconstruction procedure is not an exhaustive search
over rational coefficient tuples within those bounds.

Each JSON ledger retains the full fitted coefficients, residual, solver
status and evaluation count. No fitted conic was promoted to a rational
curve or to an integer candidate. The known real seed and its exact
irrationality obstruction are in the independent-geometry report.

Dependencies: Python3, NumPy, SciPy, SymPy. SymPy was installed in this
execution pass. No generic integer frontier was extended.
