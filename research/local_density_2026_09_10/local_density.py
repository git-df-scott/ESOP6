"""Exact primitive congruence counts for x1^6+...+x5^6=x6^6."""
from fractions import Fraction
from pathlib import Path
import json

def histogram(q, p=None):
    h=[0]*q
    for x in range(q):
        if p is None or x%p==0: h[pow(x,6,q)]+=1
    return h

def convolution(a,b,q):
    out=[0]*q
    aa=[(i,v) for i,v in enumerate(a) if v]
    bb=[(i,v) for i,v in enumerate(b) if v]
    for i,v in aa:
        for j,w in bb: out[(i+j)%q]+=v*w
    return out

def count(h,q):
    sums=h
    for _ in range(4): sums=convolution(sums,h,q)
    assert sum(sums)==sum(h)**5
    return sum(a*b for a,b in zip(sums,h))

def row(p,k):
    q=p**k
    total=count(histogram(q),q)
    nonprimitive=count(histogram(q,p),q)
    primitive=total-nonprimitive
    delta=Fraction(primitive,q**5)
    return dict(p=p,k=k,modulus=q,total=total,nonprimitive=nonprimitive,
                primitive=primitive,delta_numerator=delta.numerator,
                delta_denominator=delta.denominator,delta_decimal=float(delta))

def primes(n):
    return [p for p in range(2,n+1) if all(p%d for d in range(2,int(p**0.5)+1))]

good=[row(p,1) for p in primes(199) if p not in (2,3)]
bad=[row(p,k) for p,maxk in ((2,10),(3,6)) for k in range(1,maxk+1)]
for r in bad:
    if r['p']==2 and r['k']>=3:
        assert Fraction(r['primitive'],r['modulus']**5)==Fraction(5,8)
    if r['p']==3 and r['k']>=2:
        assert Fraction(r['primitive'],r['modulus']**5)==Fraction(20,81)
    if r['k']>=7:
        earlier=next(a for a in bad if a['p']==r['p'] and a['k']==r['k']-6)
        assert r['nonprimitive']==r['p']**30*earlier['total']
for r in good:
    assert r['nonprimitive']==1
    if r['p']%6==5:
        # Sixth powers are squares with identical multiplicities; the split
        # six-variable quadratic has p^5+(p-1)p^2 zeros in this case.
        p=r['p']; assert r['total']==p**5+(p-1)*p**2

# Independent brute-force control at q=2,3,4. No convolution reuse.
import itertools
for q,p,k in ((2,2,1),(3,3,1),(4,2,2)):
    brute=sum(1 for xs in itertools.product(range(q),repeat=6)
              if (sum(x**6 for x in xs[:5])-xs[5]**6)%q==0
              and any(x%p for x in xs))
    assert brute==next(r['primitive'] for r in bad if r['p']==p and r['k']==k)

# Good-prime lift checks, calculated afresh at p^2.
lift_checks=[row(p,2) for p in (5,7,11,13)]
for r in lift_checks:
    base=next(a for a in good if a['p']==r['p'])
    assert Fraction(r['primitive'],r['modulus']**5)==Fraction(base['primitive'],base['modulus']**5)

product=Fraction(1)
for p in (2,3):
    r=max((r for r in bad if r['p']==p),key=lambda r:r['k'])
    product*=Fraction(r['delta_numerator'],r['delta_denominator'])
for r in good: product*=Fraction(r['delta_numerator'],r['delta_denominator'])
data=dict(equation='x1^6+x2^6+x3^6+x4^6+x5^6=x6^6',
    definition='delta(p,k)=primitive_solution_count_mod_p^k / p^(5k)',
    good_primes=good,bad_prime_powers=bad,good_prime_lift_checks=lift_checks,
    finite_product_through_199=float(product),
    warning='Finite local product only: neither rational existence nor a height forecast.')
Path(__file__).with_name('results.json').write_text(json.dumps(data,indent=2)+'\n')
for r in bad: print(r)
print('good prime rows:',len(good),'finite product:',float(product))
print('all independent brute-force and good-prime lift controls passed')
