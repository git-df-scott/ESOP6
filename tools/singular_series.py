#!/usr/bin/env python3
"""Hardy-Littlewood heuristic constant for  a_1^k + ... + a_n^k = f^k.

For the (k,1,n) equation the heuristic count of PRIMITIVE solutions with
f <= F is  N(F) ~ C log F,  with

    C = (rho / n!) * prod_p sigma_p,
    rho = (n/k) * Gamma(1+1/k)^n / Gamma(1+n/k)      (real density),
    sigma_p = lim_j  N*_j(p) / p^{(n)j}                (p-adic density),

where N*_j(p) counts PRIMITIVE (not all divisible by p) solutions of the
congruence modulo p^j.  Primitivity matters: a single primitive solution
yields F/f imprimitive multiples, so the log law is a statement about
primitive tuples only.

For p not dividing k the congruence is smooth on primitive tuples and
j=1 is exact.  For p | k we compute increasing j until the ratio is
stationary and report the sequence so the stabilisation is visible.

Distributions are convolved exactly with Python integers (no floats
inside any count).  Only the final product is a float.
"""
import argparse, json, math, sys
from collections import Counter

def primes_upto(n):
    s = bytearray([1])*(n+1); s[0]=s[1]=0
    for i in range(2,int(n**.5)+1):
        if s[i]: s[i*i::i]=bytearray(len(s[i*i::i]))
    return [i for i in range(n+1) if s[i]]

def power_dist(m,k):
    """Counter r -> #{x mod m : x^k = r}."""
    c=Counter()
    for x in range(m): c[pow(x,k,m)]+=1
    return c

def convolve(a,b,m):
    """Exact cyclic convolution of two Counters modulo m (list output)."""
    out=[0]*m
    bl=[0]*m
    for r,v in b.items(): bl[r]=v
    for r,v in a.items():
        if v==0: continue
        # out[(r+s)%m] += v*bl[s]
        for s in range(m):
            w=bl[s]
            if w: out[(r+s)%m]+=v*w
    return out

def count_solutions(m,k,n):
    """N_j = #{(a_1..a_n,f) mod m : sum a_i^k = f^k}, exact integer."""
    d=power_dist(m,k)
    acc=Counter(d)
    for _ in range(n-1):
        acc=Counter({r:v for r,v in enumerate(convolve(acc,d,m)) if v})
    return sum(acc[r]*d[r] for r in d)

def primitive_count(p,j,k,n,memo):
    """Solutions mod p^j with not all coordinates divisible by p."""
    m=p**j
    tot=count_solutions(m,k,n)
    # tuples with all n+1 coordinates = 0 mod p: x=p*y, y mod p^{j-1}
    if j<=k:
        allzero=p**((n+1)*(j-1))          # equation automatic mod p^j
    else:
        # sum (p y_i)^k = (p g)^k mod p^j  <=>  sum y_i^k = g^k mod p^{j-k}
        # y ranges mod p^{j-1}: each residue mod p^{j-k} lifts p^{(k-1)(n+1)} ways
        allzero=count_solutions(p**(j-k),k,n)*p**((k-1)*(n+1))
    return tot-allzero

def sigma_p(p,k,n,maxj):
    seq=[]
    for j in range(1,maxj+1):
        m=p**j
        if m>3000: break
        N=primitive_count(p,j,k,n,None)
        seq.append((j, N, N/(p**(n*j))))
    return seq

def sigma_smooth_fft(p,k,n):
    """For p not dividing k the congruence is smooth on primitive tuples,
    so j=1 is exact.  Counts are computed with float FFT: every
    intermediate is below 2^53 for p<=2^10 and relative error is ~1e-13
    beyond, far below the accuracy we report."""
    import numpy as np
    d=np.zeros(p)
    x=np.arange(p,dtype=object)
    for r in (pow(int(v),k,p) for v in range(p)): d[r]+=1
    D=np.fft.rfft(d)
    tot=np.fft.irfft(D**n,p)
    N=float(np.rint(tot@d))
    N-=1.0   # the all-zero tuple is the only imprimitive one mod p
    return N/p**n

def real_density(k,n):
    return (n/k)*math.gamma(1+1/k)**n/math.gamma(1+n/k)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--k',type=int,default=6)
    ap.add_argument('--n',type=int,default=5)
    ap.add_argument('--pmax',type=int,default=2000)
    ap.add_argument('--output')
    a=ap.parse_args()
    k,n=a.k,a.n
    rho=real_density(k,n)
    report={'k':k,'n':n,'rho':rho,'primes':{}}
    prod=1.0
    prev_p=0
    for p in primes_upto(a.pmax):
        if k%p==0:
            seq=sigma_p(p,k,n,12)
            # take the last two and require agreement to 1e-9 relative
            vals=[s[2] for s in seq]
            stable = len(vals)>=2 and abs(vals[-1]-vals[-2])<=1e-9*max(vals[-1],1e-300)
            sig=vals[-1]
            report['primes'][p]={'sequence':[(j,str(N),v) for j,N,v in seq],'sigma':sig,'stable':stable}
        else:
            sig=sigma_smooth_fft(p,k,n)
            report['primes'][p]={'sigma':sig}
        prod*=sig
        if p in (1000,) or (p>1000 and prev_p<=1000) or (p>5000 and prev_p<=5000):
            report.setdefault('partial_products',{})[str(prev_p)]=prod/sig
        prev_p=p
    # tail bound: |sigma_p - 1| <= (k-1)^(n+1) p^{-(n-1)/2} for p not | k (Weil),
    # sum over p > pmax integrated
    P=a.pmax
    # Weil-type bound |sigma_p-1| <= (k-1)^(n+1) p^{-(n-1)/2} for p not | k,
    # summed over p>P by the integral  (k-1)^(n+1) * 2/(n-3) * P^{-(n-3)/2}/log P.
    # Weil bound for the diagonal hypersurface: for p not dividing k,
    # |N*(p) - p^n| <= (p-1) (k-1)^(n+1)/k * p^((n-1)/2), hence
    # |sigma_p - 1| <= (k-1)^(n+1)/k * p^((1-n)/2).  Summing p > P with the
    # prime number theorem gives the log-bound below on the tail product.
    tail=(k-1)**(n+1)/k * 2/(n-3) * P**(-(n-3)/2)/math.log(P)
    C=rho*prod/math.factorial(n)
    e=n-k+1   # growth exponent: N(F) ~ C F^e/e  (e>0)  or  C log F  (e==0)
    def expected(F): return C*math.log(F) if e==0 else C*F**e/e
    report.update({'product_sigma_p_le_pmax':prod,'C_primitive_unordered':C,
                   'growth_exponent':e,'tail_log_bound':tail,
                   'expected_primitive_solutions':{str(F):expected(F) for F in
                       [144,1141,85359,730000,4300000,10**7,10**8,10**12,10**20,10**100]}})
    js=json.dumps(report,indent=1)
    print(js if not a.output else f"k={k} n={n} C={C:.6e}  prod={prod:.6e}  growth=F^{e}  tail<= x{math.exp(tail):.4f}")
    if a.output:
        open(a.output,'w').write(js)

if __name__=='__main__': main()
