#!/usr/bin/env python3
"""Sparse, bounded, exact search for 2X^6+2Y^6+Z^6=W^6.

Primitive solutions have 42|X,Y, so W^6-Z^6 = 2*42^6*(x^6+y^6).
Generate specified CRT-root targets, factor the four cyclotomic factors of
W^6-Z^6, and solve the two-sixth-power problem by divisor/discriminant tests.
All factor primes receive recursive Lucas certificates. No Bloom filter,
floating point, generic height sweep, or probabilistic acceptance is used.

Requires python-flint for factorization, or falls back to SymPy. The output
can be replayed without either package by replay_surface_divisor.py.
"""
import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import random
import time

try:
    from flint import fmpz
    BACKEND = 'python-flint'
    def raw_factor(n):
        return {int(p):int(e) for p,e in fmpz(n).factor()}
except ImportError:
    from sympy import factorint
    BACKEND = 'sympy'
    def raw_factor(n):
        return {int(p):int(e) for p,e in factorint(n).items()}

MODULUS = 2*42**6
PAIR_MODS = (13,19,31,37,43)
PRIME_CERTS = {}


def factor_product(factors):
    return math.prod(p**e for p,e in factors.items())


def certify_prime(n):
    """Lucas converse to Fermat with the complete factorization of n-1."""
    if str(n) in PRIME_CERTS:
        return
    if n <= 1000:
        assert n >= 2 and all(n%d for d in range(2,math.isqrt(n)+1))
        PRIME_CERTS[str(n)] = {'method':'trial_division'}
        return
    fac = raw_factor(n-1)
    assert factor_product(fac) == n-1
    witnesses = {}
    for q in fac:
        certify_prime(q)
        for a in range(2,1025):
            if (pow(a,n-1,n) == 1 and
                math.gcd(pow(a,(n-1)//q,n)-1,n) == 1):
                witnesses[str(q)] = str(a)
                break
        else:
            raise ArithmeticError(f'No Lucas certificate for {n}, factor {q}')
    PRIME_CERTS[str(n)] = {'method':'lucas',
        'n_minus_one':{str(p):e for p,e in sorted(fac.items())},
        'witnesses':witnesses}


def certified_factor(n):
    fac = raw_factor(n)
    assert factor_product(fac) == n
    for p in fac:
        certify_prime(p)
    return fac


def iroot(n,k):
    if n < 0:
        raise ValueError('negative radicand')
    lo,hi=0,1 << ((n.bit_length()+k-1)//k)
    while lo < hi:
        mid=(lo+hi+1)//2
        if mid**k <= n:
            lo=mid
        else:
            hi=mid-1
    return lo


def pair_from_divisors(n,fac):
    """All unordered positive representations, with x>=y."""
    assert factor_product(fac) == n
    lo=iroot(n,3)+1  # h^3>n when both terms are positive
    hi=iroot(4*n,3)
    ds=[1]
    for p,e in sorted(fac.items()):
        old=ds
        ds=[]
        pe=1
        for _ in range(e+1):
            ds.extend(d*pe for d in old if d*pe <= hi)
            pe*=p
    pairs=[]
    tested=0
    for h in sorted(d for d in ds if d>=lo):
        tested += 1
        disc3=4*(n//h)-h*h
        if disc3<0 or disc3%3:
            continue
        disc=math.isqrt(disc3//3)
        if 3*disc*disc != disc3 or (h+disc)%2:
            continue
        x2,y2=(h+disc)//2,(h-disc)//2
        if y2<=0:
            continue
        x,y=math.isqrt(x2),math.isqrt(y2)
        if x*x==x2 and y*y==y2:
            assert x**6+y**6==n
            pairs.append([x,y])
    return pairs,tested


def monotone_pair(n):
    """Independent exact oracle for small controls, not divisor-based."""
    x,y=1,iroot(n,6)
    ans=[]
    while x<=y:
        total=x**6+y**6
        if total==n:
            ans.append([y,x]);x+=1;y-=1
        elif total<n:
            x+=1
        else:
            y-=1
    return sorted(ans)


def controls():
    rng=random.Random(660601)
    targets={x**6+y**6 for x in range(1,31) for y in range(1,x+1)}
    planted=len(targets)
    targets.update(rng.randrange(1,2*30**6+1) for _ in range(384))
    for n in sorted(targets):
        got,_=pair_from_divisors(n,certified_factor(n))
        assert sorted(got)==monotone_pair(n),(n,got,monotone_pair(n))
    return {'planted_targets':planted,'total_targets':len(targets),
            'independent_oracle':'monotone exact sixth-power search','result':'PASS'}


def roots():
    parts=(128,729,117649)
    local=[[r for r in range(1,m) if pow(r,6,m)==1] for m in parts]
    answer=set()
    for a in local[0]:
        for b in local[1]:
            for c in local[2]:
                r=sum(v*(MODULUS//m)*pow(MODULUS//m,-1,m)
                      for v,m in zip((a,b,c),parts))%MODULUS
                answer.add(r)
    assert len(answer)==144 and all(pow(r,6,MODULUS)==1 for r in answer)
    return sorted(answer)


def candidate_domain(per_band):
    rng=random.Random(660604)
    bands=[(4300001,10000000),(10000001,100000000),
           (100000001,1000000000),(1000000001,10000000000),
           (10000000001,100000000000),(100000000001,1000000000000)]
    fs=[]
    for lo,hi in bands:
        selected=set()
        while len(selected)<per_band:
            f=rng.randint(lo,hi)
            if math.gcd(f,42)==1:
                selected.add(f)
        fs.extend(sorted(selected))
    candidates=set()
    for f in fs:
        for r in roots():
            first=r*f%MODULUS
            if first>=f:
                continue
            last=first+((f-1-first)//MODULUS)*MODULUS
            for t in (first,last):
                assert 0<t<f and (pow(f,6,MODULUS)-pow(t,6,MODULUS))%MODULUS==0
                candidates.add((f,t))
    return bands,fs,sorted(candidates)


def valuation_reject(n):
    for p,allowed in ((2,{0,1}),(3,{0}),(7,{0})):
        v=0
        while n%p==0:
            n//=p;v+=1
        if v%6 not in allowed:
            return p
    return None


def factored_target(f,t):
    out=Counter()
    for n in (f-t,f+t,f*f+f*t+t*t,f*f-f*t+t*t):
        out.update(certified_factor(n))
    out.subtract({2:7,3:6,7:6})
    assert all(e>=0 for e in out.values())
    fac={p:e for p,e in out.items() if e}
    assert factor_product(fac)==(f**6-t**6)//MODULUS
    return fac


def run(args):
    start=time.monotonic()
    destination=Path(args.output)
    destination.mkdir(parents=True,exist_ok=True)
    control=controls()
    print(json.dumps({'controls':control,'backend':BACKEND}),flush=True)
    bands,fs,candidates=candidate_domain(args.per_band)
    digest=hashlib.sha256(''.join(f'{f},{t}\n' for f,t in candidates).encode()).hexdigest()
    domain={'seed':660604,'bands':bands,'f_values':fs,'per_band':args.per_band,
            'modulus':MODULUS,'roots':roots(),'selection':'first and last positive t<f on each root ray',
            'candidate_count':len(candidates),'candidate_sha256':digest}
    (destination/'divisor_domain.json').write_text(json.dumps(domain,indent=2)+'\n')
    masks={m:{(pow(x,6,m)+pow(y,6,m))%m for x in range(m) for y in range(m)}
           for m in PAIR_MODS}
    counts=Counter()
    solution=None
    with (destination/'divisor_targets.jsonl').open('w') as log:
        for index,(f,t) in enumerate(candidates):
            n=(f**6-t**6)//MODULUS
            row={'f':str(f),'t':str(t),'target':str(n)}
            bad=valuation_reject(n)
            if bad:
                row.update(stage='valuation_reject',prime=bad)
            else:
                bad=next((m for m in PAIR_MODS if n%m not in masks[m]),None)
                if bad:
                    row.update(stage='residue_reject',modulus=bad)
                else:
                    fac=factored_target(f,t)
                    row['factors']={str(p):e for p,e in sorted(fac.items())}
                    bad=next((p for p,e in sorted(fac.items()) if p%4==3 and e%6),None)
                    if bad:
                        row.update(stage='inert_prime_reject',prime=str(bad))
                    else:
                        pairs,tested=pair_from_divisors(n,fac)
                        row.update(stage='divisor_test',divisors_in_window=tested,
                                   representations=[[str(x),str(y)] for x,y in pairs])
                        counts['divisors_in_window']+=tested
                        if pairs:
                            x,y=pairs[0]
                            values=[42*x,42*x,42*y,42*y,t,f]
                            g=math.gcd(*values)
                            values=sorted(v//g for v in values[:5])+[f//g]
                            assert min(values)>0 and sum(v**6 for v in values[:5])==values[5]**6
                            solution=values
                            row['solution']=[str(v) for v in values]
            counts[row['stage']]+=1
            counts['processed']+=1
            log.write(json.dumps(row,sort_keys=True)+'\n')
            if solution:
                log.flush()
                print('SOLUTION '+json.dumps([str(v) for v in solution]),flush=True)
                break
            if (index+1)%1000==0:
                log.flush()
                print(json.dumps({'processed':index+1,'total':len(candidates),
                                  'seconds':round(time.monotonic()-start,3),
                                  'counts':dict(counts)}),flush=True)
    (destination/'prime_certificates.json').write_text(json.dumps(PRIME_CERTS,sort_keys=True,indent=2)+'\n')
    summary={'domain':domain,'controls':control,'backend':BACKEND,
             'counts':dict(counts),'solution':None if solution is None else [str(v) for v in solution],
             'seconds':time.monotonic()-start,'prime_certificates':len(PRIME_CERTS)}
    (destination/'divisor_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({'summary':summary},sort_keys=True),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--per-band',type=int,default=16)
    parser.add_argument('--output',default='results/astra_direct_2026_09_04')
    arguments=parser.parse_args()
    if not 1<=arguments.per_band<=128:
        parser.error('per-band must lie in [1,128]')
    run(arguments)
