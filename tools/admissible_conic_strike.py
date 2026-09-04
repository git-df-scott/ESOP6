#!/usr/bin/env python3
"""Search conics selected by the proved 2/3/7 valuation constraints.

Seed points below are checked by exact substitution; the four missing
conics in the seed-discovery record are NOT claimed independently empty.
First test projective square lifting on every point of P^1(Z/p^k Z).
If a conic passes all these necessary checks, search the declared rational
parameter box. No local or finite-search result is a global surface result.
"""
import argparse
from collections import Counter
import hashlib
import itertools
import json
import math
from pathlib import Path
import time

SEEDS=[(864,931,938448,607860,1084489),
       (864,1225,1543500,373248,1500625),
       (864,1519,109368,13428,115601),
       (864,1813,456876,257400,744485)]


def quotient(seed,u,v):
    p,q,a,b,d=seed
    plus=2*q*u*u+p*v*v
    A=a*(2*q*u*u-p*v*v)
    B=-b*plus+2*p*d*u*v
    D=d*plus-4*q*b*u*v
    return [q*(A+B),q*(A-B),q*D-p*A,q*D+p*A]


def projective_square_possible(values,p,k):
    """Necessary and sufficient for a nonzero residue vector after removing
    its known common p-power to be unit-scaled to squares at that precision.
    If all residues vanish, return True (unknown at this precision).
    """
    modulus=p**k
    values=[x%modulus for x in values]
    if not any(values):
        return True
    def val(x):
        if not x:
            return k
        v=0
        while x%p==0:
            x//=p;v+=1
        return v
    vs=[val(x) for x in values]
    common=min(vs)
    unit_conditions=[]
    for x,v in zip(values,vs):
        if not x:
            continue
        if (v-common)%2:
            return False
        unit=x//(p**v)
        if p==2:
            bits=min(3,k-v)
            unit_conditions.append((unit%(2**bits),bits))
        else:
            unit_conditions.append(pow(unit%p,(p-1)//2,p))
    if p!=2:
        return len(set(unit_conditions))<=1
    return all((a-b)%(2**min(ka,kb))==0
               for a,ka in unit_conditions for b,kb in unit_conditions)


def controls():
    # Exhaust all vectors modulo 8, 9, 7, and compare with direct scaling
    # after removing the common p-adic valuation. This is independent of
    # the parity/Legendre implementation above.
    tests=0
    for p,k in [(2,3),(3,2),(7,1)]:
        mod=p**k
        for vector in itertools.product(range(mod),repeat=4):
            if not any(vector):
                expected=True
            else:
                reduced=list(vector);precision=mod
                while all(x%p==0 for x in reduced):
                    reduced=[x//p for x in reduced];precision//=p
                squares={i*i%precision for i in range(precision)}
                expected=any(all((c*x)%precision in squares for x in reduced)
                             for c in range(1,precision) if c%p)
            assert projective_square_possible(vector,p,k)==expected,(vector,p,k)
            tests+=1
    return tests


def local_scan(seed,p,k):
    mod=p**k
    allowed=[]
    for u in range(mod):
        if projective_square_possible(quotient(seed,u,1),p,k):
            allowed.append((u,1))
    for v in range(0,mod,p):
        if projective_square_possible(quotient(seed,1,v),p,k):
            allowed.append((1,v))
    return {'prime':p,'exponent':k,'projective_parameters':mod+mod//p,
            'allowed':len(allowed),
            'allowed_sha256':hashlib.sha256(''.join(f'{u},{v}\n' for u,v in allowed).encode()).hexdigest()}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--height',type=int,default=2000)
    parser.add_argument('--output',type=Path,default=Path('results/astra_direct_2026_09_04'))
    args=parser.parse_args();start=time.monotonic()
    tested=controls();print(f'Local square-lift controls PASS: {tested}',flush=True)
    rows=[];solution=None
    for seed in SEEDS:
        p,q,a,b,d=seed
        assert (2*q**3-p**3)*a*a+6*q**3*b*b-3*p*q*q*d*d==0
        # Check representative secants including both charts and infinity.
        for u,v in [(1,0),(0,1),(1,1),(-3,2),(7,11)]:
            U,V,R,T=quotient(seed,u,v)
            assert 2*U**3+2*V**3+R**3==T**3
            assert q*(T-R)==p*(U+V)
        row={'lambda':[p,q],'base_point':[a,b,d], 'local_checks':[]}
        for pp,kk in [(2,8),(2,12),(3,5),(3,7),(7,3),(7,4),
                      (13,2),(19,2),(31,2)]:
            check=local_scan(seed,pp,kk)
            row['local_checks'].append(check)
            if check['allowed']==0:
                row['result']='NO LOCAL SQUARE LIFT'
                break
        if row.get('result'):
            print(json.dumps(row),flush=True);rows.append(row);continue
        counts=Counter();digest=hashlib.sha256()
        params=itertools.chain([(1,0)],((u,v) for v in range(1,args.height+1)
                 for u in range(-args.height,args.height+1) if math.gcd(u,v)==1))
        for u,v in params:
            counts['coprime_parameters']+=1
            Q=quotient(seed,u,v)
            if max(Q)<0:
                Q=[-x for x in Q]
            if min(Q)<=0:
                continue
            counts['positive_quotient_points']+=1
            digest.update(f'{u},{v}\n'.encode())
            g=math.gcd(*Q);primitive=[x//g for x in Q]
            rr=[math.isqrt(x) for x in primitive]
            if all(x*x==y for x,y in zip(rr,primitive)):
                x,y,z,w=rr
                assert 2*x**6+2*y**6+z**6==w**6
                solution=[x,x,y,y,z,w]
                counts['square_lifts']+=1
                print('SOLUTION '+json.dumps([str(x) for x in solution]),flush=True)
                break
        row.update(result='SOLUTION' if solution else 'NO LIFT IN PARAMETER BOX',
            domain={'u':[-args.height,args.height],'v':[1,args.height],
                    'gcd_uv':1,'also_infinity':[1,0]},counts=dict(counts),
            positive_parameter_sha256=digest.hexdigest())
        rows.append(row);print(json.dumps(row),flush=True)
        if solution:
            break
    output={'local_control_vectors':tested,'conics':rows,
            'solution':None if solution is None else [str(x) for x in solution],
            'seconds':time.monotonic()-start}
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/'admissible_conic_search.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({'seconds':output['seconds'],'solution':output['solution']}),flush=True)


if __name__=='__main__':
    main()
