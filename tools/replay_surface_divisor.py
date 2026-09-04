#!/usr/bin/env python3
"""Independent replay of the retained sparse surface search.

Python standard library only; does not import the search, SymPy, or FLINT.
Proves all used factors prime from Lucas certificates, reconstructs the
specified CRT domain, and independently checks each rejection and each
divisor in the complete representation window.
"""
import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('directory',type=Path)
    directory=parser.parse_args().directory
    certs=json.loads((directory/'prime_certificates.json').read_text())
    proved=set()
    def prime(n):
        if n in proved:
            return
        c=certs[str(n)]
        if c['method']=='trial_division':
            assert 2<=n<=1000 and all(n%d for d in range(2,math.isqrt(n)+1))
        else:
            assert c['method']=='lucas'
            factors={int(p):e for p,e in c['n_minus_one'].items()}
            assert all(e>=1 for e in factors.values())
            assert math.prod(p**e for p,e in factors.items())==n-1
            for p in factors:
                assert p<n
                prime(p)
                a=int(c['witnesses'][str(p)])
                assert pow(a,n-1,n)==1
                assert math.gcd(pow(a,(n-1)//p,n)-1,n)==1
        proved.add(n)
    for n in certs:
        prime(int(n))
    domain=json.loads((directory/'divisor_domain.json').read_text())
    modulus=2*42**6
    assert domain['modulus']==modulus
    roots=domain['roots']
    assert len(roots)==len(set(roots))==144
    assert all(0<r<modulus and pow(r,6,modulus)==1 for r in roots)
    # Independently count all local roots. CRT gives exactly their product.
    assert math.prod(sum(pow(r,6,m)==1 for r in range(m))
                     for m in (128,729,117649))==144
    pairs=set()
    for f in domain['f_values']:
        assert math.gcd(f,42)==1
        for r in roots:
            t=r*f%modulus
            if 0<t<f:
                pairs.add((f,t))
                pairs.add((f,t+modulus*((f-1-t)//modulus)))
    pairs=sorted(pairs)
    digest=hashlib.sha256(''.join(f'{f},{t}\n' for f,t in pairs).encode()).hexdigest()
    assert len(pairs)==domain['candidate_count'] and digest==domain['candidate_sha256']
    counts=Counter()
    solutions=[]
    with (directory/'divisor_targets.jsonl').open() as inp:
        rows=[json.loads(line) for line in inp]
    assert len(rows)==len(pairs), 'Incomplete search cannot certify the full domain'
    for row,(f,t) in zip(rows,pairs):
        assert (int(row['f']),int(row['t']))==(f,t)
        assert (f**6-t**6)%modulus==0
        n=(f**6-t**6)//modulus
        assert int(row['target'])==n>0
        stage=row['stage']
        if stage=='valuation_reject':
            p=row['prime'];v=0;z=n
            while z%p==0:
                z//=p;v+=1
            assert p in (2,3,7)
            assert v%6 not in ({0,1} if p==2 else {0})
        elif stage=='residue_reject':
            m=row['modulus']
            sixths={pow(x,6,m) for x in range(m)}
            assert n%m not in {(x+y)%m for x in sixths for y in sixths}
        else:
            factors={int(p):e for p,e in row['factors'].items()}
            assert all(p in proved and e>=1 for p,e in factors.items())
            assert math.prod(p**e for p,e in factors.items())==n
            if stage=='inert_prime_reject':
                p=int(row['prime'])
                assert p%4==3 and factors[p]%6!=0
            else:
                assert stage=='divisor_test'
                # Different implementation: generate ALL divisors first,
                # then use cubic inequalities and the product x^2*y^2.
                divisors=[1]
                for p,e in factors.items():
                    divisors=[d*p**j for d in divisors for j in range(e+1)]
                tested=0;representations=[]
                for h in sorted(divisors):
                    if not n<h**3<=4*n:
                        continue
                    tested+=1
                    if (h*h-n//h)%3:
                        continue
                    product=(h*h-n//h)//3
                    disc=h*h-4*product
                    if product<=0 or disc<0:
                        continue
                    delta=math.isqrt(disc)
                    if delta*delta!=disc or (h+delta)%2:
                        continue
                    a,b=(h+delta)//2,(h-delta)//2
                    x,y=math.isqrt(a),math.isqrt(b)
                    if x*x==a and y*y==b:
                        assert x>0 and y>0 and x**6+y**6==n
                        representations.append([str(x),str(y)])
                assert tested==row['divisors_in_window']
                assert representations==row['representations']
                counts['divisors_in_window']+=tested
                solutions.extend(row.get('solution',[]))
        counts[stage]+=1
        counts['processed']+=1
    summary=json.loads((directory/'divisor_summary.json').read_text())
    assert dict(counts)==summary['counts']
    assert not solutions and summary['solution'] is None
    print(json.dumps({'result':'PASS','targets_replayed':len(rows),
                      'primes_proved':len(proved),'counts':dict(counts),
                      'candidate_sha256':digest,'solutions':0},indent=2))


if __name__=='__main__':
    main()
