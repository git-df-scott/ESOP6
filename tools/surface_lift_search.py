#!/usr/bin/env python3
"""Bounded square-lift search on an explicitly parametrized quotient conic.

The conic lies on 2U^3+2V^3+R^3=T^3 at
2(T-R)=U+V. A projective point lifts iff its positive primitive integer
coordinate vector consists of four squares. No floating point is used.
"""
import argparse
import hashlib
import json
import math
import time


def square_lift(values):
    if min(values) <= 0:
        return None
    g = math.gcd(*values)
    primitive = [v//g for v in values]
    roots = [math.isqrt(v) for v in primitive]
    return roots if all(r*r == v for r,v in zip(roots, primitive)) else None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--height', type=int, default=2000)
    args = parser.parse_args()
    assert square_lift([37,148,333,592]) == [1,2,3,4]
    assert square_lift([1,4,9,17]) is None
    assert square_lift([0,4,9,16]) is None
    start = time.monotonic()
    counts = {'coprime_parameters':0,'positive_quotient_points':0,
              'square_lifts':0}
    digest = hashlib.sha256()
    for v in range(1,args.height+1):
        for u in range(1,args.height+1):
            if math.gcd(u,v) != 1:
                continue
            counts['coprime_parameters'] += 1
            Q = [u*u+8*u*v-5*v*v,-u*u+8*u*v+5*v*v,
                 2*(u*u-2*u*v+5*v*v),2*(u*u+2*u*v+5*v*v)]
            if min(Q) <= 0:
                continue
            counts['positive_quotient_points'] += 1
            digest.update(f'{u},{v}\n'.encode())
            roots = square_lift(Q)
            if roots:
                x,y,z,w=roots
                assert 2*x**6+2*y**6+z**6 == w**6
                counts['square_lifts'] += 1
                print(json.dumps({'SOLUTION':[x,x,y,y,z,w],
                                  'parameter':[u,v]}),flush=True)
                return
    print(json.dumps({'domain':{'lambda':'1/2','u':[1,args.height],
                     'v':[1,args.height],'gcd_uv':1},
                     'counts':counts, 'positive_parameter_sha256':digest.hexdigest(),
                     'seconds':time.monotonic()-start},indent=2))


if __name__ == '__main__':
    main()
