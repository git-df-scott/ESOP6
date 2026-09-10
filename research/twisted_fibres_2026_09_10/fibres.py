#!/usr/bin/env python3
"""Exact fibre normalization, elliptic maps, and certificate-rooted search.

Standard library only. All arithmetic and lifting are exact. A finite point
search never closes a fibre; only a certified rank upper bound of zero plus
complete torsion checking does. Local conditions describe charts, not fibres.
"""
import argparse
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations_with_replacement
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
KINDS = ('minus_cubed', 'plus_cubed', 'plus', 'minus', 'square', 'difference')


def source_coefficient(c, kind):
    return {'plus':c, 'minus':-c, 'square':c*c, 'difference':-432*c*c,
            'minus_cubed':-c**3, 'plus_cubed':c**3}[kind]


def factor(n):
    """Complete trial division; deliberately bounded to the small pilot."""
    assert n > 0
    out = {}
    p = 2
    while p*p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0)+1
            n //= p
        p = 3 if p == 2 else p+2
    if n > 1:
        out[n] = out.get(n, 0)+1
    return out


def extract_power(fac, k):
    scale = math.prod(p**(e//k) for p, e in fac.items())
    residue = math.prod(p**(e % k) for p, e in fac.items())
    return residue, scale


def integer_root(n, k):
    if n < 0:
        if k % 2 == 0:
            return None
        t = integer_root(-n, k)
        return None if t is None else -t
    lo, hi = 0, 1 << ((n.bit_length()+k-1)//k)
    while lo <= hi:
        mid = (lo+hi)//2
        v = mid**k
        if v == n:
            return mid
        if v < n:
            lo = mid+1
        else:
            hi = mid-1
    return None


def rational_root(q, k):
    q = Q(q)
    a, b = integer_root(q.numerator, k), integer_root(q.denominator, k)
    return Q(a, b) if a is not None and b is not None else None


def models(c, fac):
    specs = {}
    for kind in KINDS:
        power = 3 if kind.endswith('_cubed') else (1 if kind in ('plus','minus') else 2)
        ff = {p:power*e for p,e in fac.items()}
        sign = -1 if kind in ('minus', 'difference', 'minus_cubed') else 1
        if kind == 'difference':
            ff[2] = ff.get(2, 0)+4
            ff[3] = ff.get(3, 0)+3
        b, h = extract_power(ff, 6)
        source_b = source_coefficient(c,kind)
        assert sign*b*h**6 == source_b
        specs[kind] = dict(b=sign*b, scale=h, source_b=source_b)
    return specs


def decode(c, kind, point, scale=1):
    """Return positive rational (A,S,B), or an exact failed lift reason."""
    if point is None or len(point) == 1:
        return None, 'infinity_boundary'
    x, y = Q(point[0])*scale**2, Q(point[1])*scale**3
    source_b = source_coefficient(c,kind)
    if y*y != x*x*x+source_b:
        raise ArithmeticError('invalid elliptic point')
    if kind == 'plus':
        A, B, S = rational_root(x, 2), rational_root(y, 3), Q(1)
    elif kind == 'minus':
        B, A, S = rational_root(x, 2), rational_root(y, 3), Q(1)
    elif kind == 'square':
        S, B, A = rational_root(x/c, 2), rational_root(y/c, 3), Q(1)
    elif kind == 'minus_cubed':
        B,S,A = rational_root(x/c,2), rational_root(y/c**2,3), Q(1)
    elif kind == 'plus_cubed':
        A,S,B = rational_root(-x/c,2), rational_root(y/c**2,3), Q(1)
    else:
        if x == 0:
            return None, 'zero_denominator'
        u, v = (36*c+y)/(6*x), (y-36*c)/(6*x)
        assert u**3-v**3 == c
        B, A, S = rational_root(u, 2), rational_root(v, 2), Q(1)
    if any(t is None for t in (A,S,B)):
        return None, 'power_class_failure'
    if not all((A,S,B)):
        return None, 'zero_coordinate_boundary'
    A,S,B = abs(A),abs(S),abs(B)
    assert A**6+c*S**6 == B**6
    return (A,S,B), 'exact_fibre_point'


def reconstruct(c, representation, point):
    aa = [Q(s) for s in representation]
    assert len(aa) == 4 and all(a > 0 for a in aa)
    assert sum(a**6 for a in aa) == c
    A,S,B = point
    rr = [a*S for a in aa]+[A,B]
    d = math.lcm(*(q.denominator for q in rr))
    vv = [int(q*d) for q in rr]
    g = math.gcd(*vv)
    vv = [v//g for v in vv]
    assert all(v > 0 for v in vv) and sum(v**6 for v in vv[:5]) == vv[5]**6
    return vv


def independent_verify(values):
    commands = [[sys.executable, str(ROOT/'tools/verify_esop6.py')],
                ['node', str(ROOT/'tools/verify_esop6.mjs')]]
    reports = []
    for cmd in commands:
        p = subprocess.run(cmd+list(map(str, values)), capture_output=True, text=True, check=True)
        r = json.loads(p.stdout)
        assert r['solution'] is True
        reports.append(r)
    assert reports[0] == reports[1]
    return reports


def local_data(c, original):
    # For a primitive fourtuple the unit counts are nonzero. Scaling by a
    # sixth power does not change these p=2,3,7 charts, since v_p(sum)<6.
    data = {}
    for p,m in ((2,8),(3,9),(7,7)):
        r = sum(a % p != 0 for a in original)
        assert 1 <= r <= 4 and c % m == r
        data[str(p)] = dict(unit_count=r, modulus=m, c_residue=c % m,
                            B='unit', charts=['S_divisible_A_unit']+
                            (['A_divisible_S_unit'] if r == 1 else []))
    return {'all_places_nonzero_soluble': True,
            'proof': 'smooth boundary; ARITHMETIC.md', 'charts':data}


def enumerate_fibres(height, include_y=True):
    ledger, counts = {}, Counter()
    tuples = []
    for aa in combinations_with_replacement(range(1, height+1),4):
        counts['sorted_tuples'] += 1
        if math.gcd(*aa) != 1:
            counts['nonprimitive_tuples_removed'] += 1
            continue
        tuples.append((aa, {'source':'box', 'height':height}))
    if include_y:
        p = ROOT/'research/cube_ansatz_2026_09_10/y_points_structured_N360.json'
        for s in json.loads(p.read_text()):
            aa = tuple(sorted(s['x']))
            assert sum(a**6 for a in aa)+s['T']**3 == s['x6']**6
            tuples.append((aa,dict(source=str(p.relative_to(ROOT)),u=s['x6']**2,v=s['T'])))
    for aa, source in tuples:
        raw_c = sum(a**6 for a in aa)
        raw_fac = factor(raw_c)
        c,h = extract_power(raw_fac,6)
        rep = [str(Q(a,h)) for a in aa]
        if c not in ledger:
            fac = {p:e % 6 for p,e in raw_fac.items() if e % 6}
            cube_n,cube_scale = extract_power(fac,3)
            ledger[c] = dict(c=c, factorization={str(p):e for p,e in fac.items()},
                cube_class=cube_n, cube_scale=cube_scale, representations=[],
                models=models(c,fac), local_data=local_data(c,aa),
                elliptic_data={}, square_cover_status='unresolved',
                genus=10, rejection=None)
        if rep not in [r['normalized'] for r in ledger[c]['representations']]:
            ledger[c]['representations'].append(dict(normalized=rep,original=list(aa),
                original_c=raw_c, sixth_scale=h, source=source))
    counts['distinct_sixth_power_classes'] = len(ledger)
    counts['distinct_cube_classes'] = len(set(x['cube_class'] for x in ledger.values()))
    counts['distinct_mordell_models'] = len(set(m['b'] for x in ledger.values() for m in x['models'].values()))
    counts['stored_representations'] = sum(len(x['representations']) for x in ledger.values())
    return ledger, dict(counts)


def bounded_cube_difference(c, denominator_bound, numerator_bound):
    """Exact optional oracle in a declared finite (u,v) rational box.

    This is not used as a rejection gate: absence in the box means nothing
    about global rank or existence. Shared denominator D is primitive.
    """
    for d in range(1,denominator_bound+1):
        for b in range(-numerator_bound,numerator_bound+1):
            a = integer_root(b**3+c*d**3,3)
            if a is not None and abs(a) <= numerator_bound and math.gcd(a,b,d) == 1:
                yield Q(a,d),Q(b,d)


def save_json(path, obj):
    path = Path(path)
    tmp = path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n')
    tmp.replace(path)


def selftest():
    # Generic diagonal-sextic positive control, not an ESOP6 certificate.
    c = 63
    points = {'plus':(1,8),'minus':(4,1),'square':(63,504),'difference':(252,3780),
              'minus_cubed':(252,3969),'plus_cubed':(Q(-63,4),Q(3969,8))}
    for kind,p in points.items():
        triple, why = decode(c,kind,p)
        A,S,B = triple
        assert (A/S,B/S) == (Q(1),Q(2)), (kind,triple,why)
        assert decode(c,kind,[0])[0] is None
    assert rational_root(Q(4,9),2) == Q(2,3)
    assert rational_root(Q(-8,27),3) == Q(-2,3)
    assert rational_root(Q(-4,9),2) is None
    assert rational_root(Q(17,16),2) is None
    assert rational_root(Q(10**120+1),2) is None
    assert extract_power(factor(4*7**6),6) == (4,7)
    for c in range(1,150):
        for m in models(c,factor(c)).values():
            assert m['b']*m['scale']**6 == m['source_b']
    # c=1 is rank-zero with only boundary cubic-difference points.
    for p in ((12,36),(12,-36)):
        assert decode(1,'difference',p)[0] is None
    # c=2: nonboundary cubic point u=1,v=-1 is not a square lift.
    assert decode(2,'difference',(12,0))[0] is None
    print(json.dumps({'result':'PASS','maps':6,'generic_fibre_positive_control':63,
                      'positive_ESOP6_fixture':False}))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--height',type=int,default=10)
    ap.add_argument('--output',default=str(HERE/'ledger.json'))
    ap.add_argument('--selftest',action='store_true')
    args = ap.parse_args()
    if args.selftest:
        selftest()
    else:
        start=time.time()
        ledger,counts=enumerate_fibres(args.height)
        out={'schema':1,'height':args.height,'counts':counts,'fibres':list(ledger.values()),
             'elapsed_seconds':time.time()-start,'scope':'all sorted primitive fourtuple rays in box, plus two retained Y seeds'}
        save_json(args.output,out)
        print(json.dumps(counts,sort_keys=True))
