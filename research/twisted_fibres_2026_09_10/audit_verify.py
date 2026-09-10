#!/usr/bin/env python3
"""Independent exact replay of the twisted-fibre ledger.

Uses standard-library integer and Fraction arithmetic. Does not recompute
elliptic ranks: it checks retained PARI transcripts and their ledger use.
Torsion point groups and good-reduction order bounds are checked separately.
The production module is imported only for differential synthetic controls.
"""
import argparse
from collections import Counter
from fractions import Fraction
import importlib.util
from itertools import combinations_with_replacement, product
import json
from math import comb, gcd, isqrt, lcm, prod
from pathlib import Path

HERE = Path(__file__).resolve().parent
KINDS = ('plus', 'minus', 'square', 'difference', 'minus_cubed', 'plus_cubed',
         'plus_four', 'minus_four', 'four_fourth')


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def prime64(n):
    """Deterministic Miller--Rabin for the explicitly enforced n<2^64."""
    require(0 < n < 2**64, f'prime certificate outside supported range: {n}')
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, s = n-1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 325, 9375, 28178, 450775, 9780504, 1795265022):
        if a % n == 0:
            continue
        z = pow(a, d, n)
        if z in (1, n-1):
            continue
        for _ in range(s-1):
            z = z*z % n
            if z == n-1:
                break
        else:
            return False
    return True


def root_integer(n, degree):
    if n < 0:
        if degree % 2 == 0:
            return None
        r = root_integer(-n, degree)
        return None if r is None else -r
    if degree == 2:
        r = isqrt(n)
        return r if r*r == n else None
    # Newton iteration, independently of production's binary search.
    if n == 0:
        return 0
    r = 1 << ((n.bit_length()+degree-1)//degree)
    while True:
        nxt = ((degree-1)*r+n//r**(degree-1))//degree
        if nxt >= r:
            break
        r = nxt
    return r if r**degree == n else None


def root_rational(value, degree):
    value = Fraction(value)
    a = root_integer(value.numerator, degree)
    b = root_integer(value.denominator, degree)
    return None if a is None or b is None else Fraction(a, b)


def coefficient(c, kind):
    if kind == 'plus': return c
    if kind == 'minus': return -c
    if kind == 'square': return c*c
    if kind == 'difference': return -432*c*c
    if kind == 'minus_cubed': return -c*c*c
    if kind == 'plus_cubed': return c*c*c
    if kind == 'plus_four': return 4*c
    if kind == 'minus_four': return -4*c
    if kind == 'four_fourth': return 4*c**4
    raise ValueError(kind)


def point(raw):
    if raw is None or raw in ([0], ['0']):
        return None
    require(len(raw) == 2, f'invalid point representation: {raw}')
    return tuple(map(Fraction, raw))


def lift(c, kind, raw, scale=1):
    """Independent exact inverse; None means this point has no positive lift."""
    p = point(raw)
    if p is None:
        return None
    x, y = p[0]*scale**2, p[1]*scale**3
    require(y*y == x*x*x+coefficient(c, kind), f'point equation failed: {kind}')
    if kind == 'difference':
        if x == 0:
            return None
        # First recover ratios to B, unlike production's u,v square tests.
        u = (y+36*c)/(6*x)
        v = (y-36*c)/(6*x)
        if u == 0:
            return None
        A, S, B = root_rational(v/u, 2), root_rational(1/u, 2), Fraction(1)
    elif kind == 'plus':
        A, S, B = root_rational(x, 2), Fraction(1), root_rational(y, 3)
    elif kind == 'minus':
        A, S, B = root_rational(y, 3), Fraction(1), root_rational(x, 2)
    elif kind == 'square':
        A, S, B = Fraction(1), root_rational(x/c, 2), root_rational(y/c, 3)
    elif kind == 'minus_cubed':
        A, S, B = Fraction(1), root_rational(y/c**2, 3), root_rational(x/c, 2)
    elif kind == 'plus_cubed':
        A, S, B = root_rational(-x/c, 2), root_rational(y/c**2, 3), Fraction(1)
    elif kind in ('plus_four', 'minus_four', 'four_fourth'):
        r = root_rational(x, 2)
        if r is None or r <= 0:
            return None
        if kind in ('plus_four', 'four_fourth') and y <= 0:
            return None
        # Test the reconstructed sixth powers directly, independently of
        # production's cube-root-then-square-root implementation.
        if kind == 'plus_four':
            t = (y+r**3)/2
            A, S, B = root_rational(t*t-c, 6), Fraction(1), root_rational(t, 3)
        elif kind == 'minus_four':
            t = (r**3-y)/2
            A, S, B = root_rational(t, 3), Fraction(1), root_rational(t*t+c, 6)
        else:
            t = (y+r**3)/(2*c*c)
            A, S, B = Fraction(1), root_rational((t*t-1)/c, 6), root_rational(t, 3)
    else:
        raise ValueError(kind)
    if any(z is None or z == 0 for z in (A, S, B)):
        return None
    A, S, B = map(abs, (A, S, B))
    require(A**6+c*S**6 == B**6, f'inverse identity failed: {kind}')
    return A/S, Fraction(1), B/S


def group_add(p, q):
    if p is None: return q
    if q is None: return p
    x1, y1 = p
    x2, y2 = q
    if x1 == x2 and y1 == -y2:
        return None
    slope = (3*x1*x1/(2*y1)) if p == q else (y2-y1)/(x2-x1)
    x3 = slope*slope-x1-x2
    return x3, slope*(x1-x3)-y1


def multiple(p, n):
    result = None
    while n:
        if n & 1:
            result = group_add(result, p)
        p = group_add(p, p)
        n //= 2
    return result


def torsion_reduction_bound(b, expected):
    g, used = 0, []
    for p in range(5, 200):
        if not prime64(p) or b % p == 0:
            continue
        order = 1
        for x in range(p):
            v = (x*x*x+b) % p
            if v == 0:
                order += 1
            elif pow(v, (p-1)//2, p) == 1:
                order += 2
        g = gcd(g, order)
        used.append((p, order))
        require(g % expected == 0, f'torsion reduction contradiction for b={b}')
        if g == expected:
            return True, used
    return False, used


def verify_model(record, filename_b):
    require(record['b'] == filename_b, 'model filename mismatch')
    if record['status'] != 'certified':
        require(record['status'] in ('timeout', 'error'), 'unknown model status')
        return dict(status=record['status'])
    b = filename_b
    lines = record['stdout'].splitlines()
    require(lines[-1] == 'COMPLETE' and 'VERSION:[2, 17, 4]' in lines,
            'missing retained GP completion/version')
    require(record['exit_code'] == 0 and '***' not in record['stderr'], 'GP error transcript')
    fields = {}
    for line in lines:
        if ':' in line:
            key, value = line.split(':', 1)
            if key != 'FIELD':
                require(key not in fields, f'duplicate GP field: {key}')
                fields[key] = value
    ranks = json.loads(fields['RANK'])
    tors = json.loads(fields['TORSION'])
    torspoints = json.loads(fields['TORSION_POINTS'])
    certs = json.loads(fields['CERTS'])
    require(certs == record['bnf_certificates'] and all(int(v) == 1 for v in certs),
            'class-group certificate mismatch')
    polynomials = [line[6:] for line in lines if line.startswith('FIELD:')]
    require(polynomials == record['bnf_polynomials'] and len(certs) == len(polynomials),
            'class-group field/certificate count mismatch')
    require([int(v) for v in ranks[:3]] == [record['rank_lower'], record['rank_upper'],
                record['sha2_lower_dimension']], 'rank transcript mismatch')
    require(0 <= record['rank_lower'] <= record['rank_upper'], 'impossible rank interval')
    require(ranks[3] == record['known_independent_points'], 'point transcript mismatch')
    require(tors == record['torsion'] and torspoints == record['torsion_points'],
            'torsion transcript mismatch')
    for raw in ranks[3]+torspoints:
        p = point(raw)
        require(p is None or p[1]**2 == p[0]**3+b, f'bad model point for {b}')
    order, invariants, generators = int(tors[0]), list(map(int, tors[1])), list(map(point, tors[2]))
    require(prod(invariants) == order and len(invariants) == len(generators), 'torsion invariants')
    listed = set(map(point, torspoints))
    require(len(listed) == len(torspoints) == order, 'duplicate/missing torsion points')
    generated = {None}
    for n, generator in zip(invariants, generators):
        require(multiple(generator, n) is None, 'torsion generator order does not divide invariant')
        cyclic = [multiple(generator, k) for k in range(n)]
        require(len(set(cyclic)) == n, 'torsion generator has smaller order')
        generated = {group_add(p, q) for p in generated for q in cyclic}
    require(generated == listed, 'torsion generators and point list disagree')
    full_bound, used = torsion_reduction_bound(b, order)
    method = 'reduction_orders' if full_bound else 'pari_elltors_only'
    if not full_bound:
        reduction_bound = gcd(*(n for _,n in used))
        # For y²=x³+b, psi_3=3x(x³+4b). A rational root of the
        # monic cubic x³+4b is integral. Check both x branches and
        # the y equation exactly. In the retained unresolved-by-order
        # cases, gcd=3 and neither branch gives rational 3-torsion.
        if reduction_bound == 3 and order == 1:
            other_x = root_integer(-4*b,3)
            has_three_torsion = (root_integer(b,2) is not None or
                (other_x is not None and root_integer(-3*b,2) is not None))
            require(not has_three_torsion, f'claimed trivial torsion misses rational 3-torsion: {b}')
            full_bound = True
            method = 'reduction_orders_and_3_division_polynomial'
    return dict(status='certified', independent_torsion_bound=full_bound,
                torsion_completeness_method=method,
                reduction_primes=len(used), exact_points=len(ranks[3])+len(torspoints))


def synthetic_controls():
    spec = importlib.util.spec_from_file_location('production_fibres_audit', HERE/'fibres.py')
    production = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(production)
    checked = 0
    for h in (1, 2, 3, 7):
        c = 63*h**6
        A, S, B = Fraction(1), Fraction(1, h), Fraction(2)
        u, v = (B/S)**2, (A/S)**2
        source = dict(plus=((A/S)**2, (B/S)**3), minus=((B/S)**2, (A/S)**3),
            square=(c*(S/A)**2, c*(B/A)**3),
            minus_cubed=(c*(B/A)**2, c*c*(S/A)**3),
            plus_cubed=(-c*(A/B)**2, c*c*(S/B)**3),
            plus_four=((A*A/(B*S))**2, (B**6+c*S**6)/(B**3*S**3)),
            minus_four=((B*B/(A*S))**2, (c*S**6-A**6)/(A**3*S**3)),
            four_fourth=((c*S*S/(A*B))**2, c*c*(B**6+A**6)/(A**3*B**3)),
            difference=(12*c/(u-v), 36*c*(u+v)/(u-v)))
        models = production.models(c, production.factor(c))
        for kind, (x, y) in source.items():
            scale = models[kind]['scale']
            for sign in (1, -1):
                normalized = (x/scale**2, sign*y/scale**3)
                independent = lift(c, kind, normalized, scale)
                actual, _ = production.decode(c, kind, normalized, scale)
                if kind in ('difference','plus_four','minus_four','four_fourth') and sign == -1:
                    require(independent is None and actual is None, 'negative difference branch')
                else:
                    require(independent == (A/S, Fraction(1), B/S), 'synthetic inverse failed')
                    require(actual and (actual[0]/actual[1], Fraction(1), actual[2]/actual[1]) == independent,
                            'production/independent inverse mismatch')
                checked += 1
    # Retained Y seeds are exact cube differences but fail the square lift.
    yfile = HERE.parent/'cube_ansatz_2026_09_10/y_points_structured_N360.json'
    seeds = json.loads(yfile.read_text())
    for seed in seeds:
        c = sum(a**6 for a in seed['x'])
        u, v = Fraction(seed['x6']**2), Fraction(seed['T'])
        require(u**3-v**3 == c, 'retained Y-seed equation failed')
        p = (12*c/(u-v), 36*c*(u+v)/(u-v))
        require(lift(c, 'difference', p) is None, 'retained Y seed unexpectedly lifts')
        require(production.decode(c, 'difference', p)[0] is None, 'production accepted retained Y seed')
    # Exact quadratic-field arithmetic: both algebraic seeds remain controls only.
    def mul_quad(z, w, d):
        return z[0]*w[0]+d*z[1]*w[1], z[0]*w[1]+z[1]*w[0]
    for a, d, tail, target in ((7, -121, (8,12,15), 17), (9, -249, (14,18,0), 22)):
        z = (1, 0)
        for _ in range(6): z = mul_quad(z, (a,1), d)
        require(2*z[0]+sum(t**6 for t in tail) == target**6, 'algebraic seed identity failed')
        require(root_integer(d, 2) is None, 'algebraic seed accidentally rational')
    return dict(synthetic_map_sign_scale_checks=checked, retained_Y_seeds=len(seeds),
                quadratic_field_seed_identities=2, ESOP6_positive_fixture=False,
                jacobian_differential_checks=differential_controls())


def differential_controls():
    """Formal Laurent-polynomial map/differential checks modulo B^6=A^6+c."""
    def mon(k=1,a=0,b=0,e=0): return {(a,b,e):Fraction(k)}
    def add(*ps):
        out = {}
        for p in ps:
            for m,v in p.items(): out[m] = out.get(m,0)+v
        return {m:v for m,v in out.items() if v}
    def mul(p,q):
        out = {}
        for m,v in p.items():
            for n,w in q.items():
                k = tuple(a+b for a,b in zip(m,n))
                out[k] = out.get(k,0)+v*w
        return {m:v for m,v in out.items() if v}
    def powp(p,n):
        r=mon()
        for _ in range(n): r=mul(r,p)
        return r
    def neg(p): return {m:-v for m,v in p.items()}
    def deriv(p):
        out=[]
        for (a,b,e),v in p.items():
            if a: out.append(mon(a*v,a-1,b,e))
            if b: out.append(mon(b*v,a+5,b-6,e))
        return add(*out)
    def zero_on_curve(p):
        if not p: return True
        shift=[-min(0,min(m[i] for m in p)) for i in range(3)]
        out={}
        for (a,b,e),v in p.items():
            a,b,e=a+shift[0],b+shift[1],e+shift[2]
            q,r=divmod(b,6)
            for k in range(q+1):
                m=(a+6*k,r,e+q-k)
                out[m]=out.get(m,0)+v*comb(q,k)
        return not any(out.values())
    # name, x, y, Weierstrass constant, differential scalar, (a,b)
    rows=[
      ('plus',mon(a=2),mon(b=3),mon(e=1),2,(2,3)),
      ('minus',mon(b=2),mon(a=3),mon(-1,e=1),2,(3,2)),
      ('square',mon(a=-2,e=1),mon(a=-3,b=3,e=1),mon(e=2),-2,(1,3)),
      ('minus_cubed',mon(a=-2,b=2,e=1),mon(a=-3,e=2),mon(-1,e=3),-2,(1,2)),
      ('plus_cubed',mon(-1,a=2,b=-2,e=1),mon(b=-3,e=2),mon(e=3),-2,(2,1)),
      ('plus_four',mon(a=4,b=-2),add(mon(b=3),mon(b=-3,e=1)),mon(4,e=1),2,(4,1)),
      ('minus_four',mon(a=-2,b=4),add(mon(a=-3,e=1),mon(-1,a=3)),mon(-4,e=1),-2,(1,4)),
      ('four_fourth',mon(a=-2,b=-2,e=2),add(mon(a=-3,b=3,e=2),mon(a=3,b=-3,e=2)),mon(4,e=4),-2,(1,1)),
      ('square_second',mon(-1,b=-2,e=1),mon(a=3,b=-3,e=1),mon(e=2),2,(3,1)),
    ]
    pairs=set()
    for name,x,y,k,scalar,(a,b) in rows:
        require(zero_on_curve(add(powp(y,2),neg(powp(x,3)),neg(k))),f'formal equation {name}')
        target=mon(scalar,a-1,b-6)
        require(zero_on_curve(add(deriv(x),neg(mul(target,y)))),f'formal differential {name}')
        pairs.add((a,b))
    # The cubic-difference x has a nonmonomial denominator D=B²−A².
    D=add(mon(b=2),mon(-1,a=2))
    N=add(mon(b=2),mon(a=2))
    target=mon(Fraction(2,3),a=1,b=-4)
    # dx/y = -D'/(3DN), since x=12c/D and y=36cN/D.
    require(zero_on_curve(add(neg(deriv(D)),neg(mul(mon(3),mul(target,mul(D,N)))))),
            'formal differential difference')
    pairs.add((2,2))
    expected={(a,b) for a in range(1,5) for b in range(1,5) if a+b<=5}
    require(pairs==expected and len(pairs)==10, 'differentials fail to span genus-ten basis')
    return dict(monomial_quotient_equations=9, pullback_differentials=10, independent_basis=True)


def verify_ledger(path):
    data = json.loads(path.read_text())
    fibres = data['fibres']
    require(len({f['c'] for f in fibres}) == len(fibres), 'duplicate normalized fibres')
    summary = Counter()
    cache = {}
    model_audits = {}
    for p in sorted((path.parent/'models').glob('*.json')):
        record = json.loads(p.read_text())
        b = int(p.stem)
        cache[b] = record
        model_audits[b] = verify_model(record, b)
        summary['model_'+record['status']] += 1
        if record['status'] == 'certified':
            summary['model_torsion_independently_bounded'] += model_audits[b]['independent_torsion_bound']
            summary['model_torsion_'+model_audits[b]['torsion_completeness_method']] += 1
            summary['exact_model_points_checked'] += model_audits[b]['exact_points']
    box_seen = set()
    rep_seen = set()
    models_seen = set()
    dispositions = Counter()
    rejection_models = set()
    for f in fibres:
        c = f['c']
        fac = {int(p):e for p,e in f['factorization'].items()}
        require(c > 0 and prod(p**e for p,e in fac.items()) == c, 'coefficient factorization')
        require(all(prime64(p) and 1 <= e <= 5 for p,e in fac.items()), 'noncanonical sixth class')
        require(f['cube_class']*f['cube_scale']**3 == c, 'cube-class scale identity')
        require(all(f['cube_class'] % p**3 != 0 for p in fac), 'cube class not cube-free')
        require(f['genus'] == 10 and f['local_data']['all_places_nonzero_soluble'] is True,
                'genus/local-data contradiction')
        for representation in f['representations']:
            aa = representation['original']
            norm = tuple(map(Fraction, representation['normalized']))
            h = representation['sixth_scale']
            require(len(aa) == 4 and all(a > 0 for a in aa) and aa == sorted(aa), 'original tuple shape')
            require(gcd(*aa) == 1 and h > 0, 'primitive tuple/scale')
            require(sum(a**6 for a in aa) == representation['original_c'] == c*h**6,
                    'raw sixth-scale identity')
            require(norm == tuple(Fraction(a,h) for a in aa) and sum(a**6 for a in norm) == c,
                    'rational representation identity')
            require((c,norm) not in rep_seen, 'duplicate stored representation')
            rep_seen.add((c,norm))
            if representation['source']['source'] == 'box':
                require(representation['source']['height'] == data['height'], 'box height provenance')
                box_seen.add(tuple(aa))
            for p,m in ((2,8),(3,9),(7,7)):
                count = sum(a % p != 0 for a in aa)
                local = f['local_data']['charts'][str(p)]
                require(1 <= count <= 4 and count == c % m == local['unit_count'] == local['c_residue'],
                        'local unit count')
                require(local['modulus'] == m and local['B'] == 'unit', 'local metadata')
                expected = ['S_divisible_A_unit']+(['A_divisible_S_unit'] if count == 1 else [])
                require(local['charts'] == expected, 'local valuation chart list')
        require(set(f['models']) == set(KINDS), 'missing elliptic quotient')
        for kind, model in f['models'].items():
            original_b = coefficient(c, kind)
            reduced_b, scale = original_b, 1
            # Normalize numerically, independent of production's exponent formulas.
            for p in set(fac) | {2,3}:
                while reduced_b % p**6 == 0:
                    reduced_b //= p**6
                    scale *= p
            require(model == dict(b=reduced_b, scale=scale, source_b=original_b), 'model normalization')
            models_seen.add(reduced_b)
        for kind, entry in f['elliptic_data'].items():
            model = f['models'][kind]
            b = model['b']
            require(entry['model_b'] == b and entry['model_cache_file'] == f'models/{b}.json', 'model link')
            require(b in cache, 'missing retained arithmetic model')
            r = cache[b]
            require(entry['status'] == r['status'], 'ledger/cache arithmetic status')
            if r['status'] == 'certified':
                require(entry['rank_lower'] == r['rank_lower'] and entry['rank_upper'] == r['rank_upper']
                        and entry['torsion_order'] == int(r['torsion'][0]), 'ledger/cache rank or torsion')
                for raw in r['torsion_points']+r['known_independent_points']:
                    p = point(raw)
                    signs = [raw] if p is None or p[1] == 0 else [raw, [p[0], -p[1]]]
                    for signed in signs:
                        require(lift(c,kind,signed,model['scale']) is None,
                                f'UNPROMOTED EXACT POSITIVE FIBRE POINT: c={c}, kind={kind}')
                        summary['exact_lift_attempts_replayed'] += 1
        rejection = f['rejection']
        if rejection is None:
            require(f['square_cover_status'] == 'unresolved', 'unresolved status mismatch')
            dispositions['unresolved'] += 1
        else:
            kind = rejection['kind']
            b = f['models'][kind]['b']
            r = cache[b]
            require(f['square_cover_status'] == 'certified_no_positive_point', 'rejection status mismatch')
            require(rejection['reason'] == 'certified_rank_zero_all_torsion_failed_lifts', 'rejection reason')
            require(rejection['model_b'] == b and r['status'] == 'certified' and r['rank_upper'] == 0,
                    'rejection lacks saved certified rank-zero upper bound')
            require(rejection['torsion_points'] == r['torsion_points'], 'rejection torsion list mismatch')
            require(kind in f['elliptic_data'], 'rejection missing elliptic data')
            dispositions[kind] += 1
            summary['certified_rejections_replayed'] += 1
            summary['rejecting_fibres_torsion_independently_complete'] += model_audits[b]['independent_torsion_bound']
            rejection_models.add(b)
    all_box = set(combinations_with_replacement(range(1,data['height']+1),4))
    expected_box = {aa for aa in all_box if gcd(*aa) == 1}
    require(box_seen == expected_box, 'bounded tuple coverage mismatch')
    counts = data['counts']
    expected_counts = dict(sorted_tuples=len(all_box), nonprimitive_tuples_removed=len(all_box-expected_box),
        distinct_sixth_power_classes=len(fibres), distinct_cube_classes=len({f['cube_class'] for f in fibres}),
        distinct_mordell_models=len(models_seen), stored_representations=len(rep_seen))
    require(counts == expected_counts, 'enumeration counts mismatch')
    # During a live run the optional stage summary can lag current entries.
    summary.update(fibres=len(fibres), representations=len(rep_seen), box_rays=len(box_seen))
    summary['distinct_rejecting_models'] = len(rejection_models)
    summary['rejecting_models_torsion_independently_complete'] = sum(
        model_audits[b]['independent_torsion_bound'] for b in rejection_models)
    return dict(result='PASS', scope='exact independent replay; saved PARI rank bounds not recomputed',
                checks=dict(summary), dispositions=dict(dispositions),
                synthetic_controls=synthetic_controls(), rank_targets=rank_targets(fibres, cache))


def rank_targets(fibres, cache):
    """Eligibility only: no p-adic heights or rational-point classification."""
    pairs = {
        't6_plus_c': ('plus','square'),
        't6_minus_c': ('minus','square'),
        't6_plus_4c': ('plus_four','difference'),
        't6_minus_4c': ('minus_four','difference'),
        't6_plus_c2_over4': ('difference','four_fourth'),
    }
    pair_counts=Counter()
    complete=[]
    eligible=[]
    genus2_fibres=[]
    for f in fibres:
        if f['rejection'] is not None:
            continue
        ranks={}
        for kind, spec in f['models'].items():
            r=cache.get(spec['b'])
            if r and r['status']=='certified':
                ranks[kind]=(r['rank_lower'],r['rank_upper'])
        available=[]
        for label,(k1,k2) in pairs.items():
            if ranks.get(k1)==(1,1) and ranks.get(k2)==(1,1):
                pair_counts[label]+=1
                available.append(label)
        if available:
            genus2_fibres.append(dict(c=f['c'],rank_one_pairs=available))
        if set(ranks)==set(KINDS):
            lower=sum(v[0] for v in ranks.values())+ranks['square'][0]
            upper=sum(v[1] for v in ranks.values())+ranks['square'][1]
            item=dict(c=f['c'],jacobian_rank_lower=lower,jacobian_rank_upper=upper)
            complete.append(item)
            if upper<=19:
                eligible.append(item)
    return dict(rational_neron_severi_lower_bound=11,
        all_nine_rank_bounds_available=len(complete),
        quadratic_chabauty_upper_bound_at_most_19=len(eligible),
        genus2_rank_one_pair_counts=dict(pair_counts),
        distinct_fibres_with_rank_one_genus2_route=len(genus2_fibres),
        first_genus2_targets=sorted(genus2_fibres,key=lambda z:z['c'])[:10],
        first_genus10_targets=sorted(eligible,key=lambda z:z['c'])[:10],
        heights_integrals_and_rational_point_classification_executed=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--ledger',type=Path,default=HERE/'ledger.json')
    parser.add_argument('--selftest',action='store_true')
    args = parser.parse_args()
    result = synthetic_controls() if args.selftest else verify_ledger(args.ledger)
    print(json.dumps(result, indent=2, sort_keys=True))
