#!/usr/bin/env python3
"""Finite exact subgroup hunt; bounded misses NEVER reject a fibre.

Consumes completed arithmetic caches without computing ranks again. The
generated subgroup is not claimed to be the full Mordell-Weil group.
Good-prime power tests discard coefficient tuples before Fraction arithmetic.
Every materialized elliptic point is checked exactly before the common lift,
reconstruction and two independent integer verification roots are called.
"""
import argparse
from collections import Counter
from fractions import Fraction as Q
import hashlib
from itertools import product
import json
from pathlib import Path
import sys
import time

from fibres import HERE, decode, independent_verify, rational_root, reconstruct, save_json

sys.set_int_max_str_digits(0)

PRIMES = (7,13,19,31,37,43,61,67,73,79,97,103,109,127,139,151)
SUPPORTED_KINDS = ('minus_cubed','plus_cubed','plus','minus','square','difference',
                   'plus_four','minus_four','four_fourth')
MOD_ROOTS = {}


def unpack(point):
    return None if point is None or len(point)==1 else tuple(map(Q,point))


def pack(point):
    return ['0'] if point is None else list(map(str,point))


def check(point,b):
    return point is None or point[1]**2 == point[0]**3+b


def add(P,R,p=None):
    if P is None:
        return R
    if R is None:
        return P
    x,y=P;u,v=R
    if p:
        if x==u and (y+v)%p==0:
            return None
        m=((3*x*x)*pow(2*y,-1,p) if x==u else (v-y)*pow(u-x,-1,p))%p
        z=(m*m-x-u)%p
        return z,(m*(x-z)-y)%p
    if x==u and y+v==0:
        return None
    m=3*x*x/(2*y) if x==u else (v-y)/(u-x)
    z=m*m-x-u
    return z,m*(x-z)-y


def mul(P,n,p=None):
    if n<0:
        P=None if P is None else (P[0],(-P[1])%p if p else -P[1])
        n=-n
    R=None
    while n:
        if n&1:
            R=add(R,P,p)
        P=add(P,P,p)
        n>>=1
    return R


def reduce_point(P,p):
    if P is None or any(v.denominator%p==0 for v in P):
        # At good reduction a rational point with a nonintegral coordinate
        # reduces to infinity. Infinity is conservatively not sieved out.
        return None
    return tuple(v.numerator*pow(v.denominator,-1,p)%p for v in P)


def possible_mod(c,kind,P,scale,p):
    assert kind in SUPPORTED_KINDS
    if P is None:
        return True
    x=P[0]*pow(scale,2,p)%p;y=P[1]*pow(scale,3,p)%p
    square=lambda z:z%p==0 or pow(z%p,(p-1)//2,p)==1
    cube=lambda z:z%p==0 or pow(z%p,(p-1)//3,p)==1
    if kind in ('plus','minus'):
        return square(x) and cube(y)
    if kind=='square':
        ci=pow(c,-1,p)
        return square(x*ci) and cube(y*ci)
    if kind in ('minus_cubed','plus_cubed'):
        ci=pow(c,-1,p)
        return square(x*ci*(-1 if kind=='plus_cubed' else 1)) and cube(y*ci*ci)
    if kind in ('plus_four','minus_four','four_fourth'):
        if p not in MOD_ROOTS:
            rr={2:{},3:{}}
            for k in (2,3):
                for z in range(p):
                    rr[k].setdefault(pow(z,k,p),[]).append(z)
            MOD_ROOTS[p]=rr
        rr=MOD_ROOTS[p];ci=pow(c,-1,p);half=pow(2,-1,p)
        for r in rr[2].get(x,[]):
            if kind=='minus_four':
                target=(r**3-y)*half%p
            else:
                target=(y+r**3)*half*(ci*ci if kind=='four_fourth' else 1)%p
            for t in rr[3].get(target,[]):
                if square(r*t*(ci if kind=='four_fourth' else 1)):
                    return True
        return False
    if x==0:
        return True
    inv=pow(6*x,-1,p)
    return square((36*c+y)*inv) and square((y-36*c)*inv)


def point_height(P):
    return max(max(abs(v.numerator).bit_length(),v.denominator.bit_length()) for v in P)


def torsion_fallback(b):
    """Exact torsion data only from cache, or the difference-model theorem.

    This function deliberately does not guess an empty torsion group. At
    good primes greater than two, all rational torsion injects into the
    finite group. A gcd bound dividing six plus the exact 2- and 3-division
    polynomials determines the entire group. The latter step matters here:
    difference curves are 3-isogenous to curves with rational 3-torsion, so
    finite-field group orders alone generally retain a spurious factor3.
    """
    orders=[]
    from math import gcd
    bound=0
    for p in (5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61):
        if b%p==0:
            continue
        n=1+sum(1 if (v:=(x*x*x+b)%p)==0 else 2 if pow(v,(p-1)//2,p)==1 else 0
                for x in range(p))
        orders.append(dict(prime=p,cardinality=n))
        bound=gcd(bound,n)
        if bound==1:
            return [None],dict(method='good_reduction_cardinality_gcd',orders=orders,bound=1)
        if bound in (2,3,6) and len(orders)>=3:
            generators=[];division=[]
            if bound%2==0:
                x=rational_root(-b,3)
                if x is not None:
                    generators.append((x,Q(0)))
                division.append(dict(order=2,polynomial='x^3+b',rational_x=None if x is None else str(x)))
            if bound%3==0:
                xx=[Q(0)];x=rational_root(-4*b,3)
                if x is not None:
                    xx.append(x)
                pp=[]
                for x in xx:
                    y=rational_root(x**3+b,2)
                    if y is not None:
                        pp.extend([(x,y),(x,-y)])
                generators+=pp
                division.append(dict(order=3,polynomial='3*x*(x^3+4*b)',
                                     rational_x_candidates=list(map(str,xx)),points=list(map(pack,pp))))
            tt={None}
            for P in generators:
                tt={add(T,mul(P,n)) for T in tt for n in range(bound)}
            assert all(check(P,b) for P in tt) and bound%len(tt)==0
            torsion=sorted(tt,key=lambda P:tuple(pack(P)))
            return torsion,dict(method='good_reduction_bound_and_exact_division_polynomials',
                               orders=orders,bound=bound,division_checks=division,
                               exact_torsion_order=len(torsion))
    return None,dict(method='good_reduction_cardinality_gcd',orders=orders,bound=bound,
                     status='insufficient_bound')


def cache_model(f,kind):
    spec=f['models'][kind]
    p=HERE/'models'/f'{spec["b"]}.json'
    if not p.exists():
        return None
    r=json.loads(p.read_text())
    return r if r.get('status')=='certified' else None


def choose_model(f):
    choices=[]
    for kind,spec in f['models'].items():
        if kind not in SUPPORTED_KINDS:
            continue
        r=cache_model(f,kind)
        if not r:
            continue
        points=[unpack(P) for P in r['known_independent_points']]
        if not points:
            continue
        points=sorted(points,key=point_height)[:3]
        # Coordinate growth tracks coefficient-square times generator height.
        estimate=max(map(point_height,points))*(144 if len(points)==1 else 12)
        choices.append((estimate,kind,r,points))
    return min(choices,key=lambda z:(z[0],z[1])) if choices else None


def hunt(f,kind,basis,torsion,provenance,rank1_bound=12,box_bound=2):
    c=f['c'];spec=f['models'][kind];b=spec['b'];scale=spec['scale']
    assert basis and all(check(P,b) for P in basis+torsion)
    rank=len(basis);bound=rank1_bound if rank==1 else box_bound
    primes=[p for p in PRIMES if (6*c*b*scale)%p][:8]
    tables={}
    for p in primes:
        tables[p]=([[mul(reduce_point(P,p),n,p) for n in range(-bound,bound+1)]
                    for P in basis],[reduce_point(T,p) for T in torsion])
    stats=Counter();rejected=Counter();outcomes=Counter();exact=[];seen=set()
    exact_multiples={};certificates=[]
    start=time.time()
    for coeffs in product(range(-bound,bound+1),repeat=rank):
        for ti,T in enumerate(torsion):
            stats['coefficient_torsion_tuples']+=1
            failed=None
            for p in primes:
                pts,tt=tables[p];P=tt[ti]
                for i,n in enumerate(coeffs):
                    P=add(P,pts[i][n+bound],p)
                if not possible_mod(c,kind,P,scale,p):
                    failed=p;break
            if failed:
                stats['modular_rejected']+=1;rejected[failed]+=1
                continue
            stats['modular_survivors']+=1
            P=T
            for i,n in enumerate(coeffs):
                key=(i,n)
                if key not in exact_multiples:
                    exact_multiples[key]=mul(basis[i],n)
                P=add(P,exact_multiples[key])
            assert check(P,b),'exact group law produced invalid elliptic point'
            key=tuple(pack(P))
            if key in seen:
                stats['duplicate_materialized_point']+=1
                continue
            seen.add(key);stats['distinct_exact_elliptic_points']+=1
            triple,why=decode(c,kind,pack(P),scale)
            outcomes[why]+=1
            exact.append(dict(coefficients=list(coeffs),torsion_index=ti,point=pack(P),lift=why))
            if triple:
                vals=reconstruct(c,f['representations'][0]['normalized'],triple)
                checks=independent_verify(vals)
                certificates.append(dict(integers=vals,certificates=checks,c=c,kind=kind,
                    elliptic_point=pack(P),coefficients=list(coeffs),torsion_index=ti))
    assert stats['coefficient_torsion_tuples']==(2*bound+1)**rank*len(torsion)
    return dict(c=c,kind=kind,model_b=b,model_scale=scale,provenance=provenance,
       input_basis=[pack(P) for P in basis],basis_claim=False,
       full_torsion_points=[pack(P) for P in torsion],
       coefficient_box=dict(min=-bound,max=bound,dimension=rank,inclusive=True),
       primes=primes,finite_field_method='necessary square/cube images, infinity skipped',
       counts=dict(stats),first_rejection_prime_counts=dict(rejected),lift_outcomes=dict(outcomes),
       materialized_points=exact,certificates=certificates,
       status='verified_counterexample' if certificates else 'bounded_subgroup_miss',
       fibre_rejection=False,wall_seconds=round(time.time()-start,4))


def y_seed_task(f,rep):
    h=rep['sixth_scale'];c=f['c'];src=rep['source']
    u=Q(src['u'],h*h);v=Q(src['v'],h*h)
    assert u**3-v**3==c and u!=v
    spec=f['models']['difference'];s=spec['scale']
    P=(12*c/(u-v)/s**2,36*c*(u+v)/(u-v)/s**3)
    assert check(P,spec['b'])
    r=cache_model(f,'difference')
    if r:
        torsion=[unpack(T) for T in r['torsion_points']]
        torsion_proof=dict(method='cached_elltors',cache=f'models/{spec["b"]}.json')
    else:
        torsion,torsion_proof=torsion_fallback(spec['b'])
    info=dict(source=src,u=str(u),v=str(v),sixth_scale=h,
              already_certified_rejected=bool(f['rejection']),torsion_proof=torsion_proof)
    if torsion is None:
        return dict(c=c,kind='difference',status='seed_imported_torsion_not_certified',
                    input_basis=[pack(P)],provenance=info,fibre_rejection=False)
    return hunt(f,'difference',[P],torsion,info)


def selftest():
    # Verify finite-field commutation and sieve soundness against exact lifted
    # positive points and arbitrary small subgroup points on all six models.
    from fibres import factor,models
    import subprocess
    gp='/workspace/scratch/3b44542ad661/runtime/pari-2.17.4/gp'
    E=4;P=(Q(0),Q(2));assert mul(P,3) is None
    P=(Q(2),Q(2));E=-4
    for n in range(-12,13):
        R=mul(P,n);assert check(R,E)
        for p in PRIMES:
            if E%p:
                assert reduce_point(R,p)==mul(reduce_point(P,p),n,p)
    controls={'plus':(1,8),'minus':(4,1),'square':(63,504),
       'difference':(252,3780),'minus_cubed':(252,3969),
       'plus_cubed':(Q(-63,4),Q(3969,8)),
       'plus_four':(Q(1,4),Q(127,8)),'minus_four':(16,62),
       'four_fourth':(Q(3969,4),Q(257985,8))}
    for kind,P in controls.items():
        s=models(63,factor(63))[kind];P=tuple(map(Q,P));P=(P[0]/s['scale']**2,P[1]/s['scale']**3)
        assert check(P,s['b'])
        for p in PRIMES:
            if (6*63*s['b']*s['scale'])%p:
                assert possible_mod(63,kind,reduce_point(P,p),s['scale'],p)
    generated_controls=0
    for A in range(1,4):
        for B in range(A+1,6):
            c=B**6-A**6
            raw={'plus':(Q(A*A),Q(B**3)), 'minus':(Q(B*B),Q(A**3)),
                'square':(Q(c,A*A),Q(c*B**3,A**3)),
                'minus_cubed':(Q(c*B*B,A*A),Q(c*c,A**3)),
                'plus_cubed':(Q(-c*A*A,B*B),Q(c*c,B**3)),
                'difference':(Q(12*c,B*B-A*A),Q(36*c*(B*B+A*A),B*B-A*A)),
                'plus_four':(Q(A**4,B*B),Q(B**6+c,B**3)),
                'minus_four':(Q(B**4,A*A),Q(c-A**6,A**3)),
                'four_fourth':(Q(c*c,A*A*B*B),Q(c*c*(B**6+A**6),A**3*B**3))}
            mm=models(c,factor(c))
            for kind,P in raw.items():
                s=mm[kind];P=(P[0]/s['scale']**2,P[1]/s['scale']**3)
                assert check(P,s['b'])
                t,reason=decode(c,kind,pack(P),s['scale'])
                assert t is not None and (t[0]/t[1],t[2]/t[1])==(A,B)
                for p in PRIMES:
                    if (6*c*s['b']*s['scale'])%p:
                        assert possible_mod(c,kind,reduce_point(P,p),s['scale'],p)
                generated_controls+=1
    # Independent PARI arithmetic comparison for nontrivial rational sums.
    code='E=ellinit([0,-4]);for(n=-12,12,print(ellmul(E,[2,2],n)));quit;\n'
    out=subprocess.run([gp,'-fq'],input=code,text=True,capture_output=True,check=True)
    assert not out.stderr
    for n,line in zip(range(-12,13),out.stdout.splitlines()):
        text=line.strip()[1:-1].split(',')
        assert unpack(text)==mul((Q(2),Q(2)),n)
    bs=(-108,-432,-1728,-16,-4,4,1)
    code=''.join(f'print(elltors(ellinit([0,{b}]))[1]);' for b in bs)+'quit;\n'
    out=subprocess.run([gp,'-fq'],input=code,text=True,capture_output=True,check=True)
    assert not out.stderr
    for b,line in zip(bs,out.stdout.splitlines()):
        tt,proof=torsion_fallback(b)
        assert tt is not None and len(tt)==int(line)
    print(json.dumps(dict(selftest='PASS',pari_group_comparisons=25,
      generic_fibre_controls=len(controls),independent_torsion_comparisons=len(bs),
      generated_positive_fibre_controls=generated_controls,
      modular_commutation_primes=len(PRIMES))),flush=True)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--ledger',type=Path,default=HERE/'ledger.json')
    ap.add_argument('--output',type=Path,default=HERE/'point_hunt_results.json')
    ap.add_argument('--selftest',action='store_true')
    args=ap.parse_args()
    if args.selftest:
        selftest();return
    raw=args.ledger.read_bytes();data=json.loads(raw);ff=data['fibres']
    result=dict(schema='esop6-bounded-point-hunt-v1',ledger_sha256=hashlib.sha256(raw).hexdigest(),
       point_hunt_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
       fibres_script_sha256=hashlib.sha256((HERE/'fibres.py').read_bytes()).hexdigest(),
       ledger_path=str(args.ledger.relative_to(HERE)),ledger_fibres=len(ff),
       scope='one lowest estimated growth model with known points per unresolved fibre; both retained Y seeds',
       supported_model_kinds=SUPPORTED_KINDS,
       rank1_coefficient_bound=12,rank2_rank3_coefficient_bound=2,max_selected_generators=3,
       unselected_generators_policy='lowest coordinate bit-height three; others explicitly omitted',
       ranks_recomputed=False,full_mordell_weil_basis_claim=False,bounded_misses_reject_fibres=False,
       tasks=[],skipped=[],counterexamples=[])
    start=time.time()
    for f in ff:
        if f['rejection']:
            continue
        choice=choose_model(f)
        if choice is None:
            result['skipped'].append(dict(c=f['c'],reason='no_certified_model_with_known_nontorsion_point'))
            continue
        estimate,kind,r,basis=choice
        t=hunt(f,kind,basis,[unpack(P) for P in r['torsion_points']],
          dict(lane='unresolved_fibre',model_cache=f'models/{r["b"]}.json',
               available_known_points=len(r['known_independent_points']),selected_known_points=len(basis),
               rank_lower=r['rank_lower'],rank_upper=r['rank_upper'],estimated_growth=estimate))
        result['tasks'].append(t);result['counterexamples']+=t['certificates']
        print(json.dumps(dict(c=t['c'],kind=kind,counts=t['counts'],outcomes=t['lift_outcomes'])),flush=True)
        save_json(args.output,result)
    for f in ff:
        for rep in f['representations']:
            if 'u' not in rep['source']:
                continue
            t=y_seed_task(f,rep);result['tasks'].append(t)
            result['counterexamples']+=t.get('certificates',[])
            print(json.dumps(dict(event='Y_seed',c=f['c'],status=t['status'],counts=t.get('counts'))),flush=True)
            save_json(args.output,result)
    total=Counter()
    for t in result['tasks']:
        total.update(t.get('counts',{}))
    result['summary']=dict(tasks=len(result['tasks']),skipped=len(result['skipped']),
       counts=dict(total),verified_counterexamples=len(result['counterexamples']),
       elapsed_seconds=round(time.time()-start,3))
    save_json(args.output,result)
    print(json.dumps(result['summary']),flush=True)


if __name__=='__main__':
    main()
