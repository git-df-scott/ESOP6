"""Deterministic sampled high-height cubic tangent search with exact integers."""
import argparse, json, math, random, time
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--samples',type=int,default=1000000);p.add_argument('--crt-seeds',action='store_true');args=p.parse_args()
out=Path('results/fourfold_routes_2026_09_05')
rng=random.Random(20260905)
primes=(16,63,65,11,13,17,19,23,31)
masks=[{x*x%p for x in range(p)} for p in primes]
count={'sampled_directions':0,'nondegenerate':0,'positive':0,'cubic_identity_checks':0,
 'square_sieve_pass':0,'hits':0,'above_730000_squared':0,'above_4300000_squared':0,'above_100000000_squared':0}
hist={};maxC=0;sample=[];best=[];bestscore=-1;start=time.monotonic()
sign=(1,1,1,1,1,-1)
for n in range(args.samples):
    if args.crt_seeds:
        a,b=[1764*rng.randint(1,8) for _ in range(2)];c=1764*rng.randint(1,8)+1
    else:
        a,b,c=[rng.randint(1,100) for _ in range(3)]
    x,y,z,r=[rng.randint(-10000,10000) for _ in range(4)]
    B=(a,-a,b,-b,c,c)
    D=(c*c*x,c*c*y,c*c*z,c*c*r,0,a*a*(x+y)+b*b*(z+r))
    gd=math.gcd(*D)
    if not gd:continue
    D=tuple(v//gd for v in D)
    count['sampled_directions']+=1
    f=sum(e*v**3 for e,v in zip(sign,D));k=sum(e*b*v*v for e,b,v in zip(sign,B,D))
    if not f or not k:continue
    count['nondegenerate']+=1
    C=tuple(f*b-3*k*v for b,v in zip(B,D))
    if C[-1]<0:C=tuple(-v for v in C)
    if min(C)<=0:continue
    count['positive']+=1
    gc=math.gcd(*C);C=tuple(v//gc for v in C)
    assert sum(v**3 for v in C[:-1])==C[-1]**3
    assert sum(e*b*b*v for e,b,v in zip(sign,B,D))==0
    count['cubic_identity_checks']+=1
    if len(sample)<1024:sample.append({'B':B,'D':D,'C':C})
    maxC=max(maxC,C[-1])
    for h in (730000,4300000,100000000):
        if C[-1]>h*h:count[f'above_{h}_squared']+=1
    score=sum(math.isqrt(v)**2==v for v in C)
    hist[score]=hist.get(score,0)+1
    if score>bestscore:bestscore=score;best=C
    if not all(v%p in mask for v in C for p,mask in zip(primes,masks)):continue
    count['square_sieve_pass']+=1
    if score==6:
        sextuple=[math.isqrt(v) for v in C]
        assert sum(v**6 for v in sextuple[:-1])==sextuple[-1]**6
        count['hits']+=1
        (out/'CANDIDATE.json').write_text(json.dumps(sextuple)+'\n')
        print('EXACT HIT',sextuple,flush=True)
        break
    if n and n%100000==0:print(n,count,flush=True)
report={'random_seed':20260905,'requested_samples':args.samples,'seed_rule':('a,b=1764*i,1764*j; c=1764*k+1; 1<=i,j,k<=8' if args.crt_seeds else '1<=a,b,c<=100'),'direction_bound':10000,
 'counts':count,'square_count_histogram':hist,'largest_primitive_cube_base':str(maxC),
 'best_square_count':bestscore,'best_cube_bases':[str(x) for x in best],
 'elapsed_seconds':time.monotonic()-start,'scope':'sampled tangent points; no integer-height interval is exhausted'}
stem='crt_tangent' if args.crt_seeds else 'large_tangent'
(out/(stem+'.json')).write_text(json.dumps(report,indent=2)+'\n')
# Decimal strings preserve exactness for independent Node verification.
(out/(stem+'_samples.json')).write_text(json.dumps([{k:[str(x) for x in v] for k,v in rec.items()} for rec in sample],indent=2)+'\n')
print(json.dumps(report,indent=2),flush=True)
