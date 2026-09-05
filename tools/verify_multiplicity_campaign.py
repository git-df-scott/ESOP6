"""Independent arithmetic, finite algebra, and bounded-engine checks."""
import itertools, json, math, random, subprocess, sys
from pathlib import Path
import sympy as S

OUT=Path('results/multiplicity_square_2026_09_05')
exe=sys.argv[1] if len(sys.argv)>1 else '/tmp/repeated_campaign_final'
def run(flag,data=''):
    return subprocess.run([exe,flag],input=data,text=True,capture_output=True,check=True).stdout

# Independently form the CRT root list, instead of stepping along 7-adic rays.
from sympy.ntheory.modular import crt
ms=(64,729,117649)
rs=[[x for x in range(m) if pow(x,6,m)==1] for m in ms]
expected=sorted(int(crt(ms,co)[0]) for co in itertools.product(*rs))
assert sorted(map(int,run('--roots').split()))==expected and len(expected)==144

rng=random.Random(20260905);M=42**6
seeds=[tuple(map(int,s.split(','))) for s in (OUT/'square_cube_ratio_seeds.csv').read_text().splitlines()]
inputs=[]
for i in range(1200):
    a,b=rng.choice(seeds);g=rng.randrange(1,100000000//a+1);f,t=g*a,g*b
    den=rng.choice((M,2*M,3*M,5*M,11*M))
    inputs.append((f,t,den))
for i in range(800):
    f=rng.randrange(1000000,100000001);t=rng.randrange(1,f);inputs.append((f,t,M))
answers=run('--quotients',''.join(f'{f} {t} {d}\n' for f,t,d in inputs)).splitlines()
assert len(answers)==len(inputs)
for (f,t,d),ans in zip(inputs,answers):
    n=f**6-t**6
    assert ans==(str(n//d) if n%d==0 else 'NONE')

pair_inputs=[(1,1),(1000000,1000000),(999999,1000000),(1999999,2000000)]
pair_inputs += [tuple(sorted((rng.randrange(1,1000001),rng.randrange(1,1000001)))) for _ in range(60)]
answers=run('--pairs',''.join(f'{a**6+b**6} {b}\n' for a,b in pair_inputs)).splitlines()
for (a,b),line in zip(pair_inputs,answers):
    ps=json.loads(line);assert [a,b] in ps
    assert all(x**6+y**6==a**6+b**6 and 1<=x<=y<=b for x,y in ps)

# Independently substitute every retained case using SymPy's polynomial ring.
t=S.symbols('t')
certificate=json.loads((OUT/'degree8_mod7.json').read_text())
seen=set()
for row in certificate['cases']:
    q,a,b=(row[k] for k in ('q','a','b'));seen.add((q,a,b))
    P=S.Poly((1-t**6)*(1+a*t),t,modulus=7)
    R=S.Poly((1-t**6)*(1+b*t),t,modulus=7)
    Q=S.Poly(1+q*t*t,t,modulus=7)
    N=S.Poly.from_list(list(reversed(row['N'])),t,modulus=7)
    tt=S.Poly(t,t,modulus=7)
    F=4*(P**6+R**6)-Q*N**5-(10*pow(27,-1,7))*tt**12*Q**3*N**3-pow(81,-1,7)*tt**24*Q**5*N
    first=next(i for i in range(43) if F.nth(i)%7)
    assert first==row['first_failure'] and int(F.nth(first))%7==row['residue'] and first>=9
assert seen==set(itertools.product(range(7),repeat=3))

# Independent direct five-term enumeration at 2 and 3, rather than MITM.
local={}
for cm,mod in ((4,8),(3,9)):
    powers=set()
    for co in itertools.product(range(cm),repeat=3):
        poly=S.Poly(sum(co[j]*t**j for j in range(3)),t)**6
        powers.add(tuple(int(poly.nth(i))%mod for i in range(13)))
    vec=sorted(powers);targets=set(vec)-{(0,)*13};solutions=0
    for ids in itertools.combinations_with_replacement(range(len(vec)),5):
        key=tuple(sum(vec[j][k] for j in ids)%mod for k in range(13))
        if key in targets:solutions+=1
    assert solutions==({8:28,9:13}[mod]);local[mod]=solutions

# Independently count the small 2111 tables, including the exact cutoff.
H=1000;counts=[0]*4
for b in range(7,H,7):
    for c in range(b,H,7):
        n2=b%2+c%2;n3=int(b%3!=0)+int(c%3!=0)
        if n2<=1 and n3<=1 and b**6+c**6+2*42**6+1<=H**6:counts[n2+2*n3]+=1
small=json.loads(subprocess.check_output([exe,'2111',str(H)],text=True,stderr=subprocess.DEVNULL))
assert small['pair_table_sizes']==counts
controls=json.loads(run('--self-test'))
report={'status':'PASS','root_list_independent_CRT':144,'python_bigint_quotients':len(inputs),
        'large_positive_pair_controls':len(pair_inputs),'degree8_sympy_cases':len(seen),
        'direct_conic_census':local,'independent_small_pair_table_sizes':counts,'engine_selftest':controls,
        'limits':'Production negative pair queries were not all rerun with a second pair algorithm; exact algorithm proof and positive/differential controls support the run.'}
(OUT/'independent_verification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
