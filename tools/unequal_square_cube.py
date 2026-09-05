"""Two forced square ratios, unequal pair sums, with CRT-admissible f/z.

This is a bounded rational construction search, not a height census.
"""
import argparse, json, math, time
from pathlib import Path
from sympy.solvers.diophantine.diophantine import diop_DN

OUT=Path('results/multiplicity_square_2026_09_05')
parser=argparse.ArgumentParser()
parser.add_argument('--replay-rejected-control',action='store_true',help='Replay the diagnostic run made before its exact norm-ratio obstruction was recognized')
args=parser.parse_args()
u,v=1,182
if not args.replay_rejected_control:
    gate={'u':u,'v':v,'prime':7,'valuation_of_u_over_v':-1,'globally_excluded':True,
          'reason':'A ratio of two nonzero sums of two rational squares has even valuation at every prime congruent to 3 modulo 4.',
          'production_search_launched':False,'diagnostic_replay_flag':'--replay-rejected-control'}
    (OUT/'unequal_square_cube_gate.json').write_text(json.dumps(gate,indent=2)+'\n')
    print(json.dumps(gate,indent=2))
    raise SystemExit(0)
sa,sb,sc=1,2,3
h=(sa*sa*u+sb*sb*v)//(sc*sc)
assert h==81
B=(sa,-sa,sb,-sb,sc,sc);sign=(1,1,1,1,1,-1)
slopes=[(1,0)]+[(m,n) for n in range(1,25) for m in range(-24,25) if math.gcd(m,n)==1]
mods=(64,9,5,7,13,19,31)
masks=[{x*x%p for x in range(p)} for p in mods]
def square_gate(z):return all(z%p in ss for p,ss in zip(mods,masks))
def square(z):return math.isqrt(z)**2==z
counts={'ratio_seeds':0,'positive_norms':0,'integer_norm_representations':0,'used_norm_representations':0,
        'directions':0,'positive_cubic_outputs':0,'square_sieve_pass':0,'exact_lifts':0}
samples=[];norm_ledger=[];best=-1;bestC=None;start=time.monotonic()
for line in (OUT/'square_cube_ratio_seeds.csv').read_text().splitlines():
    a,b=map(int,line.split(','));counts['ratio_seeds']+=1
    assert math.gcd(a,b)==1 and (a**6-b**6)%(42**6)==0
    d=a*a-b*b
    K=sc*sc*(d*d*((4*h**3-u**3-v**3)//3)+4*a*a*b*b*h**3)
    if K<=0:
        norm_ledger.append({'a':a,'b':b,'K':str(K),'reason':'negative real norm'});continue
    counts['positive_norms']+=1
    seeds=sorted(set(tuple(map(int,xy)) for xy in diop_DN(-v,K)))
    counts['integer_norm_representations']+=len(seeds)
    norm_ledger.append({'a':a,'b':b,'K':str(K),'integer_norm_representations':len(seeds),'used':seeds[:16]})
    for X,Y in seeds[:16]:
        assert X*X+v*Y*Y==K
        counts['used_norm_representations']+=1
        for sx,sy in ((1,1),(1,-1),(-1,1),(-1,-1)):
            for m,n in slopes:
                counts['directions']+=1
                E=m*m+v*n*n
                XX=(m*m-v*n*n)*sx*X-2*v*m*n*sy*Y
                YY=2*m*n*sx*X+(m*m-v*n*n)*sy*Y
                assert XX*XX+v*YY*YY==K*E*E
                # x=(XX/E-2*h*sa*b^2)/(sc*d), and similarly y.
                xx=XX-2*h*sa*b*b*E;yy=YY-2*h*sb*b*b*E;den=sc*d*E
                D=(u*den+xx,u*den-xx,v*den+yy,v*den-yy,0,2*h*den)
                g=math.gcd(*D);D=tuple(x//g for x in D)
                F=sum(e*x**3 for e,x in zip(sign,D))
                kk=sum(e*z*x*x for e,z,x in zip(sign,B,D))
                C=tuple(F*z-3*kk*x for z,x in zip(B,D))
                assert C[-1]*b*b==C[-2]*a*a
                assert v*(C[0]+C[1])==u*(C[2]+C[3])
                if C[-1]<0:C=tuple(-x for x in C)
                if min(C)<=0:continue
                counts['positive_cubic_outputs']+=1
                g=math.gcd(*C);C=tuple(x//g for x in C)
                assert sum(x**3 for x in C[:5])==C[5]**3
                if len(samples)<32:samples.append({'a':a,'b':b,'C':[str(x) for x in C]})
                # Ci/C5 is rational square iff Ci*C5 is integer square.
                # Local residue sieve is necessary and applied before square roots.
                products=[x*C[4] for x in C[:4]]
                passes=[square_gate(z) for z in products]
                if all(passes):counts['square_sieve_pass']+=1
                score=sum(ok and square(z) for ok,z in zip(passes,products))
                if score>best:best,bestC=score,C
                if score==4:
                    assert all(square(x) for x in C)
                    sol=[math.isqrt(x) for x in C]
                    assert sum(x**6 for x in sol[:5])==sol[5]**6
                    counts['exact_lifts']+=1
                    (OUT/'CANDIDATE.json').write_text(json.dumps(sol)+'\n')
                    print('EXACT COUNTEREXAMPLE',sol,flush=True)
    print(a,b,'directions',counts['directions'],'positive',counts['positive_cubic_outputs'],flush=True)
record={'seed':B,'pair_sum_ratio':[u,v],'h':h,'counts':counts,'slopes':len(slopes),
        'norm_representation_cap_per_ratio':16,'best_remaining_square_ratios':best,
        'best_C':[str(x) for x in bestC] if bestC else None,'norm_ledger':norm_ledger,
        'samples':samples,'seconds':time.monotonic()-start,
        'scope':'bounded conic rotations with unequal pair sums and locally admissible RHS/fifth ratio; no height interval exhausted'}
(OUT/'unequal_square_cube.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({k:v for k,v in record.items() if k not in ('samples','norm_ledger')},indent=2))
