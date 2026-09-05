"""Exact necessary conditions. No finite local census asserts global existence."""
import hashlib
import itertools
import json
from collections import Counter, defaultdict
from pathlib import Path

OUT = Path('results/multiplicity_square_2026_09_05')
OUT.mkdir(parents=True, exist_ok=True)

def mul(a, b, mod=None):
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i+j] += x*y
    return [x % mod for x in c] if mod else c

def power(a, k, mod=None):
    c = [1]
    for _ in range(k):
        c = mul(c, a, mod)
    return c

def add(a, b, fac=1, shift=0, mod=7):
    c = a + [0] * max(0, len(b) + shift - len(a))
    for i, x in enumerate(b):
        c[i+shift] += fac*x
    return [x % mod for x in c]

def degree8():
    rows = []
    for q, a, b in itertools.product(range(7), repeat=3):
        P, R = [1,a,0,0,0,0,-1,-a], [1,b,0,0,0,0,-1,-b]
        Q = [1,0,q]
        D = [v*4 % 7 for v in add(power(P,6,7), power(R,6,7))]
        N = [1]
        for n in range(1,9):
            rhs = mul(Q, power(N,5,7),7)
            N.append((D[n] - (rhs[n] if n < len(rhs) else 0))*3 % 7)
        rhs = add(mul(Q,power(N,5,7),7), mul(power(Q,3,7),power(N,3,7),7), 10*pow(27,-1,7),12)
        rhs = add(rhs,mul(power(Q,5,7),N,7),pow(81,-1,7),24)
        residual = add(D,rhs,-1)
        first = next((i for i,x in enumerate(residual) if x),None)
        assert first is not None and first >= 9
        rows.append({'q':q,'a':a,'b':b,'N':N,'first_failure':first,'residue':residual[first]})
    record = {'modulus':7,'all_q_cases':343,'square_q_cases':196,'survivors':0,
              'first_failure_histogram':dict(Counter(x['first_failure'] for x in rows)),
              'scope':'Normalized P,R,N,q all in Z_(7); rational coefficients with a 7 denominator not excluded.',
              'cases':rows}
    (OUT/'degree8_mod7.json').write_text(json.dumps(record,indent=2)+'\n')
    return {k:v for k,v in record.items() if k!='cases'}

def rank_mod(rows,p):
    a=[list(x) for x in rows];r=0
    for j in range(3):
        k=next((k for k in range(r,len(a)) if a[k][j]%p),None)
        if k is None: continue
        a[r],a[k]=a[k],a[r];inv=pow(a[r][j]%p,-1,p)
        a[r]=[v*inv%p for v in a[r]]
        for k in range(len(a)):
            if k!=r:
                v=a[k][j];a[k]=[(x-v*y)%p for x,y in zip(a[k],a[r])]
        r+=1
    return r

def conic_census(p, coefficient_modulus, equation_modulus):
    # Sixth powers mod 8 depend only on coefficients mod 4; mod 9 only
    # on coefficients mod 3. Duplicate polynomial powers are identified.
    reps={}
    for co in itertools.product(range(coefficient_modulus),repeat=3):
        key=tuple(power(co,6,equation_modulus))
        reps.setdefault(key,co)
    vec=list(reps);forms=list(reps.values());M=equation_modulus;n=len(vec)
    triples=defaultdict(list)
    for ids in itertools.combinations_with_replacement(range(n),3):
        triples[tuple(sum(vec[i][k] for i in ids)%M for k in range(13))].append(ids)
    solutions=set()
    for rhs in range(n):
        if p!=5 and not any(x%p for x in forms[rhs]):continue
        for i,j in itertools.combinations_with_replacement(range(n),2):
            key=tuple((vec[rhs][k]-vec[i][k]-vec[j][k])%M for k in range(13))
            for tri in triples.get(key,[]):
                if any(x%p for index in (i,j)+tri+(rhs,) for x in forms[index]):
                    solutions.add((tuple(sorted((i,j)+tri)),rhs))
    hist=Counter();examples={};rhs_forms=set()
    for ids,rhs in sorted(solutions):
        rr=rank_mod([forms[i] for i in ids]+[forms[rhs]],p)
        hist[rr]+=1;rhs_forms.add(forms[rhs])
        examples.setdefault(rr,{'left':[forms[i] for i in ids],'right':forms[rhs]})
    payload=json.dumps(sorted(solutions),separators=(',',':')).encode()
    return {'prime':p,'coefficient_modulus':coefficient_modulus,'equation_modulus':M,
            'distinct_sixth_polynomials':n,'solution_multisets':len(solutions),
            'span_rank_mod_p_histogram':dict(hist),'examples_by_rank':examples,
            'rhs_representatives':sorted(rhs_forms),'solution_ids_sha256':hashlib.sha256(payload).hexdigest(),
            'scope':'complete finite coefficient census; each prime is a separate necessary test, not a global rational-point certificate'}

def local_multiplicities():
    patterns=[(1,1,1,1,1),(2,1,1,1),(2,2,1),(3,1,1),(3,2),(4,1),(5,)]
    data=[]
    for pattern in patterns:
        row={'pattern':pattern,'singleton_slots':[i for i,m in enumerate(pattern) if m==1]}
        # At 8,9,7 every unit sixth power is 1. For a primitive solution,
        # f is a unit at the associated prime, hence the weighted count is 1.
        row['unit_supports']=[bits for bits in itertools.product((0,1),repeat=len(pattern)) if sum(a*b for a,b in zip(pattern,bits))==1]
        row['role_assignments_2_3_7']=len(row['unit_supports'])**3
        data.append(row)
    return data

if __name__=='__main__':
    record={'multiplicity_roles':local_multiplicities(),'degree8':degree8(),'conics':[]}
    print(json.dumps(record,indent=2),flush=True)
    for params in [(2,4,8),(3,3,9),(5,5,5),(7,7,7)]:
        row=conic_census(*params);record['conics'].append(row)
        print(json.dumps(row,indent=2),flush=True)
    (OUT/'algebra_local.json').write_text(json.dumps(record,indent=2)+'\n')
