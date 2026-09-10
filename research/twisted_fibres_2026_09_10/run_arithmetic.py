#!/usr/bin/env python3
"""Resumable nine-model rejection DAG, certified PARI 2.17.4 only.

Each GP task is independent and has a wall-clock limit. Failures/timeouts
remain unresolved and are cached. No jobs or old lottery files are changed.
"""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import subprocess
import time
from fibres import HERE, KINDS, decode, reconstruct, independent_verify, save_json


def gp_program(b, effort=0):
    # Exact version-pinned inspection is necessary: cubic class groups used
    # in ellrankinit are explicitly certified before rank claims are used.
    return f'''
setrand(20260910);
J(v)={{my(s="[");if(type(v)!="t_VEC",return(Str("\\\"",v,"\\\"")));for(i=1,#v,if(i>1,s=concat(s,","));s=concat(s,J(v[i])));concat(s,"]")}};
print("VERSION:",version());
E=ellinit([0,{b}]);
R=ellrankinit(E);
certs=List();
for(i=1,#R[3],if(#R[3][i]==10,print("FIELD:",R[3][i].pol);z=bnfcertify(R[3][i]);if(z!=1,error("uncertified class group"));listput(certs,z)));
print("CERTS:",J(Vec(certs)));
V=ellrank(R,{effort});
print("RANK:",J(V));
T=elltors(E);
print("TORSION:",J(T));
TT=List();listput(TT,[0]);
for(i=1,#T[2],old=Vec(TT);TT=List();for(j=1,#old,for(k=0,T[2][i]-1,listput(TT,elladd(E,old[j],ellmul(E,T[3][i],k))))));
print("TORSION_POINTS:",J(Vec(TT)));
print("COMPLETE");
quit;
'''


def rank_task(gp, b, seconds, effort):
    program = gp_program(b,effort)
    t=time.time()
    try:
        p=subprocess.run([gp,'-fq','-s','128M'],input=program,capture_output=True,
                         text=True,timeout=seconds)
    except subprocess.TimeoutExpired as e:
        return dict(b=b,status='timeout',wall_seconds=time.time()-t,limit_seconds=seconds,
                    stdout=(e.stdout or b'').decode() if isinstance(e.stdout,bytes) else e.stdout,
                    stderr=(e.stderr or b'').decode() if isinstance(e.stderr,bytes) else e.stderr,
                    effort=effort,program_sha256=hashlib.sha256(program.encode()).hexdigest())
    result=dict(b=b,status='error',wall_seconds=time.time()-t,limit_seconds=seconds,
                stdout=p.stdout,stderr=p.stderr,exit_code=p.returncode,effort=effort,
                program_sha256=hashlib.sha256(program.encode()).hexdigest())
    try:
        lines=p.stdout.splitlines()
        assert p.returncode == 0 and lines[-1] == 'COMPLETE'
        assert 'VERSION:[2, 17, 4]' in lines
        # GP can continue after an error; a COMPLETE marker alone is not enough.
        assert '***' not in p.stderr
        fields={s.split(':',1)[0]:s.split(':',1)[1] for s in lines if ':' in s}
        certs=json.loads(fields['CERTS'])
        assert all(int(v)==1 for v in certs)
        rank=json.loads(fields['RANK'])
        tors=json.loads(fields['TORSION'])
        torspts=json.loads(fields['TORSION_POINTS'])
        assert len(torspts) == int(tors[0])
        for pnt in rank[3]+torspts:
            if len(pnt)==2:
                x,y=map(Q,pnt)
                assert y*y == x*x*x+b
            else:
                assert pnt == ['0']
        result.update(status='certified',rank_lower=int(rank[0]),rank_upper=int(rank[1]),
                      sha2_lower_dimension=int(rank[2]),known_independent_points=rank[3],
                      torsion=tors,torsion_points=torspts,bnf_certificates=certs,
                      bnf_polynomials=[s[6:] for s in lines if s.startswith('FIELD:')])
    except (AssertionError,ValueError,KeyError,IndexError) as e:
        result['parse_error']=repr(e)
    return result


def process_fibre(f,kind,r,outdir):
    spec=f['models'][kind]
    f['elliptic_data'][kind]=dict(model_b=spec['b'],model_cache_file=f'models/{spec["b"]}.json',status=r['status'])
    if r['status'] != 'certified':
        return
    entry=f['elliptic_data'][kind]
    entry.update(rank_lower=r['rank_lower'],rank_upper=r['rank_upper'],torsion_order=int(r['torsion'][0]))
    reasons=Counter()
    # These are available exact points, not a complete group when rank>0.
    for p in r['torsion_points']+r['known_independent_points']:
        for signed in ([p] if len(p)==1 or Q(p[1])==0 else [p,[p[0],str(-Q(p[1]))]]):
            triple,why=decode(f['c'],kind,signed,spec['scale'])
            reasons[why]+=1
            if triple:
                vals=reconstruct(f['c'],f['representations'][0]['normalized'],triple)
                certs=independent_verify(vals)
                save_json(outdir/'EXACT_COUNTEREXAMPLE.json',dict(integers=vals,certificates=certs,fibre=f['c'],kind=kind))
                raise RuntimeError('EXACT COUNTEREXAMPLE VERIFIED; stop search')
    entry['point_lift_outcomes']=dict(reasons)
    if r['rank_upper']==0:
        f['square_cover_status']='certified_no_positive_point'
        f['rejection']=dict(kind=kind,model_b=spec['b'],reason='certified_rank_zero_all_torsion_failed_lifts',
                           torsion_points=r['torsion_points'],outcomes=dict(reasons))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--gp',required=True)
    ap.add_argument('--ledger',default=str(HERE/'ledger.json'))
    ap.add_argument('--seconds',type=float,default=10)
    ap.add_argument('--workers',type=int,default=3)
    ap.add_argument('--effort',type=int,default=0)
    ap.add_argument('--kinds',nargs='+',choices=KINDS,default=list(KINDS))
    ap.add_argument('--retry-failed',action='store_true')
    args=ap.parse_args()
    path=Path(args.ledger); outdir=path.parent; cache=outdir/'models';cache.mkdir(exist_ok=True)
    data=json.loads(path.read_text());fibres=data['fibres']
    data['arithmetic_run']=dict(version='2.17.4',seconds=args.seconds,workers=args.workers,
          effort=args.effort,seed=20260910,kinds=args.kinds,rank_basis_claim=False,
          class_groups_explicitly_certified=True)
    start=time.time()
    for kind in args.kinds:
        todo={}
        for f in fibres:
            if f['rejection'] is None:
                b=f['models'][kind]['b']
                todo.setdefault(b,[]).append(f)
        print(json.dumps(dict(event='stage_start',kind=kind,models=len(todo),fibres=sum(map(len,todo.values())))),flush=True)
        futures={}
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            for b, ff in todo.items():
                p=cache/f'{b}.json'
                r=json.loads(p.read_text()) if p.exists() else None
                if r and not (args.retry_failed and r['status']!='certified'):
                    for f in ff:process_fibre(f,kind,r,outdir)
                else:
                    futures[pool.submit(rank_task,args.gp,b,args.seconds,args.effort)]=b
            complete=0
            for future in as_completed(futures):
                b=futures[future];r=future.result();save_json(cache/f'{b}.json',r)
                for f in todo[b]:process_fibre(f,kind,r,outdir)
                complete+=1
                if complete%25==0:
                    save_json(path,data)
                    print(json.dumps(dict(event='progress',kind=kind,models_completed=complete,
                       certified_rejected=sum(f['rejection'] is not None for f in fibres),
                       elapsed_seconds=round(time.time()-start,2))),flush=True)
        stats=Counter(f['rejection']['kind'] if f['rejection'] else 'unresolved' for f in fibres)
        data['arithmetic_summary']=dict(stats)
        save_json(path,data)
        print(json.dumps(dict(event='stage_complete',kind=kind,dispositions=dict(stats),
                            elapsed_seconds=round(time.time()-start,2))),flush=True)
    models=[json.loads(p.read_text()) for p in cache.glob('*.json')]
    data['model_summary']=dict(Counter(r['status'] for r in models))
    save_json(path,data)


if __name__=='__main__':main()
