#!/usr/bin/env python3
"""Verify retained Y2 lists without repeating their search or launching jobs."""
import hashlib
import json
import math
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
SRC=ROOT/'research/cube_ansatz_2026_09_10'


def main():
    audit={}
    for p in sorted(SRC.glob('y2_*.json')):
        rows=json.loads(p.read_text())
        good=[];bad=[];seen=set();duplicates=0
        for i,r in enumerate(rows):
            if not isinstance(r,dict) or not all(k in r for k in ('x','x6','S')):
                continue
            x,f,s=r['x'],r['x6'],r['S']
            residual=sum(a**6 for a in x)+s*s-f**6
            if residual != 0 or len(x)!=4 or not all(a>0 for a in x) or f<=0 or s==0:
                bad.append(dict(index_zero_based=i,point=r,residual=str(residual)))
                continue
            key=(tuple(sorted(x)),f,s)
            if key in seen:duplicates+=1;continue
            seen.add(key);good.append(r)
        audit[p.name]=dict(source_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
                          rows=len(rows),valid_distinct=len(good),duplicates=duplicates,invalid=bad)
        if p.name=='y2_points_N500.json':
            (HERE/'y2_points_N500_verified.json').write_text(json.dumps(good,indent=1)+'\n')
            for i in range(4):
                chunk=[r for r in good if r['x6']>300][i::4]
                (HERE/f'y2_b_verified_chunk{i}.json').write_text(json.dumps(chunk,indent=1)+'\n')
    (HERE/'seed_validation.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps({p:{k:v for k,v in a.items() if k!='invalid'}|{'invalid_count':len(a['invalid'])}
                      for p,a in audit.items()},indent=2))


if __name__=='__main__':main()
