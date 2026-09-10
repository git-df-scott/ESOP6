#!/usr/bin/env python3
"""Produce reviewable fibre and quadratic-Chabauty target ledgers."""
from collections import Counter
import json
from pathlib import Path
from fibres import HERE,KINDS,save_json

ROUTES={'t6_plus_c':('plus','square'), 't6_minus_c':('minus','square'),
        't6_plus_4c':('plus_four','difference'), 't6_minus_4c':('minus_four','difference'),
        't6_plus_c2_over4':('difference','four_fourth')}


def main():
    d=json.loads((HERE/'ledger.json').read_text())
    rows=sorted(d['fibres'],key=lambda f:f['c'])
    survivors=[]
    table=['# Complete normalized fibre ledger','',
       'All 627 retained fibres. Exact representations, scales, local charts and raw arithmetic are in arithmetic_evidence.zip.',
       'A dash means the later test was unnecessary after an earlier certified rejection. Timeout means unresolved.',
       'The three local numbers are unit counts at 2,3,7; they are chart restrictions, not fibre obstructions.','',
       '| c | Fourtuple representation(s) | Local counts | Ranks in execution order | Status |',
       '|---:|---|---|---|---|']
    for f in rows:
        entries=f['elliptic_data'];cells=[]
        for kind in KINDS:
            e=entries.get(kind)
            if e is None:cells.append('-')
            elif e['status']!='certified':cells.append(e['status'])
            else:
                lo,hi=e['rank_lower'],e['rank_upper']
                cells.append(str(lo) if lo==hi else f'{lo}..{hi}')
        reps='; '.join('('+','.join(r['normalized'])+')' for r in f['representations'])
        counts=','.join(str(f['local_data']['charts'][str(p)]['unit_count']) for p in (2,3,7))
        status='excluded: '+f['rejection']['kind'] if f['rejection'] else 'unresolved'
        table.append(f'| {f["c"]} | {reps} | {counts} | {", ".join(cells)} | {status} |')
        if f['rejection']:continue
        s=dict(c=f['c'],representations=f['representations'],factorization=f['factorization'],
               elliptic_data=entries,genus2_rank_one_routes=[])
        for name,pair in ROUTES.items():
            if all(entries[k]['status']=='certified' and
                   entries[k]['rank_lower']==entries[k]['rank_upper']==1 for k in pair):
                s['genus2_rank_one_routes'].append(name)
        if all(e['status']=='certified' for e in entries.values()) and len(entries)==len(KINDS):
            s['jacobian_rank_lower']=sum(e['rank_lower'] for e in entries.values())+entries['square']['rank_lower']
            s['jacobian_rank_upper']=sum(e['rank_upper'] for e in entries.values())+entries['square']['rank_upper']
            s['genus10_qc_finiteness_criterion']=s['jacobian_rank_upper']<=19
        else:
            s['jacobian_rank_lower']=None;s['jacobian_rank_upper']=None
            s['genus10_qc_finiteness_criterion']=None
        survivors.append(s)
    table[3]+=' Execution order: '+', '.join(KINDS)+'.'
    (HERE/'FIBRE_LEDGER.md').write_text('\n'.join(table)+'\n')
    save_json(HERE/'survivors.json',dict(count=len(survivors),fibres=survivors,
            exact_point_classification_executed=False,rational_NS_lower_bound=11))
    summary=dict(fibres=len(rows),representations=d['counts']['stored_representations'],
       primitive_box_rays=d['counts']['sorted_tuples']-d['counts']['nonprimitive_tuples_removed'],
       box_height=d['height'],certified_excluded=len(rows)-len(survivors),unresolved=len(survivors),
       stage_rejections=d['arithmetic_summary'],model_status=d['model_summary'],
       rank_one_genus2_fibres=sum(bool(s['genus2_rank_one_routes']) for s in survivors),
       complete_nine_rank_intervals=sum(s['jacobian_rank_upper'] is not None for s in survivors),
       genus10_qc_eligible=sum(s['genus10_qc_finiteness_criterion'] is True for s in survivors),
       counterexample=False)
    point=HERE/'point_hunt_results.json'
    if point.exists():summary['bounded_point_search']=json.loads(point.read_text())['summary']
    save_json(HERE/'search_summary.json',summary)
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
