#!/usr/bin/env python3
"""Write results/astra_third_2026_09_11/manifest.json (sha256 of every file)."""
import datetime
import hashlib
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/astra_third_2026_09_11'
FILES = ['ASTRA_THIRD_STRIKE.md',
         'tools/degree_eight_strike.py',
         'tools/unequal_slope_six.py',
         'tools/third_manifest.py',
         'tests/third_strike.py',
         'results/astra_third_2026_09_11/degree_eight.json',
         'results/astra_third_2026_09_11/unequal_slope_six.json',
         'results/astra_third_2026_09_11/mod3_shapes.json',
         'results/astra_third_2026_09_11/third_gates.json',
         'results/astra_third_2026_09_11/inherited_gates.log']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    base = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT,
                                   text=True).strip()
    entries = {}
    for rel in FILES:
        p = ROOT / rel
        if p.exists():
            entries[rel] = {'bytes': p.stat().st_size, 'sha256': sha(p)}
    man = {
        'repository': 'ESOP6',
        'branch': 'claude/great-dirac-nd86f4',
        'base_commit': base,
        'created_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'status': {
            'esop6_solution_found': False,
            'positive_rational_surface_point': False,
            'degree_eight_identity_found': False,
            'degree_eight_obstruction_obtained': False,
            'degree_eight_top_equation_equivalent_to_esop6': True,
            'degree_eight_3_integral_case': 'IMPOSSIBLE',
            'degree_eight_remaining_branch': 'LIVE: deg N = 2',
            'unequal_slope_six_degM0': 'IMPOSSIBLE for every rational rho>0',
            'unequal_slope_six_degM6': 'OPEN (3-integral case impossible)',
            'new_integer_search_run': False,
        },
        'scope': 'RATIONAL_CURVE_ATTEMPT.md equation (5); and the unequal '
                 'leading-slope degree-six sibling',
        'files': dict(sorted(entries.items())),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'manifest.json').write_text(json.dumps(man, indent=2) + '\n')
    print(json.dumps({'result': 'PASS', 'files': len(entries)}, indent=2))


if __name__ == '__main__':
    main()
