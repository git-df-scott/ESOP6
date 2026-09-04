#!/usr/bin/env python3
"""Cross-check the standalone Python and JavaScript integer certificates.

No positive ESOP6 solution is available as a fixture. The controls include
arithmetically equal boundary inputs, positivity failures, large integers,
common factors, and independent exact expected values.
"""
import json
import math
from pathlib import Path
import random
import subprocess
import sys

root=Path(__file__).resolve().parents[1]
cases=[[0,0,0,0,1,1],[0,0,0,0,0,0],[-1,0,0,0,0,1],
       [1,1,1,1,1,1],[2,4,6,8,10,12],
       [10**100+i for i in range(6)]]
rng=random.Random(660612)
cases.extend([[rng.randrange(1,10**60) for _ in range(6)] for _ in range(12)])
for values in cases:
    arguments=list(map(str,values))
    py=subprocess.run([sys.executable,str(root/'tools/verify_esop6.py'),*arguments],
                      text=True,capture_output=True)
    js=subprocess.run(['node',str(root/'tools/verify_esop6.mjs'),*arguments],
                      text=True,capture_output=True)
    a,b=json.loads(py.stdout),json.loads(js.stdout)
    assert a==b,(values,a,b)
    assert a['sixth_powers']==[str(n**6) for n in values]
    assert a['positivity']==all(n>0 for n in values)
    assert a['equality']==(sum(n**6 for n in values[:5])==values[5]**6)
    assert a['gcd']==str(math.gcd(*values))
    assert py.returncode==js.returncode==(0 if a['solution'] else 1)
print(json.dumps({'result':'PASS','cases':len(cases),
                  'implementations':['Python integers','JavaScript BigInt'],
                  'known_positive_ESOP6_fixture':False},indent=2))
