#!/usr/bin/env python3
"""Regression and independent verification of the second-strike proof."""
import json
from fractions import Fraction
from pathlib import Path
import random
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
import boundary_contact as bc
import boundary_symbolics as symbolic


def main():
    exact=symbolic.checks()
    assert exact['result']=='PASS'
    # An earlier development-only zero-padding path created floats. Guard
    # the entire formal jet, including zero coefficients, against recurrence.
    for a,b,n in [([1],[1],72),([1,1],[1,-1],30)]:
        assert all(isinstance(c,Fraction) for c in bc.formal_midpoint(a,b,n))
    space=bc.LiftSpace();rng=random.Random(660602)
    for n in range(1,4):
        place=3**n;modulus=9*place
        for _ in range(8):
            q=[a+3*rng.randrange(27) for a in bc.BASE]
            h=[rng.randrange(3) for _ in range(9)]
            before=bc.residual(q,modulus)
            after=bc.residual([x+place*y for x,y in zip(q,h)],modulus)
            delta=[(x-y)%modulus for x,y in zip(after,before)]
            expected=[3*place*sum(a*b for a,b in zip(row,h))%modulus for row in space.matrix]
            assert delta==expected
    py=bc.certificate()
    cpp=json.loads(subprocess.check_output([str(ROOT/'bin/boundary_contact_verify')],text=True))
    assert py['result']==cpp['result']
    for a,b in zip(py['levels'],cpp['levels']):
        for key in ['modulus','prefixes','liftable','children','prefix_fnv1a64']:
            assert a[key]==b[key],(key,a,b)
    retained=json.loads((ROOT/'results/astra_second_2026_09_04/coefficient_obstruction.json').read_text())
    assert py['levels']==retained['levels']
    assert py['jacobian_divided_by_3_mod_3']==retained['jacobian_divided_by_3_mod_3']
    print(json.dumps({'result':'PASS','exact_symbolic_identities':'PASS',
        'formal_fraction_type_guards':'PASS','formal_jet_orders':[72,30],
        'linear_lift_controls':24,'independent_implementations':['Python','C++17'],
        'prefix_checks_per_implementation':py['prefix_checks_including_initial'],
        'obstruction_modulus':729,'retained_certificate_match':True},indent=2))


if __name__=='__main__':main()
