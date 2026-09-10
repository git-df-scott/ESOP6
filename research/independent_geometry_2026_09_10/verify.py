#!/usr/bin/env python3
"""Exact finite checks supporting the 7-adic centrally symmetric conic obstruction.
The complete proof, including passage from finite residues to Q_7, is in REPORT.md.
"""
import itertools, json
from pathlib import Path
p=7
sixth = {pow(x,6,p) for x in range(p)}
cubes = {pow(x,3,p) for x in range(p)}
squares = {pow(x,2,p) for x in range(1,p)}
assert sixth == {0,1}
assert cubes == {0,1,6}
assert squares == {1,2,4}
counts = sorted({sum(pow(x,6,p) for x in row)%p
                 for row in itertools.product(range(p),repeat=5)
                 if any(row)})
assert counts == [1,2,3,4,5]
assert set(counts) & cubes == {1}
assert {x for x in range(1,p) if pow(x,3,p)==1} == squares
assert (1**2+2**2)%7 == 5 and 5 not in squares
# Independent complete residue identity check for each nonzero (A,B).
# A nondegenerate diagonal quadratic cannot take only zero/square residues.
# This is the unit-coefficient case of the diagonalization argument.
bad_diagonals=[]
for A,B in itertools.product(range(1,p), repeat=2):
    values={(A*u*u+B*v*v)%p for u,v in itertools.product(range(p),repeat=2)}
    if values <= squares|{0}: bad_diagonals.append([A,B])
assert not bad_diagonals
result={"sixth_residues": sorted(sixth), "cube_residues": sorted(cubes),
        "unit_square_residues": sorted(squares),
        "primitive_five_term_sum_residues": counts,
        "cube_intersection": sorted(set(counts)&cubes),
        "ordered_residue_tuples_checked": 7**5,
        "nonzero_diagonal_pairs_checked":36,
        "diagonal_pairs_taking_only_square_or_zero_values":bad_diagonals,
        "counterexample_found":False,
        "scope":"Finite checks support the proof; not a search for integer ESOP6 points."}
path=Path(__file__).with_name('verification.json')
path.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
