#!/usr/bin/env python3
"""Replay finite local certificates without importing the generating code."""
import json
from math import isqrt
from pathlib import Path


def main():
    root = Path(__file__).parent
    result = json.loads((root / "exact_results.json").read_text())
    cert = result["ternary_fermat_local_certificate"]
    bad = []
    for p in range(5, 398):
        if any(p % d == 0 for d in range(2, isqrt(p)+1)):
            continue
        # Direct variable enumeration, not a sumset-based projective search.
        has = any((pow(x, 6, p)+pow(y, 6, p)+1) % p == 0
                  for x in range(p) for y in range(p))
        if not has:
            bad.append(p)
    assert bad == cert["anisotropic_primes_excluding_2_3"]
    for p, modulus, unit_modulus in [(2, 256, 8), (3, 729, 9)]:
        sixth = {pow(a, 6, modulus) for a in range(modulus) if a % p}
        expected = {a for a in range(modulus) if a % unit_modulus == 1}
        assert sixth == expected
        assert all((3*k*k-2) % modulus in sixth
                   for k in range(modulus) if k % p)
    # Check each retained Pythagorean triple and every claimed finite gate.
    survivors = []
    for row in result["pythagorean_enumeration"]["survivor_ledger"]:
        a, b, c = row["pythagorean_triple"]
        assert a*a+b*b == c*c and a*b*c == row["K"]
        statuses = []
        for p in [2, 3]+bad:
            unit, v = row["K"], 0
            while unit % p == 0:
                v += 1
                unit //= p
            allowed = v % 3 == 0
            if p not in (2, 3):
                residues = [pow(x, 6, p) for x in range(p)]
                pair = {(x+y) % p for x in residues for y in residues}
                allowed = allowed and any((3*unit*unit-z) % p in pair for z in residues)
            assert allowed == row["local_data"][str(p)]["locally_soluble"]
            statuses.append(allowed)
        assert all(statuses) == row["all_places_locally_soluble"]
        if all(statuses):
            survivors.append([row["m"], row["n"]])
    output = {"status": "PASS", "direct_prime_enumeration_bound": 397,
              "sixth_unit_moduli": [256, 729], "retained_all_place_survivors": survivors,
              "scope": "independent finite certificates; all-place theorem also uses proof and Weil bound in REPORT.md"}
    print(json.dumps(output, indent=2))
    (root / "independent_finite_verification.json").write_text(json.dumps(output, indent=2)+"\n")


if __name__ == "__main__":
    main()
