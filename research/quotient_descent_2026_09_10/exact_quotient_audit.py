#!/usr/bin/env python3
"""Exact quotient identities, seed-slice obstructions, and Pythagorean local ledger.

Standard library only. No floating point, rational recognition, or point-search
claim. The all-place local theorem used by the finite classifier is in REPORT.md.
"""
import argparse
from collections import Counter
from hashlib import sha256
import json
from math import gcd, isqrt
from pathlib import Path


def valuation(n, p):
    assert n
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v, n


def primes(bound):
    return [p for p in range(2, bound + 1)
            if all(p % d for d in range(2, isqrt(p) + 1))]


# Sparse polynomial arithmetic in four formal variables.
ONE = {(0, 0, 0, 0): 1}


def variable(i):
    e = [0] * 4
    e[i] = 1
    return {tuple(e): 1}


def add(*polys):
    out = Counter()
    for poly in polys:
        for e, c in poly.items():
            out[e] += c
    return {e: c for e, c in out.items() if c}


def scale(poly, c):
    return {e: a*c for e, a in poly.items() if a*c}


def mul(a, b):
    out = Counter()
    for e, x in a.items():
        for f, y in b.items():
            out[tuple(u+v for u, v in zip(e, f))] += x*y
    return {e: c for e, c in out.items() if c}


def power(a, n):
    out = ONE
    for _ in range(n):
        out = mul(out, a)
    return out


def phi(a, om):
    return add(power(a, 6), scale(mul(power(a, 4), om), 15),
               scale(mul(power(a, 2), power(om, 2)), 15), power(om, 3))


def symbolic_checks():
    a, b, c, d = [variable(i) for i in range(4)]
    e1, e2 = add(a, b), mul(a, b)
    quotient = add(power(e1, 6), scale(mul(power(e1, 4), e2), -6),
                   scale(mul(power(e1, 2), power(e2, 2)), 9),
                   scale(power(e2, 3), -2))
    assert quotient == add(power(a, 6), power(b, 6))
    assert add(power(add(a, b), 6), power(add(a, scale(b, -1)), 6)) == scale(phi(a, power(b, 2)), 2)
    # T = Omega + 5A^2; Phi = T^3 - 60A^4 T + 176A^6.
    t = add(b, scale(power(a, 2), 5))
    assert phi(a, b) == add(power(t, 3), scale(mul(power(a, 4), t), -60), scale(power(a, 6), 176))
    # If H^6 + 2 Phi(A,Omega) = C then x=-2T,y=2H^3
    # has y^2=x^3-240 A^4 x+4C-1408 A^6.
    x = scale(t, -2)
    ec_rhs = add(power(x, 3), scale(mul(power(a, 4), x), -240),
                 scale(c, 4), scale(power(a, 6), -1408))
    assert ec_rhs == scale(add(c, scale(phi(a, b), -2)), 4)
    # Pythagorean compression with m=a,n=b.
    leg1 = scale(mul(a, b), 2)
    leg2 = add(power(a, 2), scale(power(b, 2), -1))
    hyp = add(power(a, 2), power(b, 2))
    assert add(power(hyp, 6), scale(power(leg1, 6), -1), scale(power(leg2, 6), -1)) == scale(power(mul(mul(leg1, leg2), hyp), 2), 3)
    return {name: "PASS" for name in ["quotient_symmetric_identity", "conjugate_pair_identity", "depressed_cubic_identity", "elliptic_compression_identity", "pythagorean_compression_identity"]}


def quadratic_field_power(a, b, disc, exponent):
    x, y = 1, 0
    for _ in range(exponent):
        x, y = x*a+y*b*disc, x*b+y*a
    return x, y


def seed_audit():
    seeds = [("gaussian", 7, -121, [8, 12, 15], 17),
             ("boundary", 9, -249, [14, 18, 0], 22)]
    output = []
    for name, a, om, tail, f in seeds:
        plus = quadratic_field_power(a, 1, om, 6)
        minus = quadratic_field_power(a, -1, om, 6)
        lhs = plus[0]+minus[0]+sum(t**6 for t in tail)
        assert plus[1]+minus[1] == 0 and lhs == f**6
        slices = []
        for free_index, free_value in enumerate(tail):
            fixed = [v for i, v in enumerate(tail) if i != free_index]
            rhs = f**6-sum(v**6 for v in fixed)
            v7, u7 = valuation(rhs, 7)
            v3, u3 = valuation(rhs, 3)
            if v7 % 6 or u7 % 7 not in (1, 2, 3):
                obstruction = {"p": 7, "valuation": v7, "unit_residue": u7 % 7,
                               "allowed_valuation_mod_6": [0], "allowed_unit_residues": [1, 2, 3]}
            else:
                assert v3 % 6 not in (0, 1)
                obstruction = {"p": 3, "valuation": v3, "unit_residue": u3 % 3,
                               "allowed_valuation_mod_6": [0, 1]}
            slices.append({"free_tail_index": free_index+3, "seed_free_value": free_value,
                           "fixed_tail_values": fixed, "rhs_sum_three_sixths": rhs,
                           "obstruction": obstruction})
        output.append({"name": name, "A": a, "Omega": om, "tail": tail, "f": f,
                       "exact_algebraic_identity_lhs": lhs, "rhs": f**6,
                       "e1": 2*a, "e2": a*a-om, "discriminant": 4*om,
                       "rational_counterexample": False, "coordinate_slices": slices})
    return output


def fermat_local_data():
    witnesses, anisotropic = {}, []
    for p in primes(397):
        if p in (2, 3):
            continue
        roots = {pow(a, 6, p): a for a in range(p)}
        witness = next(((roots[a], roots[b], 1)
                        for a in roots for b in roots if (a+b+1) % p == 0), None)
        if witness:
            assert sum(pow(a, 6, p) for a in witness) % p == 0
            witnesses[str(p)] = list(witness)
        else:
            anisotropic.append(p)
    assert anisotropic == [7, 31, 67, 79, 139, 223]
    sets = {}
    for p in anisotropic:
        r6 = {pow(a, 6, p) for a in range(p)}
        sets[p] = {(a+b+c) % p for a in r6 for b in r6 for c in r6}
    return witnesses, anisotropic, sets


def local_status(k, anisotropic, sumsets):
    data = {}
    for p in [2, 3]+anisotropic:
        v, unit = valuation(k, p)
        residue = 3*unit*unit % p
        soluble = v % 3 == 0 and (p in (2, 3) or residue in sumsets[p])
        data[str(p)] = {"v_K": v, "normalized_rhs_mod_p": residue, "locally_soluble": soluble}
    return data


def enumerate_pythagorean(bound, anisotropic, sumsets):
    counts = Counter()
    ledger = []
    dispositions = sha256()
    for m in range(2, bound+1):
        for n in range(1, m):
            if (m-n) % 2 == 0 or gcd(m, n) != 1:
                continue
            counts["primitive_parameters"] += 1
            leg1, leg2, hyp = 2*m*n, m*m-n*n, m*m+n*n
            k = leg1*leg2*hyp
            data = local_status(k, anisotropic, sumsets)
            first_reject = next((p for p, row in data.items() if not row["locally_soluble"]), None)
            counts["first_reject_"+str(first_reject)] += 1
            dispositions.update(f"{m},{n}:{first_reject}\n".encode())
            if all(data[str(p)]["locally_soluble"] for p in (2, 3, 7)):
                counts["survives_2_3_7"] += 1
                ledger.append({"m": m, "n": n, "pythagorean_triple": [leg1, leg2, hyp],
                               "K": k, "rhs_sum_three_sixths": 3*k*k,
                               "local_data": data, "all_places_locally_soluble": first_reject is None,
                               "rational_point_status": "NOT SEARCHED / UNKNOWN"})
    counts["all_places_locally_soluble"] = sum(row["all_places_locally_soluble"] for row in ledger)
    return {"domain": {"m_min": 2, "m_max": bound, "n_min": 1, "n_max": "m-1", "gcd_mn": 1, "opposite_parity": True},
            "counts": dict(counts), "disposition_sha256": dispositions.hexdigest(), "survivor_ledger": ledger}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bound", type=int, default=1000)
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("exact_results.json"))
    args = parser.parse_args()
    witnesses, anisotropic, sumsets = fermat_local_data()
    result = {"arithmetic": "Python arbitrary-precision integers; symbolic coefficient dictionaries",
              "symbolic_checks": symbolic_checks(), "seeds": seed_audit(),
              "ternary_fermat_local_certificate": {"finite_prime_bound": 397, "genus": 10,
                  "large_prime_argument": "p>=401: (p+1)^2>400p, hence p+1-20sqrt(p)>0 by Weil",
                  "anisotropic_primes_excluding_2_3": anisotropic, "isotropic_witnesses": witnesses,
                  "sum_three_sixth_residue_sets": {str(p): sorted(values) for p, values in sumsets.items()}},
              "pythagorean_enumeration": enumerate_pythagorean(args.bound, anisotropic, sumsets),
              "counterexample_status": "NONE; no global rational-point existence claim"}
    args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"symbolic_checks": result["symbolic_checks"], "seed_slice_obstructions": 6,
                      "anisotropic_primes_excluding_2_3": anisotropic,
                      "enumeration_counts": result["pythagorean_enumeration"]["counts"],
                      "output": str(args.output)}, indent=2))


if __name__ == "__main__":
    main()
