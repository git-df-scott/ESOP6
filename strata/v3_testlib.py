#!/usr/bin/env python3
"""Helpers for strata/test_v3.sh: case generation, brute-force two-sum oracle,
pattern classification.  Pure Python big integers; knows nothing about the
engine's internals."""
import random, sys, os

def sieve(n):
    s = bytearray([1]) * (n + 1)
    s[0:2] = b"\0\0"
    for i in range(2, int(n ** .5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(range(i * i, n + 1, i)))
    return s

def primes_upto(n):
    s = sieve(n)
    return [i for i in range(2, n + 1) if s[i]]

_CACHE = {}
def _pows(B):
    if B not in _CACHE:
        _CACHE[B] = (set(x ** 6 for x in range(1, B + 1)), [x ** 6 for x in range(0, B + 1)])
    return _CACHE[B]

def two_sum_oracle(R, B):
    """Brute force: is R = x^6 + y^6 with 1 <= y <= x <= B?  Returns (x,y) or None."""
    if R < 2:
        return None
    S, P = _pows(B)
    x = int(round(R ** (1.0 / 6)))
    while x ** 6 > R:
        x -= 1
    while (x + 1) ** 6 <= R:
        x += 1
    if x > B:
        x = B
    while x >= 1:
        p = P[x]
        if 2 * p < R:
            break
        rest = R - p
        if rest >= 1 and rest in S:
            y = int(round(rest ** (1.0 / 6)))
            while y > 0 and y ** 6 > rest:
                y -= 1
            while (y + 1) ** 6 <= rest:
                y += 1
            if y >= 1 and y <= x and y ** 6 == rest:
                return (x, y)
        x -= 1
    return None

def pattern(a, b, pcut):
    """Classify the divisibility pattern of the pair (a,b) w.r.t. primes <= pcut."""
    def fac(n):
        f = {}
        d = 2
        while d * d <= n:
            while n % d == 0:
                f[d] = f.get(d, 0) + 1
                n //= d
            d += 1
        if n > 1:
            f[n] = f.get(n, 0) + 1
        return f
    fa, fb = fac(a), fac(b)
    common = [p for p in fa if p in fb and p <= pcut]
    tags = []
    if common:
        depth = max(min(fa[p], fb[p]) for p in common)
        tags.append("both-div-p")
        if depth >= 2 or len(common) >= 2:
            tags.append("recursion-depth>=2")
        else:
            tags.append("recursion-depth=1")
    one = [p for p in set(list(fa) + list(fb)) if p <= pcut and (p in fa) != (p in fb)]
    if one:
        tags.append("exactly-one-div-p")
        if any(p ** (6 * (fa.get(p, 0) + fb.get(p, 0))) > 0 for p in one):
            pass
    if all(p > pcut for p in fa) and all(p > pcut for p in fb):
        tags.append("rough-rough")
    if (a == 1 or (len(fa) == 1 and list(fa.values())[0] == 1 and list(fa)[0] > pcut)) and \
       (b == 1 or (len(fb) == 1 and list(fb.values())[0] == 1 and list(fb)[0] > pcut)):
        tags.append("spec-ROUGH-pair")
    return tags

def gen_cases(B, npos, nneg, seed, pcut):
    """Return (cases, expected, tags) where cases is a list of (R, bound)."""
    random.seed(seed)
    P = primes_upto(B)
    big = [p for p in P if p > 700] or [p for p in P if p > pcut]
    small = [p for p in P if p <= 60]
    cases, exp, tags = [], [], []

    def addpos(a, b, note):
        if a < 1 or b < 1 or a > B or b > B:
            return
        cases.append((a ** 6 + b ** 6, B)); exp.append(1)
        tags.append(pattern(a, b, pcut) + [note])

    # structured families: every divisibility pattern the engine distinguishes
    fams = 0
    while fams < npos // 2:
        fams += 1
        k = fams % 10
        if k == 0 and big:                      # rough x rough (both primes > 700)
            addpos(random.choice(big), random.choice(big), "rough-rough")
        elif k == 1 and big:                    # e*q x q' : the class SPEC_V3 misses
            e = random.choice([2, 3, 4, 5, 6, 7, 8, 10, 12])
            addpos(e * random.choice(big), random.choice(big), "smooth*rough")
        elif k == 2 and big:                    # common large prime
            q = random.choice(big)
            addpos(random.choice([1, 2, 3]) * q, random.choice([1, 2, 3]) * q, "common-large-prime")
        elif k == 3:                            # both divisible by p, depth >= 2
            p = random.choice([2, 3, 5, 7])
            g = p * p
            addpos(g * random.randint(1, B // g), g * random.randint(1, B // g), "both-div-p^2")
        elif k == 4:                            # both divisible by two distinct primes
            p, q = random.sample([2, 3, 5, 7, 11, 13], 2)
            g = p * q
            if B // g >= 1:
                addpos(g * random.randint(1, B // g), g * random.randint(1, B // g), "both-div-pq")
        elif k == 5:                            # exactly one divisible by a mid prime
            p = random.choice([11, 13, 17, 19, 23, 29, 31, 37, 41, 43])
            a = p * random.randint(1, max(1, B // p))
            b = random.randint(1, B)
            while b % p == 0:
                b = random.randint(1, B)
            addpos(a, b, "one-div-midprime")
        elif k == 6:                            # exactly one even / one div by 3
            a = 2 * random.randint(1, B // 2)
            b = random.randrange(1, B, 2)
            addpos(a, b, "one-even")
        elif k == 7:                            # small smooth bases
            addpos(random.randint(1, 60), random.randint(1, 60), "small")
        elif k == 8:                            # one base equal to 1
            addpos(1, random.randint(1, B), "one-is-1")
        else:
            addpos(random.randint(1, B), random.randint(1, B), "uniform")
    while len(cases) < npos:
        addpos(random.randint(1, B), random.randint(1, B), "uniform")

    # negatives: perturbed sums and random 128-bit values
    tried = 0
    while len(cases) < npos + nneg and tried < 20 * nneg:
        tried += 1
        if random.random() < 0.7:
            a = random.randint(1, B); b = random.randint(1, B)
            R = a ** 6 + b ** 6 + random.choice([1, 2, -1, 42 ** 6, -(42 ** 6), 3])
        else:
            R = random.getrandbits(random.choice([40, 60, 80, 100, 120]))
        if R < 2:
            continue
        if two_sum_oracle(R, B) is not None:
            continue
        cases.append((R, B)); exp.append(0); tags.append(["negative"])
    return cases, exp, tags

if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "gen":
        B, npos, nneg, seed, pcut, out = (int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]),
                                          int(sys.argv[5]), int(sys.argv[6]), sys.argv[7])
        cases, exp, tags = gen_cases(B, npos, nneg, seed, pcut)
        with open(out, "w") as f:
            for R, b in cases:
                f.write("%d %d\n" % (R, b))
        with open(out + ".exp", "w") as f:
            for e, t in zip(exp, tags):
                f.write("%d %s\n" % (e, ",".join(t)))
        print("generated %d cases (%d positive, %d negative) into %s" %
              (len(cases), sum(exp), len(exp) - sum(exp), out))
    elif cmd == "check":
        casef, outf = sys.argv[2], sys.argv[3]
        cases = [l.split() for l in open(casef)]
        expl = [l.split(None, 1) for l in open(casef + ".exp")]
        got = [l.split() for l in open(outf)]
        assert len(got) == len(cases), "engine produced %d lines for %d cases" % (len(got), len(cases))
        bad = 0
        from collections import Counter
        cov = Counter()
        for i, (c, e, g) in enumerate(zip(cases, expl, got)):
            R = int(c[0]); B = int(c[1]); want = int(e[0])
            for t in e[1].strip().split(","):
                cov[t] += 1
            yes = (g[1] == "YES")
            if yes:
                x, y = int(g[2]), int(g[3])
                if x ** 6 + y ** 6 != R or not (1 <= x <= B and 1 <= y <= B):
                    print("BAD PAIR at %d: R=%d pair=(%d,%d)" % (i, R, x, y)); bad += 1; continue
            if int(yes) != want:
                bad += 1
                if bad <= 10:
                    print("MISMATCH at %d: R=%d B=%d expected=%d got=%s (%s)" % (i, R, B, want, g[1], e[1].strip()))
        print("pattern coverage:")
        for k, v in sorted(cov.items()):
            print("    %-22s %d" % (k, v))
        print("checked %d cases, %d mismatches" % (len(cases), bad))
        sys.exit(1 if bad else 0)
    elif cmd == "candcount":
        # independent Python candidate count for the band (fmin,fmax]
        fmin, fmax = int(sys.argv[2]), int(sys.argv[3])
        M = 42 ** 6
        roots = [x for x in range(1, M) if pow(x, 6, M) == 1] if M < 10 ** 6 else None
        # build the 144 roots by CRT (M is too large to scan)
        r2 = [x for x in range(1, 64, 2) if pow(x, 6, 64) == 1]
        r3 = [x for x in range(1, 729) if x % 3 and pow(x, 6, 729) == 1]
        r7 = [x for x in range(1, 117649) if x % 7 and pow(x, 6, 117649) == 1]
        roots = []
        for a in r2:
            for b in r3:
                for c in r7:
                    x = a
                    while x % 729 != b:
                        x += 64
                    while x % 117649 != c:
                        x += 64 * 729
                    roots.append(x)
        assert len(roots) == 144
        n = 0
        for f in range(fmin + 1, fmax + 1):
            if f % 2 == 0 or f % 3 == 0 or f % 7 == 0:
                continue
            for z in roots:
                t = z * f % M
                if 0 < t < f:
                    n += 1
        print(n)
    else:
        sys.exit("unknown command " + cmd)
