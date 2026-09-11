#!/usr/bin/env python3
"""Rewrite the COVERAGE section of LOG716.md from runs/coverage.txt.

Called by run716.sh after every completed chunk, so LOG716.md always reflects
the last chunk that was fully covered by ALL Q residue passes.
"""
import os, re
d = os.path.dirname(os.path.abspath(__file__))
rows = []
for ln in open(os.path.join(d, "runs/coverage.txt")):
    p = ln.split()
    if not p or p[0] != "COVERED":
        continue
    rows.append((int(p[1]), int(p[2]),
                 int(p[3].split("=")[1]), int(p[4].split("=")[1]),
                 int(p[5].split("=")[1].rstrip("s")), p[6].split("=")[1]))
rows.sort()
out = ["<!--COVERAGE-->", "## COVERAGE (auto-generated from runs/coverage.txt; do not hand-edit)", ""]
out.append("Prior session, single pass, same engine and same host:")
out.append("")
out.append("    (7,1,6): 2 < f <= 1500, 2849012361 leaves, 0 solutions, 4-core host, 238 s")
out.append("")
if rows:
    out.append("This session, Q=7 residue-class bucketing, chunked so that each line below")
    out.append("is fully covered by all 7 passes:")
    out.append("")
    for lo, hi, lv, sl, w, q in rows:
        out.append("    (7,1,6): %d < f <= %d, %d leaves, %d solutions, 4-core host, %d s"
                   % (lo, hi, lv, sl, w))
    # contiguous merge upward from the lowest lo
    lo0, hi0 = rows[0][0], rows[0][1]
    tot, wall = rows[0][2], rows[0][4]
    for lo, hi, lv, sl, w, q in rows[1:]:
        if lo != hi0:
            break
        hi0, tot, wall = hi, tot + lv, wall + w
    out.append("")
    out.append("Aggregate for this session (contiguous):")
    out.append("")
    out.append("    (7,1,6): %d < f <= %d, %d leaves, %d solutions, 4-core host, %d s"
               % (lo0, hi0, tot, sum(r[3] for r in rows), wall))
    out.append("")
    out.append("Total coverage to date (contiguous from f=2):")
    out.append("")
    out.append("    (7,1,6): 2 < f <= %d, %d leaves, 0 solutions, 4-core host, %d s"
               % (hi0, 2849012361 + tot, 238 + wall))
out.append("")
out.append("The EulerNet search bound for k=7 is **unknown to us** (euler.free.fr is")
out.append("unreachable), so no claim of novel coverage is made -- these lines record")
out.append("only what this run itself covered.")
out.append("<!--/COVERAGE-->")
block = "\n".join(out) + "\n"
p = os.path.join(d, "LOG716.md")
s = open(p).read()
if "<!--COVERAGE-->" in s:
    s = re.sub(r"<!--COVERAGE-->.*?<!--/COVERAGE-->\n?", block, s, flags=re.S)
else:
    s = s.rstrip() + "\n\n" + block
open(p, "w").write(s)
