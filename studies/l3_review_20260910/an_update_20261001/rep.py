"""Exact-string replacement helper for the AN update. usage: rep.py <tex> <pairs.py>
pairs.py defines PAIRS = [(old, new), ...]; each old must occur exactly once."""
import sys, runpy
tex, pf = sys.argv[1], sys.argv[2]
P = runpy.run_path(pf)["PAIRS"]
s = open(tex).read()
for i, (o, n) in enumerate(P):
    c = s.count(o)
    if c != 1:
        sys.exit(f"[rep] {tex}: pair {i} found {c} times: {o[:80]!r}")
    s = s.replace(o, n)
open(tex, "w").write(s)
print(f"[rep] {tex}: {len(P)} replacements")
