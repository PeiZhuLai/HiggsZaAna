#!/usr/bin/env python3
"""Smoke check of the condor re-scoring: the condor output (packed conda env, staged input) must equal, branch by
branch and bit by bit, the output of the same file re-scored locally (EOS conda env, FUSE input).
Usage: compare_smoke.py <file A> <file B>"""
import sys
import numpy as np
import uproot
a, b = uproot.open(sys.argv[1]), uproot.open(sys.argv[2])
ta = sorted(k for k in a.keys(cycle=False) if a[k].classname == "TTree")
tb = sorted(k for k in b.keys(cycle=False) if b[k].classname == "TTree")
assert ta == tb, (ta, tb)
nb = 0; bad = []
for t in ta:
    assert a[t].num_entries == b[t].num_entries and set(a[t].keys()) == set(b[t].keys()), t
    for br in a[t].keys():
        x = a[t][br].array(library="np"); y = b[t][br].array(library="np"); nb += 1
        same = np.array_equal(x, y, equal_nan=True) if x.dtype.kind == "f" else (x.dtype == object or np.array_equal(x, y))
        if x.dtype == object:
            same = all(np.array_equal(u, v) for u, v in zip(x, y))
        if not same: bad.append("%s:%s" % (t, br))
print("trees %d branches compared %d differing %d %s" % (len(ta), nb, len(bad), bad[:5]))
sys.exit(1 if bad else 0)
