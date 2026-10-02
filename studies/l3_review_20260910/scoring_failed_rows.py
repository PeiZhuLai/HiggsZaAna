#!/usr/bin/env python3
"""Write the joblist rows whose scored ROOT fails GATE 5 (same check as gate5_scored.py).

Usage: scoring_failed_rows.py <joblist.tsv> <since_epoch> <out.tsv>
Prints the number of failed rows; exit 0 always (the caller decides).
"""
import sys
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910")
from gate5_scored import check  # noqa: E402

jl, since, out = sys.argv[1], float(sys.argv[2]), sys.argv[3]
rows = [l.rstrip("\n").split("\t") for l in open(jl) if l.strip()]
with ThreadPoolExecutor(16) as ex:
    res = list(ex.map(lambda r: check(r, since), rows))
bad = [r for r, x in zip(rows, res) if x[1] != "OK"]
with open(out, "w") as f:
    for r in bad:
        f.write("\t".join(r) + "\n")
print(len(bad))
