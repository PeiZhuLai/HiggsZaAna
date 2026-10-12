#!/usr/bin/env python3
"""Compare two R1 outputs: every parquet column bit by bit (NaN == NaN), schema, and summary numbers.
Usage: compare_r1.py <dirA> <dirB>"""
import glob, json, os, sys
import numpy as np
import pyarrow.parquet as pq
a, b = sys.argv[1:3]
fa = sorted(os.path.relpath(p, a) for p in glob.glob(a + "/**/*.parquet", recursive=True))
fb = sorted(os.path.relpath(p, b) for p in glob.glob(b + "/**/*.parquet", recursive=True))
bad = 0
if fa != fb:
    print("file lists differ:\n  only A:", sorted(set(fa) - set(fb)), "\n  only B:", sorted(set(fb) - set(fa))); bad += 1
for rel in sorted(set(fa) & set(fb)):
    ta, tb = pq.read_table(os.path.join(a, rel)), pq.read_table(os.path.join(b, rel))
    if ta.schema != tb.schema:
        print("schema differs:", rel, set(ta.column_names) ^ set(tb.column_names)); bad += 1; continue
    if ta.num_rows != tb.num_rows:
        print("rows differ:", rel, ta.num_rows, tb.num_rows); bad += 1; continue
    for col in ta.column_names:
        x, y = ta[col].to_numpy(zero_copy_only=False), tb[col].to_numpy(zero_copy_only=False)
        try:
            same = np.array_equal(x, y, equal_nan=True)
        except TypeError:
            same = all(np.array_equal(np.asarray(u), np.asarray(v)) for u, v in zip(x, y))
        if not same:
            print("column differs:", rel, col); bad += 1
def norm(d, root):
    # timing fields and output paths legitimately differ between runs
    d = {k: v for k, v in d.items() if not k.startswith("time") and k not in ("config", "outputs")}
    return json.loads(json.dumps(d).replace(root, "<OUT>"))
sa, sb = (norm(json.load(open(os.path.join(d, "summary.json"))), d) for d in (a, b))
if sa != sb:
    print("summary differs:", sorted(k for k in set(sa) | set(sb) if sa.get(k) != sb.get(k))); bad += 1
print("parquet files compared:", len(set(fa) & set(fb)), "| differences:", bad)
sys.exit(1 if bad else 0)
