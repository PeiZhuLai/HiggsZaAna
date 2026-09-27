#!/usr/bin/env python
"""
Build a resubmission joblist from the rows of the full joblist whose scored ROOT file is
missing, unreadable, short, or has all-NaN MVA scores.

Same acceptance test as GATE 5, so a row that survives here is a row GATE 5 will pass.
"""
import os, sys
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pyarrow.parquet as pq
import uproot

JOBLIST = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Parquet2Rootfile/Condor/joblist.tsv"
OUT     = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Parquet2Rootfile/Condor/joblist_missing.tsv"
PROBE   = [1, 4, 30]


def bad(row):
    pqf, rootf = row[0], row[1]
    if not os.path.exists(rootf):
        return True
    try:
        rows = pq.ParquetFile(pqf).metadata.num_rows
        with uproot.open(rootf) as f:
            t = f["inclusive"]
            if t.num_entries != rows:
                return True
            for m in PROBE:
                b = "MVA_Score_mA_M%d" % m
                if b not in t:
                    return True
                a = t[b].array(library="np")
                if a.size and not np.isfinite(a).any():
                    return True
    except Exception:
        return True
    return False


rows = [l.rstrip("\n").split("\t") for l in open(JOBLIST) if l.strip()]
with ThreadPoolExecutor(max_workers=6) as ex:
    flags = list(ex.map(bad, rows))
missing = [r for r, f in zip(rows, flags) if f]
with open(OUT, "w") as fh:
    for r in missing:
        fh.write("\t".join(r) + "\n")
print("total %d   good %d   to resubmit %d" % (len(rows), len(rows) - len(missing), len(missing)))
print("wrote", OUT)
