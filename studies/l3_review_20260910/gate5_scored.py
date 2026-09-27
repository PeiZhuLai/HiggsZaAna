#!/usr/bin/env python
"""
GATE 5 -- validate BDT-scored ROOT files.

Two things must hold, and the second is the one that a file count or an exit code
cannot see:
  1. entries match the merged parquet the job read (nothing dropped in conversion);
  2. the MVA_Score_mA_M* branches exist AND are not all-NaN.

(2) matters because add_mva_scores() catches a model-loading failure and fills every
score with NaN, prints a WARNING, and returns normally. The job then exits 0 and writes
a perfectly valid ROOT file with no usable scores in it. Scoring with the wrong or a
missing model looks exactly like success everywhere except in the score values.

Usage: python gate5_scored.py [--joblist PATH] [--since EPOCH]
"""
import os, sys, time
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pyarrow.parquet as pq
import uproot

DEF_JOBLIST = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Parquet2Rootfile/Condor/joblist.tsv"
REPORT = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/logs_fsrfix/gate5_report.txt"
PROBE_MASSES = [1, 2, 3, 4, 10, 30]   # 1-3 -> low-mass model, >=4 -> high-mass model


def check(row, since):
    pqf, rootf = row[0], row[1]
    label = "/".join(rootf.rsplit("/", 2)[-2:])
    if not os.path.exists(rootf):
        return (label, "MISSING_ROOT", None, None, "")
    mt = os.path.getmtime(rootf)
    if since and mt < since:
        return (label, "STALE_ROOT", None, None,
                time.strftime('%F %T', time.localtime(mt)))
    try:
        rows = pq.ParquetFile(pqf).metadata.num_rows
    except Exception as e:
        return (label, "PARQUET_UNREADABLE", None, None, str(e)[:70])
    try:
        with uproot.open(rootf) as f:
            t = f["inclusive"]
            n = t.num_entries
            missing, allnan = [], []
            for m in PROBE_MASSES:
                b = "MVA_Score_mA_M%d" % m
                if b not in t:
                    missing.append(b); continue
                a = t[b].array(library="np")
                if a.size and not np.isfinite(a).any():
                    allnan.append(b)
    except Exception as e:
        return (label, "ROOT_UNREADABLE", rows, None, str(e)[:70])
    if missing:
        return (label, "NO_SCORE_BRANCH", rows, n, ",".join(missing)[:60])
    if allnan:
        return (label, "SCORES_ALL_NAN", rows, n, ",".join(allnan)[:60])
    if rows != n:
        return (label, "ROW_MISMATCH", rows, n, "diff %+d" % (n - rows))
    return (label, "OK", rows, n, "")


def main():
    jl, since = DEF_JOBLIST, 0.0
    if "--joblist" in sys.argv:
        jl = sys.argv[sys.argv.index("--joblist") + 1]
    if "--since" in sys.argv:
        since = float(sys.argv[sys.argv.index("--since") + 1])
    rows = [l.rstrip("\n").split("\t") for l in open(jl) if l.strip()]
    with ThreadPoolExecutor(max_workers=6) as ex:
        res = list(ex.map(lambda r: check(r, since), rows))
    bad = [r for r in res if r[1] != "OK"]
    L = ["GATE 5 -- scored ROOT: entries match parquet AND MVA scores are real",
         "generated %s" % time.strftime('%F %T'),
         "joblist: %s" % jl,
         "checked: %d   OK: %d   FAILED: %d" % (len(res), len(res) - len(bad), len(bad)), ""]
    if bad:
        by = {}
        for r in bad:
            by.setdefault(r[1], []).append(r)
        for st in sorted(by):
            L.append("  %s : %d" % (st, len(by[st])))
            for r in sorted(by[st])[:30]:
                L.append("      %-42s %10s %10s  %s" % (r[0], r[2], r[3], r[4]))
            if len(by[st]) > 30:
                L.append("      ... and %d more" % (len(by[st]) - 30))
        L += ["", "VERDICT: FAILED"]
    else:
        L += ["total inclusive entries: %d" % sum(r[3] for r in res), "VERDICT: PASSED"]
    txt = "\n".join(L)
    open(REPORT, "w").write(txt + "\n")
    print(txt)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
