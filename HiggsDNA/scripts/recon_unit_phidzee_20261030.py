#!/usr/bin/env python3
"""Completeness and reconciliation check of one unit of the photon-ID (Z->ee) re-production.

For every task directory (<sample>_<era>) under UNIT_DIR:
  * every job_N has a summary json with "successful": true (completion marker = summary, not parquet;
    with fpo=1 a job whose events all fail the selection writes no parquet, which is normal);
  * for every merged_<syst>.parquet: sum of chunk rows (output_job_*_<syst>.parquet, metadata only,
    merged files excluded) == merged rows;
  * merged_<syst>.parquet is not older than any of its chunks (merge-before-rechunk race);
  * --expect a,b,c: these task names must exist (stage 0: a missing sample would pass every other check).
Exit code 1 on any failure.

usage: recon_unit_phidzee_20261030.py UNIT_DIR [--expect task1,task2,...]
"""
import glob
import json
import os
import sys

import pyarrow.parquet as pq


def rows(path):
    return pq.ParquetFile(path).metadata.num_rows


def main():
    unit = sys.argv[1]
    expect = []
    if "--expect" in sys.argv:
        expect = sys.argv[sys.argv.index("--expect") + 1].split(",")
    tasks = sorted(t for t in os.listdir(unit) if os.path.isdir(os.path.join(unit, t)) and glob.glob(os.path.join(unit, t, "job_*")))
    bad = 0
    for t in expect:
        if t not in tasks:
            print("[FAIL] %-34s no output dir" % t)
            bad += 1
    total = 0
    for t in tasks:
        d = os.path.join(unit, t)
        jobs = glob.glob(os.path.join(d, "job_*"))
        summ = [glob.glob(os.path.join(j, "*_summary_job*.json")) for j in jobs]
        ok = sum(1 for s in summ if s and json.load(open(s[0])).get("successful") is True)
        msg = []
        if ok != len(jobs):
            msg.append("summaries ok %d/%d" % (ok, len(jobs)))
        merged = sorted(glob.glob(os.path.join(d, "merged_*.parquet")))
        if not merged:
            msg.append("no merged parquet")
        n_nom = -1
        for m in merged:
            syst = os.path.basename(m)[len("merged_"):-len(".parquet")]
            chunks = glob.glob(os.path.join(d, "job_*", "output_job_*_%s.parquet" % syst))
            n_chunk = sum(rows(c) for c in chunks)
            n_merged = rows(m)
            if n_chunk != n_merged:
                msg.append("%s chunk %d != merged %d" % (syst, n_chunk, n_merged))
            if chunks and os.path.getmtime(m) < max(os.path.getmtime(c) for c in chunks):
                msg.append("%s merged older than a chunk" % syst)
            if syst == "nominal":
                n_nom = n_merged
        total += max(n_nom, 0)
        bad += bool(msg)
        print("[%s] %-34s jobs %5d  merged %2d  nominal rows %9d  %s" % ("FAIL" if msg else "ok", t, len(jobs), len(merged), n_nom, "; ".join(msg)))
    print("SUMMARY tasks=%d nominal_rows=%d bad=%d" % (len(tasks), total, bad))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
