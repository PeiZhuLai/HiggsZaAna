#!/usr/bin/env python
"""Where are the missing fpo=1 chunks? Per sample: job dirs vs chunks on disk."""
import os, glob, sys
B = "/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix_fpo1"
tot_j = tot_c = 0
for stage in ("Bkg_MC", "Data"):
    base = os.path.join(B, stage)
    if not os.path.isdir(base):
        continue
    print("--- %s" % stage)
    for smp in sorted(os.listdir(base)):
        d = os.path.join(base, smp)
        if not os.path.isdir(d):
            continue
        jobs = [x for x in os.listdir(d) if x.startswith("job_")]
        miss = [j for j in jobs
                if not os.path.exists(os.path.join(d, j, "output_%s_nominal.parquet" % j))]
        tot_j += len(jobs); tot_c += len(jobs) - len(miss)
        pct = 100.0 * len(miss) / len(jobs) if jobs else 0
        flag = "" if not miss else "   MISSING %d (%.1f%%)" % (len(miss), pct)
        print("  %-34s jobs=%6d done=%6d%s" % (smp, len(jobs), len(jobs) - len(miss), flag))
        if miss:
            with open("/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/logs_fsrfix/missing_jobs.txt", "a") as f:
                for j in sorted(miss):
                    f.write("%s %s %s\n" % (stage, smp, j))
print("TOTAL jobs=%d done=%d missing=%d" % (tot_j, tot_c, tot_j - tot_c))
