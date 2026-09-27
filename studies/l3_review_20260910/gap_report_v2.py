#!/usr/bin/env python
"""Real completeness for the fpo=1 production.

The completion marker is the SUMMARY json, not the parquet. At fpo=1 a job
covers one input file and "zero events selected" is routine: HiggsDNA then
writes the summary with n_events_selected=0 and NO parquet, and the job is
successful. Counting parquets called 2286 healthy jobs "missing" on 2026-09-16
and sent 200 needless resubmissions. (ref_progress_count_needs_completion_marker)
"""
import os, sys, json, glob
B = "/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix_fpo1"
OUT = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/logs_fsrfix/missing_jobs_v2.txt"
open(OUT, "w").close()
tj = ts = tp = tz = 0
for stage in ("Bkg_MC", "Data"):
    base = os.path.join(B, stage)
    if not os.path.isdir(base): continue
    print("--- %s" % stage)
    for smp in sorted(os.listdir(base)):
        d = os.path.join(base, smp)
        if not os.path.isdir(d): continue
        jobs = sorted(x for x in os.listdir(d) if x.startswith("job_"))
        nosum, npq, nzero = [], 0, 0
        for j in jobs:
            n = j[4:]
            jd = os.path.join(d, j)
            s = glob.glob(os.path.join(jd, "*_summary_job%s.json" % n))
            if not s:
                nosum.append(j); continue
            if os.path.exists(os.path.join(jd, "output_job_%s_nominal.parquet" % n)):
                npq += 1
            else:
                nzero += 1
        tj += len(jobs); ts += len(jobs) - len(nosum); tp += npq; tz += nzero
        flag = "" if not nosum else "   NOT RUN %d" % len(nosum)
        print("  %-34s jobs=%6d ran=%6d parquet=%6d zero-sel=%5d%s"
              % (smp, len(jobs), len(jobs) - len(nosum), npq, nzero, flag))
        if nosum:
            with open(OUT, "a") as f:
                for j in nosum: f.write("%s %s %s\n" % (stage, smp, j))
print("TOTAL jobs=%d ran=%d parquet=%d zero-selection=%d NOT-RUN=%d"
      % (tj, ts, tp, tz, tj - ts))
