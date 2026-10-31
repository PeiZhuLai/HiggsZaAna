#!/usr/bin/env python3
"""Completeness and reconciliation check of the SF-uncertainty signal re-production.

For every <mA>_<era> directory of the new production:
  * every job_N has a summary json with "successful": true (completion marker = summary, not parquet);
  * for each of the 17 variations: sum of chunk rows (output_job_*_<syst>.parquet, read from
    metadata only) == rows of merged_<syst>.parquet;
  * merged_<syst>.parquet is not older than any of its chunks (merge-before-rechunk race);
  * merged_nominal rows equal the previous production (same selection, same inputs).
Exit code 1 on any failure.

usage: recon_sig_sfunc_20261024.py [new_Sig_MC_dir]
"""
import glob
import json
import os
import sys

import pyarrow.parquet as pq

NEW = sys.argv[1] if len(sys.argv) > 1 else "/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_sfunc/Sig_MC"
OLD = "/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix/Sig_MC"
MASSES = ["mA_M%d" % m for m in (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30)]
ERAS = os.environ.get("RECON_ERAS", "2022preEE,2022postEE,2023preBPix,2023postBPix,2024").split(",")
SYSTS = ["nominal"] + ["%s_%s" % (s, ud) for s in ("FNUF", "Material", "Electron_scale", "Electron_smear",
                                                    "Muon_scale", "Muon_smear", "Photon_scale", "Photon_smear")
                       for ud in ("up", "down")]


def rows(path):
    return pq.ParquetFile(path).metadata.num_rows


def main():
    bad = 0
    for m in MASSES:
        for e in ERAS:
            sy = "%s_%s" % (m, e)
            d = os.path.join(NEW, sy)
            if not os.path.isdir(d):
                print("!! %s: no output dir" % sy)
                bad += 1
                continue
            jobs = sorted(glob.glob(os.path.join(d, "job_*")))
            summ = [glob.glob(os.path.join(j, "*_summary_job*.json")) for j in jobs]
            ok = sum(1 for s in summ if s and json.load(open(s[0])).get("successful") is True)
            msg = []
            n_nom = -1
            if ok != len(jobs):
                msg.append("summaries ok %d/%d" % (ok, len(jobs)))
            for syst in SYSTS:
                merged = os.path.join(d, "merged_%s.parquet" % syst)
                chunks = glob.glob(os.path.join(d, "job_*", "output_job_*_%s.parquet" % syst))
                if not os.path.exists(merged):
                    msg.append("missing merged_%s" % syst)
                    continue
                n_chunk = sum(rows(c) for c in chunks)
                n_merged = rows(merged)
                if n_chunk != n_merged:
                    msg.append("%s chunk %d != merged %d" % (syst, n_chunk, n_merged))
                if chunks and os.path.getmtime(merged) < max(os.path.getmtime(c) for c in chunks):
                    msg.append("%s merged older than a chunk" % syst)
                if syst == "nominal":
                    n_nom = n_merged
                    n_old = rows(os.path.join(OLD, sy, "merged_nominal.parquet"))
                    if n_old != n_merged:
                        msg.append("nominal rows new %d != old %d" % (n_merged, n_old))
            status = "FAIL" if msg else "ok"
            bad += bool(msg)
            print("[%s] %-22s jobs %3d  nominal rows %6d  %s" % (status, sy, len(jobs), n_nom, "; ".join(msg)))
    print("SUMMARY bad=%d" % bad)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
