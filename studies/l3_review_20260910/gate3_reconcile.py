#!/usr/bin/env python
"""
GATE 3 for the FSR-fix production: merged parquet rows == ROOT `inclusive` entries,
for every one of the 1209 sample/era/systematic conversions.

Why this and not the exit code: 1_run_P2Root.sh retries a failed conversion twice,
prints "Sample X completed successfully" regardless, and exits 0. On 2026-09-20 all
1209 conversions died on an import error and it still reported RC=0. The only
trustworthy signal is the product.

Checks per entry:
  merged_rows == inclusive       (nothing dropped in conversion)
  train+validation+test == incl  (the --split partition is exhaustive, internal check)
  ROOT mtime is from this round  (not a leftover from an earlier production)

Reads only parquet footers and ROOT headers, so it is I/O-latency bound, not
throughput bound -- hence the small thread pool rather than a serial loop.

Usage:  python gate3_reconcile.py [--since EPOCH]
Output: logs_fsrfix/gate3_report.txt   (and a non-zero exit code if anything fails)
"""
import os, sys, time
from concurrent.futures import ThreadPoolExecutor

import pyarrow.parquet as pq
import uproot

SIG_IN  = "/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix/Sig_MC"
BKG_IN  = "/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix_fpo1"
OUT     = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_inputs_fsrfix"
REPORT  = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/logs_fsrfix/gate3_report.txt"

YEARS  = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]
SIGS   = ["mA_M%s" % m for m in [1,2,3,4,5,6,7,8,9,10,15,20,25,30]]
SYSTS  = ["FNUF","Material","Electron_scale","Electron_smear",
          "Muon_scale","Muon_smear","Photon_scale","Photon_smear"]
# 2026-09-21: 2022 DYGto2LG comes from the inclusive PTG-10to100 sample, not the
# retired 10to50 + 50to100 slices -- that is what the FSR-fix production actually
# contains. Expected conversions therefore drop from 1209 to 1207.
BKG    = ([("DYGto2LG_10to100", y) for y in YEARS] +
          [("DYJetsToLL", y) for y in ["2022preEE","2022postEE","2023preBPix","2023postBPix"]] +
          [("DYJetsTo2E","2024"), ("DYJetsTo2Mu","2024"), ("DYJetsTo2Tau","2024")])


def jobs():
    for s in SIGS:
        for y in YEARS:
            yield ("%s/%s" % (s, y),
                   "%s/%s_%s/merged_nominal.parquet" % (SIG_IN, s, y),
                   "%s/%s/%s.root" % (OUT, s, y))
            for sy in SYSTS:
                for ud in ("up", "down"):
                    yield ("%s_%s_%s/%s" % (s, sy, ud, y),
                           "%s/%s_%s/merged_%s_%s.parquet" % (SIG_IN, s, y, sy, ud),
                           "%s/%s_%s_%s/%s.root" % (OUT, s, sy, ud, y))
    for s, y in BKG:
        yield ("%s/%s" % (s, y),
               "%s/Bkg_MC/%s_%s/merged_nominal.parquet" % (BKG_IN, s, y),
               "%s/%s/%s.root" % (OUT, s, y))
    for y in YEARS:
        yield ("Data/%s" % y,
               "%s/Data/Data_%s/merged_nominal.parquet" % (BKG_IN, y),
               "%s/Data/%s.root" % (OUT, y))


def check(job, since):
    name, pqf, rf = job
    if not os.path.exists(rf):
        return (name, "MISSING_ROOT", None, None, "")
    if not os.path.exists(pqf):
        return (name, "MISSING_PARQUET", None, None, "")
    try:
        rows = pq.ParquetFile(pqf).metadata.num_rows
    except Exception as e:
        return (name, "PARQUET_UNREADABLE", None, None, str(e)[:120])
    try:
        f = uproot.open(rf)
        incl = f["inclusive"].num_entries
        parts = sum(f[k].num_entries for k in ("train", "validation", "test") if k in f)
        f.close()
    except Exception as e:
        return (name, "ROOT_UNREADABLE", rows, None, str(e)[:120])
    mt = os.path.getmtime(rf)
    if since and mt < since:
        return (name, "STALE_ROOT", rows, incl,
                "mtime %s" % time.strftime('%F %T', time.localtime(mt)))
    if rows != incl:
        return (name, "ROW_MISMATCH", rows, incl, "diff %+d" % (incl - rows))
    if parts != incl:
        return (name, "SPLIT_MISMATCH", rows, incl, "train+val+test=%d" % parts)
    return (name, "OK", rows, incl, "")


def main():
    since = 0.0
    if "--since" in sys.argv:
        since = float(sys.argv[sys.argv.index("--since") + 1])
    js = list(jobs())
    with ThreadPoolExecutor(max_workers=6) as ex:
        res = list(ex.map(lambda j: check(j, since), js))

    bad = [r for r in res if r[1] != "OK"]
    lines = ["GATE 3 -- merged parquet rows vs ROOT inclusive entries",
             "generated %s" % time.strftime('%F %T'),
             "expected conversions: %d   checked: %d   OK: %d   FAILED: %d"
             % (len(js), len(res), len(res) - len(bad), len(bad)), ""]
    if bad:
        lines.append("FAILURES (status, merged_rows, inclusive, note):")
        by = {}
        for r in bad:
            by.setdefault(r[1], []).append(r)
        for st in sorted(by):
            lines.append("  %s : %d" % (st, len(by[st])))
            for r in sorted(by[st])[:40]:
                lines.append("      %-40s %10s %10s  %s"
                             % (r[0], r[2], r[3], r[4]))
            if len(by[st]) > 40:
                lines.append("      ... and %d more" % (len(by[st]) - 40))
        lines.append("")
        lines.append("VERDICT: FAILED")
    else:
        tot = sum(r[3] for r in res)
        lines.append("total inclusive entries across all %d files: %d" % (len(res), tot))
        lines.append("VERDICT: PASSED")
    txt = "\n".join(lines)
    open(REPORT, "w").write(txt + "\n")
    print(txt)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
