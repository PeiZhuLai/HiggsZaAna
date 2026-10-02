#!/usr/bin/env python3
"""Category 3 of cleanup_old_20260929.sh: delete chunk parquets of the current production.

A sample dir is eligible only if, for EVERY merged_<syst>.parquet at its top level,
  sum(num_rows of job_*/output_job_*_<syst>.parquet) == num_rows(merged_<syst>.parquet)
  and mtime(merged) >= max(mtime of those chunks)   (merge-before-rechunk race)
and every chunk syst has a merged file. merged_*.parquet, summaries and configs are kept.
(merged == ROOT was checked by GATE 5 / stage 3 of the scoring and p2root chains.)

Usage: cleanup_chunks.py <APPLY 0|1>
"""
import os, re, subprocess, sys, time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
import pyarrow.parquet as pq

APPLY = len(sys.argv) > 1 and sys.argv[1] == "1"
MGM = "root://eosproject-h.cern.ch"
P = "/eos/project/h/htozg-dy-privatemc/pelai/HZa"
PARENTS = [f"{P}/parquet_DNA_tmp_fsrfix/Sig_MC", f"{P}/parquet_DNA_tmp_fsrfix_fpo1/Bkg_MC",
           f"{P}/parquet_DNA_tmp_fsrfix_fpo1/Data", f"{P}/parquet_DNA_tmp_fsrfix_extraBkg/Bkg_MC_extraBkg2024"]
REPORT = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/logs_fsrfix/cleanup_chunks_report.txt"
WORK = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/logs_fsrfix/cleanup_chunks_work"
CH = re.compile(r"/job_\d+/output_job_\d+_(.+)\.parquet$")
env = dict(os.environ, EOS_MGM_URL=MGM)


def eos_find(d):
    r = subprocess.run(["eos", "find", "-f", d], capture_output=True, text=True, env=env, timeout=600)
    if r.returncode != 0:
        raise RuntimeError(f"eos find {d}: {r.stderr[:200]}")
    return [l.strip() for l in r.stdout.splitlines() if l.strip()]


def meta(f):
    return pq.ParquetFile(f).metadata.num_rows, os.path.getmtime(f)


def check(sd):
    files = eos_find(sd)
    merged = {os.path.basename(f)[7:-8]: f for f in files
              if os.path.dirname(f) == sd and os.path.basename(f).startswith("merged_") and f.endswith(".parquet")}
    chunks = defaultdict(list)
    for f in files:
        m = CH.search(f)
        if m:
            chunks[m.group(1)].append(f)
    if not chunks:
        return sd, "NO_CHUNKS", [], ""
    why = []
    if set(chunks) - set(merged):
        why.append("chunk systs without merged: %s" % sorted(set(chunks) - set(merged)))
    for s, mf in sorted(merged.items()):
        mr, mt = meta(mf)
        with ThreadPoolExecutor(12) as ex:
            res = list(ex.map(meta, chunks.get(s, [])))
        cr = sum(r for r, _ in res); cmax = max((t for _, t in res), default=0)
        if cr != mr:
            why.append(f"{s}: chunk_rows {cr} != merged_rows {mr}")
        if cmax > mt:
            why.append(f"{s}: late chunk (merge race)")
    allc = [f for v in chunks.values() for f in v]
    return sd, ("OK" if not why else "FAIL"), allc, "; ".join(why)


def main():
    sds = []
    for par in PARENTS:
        sds += [f"{par}/{d}" for d in sorted(os.listdir(par)) if os.path.isdir(f"{par}/{d}")]
    out, todel = [], []
    with ThreadPoolExecutor(4) as ex:
        for sd, st, allc, why in ex.map(check, sds):
            out.append(f"{st:9s} chunks={len(allc):6d} {sd} {why}")
            print(out[-1], flush=True)
            if st == "OK":
                todel += allc
    nfail = sum(1 for l in out if l.startswith("FAIL"))
    with open(REPORT, "w") as f:
        f.write("generated %s APPLY=%s\n" % (time.strftime("%F %T"), APPLY))
        f.write("\n".join(out) + "\n")
        f.write(f"sample dirs {len(sds)}  FAIL {nfail}  chunk files eligible {len(todel)}\n")
    print(f"sample dirs {len(sds)} FAIL {nfail} eligible chunk files {len(todel)}", flush=True)
    if not APPLY or not todel:
        return 0
    os.makedirs(WORK, exist_ok=True)
    parts = [todel[i::6] for i in range(6)]
    procs = []
    for i, part in enumerate(parts):
        fn = f"{WORK}/rm_{i}.txt"
        with open(fn, "w") as f:
            f.write("".join(f'rm "{p}"\n' for p in part))
        procs.append(subprocess.Popen(f"eos < {fn} > {fn}.log 2>&1", shell=True, env=env))
    for p in procs:
        p.wait()
    left = [p for p in todel if os.path.exists(p)]
    with open(REPORT, "a") as f:
        f.write(f"deleted {len(todel) - len(left)} / {len(todel)}; left {len(left)}\n")
    print(f"deleted {len(todel) - len(left)} / {len(todel)}; left {len(left)}", flush=True)
    return 1 if left or nfail else 0


if __name__ == "__main__":
    sys.exit(main())
