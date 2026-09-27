#!/usr/bin/env python
"""Sum genEventSumw and LHEScaleSumw over every signal input file, per sample.

Needed as the DENOMINATOR of the acceptance ratio: without it the scale
variation of the cross section is folded into what looks like an acceptance
uncertainty, and would be counted twice against the separate xs uncertainty.

Reads only the Runs tree (a handful of values per file), so it is cheap even
though it touches every input file.
"""
import json, glob, os, sys
from concurrent.futures import ThreadPoolExecutor
import numpy as np, uproot

SIG = "/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix/Sig_MC"
OUT = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/logs_fsrfix/gen_scale_sums.json"

def one(u):
    try:
        r = uproot.open(u, timeout=120)["Runs"]
        g = float(np.asarray(r["genEventSumw"].array(library="np")).sum())
        s = np.asarray(r["LHEScaleSumw"].array(library="np"))
        s = np.vstack([np.asarray(x, dtype=float) for x in s]) if s.dtype == object else np.atleast_2d(s)
        gs = np.asarray(r["genEventSumw"].array(library="np"), dtype=float)
        return g, (s * gs[:, None]).sum(axis=0)
    except Exception as e:
        print("FAIL", os.path.basename(u), type(e).__name__, str(e)[:70], flush=True)
        return None

res = {}
for d in sorted(glob.glob(os.path.join(SIG, "*"))):
    smp = os.path.basename(d)
    files = []
    for cfg in glob.glob(os.path.join(d, "job_*", "*_config_job*.json")):
        files += json.load(open(cfg)).get("files") or []
    files = sorted(set(files))
    if not files:
        continue
    with ThreadPoolExecutor(max_workers=6) as ex:
        out = [x for x in ex.map(one, files) if x]
    if not out:
        continue
    gtot = sum(o[0] for o in out)
    stot = np.sum([o[1] for o in out], axis=0)
    res[smp] = (stot / gtot).tolist()
    print("%-28s files=%4d/%4d  gen=%.4g  ratios=%s" %
          (smp, len(out), len(files), gtot,
           " ".join("%.4f" % v for v in (stot / gtot))), flush=True)
json.dump(res, open(OUT, "w"), indent=1)
print("wrote", OUT)
