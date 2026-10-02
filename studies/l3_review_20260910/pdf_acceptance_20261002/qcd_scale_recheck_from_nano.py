#!/usr/bin/env python
"""Recompute the QCD-scale acceptance of ../qcd_scale_acceptance.py for given samples with the
LHEScaleWeight read back from the input NanoAOD instead of the parquet.

Reason: in mA_M7_2024 the parquet carries LHEScaleWeight_* = 1 and a zero PDF band for all
11406 preselected events of one input file (4bcc6de5-...; HiggsDNA attach_lhe_weights fell back
to unit weights), i.e. 42% of the sample. That sample is the 2.92% "worst case" of the
QCD-scale study. Same method, same selection (preselection, weight_central), same denominator
(logs_fsrfix/gen_scale_sums.json).
Per-file cache in cache_scale/<sample>/; xrootd reads in worker processes with a hard timeout
(single reads hung > 20 min despite uproot's timeout), up to 6 rounds with rotating redirectors.
Usage: python qcd_scale_recheck_from_nano.py mA_M7_2024 [more samples]
"""
import sys, os, json, time, numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pdf_acceptance as pa
import uproot, awkward as ak, pyarrow.parquet as pq, fsspec_xrootd, XRootD.client  # noqa: F401

GEN = json.load(open("/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/logs_fsrfix/gen_scale_sums.json"))
NAMES = ["Zero", "One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight"]
SEVEN = [0, 1, 3, 4, 5, 7, 8]
REDIRS = ["root://xrootd-cms.infn.it/", "root://cms-xrd-global.cern.ch/", "root://cmsxrootd.fnal.gov/"]
CS = os.path.join(pa.WORK, "cache_scale")


def cpath(s, u):
    return pa.cache_path(s, u).replace(pa.CACHE, CS)


def one(args):
    s, u, rnd = args
    lfn = "/store" + u.split("/store", 1)[1]
    with uproot.open(REDIRS[rnd % len(REDIRS)] + lfn, timeout=300) as f:
        parts = list(f["Events"].iterate(["run", "luminosityBlock", "event", "LHEScaleWeight"], step_size=50000))
    a = ak.concatenate(parts)
    out = cpath(s, u)
    np.savez_compressed(out + ".part", key=pa.keyof(a["run"], a["luminosityBlock"], a["event"]),
                        scale=ak.to_numpy(a["LHEScaleWeight"]).astype(np.float32))
    os.replace(out + ".part.npz", out)
    return u


def fetch(s):
    os.makedirs(os.path.join(CS, s), exist_ok=True)
    for rnd in range(6):
        todo = [u for u in pa.input_files(s) if not os.path.exists(cpath(s, u))]
        if not todo:
            return
        print("[round %d] %s: %d files to read via %s" % (rnd, s, len(todo), REDIRS[rnd % len(REDIRS)]), flush=True)
        pool = Pool(min(6, len(todo)))
        jobs = [(u, pool.apply_async(one, ((s, u, rnd),))) for u in todo]
        deadline = time.time() + 900
        for u, j in jobs:
            try:
                j.get(timeout=max(1.0, deadline - time.time())); print("  ok", os.path.basename(u), flush=True)
            except Exception as e:  # noqa: BLE001
                print("  FAIL", os.path.basename(u), type(e).__name__, str(e)[:100], flush=True)
        pool.terminate(); pool.join()
    sys.exit("could not read all files of " + s)


out = {}
for s in sys.argv[1:]:
    fetch(s)
    t = pq.read_table(os.path.join(pa.SIG, s, "merged_nominal.parquet"),
                      columns=["run", "luminosityBlock", "event", "weight_central"] + ["LHEScaleWeight_" + n for n in NAMES])
    k = pa.keyof(np.asarray(t["run"]), np.asarray(t["luminosityBlock"]), np.asarray(t["event"]))
    w = np.asarray(t["weight_central"], dtype=float)
    pq_sc = np.stack([np.asarray(t["LHEScaleWeight_" + n]) for n in NAMES], axis=1).astype(float)
    sc = np.full((len(k), 9), np.nan)
    for u in pa.input_files(s):
        z = np.load(cpath(s, u)); kk = z["key"]
        m = np.isin(k, kk)
        order = np.argsort(kk); idx = order[np.searchsorted(kk, k[m], sorter=order)]
        sc[m] = z["scale"][idx]
    assert not np.isnan(sc).any(), "unmatched events"
    same = np.all(np.abs(sc - pq_sc) < 1e-5, axis=1)
    unit = np.all(pq_sc == 1.0, axis=1)
    gen = np.asarray(GEN[s], dtype=float); gen /= gen[4]
    sel = (w[:, None] * sc).sum(axis=0); sel /= sel[4]
    sel_pq = (w[:, None] * pq_sc).sum(axis=0); sel_pq /= sel_pq[4]
    acc = sel / gen; acc_pq = sel_pq / gen
    r = {"n": int(len(k)), "n_parquet_equal_nano": int(same.sum()), "n_parquet_unit": int(unit.sum()),
         "n_mismatch_not_unit": int((~same & ~unit).sum()),
         "acc_up": float(acc[SEVEN].max() - 1), "acc_dn": float(acc[SEVEN].min() - 1),
         "yield_up": float(sel[SEVEN].max() - 1), "yield_dn": float(sel[SEVEN].min() - 1),
         "parquet_acc_up": float(acc_pq[SEVEN].max() - 1), "parquet_acc_dn": float(acc_pq[SEVEN].min() - 1),
         "parquet_yield_up": float(sel_pq[SEVEN].max() - 1), "parquet_yield_dn": float(sel_pq[SEVEN].min() - 1)}
    out[s] = r
    print(json.dumps({s: r}, indent=1), flush=True)
os.makedirs(pa.RES, exist_ok=True)
json.dump(out, open(os.path.join(pa.RES, "qcd_scale_recheck_from_nano.json"), "w"), indent=1)
