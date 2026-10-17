#!/usr/bin/env python3
"""Retrain the two BDTs selected by the input-variable search (bdt_varsearch, 2026-10-14) in float64.

The search jobs trained on float32 slim copies (bdt_varsearch/build_slim.py). Here the training sets are read
straight from run3_bdt_inputs_fsrfix 'train' trees in float64, with the event selection of hza_features._load
(95 < m_llgg < 180, finite factor; background 600k events, default_rng(0)) and the production weights
(hza_features.signal_weight / _sideband_reweight). Feature definitions are imported from search_job.column,
so the inputs are the ones the search evaluated. Recipe = search_job.train = train_{low,high}mass_final.py
(XGB depth 4, 300 trees, lr 0.1, subsample 0.8, min_child_weight 5); m_a = 1 signal restricted to 115-133 GeV
for the low-mass model (weight rescaled to the full m_a = 1 weight).

Check: the selected events must equal the slim training sets used by the search (count and m_llgg).
Env: higgs-alp-ana. Output: models/model_<low|high>.{pkl,meta.json}
Usage: train_selected.py low|high [...]
"""
import json
import os
import pickle
import sys
import time

import numpy as np
import uproot
import xgboost as xgb

HERE = os.path.dirname(os.path.abspath(__file__))
SEARCH = os.path.join(os.path.dirname(HERE), "bdt_varsearch")
sys.path.insert(0, SEARCH)
sys.path.insert(0, "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts")
import search_job as sj  # noqa: E402
import hza_features as hf  # noqa: E402

INP = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_inputs_fsrfix"
SLIM = "/eos/home-p/pelai/HZa/bdt_varsearch_20261013/samples"
SEL = json.load(open(os.path.join(HERE, "selected.json")))


def raw_columns(feats):
    cols = set(hf.READ) | set(hf.REWEIGHT_EXTRA) | {"H_m", "ALP_m", "factor", "pho1Pt", "pho2Pt"}
    for f in feats:
        if f in ("param_rel", "param", "pho_pt_asym", "pho_pt_asym_oHm"):
            continue
        if f in ("pho1_relEerr", "pho2_relEerr"):
            w = "lead" if f.startswith("pho1") else "sublead"
            cols |= {"ALP_%s_photon_energyErr" % w, "ALP_%s_photon_eta" % w}
            continue
        cols.add(f)
    return sorted(cols)


def load(path, tree, nmax, cols):
    """hza_features._load selection (mask, then default_rng(0) subsample) with any columns, float64."""
    a = uproot.open(path)[tree].arrays(cols, library="np")
    m = (a["H_m"] > 95) & (a["H_m"] < 180) & np.isfinite(a["factor"])
    a = {k: v[m].astype(np.float64) for k, v in a.items()}
    if nmax and len(a["H_m"]) > nmax:
        idx = np.random.default_rng(0).choice(len(a["H_m"]), nmax, replace=False)
        a = {k: v[idx] for k, v in a.items()}
    return a


def check_vs_slim(a, name):
    s = np.load(os.path.join(SLIM, name + ".npz"))
    same = len(s["H_m"]) == len(a["H_m"]) and np.array_equal(s["H_m"], a["H_m"].astype(np.float32))
    print("  [check] %-14s %d events, equal to the search slim set: %s" % (name, len(a["H_m"]), same), flush=True)
    if not same:
        sys.exit("selection differs from the slim set %s" % name)


def build(model):
    cfg = SEL[model]
    feats, masses = cfg["feats"], cfg["masses"]
    cols = raw_columns(feats)
    Xs, ws = [], []
    for m in masses:
        a = load("%s/mA_M%d/run3.root" % (INP, m), "train", None, cols)
        check_vs_slim(a, "train_sig_m%d" % m)
        w = np.clip(hf.signal_weight({k: a[k] for k in hf.READ + hf.REWEIGHT_EXTRA}, m), 0, None)
        x = sj.matrix(a, feats, float(m))
        if (cfg["win_mA1"] and m == 1) or m in cfg.get("win_masses", []):   # win_masses: 10-15, window training at m_a 2, 3 too
            keep = sj.inr(a["H_m"], *sj.WIN)
            print("  mA%d window" % m + " %s: keep %d/%d, %.3f of the weight" % (sj.WIN, keep.sum(), len(keep), w[keep].sum() / w.sum()))
            x, w = x[keep], w[keep] * (w.sum() / w[keep].sum())
        Xs.append(x); ws.append(w)
    Xs = np.vstack(Xs); ws = np.concatenate(ws)
    b = load("%s/All_Bkg/run3.root" % INP, "train", 600000, cols)
    check_vs_slim(b, "train_bkg")
    bhyp = np.random.default_rng(1).choice(masses, size=len(b["H_m"])).astype(np.float64)
    Xb = sj.matrix(b, feats, bhyp)
    wb = np.clip(b["factor"] * hf._sideband_reweight({k: b[k] for k in hf.READ + hf.REWEIGHT_EXTRA}), 0, None)
    X = np.vstack([Xs, Xb]); y = np.concatenate([np.ones(len(Xs)), np.zeros(len(Xb))])
    w = np.concatenate([ws * (wb.sum() / ws.sum()), wb])
    return X, y, w


def main():
    for model in sys.argv[1:] or ["low", "high"]:
        t0 = time.time()
        cfg = SEL[model]
        print("===== %s: %s" % (model, cfg["feats"]), flush=True)
        X, y, w = build(model)
        print("  X %s dtype %s, %d signal / %d background" % (X.shape, X.dtype, int(y.sum()), int((1 - y).sum())), flush=True)
        clf = xgb.XGBClassifier(max_depth=4, n_estimators=300, learning_rate=0.1, subsample=0.8, min_child_weight=5,
                                eval_metric="logloss", n_jobs=4, verbosity=0)
        clf.fit(X, y, sample_weight=w)
        out = os.path.join(HERE, "models", "model_%s" % model)
        pickle.dump(clf, open(out + ".pkl.partial", "wb")); os.replace(out + ".pkl.partial", out + ".pkl")
        json.dump({"model": model, "features": cfg["feats"], "masses": cfg["masses"], "mA1_window": sj.WIN if cfg["win_mA1"] else None,
                   "source": cfg["source"], "inputs": INP, "dtype": "float64", "xgboost": xgb.__version__,
                   "recipe": "XGB depth4/300/lr0.1/subsample0.8/mcw5 (search_job.train = train_*_final.py)"},
                  open(out + ".meta.json", "w"), indent=1)
        print("  [save] %s.pkl (%.0f s)" % (out, time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
