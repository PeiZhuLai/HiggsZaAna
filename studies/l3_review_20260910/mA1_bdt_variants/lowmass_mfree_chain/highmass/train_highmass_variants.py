#!/usr/bin/env python3
"""Sandbox retraining of the high-mass (m_a = 4-30) BDT with the input definitions of the new low-mass
model, for the inputs the two models share (user rule 2026-10-13: a variable used by both models must
have the same definition in both):
  pho1Pt_oHm -> pho1Pt = pho1Pt_oHm * m_llgg (GeV)
  param      -> param_rel = (m_gg - m_a,hyp) / m_a,hyp
The other 14 inputs and the recipe are copied from HZaMVA/scripts/train_highmass_final.py
(same inputs, sampling, weights, hyper-parameters). No signal-window restriction (that is a choice of
training sample for m_a = 1, not a variable definition).

  base    production definitions: must reproduce the deployed model_Za_BDT_highmass_run3.pkl
  shared  the two shared inputs redefined as above
Env: higgs-alp-ana. Output: model_hm_<variant>.{pkl,meta.json} next to this script.
"""
import json
import os
import pickle
import sys

import numpy as np
import xgboost as xgb

sys.path.insert(0, "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts")
from hza_features import _load, _feat, BASE_VARS, _sideband_reweight, signal_weight, ROOT_DIR  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
MASSES = [4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30]
FULL16 = BASE_VARS + ["pho_pt_asym", "param"]
FEATS = ["pho1Pt_oHm", "pho1R9", "pho1IetaIeta55", "pho1PIso_noCorr", "pho2Pt_oHm", "pho2R9", "pho2IetaIeta55",
         "pho2PIso_noCorr", "ALP_calculatedPhotonIso", "var_dR_Za", "var_dR_g1g2", "var_dR_g1Z", "var_PtaOverMh",
         "H_pt_oHm", "pho_pt_asym", "param"]
IDX = [FULL16.index(n) for n in FEATS]
NAMES = {"base": FEATS, "shared": [{"pho1Pt_oHm": "pho1Pt", "param": "param_rel"}.get(f, f) for f in FEATS]}


def matrix(a, mhyp, variant):
    """16-column matrix in FEATS order; mhyp scalar or per-event array."""
    X = _feat(a, 0)
    X[:, FULL16.index("param")] = (a["ALP_m"] - mhyp) / a["H_m"]
    X = X[:, IDX]
    if variant == "shared":
        X[:, FEATS.index("pho1Pt_oHm")] = a["pho1Pt_oHm"] * a["H_m"]
        X[:, FEATS.index("param")] = (a["ALP_m"] - mhyp) / mhyp
    return X


def build(tree, variant, n_sig=None, n_bkg=600000):
    sig = {m: _load(f"{ROOT_DIR}/mA_M{m}/run3.root", tree, n_sig, extra=["Z_m", "event"]) for m in MASSES}
    bk = _load(f"{ROOT_DIR}/All_Bkg/run3.root", tree, n_bkg, extra=["Z_m", "event"])
    Xs = np.vstack([matrix(sig[m], float(m), variant) for m in MASSES])
    ws = np.concatenate([np.clip(signal_weight(sig[m], m), 0, None) for m in MASSES])
    bhyp = np.random.default_rng(1).choice(MASSES, size=len(bk["H_m"])).astype(np.float64)
    Xb = matrix(bk, bhyp, variant)
    wb = np.clip(bk["factor"] * _sideband_reweight(bk), 0, None)
    X = np.vstack([Xs, Xb]); y = np.concatenate([np.ones(len(Xs)), np.zeros(len(Xb))])
    w = np.concatenate([ws * (wb.sum() / ws.sum()), wb])
    return X, y, w


def main():
    for v in sys.argv[1:] or ["base", "shared"]:
        print(f"===== high-mass variant {v}: {NAMES[v]}", flush=True)
        X, y, w = build("train", v)
        clf = xgb.XGBClassifier(max_depth=4, n_estimators=300, learning_rate=0.1, subsample=0.8,
                                min_child_weight=5, eval_metric="logloss", n_jobs=4, verbosity=0)
        clf.fit(X, y, sample_weight=w)
        out = os.path.join(HERE, f"model_hm_{v}")
        pickle.dump(clf, open(out + ".pkl", "wb"))
        json.dump({"variant": v, "features": NAMES[v], "masses": MASSES, "root_dir": ROOT_DIR,
                   "model": "XGB depth4/300/lr0.1 (train_highmass_final.py recipe)"}, open(out + ".meta.json", "w"), indent=2)
        print(f"[save] {out}.pkl", flush=True)


if __name__ == "__main__":
    main()
