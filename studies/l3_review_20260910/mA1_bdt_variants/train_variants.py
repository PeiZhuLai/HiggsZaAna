#!/usr/bin/env python3
"""Sandbox retraining of the low-mass (m_a = 1-3 GeV) BDT to test whether a change of the
training removes the m_a = 1 GeV second m_llgg population (overlapping showers whose photon
energy is overestimated, see doc/HZa/hza_mA1_overlapping_showers_run2_vs_run3_2026-10-13.md).

The recipe is a copy of HZaMVA/scripts/train_lowmass_final.py (same inputs, sampling, weights,
hyper-parameters); only the variant-specific change differs:

  base      production recipe, used to check that the retraining reproduces the deployed model
  win       (1) m_a = 1 signal training events restricted to the signal-model window
            115 < m_llgg < 133 GeV (Signal/scripts/signalFit.py), rescaled to the same total
            weight so the m_a = 1 share in the training is unchanged; m_a = 2, 3 unchanged
            (their signal models use the full range)
  nopt      (2) pho1Pt_oHm removed from the inputs; the boost information is left to var_dR_g1g2
  win_nopt  (1) + (2)
  win_mfree (1) with every input free of m_llgg: pho1Pt_oHm -> pho1Pt (GeV) and
            param -> (m_gg - m_a,hyp) / m_a,hyp. Diagnostic for the sculpting seen with (1):
            pho1Pt_oHm and param both divide by m_llgg, so a window-restricted signal can teach
            the BDT m_llgg directly.

Nothing in HZaMVA/using or the production chain is touched; models are written next to this
script as model_<variant>.{pkl,meta.json}.
Env: higgs-alp-ana. Usage: train_variants.py <variant> [<variant> ...]
"""
import json
import os
import pickle
import sys

import numpy as np
import xgboost as xgb

sys.path.insert(0, "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts")
from hza_features import _load, BASE_VARS, _sideband_reweight, signal_weight, ROOT_DIR  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
LOWM = [1, 2, 3]
FULL16 = BASE_VARS + ["pho_pt_asym", "param"]
WINDOW_MA1 = (115.0, 133.0)
VARIANTS = {
    "base":     {"feats": ["ALP_calculatedPhotonIso", "var_dR_g1g2", "pho1R9", "pho1Pt_oHm", "param"], "win": False},
    "win":      {"feats": ["ALP_calculatedPhotonIso", "var_dR_g1g2", "pho1R9", "pho1Pt_oHm", "param"], "win": True},
    "nopt":     {"feats": ["ALP_calculatedPhotonIso", "var_dR_g1g2", "pho1R9", "param"], "win": False},
    "win_nopt": {"feats": ["ALP_calculatedPhotonIso", "var_dR_g1g2", "pho1R9", "param"], "win": True},
    "win_mfree": {"feats": ["ALP_calculatedPhotonIso", "var_dR_g1g2", "pho1R9", "pho1Pt", "param_rel"], "win": True},
}
DERIVED = ("pho1Pt", "param_rel")


def matrix(a, mhyp, feats):
    """Feature matrix; mhyp may be a scalar or a per-event array (background hypotheses)."""
    cols = []
    for f in feats:
        if f == "param":
            cols.append((a["ALP_m"] - mhyp) / a["H_m"])
        elif f == "param_rel":
            cols.append((a["ALP_m"] - mhyp) / mhyp)
        elif f == "pho1Pt":
            cols.append(a["pho1Pt_oHm"] * a["H_m"])
        else:
            cols.append(a[f])
    return np.column_stack(cols)


def build(tree, feats, win, n_sig=None, n_bkg=600000):
    sig = {m: _load(f"{ROOT_DIR}/mA_M{m}/run3.root", tree, n_sig, extra=["Z_m", "event"]) for m in LOWM}
    Xs, ws = [], []
    for m in LOWM:
        x = matrix(sig[m], m, feats)
        w = np.clip(signal_weight(sig[m], m), 0, None)
        if win and m == 1:
            keep = (sig[m]["H_m"] > WINDOW_MA1[0]) & (sig[m]["H_m"] < WINDOW_MA1[1])
            print(f"[build:{tree}] mA1 window {WINDOW_MA1}: keep {keep.sum()}/{len(keep)} events, "
                  f"{w[keep].sum() / w.sum():.3f} of the weight (rescaled back to the full weight)")
            x, w = x[keep], w[keep] * (w.sum() / w[keep].sum())
        Xs.append(x); ws.append(w)
    Xs = np.vstack(Xs); ws = np.concatenate(ws)
    bk = _load(f"{ROOT_DIR}/All_Bkg/run3.root", tree, n_bkg, extra=["Z_m", "event"])
    bhyp = np.random.default_rng(1).choice(LOWM, size=len(bk["H_m"]))
    Xb = matrix(bk, bhyp.astype(np.float64), feats)
    wb = np.clip(bk["factor"] * _sideband_reweight(bk), 0, None)
    X = np.vstack([Xs, Xb]); y = np.concatenate([np.ones(len(Xs)), np.zeros(len(Xb))])
    w = np.concatenate([ws * (wb.sum() / ws.sum()), wb])
    return X, y, w


def main():
    names = sys.argv[1:] or list(VARIANTS)
    for name in names:
        cfg = VARIANTS[name]
        print(f"\n===== variant {name}: features {cfg['feats']}, mA1 window {cfg['win']}", flush=True)
        X, y, w = build("train", cfg["feats"], cfg["win"])
        clf = xgb.XGBClassifier(max_depth=4, n_estimators=300, learning_rate=0.1, subsample=0.8,
                                min_child_weight=5, eval_metric="logloss", n_jobs=4, verbosity=0)
        clf.fit(X, y, sample_weight=w)
        out = os.path.join(HERE, f"model_{name}")
        pickle.dump(clf, open(out + ".pkl", "wb"))
        json.dump({"variant": name, "features": cfg["feats"], "mA1_window": WINDOW_MA1 if cfg["win"] else None,
                   "masses": LOWM, "root_dir": ROOT_DIR, "model": "XGB depth4/300/lr0.1 (train_lowmass_final.py recipe)"},
                  open(out + ".meta.json", "w"), indent=2)
        print(f"[save] {out}.pkl", flush=True)


if __name__ == "__main__":
    main()
