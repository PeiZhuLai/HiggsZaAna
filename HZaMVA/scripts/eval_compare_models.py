#!/usr/bin/env python3
"""Compare two pickled main-BDT models on the SAME test set: AUC and R(mA).
  OLD = model_Za_BDT_run3.pkl.bak_preSignalReweight_20260705  (Jun-3 production, no signal reweight)
  NEW = model_Za_BDT_run3.pkl                                  (2026-07-05 retrain, signal reweight ON)
NOTE: this is a deployed-model vs deployed-model comparison, so the delta mixes (a) the signal
sideband reweight and (b) the fresh Optuna hyper-parameters. For the isolated reweight effect see
ab_signal_reweight.py (fixed HP).

Both models take the 16-feature matrix in order `variables + [param]` == hza_features FEATURES, so
hza_features._feat scores either pickle directly. Metrics mirror ab_signal_reweight.py:
  - R(mA): background-only sculpting ratio on inclusive bkg (top-0.2% score slice, peak 120-130
    fraction / inclusive peak fraction). R=1 => no sculpting.
  - AUC : pooled weighted test AUC over all masses, NOMINAL signal weights for both models.
CPU only (higgs-alp-ana).
"""
import pickle, json
import numpy as np
from hza_features import _load, _feat, FEATURES, MASSES, ROOT_DIR

OLD_PKL = "model_Za_BDT_run3.pkl.bak_preSignalReweight_20260705"
NEW_PKL = "model_Za_BDT_run3.pkl"
OUT_JSON = "eval_compare_models_results.json"

def weighted_auc(vals, y, w):
    o = np.argsort(vals, kind="mergesort"); y = y[o]; w = w[o]; is_b = y < 0.5
    below_b = np.cumsum(w * is_b) - w * is_b
    return (w[~is_b] * below_b[~is_b]).sum() / (w[~is_b].sum() * w[is_b].sum() + 1e-12)

def wq(v, wt, q):
    o = np.argsort(v); c = np.cumsum(wt[o]) / wt.sum(); return np.interp(q, c, v[o])

print("[load] models ...", flush=True)
old = pickle.load(open(OLD_PKL, "rb")); new = pickle.load(open(NEW_PKL, "rb"))

print("[load] signal test (all masses) + background test/inclusive ...", flush=True)
sig_te = {m: _load(f"{ROOT_DIR}/mA_M{m}/run3.root", "test", None) for m in MASSES}
bk_te  = _load(f"{ROOT_DIR}/All_Bkg/run3.root", "test", None)
bk_in  = _load(f"{ROOT_DIR}/All_Bkg/run3.root", "inclusive")

# precompute feature matrices (param set per mass); scores reused per model
Xs_te = {m: _feat(sig_te[m], m) for m in MASSES}
Xb_te = {m: _feat(bk_te,  m) for m in MASSES}
Xincl = {m: _feat(bk_in,  m) for m in MASSES}
ws_te = {m: np.clip(sig_te[m]["factor"], 0, None) for m in MASSES}   # NOMINAL signal weights
wb_te = np.clip(bk_te["factor"], 0, None)
Hmb = bk_in["H_m"]; wBe = np.clip(bk_in["factor"], 0, None)
peak = (Hmb > 120) & (Hmb < 130); frac = wBe[peak].sum() / wBe.sum()

def evaluate(clf, tag):
    print(f"[{tag}] scoring ...", flush=True)
    R = {}
    for m in MASSES:
        s = clf.predict_proba(Xincl[m])[:, 1]; thr = wq(s, wBe, 1 - 0.002); p = s > thr
        R[m] = float((wBe[p & peak].sum() / wBe[p].sum()) / frac) if wBe[p].sum() > 0 else float("nan")
    sv, yv, wv = [], [], []
    for m in MASSES:
        ss = clf.predict_proba(Xs_te[m])[:, 1]; sb = clf.predict_proba(Xb_te[m])[:, 1]
        wsm = ws_te[m]
        sv += [ss, sb]; yv += [np.ones(len(ss)), np.zeros(len(sb))]
        wv += [wsm * (wb_te.sum() / wsm.sum()), wb_te]
    auc = float(weighted_auc(np.concatenate(sv), np.concatenate(yv), np.concatenate(wv)))
    return {"auc": auc, "R": R,
            "maxR": float(np.nanmax(list(R.values()))),
            "mean_absR_minus1": float(np.nanmean([abs(R[m] - 1) for m in MASSES]))}

res_old = evaluate(old, "OLD")
res_new = evaluate(new, "NEW")

print("\n" + "=" * 74)
print("Deployed main HZa BDT: OLD (Jun-3, no sig reweight) vs NEW (2026-07-05, sig reweight)")
print("(delta mixes signal reweight + fresh Optuna HP; test set + metrics identical)")
print("=" * 74)
print(f"pooled test AUC :  OLD={res_old['auc']:.4f}   NEW={res_new['auc']:.4f}   dAUC={res_new['auc']-res_old['auc']:+.4f}")
print(f"max R(mA)       :  OLD={res_old['maxR']:.3f}    NEW={res_new['maxR']:.3f}    dmaxR={res_new['maxR']-res_old['maxR']:+.3f}")
print(f"mean|R-1|       :  OLD={res_old['mean_absR_minus1']:.3f}    NEW={res_new['mean_absR_minus1']:.3f}    d={res_new['mean_absR_minus1']-res_old['mean_absR_minus1']:+.3f}")
print("-" * 74)
print(f"{'mA':>4} {'R_OLD':>8} {'R_NEW':>8} {'dR':>8}")
for m in MASSES:
    print(f"{m:>4} {res_old['R'][m]:>8.3f} {res_new['R'][m]:>8.3f} {res_new['R'][m]-res_old['R'][m]:>+8.3f}")
print("=" * 74)

json.dump({"OLD_pkl": OLD_PKL, "NEW_pkl": NEW_PKL, "masses": MASSES,
           "OLD": res_old, "NEW": res_new}, open(OUT_JSON, "w"), indent=2)
print(f"[save] {OUT_JSON}")
