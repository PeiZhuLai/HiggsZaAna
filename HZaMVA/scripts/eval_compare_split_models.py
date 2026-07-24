#!/usr/bin/env python3
"""Compare OLD vs NEW *split* BDT models (the ones production scoring actually uses via
Parque2Root_BDT.add_mva_scores) on the SAME test set: AUC + R(mA).
  low  : mA 1-3, 5 features  (LOW_FEATURE_COLUMNS)
  high : mA 4-30, 16 features (HIGH_FEATURE_COLUMNS == hza_features FULL16)
OLD = *.bak_preSignalReweight_20260706 (Jun, no signal reweight)
NEW = model_Za_BDT_{low,high}mass_run3.pkl (2026-07-06 retrain, signal reweight ON)
Metrics mirror train_{low,high}mass_final.py: R(mA)=top-0.2% inclusive-bkg sculpting ratio
(peak 120-130 fraction / inclusive peak fraction); AUC=pooled weighted test AUC, nominal
signal weights for both models. CPU (higgs-alp-ana).
"""
import pickle, json
import numpy as np
from hza_features import _load, _feat, BASE_VARS, ROOT_DIR

FULL16 = BASE_VARS + ["pho_pt_asym", "param"]
LOW_MASSES  = [1, 2, 3]
HIGH_MASSES = [4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30]
LOW_FEATS = ["ALP_calculatedPhotonIso", "var_dR_g1g2", "pho1R9", "pho1Pt_oHm", "param"]
LOW_IDX   = [FULL16.index(n) for n in LOW_FEATS]
HIGH_IDX  = list(range(len(FULL16)))   # all 16

MODELS = {
    "low":  {"masses": LOW_MASSES,  "idx": LOW_IDX,
             "OLD": "model_Za_BDT_lowmass_run3.pkl.bak_preSignalReweight_20260706",
             "NEW": "model_Za_BDT_lowmass_run3.pkl"},
    "high": {"masses": HIGH_MASSES, "idx": HIGH_IDX,
             "OLD": "model_Za_BDT_highmass_run3.pkl.bak_preSignalReweight_20260706",
             "NEW": "model_Za_BDT_highmass_run3.pkl"},
}

def weighted_auc(vals, y, w):
    o = np.argsort(vals, kind="mergesort"); y = y[o]; w = w[o]; is_b = y < 0.5
    below_b = np.cumsum(w * is_b) - w * is_b
    return (w[~is_b] * below_b[~is_b]).sum() / (w[~is_b].sum() * w[is_b].sum() + 1e-12)

def wq(v, wt, q):
    o = np.argsort(v); c = np.cumsum(wt[o]) / wt.sum(); return np.interp(q, c, v[o])

print("[load] background test + inclusive ...", flush=True)
bk_te = _load(f"{ROOT_DIR}/All_Bkg/run3.root", "test", None)
bk_in = _load(f"{ROOT_DIR}/All_Bkg/run3.root", "inclusive")
Hmb = bk_in["H_m"]; wBe = np.clip(bk_in["factor"], 0, None)
peak = (Hmb > 120) & (Hmb < 130); frac = wBe[peak].sum() / wBe.sum()
wb_te = np.clip(bk_te["factor"], 0, None)

results = {}
for kind, cfg in MODELS.items():
    masses, idx = cfg["masses"], cfg["idx"]
    print(f"[load] {kind} signal test (mA {masses[0]}..{masses[-1]}) ...", flush=True)
    sig_te = {m: _load(f"{ROOT_DIR}/mA_M{m}/run3.root", "test", None) for m in masses}
    Xs_te  = {m: _feat(sig_te[m], m)[:, idx] for m in masses}
    ws_te  = {m: np.clip(sig_te[m]["factor"], 0, None) for m in masses}
    Xb_te  = {m: _feat(bk_te, m)[:, idx] for m in masses}
    Xincl  = {m: _feat(bk_in, m)[:, idx] for m in masses}

    def evaluate(pkl):
        clf = pickle.load(open(pkl, "rb"))
        R = {}
        for m in masses:
            s = clf.predict_proba(Xincl[m])[:, 1]; thr = wq(s, wBe, 1 - 0.002); p = s > thr
            R[m] = float((wBe[p & peak].sum() / wBe[p].sum()) / frac) if wBe[p].sum() > 0 else float("nan")
        sv, yv, wv = [], [], []
        for m in masses:
            ss = clf.predict_proba(Xs_te[m])[:, 1]; sb = clf.predict_proba(Xb_te[m])[:, 1]
            wsm = ws_te[m]
            sv += [ss, sb]; yv += [np.ones(len(ss)), np.zeros(len(sb))]
            wv += [wsm * (wb_te.sum() / wsm.sum()), wb_te]
        auc = float(weighted_auc(np.concatenate(sv), np.concatenate(yv), np.concatenate(wv)))
        return {"auc": auc, "R": R,
                "maxR": float(np.nanmax(list(R.values()))),
                "mean_absR_minus1": float(np.nanmean([abs(R[m] - 1) for m in masses]))}

    results[kind] = {"masses": masses, "OLD": evaluate(cfg["OLD"]), "NEW": evaluate(cfg["NEW"])}

print("\n" + "=" * 78)
print("DEPLOYED split BDT models: OLD (no sig reweight) vs NEW (sig reweight), same test set")
print("=" * 78)
for kind in ("low", "high"):
    r = results[kind]; o, n = r["OLD"], r["NEW"]
    print(f"\n### {kind}-mass (mA {r['masses'][0]}-{r['masses'][-1]}) ###")
    print(f"pooled test AUC :  OLD={o['auc']:.4f}   NEW={n['auc']:.4f}   dAUC={n['auc']-o['auc']:+.4f}")
    print(f"max R(mA)       :  OLD={o['maxR']:.3f}    NEW={n['maxR']:.3f}    dmaxR={n['maxR']-o['maxR']:+.3f}")
    print(f"mean|R-1|       :  OLD={o['mean_absR_minus1']:.3f}    NEW={n['mean_absR_minus1']:.3f}    d={n['mean_absR_minus1']-o['mean_absR_minus1']:+.3f}")
    print(f"{'mA':>4} {'R_OLD':>8} {'R_NEW':>8} {'dR':>8}")
    for m in r["masses"]:
        print(f"{m:>4} {o['R'][m]:>8.3f} {n['R'][m]:>8.3f} {n['R'][m]-o['R'][m]:>+8.3f}")
print("=" * 78)
json.dump(results, open("eval_compare_split_models_results.json", "w"), indent=2)
print("[save] eval_compare_split_models_results.json")
