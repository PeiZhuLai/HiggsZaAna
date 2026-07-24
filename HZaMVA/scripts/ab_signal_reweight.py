#!/usr/bin/env python3
"""Controlled A/B for the MAIN HZa BDT: does applying the sideband reweight to SIGNAL
(before training) change AUC / R(mA)?  Everything is held fixed EXCEPT the signal training
weight:
  OLD arm : signal weight = nominal factor          (production before 2026-07-05)
  NEW arm : signal weight = factor x sideband reweight (true-mass param)
Background is sideband-reweighted in BOTH arms (unchanged). Same 16 features, same train/test
trees, same regularized hyper-parameters, same seed. Reuses hza_features so the feature matrix
and reweight are byte-identical to run3_Za_BDT.py.

Metrics (mirrors the production 'final' trainers):
  - R(mA): background-only sculpting ratio. On inclusive background, take the top-0.2% BDT score
    slice (per mass hypothesis), R = (peak 120-130 fraction inside the slice) / (inclusive peak
    fraction). R=1 => no sculpting. R depends on the trained model only, NOT on signal weights.
  - AUC : pooled weighted test AUC over all masses, evaluated with NOMINAL signal weights in BOTH
    arms so the number reflects the model's discrimination, not the weight convention.

CPU only (higgs-alp-ana). Prints an OLD-vs-NEW table and writes a JSON next to this script.
"""
import json
import numpy as np, xgboost as xgb
from hza_features import _load, _feat, FEATURES, MASSES, REWEIGHT_EXTRA, ROOT_DIR, \
                         _sideband_reweight, signal_weight

PARAM_IDX = FEATURES.index("param")
SEED = 42
N_SIG_PER = 40000      # per-mass signal cap for training (plenty; test/inclusive uncapped)
N_BKG_TR  = 400000     # background training cap (matches the aux 'final' trainers)
# Production regularized base config (run3_Za_BDT.py regularized_xgb_candidates, lines ~1081-1091)
XGB_PARAMS = dict(max_depth=3, learning_rate=0.06, n_estimators=800, min_child_weight=20,
                  gamma=1.0, reg_lambda=10.0, reg_alpha=0.1, subsample=0.8, colsample_bytree=0.8,
                  eval_metric="logloss", tree_method="hist", n_jobs=8, verbosity=0,
                  random_state=SEED)
OUT_JSON = "ab_signal_reweight_results.json"

def weighted_auc(vals, y, w):
    o = np.argsort(vals, kind="mergesort"); y = y[o]; w = w[o]; is_b = y < 0.5
    below_b = np.cumsum(w * is_b) - w * is_b
    return (w[~is_b] * below_b[~is_b]).sum() / (w[~is_b].sum() * w[is_b].sum() + 1e-12)

def wq(v, wt, q):
    o = np.argsort(v); c = np.cumsum(wt[o]) / wt.sum(); return np.interp(q, c, v[o])

# ---------------- load once (shared across both arms) ----------------
print("[load] signal (train+test) all masses + background (train+test+inclusive) ...", flush=True)
sig_tr = {m: _load(f"{ROOT_DIR}/mA_M{m}/run3.root", "train", N_SIG_PER, extra=REWEIGHT_EXTRA) for m in MASSES}
sig_te = {m: _load(f"{ROOT_DIR}/mA_M{m}/run3.root", "test",  None,      extra=REWEIGHT_EXTRA) for m in MASSES}
bk_tr  = _load(f"{ROOT_DIR}/All_Bkg/run3.root", "train", N_BKG_TR, extra=REWEIGHT_EXTRA)
bk_te  = _load(f"{ROOT_DIR}/All_Bkg/run3.root", "test",  None,     extra=REWEIGHT_EXTRA)
bk_in  = _load(f"{ROOT_DIR}/All_Bkg/run3.root", "inclusive")

# ---- shared feature matrices + shared background weight (identical in both arms) ----
Xs = np.vstack([_feat(sig_tr[m], m) for m in MASSES])                 # signal features (16), true param
rng = np.random.default_rng(SEED)
bhyp = rng.choice(MASSES, size=len(bk_tr["H_m"]))                     # per-event bkg mass hypothesis
Xb = _feat(bk_tr, 0); Xb[:, PARAM_IDX] = (bk_tr["ALP_m"] - bhyp) / bk_tr["H_m"]
X = np.vstack([Xs, Xb])
y = np.concatenate([np.ones(len(Xs)), np.zeros(len(Xb))])
wb = np.clip(bk_tr["factor"] * _sideband_reweight(bk_tr), 0, None)    # bkg reweighted (both arms)

# signal training weights: OLD = factor ; NEW = factor x sideband reweight
ws_old = np.concatenate([np.clip(sig_tr[m]["factor"], 0, None)       for m in MASSES])
ws_new = np.concatenate([np.clip(signal_weight(sig_tr[m], m), 0, None) for m in MASSES])

def balanced_w(ws):
    return np.concatenate([ws * (wb.sum() / ws.sum()), wb])          # sig scaled to bkg total yield

# ---- inclusive-background R(mA) setup (background-only sculpting metric) ----
Hmb = bk_in["H_m"]; wBe = np.clip(bk_in["factor"], 0, None)
peak = (Hmb > 120) & (Hmb < 130); frac = wBe[peak].sum() / wBe.sum()
Xincl = {m: _feat(bk_in, m) for m in MASSES}                         # inclusive bkg, param per mass

# ---- test set for pooled AUC (nominal signal weights in BOTH arms) ----
Xs_te = {m: _feat(sig_te[m], m) for m in MASSES}
Xb_te = {m: _feat(bk_te,  m) for m in MASSES}
ws_te = {m: np.clip(sig_te[m]["factor"], 0, None) for m in MASSES}   # NOMINAL -> isolates the model
wb_te = np.clip(bk_te["factor"], 0, None)

def train_and_eval(ws, tag):
    print(f"\n[{tag}] fit XGB {XGB_PARAMS['max_depth']}/{XGB_PARAMS['n_estimators']}/lr{XGB_PARAMS['learning_rate']} ...", flush=True)
    clf = xgb.XGBClassifier(**XGB_PARAMS)
    clf.fit(X, y, sample_weight=balanced_w(ws))
    # R(mA): top-0.2% score slice on inclusive background, peak fraction ratio
    R = {}
    for m in MASSES:
        s = clf.predict_proba(Xincl[m])[:, 1]; thr = wq(s, wBe, 1 - 0.002); p = s > thr
        R[m] = float((wBe[p & peak].sum() / wBe[p].sum()) / frac) if wBe[p].sum() > 0 else float("nan")
    # pooled test AUC over all masses (nominal signal weights)
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

res_old = train_and_eval(ws_old, "OLD (signal = factor)")
res_new = train_and_eval(ws_new, "NEW (signal = factor x reweight)")

# ---------------- report ----------------
print("\n" + "=" * 74)
print("A/B: main HZa BDT, signal sideband reweight OFF (OLD) vs ON (NEW)")
print("regularized HP depth3/n_est800/lr0.06 ; features=16 ; seed=%d" % SEED)
print("=" * 74)
print(f"pooled test AUC :  OLD={res_old['auc']:.4f}   NEW={res_new['auc']:.4f}   dAUC={res_new['auc']-res_old['auc']:+.4f}")
print(f"max R(mA)       :  OLD={res_old['maxR']:.3f}    NEW={res_new['maxR']:.3f}    dmaxR={res_new['maxR']-res_old['maxR']:+.3f}")
print(f"mean|R-1|       :  OLD={res_old['mean_absR_minus1']:.3f}    NEW={res_new['mean_absR_minus1']:.3f}    d={res_new['mean_absR_minus1']-res_old['mean_absR_minus1']:+.3f}")
print("-" * 74)
print(f"{'mA':>4} {'R_OLD':>8} {'R_NEW':>8} {'dR':>8}")
for m in MASSES:
    print(f"{m:>4} {res_old['R'][m]:>8.3f} {res_new['R'][m]:>8.3f} {res_new['R'][m]-res_old['R'][m]:>+8.3f}")
print("=" * 74)

json.dump({"config": {"xgb": {k: v for k, v in XGB_PARAMS.items() if k != "random_state"},
                      "seed": SEED, "n_sig_per": N_SIG_PER, "n_bkg_tr": N_BKG_TR,
                      "masses": MASSES, "features": FEATURES},
           "OLD": res_old, "NEW": res_new},
          open(OUT_JSON, "w"), indent=2)
print(f"[save] {OUT_JSON}")
