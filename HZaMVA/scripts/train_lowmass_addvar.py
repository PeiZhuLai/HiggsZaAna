#!/usr/bin/env python3
"""Extend the low-mass BDT to mA=1-3 (current scope) and find the single most worthwhile variable to ADD:
maximize AUC gain while keeping sculpting (R) lowest. CPU (higgs-alp-ana); data loaded once.
Baseline = the current 5-feature low-mass set; candidates = every other input variable.
"""
import json, pickle, itertools
import numpy as np, xgboost as xgb
from hza_features import _load, _feat, BASE_VARS, _sideband_reweight, signal_weight, ROOT_DIR

LOWM = [1, 2, 3]
FULL16 = BASE_VARS + ["pho_pt_asym", "param"]
BASE = ["ALP_calculatedPhotonIso", "var_dR_g1g2", "pho1R9", "pho1Pt_oHm"]   # current low-mass set (+param)
POOL = [v for v in BASE_VARS + ["pho_pt_asym"] if v not in BASE]            # candidates to add

def weighted_auc(vals, y, w):
    o = np.argsort(vals, kind="mergesort"); y = y[o]; w = w[o]; is_b = y < 0.5
    below_b = np.cumsum(w * is_b) - w * is_b
    return (w[~is_b] * below_b[~is_b]).sum() / (w[~is_b].sum() * w[is_b].sum() + 1e-12)

def wq(v, wt, q): o = np.argsort(v); c = np.cumsum(wt[o]) / wt.sum(); return np.interp(q, c, v[o])

print("[load] signal mA=1,2,3 + All_Bkg ...")
sig_tr = {m: _load(f"{ROOT_DIR}/mA_M{m}/run3.root", "train", 40000, extra=["Z_m", "event"]) for m in LOWM}
sig_in = {m: _load(f"{ROOT_DIR}/mA_M{m}/run3.root", "inclusive", extra=["Z_m", "event"]) for m in LOWM}
bk_tr = _load(f"{ROOT_DIR}/All_Bkg/run3.root", "train", 400000, extra=["Z_m", "event"])
bk_in = _load(f"{ROOT_DIR}/All_Bkg/run3.root", "inclusive")

Xs16 = np.vstack([_feat(sig_tr[m], m) for m in LOWM])
ws = np.concatenate([np.clip(signal_weight(sig_tr[m], m), 0, None) for m in LOWM])  # signal: factor x sideband reweight
bhyp = np.random.default_rng(1).choice(LOWM, size=len(bk_tr["H_m"]))
Xb16 = _feat(bk_tr, 0); Xb16[:, FULL16.index("param")] = (bk_tr["ALP_m"] - bhyp) / bk_tr["H_m"]
wb = np.clip(bk_tr["factor"] * _sideband_reweight(bk_tr), 0, None)
y = np.concatenate([np.ones(len(Xs16)), np.zeros(len(Xb16))]); w = np.concatenate([ws * (wb.sum()/ws.sum()), wb])
X16 = np.vstack([Xs16, Xb16])
Hmb = bk_in["H_m"]; wBe = np.clip(bk_in["factor"], 0, None); peak = (Hmb > 120) & (Hmb < 130); frac = wBe[peak].sum()/wBe.sum()
Xb_in16 = {m: _feat(bk_in, m) for m in LOWM}; Xs_in16 = {m: _feat(sig_in[m], m) for m in LOWM}

def evaluate(keep):
    idx = [FULL16.index(n) for n in keep + ["param"]]
    clf = xgb.XGBClassifier(max_depth=4, n_estimators=300, learning_rate=0.1, subsample=0.8,
                            min_child_weight=5, eval_metric="logloss", n_jobs=4, verbosity=0)
    clf.fit(X16[:, idx], y, sample_weight=w)
    R = {}
    for m in LOWM:
        s = clf.predict_proba(Xb_in16[m][:, idx])[:, 1]; thr = wq(s, wBe, 1 - 0.002); p = s > thr
        R[m] = (wBe[p & peak].sum() / wBe[p].sum()) / frac if wBe[p].sum() > 0 else np.nan
    sv, yv, wv = [], [], []
    for m in LOWM:
        ss = clf.predict_proba(Xs_in16[m][:, idx])[:, 1]; sb = clf.predict_proba(Xb_in16[m][:, idx])[:, 1]
        wss = np.clip(signal_weight(sig_in[m], m), 0, None)
        sv += [ss, sb]; yv += [np.ones(len(ss)), np.zeros(len(sb))]; wv += [wss * (wBe.sum()/wss.sum()), wBe]
    auc = weighted_auc(np.concatenate(sv), np.concatenate(yv), np.concatenate(wv))
    return auc, R

aucB, RB = evaluate(BASE)
print(f"\n[baseline mA=1-3] feats={BASE}+param  AUC={aucB:.4f}  "
      f"R=" + "/".join(f"{RB[m]:.2f}" for m in LOWM) + f"  maxR={max(RB.values()):.2f}")
rh = "".join(f"{'R'+str(m):>6}" for m in LOWM)
print(f"\n{'+variable':<26}{'AUC':>8}{'dAUC':>8}{rh}{'maxR':>7}{'dmaxR':>7}")
rows = []
for cand in POOL:
    auc, R = evaluate(BASE + [cand])
    mR = max(R.values()); rows.append((cand, auc, auc - aucB, R, mR, mR - max(RB.values())))
for cand, auc, dauc, R, mR, dmR in sorted(rows, key=lambda r: -r[2]):
    rr = "".join(f"{R[m]:>6.2f}" for m in LOWM)
    print(f"{cand:<26}{auc:>8.4f}{dauc:>+8.4f}{rr}{mR:>7.2f}{dmR:>+7.2f}")

# recommend: best AUC gain among those that keep maxR low (<= baseline maxR + 0.15)
tol = max(RB.values()) + 0.15
cand_ok = [r for r in rows if r[4] <= tol]
best = max(cand_ok or rows, key=lambda r: r[2])
print(f"\n>>> baseline maxR={max(RB.values()):.2f}; keep-low-sculpting tol(maxR<={tol:.2f})")
print(f">>> BEST to add: +{best[0]}  AUC {aucB:.4f}->{best[1]:.4f} (dAUC {best[2]:+.4f})  "
      f"maxR {max(RB.values()):.2f}->{best[4]:.2f}")
