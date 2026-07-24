#!/usr/bin/env python3
"""Push the dedicated low-mA (mA=1,2,3) XGBoost to MINIMAL variable sets to drive R(mA) down.
CPU only (higgs-alp-ana); data loaded ONCE, reused. KEEP-based variants (param always added).

Diagnostic ranking (sep AUC@mA1 / |corr(m_llgg)|):
  ALP_calculatedPhotonIso 0.996/0.010  var_dR_g1g2 0.773/0.068  var_PtaOverMh 0.746/0.025
  H_pt_oHm 0.650/0.093  pho2R9 0.639/0.072  pho1R9 0.634/0.065  pho1Pt_oHm 0.625/0.042
Key idea: ID variables (iso, R9) don't encode kinematics -> can't rebuild m_llgg -> least sculpt;
kinematic vars (dR, pT, PtaOverMh) discriminate but can combine into m_llgg -> more sculpt.
"""
import numpy as np, xgboost as xgb
from hza_features import _load, _feat, BASE_VARS, _sideband_reweight, signal_weight, ROOT_DIR

LOWM = [1, 2, 3]
FULL16 = BASE_VARS + ["pho_pt_asym", "param"]
VARIANTS = {  # base-var names to KEEP (param appended automatically)
    "clean8":  ["ALP_calculatedPhotonIso", "var_dR_g1g2", "var_PtaOverMh", "H_pt_oHm", "pho2R9", "pho1R9", "pho1Pt_oHm"],
    "min6":    ["ALP_calculatedPhotonIso", "var_dR_g1g2", "var_PtaOverMh", "pho1Pt_oHm", "pho1R9"],
    "min4":    ["ALP_calculatedPhotonIso", "var_dR_g1g2", "var_PtaOverMh", "pho1R9"],
    "kin3":    ["ALP_calculatedPhotonIso", "var_dR_g1g2", "var_PtaOverMh"],
    "id3":     ["ALP_calculatedPhotonIso", "pho1R9", "pho2R9"],          # ID-only, no kinematics
    "iso+pta": ["ALP_calculatedPhotonIso", "var_PtaOverMh"],
    "iso1":    ["ALP_calculatedPhotonIso"],
}

def cols_for(keep):
    names = keep + ["param"]
    return names, [FULL16.index(n) for n in names]

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
Hmb = bk_in["H_m"]; wBe = np.clip(bk_in["factor"], 0, None)
peak = (Hmb > 120) & (Hmb < 130); frac = wBe[peak].sum() / wBe.sum()
Xb_in16 = {m: _feat(bk_in, m) for m in LOWM}; Xs_in16 = {m: _feat(sig_in[m], m) for m in LOWM}

print(f"\n{'variant':<10}{'nvar':>5}{'AUC(1-3)':>10}{'R(mA1)':>8}{'R(mA2)':>8}{'R(mA3)':>8}  keep")
for name, keep in VARIANTS.items():
    names, idx = cols_for(keep)
    clf = xgb.XGBClassifier(max_depth=4, n_estimators=300, learning_rate=0.1, subsample=0.8,
                            min_child_weight=5, eval_metric="logloss", n_jobs=4, verbosity=0)
    clf.fit(X16[:, idx], y, sample_weight=w)
    R = {}
    for m in LOWM:
        s = clf.predict_proba(Xb_in16[m][:, idx])[:, 1]
        thr = wq(s, wBe, 1 - 0.002); p = s > thr
        R[m] = (wBe[p & peak].sum() / wBe[p].sum()) / frac if wBe[p].sum() > 0 else float("nan")
    sv, yv, wv = [], [], []
    for m in LOWM:
        ss = clf.predict_proba(Xs_in16[m][:, idx])[:, 1]; sb = clf.predict_proba(Xb_in16[m][:, idx])[:, 1]
        wss = np.clip(signal_weight(sig_in[m], m), 0, None)
        sv += [ss, sb]; yv += [np.ones(len(ss)), np.zeros(len(sb))]; wv += [wss * (wBe.sum()/wss.sum()), wBe]
    auc = weighted_auc(np.concatenate(sv), np.concatenate(yv), np.concatenate(wv))
    print(f"{name:<10}{len(names):>5}{auc:>10.4f}{R[1]:>8.2f}{R[2]:>8.2f}{R[3]:>8.2f}  {'+'.join(keep)}")
