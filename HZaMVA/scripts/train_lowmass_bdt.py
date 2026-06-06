#!/usr/bin/env python3
"""Dedicated low-mA (mA=1,2,3) XGBoost with reduced input variables, to test whether dropping
the diagnosed sculpting-driver variables lowers R(mA) without losing much AUC.
CPU only (runs in higgs-alp-ana) -- avoids the GPU/condor EOS flakiness. Data is loaded ONCE
and reused across variants. Compares full(16) vs progressively reduced feature sets.

variants (diagnosed in diag_lowma_vars.py):
  full         : all 16
  conservative : drop 3 pure liabilities (sep~0.5 at mA1 but corr w/ m_llgg)
  middle       : + drop 3 sculpting engines (high corr, moderate sep)
  clean        : + drop 2 high-corr photon isolations
"""
import numpy as np, xgboost as xgb
from hza_features import _load, _feat, BASE_VARS, _sideband_reweight, ROOT_DIR

LOWM = [1, 2, 3]
FULL16 = BASE_VARS + ["pho_pt_asym", "param"]
DROP = {
    "full":         [],
    "conservative": ["pho1IetaIeta55", "pho2IetaIeta55", "pho_pt_asym"],
    "middle":       ["pho1IetaIeta55", "pho2IetaIeta55", "pho_pt_asym", "var_dR_Za", "var_dR_g1Z", "pho2Pt_oHm"],
    "clean":        ["pho1IetaIeta55", "pho2IetaIeta55", "pho_pt_asym", "var_dR_Za", "var_dR_g1Z", "pho2Pt_oHm",
                     "pho1PIso_noCorr", "pho2PIso_noCorr"],
}

def cols_for(drop):
    keep = [n for n in FULL16 if n not in drop]
    assert "param" in keep
    return keep, [FULL16.index(n) for n in keep]

def weighted_auc(vals, y, w):
    o = np.argsort(vals, kind="mergesort"); y = y[o]; w = w[o]; is_b = y < 0.5
    below_b = np.cumsum(w * is_b) - w * is_b
    return (w[~is_b] * below_b[~is_b]).sum() / (w[~is_b].sum() * w[is_b].sum() + 1e-12)

def wq(v, wt, q): o = np.argsort(v); c = np.cumsum(wt[o]) / wt.sum(); return np.interp(q, c, v[o])

# ---- load once: low-mass signal (train + inclusive), background (train + inclusive) ----
print("[load] signal mA=1,2,3 + All_Bkg ...")
sig_tr = {m: _load(f"{ROOT_DIR}/mA_M{m}/run3.root", "train", 40000) for m in LOWM}
sig_in = {m: _load(f"{ROOT_DIR}/mA_M{m}/run3.root", "inclusive") for m in LOWM}
bk_tr = _load(f"{ROOT_DIR}/All_Bkg/run3.root", "train", 400000, extra=["Z_m", "event"])
bk_in = _load(f"{ROOT_DIR}/All_Bkg/run3.root", "inclusive")

Xs16 = np.vstack([_feat(sig_tr[m], m) for m in LOWM])
ws = np.concatenate([np.clip(sig_tr[m]["factor"], 0, None) for m in LOWM])
bhyp = np.random.default_rng(1).choice(LOWM, size=len(bk_tr["H_m"]))
Xb16 = _feat(bk_tr, 0); Xb16[:, FULL16.index("param")] = (bk_tr["ALP_m"] - bhyp) / bk_tr["H_m"]
wb = np.clip(bk_tr["factor"] * _sideband_reweight(bk_tr), 0, None)
y = np.concatenate([np.ones(len(Xs16)), np.zeros(len(Xb16))])
wsig = ws * (wb.sum() / ws.sum())
w = np.concatenate([wsig, wb])
X16 = np.vstack([Xs16, Xb16])

Hmb = bk_in["H_m"]; wBe = np.clip(bk_in["factor"], 0, None)
peak = (Hmb > 120) & (Hmb < 130); frac = wBe[peak].sum() / wBe.sum()
Xb_in16 = {m: _feat(bk_in, m) for m in LOWM}
# AUC eval: inclusive signal@true vs bkg@same hypothesis, pooled over masses
Xs_in16 = {m: _feat(sig_in[m], m) for m in LOWM}

print(f"\n{'variant':<14}{'nvar':>5}{'AUC(1-3)':>10}{'R(mA1)':>8}{'R(mA2)':>8}{'R(mA3)':>8}")
for name, drop in DROP.items():
    keep, idx = cols_for(drop)
    clf = xgb.XGBClassifier(max_depth=4, n_estimators=300, learning_rate=0.1, subsample=0.8,
                            min_child_weight=5, eval_metric="logloss", n_jobs=4, verbosity=0)
    clf.fit(X16[:, idx], y, sample_weight=w)
    # R(mA) on inclusive background
    R = {}
    for m in LOWM:
        s = clf.predict_proba(Xb_in16[m][:, idx])[:, 1]
        thr = wq(s, wBe, 1 - 0.002); p = s > thr
        R[m] = (wBe[p & peak].sum() / wBe[p].sum()) / frac if wBe[p].sum() > 0 else float("nan")
    # pooled AUC over masses
    sv, yv, wv = [], [], []
    for m in LOWM:
        ss = clf.predict_proba(Xs_in16[m][:, idx])[:, 1]; sb = clf.predict_proba(Xb_in16[m][:, idx])[:, 1]
        wss = np.clip(sig_in[m]["factor"], 0, None)
        sv += [ss, sb]; yv += [np.ones(len(ss)), np.zeros(len(sb))]
        wv += [wss * (wBe.sum() / wss.sum()), wBe]
    auc = weighted_auc(np.concatenate(sv), np.concatenate(yv), np.concatenate(wv))
    print(f"{name:<14}{len(keep):>5}{auc:>10.4f}{R[1]:>8.2f}{R[2]:>8.2f}{R[3]:>8.2f}")
