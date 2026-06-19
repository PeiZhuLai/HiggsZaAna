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
from hza_features import (_load, build_matrix, BASE_VARS, HZG_GGF_VARS,
                          _sideband_reweight, ROOT_DIR)

LOWM = [1, 2, 3]
FULL16 = BASE_VARS + ["pho_pt_asym", "param"]
# Current production low-mass set (model_Za_BDT_lowmass_run3.meta.json)
MIN5   = ["ALP_calculatedPhotonIso", "var_dR_g1g2", "pho1R9", "pho1Pt_oHm", "param"]
# HZgamma ggF-style groups (orthogonal to m_llgammagamma); jets excluded for HZa
ANGLES = ["Z_cos_theta", "lep_cos_theta", "lep_phi"]
DRGL   = ["min_lg_deltaR", "max_lg_deltaR"]
ETAS   = ["lead_lepton_eta", "sublead_lepton_eta", "gamma_eta_raw"]
ORTH   = ANGLES + DRGL + ETAS + ["H_ptt", "gamma_ptRelErr"]   # the sculpting-safe ggF additions

# Each variant is an explicit keep-list ('param' must be present; it is the per-event mass hyp).
VARIANTS = {
    "min5":            MIN5,
    "min5+angles":     MIN5[:-1] + ANGLES + ["param"],
    "min5+ggf_orth":   MIN5[:-1] + ORTH + ["param"],
    "full16":          FULL16,
    "full16+ggf_orth": BASE_VARS + ["pho_pt_asym"] + ORTH + ["param"],
    "full16+ggf_all":  BASE_VARS + ["pho_pt_asym"] + HZG_GGF_VARS + ["param"],
}
# Extra branches to read from the ROOT trees (the new ggF candidates).
GGF_READ = list(HZG_GGF_VARS)

def weighted_auc(vals, y, w):
    o = np.argsort(vals, kind="mergesort"); y = y[o]; w = w[o]; is_b = y < 0.5
    below_b = np.cumsum(w * is_b) - w * is_b
    return (w[~is_b] * below_b[~is_b]).sum() / (w[~is_b].sum() * w[is_b].sum() + 1e-12)

def wq(v, wt, q): o = np.argsort(v); c = np.cumsum(wt[o]) / wt.sum(); return np.interp(q, c, v[o])

# ---- load once: low-mass signal (train + inclusive), background (train + inclusive) ----
# extra=GGF_READ pulls in the new HZgamma ggF-style branches (require a Parque2Root_BDT.py re-run).
print("[load] signal mA=1,2,3 + All_Bkg (with ggF candidates) ...")
sig_tr = {m: _load(f"{ROOT_DIR}/mA_M{m}/run3.root", "train", 40000, extra=GGF_READ) for m in LOWM}
sig_in = {m: _load(f"{ROOT_DIR}/mA_M{m}/run3.root", "inclusive", extra=GGF_READ) for m in LOWM}
bk_tr = _load(f"{ROOT_DIR}/All_Bkg/run3.root", "train", 400000, extra=GGF_READ + ["Z_m", "event"])
bk_in = _load(f"{ROOT_DIR}/All_Bkg/run3.root", "inclusive", extra=GGF_READ)

ws = np.concatenate([np.clip(sig_tr[m]["factor"], 0, None) for m in LOWM])
bhyp = np.random.default_rng(1).choice(LOWM, size=len(bk_tr["H_m"]))   # per-event bkg mass hyp
wb = np.clip(bk_tr["factor"] * _sideband_reweight(bk_tr), 0, None)
y = np.concatenate([np.ones(len(ws)), np.zeros(len(wb))])
wsig = ws * (wb.sum() / ws.sum())
w = np.concatenate([wsig, wb])

Hmb = bk_in["H_m"]; wBe = np.clip(bk_in["factor"], 0, None)
peak = (Hmb > 120) & (Hmb < 130); frac = wBe[peak].sum() / wBe.sum()

print(f"\n{'variant':<18}{'nvar':>5}{'AUC(1-3)':>10}{'R(mA1)':>8}{'R(mA2)':>8}{'R(mA3)':>8}")
for name, keep in VARIANTS.items():
    Xs = np.vstack([build_matrix(sig_tr[m], m, keep) for m in LOWM])
    Xb = build_matrix(bk_tr, bhyp, keep)              # bhyp array -> per-event param column
    X = np.vstack([Xs, Xb])
    clf = xgb.XGBClassifier(max_depth=4, n_estimators=300, learning_rate=0.1, subsample=0.8,
                            min_child_weight=5, eval_metric="logloss", n_jobs=4, verbosity=0)
    clf.fit(X, y, sample_weight=w)
    # R(mA) on inclusive background
    R = {}
    for m in LOWM:
        s = clf.predict_proba(build_matrix(bk_in, m, keep))[:, 1]
        thr = wq(s, wBe, 1 - 0.002); p = s > thr
        R[m] = (wBe[p & peak].sum() / wBe[p].sum()) / frac if wBe[p].sum() > 0 else float("nan")
    # pooled AUC over masses
    sv, yv, wv = [], [], []
    for m in LOWM:
        ss = clf.predict_proba(build_matrix(sig_in[m], m, keep))[:, 1]
        sb = clf.predict_proba(build_matrix(bk_in, m, keep))[:, 1]
        wss = np.clip(sig_in[m]["factor"], 0, None)
        sv += [ss, sb]; yv += [np.ones(len(ss)), np.zeros(len(sb))]
        wv += [wss * (wBe.sum() / wss.sum()), wBe]
    auc = weighted_auc(np.concatenate(sv), np.concatenate(yv), np.concatenate(wv))
    print(f"{name:<18}{len(keep):>5}{auc:>10.4f}{R[1]:>8.2f}{R[2]:>8.2f}{R[3]:>8.2f}")
