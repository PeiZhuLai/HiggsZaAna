#!/usr/bin/env python3
"""A/B test: does dropping `pho_pt_asym` from the HIGH-mass model lower R(m_a) (esp. m_a=30=1.21)
without costing AUC? Trains BASELINE (full 16, = deployed) and REDUCED (drop pho_pt_asym) under
identical data/config/seed, prints R(m_a) for all high masses + test AUC side by side.
Does NOT touch the deployed model. Runs in higgs-alp-ana.
"""
import numpy as np, xgboost as xgb
from hza_features import _load, _feat, BASE_VARS, _sideband_reweight, ROOT_DIR

HM = [4,5,6,7,8,9,10,15,20,25,30]
FULL16 = BASE_VARS + ["pho_pt_asym", "param"]
BASE_FEATS = [f for f in FULL16]                                  # full 16 (deployed)
RED_FEATS  = [f for f in FULL16 if f != "pho_pt_asym"]            # 15 (drop pho_pt_asym)

def weighted_auc(vals, y, w):
    o = np.argsort(vals, kind="mergesort"); y = y[o]; w = w[o]; is_b = y < 0.5
    below_b = np.cumsum(w * is_b) - w * is_b
    return (w[~is_b] * below_b[~is_b]).sum() / (w[~is_b].sum() * w[is_b].sum() + 1e-12)

def wq(v, wt, q): o = np.argsort(v); c = np.cumsum(wt[o]) / wt.sum(); return np.interp(q, c, v[o])

# build once with FULL feature matrix; slice per feature-set via index lists
def build(tree, n_bkg=600000):
    sig = {m: _load(f"{ROOT_DIR}/mA_M{m}/run3.root", tree, None) for m in HM}
    bk  = _load(f"{ROOT_DIR}/All_Bkg/run3.root", tree, n_bkg, extra=["Z_m", "event"])
    Xs = np.vstack([_feat(sig[m], m) for m in HM])
    ws = np.concatenate([np.clip(sig[m]["factor"], 0, None) for m in HM])
    bhyp = np.random.default_rng(1).choice(HM, size=len(bk["H_m"]))
    Xb = _feat(bk, 0); Xb[:, FULL16.index("param")] = (bk["ALP_m"] - bhyp) / bk["H_m"]
    wb = np.clip(bk["factor"] * _sideband_reweight(bk), 0, None)
    X = np.vstack([Xs, Xb]); y = np.concatenate([np.ones(len(Xs)), np.zeros(len(Xb))])
    w = np.concatenate([ws * (wb.sum()/ws.sum()), wb])
    return X, y, w

print("[build] train / test (shared) ...")
Xtr, ytr, wtr = build("train"); Xte, yte, wte = build("test")
bk_in = _load(f"{ROOT_DIR}/All_Bkg/run3.root", "inclusive")
Hmb = bk_in["H_m"]; wBe = np.clip(bk_in["factor"], 0, None)
peak = (Hmb > 120) & (Hmb < 130); frac = wBe[peak].sum()/wBe.sum()
Xb_in_full = _feat(bk_in, 0)
# precompute per-mass inclusive feature matrices (param overwritten per mass)
def incl_X(m):
    Xm = Xb_in_full.copy(); Xm[:, FULL16.index("param")] = (bk_in["ALP_m"] - m) / bk_in["H_m"]; return Xm

def run(feats, tag):
    idx = [FULL16.index(f) for f in feats]
    clf = xgb.XGBClassifier(max_depth=4, n_estimators=300, learning_rate=0.1, subsample=0.8,
                            min_child_weight=5, eval_metric="logloss", n_jobs=4, verbosity=0)
    clf.fit(Xtr[:, idx], ytr, sample_weight=wtr)
    auc = weighted_auc(clf.predict_proba(Xte[:, idx])[:, 1], yte, wte)
    R = {}
    for m in HM:
        s = clf.predict_proba(incl_X(m)[:, idx])[:, 1]; thr = wq(s, wBe, 1-0.002); p = s > thr
        R[m] = float((wBe[p&peak].sum()/wBe[p].sum())/frac) if wBe[p].sum()>0 else float("nan")
    print(f"[{tag}] n_feat={len(feats)}  AUC(test)={auc:.4f}")
    return R, auc

print("\n=== training BASELINE (full 16) ===");  Rb, ab = run(BASE_FEATS, "baseline")
print("=== training REDUCED (drop pho_pt_asym, 15) ===");  Rr, ar = run(RED_FEATS, "reduced")

print("\n############ A/B RESULT: R(m_a) baseline vs reduced ############")
print(f"{'m_a':>4} {'R_base':>8} {'R_red':>8} {'dR':>7}")
for m in HM:
    print(f"{m:>4} {Rb[m]:>8.3f} {Rr[m]:>8.3f} {Rr[m]-Rb[m]:>+7.3f}")
print(f"\nAUC:  baseline={ab:.4f}  reduced={ar:.4f}  dAUC={ar-ab:+.4f}")
maxb=max(Rb.values()); maxr=max(Rr.values())
print(f"max R: baseline={maxb:.3f} (mA{max(Rb,key=Rb.get)})  reduced={maxr:.3f} (mA{max(Rr,key=Rr.get)})")
print(f"R(mA30): baseline={Rb[30]:.3f}  reduced={Rr[30]:.3f}")
