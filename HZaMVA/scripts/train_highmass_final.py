#!/usr/bin/env python3
"""Formal dedicated low-mA (mA=1,2,3) XGBoost with the selected 5-feature set
(combinatorial-search winner: R(mA1,2,3)=0.99/0.95/1.02, AUC 0.9434).
CPU only (higgs-alp-ana). Trains on the 'train' tree, checks overtraining on the 'test' tree,
reports R(mA)/AUC, saves model (pkl + json) + feature metadata, and writes an overtraining PDF.

  ALP_calculatedPhotonIso + var_dR_g1g2 + pho1R9 + pho1Pt_oHm + param
"""
import json, pickle
import numpy as np, xgboost as xgb
from hza_features import _load, _feat, BASE_VARS, _sideband_reweight, signal_weight, ROOT_DIR

LOWM = [4,5,6,7,8,9,10,15,20,25,30]
FULL16 = BASE_VARS + ["pho_pt_asym", "param"]
LOWVARS = ["pho1Pt_oHm","pho1R9","pho1IetaIeta55","pho1PIso_noCorr","pho2Pt_oHm","pho2R9","pho2IetaIeta55","pho2PIso_noCorr","ALP_calculatedPhotonIso","var_dR_Za","var_dR_g1g2","var_dR_g1Z","var_PtaOverMh","H_pt_oHm","pho_pt_asym"]
FEATS = LOWVARS + ["param"]
IDX = [FULL16.index(n) for n in FEATS]
OUT = "model_Za_BDT_highmass_run3"
PLOT = "../plots_MVA/run3_highmass/overtrain_highmass_run3.pdf"

def weighted_auc(vals, y, w):
    o = np.argsort(vals, kind="mergesort"); y = y[o]; w = w[o]; is_b = y < 0.5
    below_b = np.cumsum(w * is_b) - w * is_b
    return (w[~is_b] * below_b[~is_b]).sum() / (w[~is_b].sum() * w[is_b].sum() + 1e-12)

def wq(v, wt, q): o = np.argsort(v); c = np.cumsum(wt[o]) / wt.sum(); return np.interp(q, c, v[o])

def weighted_ks(x1, w1, x2, w2):
    x1=np.asarray(x1,float); x2=np.asarray(x2,float); w1=np.abs(w1); w2=np.abs(w2)
    o1=np.argsort(x1); o2=np.argsort(x2); x1,w1=x1[o1],w1[o1]; x2,w2=x2[o2],w2[o2]
    c1=np.cumsum(w1)/w1.sum(); c2=np.cumsum(w2)/w2.sum(); xa=np.sort(np.concatenate([x1,x2]))
    i1=np.searchsorted(x1,xa,side="right")-1; i2=np.searchsorted(x2,xa,side="right")-1
    s=np.max(np.abs(np.where(i1>=0,c1[np.maximum(i1,0)],0)-np.where(i2>=0,c2[np.maximum(i2,0)],0)))
    n1=w1.sum()**2/(w1**2).sum(); n2=w2.sum()**2/(w2**2).sum(); en=np.sqrt(n1*n2/(n1+n2))
    from math import exp
    return s, 2*sum((-1)**(k-1)*exp(-2*k*k*(en*s)**2) for k in range(1,101))

def ks_capped(x1,w1,x2,w2,cap=4000,reps=25,seed=0):
    rng=np.random.default_rng(seed); ps=[]
    for _ in range(reps):
        i1=rng.choice(len(x1),min(cap,len(x1)),replace=False); i2=rng.choice(len(x2),min(cap,len(x2)),replace=False)
        ps.append(weighted_ks(x1[i1],w1[i1],x2[i2],w2[i2])[1])
    return float(np.median(ps))

def build(tree, n_sig=None, n_bkg=600000):
    sig = {m: _load(f"{ROOT_DIR}/mA_M{m}/run3.root", tree, n_sig, extra=["Z_m", "event"]) for m in LOWM}
    bk  = _load(f"{ROOT_DIR}/All_Bkg/run3.root", tree, n_bkg, extra=["Z_m", "event"])
    Xs = np.vstack([_feat(sig[m], m)[:, IDX] for m in LOWM])
    ws = np.concatenate([np.clip(signal_weight(sig[m], m), 0, None) for m in LOWM])  # signal: factor x sideband reweight
    bhyp = np.random.default_rng(1).choice(LOWM, size=len(bk["H_m"]))
    Xb = _feat(bk, 0); Xb[:, FULL16.index("param")] = (bk["ALP_m"] - bhyp) / bk["H_m"]; Xb = Xb[:, IDX]
    wb = np.clip(bk["factor"] * _sideband_reweight(bk), 0, None)
    X = np.vstack([Xs, Xb]); y = np.concatenate([np.ones(len(Xs)), np.zeros(len(Xb))])
    w = np.concatenate([ws * (wb.sum()/ws.sum()), wb])
    return X, y, w

print("[build] train ...")
Xtr, ytr, wtr = build("train")
clf = xgb.XGBClassifier(max_depth=4, n_estimators=300, learning_rate=0.1, subsample=0.8,
                        min_child_weight=5, eval_metric="logloss", n_jobs=4, verbosity=0)
clf.fit(Xtr, ytr, sample_weight=wtr)

# --- R(mA) on inclusive background, AUC + overtraining on test tree ---
bk_in = _load(f"{ROOT_DIR}/All_Bkg/run3.root", "inclusive")
Hmb = bk_in["H_m"]; wBe = np.clip(bk_in["factor"], 0, None)
peak = (Hmb > 120) & (Hmb < 130); frac = wBe[peak].sum()/wBe.sum()
R = {}
for m in LOWM:
    s = clf.predict_proba(_feat(bk_in, m)[:, IDX])[:, 1]; thr = wq(s, wBe, 1-0.002); p = s > thr
    R[m] = float((wBe[p&peak].sum()/wBe[p].sum())/frac) if wBe[p].sum()>0 else float("nan")

print("[build] test (overtraining) ...")
Xte, yte, wte = build("test")
str_ = clf.predict_proba(Xtr)[:, 1]; ste = clf.predict_proba(Xte)[:, 1]
# pooled test AUC
auc = weighted_auc(ste, yte, wte)
sig_tr, bk_tr = str_[ytr>0.5], str_[ytr<0.5]; sig_te, bk_te = ste[yte>0.5], ste[yte<0.5]
wst, wbt = wtr[ytr>0.5], wtr[ytr<0.5]; wse, wbe = wte[yte>0.5], wte[yte<0.5]
p_sig = ks_capped(sig_tr, wst, sig_te, wse, seed=1); p_bk = ks_capped(bk_tr, wbt, bk_te, wbe, seed=2)

print(f"\n[result] features={FEATS} (n={len(FEATS)})")
print("[result] R: " + "  ".join(f"mA{m}={R[m]:.2f}" for m in LOWM) + f"   AUC(test)={auc:.4f}")
print(f"[result] overtraining capped-KS p: sig={p_sig:.3g}  bkg={p_bk:.3g}")

# --- save model + metadata ---
pickle.dump(clf, open(f"{OUT}.pkl", "wb"))
clf.get_booster().save_model(f"{OUT}.json")
json.dump({"features": FEATS, "masses": LOWM, "R": R, "auc_test": auc,
           "ks_sig": p_sig, "ks_bkg": p_bk, "model": "XGB depth4/300/lr0.1"},
          open(f"{OUT}.meta.json", "w"), indent=2)
print(f"[save] {OUT}.pkl / .json / .meta.json")

# --- overtraining plot (best-effort; matplotlib/scipy import can be flaky on EOS) ---
try:
    import os; os.makedirs("../plots_MVA/run3_highmass", exist_ok=True)
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(8, 6)); rng = (0, 1); bins = 40
    plt.plot([], [], " ", label=f"Sig K-S p = {p_sig:.3g}")
    plt.plot([], [], " ", label=f"Bkg K-S p = {p_bk:.3g}")
    plt.hist(sig_tr, color="r", alpha=0.5, range=rng, bins=bins, density=True, weights=wst, histtype="stepfilled", label="Sig (Train)")
    plt.hist(bk_tr,  color="b", alpha=0.5, range=rng, bins=bins, density=True, weights=wbt, histtype="stepfilled", label="Bkg (Train)")
    for v, ww, c, lab in ((sig_te, wse, "r", "Sig (Test)"), (bk_te, wbe, "b", "Bkg (Test)")):
        h, e = np.histogram(v, bins=bins, range=rng, density=True, weights=ww); ctr = (e[:-1]+e[1:])/2
        err = np.sqrt(np.clip(h*len(v)/max(h.sum(),1e-9), 0, None))/(len(v)/max(h.sum(),1e-9))
        plt.errorbar(ctr, h, yerr=err, fmt=".", c=c, label=lab, markersize=8)
    plt.xlim(-0.1, 1.1)
    plt.xlabel(r"$m_{a} \geq 4$ GeV BDT score", fontsize=16); plt.ylabel("Arbitrary Units", fontsize=16)
    plt.legend(loc="upper center", fontsize=12, frameon=False); plt.tight_layout()
    fig.savefig(PLOT, dpi=200); print(f"[plot] {PLOT}")
except Exception as e:
    print("[plot] skipped (matplotlib/EOS):", e)
