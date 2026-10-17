#!/usr/bin/env python3
"""Formal dedicated low-mA (mA=1,2,3) XGBoost = low_relE2_w123: the low_swapPt_relE2 inputs (input-variable search
winner), with the signal of every low mass m_a = 1, 2, 3 trained inside the window (one rule for the three points).
CPU only (higgs-alp-ana). Trains on the 'train' tree, checks overtraining on the 'test' tree,
reports R(mA)/AUC, saves model (pkl + json) + feature metadata, and writes an overtraining PDF.

  ALP_calculatedPhotonIso + var_dR_g1g2 + pho1R9 + pho2_relEerr + param
  pho2_relEerr = sigma_E/E of the sub-leading photon (hza_features.feature_matrix); param = (m_gg - m_a)/m_llgg

Training set (same events, weights and recipe as the search model, varsearch_validate/train_selected.build with
win_masses 2, 3, model low_relE2_w123):
  - signal m_a = 1,2,3: 'train' events inside WINDOW (115 < m_llgg < 133) only, weight = factor x sideband reweight
    (true-mass param), clipped at 0, rescaled per mass to the full weight of that mass (WINDOW_MASSES = 1, 2, 3)
  - background: 600k 'train' events (default_rng(0)), weight = factor x sideband reweight, mass hypothesis
    default_rng(1) over 1,2,3; signal normalized to the background sum
  - XGB depth 4, 300 trees, lr 0.1, subsample 0.8, min_child_weight 5, n_jobs 1 (as for the deployed model)
After any change here, retrain to a side path (HZA_LOWMODEL_OUT) and check that the booster equals the installed
one before installing. FEATS must equal LOW_FEATURE_COLUMNS in Parquet2Rootfile/Parque2Root_BDT.py:
install_lowhigh_models.sh and the scorer refuse a model whose meta features differ.
The AUC/KS quoted in the AN come from plot_mva_diagnostics.py / plot_overtrain_lowmass.py on the installed set;
the numbers printed here are a quick check on the training definition (m_a 1-3 windowed in train and test).
Output: <stem>.pkl / .json / .meta.json with stem = HZA_LOWMODEL_OUT (default model_Za_BDT_lowmass_run3 in the cwd),
        overtraining PDF in <dir of stem>/plots_lowmass/
"""
import json, os, pickle
import numpy as np, xgboost as xgb
from hza_features import (_load, _sideband_reweight, signal_weight, feature_matrix, relerr_extra,
                          READ, REWEIGHT_EXTRA, ROOT_DIR)

LOWM = [1, 2, 3]
FEATS = ["ALP_calculatedPhotonIso", "var_dR_g1g2", "pho1R9", "pho2_relEerr", "param"]
WINDOW = (115.0, 133.0)
WINDOW_MASSES = [1, 2, 3]
EXTRA = REWEIGHT_EXTRA + relerr_extra(FEATS)
OUT = os.environ.get("HZA_LOWMODEL_OUT", "model_Za_BDT_lowmass_run3")
PLOT = os.path.join(os.path.dirname(OUT) or ".", "plots_lowmass", "overtrain_lowmass_run3.pdf")

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

def reweight_view(a):
    """Only the branches the sideband reweighter uses (same frame as the search training)."""
    return {k: a[k] for k in READ + REWEIGHT_EXTRA}

def build(tree, n_bkg=600000):
    Xs, ws = [], []
    for m in LOWM:
        a = _load(f"{ROOT_DIR}/mA_M{m}/run3.root", tree, None, extra=EXTRA)
        w = np.clip(signal_weight(reweight_view(a), m), 0, None)   # signal: factor x sideband reweight
        x = feature_matrix(a, float(m), FEATS)
        if m in WINDOW_MASSES:
            keep = (a["H_m"] > WINDOW[0]) & (a["H_m"] < WINDOW[1])
            print(f"  [{tree}] mA{m} window {WINDOW}: keep {keep.sum()}/{len(keep)}, {w[keep].sum()/w.sum():.3f} of the weight")
            x, w = x[keep], w[keep] * (w.sum() / w[keep].sum())
        Xs.append(x); ws.append(w)
    Xs = np.vstack(Xs); ws = np.concatenate(ws)
    bk = _load(f"{ROOT_DIR}/All_Bkg/run3.root", tree, n_bkg, extra=EXTRA)
    bhyp = np.random.default_rng(1).choice(LOWM, size=len(bk["H_m"])).astype(np.float64)
    Xb = feature_matrix(bk, bhyp, FEATS)
    wb = np.clip(bk["factor"] * _sideband_reweight(reweight_view(bk)), 0, None)
    X = np.vstack([Xs, Xb]); y = np.concatenate([np.ones(len(Xs)), np.zeros(len(Xb))])
    w = np.concatenate([ws * (wb.sum()/ws.sum()), wb])
    return X, y, w

print(f"[cfg] features={FEATS}  window={WINDOW} for mA {WINDOW_MASSES}  inputs={ROOT_DIR}")
print("[build] train ...")
Xtr, ytr, wtr = build("train")
print(f"  X {Xtr.shape} dtype {Xtr.dtype}, {int(ytr.sum())} signal / {int((1-ytr).sum())} background")
clf = xgb.XGBClassifier(max_depth=4, n_estimators=300, learning_rate=0.1, subsample=0.8,
                        min_child_weight=5, eval_metric="logloss", n_jobs=1, verbosity=0)
clf.fit(Xtr, ytr, sample_weight=wtr)

# --- R(mA) on inclusive background, AUC + overtraining on test tree ---
bk_in = _load(f"{ROOT_DIR}/All_Bkg/run3.root", "inclusive", extra=relerr_extra(FEATS))
Hmb = bk_in["H_m"]; wBe = np.clip(bk_in["factor"], 0, None)
peak = (Hmb > 120) & (Hmb < 130); frac = wBe[peak].sum()/wBe.sum()
R = {}
for m in LOWM:
    s = clf.predict_proba(feature_matrix(bk_in, float(m), FEATS))[:, 1]; thr = wq(s, wBe, 1-0.002); p = s > thr
    R[m] = float((wBe[p&peak].sum()/wBe[p].sum())/frac) if wBe[p].sum()>0 else float("nan")

print("[build] test (overtraining) ...")
Xte, yte, wte = build("test")
str_ = clf.predict_proba(Xtr)[:, 1]; ste = clf.predict_proba(Xte)[:, 1]
# pooled test AUC
auc = weighted_auc(ste, yte, wte)
sig_tr, bk_tr = str_[ytr>0.5], str_[ytr<0.5]; sig_te, bk_te = ste[yte>0.5], ste[yte<0.5]
wst, wbt = wtr[ytr>0.5], wtr[ytr<0.5]; wse, wbe = wte[yte>0.5], wte[yte<0.5]
p_sig = ks_capped(sig_tr, wst, sig_te, wse, seed=1); p_bk = ks_capped(bk_tr, wbt, bk_te, wbe, seed=2)

print(f"\n[result] features={FEATS}")
print(f"[result] R(mA1)={R[1]:.2f} R(mA2)={R[2]:.2f} R(mA3)={R[3]:.2f}  AUC(test,1-3)={auc:.4f}")
print(f"[result] overtraining capped-KS p: sig={p_sig:.3g}  bkg={p_bk:.3g}")

# --- save model + metadata (the three files are written together: install_lowhigh_models.sh checks their mtimes) ---
pickle.dump(clf, open(f"{OUT}.pkl", "wb"))
clf.get_booster().save_model(f"{OUT}.json")
json.dump({"features": FEATS, "masses": LOWM, "window": list(WINDOW), "window_masses": WINDOW_MASSES,
           "mA1_window": list(WINDOW), "name": "low_relE2_w123",
           # auc_test / ks_sig / ks_bkg (AN numbers) are added from plot_mva_diagnostics.py / plot_overtrain_lowmass.py
           "R": R, "check_auc_test_window": auc, "check_ks_sig_window": p_sig, "check_ks_bkg_window": p_bk,
           "model": "XGB depth4/300/lr0.1/subsample0.8/mcw5", "inputs": ROOT_DIR, "dtype": "float64",
           "xgboost": xgb.__version__, "training": "HZaMVA/scripts/train_lowmass_final.py"},
          open(f"{OUT}.meta.json", "w"), indent=2)
print(f"[save] {OUT}.pkl / .json / .meta.json")

# --- overtraining plot (best-effort; matplotlib/scipy import can be flaky on EOS) ---
try:
    os.makedirs(os.path.dirname(PLOT), exist_ok=True)
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(8, 6)); rng = (0, 1); bins = 40
    plt.hist(0, color="w", label=f"Sig K-S p = {p_sig:.3g} (n=4000)")
    plt.hist(0, color="w", label=f"Bkg K-S p = {p_bk:.3g} (n=4000)")
    plt.hist(sig_tr, color="r", alpha=0.5, range=rng, bins=bins, density=True, weights=wst, histtype="stepfilled", label="Sig (Train)")
    plt.hist(bk_tr,  color="b", alpha=0.5, range=rng, bins=bins, density=True, weights=wbt, histtype="stepfilled", label="Bkg (Train)")
    for v, ww, c, lab in ((sig_te, wse, "r", "Sig (Test)"), (bk_te, wbe, "b", "Bkg (Test)")):
        h, e = np.histogram(v, bins=bins, range=rng, density=True, weights=ww); ctr = (e[:-1]+e[1:])/2
        err = np.sqrt(np.clip(h*len(v)/max(h.sum(),1e-9), 0, None))/(len(v)/max(h.sum(),1e-9))
        plt.errorbar(ctr, h, yerr=err, fmt=".", c=c, label=lab, markersize=8)
    plt.xlabel("low-mass BDT score", fontsize=16); plt.ylabel("Arbitrary Units", fontsize=16)
    plt.legend(loc="upper center", fontsize=12, frameon=False); plt.tight_layout()
    fig.savefig(PLOT, dpi=200); print(f"[plot] {PLOT}")
except Exception as e:
    print("[plot] skipped (matplotlib/EOS):", e)
