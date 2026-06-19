#!/usr/bin/env python3
"""Regenerate ONLY the low-mass overtraining plot (train/test BDT score + capped-KS) from the
DEPLOYED model pkl -- NO retraining, so it matches bit-exact the model that produced the limits.
Mirrors the build()/ks/plot logic of train_lowmass_final.py. CPU (higgs-alp-ana).
"""
import json, pickle
import numpy as np
from hza_features import _load, _feat, BASE_VARS, _sideband_reweight, ROOT_DIR

_HZAMVA = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA"
MODEL = f"{_HZAMVA}/using/model_Za_BDT_lowmass_run3.pkl"
META  = f"{_HZAMVA}/using/model_Za_BDT_lowmass_run3.meta.json"
PLOT  = f"{_HZAMVA}/plots_MVA/run3_lowmass/overtrain_lowmass_run3.pdf"
LOWM  = [1, 2, 3]
FULL16 = BASE_VARS + ["pho_pt_asym", "param"]

meta = json.load(open(META)); FEATS = meta["features"]; IDX = [FULL16.index(n) for n in FEATS]
clf = pickle.load(open(MODEL, "rb"))
print(f"[cfg] features={FEATS}")

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
    sig = {m: _load(f"{ROOT_DIR}/mA_M{m}/run3.root", tree, n_sig) for m in LOWM}
    bk  = _load(f"{ROOT_DIR}/All_Bkg/run3.root", tree, n_bkg, extra=["Z_m", "event"])
    Xs = np.vstack([_feat(sig[m], m)[:, IDX] for m in LOWM])
    ws = np.concatenate([np.clip(sig[m]["factor"], 0, None) for m in LOWM])
    bhyp = np.random.default_rng(1).choice(LOWM, size=len(bk["H_m"]))
    Xb = _feat(bk, 0); Xb[:, FULL16.index("param")] = (bk["ALP_m"] - bhyp) / bk["H_m"]; Xb = Xb[:, IDX]
    wb = np.clip(bk["factor"] * _sideband_reweight(bk), 0, None)
    X = np.vstack([Xs, Xb]); y = np.concatenate([np.ones(len(Xs)), np.zeros(len(Xb))])
    w = np.concatenate([ws * (wb.sum()/ws.sum()), wb])
    return X, y, w

print("[build] train / test ...")
Xtr, ytr, wtr = build("train"); Xte, yte, wte = build("test")
str_ = clf.predict_proba(Xtr)[:, 1]; ste = clf.predict_proba(Xte)[:, 1]
sig_tr, bk_tr = str_[ytr>0.5], str_[ytr<0.5]; sig_te, bk_te = ste[yte>0.5], ste[yte<0.5]
wst, wbt = wtr[ytr>0.5], wtr[ytr<0.5]; wse, wbe = wte[yte>0.5], wte[yte<0.5]
p_sig = ks_capped(sig_tr, wst, sig_te, wse, seed=1); p_bk = ks_capped(bk_tr, wbt, bk_te, wbe, seed=2)
print(f"[result] capped-KS p: sig={p_sig:.3g}  bkg={p_bk:.3g}")

import os; os.makedirs(os.path.dirname(PLOT), exist_ok=True)
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

# ---- CMS-Preliminary label, copied verbatim from plot_mva_diagnostics.py so the style
# ---- (bold "CMS" + italic "Preliminary" + lumi) matches ROC.pdf exactly. ----
def add_cms_preliminary(fig, cms_fontsize=36):
    if fig.subplotpars.top > 0.88:
        fig.subplots_adjust(top=0.88)
    left = max(fig.subplotpars.left, 0.02)
    right = min(fig.subplotpars.right, 0.98)
    prelim_x = min(left + 2.7 * cms_fontsize / 72.0 / fig.get_figwidth(), right - 0.20)
    fig.text(left, 0.965, "CMS", fontsize=cms_fontsize, fontweight="bold", ha="left", va="top")
    fig.text(prelim_x, 0.965, "Preliminary", fontsize=cms_fontsize - 4, fontstyle="italic", ha="left", va="top")
    fig.text(right, 0.965, r"$172.13\ \mathrm{fb}^{-1}\ (13.6\ \mathrm{TeV})$",
             fontsize=cms_fontsize - 4, ha="right", va="top")

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
plt.xlabel(r"$m_{a} = 1, 2, 3$ GeV BDT score", fontsize=16); plt.ylabel("Arbitrary Units", fontsize=16)
plt.legend(loc="upper center", fontsize=12, frameon=False)
fig.tight_layout(rect=[0, 0, 1, 0.93])
add_cms_preliminary(fig, cms_fontsize=18)   # identical CMS-Preliminary style to ROC.pdf
fig.savefig(PLOT, dpi=200); print(f"[plot] {PLOT}")
