#!/usr/bin/env python3
"""Regenerate the BDT overtraining PDF from a SAVED model (no retraining).

Loads model_*.pkl + its .meta.json, rebuilds the train/test sets via hza_features, scores them,
and replots the train-vs-test BDT score overlay. KS p-values are read straight from the meta so
the plot stays consistent with the stored model. Cosmetic fixes vs the in-training plot:
  - x-axis starts at -0.1 (the old plt.hist(0,...) legend dummies pulled it to ~-0.5)
  - no "(n=4000)" in the legend
  - x-axis label is passed in via --xlabel (e.g. per-mass-range wording)
Runs in higgs-alp-ana (xgboost).

  python3 replot_overtrain_bdt.py --model model_Za_BDT_highmass_run3.pkl \
      --meta model_Za_BDT_highmass_run3.meta.json \
      --out ../plots_MVA/run3_highmass/overtrain_highmass_run3.pdf \
      --xlabel '$m_{a} \geq 4$ GeV BDT score'
"""
import argparse, json, pickle
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from hza_features import _load, _feat, BASE_VARS, _sideband_reweight, ROOT_DIR

FULL16 = BASE_VARS + ["pho_pt_asym", "param"]

# ---- CMS-Preliminary label, mirroring plot_mva_diagnostics.add_cms_preliminary ----
def add_cms_preliminary(fig, cms_fontsize=18):
    if fig.subplotpars.top > 0.88:
        fig.subplots_adjust(top=0.88)
    left = max(fig.subplotpars.left, 0.02)
    right = min(fig.subplotpars.right, 0.98)
    prelim_x = min(left + 2.7 * cms_fontsize / 72.0 / fig.get_figwidth(), right - 0.20)
    fig.text(left, 0.965, "CMS", fontsize=cms_fontsize, fontweight="bold", ha="left", va="top")
    fig.text(prelim_x, 0.965, "Preliminary", fontsize=cms_fontsize - 4, fontstyle="italic", ha="left", va="top")
    fig.text(right, 0.965, r"$172.13\ \mathrm{fb}^{-1}\ (13.6\ \mathrm{TeV})$",
             fontsize=cms_fontsize - 4, ha="right", va="top")

def build(masses, IDX, tree, n_sig=None, n_bkg=600000):
    """Mirror train_{high,low}mass_final.build(): balanced sig/bkg feature set for one split."""
    sig = {m: _load(f"{ROOT_DIR}/mA_M{m}/run3.root", tree, n_sig) for m in masses}
    bk = _load(f"{ROOT_DIR}/All_Bkg/run3.root", tree, n_bkg, extra=["Z_m", "event"])
    Xs = np.vstack([_feat(sig[m], m)[:, IDX] for m in masses])
    ws = np.concatenate([np.clip(sig[m]["factor"], 0, None) for m in masses])
    bhyp = np.random.default_rng(1).choice(masses, size=len(bk["H_m"]))
    Xb = _feat(bk, 0); Xb[:, FULL16.index("param")] = (bk["ALP_m"] - bhyp) / bk["H_m"]; Xb = Xb[:, IDX]
    wb = np.clip(bk["factor"] * _sideband_reweight(bk), 0, None)
    X = np.vstack([Xs, Xb]); y = np.concatenate([np.ones(len(Xs)), np.zeros(len(Xb))])
    w = np.concatenate([ws * (wb.sum() / ws.sum()), wb])
    return X, y, w

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True)
    ap.add_argument("--meta", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--xlabel", required=True)
    ap.add_argument("--n-bkg", type=int, default=600000)
    a = ap.parse_args()

    clf = pickle.load(open(a.model, "rb"))
    meta = json.load(open(a.meta))
    masses = meta["masses"]; feats = meta["features"]
    IDX = [FULL16.index(f) for f in feats]
    p_sig = meta.get("ks_sig", float("nan")); p_bk = meta.get("ks_bkg", float("nan"))
    print(f"[cfg] model={a.model} masses={masses} feats={feats}")

    Xtr, ytr, wtr = build(masses, IDX, "train", n_bkg=a.n_bkg)
    Xte, yte, wte = build(masses, IDX, "test", n_bkg=a.n_bkg)
    str_ = clf.predict_proba(Xtr)[:, 1]; ste = clf.predict_proba(Xte)[:, 1]
    sig_tr, bk_tr = str_[ytr > 0.5], str_[ytr < 0.5]; sig_te, bk_te = ste[yte > 0.5], ste[yte < 0.5]
    wst, wbt = wtr[ytr > 0.5], wtr[ytr < 0.5]; wse, wbe = wte[yte > 0.5], wte[yte < 0.5]

    fig = plt.figure(figsize=(8, 6)); rng = (0, 1); bins = 40
    # legend-only KS text (no dummy histogram, so the x-range stays clean)
    plt.plot([], [], " ", label=f"Sig K-S p = {p_sig:.3g}")
    plt.plot([], [], " ", label=f"Bkg K-S p = {p_bk:.3g}")
    plt.hist(sig_tr, color="r", alpha=0.5, range=rng, bins=bins, density=True, weights=wst, histtype="stepfilled", label="Sig (Train)")
    plt.hist(bk_tr,  color="b", alpha=0.5, range=rng, bins=bins, density=True, weights=wbt, histtype="stepfilled", label="Bkg (Train)")
    for v, ww, c, lab in ((sig_te, wse, "r", "Sig (Test)"), (bk_te, wbe, "b", "Bkg (Test)")):
        h, e = np.histogram(v, bins=bins, range=rng, density=True, weights=ww); ctr = (e[:-1] + e[1:]) / 2
        err = np.sqrt(np.clip(h * len(v) / max(h.sum(), 1e-9), 0, None)) / (len(v) / max(h.sum(), 1e-9))
        plt.errorbar(ctr, h, yerr=err, fmt=".", c=c, label=lab, markersize=8)
    plt.xlim(-0.1, 1.1)
    plt.xlabel(a.xlabel, fontsize=16); plt.ylabel("Arbitrary Units", fontsize=16)
    plt.legend(loc="upper center", fontsize=12, frameon=False)
    fig.tight_layout(rect=[0, 0, 1, 0.93]); add_cms_preliminary(fig, cms_fontsize=18)
    fig.savefig(a.out, dpi=200)
    print(f"[plot] {a.out}  (KS p from meta: sig={p_sig:.3g} bkg={p_bk:.3g})")

if __name__ == "__main__":
    main()
