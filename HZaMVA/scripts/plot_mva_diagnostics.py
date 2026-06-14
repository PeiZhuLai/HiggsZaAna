#!/usr/bin/env python3
"""Standard MVA diagnostic plots for the dedicated low/high-mass BDTs, matching the plot types in
plots_MVA/run3 (the ones that come from a model, not the analysis pipeline):
  ROC, feature_importance, per-feature sig/bkg distributions, correlation matrices, sig_bkg_BDT score.
CPU (higgs-alp-ana). Loads a model .pkl + its .meta.json (features, masses) and writes PDFs to --outdir.
NOT generated (analysis/HPO-specific): optuna_*, cat*_H_mass, bkg_mass_shapes/smoothing, loss_vs_nEstimators.
"""
import argparse, json, os, pickle
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from hza_features import _load, _feat, BASE_VARS, ROOT_DIR

FULL16 = BASE_VARS + ["pho_pt_asym", "param"]

# ---- CMS-Preliminary label + axis formatting, mirroring run3_Za_BDT.py corr style ----
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

def format_corr_axis(ax, labels, tick_size, show_ylabels=True):
    pos = np.arange(len(labels)) + 0.5
    ax.set_xticks(pos); ax.set_xticklabels(labels, rotation=30, ha="right", rotation_mode="anchor")
    ax.set_yticks(pos)
    # force HORIZONTAL y tick labels (seaborn/MPL otherwise rotate the long feature names
    # to vertical and pile them on top of each other)
    ax.set_yticklabels(labels if show_ylabels else [], rotation=0, va="center")
    ax.set_xlabel(""); ax.set_ylabel("")
    ax.tick_params(axis="both", which="major", labelsize=tick_size)

def weighted_auc_curve(score, y, w):
    o = np.argsort(-score); y = y[o]; w = w[o]
    P = w[y > 0.5].sum(); N = w[y < 0.5].sum()
    tp = np.cumsum(np.where(y > 0.5, w, 0)) / P
    fp = np.cumsum(np.where(y < 0.5, w, 0)) / N
    auc = np.trapz(tp, fp)
    return fp, tp, auc

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True)
    ap.add_argument("--meta", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--n-sig-per", type=int, default=20000)
    ap.add_argument("--n-bkg", type=int, default=200000)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    clf = pickle.load(open(a.model, "rb"))
    meta = json.load(open(a.meta))
    feats = meta["features"]; masses = meta["masses"]
    idx = [FULL16.index(f) for f in feats]
    print(f"[cfg] model={a.model} feats={feats} masses={masses} -> {a.outdir}")

    # ---- build a balanced sig/bkg eval set for a given split (signal@true mass, bkg@random mass) ----
    def build_eval(tree):
        sig = {m: _load(f"{ROOT_DIR}/mA_M{m}/run3.root", tree, a.n_sig_per) for m in masses}
        bk = _load(f"{ROOT_DIR}/All_Bkg/run3.root", tree, a.n_bkg)
        Xs = np.vstack([_feat(sig[m], m)[:, idx] for m in masses])
        ws = np.concatenate([np.clip(sig[m]["factor"], 0, None) for m in masses])
        bhyp = np.random.default_rng(1).choice(masses, size=len(bk["H_m"]))
        Xb = _feat(bk, 0); Xb[:, FULL16.index("param")] = (bk["ALP_m"] - bhyp) / bk["H_m"]; Xb = Xb[:, idx]
        wb = np.clip(bk["factor"], 0, None)
        ss = clf.predict_proba(Xs)[:, 1]; sb = clf.predict_proba(Xb)[:, 1]
        wsb = ws * (wb.sum() / ws.sum())   # balance for ROC/score shapes
        return dict(sig=sig, bk=bk, Xs=Xs, Xb=Xb, ws=ws, wb=wb, ss=ss, sb=sb, wsb=wsb)

    def roc_of(ev):
        score = np.concatenate([ev["ss"], ev["sb"]])
        yy = np.concatenate([np.ones(len(ev["ss"])), np.zeros(len(ev["sb"]))])
        ww = np.concatenate([ev["wsb"], ev["wb"]])
        return weighted_auc_curve(score, yy, ww)

    ev_test = build_eval("test")
    ev_val = build_eval("validation")
    # downstream distributions / correlation / score plots use the test split
    sig, bk = ev_test["sig"], ev_test["bk"]
    Xs, Xb, ws, wb = ev_test["Xs"], ev_test["Xb"], ev_test["ws"], ev_test["wb"]
    ss, sb, wsb = ev_test["ss"], ev_test["sb"], ev_test["wsb"]

    # ---- ROC: Test + Validation curves ----
    # Match run3_Za_BDT.py convention: x = Signal Efficiency (TPR), y = Background Rejection (1-FPR).
    fp_t, tp_t, auc = roc_of(ev_test)
    fp_v, tp_v, auc_v = roc_of(ev_val)
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(tp_t, 1 - fp_t, color="darkgreen", lw=2, label=f"Test (AUC = {auc:.4f})")
    ax.plot(tp_v, 1 - fp_v, color="darkorange", lw=2, label=f"Validation (AUC = {auc_v:.4f})")
    ax.plot([1, 0], [0, 1], color="navy", lw=2, ls="--", label="Random Guess")
    ax.set_xlim(0.01, 1.1); ax.set_ylim(0.0, 1.1)
    ax.set_xlabel("Signal Efficiency", fontsize=20, labelpad=12)
    ax.set_ylabel("Background Rejection", fontsize=20, labelpad=12)
    ax.tick_params(axis="both", labelsize=16, direction="in", top=True, right=True, length=11)
    ax.legend(loc="lower left", fontsize=15, frameon=False)
    fig.subplots_adjust(left=0.15, right=0.96, top=0.91, bottom=0.15)
    add_cms_preliminary(fig, cms_fontsize=18)
    fig.savefig(f"{a.outdir}/ROC.pdf"); plt.close(fig)
    print(f"[plot] ROC.pdf (Test AUC={auc:.4f}, Validation AUC={auc_v:.4f})")

    # ---- feature_importance ----
    imp = clf.feature_importances_
    order = np.argsort(imp)
    fig, ax = plt.subplots(figsize=(8, 0.4 * len(feats) + 2))
    ax.barh(np.array(feats)[order], imp[order], color="tab:blue")
    ax.set_xlabel("XGBoost importance (gain-weighted)", fontsize=13); ax.tick_params(labelsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.93]); add_cms_preliminary(fig, cms_fontsize=18)
    fig.savefig(f"{a.outdir}/feature_importance.pdf"); plt.close(fig); print("[plot] feature_importance.pdf")

    # ---- sig_bkg_BDT (score distribution) ----
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.hist(ss, bins=50, range=(0, 1), weights=wsb, density=True, histtype="step", color="r", lw=2, label="Signal")
    ax.hist(sb, bins=50, range=(0, 1), weights=wb, density=True, histtype="step", color="b", lw=2, label="Background")
    ax.set_xlabel("BDT score", fontsize=14); ax.set_ylabel("Arbitrary units", fontsize=14)
    ax.set_yscale("log"); ax.legend(fontsize=13, frameon=False); ax.tick_params(direction="in", top=True, right=True)
    fig.tight_layout(); fig.savefig(f"{a.outdir}/sig_bkg_BDT.pdf"); plt.close(fig); print("[plot] sig_bkg_BDT.pdf")

    # ---- per-feature sig vs bkg distributions ----
    for j, fname in enumerate(feats):
        xs = Xs[:, j]; xb = Xb[:, j]
        lo, hi = np.percentile(np.concatenate([xs, xb]), [0.5, 99.5])
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.hist(xs, bins=50, range=(lo, hi), weights=ws, density=True, histtype="step", color="r", lw=2, label="Sig")
        ax.hist(xb, bins=50, range=(lo, hi), weights=wb, density=True, histtype="step", color="b", lw=2, ls="-.", label="Bkg")
        ax.set_xlabel(fname, fontsize=13); ax.set_ylabel("a.u.", fontsize=12)
        ax.legend(fontsize=12, frameon=False); ax.tick_params(direction="in", top=True, right=True)
        fig.tight_layout(); fig.savefig(f"{a.outdir}/{fname}.pdf"); plt.close(fig)
    print(f"[plot] {len(feats)} per-feature distributions")

    # ---- correlation matrices (run3_Za_BDT.py style; H_m / m_H appended as the LAST row/col) ----
    Hm_s = np.concatenate([sig[m]["H_m"] for m in masses])   # row order matches Xs (vstack over masses)
    Hm_b = bk["H_m"]
    df_sig = pd.DataFrame(Xs, columns=feats); df_sig["H_m"] = Hm_s
    df_bkg = pd.DataFrame(Xb, columns=feats); df_bkg["H_m"] = Hm_b
    corr_cols   = feats + ["H_m"]
    corr_labels = feats + [r"$m_{H}$"]
    n = len(corr_cols)
    ann = int(np.clip(360.0 / n, 12, 29))    # annotation font scales with matrix size
    tk  = int(np.clip(380.0 / n, 12, 30))    # tick-label font

    def draw_single(df, title, out):
        fig, ax = plt.subplots(figsize=(1.25 * n + 4, 1.15 * n + 3))
        sns.heatmap(df[corr_cols].corr(), annot=True, fmt=".2f", ax=ax, cmap="PiYG",
                    annot_kws={"size": ann}, cbar_kws={"pad": 0.01}, vmin=-1, vmax=1)
        ax.set_title(title, fontsize=int(ann * 1.25), fontstyle="italic", pad=18)
        format_corr_axis(ax, corr_labels, tk)
        cbar = ax.collections[0].colorbar
        if cbar is not None:
            cbar.ax.tick_params(labelsize=ann)
        fig.subplots_adjust(left=0.18, right=1.02, top=0.93, bottom=0.16)
        add_cms_preliminary(fig)
        fig.savefig(out); plt.close(fig)

    draw_single(df_sig, "Signal", f"{a.outdir}/corr_Sig.pdf")
    draw_single(df_bkg, "Background", f"{a.outdir}/corr_Bkg.pdf")

    # side-by-side Signal | Background (matches run3 corr_Sig_Bkg.pdf layout).
    # square=True + equal-size panels + ONE shared colorbar => both matrices have the
    # exact same aspect ratio; the colorbar steals equally from both so neither is skewed.
    # Fonts are a touch smaller than the single plots so the per-cell text never collides,
    # and scale with n so the 17-var (high-mass) and 7-var (low-mass) plots both read well.
    maxlen = max(len(l) for l in corr_labels)
    ann_c = int(np.clip(300.0 / n, 10, 22))   # cell annotation
    tk_c  = int(np.clip(320.0 / n, 10, 24))   # tick labels
    cell  = 1.05                               # inches per square cell
    left_in = 0.60 * maxlen * tk_c / 72.0 + 0.4    # room for the longest horizontal y label
    figw = left_in + 2.0 * cell * n + 0.6 + 3.0    # labels + 2 panels + gap + colorbar strip
    figh = cell * n + 4.6                          # panel + CMS/title top + rotated-xlabel bottom
    f, (ax1, ax2) = plt.subplots(1, 2, figsize=(figw, figh))
    for ax, df, title in ((ax1, df_sig, "Signal"), (ax2, df_bkg, "Background")):
        sns.heatmap(df[corr_cols].corr(), annot=True, fmt=".2f", ax=ax, cmap="PiYG",
                    annot_kws={"size": ann_c}, vmin=-1, vmax=1, cbar=False, square=True)
        ax.set_title(title, fontsize=int(ann_c * 1.5), fontstyle="italic", pad=14)
    format_corr_axis(ax1, corr_labels, tk_c)
    format_corr_axis(ax2, corr_labels, tk_c, show_ylabels=False)
    # panels occupy [left, 0.85]; the colorbar lives in its own axes in the reserved
    # right strip so its tick labels never get clipped and the two panels stay identical.
    f.subplots_adjust(left=min(0.36, left_in / figw), right=0.85,
                      top=0.88, bottom=0.16, wspace=0.06)
    sm = plt.cm.ScalarMappable(cmap="PiYG", norm=plt.Normalize(-1, 1)); sm.set_array([])
    cax = f.add_axes([0.875, 0.24, 0.016, 0.46])   # [x, y, w, h] in figure fraction
    cbar = f.colorbar(sm, cax=cax)
    cbar.ax.tick_params(labelsize=tk_c)
    add_cms_preliminary(f, cms_fontsize=int(ann_c * 1.6))
    f.savefig(f"{a.outdir}/corr_Sig_Bkg.pdf"); plt.close(f)
    print("[plot] corr_Sig/Bkg/Sig_Bkg.pdf (incl. m_H as last row, run3 style)")
    print(f"[done] AUC={auc:.4f}  -> {a.outdir}")

if __name__ == "__main__":
    main()
