#!/usr/bin/env python3
"""Compare the per-mass-DisCo NN against the production XGBoost BDT on the SAME events.

Runs in the BDT's native env (higgs-alp-ana, xgboost 2.1.1) -- NO torch, does NOT touch hza_NN.
  - NN scores  : read from the step-(b) slim files (run3_nn_scored_inputs_nominal/<sample>),
                 produced earlier in hza_NN. (So this also validates that pipeline.)
  - BDT scores : computed here from the bdt_inputs features via model_Za_BDT_run3.pkl.
Both files come from the SAME inclusive trees in the same order; a single H_m in [95,180]
mask is applied to both so the events line up exactly.

Metrics: weighted AUC (parametrized, signal at true hypothesis vs bkg) and sculpting
R(mA) = post-cut 125-peak bkg fraction / pre-cut at fixed bkg efficiency. Writes a PDF + table.
"""
import argparse, pickle
import numpy as np, uproot
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from hza_features import BASE_VARS, _feat, ROOT_DIR, MASSES

def roc_auc_score(y, vals, sample_weight):
    """Weighted AUC via rank formula (no sklearn; EOS .so imports are flaky tonight)."""
    o = np.argsort(vals, kind="mergesort"); y = np.asarray(y)[o]; w = np.asarray(sample_weight)[o]
    is_b = y < 0.5; below_b = np.cumsum(w * is_b) - w * is_b
    return (w[~is_b] * below_b[~is_b]).sum() / (w[~is_b].sum() * w[is_b].sum() + 1e-12)

NN_SCORED   = "/eos/home-p/pelai/HZa/root_P2Root/run3_nn_scored_inputs_nominal"
EVAL_MASSES = [1, 2, 3, 4, 5, 10, 20]
FEAT_READ   = BASE_VARS + ["H_m", "ALP_m", "factor"]

def load_sample(sample):
    """Return (feat_dict, H_m, factor, nn_tree) for one sample's inclusive tree, with a common
    H_m in [95,180] & finite-factor mask applied to BOTH the bdt_inputs features and the NN file."""
    bdf = uproot.open(f"{ROOT_DIR}/{sample}/run3.root")["inclusive"].arrays(FEAT_READ, library="np")
    nnf = uproot.open(f"{NN_SCORED}/{sample}/run3.root")["inclusive"]
    Hm = bdf["H_m"]; fac = bdf["factor"]
    mask = (Hm > 95) & (Hm < 180) & np.isfinite(fac)
    feat = {b: bdf[b][mask] for b in FEAT_READ}
    nn = {b: nnf[b].array(library="np")[mask] for b in nnf.keys() if b.startswith("NN_Score")}
    return feat, Hm[mask], np.clip(fac[mask], 0, None), nn

def wq(v, wt, q):
    o = np.argsort(v); c = np.cumsum(wt[o]) / wt.sum(); return np.interp(q, c, v[o])

def R_at_eff(score, wB, peak, frac, eff=0.002):
    thr = wq(score, wB, 1 - eff); p = score > thr
    return (wB[p & peak].sum() / wB[p].sum()) / frac if wB[p].sum() > 0 else float("nan")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bdt", default="model_Za_BDT_run3.pkl")
    ap.add_argument("--out", default="plots_nn/compare_nn_vs_bdt.pdf")
    a = ap.parse_args()
    bdt = pickle.load(open(a.bdt, "rb"))
    print(f"[cfg] bdt={a.bdt}  NN scores from {NN_SCORED}")

    # ---- background: sculpting R(mA) for both models ----
    fb, Hmb, wB, nnb = load_sample("All_Bkg")
    peak = (Hmb > 120) & (Hmb < 130); frac = wB[peak].sum() / wB.sum()
    rows = []
    for mh in EVAL_MASSES:
        s_bdt = bdt.predict_proba(_feat(fb, mh))[:, 1]
        s_nn = nnb[f"NN_Score_mA_M{mh}"]
        rows.append((mh, R_at_eff(s_bdt, wB, peak, frac), R_at_eff(s_nn, wB, peak, frac)))

    # ---- parametrized AUC: signal@true-hypothesis vs bkg@same-hypothesis, pooled over masses ----
    sb_bdt, sn_bdt, sb_nn, sn_nn, yb, wb = [], [], [], [], [], []
    for m in MASSES:
        fs, _, ws, nns = load_sample(f"mA_M{m}")
        sb_bdt += [bdt.predict_proba(_feat(fs, m))[:, 1], bdt.predict_proba(_feat(fb, m))[:, 1]]
        sb_nn  += [nns[f"NN_Score_mA_M{m}"],               nnb[f"NN_Score_mA_M{m}"]]
        yb += [np.ones(len(ws)), np.zeros(len(wB))]; wb += [ws, wB]
    y = np.concatenate(yb); w = np.concatenate(wb)
    auc_bdt = roc_auc_score(y, np.concatenate(sb_bdt), sample_weight=w)
    auc_nn  = roc_auc_score(y, np.concatenate(sb_nn),  sample_weight=w)

    print(f"\n{'mA':>4} {'R_BDT':>8} {'R_NN':>8} {'NN-BDT':>8}")
    for mh, rb, rn in rows:
        print(f"{mh:>4} {rb:>8.2f} {rn:>8.2f} {rn-rb:>+8.2f}")
    print(f"\nAUC   BDT={auc_bdt:.4f}   NN={auc_nn:.4f}   (NN-BDT={auc_nn-auc_bdt:+.4f})")

    mA = [r[0] for r in rows]
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.axhline(1.0, color="grey", ls=":", lw=1, label="no sculpting (R=1)")
    ax.plot(mA, [r[1] for r in rows], "o-", color="tab:orange", lw=2, label=f"BDT (AUC {auc_bdt:.3f})")
    ax.plot(mA, [r[2] for r in rows], "s-", color="tab:blue",
            lw=2, label=f"NN per-mass $\\lambda$100 (AUC {auc_nn:.3f})")
    ax.set_xlabel("$m_A$ hypothesis [GeV]", fontsize=15)
    ax.set_ylabel("Sculpting R (125-peak, eff=0.002)", fontsize=15)
    ax.tick_params(labelsize=12, direction="in", top=True, right=True); ax.legend(fontsize=12, frameon=False)
    fig.tight_layout(); fig.savefig(a.out, dpi=200); print("[plot]", a.out)

if __name__ == "__main__":
    main()
