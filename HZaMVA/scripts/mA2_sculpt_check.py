#!/usr/bin/env python3
"""Diagnose the mA=2 limit anomaly: is the background m_llgammagamma SCULPTED to peak at 125
at the actual working point (low-mass BDT score > 0.980)? If the post-cut background develops a
125 peak like the signal, signal and background become degenerate -> weak (~1) limit.
Plots normalized m_llgg for: bkg pre-cut, bkg post-cut(0.980), signal mA=2. CPU (higgs-alp-ana).
"""
import json, pickle
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from hza_features import _load, _feat, BASE_VARS, ROOT_DIR

FULL16 = BASE_VARS + ["pho_pt_asym", "param"]
MA = 2; CUT = 0.980; LUMI = 172.0
clf = pickle.load(open("model_Za_BDT_lowmass_run3.pkl", "rb"))
feats = json.load(open("model_Za_BDT_lowmass_run3.meta.json"))["features"]
idx = [FULL16.index(f) for f in feats]

bk = _load(f"{ROOT_DIR}/All_Bkg/run3.root", "inclusive")
sg = _load(f"{ROOT_DIR}/mA_M{MA}/run3.root", "inclusive")
Hmb = bk["H_m"]; wb = np.clip(bk["factor"], 0, None)
Hms = sg["H_m"]; wsig = np.clip(sg["factor"], 0, None)
scb = clf.predict_proba(_feat(bk, MA)[:, idx])[:, 1]
pas = scb > CUT
# R at this working point: post-cut 120-130 fraction / pre-cut
peak = lambda H: (H > 120) & (H < 130)
frac_pre = wb[peak(Hmb)].sum() / wb.sum()
frac_post = wb[pas & peak(Hmb)].sum() / wb[pas].sum()
R_wp = frac_post / frac_pre
bkg_eff = wb[pas].sum() / wb.sum()
print(f"[mA=2 cut={CUT}] bkg eff={bkg_eff:.4f}  pre-cut 120-130 frac={frac_pre:.3f}  "
      f"post-cut frac={frac_post:.3f}  R(working point)={R_wp:.2f}")

rng = (100, 180); bins = 80
fig, ax = plt.subplots(figsize=(8, 6))
ax.hist(Hmb, bins=bins, range=rng, weights=wb, density=True, histtype="step", color="0.6", lw=2, label="Bkg (pre-cut)")
ax.hist(Hmb[pas], bins=bins, range=rng, weights=wb[pas], density=True, histtype="stepfilled",
        color="tab:red", alpha=0.45, label=f"Bkg (score>{CUT})")
ax.hist(Hms, bins=bins, range=rng, weights=wsig, density=True, histtype="step", color="tab:blue", lw=2, label="Signal $m_A$=2")
ax.axvspan(120, 130, color="orange", alpha=0.12)
ax.axvline(125, color="k", ls=":", lw=1)
ax.set_xlabel(r"$m_{\ell\ell\gamma\gamma}$ [GeV]", fontsize=14)
ax.set_ylabel("a.u. (normalized)", fontsize=14)
ax.legend(fontsize=12, frameon=False, loc="upper right")
ax.tick_params(direction="in", top=True, right=True)
ax.text(0.0, 1.02, "CMS", transform=ax.transAxes, fontsize=15, fontweight="bold")
ax.text(0.085, 1.02, "Preliminary", transform=ax.transAxes, fontsize=13, style="italic")
ax.text(1.0, 1.02, f"{LUMI:.0f} fb$^{{-1}}$ (13.6 TeV)", transform=ax.transAxes, fontsize=13, ha="right")
ax.text(0.03, 0.80, f"low-mass BDT, score>{CUT}\nR(working point)={R_wp:.2f}", transform=ax.transAxes, fontsize=11)
out = "../plots_MVA/run3_lowmass/mA2_sculpt_check.pdf"
fig.tight_layout(); fig.savefig(out, dpi=200); print("[plot]", out)
