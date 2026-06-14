#!/usr/bin/env python3
"""Scan the mA=2 BDT working point vs expected limit (anchored approximation).
For each cut c: S(c), B(c) = signal/background weighted yield in the 120-130 SR after score>c
(low-mass model, mA=2 hypothesis). Expected 95% CL limit ~ 1.64*sqrt(B)/S (single-bin Asimov);
ANCHORED to the real combine expected limit at the production cut (0.980 -> 0.87) so the absolute
scale (lumi*xs*br, ele+mu cats, systematics) is folded in via the anchor. Identifies a better cut.
"""
import json, pickle
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from hza_features import _load, _feat, BASE_VARS, ROOT_DIR

FULL16 = BASE_VARS + ["pho_pt_asym", "param"]
MA = 2; SR = (120, 130); ANCHOR_CUT = 0.980; ANCHOR_LIM = 0.871; LUMI = 172.0
clf = pickle.load(open("model_Za_BDT_lowmass_run3.pkl", "rb"))
feats = json.load(open("model_Za_BDT_lowmass_run3.meta.json"))["features"]
idx = [FULL16.index(f) for f in feats]

bk = _load(f"{ROOT_DIR}/All_Bkg/run3.root", "inclusive")
sg = _load(f"{ROOT_DIR}/mA_M{MA}/run3.root", "inclusive")
scb = clf.predict_proba(_feat(bk, MA)[:, idx])[:, 1]; wb = np.clip(bk["factor"], 0, None); Hmb = bk["H_m"]
scs = clf.predict_proba(_feat(sg, MA)[:, idx])[:, 1]; ws = np.clip(sg["factor"], 0, None); Hms = sg["H_m"]
srb = (Hmb > SR[0]) & (Hmb < SR[1]); srs = (Hms > SR[0]) & (Hms < SR[1])

cuts = [0.90, 0.93, 0.95, 0.96, 0.97, 0.975, 0.980, 0.985, 0.990, 0.993, 0.995, 0.997, 0.998, 0.999]
def fom(c):
    S = ws[(scs > c) & srs].sum(); B = wb[(scb > c) & srb].sum()
    nB = int(((scb > c) & srb).sum())   # raw bkg events (reliability)
    return S, B, nB
S0, B0, _ = fom(ANCHOR_CUT); f0 = np.sqrt(B0) / S0
print(f"anchor: cut={ANCHOR_CUT} S={S0:.4g} B={B0:.4g} sqrtB/S={f0:.4g} -> combine expected={ANCHOR_LIM}")
print(f"\n{'cut':>6} {'S(rel)':>9} {'B(rel)':>9} {'nBkg':>6} {'sqrtB/S':>9} {'~exp.lim':>9}")
rows = []
for c in cuts:
    S, B, nB = fom(c)
    lim = ANCHOR_LIM * (np.sqrt(B) / S) / f0 if S > 0 and B > 0 else np.nan
    rows.append((c, S, B, nB, lim))
    print(f"{c:>6.3f} {S:>9.4g} {B:>9.4g} {nB:>6d} {np.sqrt(B)/S:>9.4g} {lim:>9.3f}")
best = min((r for r in rows if r[3] >= 20 and r[4] == r[4]), key=lambda r: r[4])
print(f"\n>>> best (with >=20 bkg events for reliability): cut={best[0]:.3f}  ~expected limit={best[4]:.3f}  "
      f"(vs {ANCHOR_LIM} at 0.980)")

fig, ax = plt.subplots(figsize=(8, 6))
cc = [r[0] for r in rows]; ll = [r[4] for r in rows]
ax.plot(cc, ll, "o-", color="tab:blue", lw=2)
ax.axhline(ANCHOR_LIM, color="grey", ls=":", label=f"current (cut 0.980) = {ANCHOR_LIM}")
ax.axvline(best[0], color="tab:red", ls="--", label=f"best cut {best[0]:.3f} -> {best[4]:.3f}")
ax.set_xlabel("low-mass BDT cut", fontsize=14); ax.set_ylabel("approx. expected limit on r (mA=2)", fontsize=14)
ax.legend(fontsize=11, frameon=False); ax.tick_params(direction="in", top=True, right=True)
ax.text(0.0, 1.02, "CMS", transform=ax.transAxes, fontsize=15, fontweight="bold")
ax.text(0.085, 1.02, "Preliminary", transform=ax.transAxes, fontsize=13, style="italic")
ax.text(1.0, 1.02, f"{LUMI:.0f} fb$^{{-1}}$ (13.6 TeV)", transform=ax.transAxes, fontsize=13, ha="right")
out = "../plots_MVA/run3_lowmass/mA2_cut_vs_limit.pdf"
fig.tight_layout(); fig.savefig(out, dpi=200); print("[plot]", out)
