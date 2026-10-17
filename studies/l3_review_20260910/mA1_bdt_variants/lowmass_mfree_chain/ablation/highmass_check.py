#!/usr/bin/env python3
"""Does the production high-mass BDT (mA4-30) enrich a second m_llgg population, as the low-mass one
did at m_a = 1? Per anchor mass point, signal test tree (30%, weight factor/0.3, 5 Run 3 eras,
95-180 GeV), stored production score MVA_Score_mA_M<m> at the working point of
Plot/output/MVAcut_points_run3.json, preselection vs after the cut:
  overlap share     gen dR(gamma,gamma) < 0.05
  sec               S(135-160) / S(115-133)   (second, high-side population)
  window fraction   S(115-133) / S(95-180)
  sigma_eff         smallest interval holding 68.3% of the 115-135 core
A BDT-driven second peak shows up as sec / overlap rising after the cut (mA1 production: overlap
0.24 -> 0.68, sec 0.30 -> 1.25). Env: higgs-alp-ana. Output: highmass_check.txt
"""
import json
import os

import numpy as np
import uproot

HERE = os.path.dirname(os.path.abspath(__file__))
SB = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix"
ERAS = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]
WPJ = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/MVAcut_points_run3.json"


def seff(x, w, frac=0.683):
    o = np.argsort(x); x = x[o]; c = np.cumsum(w[o]); c /= c[-1]; best = np.inf; j = 0
    for i in range(len(x)):
        t = (c[i - 1] if i else 0) + frac
        j = max(j, i)
        while j < len(x) and c[j] < t: j += 1
        if j >= len(x): break
        best = min(best, x[j] - x[i])
    return 0.5 * best


def stats(m, w, g):
    S = w.sum(); Sw = w[(m > 115) & (m < 133)].sum()
    core = (m > 115) & (m < 135)
    return dict(ovl=w[g < 0.05].sum() / S, sec=w[(m > 135) & (m < 160)].sum() / Sw, win=Sw / S,
                seff=seff(m[core], w[core]))


def main():
    wp = {int(e["mA"]): float(e["MVAcut"]) for e in json.load(open(WPJ))["results"]}
    lines = ["mA   cut    | presel: overlap  sec   winfrac seff | after WP: overlap  sec   winfrac seff | eff"]
    for mA in (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30):
        cols = ["H_mass", "factor", f"MVA_Score_mA_M{mA}", "GenALP_dR_gg"]
        acc = {k: [] for k in cols}
        for e in ERAS:
            a = uproot.open(f"{SB}/mA_M{mA}/{e}.root")["test"].arrays(cols, library="np")
            for k in cols: acc[k].append(a[k])
        a = {k: np.concatenate(v) for k, v in acc.items()}
        g = np.array([x[0] if len(x) else np.nan for x in a["GenALP_dR_gg"]], float)
        m = a["H_mass"]; w = a["factor"].astype(float) / 0.3; s = a[f"MVA_Score_mA_M{mA}"]
        ok = (m > 95) & (m < 180) & np.isfinite(w) & np.isfinite(s)
        m, w, s, g = m[ok], w[ok], s[ok], g[ok]
        p = s > wp[mA]
        A, B = stats(m, w, g), stats(m[p], w[p], g[p])
        lines.append(f"{mA:<4d} {wp[mA]:.3f} |        {A['ovl']:.3f}  {A['sec']:.3f}  {A['win']:.3f}  {A['seff']:.2f} |"
                     f"           {B['ovl']:.3f}  {B['sec']:.3f}  {B['win']:.3f}  {B['seff']:.2f} | {w[p].sum() / w.sum():.3f}")
        print(lines[-1], flush=True)
    open(os.path.join(HERE, "highmass_check.txt"), "w").write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
