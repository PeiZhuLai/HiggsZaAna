#!/usr/bin/env python3
"""
V3 / V8 reviewer follow-up: yields of ttbar (TTto2L2Nu, 10% file subset), ttbar+gamma (TTG),
ttbar+gammagamma (TTGG) and Z+gammagamma (ZGG) after the full HZa selection, at the adopted BDT
working point of every mass hypothesis m_a = 1..30 GeV, in 115 < m_llgg < 135 GeV, compared with
DY+gamma (DYGto2LG_10to100) and DY+jets from run3_bdt_scored_fsrfix.

Conventions copied from Plot/scripts/make_sculpt_R_table_from_wp.py (the script that consumes the
same scored files and the same WP json):
  tree 'inclusive', mass 'H_mass', weight 'factor' (== weight_central, xs*lumi/sum(genWeight)
  applied at HiggsDNA merge), score MVA_Score_mA_M{m}, pass = score > cut.
Working points: Plot/output/MVAcut_points_run3.json (14 anchors); an interpolated mass takes the
cut of the nearest anchor (ties -> lower anchor; none occur for integer masses).

The new samples exist only for 2024 (109.82 fb^-1). Full Run 3 (172.13 fb^-1) is obtained by
scaling their 2024 yield by 172.13/109.82; DY is quoted both for 2024 and summed over all eras.

TTto2L2Nu uses 147 of 1466 files. Its 'factor' is normalized with the sum of genWeight of the
files actually processed, so the yield already corresponds to the full cross section; no extra
1/fraction factor is applied (applying one would double count).

Usage:  python compute_extrabkg_yields.py [--new-base DIR] [--out-prefix PATH]
"""
import argparse
import json
import os

import numpy as np
import uproot

DY_BASE = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix"
NEW_BASE = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix_extraBkg"
WP_JSON = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/MVAcut_points_run3.json"
LUMI = {"2022preEE": 7.99, "2022postEE": 26.68, "2023preBPix": 17.96, "2023postBPix": 9.68, "2024": 109.82}
LUMI_RUN3 = sum(LUMI.values())
EXTRAP = LUMI_RUN3 / LUMI["2024"]
WINDOW = (115.0, 135.0)
MASSES = list(range(1, 31))

DY_GAMMA = {e: ["DYGto2LG_10to100"] for e in LUMI}
DY_JETS = {e: ["DYJetsToLL"] for e in LUMI if e != "2024"}
DY_JETS["2024"] = ["DYJetsTo2E", "DYJetsTo2Mu", "DYJetsTo2Tau"]
NEW = ["TTto2L2Nu", "TTG", "TTGG_Run3", "ZGG"]


def load_wp():
    anchors = {int(e["mA"]): float(e["MVAcut"]) for e in json.load(open(WP_JSON))["results"]
               if e.get("MVAcut") is not None}
    wp = {}
    for m in MASSES:
        a = min(sorted(anchors), key=lambda x: (abs(x - m), x))
        wp[m] = (a, anchors[a])
    return wp


def load(path):
    t = uproot.open(path)["inclusive"]
    br = ["H_mass", "factor", "n_iso_photons"] + ["MVA_Score_mA_M%d" % m for m in MASSES]
    br = [b for b in br if b in t.keys()]
    return t.arrays(br, library="np")


def cat(arrs):
    keys = set.intersection(*[set(a) for a in arrs]) if arrs else set()
    return {k: np.concatenate([a[k] for a in arrs]) for k in keys}


def yields(arr, wp, extra_mask=None):
    """-> {m: (sumw, sqrt(sumw2), n_raw)} after window + WP; plus the window-only preselection."""
    m, w = arr["H_mass"], arr["factor"]
    base = (m > WINDOW[0]) & (m < WINDOW[1]) & np.isfinite(w)
    if extra_mask is not None:
        base &= extra_mask
    out = {"presel_window": (float(w[base].sum()), float(np.sqrt((w[base] ** 2).sum())), int(base.sum()))}
    for mm in MASSES:
        s = arr["MVA_Score_mA_M%d" % mm]
        sel = base & np.isfinite(s) & (s > wp[mm][1])
        out[mm] = (float(w[sel].sum()), float(np.sqrt((w[sel] ** 2).sum())), int(sel.sum()))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--new-base", default=NEW_BASE)
    ap.add_argument("--out-prefix", default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                           "yields_extrabkg"))
    a = ap.parse_args()
    wp = load_wp()

    res = {"window": WINDOW, "lumi_2024": LUMI["2024"], "lumi_run3": LUMI_RUN3, "extrap_2024_to_run3": EXTRAP,
           "wp": {m: {"anchor": wp[m][0], "cut": wp[m][1]} for m in MASSES}, "yields": {}}

    # DY, 2024 and full Run 3 (all eras)
    for label, table in (("DYgamma", DY_GAMMA), ("DYjets", DY_JETS)):
        a24 = cat([load(os.path.join(DY_BASE, s, "2024.root")) for s in table["2024"]])
        aall = cat([load(os.path.join(DY_BASE, s, e + ".root")) for e in table for s in table[e]])
        res["yields"][label + "_2024"] = yields(a24, wp)
        res["yields"][label + "_run3"] = yields(aall, wp)
        if "n_iso_photons" in a24:
            # informational: DY+jets with a generator-isolated photon (overlap with DY+gamma)
            res["yields"][label + "_2024_nIsoPho0"] = yields(a24, wp, a24["n_iso_photons"] == 0)

    missing = []
    for s in NEW:
        p = os.path.join(a.new_base, s, "2024.root")
        if not os.path.exists(p):
            # e.g. every job selected zero events -> no merged parquet -> no ROOT file.
            # Counted as 0 in the sums and listed at the top of the table; check the
            # reconcile log before believing a 0.
            print("MISSING", p)
            missing.append(s)
            continue
        arr = load(p)
        res["yields"][s + "_2024"] = yields(arr, wp)
        if s == "TTto2L2Nu":   # overlap removal w.r.t. TTG, same convention as MCOverlapTagger
            res["yields"][s + "_2024_nIsoPho0"] = yields(arr, wp, arr["n_iso_photons"] == 0)
        if s == "TTG":
            res["yields"][s + "_2024_nIsoPhoGe1"] = yields(arr, wp, arr["n_iso_photons"] > 0)

    json.dump(res, open(a.out_prefix + ".json", "w"), indent=1)

    Y = res["yields"]
    res["missing_new_samples"] = missing
    json.dump(res, open(a.out_prefix + ".json", "w"), indent=1)
    def g(k, m):
        return Y[k][m] if k in Y else (0.0, 0.0, 0)
    cols = [("TTto2L2Nu_2024", "tt(2l2nu)"), ("TTG_2024", "ttG"), ("TTGG_Run3_2024", "ttGG"),
            ("ZGG_2024", "ZGG"), ("DYgamma_2024", "DY+G"), ("DYjets_2024", "DY+jets")]
    lines = ["# Extra-background yields, 2024 (%.2f fb^-1), %g<m_llgg<%g GeV, weight=factor" % (LUMI["2024"], *WINDOW), "",
             "Missing new-sample ROOT files (counted as 0): %s" % (", ".join(missing) or "none"), "",
             "Cells: yield ± stat [raw MC events]", "",
             "| m_a | WP (anchor) | " + " | ".join(c[1] for c in cols) + " | (tt+ttG+ttGG)/DY | ZGG/DY |",
             "|---|---|" + "---|" * len(cols) + "---|---|"]
    for m in ["presel_window"] + MASSES:
        row = [str(m), "-" if m == "presel_window" else "%.3f (%d)" % (wp[m][1], wp[m][0])]
        for k, _ in cols:
            y, e, n = g(k, m)
            row.append("%.3g ± %.2g [%d]" % (y, e, n))
        dy = g("DYgamma_2024", m)[0] + g("DYjets_2024", m)[0]
        tt = sum(g(k, m)[0] for k in ("TTto2L2Nu_2024", "TTG_2024", "TTGG_Run3_2024"))
        z = g("ZGG_2024", m)[0]
        row.append("%.3g" % (tt / dy) if dy > 0 else "n/a")
        row.append("%.3g" % (z / dy) if dy > 0 else "n/a")
        lines.append("| " + " | ".join(row) + " |")
    lines += ["", "Full Run 3 (%.2f fb^-1): new samples x %.4f; DY summed over eras." % (LUMI_RUN3, EXTRAP), "",
              "| m_a | tt+ttG+ttGG | ZGG | DY+G | DY+jets | (tt+ttG+ttGG)/DY | ZGG/DY |", "|---|---|---|---|---|---|---|"]
    for m in ["presel_window"] + MASSES:
        tt = EXTRAP * sum(g(k, m)[0] for k in ("TTto2L2Nu_2024", "TTG_2024", "TTGG_Run3_2024"))
        z = EXTRAP * g("ZGG_2024", m)[0]
        dg, dj = g("DYgamma_run3", m)[0], g("DYjets_run3", m)[0]
        dy = dg + dj
        lines.append("| %s | %.3g | %.3g | %.3g | %.3g | %s | %s |" % (
            m, tt, z, dg, dj, ("%.3g" % (tt / dy)) if dy > 0 else "n/a", ("%.3g" % (z / dy)) if dy > 0 else "n/a"))
    open(a.out_prefix + ".md", "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
