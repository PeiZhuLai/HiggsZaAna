#!/usr/bin/env python3
"""Evaluate the sandbox low-mass BDT variants (train_variants.py) at m_a = 1 GeV.

Samples (run3_bdt_scored_fsrfix, 5 Run 3 eras):
  signal      mA_M1/<era>.root tree "test" (the 30% the signal model uses), weight factor / 0.3
  background  combined DY recipe of make_sculpt_R_table_from_wp.py, tree "inclusive", raw factor
  data        Data/<era>.root, only the sidebands 105-115 and 133-180 GeV are used (the signal
              window 115-133 is never looked at)

Working points per variant:
  R1      production rule for m_a = 1: tightest R = 1 crossing on the scan grid
          (scan_score_R_significance.py: 45 cuts in 0.70-0.99, >= 20 bkg events in 120-130,
          R with raw weights)
  eqdata  cut giving the same number of data events in the sidebands as production
          (stored score > 0.963), i.e. the same background statistics for the fit

Reported per working point: R, signal yield in 115-133 and in 95-180, overlap share
(gen dR(gg) < 0.05) after the cut and in the window, data sideband count, background
expected in 115-133 (data sidebands x DY window/sideband ratio, and a linear interpolation
of the data density) and an Asimov Z proxy. The final judgment is the expected limit.
Env: higgs-alp-ana. Usage: eval_variants.py [variant ...]
"""
import json
import os
import pickle
import sys

import numpy as np
import uproot

HERE = os.path.dirname(os.path.abspath(__file__))
SB = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix"
ERAS = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]
BKG = {"2022preEE": ["DYGto2LG_10to100", "DYJetsToLL"], "2022postEE": ["DYGto2LG_10to100", "DYJetsToLL"],
       "2023preBPix": ["DYGto2LG_10to100", "DYJetsToLL"], "2023postBPix": ["DYGto2LG_10to100", "DYJetsToLL"],
       "2024": ["DYGto2LG_10to100", "DYJetsTo2E", "DYJetsTo2Mu", "DYJetsTo2Tau"]}
RAW = ["ALP_calculatedPhotonIso", "var_dR_g1g2", "pho1R9", "pho1Pt_oHm", "ALP_mass", "H_mass", "factor", "MVA_Score_mA_M1"]
TEST_FRAC = 0.3
WIN = (115.0, 133.0)
SIDE = ((105.0, 115.0), (133.0, 180.0))
PROD_WP = 0.963
MHYP = 1.0


def load(paths, tree, extra=()):
    out = {}
    for p in paths:
        a = uproot.open(p)[tree].arrays(RAW + list(extra), library="np")
        for k, v in a.items():
            out.setdefault(k, []).append(v)
    a = {k: np.concatenate(v).astype(np.float64) for k, v in out.items()}
    sel = (a["H_mass"] > 95) & (a["H_mass"] < 180) & np.isfinite(a["factor"])
    return {k: v[sel] for k, v in a.items()}


def score(clf, feats, a, mhyp=MHYP):
    der = {"param": lambda: (a["ALP_mass"] - mhyp) / a["H_mass"],
           "param_rel": lambda: (a["ALP_mass"] - mhyp) / mhyp,
           "pho1Pt": lambda: a["pho1Pt_oHm"] * a["H_mass"]}
    cols = [der[f]() if f in der else a[f] for f in feats]
    return clf.predict_proba(np.column_stack(cols))[:, 1]


def in_range(m, lo, hi):
    return (m > lo) & (m < hi)


def side(m):
    return in_range(m, *SIDE[0]) | in_range(m, *SIDE[1])


def R_at(bk, s, cut):
    w, m = bk["factor"], bk["H_mass"]
    pk = in_range(m, 120, 130); frac = w[pk].sum() / w.sum()
    p = s > cut
    return float(w[p & pk].sum() / w[p].sum() / frac) if w[p].sum() > 0 else float("nan")


def r1_cut(bk, s):
    w, m = bk["factor"], bk["H_mass"]
    pk = in_range(m, 120, 130); frac = w[pk].sum() / w.sum()
    pts = []
    for thr in np.linspace(0.70, 0.99, 45):
        p = s > thr
        if w[p & pk].sum() <= 0 or (p & pk).sum() < 20:
            continue
        pts.append((thr, w[p & pk].sum() / w[p].sum() / frac))
    pts = np.array(pts); best = None
    for i in range(len(pts) - 1):
        r0, r1 = pts[i, 1], pts[i + 1, 1]
        if (r0 - 1) * (r1 - 1) <= 0 and r0 != r1:
            best = pts[i, 0] + (1 - r0) / (r1 - r0) * (pts[i + 1, 0] - pts[i, 0])
    if best is None:
        best = pts[int(np.argmin(np.abs(pts[:, 1] - 1))), 0]
    return float(best), pts


def eqdata_cut(m_data, s_data, n_target):
    s = np.sort(s_data[side(m_data)])[::-1]
    return float(0.5 * (s[n_target - 1] + s[n_target])) if len(s) > n_target else float(s[-1])


def zasimov(s, b):
    return float(np.sqrt(2 * ((s + b) * np.log(1 + s / b) - s))) if b > 0 and s > 0 else float("nan")


def metrics(cut, sg, ss, bk, sb, da, sd):
    w, m, ov = sg["factor"] / TEST_FRAC, sg["H_mass"], sg["GenALP_dR_gg"] < 0.05
    p = ss > cut; win = in_range(m, *WIN)
    wb, mb = bk["factor"], bk["H_mass"]; pb = sb > cut
    md = da["H_mass"]; pd_ = sd > cut
    n_side = int((pd_ & side(md)).sum())
    mc_ratio = wb[pb & in_range(mb, *WIN)].sum() / wb[pb & side(mb)].sum()
    b_mc = n_side * mc_ratio
    lo = (pd_ & in_range(md, 105, 115)).sum() / 10.0; hi = (pd_ & in_range(md, 133, 143)).sum() / 10.0
    b_lin = 0.5 * (lo + hi) * (WIN[1] - WIN[0])
    s_win = w[p & win].sum()
    return {
        "cut": cut, "R": R_at(bk, sb, cut),
        "S_window": s_win, "S_95_180": w[p].sum(), "eff_window": s_win / w[win].sum(),
        "overlap_share_all": w[p & ov].sum() / w[p].sum(),
        "overlap_share_window": w[p & ov & win].sum() / w[p & win].sum(),
        "window_fraction": s_win / w[p].sum(),
        "data_side": n_side, "data_105_115": int((pd_ & in_range(md, 105, 115)).sum()),
        "B_window_mcratio": b_mc, "B_window_lininterp": b_lin,
        "Z_mcratio": zasimov(s_win, b_mc), "Z_lininterp": zasimov(s_win, b_lin),
    }


def main():
    names = sys.argv[1:] or ["base", "win", "nopt", "win_nopt"]
    sg = load([f"{SB}/mA_M1/{e}.root" for e in ERAS], "test", extra=("GenALP_dR_gg",))
    bk = load([f"{SB}/{s}/{e}.root" for e in ERAS for s in BKG[e]], "inclusive")
    da = load([f"{SB}/Data/{e}.root" for e in ERAS], "inclusive")
    print(f"loaded signal(test) {len(sg['H_mass'])}, DY {len(bk['H_mass'])}, data {len(da['H_mass'])}")

    # production reference: stored scores
    ref = {"production": metrics(PROD_WP, sg, sg["MVA_Score_mA_M1"], bk, bk["MVA_Score_mA_M1"], da, da["MVA_Score_mA_M1"])}
    n_target = ref["production"]["data_side"]
    results = {"production_stored_score": ref}
    for name in names:
        meta = json.load(open(os.path.join(HERE, f"model_{name}.meta.json")))
        clf = pickle.load(open(os.path.join(HERE, f"model_{name}.pkl"), "rb"))
        feats = meta["features"]
        ss, sb, sd = score(clf, feats, sg), score(clf, feats, bk), score(clf, feats, da)
        res = {"features": feats}
        if name == "base":
            res["max_abs_diff_vs_stored"] = {k: float(np.abs(x - a["MVA_Score_mA_M1"]).max())
                                             for k, x, a in (("signal", ss, sg), ("DY", sb, bk), ("data", sd, da))}
            print("base vs deployed (stored score) max |diff|:", res["max_abs_diff_vs_stored"])
        c_r1, pts = r1_cut(bk, sb)
        res["R1"] = metrics(c_r1, sg, ss, bk, sb, da, sd)
        res["eqdata"] = metrics(eqdata_cut(da["H_mass"], sd, n_target), sg, ss, bk, sb, da, sd)
        res["R_scan"] = pts.tolist()
        results[name] = res
        np.savez_compressed(os.path.join(HERE, f"scores_{name}.npz"), sig=ss, bkg=sb, data=sd)
    np.savez_compressed(os.path.join(HERE, "samples_meta.npz"),
                        sig_m=sg["H_mass"], sig_w=sg["factor"] / TEST_FRAC, sig_dr=sg["GenALP_dR_gg"],
                        sig_stored=sg["MVA_Score_mA_M1"], bkg_m=bk["H_mass"], bkg_w=bk["factor"],
                        bkg_stored=bk["MVA_Score_mA_M1"], data_m=da["H_mass"], data_stored=da["MVA_Score_mA_M1"])
    json.dump(results, open(os.path.join(HERE, "eval_variants.json"), "w"), indent=1)

    keys = ["cut", "R", "S_window", "S_95_180", "eff_window", "window_fraction", "overlap_share_all",
            "overlap_share_window", "data_side", "B_window_mcratio", "B_window_lininterp", "Z_mcratio", "Z_lininterp"]
    rows = [("production", "stored", ref["production"])]
    rows += [(n, wp, results[n][wp]) for n in names for wp in ("R1", "eqdata")]
    with open(os.path.join(HERE, "eval_variants.txt"), "w") as f:
        hdr = "%-10s %-7s " % ("variant", "WP") + " ".join("%10s" % k[:10] for k in keys)
        print(hdr); f.write(hdr + "\n")
        for n, wp, r in rows:
            line = "%-10s %-7s " % (n, wp) + " ".join("%10.4g" % r[k] for k in keys)
            print(line); f.write(line + "\n")


if __name__ == "__main__":
    main()
