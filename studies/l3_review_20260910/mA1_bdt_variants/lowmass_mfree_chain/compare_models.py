#!/usr/bin/env python3
"""Production low-mass BDT vs win_mfree at m_a = 1, 2, 3 GeV, before investing in the full chain.

Both models are scored on the same events: production scores from run3_bdt_scored_fsrfix, new
scores from the sandbox copy (rescore_lowmass.py; identical files except MVA_Score_mA_M{1,2,3}).

Per model, mass point and cut:
  sculpting  DY MC (sculpt-R recipe, raw factor, tree inclusive, 95-180):
               R      = f(120-130)^pass / f(120-130)^incl          (+- 200 bootstrap replicas)
               R_lo   = same for 105-115, R_hi = same for 135-145  (turn-on / shoulder probes)
               chi2   = pass vs inclusive normalized shape, 5 GeV bins, MC stat errors
             data, sidebands only (95-115, 135-180; the 115-135 window is never read):
               shape of the passing events / shape of all preselected data, per band, and chi2
  signal     mA_M<m> test tree (30%), weight factor/0.3, no signal reweight (proxy only):
               yield in 95-180, in 115-133, window fraction, overlap share (gen dR(gg) < 0.05),
               second-population fraction f(135-160)/f(115-133)
  proxy Z    S(115-133) / sqrt(B(115-133)), B = data sidebands x DY window/sideband ratio
Cuts: the production working point (0.963 / 0.975 / 0.988) for the production model, and for the
new model the cut with the same data sideband count (eqdata) plus the R = 1 crossing.
Env: higgs-alp-ana. Output: compare_models.{txt,json}, compare_scan_<model>_mA<m>.csv
"""
import json
import os

import numpy as np
import uproot

HERE = os.path.dirname(os.path.abspath(__file__))
PROD = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix"
NEW = "/eos/home-p/pelai/HZa/sandbox_lowmass_mfree/run3_bdt_scored"
ERAS = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]
BKG = {"2022preEE": ["DYGto2LG_10to100", "DYJetsToLL"], "2022postEE": ["DYGto2LG_10to100", "DYJetsToLL"],
       "2023preBPix": ["DYGto2LG_10to100", "DYJetsToLL"], "2023postBPix": ["DYGto2LG_10to100", "DYJetsToLL"],
       "2024": ["DYGto2LG_10to100", "DYJetsTo2E", "DYJetsTo2Mu", "DYJetsTo2Tau"]}
PROD_WP = {1: 0.963, 2: 0.975, 3: 0.988}   # = Plot/output/MVAcut_points_run3.json
MASSES = (1, 2, 3)
WIN = (115.0, 133.0)
SB_BANDS = [(95, 105), (105, 115), (135, 145), (145, 160), (160, 180)]
SHAPE_EDGES = np.arange(95, 181, 5.0)
TEST_FRAC = 0.3
GRID = np.unique(np.concatenate([np.linspace(0.70, 0.95, 26), np.linspace(0.95, 0.999, 99)]))


def load(sample_eras, tree, extra=()):
    br = ["H_mass", "factor"] + [f"MVA_Score_mA_M{m}" for m in MASSES] + list(extra)
    out = {"prod": {}, "new": {}}
    for s, e in sample_eras:
        for tag, base in (("prod", PROD), ("new", NEW)):
            a = uproot.open(f"{base}/{s}/{e}.root")[tree].arrays(br, library="np")
            for k, v in a.items():
                out[tag].setdefault(k, []).append(v)
    res = {}
    for tag in out:
        res[tag] = {k: (np.concatenate(v) if not isinstance(v[0], np.ndarray) or v[0].dtype != object
                        else np.concatenate(v)) for k, v in out[tag].items()}
    # the two copies must hold the same events in the same order
    assert np.array_equal(res["prod"]["H_mass"], res["new"]["H_mass"], equal_nan=True)
    return res


def inr(m, lo, hi):
    return (m > lo) & (m < hi)


def is_sb(m):
    return inr(m, 95, 115) | inr(m, 135, 180)


def frac(w, m, lo, hi, sel):
    t = w[sel].sum()
    return w[sel & inr(m, lo, hi)].sum() / t if t > 0 else np.nan


def ratio_band(w, m, s, cut, lo, hi):
    base = np.ones_like(m, bool)
    p = s > cut
    f0 = frac(w, m, lo, hi, base); f1 = frac(w, m, lo, hi, p)
    return f1 / f0 if f0 > 0 else np.nan


def shape_chi2(w, m, s, cut, edges, sel_extra=None):
    """chi2/ndf of passing vs inclusive normalized shapes (MC-stat errors of the passing set)."""
    sel = np.ones_like(m, bool) if sel_extra is None else sel_extra
    p = sel & (s > cut)
    h0, _ = np.histogram(m[sel], edges, weights=w[sel]); h1, _ = np.histogram(m[p], edges, weights=w[p])
    e1, _ = np.histogram(m[p], edges, weights=w[p] ** 2)
    if h1.sum() <= 0:
        return np.nan, 0
    n0, n1 = h0 / h0.sum(), h1 / h1.sum()
    var = e1 / h1.sum() ** 2
    ok = (var > 0) & (h0 > 0)
    c2 = np.sum((n1[ok] - n0[ok]) ** 2 / var[ok])
    return float(c2), int(ok.sum() - 1)


def main():
    dy = load([(s, e) for e in ERAS for s in BKG[e]], "inclusive")
    da = load([("Data", e) for e in ERAS], "inclusive")
    rng = np.random.default_rng(13)
    summary = {}
    lines = []
    for mA in MASSES:
        sig = load([(f"mA_M{mA}", e) for e in ERAS], "test", extra=["GenALP_dR_gg"])
        for tag in ("prod", "new"):
            sc = f"MVA_Score_mA_M{mA}"
            b, d, s = dy[tag], da[tag], sig[tag]
            bm, bw, bs = b["H_mass"], b["factor"].astype(float), b[sc]
            okb = inr(bm, 95, 180) & np.isfinite(bw) & np.isfinite(bs)
            bm, bw, bs = bm[okb], bw[okb], bs[okb]
            dm, ds = d["H_mass"], d[sc]
            okd = inr(dm, 95, 180) & is_sb(dm) & np.isfinite(ds)     # sidebands only
            dm, ds = dm[okd], ds[okd]
            dw = np.ones_like(dm)
            sm, sw, ss = s["H_mass"], s["factor"].astype(float) / TEST_FRAC, s[sc]
            gdr = np.array([x[0] if len(x) else np.nan for x in s["GenALP_dR_gg"]], float)
            oks = inr(sm, 95, 180) & np.isfinite(ss)
            sm, sw, ss, gdr = sm[oks], sw[oks], ss[oks], gdr[oks]
            pk_incl = frac(bw, bm, 120, 130, np.ones_like(bm, bool))
            rows = []
            boots = rng.poisson(1.0, size=(200, len(bm))).astype(float)
            for c in GRID:
                p = bs > c
                if bw[p].sum() <= 0:
                    continue
                n_dy120 = int((p & inr(bm, 120, 130)).sum())
                R = frac(bw, bm, 120, 130, p) / pk_incl
                nd = int((ds > c).sum())
                ps = ss > c
                S_all = sw[ps].sum(); S_win = sw[ps & inr(sm, *WIN)].sum()
                sec = sw[ps & inr(sm, 135, 160)].sum()
                ov = sw[ps & (gdr < 0.05)].sum() / S_all if S_all > 0 else np.nan
                ov_win = (sw[ps & inr(sm, *WIN) & (gdr < 0.05)].sum() / S_win) if S_win > 0 else np.nan
                dy_win = bw[p & inr(bm, *WIN)].sum(); dy_sb = bw[p & is_sb(bm)].sum()
                B_win = nd * dy_win / dy_sb if dy_sb > 0 else np.nan
                rows.append(dict(cut=float(c), R=float(R), n_dy120=n_dy120,
                                 R_lo=float(ratio_band(bw, bm, bs, c, 105, 115)),
                                 R_hi=float(ratio_band(bw, bm, bs, c, 135, 145)),
                                 data_sb=nd, S_all=float(S_all), S_win=float(S_win),
                                 win_frac=float(S_win / S_all) if S_all > 0 else np.nan,
                                 sec_frac=float(sec / S_win) if S_win > 0 else np.nan,
                                 overlap=float(ov), overlap_win=float(ov_win), B_win=float(B_win),
                                 Z=float(S_win / np.sqrt(B_win)) if B_win > 0 else np.nan))
            # preselection reference (no cut)
            S0 = sw.sum(); S0w = sw[inr(sm, *WIN)].sum()
            pre = dict(win_frac=float(S0w / S0), overlap=float(sw[gdr < 0.05].sum() / S0),
                       sec_frac=float(sw[inr(sm, 135, 160)].sum() / S0w))
            import csv
            with open(os.path.join(HERE, f"compare_scan_{tag}_mA{mA}.csv"), "w") as fo:
                wr = csv.DictWriter(fo, fieldnames=list(rows[0])); wr.writeheader(); wr.writerows(rows)
            # working points
            nd_prod = int((da["prod"][f"MVA_Score_mA_M{mA}"][ (inr(da["prod"]["H_mass"], 95, 180)) & is_sb(da["prod"]["H_mass"])] > PROD_WP[mA]).sum())
            wps = {}
            if tag == "prod":
                wps["prod_wp"] = PROD_WP[mA]
            else:
                # eqdata: smallest cut whose data sideband count <= production's
                cands = [r for r in rows if r["data_sb"] <= nd_prod]
                wps["eqdata"] = cands[0]["cut"] if cands else np.nan
                # R=1 crossing: tightest grid cut with R >= 1 above which R < 1, with >= 20 DY events in 120-130
                good = [r for r in rows if r["n_dy120"] >= 20]
                cr = [good[i]["cut"] for i in range(len(good) - 1) if (good[i]["R"] - 1) * (good[i + 1]["R"] - 1) <= 0]
                wps["R1"] = cr[-1] if cr else np.nan
            for name, c in wps.items():
                if not np.isfinite(c):
                    continue
                r = min(rows, key=lambda x: abs(x["cut"] - c))
                p = bs > c
                Rb = []
                for bb in boots[:200]:
                    ww = bw * bb
                    Rb.append(frac(ww, bm, 120, 130, p) / frac(ww, bm, 120, 130, np.ones_like(bm, bool)))
                c2, ndf = shape_chi2(bw, bm, bs, c, SHAPE_EDGES)
                # data sideband shape: passing / inclusive per band
                dsel = ds > c
                band = {f"{lo}-{hi}": float((inr(dm[dsel], lo, hi).sum() / max(dsel.sum(), 1)) /
                                             (inr(dm, lo, hi).sum() / len(dm))) for lo, hi in SB_BANDS}
                ed = np.array([95, 105, 115, 135, 145, 160, 180], float)
                sbmask = np.ones_like(dm, bool)
                c2d, ndfd = shape_chi2(dw, dm, ds, c, ed, sbmask)
                rec = dict(r, R_err=float(np.nanstd(Rb)), dy_shape_chi2=c2, dy_shape_ndf=ndf,
                           data_sb_band_ratio=band, data_sb_chi2=c2d, data_sb_ndf=ndfd, wp=name, pre=pre)
                summary[f"{tag}_mA{mA}_{name}"] = rec
                lines.append(f"{tag:4s} mA{mA} {name:8s} cut={r['cut']:.4f} R={r['R']:.2f}+-{rec['R_err']:.2f} "
                             f"R_lo={r['R_lo']:.2f} R_hi={r['R_hi']:.2f} DYchi2={c2:.1f}/{ndf} "
                             f"dataSB={r['data_sb']} dataSBchi2={c2d:.1f}/{ndfd} "
                             f"band[{' '.join('%.2f' % v for v in band.values())}] | "
                             f"S_win={r['S_win']:.1f} winfrac={r['win_frac']:.2f} sec={r['sec_frac']:.2f} "
                             f"ovl={r['overlap']:.2f} ovl_win={r['overlap_win']:.2f} Z={r['Z']:.2f} "
                             f"| presel winfrac={pre['win_frac']:.2f} ovl={pre['overlap']:.2f} sec={pre['sec_frac']:.2f}")
                print(lines[-1], flush=True)
    open(os.path.join(HERE, "compare_models.txt"), "w").write("\n".join(lines) + "\n")
    json.dump(summary, open(os.path.join(HERE, "compare_models.json"), "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
