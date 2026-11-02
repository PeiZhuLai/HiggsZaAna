#!/usr/bin/env python3
"""Diagnose the data/MC excess at high BDT score for m_a = 3, 4, 5 GeV (AN Fig. 104).

Same inputs and weights as the production dataVmc (1_prepare_dataVmc.py): the scored ntuples of
run3_bdt_scored_fsrfix, 95 < H_m < 180, background weight = weight_sideband_rwgt, data weight 1.
Read-only. Output: one JSON with all numbers and a text summary next to this script.

A  data/MC per 0.05 score bin and chi2 over score > 0.5, inclusive (data blinded >= 0.9) and
   sideband (H_m outside 115-135, not blinded); total, per era, per lepton channel.
B  sideband score template fit with free Z+gamma and Z+jets normalizations.
C  shape of the BDT inputs in the sideband, score > 0.5 vs score < 0.5 (data vs MC, normalized).
"""
import json, os, sys
import numpy as np
import uproot
from scipy import stats

B = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix"
ERAS = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]
ZJ = {e: ["DYJetsToLL"] for e in ERAS[:4]}
ZJ["2024"] = ["DYJetsTo2E", "DYJetsTo2Mu", "DYJetsTo2Tau"]
MASSES = [3, 4, 5]
LOW = ["ALP_calculatedPhotonIso", "var_dR_g1g2", "pho1R9", "pho2_relEerr"]
HIGH = ["pho1Pt_oHm", "pho1R9", "pho1IetaIeta55", "pho1PIso_noCorr", "pho2Pt_oHm", "pho2R9",
        "pho2IetaIeta55", "pho2PIso_noCorr", "ALP_calculatedPhotonIso", "var_dR_Za", "var_dR_g1g2",
        "var_dR_g1Z", "var_PtaOverMh", "H_pt_oHm", "pho_pt_asym"]
FEATS = {3: LOW, 4: HIGH, 5: HIGH}
EDGES = np.linspace(0, 1, 21)
HI = 0.5
OUT = os.path.dirname(os.path.abspath(__file__))


def load(path, is_data):
    t = uproot.open(path)["inclusive"]
    vars_ = sorted(set(LOW + HIGH) - {"pho_pt_asym", "pho2_relEerr"}) + [
        "ALP_lead_photon_pt", "ALP_sublead_photon_pt", "ALP_sublead_photon_energyErr", "ALP_sublead_photon_eta"]
    br = ["H_m", "z_mumu", "z_ee"] + ["MVA_Score_mA_M%d" % m for m in MASSES] + vars_
    br += [] if is_data else ["weight_sideband_rwgt"]
    a = t.arrays([b for b in br if b in t.keys()], library="np")
    a["w"] = np.ones(len(a["H_m"])) if is_data else a["weight_sideband_rwgt"].astype(float)
    # derived as in 1_prepare_dataVmc.py (the data ntuples carry no pho2_relEerr branch)
    pt1, pt2 = a["ALP_lead_photon_pt"].astype(float), a["ALP_sublead_photon_pt"].astype(float)
    a["pho_pt_asym"] = (pt1 - pt2) / (pt1 + pt2 + 1e-6)
    a["pho2_relEerr"] = a["ALP_sublead_photon_energyErr"] / (pt2 * np.cosh(a["ALP_sublead_photon_eta"].astype(float)))
    sel = (a["H_m"] > 95) & (a["H_m"] < 180)
    return {k: v[sel] for k, v in a.items()}


def cat(parts):
    return {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}


samples = {}
for e in ERAS:
    samples[("Data", e)] = load("%s/Data/%s.root" % (B, e), True)
    samples[("Zg", e)] = load("%s/DYGto2LG_10to100/%s.root" % (B, e), False)
    samples[("Zj", e)] = cat([load("%s/%s/%s.root" % (B, s, e), False) for s in ZJ[e]])
    print("loaded", e, {k[0]: len(v["w"]) for k, v in samples.items() if k[1] == e}, flush=True)


def mask(a, era=None, ch=None, region="incl"):
    m = np.ones(len(a["w"]), bool)
    if ch == "ee":
        m &= np.abs(a["z_mumu"]) != 1
    if ch == "mumu":
        m &= np.abs(a["z_ee"]) != 1
    if region == "sb":
        m &= (a["H_m"] < 115) | (a["H_m"] > 135)
    return m


def hist(proc, mA, eras, ch, region, var=None, edges=EDGES, extra=None):
    h = np.zeros(len(edges) - 1); s2 = np.zeros(len(edges) - 1)
    for e in eras:
        a = samples[(proc, e)]
        m = mask(a, e, ch, region)
        sc = a["MVA_Score_mA_M%d" % mA]
        if extra is not None:
            m &= extra(sc)
        if proc == "Data" and region == "incl":
            m &= sc < 0.9
        x = sc if var is None else a[var]
        hh, _ = np.histogram(x[m], edges, weights=a["w"][m]); h += hh
        ss, _ = np.histogram(x[m], edges, weights=a["w"][m] ** 2); s2 += ss
    return h, s2


def chi2(d, m, s2, sel):
    v = np.where(d > 0, d, m) + s2
    ok = sel & (v > 0)
    c = float(np.sum((d[ok] - m[ok]) ** 2 / v[ok])); n = int(ok.sum())
    return c, n, float(stats.chi2.sf(c, n)) if n else float("nan")


res = {"A": {}, "B": {}, "C": {}}
centers = 0.5 * (EDGES[1:] + EDGES[:-1])
splits = [("all", ERAS, None)] + [(e, [e], None) for e in ERAS] + [("ee", ERAS, "ee"), ("mumu", ERAS, "mumu")]
for mA in MASSES:
    for region in ("incl", "sb"):
        for name, eras, ch in splits:
            d, _ = hist("Data", mA, eras, ch, region)
            zg, zg2 = hist("Zg", mA, eras, ch, region)
            zj, zj2 = hist("Zj", mA, eras, ch, region)
            mc, s2 = zg + zj, zg2 + zj2
            top = 0.9 if region == "incl" else 1.0
            hi = (centers > HI) & (centers < top)
            c, n, p = chi2(d, mc, s2, hi)
            lo = centers < HI
            res["A"]["M%d/%s/%s" % (mA, region, name)] = {
                "norm_lo": float(d[lo].sum() / mc[lo].sum()),
                "ratio_hi": float(d[hi].sum() / mc[hi].sum()), "n_data_hi": float(d[hi].sum()),
                "mc_hi": float(mc[hi].sum()), "mc_hi_err": float(np.sqrt(s2[hi].sum())),
                "zg_frac_hi": float(zg[hi].sum() / mc[hi].sum()), "zg_frac_lo": float(zg[lo].sum() / mc[lo].sum()),
                "chi2_hi": c, "ndf_hi": n, "p_hi": p,
                "bins": {"%.3f" % x: [float(a), float(b), float(np.sqrt(e))] for x, a, b, e in zip(centers, d, mc, s2)}}
    # B: two-normalization fit on the sideband score distribution
    d, _ = hist("Data", mA, ERAS, None, "sb")
    zg, zg2 = hist("Zg", mA, ERAS, None, "sb"); zj, zj2 = hist("Zj", mA, ERAS, None, "sb")
    v = np.where(d > 0, d, 1) + zg2 + zj2
    X = np.vstack([zg, zj]).T / np.sqrt(v)[:, None]; y = d / np.sqrt(v)
    (a_g, a_j), *_ = np.linalg.lstsq(X, y, rcond=None)
    fit = a_g * zg + a_j * zj
    hi = centers > HI
    c, n, p = chi2(d, fit, a_g ** 2 * zg2 + a_j ** 2 * zj2, hi)
    c_all, n_all, p_all = chi2(d, fit, a_g ** 2 * zg2 + a_j ** 2 * zj2, np.ones_like(hi))
    res["B"]["M%d" % mA] = {"k_Zg": float(a_g), "k_Zj": float(a_j), "chi2_hi_after": c, "ndf_hi": n, "p_hi_after": p,
                            "chi2_all_after": c_all, "ndf_all": n_all - 2, "ratio_hi_after": float(d[hi].sum() / fit[hi].sum())}
    # C: input shapes in the sideband, score > 0.5 vs < 0.5
    for var in FEATS[mA]:
        allx = np.concatenate([samples[("Data", e)][var] for e in ERAS])
        allx = allx[np.isfinite(allx)]
        edges = np.unique(np.quantile(allx, np.linspace(0, 1, 16)))
        if len(edges) < 4:
            continue
        row = {}
        for tag, ext in (("hi", lambda s: s > HI), ("lo", lambda s: s < HI)):
            dd, _ = hist("Data", mA, ERAS, None, "sb", var, edges, ext)
            g, g2 = hist("Zg", mA, ERAS, None, "sb", var, edges, ext)
            j, j2 = hist("Zj", mA, ERAS, None, "sb", var, edges, ext)
            mc = g + j; k = dd.sum() / mc.sum()
            cc, nn, pp = chi2(dd, k * mc, k * k * (g2 + j2), np.ones(len(dd), bool))
            row[tag] = {"chi2": cc, "ndf": nn - 1, "p": pp, "norm": float(k), "n_data": float(dd.sum()),
                        "ratio_by_bin": [float(x) for x in dd / np.where(mc > 0, k * mc, np.nan)]}
        row["edges"] = [float(x) for x in edges]
        res["C"]["M%d/%s" % (mA, var)] = row

json.dump(res, open(os.path.join(OUT, "diag_highscore.json"), "w"), indent=1)
L = []
for mA in MASSES:
    L.append("=== m_a = %d GeV (%s model)" % (mA, "low-mass" if mA <= 3 else "high-mass"))
    for region in ("incl", "sb"):
        L.append("  [%s] split   norm(<0.5)  data/MC(>0.5)   n_data  MC+-stat       Zg frac lo/hi  chi2/ndf  p" % region)
        for name, _, _ in splits:
            r = res["A"]["M%d/%s/%s" % (mA, region, name)]
            L.append("   %-13s %6.3f      %6.3f        %6.0f  %7.1f+-%-5.1f   %4.2f/%4.2f    %5.1f/%-3d %.3g" % (
                name, r["norm_lo"], r["ratio_hi"], r["n_data_hi"], r["mc_hi"], r["mc_hi_err"], r["zg_frac_lo"],
                r["zg_frac_hi"], r["chi2_hi"], r["ndf_hi"], r["p_hi"]))
    b = res["B"]["M%d" % mA]
    L.append("  [fit sb] k(Zg)=%.3f k(Zj)=%.3f  -> data/fit(>0.5)=%.3f  chi2(>0.5)=%.1f/%d p=%.3g  chi2(all)=%.1f/%d" % (
        b["k_Zg"], b["k_Zj"], b["ratio_hi_after"], b["chi2_hi_after"], b["ndf_hi"], b["p_hi_after"], b["chi2_all_after"], b["ndf_all"]))
    L.append("  [inputs sb] var  chi2/ndf (score>0.5)  p  |  chi2/ndf (score<0.5)  p")
    rows = sorted(((v, res["C"][k]) for k, v in [(k, k.split("/")[1]) for k in res["C"] if k.startswith("M%d/" % mA)]),
                  key=lambda x: x[1]["hi"]["p"])
    for v, r in rows:
        L.append("   %-24s %6.1f/%-3d %.2g   |  %7.1f/%-3d %.2g" % (v, r["hi"]["chi2"], r["hi"]["ndf"], r["hi"]["p"],
                                                               r["lo"]["chi2"], r["lo"]["ndf"], r["lo"]["p"]))
open(os.path.join(OUT, "diag_highscore.txt"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
