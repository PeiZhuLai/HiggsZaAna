#!/usr/bin/env python3
"""m_a = 3 GeV high-score excess concentrated in the ee channel: which variables differ?
Sideband (H_m outside 115-135, not blinded), BDT score for m_a = 3 above / below 0.5, background MC
without the BDT training events (bkg_split friend from make_bkg_split_friends.py, weight x2), MC
scaled to data in each (channel, score) slice. chi2 of the shape for the low-mass inputs and for
lepton / photon kinematics, ee vs mumu. Read-only; writes diag_mA3_ee.txt next to this script.
"""
import os
import numpy as np
import uproot
from scipy import stats

B = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix"
F = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix_bkgsplit"
ERAS = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]
BKG = {e: ["DYGto2LG_10to100", "DYJetsToLL"] for e in ERAS[:4]}
BKG["2024"] = ["DYGto2LG_10to100", "DYJetsTo2E", "DYJetsTo2Mu", "DYJetsTo2Tau"]
RAW = ["H_m", "z_mumu", "z_ee", "MVA_Score_mA_M3", "ALP_calculatedPhotonIso", "var_dR_g1g2", "pho1R9",
       "ALP_sublead_photon_energyErr", "ALP_sublead_photon_pt", "ALP_sublead_photon_eta", "ALP_lead_photon_pt",
       "ALP_lead_photon_eta", "Z_lead_lepton_pt", "Z_lead_lepton_eta", "Z_sublead_lepton_pt", "Z_sublead_lepton_eta",
       "Z_m", "pho2R9", "pho1IetaIeta55", "pho2IetaIeta55", "pho1PIso_noCorr", "pho2PIso_noCorr",
       "ALP_lead_photon_electronVeto", "ALP_sublead_photon_electronVeto", "ALP_lead_photon_lep_near_dR",
       "ALP_sublead_photon_lep_near_dR", "var_dR_Za", "H_pt_oHm"]
VARS = ["ALP_calculatedPhotonIso", "var_dR_g1g2", "pho1R9", "pho2_relEerr",
        "ALP_lead_photon_pt", "ALP_sublead_photon_pt", "ALP_lead_photon_eta", "ALP_sublead_photon_eta",
        "Z_lead_lepton_pt", "Z_sublead_lepton_pt", "Z_lead_lepton_eta", "Z_sublead_lepton_eta", "Z_m",
        "pho2R9", "pho1IetaIeta55", "pho2IetaIeta55", "pho1PIso_noCorr", "pho2PIso_noCorr",
        "ALP_lead_photon_electronVeto", "ALP_sublead_photon_electronVeto", "ALP_lead_photon_lep_near_dR",
        "ALP_sublead_photon_lep_near_dR", "var_dR_Za", "H_pt_oHm"]
OUT = os.path.dirname(os.path.abspath(__file__))


def load(path, friend=None):
    t = uproot.open(path)["inclusive"]
    a = t.arrays([b for b in RAW if b in t.keys()], library="np")
    a["w"] = np.ones(len(a["H_m"]))
    if friend:
        a["w"] = t["weight_sideband_rwgt"].array(library="np").astype(float)
        sp = uproot.open(friend)["bkgsplit"]["bkg_split"].array(library="np")
        assert len(sp) == len(a["w"])
        a["w"] = a["w"] * np.where((sp == 1) | (sp == 2), 2.0, 0.0)
    a["pho2_relEerr"] = a["ALP_sublead_photon_energyErr"] / (a["ALP_sublead_photon_pt"] * np.cosh(a["ALP_sublead_photon_eta"]))
    m = (a["H_m"] > 95) & (a["H_m"] < 180) & ((a["H_m"] < 115) | (a["H_m"] > 135))
    return {k: v[m] for k, v in a.items()}


def cat(xs):
    return {k: np.concatenate([x[k] for x in xs]) for k in xs[0]}


data = cat([load("%s/Data/%s.root" % (B, e)) for e in ERAS])
mc = cat([load("%s/%s/%s.root" % (B, s, e), "%s/%s/%s.root" % (F, s, e)) for e in ERAS for s in BKG[e]])


def chan(a, ch):
    return (np.abs(a["z_mumu"]) != 1) if ch == "ee" else (np.abs(a["z_ee"]) != 1)


L = ["mA3, sideband, MC = non-training events x2, MC scaled to data in each slice; shape chi2/ndf and p",
     "%-32s | %-22s %-22s | %-22s %-22s" % ("variable", "ee score>0.5", "ee score<0.5", "mumu score>0.5", "mumu score<0.5")]
rows = []
for v in VARS:
    x_all = data[v][np.isfinite(data[v])]
    edges = np.unique(np.quantile(x_all, np.linspace(0, 1, 13)))
    if len(edges) < 4:
        edges = np.unique(np.concatenate([edges, [x_all.min() - 1e-6, x_all.max() + 1e-6]]))
    cells = []; worst = 1.0
    for ch in ("ee", "mumu"):
        for tag, cut in (("hi", lambda s: s > 0.5), ("lo", lambda s: s < 0.5)):
            dm = chan(data, ch) & cut(data["MVA_Score_mA_M3"]); mm = chan(mc, ch) & cut(mc["MVA_Score_mA_M3"])
            d, _ = np.histogram(data[v][dm], edges)
            m, _ = np.histogram(mc[v][mm], edges, weights=mc["w"][mm]); m2, _ = np.histogram(mc[v][mm], edges, weights=mc["w"][mm] ** 2)
            k = d.sum() / m.sum(); m, m2 = m * k, m2 * k * k
            var = np.where(d > 0, d, m) + m2; ok = var > 0
            c = float(np.sum((d[ok] - m[ok]) ** 2 / var[ok])); n = int(ok.sum()) - 1
            p = float(stats.chi2.sf(c, n)) if n > 0 else float("nan")
            cells.append("%6.1f/%-3d p=%-8.2g" % (c, n, p))
            if ch == "ee" and tag == "hi":
                worst = p
    rows.append((worst, "%-32s | %s %s | %s %s" % (v, *cells)))
for _, r in sorted(rows):
    L.append(r)
nd = {ch: int((chan(data, ch) & (data["MVA_Score_mA_M3"] > 0.5)).sum()) for ch in ("ee", "mumu")}
L.insert(1, "data events with score>0.5 in the sideband: ee %d, mumu %d" % (nd["ee"], nd["mumu"]))
open(os.path.join(OUT, "diag_mA3_ee.txt"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
