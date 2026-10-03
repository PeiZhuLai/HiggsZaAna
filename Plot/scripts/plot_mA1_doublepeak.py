#!/usr/bin/env python3
"""Origin of the double m_gammagamma peak of the m_a = 1 GeV signal (reviewer question, v1 Fig. 33).

The reconstructed m_gammagamma of the m_a = 1 GeV signal shows a main peak at 1 GeV and a second
one at about 1.3-1.5 GeV, and m_llgammagamma a second population at 134-137 GeV. This script
splits the selected signal events (FSR-fix production, all Run 3 eras, adopted BDT working point
of m_a = 1 GeV) by the generator-level opening angle of the two ALP photons and measures the
photon-pair energy response, to test whether the second peak comes from overlapping showers of
correctly matched photons.

Input : run3_bdt_scored_fsrfix/mA_M1/<era>.root (tree inclusive), weight = factor x normalized signal reweight
        working point from Plot/output/MVAcut_points_run3.json
Output: Plot/plots/mA1_doublepeak/{mgg_by_gendR,mllgg_by_gendR,response_vs_gendR}.pdf/.png
        and summary.txt (fractions and medians quoted in the AN / review answer)
Env   : higgs-alp-ana python (uproot + PyROOT)
"""
import ast
import json
import os

import numpy as np
import ROOT
import uproot

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

BASE = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix/mA_M1"
ERAS = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]
WP_JSON = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/MVAcut_points_run3.json"
OUT = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/mA1_doublepeak"
CONST = "/afs/cern.ch/work/p/pelai/HZgamma/higgsdna-hzg-run3/higgs_dna/metaconditions/corrections/constants.py"
DR_SPLIT = 0.05   # ~ three ECAL crystals: below it the two photon showers overlap

def _read_lumi(path):
    """LUMI dict from the authoritative constants.py, read with ast (the module itself does not
    import under this python)."""
    for node in ast.parse(open(path).read()).body:
        if isinstance(node, ast.Assign) and any(getattr(t, "id", "") == "LUMI" for t in node.targets):
            return ast.literal_eval(node.value)
    raise RuntimeError("LUMI not found in %s" % path)


_L = _read_lumi(CONST)
LUMI = _L["2022"] + _L["2023"] + _L["2024"]
LUMI_LABEL = "%.2f fb^{-1} (13.6 TeV)" % LUMI

BR = ["H_mass", "ALP_mass", "factor", "MVA_Score_mA_M1", "GenALP_dR_gg",
      "ALP_lead_photon_pt", "ALP_lead_photon_eta", "ALP_lead_photon_phi",
      "ALP_sublead_photon_pt", "ALP_sublead_photon_eta", "ALP_sublead_photon_phi",
      "GenALPLeadPho_pt", "GenALPLeadPho_eta", "GenALPLeadPho_phi",
      "GenALPSubleadPho_pt", "GenALPSubleadPho_eta", "GenALPSubleadPho_phi"]


def dr(eta1, phi1, eta2, phi2):
    dphi = np.mod(phi1 - phi2 + np.pi, 2 * np.pi) - np.pi
    return np.hypot(eta1 - eta2, dphi)


# [PZ 2026-10-03] The signal is reweighted with the nominal sideband reweight, as in apply_bdt_sig.py and
# signal_eff_sumw.py: true-mass param (ALP_m - m_a)/H_m, normalized per (era, channel) so that the
# preselected yield is unchanged. HZA_SIGNAL_REWEIGHT=0 restores the unreweighted signal.
SIGNAL_RW_JSON = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/reweights/sideband_run3_iterative.json"
_RW = None


def signal_rw(path, tree, ma):
    """Normalized per-event signal reweight, aligned with the entries of <tree> in <path>."""
    global _RW
    if os.environ.get("HZA_SIGNAL_REWEIGHT", "1") == "0":
        return None
    if _RW is None:
        import sys
        sys.path.insert(0, "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts")
        from sideband_reweight import SidebandReweighter
        _RW = SidebandReweighter.from_json(SIGNAL_RW_JSON)
    fr = uproot.open(path)[tree].arrays(library="pd")
    for logical, cands in (("pho1ECALIso", ("pho1PIso_noCorr",)), ("pho2ECALIso", ("pho2PIso_noCorr",)),
                           ("H_m", ("H_mass",)), ("ALP_m", ("ALP_mass",))):
        if logical not in fr.columns:
            for c in cands:
                if c in fr.columns:
                    fr[logical] = fr[c]; break
    fr["param"] = (fr["ALP_m"].to_numpy(dtype=float) - float(ma)) / fr["H_m"].to_numpy(dtype=float)
    r = np.asarray(_RW.weights_for_dataframe(fr), dtype=float)
    w = fr["factor"].to_numpy(dtype=float)
    out = np.ones(len(fr), dtype=float)
    for col in ("n_electrons", "n_muons"):
        sel = fr[col].to_numpy() == 2
        den = float(np.sum(w[sel] * r[sel]))
        out[sel] = r[sel] * (float(np.sum(w[sel])) / den if den > 0 else 1.0)
    return out


def load(cut):
    parts = []
    for era in ERAS:
        path = "%s/%s.root" % (BASE, era)
        a = uproot.open(path)["inclusive"].arrays(BR, library="np")
        rw = signal_rw(path, "inclusive", 1)
        if rw is not None:
            a["factor"] = a["factor"] * rw
        sel = (a["MVA_Score_mA_M1"] > cut) & (a["H_mass"] > 95) & (a["H_mass"] < 180)
        parts.append({k: v[sel] for k, v in a.items()})
    return {k: np.concatenate([p[k] for p in parts]) for k in BR}


def wmedian(x, w):
    o = np.argsort(x); x, w = x[o], w[o]; c = np.cumsum(w)
    return float(x[np.searchsorted(c, 0.5 * c[-1])])


def style(h, col):
    h.SetLineColor(col); h.SetLineWidth(2); h.SetMarkerColor(col)
    for ax in (h.GetXaxis(), h.GetYaxis()):
        ax.SetTitleSize(0.055); ax.SetLabelSize(0.050)
    h.GetYaxis().SetTitleOffset(1.25)


def labels():
    t = ROOT.TLatex(); t.SetNDC(); t.SetTextSize(0.045)
    t.DrawLatex(0.16, 0.92, "#bf{CMS} #it{Preliminary}")
    t.SetTextAlign(31); t.DrawLatex(0.95, 0.92, LUMI_LABEL)
    return t


def overlay(name, xtitle, xs, ws, masks, nb, lo, hi):
    c = ROOT.TCanvas(name, "", 800, 700)
    c.SetLeftMargin(0.16); c.SetBottomMargin(0.14); c.SetTopMargin(0.09)
    hs = []
    specs = [("all selected events", ROOT.kBlack, np.ones_like(xs, bool)),
             ("gen #DeltaR(#gamma,#gamma) < %.2f" % DR_SPLIT, ROOT.kRed + 1, masks[0]),
             ("gen #DeltaR(#gamma,#gamma) #geq %.2f" % DR_SPLIT, ROOT.kAzure + 2, masks[1])]
    for i, (lab, col, m) in enumerate(specs):
        h = ROOT.TH1D("%s_%d" % (name, i), "", nb, lo, hi)
        for x, w in zip(xs[m], ws[m]):
            h.Fill(x, w)
        style(h, col); h.GetXaxis().SetTitle(xtitle); h.GetYaxis().SetTitle("Events / %.3g GeV" % ((hi - lo) / nb))
        hs.append((h, lab))
    ymax = max(h.GetMaximum() for h, _ in hs)
    for i, (h, _) in enumerate(hs):
        h.SetMaximum(1.35 * ymax); h.Draw("hist" if i == 0 else "hist same")
    leg = ROOT.TLegend(0.50, 0.66, 0.93, 0.87); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.037)
    leg.SetHeader("m_{a} = 1 GeV, BDT working point")
    for h, lab in hs:
        leg.AddEntry(h, lab, "l")
    leg.Draw(); t = labels()
    for ext in ("pdf", "png"):
        c.SaveAs("%s/%s.%s" % (OUT, name, ext))


def main():
    os.makedirs(OUT, exist_ok=True)
    cut = {e["mA"]: e["MVAcut"] for e in json.load(open(WP_JSON))["results"]}[1]
    a = load(cut)
    w = a["factor"]
    gdr = a["GenALP_dR_gg"]
    reco = a["ALP_lead_photon_pt"] + a["ALP_sublead_photon_pt"]
    gen = a["GenALPLeadPho_pt"] + a["GenALPSubleadPho_pt"]
    resp = np.where(gen > 0, reco / np.where(gen > 0, gen, 1), np.nan)
    # one-to-one matching of the two reco photons to the two gen ALP photons (either assignment)
    d11 = dr(a["ALP_lead_photon_eta"], a["ALP_lead_photon_phi"], a["GenALPLeadPho_eta"], a["GenALPLeadPho_phi"])
    d22 = dr(a["ALP_sublead_photon_eta"], a["ALP_sublead_photon_phi"], a["GenALPSubleadPho_eta"], a["GenALPSubleadPho_phi"])
    d12 = dr(a["ALP_lead_photon_eta"], a["ALP_lead_photon_phi"], a["GenALPSubleadPho_eta"], a["GenALPSubleadPho_phi"])
    d21 = dr(a["ALP_sublead_photon_eta"], a["ALP_sublead_photon_phi"], a["GenALPLeadPho_eta"], a["GenALPLeadPho_phi"])
    matched = (np.maximum(d11, d22) < 0.02) | (np.maximum(d12, d21) < 0.02)   # tighter than the photon separation
    ok = np.isfinite(resp) & (gdr > 0)
    close, wide = ok & (gdr < DR_SPLIT), ok & (gdr >= DR_SPLIT)
    main_pk, second_pk = ok & (a["ALP_mass"] < 1.15), ok & (a["ALP_mass"] >= 1.15)
    core, tail = ok & (a["H_mass"] > 120) & (a["H_mass"] < 130), ok & (a["H_mass"] > 134) & (a["H_mass"] < 137)

    lines = ["m_a = 1 GeV signal, FSR-fix production, eras %s, BDT > %.3f, 95 < m_llgg < 180" % (",".join(ERAS), cut),
             "weighted events: %.2f (raw %d)" % (w.sum(), w.size), ""]
    def frac(m, sub):
        return w[m & sub].sum() / w[m].sum()
    for lab, m in (("all", ok), ("m_gg < 1.15 (main peak)", main_pk), ("m_gg >= 1.15 (second peak)", second_pk),
                   ("m_llgg 120-130 (core)", core), ("m_llgg 134-137 (second population)", tail)):
        lines.append("%-36s w=%7.2f  frac(gen dR<%.2f)=%.2f  frac(1-to-1 matched)=%.2f  median response=%.3f  median gen dR=%.3f"
                     % (lab, w[m].sum(), DR_SPLIT, frac(m, close), frac(m, matched), wmedian(resp[m], w[m]), wmedian(gdr[m], w[m])))
    lines.append("")
    for lab, m in (("gen dR < %.2f" % DR_SPLIT, close), ("gen dR >= %.2f" % DR_SPLIT, wide)):
        lines.append("%-20s frac(m_gg >= 1.15)=%.2f  median response=%.3f  median m_gg=%.3f"
                     % (lab, frac(m, second_pk), wmedian(resp[m], w[m]), wmedian(a["ALP_mass"][m], w[m])))
    open("%s/summary.txt" % OUT, "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))

    overlay("mgg_by_gendR", "m_{#gamma#gamma} [GeV]", a["ALP_mass"][ok], w[ok], (close[ok], wide[ok]), 40, 0.5, 2.5)
    overlay("mllgg_by_gendR", "m_{ll#gamma#gamma} [GeV]", a["H_mass"][ok], w[ok], (close[ok], wide[ok]), 40, 110, 145)
    overlay("response_vs_gendR", "p_{T}^{reco}(#gamma#gamma) / p_{T}^{gen}(#gamma#gamma)", resp[ok], w[ok], (close[ok], wide[ok]), 40, 0.6, 2.2)


if __name__ == "__main__":
    main()
