#!/usr/bin/env python3
"""Signal BDT input variables before and after the sideband reweighting, per mass point.

Companion of plot_signal_reweight_mllgg_bdtwp.py (same events: the MVA-cut nominal trees that
enter the datacards, five Run 3 eras and both lepton channels summed, unit-area normalization so
that only the shape is compared; "before" = weight / weight_sideband_rwgt). The 16 variables (+ pho2_relEerr, low-mass input),
their ranges and axis titles are those of the data/simulation comparison in App. F
(1_prepare_dataVmc.py, Plot_Configs.py); param uses the true ALP mass, as in the training.

Input : /eos/home-p/pelai/HZa/root_MVAcut/sig/mA_M<ma>/output_<era>.root
Output: Plot/plots/signalReweight_bdtwp/trainvars_mA_M<ma>/<n>_<var>.pdf (+ .png)
        Plot/plots/signalReweight_bdtwp/trainvars_summary.txt

Usage : python plot_signal_reweight_trainvars_bdtwp.py [--ma 1,5,30]   (default: all 14 mass points)
"""

from __future__ import annotations

import argparse
import ast
import os
import sys
from pathlib import Path

import numpy as np
import uproot

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import ROOT  # noqa: E402
from _cms_style import LINE_WIDTH, cms_label, format_axes, make_legend  # noqa: E402

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

IN_BASE = Path("/eos/home-p/pelai/HZa/root_MVAcut/sig")
OUT_BASE = Path("/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/signalReweight_bdtwp")
CONSTANTS = ("/afs/cern.ch/work/p/pelai/HZgamma/higgsdna-hzg-run3/higgs_dna/"
             "metaconditions/corrections/constants.py")
ERAS = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]
CHANNELS = ["ele", "mu"]

# (name, branch or derived, nbins, lo, hi, x title) -- binning as in 1_prepare_dataVmc.py
VARS = [
    ("pho1Pt_oHm", "pho1Pt_oHm", 30, 0.0, 0.6, "p_{T,#gamma1}/m_{ll#gamma#gamma}"),
    ("pho1R9", "pho1R9", 25, 0.1, 1.0, "R_{9,#gamma1}"),
    ("pho1IetaIeta55", "pho1IetaIeta55", 25, 0.0, 0.08, "#sigma_{i#etai#eta, #gamma1}"),
    ("pho2Pt_oHm", "pho2Pt_oHm", 25, 0.0, 0.5, "p_{T,#gamma2}/m_{ll#gamma#gamma}"),
    ("pho2R9", "pho2R9", 25, 0.1, 1.0, "R_{9,#gamma2}"),
    ("pho2IetaIeta55", "pho2IetaIeta55", 25, 0.0, 0.08, "#sigma_{i#etai#eta, #gamma2}"),
    ("pho1ECALIso", "ALP_lead_photon_ecalPFClusterIso", 50, 0.0, 50.0, "I_{ECAL,#gamma1}"),
    ("pho2ECALIso", "ALP_sublead_photon_ecalPFClusterIso", 50, 0.0, 50.0, "I_{ECAL,#gamma2}"),
    ("ALP_calculatedPhotonIso", "ALP_calculatedPhotonIso", 25, 0.0, 125.0, "I_{#gamma,ALPs}"),
    ("var_dR_Za", "var_dR_Za", 35, 0.0, 7.0, "#DeltaR(Z,a)"),
    ("var_dR_g1g2", "var_dR_g1g2", 25, 0.0, 5.0, "#DeltaR(#gamma1,#gamma2)"),
    ("var_dR_g1Z", "var_dR_g1Z", 35, 0.0, 7.0, "#DeltaR(#gamma1,Z)"),
    ("var_PtaOverMh", "var_PtaOverMh", 25, 0.0, 0.75, "p_{T,a}/m_{ll#gamma#gamma}"),
    ("H_pt_oHm", "H_pt_oHm", 50, 0.0, 2.5, "p_{T,H}/m_{ll#gamma#gamma}"),
    ("pho_pt_asym", None, 40, 0.0, 1.0,
     "(p_{T,#gamma1}-p_{T,#gamma2})/(p_{T,#gamma1}+p_{T,#gamma2})"),
    # range set per mass point in main() (param_range): the signal is much narrower than in App. F
    ("param", None, 40, -0.01, 0.01, "(m_{a} - m_{a,hyp})/m_{ll#gamma#gamma}"),
    # low-mass BDT input (m_a = 1-3): sub-leading-photon sigma_E/E = energyErr / (pT cosh eta)
    ("pho2_relEerr", None, 55, 0.0, 0.55, "#sigma_{E}/E_{#gamma2}"),
]
BRANCHES = sorted({v[1] for v in VARS if v[1]} | {"H_m", "ALP_m", "weight", "weight_sideband_rwgt",
                                                  "ALP_sublead_photon_energyErr", "ALP_sublead_photon_pt", "ALP_sublead_photon_eta"})


def run3_lumi() -> float:
    # The file cannot be imported (its ERA dict uses sets as keys): parse only LUMI.
    tree = ast.parse(Path(CONSTANTS).read_text())
    for node in tree.body:
        if isinstance(node, ast.Assign) and getattr(node.targets[0], "id", None) == "LUMI":
            lumi = ast.literal_eval(node.value)
            return sum(lumi[e] for e in ERAS)
    raise RuntimeError(f"LUMI not found in {CONSTANTS}")


def load(ma: int):
    parts = []
    for era in ERAS:
        f = uproot.open(IN_BASE / f"mA_M{ma}" / f"output_{era}.root")
        for ch in CHANNELS:
            parts.append(f[f"DiphotonTree/ggh_125_Za_{ch}_13p6TeV_cat0"].arrays(BRANCHES, library="np"))
    d = {k: np.concatenate([p[k] for p in parts]).astype(np.float64) for k in BRANCHES}
    p1, p2 = d["pho1Pt_oHm"], d["pho2Pt_oHm"]
    d["pho_pt_asym"] = (p1 - p2) / (p1 + p2 + 1e-6)
    d["param"] = (d["ALP_m"] - float(ma)) / d["H_m"]
    d["pho2_relEerr"] = d["ALP_sublead_photon_energyErr"] / (d["ALP_sublead_photon_pt"] * np.cosh(d["ALP_sublead_photon_eta"]))
    ok = (d["weight_sideband_rwgt"] != 0) & (d["weight"] != 0)  # zero-weight events carry no information
    wa = d["weight"][ok]
    wb = wa / d["weight_sideband_rwgt"][ok]
    return {k: v[ok] for k, v in d.items()}, wa, wb


def hist(name, x, w, edges):
    h = ROOT.TH1D(name, "", len(edges) - 1, edges)
    h.Sumw2()
    xc = np.clip(x, edges[0], edges[-1] - 1e-9 * (edges[-1] - edges[0]))  # overflow into last bin
    h.FillN(len(xc), np.ascontiguousarray(xc), np.ascontiguousarray(w))
    h.SetDirectory(0)
    h.Scale(1.0 / h.Integral())
    return h


def ratio(name, x, wa, wb, edges, min_frac=1e-3):
    """After/before per bin (both unit area): the wb-weighted mean of r = wa/wb in the bin over
    the global one, with the standard error of that weighted mean; bins with less than
    min_frac of the signal are left empty."""
    h = ROOT.TH1D(name, "", len(edges) - 1, edges)
    h.SetDirectory(0)
    norm = wb.sum() / wa.sum()
    xc = np.clip(x, edges[0], edges[-1] - 1e-9 * (edges[-1] - edges[0]))
    idx = np.digitize(xc, edges) - 1
    r = wa / wb
    for i in range(len(edges) - 1):
        sel = idx == i
        sw = wb[sel].sum()
        if sel.sum() < 2 or sw <= 0 or sw < min_frac * wb.sum():
            continue
        mean = np.sum(wb[sel] * r[sel]) / sw
        var = np.sum(wb[sel] ** 2 * (r[sel] - mean) ** 2) / sw ** 2
        h.SetBinContent(i + 1, mean * norm)
        h.SetBinError(i + 1, np.sqrt(var) * norm)
    return h


def draw(ma, tag, xtitle, hb, ha, hr, lumi, out_path: Path):
    canvas = ROOT.TCanvas(f"c_{tag}", "", 800, 800)
    pad1 = ROOT.TPad("pad1", "", 0, 0.30, 1, 1)
    pad2 = ROOT.TPad("pad2", "", 0, 0.0, 1, 0.30)
    for p in (pad1, pad2):
        p.SetLeftMargin(0.17)
    pad1.SetBottomMargin(0.02)
    pad2.SetTopMargin(0.03)
    pad2.SetBottomMargin(0.40)
    pad1.Draw()
    pad2.Draw()

    pad1.cd()
    for h, color in ((hb, ROOT.kBlack), (ha, ROOT.kRed + 1)):
        h.SetLineColor(color)
        h.SetLineWidth(LINE_WIDTH)
    hb.GetYaxis().SetTitle("Fraction of events")
    format_axes(hb, title_size=0.055 / 0.7, label_size=0.05 / 0.7)
    hb.GetXaxis().SetTitleSize(0)
    hb.GetXaxis().SetLabelSize(0)
    hb.GetYaxis().SetTitleOffset(1.2)
    hb.SetMaximum(1.9 * max(hb.GetMaximum(), ha.GetMaximum()))
    hb.SetMinimum(0.0)
    hb.Draw("HIST")
    ha.Draw("HIST SAME")

    leg = make_legend(0.44, 0.58, 0.82, 0.72, text_size=0.045 / 0.7)
    leg.AddEntry(hb, "Before reweighting", "l")
    leg.AddEntry(ha, "After reweighting", "l")
    leg.Draw()
    txt = ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextFont(42)
    txt.SetTextSize(0.045 / 0.7)
    txt.DrawLatex(0.21, 0.84, f"H#rightarrowZa, m_{{a}} = {ma} GeV")
    txt.DrawLatex(0.21, 0.76, "After the BDT selection")
    labels = cms_label(lumi=lumi, x_left=0.17, cms_text="CMS Simulation",
                       text_size=0.045 / 0.7, y=0.94)

    pad2.cd()
    hr.SetLineColor(ROOT.kRed + 1)
    hr.SetMarkerColor(ROOT.kRed + 1)
    hr.SetMarkerStyle(20)
    hr.SetMarkerSize(0.8)
    hr.SetLineWidth(LINE_WIDTH)
    hr.GetYaxis().SetTitle("After / before")
    hr.GetXaxis().SetTitle(xtitle)
    format_axes(hr, title_size=0.055 / 0.3, label_size=0.05 / 0.3)
    hr.GetYaxis().SetTitleOffset(0.52)
    hr.GetXaxis().SetTitleOffset(0.95)
    hr.GetYaxis().SetNdivisions(505)
    hr.GetXaxis().SetNdivisions(505)
    hr.SetMinimum(0.5)
    hr.SetMaximum(1.5)
    hr.Draw("E1")
    edges = hr.GetXaxis()
    line = ROOT.TLine(edges.GetXmin(), 1.0, edges.GetXmax(), 1.0)
    line.SetLineStyle(2)
    line.Draw("SAME")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    canvas.SaveAs(str(out_path))
    canvas.SaveAs(str(out_path.with_suffix(".png")))
    canvas.Close()
    del labels, leg, line


MASSES = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30]
# variables whose mean is ~0 for the signal: a relative mean shift is not meaningful
NO_REL_MEAN = {"ALP_calculatedPhotonIso", "param"}


def param_range(x, w):
    """Symmetric range covering 99% of the signal: param (true mass) is narrow and its width
    grows with m_a, so the App. F range would put the whole signal into one or two bins."""
    o = np.argsort(x)
    c = np.cumsum(w[o]) / w.sum()
    lo, hi = np.interp([0.005, 0.995], c, x[o])
    lim = np.ceil(1.2 * max(abs(lo), abs(hi)) * 1000.0) / 1000.0
    return -lim, lim


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ma", default="all")
    parser.add_argument("--out", default=str(OUT_BASE))
    args = parser.parse_args()
    masses = MASSES if args.ma == "all" else [int(v) for v in args.ma.split(",")]
    lumi = run3_lumi()
    lines = ["ma  variable                  range              mean_shift[%]  max|ratio-1| (bins >1% of signal)"]
    for ma in masses:
        d, wa, wb = load(ma)
        out_dir = Path(args.out) / f"trainvars_mA_M{ma}"
        print(f"m_a = {ma} GeV: {len(wa)} events, yield ratio after/before {wa.sum() / wb.sum():.4f}")
        for n, (name, branch, nb, lo, hi, xtitle) in enumerate(VARS, start=1):
            x = d[branch or name]
            if name == "param":
                lo, hi = param_range(x, wb)
            edges = np.linspace(lo, hi, nb + 1)
            hb = hist(f"hb_{ma}_{name}", x, wb, edges)
            ha = hist(f"ha_{ma}_{name}", x, wa, edges)
            hr = ratio(f"hr_{ma}_{name}", x, wa, wb, edges)
            big = [abs(hr.GetBinContent(i) - 1) for i in range(1, nb + 1)
                   if hb.GetBinContent(i) > 0.01 and hr.GetBinContent(i) > 0]
            mb, mafter = np.average(x, weights=wb), np.average(x, weights=wa)
            rel = "%+.2f" % (100.0 * (mafter / mb - 1.0)) if name not in NO_REL_MEAN and mb else "-"
            lines.append("%2d  %-24s [%7.3f,%7.3f]  %13s  %.3f" % (ma, name, lo, hi, rel, max(big) if big else 0.0))
            draw(ma, f"{ma}_{name}", xtitle, hb, ha, hr, lumi, out_dir / f"{n}_{name}.pdf")
        print("wrote", out_dir)
    out = Path(args.out) / "trainvars_summary.txt"
    out.write_text("\n".join(lines) + "\n")
    print("wrote", out)


if __name__ == "__main__":
    main()
