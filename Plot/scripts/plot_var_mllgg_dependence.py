#!/usr/bin/env python3
"""m_llgg dependence of BDT input candidates in the simulated background (AN appendix on the low-mass input choice).

Two figures from the background simulation (All_Bkg, 'test' tree of the BDT inputs, nominal weights clipped at 0,
95 < m_llgg < 180 GeV):
  corr_matrix_bkg.pdf      weighted Pearson correlation matrix of m_llgg, the inputs rejected for the low-mass model
                           (sculpting) and the inputs of the current low-mass model
  median_vs_mllgg_bkg.pdf  weighted median of each variable in 5 GeV bins of m_llgg divided by its median over the
                           whole range: rejected inputs solid, current low-mass inputs dashed
A linear correlation can stay small while the median moves with m_llgg (the leading-photon sigma_E/E rises by about
half from 95 to 180 GeV for a correlation of 0.15), so both views are shown. param = (m_gg - m_a)/m_llgg is not drawn:
it is divided by m_llgg by construction and depends on the mass hypothesis.
Numbers: median_vs_mllgg_bkg.txt. Env: higgs-alp-ana (PyROOT). Output: Plot/plots/var_mllgg_dependence/
"""
import ast
import os
import sys

import numpy as np
import ROOT

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "lib"))
sys.path.insert(0, "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts")
from hza_features import ROOT_DIR, _load, feature_matrix, relerr_extra  # noqa: E402
from Plot_Helper import MakeLumiLabel  # noqa: E402

# luminosities from the authoritative HZgamma constants.py; only the LUMI literal is parsed (the module itself does not
# import: its ERA dict has set keys)
CONSTANTS = "/afs/cern.ch/work/p/pelai/HZgamma/higgsdna-hzg-run3/higgs_dna/metaconditions/corrections/constants.py"
LUMI = next(ast.literal_eval(n.value) for n in ast.parse(open(CONSTANTS).read()).body
            if isinstance(n, ast.Assign) and getattr(n.targets[0], "id", "") == "LUMI")

OUT = os.path.join(os.path.dirname(HERE), "plots", "var_mllgg_dependence")
ERAS = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]
REJECTED = [("pho1_relEerr", "#sigma_{E}/E(#gamma_{1})"), ("pho1Pt_oHm", "p_{T}^{#gamma_{1}}/m_{ll#gamma#gamma}"),
            ("pho2Pt_oHm", "p_{T}^{#gamma_{2}}/m_{ll#gamma#gamma}"), ("H_pt_oHm", "p_{T}^{ll#gamma#gamma}/m_{ll#gamma#gamma}"),
            ("var_PtaOverMh", "p_{T}^{a}/m_{ll#gamma#gamma}"), ("pho_pt_asym", "p_{T} asymmetry"),
            ("ALP_lead_photon_etaWidth", "#eta width(#gamma_{1})")]
CURRENT = [("pho2_relEerr", "#sigma_{E}/E(#gamma_{2})"), ("ALP_calculatedPhotonIso", "photon isolation"),
           ("var_dR_g1g2", "#DeltaR(#gamma_{1},#gamma_{2})"), ("pho1R9", "R_{9}(#gamma_{1})")]
EDGES = np.arange(95.0, 180.1, 5.0)


def wmedian(x, w):
    o = np.argsort(x)
    c = np.cumsum(w[o])
    return float(x[o][np.searchsorted(c, 0.5 * c[-1])])


def wcorr(x, y, w):
    mx, my = np.average(x, weights=w), np.average(y, weights=w)
    return float(np.average((x - mx) * (y - my), weights=w) / np.sqrt(np.average((x - mx) ** 2, weights=w) * np.average((y - my) ** 2, weights=w)))


def cms_label(pad):
    t = ROOT.TLatex()
    t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.045); t.SetTextAlign(11)
    t.DrawLatex(pad.GetLeftMargin(), 0.91, "#bf{CMS} #it{Preliminary}")
    return t


def axis_style(h):
    for ax in (h.GetXaxis(), h.GetYaxis()):
        ax.SetTitleSize(0.055); ax.SetLabelSize(0.050)


def main():
    ROOT.gROOT.SetBatch(True); ROOT.gStyle.SetOptStat(0); ROOT.gStyle.SetOptTitle(0)
    os.makedirs(OUT, exist_ok=True)
    feats = [f for f, _ in REJECTED + CURRENT]
    raw = ["ALP_lead_photon_etaWidth", "ALP_calculatedPhotonIso", "var_dR_g1g2", "pho1R9", "H_pt_oHm", "var_PtaOverMh"]
    a = _load("%s/All_Bkg/run3.root" % ROOT_DIR, "test", None, extra=list(dict.fromkeys(relerr_extra(feats) + raw)))
    X = feature_matrix(a, 0.0, feats)
    m = a["H_m"]; w = np.clip(a["factor"], 0, None)
    ok = np.isfinite(X).all(axis=1) & (w > 0)
    X, m, w = X[ok], m[ok], w[ok]
    lumi = sum(LUMI[e] for e in ERAS)

    # correlation matrix: m_llgg first, then rejected, then current inputs
    names = ["m_{ll#gamma#gamma}"] + [l for _, l in REJECTED + CURRENT]
    cols = [m] + [X[:, i] for i in range(X.shape[1])]
    n = len(cols)
    h = ROOT.TH2D("corr", "", n, 0, n, n, 0, n)
    for i in range(n):
        for j in range(n):
            h.SetBinContent(i + 1, n - j, round(wcorr(cols[i], cols[j], w), 2))
        h.GetXaxis().SetBinLabel(i + 1, names[i]); h.GetYaxis().SetBinLabel(n - i, names[i])
    c = ROOT.TCanvas("c_corr", "", 1000, 900)
    c.SetLeftMargin(0.24); c.SetBottomMargin(0.24); c.SetRightMargin(0.13); c.SetTopMargin(0.10)
    ROOT.gStyle.SetPaintTextFormat(".2f"); h.SetMarkerSize(0.9)
    h.GetZaxis().SetRangeUser(-1, 1); axis_style(h)
    h.GetXaxis().SetLabelSize(0.030); h.GetYaxis().SetLabelSize(0.030); h.GetXaxis().LabelsOption("v")
    h.Draw("COLZ TEXT")
    keep = [cms_label(c), MakeLumiLabel(round(lumi, 2))]
    c.SaveAs(os.path.join(OUT, "corr_matrix_bkg.pdf"))

    # median vs m_llgg, normalized to the median over the whole range
    ctr = 0.5 * (EDGES[:-1] + EDGES[1:])
    lines = ["m_llgg bin centers: " + " ".join("%.1f" % x for x in ctr)]
    c2 = ROOT.TCanvas("c_med", "", 900, 800)
    c2.SetLeftMargin(0.16); c2.SetBottomMargin(0.14); c2.SetRightMargin(0.04); c2.SetTopMargin(0.10)
    frame = ROOT.TH1D("frame", ";m_{ll#gamma#gamma} (GeV);median / median over 95-180 GeV", 1, 95, 180)
    frame.SetMinimum(0.5); frame.SetMaximum(1.9); axis_style(frame); frame.GetYaxis().SetTitleOffset(1.3)
    frame.Draw()
    leg = ROOT.TLegend(0.19, 0.62, 0.95, 0.88); leg.SetNColumns(2); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.030)
    greys = [ROOT.kBlack, ROOT.kGray + 3, ROOT.kGray + 2, ROOT.kGray + 1]
    colors = [ROOT.kRed + 1, ROOT.kOrange + 7, ROOT.kMagenta + 1, ROOT.kBlue + 1, ROOT.kCyan + 2, ROOT.kGreen + 2, ROOT.kViolet + 2]
    graphs = []
    for k, (f, label) in enumerate(REJECTED + CURRENT):
        x = X[:, feats.index(f)]
        tot = wmedian(x, w)
        y = []
        for lo, hi in zip(EDGES[:-1], EDGES[1:]):
            s = (m >= lo) & (m < hi)
            y.append(wmedian(x[s], w[s]) / tot if s.any() and tot != 0 else float("nan"))
        rej = k < len(REJECTED)
        g = ROOT.TGraph(len(ctr), np.array(ctr, dtype=float), np.array(y, dtype=float))
        g.SetLineWidth(3 if rej else 2); g.SetLineStyle(1 if rej else 2)
        g.SetLineColor(colors[k] if rej else greys[k - len(REJECTED)])
        g.SetMarkerStyle(20 if rej else 24); g.SetMarkerSize(0.8); g.SetMarkerColor(g.GetLineColor())
        g.Draw("LP SAME"); graphs.append(g)
        leg.AddEntry(g, label + ("" if rej else " (current)"), "lp")
        lines.append("%-26s %s  median %.4g  first/last bin %.3f / %.3f  corr(m_llgg) %+.3f" % (
            f, "rejected" if rej else "current ", tot, y[0], y[-1], wcorr(m, x, w)))
    leg.Draw()
    keep += [cms_label(c2), MakeLumiLabel(round(lumi, 2))]
    c2.SaveAs(os.path.join(OUT, "median_vs_mllgg_bkg.pdf"))
    open(os.path.join(OUT, "median_vs_mllgg_bkg.txt"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
