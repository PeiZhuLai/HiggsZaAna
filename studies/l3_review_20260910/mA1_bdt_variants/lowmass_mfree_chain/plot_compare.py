#!/usr/bin/env python3
"""m_llgg shapes, production low-mass BDT vs win_mfree, m_a = 1, 2, 3 GeV, at the working points of
compare_models.json (production WP; new model at equal data-sideband count and at the R = 1 crossing).
The preselection shape (no BDT cut) is the dashed reference. PyROOT, env higgs-alp-ana.

  plots/cmp_mA<m>_signal.{pdf,png}   signal (test tree), unit area
  plots/cmp_mA<m>_dy.{pdf,png}       DY simulation, unit area (sculpting)
  plots/cmp_mA<m>_data.{pdf,png}     data sidebands only (115-135 GeV never filled), unit area
"""
import ast
import json
import os

import numpy as np
import ROOT

from compare_models import load, ERAS, BKG, inr, is_sb, TEST_FRAC

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "plots")
CONST = "/afs/cern.ch/work/p/pelai/HZgamma/higgsdna-hzg-run3/higgs_dna/metaconditions/corrections/constants.py"


def _read_lumi(path):
    for node in ast.parse(open(path).read()).body:
        if isinstance(node, ast.Assign) and any(getattr(t, "id", "") == "LUMI" for t in node.targets):
            return ast.literal_eval(node.value)
    raise RuntimeError("LUMI not found in %s" % path)


_L = _read_lumi(CONST)
LUMI_LABEL = "%.2f fb^{-1} (13.6 TeV)" % (_L["2022"] + _L["2023"] + _L["2024"])
STYLE = {"presel": (ROOT.kGray + 2, 2, "preselection (no BDT cut)"),
         "prod_wp": (ROOT.kBlack, 1, "production BDT, WP %.4f"),
         "eqdata": (ROOT.kAzure + 2, 1, "new BDT, same data count, %.4f"),
         "R1": (ROOT.kRed + 1, 1, "new BDT, R = 1 crossing, %.4f")}


def hist(name, x, w, edges, key, skip=None):
    h = ROOT.TH1D(name, "", len(edges) - 1, np.asarray(edges, float)); h.Sumw2()
    for xi, wi in zip(x, w):
        if skip and skip[0] < xi < skip[1]:
            continue
        h.Fill(xi, wi)
    if h.Integral() > 0:
        h.Scale(1.0 / h.Integral())
    col, ls, _ = STYLE[key]
    h.SetLineColor(col); h.SetLineStyle(ls); h.SetLineWidth(2)
    for ax in (h.GetXaxis(), h.GetYaxis()):
        ax.SetTitleSize(0.055); ax.SetLabelSize(0.050)
    h.GetYaxis().SetTitleOffset(1.35)
    return h


def draw(name, hs, ytitle, header):
    c = ROOT.TCanvas(name, "", 800, 700)
    c.SetLeftMargin(0.17); c.SetBottomMargin(0.14); c.SetTopMargin(0.09)
    ymax = max(h.GetMaximum() for h, _ in hs)
    for i, (h, _) in enumerate(hs):
        h.SetMaximum(1.55 * ymax); h.SetMinimum(0)
        h.GetXaxis().SetTitle("m_{ll#gamma#gamma} (GeV)"); h.GetYaxis().SetTitle(ytitle)
        h.Draw("hist" if i == 0 else "hist same")
    leg = ROOT.TLegend(0.40, 0.66, 0.94, 0.89); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.029)
    leg.SetHeader(header)
    for h, lab in hs:
        leg.AddEntry(h, lab, "l")
    leg.Draw()
    t = ROOT.TLatex(); t.SetNDC(); t.SetTextSize(0.045)
    t.DrawLatex(0.17, 0.92, "#bf{CMS} #it{Preliminary}")
    t.SetTextAlign(31); t.DrawLatex(0.95, 0.92, LUMI_LABEL)
    for ext in ("pdf", "png"):
        c.SaveAs(os.path.join(OUT, "%s.%s" % (name, ext)))


def main():
    os.makedirs(OUT, exist_ok=True)
    S = json.load(open(os.path.join(HERE, "compare_models.json")))
    dy = load([(s, e) for e in ERAS for s in BKG[e]], "inclusive")
    da = load([("Data", e) for e in ERAS], "inclusive")
    for mA in (1, 2, 3):
        sig = load([(f"mA_M{mA}", e) for e in ERAS], "test")
        sc = f"MVA_Score_mA_M{mA}"
        wps = [("prod", "prod_wp")] + [("new", k) for k in ("eqdata", "R1") if f"new_mA{mA}_{k}" in S]
        out = {"signal": [], "dy": [], "data": []}
        for kind, src, edges, skip in (("signal", sig, np.arange(95, 181, 1.0), None),
                                        ("dy", dy, np.arange(95, 181, 5.0), None),
                                        ("data", da, np.arange(95, 181, 5.0), (115, 135))):
            base = src["prod"]
            m = base["H_mass"]
            w = (base["factor"].astype(float) / TEST_FRAC) if kind == "signal" else \
                (base["factor"].astype(float) if kind == "dy" else np.ones_like(m))
            ok = inr(m, 95, 180) & np.isfinite(w)
            out[kind].append((hist(f"{kind}{mA}_pre", m[ok], w[ok], edges, "presel", skip), STYLE["presel"][2]))
            for tag, key in wps:
                c = S[f"{tag}_mA{mA}_{key}"]["cut"]
                p = ok & (src[tag][sc] > c)
                out[kind].append((hist(f"{kind}{mA}_{key}", m[p], w[p], edges, key, skip), STYLE[key][2] % c))
        draw(f"cmp_mA{mA}_signal", out["signal"], "Fraction / 1 GeV", f"m_{{a}} = {mA} GeV signal")
        draw(f"cmp_mA{mA}_dy", out["dy"], "Fraction / 5 GeV", f"DY simulation, m_{{a}} = {mA} GeV score")
        draw(f"cmp_mA{mA}_data", out["data"], "Fraction / 5 GeV", f"Data sidebands, m_{{a}} = {mA} GeV score")


if __name__ == "__main__":
    main()
