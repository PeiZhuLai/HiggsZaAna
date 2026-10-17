#!/usr/bin/env python3
"""Plots and R bootstrap for the sandbox BDT variants at their equal-data working points
(eval_variants.py must have run). PyROOT, env higgs-alp-ana.

  mllgg_signal.pdf   m_a = 1 signal m_llgg after the cut (absolute yield)
  mllgg_dy.pdf       DY simulation m_llgg after the cut, unit area (sculpting)
  mllgg_data.pdf     data m_llgg after the cut, 115-133 GeV not shown, unit area of the sidebands
  R_bootstrap.txt    R at the cut with its bootstrap spread (200 resamplings of the DY events)
"""
import ast
import json
import os

import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "plots")
CONST = "/afs/cern.ch/work/p/pelai/HZgamma/higgsdna-hzg-run3/higgs_dna/metaconditions/corrections/constants.py"
SHOW = [("production", ROOT.kBlack, "production (stored score > 0.963)"),
        ("win", ROOT.kRed + 1, "(1) signal window only"),
        ("nopt", ROOT.kGreen + 2, "(2) without p_{T}^{#gamma1}/m_{ll#gamma#gamma}"),
        ("win_mfree", ROOT.kAzure + 2, "(1) + inputs without m_{ll#gamma#gamma}")]


def _read_lumi(path):
    for node in ast.parse(open(path).read()).body:
        if isinstance(node, ast.Assign) and any(getattr(t, "id", "") == "LUMI" for t in node.targets):
            return ast.literal_eval(node.value)
    raise RuntimeError("LUMI not found in %s" % path)


_L = _read_lumi(CONST)
LUMI_LABEL = "%.2f fb^{-1} (13.6 TeV)" % (_L["2022"] + _L["2023"] + _L["2024"])


def labels():
    t = ROOT.TLatex(); t.SetNDC(); t.SetTextSize(0.045)
    t.DrawLatex(0.16, 0.92, "#bf{CMS} #it{Preliminary}")
    t.SetTextAlign(31); t.DrawLatex(0.95, 0.92, LUMI_LABEL)
    return t


def hist(name, x, w, nb, lo, hi, col, norm=None):
    h = ROOT.TH1D(name, "", nb, lo, hi); h.Sumw2()
    for xi, wi in zip(x, w):
        h.Fill(xi, wi)
    if norm == "area" and h.Integral() > 0:
        h.Scale(1.0 / h.Integral())
    h.SetLineColor(col); h.SetMarkerColor(col); h.SetLineWidth(2)
    for ax in (h.GetXaxis(), h.GetYaxis()):
        ax.SetTitleSize(0.055); ax.SetLabelSize(0.050)
    h.GetYaxis().SetTitleOffset(1.3)
    return h


def draw(name, hs, ytitle, header, opt="hist"):
    c = ROOT.TCanvas(name, "", 800, 700)
    c.SetLeftMargin(0.16); c.SetBottomMargin(0.14); c.SetTopMargin(0.09)
    ymax = max(h.GetMaximum() for h, _ in hs)
    for i, (h, _) in enumerate(hs):
        h.SetMaximum(1.45 * ymax); h.SetMinimum(0)
        h.GetXaxis().SetTitle("m_{ll#gamma#gamma} (GeV)"); h.GetYaxis().SetTitle(ytitle)
        h.Draw(opt if i == 0 else opt + " same")
    leg = ROOT.TLegend(0.36, 0.64, 0.94, 0.88); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.029)
    leg.SetHeader(header)
    for h, lab in hs:
        leg.AddEntry(h, lab, "l")
    leg.Draw(); t = labels()
    for ext in ("pdf", "png"):
        c.SaveAs(os.path.join(OUT, "%s.%s" % (name, ext)))


def main():
    os.makedirs(OUT, exist_ok=True)
    ev = json.load(open(os.path.join(HERE, "eval_variants.json")))
    m = np.load(os.path.join(HERE, "samples_meta.npz"))
    sets = {}
    for name, col, lab in SHOW:
        if name == "production":
            cut = ev["production_stored_score"]["production"]["cut"]
            sc = {"sig": m["sig_stored"], "bkg": m["bkg_stored"], "data": m["data_stored"]}
        else:
            cut = ev[name]["eqdata"]["cut"]
            sc = np.load(os.path.join(HERE, "scores_%s.npz" % name))
        sets[name] = (cut, sc, col, lab)

    hs_s, hs_b, hs_d = [], [], []
    rng = np.random.default_rng(7)
    lines = ["%-11s %7s %7s %7s  (200 bootstrap resamplings of the DY events, raw weights)" % ("variant", "cut", "R", "sd(R)")]
    bm, bw = m["bkg_m"], m["bkg_w"]
    sel = (bm > 95) & (bm < 180); pk = (bm > 120) & (bm < 130)
    for name, (cut, sc, col, lab) in sets.items():
        ps, pb, pd_ = sc["sig"] > cut, sc["bkg"] > cut, sc["data"] > cut
        hs_s.append((hist("s_" + name, m["sig_m"][ps], m["sig_w"][ps], 85, 95, 180, col), "%s, %.4f" % (lab, cut)))
        hs_b.append((hist("b_" + name, bm[pb], bw[pb], 17, 95, 180, col, "area"), lab))
        dm = m["data_m"][pd_]; dm = dm[(dm > 105) & (dm < 180) & ~((dm > 115) & (dm < 133))]
        hs_d.append((hist("d_" + name, dm, np.ones_like(dm), 15, 105, 180, col, "area"), "%s, N=%d" % (lab, len(dm))))
        idx = np.flatnonzero(sel)
        Rs = []
        for _ in range(200):
            j = rng.choice(idx, len(idx), replace=True)
            w, p, k = bw[j], pb[j], pk[j]
            frac = w[k].sum() / w.sum()
            Rs.append(w[p & k].sum() / w[p].sum() / frac)
        w, p, k = bw[sel], pb[sel], pk[sel]
        R = w[p & k].sum() / w[p].sum() / (w[k].sum() / w.sum())
        neff = w[p & k].sum() ** 2 / (w[p & k] ** 2).sum()
        lines.append("%-11s %7.4f %7.3f %7.3f   n_eff(120-130, pass) = %.0f" % (name, cut, R, np.std(Rs), neff))
    draw("mllgg_signal", hs_s, "Events / 1 GeV", "m_{a} = 1 GeV signal, same data sideband count")
    draw("mllgg_dy", hs_b, "Fraction / 5 GeV", "DY simulation")
    draw("mllgg_data", hs_d, "Fraction / 5 GeV", "Data, 115-133 GeV not shown", opt="hist")
    open(os.path.join(HERE, "R_bootstrap.txt"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
