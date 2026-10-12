#!/usr/bin/env python3
"""Signal m(llgg) before and after the sideband reweighting, after the BDT working point.

Answers the L3 question (Sec. 9.2): is the after/before ratio of the signal m(llgg)
flat, so that the reweighting can be treated as a normalization (lnN) uncertainty?
Same selection as studies/l3_review_20260910/an_update_sigrw_20261003/sig_mllgg_rw_shape.py:
the events entering the datacards (MVA-cut nominal trees), the five Run 3 eras summed,
electron and muon channels separately, both histograms normalized to unit area so that
only the shape is compared. "Before" = weight / weight_sideband_rwgt.

The ratio uncertainty is the statistical spread of the reweighting factor within the bin
(before and after are the same events, so the two histograms are fully correlated).

Input : /eos/home-p/pelai/HZa/root_MVAcut/sig/mA_M<ma>/output_<era>.root
Output: Plot/plots/signalReweight_bdtwp/sigRwgt_bdtwp_mllgg_mA_M<ma>_<ch>.pdf
        Plot/plots/signalReweight_bdtwp/sigRwgt_bdtwp_summary.txt

Usage : python plot_signal_reweight_mllgg_bdtwp.py [--ma 1,5,15,30]
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
MASSES = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30]
CHANNELS = {"ele": "ee", "mu": "#mu#mu"}
# 2.5 GeV bins over the full m(llgg) fit range; 115-130 GeV is an integer number of bins.
EDGES = np.arange(95.0, 180.0 + 1e-6, 2.5)


def run3_lumi() -> float:
    # The file cannot be imported (its ERA dict uses sets as keys): parse only LUMI.
    tree = ast.parse(Path(CONSTANTS).read_text())
    for node in tree.body:
        if isinstance(node, ast.Assign) and getattr(node.targets[0], "id", None) == "LUMI":
            lumi = ast.literal_eval(node.value)
            return sum(lumi[e] for e in ERAS)
    raise RuntimeError(f"LUMI not found in {CONSTANTS}")


def load(ma: int, ch: str):
    xs, wa, wb = [], [], []
    for era in ERAS:
        path = IN_BASE / f"mA_M{ma}" / f"output_{era}.root"
        tree = uproot.open(path)[f"DiphotonTree/ggh_125_Za_{ch}_13p6TeV_cat0"]
        a = tree.arrays(["CMS_hza_mass", "weight", "weight_sideband_rwgt"], library="np")
        xs.append(a["CMS_hza_mass"])
        wa.append(a["weight"])
        wb.append(a["weight"] / a["weight_sideband_rwgt"])
    return np.concatenate(xs), np.concatenate(wa), np.concatenate(wb)


def to_hist(name, x, w):
    h = ROOT.TH1D(name, "", len(EDGES) - 1, EDGES)
    h.Sumw2()
    h.FillN(len(x), np.ascontiguousarray(x, dtype=np.float64),
            np.ascontiguousarray(w, dtype=np.float64))
    h.SetDirectory(0)
    h.Scale(1.0 / h.Integral())
    return h


def ratio_hist(name, x, wa, wb, min_frac=1e-3):
    """After/before per bin, both normalized to unit area. With r = wa/wb the per-event
    reweighting factor, the bin ratio is the wb-weighted mean of r divided by the global
    one; its uncertainty is the standard error of that weighted mean. Bins with less than
    min_frac of the signal are left empty."""
    h = ROOT.TH1D(name, "", len(EDGES) - 1, EDGES)
    h.SetDirectory(0)
    norm = wb.sum() / wa.sum()
    idx = np.digitize(x, EDGES) - 1
    # events with zero weight carry no information (and would give 0/0)
    r = np.divide(wa, wb, out=np.ones_like(wa), where=wb != 0)
    for i in range(len(EDGES) - 1):
        sel = idx == i
        sw = wb[sel].sum()
        # bins holding less than 0.1% of the signal are not shown
        if sel.sum() < 2 or sw <= 0 or sw < min_frac * wb.sum():
            continue
        mean = np.sum(wb[sel] * r[sel]) / sw
        var = np.sum(wb[sel] ** 2 * (r[sel] - mean) ** 2) / sw ** 2
        h.SetBinContent(i + 1, mean * norm)
        h.SetBinError(i + 1, np.sqrt(var) * norm)
    return h


def width68(x, w):
    o = np.argsort(x)
    c = np.cumsum(w[o]) / w.sum()
    lo, hi = np.interp([0.16, 0.84], c, x[o])
    return hi - lo


def draw(ma, ch, hb, ha, hr, lumi, out_path: Path):
    canvas = ROOT.TCanvas(f"c_{ma}_{ch}", "", 800, 800)
    pad1 = ROOT.TPad("pad1", "", 0, 0.30, 1, 1)
    pad2 = ROOT.TPad("pad2", "", 0, 0.0, 1, 0.30)
    pad1.SetBottomMargin(0.02)
    pad1.SetLeftMargin(0.17)
    pad2.SetTopMargin(0.03)
    pad2.SetBottomMargin(0.40)
    pad2.SetLeftMargin(0.17)
    pad1.Draw()
    pad2.Draw()

    pad1.cd()
    for h, color in ((hb, ROOT.kBlack), (ha, ROOT.kRed + 1)):
        h.SetLineColor(color)
        h.SetLineWidth(LINE_WIDTH)
        h.SetMarkerColor(color)
    hb.GetYaxis().SetTitle("Fraction / 2.5 GeV")
    format_axes(hb, title_size=0.055 / 0.7, label_size=0.05 / 0.7)
    hb.GetXaxis().SetTitleSize(0)
    hb.GetXaxis().SetLabelSize(0)
    hb.GetYaxis().SetTitleOffset(1.2)
    hb.SetMaximum(1.9 * max(hb.GetMaximum(), ha.GetMaximum()))
    hb.SetMinimum(0.0)
    # zoom on the populated range (m_a = 1 GeV has a second, high-mass peak)
    xlo, xhi = (95.0, 170.0) if ma == 1 else (100.0, 150.0)
    hb.GetXaxis().SetRangeUser(xlo, xhi)
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
    txt.DrawLatex(0.21, 0.84, f"H#rightarrowZa#rightarrow{CHANNELS[ch]}#gamma#gamma, m_{{a}} = {ma} GeV")
    txt.DrawLatex(0.21, 0.76, "After the BDT selection")
    labels = cms_label(lumi=lumi, x_left=0.17, cms_text="CMS Simulation", text_size=0.045 / 0.7, y=0.94)

    pad2.cd()
    hr.SetLineColor(ROOT.kRed + 1)
    hr.SetMarkerColor(ROOT.kRed + 1)
    hr.SetMarkerStyle(20)
    hr.SetMarkerSize(0.8)
    hr.SetLineWidth(LINE_WIDTH)
    hr.GetYaxis().SetTitle("After / before")
    hr.GetXaxis().SetTitle("m_{#font[12]{ll}#gamma#gamma} [GeV]")
    format_axes(hr, title_size=0.055 / 0.3, label_size=0.05 / 0.3)
    hr.GetYaxis().SetTitleOffset(0.52)
    hr.GetXaxis().SetTitleOffset(0.95)
    hr.GetYaxis().SetNdivisions(505)
    hr.SetMinimum(0.7)
    hr.SetMaximum(1.3)
    hr.GetXaxis().SetRangeUser(xlo, xhi)
    hr.GetXaxis().SetNdivisions(505)
    hr.Draw("E1")
    keep = []
    line = ROOT.TLine(xlo, 1.0, xhi, 1.0)
    line.SetLineStyle(2)
    line.Draw("SAME")
    # Inside 115-135 GeV the m(llgg) reweight factors are not defined (signal region
    # excluded from the derivation), so any structure there comes from the other variables.
    for edge in (115.0, 135.0):
        m = ROOT.TLine(edge, 0.7, edge, 1.3)
        m.SetLineStyle(3)
        m.SetLineColor(ROOT.kGray + 2)
        m.Draw("SAME")
        keep.append(m)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    canvas.SaveAs(str(out_path))
    canvas.SaveAs(str(out_path.with_suffix(".png")))
    canvas.Close()
    del labels, leg, line, keep


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ma", default="all")
    parser.add_argument("--out", default=str(OUT_BASE))
    args = parser.parse_args()
    masses = MASSES if args.ma == "all" else [int(v) for v in args.ma.split(",")]
    out_base = Path(args.out)
    lumi = run3_lumi()
    print(f"lumi (2022-2024) = {lumi:.2f} fb-1")

    lines = ["ma ch   yield_ratio mean_shift[GeV] width68_rel  max|ratio-1| 115-130 / 110-140 "
             "(5 GeV bins, shape only)"]
    for ma in masses:
        for ch in CHANNELS:
            x, wa, wb = load(ma, ch)
            hb = to_hist(f"hb_{ma}_{ch}", x, wb)
            ha = to_hist(f"ha_{ma}_{ch}", x, wa)
            hr = ratio_hist(f"hr_{ma}_{ch}", x, wa, wb)
            # cross-check with the 5 GeV numbers quoted in the TWiki answer
            e5 = np.arange(110, 141, 5)
            a5, _ = np.histogram(x, e5, weights=wa)
            b5, _ = np.histogram(x, e5, weights=wb)
            r5 = (a5 / a5.sum()) / (b5 / b5.sum())
            core = np.abs(r5[1:4] - 1).max()
            outer = np.abs(r5 - 1).max()
            lines.append("%2d %-4s %9.3f %13.3f %12.3f   %.3f / %.3f" % (
                ma, ch, wa.sum() / wb.sum(), np.average(x, weights=wa) - np.average(x, weights=wb),
                width68(x, wa) / width68(x, wb) - 1, core, outer))
            print(lines[-1])
            draw(ma, ch, hb, ha, hr, lumi, out_base / f"sigRwgt_bdtwp_mllgg_mA_M{ma}_{ch}.pdf")
    out_base.mkdir(parents=True, exist_ok=True)
    (out_base / "sigRwgt_bdtwp_summary.txt").write_text("\n".join(lines) + "\n")
    print("wrote", out_base / "sigRwgt_bdtwp_summary.txt")


if __name__ == "__main__":
    main()
