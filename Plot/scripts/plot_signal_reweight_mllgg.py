#!/usr/bin/env python3
"""Signal m(llgg) before and after the sideband shape reweighting.

Answers the L3 question (Section 9.2): show the m(llgg) distribution of the
simulated signal with and without the BDT-training sideband reweighting, and
show whether the after/before ratio is flat -- which is what justifies quoting
the reweighting uncertainty as a flat lnN on the yield rather than as a shape.

Input : /eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix/mA_M<ma>/<era>.root
        (tree "inclusive"; branch "factor" is the training weight)
        HZaMVA/reweights/sideband_run3_iterative.json
Output: Plot/plots/signalReweight/<era>/sigRwgt_mllgg_mA_M<ma>_<era>.pdf
        Plot/plots/signalReweight/sigRwgt_summary.txt

Usage : python plot_signal_reweight_mllgg.py [--ma 1,5,15,30] [--year all]
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import uproot

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
sys.path.append("/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts")

import ROOT  # noqa: E402
from _cms_style import (  # noqa: E402
    LINE_WIDTH,
    cms_label,
    format_axes,
    make_legend,
)
from sideband_reweight import load_sideband_reweighter  # noqa: E402

_KEEP = []

ROOT_BASE = Path("/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix")
OUT_BASE = Path("/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/signalReweight")
ERAS = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]
MASSES = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30]

# m(llgg) histogram: full analysis range, 2 GeV bins.
MLLGG_MIN, MLLGG_MAX, MLLGG_NBINS = 95.0, 180.0, 34

# Columns the reweighter needs that are stored under these names in the ROOT files.
NEEDED = [
    "var_dR_Za", "var_dR_g1g2", "var_dR_g1Z",
    "pho1IetaIeta55", "pho2IetaIeta55",
    # ECAL isolation is stored under its ALP_* name; sideband_reweight.VAR_ALIASES
    # maps pho1ECALIso/pho2ECALIso onto these.
    "ALP_lead_photon_ecalPFClusterIso", "ALP_sublead_photon_ecalPFClusterIso",
    "ALP_calculatedPhotonIso",
    "pho1R9", "pho2R9",
    "var_PtaOverMh", "pho1Pt_oHm", "pho2Pt_oHm", "H_pt_oHm",
    "Z_m", "H_m", "ALP_m", "factor",
]


def load_frame(path: Path, ma: int) -> pd.DataFrame | None:
    if not path.exists():
        return None
    with uproot.open(f"{path}:inclusive") as tree:
        have = [b for b in NEEDED if b in tree.keys()]
        missing = [b for b in NEEDED if b not in tree.keys()]
        if missing:
            print(f"  [warn] {path.name}: missing {missing}")
        arrays = tree.arrays(have, library="np")
    frame = pd.DataFrame({k: np.asarray(v, dtype=float) for k, v in arrays.items()})
    # For signal the mass-hypothesis variable uses the TRUE generated ALP mass,
    # exactly as in run3_Za_BDT.py.
    frame["param"] = (frame["ALP_m"] - float(ma)) / frame["H_m"]
    return frame


def fill_hist(name, values, weights):
    hist = ROOT.TH1D(name, "", MLLGG_NBINS, MLLGG_MIN, MLLGG_MAX)
    hist.Sumw2()
    values = np.ascontiguousarray(np.asarray(values, dtype=np.float64))
    weights = np.ascontiguousarray(np.asarray(weights, dtype=np.float64))
    if len(values):
        hist.FillN(len(values), values, weights)
    hist.SetDirectory(0)
    return hist


def draw(ma, era, h_before, h_after, out_path: Path, lumi_key: str):
    canvas = ROOT.TCanvas(f"c_{ma}_{era}", "", 800, 800)
    pad1 = ROOT.TPad("pad1", "", 0, 0.30, 1, 1)
    pad2 = ROOT.TPad("pad2", "", 0, 0.0, 1, 0.30)
    pad1.SetBottomMargin(0.02)
    pad1.SetLeftMargin(0.14)
    pad2.SetTopMargin(0.03)
    pad2.SetBottomMargin(0.35)
    pad2.SetLeftMargin(0.14)
    pad1.Draw()
    pad2.Draw()

    pad1.cd()
    for hist, color in ((h_before, ROOT.kBlack), (h_after, ROOT.kRed + 1)):
        hist.SetLineColor(color)
        hist.SetLineWidth(LINE_WIDTH)
        hist.SetMarkerColor(color)
    h_before.GetYaxis().SetTitle("A.U.")
    format_axes(h_before, title_size=0.055 / 0.7, label_size=0.05 / 0.7)
    h_before.GetXaxis().SetTitleSize(0)
    h_before.GetXaxis().SetLabelSize(0)
    h_before.GetYaxis().SetTitleOffset(1.05)
    h_before.SetMaximum(1.45 * max(h_before.GetMaximum(), h_after.GetMaximum()))
    h_before.Draw("HIST")
    h_after.Draw("HIST SAME")

    leg = make_legend(0.55, 0.66, 0.92, 0.84, text_size=0.045 / 0.7)
    leg.AddEntry(h_before, "Before reweighting", "l")
    leg.AddEntry(h_after, "After reweighting", "l")
    leg.Draw()

    txt = ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextFont(42)
    txt.SetTextSize(0.045 / 0.7)
    txt.DrawLatex(0.18, 0.84, f"H#rightarrowZa, m_{{a}} = {ma} GeV")
    keep_labels = cms_label(year=lumi_key, cms_text="CMS Simulation",
                            text_size=0.045 / 0.7, y=0.94)

    pad2.cd()
    ratio = h_after.Clone(f"ratio_{ma}_{era}")
    ratio.Divide(h_before)
    ratio.SetLineColor(ROOT.kRed + 1)
    ratio.SetLineWidth(LINE_WIDTH)
    ratio.GetYaxis().SetTitle("After / before")
    ratio.GetXaxis().SetTitle("m_{#font[12]{ll}#gamma#gamma} [GeV]")
    format_axes(ratio, title_size=0.055 / 0.3, label_size=0.05 / 0.3)
    ratio.GetYaxis().SetTitleOffset(0.45)
    ratio.GetXaxis().SetTitleOffset(1.0)
    ratio.GetYaxis().SetNdivisions(505)
    ratio.SetMinimum(0.5)
    ratio.SetMaximum(1.5)
    ratio.Draw("HIST")
    # Draw the unity line inside the frame only (a TLine spanning the full user
    # range leaks into the pad margin).
    line = ROOT.TLine(ratio.GetXaxis().GetXmin(), 1.0,
                      ratio.GetXaxis().GetXmax(), 1.0)
    line.SetLineStyle(2)
    line.Draw("SAME")
    # Mark the sideband boundaries: inside 115--135 GeV the m(llgg) reweight
    # factors are not defined, so any structure there comes from the other
    # training variables.
    for edge in (115.0, 135.0):
        marker = ROOT.TLine(edge, 0.5, edge, 1.5)
        marker.SetLineStyle(3)
        marker.SetLineColor(ROOT.kGray + 2)
        marker.Draw("SAME")
        _KEEP.append(marker)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    canvas.SaveAs(str(out_path))
    canvas.Close()
    del keep_labels, leg, line, ratio


def spread(h_before, h_after, x_lo=None, x_hi=None):
    """Yield ratio and shape non-flatness of after/before in [x_lo, x_hi].

    The shape numbers are computed after renormalizing the two histograms to
    the same integral, so they measure only the departure from a flat lnN.
    """
    n_before = h_before.Integral()
    n_after = h_after.Integral()
    devs = []
    for i in range(1, h_before.GetNbinsX() + 1):
        center = h_before.GetBinCenter(i)
        if x_lo is not None and center < x_lo:
            continue
        if x_hi is not None and center > x_hi:
            continue
        b = h_before.GetBinContent(i)
        if b <= 0 or b < 0.01 * h_before.GetMaximum():
            continue
        a = h_after.GetBinContent(i)
        devs.append(a / b * (n_before / n_after) - 1.0)
    devs = np.asarray(devs) if devs else np.zeros(1)
    return (n_after / n_before if n_before else float("nan"),
            float(np.max(np.abs(devs))), float(np.sqrt(np.mean(devs ** 2))))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ma", default="all")
    parser.add_argument("--year", default="all")
    parser.add_argument("--out", default=str(OUT_BASE))
    args = parser.parse_args()

    masses = MASSES if args.ma == "all" else [int(x) for x in args.ma.split(",")]
    eras = ERAS if args.year == "all" else args.year.split(",")

    reweighter = load_sideband_reweighter(required=True)
    print("reweight JSON:", reweighter.source_path)

    out_base = Path(args.out)
    lines = ["ma  era            yield_ratio  full range: max|dev| rms(dev)   "
             "115-135 GeV: max|dev| rms(dev)"]
    for ma in masses:
        for era in eras:
            path = ROOT_BASE / f"mA_M{ma}" / f"{era}.root"
            frame = load_frame(path, ma)
            if frame is None or frame.empty:
                print(f"  [skip] {path}")
                continue
            rwgt = reweighter.weights_for_dataframe(frame)
            w0 = frame["factor"].to_numpy(dtype=float)
            mllgg = frame["H_m"].to_numpy(dtype=float)
            h_before = fill_hist(f"hb_{ma}_{era}", mllgg, w0)
            h_after = fill_hist(f"ha_{ma}_{era}", mllgg, w0 * rwgt)
            y_ratio, dev_max, dev_rms = spread(h_before, h_after)
            _, sr_max, sr_rms = spread(h_before, h_after, 115.0, 135.0)
            lines.append(f"{ma:>2}  {era:<14} {y_ratio:10.4f}  "
                         f"{dev_max:9.4f} {dev_rms:9.4f}   {sr_max:9.4f} {sr_rms:9.4f}")
            print(lines[-1])
            draw(ma, era, h_before, h_after,
                 out_base / era / f"sigRwgt_mllgg_mA_M{ma}_{era}.pdf", era)

    out_base.mkdir(parents=True, exist_ok=True)
    (out_base / "sigRwgt_summary.txt").write_text("\n".join(lines) + "\n")
    print("wrote", out_base / "sigRwgt_summary.txt")


if __name__ == "__main__":
    main()
