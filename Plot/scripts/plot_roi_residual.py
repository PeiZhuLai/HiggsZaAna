"""
PyROOT version of the m_Gamma ROI residual plot (AN-20-142 Fig 9 analogue).

Reads the merged_nominal parquet produced by HiggsDNA on the MLNanoAOD signal
sample, and draws two TH1F overlays per mass point on one canvas:
  Left  : m_Gamma^regressed distribution + dashed line at gen m_a
  Right : residual = MergedML_residual_mass = m_Gamma - m_a^gen
          (the "ROI" of AN §5.4)

Input parquet typically lives under
    /eos/.../parquet_merged_DNA_tmp/Sig_MC_MLNANO{,_all}/<sample>_2024/merged_nominal.parquet

Output:
    /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/roi_residual/<sample>.{png,pdf,root}
"""

from __future__ import annotations

import argparse
import os
import sys

import ROOT

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from _cms_style import (
    apply_line_width, format_axes, make_legend, LINE_WIDTH,
)

# 2024 integrated luminosity, matching plot_rechits_2d.py's DEFAULT_LUMI.
DEFAULT_LUMI = "109.82 fb^{-1}"


def _format_axes(h: ROOT.TH1F) -> None:
    format_axes(h, y_title_offset=1.60)


def _draw_header(left_x: float = 0.15, right_x: float = 0.95, y: float = 0.92,
                 lumi_label: str = DEFAULT_LUMI):
    """CMS Simulation (top-left) + '<lumi> (13.6 TeV)' (top-right), drawn above
    the frame in the same style as plot_rechits_2d.py. Returns the TLatex
    objects so the caller keeps them alive until the canvas is saved.
    """
    txt_left = ROOT.TLatex()
    txt_left.SetNDC(True)
    txt_left.SetTextFont(62)
    txt_left.SetTextSize(0.045)
    txt_left.SetTextAlign(11)
    txt_left.DrawLatex(left_x, y, "CMS Simulation")

    lumi_text = "%s (13.6 TeV)" % lumi_label if lumi_label else "(13.6 TeV)"
    txt_right = ROOT.TLatex()
    txt_right.SetNDC(True)
    txt_right.SetTextFont(42)
    txt_right.SetTextSize(0.045)
    txt_right.SetTextAlign(31)
    txt_right.DrawLatex(right_x, y, lumi_text)
    return txt_left, txt_right


def _load(parquet_path: str):
    import pandas as pd
    return pd.read_parquet(parquet_path)


def make_plot(parquet_path: str, gen_ma: float, out_stub: str) -> None:
    df = _load(parquet_path)
    res = df["MergedML_residual_mass"].to_numpy()
    dr = df["MergedML_dR_to_genA"].to_numpy()
    mass = df["MergedML_mass"].to_numpy()

    import numpy as np
    ok = res > -100.0
    truth = ok & (dr < 0.05)
    sel = truth & df["pass_allcuts_merged_ML"].astype(bool).to_numpy()

    # left: m_Gamma
    h_m_all = ROOT.TH1F("h_m_all",   "all (#DeltaR<3);m_{#Gamma}^{regressed} [GeV];Events / 25 MeV", 48, 0.0, 1.2)
    h_m_tm  = ROOT.TH1F("h_m_tm",    "truth-match #DeltaR<0.05",                                       48, 0.0, 1.2)

    # right: residual
    h_r_all = ROOT.TH1F("h_r_all",   "all;m_{#Gamma}^{regressed} - m_{a}^{gen} [GeV];Events / 25 MeV", 60, -0.5, 1.0)
    h_r_tm  = ROOT.TH1F("h_r_tm",    "truth-match #DeltaR<0.05",                                       60, -0.5, 1.0)
    h_r_pass = ROOT.TH1F("h_r_pass", "+ pass_allcuts_merged_ML",                                       60, -0.5, 1.0)

    for v in mass[ok & (dr < 3.0)]:
        h_m_all.Fill(v)
    for v in mass[truth]:
        h_m_tm.Fill(v)
    for v in res[ok & (dr < 3.0)]:
        h_r_all.Fill(v)
    for v in res[truth]:
        h_r_tm.Fill(v)
    for v in res[sel]:
        h_r_pass.Fill(v)

    for h, color in [
        (h_m_all,  ROOT.kAzure + 2),
        (h_m_tm,   ROOT.kOrange + 1),
        (h_r_all,  ROOT.kAzure + 2),
        (h_r_tm,   ROOT.kOrange + 1),
        (h_r_pass, ROOT.kGreen + 2),
    ]:
        h.SetLineColor(color)
        apply_line_width(h, LINE_WIDTH)
        _format_axes(h)

    # The blue "all" curves are drawn dashed.
    h_m_all.SetLineStyle(2)
    h_r_all.SetLineStyle(2)

    os.makedirs(os.path.dirname(out_stub), exist_ok=True)

    # =================== Canvas 1: m_Gamma ===================
    c1 = ROOT.TCanvas("c_mGamma", "m_Gamma regressed", 800, 700)
    c1.SetLeftMargin(0.19)
    c1.SetRightMargin(0.05)
    c1.SetBottomMargin(0.16)
    c1.SetTopMargin(0.10)
    hmax = max(h_m_all.GetMaximum(), h_m_tm.GetMaximum())
    h_m_all.SetMaximum(hmax * 1.6)   # y-axis range = 1.6x the histogram peak
    h_m_all.Draw("HIST")
    h_m_tm.Draw("HIST SAME")

    # vertical line only as tall as the histogram peak (not the full y-axis)
    line_gen = ROOT.TLine(gen_ma, 0.0, gen_ma, hmax)
    line_gen.SetLineColor(ROOT.kRed + 1)
    line_gen.SetLineStyle(2)
    line_gen.SetLineWidth(LINE_WIDTH)
    line_gen.Draw()

    leg1 = make_legend(0.32, 0.73, 0.80, 0.87, text_size=0.033)
    leg1.AddEntry(h_m_all, f"all (#DeltaR<3) N={int(h_m_all.GetEntries())}", "l")
    leg1.AddEntry(h_m_tm,  f"truth-match #DeltaR<0.05 N={int(h_m_tm.GetEntries())}", "l")
    leg1.AddEntry(line_gen, f"gen m_{{a}}={gen_ma:.1f} GeV", "l")
    leg1.Draw()

    cms_lab1 = _draw_header()

    out1 = out_stub + "_mGamma"
    for ext in ("png", "pdf"):
        c1.SaveAs(f"{out1}.{ext}")
        print(f"[plot_roi_residual] saved {out1}.{ext}")

    # =================== Canvas 2: residual ===================
    c2 = ROOT.TCanvas("c_residual", "ROI residual", 800, 700)
    c2.SetLeftMargin(0.19)
    c2.SetRightMargin(0.05)
    c2.SetBottomMargin(0.16)
    c2.SetTopMargin(0.10)
    hmax2 = max(h_r_all.GetMaximum(), h_r_tm.GetMaximum(), h_r_pass.GetMaximum())
    h_r_all.SetMaximum(hmax2 * 1.6)   # y-axis range = 1.6x the histogram peak
    h_r_all.Draw("HIST")
    h_r_tm.Draw("HIST SAME")
    h_r_pass.Draw("HIST SAME")

    # vertical line only as tall as the histogram peak (not the full y-axis)
    line_zero = ROOT.TLine(0.0, 0.0, 0.0, hmax2)
    line_zero.SetLineColor(ROOT.kRed + 1)
    line_zero.SetLineStyle(2)
    line_zero.SetLineWidth(LINE_WIDTH)
    line_zero.Draw()

    leg2 = make_legend(0.32, 0.69, 0.80, 0.87, text_size=0.033)
    leg2.AddEntry(h_r_all,  f"all (#DeltaR<3) N={int(h_r_all.GetEntries())}", "l")
    leg2.AddEntry(h_r_tm,   f"truth-match #DeltaR<0.05 N={int(h_r_tm.GetEntries())}", "l")
    leg2.AddEntry(h_r_pass, f"+ pass_allcuts_merged_ML N={int(h_r_pass.GetEntries())}", "l")
    leg2.AddEntry(line_zero, "ideal (residual=0)", "l")
    leg2.Draw()

    cms_lab2 = _draw_header()

    out2 = out_stub + "_residual"
    for ext in ("png", "pdf"):
        c2.SaveAs(f"{out2}.{ext}")
        print(f"[plot_roi_residual] saved {out2}.{ext}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("parquet", help="merged_nominal.parquet path")
    ap.add_argument("--gen-ma", type=float, required=True,
                    help="generated m_a for this sample (GeV)")
    ap.add_argument(
        "--output-stub", "-o",
        default=None,
        help="output stub (no ext). default: Plot/plots/roi_residual/mA_<ma>")
    args = ap.parse_args()

    if args.output_stub is None:
        tag = f"mA_{args.gen_ma:.3f}".replace(".", "p")
        args.output_stub = f"/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/roi_residual/{tag}"
    make_plot(args.parquet, args.gen_ma, args.output_stub)


if __name__ == "__main__":
    main()
