#!/usr/bin/env python3
"""Background m_llgg shape in the signal-region BDT slice vs the adjacent slice of equal width.

Reviewer follow-up (NPS-25-014 TWiki, v1 Fig. 23 / current fig:bdt_mass_shapes): instead of
0.1-wide BDT-score slices, compare for every mass hypothesis m_a with adopted working point c
  SR slice       : c      <= score <= 1
  adjacent slice : 2c - 1 <= score <  c        (same width, 1 - c)
of the combined simulated DY background (DY+gamma and DY+jets), normalized to unit area, with
the SR/adjacent ratio. The score is the analysis's stored MVA_Score_mA_M{m} (the parametric BDT
evaluated at m); the cut c follows apply_bdt_data.py: anchors from MVAcut_points_run3.json, and
interpolated masses take the NEAREST anchor's cut (find_nearest_cut; not the linear
interpolation that plot_bkgmcSculptingCheck._complete_mva_cuts uses).

Samples/weights are the ones of the main (Nominal) suite of plot_bkgmcSculptingCheck.py:
Analyzer_Config('run3').sample_loc (run3_bdt_scored_fsrfix), DY+jets 2024 = flavor-split samples
(Plot_Helper._run3_sources_for_sample), branch H_m, raw 'weight' (negative weights kept).

Per mass it also reports (stdout + JSON):
  chi2/ndf and p-value of the two normalized 5 GeV shapes in 95-180 GeV (bin variances from sumw2),
  weighted unbinned KS distance and p-value (n_eff), and the peak-window fraction
  f = w(120-130)/w(95-180) in each slice with its uncertainty and their ratio.

Outputs (new subdirectory; nothing existing is overwritten):
  Plot/plots/bkgmcScupltingCheck/srAdjacentSlices/bkg_mass_shape_srAdjacent_mA{NN}.{pdf,png}
  Plot/plots/bkgmcScupltingCheck/srAdjacentSlices/srAdjacent_peakFraction_summary.{pdf,png}
  --json <path>  per-mass numbers

Env: higgs-alp-ana (ROOT 6.24 + uproot).
"""
import argparse
import json
import math
import os
import sys
from array import array
from pathlib import Path

import numpy as np
import uproot
import ROOT

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
sys.path.insert(0, str(SCRIPT_DIR.parent / "lib"))

import plot_bkgmcSculptingCheck as SC            # reuse constants / cut parsing / anchors
from Analyzer_Configs import Analyzer_Config
from Plot_Helper import _run3_sources_for_sample, SaveCanvPic

OUT_DEFAULT = SC.LOCAL_OUTPUT_DIR / "srAdjacentSlices"
PEAK = (120.0, 130.0)
XMIN, XMAX = SC.H_M_XMIN, SC.H_M_XMAX
BINW = SC.BDT_SHAPE_MASS_BIN_WIDTH                 # 5 GeV, as fig:bdt_mass_shapes
LUMI_TEXT = SC.BDT_SHAPE_LUMI_TEXT                 # "172.13 fb^{-1} (13.6 TeV)"
COL_SR, COL_ADJ = ROOT.kRed + 1, ROOT.kBlack       # new/SR red, reference/adjacent black

# house style (hzg-plot-style): numbers are single-pad ROOT values; divide by pad height
TITLE, LABEL, CMS_SIZE, NOTE, LEG = 0.055, 0.050, 0.038, 0.028, 0.032
CANVAS_W, CANVAS_H, SPLIT = 900, 950, 0.30
Y_TITLE_OFFSET, RATIO_Y_TITLE_OFFSET, X_TITLE_OFFSET = 1.10, 0.68, 1.05


def sz(house, pad):
    return house / float(pad.GetAbsHNDC())


def left_margin_for(w, h):
    return (2.1 * Y_TITLE_OFFSET + 1.3) * 0.055 * h / w + 0.012


def nearest_cuts(json_path):
    anchors = SC._parse_mva_cuts(str(json_path))
    cuts = {}
    for m in SC.TARGET_MASSES:
        a = m if m in anchors else min(sorted(anchors), key=lambda k: abs(k - m))
        cuts[m] = (a, float(anchors[a]))
    return cuts


def load_background(masses):
    cfg = Analyzer_Config("inclusive", "run3", 0, True)
    paths = []
    for sample in cfg.bkg_names:
        for directory, year in _run3_sources_for_sample(sample, cfg):
            p = os.path.join(cfg.sample_loc, directory, f"{year}.root")
            if not os.path.exists(p):
                raise FileNotFoundError(p)
            paths.append(p)
    brs = ["H_m", "weight"] + [f"MVA_Score_mA_M{m}" for m in masses]
    parts = {b: [] for b in brs}
    for p in paths:
        a = uproot.open(p)["inclusive"].arrays(brs, library="np")
        for b in brs:
            parts[b].append(a[b])
    d = {b: np.concatenate(v) for b, v in parts.items()}
    sel = (d["H_m"] >= XMIN) & (d["H_m"] <= XMAX) & np.isfinite(d["weight"])
    print(f"[input] {len(paths)} files from {cfg.sample_loc}; {sel.sum()} events in {XMIN:.0f}-{XMAX:.0f} GeV")
    for p in paths:
        print("   ", p)
    return {b: v[sel] for b, v in d.items()}


def make_hist(name, m, w):
    nb = int(round((XMAX - XMIN) / BINW))
    h = ROOT.TH1D(name, "", nb, XMIN, XMAX)
    h.Sumw2(); h.SetDirectory(0)
    if len(m):
        h.FillN(len(m), array("d", m), array("d", w))
    return h


def peak_fraction(m, w):
    W = w.sum()
    pk = (m > PEAK[0]) & (m < PEAK[1])
    f = w[pk].sum() / W
    var = np.sum(w ** 2 * (pk.astype(float) - f) ** 2) / W ** 2
    return float(f), float(math.sqrt(var))


def chi2_shapes(h1, h2):
    W1, W2 = h1.Integral(), h2.Integral()
    chi2, nb = 0.0, 0
    for i in range(1, h1.GetNbinsX() + 1):
        p1, p2 = h1.GetBinContent(i) / W1, h2.GetBinContent(i) / W2
        v = (h1.GetBinError(i) / W1) ** 2 + (h2.GetBinError(i) / W2) ** 2
        if v > 0:
            chi2 += (p1 - p2) ** 2 / v; nb += 1
    ndf = max(nb - 1, 1)
    return chi2, ndf, ROOT.TMath.Prob(chi2, ndf)


def weighted_ks(m1, w1, m2, w2):
    grid = np.sort(np.concatenate([m1, m2]))
    def ecdf(m, w):
        o = np.argsort(m); cw = np.cumsum(w[o]) / w.sum()
        idx = np.searchsorted(m[o], grid, side="right")
        return np.where(idx > 0, cw[np.maximum(idx - 1, 0)], 0.0)
    D = float(np.max(np.abs(ecdf(m1, w1) - ecdf(m2, w2))))
    n1 = w1.sum() ** 2 / np.sum(w1 ** 2); n2 = w2.sum() ** 2 / np.sum(w2 ** 2)
    ne = n1 * n2 / (n1 + n2)
    return D, float(ROOT.TMath.KolmogorovProb(D * math.sqrt(ne)))


def neff(w):
    return float(w.sum() ** 2 / np.sum(w ** 2)) if len(w) else 0.0


keep = []


def draw_mass(mass, anchor, cut, h_sr, h_adj, res, outdir):
    lo = 2 * cut - 1
    hs = h_sr.Clone(f"hs_{mass}"); ha = h_adj.Clone(f"ha_{mass}")
    for h in (hs, ha):
        h.SetDirectory(0); h.Scale(1.0 / h.Integral("width"))
    ratio = hs.Clone(f"r_{mass}"); ratio.SetDirectory(0); ratio.Divide(ha)

    c = ROOT.TCanvas(f"c_{mass}", "", CANVAS_W, CANVAS_H)
    lm = left_margin_for(CANVAS_W, CANVAS_H)
    up = ROOT.TPad(f"up_{mass}", "", 0, SPLIT, 1, 1)
    up.SetLeftMargin(lm); up.SetRightMargin(0.05)
    up.SetTopMargin(0.09 / (1 - SPLIT)); up.SetBottomMargin(0.02)
    up.SetTickx(1); up.SetTicky(1); up.Draw()
    dn = ROOT.TPad(f"dn_{mass}", "", 0, 0, 1, SPLIT)
    dn.SetLeftMargin(lm); dn.SetRightMargin(0.05); dn.SetTopMargin(0.03); dn.SetBottomMargin(0.42)
    dn.SetTickx(1); dn.SetTicky(1); dn.SetGridy(); dn.Draw()

    up.cd()
    ymax = max(hs.GetMaximum() + hs.GetBinError(hs.GetMaximumBin()),
               ha.GetMaximum() + ha.GetBinError(ha.GetMaximumBin()))
    frame = ROOT.TH1D(f"fr_{mass}", "", 1, XMIN, XMAX); frame.SetDirectory(0)
    HEAD = 2.3   # headroom: 3 text lines + 2-entry legend above the curves
    frame.SetMaximum(HEAD * ymax); frame.SetMinimum(-0.02 * ymax)
    fx, fy = frame.GetXaxis(), frame.GetYaxis()
    fx.SetLabelSize(0); fx.SetTitleSize(0)
    fy.SetTitle(f"A.U. / {BINW:.0f} GeV"); fy.SetTitleSize(sz(TITLE, up)); fy.SetLabelSize(sz(LABEL, up))
    fy.SetTitleOffset(Y_TITLE_OFFSET); fy.SetNdivisions(505)
    frame.Draw("AXIS")
    for h, col, st in ((ha, COL_ADJ, 1), (hs, COL_SR, 1)):
        h.SetLineColor(col); h.SetLineStyle(st); h.SetLineWidth(3); h.SetMarkerSize(0); h.SetFillStyle(0)
        h.Draw("HIST SAME"); h.Draw("E1 X0 SAME")
    lines = []
    for x in PEAK:
        ln = ROOT.TLine(x, frame.GetMinimum(), x, 1.1 * ymax)
        ln.SetLineColor(ROOT.kGray + 2); ln.SetLineStyle(2); ln.SetLineWidth(2); ln.Draw(); lines.append(ln)

    top = 1 - up.GetTopMargin()
    t = ROOT.TLatex(); t.SetNDC(); t.SetTextFont(42)
    t.SetTextSize(sz(CMS_SIZE, up)); t.SetTextAlign(11)
    t.DrawLatex(lm + 0.005, top + 0.020 / up.GetAbsHNDC(), "#bf{CMS} #it{Preliminary}")
    t.SetTextAlign(31); t.DrawLatex(0.95, top + 0.020 / up.GetAbsHNDC(), LUMI_TEXT)
    t.SetTextAlign(13); t.SetTextSize(sz(0.034, up))
    ylab = top - 0.035 / up.GetAbsHNDC()
    t.DrawLatex(lm + 0.03, ylab, f"m_{{a}} = {mass} GeV,  BDT cut c = {cut:.3f}")
    t.SetTextSize(sz(NOTE, up))
    dy = 0.050 / up.GetAbsHNDC()
    t.DrawLatex(lm + 0.03, ylab - 1.1 * dy,
                f"#chi^{{2}}/ndf = {res['chi2']:.1f}/{res['ndf']} (p = {res['chi2_p']:.2f}),  KS p = {res['ks_p']:.2f}")
    t.DrawLatex(lm + 0.03, ylab - 1.9 * dy,
                f"f_{{120-130}}:  SR {res['f_sr']:.3f} #pm {res['f_sr_err']:.3f},  adj. {res['f_adj']:.3f} #pm {res['f_adj_err']:.3f}")

    nl = 2.7
    if anchor != mass:
        t.DrawLatex(lm + 0.03, ylab - 2.7 * dy, f"(c of the nearest simulated m_{{a}} = {anchor} GeV)")
        nl = 3.5
    ytop_leg = ylab - nl * dy
    leg = ROOT.TLegend(0.52, ytop_leg - 2 * 0.060 / up.GetAbsHNDC(), 0.93, ytop_leg)
    leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextFont(42); leg.SetTextSize(sz(LEG, up))
    leg.AddEntry(hs, f"SR [{cut:.3f}, 1]", "l")
    leg.AddEntry(ha, f"Adjacent [{lo:.3f}, {cut:.3f})", "l")
    leg.Draw()

    dn.cd()
    rmax = 0.0
    for i in range(1, ratio.GetNbinsX() + 1):
        if ratio.GetBinContent(i) != 0:
            rmax = max(rmax, ratio.GetBinContent(i) + ratio.GetBinError(i))
    rhi = min(max(2.0, 1.1 * rmax), 4.0)
    rf = ROOT.TH1D(f"rf_{mass}", "", 1, XMIN, XMAX); rf.SetDirectory(0)
    rf.SetMinimum(0.0); rf.SetMaximum(rhi)
    rx, ry = rf.GetXaxis(), rf.GetYaxis()
    rx.SetTitle("m_{ll#gamma#gamma} [GeV]"); rx.SetTitleSize(sz(TITLE, dn)); rx.SetLabelSize(sz(LABEL, dn))
    rx.SetTitleOffset(X_TITLE_OFFSET); rx.SetLabelOffset(0.010); rx.SetTickLength(0.03 / SPLIT)
    ry.SetTitle("SR / adj."); ry.SetTitleSize(sz(0.038, dn)); ry.SetLabelSize(sz(0.036, dn))
    ry.SetTitleOffset(RATIO_Y_TITLE_OFFSET * 0.9); ry.SetNdivisions(505); ry.CenterTitle(True)
    rf.Draw("AXIS")
    one = ROOT.TLine(XMIN, 1.0, XMAX, 1.0); one.SetLineStyle(2); one.SetLineColor(ROOT.kGray + 2); one.Draw()
    ratio.SetLineColor(COL_SR); ratio.SetMarkerColor(COL_SR); ratio.SetMarkerStyle(20); ratio.SetMarkerSize(1.2)
    ratio.SetLineWidth(2); ratio.Draw("E1 SAME")
    keep.extend([hs, ha, ratio, frame, rf, lines, one, leg, up, dn])

    name = f"bkg_mass_shape_srAdjacent_mA{mass:02d}"
    c.cd(); c.SaveAs(str(outdir / f"{name}.png"))
    SaveCanvPic(c, str(outdir), name)


def draw_summary(results, outdir):
    ms = sorted(results)
    g = ROOT.TGraphErrors(len(ms))
    for i, m in enumerate(ms):
        r = results[m]
        g.SetPoint(i, m, r["f_ratio"]); g.SetPointError(i, 0, r["f_ratio_err"])
    c = ROOT.TCanvas("c_sum", "", 800, 600)
    c.SetLeftMargin(left_margin_for(800, 600) + 0.01); c.SetRightMargin(0.05); c.SetTopMargin(0.09); c.SetBottomMargin(0.15)
    c.SetTickx(1); c.SetTicky(1)
    ylo = min(0.0, min(results[m]["f_ratio"] - results[m]["f_ratio_err"] for m in ms))
    yhi = max(2.0, max(results[m]["f_ratio"] + results[m]["f_ratio_err"] for m in ms) * 1.15)
    fr = ROOT.TH1D("fr_sum", "", 1, 0, 31); fr.SetDirectory(0); fr.SetMinimum(ylo); fr.SetMaximum(yhi)
    fr.GetXaxis().SetTitle("m_{a} [GeV]"); fr.GetYaxis().SetTitle("f_{120-130}^{SR} / f_{120-130}^{adj.}")
    for ax in (fr.GetXaxis(), fr.GetYaxis()):
        ax.SetTitleSize(TITLE); ax.SetLabelSize(LABEL)
    fr.GetYaxis().SetTitleOffset(Y_TITLE_OFFSET); fr.GetXaxis().SetTitleOffset(X_TITLE_OFFSET)
    fr.GetXaxis().SetLabelOffset(0.010); fr.GetYaxis().SetNdivisions(505)
    fr.Draw("AXIS")
    one = ROOT.TLine(0, 1, 31, 1); one.SetLineStyle(2); one.SetLineColor(ROOT.kGray + 2); one.Draw()
    g.SetMarkerStyle(20); g.SetMarkerSize(1.4); g.SetMarkerColor(COL_SR); g.SetLineColor(COL_SR); g.SetLineWidth(2)
    g.Draw("PZ SAME")
    t = ROOT.TLatex(); t.SetNDC(); t.SetTextFont(42); t.SetTextSize(CMS_SIZE)
    t.SetTextAlign(11); t.DrawLatex(c.GetLeftMargin() + 0.005, 1 - 0.09 + 0.02, "#bf{CMS} #it{Preliminary}")
    t.SetTextAlign(31); t.DrawLatex(0.95, 1 - 0.09 + 0.02, LUMI_TEXT)
    t.SetTextAlign(13); t.SetTextSize(NOTE)
    t.DrawLatex(c.GetLeftMargin() + 0.03, 1 - 0.09 - 0.04,
                "Simulated DY background, SR slice [c, 1] vs adjacent [2c-1, c)")
    keep.extend([g, fr, one])
    c.SaveAs(str(outdir / "srAdjacent_peakFraction_summary.png"))
    SaveCanvPic(c, str(outdir), "srAdjacent_peakFraction_summary")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mva-cut-json", default=SC.DEFAULT_MVA_CUT_JSON)
    ap.add_argument("--output-dir", default=str(OUT_DEFAULT))
    ap.add_argument("--json", default=None, help="write per-mass numbers to this JSON")
    ap.add_argument("--masses", default=None, help="comma list (default: 1..30)")
    args = ap.parse_args()

    ROOT.gROOT.SetBatch(True); ROOT.gStyle.SetOptStat(0); ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetEndErrorSize(0); ROOT.gStyle.SetErrorX(0.)
    outdir = Path(args.output_dir); outdir.mkdir(parents=True, exist_ok=True)
    masses = [int(x) for x in args.masses.split(",")] if args.masses else list(SC.TARGET_MASSES)
    cuts = nearest_cuts(args.mva_cut_json)
    d = load_background(masses)
    m_all, w_all = d["H_m"], d["weight"]

    results = {}
    hdr = (f"{'mA':>3} {'anc':>3} {'cut':>6} {'adj.lo':>6} | {'nSR':>6} {'neffSR':>6} {'nAdj':>6} {'neffAdj':>7} | "
           f"{'chi2/ndf':>10} {'p':>5} {'KS D':>6} {'KS p':>5} | {'fSR':>13} {'fAdj':>13} {'ratio':>12} {'pull':>5}")
    print(hdr)
    for mass in masses:
        anchor, cut = cuts[mass]
        s = d[f"MVA_Score_mA_M{mass}"]
        sr = (s >= cut) & (s <= 1.0)
        adj = (s >= 2 * cut - 1) & (s < cut)
        h_sr = make_hist(f"hsr_{mass}", m_all[sr], w_all[sr])
        h_adj = make_hist(f"hadj_{mass}", m_all[adj], w_all[adj])
        chi2, ndf, p = chi2_shapes(h_sr, h_adj)
        D, ksp = weighted_ks(m_all[sr], w_all[sr], m_all[adj], w_all[adj])
        fs, fse = peak_fraction(m_all[sr], w_all[sr])
        fa, fae = peak_fraction(m_all[adj], w_all[adj])
        rat = fs / fa if fa > 0 else float("nan")
        rat_err = rat * math.sqrt((fse / fs) ** 2 + (fae / fa) ** 2) if fs > 0 and fa > 0 else float("nan")
        pull = (fs - fa) / math.sqrt(fse ** 2 + fae ** 2)
        res = dict(anchor=anchor, cut=cut, adj_low=2 * cut - 1,
                   n_sr=int(sr.sum()), neff_sr=neff(w_all[sr]), n_adj=int(adj.sum()), neff_adj=neff(w_all[adj]),
                   n_sr_peak=int((sr & (m_all > PEAK[0]) & (m_all < PEAK[1])).sum()),
                   n_adj_peak=int((adj & (m_all > PEAK[0]) & (m_all < PEAK[1])).sum()),
                   chi2=chi2, ndf=ndf, chi2_p=p, ks_D=D, ks_p=ksp,
                   f_sr=fs, f_sr_err=fse, f_adj=fa, f_adj_err=fae, f_ratio=rat, f_ratio_err=rat_err, f_pull=pull)
        results[mass] = res
        print(f"{mass:>3} {anchor:>3} {cut:>6.3f} {2*cut-1:>6.3f} | {res['n_sr']:>6} {res['neff_sr']:>6.0f} {res['n_adj']:>6} {res['neff_adj']:>7.0f} | "
              f"{chi2:>5.1f}/{ndf:<4} {p:>5.2f} {D:>6.3f} {ksp:>5.2f} | {fs:.3f}+-{fse:.3f} {fa:.3f}+-{fae:.3f} {rat:>5.2f}+-{rat_err:.2f} {pull:>5.1f}")
        draw_mass(mass, anchor, cut, h_sr, h_adj, res, outdir)

    if len(results) > 1:
        draw_summary(results, outdir)
    if args.json:
        with open(args.json, "w") as f:
            json.dump({str(k): v for k, v in results.items()}, f, indent=1)
        print("[json] wrote", args.json)
    print("[output]", outdir)


if __name__ == "__main__":
    main()
